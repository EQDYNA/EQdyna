"""
JAX rewrite of EQdyna's tpv104 (friclaw=4, rate-and-state with strong rate
weakening) time loop -- phase 3 of the feasibility spike. Elastic kernels
identical to port_jax.py (tpv8); `faulting` replaced with the RSF path
(port_rsf.py's NumPy reference, ported statement-for-statement), with the
Newton-Raphson solve expressed as a fixed `jax.lax.fori_loop` (20 iterations)
over a per-node "converged" mask -- see port_rsf.py's module docstring for
the proof that masking-with-idempotent-recompute reproduces the Fortran
`do`-loop-with-`exit` exactly, including the "never converged in 20 tries"
fallback (that node's v_trial DOES get the 20th update applied, matching
Fortran's `do` loop running its full body on iv==ivmax).

float64 enabled first, before any other jax import.

Milestone 10 (drv.a6, TPV==2802, C_elastic==0): faulting.f90:363-374's
rsfNucleation TPV==2802 branch mirrors port_rsf.py's nt==1 special case
(see that module's docstring/step function for the full derivation) --
detected here via `timeElapsed == dt` since this scan's carry has no
explicit step counter.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from loading import load, region_damp  # noqa: E402
from port_jax import build, enable_compilation_cache, time_loop  # elastic-kernel invariants: identical for any friclaw
import kernels_jax


def make_step(inv, S):
    dt = inv['dt']; rdampk = inv['rdampk']
    nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    un = jnp.asarray(S['un']); us = jnp.asarray(S['us']); ud = jnp.asarray(S['ud']); arn = jnp.asarray(S['arn'])
    massSlave = jnp.asarray(S['fnms'][nsmp1]); massMaster = jnp.asarray(S['fnms'][nsmp2])
    mr = massMaster * massSlave / (massMaster + massSlave)
    T_coeff_const = arn * dt / mr
    xsource, ysource, zsource = S['xsource'], S['ysource'], S['zsource']
    nucR, nucT, nucdtau0, TPV = S['nucR'], S['nucT'], S['nucdtau0'], float(S['TPV'])
    coor_s = jnp.asarray(S['meshCoor'][nsmp1])
    radius = jnp.sqrt((coor_s[:, 0] - xsource) ** 2 + (coor_s[:, 1] - ysource) ** 2 +
                       (coor_s[:, 2] - zsource) ** 2)
    slipRateThres = S['slipRateThres']; C_elastic = S['C_elastic']
    # faulting.f90 solveRSF min_norm/max_norm clamp (see port_rsf.py's
    # comment) -- insertFaultType/C_elastic are static (Python int) at
    # trace time, so this `if` is resolved once, not per-node.
    insertFaultType = S.get('insertFaultType', 0)
    min_norm = -10.0e6
    max_norm = -40.0e6
    idxF_s = [jnp.asarray(S['eq_ids'][nsmp1, d]) for d in range(3)]
    idxF_m = [jnp.asarray(S['eq_ids'][nsmp2, d]) for d in range(3)]

    def rsf_slip_law(V2, psi, a, b, Dc, f0, V0, fw, Vw):
        tmpc = 1.0 / (2.0 * V0) * jnp.exp(psi / a)
        tmp = (V2 + 1.0e-30) * tmpc
        xmu = a * jnp.log(tmp + jnp.sqrt(tmp ** 2 + 1.0))
        dxmudv = a * tmpc / jnp.sqrt(1.0 + tmp ** 2)
        fLV = f0 - (b - a) * jnp.log(V2 / V0)
        fss = fw + (fLV - fw) / ((1.0 + (V2 / Vw) ** 8) ** 0.125)
        fssa = fss / a
        psiss = a * jnp.log(2.0 * V0 / V2 * (jnp.exp(fssa) - jnp.exp(-fssa)) / 2.0)
        psi_new = psiss + (psi - psiss) * jnp.exp(-V2 * dt / Dc)
        return xmu, dxmudv, psi_new

    def rsf_normal_stress(V2, theta_pc, tnrm, L):
        theta_pc_dot = -V2 / L * (theta_pc - jnp.abs(tnrm))
        return theta_pc + theta_pc_dot * dt

    def step(carry, _):
        v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
        timeElapsed = timeElapsed + dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_jax.elastic_step(
            inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk)

        # ---- faulting: getNsdSlipSliprateTraction (shared prelude) ----
        fx_s = force[idxF_s[0]]; fy_s = force[idxF_s[1]]; fz_s = force[idxF_s[2]]
        fx_m = force[idxF_m[0]]; fy_m = force[idxF_m[1]]; fz_m = force[idxF_m[2]]
        vx_s, vy_s, vz_s = velArr[nsmp1, 0], velArr[nsmp1, 1], velArr[nsmp1, 2]
        vx_m, vy_m, vz_m = velArr[nsmp2, 0], velArr[nsmp2, 1], velArr[nsmp2, 2]
        dx_s, dy_s, dz_s = dispArr[nsmp1, 0], dispArr[nsmp1, 1], dispArr[nsmp1, 2]
        dx_m, dy_m, dz_m = dispArr[nsmp2, 0], dispArr[nsmp2, 1], dispArr[nsmp2, 2]

        def rot(vx, vy, vz, dirvec):
            return vx * dirvec[:, 0] + vy * dirvec[:, 1] + vz * dirvec[:, 2]

        fN_s, fS_s, fD_s = rot(fx_s, fy_s, fz_s, un), rot(fx_s, fy_s, fz_s, us), rot(fx_s, fy_s, fz_s, ud)
        fN_m, fS_m, fD_m = rot(fx_m, fy_m, fz_m, un), rot(fx_m, fy_m, fz_m, us), rot(fx_m, fy_m, fz_m, ud)
        vN_s, vS_s, vD_s = rot(vx_s, vy_s, vz_s, un), rot(vx_s, vy_s, vz_s, us), rot(vx_s, vy_s, vz_s, ud)
        vN_m, vS_m, vD_m = rot(vx_m, vy_m, vz_m, un), rot(vx_m, vy_m, vz_m, us), rot(vx_m, vy_m, vz_m, ud)
        dN_s, dS_s, dD_s = rot(dx_s, dy_s, dz_s, un), rot(dx_s, dy_s, dz_s, us), rot(dx_s, dy_s, dz_s, ud)
        dN_m, dS_m, dD_m = rot(dx_m, dy_m, dz_m, un), rot(dx_m, dy_m, dz_m, us), rot(dx_m, dy_m, dz_m, ud)

        slipN0 = dN_m - dN_s; slipS0 = dS_m - dS_s; slipD0 = dD_m - dD_s
        srN0 = vN_m - vN_s; srS0 = vS_m - vS_s; srD0 = vD_m - vD_s
        srMag_full = jnp.sqrt(srN0 ** 2 + srS0 ** 2 + srD0 ** 2)

        fric = fric.at[:, 70].set(slipS0); fric = fric.at[:, 71].set(slipD0); fric = fric.at[:, 72].set(slipN0)
        fric = fric.at[:, 73].set(srS0); fric = fric.at[:, 74].set(srD0)
        fric = fric.at[:, 75].set(jnp.maximum(fric[:, 75], srMag_full))
        fric = fric.at[:, 76].add(srMag_full * dt)

        totalMass = (massSlave + massMaster) * arn
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic

        # ---- rsfNucleation ----
        F = jnp.where(radius < nucR, jnp.exp(radius ** 2 / (radius ** 2 - nucR ** 2)), 0.0)
        G = jnp.where(timeElapsed <= nucT,
                      jnp.exp((timeElapsed - nucT) ** 2 / (timeElapsed * (timeElapsed - 2.0 * nucT))), 1.0)
        if TPV == 104.0 or TPV == 105.0:
            dtau_nuc = nucdtau0 * F * G
        elif TPV == 2802.0:
            # faulting.f90:363-374 (drv.a6, C_elastic==0/plastic) -- mirrors
            # port_rsf.py's nt==1 special case exactly (see that module's
            # docstring for the derivation/why-it-matters), except `nt==1`
            # is detected via `timeElapsed == dt` (a traced comparison,
            # since there is no explicit step counter in this scan's carry)
            # -- `0.0 + dt == dt` exactly in IEEE754, so this is bit-exact,
            # not an approximation. `TPV==2802.0` is resolved at TRACE time
            # (a static Python float), so this whole branch (including the
            # is_step1 jnp.where machinery) is compiled in ONLY for this
            # case -- zero-cost/absent for TPV==104/105 (port_rsf_jax.py's
            # only other caller).
            is_step1 = (timeElapsed == dt)
            ttao = jnp.sqrt(Ts ** 2 + Td ** 2)
            backSliprate = jnp.sqrt((srS0 + fric[:, 25]) ** 2 + (srD0 + fric[:, 26]) ** 2)
            rsf_a0 = fric[:, 8]; rsf_v0_0 = fric[:, 11]
            state_step1 = rsf_a0 * jnp.log(2.0 * rsf_v0_0 / backSliprate
                                             * jnp.sinh(ttao / jnp.abs(Tn) / rsf_a0))
            thetapc_step1 = jnp.abs(Tn)
            fric = fric.at[:, 19].set(jnp.where(is_step1, state_step1, fric[:, 19]))
            fric = fric.at[:, 22].set(jnp.where(is_step1, thetapc_step1, fric[:, 22]))
            fric = fric.at[:, 80].set(jnp.where(is_step1, nucdtau0, fric[:, 80]))
            dtau_nuc = fric[:, 80] * F * G
        else:
            dtau_nuc = 0.0 * F
        Ts = Ts + dtau_nuc

        # ---- solveRSF ----
        tnrm = Tn + fric[:, 5]
        if insertFaultType > 0 and C_elastic == 1:
            tnrm = jnp.where(tnrm >= min_norm, min_norm,
                              jnp.where(tnrm <= max_norm, max_norm, tnrm))
        tnrm = jnp.where(tnrm > 0.0, 0.0, tnrm)

        slipN = slipN0 + fric[:, 24] * timeElapsed
        srN = srN0 + fric[:, 24]; srS = srS0 + fric[:, 25]; srD = srD0 + fric[:, 26]
        srMag_sd = jnp.sqrt(srS ** 2 + srD ** 2)

        v_trial0 = srMag_sd
        a = fric[:, 8]; b = fric[:, 9]; Dc = fric[:, 10]; V0 = fric[:, 11]; f0 = fric[:, 12]
        fw = fric[:, 13]; Vw = fric[:, 14]
        theta_pc_tmp = fric[:, 22]
        state0 = fric[:, 19]
        xmu_pre, _, _ = rsf_slip_law(v_trial0, state0, a, b, Dc, f0, V0, fw, Vw)
        taoc_old = xmu_pre * theta_pc_tmp

        T_coeff = T_coeff_const
        trialS = Ts - taoc_old * 0.5 * (srS / v_trial0) + fric[:, 25] / T_coeff
        trialD = Td - taoc_old * 0.5 * (srD / v_trial0) + fric[:, 26] / T_coeff
        trialMag = jnp.sqrt(trialS ** 2 + trialD ** 2)

        def newton_body(i, nc):
            v_trial, done, state_out, thetaPc_out, taoc_new = nc
            xmu, dxmudv, state_this = rsf_slip_law(v_trial, state0, a, b, Dc, f0, V0, fw, Vw)
            thetaPc_this = rsf_normal_stress(v_trial, theta_pc_tmp, tnrm, Dc)
            taoc_this = xmu * thetaPc_this
            rsfeq = v_trial + T_coeff * (taoc_this * 0.5 - trialMag)
            drsfeqdv = 1.0 + T_coeff * (dxmudv * thetaPc_this) * 0.5
            exit_now = (jnp.abs(rsfeq / drsfeqdv) < 1.0e-14 * jnp.abs(v_trial)) & \
                       (jnp.abs(rsfeq) < 1.0e-6 * jnp.abs(v_trial)) & (~done)
            done_after = done | exit_now
            newSliprate = v_trial - rsfeq / drsfeqdv
            v_next = jnp.where(newSliprate <= 0.0, v_trial / 2.0, newSliprate)
            v_trial_new = jnp.where(done_after, v_trial, v_next)
            return (v_trial_new, done_after, state_this, thetaPc_this, taoc_this)

        done0 = jnp.zeros_like(v_trial0, dtype=bool)
        nc0 = (v_trial0, done0, state0, theta_pc_tmp, jnp.zeros_like(v_trial0))
        v_trial, done, state_out, thetaPc_out, taoc_new = jax.lax.fori_loop(0, 20, newton_body, nc0)

        fric = fric.at[:, 19].set(state_out); fric = fric.at[:, 22].set(thetaPc_out)

        Ts = taoc_old * 0.5 * (srS / srMag_sd) + taoc_new * 0.5 * (trialS / trialMag)
        Td = taoc_old * 0.5 * (srD / srMag_sd) + taoc_new * 0.5 * (trialD / trialMag)

        fric = fric.at[:, 77].set(tnrm); fric = fric.at[:, 78].set(Ts); fric = fric.at[:, 79].set(Td)
        fric = fric.at[:, 46].set(v_trial)
        fric = fric.at[:, 47].set(jnp.sqrt(Ts ** 2 + Td ** 2))

        accN = -srN / dt - slipN / dt / dt
        accS = (v_trial * (trialS / trialMag) - srS) / dt
        accD = (v_trial * (trialD / trialMag) - srD) / dt
        xAcc = accN * un[:, 0] + accS * us[:, 0] + accD * ud[:, 0]
        yAcc = accN * un[:, 1] + accS * us[:, 1] + accD * ud[:, 1]
        zAcc = accN * un[:, 2] + accS * us[:, 2] + accD * ud[:, 2]
        xR = fx_s + fx_m; yR = fy_s + fy_m; zR = fz_s + fz_m

        force = force.at[idxF_s[0]].set((-xAcc + xR / massMaster) * mr)
        force = force.at[idxF_s[1]].set((-yAcc + yR / massMaster) * mr)
        force = force.at[idxF_s[2]].set((-zAcc + zR / massMaster) * mr)
        force = force.at[idxF_m[0]].set((xAcc + xR / massSlave) * mr)
        force = force.at[idxF_m[1]].set((yAcc + yR / massSlave) * mr)
        force = force.at[idxF_m[2]].set((zAcc + zR / massSlave) * mr)

        fric = fric.at[:, 30].set(vx_m + (xAcc + xR / massSlave) * dt)
        fric = fric.at[:, 31].set(vy_m + (yAcc + yR / massSlave) * dt)
        fric = fric.at[:, 32].set(vz_m + (zAcc + zR / massSlave) * dt)
        fric = fric.at[:, 33].set(vx_s + (-xAcc + xR / massMaster) * dt)
        fric = fric.at[:, 34].set(vy_s + (-yAcc + yR / massMaster) * dt)
        fric = fric.at[:, 35].set(vz_s + (-zAcc + zR / massMaster) * dt)

        force = force.at[0].set(0.0)

        need = fnft > 5000.0
        fnft = jnp.where(need & (srMag_sd >= slipRateThres), timeElapsed, fnft)

        force = force * inv['inv_mass_full']

        new_carry = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
        return new_carry, None

    return step


def run(S, nsteps=None, verbose=True):
    nsteps = nsteps or S['nstep']
    enable_compilation_cache()
    inv = build(S)
    conn = S['conn']; elemType = S['elemType']
    E_pml = np.nonzero(elemType == 2)[0]
    inv['conn_p'] = jnp.asarray(conn[E_pml])

    N = S['N']; NEQ = S['NEQ']; nftnd = S['nftnd']
    v1 = jnp.zeros(NEQ + 1, dtype=jnp.float64)
    velArr = jnp.zeros((N, 3), dtype=jnp.float64)
    dispArr = jnp.zeros((N, 3), dtype=jnp.float64)
    force = jnp.zeros(NEQ + 1, dtype=jnp.float64)
    stress_i = inv['stress_i0']  # zeros for C_elastic==1; lithostatic pre-stress otherwise
    # (Milestone 10, drv.a6 -- see port_jax.py's build(), shared by this module).
    s_p = jnp.zeros((inv['Ep'], 15), dtype=jnp.float64)
    fric = jnp.asarray(S['fric_init'].copy())
    fnft = jnp.full(nftnd, 99999.0, dtype=jnp.float64)
    timeElapsed = jnp.asarray(0.0, dtype=jnp.float64)

    carry0 = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
    # `inv` goes in as a jit ARGUMENT, not a closure constant -- hence the
    # factory (lambda i: make_step(i, S)) rather than a pre-built step; see
    # port_jax.split_inv() for the measured peak-RSS reason. nsteps traced.
    carry = time_loop(lambda i: make_step(i, S), inv, carry0, nsteps)
    jax.block_until_ready(carry)
    v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr), fnft=np.asarray(fnft),
                fric=np.asarray(fric), force=np.asarray(force))
