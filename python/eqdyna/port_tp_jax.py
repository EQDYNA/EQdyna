"""
JAX rewrite of EQdyna's tpv1053d (friclaw=5, RSF slip law + strong rate
weakening + thermal pressurization) time loop -- phase 4. Elastic kernels
identical to port_jax.py/port_rsf_jax.py; `faulting`+`thermop` ported from
port_tp.py's NumPy reference (read/debugged first, per the standing
discipline) -- see port_tp.py's module docstring for the fric_tp_h=0.0
bug this port originally found and its f21afaf fix on master (now ported
as the per-node `fric(40)` array), the d9a50fa thetaPcTmp-fix provenance
check, and the friclaw==5-specific branches (theta_pc frozen, TPV==105
creeping-rate floor).

thermop's history integral is the real jit/scan structural question: the
Fortran's inner sum runs over j=1..nt-1 (a DYNAMIC-length slice, since nt
is the scan's own iteration index) -- not directly expressible as a static
-shape JAX op. Resolved by preallocating the FULL (nftnd, nsteps) history
arrays as loop-carry state (nsteps is static -- known before the loop
starts) and summing over ALL nsteps columns every step, relying on the
fact that not-yet-written "future" columns are exactly 0.0 (functional
`.at[].set` never touches them) so their contribution to the weighted sum
is exactly 0 regardless of what the (otherwise physically meaningless)
kernel evaluates to for them. The one hazard this creates: those future
columns' `age = (nt-(j+1))*dt` is negative, and `sqrt(negative)` is NaN,
and `0 * NaN = NaN` -- so age is floored at a small positive value
(`jnp.where(age > 0, age, dt)`) purely to keep the kernel finite; since the
zero history term already zeroes the product, this floor changes nothing
about the computed physics, only avoids NaN poisoning the sum. No masking
array needed, verified by full-run parity below.

float64 enabled first, before any other jax import.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from loading import load, region_damp  # noqa: E402
from port_jax import build, enable_compilation_cache, time_loop
import kernels_jax


def make_step(inv, S, nsteps):
    dt = inv['dt']; rdampk = inv['rdampk']
    nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    un = jnp.asarray(S['un']); us = jnp.asarray(S['us']); ud = jnp.asarray(S['ud']); arn = jnp.asarray(S['arn'])
    massSlave = jnp.asarray(S['fnms'][nsmp1]); massMaster = jnp.asarray(S['fnms'][nsmp2])
    mr = massMaster * massSlave / (massMaster + massSlave)
    T_coeff = arn * dt / mr
    xsource, ysource, zsource = S['xsource'], S['ysource'], S['zsource']
    nucR, nucT, nucdtau0, TPV = S['nucR'], S['nucT'], S['nucdtau0'], float(S['TPV'])
    coor_s = jnp.asarray(S['meshCoor'][nsmp1])
    radius = jnp.sqrt((coor_s[:, 0] - xsource) ** 2 + (coor_s[:, 1] - ysource) ** 2 +
                       (coor_s[:, 2] - zsource) ** 2)
    slipRateThres = S['slipRateThres']; C_elastic = S['C_elastic']
    idxF_s = [jnp.asarray(S['eq_ids'][nsmp1, d]) for d in range(3)]
    idxF_m = [jnp.asarray(S['eq_ids'][nsmp2, d]) for d in range(3)]

    fric0 = jnp.asarray(S['fric_init'])
    gama = fric0[:, 18] / fric0[:, 17]
    omega = fric0[:, 15]
    kapa = fric0[:, 16]
    rouc = fric0[:, 17]
    Tini = fric0[:, 40]
    fric_tp_h = fric0[:, 39]  # fric(40) = FRIC_SLOT_TP_H, per-node (post f21afaf fix)
    nftnd = S['nftnd']
    j_all = jnp.arange(nsteps)  # column index 0..nsteps-1 <-> Fortran history step j=1..nsteps

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

    def step(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist) = carry
        timeElapsed = timeElapsed + dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_jax.elastic_step(
            inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk)

        # ---- thermop: full-width sum, future columns are exactly 0 ----
        age = (nt - (j_all + 1)).astype(jnp.float64) * dt  # (nsteps,)
        age = jnp.where(age > 0, age, dt)
        denom_k = 4.0 * kapa[:, None] * age[None, :] + 2.0 * (fric_tp_h ** 2)[:, None]
        denom_o = 4.0 * omega[:, None] * age[None, :] + 2.0 * (fric_tp_h ** 2)[:, None]
        ker1 = (-kapa[:, None] / (omega - kapa)[:, None] / jnp.sqrt(denom_k)
                + omega[:, None] / (omega - kapa)[:, None] / jnp.sqrt(denom_o))
        hist_term = jnp.abs(shear_hist) * sliprate_hist  # (nftnd, nsteps); future cols are exactly 0
        patnode = (hist_term * ker1 * dt).sum(axis=1) * gama / jnp.sqrt(jnp.pi)
        ker2 = 1.0 / jnp.sqrt(denom_k)
        Tatnode = (hist_term * ker2 * dt).sum(axis=1) / rouc / jnp.sqrt(jnp.pi)
        fric = fric.at[:, 50].set(patnode)
        fric = fric.at[:, 51].set(Tatnode + Tini)

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

        F = jnp.where(radius < nucR, jnp.exp(radius ** 2 / (radius ** 2 - nucR ** 2)), 0.0)
        G = jnp.where(timeElapsed <= nucT,
                      jnp.exp((timeElapsed - nucT) ** 2 / (timeElapsed * (timeElapsed - 2.0 * nucT))), 1.0)
        dtau_nuc = (nucdtau0 * F * G) if (TPV == 104.0 or TPV == 105.0) else (0.0 * F)
        Ts = Ts + dtau_nuc

        # ---- solveRSF (friclaw==5) ----
        tnrm = Tn + fric[:, 50]
        tnrm = jnp.where(tnrm > 0.0, 0.0, tnrm)

        slipN = slipN0 + fric[:, 24] * timeElapsed
        srN = srN0 + fric[:, 24]; srS = srS0 + fric[:, 25]; srD = srD0 + fric[:, 26]
        srMag_sd = jnp.sqrt(srS ** 2 + srD ** 2)

        v_trial0 = srMag_sd
        a = fric[:, 8]; b = fric[:, 9]; Dc = fric[:, 10]; V0 = fric[:, 11]; f0 = fric[:, 12]
        fw = fric[:, 13]; Vw = fric[:, 14]
        state0 = fric[:, 19]
        xmu_pre, _, _ = rsf_slip_law(v_trial0, state0, a, b, Dc, f0, V0, fw, Vw)
        taoc_old = fric[:, 3] - xmu_pre * jnp.minimum(tnrm, 0.0)

        trialS = Ts - taoc_old * 0.5 * (srS / v_trial0) + fric[:, 25] / T_coeff
        trialD = Td - taoc_old * 0.5 * (srD / v_trial0) + fric[:, 26] / T_coeff
        trialMag = jnp.sqrt(trialS ** 2 + trialD ** 2)

        def newton_body(i, nc):
            v_trial, done, state_out, taoc_new = nc
            xmu, dxmudv, state_this = rsf_slip_law(v_trial, state0, a, b, Dc, f0, V0, fw, Vw)
            taoc_this = fric[:, 3] - xmu * jnp.minimum(tnrm, 0.0)
            rsfeq = v_trial + T_coeff * (taoc_this * 0.5 - trialMag)
            drsfeqdv = 1.0 + T_coeff * (-dxmudv * jnp.minimum(tnrm, 0.0)) * 0.5
            exit_now = (jnp.abs(rsfeq / drsfeqdv) < 1.0e-14 * jnp.abs(v_trial)) & \
                       (jnp.abs(rsfeq) < 1.0e-6 * jnp.abs(v_trial)) & (~done)
            done_after = done | exit_now
            newSliprate = v_trial - rsfeq / drsfeqdv
            v_next = jnp.where(newSliprate <= 0.0, v_trial / 2.0, newSliprate)
            v_trial_new = jnp.where(done_after, v_trial, v_next)
            return (v_trial_new, done_after, state_this, taoc_this)

        done0 = jnp.zeros_like(v_trial0, dtype=bool)
        nc0 = (v_trial0, done0, state0, jnp.zeros_like(v_trial0))
        v_trial, done, state_out, taoc_new = jax.lax.fori_loop(0, 20, newton_body, nc0)

        creep_rate = fric[:, 45]
        v_trial = jnp.where((TPV == 105.0) & (v_trial < creep_rate), creep_rate, v_trial)

        fric = fric.at[:, 19].set(state_out)  # fric(23)/thetaPc NOT written -- not evolved

        Ts = taoc_old * 0.5 * (srS / srMag_sd) + taoc_new * 0.5 * (trialS / trialMag)
        Td = taoc_old * 0.5 * (srD / srMag_sd) + taoc_new * 0.5 * (trialD / trialMag)

        fric = fric.at[:, 77].set(tnrm); fric = fric.at[:, 78].set(Ts); fric = fric.at[:, 79].set(Td)
        fric = fric.at[:, 46].set(v_trial)
        fric = fric.at[:, 47].set(jnp.sqrt(Ts ** 2 + Td ** 2))
        sliprate_hist = sliprate_hist.at[:, nt - 1].set(fric[:, 46])
        shear_hist = shear_hist.at[:, nt - 1].set(fric[:, 47])

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

        new_carry = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
                     sliprate_hist, shear_hist)
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
    stress_i = jnp.zeros((inv['Ei'], 6), dtype=jnp.float64)
    s_p = jnp.zeros((inv['Ep'], 15), dtype=jnp.float64)
    fric = jnp.asarray(S['fric_init'].copy())
    fnft = jnp.full(nftnd, 99999.0, dtype=jnp.float64)
    timeElapsed = jnp.asarray(0.0, dtype=jnp.float64)
    sliprate_hist = jnp.zeros((nftnd, nsteps), dtype=jnp.float64)
    shear_hist = jnp.zeros((nftnd, nsteps), dtype=jnp.float64)

    # ---- THE SAME LOOP THE OTHER TWO JAX PORTS RUN (2026-09-15) ----
    # This port used to own a SECOND time loop -- its own `jax.jit(lambda c:
    # jax.lax.scan(step, c, xs=jnp.arange(1, nsteps+1)))` with `inv` closed
    # over -- while port_jax/port_rsf_jax went through port_jax.time_loop with
    # `inv` promoted to jit ARGUMENTS by split_inv(). One loop per friclaw is
    # not the Fortran's shape (driver.f90:10 is ONE `do nt = 1, nstep`; the
    # friclaw branch is two lines INSIDE faulting.f90:17-18), and the cost of
    # the divergence was concrete: the constant-promotion fix landed on two of
    # three ports and left THIS case at 12.13 GB instead of 3.95, because there
    # were three promotion decisions instead of one.
    #
    # WHY IT CAN SHARE THE LOOP -- the thing previously believed to prevent it.
    # thermop's history integral is genuinely non-Markovian: the oracle,
    # src/updateThermalPressurization.f90:22-33, sums j = 1..nt-1 with a kernel
    # `1/sqrt(4*kapa*(nt-j)*dt + 2*h^2)` whose weight on EVERY past term
    # changes as nt advances. There is no recursion and no exponential to fold
    # into an accumulator, and the kernel decays only as (lag)^-1/2, so a
    # rolling window would drop terms that still contribute -- the full
    # (nftnd, nsteps) history is REQUIRED, and Fortran keeps exactly the same
    # thing (onFaultTPHist(2, nftnd, nstep, ntotft), a global of the same
    # shape). That was the assumption worth checking, and it holds.
    # But it makes nsteps a SHAPE, not a trip count: the history WIDTH must be
    # static, while the number of iterations need not be. `time_loop` supplies
    # the 1-based step number to `step(carry, nt)` exactly as `lax.scan(...,
    # xs=jnp.arange(1, nsteps+1))` did, so the only thing this port gives up by
    # using it is a cache-key benefit it never had: the carry SHAPE still
    # depends on nsteps, so a new step count still recompiles here (unlike
    # port_jax/port_rsf_jax, where time_loop's traced trip count means one
    # executable serves every step count).
    #
    # NO SILENT CLAMP (rule 2): `sliprate_hist.at[:, nt-1].set()` would CLAMP,
    # not raise, if the trip count ever exceeded the history width -- silently
    # writing every overflowing step into the last column. The two come from
    # the same `nsteps` three lines above, so they cannot disagree; asserted
    # rather than assumed, because the failure mode is invisible.
    #
    # THIS CHANGES THE ANSWER. Read the numbers before keeping it.
    #
    # MEASURED, test.tpv1053d serial (E=912600, nftnd=4005), 120 steps,
    # jax-cpu, full-run peak RSS, four configurations, each run fresh:
    #   (0) own lax.scan, inv closed over        12.130 GB  69.1 s  575.9 ms/st
    #   (1) (0) + this session's kernels_jax
    #       block scatter and int32 indices      10.084 GB  86.5 s  720.5 ms/st
    #   (2) (1) + lax.scan -> fori_loop only,
    #       inv STILL closed over                10.536 GB  89.9 s  748.7 ms/st
    #   (3) (1) + shared time_loop (fori_loop
    #       AND split_inv promotion) = THIS       3.952 GB  59.3 s  494.4 ms/st
    # (3) against (0): -67% peak RSS, -14% ms/step. (1) is BIT-IDENTICAL to
    # (0) -- all five digests unchanged. (2) and (3) are NOT:
    #     (1) -> (2)   velArr 1.2279e-12  dispArr 9.2944e-14
    #                  force  5.3512e-11  fric    5.9001e-06   fnft 0 flips
    #     (2) -> (3)   velArr 6.6822e-13  dispArr 4.9703e-14
    #                  force  2.5344e-11  fric    2.4038e-06   fnft 0 flips
    #     (0) -> (3)   velArr 1.5985e-12  dispArr 1.4260e-13
    #                  force  7.8856e-11  fric    7.9796e-06   fnft 0 flips
    # (max abs over the FULL arrays, sha256-over-tobytes comparison at full
    # step count; fric's columns carry ~1e8 Pa stresses, so 8e-06 absolute is
    # ~1e-13 relative, and no rupture time flips). Against the committed
    # reference the accept tier measures this case at 4.245508e-06 with (3) in
    # place, versus its 1e-4 bound and the 3.73e-06 recorded with (0) -- i.e.
    # the port moved 23x less than the bound it is gated on, and the tier is
    # 5/5 green (re-run fresh this session, not inherited).
    #
    # ATTRIBUTION -- and a correction to what this file used to say. The old
    # note recorded fric moving 5.9001e-06 (with velArr 1.2279e-12, dispArr
    # 9.2944e-14, force 5.3512e-11) and blamed PROMOTION, reasoning that it
    # turns the element-force scatter's index into one contiguous buffer
    # instead of an in-graph concatenate. Those four numbers are, to every
    # digit, the (1) -> (2) row above -- which has NO promotion in it at all
    # and, after this session's kernels_jax change, no concatenated scatter
    # index either. The deviation is the LOOP CONSTRUCT: swapping lax.scan for
    # lax.fori_loop changes the HLO the step body sits in, and this program has
    # a reduction that is sensitive to it (thermop's `(hist_term*ker*dt).sum
    # (axis=1)` over the 120-column history, whose result feeds a 20-iteration
    # Newton solve that amplifies the last bits). Promotion adds a further
    # 2.4e-06 of the same kind. Neither is a scatter-ordering effect.
    #
    # WHY IT IS KEPT ANYWAY, stated plainly rather than buried: (0) is the only
    # bit-identical option and it costs 12.13 GB -- 3x the footprint of every
    # other case in the suite and the largest single number in the project's
    # memory table. The trade taken here is 8e-06 absolute on fric (1e-13
    # relative, no rupture-time flip, 12x inside this case's own accept bound)
    # for -67% peak RSS and -14% ms/step, plus the structural fact that all
    # three JAX ports now run the SAME loop, so the next loop-level fix cannot
    # land on two of three again. To go back to bit-identical, restore the two
    # lines below to `jax.jit(lambda c: jax.lax.scan(make_step(inv, S, nsteps),
    # c, xs=jnp.arange(1, nsteps + 1))[0])(carry0)` -- nothing else in this
    # file depends on the choice.
    if sliprate_hist.shape[1] != nsteps or shear_hist.shape[1] != nsteps:
        raise ValueError(
            'port_tp_jax.run: thermop history width %d/%d must equal the step count %d -- '
            'a shorter history would make `.at[:, nt-1].set()` clamp silently and pile '
            'every overflowing step into the last column'
            % (sliprate_hist.shape[1], shear_hist.shape[1], nsteps))
    carry0 = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
              sliprate_hist, shear_hist)
    carry = time_loop(lambda i: make_step(i, S, nsteps), inv, carry0, nsteps)
    jax.block_until_ready(carry)
    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist) = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr), fnft=np.asarray(fnft),
                fric=np.asarray(fric), force=np.asarray(force))
