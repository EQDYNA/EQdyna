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
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from port import load, region_damp  # noqa: E402
from port_jax import build  # elastic-kernel invariants: identical for any friclaw


def _c(dN, v):
    return (dN * v).sum(axis=1)


def make_step(inv, S):
    dt = inv['dt']; rdampk = inv['rdampk']; NEQ1 = inv['NEQ1']
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

        # ---- velDispUpdate (identical to port_jax.py) ----
        idx3_v = inv['idx3_v']
        accel3 = force[idx3_v]
        v1 = v1.at[idx3_v].set(v1[idx3_v] + accel3 * dt)
        velArr = velArr.at[inv['int_nodes_idx']].set(v1[idx3_v])
        dispArr = dispArr.at[inv['int_nodes_idx']].add(v1[idx3_v] * dt)

        idx12_v = inv['idx12_v']
        f9 = force[idx12_v[:, 0:9]]
        vold9 = v1[idx12_v[:, 0:9]]
        v1 = v1.at[idx12_v[:, 0:9]].set((f9 + vold9 * inv['a9']) / inv['b9'])
        f3v9 = force[idx12_v[:, 9:12]]
        v1 = v1.at[idx12_v[:, 9:12]].set(v1[idx12_v[:, 9:12]] + f3v9 * dt)
        v1 = v1.at[0].set(0.0)
        vA = v1[idx12_v[:, 0]] + v1[idx12_v[:, 1]] + v1[idx12_v[:, 2]] + v1[idx12_v[:, 9]]
        vB = v1[idx12_v[:, 3]] + v1[idx12_v[:, 4]] + v1[idx12_v[:, 5]] + v1[idx12_v[:, 10]]
        vC = v1[idx12_v[:, 6]] + v1[idx12_v[:, 7]] + v1[idx12_v[:, 8]] + v1[idx12_v[:, 11]]
        has_eq = idx12_v[:, 0] > 0
        pnodes = inv['pml_nodes_idx']
        velArr = velArr.at[pnodes, 0].set(jnp.where(has_eq, vA, 0.0))
        velArr = velArr.at[pnodes, 1].set(jnp.where(has_eq, vB, 0.0))
        velArr = velArr.at[pnodes, 2].set(jnp.where(has_eq, vC, 0.0))
        dispArr = dispArr.at[pnodes, 0].set(jnp.where(has_eq, dispArr[pnodes, 0] + vA * dt, 0.0))
        dispArr = dispArr.at[pnodes, 1].set(jnp.where(has_eq, dispArr[pnodes, 1] + vB * dt, 0.0))
        dispArr = dispArr.at[pnodes, 2].set(jnp.where(has_eq, dispArr[pnodes, 2] + vC * dt, 0.0))

        scat_idx = []; scat_val = []

        # ---- assembleGlobalKU interior ----
        conn_i = inv['conn_i']
        vl = velArr[conn_i]
        dNx_i, dNy_i, dNz_i = inv['dNx_i'], inv['dNy_i'], inv['dNz_i']
        lam_i, miu_i, constk_w_i = inv['lam_i'], inv['miu_i'], inv['constk_w_i']
        sr1 = _c(dNx_i, vl[:, :, 0]); sr2 = _c(dNy_i, vl[:, :, 1]); sr3 = _c(dNz_i, vl[:, :, 2])
        sr4 = _c(dNz_i, vl[:, :, 1]) + _c(dNy_i, vl[:, :, 2])
        sr5 = _c(dNz_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 2])
        sr6 = _c(dNy_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 1])
        vol_sr = sr1 + sr2 + sr3
        strr1 = lam_i * vol_sr + 2 * miu_i * sr1
        strr2 = lam_i * vol_sr + 2 * miu_i * sr2
        strr3 = lam_i * vol_sr + 2 * miu_i * sr3
        strr4 = miu_i * sr4; strr5 = miu_i * sr5; strr6 = miu_i * sr6
        stress_i = stress_i.at[:, 0].add(strr1 * dt); stress_i = stress_i.at[:, 1].add(strr2 * dt)
        stress_i = stress_i.at[:, 2].add(strr3 * dt); stress_i = stress_i.at[:, 3].add(strr4 * dt)
        stress_i = stress_i.at[:, 4].add(strr5 * dt); stress_i = stress_i.at[:, 5].add(strr6 * dt)
        st1 = constk_w_i * (stress_i[:, 0] + rdampk * strr1)
        st2 = constk_w_i * (stress_i[:, 1] + rdampk * strr2)
        st3 = constk_w_i * (stress_i[:, 2] + rdampk * strr3)
        st4 = constk_w_i * (stress_i[:, 3] + rdampk * strr4)
        st5 = constk_w_i * (stress_i[:, 4] + rdampk * strr5)
        st6 = constk_w_i * (stress_i[:, 5] + rdampk * strr6)
        Fx = dNx_i * st1[:, None] + dNz_i * st5[:, None] + dNy_i * st6[:, None]
        Fy = dNy_i * st2[:, None] + dNz_i * st4[:, None] + dNx_i * st6[:, None]
        Fz = dNz_i * st3[:, None] + dNy_i * st4[:, None] + dNx_i * st5[:, None]
        scat_idx += [inv['idxIx'], inv['idxIy'], inv['idxIz']]
        scat_val += [Fx.ravel(), Fy.ravel(), Fz.ravel()]

        # ---- assembleGlobalKU PML ----
        if inv['Ep'] > 0:
            conn_p = inv['conn_p']
            vlp = velArr[conn_p]
            dNx_p, dNy_p, dNz_p = inv['dNx_p'], inv['dNy_p'], inv['dNz_p']
            lam_p, miu_p, det_w_p = inv['lam_p'], inv['miu_p'], inv['det_w_p']
            srp1 = _c(dNx_p, vlp[:, :, 0]); srp2 = _c(dNy_p, vlp[:, :, 1]); srp3 = _c(dNz_p, vlp[:, :, 2])
            srp4 = _c(dNz_p, vlp[:, :, 1]) + _c(dNy_p, vlp[:, :, 2])
            srp5 = _c(dNz_p, vlp[:, :, 0]) + _c(dNx_p, vlp[:, :, 2])
            srp6 = _c(dNy_p, vlp[:, :, 0]) + _c(dNx_p, vlp[:, :, 1])
            volp = srp1 + srp2 + srp3
            srate = [lam_p * volp + 2 * miu_p * srp1, lam_p * volp + 2 * miu_p * srp2,
                     lam_p * volp + 2 * miu_p * srp3, miu_p * srp4, miu_p * srp5, miu_p * srp6]
            Dxvx = _c(dNx_p, vlp[:, :, 0]); Dyvy = _c(dNy_p, vlp[:, :, 1]); Dzvz = _c(dNz_p, vlp[:, :, 2])
            Dxvy = _c(dNx_p, vlp[:, :, 1]); Dyvx = _c(dNy_p, vlp[:, :, 0])
            Dxvz = _c(dNx_p, vlp[:, :, 2]); Dzvx = _c(dNz_p, vlp[:, :, 0])
            Dyvz = _c(dNy_p, vlp[:, :, 2]); Dzvy = _c(dNz_p, vlp[:, :, 1])
            a1, b1, a2, b2, a3, b3 = inv['a1'], inv['b1'], inv['a2'], inv['b2'], inv['a3'], inv['b3']
            sp = s_p
            sp = sp.at[:, 0].set(((lam_p + 2 * miu_p) * Dxvx + a1 * sp[:, 0]) / b1)
            sp = sp.at[:, 1].set((lam_p * Dyvy + a2 * sp[:, 1]) / b2)
            sp = sp.at[:, 2].set((lam_p * Dzvz + a3 * sp[:, 2]) / b3)
            sp = sp.at[:, 3].set((lam_p * Dxvx + a1 * sp[:, 3]) / b1)
            sp = sp.at[:, 4].set(((lam_p + 2 * miu_p) * Dyvy + a2 * sp[:, 4]) / b2)
            sp = sp.at[:, 5].set((lam_p * Dzvz + a3 * sp[:, 5]) / b3)
            sp = sp.at[:, 6].set((lam_p * Dxvx + a1 * sp[:, 6]) / b1)
            sp = sp.at[:, 7].set((lam_p * Dyvy + a2 * sp[:, 7]) / b2)
            sp = sp.at[:, 8].set(((lam_p + 2 * miu_p) * Dzvz + a3 * sp[:, 8]) / b3)
            sp = sp.at[:, 9].set((miu_p * Dxvy + a1 * sp[:, 9]) / b1)
            sp = sp.at[:, 10].set((miu_p * Dyvx + a2 * sp[:, 10]) / b2)
            sp = sp.at[:, 11].set((miu_p * Dxvz + a1 * sp[:, 11]) / b1)
            sp = sp.at[:, 12].set((miu_p * Dzvx + a3 * sp[:, 12]) / b3)
            sp = sp.at[:, 13].set((miu_p * Dyvz + a2 * sp[:, 13]) / b2)
            sp = sp.at[:, 14].set((miu_p * Dzvy + a3 * sp[:, 14]) / b3)
            s_p = sp
            sxx = s_p[:, 0] + s_p[:, 1] + s_p[:, 2]
            syy = s_p[:, 3] + s_p[:, 4] + s_p[:, 5]
            szz = s_p[:, 6] + s_p[:, 7] + s_p[:, 8]
            sxy = s_p[:, 9] + s_p[:, 10]; sxz = s_p[:, 11] + s_p[:, 12]; syz = s_p[:, 13] + s_p[:, 14]
            s0 = [rdampk * srate[k] for k in range(6)]
            f1 = -det_w_p[:, None] * dNx_p * sxx[:, None]
            f2 = -det_w_p[:, None] * dNy_p * sxy[:, None]
            f3v = -det_w_p[:, None] * dNz_p * sxz[:, None]
            f4 = -det_w_p[:, None] * dNx_p * sxy[:, None]
            f5 = -det_w_p[:, None] * dNy_p * syy[:, None]
            f6 = -det_w_p[:, None] * dNz_p * syz[:, None]
            f7 = -det_w_p[:, None] * dNx_p * sxz[:, None]
            f8 = -det_w_p[:, None] * dNy_p * syz[:, None]
            f9v = -det_w_p[:, None] * dNz_p * szz[:, None]
            f10 = -det_w_p[:, None] * (dNx_p * s0[0][:, None] + dNz_p * s0[4][:, None] + dNy_p * s0[5][:, None])
            f11 = -det_w_p[:, None] * (dNy_p * s0[1][:, None] + dNz_p * s0[3][:, None] + dNx_p * s0[5][:, None])
            f12 = -det_w_p[:, None] * (dNz_p * s0[2][:, None] + dNy_p * s0[3][:, None] + dNx_p * s0[4][:, None])
            efPML12 = [f1, f2, f3v, f4, f5, f6, f7, f8, f9v, f10, f11, f12]
            for j in range(12):
                scat_idx.append(inv['idxP12'][j]); scat_val.append(efPML12[j].ravel())
            grp_sums = [efPML12[0] + efPML12[1] + efPML12[2] + efPML12[9],
                        efPML12[3] + efPML12[4] + efPML12[5] + efPML12[10],
                        efPML12[6] + efPML12[7] + efPML12[8] + efPML12[11]]
            for g in range(3):
                scat_idx.append(inv['idxP3'][g]); scat_val.append(grp_sums[g].ravel())

        # ---- hrglss ----
        conn = inv['conn']; phi = inv['phi']; ss = inv['ss']
        dl_all = dispArr[conn] + rdampk * velArr[conn]
        for m in range(4):
            phid = (phi[:, m, :, None] * dl_all).sum(axis=1)
            r0 = ss[:, 0] * phid[:, 0] + ss[:, 1] * phid[:, 1] + ss[:, 2] * phid[:, 2]
            r1 = ss[:, 1] * phid[:, 0] + ss[:, 3] * phid[:, 1] + ss[:, 4] * phid[:, 2]
            r2 = ss[:, 2] * phid[:, 0] + ss[:, 4] * phid[:, 1] + ss[:, 5] * phid[:, 2]
            fh0 = -(phi[:, m, :] * r0[:, None]).ravel()
            fh1 = -(phi[:, m, :] * r1[:, None]).ravel()
            fh2 = -(phi[:, m, :] * r2[:, None]).ravel()
            scat_idx += [inv['idxH0'], inv['idxH1'], inv['idxH2']]
            scat_val += [fh0, fh1, fh2]

        all_idx = jnp.concatenate(scat_idx)
        all_val = jnp.concatenate(scat_val)
        force = jnp.zeros(NEQ1, dtype=jnp.float64).at[all_idx].add(all_val)
        force = force.at[0].set(0.0)

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
        dtau_nuc = (nucdtau0 * F * G) if (TPV == 104.0 or TPV == 105.0) else (0.0 * F)
        Ts = Ts + dtau_nuc

        # ---- solveRSF ----
        tnrm = Tn + fric[:, 5]
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

    step = make_step(inv, S)
    scan_fn = jax.jit(lambda c: jax.lax.scan(step, c, xs=None, length=nsteps)[0])

    carry0 = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
    carry = scan_fn(carry0)
    jax.block_until_ready(carry)
    v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr), fnft=np.asarray(fnft),
                fric=np.asarray(fric), force=np.asarray(force))
