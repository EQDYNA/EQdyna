"""
NumPy reference port of EQdyna's tpv104 (friclaw=4, rate-and-state with
strong rate weakening) time loop. Elastic kernels (velDispUpdate,
assembleGlobalKU interior+PML, hrglss) are IDENTICAL to python/eqdyna/port.py
(tpv8, friclaw=1) -- only `faulting` differs, so this module duplicates
those kernels verbatim rather than refactor port.py, per rule 1 ("minimal
changes"): the tpv8 port stays untouched and this file is self-contained.

RSF-specific subtleties reproduced exactly (see faulting.f90/fric.f90):
  - solveRSF (friclaw>=3) redefines nsdSliprateVector(4) as the S,D-ONLY
    magnitude (sqrt(s^2+d^2), dropping the normal component) -- and this
    mutated value, not the original 3-component magnitude, is what
    storeRuptureTime sees. fric(71-77) are written from the ORIGINAL
    3-component magnitude in getNsdSlipSliprateTraction, before this
    mutation, so they are unaffected.
  - The two "pre-loop" calls (rate_state_slip_law and
    rate_state_normal_stress, called once before NewtonRaphson with the
    raw measured sliprate as the trial value) mutate fric(20)/fric(23) as
    a side effect in the Fortran, but those mutations are provably
    discarded: fric(20)/fric(23) are overwritten again immediately after
    NewtonRaphson returns, using NewtonRaphson's OWN internal state that
    was seeded from a snapshot taken BEFORE the pre-loop mutation. Only
    the pre-loop xmu (for taoc_old) has a lasting effect; theta_pc_dot
    (fric(24)) is fully dead (never read downstream, not in frt.txt) and
    is not computed here at all.
  - fric(31-36) (vxm/vym/vzm/vxs/vys/vzs) ARE written for friclaw>=3
    (unlike friclaw<=2, where they stay 0 -- see port.py's comment on the
    same columns) -- this is the "burned people" column faulting.f90/
    pathway_forward.md flags; get the branch (friclaw==5 vs <5) right:
    tpv104 is friclaw==4, so the theta_pc/NewtonRaphson friclaw<5 branch
    IS exercised (taoc_old = xmu*theta_pc, not fric(4)-xmu*tnrm).
  - solveRSF REPLACES nodalForceArr at the fault node pair (not += like
    solveSWTW) with an acceleration-consistent value derived from the
    Newton-converged slip rate.

Newton-Raphson: reproduced as a fixed 20-iteration loop with a per-node
"converged" mask. Because state/theta_pc/xmu/taoc_new/rsfeq/drsfeqdv are
PURE functions of (v_trial, the frozen state0/thetaPc0 baseline, and the
fixed inputs) -- not accumulated across iterations -- freezing v_trial for
a node that has converged makes every subsequent iteration's recompute for
that node idempotent (same inputs -> bit-identical outputs). So "always
recompute for all nodes, then mask the v_trial update by the *converged*
flag (using the pre-this-iteration flag so a newly-converging node does
NOT get this iteration's update applied, matching Fortran's `exit` before
the update)" reproduces the scalar early-exit Fortran exactly, including
the fallback for nodes that never converge within 20 iterations (their
v_trial DOES get the 20th update applied, because Fortran's `do` loop runs
its full body on iv==ivmax too when it doesn't exit).
"""
import numpy as np
from loading import load, region_damp  # noqa: F401
import kernels_numpy


def run(S, nsteps=None, verbose=True):
    dt = S['dt']; rdampk = S['rdampk']
    eq_ids = S['eq_ids']
    nsteps = nsteps or S['nstep']

    inv = kernels_numpy.build(S)
    v1, velArr, dispArr, force, stress_i, s_p = kernels_numpy.init_state(
        inv, inv['E_int'].shape[0], inv['E_pml'].shape[0])
    mass_full = np.concatenate(([1.0], S['nodalMassArr']))
    inv_mass_full = np.where(mass_full > 0, 1.0 / np.where(mass_full > 0, mass_full, 1.0), 0.0)

    # faulting.f90 solveRSF, line ~201: "if (insertFaultType>0 .and.
    # C_elastic==1)" clamp on the normal traction for non-planar/elastic
    # faults. min_norm/max_norm are globalvar.f90 module-level constants
    # (-10.0d6/-40.0d6 Pa), NOT case parameters -- reproduced verbatim,
    # including the non-clamping gap between them (min_norm > max_norm, so
    # the two branches are mutually exclusive and values strictly between
    # max_norm and min_norm pass through untouched, exactly like the
    # Fortran if/elseif).
    insertFaultType = S.get('insertFaultType', 0)
    min_norm = -10.0e6
    max_norm = -40.0e6

    nftnd = S['nftnd']; nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    un = S['un']; us = S['us']; ud = S['ud']; arn = S['arn']
    fric = S['fric_init'].copy()
    fnft = np.full(nftnd, 99999.0)
    timeElapsed = 0.0
    slipRateThres = S['slipRateThres']; C_elastic = S['C_elastic']
    idxF_s = [eq_ids[nsmp1, d] for d in range(3)]
    idxF_m = [eq_ids[nsmp2, d] for d in range(3)]
    massSlave = S['fnms'][nsmp1]; massMaster = S['fnms'][nsmp2]
    mr = massMaster * massSlave / (massMaster + massSlave)

    xsource, ysource, zsource = S['xsource'], S['ysource'], S['zsource']
    nucR, nucT, nucdtau0, TPV = S['nucR'], S['nucT'], S['nucdtau0'], S['TPV']
    nsmp1_coor = S['meshCoor'][nsmp1]
    radius = np.sqrt((nsmp1_coor[:, 0] - xsource) ** 2 + (nsmp1_coor[:, 1] - ysource) ** 2 +
                      (nsmp1_coor[:, 2] - zsource) ** 2)
    C_nuclea = S['C_nuclea']; nucfault = S['nucfault']  # ift is always 1 here (ntotft==1)

    def rsf_slip_law(V2, psi, a, b, Dc, f0, V0, fw, Vw):
        tmpc = 1.0 / (2.0 * V0) * np.exp(psi / a)
        tmp = (V2 + 1.0e-30) * tmpc
        xmu = a * np.log(tmp + np.sqrt(tmp ** 2 + 1.0))
        dxmudv = a * tmpc / np.sqrt(1.0 + tmp ** 2)
        fLV = f0 - (b - a) * np.log(V2 / V0)
        fss = fw + (fLV - fw) / ((1.0 + (V2 / Vw) ** 8) ** 0.125)
        fssa = fss / a
        psiss = a * np.log(2.0 * V0 / V2 * (np.exp(fssa) - np.exp(-fssa)) / 2.0)
        psi_new = psiss + (psi - psiss) * np.exp(-V2 * dt / Dc)
        return xmu, dxmudv, psi_new

    def rsf_normal_stress(V2, theta_pc, tnrm, L):
        theta_pc_dot = -V2 / L * (theta_pc - np.abs(tnrm))
        return theta_pc + theta_pc_dot * dt

    for nt in range(1, nsteps + 1):
        timeElapsed += dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_numpy.elastic_step(
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
        srMag_full = np.sqrt(srN0 ** 2 + srS0 ** 2 + srD0 ** 2)  # 3-component, for fric(71-77)

        fric[:, 70] = slipS0; fric[:, 71] = slipD0; fric[:, 72] = slipN0
        fric[:, 73] = srS0; fric[:, 74] = srD0
        fric[:, 75] = np.maximum(fric[:, 75], srMag_full)
        fric[:, 76] = fric[:, 76] + srMag_full * dt

        totalMass = (massSlave + massMaster) * arn
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic

        # ---- rsfNucleation (friclaw>=3, C_nuclea==1, ift==nucfault==1 always here) ----
        F = np.where(radius < nucR, np.exp(radius ** 2 / (radius ** 2 - nucR ** 2)), 0.0)
        G = np.where(timeElapsed <= nucT,
                     np.exp((timeElapsed - nucT) ** 2 / (timeElapsed * (timeElapsed - 2.0 * nucT))), 1.0)
        dtau_nuc = nucdtau0 * F * G if (TPV == 104 or TPV == 105) else 0.0 * F
        Ts = Ts + dtau_nuc

        Tmag = np.sqrt(Ts ** 2 + Td ** 2)  # shear magnitude AFTER nucleation perturbation

        # ---- solveRSF ----
        tnrm = Tn + fric[:, 5]  # friclaw==4<5 -> += fric(6)
        if insertFaultType > 0 and C_elastic == 1:
            tnrm = np.where(tnrm >= min_norm, min_norm,
                            np.where(tnrm <= max_norm, max_norm, tnrm))
        tnrm = np.where(tnrm > 0.0, 0.0, tnrm)

        slipN = slipN0 + fric[:, 24] * timeElapsed
        slipS = slipS0 + fric[:, 25] * timeElapsed
        slipD = slipD0 + fric[:, 26] * timeElapsed
        srN = srN0 + fric[:, 24]; srS = srS0 + fric[:, 25]; srD = srD0 + fric[:, 26]
        srMag_sd = np.sqrt(srS ** 2 + srD ** 2)  # RSF-redefined magnitude (S,D only) -- what storeRuptureTime sees

        v_trial = srMag_sd
        a = fric[:, 8]; b = fric[:, 9]; Dc = fric[:, 10]; V0 = fric[:, 11]; f0 = fric[:, 12]
        fw = fric[:, 13]; Vw = fric[:, 14]
        theta_pc_tmp = fric[:, 22].copy()  # snapshot BEFORE the (discarded) pre-loop mutation
        state0 = fric[:, 19].copy()  # snapshot BEFORE the (discarded) pre-loop mutation
        xmu_pre, _, _ = rsf_slip_law(v_trial, state0, a, b, Dc, f0, V0, fw, Vw)
        taoc_old = xmu_pre * theta_pc_tmp

        T_coeff = arn * dt / mr
        trialS = Ts - taoc_old * 0.5 * (srS / v_trial) + fric[:, 25] / T_coeff
        trialD = Td - taoc_old * 0.5 * (srD / v_trial) + fric[:, 26] / T_coeff
        trialMag = np.sqrt(trialS ** 2 + trialD ** 2)

        # ---- NewtonRaphson: fixed 20 iterations, per-node convergence mask ----
        done = np.zeros_like(v_trial, dtype=bool)
        state_out = state0.copy(); thetaPc_out = theta_pc_tmp.copy(); taoc_new = np.zeros_like(v_trial)
        for _ in range(20):
            xmu, dxmudv, state_this = rsf_slip_law(v_trial, state0, a, b, Dc, f0, V0, fw, Vw)
            thetaPc_this = rsf_normal_stress(v_trial, theta_pc_tmp, tnrm, Dc)
            taoc_this = xmu * thetaPc_this
            rsfeq = v_trial + T_coeff * (taoc_this * 0.5 - trialMag)
            drsfeqdv = 1.0 + T_coeff * (dxmudv * thetaPc_this) * 0.5
            exit_now = (np.abs(rsfeq / drsfeqdv) < 1.0e-14 * np.abs(v_trial)) & \
                       (np.abs(rsfeq) < 1.0e-6 * np.abs(v_trial)) & (~done)
            done_after = done | exit_now
            newSliprate = v_trial - rsfeq / drsfeqdv
            v_next = np.where(newSliprate <= 0.0, v_trial / 2.0, newSliprate)
            v_trial = np.where(done_after, v_trial, v_next)
            state_out = state_this; thetaPc_out = thetaPc_this; taoc_new = taoc_this
            done = done_after

        fric[:, 19] = state_out; fric[:, 22] = thetaPc_out

        Ts = taoc_old * 0.5 * (srS / srMag_sd) + taoc_new * 0.5 * (trialS / trialMag)
        Td = taoc_old * 0.5 * (srD / srMag_sd) + taoc_new * 0.5 * (trialD / trialMag)
        Tn_final = tnrm

        fric[:, 77] = Tn_final; fric[:, 78] = Ts; fric[:, 79] = Td
        fric[:, 46] = v_trial
        fric[:, 47] = np.sqrt(Ts ** 2 + Td ** 2)

        accN = -srN / dt - slipN / dt / dt
        accS = (v_trial * (trialS / trialMag) - srS) / dt
        accD = (v_trial * (trialD / trialMag) - srD) / dt
        xAcc = accN * un[:, 0] + accS * us[:, 0] + accD * ud[:, 0]
        yAcc = accN * un[:, 1] + accS * us[:, 1] + accD * ud[:, 1]
        zAcc = accN * un[:, 2] + accS * us[:, 2] + accD * ud[:, 2]
        xR = fx_s + fx_m; yR = fy_s + fy_m; zR = fz_s + fz_m

        force[idxF_s[0]] = (-xAcc + xR / massMaster) * mr
        force[idxF_s[1]] = (-yAcc + yR / massMaster) * mr
        force[idxF_s[2]] = (-zAcc + zR / massMaster) * mr
        force[idxF_m[0]] = (xAcc + xR / massSlave) * mr
        force[idxF_m[1]] = (yAcc + yR / massSlave) * mr
        force[idxF_m[2]] = (zAcc + zR / massSlave) * mr

        fric[:, 30] = vx_m + (xAcc + xR / massSlave) * dt
        fric[:, 31] = vy_m + (yAcc + yR / massSlave) * dt
        fric[:, 32] = vz_m + (zAcc + zR / massSlave) * dt
        fric[:, 33] = vx_s + (-xAcc + xR / massMaster) * dt
        fric[:, 34] = vy_s + (-yAcc + yR / massMaster) * dt
        fric[:, 35] = vz_s + (-zAcc + zR / massMaster) * dt

        force[0] = 0.0

        need = fnft > 5000.0
        fnft = np.where(need & (srMag_sd >= slipRateThres), timeElapsed, fnft)

        force[1:] = force[1:] * inv_mass_full[1:]
        if verbose and nt % 10 == 0:
            print('step', nt, '/', nsteps, 'max slip rate', v_trial.max())

    return dict(velArr=velArr, dispArr=dispArr, fnft=fnft, fric=fric, force=force)
