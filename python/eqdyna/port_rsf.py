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
from port import load, region_damp  # noqa: F401 (re-exported for callers)


def _c(dN, v):
    return (dN * v).sum(axis=1)


def run(S, nsteps=None, verbose=True):
    N = S['N']; dt = S['dt']; w = S['w']; rdampk = S['rdampk']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']
    nsteps = nsteps or S['nstep']

    is_int = elemType == 1
    E_int = np.nonzero(is_int)[0]
    E_pml = np.nonzero(elemType == 2)[0]

    velArr = np.zeros((N, 3)); dispArr = np.zeros((N, 3))
    v1 = np.zeros(NEQ + 1)
    force = np.zeros(NEQ + 1)
    mass_full = np.concatenate(([1.0], S['nodalMassArr']))
    inv_mass_full = np.where(mass_full > 0, 1.0 / np.where(mass_full > 0, mass_full, 1.0), 0.0)

    int_nodes = ndof == 3
    idx3_v = eq_ids[int_nodes, 0:3]
    pml_nodes = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes, 0], S['meshCoor'][pml_nodes, 1],
                                 S['meshCoor'][pml_nodes, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'], True)
    dampv_pml = np.zeros((pml_nodes.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes, :]
    a9 = 1.0 / dt - dampv_pml / 2.0; b9 = 1.0 / dt + dampv_pml / 2.0

    lam_i = mat[E_int, 3]; miu_i = mat[E_int, 4]
    dN_i = eleshp[E_int]
    dNx_i, dNy_i, dNz_i = dN_i[:, :, 0], dN_i[:, :, 1], dN_i[:, :, 2]
    constk_w_i = (-eledet[E_int]) * w
    stress_i = np.zeros((E_int.shape[0], 6))
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel(); idxIy = eq_ids[conn_i, 1].ravel(); idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p, dNy_p, dNz_p = dN_p[:, :, 0], dN_p[:, :, 1], dN_p[:, :, 2]
    det_w_p = eledet[E_pml] * w
    conn_p = conn[E_pml]
    xc_p = S['meshCoor'][conn_p].mean(axis=1)
    d1p, d2p, d3p = region_damp(xc_p[:, 0], xc_p[:, 1], xc_p[:, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'], False)
    a1 = 1.0 / dt - d1p / 2.0; b1 = 1.0 / dt + d1p / 2.0
    a2 = 1.0 / dt - d2p / 2.0; b2 = 1.0 / dt + d2p / 2.0
    a3 = 1.0 / dt - d3p / 2.0; b3 = 1.0 / dt + d3p / 2.0
    s_p = np.zeros((E_pml.shape[0], 15))
    nd_p = ndof[conn_p]; is12_p = nd_p == 12; is3_p = nd_p == 3
    idxP12 = [np.where(is12_p, eq_ids[conn_p, j], 0).ravel() for j in range(12)]
    idxP3 = [np.where(is3_p, eq_ids[conn_p, g], 0).ravel() for g in range(3)]

    nd_all = ndof[conn]
    slot0 = np.where(nd_all == 3, 0, 9); slot1 = np.where(nd_all == 3, 1, 10); slot2 = np.where(nd_all == 3, 2, 11)
    idxH0 = np.take_along_axis(eq_ids[conn], slot0[:, :, None], axis=2)[:, :, 0].ravel()
    idxH1 = np.take_along_axis(eq_ids[conn], slot1[:, :, None], axis=2)[:, :, 0].ravel()
    idxH2 = np.take_along_axis(eq_ids[conn], slot2[:, :, None], axis=2)[:, :, 0].ravel()

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
    NEQ1 = NEQ + 1

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

        # ---- velDispUpdate (identical to port.py) ----
        accel3 = force[idx3_v]
        v1[idx3_v] = v1[idx3_v] + accel3 * dt
        velArr[int_nodes] = v1[idx3_v]
        dispArr[int_nodes] += v1[idx3_v] * dt
        if pml_nodes.size:
            f9 = force[idx12_v[:, 0:9]]
            vold9 = v1[idx12_v[:, 0:9]]
            v1[idx12_v[:, 0:9]] = (f9 + vold9 * a9) / b9
            f3v9 = force[idx12_v[:, 9:12]]
            v1[idx12_v[:, 9:12]] = v1[idx12_v[:, 9:12]] + f3v9 * dt
            v1[0] = 0.0
            vA = v1[idx12_v[:, 0]] + v1[idx12_v[:, 1]] + v1[idx12_v[:, 2]] + v1[idx12_v[:, 9]]
            vB = v1[idx12_v[:, 3]] + v1[idx12_v[:, 4]] + v1[idx12_v[:, 5]] + v1[idx12_v[:, 10]]
            vC = v1[idx12_v[:, 6]] + v1[idx12_v[:, 7]] + v1[idx12_v[:, 8]] + v1[idx12_v[:, 11]]
            has_eq = idx12_v[:, 0] > 0
            velArr[pml_nodes, 0] = np.where(has_eq, vA, 0.0)
            velArr[pml_nodes, 1] = np.where(has_eq, vB, 0.0)
            velArr[pml_nodes, 2] = np.where(has_eq, vC, 0.0)
            dispArr[pml_nodes, 0] = np.where(has_eq, dispArr[pml_nodes, 0] + vA * dt, 0.0)
            dispArr[pml_nodes, 1] = np.where(has_eq, dispArr[pml_nodes, 1] + vB * dt, 0.0)
            dispArr[pml_nodes, 2] = np.where(has_eq, dispArr[pml_nodes, 2] + vC * dt, 0.0)

        force[:] = 0.0
        scat_idx = []; scat_val = []

        # ---- assembleGlobalKU interior (identical to port.py) ----
        vl = velArr[conn_i]
        sr1 = _c(dNx_i, vl[:, :, 0]); sr2 = _c(dNy_i, vl[:, :, 1]); sr3 = _c(dNz_i, vl[:, :, 2])
        sr4 = _c(dNz_i, vl[:, :, 1]) + _c(dNy_i, vl[:, :, 2])
        sr5 = _c(dNz_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 2])
        sr6 = _c(dNy_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 1])
        vol_sr = sr1 + sr2 + sr3
        strr1 = lam_i * vol_sr + 2 * miu_i * sr1
        strr2 = lam_i * vol_sr + 2 * miu_i * sr2
        strr3 = lam_i * vol_sr + 2 * miu_i * sr3
        strr4 = miu_i * sr4; strr5 = miu_i * sr5; strr6 = miu_i * sr6
        stress_i[:, 0] += strr1 * dt; stress_i[:, 1] += strr2 * dt; stress_i[:, 2] += strr3 * dt
        stress_i[:, 3] += strr4 * dt; stress_i[:, 4] += strr5 * dt; stress_i[:, 5] += strr6 * dt
        st1 = constk_w_i * (stress_i[:, 0] + rdampk * strr1)
        st2 = constk_w_i * (stress_i[:, 1] + rdampk * strr2)
        st3 = constk_w_i * (stress_i[:, 2] + rdampk * strr3)
        st4 = constk_w_i * (stress_i[:, 3] + rdampk * strr4)
        st5 = constk_w_i * (stress_i[:, 4] + rdampk * strr5)
        st6 = constk_w_i * (stress_i[:, 5] + rdampk * strr6)
        Fx = dNx_i * st1[:, None] + dNz_i * st5[:, None] + dNy_i * st6[:, None]
        Fy = dNy_i * st2[:, None] + dNz_i * st4[:, None] + dNx_i * st6[:, None]
        Fz = dNz_i * st3[:, None] + dNy_i * st4[:, None] + dNx_i * st5[:, None]
        scat_idx += [idxIx, idxIy, idxIz]; scat_val += [Fx.ravel(), Fy.ravel(), Fz.ravel()]

        # ---- assembleGlobalKU PML (identical to port.py) ----
        if E_pml.shape[0]:
            vlp = velArr[conn_p]
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
            s_p[:, 0] = ((lam_p + 2 * miu_p) * Dxvx + a1 * s_p[:, 0]) / b1
            s_p[:, 1] = (lam_p * Dyvy + a2 * s_p[:, 1]) / b2
            s_p[:, 2] = (lam_p * Dzvz + a3 * s_p[:, 2]) / b3
            s_p[:, 3] = (lam_p * Dxvx + a1 * s_p[:, 3]) / b1
            s_p[:, 4] = ((lam_p + 2 * miu_p) * Dyvy + a2 * s_p[:, 4]) / b2
            s_p[:, 5] = (lam_p * Dzvz + a3 * s_p[:, 5]) / b3
            s_p[:, 6] = (lam_p * Dxvx + a1 * s_p[:, 6]) / b1
            s_p[:, 7] = (lam_p * Dyvy + a2 * s_p[:, 7]) / b2
            s_p[:, 8] = ((lam_p + 2 * miu_p) * Dzvz + a3 * s_p[:, 8]) / b3
            s_p[:, 9] = (miu_p * Dxvy + a1 * s_p[:, 9]) / b1
            s_p[:, 10] = (miu_p * Dyvx + a2 * s_p[:, 10]) / b2
            s_p[:, 11] = (miu_p * Dxvz + a1 * s_p[:, 11]) / b1
            s_p[:, 12] = (miu_p * Dzvx + a3 * s_p[:, 12]) / b3
            s_p[:, 13] = (miu_p * Dyvz + a2 * s_p[:, 13]) / b2
            s_p[:, 14] = (miu_p * Dzvy + a3 * s_p[:, 14]) / b3
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
                scat_idx.append(idxP12[j]); scat_val.append(efPML12[j].ravel())
            grp_sums = [efPML12[0] + efPML12[1] + efPML12[2] + efPML12[9],
                        efPML12[3] + efPML12[4] + efPML12[5] + efPML12[10],
                        efPML12[6] + efPML12[7] + efPML12[8] + efPML12[11]]
            for g in range(3):
                scat_idx.append(idxP3[g]); scat_val.append(grp_sums[g].ravel())

        # ---- hrglss (identical to port.py) ----
        dl_all = dispArr[conn] + rdampk * velArr[conn]
        for m in range(4):
            phid = (phi[:, m, :, None] * dl_all).sum(axis=1)
            r0 = ss[:, 0] * phid[:, 0] + ss[:, 1] * phid[:, 1] + ss[:, 2] * phid[:, 2]
            r1 = ss[:, 1] * phid[:, 0] + ss[:, 3] * phid[:, 1] + ss[:, 4] * phid[:, 2]
            r2 = ss[:, 2] * phid[:, 0] + ss[:, 4] * phid[:, 1] + ss[:, 5] * phid[:, 2]
            fh0 = -(phi[:, m, :] * r0[:, None]).ravel()
            fh1 = -(phi[:, m, :] * r1[:, None]).ravel()
            fh2 = -(phi[:, m, :] * r2[:, None]).ravel()
            scat_idx += [idxH0, idxH1, idxH2]; scat_val += [fh0, fh1, fh2]

        all_idx = np.concatenate(scat_idx); all_val = np.concatenate(scat_val)
        force += np.bincount(all_idx, weights=all_val, minlength=NEQ1)
        force[0] = 0.0

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
