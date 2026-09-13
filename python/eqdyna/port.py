"""
Vectorized NumPy port of EQdyna's explicit time loop for tpv8 (friclaw=1,
serial/np=1). NOT a from-scratch mesh generator: mesh/mass/shape-function
state is loaded from Fortran `pydump_*` files written by `src/pydump.f90`
(instrumentation added for this spike, called once before the time loop
starts). This ports only the per-step kernels: velDispUpdate, assembleGlobalKU
(interior calcElemKU + PML calcPMLElemKU), hrglss (C_hg==1 KF78), faulting
(getNsdSlipSliprateTraction + solveSWTW friclaw==1 + swtwNucleation +
storeRuptureTime). calcElemMass's contribution is analytically zero for this
case (rdampm=0.0d0, C_elastic=1) and is not computed at all -- see README.

UPDATE (2026-09-13 refactor): `load()` and `region_damp()` moved to
eqdyna/loading.py (shared by every friction-family port, not just this
one) and re-exported here for backward compatibility with existing
callers (`from port import load, region_damp`). `region_damp()` also lost
its `bound_inclusive` parameter there: src/func_lib.f90's
`pmlRegionDistance` was standardized to INCLUSIVE (>=/<=) for both callers
(pathway_forward.md item 13, 2026-09-13), so there is no longer a
per-caller distinction to reproduce -- see loading.py's docstring for the
full rationale and the re-run-parity note.
"""
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from loading import load, region_damp  # noqa: F401 (re-exported for existing callers)
import kernels_numpy


def run(S, nsteps=None, verbose=True):
    dt = S['dt']; rdampk = S['rdampk']; NEQ = S['NEQ']
    eq_ids = S['eq_ids']
    nsteps = nsteps or S['nstep']
    NEQ1 = NEQ + 1

    inv = kernels_numpy.build(S)
    v1, velArr, dispArr, force, stress_i, s_p = kernels_numpy.init_state(
        inv, inv['E_int'].shape[0], inv['E_pml'].shape[0])
    mass = np.concatenate(([1.0], S['nodalMassArr']))  # index0 unused (never divided)

    nftnd = S['nftnd']; nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    un = S['un']; us = S['us']; ud = S['ud']; arn = S['arn']
    fric = S['fric_init'].copy()
    fnft = np.full(nftnd, 99999.0)
    timeElapsed = 0.0
    slipRateThres = S['slipRateThres']; C_elastic = S['C_elastic']
    idxF_s = [eq_ids[nsmp1, d] for d in range(3)]
    idxF_m = [eq_ids[nsmp2, d] for d in range(3)]

    for nt in range(1, nsteps + 1):
        timeElapsed += dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_numpy.elastic_step(
            inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk)

        # ---- faulting (friclaw==1: solveSWTW), C_nuclea==1 nucfault==ift==1 ----
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

        slipN = dN_m - dN_s; slipS = dS_m - dS_s; slipD = dD_m - dD_s
        srN = vN_m - vN_s; srS = vS_m - vS_s; srD = vD_m - vD_s
        srMag = np.sqrt(srN ** 2 + srS ** 2 + srD ** 2)

        fric[:, 70] = slipS; fric[:, 71] = slipD; fric[:, 72] = slipN
        fric[:, 73] = srS; fric[:, 74] = srD
        fric[:, 75] = np.maximum(fric[:, 75], srMag)
        fric[:, 76] = fric[:, 76] + srMag * dt

        massSlave = S['fnms'][nsmp1]; massMaster = S['fnms'][nsmp2]
        totalMass = (massSlave + massMaster) * arn
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic
        Tmag = np.sqrt(Ts ** 2 + Td ** 2)

        fs = fric[:, 0]; fd = fric[:, 1]; D0 = fric[:, 2]
        slip = fric[:, 76]
        fricCoeff = np.where(np.abs(slip) < 1.0e-10, fs, fs - (fs - fd) * slip / D0)
        fricCoeff = np.where(slip >= D0, fd, fricCoeff)

        fricCoeff = np.minimum(fs, fricCoeff)  # swtwNucleation no-op for TPV==8 (tr default 1e9)

        effNorm = np.where((Tn + fric[:, 5]) > 0.0, 0.0, Tn + fric[:, 5])
        trialShear = fric[:, 3] - fricCoeff * effNorm
        over = Tmag > trialShear
        scale = np.where(over, trialShear / np.where(Tmag == 0, 1.0, Tmag), 1.0)
        Ts = Ts * scale; Td = Td * scale

        xN = Tn * un[:, 0] + Ts * us[:, 0] + Td * ud[:, 0]
        xN_y = Tn * un[:, 1] + Ts * us[:, 1] + Td * ud[:, 1]
        xN_z = Tn * un[:, 2] + Ts * us[:, 2] + Td * ud[:, 2]
        xTrac = np.stack([xN, xN_y, xN_z], axis=1) * arn[:, None]
        initN, initS, initD = fric[:, 6], fric[:, 7], fric[:, 48]
        xInit = np.stack([initN * un[:, 0] + initS * us[:, 0] + initD * ud[:, 0],
                           initN * un[:, 1] + initS * us[:, 1] + initD * ud[:, 1],
                           initN * un[:, 2] + initS * us[:, 2] + initD * ud[:, 2]], axis=1) * arn[:, None]
        delta_f = xTrac - xInit * C_elastic
        flt_idx = np.concatenate([idxF_s[0], idxF_m[0], idxF_s[1], idxF_m[1], idxF_s[2], idxF_m[2]])
        flt_val = np.concatenate([delta_f[:, 0], -delta_f[:, 0], delta_f[:, 1], -delta_f[:, 1],
                                   delta_f[:, 2], -delta_f[:, 2]])
        force += np.bincount(flt_idx, weights=flt_val, minlength=NEQ1)
        force[0] = 0.0

        fric[:, 77] = Tn; fric[:, 78] = Ts; fric[:, 79] = Td

        need = fnft > 5000.0
        fnft = np.where(need & (srMag >= slipRateThres), timeElapsed, fnft)

        force[0] = 0.0
        force[1:] = force[1:] / mass[1:]
        if verbose and nt % 20 == 0:
            print('step', nt, '/', nsteps, 'max|v|', np.abs(velArr).max())

    return dict(velArr=velArr, dispArr=dispArr, fnft=fnft, fric=fric, force=force)

