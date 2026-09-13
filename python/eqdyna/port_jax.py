"""
JAX rewrite of the EQdyna tpv8 (friclaw=1) time loop -- phase 2 of the
feasibility spike. Mirrors python/eqdyna/port.py's NumPy `run()` exactly
(same formulas, same order, same known-bug reproduction for the PML
region cascade), but structured as one jitted, scanned step function
instead of a per-step Python-dispatched loop.

float64 is enabled FIRST, before any other jax import touches an array,
per the phase-2 requirement -- silent float32 would fake a speedup and
break parity against the double-precision Fortran/NumPy reference.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from loading import load, region_damp  # noqa: E402  (must follow the x64 config line)
import kernels_jax


def build(S):
    """Precompute every loop-invariant array ONCE (Python/NumPy side, static
    shapes), then hand them to a closure that lax.scan iterates. Nothing in
    the scanned step function does a nonzero()/boolean-mask reshape -- all
    such shape-determining operations happen here, outside jit, exactly like
    the NumPy port's pre-loop setup."""
    N = S['N']; dt = S['dt']; w = S['w']; rdampk = S['rdampk']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']

    is_int = elemType == 1
    E_int = np.nonzero(is_int)[0]
    E_pml = np.nonzero(elemType == 2)[0]

    int_nodes_idx = np.nonzero(ndof == 3)[0]
    idx3_v = eq_ids[int_nodes_idx, 0:3]

    pml_nodes_idx = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes_idx, 0], S['meshCoor'][pml_nodes_idx, 1],
                                 S['meshCoor'][pml_nodes_idx, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'])
    dampv_pml = np.zeros((pml_nodes_idx.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes_idx, :]
    a9 = 1.0 / dt - dampv_pml / 2.0; b9 = 1.0 / dt + dampv_pml / 2.0

    lam_i = mat[E_int, 3]; miu_i = mat[E_int, 4]
    dN_i = eleshp[E_int]
    dNx_i, dNy_i, dNz_i = dN_i[:, :, 0], dN_i[:, :, 1], dN_i[:, :, 2]
    constk_w_i = (-eledet[E_int]) * w
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel(); idxIy = eq_ids[conn_i, 1].ravel(); idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p, dNy_p, dNz_p = dN_p[:, :, 0], dN_p[:, :, 1], dN_p[:, :, 2]
    det_w_p = eledet[E_pml] * w
    conn_p = conn[E_pml]
    xc_p = S['meshCoor'][conn_p].mean(axis=1)
    d1p, d2p, d3p = region_damp(xc_p[:, 0], xc_p[:, 1], xc_p[:, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'])
    a1 = 1.0 / dt - d1p / 2.0; b1 = 1.0 / dt + d1p / 2.0
    a2 = 1.0 / dt - d2p / 2.0; b2 = 1.0 / dt + d2p / 2.0
    a3 = 1.0 / dt - d3p / 2.0; b3 = 1.0 / dt + d3p / 2.0
    nd_p = ndof[conn_p]; is12_p = nd_p == 12; is3_p = nd_p == 3
    idxP12 = [np.where(is12_p, eq_ids[conn_p, j], 0).ravel() for j in range(12)]
    idxP3 = [np.where(is3_p, eq_ids[conn_p, g], 0).ravel() for g in range(3)]

    nd_all = ndof[conn]
    slot0 = np.where(nd_all == 3, 0, 9); slot1 = np.where(nd_all == 3, 1, 10); slot2 = np.where(nd_all == 3, 2, 11)
    idxH0 = np.take_along_axis(eq_ids[conn], slot0[:, :, None], axis=2)[:, :, 0].ravel()
    idxH1 = np.take_along_axis(eq_ids[conn], slot1[:, :, None], axis=2)[:, :, 0].ravel()
    idxH2 = np.take_along_axis(eq_ids[conn], slot2[:, :, None], axis=2)[:, :, 0].ravel()

    nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    idxF_s = [eq_ids[nsmp1, d] for d in range(3)]
    idxF_m = [eq_ids[nsmp2, d] for d in range(3)]

    mass_full = np.concatenate(([1.0], S['nodalMassArr']))
    inv_mass_full = np.where(mass_full > 0, 1.0 / np.where(mass_full > 0, mass_full, 1.0), 0.0)

    j = jnp.asarray
    inv = dict(
        N=N, dt=dt, rdampk=rdampk, NEQ1=NEQ + 1,
        conn=j(conn), phi=j(phi), ss=j(ss),
        int_nodes_idx=j(int_nodes_idx), idx3_v=j(idx3_v),
        pml_nodes_idx=j(pml_nodes_idx), idx12_v=j(idx12_v), a9=j(a9), b9=j(b9),
        lam_i=j(lam_i), miu_i=j(miu_i), dNx_i=j(dNx_i), dNy_i=j(dNy_i), dNz_i=j(dNz_i),
        constk_w_i=j(constk_w_i), conn_i=j(conn_i), idxIx=j(idxIx), idxIy=j(idxIy), idxIz=j(idxIz),
        lam_p=j(lam_p), miu_p=j(miu_p), dNx_p=j(dNx_p), dNy_p=j(dNy_p), dNz_p=j(dNz_p),
        det_w_p=j(det_w_p), a1=j(a1), b1=j(b1), a2=j(a2), b2=j(b2), a3=j(a3), b3=j(b3),
        idxP12=[j(x) for x in idxP12], idxP3=[j(x) for x in idxP3],
        idxH0=j(idxH0), idxH1=j(idxH1), idxH2=j(idxH2),
        nsmp1=j(nsmp1), nsmp2=j(nsmp2), un=j(S['un']), us=j(S['us']), ud=j(S['ud']), arn=j(S['arn']),
        idxF_s=[j(x) for x in idxF_s], idxF_m=[j(x) for x in idxF_m],
        inv_mass_full=j(inv_mass_full), slipRateThres=S['slipRateThres'], C_elastic=S['C_elastic'],
        Ei=E_int.shape[0], Ep=E_pml.shape[0], E=conn.shape[0],
    )
    return inv


def make_step(inv):
    dt = inv['dt']; rdampk = inv['rdampk']

    def step(carry, _):
        v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
        timeElapsed = timeElapsed + dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_jax.elastic_step(
            inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk)

        # ---- faulting (friclaw==1: solveSWTW; no Newton-Raphson -- that
        # branch (RSF, friclaw>=3) is not part of tpv8's executed path, so
        # nothing here resisted jitting) ----
        nsmp1, nsmp2 = inv['nsmp1'], inv['nsmp2']
        un, us, ud, arn = inv['un'], inv['us'], inv['ud'], inv['arn']
        idxF_s, idxF_m = inv['idxF_s'], inv['idxF_m']
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
        srMag = jnp.sqrt(srN ** 2 + srS ** 2 + srD ** 2)

        fric = fric.at[:, 70].set(slipS); fric = fric.at[:, 71].set(slipD); fric = fric.at[:, 72].set(slipN)
        fric = fric.at[:, 73].set(srS); fric = fric.at[:, 74].set(srD)
        fric = fric.at[:, 75].set(jnp.maximum(fric[:, 75], srMag))
        fric = fric.at[:, 76].add(srMag * dt)

        massSlave = inv['fnms'][nsmp1] if False else inv['massSlave']
        massMaster = inv['massMaster']
        totalMass = (massSlave + massMaster) * arn
        C_elastic = inv['C_elastic']
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic
        Tmag = jnp.sqrt(Ts ** 2 + Td ** 2)

        fs = fric[:, 0]; fd = fric[:, 1]; D0 = fric[:, 2]
        slip = fric[:, 76]
        fricCoeff = jnp.where(jnp.abs(slip) < 1.0e-10, fs, fs - (fs - fd) * slip / D0)
        fricCoeff = jnp.where(slip >= D0, fd, fricCoeff)
        fricCoeff = jnp.minimum(fs, fricCoeff)  # swtwNucleation no-op for TPV==8

        effNorm = jnp.where((Tn + fric[:, 5]) > 0.0, 0.0, Tn + fric[:, 5])
        trialShear = fric[:, 3] - fricCoeff * effNorm
        over = Tmag > trialShear
        scale = jnp.where(over, trialShear / jnp.where(Tmag == 0, 1.0, Tmag), 1.0)
        Ts = Ts * scale; Td = Td * scale

        xN = Tn * un[:, 0] + Ts * us[:, 0] + Td * ud[:, 0]
        xNy = Tn * un[:, 1] + Ts * us[:, 1] + Td * ud[:, 1]
        xNz = Tn * un[:, 2] + Ts * us[:, 2] + Td * ud[:, 2]
        xTrac = jnp.stack([xN, xNy, xNz], axis=1) * arn[:, None]
        initN, initS, initD = fric[:, 6], fric[:, 7], fric[:, 48]
        xInit = jnp.stack([initN * un[:, 0] + initS * us[:, 0] + initD * ud[:, 0],
                            initN * un[:, 1] + initS * us[:, 1] + initD * ud[:, 1],
                            initN * un[:, 2] + initS * us[:, 2] + initD * ud[:, 2]], axis=1) * arn[:, None]
        delta_f = xTrac - xInit * C_elastic

        flt_idx = jnp.concatenate([idxF_s[0], idxF_m[0], idxF_s[1], idxF_m[1], idxF_s[2], idxF_m[2]])
        flt_val = jnp.concatenate([delta_f[:, 0], -delta_f[:, 0], delta_f[:, 1], -delta_f[:, 1],
                                    delta_f[:, 2], -delta_f[:, 2]])
        force = force.at[flt_idx].add(flt_val)
        force = force.at[0].set(0.0)

        fric = fric.at[:, 77].set(Tn); fric = fric.at[:, 78].set(Ts); fric = fric.at[:, 79].set(Td)

        need = fnft > 5000.0
        fnft = jnp.where(need & (srMag >= inv['slipRateThres']), timeElapsed, fnft)

        force = force * inv['inv_mass_full']

        new_carry = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
        return new_carry, None

    return step


def run(S, nsteps=None, verbose=True):
    nsteps = nsteps or S['nstep']
    inv = build(S)
    # a few arrays didn't fit the flat dict-literal above; attach them here.
    conn = S['conn']; elemType = S['elemType']
    E_pml = np.nonzero(elemType == 2)[0]
    inv['conn_p'] = jnp.asarray(conn[E_pml])
    inv['massSlave'] = jnp.asarray(S['fnms'][S['nsmp1']])
    inv['massMaster'] = jnp.asarray(S['fnms'][S['nsmp2']])

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

    step = make_step(inv)
    scan_fn = jax.jit(lambda c: jax.lax.scan(step, c, xs=None, length=nsteps)[0])

    carry0 = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
    carry = scan_fn(carry0)
    jax.block_until_ready(carry)
    v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr), fnft=np.asarray(fnft),
                fric=np.asarray(fric), force=np.asarray(force))
