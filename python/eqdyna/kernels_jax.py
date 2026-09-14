"""
JAX elastic time-step kernel shared by every friclaw's JAX port
(port_jax.py/port_rsf_jax.py/port_tp_jax.py). Before this module existed,
the same `velDispUpdate` + `assembleGlobalKU` (interior + PML) + `hrglss`
block (functionally identical to kernels_numpy.py's NumPy version, just
expressed with `.at[].set()`/`.at[].add()` instead of in-place mutation)
was duplicated verbatim in all three files.

Risk-list note (mira's port risk #1): this module does NOT touch any
scan-carry tuple. Each port's `step(carry, _)` still unpacks/repacks its
OWN carry tuple (base 9-tuple for port_jax.py/port_rsf_jax.py; +2 history
arrays for port_tp_jax.py) in exactly the order it always did -- this
function only replaces the elastic sub-block's body, taking/returning the
6 arrays it touches (v1, velArr, dispArr, force, stress_i, s_p) as plain
arguments, never the carry tuple itself. Stops right after the
post-hourglass force scrub (`force.at[0].set(0.0)`), matching
kernels_numpy.py's elastic_step -- friction-specific force overwrite and
the final mass multiply stay in each port's own step().

Milestone 10 (drv.a6, C_elastic==0, friclaw==4): mirrors kernels_numpy.py's
gravity body-force term and Drucker-Prager return-mapping EXACTLY, statement
for statement (see that module's top docstring for the full derivation/
no-op proof) -- `inv['C_elastic']` is a static Python int at trace time
(never a jnp array, see port_jax.py's build()), so `if inv['C_elastic']==0:`
is resolved once at trace time, identical to kernels_numpy.py's runtime
branch: skipped entirely, not merely masked, for every C_elastic==1 caller.
"""
import jax.numpy as jnp


def _c(dN, v):
    return (dN * v).sum(axis=1)


def elastic_step(inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk):
    """velDispUpdate + assembleGlobalKU (interior + PML) + hrglss (C_hg==1,
    all elements), functional/JAX style. `inv` is the dict built by
    port_jax.py's `build(S)` (shared, unchanged, by all three JAX ports)."""
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

    # ---- assembleGlobalKU: interior elements ----
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

    if inv['C_elastic'] == 0:
        # calcElemKU.f90:127-161, mirrors kernels_numpy.py's elastic_step
        # statement for statement (see that module's docstring).
        ccosphi = inv['ccosphi']; sinphi = inv['sinphi']; tv = inv['tv']
        strmea = (stress_i[:, 0] + stress_i[:, 1] + stress_i[:, 2]) / 3.0
        strdev0 = stress_i[:, 0] - strmea
        strdev1 = stress_i[:, 1] - strmea
        strdev2 = stress_i[:, 2] - strmea
        strdev3 = stress_i[:, 3]; strdev4 = stress_i[:, 4]; strdev5 = stress_i[:, 5]
        taomax = jnp.sqrt(0.5 * (strdev0 ** 2 + strdev1 ** 2 + strdev2 ** 2)
                           + strdev3 ** 2 + strdev4 ** 2 + strdev5 ** 2)
        porep = 0.0  # eleporep is provably always exactly 0.0 -- see kernels_numpy.py.
        yld = jnp.maximum(ccosphi - sinphi * (strmea + porep), 0.0)
        mask = taomax > yld
        safe_tao = jnp.where(mask, taomax, 1.0)
        ratio = yld / safe_tao
        rjust = jnp.where(mask, ratio + (1.0 - ratio) * jnp.exp(-dt / tv), 1.0)
        stress_i = stress_i.at[:, 0].set(strdev0 * rjust + strmea)
        stress_i = stress_i.at[:, 1].set(strdev1 * rjust + strmea)
        stress_i = stress_i.at[:, 2].set(strdev2 * rjust + strmea)
        stress_i = stress_i.at[:, 3].set(strdev3 * rjust)
        stress_i = stress_i.at[:, 4].set(strdev4 * rjust)
        stress_i = stress_i.at[:, 5].set(strdev5 * rjust)

    st1 = constk_w_i * (stress_i[:, 0] + rdampk * strr1)
    st2 = constk_w_i * (stress_i[:, 1] + rdampk * strr2)
    st3 = constk_w_i * (stress_i[:, 2] + rdampk * strr3)
    st4 = constk_w_i * (stress_i[:, 3] + rdampk * strr4)
    st5 = constk_w_i * (stress_i[:, 4] + rdampk * strr5)
    st6 = constk_w_i * (stress_i[:, 5] + rdampk * strr6)
    Fx = dNx_i * st1[:, None] + dNz_i * st5[:, None] + dNy_i * st6[:, None]
    Fy = dNy_i * st2[:, None] + dNz_i * st4[:, None] + dNx_i * st6[:, None]
    Fz = dNz_i * st3[:, None] + dNy_i * st4[:, None] + dNx_i * st5[:, None]
    Fz = Fz - (inv['m_e_i'] * inv['grav_const'])[:, None]
    scat_idx += [inv['idxIx'], inv['idxIy'], inv['idxIz']]
    scat_val += [Fx.ravel(), Fy.ravel(), Fz.ravel()]

    # ---- assembleGlobalKU: PML elements ----
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
        s0 = [rdampk * srate[k] + inv['pml_init6'][:, k] for k in range(6)]
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
        grav_p = (inv['m_e_p'] * inv['grav_const'])[:, None]
        grp_sums = [efPML12[0] + efPML12[1] + efPML12[2] + efPML12[9],
                    efPML12[3] + efPML12[4] + efPML12[5] + efPML12[10],
                    efPML12[6] + efPML12[7] + efPML12[8] + efPML12[11] - grav_p]
        for g in range(3):
            scat_idx.append(inv['idxP3'][g]); scat_val.append(grp_sums[g].ravel())

    # ---- hrglss (C_hg==1, all elements) ----
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
    force = jnp.zeros(inv['NEQ1'], dtype=jnp.float64).at[all_idx].add(all_val)
    force = force.at[0].set(0.0)

    return v1, velArr, dispArr, force, stress_i, s_p
