"""
NumPy elastic time-step kernel shared by every friclaw's NumPy port
(port.py/port_rsf.py/port_tp.py). Before this module existed, `velDispUpdate`
+ `assembleGlobalKU` (interior + PML) + `hrglss` (C_hg==1 KF78) were
duplicated verbatim in all three files (only `faulting`/`thermop` differ
per friclaw) -- this extracts that identical block once.

Stops right after the post-hourglass force scrub (`force[0] = 0.0`):
friction-specific code (the fault-node force overwrite/add, and the final
mass divide) stays in each port's own `run()`, because port.py divides
`force` by `mass` directly while port_rsf.py/port_tp.py multiply by a
precomputed reciprocal (`inv_mass_full`) -- NOT bit-identical in floating
point, and a genuine pre-existing difference between the ports that this
refactor does not unify.

Milestone 10 (drv.a6, C_elastic==0, friclaw==4) adds two C_elastic-gated
terms, both PROVABLE no-ops for every C_elastic==1 caller (port.py/
port_tp.py, and port_rsf.py for tpv104/tpv10):

  - Gravity body force (src/assembleGlobalKU.f90:16 al(3,:) fed through
    src/calcElemMass.f90 into elresf, BEFORE calcElemKU's stress-based
    `work` is added on top of it -- floating-point addition being
    commutative, adding this constant term into Fz/the PML z-group BEFORE
    the same single per-element scatter this kernel already does
    reproduces the combined elresf, at worst a few-ULP reduction-order
    difference from Fortran's own per-element grouping, well inside the
    roundoff already tolerated by every existing per-case ABS_BOUND).
    `grav_const = (1-C_elastic)*grav*(roumax-(gamar+1)*rhow)/roumax` is
    EXACTLY 0.0 (bitwise) when C_elastic==1, so the term added is `x - 0.0
    == x` exactly -- not merely small, IDENTICAL.
    Derivation of `-m_e*grav_const` (not a per-element/per-node loop):
    assembleGlobalMass.f90's assembleElementMassDetShg adds the SAME
    per-element lumped mass `elementMass(3*(i-1)+ixyz)` into nodalMassArr
    at ALL 12 dof slots for a PML node (including the vhg group, ixyz+9)
    or all 3 for an interior node -- so nodalMassArr at any node's
    z-type equation slot already equals sum-over-touching-elements of
    that same per-element mass `m_e` (`contm`'s `constm*det`, identical
    for every node/direction of an element) -- meaning the gravity
    contribution at that slot is just `-nodalMassArr[slot]*grav_const`;
    this kernel instead adds the PER-ELEMENT `-m_e*grav_const` before the
    scatter (matching Fortran's per-element accumulation order), which
    sums to the identical total.
  - Drucker-Prager viscoplastic return-mapping (src/calcElemKU.f90:127-161,
    interior elements only -- PML elements never go through calcElemKU),
    applied to `stress_i` right after the elastic trial-stress increment,
    gated by a plain Python `if C_elastic == 0:` (a static per-run scalar,
    never a per-element array) -- for C_elastic==1 this entire block does
    not execute, not merely masked to a no-op array-wise, so it cannot
    perturb any previously-verified case even at the ULP level.
    `porep` (pore pressure) is hardcoded 0.0: src/meshgen.f90's
    setPlasticStress and src/eqdyna3d.f90's zero-init are the ONLY writes
    to Fortran's `eleporep` anywhere in src/*.f90 (grep-verified) -- the
    lithostatic pore-pressure formula is commented out in the Fortran
    itself, so eleporep is provably always exactly 0.0, not assumed.
    `pstrain`/`pstrinc` (the plastic-strain accumulator) is NOT computed:
    it is write-only (never read back by any downstream physics) and
    test.reference.results/test.drv.a6 has no pstr.txt* output to gate
    against (only frt.txt1/frt.txt3) -- computing it would be unverifiable
    scope creep, flagged here rather than silently added.
"""
import numpy as np


def _c(dN, v):
    """Specialized replacement for np.einsum('ei,ei->e', dN, v): a plain
    broadcast-multiply + axis-sum."""
    return (dN * v).sum(axis=1)


def build(S):
    """Precompute every loop-invariant array ONCE: mesh connectivity split
    by interior/PML element, equation-index gathers for the velocity/force
    scatter-gathers, and the static PML damping factors from region_damp.
    Returns a dict consumed by `elastic_step`."""
    from loading import region_damp

    N = S['N']; dt = S['dt']; w = S['w']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']

    is_int = elemType == 1
    is_pml = elemType == 2
    E_int = np.nonzero(is_int)[0]
    E_pml = np.nonzero(is_pml)[0]

    int_nodes = ndof == 3
    idx3_v = eq_ids[int_nodes, 0:3]

    pml_nodes = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes, 0], S['meshCoor'][pml_nodes, 1],
                                 S['meshCoor'][pml_nodes, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'])
    dampv_pml = np.zeros((pml_nodes.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes, :]
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

    # ---- Milestone 10 (drv.a6, C_elastic==0): gravity + Drucker-Prager ----
    # grav_const is EXACTLY 0.0 (bitwise) whenever C_elastic==1 -- see this
    # module's top docstring for the no-op proof and the nodalMassArr
    # derivation.
    C_elastic = S['C_elastic']
    grav_const = ((1.0 - C_elastic) * S['grav']
                  * (S['roumax'] - (S['gamar'] + 1.0) * S['rhow']) / S['roumax'])
    m_e_all = mat[:, 2] * eledet  # contm's per-element lumped mass (same for every node/direction)
    m_e_i = m_e_all[E_int]
    m_e_p = m_e_all[E_pml]

    if C_elastic == 0:
        init_stress = S['init_stress']  # (E,6), main.py's setPlasticStress port -- raises
        stress_i0 = init_stress[E_int].copy()  # (KeyError) loudly if missing, never a silent zero.
        pml_init6 = init_stress[E_pml].copy()
        ccosphi = S['ccosphi']; sinphi = S['sinphi']; tv = S['tv']
    else:
        stress_i0 = np.zeros((E_int.shape[0], 6))
        pml_init6 = np.zeros((E_pml.shape[0], 6))
        ccosphi = sinphi = tv = 0.0

    return dict(
        N=N, dt=dt, w=w, NEQ=NEQ, NEQ1=NEQ + 1, conn=conn, phi=phi, ss=ss,
        int_nodes=int_nodes, idx3_v=idx3_v, pml_nodes=pml_nodes, idx12_v=idx12_v,
        a9=a9, b9=b9,
        E_int=E_int, lam_i=lam_i, miu_i=miu_i, dNx_i=dNx_i, dNy_i=dNy_i, dNz_i=dNz_i,
        constk_w_i=constk_w_i, conn_i=conn_i, idxIx=idxIx, idxIy=idxIy, idxIz=idxIz,
        E_pml=E_pml, lam_p=lam_p, miu_p=miu_p, dNx_p=dNx_p, dNy_p=dNy_p, dNz_p=dNz_p,
        det_w_p=det_w_p, conn_p=conn_p, a1=a1, b1=b1, a2=a2, b2=b2, a3=a3, b3=b3,
        idxP12=idxP12, idxP3=idxP3, idxH0=idxH0, idxH1=idxH1, idxH2=idxH2,
        C_elastic=C_elastic, grav_const=grav_const, m_e_i=m_e_i, m_e_p=m_e_p,
        stress_i0=stress_i0, pml_init6=pml_init6, ccosphi=ccosphi, sinphi=sinphi, tv=tv,
    )


def init_state(inv, E_int_count, E_pml_count):
    """Fresh step state for a new run: velArr/dispArr/v1/force/s_p start at
    zero; `stress_i` starts at `inv['stress_i0']` -- exactly zeros for
    C_elastic==1 (unchanged behavior for every existing caller), or
    meshgen.f90's setPlasticStress lithostatic pre-stress for C_elastic==0
    (Milestone 10, see kernels_numpy.py's build())."""
    N = inv['N']; NEQ1 = inv['NEQ1']
    velArr = np.zeros((N, 3)); dispArr = np.zeros((N, 3))
    v1 = np.zeros(NEQ1)
    force = np.zeros(NEQ1)
    stress_i = inv['stress_i0'].copy()
    s_p = np.zeros((E_pml_count, 15))
    return v1, velArr, dispArr, force, stress_i, s_p


def elastic_step(inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk):
    """velDispUpdate + assembleGlobalKU (interior + PML) + hrglss (C_hg==1,
    all elements). Mutates v1/velArr/dispArr/force/stress_i/s_p IN PLACE
    (matching the per-file loops this replaces) and returns them. Ends
    right after the post-hourglass force scrub -- see module docstring for
    why the fault-force/mass-divide step is NOT included here."""
    idx3_v = inv['idx3_v']; int_nodes = inv['int_nodes']
    idx12_v = inv['idx12_v']; pml_nodes = inv['pml_nodes']
    a9 = inv['a9']; b9 = inv['b9']

    # ---- velDispUpdate ----
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

    # ---- assembleGlobalKU: interior elements ----
    conn_i = inv['conn_i']; dNx_i, dNy_i, dNz_i = inv['dNx_i'], inv['dNy_i'], inv['dNz_i']
    lam_i, miu_i, constk_w_i = inv['lam_i'], inv['miu_i'], inv['constk_w_i']
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

    if inv['C_elastic'] == 0:
        # calcElemKU.f90:127-161 Drucker-Prager viscoplastic return-mapping.
        # Static Python `if` (not a per-element mask): this entire block is
        # SKIPPED, not merely a no-op array-wise, for every C_elastic==1
        # caller -- see build()'s docstring note.
        ccosphi = inv['ccosphi']; sinphi = inv['sinphi']; tv = inv['tv']
        strmea = (stress_i[:, 0] + stress_i[:, 1] + stress_i[:, 2]) / 3.0
        strdev0 = stress_i[:, 0] - strmea
        strdev1 = stress_i[:, 1] - strmea
        strdev2 = stress_i[:, 2] - strmea
        strdev3 = stress_i[:, 3]; strdev4 = stress_i[:, 4]; strdev5 = stress_i[:, 5]
        taomax = np.sqrt(0.5 * (strdev0 ** 2 + strdev1 ** 2 + strdev2 ** 2)
                          + strdev3 ** 2 + strdev4 ** 2 + strdev5 ** 2)
        porep = 0.0  # eleporep is provably always exactly 0.0 -- see build()'s docstring.
        yld = ccosphi - sinphi * (strmea + porep)
        yld = np.maximum(yld, 0.0)
        mask = taomax > yld
        safe_tao = np.where(mask, taomax, 1.0)  # avoid 0/0 on non-yielding elements
        ratio = yld / safe_tao
        # rjust==1.0 where not yielding reproduces Fortran's "stress
        # unchanged" branch EXACTLY (strdev*1+strmea == the pre-correction
        # stress for i<=3; strdev*1 == it for i=4,5), so no separate
        # np.where is needed on the stress-assignment lines below.
        rjust = np.where(mask, ratio + (1.0 - ratio) * np.exp(-dt / tv), 1.0)
        stress_i[:, 0] = strdev0 * rjust + strmea
        stress_i[:, 1] = strdev1 * rjust + strmea
        stress_i[:, 2] = strdev2 * rjust + strmea
        stress_i[:, 3] = strdev3 * rjust
        stress_i[:, 4] = strdev4 * rjust
        stress_i[:, 5] = strdev5 * rjust

    st1 = constk_w_i * (stress_i[:, 0] + rdampk * strr1)
    st2 = constk_w_i * (stress_i[:, 1] + rdampk * strr2)
    st3 = constk_w_i * (stress_i[:, 2] + rdampk * strr3)
    st4 = constk_w_i * (stress_i[:, 3] + rdampk * strr4)
    st5 = constk_w_i * (stress_i[:, 4] + rdampk * strr5)
    st6 = constk_w_i * (stress_i[:, 5] + rdampk * strr6)
    Fx = dNx_i * st1[:, None] + dNz_i * st5[:, None] + dNy_i * st6[:, None]
    Fy = dNy_i * st2[:, None] + dNz_i * st4[:, None] + dNx_i * st6[:, None]
    Fz = dNz_i * st3[:, None] + dNy_i * st4[:, None] + dNx_i * st5[:, None]
    # calcElemMass's gravity body force (al(3,:), z-direction only) -- see
    # build()'s docstring for the -m_e*grav_const derivation; EXACTLY 0.0
    # when C_elastic==1.
    Fz = Fz - (inv['m_e_i'] * inv['grav_const'])[:, None]
    scat_idx += [inv['idxIx'], inv['idxIy'], inv['idxIz']]
    scat_val += [Fx.ravel(), Fy.ravel(), Fz.ravel()]

    # ---- assembleGlobalKU: PML elements ----
    E_pml = inv['E_pml']
    if E_pml.shape[0]:
        conn_p = inv['conn_p']; dNx_p, dNy_p, dNz_p = inv['dNx_p'], inv['dNy_p'], inv['dNz_p']
        lam_p, miu_p, det_w_p = inv['lam_p'], inv['miu_p'], inv['det_w_p']
        a1, b1, a2, b2, a3, b3 = inv['a1'], inv['b1'], inv['a2'], inv['b2'], inv['a3'], inv['b3']
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
        sxy = s_p[:, 9] + s_p[:, 10]
        sxz = s_p[:, 11] + s_p[:, 12]
        syz = s_p[:, 13] + s_p[:, 14]
        # s0(1:6) == Fortran's s(16:21): setPlasticStress's lithostatic
        # pre-stress, a CONSTANT never reassigned inside calcPMLElemKU
        # itself (confirmed by reading assembleGlobalKU.f90's
        # calcPMLElemKU -- it only ever READS s(16:21), added to
        # rdampk*stressrate) -- exactly 0 for C_elastic==1 (pml_init6 is
        # all-zero then, see build()).
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
        # calcElemMass's gravity body force lands in the z ("vhg_z") slot
        # (efPML12[11]) BEFORE calcPMLElemKU subtracts from it -- added
        # here to the z group sum instead (addition is commutative);
        # EXACTLY 0.0 when C_elastic==1.
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

    all_idx = np.concatenate(scat_idx)
    all_val = np.concatenate(scat_val)
    force += np.bincount(all_idx, weights=all_val, minlength=inv['NEQ1'])
    force[0] = 0.0  # scrub sink again (bincount may have accumulated masked-out contributions there)

    return v1, velArr, dispArr, force, stress_i, s_p
