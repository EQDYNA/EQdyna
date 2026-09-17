"""assembleGlobalKU.py <- src/assembleGlobalKU.f90 (+ calcElemKU.f90,
calcHourglassResist.f90, calcElemMass.f90's gravity body force).

ONE implementation, both backends. This file replaces kernels_numpy.py
(592 lines) and kernels_jax.py (268 lines), which were the same physics
written twice -- and where the JAX copy silently accumulated two different
numerical choices that nobody could see side by side. Those two choices did
not disappear; they are the two named predicates below (_fuse_modes,
_contract_modes), which live HERE, next to the arithmetic they change,
because they are kernel decisions and not array-library spelling. Everything
else in this file is literally one expression serving both backends.

WHAT THE ELEMENT FORCE SCATTER COSTS, AND WHY IT IS BLOCK-AT-A-TIME
There are 148 scatter entries per element (interior force 3x8, hourglass
4 modes x 3 x 8, PML 15x8 on PML elements). Materialising that as one
concatenated index + value pair is 874 MB + 874 MB on test.tpv104. Each
block is therefore built and scattered immediately, and only the arrays
that must be LIVE AT THE SAME TIME exist: the three PML group sums (each
accumulates four f-blocks) plus one in-flight block.

BIT-IDENTITY OF THAT SHAPE: numpy's add.at accumulates in index-array
order and XLA's CPU scatter applies a duplicate index list in array order,
so calling either once per block, in the SAME fixed block order, performs
exactly the same additions on the same partial sums. Block order is
therefore FIXED and is not an implementation detail: interior Fx/Fy/Fz,
then (PML only) f1..f12 and the three PML group sums, then hourglass.
force starts at exactly 0.0 and 0.0 + x == x.

C_elastic==0 (test.drv.a6) adds two terms, both PROVABLE no-ops when
C_elastic==1: the gravity body force, whose (1-C_elastic) factor makes
grav_const EXACTLY 0.0 (so the term added is x - 0.0 == x, identical, not
merely small), and the Drucker-Prager return-mapping, which is skipped by a
static Python if rather than masked array-wise.
"""
import numpy as np

from . import backend as B


def _fuse_modes(xp):
    """Whether hourglass() sums the 4 modes per (element, node) BEFORE the
    scatter, instead of scattering each mode separately.

    A REAL NUMERICAL DIVERGENCE, not a spelling difference, and the one thing
    in this kernel that does not collapse to a single answer. Both orderings
    compute the same sum; they associate it differently, and float addition
    is not associative.

      numpy (False): mode-major -- mode 0's contribution to an equation, then
        mode 1's, ... This is the association calcHourglassResist.f90
        performs, so it is the bit-faithful one, and it is the reference the
        jax column is measured against.

      jax (True): element-major -- 4 modes summed per (element, node) first.
        Takes the scatter from 96 to 24 entries per element; measured -16.0%
        (12.903 -> 10.844 ms/step, tpv8, A100). Gated on magnitude rather
        than bits, against this port's OWN nondeterminism floor: it moves the
        answer LESS than re-running the unmodified jax code does (fric
        2.086e-07 vs 2.235e-07 for base-vs-base-rerun).
    """
    return B.is_jax(xp)


def _contract_modes(xp, phi_m, dl_all):
    """Sum over the 8 element nodes: (E,8) x (E,8,3) -> (E,3).

    A second per-backend choice, and a MEMORY one rather than a numerical
    one: the two spellings were verified bitwise equal on numpy for this
    contraction, so what differs is only what gets materialised.

      numpy: np.einsum('ei,eij->ej', ...) reduces without building the
        (E,8,3) product (45 MB per step on test.tpv104). NOT replaced by a
        stacked matmul: matmul is free to contract with FMA, which is exactly
        how an earlier change here looked identical across 2.8M values and
        then diverged at step 5.

      jax: (phi_m[:,:,None] * dl_all).sum(axis=1) -- XLA fuses the product
        into the reduction so nothing is materialised anyway, and this is the
        form the jax column's committed digests were produced with.
    """
    if B.is_jax(xp):
        return (phi_m[:, :, None] * dl_all).sum(axis=1)
    return xp.einsum('ei,eij->ej', phi_m, dl_all)


def _c(xp, dN, v):
    """Contract an (E,8) shape-derivative block against an (E,8) nodal block
    over the 8 element nodes. (dN*v).sum(axis=1), NOT einsum: the product
    is a contiguous (E,8) temporary so numpy's .sum(axis=1) runs pairwise
    summation over the 8 terms, while einsum accumulates sequentially. They
    disagree in the last bits and this path is gated on bit-identity."""
    return (dN * v).sum(axis=1)


def build(S):
    """Every loop-invariant array, computed ONCE on the host in numpy.

    Host-side for both backends: this is the Fortran's pre-time-loop setup,
    it does the shape-determining nonzero()/boolean-mask work that must not
    happen inside a jitted step, and computing it with jax would buy nothing.
    backend.to_device moves/narrows the result for the jax column.
    """
    from .func_lib import region_damp

    dt = S['dt']; w = S['w']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']

    # assembleGlobalKU.f90:25's dispatch is `elemTypeArr(nel)==1 .or.
    # elemTypeArr(nel)>10` -> calcElemKU (interior kernel), matching
    # `elseif elemTypeArr(nel)==2` -> calcPMLElemKU (PML). elemType 11/12
    # (wedge) and 13 (plain non-degenerate brick, meshgen.py's own
    # docstring) all route through calcElemKU exactly like elemType 1 --
    # confirmed by reading assembleGlobalKU.f90 directly: calcElemKU takes
    # `eleshp(1,1,nel)`/`eledet(nel)` as opaque precomputed arrays and has
    # no elemTypeArr branch of its own (grepped calcElemKU.f90,
    # calcHourglassResist.f90, calcElemMass.f90 -- none reference
    # elemTypeArr), so once those arrays are correct (assembleGlobalMass.py's
    # compute_element_shape) this dispatch is the only place that needs to
    # know about elemType>10 at all. Verified against Fortran directly by
    # testsys/parity/evidence_wedge_kernel.py.
    E_int = np.nonzero((elemType == 1) | (elemType > 10))[0]
    E_pml = np.nonzero(elemType == 2)[0]

    int_nodes_idx = np.nonzero(ndof == 3)[0]
    idx3_v = eq_ids[int_nodes_idx, 0:3]

    pml_nodes_idx = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes_idx, 0],
                                S['meshCoor'][pml_nodes_idx, 1],
                                S['meshCoor'][pml_nodes_idx, 2],
                                S['PMLb'], S['nPML'], S['vmaxPML'], S['R'])
    dampv_pml = np.zeros((pml_nodes_idx.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes_idx, :]
    a9 = 1.0 / dt - dampv_pml / 2.0
    b9 = 1.0 / dt + dampv_pml / 2.0

    # eleshp is a transposed VIEW, so dN[:,:,k] has a 24-byte inner stride --
    # one useful double per cache line in _c's multiply. Loop-invariant, so the
    # contiguous copy is paid once here. A copy is exact, and _c's reduction is
    # unaffected (dN*v allocates a contiguous temporary either way).
    lam_i = mat[E_int, 3]; miu_i = mat[E_int, 4]
    dN_i = eleshp[E_int]
    dNx_i = np.ascontiguousarray(dN_i[:, :, 0])
    dNy_i = np.ascontiguousarray(dN_i[:, :, 1])
    dNz_i = np.ascontiguousarray(dN_i[:, :, 2])
    constk_w_i = (-eledet[E_int]) * w
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel()
    idxIy = eq_ids[conn_i, 1].ravel()
    idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p = np.ascontiguousarray(dN_p[:, :, 0])
    dNy_p = np.ascontiguousarray(dN_p[:, :, 1])
    dNz_p = np.ascontiguousarray(dN_p[:, :, 2])
    det_w_p = eledet[E_pml] * w
    # f1..f9 all have the shape -det_w_p[:,None] * dN?_p * s[:,None], which
    # Python evaluates strictly left to right, so the left factor is loop-
    # invariant and hoisted. Same operands, same association, same roundings.
    wx_p = -det_w_p[:, None] * dNx_p
    wy_p = -det_w_p[:, None] * dNy_p
    wz_p = -det_w_p[:, None] * dNz_p
    neg_det_w_p = -det_w_p[:, None]
    conn_p = conn[E_pml]
    xc_p = S['meshCoor'][conn_p].mean(axis=1)
    d1p, d2p, d3p = region_damp(xc_p[:, 0], xc_p[:, 1], xc_p[:, 2], S['PMLb'],
                                S['nPML'], S['vmaxPML'], S['R'])
    a1 = 1.0 / dt - d1p / 2.0; b1 = 1.0 / dt + d1p / 2.0
    a2 = 1.0 / dt - d2p / 2.0; b2 = 1.0 / dt + d2p / 2.0
    a3 = 1.0 / dt - d3p / 2.0; b3 = 1.0 / dt + d3p / 2.0
    nd_p = ndof[conn_p]; is12_p = nd_p == 12; is3_p = nd_p == 3
    idxP12 = [np.where(is12_p, eq_ids[conn_p, jj], 0).ravel() for jj in range(12)]
    idxP3 = [np.where(is3_p, eq_ids[conn_p, gg], 0).ravel() for gg in range(3)]

    nd_all = ndof[conn]
    slot0 = np.where(nd_all == 3, 0, 9)
    slot1 = np.where(nd_all == 3, 1, 10)
    slot2 = np.where(nd_all == 3, 2, 11)
    idxH0 = np.take_along_axis(eq_ids[conn], slot0[:, :, None], axis=2)[:, :, 0].ravel()
    idxH1 = np.take_along_axis(eq_ids[conn], slot1[:, :, None], axis=2)[:, :, 0].ravel()
    idxH2 = np.take_along_axis(eq_ids[conn], slot2[:, :, None], axis=2)[:, :, 0].ravel()

    # --- C_elastic==0 terms; exact no-ops otherwise (see module docstring) ---
    C_elastic = S['C_elastic']
    grav_const = ((1.0 - C_elastic) * S['grav']
                  * (S['roumax'] - (S['gamar'] + 1.0) * S['rhow']) / S['roumax'])
    # contm's per-element lumped mass -- the uniform (non-wedge) closed
    # form ONLY (assembleGlobalMass.py's `contm`), not
    # `_contm_wedge_node_mass`'s per-node wedge value. This is a latent gap
    # for wedge elements, but a PROVEN no-op in this port's scope: `m_e_all`
    # is used only via `grav_const` below, which is EXACTLY 0.0 whenever
    # C_elastic==1 -- and eqdyna3d.py's build_solver_state already refuses
    # C_degen>3 (the only source of wedge elements) combined with
    # C_elastic==0, so a wedge run can never reach the C_elastic==0 branch
    # that would need this value to be right. If that refusal is ever
    # lifted, `m_e_all` must be recomputed per-node (elemType 11/12) here.
    m_e_all = mat[:, 2] * eledet
    if C_elastic == 0:
        init_stress = S['init_stress']   # raises loudly (KeyError) if missing;
        stress_i0 = init_stress[E_int].copy()    # never a silent zero.
        pml_init6 = init_stress[E_pml].copy()
        ccosphi = S['ccosphi']; sinphi = S['sinphi']; tv = S['tv']
    else:
        stress_i0 = np.zeros((E_int.shape[0], 6))
        pml_init6 = np.zeros((E_pml.shape[0], 6))
        ccosphi = sinphi = tv = 0.0

    return dict(
        N=S['N'], dt=dt, rdampk=S['rdampk'], NEQ=NEQ, NEQ1=NEQ + 1,
        conn=conn, phi=phi, ss=ss,
        int_nodes_idx=int_nodes_idx, idx3_v=idx3_v,
        pml_nodes_idx=pml_nodes_idx, idx12_v=idx12_v, a9=a9, b9=b9,
        Ei=E_int.shape[0], Ep=E_pml.shape[0], E=conn.shape[0],
        lam_i=lam_i, miu_i=miu_i, dNx_i=dNx_i, dNy_i=dNy_i, dNz_i=dNz_i,
        constk_w_i=constk_w_i, conn_i=conn_i,
        idxIx=idxIx, idxIy=idxIy, idxIz=idxIz,
        lam_p=lam_p, miu_p=miu_p, dNx_p=dNx_p, dNy_p=dNy_p, dNz_p=dNz_p,
        det_w_p=det_w_p, wx_p=wx_p, wy_p=wy_p, wz_p=wz_p,
        neg_det_w_p=neg_det_w_p, conn_p=conn_p,
        a1=a1, b1=b1, a2=a2, b2=b2, a3=a3, b3=b3,
        idxP12=idxP12, idxP3=idxP3, idxH0=idxH0, idxH1=idxH1, idxH2=idxH2,
        C_elastic=C_elastic, grav_const=grav_const,
        m_e_i=m_e_all[E_int], m_e_p=m_e_all[E_pml],
        stress_i0=stress_i0, pml_init6=pml_init6,
        ccosphi=ccosphi, sinphi=sinphi, tv=tv,
    )


def alloc_scratch(xp, inv):
    """Reusable element-force block buffers, or None on jax.

    The kernel builds ~30 blocks of shape (E, 8) per step and scatters each
    one the moment it is complete, so at most a handful are live at once and
    one buffer per role suffices. numpy needs them: allocating each block
    fresh is ~250 MB of mmap churn per step on test.tpv8 (see
    backend.mul_into). jax does not: arrays are immutable and XLA allocates,
    so returning None keeps these off the device entirely rather than
    shipping dead buffers through to_device.

    `stage_p` is (4, Ep, 8): three PML group accumulators plus one in-flight
    block. The in-flight buffer MUST be distinct from the accumulators --
    seeding a group by aliasing the block would let the next block overwrite
    the group. That is what store_into guards.
    """
    if B.is_jax(xp):
        return None
    return dict(stage=np.empty((inv['E'], 8)),
                stage_p=np.empty((4, inv['Ep'], 8)))


def assembleGlobalKU(xp, inv, velArr, force, stress_i, s_p, dt, rdampk, scratch):
    """assembleGlobalKU.f90 -- interior (calcElemKU) then PML
    (calcPMLElemKU) element forces, scattered block by block into force.

    force must already be zeroed by the caller (driver.f90:23). Returns
    (force, stress_i, s_p) -- the numpy backend mutates force in place and
    hands back the same buffer, the jax backend returns new arrays; callers
    must use the return value and assume neither.
    """
    conn_i = inv['conn_i']
    dNx_i, dNy_i, dNz_i = inv['dNx_i'], inv['dNy_i'], inv['dNz_i']
    lam_i, miu_i, constk_w_i = inv['lam_i'], inv['miu_i'], inv['constk_w_i']
    # Gather each velocity component separately. velArr[conn_i] would build an
    # (E,8,3) block whose [:,:,k] slices are 24-byte-strided and which _c then
    # reads nine times; this gathers the SAME values (a gather is exact)
    # straight into contiguous (E,8) arrays, measured faster on both ends
    # (each _c 6.40 -> 3.99 ms).
    #
    # This was briefly a per-backend predicate, on the theory that handing
    # XLA a differently SHAPED operand for the following reduction would
    # change its lowering. MEASURED: it does not. Switching the jax column to
    # the blocked-gather spelling left its divergence from the previous jax
    # port bit-for-bit unchanged (10 of 1,425,874 force entries, max 2.209e-14
    # at step 2 of test.tpv8, identical either way). A divergence that
    # measurement cannot detect is not a divergence, so the predicate is gone
    # and both backends use this one spelling.
    vlx = velArr[conn_i, 0]; vly = velArr[conn_i, 1]; vlz = velArr[conn_i, 2]
    sr1 = _c(xp, dNx_i, vlx); sr2 = _c(xp, dNy_i, vly); sr3 = _c(xp, dNz_i, vlz)
    sr4 = _c(xp, dNz_i, vly) + _c(xp, dNy_i, vlz)
    sr5 = _c(xp, dNz_i, vlx) + _c(xp, dNx_i, vlz)
    sr6 = _c(xp, dNy_i, vlx) + _c(xp, dNx_i, vly)
    vol_sr = sr1 + sr2 + sr3
    strr = (lam_i * vol_sr + 2 * miu_i * sr1,
            lam_i * vol_sr + 2 * miu_i * sr2,
            lam_i * vol_sr + 2 * miu_i * sr3,
            miu_i * sr4, miu_i * sr5, miu_i * sr6)

    for k in range(6):
        stress_i = B.addat(xp, stress_i, (slice(None), k), strr[k] * dt)

    if inv['C_elastic'] == 0:
        stress_i = _drucker_prager(xp, inv, stress_i, dt)

    st = [constk_w_i * (stress_i[:, k] + rdampk * strr[k]) for k in range(6)]

    # Fx, Fy, Fz -- each block is scattered the moment it is complete, so ONE
    # (Ei,8) buffer serves all three on numpy: each block's first product is an
    # `out=` overwrite, never a read-modify-write, so no part of Fx can survive
    # into Fy. `a + b + c` is left-associative, so seeding with the first
    # product and then `+=`-ing performs the same two additions on the same two
    # roundings as the functional form -- bit-identical, just not reallocated.
    buf = None if scratch is None else scratch['stage'][:conn_i.shape[0]]
    F = B.mul_into(xp, buf, dNx_i, st[0][:, None])
    F = B.iadd(xp, F, dNz_i * st[4][:, None])
    F = B.iadd(xp, F, dNy_i * st[5][:, None])
    force = B.scatter_add(xp, force, inv['idxIx'], F.ravel())
    F = B.mul_into(xp, buf, dNy_i, st[1][:, None])
    F = B.iadd(xp, F, dNz_i * st[3][:, None])
    F = B.iadd(xp, F, dNx_i * st[5][:, None])
    force = B.scatter_add(xp, force, inv['idxIy'], F.ravel())
    F = B.mul_into(xp, buf, dNz_i, st[2][:, None])
    F = B.iadd(xp, F, dNy_i * st[3][:, None])
    F = B.iadd(xp, F, dNx_i * st[4][:, None])
    # calcElemMass's gravity body force (al(3,:), z only) -- EXACTLY 0.0 when
    # C_elastic==1, so this is x - 0.0 == x, identical.
    F = B.iadd(xp, F, -(inv['m_e_i'] * inv['grav_const'])[:, None])
    force = B.scatter_add(xp, force, inv['idxIz'], F.ravel())

    if inv['Ep'] > 0:
        force, s_p = _pml(xp, inv, velArr, force, s_p, rdampk, scratch)

    return force, stress_i, s_p


def _drucker_prager(xp, inv, stress_i, dt):
    """calcElemKU.f90:127-161 -- viscoplastic return-mapping, interior
    elements only (PML elements never go through calcElemKU).

    porep is hardcoded 0.0: setPlasticStress and eqdyna3d.f90's zero-init
    are the ONLY writes to Fortran's eleporep anywhere in src/*.f90 -- the
    lithostatic pore-pressure formula is commented out in the Fortran
    itself, so it is provably 0.0, not assumed.

    pstrain is deliberately NOT computed: it is write-only (never read back
    by any downstream physics) and test.reference.results/test.drv.a6 has no
    pstr.txt to gate it against, so computing it would be unverifiable.
    Flagged, not silently included.
    """
    ccosphi = inv['ccosphi']; sinphi = inv['sinphi']; tv = inv['tv']
    strmea = (stress_i[:, 0] + stress_i[:, 1] + stress_i[:, 2]) / 3.0
    sd0 = stress_i[:, 0] - strmea
    sd1 = stress_i[:, 1] - strmea
    sd2 = stress_i[:, 2] - strmea
    sd3 = stress_i[:, 3]; sd4 = stress_i[:, 4]; sd5 = stress_i[:, 5]
    taomax = xp.sqrt(0.5 * (sd0 ** 2 + sd1 ** 2 + sd2 ** 2)
                     + sd3 ** 2 + sd4 ** 2 + sd5 ** 2)
    porep = 0.0
    yld = xp.maximum(ccosphi - sinphi * (strmea + porep), 0.0)
    mask = taomax > yld
    safe_tao = xp.where(mask, taomax, 1.0)    # avoid 0/0 on non-yielding elements
    ratio = yld / safe_tao
    # rjust == 1.0 where not yielding reproduces Fortran's "stress unchanged"
    # branch EXACTLY, so the assignments below need no second where.
    rjust = xp.where(mask, ratio + (1.0 - ratio) * xp.exp(-dt / tv), 1.0)
    for k, v in enumerate((sd0 * rjust + strmea, sd1 * rjust + strmea,
                           sd2 * rjust + strmea, sd3 * rjust,
                           sd4 * rjust, sd5 * rjust)):
        stress_i = B.setat(xp, stress_i, (slice(None), k), v)
    return stress_i


def _pml(xp, inv, velArr, force, s_p, rdampk, scratch):
    """calcPMLElemKU -- the 15 split stress components and the 12 f-blocks.

    Only the three group sums stay live; each of the twelve f-blocks is
    scattered the moment it is built. Accumulating a group incrementally
    (g = f1; g += f2; g += f3v; g += f10) is the SAME left-associative sum
    on the same two-at-a-time roundings as summing the four at once.
    """
    conn_p = inv['conn_p']
    dNx_p, dNy_p, dNz_p = inv['dNx_p'], inv['dNy_p'], inv['dNz_p']
    lam_p, miu_p = inv['lam_p'], inv['miu_p']
    a1, b1, a2, b2 = inv['a1'], inv['b1'], inv['a2'], inv['b2']
    a3, b3 = inv['a3'], inv['b3']

    vpx = velArr[conn_p, 0]; vpy = velArr[conn_p, 1]; vpz = velArr[conn_p, 2]
    # The nine velocity-gradient contractions are computed ONCE and the six
    # strain rates read off them: srp1 == Dxvx, srp2 == Dyvy, srp3 == Dzvz,
    # srp4 == Dzvy + Dyvz, srp5 == Dzvx + Dxvz, srp6 == Dyvx + Dxvy, operand
    # for operand -- so every srp is the same float64, only computed once.
    Dxvx = _c(xp, dNx_p, vpx); Dyvy = _c(xp, dNy_p, vpy); Dzvz = _c(xp, dNz_p, vpz)
    Dxvy = _c(xp, dNx_p, vpy); Dyvx = _c(xp, dNy_p, vpx)
    Dxvz = _c(xp, dNx_p, vpz); Dzvx = _c(xp, dNz_p, vpx)
    Dyvz = _c(xp, dNy_p, vpz); Dzvy = _c(xp, dNz_p, vpy)
    srp4 = Dzvy + Dyvz; srp5 = Dzvx + Dxvz; srp6 = Dyvx + Dxvy
    volp = Dxvx + Dyvy + Dzvz
    srate = [lam_p * volp + 2 * miu_p * Dxvx, lam_p * volp + 2 * miu_p * Dyvy,
             lam_p * volp + 2 * miu_p * Dzvz,
             miu_p * srp4, miu_p * srp5, miu_p * srp6]

    news = [((lam_p + 2 * miu_p) * Dxvx + a1 * s_p[:, 0]) / b1,
            (lam_p * Dyvy + a2 * s_p[:, 1]) / b2,
            (lam_p * Dzvz + a3 * s_p[:, 2]) / b3,
            (lam_p * Dxvx + a1 * s_p[:, 3]) / b1,
            ((lam_p + 2 * miu_p) * Dyvy + a2 * s_p[:, 4]) / b2,
            (lam_p * Dzvz + a3 * s_p[:, 5]) / b3,
            (lam_p * Dxvx + a1 * s_p[:, 6]) / b1,
            (lam_p * Dyvy + a2 * s_p[:, 7]) / b2,
            ((lam_p + 2 * miu_p) * Dzvz + a3 * s_p[:, 8]) / b3,
            (miu_p * Dxvy + a1 * s_p[:, 9]) / b1,
            (miu_p * Dyvx + a2 * s_p[:, 10]) / b2,
            (miu_p * Dxvz + a1 * s_p[:, 11]) / b1,
            (miu_p * Dzvx + a3 * s_p[:, 12]) / b3,
            (miu_p * Dyvz + a2 * s_p[:, 13]) / b2,
            (miu_p * Dzvy + a3 * s_p[:, 14]) / b3]
    for k in range(15):
        s_p = B.setat(xp, s_p, (slice(None), k), news[k])

    sxx = s_p[:, 0] + s_p[:, 1] + s_p[:, 2]
    syy = s_p[:, 3] + s_p[:, 4] + s_p[:, 5]
    szz = s_p[:, 6] + s_p[:, 7] + s_p[:, 8]
    sxy = s_p[:, 9] + s_p[:, 10]
    sxz = s_p[:, 11] + s_p[:, 12]
    syz = s_p[:, 13] + s_p[:, 14]
    # s0(1:6) == Fortran's s(16:21): setPlasticStress's lithostatic pre-stress,
    # only ever READ inside calcPMLElemKU. Exactly 0 for C_elastic==1.
    s0 = [rdampk * srate[k] + inv['pml_init6'][:, k] for k in range(6)]

    wx_p, wy_p, wz_p = inv['wx_p'], inv['wy_p'], inv['wz_p']
    ndw = inv['neg_det_w_p']
    ip12 = inv['idxP12']; ip3 = inv['idxP3']

    sp_ = (None, None, None, None) if scratch is None else scratch['stage_p']
    groups = [None, None, None]
    # f1..f9: block, scatter, then fold into its group. g = t aliases the
    # block on numpy, which is safe because the scatter has already consumed
    # it and nothing else references it.
    blocks9 = ((wx_p, sxx, 0), (wy_p, sxy, 0), (wz_p, sxz, 0),
               (wx_p, sxy, 1), (wy_p, syy, 1), (wz_p, syz, 1),
               (wx_p, sxz, 2), (wy_p, syz, 2), (wz_p, szz, 2))
    for jj, (wgt, sv, grp) in enumerate(blocks9):
        t = B.mul_into(xp, sp_[3], wgt, sv[:, None])
        force = B.scatter_add(xp, force, ip12[jj], t.ravel())
        # A group is SEEDED with a copy, not an alias: `t` is the shared
        # in-flight buffer on numpy, and the next block overwrites it.
        groups[grp] = (B.store_into(xp, sp_[grp], t) if groups[grp] is None
                       else B.iadd(xp, groups[grp], t))

    # f10..f12 multiply det_w_p into a SUM of three terms -- a different
    # association that wx_p/wy_p/wz_p cannot express, so they are rebuilt.
    blocks3 = ((dNx_p, 0, dNz_p, 4, dNy_p, 5, 0),
               (dNy_p, 1, dNz_p, 3, dNx_p, 5, 1),
               (dNz_p, 2, dNy_p, 3, dNx_p, 4, 2))
    for jj, (dA, kA, dB2, kB, dC, kC, grp) in enumerate(blocks3):
        t = B.mul_into(xp, sp_[3], dA, s0[kA][:, None])
        t = B.iadd(xp, t, dB2 * s0[kB][:, None])
        t = B.iadd(xp, t, dC * s0[kC][:, None])
        t = B.mul_into(xp, sp_[3], t, ndw)
        force = B.scatter_add(xp, force, ip12[9 + jj], t.ravel())
        groups[grp] = B.iadd(xp, groups[grp], t)

    # calcElemMass's gravity lands in the z (vhg_z) slot before calcPMLElemKU
    # subtracts from it; subtracted from the z group instead (addition is
    # commutative). EXACTLY 0.0 when C_elastic==1.
    groups[2] = B.iadd(xp, groups[2], -(inv['m_e_p'] * inv['grav_const'])[:, None])
    for g in range(3):
        force = B.scatter_add(xp, force, ip3[g], groups[g].ravel())
    return force, s_p


def calcHourglassResist(xp, inv, dispArr, velArr, force, rdampk):
    """calcHourglassResist.f90 -- KF78 hourglass control, C_hg==1, all
    elements.

    dispArr + rdampk*velArr is combined at NODE level and gathered once:
    gathering both (N,3) arrays up to (E,8,3) first would recompute each
    node's value once per element that touches it (8x over) on 141 MB arrays
    instead of 18 MB ones. A gather reproduces values exactly.

    The four modes' accumulation order is _fuse_modes's declared divergence
    and the contraction spelling is _contract_modes's -- see both.
    """
    conn = inv['conn']; phi = inv['phi']; ss = inv['ss']
    idxH = (inv['idxH0'], inv['idxH1'], inv['idxH2'])
    dl_all = (dispArr + rdampk * velArr)[conn]

    fuse = _fuse_modes(xp)
    acc = [None, None, None]
    for m in range(4):
        phid = _contract_modes(xp, phi[:, m, :], dl_all)
        r = (ss[:, 0] * phid[:, 0] + ss[:, 1] * phid[:, 1] + ss[:, 2] * phid[:, 2],
             ss[:, 1] * phid[:, 0] + ss[:, 3] * phid[:, 1] + ss[:, 4] * phid[:, 2],
             ss[:, 2] * phid[:, 0] + ss[:, 4] * phid[:, 1] + ss[:, 5] * phid[:, 2])
        for d in range(3):
            # -(phi*r), not a hoisted (-phi)*r. The two are identical in
            # scalar IEEE (negation only flips the sign bit), and an earlier
            # version of this kernel hoisted `neg_phi = -phi` so numpy could
            # write the product straight through an `out=` buffer with no
            # negation pass. That buffer is gone, so the hoist bought nothing
            # and cost an (E,4,8) array -- 188 MB on test.tpv104, on the host
            # AND, once promoted, on the device. Removing it is BIT-IDENTICAL
            # on numpy (verified, full-output digests at steps 1/5/114) and
            # restores the form the jax column was originally written with;
            # under XLA the hoisted form was measurably NOT equivalent.
            blk = -(phi[:, m, :] * r[d][:, None])
            if fuse:
                acc[d] = blk if acc[d] is None else acc[d] + blk
            else:
                force = B.scatter_add(xp, force, idxH[d], blk.ravel())
    if fuse:
        for d in range(3):
            force = B.scatter_add(xp, force, idxH[d], acc[d].ravel())

    # Scrub the sink: the scatter accumulated masked-out contributions there.
    return B.setat(xp, force, 0, 0.0)
