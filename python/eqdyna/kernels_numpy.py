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
    broadcast-multiply + axis-sum.

    NOT interchangeable with `np.einsum('ei,ei->e', ...)`, despite the name:
    the product here is a contiguous (E,8) temporary, so `.sum(axis=1)` runs
    NumPy's pairwise summation over the 8 terms, while einsum accumulates them
    sequentially. The two disagree in the last bits, and this port is gated on
    bit-identical output, so the einsum form is not an available substitution.
    """
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
    # PERF (2026-09-14): `eleshp` is a transposed view, so `dN_i[:,:,k]` is an
    # (E,8) array with a 24-byte inner stride -- one useful double per cache
    # line in `_c`'s multiply. These are loop-invariant, so the contiguous copy
    # is paid once at build time. Values are untouched (a copy is exact), and
    # `_c`'s reduction is unaffected: `dN*v` allocates a CONTIGUOUS temporary
    # either way, so `.sum(axis=1)` runs the same pairwise summation over the
    # same 8 floats in the same order -- bit-identical, verified end-to-end.
    dNx_i = np.ascontiguousarray(dN_i[:, :, 0])
    dNy_i = np.ascontiguousarray(dN_i[:, :, 1])
    dNz_i = np.ascontiguousarray(dN_i[:, :, 2])
    constk_w_i = (-eledet[E_int]) * w
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel(); idxIy = eq_ids[conn_i, 1].ravel(); idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p = np.ascontiguousarray(dN_p[:, :, 0])  # same one-time contiguous copy
    dNy_p = np.ascontiguousarray(dN_p[:, :, 1])  # as dNx_i/dNy_i/dNz_i above
    dNz_p = np.ascontiguousarray(dN_p[:, :, 2])
    det_w_p = eledet[E_pml] * w
    # PERF (2026-09-14): f1..f9 below all have the shape
    # `-det_w_p[:,None] * dN?_p * s[:,None]`, which Python evaluates strictly
    # left to right as `((-det_w_p[:,None]) * dN?_p) * s[:,None]`. The left
    # factor is loop-invariant and was being rebuilt nine times per step, so it
    # is hoisted here. Same operands, same association, same roundings -- the
    # products are bit-identical, only recomputed once instead of every step.
    # (f10..f12 keep their own form: they multiply det_w_p into a SUM of three
    # terms, a different association that these arrays cannot express.)
    wx_p = -det_w_p[:, None] * dNx_p
    wy_p = -det_w_p[:, None] * dNy_p
    wz_p = -det_w_p[:, None] * dNz_p
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

    # ---- PERF (2026-09-15): in-place block scatter, no staged function space ----
    # `elastic_step` accumulates every element force contribution into `force`.
    # Two earlier shapes of this, both replaced here:
    #   (a) build the index array AND the value array with a fresh
    #       `np.concatenate` every step, then one `np.bincount`;
    #   (b) (2026-09-14) hoist the concatenated index array to build() and
    #       stage the values in one persistent buffer written through with
    #       `out=`, still one `np.bincount`.
    # Both MATERIALISE THE WHOLE FUNCTION SPACE: 148 scatter entries per
    # element (interior force 3x8, hourglass 4 modes x 3 x 8, PML 15 x 8 on
    # PML elements only), i.e. on test.tpv104 (735,000 elements) 109,214,784
    # int64 indices (874 MB) plus 109,214,784 float64 values (874 MB) -- and
    # the index concatenation is pure DUPLICATION, because the individual
    # blocks it concatenates (idxIx/idxIy/idxIz, idxP12, idxP3, idxH0/1/2) are
    # kept anyway, and the four hourglass modes all reuse the same three.
    #
    # `np.add.at` accumulates into `force` IN PLACE, one block at a time, so
    # neither array is needed: no concatenated index, no staged value space.
    # The only value buffers that must exist are the ones that have to be LIVE
    # AT THE SAME TIME -- the three PML group sums (they accumulate across all
    # twelve PML f-blocks) plus one in-flight block -- because every other
    # block is consumed by its scatter the moment it is built. What is left is
    # `stage` (one (E,8) block, reused by the interior and hourglass blocks in
    # turn) and `stage_p` (4 x (E_pml,8)).
    #
    # MEASURED, full-length serial runs, peak RSS (/usr/bin/time -v's number):
    #     test.tpv8    (E=235008, 114 steps)   1.389 GB -> 0.855 GB  (-38%)
    #     test.tpv104  (E=735000, 120 steps)   3.968 GB -> 2.406 GB  (-39%)
    # and ms/step 543.40 -> 527.60 and 1691.42 -> 1660.95 respectively, i.e.
    # the footprint is not bought with time -- it is slightly cheaper too.
    #
    # BIT-IDENTITY -- the reason this is a legal substitution and not merely a
    # cheaper one. `np.bincount` accumulates in ARRAY order: block by block in
    # concatenation order and, within a block, in (element, node) order.
    # `np.add.at` accumulates in INDEX-ARRAY order. Calling it once per block,
    # in the SAME block order, therefore performs exactly the same additions
    # on exactly the same partial sums for every target equation. That
    # equivalence is load-bearing rather than cosmetic: float addition is not
    # associative, so ANY regrouping (per-block `np.bincount` calls, summing
    # the four hourglass modes before the scatter, reordering the blocks)
    # changes the answer in the last bits. Block order is therefore fixed:
    # interior Fx/Fy/Fz, then (PML elements only) f1..f12 and the three PML
    # group sums, then hourglass mode 0..3 x direction 0..2. `force` is zeroed
    # before the first block and `0.0 + x == x` exactly (including x == -0.0,
    # which bincount's own zero-initialised accumulator also renders +0.0), so
    # dropping bincount's separate output array changes nothing either.
    # Verified end-to-end, sha256 over tobytes() of the full velArr/dispArr/
    # force/fric/fnft at full step count: test.tpv8 (114 steps, 235,008
    # elements) and test.tpv104 (120 steps, 735,000 elements) digests
    # UNCHANGED against the bincount form.
    #
    # `np.add.at` used to be the slow path it is still remembered as (an
    # unbuffered per-element ufunc loop); NumPy 1.25's ufunc.at fast path
    # removed that, and on this box (NumPy 2.2.6) it is also FASTER than the
    # bincount form it replaces. Timed in ISOLATION on test.tpv8's real index
    # arrays (30 blocks, 37,552,896 entries, one core), best of 3:
    #     np.add.at, block at a time                    91.5 ms
    #     prebuilt index + staged values + bincount    136.8 ms
    #     concatenate index AND values + bincount      179.5 ms
    # -- and byte-identical output across all three. A previous measurement
    # reporting add.at 1.35x SLOWER does not reproduce on this NumPy; the
    # cross-check is testsys/unit/test_kernels_numpy_scatter.py, which pins the
    # ORDERING (not the speed), since that is what bit-identity depends on.
    #
    # THE ENTRY COUNT is the lever this does NOT pull, priced here so the next
    # reader does not have to re-derive it. 96 of the 148 entries per element
    # are hourglass (4 modes x 3 directions x 8 nodes), and all four modes
    # scatter through the SAME three index arrays, so summing the modes per
    # (element, node) first would take the whole scatter from 148 to 76 entries
    # per element. Timed the same way on test.tpv8: 91.5 ms -> 59.1 ms, i.e.
    # -32 ms of a 528 ms step (-6%), and now that nothing is staged it costs no
    # memory to do. It is still not done, because it is NOT bit-identical: for
    # a given equation the contributions currently arrive mode-major
    # (mode 0's eight elements, then mode 1's, ...) and fusing makes them
    # element-major, which is a different association of the same sum. The JAX
    # kernel DOES fuse the modes (kernels_jax.py's hourglass note, where the
    # accepted deviation is documented and gated against that port's own
    # nondeterminism floor); this one is the reference the JAX port is measured
    # against, so it keeps the Fortran-order association.
    stage = np.empty((conn.shape[0], 8))        # interior blocks use rows [:E_int]
    stage_p = np.empty((4, E_pml.shape[0], 8))  # PML group sums g0/g1/g2 + in-flight f

    # -(phi*r) == (-phi)*r exactly in IEEE 754 (negation only flips the sign
    # bit; multiplication's rounding is sign-symmetric), so folding the sign
    # into a hoisted copy lets the hourglass force write straight through `out=`
    # with no negation pass and no temporary. Bit-identical, not approximate.
    neg_phi = -phi
    neg_det_w_p = -det_w_p[:, None]
    # MEASURED AND DROPPED (2026-09-15), recorded so it is not re-attempted on
    # the strength of the same plausible-but-wrong hypothesis: `phi` as handed
    # in by main.py/loading.py is `np.transpose(phi48,(0,2,1))`, a VIEW, so the
    # per-mode einsum operand `phi[:, m, :]` has an inner stride of 4 doubles
    # -- one useful double per 32 bytes, exactly the pattern the
    # dNx_i/dNy_i/dNz_i contiguous copies above exist to avoid -- and that
    # einsum is 89.5 ms/step of test.tpv8's 528 (second only to the scatter, by
    # cProfile). Contracting against `neg_phi` instead (a fresh array, hence
    # C-contiguous) with the sign carried by a negated `ss` is bit-identical --
    # verified, full-length digests unchanged on test.tpv8 and test.tpv104,
    # since (-ss).(-phid) == +(ss.phid) term by term in IEEE -- and it is
    # WORTHLESS: 528.81 vs 527.91 ms/step on test.tpv8, 1656.30 vs 1660.95 on
    # test.tpv104, i.e. inside run-to-run spread, while costing the extra (E,6)
    # `neg_ss` array (+33 MB on test.tpv104). The stride is not what that
    # einsum is spending its time on; np.einsum's 'ei,eij->ej' inner loop is
    # scalar either way (~250 M mul-add/s here, far below this core's
    # bandwidth). Reverted. A real win there needs a different reduction, and
    # every faster reduction tried (stacked matmul, contiguous-axis reduce) is
    # a DIFFERENT accumulation order -- see the matmul/FMA note at the einsum
    # itself.

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
        det_w_p=det_w_p, wx_p=wx_p, wy_p=wy_p, wz_p=wz_p,
        conn_p=conn_p, a1=a1, b1=b1, a2=a2, b2=b2, a3=a3, b3=b3,
        idxP12=idxP12, idxP3=idxP3, idxH0=idxH0, idxH1=idxH1, idxH2=idxH2,
        stage=stage, stage_p=stage_p,
        neg_phi=neg_phi, neg_det_w_p=neg_det_w_p,
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
    # Read back ONCE, after the store, and reuse. The read-back must stay AFTER
    # the store and must not be replaced by the stored expression: `idx3_v`
    # repeats index 0 (main.py maps the -1 fixed-boundary sentinel onto the
    # sink slot), so for those entries the store is last-write-wins and the
    # value that comes back is NOT the value that went in. Preserving that is
    # the point -- this only removes a duplicate gather.
    v3 = v1[idx3_v]
    velArr[int_nodes] = v3
    dispArr[int_nodes] += v3 * dt

    if pml_nodes.size:
        f9 = force[idx12_v[:, 0:9]]
        vold9 = v1[idx12_v[:, 0:9]]
        v1[idx12_v[:, 0:9]] = (f9 + vold9 * a9) / b9
        f3v9 = force[idx12_v[:, 9:12]]
        v1[idx12_v[:, 9:12]] = v1[idx12_v[:, 9:12]] + f3v9 * dt
        v1[0] = 0.0
        # One (n_pml,12) gather instead of twelve (n_pml,) gathers of the same
        # `v1`. No store happens between them, so `V[:,k]` is element-for-
        # element the array `v1[idx12_v[:,k]]` used to return.
        V = v1[idx12_v]
        vA = V[:, 0] + V[:, 1] + V[:, 2] + V[:, 9]
        vB = V[:, 3] + V[:, 4] + V[:, 5] + V[:, 10]
        vC = V[:, 6] + V[:, 7] + V[:, 8] + V[:, 11]
        has_eq = idx12_v[:, 0] > 0
        velArr[pml_nodes, 0] = np.where(has_eq, vA, 0.0)
        velArr[pml_nodes, 1] = np.where(has_eq, vB, 0.0)
        velArr[pml_nodes, 2] = np.where(has_eq, vC, 0.0)
        dispArr[pml_nodes, 0] = np.where(has_eq, dispArr[pml_nodes, 0] + vA * dt, 0.0)
        dispArr[pml_nodes, 1] = np.where(has_eq, dispArr[pml_nodes, 1] + vB * dt, 0.0)
        dispArr[pml_nodes, 2] = np.where(has_eq, dispArr[pml_nodes, 2] + vC * dt, 0.0)

    # `force` is accumulated BLOCK BY BLOCK from here on, with `np.add.at` and
    # in the fixed block order build() documents -- zeroed first, so the first
    # contribution to any equation lands on an exact 0.0 exactly as
    # `np.bincount`'s own accumulator did. `blk` is the single reusable
    # (E,8) value buffer; nothing stages the whole 148-entries-per-element
    # function space any more (see build()).
    force[:] = 0.0
    blk = inv['stage']

    # ---- assembleGlobalKU: interior elements ----
    conn_i = inv['conn_i']; dNx_i, dNy_i, dNz_i = inv['dNx_i'], inv['dNy_i'], inv['dNz_i']
    lam_i, miu_i, constk_w_i = inv['lam_i'], inv['miu_i'], inv['constk_w_i']
    # PERF (2026-09-14): gather each velocity component separately. `velArr[
    # conn_i]` builds an (E,8,3) block whose `[:,:,k]` slices are 24-byte-
    # strided, which `_c` then reads nine times; `velArr[conn_i, k]` gathers
    # the SAME values (a gather with a constant last index -- verified equal
    # element-for-element) straight into a contiguous (E,8) array. Measured
    # faster on BOTH ends: the three gathers cost less than the one blocked
    # gather, and each `_c` drops from 6.40 ms to 3.99 ms.
    vlx = velArr[conn_i, 0]; vly = velArr[conn_i, 1]; vlz = velArr[conn_i, 2]
    sr1 = _c(dNx_i, vlx); sr2 = _c(dNy_i, vly); sr3 = _c(dNz_i, vlz)
    sr4 = _c(dNz_i, vly) + _c(dNy_i, vlz)
    sr5 = _c(dNz_i, vlx) + _c(dNx_i, vlz)
    sr6 = _c(dNy_i, vlx) + _c(dNx_i, vly)
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
    # `a + b + c` is left-associative, so `x = a; x += b; x += c` performs the
    # same two additions on the same two roundings -- the only thing that
    # changes is that the sum lands in the reusable block buffer.
    # Each block is scattered the moment it is complete, so ONE (E_int,8)
    # buffer serves all three: Fy's store is an `out=` overwrite, never a
    # read-modify-write, so no part of Fx can survive into it.
    F = blk[:conn_i.shape[0]]
    np.multiply(dNx_i, st1[:, None], out=F); F += dNz_i * st5[:, None]; F += dNy_i * st6[:, None]
    np.add.at(force, inv['idxIx'], F.ravel())
    np.multiply(dNy_i, st2[:, None], out=F); F += dNz_i * st4[:, None]; F += dNx_i * st6[:, None]
    np.add.at(force, inv['idxIy'], F.ravel())
    np.multiply(dNz_i, st3[:, None], out=F); F += dNy_i * st4[:, None]; F += dNx_i * st5[:, None]
    # calcElemMass's gravity body force (al(3,:), z-direction only) -- see
    # build()'s docstring for the -m_e*grav_const derivation; EXACTLY 0.0
    # when C_elastic==1.
    F -= (inv['m_e_i'] * inv['grav_const'])[:, None]
    np.add.at(force, inv['idxIz'], F.ravel())

    # ---- assembleGlobalKU: PML elements ----
    E_pml = inv['E_pml']
    if E_pml.shape[0]:
        conn_p = inv['conn_p']; dNx_p, dNy_p, dNz_p = inv['dNx_p'], inv['dNy_p'], inv['dNz_p']
        lam_p, miu_p, det_w_p = inv['lam_p'], inv['miu_p'], inv['det_w_p']
        a1, b1, a2, b2, a3, b3 = inv['a1'], inv['b1'], inv['a2'], inv['b2'], inv['a3'], inv['b3']
        vpx = velArr[conn_p, 0]; vpy = velArr[conn_p, 1]; vpz = velArr[conn_p, 2]
        # PERF (2026-09-14): the nine velocity-gradient contractions are
        # computed ONCE and the six strain rates are read off them. This is
        # exact re-use, not an approximation: the previous code called `_c`
        # 15 times, and six of those calls had operands character-for-character
        # identical to one of the nine D-calls below
        #   srp1 == Dxvx, srp2 == Dyvy, srp3 == Dzvz,
        #   srp4 == Dzvy + Dyvz, srp5 == Dzvx + Dxvz, srp6 == Dyvx + Dxvy
        # so every srp below is the SAME float64 as before, bit for bit; only
        # the duplicate evaluation is gone.
        Dxvx = _c(dNx_p, vpx); Dyvy = _c(dNy_p, vpy); Dzvz = _c(dNz_p, vpz)
        Dxvy = _c(dNx_p, vpy); Dyvx = _c(dNy_p, vpx)
        Dxvz = _c(dNx_p, vpz); Dzvx = _c(dNz_p, vpx)
        Dyvz = _c(dNy_p, vpz); Dzvy = _c(dNz_p, vpy)

        srp1 = Dxvx; srp2 = Dyvy; srp3 = Dzvz
        srp4 = Dzvy + Dyvz
        srp5 = Dzvx + Dxvz
        srp6 = Dyvx + Dxvy
        volp = srp1 + srp2 + srp3
        srate = [lam_p * volp + 2 * miu_p * srp1, lam_p * volp + 2 * miu_p * srp2,
                 lam_p * volp + 2 * miu_p * srp3, miu_p * srp4, miu_p * srp5, miu_p * srp6]

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

        wx_p, wy_p, wz_p = inv['wx_p'], inv['wy_p'], inv['wz_p']  # == -det_w_p[:,None]*dN?_p, hoisted
        # The twelve PML f-blocks are built ONE AT A TIME in `t` and scattered
        # immediately; only the three group sums g0/g1/g2 have to stay live,
        # because each accumulates four of the f-blocks (f1+f2+f3v+f10,
        # f4+f5+f6+f11, f7+f8+f9v+f12). Accumulating them incrementally --
        # `g = f1; g += f2; g += f3v; g += f10` -- is the SAME left-associative
        # sum `np.add(f1,f2,out=g); g += f3v; g += f10` performed, on the same
        # two-at-a-time roundings, so the group sums are bit-identical while
        # only 4 blocks of (E_pml,8) exist instead of 15.
        sp_ = inv['stage_p']; g0 = sp_[0]; g1 = sp_[1]; g2 = sp_[2]; t = sp_[3]
        ip12 = inv['idxP12']; ip3 = inv['idxP3']
        np.multiply(wx_p, sxx[:, None], out=t); np.copyto(g0, t); np.add.at(force, ip12[0], t.ravel())
        np.multiply(wy_p, sxy[:, None], out=t); g0 += t; np.add.at(force, ip12[1], t.ravel())
        np.multiply(wz_p, sxz[:, None], out=t); g0 += t; np.add.at(force, ip12[2], t.ravel())
        np.multiply(wx_p, sxy[:, None], out=t); np.copyto(g1, t); np.add.at(force, ip12[3], t.ravel())
        np.multiply(wy_p, syy[:, None], out=t); g1 += t; np.add.at(force, ip12[4], t.ravel())
        np.multiply(wz_p, syz[:, None], out=t); g1 += t; np.add.at(force, ip12[5], t.ravel())
        np.multiply(wx_p, sxz[:, None], out=t); np.copyto(g2, t); np.add.at(force, ip12[6], t.ravel())
        np.multiply(wy_p, syz[:, None], out=t); g2 += t; np.add.at(force, ip12[7], t.ravel())
        np.multiply(wz_p, szz[:, None], out=t); g2 += t; np.add.at(force, ip12[8], t.ravel())
        # f10..f12 multiply det_w_p into a SUM of three terms, so they cannot
        # reuse wx_p/wy_p/wz_p. `t *= neg_det_w_p` instead of `neg_det_w_p * t`
        # is the same product -- IEEE multiplication is commutative exactly.
        ndw = inv['neg_det_w_p']
        np.multiply(dNx_p, s0[0][:, None], out=t); t += dNz_p * s0[4][:, None]
        t += dNy_p * s0[5][:, None]; t *= ndw
        g0 += t; np.add.at(force, ip12[9], t.ravel())
        np.multiply(dNy_p, s0[1][:, None], out=t); t += dNz_p * s0[3][:, None]
        t += dNx_p * s0[5][:, None]; t *= ndw
        g1 += t; np.add.at(force, ip12[10], t.ravel())
        np.multiply(dNz_p, s0[2][:, None], out=t); t += dNy_p * s0[3][:, None]
        t += dNx_p * s0[4][:, None]; t *= ndw
        g2 += t; np.add.at(force, ip12[11], t.ravel())

        # calcElemMass's gravity body force lands in the z ("vhg_z") slot
        # (efPML12[11]) BEFORE calcPMLElemKU subtracts from it -- subtracted
        # here from the z group sum instead (addition is commutative);
        # EXACTLY 0.0 when C_elastic==1.
        g2 -= (inv['m_e_p'] * inv['grav_const'])[:, None]
        np.add.at(force, ip3[0], g0.ravel())
        np.add.at(force, ip3[1], g1.ravel())
        np.add.at(force, ip3[2], g2.ravel())

    # ---- hrglss (C_hg==1, all elements) ----
    conn = inv['conn']; phi = inv['phi']; ss = inv['ss']
    neg_phi = inv['neg_phi']
    idxH0 = inv['idxH0']; idxH1 = inv['idxH1']; idxH2 = inv['idxH2']
    # PERF (2026-09-14): combine at NODE level, then gather once. The old form
    # gathered both (N,3) arrays up to (E,8,3) -- 45 MB each here -- and did the
    # arithmetic on the expanded copies, recomputing the same node value once
    # per element that touches the node (8x over). `dispArr + rdampk*velArr` is
    # the identical expression evaluated on the 6 MB node arrays instead, and a
    # gather reproduces values exactly, so `dl_all` is bit-identical.
    dl_all = (dispArr + rdampk * velArr)[conn]
    # PERF (2026-09-14): each mode's contraction against dl_all is ONE
    # `np.einsum('ei,eij->ej', ...)` -- a reduce over the 8 nodes -- rather than
    # `(phi[:,m,:,None]*dl_all).sum(axis=1)`, which materialised a full (E,8,3)
    # product temporary (45 MB here) only to reduce it away. NOT replaced by a
    # stacked (E,4,8)@(E,8,3) `matmul`: matmul is free to contract with FMA,
    # which changes the answer in the last bits (an FMA contraction inside a
    # matmul is exactly how an earlier change here looked bit-identical across
    # 2.8M values and then diverged at step 5).
    # `phi[:, m, :]` is a STRIDED slice (phi is a transposed view); switching
    # the operand to the contiguous `neg_phi` was tried, is bit-identical, and
    # buys nothing -- see build()'s "MEASURED AND DROPPED" note.
    for m in range(4):
        phid = np.einsum('ei,eij->ej', phi[:, m, :], dl_all)
        r0 = ss[:, 0] * phid[:, 0] + ss[:, 1] * phid[:, 1] + ss[:, 2] * phid[:, 2]
        r1 = ss[:, 1] * phid[:, 0] + ss[:, 3] * phid[:, 1] + ss[:, 4] * phid[:, 2]
        r2 = ss[:, 2] * phid[:, 0] + ss[:, 4] * phid[:, 1] + ss[:, 5] * phid[:, 2]
        # Each mode's three direction blocks are scattered as they are built,
        # through the same reusable buffer. The four modes are NOT summed per
        # (element, node) first: that would be a cheaper scatter (24 entries
        # per element instead of 96) but a DIFFERENT association -- see the
        # measured tradeoff in build()'s note -- and this path is gated on
        # bit-identity.
        np.multiply(neg_phi[:, m, :], r0[:, None], out=blk); np.add.at(force, idxH0, blk.ravel())
        np.multiply(neg_phi[:, m, :], r1[:, None], out=blk); np.add.at(force, idxH1, blk.ravel())
        np.multiply(neg_phi[:, m, :], r2[:, None], out=blk); np.add.at(force, idxH2, blk.ravel())

    force[0] = 0.0  # scrub sink again (the scatter accumulated masked-out contributions there)

    return v1, velArr, dispArr, force, stress_i, s_p
