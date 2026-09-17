"""
Milestone 6.5: standalone port of EQdyna's lumped-mass assembly --
src/assembleGlobalMass.f90's element loop (calcGlobalShapeFunc's Jacobian-
determinant computation only, NOT its wedge-degeneration branch or the
full spatial-derivative/hourglass output -- see scope note below) +
contm (per-element lumped mass) + assembleElementMassDetShg's
nodalMassArr/fnms scatter.

SCOPE: this milestone ports `nodalMassArr` (the equation-numbered lumped
mass vector the time-stepper divides nodal force by) and `fnms` (the
per-node total lumped mass, keyed by node id, used by faulting.f90's
split-node mass pairing) -- named "lumped-mass integration" by the
coordinator, and nothing else assembleGlobalMass.f90 touches. In
particular:
  - `eledet`/`eleshp` (calcGlobalShapeFunc's full output: spatial shape-
    function derivatives + the transformed Jacobian-inverse `xs`) and
    `ss`/`phi` (calcSSPhi4Hrgls's hourglass-control tensors) are NOT ported
    here -- out of scope for mass assembly, needed instead by whichever
    future milestone ports assembleGlobalKU (the FEM stiffness kernel).
    `pydump_elemgeo.txt` remains available as that milestone's oracle.
  - MPI4NodalQuant's neighbor-exchange is entirely guarded behind
    `if (npxyz(ixyz)>1)` for each of the 3 axes -- for the serial
    (npx=npy=npz=1) case this milestone targets, confirmed (by reading the
    guard directly) to be a no-op beyond an MPI barrier/timer; nothing to
    port.
  - C_degen==0 (no wedge degeneration): elemTypeArr never takes the 11/12
    values calcGlobalShapeFunc's degeneration branch checks for in this
    milestone's scope (matching every other meshgen milestone's scope
    note) -- `globalShapeFunc == localShapeFunc` going into the Jacobian
    for every element, exactly Fortran's own unconditional
    `globalShapeFunc = localShapeFunc` copy for this case.
    STILL TRUE as of the C_degen>3 (tpv36/tpv37) mesh port (meshgen.py's
    build_elements/library_degeneration.f90 wedge()/reorder()): this
    module's `compute_element_shape`/`compute_hourglass`/`contm` were NOT
    extended with calcGlobalShapeFunc.f90:22-28's elemTypeArr==11/12
    shape-function-merge branch -- left explicitly REFUSING, not silently
    wrong: eqdyna3d.py's `build_solver_state` raises NotImplementedError
    before calling into this module whenever the mesh contains any
    elemTypeArr 11/12 element, rather than feeding them through this
    module's generic (non-degenerate) formulas.

Verified against `testsys/parity/fixtures/test_tpv8_serial/
pydump_nodalmass.txt` / `pydump_fnms.txt` via
the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15) -- these fixtures are the
PARITY ORACLE ONLY. Per the explicit provenance decision recorded for M8
(meshgen.py's top docstring): the standalone solver
path this milestone eventually feeds must call THIS port at runtime, and
must never read pydump_nodalmass.txt/pydump_fnms.txt in production --
those files only exist because a Fortran binary with pydump instrumentation
was run to generate them, which the standalone path is explicitly meant to
not require.

Milestone 7.5 (closes the M6.5 scope gap flagged above): adds
`compute_element_shape` (calcGlobalShapeFunc.f90's FULL output -- the
Jacobian-inverse-transformed spatial shape-function derivatives `eleshp`
plus the cofactor/det matrix `xs`, for BOTH interior (elemTypeArr==1) and
PML (elemTypeArr==2) elements -- the C_degen==0 scope this whole module
targets means neither type ever takes the wedge-degeneration (elemTypeArr
11/12) branch, so one code path serves both element types, exactly as the
Fortran's single `calcGlobalShapeFunc` call from `assembleGlobalMass`'s
element loop does), `compute_hourglass` (`calcSSPhi4Hrgls`'s `ss`/`phi`
hourglass-control tensors, including a port of `func_lib.f90`'s `vlm`
Belytschko-formula element volume), and `init_vel` (`eqdyna3d.f90`'s
fault-node velocity seeding of the 1D solver-state vector `v1` from
`fric(31:36)` via the equation-number map -- mode==1 (non-restart) only,
matching this port's case scope; mode==2 restart reads `fric(31:36)` from
a restart netCDF this port does not read, out of scope).

Verified against `testsys/parity/fixtures/test_tpv8_serial/
pydump_elemgeo.txt` (eledet/eleshp/ss/phi per element) and
`pydump_v1.txt` (the full v1 array) via
the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15) -- ORACLE ONLY, same runtime-vs-
-oracle rule as above: the standalone solver must call these ports, never
read the pydump files.
"""
import numpy as np

# calcLocalShapeFunc.f90's single reduced-Gauss-point constants, verbatim.
_CST = 1.0 / 8.0
_ACOOR = np.array([
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
]).T  # (3,8), acoor(j,i) in Fortran
_LOCAL_DERIV = _CST * _ACOOR  # (3,8): localShapeFunc(1:3,:)
_W = 8.0  # calcLocalShapeFunc's `w = 8.0d0`, the 1-point Gaussian weight


def compute_element_det(xl):
    """Port of calcGlobalShapeFunc's Jacobian-determinant computation
    ONLY (C_degen==0 branch, no wedge degeneration -- see module docstring).

    xl: (E,8,3) physical coordinates of each element's 8 corner nodes, in
    Fortran nodeElemIdRelation column order (conn from build_elements).

    Returns det: (E,) float array, the Jacobian determinant at the single
    reduced-Gauss point. Raises ValueError on any non-positive determinant
    (Fortran's own `if (det <= 0.0d0) ... stop`, made a loud Python
    exception instead of a Fortran STOP).
    """
    # xs[e,j,i] = sum_k localDeriv[i,k] * xl[e,k,j]  (Fortran's xs(j,i))
    xs = np.einsum('ik,ekj->eji', _LOCAL_DERIV, xl)

    cof11 = xs[:, 1, 1] * xs[:, 2, 2] - xs[:, 1, 2] * xs[:, 2, 1]
    cof12 = xs[:, 1, 2] * xs[:, 2, 0] - xs[:, 1, 0] * xs[:, 2, 2]
    cof13 = xs[:, 1, 0] * xs[:, 2, 1] - xs[:, 1, 1] * xs[:, 2, 0]
    det = xs[:, 0, 0] * cof11 + xs[:, 0, 1] * cof12 + xs[:, 0, 2] * cof13

    bad = np.nonzero(det <= 0.0)[0]
    if bad.size:
        raise ValueError('compute_element_det: non-positive determinant at '
                          'element(s) %r (det=%r)' % (bad[:5].tolist(), det[bad[:5]].tolist()))
    return det


def contm(det, constm):
    """Port of contm's lumped-mass formula (C_degen==0, single reduced-
    Gauss-point, homogeneous 8-node brick).

    Closed form, not a re-implementation of Fortran's 8-way accumulation
    loop: contm's per-node shape-function VALUE (row 4 of globalShapeFunc)
    is `localShapeFunc(4,j) = 1/8` for every node j, in every element,
    unconditionally (calcLocalShapeFunc.f90's hardcoded constant, never
    recomputed per element) -- so `temp2 = totmas*(1/8)**2` is IDENTICAL
    for all 8 nodes, `dsum = 8*temp2`, `temp1 = totmas/dsum = 8.0`
    (an exact algebraic identity given the hardcoded 1/8, not an
    approximation), and every node's lumped mass is
    `temp1*temp2 = totmas/8 = constm*w*det/8 = constm*det` (using
    calcLocalShapeFunc's `w=8.0d0`). Verified against the fixture below at
    the same tolerance every other milestone in this port uses for
    multi-element-accumulated quantities.

    Returns m_e: (E,) float array, the (identical-for-x/y/z, identical-
    for-all-8-nodes) per-element-per-node lumped mass contribution.
    """
    return constm * det


def assemble_mass(conn, mat, det, num_dof, eq_start, eq_nums, n_equations, n_nodes):
    """Port of assembleElementMassDetShg's nodalMassArr/fnms scatter
    (mass-only slice -- eledet/eleshp writes are out of scope, see module
    docstring), in the verbatim element-major/local-node-minor accumulation
    ORDER of the Fortran (and of this function's original scalar loop).

    PERFORMANCE NOTE (order-preserving vectorization, replacing the original
    triple scalar loop -- the scalar version is kept below as
    `_assemble_mass_scalar` and is the bit-for-bit oracle the unit test
    `testsys/unit/test_assembleGlobalMass_vectorized.py` checks this against):
    the original docstring's constraint -- "NOT vectorizing a scatter whose
    accumulation order affects bit-for-bit output" -- is respected here
    exactly, NOT relaxed. `np.add.at`/`np.bincount` are deliberately NOT
    used: both would reassociate the per-target summation. Instead:

      1. Every (element, local-node) contribution is bucketed by its target
         node id with a STABLE sort, so within each node the contributions
         stay in increasing `e*8 + k` order -- the same order the scalar
         loop visits them.
      2. Those buckets are laid out as a dense (n_nodes+1, kmax) matrix,
         zero-padded AT THE TAIL, and summed by `kmax` sequential
         whole-column adds. Column j adds exactly the j-th contribution of
         every node, so each node sees the identical left-to-right addition
         chain the scalar loop performs; the tail padding contributes
         `x + 0.0`, which is exact in IEEE-754 for every x this function
         can produce (all contributions are strictly positive -- `constm`
         and `det` are both > 0, det checked by compute_element_det).
      3. `nodalMassArr` is then NOT accumulated separately: the scalar loop
         adds the SAME value, in the SAME order, to every positive equation
         slot of a node as it adds to `fnms[node]`, so
         `nodalMassArr[eq] == fnms[node]` bit-for-bit by construction, for
         every eq>0 in eq_nums[node]. This is a re-derivation of the scalar
         loop's own invariant, not an approximation of it.

    Both the numOfDofPerNodeArr==12 (PML) and ==3 (interior/exterior)
    Fortran branches turn out to add the SAME scalar `m_e[e]` to every
    valid (`eq>0`) equation-number slot for that (element,node) pair --
    confirmed by tracing contm's output (all 8 nodes' 3 components equal)
    through both branches, not assumed -- so this port uses one unified
    loop instead of reproducing the two syntactically-different Fortran
    branches; the two are provably identical in the VALUES they add, only
    Fortran's 12-dof branch has the `eq>0` guard the 3-dof branch omits
    (safe by construction there: 3-dof nodes are never at the model's
    fixed outer boundary, see readInputFiles.py/meshgen.py's PML-vs-interior
    classification) -- this port applies the `eq>0` guard universally,
    which is a no-op for 3-dof nodes and correct for 12-dof nodes.

    conn: (E,8) 1-indexed node ids (from build_elements).
    mat: (E,5) [vp,vs,rho,lambda,mu] (from build_elements); density is
        column 2 (0-indexed), matching Fortran's `mat(nel,3)`.
    det: (E,) from compute_element_det.
    num_dof, eq_start, eq_nums: from build_equation_numbers (1-indexed by
        node id; eq_nums[node] is that node's array of equation numbers,
        positive or -1 for the fixed-boundary sentinel).
    n_equations: totalNumOfEquations (nodalMassArr size).
    n_nodes: totalNumOfNodes (fnms size).

    Returns (nodalMassArr, fnms): 1-indexed float arrays (row 0 unused),
    sizes (n_equations+1,) and (n_nodes+1,).
    """
    m_e = contm(det, mat[:, 2])

    # (1) one contribution per (element, local node), in the scalar loop's
    # own visit order: flat index p = e*8 + k.
    flat = np.ascontiguousarray(conn).ravel()
    if flat.size and (flat.min() < 1 or flat.max() > n_nodes):
        raise ValueError('assemble_mass: conn references node id(s) outside '
                          '1..%d (min=%d, max=%d)'
                          % (n_nodes, int(flat.min()), int(flat.max())))
    vals = np.repeat(m_e, 8)

    # (2) stable bucket-by-node, then dense (node, occurrence) layout.
    order = np.argsort(flat, kind='stable')
    sorted_nodes = flat[order]
    counts = np.bincount(flat, minlength=n_nodes + 1)
    kmax = int(counts.max()) if counts.size else 0
    starts = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(counts[:-1], out=starts[1:])
    rank = np.arange(sorted_nodes.size, dtype=np.int64) - starts[sorted_nodes]

    contrib = np.zeros((n_nodes + 1, kmax))
    contrib[sorted_nodes, rank] = vals[order]

    # (3) sequential left-to-right accumulation, one column at a time.
    fnms = np.zeros(n_nodes + 1)
    for j in range(kmax):
        fnms += contrib[:, j]

    # (4) nodalMassArr[eq] == fnms[node] for every positive eq of that node
    # (identical value, identical order -- see docstring).
    nd = np.asarray(num_dof[1:], dtype=np.int64)
    eq_flat = np.concatenate(eq_nums[1:]).astype(np.int64, copy=False)
    if eq_flat.size != int(nd.sum()):
        raise ValueError('assemble_mass: eq_nums slot count %d disagrees with '
                          'num_dof total %d' % (eq_flat.size, int(nd.sum())))
    node_of_slot = np.repeat(np.arange(1, n_nodes + 1, dtype=np.int64), nd)
    live = eq_flat > 0
    nodalMassArr = np.zeros(n_equations + 1)
    nodalMassArr[eq_flat[live]] = fnms[node_of_slot[live]]

    return nodalMassArr, fnms


def _assemble_mass_scalar(conn, mat, det, num_dof, eq_start, eq_nums, n_equations, n_nodes):
    """The original verbatim scalar port of assembleElementMassDetShg's
    scatter, kept as the bit-for-bit ORACLE for `assemble_mass`'s
    order-preserving vectorization (see that function's docstring). Not
    used on any production path -- `testsys/unit/test_assembleGlobalMass_vectorized.py`
    asserts byte-equality of the two on a real mesh.
    """
    E = conn.shape[0]
    m_e = contm(det, mat[:, 2])
    nodalMassArr = np.zeros(n_equations + 1)
    fnms = np.zeros(n_nodes + 1)

    for e in range(E):
        me_val = m_e[e]
        for k in range(8):
            nodeID = int(conn[e, k])
            for eq in eq_nums[nodeID]:
                if eq > 0:
                    nodalMassArr[eq] += me_val
            fnms[nodeID] += me_val

    return nodalMassArr, fnms


def compute_element_shape(xl):
    """Milestone 7.5: port of calcGlobalShapeFunc's FULL output (C_degen==0,
    no wedge-degeneration branch -- see module docstring), for both interior
    and PML elements (both take the identical code path in this scope).

    xl: (E,8,3) physical coordinates of each element's 8 corner nodes, in
    Fortran nodeElemIdRelation column order (conn from build_elements) --
    same convention as compute_element_det.

    Returns (det, eleshp, xs):
      det: (E,) Jacobian determinant (bit-identical to compute_element_det's
        return -- same formula, kept as a separate function per that
        function's own docstring/scope, not recomputed differently here).
      eleshp: (E,3,8) spatial shape-function derivatives dNx/dNy/dNz per
        node (Fortran's eleshp(1:nrowsh-1,1:nen,elemID) -- row 4, the
        shape-function VALUE, is never transformed by this subroutine, per
        the Fortran's own tmpGlobalShapeFunc(4,i) never being reassigned,
        and is not part of eleshp -- contm's closed form already accounts
        for it, see compute_element_det/contm's docstrings).
      xs: (E,3,3) the cofactor-matrix-over-determinant Fortran RETURNS in
        its own `xs` argument (reused variable name, NOT the same `xs` as
        the Jacobian computed mid-subroutine) -- this is what
        calcSSPhi4Hrgls's `xs` parameter actually receives.
    """
    xs_jac = np.einsum('ik,ekj->eji', _LOCAL_DERIV, xl)  # xs(j,i) mid-subroutine

    cof11 = xs_jac[:, 1, 1] * xs_jac[:, 2, 2] - xs_jac[:, 1, 2] * xs_jac[:, 2, 1]
    cof12 = xs_jac[:, 1, 2] * xs_jac[:, 2, 0] - xs_jac[:, 1, 0] * xs_jac[:, 2, 2]
    cof13 = xs_jac[:, 1, 0] * xs_jac[:, 2, 1] - xs_jac[:, 1, 1] * xs_jac[:, 2, 0]
    cof21 = xs_jac[:, 2, 1] * xs_jac[:, 0, 2] - xs_jac[:, 2, 2] * xs_jac[:, 0, 1]
    cof22 = xs_jac[:, 2, 2] * xs_jac[:, 0, 0] - xs_jac[:, 2, 0] * xs_jac[:, 0, 2]
    cof23 = xs_jac[:, 2, 0] * xs_jac[:, 0, 1] - xs_jac[:, 2, 1] * xs_jac[:, 0, 0]
    cof31 = xs_jac[:, 0, 1] * xs_jac[:, 1, 2] - xs_jac[:, 0, 2] * xs_jac[:, 1, 1]
    cof32 = xs_jac[:, 0, 2] * xs_jac[:, 1, 0] - xs_jac[:, 0, 0] * xs_jac[:, 1, 2]
    cof33 = xs_jac[:, 0, 0] * xs_jac[:, 1, 1] - xs_jac[:, 0, 1] * xs_jac[:, 1, 0]

    det = xs_jac[:, 0, 0] * cof11 + xs_jac[:, 0, 1] * cof12 + xs_jac[:, 0, 2] * cof13
    bad = np.nonzero(det <= 0.0)[0]
    if bad.size:
        raise ValueError('compute_element_shape: non-positive determinant at '
                          'element(s) %r (det=%r)' % (bad[:5].tolist(), det[bad[:5]].tolist()))

    # eleshp(row,node) = (tmpG(1,node)*cof_row1 + tmpG(2,node)*cof_row2
    #                      + tmpG(3,node)*cof_row3) / det, tmpG == localDeriv
    # (globalShapeFunc == localShapeFunc going in, C_degen==0, no
    # degeneration overwrite).
    cof_rows = np.stack([
        np.stack([cof11, cof12, cof13], axis=1),
        np.stack([cof21, cof22, cof23], axis=1),
        np.stack([cof31, cof32, cof33], axis=1),
    ], axis=1)  # (E,3,3): cof_rows[e,row,col]
    eleshp = np.einsum('erc,ck->erk', cof_rows, _LOCAL_DERIV) / det[:, None, None]

    # xs (Fortran's OUTPUT xs, reshape((/cof11,cof12,cof13,cof21,...,cof33/),(3,3))/det):
    # column-major reshape -> xs(:,1)=[cof11,cof12,cof13], xs(:,2)=[cof21,cof22,cof23],
    # xs(:,3)=[cof31,cof32,cof33] (each triple is one COLUMN, not one row --
    # a previous version of this line stacked them as rows, transposing the
    # result; invisible on every axis-aligned mesh tried so far because the
    # off-diagonal cofactors are then exactly 0 on both sides, and only
    # exposed by tpv10's genuinely skewed dipping-fault elements, where
    # ss(2)/(3)/(5) showed a real ~25% relative miss, not roundoff).
    xs_out = np.stack([
        np.stack([cof11, cof12, cof13], axis=1),  # column 0
        np.stack([cof21, cof22, cof23], axis=1),  # column 1
        np.stack([cof31, cof32, cof33], axis=1),  # column 2
    ], axis=2) / det[:, None, None]  # (E,3,3): xs_out[e,row,col]

    return det, eleshp, xs_out


# func_lib.f90's vlm: 1-indexed node-permutation table `it(8,8)`, verbatim.
_VLM_IT = np.array([
    [1, 2, 3, 4, 5, 6, 7, 8],
    [2, 3, 4, 1, 6, 7, 8, 5],
    [3, 4, 1, 2, 7, 8, 5, 6],
    [4, 1, 2, 3, 8, 5, 6, 7],
    [5, 8, 7, 6, 1, 4, 3, 2],
    [6, 5, 8, 7, 2, 1, 4, 3],
    [7, 6, 5, 8, 3, 2, 1, 4],
    [8, 7, 6, 5, 4, 3, 2, 1],
]) - 1  # 0-indexed; _VLM_IT[k-1, i-1] == Fortran's it(k,i)


def compute_element_volume(xl):
    """Port of func_lib.f90's vlm (Belytschko et al. 1984 hexahedron
    volume formula), vectorized over elements.

    xl: (E,8,3) physical coordinates, same convention as compute_element_det.
    Returns volume: (E,) float array.
    """
    y = xl[:, :, 1]  # xl(2,k) in Fortran -> xl[:, k-1, 1]
    z = xl[:, :, 2]  # xl(3,k)
    x = xl[:, :, 0]  # xl(1,k)

    # _VLM_IT[i-1, k-1] == Fortran it(k,i) (each hardcoded row above is one
    # COLUMN of Fortran's column-major reshape) -- so for fixed k, varying
    # i=1..8, the node-id sequence is _VLM_IT[:, k-1].
    # PERFORMANCE: the formula below re-reads only k in {2,3,4,5,6,8}, but
    # 24 times; gather each (E,8) column block ONCE. Pure gather, no
    # arithmetic touched, so the sum below is bit-for-bit what it was.
    _KS = (2, 3, 4, 5, 6, 8)
    ycol = {k: y[:, _VLM_IT[:, k - 1]] for k in _KS}
    zcol = {k: z[:, _VLM_IT[:, k - 1]] for k in _KS}

    bb = (ycol[2] * (zcol[6] - zcol[3] + zcol[5] - zcol[4]) +
          ycol[3] * (zcol[2] - zcol[4]) +
          ycol[4] * (zcol[3] - zcol[8] + zcol[2] - zcol[5]) +
          ycol[5] * (zcol[8] - zcol[6] + zcol[4] - zcol[2]) +
          ycol[6] * (zcol[5] - zcol[2]) +
          ycol[8] * (zcol[4] - zcol[5]))  # (E,8): bb(i)

    volume = np.sum(x * bb, axis=1) / 12.0
    return volume


# calcSSPhi4Hrgls's ha(8,4) hourglass-mode matrix, verbatim (Fortran reshape
# is column-major: ha(:,1)=col1, ha(:,2)=col2, ...).
_HA = np.array([
    [1, 1, -1, -1, -1, -1, 1, 1],
    [1, -1, -1, 1, -1, 1, 1, -1],
    [1, -1, 1, -1, 1, -1, 1, -1],
    [-1, 1, -1, 1, 1, -1, 1, -1],
], dtype=np.float64).T  # (8,4): _HA[node-1, mode-1] == Fortran ha(node,mode)


def compute_hourglass(xl, xs, mat, eleshp):
    """Milestone 7.5: port of assembleGlobalMass.f90's calcSSPhi4Hrgls.

    xl: (E,8,3), same convention as compute_element_det.
    xs: (E,3,3), from compute_element_shape's `xs` return.
    mat: (E,5) [vp,vs,rho,lambda,mu] (from build_elements).
    eleshp: (E,3,8), from compute_element_shape.

    Returns (ss, phi): ss is (E,6) [Fortran ss(1..6)], phi is (E,8,4)
    [Fortran phi(node,mode,elemID)].
    """
    lam = mat[:, 3]
    mu = mat[:, 4]
    vol = compute_element_volume(xl)
    ce = mu * (3.0 * lam + 2.0 * mu) / (lam + mu)
    ce = 16.0 * ce / 15.0
    co = ce * vol / 48.0

    x1, x2, x3 = xs[:, :, 0], xs[:, :, 1], xs[:, :, 2]  # xs(:,1)/xs(:,2)/xs(:,3) columns
    ss = np.stack([
        co * np.sum(x1 * x1, axis=1),
        co * np.sum(x1 * x2, axis=1),
        co * np.sum(x1 * x3, axis=1),
        co * np.sum(x2 * x2, axis=1),
        co * np.sum(x2 * x3, axis=1),
        co * np.sum(x3 * x3, axis=1),
    ], axis=1)

    # phi'(j,i) = ha(j,i) - sum_k ha(k,i) * (xl(:,k) . eleshp(:,j))
    # dot(e,k,j) = sum_c xl[e,k,c]*eleshp[e,c,j]
    dot = np.einsum('ekc,ecj->ekj', xl, eleshp)  # (E,8,8): dot[e,k,j]
    inner = np.einsum('ki,ekj->eij', _HA, dot)  # (E,4,8): sum_k ha(k,i)*dot(k,j)
    phi_prime = _HA.T[None, :, :] - inner  # (E,4,8): ha(j,i) - inner(i,j)
    norm = np.sqrt(np.sum(phi_prime ** 2, axis=2) / 8.0)  # (E,4)
    phi = phi_prime / norm[:, :, None]  # (E,4,8)
    phi = np.transpose(phi, (0, 2, 1))  # -> (E,8,4): phi[e,node,mode]

    return ss, phi


def init_vel(nsmp, eq_nums, fric, n_equations):
    """Milestone 7.5: port of eqdyna3d.f90's init_vel, mode==1 (non-restart)
    case only -- matching this port's scope (mode==2 reads fric(31:36) from
    a restart netCDF this port does not read).

    Seeds the 1D solver-state velocity vector v1 at fault-node equations
    from fric(31:36) (FRIC_SLOT_VEL_MASTER_X/Y/Z=31/32/33,
    FRIC_SLOT_VEL_SLAVE_X/Y/Z=34/35/36) via the equation-number map -- for
    mode==1 these fric slots are exactly 0.0 (readInputFiles.read_on_fault_vars
    never writes them, matching Fortran's fric=0.0d0 zero-init that is never
    overwritten for mode==1), so this reproduces Fortran's actual observed
    behavior (a documented no-op for this case), not an assumption.

    nsmp: (nftnd,2) int64 [slave_id, master_id], 1-indexed.
    eq_nums: from build_equation_numbers (1-indexed by node id; eq_nums[node]
        is that node's array of equation numbers).
    fric: (nftnd+1,101) friction-state array, 1-indexed rows/cols.
    n_equations: totalNumOfEquations (v1 size).

    Returns v1: (n_equations+1,) float array (row 0 unused, 1-indexed
    equation numbers, matching every other 1-indexed array in this port).
    """
    v1 = np.zeros(n_equations + 1)
    nftnd = nsmp.shape[0]
    FRIC_SLOT_VEL_MASTER_X, FRIC_SLOT_VEL_MASTER_Y, FRIC_SLOT_VEL_MASTER_Z = 31, 32, 33
    FRIC_SLOT_VEL_SLAVE_X, FRIC_SLOT_VEL_SLAVE_Y, FRIC_SLOT_VEL_SLAVE_Z = 34, 35, 36

    for i in range(1, nftnd + 1):
        slave = int(nsmp[i - 1, 0])
        master = int(nsmp[i - 1, 1])
        eqs_slave = eq_nums[slave]  # eq_start[slave] is eqs_slave's own offset -- not needed separately
        v1[eqs_slave[0]] = fric[i, FRIC_SLOT_VEL_SLAVE_X]
        v1[eqs_slave[1]] = fric[i, FRIC_SLOT_VEL_SLAVE_Y]
        v1[eqs_slave[2]] = fric[i, FRIC_SLOT_VEL_SLAVE_Z]
        eqs_master = eq_nums[master]
        v1[eqs_master[0]] = fric[i, FRIC_SLOT_VEL_MASTER_X]
        v1[eqs_master[1]] = fric[i, FRIC_SLOT_VEL_MASTER_Y]
        v1[eqs_master[2]] = fric[i, FRIC_SLOT_VEL_MASTER_Z]

    return v1
