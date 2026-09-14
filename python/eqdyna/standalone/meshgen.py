"""
Standalone (no-Fortran-in-the-loop) port of EQdyna's mesh generation for
the SERIAL (npx=npy=npz=1) case only -- src/countMeshEntities.f90 (mesh4num)
+ src/meshgen.f90's `getLocalOneDimCoorArrAndSize` + node-coordinate loop
+ createElement/setElementMaterial/replaceSlaveWithMasterNode.

Milestone 1 (node coordinates), Milestone 2 (element connectivity +
per-element material), Milestone 3 (equation numbering), Milestone 4
(split-node fault geometry: un/us/ud unit vectors + arn nodal area),
Milestone 5 (on-/off-fault station node matching), and Milestone 6
(native case-input reading -- see native_input.py, which replaces every
hand-transcribed PARAMS/MATERIAL/station value M1-M5 previously used) of
the standalone-phase spec are ported and verified against
`testsys/parity/fixtures/test_tpv8_serial/pydump_meshCoor.txt` /
`pydump_conn.txt` / `pydump_nodeinfo.txt` / `pydump_fault.txt` /
`pydump_stations.txt` (the Fortran ground truth, dumped by
`src/pydump.f90`), via `testsys/parity/test_standalone_meshgen.py` (wired
into `testsys/run.py parity`): M1 max abs diff 3.6e-12 across 249,047
nodes; M2 zero mismatched elements (connectivity, elemType, material)
across 235,008 elements; M3 zero mismatches (numOfDofPerNodeArr,
eqNumStartIndexLoc, eqNumIndexArr incl. the -1 fixed-boundary sentinel and
the 3-vs-12 dof PML split) across all 249,047 nodes and totalNumOfEquations
bit-exact (1,425,873); M4 un/us/ud bit-exact and arn max abs diff at
roundoff across all 1,891 fault nodes; M5 station counts and every matched
(station,node) pair bit-exact, in Fortran's own match order, across all 8
on-fault and 11/15 matched off-fault stations; M6 re-runs M1-M5 against
NATIVELY-READ inputs (unchanged results, confirming the readers) plus a new
initial-fric-state check (netCDF4 read of on_fault_vars_input.nc) at max
abs diff 7.45e-9 across all 1,891 fault nodes x 100 friction slots; M6.5
(see mass_assembly.py) ports assembleGlobalMass.f90's lumped-mass
integration, verified against pydump_nodalmass.txt/pydump_fnms.txt at max
relative diff 8.2e-16 (true roundoff) across 1,425,873 equations / 249,047
nodes.

M3 found one real bug before it shipped: `setEquationNumber`'s fixed-
-boundary test compares against the model's ACTUAL computed grid extent
(`xline[0]`/`xline[-1]`/etc, i.e. Fortran's `xmin=modelBoundCoor(1,1)`
overwrite after grid-building), NOT the requested input bounds in
`params` -- the two differ because the grid-building loop stops once it
reaches or passes the requested bound. Using the input bounds silently
inflated totalNumOfEquations by ~230k (too few nodes classified fixed).
See `build_equation_numbers`'s docstring/comment for the fix.

Milestone 7.5 (see mass_assembly.py) closes the assembly gap this
docstring used to flag as "not yet ported": `compute_element_shape` now
ports eledet/eleshp (calcGlobalShapeFunc's full output -- spatial shape-
function derivatives dNx/dNy/dNz -- for both interior and PML elements,
since C_degen==0 means neither ever takes the wedge-degeneration branch),
`compute_hourglass` ports ss/phi (calcSSPhi4Hrgls's hourglass-control
tensors, including func_lib.f90's vlm volume formula), and `init_vel`
ports eqdyna3d.f90's fault-node velocity seeding of v1 from fric(31-36)
via the equation map (mode==1 only). Verified against
pydump_elemgeo.txt/pydump_v1.txt via test_standalone_meshgen.py.

NOT YET PORTED (explicitly, not silently deferred): assembleGlobalKU (the
FEM stiffness/internal-force kernel itself), hrglss's per-step hourglass
FORCE (as opposed to the ss/phi tensors M7.5 now provides), velDispUpdate,
faulting -- i.e. everything in the actual time-stepping loop. Scoped this
way deliberately -- shipping a "mostly right" mesh generator with the
fault/station logic silently stubbed would be exactly the kind of
unverified deliverable this discipline forbids. See the module's TODO
list at the bottom for what a follow-up milestone needs next, in
dependency order.

Serial-only simplification used throughout: `getLocalOneDimCoorArrAndSize`'s
MPI-partition slicing (numOfMPIXyz>1 branch) is never exercised at
npx=npy=npz=1 -- the "local" 1D coordinate array IS the global one. This
is not a shortcut that changes results for npx=npy=npz=1; it is the same
code path the Fortran takes for MPIXyzId=0, numOfMPIXyz=1
(`residualNumOfNodes=0`, `localOneDimCoorArrSize=globalOneDimCoorArrSize`).
"""
import numpy as np


def one_dim_coor_array(dim_id, dx, dy, dz, fxmin, fxmax, fzmin, fzmax,
                        dis4uniF, dis4uniB, xmin, xmax, ymin, ymax, zmin, zmax,
                        rat, nPML, np_max=1000000):
    """Port of getLocalOneDimCoorArrAndSize, serial case (numOfMPIXyz=1).
    dim_id: 1=x, 2=y, 3=z. Returns (coor_array, PMLb_partial dict)."""
    if dim_id == 1:
        n_uniform = round((fxmax - fxmin) / dx) + 1
        grid_size = dx
        front_edge, back_edge = fxmin, fxmax
        min_coor, max_coor = xmin, xmax
    elif dim_id == 2:
        n_uniform = dis4uniF + dis4uniB + 1
        grid_size = dy
        front_edge, back_edge = -dis4uniF * dy, dis4uniB * dy
        min_coor, max_coor = ymin, ymax
    elif dim_id == 3:
        n_uniform = round((fzmax - fzmin) / dz) + 1
        grid_size = dz
        front_edge, back_edge = fzmin, fzmax
        min_coor, max_coor = zmin, zmax
    else:
        raise ValueError(dim_id)

    coor_tmp = front_edge
    grid_tmp = grid_size
    i = 0
    for i in range(1, np_max + 1):
        grid_tmp *= rat
        coor_tmp -= grid_tmp
        if coor_tmp <= min_coor:
            break
    front_edge_node_id = i + nPML

    coor_tmp = back_edge
    grid_tmp = grid_size
    j = 0
    for j in range(1, np_max + 1):
        grid_tmp *= rat
        coor_tmp += grid_tmp
        if coor_tmp >= max_coor:
            break
    if dim_id == 3:
        j = -nPML

    global_size = n_uniform + front_edge_node_id + j + nPML
    arr = np.zeros(global_size)

    arr[front_edge_node_id + 1 - 1] = front_edge  # 1-indexed Fortran -> 0-indexed here
    grid_tmp = grid_size
    for i in range(front_edge_node_id, 0, -1):
        grid_tmp *= rat
        arr[i - 1] = arr[i + 1 - 1] - grid_tmp
    for i in range(front_edge_node_id + 2, front_edge_node_id + n_uniform + 1):
        arr[i - 1] = arr[i - 1 - 1] + grid_size
    if dim_id < 3:
        grid_tmp = grid_size
        for i in range(front_edge_node_id + n_uniform + 1, global_size + 1):
            grid_tmp *= rat
            arr[i - 1] = arr[i - 1 - 1] + grid_tmp

    pmlb = {}
    if dim_id == 1:
        pmlb['xmax0'] = arr[global_size - nPML - 1]
        pmlb['xmin0'] = arr[nPML + 1 - 1]
        pmlb['maxdx'] = arr[global_size - 1] - arr[global_size - 1 - 1]
        model_bound = (arr[0], arr[-1])
    elif dim_id == 2:
        pmlb['ymax0'] = arr[global_size - nPML - 1]
        pmlb['ymin0'] = arr[nPML + 1 - 1]
        pmlb['maxdy'] = arr[global_size - 1] - arr[global_size - 1 - 1]
        model_bound = (arr[0], arr[-1])
    else:
        pmlb['zmin0'] = arr[nPML + 1 - 1]
        pmlb['maxdz'] = arr[1] - arr[0]
        model_bound = (arr[0], arr[-1])
    return arr, pmlb, model_bound


def build_grid_lines(params):
    """params: dict with dx,dy,dz,fxmin,fxmax,fzmin,fzmax,dis4uniF,dis4uniB,
    xmin,xmax,ymin,ymax,zmin,zmax,rat,nPML. Returns (xline, yline, zline,
    PMLb dict, model bound coors)."""
    p = params
    xline, pmlx, xbound = one_dim_coor_array(
        1, p['dx'], p['dy'], p['dz'], p['fxmin'], p['fxmax'], p['fzmin'], p['fzmax'],
        p['dis4uniF'], p['dis4uniB'], p['xmin'], p['xmax'], p['ymin'], p['ymax'],
        p['zmin'], p['zmax'], p['rat'], p['nPML'])
    yline, pmly, ybound = one_dim_coor_array(
        2, p['dx'], p['dy'], p['dz'], p['fxmin'], p['fxmax'], p['fzmin'], p['fzmax'],
        p['dis4uniF'], p['dis4uniB'], p['xmin'], p['xmax'], p['ymin'], p['ymax'],
        p['zmin'], p['zmax'], p['rat'], p['nPML'])
    zline, pmlz, zbound = one_dim_coor_array(
        3, p['dx'], p['dy'], p['dz'], p['fxmin'], p['fxmax'], p['fzmin'], p['fzmax'],
        p['dis4uniF'], p['dis4uniB'], p['xmin'], p['xmax'], p['ymin'], p['ymax'],
        p['zmin'], p['zmax'], p['rat'], p['nPML'])
    pmlb = dict(**pmlx, **pmly, **pmlz)
    return xline, yline, zline, pmlb, (xbound, ybound, zbound)


def is_on_fault(x, y, z, fxmin, fxmax, fymin, fymax, fzmin, fzmax, tol, c_degen=0.0):
    """Port of checkIsOnFault, C_degen==0 branch only (planar fault at y=0;
    tpv8/tpv104 both use this branch -- insertFaultType>0/C_degen>3 rough-
    fault branches are NOT ported here, matching this milestone's scope)."""
    if c_degen != 0.0:
        raise NotImplementedError('only the planar (C_degen==0) fault branch is ported')
    in_box = (x >= fxmin - tol and x <= fxmax + tol and
              y >= fymin - tol and y <= fymax + tol and
              z >= fzmin - tol and z <= fzmax + tol)
    return in_box and y == 0.0


def _fortran_nint(x):
    """Fortran nint(): round-half-away-from-zero. Same formula as
    native_input.py's helper of the same purpose (kept local here to avoid
    a cross-module dependency for a two-line function)."""
    return int(np.sign(x) * np.floor(np.abs(x) + 0.5)) if x != 0 else 0


def insert_fault_interface(x, y, z, rough, dx, dz, ymin, ymax, tol):
    """Port of func_lib.f90's insertFaultInterface (Milestone 9): given a
    node's UNDISTORTED (x, y, z) -- the same coordinates checkIsOnFault
    tests against, per the Fortran calling this from meshgen's loop BEFORE
    checkIsOnFault runs on the same nodeCoor -- looks up the local fault-
    surface height `peak` and its along-strike/along-dip slopes (pfx, pfz)
    from `rough` (native_input.read_fault_rough_geometry's dict), then
    returns the MORPHED y-coordinate `ycoort` via meshgen.f90's linear
    blend between the fault surface and the model's ymin/ymax planes.

    NOTE (verbatim, not "fixed"): the column-index math (ixx, izz) uses
    the MODEL's dx/dz (bModelGeometry.txt), NOT the rough-geometry file's
    own grid spacing (rough['dx'], used only to derive rough_fx_max) --
    this is what src/func_lib.f90 actually does (bare `dx`/`dz` module
    globals inside insertFaultInterface), and is only correct because
    these benchmark cases' rough-geometry sampling grid is generated at
    the same spacing as the model mesh; reproduced as-is per this
    discipline's rule against silently "fixing" behavior it doesn't
    understand outside its ported scope.

    Returns (ycoort, pfx, pfz) -- pfx/pfz are also needed by
    build_fault_geometry's insertFaultType>0 un/us/ud branch.
    """
    fx1, fx2, fz1 = rough['fx_min'], rough['fx_max'], rough['fz_min']
    nnx, nnz, rough_geo = rough['nnx'], rough['nnz'], rough['rough_geo']

    if fx1 - tol < x < fx2 + tol and z > fz1 - tol:
        ixx = _fortran_nint((x - fx1) / dx) + 1
        izz = _fortran_nint((z - fz1) / dz) + 1
    elif x < fx1 - tol and z > fz1 - tol:
        ixx = 1
        izz = _fortran_nint((z - fz1) / dz) + 1
    elif x > fx2 + tol and z > fz1 - tol:
        ixx = nnx
        izz = _fortran_nint((z - fz1) / dz) + 1
    elif fx1 - tol < x < fx2 + tol and z < fz1 - tol:
        ixx = _fortran_nint((x - fx1) / dx) + 1
        izz = 1
    elif x < fx1 - tol and z < fz1 - tol:
        ixx, izz = 1, 1
    elif x > fx2 + tol and z < fz1 - tol:
        ixx, izz = nnx, 1
    else:
        raise ValueError('insert_fault_interface: (x=%r, z=%r) matched none of the six '
                          'Fortran if/elseif branches (exactly on a boundary?)' % (x, z))

    col = nnz * (ixx - 1) + izz  # 1-indexed, Fortran rough_geo(:, nnz*(ixx-1)+izz)
    peak, pfx, pfz = rough_geo[0, col - 1], rough_geo[1, col - 1], rough_geo[2, col - 1]

    if y > -tol:
        ycoort = y * (ymax - peak) / ymax + peak
    elif y < -tol:
        ycoort = y * (peak - ymin) / (-ymin) + peak
    else:
        raise ValueError('insert_fault_interface: y=%r is within tol of 0 without being '
                          '> -tol or < -tol (exact-tie edge case Fortran leaves undefined)' % y)

    return ycoort, pfx, pfz


def build_node_coordinates(xline, yline, zline, params):
    """Port of countMeshEntities/meshgen's node-creation loop (do ix; do iz;
    do iy) for a SINGLE planar fault (ntotft==1, C_degen==0). Returns
    meshCoor (N,3) in Fortran 1-indexed node-id order (row 0 unused,
    matching the pydump dumps' 1-indexed convention read elsewhere in
    python/eqdyna), nftnd (int), and nsmp ((nftnd,2) int64 [slave_id,
    master_id], both 1-indexed, in fault-encounter order -- the same pairs
    meshgen.f90's createMasterNode writes into nsmp(1,:,1)/nsmp(2,:,1)).

    Master nodes are NOT interleaved with their slave: meshgen.f90 assigns
    master-node ids starting at `msnode = nx*ny*nz + nftnd0(iFault)`, i.e.
    APPENDED after all nx*ny*nz regular-grid nodes, in fault-encounter
    order (still ix-outer, iz-middle, iy-inner, since that's the loop
    order that reaches them) -- NOT interleaved immediately after each
    slave row. (countMeshEntities.f90's counting-pass loop increments its
    running node COUNT at the point a fault node is found only to get the
    right TOTAL; it says nothing about final array layout -- confirmed by
    diffing this port's first attempt, which interleaved, against the
    golden `pydump_meshCoor.txt`: node counts matched exactly, per-node
    coordinates did not, until switched to this appended-at-end layout.)
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    rough = p.get('rough')
    insert_fault_type = p.get('insertFaultType', 0)
    regular = []
    master = []
    nsmp = []  # (slave_node_id, master_node_id), both 1-indexed, in fault-encounter order
    for ix in range(nx):
        xcoor = xline[ix]
        for iz in range(nz):
            zcoor = zline[iz]
            for iy in range(ny):
                ycoor = yline[iy]
                # checkIsOnFault (via is_on_fault below) always tests the
                # UNDISTORTED ycoor -- meshgen.f90 calls insertFaultInterface
                # only to get ycoort for STORAGE, never mutates nodeCoor
                # itself before the fault test runs.
                is_fault = is_on_fault(xcoor, ycoor, zcoor, p['fxmin'], p['fxmax'],
                                        p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], tol)
                y_store = ycoor
                if insert_fault_type > 0:
                    # ymin/ymax here MUST be the grid-derived bounds (yline[0]/
                    # yline[-1]), not params['ymin']/['ymax'] (the requested
                    # input bounds) -- same trap M3's build_equation_numbers
                    # docstring already flags for xmin/xmax/zmin.
                    y_store, _, _ = insert_fault_interface(
                        xcoor, ycoor, zcoor, rough, p['dx'], p['dz'], yline[0], yline[-1], tol)
                regular.append((xcoor, y_store, zcoor))
                slave_id = len(regular)  # 1-indexed nodeCount at this point
                if is_fault:
                    master.append((xcoor, y_store, zcoor))
                    master_id = nx * ny * nz + len(master)  # msnode = nx*ny*nz + nftnd0
                    nsmp.append((slave_id, master_id))
    nftnd = len(master)
    N = len(regular) + nftnd
    meshCoor = np.zeros((N + 1, 3))  # index 0 unused (1-indexed node ids)
    meshCoor[1:1 + len(regular)] = regular
    meshCoor[1 + len(regular):] = master
    nsmp = np.array(nsmp, dtype=np.int64)  # (nftnd, 2)
    return meshCoor, nftnd, nsmp

def build_elements(xline, yline, zline, params, pmlb, nsmp, material, meshCoor):
    """Milestone 2: port of meshgen.f90's createElement + setElementMaterial +
    replaceSlaveWithMasterNode, single planar fault (ntotft==1, C_degen==0),
    homogeneous or 1D-layered material (setElementMaterial's two branches).

    Verified against `testsys/parity/fixtures/test_tpv8_serial/pydump_conn.txt`
    (nodeElemIdRelation + elemTypeArr + mat(1:5) for every element, dumped by
    src/pydump.f90) -- see testsys/parity/test_standalone_meshgen.py.

    Reproduces the Fortran plane1/plane2 sliding-column bookkeeping exactly
    (scalar loop, not vectorized -- Phase B of the port workflow: verbatim
    first, vectorize only after this checkpoint passes) because the master-
    node overwrite into plane2's extra row (row index ny, 0-indexed, one row
    beyond the ny regular rows) has to land in the SAME column-shift timing
    as Fortran's `plane1 = plane2` at the end of each ix iteration, or
    elements adjacent to the fault silently get the wrong corner node.

    params: same dict as build_grid_lines.
    pmlb: dict from build_grid_lines (xmax0/xmin0/ymax0/ymin0/zmin0 keys).
    nsmp: (nftnd,2) int64 [slave_id, master_id], from build_node_coordinates.
    material: (nmat, n2mat) array (n2mat==3: homogeneous [vp,vs,rho];
        n2mat==4: 1D layered [depth_bottom,vp,vs,rho], see setElementMaterial).
    meshCoor: (N+1,3) from build_node_coordinates -- Milestone 9 addition:
        element-center classification (PML alignment/type, material depth)
        reads the ACTUAL node coordinates here (matching createElement.f90's
        `meshCoor(j,nodeElemIdRelation(i,elemCount))`), not the grid lines
        directly, so insertFaultType>0's y-morph is correctly reflected.

    Returns (conn, elem_type, mat, depth) where conn is (E,8) 1-indexed node
    ids in the Fortran nodeElemIdRelation column order (0-unused row NOT
    included -- conn is 0-indexed by element, elements 1..E map to rows
    0..E-1), elem_type is (E,) int (1=interior/hourglass-controlled, 2=PML),
    mat is (E,5) [vp,vs,rho,lambda,mu], and depth (E,) is
    `-0.5*(zline[iz]+zline[iz-1]) + 7.3215` (meshgen.f90:103's argument to
    setPlasticStress, verbatim including the 7.3215 magic-number shift --
    ALWAYS computed, harmless when C_elastic==1 since no caller reads it
    then; consumed by main.py's build_solver_state only when C_elastic==0,
    to seed each interior/PML element's lithostatic pre-stress -- Milestone
    10, drv.a6's Drucker-Prager viscoplasticity).
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    dx, dy = p['dx'], p['dy']
    fxmin, fxmax, fzmin = p['fxmin'], p['fxmax'], p['fzmin']

    # PMLb per getLocalOneDimCoorArrAndSize's mapping (build_grid_lines' pmlb dict).
    xmax0, xmin0 = pmlb['xmax0'], pmlb['xmin0']
    ymax0, ymin0 = pmlb['ymax0'], pmlb['ymin0']
    zmin0 = pmlb['zmin0']

    slave2master = {int(s): int(m) for s, m in nsmp}

    n_elem = (nx - 1) * (ny - 1) * (nz - 1)
    conn = np.zeros((n_elem, 8), dtype=np.int64)
    elem_type = np.zeros(n_elem, dtype=np.int64)

    nmat, n2mat = material.shape
    if nmat == 1 and n2mat == 3:
        vp, vs, rho = material[0, 0], material[0, 1], material[0, 2]
        mu = vs * vs * rho
        lam = vp * vp * rho - 2.0 * mu
        mat_row_homog = np.array([vp, vs, rho, lam, mu])
    elif not (nmat > 1 and n2mat == 4):
        raise NotImplementedError('setElementMaterial: only nmat==1/n2mat==3 '
                                   '(homogeneous) or nmat>1/n2mat==4 (1D '
                                   'layered) branches are ported')

    def material_for(elem_center_z):
        if nmat == 1 and n2mat == 3:
            return mat_row_homog
        depth = abs(elem_center_z)
        if depth < material[0, 0]:
            vp, vs, rho = material[0, 1], material[0, 2], material[0, 3]
        else:
            row = None
            for i in range(1, nmat):
                if depth < material[i, 0] and depth >= material[i - 1, 0]:
                    row = i
                    break
            if row is None:
                raise ValueError('setElementMaterial: depth %r not covered '
                                  'by any material layer' % depth)
            vp, vs, rho = material[row, 1], material[row, 2], material[row, 3]
        mu = vs * vs * rho
        lam = vp * vp * rho - 2.0 * mu
        return np.array([vp, vs, rho, lam, mu])

    mat = np.zeros((n_elem, 5))
    depth = np.zeros(n_elem)

    # plane1/plane2: (ny+1) x nz, row ny (0-indexed) is the ntotft==1 master row.
    plane1 = np.zeros((ny + 1, nz), dtype=np.int64)
    plane2 = np.zeros((ny + 1, nz), dtype=np.int64)

    node_count = 0
    elem_count = 0
    master_count = 0
    for ix in range(nx):
        xcoor = xline[ix]
        for iz in range(nz):
            zcoor = zline[iz]
            for iy in range(ny):
                ycoor = yline[iy]
                node_count += 1
                plane2[iy, iz] = node_count

                if is_on_fault(xcoor, ycoor, zcoor, fxmin, fxmax, p['fymin'],
                                p['fymax'], fzmin, p['fzmax'], tol):
                    master_count += 1
                    msnode = nx * ny * nz + master_count
                    plane2[ny, iz] = msnode

                if ix >= 1 and iy >= 1 and iz >= 1:
                    c = np.array([
                        plane1[iy - 1, iz - 1], plane2[iy - 1, iz - 1],
                        plane2[iy, iz - 1], plane1[iy, iz - 1],
                        plane1[iy - 1, iz], plane2[iy - 1, iz],
                        plane2[iy, iz], plane1[iy, iz],
                    ], dtype=np.int64)

                    # elementCenterCoor = mean of the 8 corner coords, from
                    # meshCoor(j, nodeElemIdRelation(i,elemCount)) per
                    # createElement.f90 -- for insertFaultType==0 this is
                    # bit-identical to the grid-line lookup (meshCoor's y IS
                    # yline's y there), but for insertFaultType>0 the actual
                    # y is the MORPHED value build_node_coordinates already
                    # stored (per-node, not per-grid-line), so this MUST come
                    # from meshCoor, not yline -- PML classification below
                    # depends on the true (warped) element center.
                    corners = meshCoor[c]  # (8,3)
                    cx, cy, cz = corners.mean(axis=0)

                    if cx == xmax0 or cx == xmin0 or cy == ymax0 or cy == ymin0 or cz == zmin0:
                        raise ValueError(
                            'checkPMLAlignment: element center exactly on a '
                            'PML bound at (%r,%r,%r)' % (cx, cy, cz))

                    etype = 1
                    if cx > xmax0 or cx < xmin0 or cy > ymax0 or cy < ymin0 or cz < zmin0:
                        etype = 2

                    # replaceSlaveWithMasterNode: only the C_degen==0,
                    # elemTypeArr==1 branch applies here (12/13 come from the
                    # wedge()/C_degen>3 branch, out of scope). Test uses the
                    # top node's coords (this ix,iy,iz), matching Fortran's
                    # `nodeCoor` at the point createElement/replaceSlave... run.
                    if etype == 1 and (xcoor > fxmin - tol and xcoor < fxmax + dx + tol
                                       and zcoor > fzmin - tol and ycoor > 0.0
                                       and abs(ycoor - dy) < tol):
                        for k in range(8):
                            nid = int(c[k])
                            if nid in slave2master:
                                c[k] = slave2master[nid]

                    conn[elem_count] = c
                    elem_type[elem_count] = etype
                    mat[elem_count] = material_for(cz)
                    # meshgen.f90:103 `setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1))
                    # + 7.3215d0, elemCount)` -- python's `iz`/`zline` here are the
                    # SAME running loop variable/array Fortran uses (both 0-indexed
                    # consistently, see build_elements' docstring math), so this is
                    # the identical expression, not a re-derivation.
                    depth[elem_count] = -0.5 * (zline[iz] + zline[iz - 1]) + 7.3215
                    elem_count += 1
        plane1 = plane2.copy()

    assert elem_count == n_elem, (elem_count, n_elem)
    return conn, elem_type, mat, depth


def build_equation_numbers(xline, yline, zline, params, pmlb):
    """Milestone 3: port of meshgen.f90's setNumDof + setEquationNumber +
    the equation-number half of createMasterNode, single planar fault
    (ntotft==1), serial (npx=npy=npz=1).

    Reproduces the Fortran's INTERLEAVED flat-array layout exactly: for each
    (ix,iz,iy) grid point, the regular node's `numOfDofPerNodeTmp` equation-
    number slots are appended to the running `eqNumIndexArr` tag FIRST (via
    setEquationNumber, called before createMasterNode in meshgen.f90's main
    loop), and if that node is on the fault, the newly-created master node's
    3 equation-number slots are appended immediately after -- NOT appended
    at the very end like meshCoor's master-node block (M1's docstring flags
    that trap for node coordinates; equation numbers are a DIFFERENT layout,
    interleaved, and must not be confused with it). `equationNumCount` (the
    actual integer values assigned) is a single running counter shared by
    both regular and master nodes in this same interleaved order, so getting
    the order right is required to reproduce the exact equation-number
    VALUES, not just the per-node dof counts.

    Verified against `testsys/parity/fixtures/test_tpv8_serial/pydump_nodeinfo.txt`
    (dumped per Fortran node id 1..totalNumOfNodes: numOfDofPerNodeArr(i),
    eqNumStartIndexLoc(i), eqNumIndexArr(start+1..start+dof)) via
    testsys/parity/test_standalone_meshgen.py.

    Returns (num_dof, eq_start, eq_nums, total_num_of_equations) where
    num_dof is (N+1,) int (3 or 12, row 0 unused), eq_start is (N+1,) int
    (0-indexed flat-array offset before this node's block, matching
    Fortran's eqNumStartIndexLoc convention of "tag value before this
    node's writes"), eq_nums is a length-(N+1) list of 1D int arrays (row 0
    unused) holding this node's assigned equation numbers (positive) or -1
    (fixed-boundary sentinel) per dof slot, and total_num_of_equations is
    the final equationNumCount (== totalNumOfEquations).
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    # NOTE: meshgen.f90 overwrites its own xmin/xmax/ymin/ymax/zmin/zmax with
    # modelBoundCoor (the ACTUAL computed grid extent) right after building
    # the grid lines -- these are xline[0]/xline[-1]/etc, NOT the requested
    # input bounds in `params` (the grid-building loop stops once it reaches
    # or exceeds the requested bound, so the two differ by design). Using
    # the input params here (an earlier version of this function did) undercounts
    # fixed-boundary nodes and silently inflates totalNumOfEquations.
    xmin, xmax = xline[0], xline[-1]
    ymin, ymax = yline[0], yline[-1]
    zmin = zline[0]
    xmax0, xmin0 = pmlb['xmax0'], pmlb['xmin0']
    ymax0, ymin0 = pmlb['ymax0'], pmlb['ymin0']
    zmin0 = pmlb['zmin0']
    ndof = 3

    N_regular = nx * ny * nz

    def is_fixed_boundary(x, y, z):
        return (abs(x - xmin) < tol or abs(x - xmax) < tol or
                abs(y - ymin) < tol or abs(y - ymax) < tol or
                abs(z - zmin) < tol)  # zmax (free surface) deliberately excluded

    def num_dof_for(x, y, z):
        if x > xmax0 or x < xmin0 or y > ymax0 or y < ymin0 or z < zmin0:
            return 12
        return 3

    # Node ids are NOT visited in id order (master ids, N_regular+k, are only
    # reached when a fault node is encountered mid-traversal, long before
    # node_count reaches N_regular) -- index by node id via dict/array
    # assignment, not list-append, to avoid corrupting the append position
    # for subsequent regular nodes.
    tag = 0
    eq_count = 0
    node_count = 0
    master_count = 0
    # nftnd is discoverable ahead of time cheaply (same predicate as below);
    # count first so num_dof/eq_start can be preallocated by node id.
    nftnd = sum(
        1 for ix in range(nx) for iz in range(nz) for iy in range(ny)
        if is_on_fault(xline[ix], yline[iy], zline[iz], p['fxmin'], p['fxmax'],
                        p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], tol))
    N_total = N_regular + nftnd
    num_dof = np.zeros(N_total + 1, dtype=np.int64)
    eq_start = np.zeros(N_total + 1, dtype=np.int64)
    eq_nums = [None] * (N_total + 1)

    for ix in range(nx):
        xcoor = xline[ix]
        for iz in range(nz):
            zcoor = zline[iz]
            for iy in range(ny):
                ycoor = yline[iy]
                node_count += 1
                ndpn = num_dof_for(xcoor, ycoor, zcoor)

                eq_start[node_count] = tag
                num_dof[node_count] = ndpn
                slots = np.empty(ndpn, dtype=np.int64)
                fixed = is_fixed_boundary(xcoor, ycoor, zcoor)
                for d in range(ndpn):
                    tag += 1
                    if fixed:
                        slots[d] = -1
                    else:
                        eq_count += 1
                        slots[d] = eq_count
                eq_nums[node_count] = slots

                if is_on_fault(xcoor, ycoor, zcoor, p['fxmin'], p['fxmax'],
                                p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], tol):
                    master_count += 1
                    msnode = N_regular + master_count
                    eq_start[msnode] = tag
                    num_dof[msnode] = ndof
                    slots2 = np.empty(ndof, dtype=np.int64)
                    for d in range(ndof):
                        tag += 1
                        eq_count += 1
                        slots2[d] = eq_count
                    eq_nums[msnode] = slots2

    assert node_count == N_regular, (node_count, N_regular)
    assert master_count == nftnd, (master_count, nftnd)
    return num_dof, eq_start, eq_nums, eq_count


def build_fault_geometry(xline, yline, zline, params, nsmp):
    """Milestone 4 (+ Milestone 9's insertFaultType>0 extension): port of
    meshgen.f90's split-node unit-vector assignment (`createMasterNode`'s
    un/us/ud writes) and the on-fault quadrilateral-area accumulation onto
    `arn` (meshgen's main-loop area block after the ix/iy/iz node loop),
    for a single (ntotft==1) fully rectangular fault -- planar
    (C_degen<=3, insertFaultType==0: angle-based un/us/ud, tpv8/tpv104's
    branch) OR dipping/rough (insertFaultType>0: pfx/pfz-derived un/us/ud
    from func_lib.f90's insertFaultInterface via meshgen.f90:804-817,
    tpv10/drv.a6's branch -- ud = us x un, matching the Fortran's explicit
    cross-product component formulas exactly). For insertFaultType>0, the
    fault surface's actual y-coordinate is `peak` (from
    insert_fault_interface), not 0 -- used for arn's corner-distance
    formula below, matching what build_node_coordinates already stored
    into meshCoor for these same nodes.

    un/us/ud depend only on the fault's strike/dip angles (`fltxyz(:,4,1)`
    in Fortran, read from bGlobal.txt as fstrike/fdip and converted:
    fltxyz(2,4,1) = 90 deg UNLESS C_degen>3, in which case it's C_degen
    itself, in degrees -- both then to radians) -- NOT the per-node
    position, for this planar branch. Verified against
    `testsys/parity/fixtures/test_tpv8_serial/pydump_fault.txt` columns
    3-11 (un,us,ud) and column 12 (arn) via
    testsys/parity/test_standalone_meshgen.py.

    `MPI4arn`'s neighbor-boundary exchange (meshgen.f90's call right after
    the area-accumulation block) is entirely guarded behind `if (npx>1)` /
    `if (npy>1)` / `if (npz>1)` -- for the serial (npx=npy=npz=1) case this
    milestone targets, it is a documented no-op, confirmed by reading the
    guard conditions directly (not assumed): nothing to port.

    NOT ported here: `fnms` (nodal fault mass, `pydump_fnms.txt`) -- that
    array is written by `assembleGlobalMass.f90` (a ~400-line FEM mass-
    integration subsystem, not mesh geometry) and is out of scope for a
    meshgen milestone; flagged for whichever future milestone ports mass
    assembly.

    params: same dict as build_grid_lines, plus 'fstrike' (degrees) and
    'C_degen' (0.0 for the planar branch this milestone ports).
    nsmp: (nftnd,2) int64 [slave_id, master_id], from build_node_coordinates,
    in fault-encounter (ix-outer, iz-middle, iy-inner) order -- the same
    order `createMasterNode`'s running fault-node counter (nftnd0) assigns,
    so row k (0-indexed) of nsmp IS fault node k+1.

    Returns (un, us, ud, arn): un/us/ud are (nftnd+1, 3) float arrays (row 0
    unused, 1-indexed fault-node rows, matching this module's other
    1-indexed conventions); arn is (nftnd+1,) float (row 0 unused).
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    nftnd = nsmp.shape[0]
    rough = p.get('rough')
    insert_fault_type = p.get('insertFaultType', 0)

    fstrike = p['fstrike'] * np.pi / 180.0
    fdip = (p['C_degen'] * np.pi / 180.0) if p['C_degen'] > 3.0 else (90.0 * np.pi / 180.0)

    un = np.zeros((nftnd + 1, 3))
    us = np.zeros((nftnd + 1, 3))
    ud = np.zeros((nftnd + 1, 3))
    un[1:] = (np.cos(fstrike) * np.sin(fdip), -np.sin(fstrike) * np.sin(fdip), np.cos(fdip))
    us[1:] = (-np.sin(fstrike), -np.cos(fstrike), 0.0)
    ud[1:] = (np.cos(fstrike) * np.cos(fdip), np.sin(fstrike) * np.cos(fdip), np.sin(fdip))

    # Reproduce createMasterNode's fltrc(ifs,ifd) bookkeeping: ixfi/izfi are
    # each set ONCE, on the very first fault-node encounter across the whole
    # traversal (not reset per ix), then ifs/ifd are offsets from those --
    # same traversal order as build_node_coordinates, so `seq` here lines up
    # 1:1 with nsmp's row order.
    grid = {}
    ixfi = izfi = None
    seq = 0
    for ix in range(nx):
        xcoor = xline[ix]
        for iz in range(nz):
            zcoor = zline[iz]
            for iy in range(ny):
                ycoor = yline[iy]
                if is_on_fault(xcoor, ycoor, zcoor, p['fxmin'], p['fxmax'],
                                p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], tol):
                    seq += 1
                    ix_f, iz_f = ix + 1, iz + 1
                    if ixfi is None:
                        ixfi = ix_f
                    if izfi is None:
                        izfi = iz_f
                    ifs = ix_f - ixfi + 1
                    ifd = iz_f - izfi + 1
                    y_geo = ycoor  # planar branch: the fault's actual y IS 0 here.
                    if insert_fault_type > 0:
                        # insertFaultType>0: createMasterNode's un/us/ud
                        # branch (meshgen.f90:804-817) OVERWRITES the
                        # angle-based un/us/ud above with pfx/pfz-derived
                        # values, per fault node -- and the fault's actual
                        # (warped) y-coordinate is `peak`, not 0, which
                        # matters for arn's corner-distance formula below
                        # (meshCoor already stores this same `peak` value,
                        # per build_node_coordinates' y-morph of fault
                        # nodes -- recomputed here rather than re-reading
                        # meshCoor, since is_on_fault's traversal order is
                        # independently re-walked in every M-builder).
                        y_geo, pfx, pfz = insert_fault_interface(
                            xcoor, ycoor, zcoor, rough, p['dx'], p['dz'],
                            yline[0], yline[-1], tol)
                        denom = (pfx ** 2 + 1.0 + pfz ** 2) ** 0.5
                        un[seq] = (-pfx / denom, 1.0 / denom, -pfz / denom)
                        us_denom = (1.0 + pfx ** 2) ** 0.5
                        us[seq] = (1.0 / us_denom, pfx / us_denom, 0.0)
                        ud[seq] = np.cross(us[seq], un[seq])
                    grid[(ifs, ifd)] = (seq, (xcoor, y_geo, zcoor))
    assert seq == nftnd, (seq, nftnd)
    ns = max(k[0] for k in grid)
    nd = max(k[1] for k in grid)
    if len(grid) != ns * nd:
        raise NotImplementedError('build_fault_geometry: only a fully '
                                   'rectangular planar fault grid is ported')

    arn = np.zeros(nftnd + 1)
    for i in range(2, nd + 1):
        for j in range(2, ns + 1):
            m1, c1 = grid[(j, i)]
            m2, c2 = grid[(j - 1, i)]
            m3, c3 = grid[(j - 1, i - 1)]
            m4, c4 = grid[(j, i - 1)]
            c1, c2, c3, c4 = (np.array(c1), np.array(c2), np.array(c3), np.array(c4))
            aa1 = np.linalg.norm(c2 - c1)
            bb1 = np.linalg.norm(c3 - c2)
            cc1 = np.linalg.norm(c4 - c3)
            dd1 = np.linalg.norm(c1 - c4)
            p1 = np.linalg.norm(c4 - c2)
            q1 = np.linalg.norm(c3 - c1)
            area = 0.25 * np.sqrt(4 * p1 * p1 * q1 * q1 -
                                   (bb1 * bb1 + dd1 * dd1 - aa1 * aa1 - cc1 * cc1) ** 2)
            area = 0.25 * area
            for m in (m1, m2, m3, m4):
                arn[m] += area
    return un, us, ud, arn


def build_station_matching(xline, yline, zline, params, xonfs, x4nds):
    """Milestone 5: port of meshgen.f90's `setSurfaceStation` (off-fault
    nearest-grid-node matching, called once per regular node, BEFORE
    `createMasterNode` in the main ix/iz/iy loop) and the on-fault station
    match embedded in `createMasterNode` ("setOnFaultStation", ntotft==1
    only -- matching this milestone's scope).

    `setSurfaceStation` has three mutually exclusive branches keyed on the
    node's position along x (`ix>1 and ix<nx` interior, `ix==1` left edge,
    `ix==nx` right edge) -- despite the "at surface only" comment in the
    Fortran, the branch guard is actually `iy` strictly interior (`iy>1 and
    iy<ny`), not a z/surface test; ported verbatim, not "fixed", per this
    discipline's rule against silently changing behavior it doesn't fully
    understand. Each branch does a nearest-node-in-x test (equality, OR
    strictly-closer-than-the-opposite-neighbor on whichever side(s) that
    branch has), then (independently) the same nearest-node-in-y test
    (always both-sided, since y is never a boundary edge case in this
    subroutine) -- ordered `if` cascades, first match wins, `n4yn` (a
    module-level "already matched" flag in Fortran, `matched` here) blocks
    rematching a station once found, exactly reproducing Fortran's `exit`
    -after-first-match / do-loop-order semantics.

    On-fault matching is a single exact-coordinate (within `tol`) test
    against `xonfs(1,:,1)` (along-strike x) and `xonfs(2,:,1)` (along-dip
    z), tried once per fault node in fault-encounter order (same order as
    `build_node_coordinates`'s nsmp), first match wins.

    Verified against `testsys/parity/fixtures/test_tpv8_serial/
    pydump_stations.txt` (dumped by src/pydump.f90's new station-matching
    block) via testsys/parity/test_standalone_meshgen.py: numOfOnFaultStCount
    and numOfOffFaultStCount match exactly, and every matched (station,node)
    pair matches Fortran's `anonfs`/`OffFaultStNodeIdIndex` in the SAME
    match order (not just as a set) -- Fortran's arrays are indexed by
    running match count, so order is part of the contract.

    params: same dict as build_grid_lines.
    xonfs: (2, nonfs) float array, columns [x_strike, z_dip] in meters
        (already *1000 from km, matching readstations2's conversion),
        1 fault only (ntotft==1).
    x4nds: (3, n_off) float array, rows [x, y, z] in meters, off-fault
        station coordinates.

    Returns (anonfs, off_fault_matches): anonfs is a list of (fault_seq,
    station_col, iFault=1) tuples, 1-indexed fault_seq/station_col, in
    match order; off_fault_matches is a list of (station_col, node_id)
    tuples, 1-indexed, in match order (node_id is the regular/slave grid
    node id, matching Fortran's `nodeCount` -- off-fault stations never
    match a master/split node since setSurfaceStation runs on the slave
    node before createMasterNode replaces anything).
    """
    p = params
    tol = p['tol']
    nx, ny, nz = len(xline), len(yline), len(zline)
    n_off = x4nds.shape[1]
    n_onf = xonfs.shape[1]
    matched = np.zeros(n_off + 1, dtype=bool)  # 1-indexed

    anonfs = []
    off_fault_matches = []
    node_count = 0
    fault_seq = 0
    for ix in range(nx):
        xcoor = xline[ix]
        for iz in range(nz):
            zcoor = zline[iz]
            for iy in range(ny):
                ycoor = yline[iy]
                node_count += 1

                def y_matches(i):
                    xs = x4nds[1, i - 1]
                    return (abs(ycoor - xs) < tol or
                            (iy >= 1 and xs > yline[iy - 1] and xs < ycoor and
                             (ycoor - xs) < (xs - yline[iy - 1])) or
                            (iy <= ny - 2 and xs > ycoor and xs < yline[iy + 1] and
                             (xs - ycoor) < (yline[iy + 1] - xs)))

                if 1 <= ix <= nx - 2 and 1 <= iy <= ny - 2:
                    for i in range(1, n_off + 1):
                        if matched[i]:
                            continue
                        xs = x4nds[0, i - 1]
                        if abs(zcoor - x4nds[2, i - 1]) >= tol:
                            continue
                        if not (abs(xcoor - xs) < tol or
                                (xs > xline[ix - 1] and xs < xcoor and
                                 (xcoor - xs) < (xs - xline[ix - 1])) or
                                (xs > xcoor and xs < xline[ix + 1] and
                                 (xs - xcoor) < (xline[ix + 1] - xs))):
                            continue
                        if y_matches(i):
                            matched[i] = True
                            off_fault_matches.append((i, node_count))
                            break
                elif ix == 0 and 1 <= iy <= ny - 2:
                    for i in range(1, n_off + 1):
                        if matched[i]:
                            continue
                        xs = x4nds[0, i - 1]
                        if abs(zcoor - x4nds[2, i - 1]) >= tol:
                            continue
                        if not (abs(xcoor - xs) < tol or
                                (xs > xcoor and xs < xline[ix + 1] and
                                 (xs - xcoor) < (xline[ix + 1] - xs))):
                            continue
                        if y_matches(i):
                            matched[i] = True
                            off_fault_matches.append((i, node_count))
                            break
                elif ix == nx - 1 and 1 <= iy <= ny - 2:
                    for i in range(1, n_off + 1):
                        if matched[i]:
                            continue
                        xs = x4nds[0, i - 1]
                        if abs(zcoor - x4nds[2, i - 1]) >= tol:
                            continue
                        if not (xs > xline[ix - 1] and xs < xcoor and
                                (xcoor - xs) < (xs - xline[ix - 1])):
                            continue
                        if y_matches(i):
                            matched[i] = True
                            off_fault_matches.append((i, node_count))
                            break

                if is_on_fault(xcoor, ycoor, zcoor, p['fxmin'], p['fxmax'],
                                p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], tol):
                    fault_seq += 1
                    for i in range(1, n_onf + 1):
                        if (abs(xcoor - xonfs[0, i - 1]) < tol and
                                abs(zcoor - xonfs[1, i - 1]) < tol):
                            anonfs.append((fault_seq, i, 1))
                            break
    return anonfs, off_fault_matches


# ---- TODO for the next milestones (explicitly not done here) ----
# - M7: frt.txt writer, '(1x,22e18.7e4)' Fortran exponent-width format.
# - eledet/eleshp/ss/phi (assembleGlobalKU's FEM-kernel inputs, beyond the
#   Jacobian determinant M6.5/mass_assembly.py already ports) -- needed
#   whenever a future milestone ports the stiffness assembly itself.
# native_input.py (M6) ports bGlobal/bModelGeometry/bFaultGeometry/
# bMaterial/bStations reading + netCDF4 on_fault_vars_input.nc -> fric.
# mass_assembly.py (M6.5) ports assembleGlobalMass's lumped-mass
# integration (nodalMassArr, fnms) -- standalone solver path must call
# this at runtime, NOT read pydump_fnms.txt/pydump_nodalmass.txt (those
# stay oracle-only, see M8 provenance note in this module's top docstring
# and README-parity.md).
