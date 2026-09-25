"""
Standalone (no-Fortran-in-the-loop) port of EQdyna's mesh generation for
the SERIAL (npx=npy=npz=1) case only -- src/countMeshEntities.f90 (mesh4num)
+ src/meshgen.f90's `getLocalOneDimCoorArrAndSize` + node-coordinate loop
+ createElement/setElementMaterial/replaceSlaveWithMasterNode.

Milestone 1 (node coordinates), Milestone 2 (element connectivity +
per-element material), Milestone 3 (equation numbering), Milestone 4
(split-node fault geometry: un/us/ud unit vectors + arn nodal area),
Milestone 5 (on-/off-fault station node matching), and Milestone 6
(native case-input reading -- see readInputFiles.py, which replaces every
hand-transcribed PARAMS/MATERIAL/station value M1-M5 previously used) of
the standalone-phase spec are ported and verified against
`testsys/parity/fixtures/test_tpv8_serial/pydump_meshCoor.txt` /
`pydump_conn.txt` / `pydump_nodeinfo.txt` / `pydump_fault.txt` /
`pydump_stations.txt` (the Fortran ground truth, dumped by
the removed pydump dumper (src/pydump.f90, deleted 2026-09-15)), via `the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15)` (wired
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
(see assembleGlobalMass.py) ports assembleGlobalMass.f90's lumped-mass
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

Milestone 7.5 (see assembleGlobalMass.py) closes the assembly gap this
docstring used to flag as "not yet ported": `compute_element_shape` now
ports eledet/eleshp (calcGlobalShapeFunc's full output -- spatial shape-
function derivatives dNx/dNy/dNz -- for both interior and PML elements,
since C_degen==0 means neither ever takes the wedge-degeneration branch),
`compute_hourglass` ports ss/phi (calcSSPhi4Hrgls's hourglass-control
tensors, including func_lib.f90's vlm volume formula), and `init_vel`
ports eqdyna3d.f90's fault-node velocity seeding of v1 from fric(31-36)
via the equation map (mode==1 only). Verified against
pydump_elemgeo.txt/pydump_v1.txt via the removed parity tier.

NOT YET PORTED (explicitly, not silently deferred): assembleGlobalKU (the
FEM stiffness/internal-force kernel itself), hrglss's per-step hourglass
FORCE (as opposed to the ss/phi tensors M7.5 now provides), velDispUpdate,
faulting -- i.e. everything in the actual time-stepping loop. Scoped this
way deliberately -- shipping a "mostly right" mesh generator with the
fault/station logic silently stubbed would be exactly the kind of
unverified deliverable this discipline forbids. See the module's TODO
list at the bottom for what a follow-up milestone needs next, in
dependency order.

RANK-LOCAL MESHES (python-jax-mpi, pathway item 64). Every builder below is
a function of `(xline, yline, zline)` exactly as the Fortran loop nest is a
function of its local lines, so a rank builds ITS OWN box by being handed its
local slices: `partition_1d`/`local_line` port getLocalOneDimCoorArrAndSize's
split (meshgen.f90:543-550, :588-596) and SLICE the globally accumulated line,
never re-accumulate it (restarting the geometric stretch rank-locally moves
coordinates by 1e-8..1e-7 m). What must NOT come from the local lines is the
MODEL boundary: Fortran overwrites xmin/.../zmax with the GLOBAL
modelBoundCoor (meshgen.f90:32-37), so the builders that test "is this node on
the model boundary" take `model_bound` explicitly; left None (the serial
path) it is read from the lines' own ends, which ARE the global line there --
the serial path's arithmetic is unchanged. `fault_boundary_lists` and
`mpi4arn` port createMasterNode's fltgm bookkeeping and MPI4arn.
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


# ---------------------------------------------------------------------------
# Rank-local decomposition: Fortran's own arithmetic, verbatim (rule 23).
# ---------------------------------------------------------------------------
def calc_xyz_mpi_id(me, npx, npy, npz):
    """meshgen.f90's calcXyzMPIId, verbatim: z fastest, then y, then x.
    Python rank r must hold the same subdomain Fortran rank r holds, or every
    per-rank number compares the wrong pair."""
    mex = me // (npy * npz)
    mey = (me - mex * npy * npz) // npz
    mez = me - mex * npy * npz - mey * npz
    return mex, mey, mez


def partition_1d(global_size, num_mpi, mpi_id):
    """getLocalOneDimCoorArrAndSize's split of one global 1D line, verbatim:
    local SIZE from meshgen.f90:543-550 (`<`), local OFFSET from :588-596
    (`<=`) -- two genuinely different comparisons in the Fortran, reproduced,
    not tidied. Returns (local_size, offset), offset 0-based into the global
    line.

    Neighbouring ranks OVERLAP by exactly one node plane (stride per-1, size
    per): node planes are shared, elements are not, because each rank makes
    elements only for its local ix,iy,iz >= 2 (meshgen.f90:86). That overlap
    is what MPI4NodalQuant and MPI4arn sum across. Checked here rather than
    trusted: a gap or a double overlap would be a partition defect, and it
    raises before any mesh is built."""
    if num_mpi < 1 or not 0 <= mpi_id < num_mpi:
        raise ValueError('partition_1d: mpi_id %r out of range for %r ranks'
                         % (mpi_id, num_mpi))
    per = (global_size + num_mpi - 1) // num_mpi
    resid = (global_size + num_mpi - 1) - per * num_mpi
    size = per if mpi_id < (num_mpi - resid) else per + 1
    if mpi_id <= (num_mpi - resid):
        off = (per - 1) * mpi_id
    else:
        off = (per - 1) * mpi_id + (mpi_id - num_mpi + resid)
    return size, off


def check_partition_1d(global_size, num_mpi):
    """Every rank's (size, offset) for one dimension, with the invariant the
    halo rests on ASSERTED: consecutive ranks share exactly one node plane,
    the first starts at 0, the last ends at global_size, and no rank is
    thinner than 2 nodes (a 1-node slab owns no element and would sit on both
    of its own faces, which fltgm's mod-10 coding cannot express). Raises
    naming the numbers rather than building a mesh on a broken split."""
    parts = [partition_1d(global_size, num_mpi, m) for m in range(num_mpi)]
    bad = []
    if parts[0][1] != 0:
        bad.append('rank 0 starts at %d' % parts[0][1])
    if parts[-1][1] + parts[-1][0] != global_size:
        bad.append('last rank ends at %d' % (parts[-1][1] + parts[-1][0]))
    for m in range(num_mpi - 1):
        if parts[m + 1][1] != parts[m][1] + parts[m][0] - 1:
            bad.append('ranks %d/%d do not share exactly one plane' % (m, m + 1))
    for m, (n, _o) in enumerate(parts):
        if n < 2:
            bad.append('rank %d holds %d node(s)' % (m, n))
    if bad:
        raise ValueError('partition_1d(%d nodes, %d ranks): %s -- %r'
                         % (global_size, num_mpi, '; '.join(bad), parts))
    return parts


def local_line(line, num_mpi, mpi_id):
    """This rank's slice of a GLOBAL grid line (meshgen.f90:588-596's copy
    out of globalOneDimCoorArr). A slice, never a re-accumulation: the
    returned doubles are the global line's own, bit for bit. Returns
    (local_line, offset)."""
    check_partition_1d(len(line), num_mpi)
    n, off = partition_1d(len(line), num_mpi, mpi_id)
    return np.array(line[off:off + n]), off


def fault_boundary_lists(nsmp, nx, ny, nz):
    """createMasterNode's fltgm bookkeeping (meshgen.f90:877-900) and
    MPI4arn's fltl..fltu fill (:271-300), for ONE fault (ntotft==1).

    fltgm codes which LOCAL faces a split node sits on (+1 ix==1, +2 ix==nx,
    +10 iy==1, +20 iy==ny, +100 iz==1, +200 iz==nz) and the six lists are
    decoded from it with the Fortran's own mod arithmetic. Returns a list of
    six int64 arrays, k = 0..5 for Fortran's k = 1..6 (x-, x+, y-, y+, z-,
    z+), each the ascending 1-based local fault-node indices on that face.
    `nsmp` is the (nftnd, 2) 1-based [slave, master] table in fault-encounter
    order, so row j IS local fault node j+1."""
    s = np.asarray(nsmp[:, 0], dtype=np.int64) - 1
    ix = s // (nz * ny)
    iz = (s % (nz * ny)) // ny
    iy = s % ny
    fltgm = ((ix == 0) * 1 + (ix == nx - 1) * 2 + (iy == 0) * 10 + (iy == ny - 1) * 20
             + (iz == 0) * 100 + (iz == nz - 1) * 200).astype(np.int64)
    masks = (fltgm % 10 == 1, fltgm % 10 == 2,
             fltgm % 100 - fltgm % 10 == 10, fltgm % 100 - fltgm % 10 == 20,
             fltgm - fltgm % 100 == 100, fltgm - fltgm % 100 == 200)
    return [np.nonzero(m)[0].astype(np.int64) + 1 for m in masks]


def flt_mpi_flags(part, flt_lists):
    """fltMPI(1:6): syncArnBoundary runs -- and sets fltMPI(k) -- exactly when
    the dimension is divided, face k is not model boundary (bnd /= 0) and this
    rank has fault nodes on it (fltnum(k) > 0), meshgen.f90:316-372."""
    return [part.dims[k // 2] > 1 and not part.at_model_edge(k // 2, k % 2)
            and flt_lists[k].size > 0 for k in range(6)]


def mpi4arn(comm, part, arn, flt_lists, params):
    """MPI4arn + syncArnBoundary (meshgen.f90:212-433), for ONE fault.

    Walks x (left, right), y (front, back), z (down, up) in the Fortran's
    order; on every face that is interior to the model AND carries fault
    nodes on THIS rank, Sendrecv's the arn of those nodes with the face
    neighbour and adds the neighbour's values back ONLY when the fault has
    non-zero nominal extent (fltxyz, i.e. params f*min/f*max) along that
    dimension -- the DIVIDE path. Along a degenerate dimension (a vertical or
    inserted fault's y) both ranks already hold the complete tributary area,
    the neighbour's value is a duplicate, and nothing is added -- but the
    exchange still happens and fltMPI(k) is still set, because
    MPI4NodalQuant's addFaultBoundaryTerm keys on fltMPI (see the long FIX
    comment at meshgen.f90:291-315 for why skipping the call was reverted).
    Evidence for both branches: testsys/regression/test_fault_mpi_boundary_arn.py
    (DUPLICATE) and test_dipping_fault_y_split.py (DIVIDE), on the Fortran.

    `arn` is the 1-indexed (nftnd+1,) array, updated IN PLACE. Returns
    fltMPI, six bools. Tags are syncArnBoundary's (tagBase+me /
    tagBase+neighbour, tagBase 1000/2000/3000)."""
    p = params
    fext = ((p['fxmin'], p['fxmax']), (p['fymin'], p['fymax']),
            (p['fzmin'], p['fzmax']))
    flt_mpi = flt_mpi_flags(part, flt_lists)
    for d in range(3):
        for ib in (0, 1):
            k = 2 * d + ib
            if not flt_mpi[k]:
                continue
            idx = flt_lists[k]
            nb = part.neighbour(d, ib)
            tag = (1000, 2000, 3000)[d]
            send = np.ascontiguousarray(arn[idx])
            recv = np.empty_like(send)
            comm.Sendrecv(send, dest=nb, sendtag=tag + part.rank, recvbuf=recv,
                          source=nb, recvtag=tag + nb)
            if fext[d][1] != fext[d][0]:
                arn[idx] += recv
    return flt_mpi


def fault_census(xline, yline, zline, params):
    """(count, key_sum) over EVERY fault node of the GLOBAL grid, key =
    0-based global regular-grid index ix*nz*ny + iz*ny + iy. Streams one
    x-plane at a time -- O(ny*nz) memory, never the (nx,nz,ny) mask -- so a
    rank-local run can check, with two integers, that the ranks' OWNED fault
    nodes are exactly the global set (a count alone would pass a node owned
    twice plus one owned never). Same predicate, elementwise, as
    on_fault_grid_mask."""
    p = params
    X, Y, Z = np.asarray(xline), np.asarray(yline), np.asarray(zline)
    ny, nz = Y.shape[0], Z.shape[0]
    dx = p['dx'] if p['C_degen'] > 3.0 else None
    count, key_sum = 0, 0
    for ix in range(X.shape[0]):
        m = _check_is_on_fault_vec(X[ix], Y[None, :], Z[:, None], p['fxmin'],
                                   p['fxmax'], p['fymin'], p['fymax'],
                                   p['fzmin'], p['fzmax'], p['tol'],
                                   p['C_degen'], dx)
        flat = np.nonzero(np.broadcast_to(m, (nz, ny)).ravel())[0]
        count += int(flat.size)
        key_sum += int(flat.sum()) + int(flat.size) * ix * nz * ny
    return count, key_sum


def equation_census(xline, yline, zline, params, pmlb, model_bound):
    """totalNumOfEquations of the GLOBAL mesh, streamed one x-plane at a time
    (countMeshEntities.f90's per-node count: 0 if on the fixed model
    boundary else numOfDof, plus 3 per master node) -- the number a
    rank-local run's OWNED equations must sum to. It is what catches a
    boundary test that used a local line's end instead of the model's."""
    p = params
    X, Y, Z = np.asarray(xline), np.asarray(yline), np.asarray(zline)
    tol = p['tol']
    (xmin, xmax), (ymin, ymax), (zmin, _zmax) = model_bound
    fy = (np.abs(Y - ymin) < tol) | (np.abs(Y - ymax) < tol)
    fz = np.abs(Z - zmin) < tol
    py = (Y > pmlb['ymax0']) | (Y < pmlb['ymin0'])
    pz = Z < pmlb['zmin0']
    dx = p['dx'] if p['C_degen'] > 3.0 else None
    total = 0
    for ix in range(X.shape[0]):
        x = X[ix]
        fx = (abs(x - xmin) < tol) or (abs(x - xmax) < tol)
        px = (x > pmlb['xmax0']) or (x < pmlb['xmin0'])
        fixed = fx | fz[:, None] | fy[None, :]
        ndpn = np.where(px | pz[:, None] | py[None, :], 12, 3)
        m = _check_is_on_fault_vec(x, Y[None, :], Z[:, None], p['fxmin'],
                                   p['fxmax'], p['fymin'], p['fymax'],
                                   p['fzmin'], p['fzmax'], tol, p['C_degen'], dx)
        total += int(np.where(fixed, 0, ndpn).sum()) + 3 * int(np.broadcast_to(m, fixed.shape).sum())
    return total


def _dip_plane_distance(y, z, c_degen):
    """checkIsOnFault/wedge()'s shared dipping-plane distance formula
    (meshgen.f90:825-826, library_degeneration.f90:10-11):
    `abs(z + y*tan(c_degen deg)) / sqrt(1+tan(c_degen deg)**2)`. `c_degen`
    is in degrees, matching the Fortran's own `C_degen/180.d0*pi` (division
    before multiplication, reproduced in that order here even though it is
    exactly commutative in double precision). Works elementwise on Python
    floats or numpy arrays."""
    tangent = np.tan((c_degen / 180.0) * np.pi)
    return np.abs(z + y * tangent) / (1.0 + tangent ** 2) ** 0.5


def _check_is_on_fault_vec(x, y, z, fxmin, fxmax, fymin, fymax, fzmin, fzmax,
                            tol, c_degen, dx=None):
    """Elementwise (broadcasting) port of checkIsOnFault -- the same two
    branches `is_on_fault` implements, written so it also works on full
    grid-shaped arrays (used by `on_fault_grid_mask`) and on point arrays
    (used by `build_elements`' type-13 retag test), not just Python floats."""
    in_box = ((x >= fxmin - tol) & (x <= fxmax + tol) &
              (y >= fymin - tol) & (y <= fymax + tol) &
              (z >= fzmin - tol) & (z <= fzmax + tol))
    if c_degen == 0.0:
        return in_box & (y == 0.0)
    if c_degen > 3.0:
        if dx is None:
            raise ValueError('_check_is_on_fault_vec: dx is required for the C_degen>3 branch')
        return in_box & (_dip_plane_distance(y, z, c_degen) < dx / 100.0)
    raise NotImplementedError(
        '_check_is_on_fault_vec: only c_degen==0 or c_degen>3 are ported (got %r); '
        'the Fortran checkIsOnFault itself takes neither if/elseif branch for '
        '0<C_degen<=3, so isOnFault stays 0 unconditionally there -- not '
        'silently mimicked here without a real case to pin it down' % c_degen)


def is_on_fault(x, y, z, fxmin, fxmax, fymin, fymax, fzmin, fzmax, tol, c_degen=0.0, dx=None):
    """Port of checkIsOnFault. C_degen==0: planar fault at y=0 (tpv8/
    tpv104). C_degen>3: distance-to-dipping-plane test (meshgen.f90:823-828;
    tpv36/tpv37's wedge-degeneration branch) -- `dx` is required then (the
    Fortran uses `dx/100.d0`, NOT the module `tol`=1e-5 the box test uses).
    Any other c_degen (0<c_degen<=3) raises: the Fortran itself falls into
    neither if/elseif branch there, so isOnFault is unconditionally 0 -- a
    behavior this port does not silently reproduce without a real case."""
    return bool(_check_is_on_fault_vec(x, y, z, fxmin, fxmax, fymin, fymax,
                                        fzmin, fzmax, tol, c_degen, dx))


def on_fault_grid_mask(xline, yline, zline, params):
    """Vectorized `is_on_fault` over the whole (ix, iz, iy) grid, in the
    traversal order every builder in this module uses (ix outer, iz middle,
    iy inner) -- now dispatching on `params['C_degen']` via the SAME
    `_check_is_on_fault_vec` helper `is_on_fault` uses (evaluated once for
    the grid instead of once per node per builder: the scalar helper was
    being called nx*ny*nz times by each of five builders).

    For C_degen==0 this is unchanged from before (same formula, just
    restructured through the shared helper) -- tpv8/tpv104/tpv10/drv.a6 all
    carry C_degen==0 (confirmed by reading their case_input/*/
    user_defined_params.py directly, not assumed), so nothing here can
    regress them. For C_degen>3 (tpv36/tpv37) the dipping-plane distance
    test genuinely depends on BOTH y and z jointly, so it is evaluated as
    one broadcast 3-D array rather than kept axis-separable the way the
    C_degen==0 y==0 test could be.

    Returns a (nx*nz*ny,) bool array in traversal order.
    """
    p = params
    tol = p['tol']
    c_degen = p['C_degen']
    X, Y, Z = np.asarray(xline), np.asarray(yline), np.asarray(zline)
    Xg = X[:, None, None]
    Yg = Y[None, None, :]
    Zg = Z[None, :, None]
    dx = p['dx'] if c_degen > 3.0 else None
    # already the full (nx,nz,ny) shape by construction: every branch above
    # combines an x-dependent, a y-dependent and a z-dependent boolean with
    # `&`, and numpy broadcasts that to the full outer-product shape before
    # this function returns.
    mask = _check_is_on_fault_vec(Xg, Yg, Zg, p['fxmin'], p['fxmax'], p['fymin'],
                                   p['fymax'], p['fzmin'], p['fzmax'], tol, c_degen, dx)
    return mask.ravel()


def _fortran_nint(x):
    """Fortran nint(): round-half-away-from-zero. Same formula as
    readInputFiles.py's helper of the same purpose (kept local here to avoid
    a cross-module dependency for a two-line function)."""
    return int(np.sign(x) * np.floor(np.abs(x) + 0.5)) if x != 0 else 0


def insert_fault_interface(x, y, z, rough, dx, dz, ymin, ymax, tol):
    """Port of func_lib.f90's insertFaultInterface (Milestone 9): given a
    node's UNDISTORTED (x, y, z) -- the same coordinates checkIsOnFault
    tests against, per the Fortran calling this from meshgen's loop BEFORE
    checkIsOnFault runs on the same nodeCoor -- looks up the local fault-
    surface height `peak` and its along-strike/along-dip slopes (pfx, pfz)
    from `rough` (readInputFiles.read_fault_rough_geometry's dict), then
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


def build_node_coordinates(xline, yline, zline, params, model_bound=None):
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
    # insertFaultInterface blends against the MODEL's ymin/ymax (the global
    # modelBoundCoor, meshgen.f90:34-35), never a rank's local line ends.
    ymin_m, ymax_m = ((yline[0], yline[-1]) if model_bound is None
                      else model_bound[1])

    if insert_fault_type == 0:
        # PERFORMANCE: insertFaultType==0 means y_store IS ycoor for every
        # node (insert_fault_interface is never called), so the whole triple
        # loop collapses to the traversal-order outer product. Bit-identical
        # by construction: the stored values are the grid-line doubles
        # themselves, copied, with no arithmetic. The insertFaultType>0 path
        # below is left as the verbatim scalar loop because its y-morph is
        # per-node and branch-heavy; both paths are covered by
        # testsys/unit/test_meshgen_vectorized.py.
        X, Y, Z = np.asarray(xline), np.asarray(yline), np.asarray(zline)
        n_reg = nx * ny * nz
        fault = on_fault_grid_mask(X, Y, Z, p)
        nftnd = int(fault.sum())
        meshCoor = np.zeros((n_reg + nftnd + 1, 3))  # index 0 unused
        meshCoor[1:n_reg + 1, 0] = np.repeat(X, nz * ny)
        meshCoor[1:n_reg + 1, 1] = np.tile(Y, nx * nz)
        meshCoor[1:n_reg + 1, 2] = np.tile(np.repeat(Z, ny), nx)
        slave_ids = np.nonzero(fault)[0] + 1  # 1-indexed, fault-encounter order
        meshCoor[n_reg + 1:] = meshCoor[slave_ids]
        master_ids = n_reg + 1 + np.arange(nftnd, dtype=np.int64)
        nsmp = np.stack([slave_ids.astype(np.int64), master_ids], axis=1)
        return meshCoor, nftnd, nsmp

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
                        xcoor, ycoor, zcoor, rough, p['dx'], p['dz'], ymin_m, ymax_m, tol)
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
    # reshape: a rank-local box that never meets the fault has NO master node,
    # and np.array([]) is shape (0,), not (0, 3) / (0, 2).
    meshCoor[1 + len(regular):] = np.asarray(master, dtype=float).reshape(-1, 3)
    nsmp = np.array(nsmp, dtype=np.int64).reshape(-1, 2)  # (nftnd, 2)
    return meshCoor, nftnd, nsmp

def build_elements(xline, yline, zline, params, pmlb, nsmp, material, meshCoor):
    """Milestone 2: port of meshgen.f90's createElement + setElementMaterial +
    replaceSlaveWithMasterNode, single planar fault (ntotft==1), homogeneous
    or 1D-layered material (setElementMaterial's two branches), for BOTH
    C_degen==0 (planar fault, tpv8/tpv104) and C_degen>3 (dipping fault
    wedge-degeneration, tpv36/tpv37) -- the latter ports
    library_degeneration.f90's wedge()/reorder() plus meshgen.f90:91-101's
    type-13 retag; see the `c_degen > 3.0` block below and
    `_wedge_material_row`/`_splice_wedge_elements`.

    C_degen==0 verified against `testsys/parity/fixtures/test_tpv8_serial/
    pydump_conn.txt` (nodeElemIdRelation + elemTypeArr + mat(1:5) for every
    element, dumped by src/pydump.f90) -- see the removed parity tier
    (test_standalone_meshgen.py, deleted 2026-09-15). C_degen>3 verified
    against a freshly-built Fortran binary on test.tpv36 by
    `testsys/parity/evidence_c_degen_port.py` (element-type-count/material/
    fault-node-count parity -- report-only, not wired into any gate).

    PERFORMANCE NOTE (vectorized; the verbatim scalar loop this replaced is
    kept as `_build_elements_scalar` below and is the bit-for-bit oracle
    `testsys/unit/test_meshgen_c_degen_wedge.py` checks this against):

    The scalar version's plane1/plane2 sliding-column bookkeeping is
    reproducible in closed form, and this was established by reading the
    loop, not assumed. Two facts make it exact:

      (a) The inner (iz, iy) loops visit EVERY (iy<ny, iz) slot of plane2 on
          every ix pass, so after `plane1 = plane2.copy()` the regular rows
          of plane1/plane2 are exactly the running node counter's value at
          that grid point:
              plane2[iy, iz] == ix*nz*ny + iz*ny + iy + 1
              plane1[iy, iz] == (ix-1)*nz*ny + iz*ny + iy + 1
          (traversal is ix outer, iz middle, iy inner, so node ids run in
          that order).
      (b) The master-node overwrite targets ROW INDEX ny -- one row past the
          ny regular rows -- while the 8-corner `c` vector only ever reads
          rows iy-1 and iy with iy <= ny-1. Row ny is therefore WRITTEN and
          NEVER READ anywhere in this function; the sliding-column timing
          the original docstring flagged as load-bearing has no observable
          effect on `conn`. Fault-adjacent corner substitution is done
          entirely by replaceSlaveWithMasterNode's slave->master lookup
          below, which is order-independent. (Kept as an explicit note so a
          future reader who adds a row-ny read knows this derivation stops
          holding.)

    Element order is unchanged: ix outer, iz middle, iy inner, over
    ix,iz,iy >= 1 -- elemIndex = (ix-1)*(nz-1)*(ny-1) + (iz-1)*(ny-1) + (iy-1).

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
    0..E-1), elem_type is (E,) int (1=interior/hourglass-controlled, 2=PML,
    and for C_degen>3 also 11=wedge below fault/12=wedge above fault/
    13=interior brick retagged adjacent-to-fault -- E itself is then LARGER
    than (nx-1)*(ny-1)*(nz-1) by the number of wedge-triggered brick
    positions, one extra element per trigger),
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
    fymin, fymax, fzmax = p['fymin'], p['fymax'], p['fzmax']
    c_degen = p['C_degen']
    if not (c_degen == 0.0 or c_degen > 3.0):
        raise NotImplementedError('build_elements: only C_degen==0 or C_degen>3 is ported '
                                   '(got %r)' % c_degen)

    # PMLb per getLocalOneDimCoorArrAndSize's mapping (build_grid_lines' pmlb dict).
    xmax0, xmin0 = pmlb['xmax0'], pmlb['xmin0']
    ymax0, ymin0 = pmlb['ymax0'], pmlb['ymin0']
    zmin0 = pmlb['zmin0']

    n_elem = (nx - 1) * (ny - 1) * (nz - 1)

    # ---- corner connectivity, from the closed form derived in the docstring ----
    # (IX,IZ,IY) in element order: ix outer, iz middle, iy inner, all >= 1.
    IX, IZ, IY = np.meshgrid(np.arange(1, nx, dtype=np.int64),
                              np.arange(1, nz, dtype=np.int64),
                              np.arange(1, ny, dtype=np.int64), indexing='ij')
    IX, IZ, IY = IX.ravel(), IZ.ravel(), IY.ravel()

    def nid(ixv, izv, iyv):
        return ixv * (nz * ny) + izv * ny + iyv + 1

    conn = np.empty((n_elem, 8), dtype=np.int64)
    conn[:, 0] = nid(IX - 1, IZ - 1, IY - 1)
    conn[:, 1] = nid(IX, IZ - 1, IY - 1)
    conn[:, 2] = nid(IX, IZ - 1, IY)
    conn[:, 3] = nid(IX - 1, IZ - 1, IY)
    conn[:, 4] = nid(IX - 1, IZ, IY - 1)
    conn[:, 5] = nid(IX, IZ, IY - 1)
    conn[:, 6] = nid(IX, IZ, IY)
    conn[:, 7] = nid(IX - 1, IZ, IY)

    # ---- element centers: mean of the 8 corner coords, from meshCoor ----
    # (for insertFaultType>0 the y stored in meshCoor is the MORPHED value,
    # so this MUST read meshCoor, not the grid lines -- see below.)
    centers = meshCoor[conn].mean(axis=1)  # (E,3)
    cx, cy, cz = centers[:, 0], centers[:, 1], centers[:, 2]

    on_bound = ((cx == xmax0) | (cx == xmin0) | (cy == ymax0) | (cy == ymin0) |
                (cz == zmin0))
    if on_bound.any():
        b = int(np.argmax(on_bound))
        raise ValueError('checkPMLAlignment: element center exactly on a '
                          'PML bound at (%r,%r,%r) (element %d)'
                          % (cx[b], cy[b], cz[b], b))

    elem_type = np.where((cx > xmax0) | (cx < xmin0) | (cy > ymax0) | (cy < ymin0) |
                          (cz < zmin0), 2, 1).astype(np.int64)

    # ---- wedge degeneration (C_degen>3): library_degeneration.f90's wedge()
    # + meshgen.f90:91-101's type-13 retag, both gathered from the ORIGINAL
    # (pre-replaceSlaveWithMasterNode) `conn` -- Fortran's wedge() reads
    # plane1/plane2 (== this `conn`) before replaceSlaveWithMasterNode ever
    # runs on that ix/iy/iz point. tpv8/tpv104/tpv10/drv.a6 all carry
    # C_degen==0 (confirmed by reading their case_input/*/
    # user_defined_params.py), so this entire block is unreachable for them.
    wedge_trigger = np.zeros(n_elem, dtype=bool)
    retag13 = np.zeros(n_elem, dtype=bool)
    wedge11_conn = wedge12_conn = wedge12_mat_row = None
    if c_degen > 3.0:
        # wedge()'s OWN box test (library_degeneration.f90:12-15): strict
        # inequalities, NO +/-tol slop on the bounds (unlike checkIsOnFault's
        # box test below), and cenz bounded ONLY from below -- reproduced
        # exactly, not "fixed" to match checkIsOnFault's box test.
        dist_wedge = _dip_plane_distance(cy, cz, c_degen)
        wedge_trigger = ((cx > fxmin) & (cx < fxmax) & (cy > fymin) & (cy < fymax) &
                          (cz > fzmin) & (dist_wedge < tol))

        # meshgen.f90:93-99's retag-to-13 test: checkIsOnFault on THIS
        # element's (unmodified, since wedge didn't fire) corner1/corner2 --
        # only meaningful where wedge did NOT fire and the element is still
        # plain interior (elemTypeArr==1); PML (2) is never retagged.
        n1 = meshCoor[conn[:, 0]]
        n2 = meshCoor[conn[:, 1]]
        onfault1 = _check_is_on_fault_vec(n1[:, 0], n1[:, 1], n1[:, 2], fxmin, fxmax,
                                           fymin, fymax, fzmin, fzmax, tol, c_degen, dx)
        onfault2 = _check_is_on_fault_vec(n2[:, 0], n2[:, 1], n2[:, 2], fxmin, fxmax,
                                           fymin, fymax, fzmin, fzmax, tol, c_degen, dx)
        retag13 = (~wedge_trigger) & (elem_type == 1) & (onfault1 | onfault2)

        # library_degeneration.f90:30-43's two `reorder` calls, vectorized:
        # neworder=(5,1,4,4,6,2,3,3) for the type-11 (below-fault) wedge,
        # (4,8,5,5,3,7,6,6) for type-12 (above-fault) -- corner k is
        # `conn[:, k-1]`. Gathered from the conn snapshot ABOVE, i.e. before
        # `replace` (built below) can mutate it in place.
        wedge11_conn = conn[:, [4, 0, 3, 3, 5, 1, 2, 2]]
        wedge12_conn = conn[:, [3, 7, 4, 4, 2, 6, 5, 5]]
        wedge12_mat_row = _wedge_material_row(material)

    elem_type = np.where(retag13, 13, elem_type)

    # ---- replaceSlaveWithMasterNode (meshgen.f90:733-736's full guard:
    # (elemTypeArr==1 AND geometric) OR elemTypeArr==12 OR elemTypeArr==13.
    # The OR-clause is unreachable dead code for C_degen==0 (elemTypeArr
    # never takes 11/12/13 then) and is exercised here for the first time;
    # elemTypeArr==12 doesn't exist at this per-brick-row granularity yet
    # (it is only created by the splice below) so it is applied directly to
    # `wedge12_conn` there instead -- UNCONDITIONALLY, matching the Fortran;
    # elemTypeArr==11 NEVER gets this substitution, also matching the
    # Fortran, and is intentionally absent from both masks below.)
    # Test uses the element's "top" node coords (this ix,iy,iz), matching
    # Fortran's `nodeCoor` at the point createElement/replaceSlave... run.
    xcoor, ycoor, zcoor = xline[IX], yline[IY], zline[IZ]
    replace = (((elem_type == 1) & (xcoor > fxmin - tol) & (xcoor < fxmax + dx + tol) &
                (zcoor > fzmin - tol) & (ycoor > 0.0) & (np.abs(ycoor - dy) < tol)) |
               (elem_type == 13))
    lut = np.arange(meshCoor.shape[0], dtype=np.int64)
    lut[nsmp[:, 0]] = nsmp[:, 1]
    if replace.any():
        conn[replace] = lut[conn[replace]]
    if wedge12_conn is not None:
        wedge12_conn = lut[wedge12_conn]

    # ---- setElementMaterial ----
    nmat, n2mat = material.shape
    mat = np.empty((n_elem, 5))
    if nmat == 1 and n2mat == 3:
        vp, vs, rho = material[0, 0], material[0, 1], material[0, 2]
        mu = vs * vs * rho
        lam = vp * vp * rho - 2.0 * mu
        mat[:] = np.array([vp, vs, rho, lam, mu])
    elif nmat > 1 and n2mat == 4:
        edges = material[:, 0]
        if np.any(np.diff(edges) <= 0.0):
            raise ValueError('setElementMaterial: bMaterial.txt layer bottoms '
                              'are not strictly ascending (%r) -- the scalar '
                              'first-match-wins layer search this replaces is '
                              'only reproducible for an ascending table' % (edges.tolist(),))
        d = np.abs(cz)
        # scalar equivalent: row 0 if d < edges[0]; else the unique i>=1 with
        # edges[i-1] <= d < edges[i]; else (d >= edges[-1]) a loud failure.
        row = np.searchsorted(edges, d, side='right')
        bad = np.nonzero(row >= nmat)[0]
        if bad.size:
            raise ValueError('setElementMaterial: depth %r (element %d) not covered '
                              'by any material layer' % (d[bad[0]], int(bad[0])))
        vp, vs, rho = material[row, 1], material[row, 2], material[row, 3]
        mu = vs * vs * rho
        lam = vp * vp * rho - 2.0 * mu
        mat[:, 0], mat[:, 1], mat[:, 2], mat[:, 3], mat[:, 4] = vp, vs, rho, lam, mu
    else:
        raise NotImplementedError('setElementMaterial: only nmat==1/n2mat==3 '
                                   '(homogeneous) or nmat>1/n2mat==4 (1D '
                                   'layered) branches are ported')

    # meshgen.f90:103 `setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0,
    # elemCount)` -- same expression, same operand order, evaluated per element.
    depth = -0.5 * (zline[IZ] + zline[IZ - 1]) + 7.3215

    if wedge_trigger.any():
        conn, elem_type, mat, depth = _splice_wedge_elements(
            conn, elem_type, mat, depth, wedge_trigger, wedge11_conn,
            wedge12_conn, wedge12_mat_row)

    return conn, elem_type, mat, depth


def _wedge_material_row(material):
    """Port of library_degeneration.f90's wedge()'s hardcoded material read
    for the type-12 (above-fault) wedge sub-element (lines 49-53):
    `mat(elemCount,1)=material(1,1)`, `,2)=material(1,2)`, `,3)=material(1,3)`,
    then `mu=vs**2*rho`, `lam=vp**2*rho-2*mu`. This reads material ROW 1,
    COLUMNS 1-3 LITERALLY, regardless of `n2mat` -- for the homogeneous
    (n2mat==3) table those columns really are [vp,vs,rho] (tpv36's case,
    matching setElementMaterial's own homogeneous branch exactly); for the
    1D-layered (n2mat==4) table, column 1 is actually the layer's DEPTH
    bound, not vp -- reproduced index-for-index as the Fortran source
    itself does, not "fixed", per this port's discipline against silently
    correcting behavior it doesn't own."""
    vp0, vs0, rho0 = material[0, 0], material[0, 1], material[0, 2]
    mu0 = vs0 * vs0 * rho0
    lam0 = vp0 * vp0 * rho0 - 2.0 * mu0
    return np.array([vp0, vs0, rho0, lam0, mu0])


def _splice_wedge_elements(conn, elem_type, mat, depth, wedge_trigger,
                            wedge11_conn, wedge12_conn, wedge12_mat_row):
    """Insert the two wedge sub-elements (types 11, 12) in place of each
    wedge_trigger brick position, shifting every later element's final
    index -- the array-level mirror of library_degeneration.f90's wedge()
    incrementing Fortran's running `elemCount` by 2 (instead of createElement's
    usual 1) for that ix/iy/iz grid point.

    `conn`/`elem_type`/`mat`/`depth` are the size-(n_brick) per-brick-position
    arrays already finalized (post material/depth/replace) for the
    non-degenerating path. The type-11 row REUSES that position's own
    mat/depth verbatim: wedge() never overwrites `mat(elemCount,:)` for the
    FIRST sub-element (only the type-12 slot gets an explicit mat write, see
    `_wedge_material_row`) and meshgen.f90:104's setPlasticStress call runs
    once per ix/iy/iz point using whatever `elemCount` is AFTER the wedge
    split (i.e. the type-12 slot) -- both facts mean this port's depth value
    for the type-11 row is not the literal Fortran value (which is simply
    never written by any call for that slot) but is harmless BY
    CONSTRUCTION: `depth` is only ever read when C_elastic==0, a combination
    this port's caller (eqdyna3d.py) refuses outright for C_degen>3.

    Returns (conn, elem_type, mat, depth) at the FINAL (post-split) element
    count.
    """
    n_brick = conn.shape[0]
    mult = np.where(wedge_trigger, 2, 1)
    starts = np.zeros(n_brick, dtype=np.int64)
    np.cumsum(mult[:-1], out=starts[1:])
    n_final = int(starts[-1] + mult[-1]) if n_brick else 0

    final_conn = np.repeat(conn, mult, axis=0)
    final_elem_type = np.repeat(elem_type, mult)
    final_mat = np.repeat(mat, mult, axis=0)
    final_depth = np.repeat(depth, mult)

    trig_idx = np.nonzero(wedge_trigger)[0]
    s = starts[trig_idx]
    final_conn[s] = wedge11_conn[trig_idx]
    final_elem_type[s] = 11
    # mat/depth at row s already correct: np.repeat carried the original
    # brick's own values forward, exactly what the type-11 sub-element keeps.
    final_conn[s + 1] = wedge12_conn[trig_idx]
    final_elem_type[s + 1] = 12
    final_mat[s + 1] = wedge12_mat_row
    # depth at row s+1: same formula, same ix/iy/iz point, already correct
    # via np.repeat too (see this function's docstring for why this is the
    # right call even though it is not the literal, never-written Fortran
    # value for the OTHER (type-11) slot).

    assert final_conn.shape[0] == n_final, (final_conn.shape[0], n_final)
    # meshgen.f90's own sanity check (assembleGlobalMass.f90:45-50): a
    # degenerate wedge's nodes 3 and 4 must coincide by construction of
    # neworder=(...,4,4,...)/(...,5,5,...) above -- verified, not assumed.
    bad11 = np.nonzero(final_conn[final_elem_type == 11][:, 2] !=
                        final_conn[final_elem_type == 11][:, 3])[0]
    if bad11.size:
        raise ValueError('_splice_wedge_elements: a type-11 wedge has unequal '
                          'node ids at corners 3/4 (%r)' % bad11[:5].tolist())
    return final_conn, final_elem_type, final_mat, final_depth


def _build_elements_scalar(xline, yline, zline, params, pmlb, nsmp, material, meshCoor):
    """The original verbatim scalar port of createElement/setElementMaterial/
    replaceSlaveWithMasterNode, kept as the bit-for-bit ORACLE for
    `build_elements`' vectorization (see that function's docstring). Not used
    on any production path -- `testsys/unit/test_meshgen_vectorized.py`
    asserts byte-equality of the two on a real mesh.
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    dx, dy = p['dx'], p['dy']
    fxmin, fxmax, fzmin = p['fxmin'], p['fxmax'], p['fzmin']
    fymin, fymax, fzmax = p['fymin'], p['fymax'], p['fzmax']
    c_degen = p['C_degen']
    if not (c_degen == 0.0 or c_degen > 3.0):
        raise NotImplementedError('_build_elements_scalar: only C_degen==0 or C_degen>3 '
                                   'is ported (got %r)' % c_degen)
    xmax0, xmin0 = pmlb['xmax0'], pmlb['xmin0']
    ymax0, ymin0 = pmlb['ymax0'], pmlb['ymin0']
    zmin0 = pmlb['zmin0']

    slave2master = {int(s): int(m) for s, m in nsmp}
    wedge12_mat_row = _wedge_material_row(material) if c_degen > 3.0 else None

    conn_rows = []       # dynamic: a wedge_trigger emits 2 rows per (ix,iy,iz)
    elem_type_rows = []
    mat_rows = []
    depth_rows = []

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

                if is_on_fault(xcoor, ycoor, zcoor, fxmin, fxmax, fymin,
                                fymax, fzmin, fzmax, tol, c_degen, dx):
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

                    corners = meshCoor[c]  # (8,3)
                    cx, cy, cz = corners.mean(axis=0)

                    if cx == xmax0 or cx == xmin0 or cy == ymax0 or cy == ymin0 or cz == zmin0:
                        raise ValueError(
                            'checkPMLAlignment: element center exactly on a '
                            'PML bound at (%r,%r,%r)' % (cx, cy, cz))

                    etype = 1
                    if cx > xmax0 or cx < xmin0 or cy > ymax0 or cy < ymin0 or cz < zmin0:
                        etype = 2

                    depth_val = -0.5 * (zline[iz] + zline[iz - 1]) + 7.3215

                    # ---- wedge degeneration (C_degen>3): verbatim scalar
                    # mirror of library_degeneration.f90's wedge() +
                    # meshgen.f90:91-101's type-13 retag. This runs BEFORE
                    # replaceSlaveWithMasterNode below, on the UNSUBSTITUTED
                    # `c`, exactly as meshgen.f90 orders the calls.
                    wedge_fired = False
                    if c_degen > 3.0:
                        dist_wedge = abs(cy * np.tan((c_degen / 180.0) * np.pi) + cz) / \
                            (1.0 + np.tan((c_degen / 180.0) * np.pi) ** 2) ** 0.5
                        if (cx > fxmin and cx < fxmax and cy > fymin and cy < fymax and
                                cz > fzmin and dist_wedge < tol):
                            wedge_fired = True
                            c11 = c[[4, 0, 3, 3, 5, 1, 2, 2]]
                            c12 = c[[3, 7, 4, 4, 2, 6, 5, 5]]
                            conn_rows.append(c11)
                            elem_type_rows.append(11)
                            mat_rows.append(material_for(cz))
                            depth_rows.append(depth_val)
                            c12 = np.array([slave2master.get(int(v), int(v)) for v in c12],
                                            dtype=np.int64)
                            conn_rows.append(c12)
                            elem_type_rows.append(12)
                            mat_rows.append(wedge12_mat_row)
                            depth_rows.append(depth_val)
                            elem_count += 2
                        elif etype == 1:
                            n1 = meshCoor[int(c[0])]
                            n2 = meshCoor[int(c[1])]
                            on1 = is_on_fault(n1[0], n1[1], n1[2], fxmin, fxmax, fymin,
                                               fymax, fzmin, fzmax, tol, c_degen, dx)
                            on2 = is_on_fault(n2[0], n2[1], n2[2], fxmin, fxmax, fymin,
                                               fymax, fzmin, fzmax, tol, c_degen, dx)
                            if on1 or on2:
                                etype = 13

                    if not wedge_fired:
                        if (etype == 1 and (xcoor > fxmin - tol and xcoor < fxmax + dx + tol
                                             and zcoor > fzmin - tol and ycoor > 0.0
                                             and abs(ycoor - dy) < tol)) or etype == 13:
                            for k in range(8):
                                nid = int(c[k])
                                if nid in slave2master:
                                    c[k] = slave2master[nid]

                        conn_rows.append(c)
                        elem_type_rows.append(etype)
                        mat_rows.append(material_for(cz))
                        depth_rows.append(depth_val)
                        elem_count += 1
        plane1 = plane2.copy()

    n_elem = (nx - 1) * (ny - 1) * (nz - 1)
    if c_degen == 0.0:
        assert elem_count == n_elem, (elem_count, n_elem)
    conn = np.array(conn_rows, dtype=np.int64)
    elem_type = np.array(elem_type_rows, dtype=np.int64)
    mat = np.array(mat_rows)
    depth = np.array(depth_rows)
    # meshgen.f90's own sanity check (assembleGlobalMass.f90:45-50): a
    # degenerate wedge's nodes 3 and 4 must coincide by construction.
    bad11 = np.nonzero(conn[elem_type == 11][:, 2] != conn[elem_type == 11][:, 3])[0]
    if bad11.size:
        raise ValueError('_build_elements_scalar: a type-11 wedge has unequal '
                          'node ids at corners 3/4 (%r)' % bad11[:5].tolist())
    return conn, elem_type, mat, depth


def build_equation_numbers(xline, yline, zline, params, pmlb, model_bound=None):
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
    the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15).

    Returns (num_dof, eq_start, eq_nums, total_num_of_equations) where
    num_dof is (N+1,) int (3 or 12, row 0 unused), eq_start is (N+1,) int
    (0-indexed flat-array offset before this node's block, matching
    Fortran's eqNumStartIndexLoc convention of "tag value before this
    node's writes"), eq_nums is a length-(N+1) list of 1D int arrays (row 0
    unused) holding this node's assigned equation numbers (positive) or -1
    (fixed-boundary sentinel) per dof slot, and total_num_of_equations is
    the final equationNumCount (== totalNumOfEquations).

    PERFORMANCE NOTE (vectorized; `_build_equation_numbers_scalar` below is
    the verbatim original and the bit-for-bit oracle
    `testsys/unit/test_meshgen_vectorized.py` checks this against): every
    quantity here is an EXACT INTEGER prefix sum over the traversal, so the
    vectorization carries no floating-point reassociation risk at all --
    `tag` and `equationNumCount` are running counters whose per-node
    increments (`ndpn` + 3*on_fault, and 0-if-fixed-else-`ndpn` + 3*on_fault
    respectively) depend only on that node's own coordinates. The
    interleaving the original docstring flags as load-bearing (regular
    node's slots, THEN its master node's 3 slots, before moving on) is
    reproduced by putting the master block inside the same per-grid-node
    block, not by appending masters at the end.
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
    # A RANK-LOCAL build passes `model_bound` (the GLOBAL modelBoundCoor,
    # countMeshEntities.f90:34-35): a local line's own ends would mark every
    # subdomain face fixed and weld the model shut.
    X, Y, Z = np.asarray(xline), np.asarray(yline), np.asarray(zline)
    if model_bound is None:
        xmin, xmax = X[0], X[-1]
        ymin, ymax = Y[0], Y[-1]
        zmin = Z[0]
    else:
        (xmin, xmax), (ymin, ymax), (zmin, _zmax) = model_bound
    xmax0, xmin0 = pmlb['xmax0'], pmlb['xmin0']
    ymax0, ymin0 = pmlb['ymax0'], pmlb['ymin0']
    zmin0 = pmlb['zmin0']
    ndof = 3

    N_regular = nx * ny * nz

    # ---- per-grid-node predicates, in traversal order (ix, iz, iy) ----
    # is_fixed_boundary: zmax (free surface) deliberately excluded.
    fx = (np.abs(X - xmin) < tol) | (np.abs(X - xmax) < tol)
    fy = (np.abs(Y - ymin) < tol) | (np.abs(Y - ymax) < tol)
    fz = np.abs(Z - zmin) < tol
    fixed = (fx[:, None, None] | fz[None, :, None] | fy[None, None, :]).ravel()

    # num_dof_for: 12 in the PML shell, 3 otherwise.
    px = (X > xmax0) | (X < xmin0)
    py = (Y > ymax0) | (Y < ymin0)
    pz = Z < zmin0
    is_pml = (px[:, None, None] | pz[None, :, None] | py[None, None, :]).ravel()
    ndpn = np.where(is_pml, 12, 3).astype(np.int64)

    fault = on_fault_grid_mask(X, Y, Z, p)
    nftnd = int(fault.sum())
    N_total = N_regular + nftnd

    # ---- running counters as exact integer prefix sums ----
    # tag advances by ndpn (regular block) + 3 (master block, if any).
    tags_per_node = ndpn + 3 * fault
    tag_before = np.zeros(N_regular, dtype=np.int64)
    np.cumsum(tags_per_node[:-1], out=tag_before[1:])
    total_tags = int(tag_before[-1] + tags_per_node[-1])

    # equationNumCount advances only for non-fixed regular slots, plus the
    # master node's 3 slots (master nodes are never fixed-boundary).
    eq_regular = np.where(fixed, 0, ndpn)
    eqs_per_node = eq_regular + 3 * fault
    eq_before = np.zeros(N_regular, dtype=np.int64)
    np.cumsum(eqs_per_node[:-1], out=eq_before[1:])
    eq_count = int(eq_before[-1] + eqs_per_node[-1])

    # ---- the flat eqNumIndexArr, laid out in tag order ----
    slot_node = np.repeat(np.arange(N_regular, dtype=np.int64), tags_per_node)
    offset = np.arange(total_tags, dtype=np.int64) - tag_before[slot_node]
    ndpn_s = ndpn[slot_node]
    is_master_slot = offset >= ndpn_s
    flat_slots = np.where(
        is_master_slot,
        eq_before[slot_node] + eq_regular[slot_node] + (offset - ndpn_s) + 1,
        np.where(fixed[slot_node], -1, eq_before[slot_node] + offset + 1))

    # ---- per-node-id views (node ids 1..N_regular regular, then masters) ----
    fault_idx = np.nonzero(fault)[0]
    num_dof = np.zeros(N_total + 1, dtype=np.int64)
    eq_start = np.zeros(N_total + 1, dtype=np.int64)
    num_dof[1:N_regular + 1] = ndpn
    num_dof[N_regular + 1:] = ndof
    eq_start[1:N_regular + 1] = tag_before
    eq_start[N_regular + 1:] = tag_before[fault_idx] + ndpn[fault_idx]

    starts = eq_start.tolist()
    counts = num_dof.tolist()
    eq_nums = [None] * (N_total + 1)
    for n in range(1, N_total + 1):
        s = starts[n]
        eq_nums[n] = flat_slots[s:s + counts[n]]

    return num_dof, eq_start, eq_nums, eq_count


def _build_equation_numbers_scalar(xline, yline, zline, params, pmlb, model_bound=None):
    """The original verbatim scalar port of setNumDof/setEquationNumber/
    createMasterNode's equation-number half, kept as the bit-for-bit ORACLE
    for `build_equation_numbers`' vectorization. Not used on any production
    path -- `testsys/unit/test_meshgen_vectorized.py` asserts equality of the
    two on a real mesh.
    """
    p = params
    nx, ny, nz = len(xline), len(yline), len(zline)
    tol = p['tol']
    if model_bound is None:
        xmin, xmax = xline[0], xline[-1]
        ymin, ymax = yline[0], yline[-1]
        zmin = zline[0]
    else:
        (xmin, xmax), (ymin, ymax), (zmin, _zmax) = model_bound
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

    tag = 0
    eq_count = 0
    node_count = 0
    master_count = 0
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


def pack_eq_ids(num_dof, eq_nums, n_nodes, ncols=12):
    """Pack `build_equation_numbers`' per-node equation-number arrays into
    the dense (n_nodes, ncols) 0-indexed table `loading.load()` expects,
    with the -1 fixed-boundary sentinel mapped to the 0 sink slot.

    Exactly the loop main.py's build_solver_state used to run per node:
        eqs = np.array(eq_nums[node]); eq_ids[node-1, :nd] = where(eqs>0, eqs, 0)
    All-integer, so there is nothing to reassociate -- this is a layout
    change only. Returns (ndof0, eq_ids), both 0-indexed by node.
    """
    nd = np.asarray(num_dof[1:n_nodes + 1], dtype=np.int64)
    if nd.size != n_nodes:
        raise ValueError('pack_eq_ids: num_dof holds %d node entries, expected %d'
                          % (nd.size, n_nodes))
    if nd.size and int(nd.max()) > ncols:
        raise ValueError('pack_eq_ids: node with %d dof exceeds ncols=%d'
                          % (int(nd.max()), ncols))
    eq_flat = np.concatenate(eq_nums[1:n_nodes + 1]).astype(np.int64, copy=False)
    if eq_flat.size != int(nd.sum()):
        raise ValueError('pack_eq_ids: eq_nums slot count %d disagrees with '
                          'num_dof total %d' % (eq_flat.size, int(nd.sum())))
    row = np.repeat(np.arange(n_nodes, dtype=np.int64), nd)
    starts = np.zeros(n_nodes, dtype=np.int64)
    np.cumsum(nd[:-1], out=starts[1:])
    col = np.arange(eq_flat.size, dtype=np.int64) - starts[row]
    eq_ids = np.zeros((n_nodes, ncols), dtype=np.int64)
    eq_ids[row, col] = np.where(eq_flat > 0, eq_flat, 0)
    return nd, eq_ids


def build_fault_geometry(xline, yline, zline, params, nsmp, model_bound=None):
    """Milestone 4 (+ Milestone 9's insertFaultType>0 extension): port of
    meshgen.f90's split-node unit-vector assignment (`createMasterNode`'s
    un/us/ud writes) and the on-fault quadrilateral-area accumulation onto
    `arn` (meshgen's main-loop area block after the ix/iy/iz node loop),
    for a single (ntotft==1) fully rectangular fault -- planar
    (C_degen<=3, insertFaultType==0: angle-based un/us/ud, tpv8/tpv104's
    branch) OR dipping/rough (insertFaultType>0: pfx/pfz-derived un/us/ud
    from func_lib.f90's insertFaultInterface via meshgen.f90:925-938,
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
    the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15).

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

    ymin_m, ymax_m = ((yline[0], yline[-1]) if model_bound is None
                      else model_bound[1])

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
    # PERFORMANCE: the original walked all nx*ny*nz grid points calling the
    # scalar `is_on_fault` just to reach the ~1e3 fault nodes. The mask gives
    # the same nodes in the same traversal order (ix outer, iz middle, iy
    # inner -- the flat index IS that order), so `seq` still lines up 1:1
    # with nsmp's rows; the per-fault-node body below is unchanged.
    fault_flat = np.nonzero(on_fault_grid_mask(xline, yline, zline, p))[0]
    ix_of = fault_flat // (nz * ny)
    iz_of = (fault_flat % (nz * ny)) // ny
    iy_of = fault_flat % ny

    grid = {}
    ixfi = izfi = None
    seq = 0
    for _f in range(fault_flat.size):
        ix, iz, iy = int(ix_of[_f]), int(iz_of[_f]), int(iy_of[_f])
        xcoor, zcoor, ycoor = xline[ix], zline[iz], yline[iy]
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
            # insertFaultType>0: createMasterNode's un/us/ud branch
            # (meshgen.f90:925-938) OVERWRITES the angle-based un/us/ud
            # above with pfx/pfz-derived values, per fault node -- and the
            # fault's actual (warped) y-coordinate is `peak`, not 0, which
            # matters for arn's corner-distance formula below (meshCoor
            # already stores this same `peak` value, per
            # build_node_coordinates' y-morph of fault nodes -- recomputed
            # here rather than re-reading meshCoor, since the fault-node
            # traversal order is independently re-walked in every M-builder).
            y_geo, pfx, pfz = insert_fault_interface(
                xcoor, ycoor, zcoor, rough, p['dx'], p['dz'],
                ymin_m, ymax_m, tol)
            denom = (pfx ** 2 + 1.0 + pfz ** 2) ** 0.5
            un[seq] = (-pfx / denom, 1.0 / denom, -pfz / denom)
            us_denom = (1.0 + pfx ** 2) ** 0.5
            us[seq] = (1.0 / us_denom, pfx / us_denom, 0.0)
            ud[seq] = np.cross(us[seq], un[seq])
        grid[(ifs, ifd)] = (seq, (xcoor, y_geo, zcoor))
    assert seq == nftnd, (seq, nftnd)
    if nftnd == 0:
        # meshgen.f90:115 `if(nftnd0(ift)>0)`: a rank whose box never meets
        # the fault has no quad grid and no area to accumulate.
        return un, us, ud, np.zeros(1)
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

    Row 94 (owner ruling 2026-09-24): the DEPTH test in each branch used to
    require `zcoor == x4nds[2,i-1]` exactly (within `tol`), so a station
    whose requested depth was not itself a grid z-plane matched no node at
    all. It now tests against `z_snap[i-1]`, the requested depth clamped to
    the nearest node of `zline` (computed once, above, before the loop) --
    x and y are unchanged, they already snap to the nearest interior node.

    On-fault matching is a single exact-coordinate (within `tol`) test
    against `xonfs(1,:,1)` (along-strike x) and `xonfs(2,:,1)` (along-dip
    z), tried once per fault node in fault-encounter order (same order as
    `build_node_coordinates`'s nsmp), first match wins.

    Verified against `testsys/parity/fixtures/test_tpv8_serial/
    pydump_stations.txt` (dumped by the removed pydump dumper's station-matching
    block) via the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15): numOfOnFaultStCount
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

    # Row 94 (owner ruling 2026-09-24): clamp each requested off-fault
    # station's depth to the nearest node of `zline` -- here the FULL
    # global z grid already (this port's serial builder holds the whole
    # domain, not a per-rank slice, unlike meshgen.f90's MPI-partitioned
    # zline), so this is a plain nearest-value search, no cross-rank
    # reduction needed. x and y are unchanged: they already snap to the
    # nearest interior node in the branches below.
    zline_arr = np.asarray(zline)
    z_snap = zline_arr[np.argmin(np.abs(zline_arr[None, :] - x4nds[2, :][:, None]), axis=1)]

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
                        if abs(zcoor - z_snap[i - 1]) >= tol:
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
                        if abs(zcoor - z_snap[i - 1]) >= tol:
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
                        if abs(zcoor - z_snap[i - 1]) >= tol:
                            continue
                        if not (xs > xline[ix - 1] and xs < xcoor and
                                (xcoor - xs) < (xs - xline[ix - 1])):
                            continue
                        if y_matches(i):
                            matched[i] = True
                            off_fault_matches.append((i, node_count))
                            break

                # row 114 fix: this call used to omit c_degen/dx entirely,
                # silently defaulting is_on_fault's c_degen=0.0 -- for
                # C_degen>3 (tpv36/tpv37's wedge-degenerate mesh) that takes
                # the WRONG branch (`y == 0.0`, never true off a dipping
                # plane) instead of the dipping-plane distance test, so
                # fault_seq never advanced and NO on-fault station ever
                # matched (numOfOnFaultStCount==0, zero faultst*.txt files,
                # against Fortran's 19). build_node_coordinates's OWN
                # C_degen>3 path (on_fault_grid_mask, used whenever
                # insertFaultType==0, true for every currently-supported
                # C_degen>3 case) already threads C_degen/dx correctly --
                # this call site is the one this milestone's own station
                # wiring newly exercises end to end, and it had not been.
                dx_for_fault = p['dx'] if p['C_degen'] > 3.0 else None
                if is_on_fault(xcoor, ycoor, zcoor, p['fxmin'], p['fxmax'],
                                p['fymin'], p['fymax'], p['fzmin'], p['fzmax'],
                                tol, p['C_degen'], dx_for_fault):
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
#   Jacobian determinant M6.5/assembleGlobalMass.py already ports) -- needed
#   whenever a future milestone ports the stiffness assembly itself.
# readInputFiles.py (M6) ports bGlobal/bModelGeometry/bFaultGeometry/
# bMaterial/bStations reading + netCDF4 on_fault_vars_input.nc -> fric.
# assembleGlobalMass.py (M6.5) ports assembleGlobalMass's lumped-mass
# integration (nodalMassArr, fnms) -- standalone solver path must call
# this at runtime, NOT read pydump_fnms.txt/pydump_nodalmass.txt (those
# stay oracle-only, see M8 provenance note in this module's top docstring
# and the removed parity tier).
