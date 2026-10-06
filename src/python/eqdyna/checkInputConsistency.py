"""checkInputConsistency.py <- src/fortran/checkInputConsistency.f90.

Refuses EXACTLY what the Fortran refuses: same condition, same message
text, same numbered exit code (src/fortran/errorCodes.f90), at the same
point in the run -- eqdyna3d.f90:104 calls checkInputConsistency right
after read_fault_rough_geometry and before calcLocalShapeFunc/meshgen (i.e.
after every bFile is read, before any mesh or solver work). eqdyna3d.py's
build_solver_state calls check() immediately after
readInputFiles.build_params returns, before its OWN scope guards
(ntotft/friclaw/npx,npy,npz/C_degen) -- none of those three checks below
need anything build_params does not already return, so this is the
earliest point on the Python side too, and it runs for BOTH run_case and
run_case_mpi since both funnel through build_solver_state.

Three checks, ported verbatim (checkInputConsistency.f90:7-19):

  1. C_elastic==0 and C_Q==1              -> ERR_CFG_Q_NEEDS_ELASTIC   (11)
  2. C_Q==1 and rat>1                      -> ERR_CFG_Q_NEEDS_UNIFORM  (12)
  3. output_plastic==1 and C_elastic!=0    -> ERR_CFG_PLASTIC_OUTPUT   (13)

C_Q is not a case-input field on EITHER side: globalvar.f90:104 declares it
`= 0` and no reader (readInputFiles.f90, scripts/case.setup) ever assigns
it from a file -- confirmed by grep, and consistent with
docs/fortran_python_correspondence.md's note on calcQAttenuationCoeff.f90
("C_Q is hardcoded 0 ... and read from no input file"). This module
hardcodes C_Q=0 for the identical reason (module-level `C_Q` below), which
makes checks 1 and 2 genuinely UNREACHABLE through any real case on BOTH
sides. They are still ported, per rule 23 ("the port refuses exactly what
Fortran refuses" -- not just the reachable subset of it), and `check()`
accepts C_Q as an explicit keyword so a future case format that DOES
surface it can exercise both checks without editing this function; see
testsys/regression/test_check_input_consistency.py for how the two
unreachable checks are pinned (direct unit calls, not a case directory) vs.
the one reachable check (output_plastic/C_elastic, real bGlobal.txt on
both binaries).
"""


import numpy as np


class InputConsistencyError(RuntimeError):
    """Raised by check(). `.code` is the matching src/fortran/errorCodes.f90
    ERR_CFG_* value, so `python3 -m eqdyna` (see eqdyna3d.main's `_abort`)
    can exit with the SAME number a Fortran run of the same bad config
    exits with (rule 23)."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


# --- 11-19 configuration and parameter consistency (errorCodes.f90:63-65) -
ERR_CFG_Q_NEEDS_ELASTIC = 11   # C_Q=1 requires C_elastic=1
ERR_CFG_Q_NEEDS_UNIFORM = 12   # C_Q=1 requires rat=1.0 (uniform elements)
ERR_CFG_PLASTIC_OUTPUT = 13    # output_plastic=1 requires C_elastic=0
ERR_CFG_NSTRESS_SIGN_INVALID = 15  # bGlobal.txt's station n-stress sign is neither +1 nor -1 (raised by readInputFiles.read_bglobal)
ERR_CFG_MATERIAL_TABLE_INVALID = 16  # bMaterial.txt's two-sided (n2mat=5) table needs one vertical planar fault, side column -1/+1, ascending per-side layer bottoms (raised by check_two_sided_material)
ERR_CFG_MATERIAL_GRID_INVALID = 17  # bMaterial.txt's 3D material grid (n2mat=6, rows x y z vp vs rho) must be a complete uniform nx*ny*nz block, every cell once, on-grid coordinates, positive vp/vs/rho (raised by build_material_grid3d)
ERR_GEOM_MULTIFAULT_Y_BAD = 32  # a fault's y-plane is not vertical/planar, coincides with another fault's, or (row 17 rebased) a fault's y bound is not an integer multiple of dy from the union-derived uniform-y belt origin (raised by meshgen.py's one_dim_coor_array)
ERR_GEOM_MULTIFAULT_XZ_BAD = 33  # (row 17 rebased) a fault's x or z bound is not an integer multiple of dx/dz from the union-derived uniform x/z belt origin -- independent per-fault x/z extents are supported, but each must land on a mesh node line (raised as InputConsistencyError by meshgen.py's one_dim_coor_array)

# C_Q is never a case-input field (see module docstring) -- hardcoded here
# exactly as globalvar.f90:104 hardcodes its default, never overridden by
# any reader on either side.
C_Q = 0


def check(C_elastic, output_plastic, rat, C_Q=C_Q):
    """checkInputConsistency.f90:7-19, verbatim, same order. Raises
    InputConsistencyError on the first failing check with the Fortran's
    message text unchanged."""
    if C_elastic == 0 and C_Q == 1:
        raise InputConsistencyError(
            ERR_CFG_Q_NEEDS_ELASTIC,
            'Q model (C_Q=1) can only work with the elastic code (C_elastic=1). '
            'Set C_Q=0 or C_elastic=1.')
    if C_Q == 1 and rat > 1:
        raise InputConsistencyError(
            ERR_CFG_Q_NEEDS_UNIFORM,
            'Q model (C_Q=1) can only work with uniform element size; rat must be 1.0.')
    if output_plastic == 1 and C_elastic != 0:
        print(' Now, C_elastic = ', C_elastic, flush=True)   # checkInputConsistency.f90:16
        raise InputConsistencyError(
            ERR_CFG_PLASTIC_OUTPUT,
            'Plastic strains are only output for C_elastic=0. Set output_plastic=0 or C_elastic=0.')


def check_multifault(faults, dy, dis4uniF, dis4uniB, C_degen, tol=1.0e-5):
    """checkInputConsistency.f90's row-17 multi-fault guards (only the
    `C_degen == 0.0` branch -- C_degen>3's dipping/wedge-degeneration
    mechanism is a different, pre-existing, single-fault-only path with a
    legitimate fymin != fymax, orthogonal to this work).
    No-op at ntotft==1 (every loop below is over a single fault, and the
    i<j distinctness loop does not execute for ntotft<2) -- bit-identical
    refusal behaviour to before this function existed.

    Row 17 REBASED (restore per-fault mesh extent, invent nothing): the "y
    outside the fixed belt" and "x/z must equal fault 1's" checks that used
    to live here are GONE, not relaxed -- they refused exactly the two
    limitations meshgen.py's one_dim_coor_array no longer has (the belt is
    now the union over every fault's own box, mirroring
    src/fortran/meshgen.f90's getLocalOneDimCoorArrAndSize/
    checkInputConsistency.f90). one_dim_coor_array itself now raises the
    commensurability guard (a fault bound not an integer multiple of
    dx/dy/dz from the union origin) at mesh-build time, same as the
    Fortran's hard refuse inside getLocalOneDimCoorArrAndSize.

    `faults`: read_bfaultgeometry's return (list of ntotft dicts with
    fxmin/fxmax/fymin/fymax/fzmin/fzmax)."""
    if C_degen != 0.0:
        return
    for i, f in enumerate(faults, start=1):
        if abs(f['fymax'] - f['fymin']) > tol:
            raise InputConsistencyError(
                ERR_GEOM_MULTIFAULT_Y_BAD,
                "checkInputConsistency: fault %d has fymin /= fymax -- only a planar, "
                "vertical fault (single y-plane) is supported; non-planar multi-fault "
                "geometry is out of scope." % i)
    for i in range(len(faults)):
        for j in range(i + 1, len(faults)):
            if abs(faults[i]['fymin'] - faults[j]['fymin']) < tol:
                raise InputConsistencyError(
                    ERR_GEOM_MULTIFAULT_Y_BAD,
                    "checkInputConsistency: faults %d and %d are both at y = %.3f -- two "
                    "faults must occupy distinct y-planes." % (i + 1, j + 1, faults[i]['fymin']))


def check_two_sided_material(material, faults, C_degen, tol=1.0e-5):
    """readInputFiles.f90's checkTwoSidedMaterialTable, same checks in the
    same order, same messages: the n2mat==5 two-sided 1D material table
    (meshgen.py build_elements' third material branch, SCEC TPV35) picks a
    layer set by the sign of (element-centre y - the fault y-plane), which
    only means something when every fault lies on ONE common vertical plane
    (`faults`: the per-fault box dicts, fymin == fymax and identical across
    faults -- the Fortran's `maxval(fltxyz(2,2,:)) - minval(fltxyz(1,2,:))`
    test); the per-side first-match lookup is only reproducible for strictly
    ascending bottoms; a side with no rows would leave every element on it
    without a material. No-op unless material has 5 columns."""
    material = np.asarray(material)
    if material.ndim != 2 or material.shape[1] != 5:
        return
    nmat = material.shape[0]
    code = ERR_CFG_MATERIAL_TABLE_INVALID
    if nmat < 2:
        raise InputConsistencyError(code, 'bMaterial.txt: a two-sided (n2mat=5) table '
                                    'needs nmat >= 2 rows (at least one per side).')
    yspan = max(f['fymax'] for f in faults) - min(f['fymin'] for f in faults)
    if yspan > tol:
        raise InputConsistencyError(code, 'bMaterial.txt: the two-sided (n2mat=5) material '
                                    'table takes sides against ONE vertical y-plane; every '
                                    'fault must have fymin == fymax and share the same y.')
    if C_degen != 0:
        raise InputConsistencyError(code, 'bMaterial.txt: the two-sided (n2mat=5) material '
                                    'table needs a vertical planar fault (C_degen=0); a '
                                    'dipping/degenerate fault has no single y-plane to take '
                                    'sides against.')
    sides = np.rint(material[:, 4]).astype(int)
    if np.any((sides != -1) & (sides != 1)):
        raise InputConsistencyError(code, 'bMaterial.txt: column 5 of a two-sided (n2mat=5) '
                                    'table must be -1 (y below the fault plane) or +1 (y above it).')
    for s in (-1, 1):
        bottoms = material[sides == s, 0]
        if bottoms.size == 0:
            raise InputConsistencyError(code, 'bMaterial.txt: a two-sided (n2mat=5) table must '
                                        'have at least one row for each side (-1 and +1).')
        if bottoms[0] <= 0.0 or np.any(np.diff(bottoms) <= 0.0):
            raise InputConsistencyError(code, 'bMaterial.txt: layer bottoms within one side of a '
                                        'two-sided (n2mat=5) table must be strictly ascending '
                                        'and positive.')


def build_material_grid3d(material, tol=1.0e-5):
    """readInputFiles.f90's buildMaterialGrid3D, same checks in the same
    order, same messages: the n2mat==6 3D structured material grid (SCEC
    TPV34, CVM-H sampled at the uniform element-centre spacing), rows
    [x y z vp vs rho] in the EQdyna frame. Self-describing: per axis the
    origin is the smallest coordinate, the spacing the smallest positive
    offset from it (1.0 for a single-plane axis), the count
    nint((max-min)/spacing)+1; nmat must equal nx*ny*nz, every cell filled
    exactly once by an on-grid row, vp/vs/rho positive. Returns None unless
    material has 6 columns; otherwise dict(origin (3,), spacing (3,),
    count (3,) int, props (3, nx, ny, nz) = vp, vs, rho). meshgen.py's
    build_elements gathers the NEAREST cell to each element centre from
    `props`, clamped -- piecewise constant, never interpolated."""
    material = np.asarray(material, dtype=float)
    if material.ndim != 2 or material.shape[1] != 6:
        return None
    nmat = material.shape[0]
    code = ERR_CFG_MATERIAL_GRID_INVALID
    if nmat < 2:
        raise InputConsistencyError(code, 'bMaterial.txt: a 3D material grid (n2mat=6) '
                                    'needs nmat >= 2 rows.')
    origin = material[:, :3].min(axis=0)
    cmax = material[:, :3].max(axis=0)
    spacing = np.empty(3)
    for k in range(3):
        off = material[:, k] - origin[k]
        pos = off[off > tol]
        spacing[k] = pos.min() if pos.size else 1.0   # single-plane axis
    count = np.rint((cmax - origin) / spacing).astype(int) + 1
    if int(np.prod(count)) != nmat:
        raise InputConsistencyError(code, 'bMaterial.txt: the 3D material grid (n2mat=6) rows '
                                    'do not form a complete uniform nx*ny*nz block '
                                    '(nmat /= nx*ny*nz).')
    off = (material[:, :3] - origin) / spacing
    idx = np.rint(off).astype(int)
    if (np.abs(off - idx).max() > 1.0e-6 or idx.min() < 0 or np.any(idx >= count)):
        raise InputConsistencyError(code, 'bMaterial.txt: a 3D material grid (n2mat=6) row has '
                                    'a coordinate that is not on the uniform grid.')
    flat = np.ravel_multi_index((idx[:, 0], idx[:, 1], idx[:, 2]), tuple(count))
    if np.unique(flat).size != nmat:
        raise InputConsistencyError(code, 'bMaterial.txt: a 3D material grid (n2mat=6) cell '
                                    'is given twice.')
    if np.any(material[:, 3:6] <= 0.0):
        raise InputConsistencyError(code, 'bMaterial.txt: a 3D material grid (n2mat=6) row has '
                                    'vp, vs or rho <= 0.')
    props = np.empty((3,) + tuple(count))
    props[:, idx[:, 0], idx[:, 1], idx[:, 2]] = material[:, 3:6].T
    return dict(origin=origin, spacing=spacing, count=count, props=props)


def material_grid3d_index(grid, cx, cy, cz):
    """The element -> grid-cell index rule of meshgen.f90 setElementMaterial's
    n2mat==6 branch: floor(off + 0.5) per axis, clamped to [0, count-1]
    (0-based here). Shared by the vectorized and scalar meshgen.py paths so
    the tie rule is written once; floor(off+0.5), not numpy.rint, so a .5
    offset rounds the same way as the Fortran."""
    c = np.stack([np.asarray(cx, dtype=float), np.asarray(cy, dtype=float),
                  np.asarray(cz, dtype=float)], axis=-1)
    off = (c - grid['origin']) / grid['spacing']
    idx = np.floor(off + 0.5).astype(int)
    return np.clip(idx, 0, grid['count'] - 1)
