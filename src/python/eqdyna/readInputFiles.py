"""
Milestone 6: native (no-Fortran-in-the-loop) readers for EQdyna's case-input
files -- ports of src/readInputFiles.f90's readglobal/readmodelgeometry/
readfaultgeometry/readmaterial/readstations1/readstations2, plus
src/netcdf_io.f90's netcdf_read_on_fault_eqdyna -- so the standalone solver
can build its geometry/friction state directly from a case directory's
bGlobal.txt/bModelGeometry.txt/bFaultGeometry.txt/bMaterial.txt/
bStations.txt/on_fault_vars_input.nc, WITHOUT running the Fortran binary
first.

Scope matches the meshgen.py milestones this module feeds: ntotft==1
(single planar fault), C_degen==0, insertFaultType==0. Multi-fault and
rough-fault (read_fault_rough_geometry) reading are explicitly NOT ported
here -- out of scope, not silently approximated.

Verified against `testsys/parity/fixtures/test_tpv8_serial/`'s
bGlobal.txt/bModelGeometry.txt/bFaultGeometry.txt/bMaterial.txt/
bStations.txt/on_fault_vars_input.nc via
the removed parity tier (test_standalone_meshgen.py, deleted 2026-09-15): the resulting PARAMS dict,
material array, xonfs/x4nds station arrays, and initial fric array are fed
straight into meshgen.py's M1-M5 builders and `read_on_fault_vars` below,
and the M1-M5 parity checks (already gated against pydump_meshCoor.txt/
pydump_conn.txt/pydump_nodeinfo.txt/pydump_fault.txt/pydump_stations.txt)
are re-run with these NATIVELY-READ inputs instead of the hand-transcribed
PARAMS dict used by earlier milestones -- same tolerances, same fixtures,
now proving the reader itself, not just the builders.
"""
import numpy as np

from .checkInputConsistency import ERR_CFG_NSTRESS_SIGN_INVALID, InputConsistencyError


def _read_records(path):
    """Yield every line of `path`, unfiltered, in file order.

    NOTHING IS SKIPPED, and that is the contract: each blank-line placeholder
    in these files is consumed by its own bare `read(1001,*)` in the Fortran,
    so the readers below advance line for line against readglobal and count
    the blanks explicitly (`next(it)  # blank`). Dropping blank or comment
    lines here -- which this docstring used to claim was happening -- would
    silently shift every field after the first placeholder."""
    with open(path) as f:
        for line in f:
            yield line


def read_bglobal(path):
    """Port of readglobal. Returns a dict of every field read from
    bGlobal.txt, keyed by its Fortran variable name."""
    lines = list(_read_records(path))
    it = iter(lines)

    def nxt():
        return next(it).split()

    g = {}
    g['mode'] = int(nxt()[0])
    g['C_elastic'] = int(nxt()[0])
    g['C_nuclea'] = int(nxt()[0])
    g['C_degen'] = float(nxt()[0])
    g['insertFaultType'] = int(nxt()[0])
    g['friclaw'] = int(nxt()[0])
    g['ntotft'] = int(nxt()[0])
    g['nucfault'] = int(nxt()[0])
    g['TPV'] = int(nxt()[0])
    g['output_plastic'] = int(nxt()[0])
    g['outputGroundMotion'] = int(nxt()[0])
    g['outputFinalSurfDisp'] = int(nxt()[0])
    next(it)  # blank
    vals = nxt()
    g['npx'], g['npy'], g['npz'] = int(vals[0]), int(vals[1]), int(vals[2])
    next(it)  # blank
    g['totalSimuTime'] = float(nxt()[0])
    g['dt'] = float(nxt()[0])
    next(it)  # blank
    vals = nxt()
    g['nmat'], g['n2mat'] = int(vals[0]), int(vals[1])
    vals = nxt()
    g['roumax'], g['rhow'], g['gamar'] = (float(v) for v in vals)
    vals = nxt()
    g['rdampk'], g['vmaxPML'] = (float(v) for v in vals)
    next(it)  # blank
    vals = nxt()
    g['xsource'], g['ysource'], g['zsource'] = (float(v) for v in vals)
    vals = nxt()
    g['nucR'], g['nucRuptVel'], g['nucdtau0'], g['nucT'] = (float(v) for v in vals)
    vals = nxt()
    g['str1ToFaultAngle'], g['devStrToStrVertRatio'] = (float(v) for v in vals)
    vals = nxt()
    g['bulk'], g['coheplas'] = (float(v) for v in vals)
    vals = nxt()
    g['fstrike'], g['fdip'] = (float(v) for v in vals)
    g['slipRateThres'] = float(nxt()[0])

    # Viscoplastic / plastic-output block. Port of the matching reads in
    # readglobal: tv (the Duvaut-Lions relaxation time, which the Fortran
    # derived as 2*dz/3464 until it got this input slot), the deviatoric
    # pre-stress depth taper, and the plastic-strain output window.
    # A file without the block was written by an older case.setup and is
    # REFUSED, not defaulted -- same verdict as the Fortran's
    # stopStaleGlobal/ERR_INPUT_FILE_STALE (PROJECT_RULES.md rule 2).
    def stale(what):
        return ValueError(
            'read_bglobal: %s ends before %s. This file was written by an '
            'older case.setup than this port; re-run case.setup in that case '
            'directory to regenerate it.' % (path, what))

    def nxtOrStale(what):
        try:
            return next(it).split()
        except StopIteration:
            raise stale(what)

    nxtOrStale('the viscoplastic/plastic-output block separator')
    try:
        g['tv'] = float(nxtOrStale('the viscoplastic relaxation time Tv')[0])
        vals = nxtOrStale('the deviatoric pre-stress taper depths')
        g['devStrTaperDepthStart'], g['devStrTaperDepthEnd'] = float(vals[0]), float(vals[1])
        vals = nxtOrStale('the plastic-strain output window half-widths')
        g['plasticOutputHalfWidth'] = tuple(float(v) for v in vals[:3])
    except IndexError:
        raise stale('a complete viscoplastic/plastic-output block')
    if len(g['plasticOutputHalfWidth']) != 3:
        raise stale('three plastic-strain output window half-widths')
    # Station n-stress sign convention (board row 22a), the file's last line:
    # +1 positive means extension, -1 positive means compression -- the case's
    # SCEC spec convention (scripts/lib.py resolveNormalStressSign). Same
    # stale/invalid verdicts as readInputFiles.f90's readglobal.
    try:
        g['nStressOutSign'] = int(nxtOrStale('the station normal-stress sign convention')[0])
    except (IndexError, ValueError):
        raise stale('a valid station normal-stress sign convention')
    if g['nStressOutSign'] not in (1, -1):
        # Same numbered exit as readInputFiles.f90's abortRun(ERR_CFG_NSTRESS_SIGN_INVALID)
        # and the same message text (rule 23); eqdyna3d.main's _abort turns it into exit 15.
        raise InputConsistencyError(
            ERR_CFG_NSTRESS_SIGN_INVALID,
            'bGlobal.txt station normal-stress sign must be +1 (extension) or -1 '
            '(compression); re-run case.setup. (read %d from %s)' % (g['nStressOutSign'], path))

    g['str1ToFaultAngle'] = g['str1ToFaultAngle'] * np.pi / 180.0
    # nstep = idnint(totalSimuTime/dt) -- Fortran round-half-away-from-zero.
    ratio = g['totalSimuTime'] / g['dt']
    g['nstep'] = int(np.sign(ratio) * np.floor(np.abs(ratio) + 0.5))
    g['rdampk'] = g['rdampk'] * g['dt']
    return g


def read_bmodelgeometry(path):
    """Port of readmodelgeometry. Returns a dict with xmin/xmax/.../dz."""
    lines = list(_read_records(path))
    it = iter(lines)
    vals = next(it).split()
    xmin, xmax = float(vals[0]), float(vals[1])
    vals = next(it).split()
    ymin, ymax = float(vals[0]), float(vals[1])
    vals = next(it).split()
    zmin, zmax = float(vals[0]), float(vals[1])
    next(it)  # blank
    vals = next(it).split()
    dis4uniF, dis4uniB = int(vals[0]), int(vals[1])
    rat = float(next(it).split()[0])
    vals = next(it).split()
    dx, dy, dz = (float(v) for v in vals)
    return dict(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                dis4uniF=dis4uniF, dis4uniB=dis4uniB, rat=rat, dx=dx, dy=dy, dz=dz)


def faultTag(ift, ntotft):
    """The per-fault filename/variable-name tag -- port of eqquasi's/
    library_output.f90's faultTag(), verbatim: '' for fault 1, 'ft<N>_' for
    fault N>=2, but only when ntotft > 1 (so a single-fault case's file/
    variable names are untouched, bit-for-bit). `ift` is 1-indexed, matching
    the Fortran convention. Mirrors scripts/lib.py's faultTag exactly (that
    one writes the case input this one reads)."""
    if ntotft > 1 and ift > 1:
        return 'ft%d_' % ift
    return ''


def read_bfaultgeometry(path, ntotft):
    """Port of readfaultgeometry. Returns a list of ntotft dicts, each with
    fxmin/fxmax/fymin/fymax/fzmin/fzmax. Row 17 (multi-fault): every fault's
    box is read structurally (always was); meshgen.py's builders now
    consume more than the first entry for ntotft>1 -- see
    checkInputConsistency's multi-fault guards (C_degen==0, planar, dy-
    aligned, inside the uniform-y belt, distinct per-fault y, same x/z
    extent as fault 1) for the scope this is actually exercised under."""
    lines = list(_read_records(path))
    it = iter(lines)
    faults = []
    for _ in range(ntotft):
        next(it)  # "For fault No. N" header line
        vals = next(it).split()
        fxmin, fxmax = float(vals[0]), float(vals[1])
        vals = next(it).split()
        fymin, fymax = float(vals[0]), float(vals[1])
        vals = next(it).split()
        fzmin, fzmax = float(vals[0]), float(vals[1])
        faults.append(dict(fxmin=fxmin, fxmax=fxmax, fymin=fymin, fymax=fymax,
                            fzmin=fzmin, fzmax=fzmax))
    return faults


def read_bmaterial(path, nmat, n2mat):
    """Port of readmaterial's file-reading half (the ccosphi/sinphi/nstep/
    rdampk/tv derived-scalar computation is done in read_bglobal/callers
    since those don't need bMaterial.txt itself). Returns (nmat, n2mat)
    float array."""
    lines = list(_read_records(path))
    material = np.zeros((nmat, n2mat))
    for i in range(nmat):
        vals = lines[i].split()
        material[i] = [float(v) for v in vals[:n2mat]]
    return material


def read_bstations(path, ntotft=1):
    """Port of readstations1+readstations2: totalNumOfOffSt; ntotft
    on-fault station COUNTS (one per fault, line 2 -- readInputFiles.f90:152
    reads `(nonfs(i), i = 1, ntotft)`); the on-fault [x,z] rows (km),
    GROUPED BY FAULT in that same order (`do i = 1, ntotft; do j = 1,
    nonfs(i)`, scripts/case.setup's resolveOnFaultStationsPerFault writes
    them this way); then totalNumOfOffSt off-fault [x,y,z] rows (km) --
    converted to m exactly like the Fortran (`xonfs=xonfs*1000.0d0`,
    `x4nds=x4nds*1000.0d0`).

    Returns (xonfs_per_fault, x4nds): xonfs_per_fault is a list of ntotft
    (2, nonfs_i) arrays (row 17; ntotft==1 callers get a 1-element list,
    identical values to the old single-array return); x4nds is (3,n_off)
    unchanged (off-fault stations are not assigned per fault)."""
    lines = list(_read_records(path))
    it = iter(lines)
    n_off = int(next(it).split()[0])
    nonfs = [int(v) for v in next(it).split()]
    if len(nonfs) != ntotft:
        raise ValueError('read_bstations: %s line 2 has %d on-fault station '
                          'count(s), expected ntotft=%d' % (path, len(nonfs), ntotft))
    next(it)  # blank
    xonfs_per_fault = []
    for n_onf in nonfs:
        xonfs = np.zeros((2, n_onf))
        for j in range(n_onf):
            vals = [float(v) for v in next(it).split()]
            xonfs[:, j] = vals[:2]
        xonfs_per_fault.append(xonfs * 1000.0)
    next(it)  # blank
    x4nds = np.zeros((3, n_off))
    for i in range(n_off):
        vals = [float(v) for v in next(it).split()]
        x4nds[:, i] = vals[:3]
    return xonfs_per_fault, x4nds * 1000.0


def read_fault_rough_geometry(path):
    """Port of readInputFiles.f90's read_fault_rough_geometry: reads
    bFault_Rough_Geometry.txt (nx, nz header; dx, fx_min, fz_min header;
    then nx*nz rows of [y, dy/dx, dy/dz], z-fastest -- i.e. row index
    nz*(ix-1)+iz, 1-indexed, matching func_lib.f90's insertFaultInterface
    lookup `rough_geo(1, nnz*(ixx-1)+izz)`).

    This file is written by the case's own `generateFaultInterface` Python
    script (scripts/generateFaultInterface) at case-setup time -- NOT by
    any Fortran binary -- so reading it here keeps the standalone solver
    path Fortran-free even for insertFaultType>0 cases (tpv10, drv.a6).

    Returns a dict: nnx, nnz (int), dx (the file's own grid spacing,
    used only to derive fx_max -- NOT the model's dx/dz, which
    insertFaultInterface's index math uses instead, per the Fortran
    source), fx_min, fz_min, fx_max, rough_geo (3, nnx*nnz) float array,
    1-indexed-by-column convention preserved via rough_geo[:, col-1].
    """
    with open(path) as f:
        lines = f.readlines()
    # Fortran's list-directed `read(unit,*) nnxTmp, nnzTmp` reads only the
    # first 2 tokens off the line, silently ignoring any extra trailing
    # value(s) -- this file's header line 1 sometimes carries a 3rd,
    # unused column; sliced here to match, not an assumption.
    nnx, nnz = (int(round(float(v))) for v in lines[0].split()[:2])
    dx_file, fx_min, fz_min = (float(v) for v in lines[1].split()[:3])
    fx_max = (nnx - 1) * dx_file + fx_min
    rough_geo = np.zeros((3, nnx * nnz))
    for i in range(nnx * nnz):
        vals = [float(v) for v in lines[2 + i].split()]
        rough_geo[:, i] = vals[:3]
    return dict(nnx=nnx, nnz=nnz, dx=dx_file, fx_min=fx_min, fz_min=fz_min,
                fx_max=fx_max, rough_geo=rough_geo)


def build_params(case_dir):
    """Convenience wrapper: reads bGlobal.txt/bModelGeometry.txt/
    bFaultGeometry.txt from `case_dir` and returns (params, globals_dict)
    where params is the dict meshgen.py's build_grid_lines/build_elements/
    build_equation_numbers/build_fault_geometry/build_station_matching
    expect (dx,dy,dz,fxmin,...,rat,nPML,tol,fstrike,C_degen), and
    globals_dict is read_bglobal's full return (for nmat/n2mat/friclaw/etc,
    needed by callers but not by the geometry builders themselves).

    Row 17 (multi-fault): ntotft>1 is now accepted for the narrow release
    scope -- two (or more) vertical, planar, parallel faults at distinct y,
    each with its OWN x/z extent (Row 17 rebased: independent per-fault x/z
    extents are supported, each must simply land on a mesh node line --
    checkInputConsistency.check_multifault, called by this function's
    caller -- eqdyna3d.py's build_solver_state -- right after this returns,
    mirroring checkInputConsistency.f90's call point, no longer requires
    every fault to share fault 1's x/z extent; meshgen.py's
    one_dim_coor_array raises the commensurability guard instead).
    `params['faults']` carries every fault's own (fxmin,...,fzmax) dict --
    every builder that needs a fault's own box reads `params['faults']` (via
    `_fault_boxes`), not the scalar fxmin/fxmax/... keys below, which remain
    fault 1's box ONLY for callers that still want "the" single box (the
    ntotft==1 convention, and any genuinely fault-1-specific need);
    `params['fault_y']` is the list of each fault's y-plane (fymin, since
    planarity requires fymin==fymax), the one per-fault quantity meshgen.py's
    builders actually need. insertFaultType>0 (rough/dipping fault) combined
    with ntotft>1 is NOT supported -- out of this release's narrow scope (the
    rough-fault y-morph is a single-fault mechanism) -- and raises rather
    than silently reading only fault 1's rough geometry for every fault.

    nPML, tol, and R are NOT present in any bFile -- nPML=6, tol=1.0e-5,
    and R=0.01d0 (theoretical PML reflection coefficient, used by
    src/func_lib.f90's pmlRegionDistance / port.py's region_damp) are
    Fortran globalvar.f90 PARAMETER constants (compile-time, not case
    input), reproduced verbatim here (see globalvar.f90's `nPML = 6`,
    `tol = 1.0d-5`, `R = 0.01d0`).

    If g['insertFaultType'] > 0 (tpv10 dipping / drv.a6 rough fault),
    also reads bFault_Rough_Geometry.txt (via read_fault_rough_geometry)
    and stores it under params['rough'] -- meshgen.py's
    build_node_coordinates/build_fault_geometry use this to apply
    func_lib.f90's insertFaultInterface y-morph. params['rough'] is None
    for insertFaultType==0 (the planar-fault milestones' unchanged path).
    """
    import os
    g = read_bglobal(os.path.join(case_dir, 'bGlobal.txt'))
    if g['ntotft'] > 1 and g['insertFaultType'] > 0:
        raise NotImplementedError(
            'build_params: ntotft=%d with insertFaultType=%d is not ported -- '
            'the rough/dipping-fault y-morph (insertFaultInterface) is a '
            'single-fault mechanism, out of this release\'s narrow multi-fault '
            'scope (two-or-more vertical PLANAR parallel faults).'
            % (g['ntotft'], g['insertFaultType']))
    mg = read_bmodelgeometry(os.path.join(case_dir, 'bModelGeometry.txt'))
    faults = read_bfaultgeometry(os.path.join(case_dir, 'bFaultGeometry.txt'), g['ntotft'])
    fg = faults[0]
    rough = None
    if g['insertFaultType'] > 0:
        rough = read_fault_rough_geometry(os.path.join(case_dir, 'bFault_Rough_Geometry.txt'))
    params = dict(
        dx=mg['dx'], dy=mg['dy'], dz=mg['dz'],
        fxmin=fg['fxmin'], fxmax=fg['fxmax'], fzmin=fg['fzmin'], fzmax=fg['fzmax'],
        fymin=fg['fymin'], fymax=fg['fymax'],
        dis4uniF=mg['dis4uniF'], dis4uniB=mg['dis4uniB'],
        xmin=mg['xmin'], xmax=mg['xmax'], ymin=mg['ymin'], ymax=mg['ymax'],
        zmin=mg['zmin'], zmax=mg['zmax'],
        rat=mg['rat'], nPML=6, tol=1.0e-5, R=0.01,
        fstrike=g['fstrike'], C_degen=g['C_degen'],
        insertFaultType=g['insertFaultType'], rough=rough,
        # Row 17 (multi-fault): every fault's own box, and the one per-fault
        # scalar (the y-plane) the C_degen==0 builders actually branch on.
        ntotft=g['ntotft'], faults=faults,
        fault_y=[f['fymin'] for f in faults],
    )
    return params, g


def read_on_fault_vars(nc_path, fxmin, fzmin, dx, dz, meshCoor, nsmp, ntotft=1, fault_of=None):
    """Port of netcdf_io.f90's netcdf_read_on_fault_eqdyna.

    Row 17 (multi-fault): `fault_of` is the (nftnd,) 1-indexed fault id for
    each row of `nsmp` (build_node_coordinates' return); each row reads its
    OWN fault's variable set, named with readInputFiles.faultTag(ift,
    ntotft) -- '' for fault 1 (bit-identical variable names at ntotft==1,
    the old default), 'ft<N>_' for fault N>=2 -- matching
    scripts/case.setup's netcdf_write_on_fault_vars and
    src/fortran/netcdf_io.f90's own per-fault variable-set read exactly.
    `fault_of=None` is the ntotft==1 shorthand (every row reads the
    untagged set), unchanged from before this fix.

    Row 17 REBASED (restore per-fault mesh extent): `fxmin`/`fzmin` are now
    PER-FAULT sequences (fxmin[ift-1], fzmin[ift-1]), not one shared scalar
    -- the Fortran side (netcdf_io.f90's `fxmin(ift)`/`fzmin(ift)`, globalvar
    arrays, never a scalar) always indexed the ii/jj offset per fault; this
    port's scalar `fxmin`/`fzmin` was the one remaining fault-1-only
    shortcut, latent as long as checkInputConsistency required every fault
    to share fault 1's x/z extent exactly. A single-fault or scalar caller
    may still pass a bare float/int: wrapped into a length-1 sequence below,
    reproducing the old ntotft==1 behaviour bit-for-bit (every row then
    reads index 0 regardless of `fault_of`, same as the old shared scalar).

    Reads on_fault_vars_input.nc directly via the netCDF4 library (the
    SAME library the Fortran-side writer, scripts/case.setup, uses to
    create the file -- not re-deriving the writer's own logic). Each
    variable is written by case.setup with dims ('dip','strike'), i.e.
    numpy shape (nfz,nfx) -- Fortran's `nf90_get_var` transposes this into
    its own (fnx,fnz,nvar) `on_fault_vars` array per netCDF's documented
    Fortran/C dimension-order reversal (see netcdf_io.f90's header
    comment); reading directly with python's netCDF4 library returns the
    array in its ORIGINAL (dip,strike)=(iz,ix) shape, so this port indexes
    `var[jj-1, ix_idx]` where Fortran indexes `on_fault_vars(ii,jj,ivar)` --
    same value, no transpose needed, just care with which index is which.

    ii/jj (Fortran 1-indexed grid-cell index along strike/dip) come from
    `nint((xcord-fxmin)/dx)+1` / `nint((zcord-fzmin)/dz)+1` -- Fortran's
    `nint` is round-half-away-from-zero, reproduced exactly (not Python's
    round-half-to-even `round()`, and not `np.round` which is also
    round-half-to-even).

    meshCoor: (N+1,3) from build_node_coordinates (1-indexed rows).
    nsmp: (nftnd,2) int64 [slave_id, master_id], from build_node_coordinates.

    Returns fric: (nftnd+1, 101) float array (row 0 and column 0 unused,
    1-indexed fault-node rows AND 1-indexed FRIC_SLOT_* columns, matching
    globalvar.f90's `fric(100,nftmx,ntotft)` slot numbering) with every
    slot netcdf_read_on_fault_eqdyna's ntotft==1 branch writes populated;
    every other slot left at exactly 0.0, matching Fortran's `fric = 0.0d0`
    zero-init (confirmed by reading src/eqdyna3d.f90 directly) that runs
    before this subroutine is called.
    """
    import netCDF4

    nftnd = nsmp.shape[0]
    fric = np.zeros((nftnd + 1, 101))
    if fault_of is None:
        fault_of = np.ones(nftnd, dtype=np.int64)
    # Row 17 rebased: accept a bare scalar (old single-box convention) or a
    # per-fault sequence indexed 0-based (fxmin_by_fault[ift-1]).
    fxmin_by_fault = [fxmin] if np.isscalar(fxmin) else list(fxmin)
    fzmin_by_fault = [fzmin] if np.isscalar(fzmin) else list(fzmin)

    varnames = ['sw_fs', 'sw_fd', 'sw_D0', 'rsf_a', 'rsf_b', 'rsf_Dc', 'rsf_v0',
                'rsf_r0', 'rsf_fw', 'rsf_vw', 'tp_a_hy', 'tp_a_th', 'tp_rouc',
                'tp_lambda', 'tp_h', 'tp_Tini', 'tp_pini', 'init_slip_rate',
                'init_strike_shear', 'init_normal_stress', 'init_state', 'tw_t0',
                'cohesion', 'init_dip_shear']

    ds = netCDF4.Dataset(nc_path, 'r')
    try:
        # Row 17: one dict of {name: array} PER FAULT TAG actually present
        # among fault_of -- read once per tag, not once per row.
        on_fault_vars_by_tag = {}
        for ift in sorted(set(int(v) for v in fault_of)):
            tag = faultTag(ift, ntotft)
            on_fault_vars_by_tag[ift] = {
                name: np.asarray(ds.variables[tag + name][:, :]) for name in varnames}
    finally:
        ds.close()

    def fortran_nint(x):
        return int(np.sign(x) * np.floor(np.abs(x) + 0.5)) if x != 0 else 0

    # FRIC_SLOT_* constants, from globalvar.f90 (verbatim, not renumbered).
    SW_FS, SW_FD, SW_D0, COHESION = 1, 2, 3, 4
    TW_T0 = 5
    INIT_NORM, INIT_STRIKE_SHEAR = 7, 8
    RSF_A, RSF_B, RSF_DC, RSF_V0, RSF_R0, RSF_FW, RSF_VW = 9, 10, 11, 12, 13, 14, 15
    TP_A_HY, TP_A_TH, TP_ROUC, TP_LAMBDA = 16, 17, 18, 19
    THETA_PC = 23
    TP_H, TP_TINI, TP_PINI = 40, 41, 42
    CREEP_VMIN, PEAK_SLIPRATE, INIT_DIP_SHEAR = 46, 47, 49
    VINI_N, VINI_X, VINI_Z = 25, 26, 27
    STATE = 20

    for i in range(1, nftnd + 1):
        slave = int(nsmp[i - 1, 0])
        xcord = meshCoor[slave, 0]
        zcord = meshCoor[slave, 2]
        ift_row = int(fault_of[i - 1])
        # Row 17 rebased: THIS row's own fault's box, not a shared one --
        # fxmin_by_fault has length 1 for the old scalar convention, so
        # ift_row-1 may be out of range there only if fault_of itself claims
        # more than one fault while the caller supplied a single box, which
        # would be a caller bug (mismatched arguments), not silently papered
        # over by clamping to index 0.
        ii = fortran_nint((xcord - fxmin_by_fault[ift_row - 1]) / dx) + 1
        jj = fortran_nint((zcord - fzmin_by_fault[ift_row - 1]) / dz) + 1
        on_fault_vars = on_fault_vars_by_tag[ift_row]

        def v(name):
            return on_fault_vars[name][jj - 1, ii - 1]

        fric[i, SW_FS] = v('sw_fs')
        fric[i, SW_FD] = v('sw_fd')
        fric[i, SW_D0] = v('sw_D0')
        fric[i, RSF_A] = v('rsf_a')
        fric[i, RSF_B] = v('rsf_b')
        fric[i, RSF_DC] = v('rsf_Dc')
        fric[i, RSF_V0] = v('rsf_v0')
        fric[i, RSF_R0] = v('rsf_r0')
        fric[i, RSF_FW] = v('rsf_fw')
        fric[i, RSF_VW] = v('rsf_vw')
        fric[i, TP_A_HY] = v('tp_a_hy')
        fric[i, TP_A_TH] = v('tp_a_th')
        fric[i, TP_ROUC] = v('tp_rouc')
        fric[i, TP_LAMBDA] = v('tp_lambda')
        fric[i, TP_H] = v('tp_h')
        fric[i, TP_TINI] = v('tp_Tini')
        fric[i, TP_PINI] = v('tp_pini')
        fric[i, CREEP_VMIN] = v('init_slip_rate')
        fric[i, INIT_STRIKE_SHEAR] = v('init_strike_shear')
        fric[i, INIT_NORM] = v('init_normal_stress')
        fric[i, STATE] = v('init_state')
        fric[i, PEAK_SLIPRATE] = fric[i, CREEP_VMIN]
        fric[i, VINI_N] = 0.0
        fric[i, VINI_X] = fric[i, CREEP_VMIN]
        fric[i, VINI_Z] = 0.0
        fric[i, TW_T0] = v('tw_t0')
        fric[i, COHESION] = v('cohesion')
        fric[i, INIT_DIP_SHEAR] = v('init_dip_shear')
        fric[i, THETA_PC] = abs(fric[i, INIT_NORM])

    return fric
