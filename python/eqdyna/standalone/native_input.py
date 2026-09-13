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
testsys/parity/test_standalone_meshgen.py: the resulting PARAMS dict,
material array, xonfs/x4nds station arrays, and initial fric array are fed
straight into meshgen.py's M1-M5 builders and `read_on_fault_vars` below,
and the M1-M5 parity checks (already gated against pydump_meshCoor.txt/
pydump_conn.txt/pydump_nodeinfo.txt/pydump_fault.txt/pydump_stations.txt)
are re-run with these NATIVELY-READ inputs instead of the hand-transcribed
PARAMS dict used by earlier milestones -- same tolerances, same fixtures,
now proving the reader itself, not just the builders.
"""
import numpy as np


def _read_records(path):
    """Yield each non-blank, non-comment-only line's whitespace-split tokens,
    in file order -- mirrors Fortran list-directed READ's own skip-blank-
    lines behavior (`read(unit,*)` with an empty format skips nothing itself,
    but every blank-line placeholder in these files is read by its own bare
    `read(1001,*)` in the Fortran, so line-for-line correspondence is exact,
    not a `skip blanks` heuristic)."""
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


def read_bfaultgeometry(path, ntotft):
    """Port of readfaultgeometry. Returns a list of ntotft dicts, each with
    fxmin/fxmax/fymin/fymax/fzmin/fzmax (this milestone: ntotft==1 only,
    matching meshgen.py's scope; ntotft>1 is read structurally but
    downstream builders raise if handed more than one)."""
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


def read_bstations(path):
    """Port of readstations1+readstations2 (ntotft==1 only, matching this
    milestone's scope): totalNumOfOffSt, nonfs(1), nonfs(1) on-fault [x,z]
    rows (km), totalNumOfOffSt off-fault [x,y,z] rows (km) -- converted to
    m exactly like the Fortran (`xonfs=xonfs*1000.0d0`, `x4nds=x4nds*1000.0d0`).
    Returns (xonfs, x4nds), shapes (2,nonfs) and (3,n_off)."""
    lines = list(_read_records(path))
    it = iter(lines)
    n_off = int(next(it).split()[0])
    n_onf = int(next(it).split()[0])
    next(it)  # blank
    xonfs = np.zeros((2, n_onf))
    for j in range(n_onf):
        vals = [float(v) for v in next(it).split()]
        xonfs[:, j] = vals[:2]
    next(it)  # blank
    x4nds = np.zeros((3, n_off))
    for i in range(n_off):
        vals = [float(v) for v in next(it).split()]
        x4nds[:, i] = vals[:3]
    return xonfs * 1000.0, x4nds * 1000.0


def build_params(case_dir):
    """Convenience wrapper: reads bGlobal.txt/bModelGeometry.txt/
    bFaultGeometry.txt from `case_dir` and returns (params, globals_dict)
    where params is the dict meshgen.py's build_grid_lines/build_elements/
    build_equation_numbers/build_fault_geometry/build_station_matching
    expect (dx,dy,dz,fxmin,...,rat,nPML,tol,fstrike,C_degen), and
    globals_dict is read_bglobal's full return (for nmat/n2mat/friclaw/etc,
    needed by callers but not by the geometry builders themselves).

    ntotft must be 1 (this milestone's scope) -- raises otherwise, not a
    silent first-fault-only truncation.

    nPML, tol, and R are NOT present in any bFile -- nPML=6, tol=1.0e-5,
    and R=0.01d0 (theoretical PML reflection coefficient, used by
    src/func_lib.f90's pmlRegionDistance / port.py's region_damp) are
    Fortran globalvar.f90 PARAMETER constants (compile-time, not case
    input), reproduced verbatim here (see globalvar.f90's `nPML = 6`,
    `tol = 1.0d-5`, `R = 0.01d0`).
    """
    import os
    g = read_bglobal(os.path.join(case_dir, 'bGlobal.txt'))
    if g['ntotft'] != 1:
        raise NotImplementedError('build_params: only ntotft==1 is ported '
                                   '(got %d)' % g['ntotft'])
    mg = read_bmodelgeometry(os.path.join(case_dir, 'bModelGeometry.txt'))
    faults = read_bfaultgeometry(os.path.join(case_dir, 'bFaultGeometry.txt'), g['ntotft'])
    fg = faults[0]
    params = dict(
        dx=mg['dx'], dy=mg['dy'], dz=mg['dz'],
        fxmin=fg['fxmin'], fxmax=fg['fxmax'], fzmin=fg['fzmin'], fzmax=fg['fzmax'],
        fymin=fg['fymin'], fymax=fg['fymax'],
        dis4uniF=mg['dis4uniF'], dis4uniB=mg['dis4uniB'],
        xmin=mg['xmin'], xmax=mg['xmax'], ymin=mg['ymin'], ymax=mg['ymax'],
        zmin=mg['zmin'], zmax=mg['zmax'],
        rat=mg['rat'], nPML=6, tol=1.0e-5, R=0.01,
        fstrike=g['fstrike'], C_degen=g['C_degen'],
    )
    return params, g


def read_on_fault_vars(nc_path, fxmin, fzmin, dx, dz, meshCoor, nsmp):
    """Port of netcdf_io.f90's netcdf_read_on_fault_eqdyna (ntotft==1 only).

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

    varnames = ['sw_fs', 'sw_fd', 'sw_D0', 'rsf_a', 'rsf_b', 'rsf_Dc', 'rsf_v0',
                'rsf_r0', 'rsf_fw', 'rsf_vw', 'tp_a_hy', 'tp_a_th', 'tp_rouc',
                'tp_lambda', 'tp_h', 'tp_Tini', 'tp_pini', 'init_slip_rate',
                'init_strike_shear', 'init_normal_stress', 'init_state', 'tw_t0',
                'cohesion', 'init_dip_shear']

    ds = netCDF4.Dataset(nc_path, 'r')
    try:
        on_fault_vars = {name: np.asarray(ds.variables[name][:, :]) for name in varnames}
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
        ii = fortran_nint((xcord - fxmin) / dx) + 1
        jj = fortran_nint((zcord - fzmin) / dz) + 1

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
