#! /usr/bin/env python3
"""
Regression guard: the n2mat==6 3D structured material grid (SCEC TPV34,
CVM-H) must REFUSE rather than silently clamp when the grid does not cover
the case's declared mesh box (bModelGeometry.txt xmin/xmax/ymin/ymax/zmin/
zmax). This was the open Medium from the TPV34 PR #93 pre-merge audit: the
invariant ("tools.materialGrid's box MUST equal the case box",
case_input/test.tpv34/user_defined_params.py) was only a comment, not a
checked assertion, on either side of the port.

The fix lives in BOTH implementations, same condition, same error code
(ERR_CFG_MATERIAL_GRID_INVALID == 17, errorCodes.f90):

  Fortran: src/fortran/readInputFiles.f90 buildMaterialGrid3D, right after
    matGridOrigin/matGridSpacing/matGridCount are derived (same place the
    existing nx*ny*nz completeness check lives).
  Python:  src/python/eqdyna/checkInputConsistency.py build_material_grid3d's
    new `domain_box` keyword (plumbed from params['xmin']/.../['zmax'] at
    both its call sites: eqdyna3d.py's early validation pass and
    meshgen.py's build_elements/_build_elements_scalar).

This guard builds ONE minimal case with a real n2mat==6 grid, twice: once
with a grid that exactly covers the declared mesh box (must run to
completion on both binaries, unaffected by this change) and once with the
SAME grid missing its outermost x+ plane (reach_x_max shrinks below xmax --
must be REFUSED by both binaries with exit code 17, not silently accepted
with the x+ elements clamped to the x=250 plane's velocity).

Cheap (rule 9): a 4x4x2-cell domain (32-row grid), dx=500m, no mesh/solve
work needed for the refused run (the check fires in readmaterial, before
meshgen). The covering run uses a 4-step term so it finishes fast on both
backends.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402  (item 95: rank-owned output past MPI_Abort)

MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE_NAME = 'test.tpv8'
ERR_CFG_MATERIAL_GRID_INVALID = 17
COVERAGE_MSG_FRAGMENT = 'does not cover the mesh box'

# A uniform domain whose element centres are an EXACT multiple of dx, so the
# "full" grid covers it exactly (reach == box, as in test.tpv34's own case).
XMIN, XMAX = -1.0e3, 1.0e3
YMIN, YMAX = -1.0e3, 1.0e3
ZMIN, ZMAX = -1.0e3, 0.0
DX = 500.0


def _full_grid_rows():
    """Every element-centre (x, y, z) of the domain above at spacing DX,
    homogeneous vp/vs/rho -- a valid, covering n2mat==6 table."""
    xs = np.arange(XMIN + DX / 2, XMAX, DX)
    ys = np.arange(YMIN + DX / 2, YMAX, DX)
    zs = np.arange(ZMIN + DX / 2, ZMAX, DX)
    rows = []
    for x in xs:
        for y in ys:
            for z in zs:
                rows.append([x, y, z, 5716.0, 3300.0, 2700.0])
    return np.array(rows)


USER_PARAMS = '''#! /usr/bin/env python3
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = %(xmin)r, %(xmax)r
par.ymin, par.ymax = %(ymin)r, %(ymax)r
par.zmin, par.zmax = %(zmin)r, %(zmax)r

par.fxmin, par.fxmax = -500.0, 500.0
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -500.0, 0.0

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -250.0

par.dx = %(dx)r
par.dy = par.dx
par.dz = par.dx

par.mat = np.array(%(mat)r)
par.nmat, par.n2mat = par.mat.shape
par.roumax = par.mat[:, 5].max()
par.vmaxPML = par.mat[:, 3].max()
par.vp, par.vs, par.rou = 5716., 3300., 2700.   # clamp floor; not used for nmat>1

par.dt = 0.5*par.dx/par.vmaxPML
par.term = 4.*par.dt
par.friclaw = 1
par.tpv = 8
par.nucR = 200.

par.nx, par.ny, par.nz = 1, 1, 1

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav = 9.8
par.fric_cohesion = 1.e6
for ix, xcoor in enumerate(par.fx):
    for iz, zcoor in enumerate(par.fz):
        par.on_fault_vars[iz, ix, 1] = par.fric_sw_fs
        if abs(abs(xcoor) - par.fxmax) < 0.01 or abs(zcoor - par.fzmin) < 0.01:
            par.on_fault_vars[iz, ix, 1] = 1000.
        par.on_fault_vars[iz, ix, 2] = par.fric_sw_fd
        par.on_fault_vars[iz, ix, 3] = par.fric_sw_D0
        par.on_fault_vars[iz, ix, 4] = par.fric_cohesion
        par.on_fault_vars[iz, ix, 7] = 7378.*zcoor
        par.on_fault_vars[iz, ix, 8] = abs(0.55*par.on_fault_vars[iz, ix, 7])
        if abs(xcoor-par.xsource) <= 200. and abs(zcoor-par.zsource) <= 200.:
            par.on_fault_vars[iz, ix, 8] = 1e6+abs(1.005*0.76*par.on_fault_vars[iz, ix, 7])

par.st_coor_on_fault = [[0.0, 0.0]]
par.st_coor_off_fault = [[0, 1, 0]]
par.n_on_fault = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
par.output_plastic = 0
'''


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(tmp, name, mat_rows):
    case_dir = os.path.join(tmp, 'case_%s' % name)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('create.newcase failed (%d):\n%s'
                           % (r.returncode, r.stdout[-2000:] + r.stderr[-2000:]))
    params = USER_PARAMS % dict(xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX,
                                 zmin=ZMIN, zmax=ZMAX, dx=DX,
                                 mat=mat_rows.tolist())
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'w') as f:
        f.write(params)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('case.setup failed for %s:\n%s' % (name, r.stdout + r.stderr))
    return case_dir


def _run_fortran(case_dir):
    binary = None
    for cand in (os.path.join(ROOT, 'bin', 'eqdyna'),
                 os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        if os.path.exists(cand):
            binary = cand
            break
    if binary is None:
        return None, None
    return mpirun_capture.run_rank_files(MPIRUN, binary, case_dir)


def _run_python(case_dir, nsteps=4):
    r = subprocess.run([sys.executable, '-u', '-m', 'eqdyna', case_dir, str(nsteps),
                        '--backend', 'numpy'],
                       cwd=ROOT, env=_env(), capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr


def check_noncovering_grid_refused_end_to_end():
    full = _full_grid_rows()
    # Drop the x+ extreme plane (the largest x value): reach_x_max shrinks
    # from XMAX to (XMAX - DX), no longer covering the declared mesh box.
    x_vals = np.unique(full[:, 0])
    bad = full[full[:, 0] != x_vals.max()]
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'noncovering', bad)

        rc_f, out_f = _run_fortran(case_dir)
        if rc_f is None:
            raise AssertionError('bin/eqdyna (or src/fortran/eqdyna) not found; '
                                 'build it first with ./install-eqdyna.sh')
        if rc_f == 0:
            raise AssertionError('Fortran ACCEPTED a material grid narrower than the '
                                 'mesh box (exit 0) -- buildMaterialGrid3D should have '
                                 'refused this')
        if rc_f != ERR_CFG_MATERIAL_GRID_INVALID:
            raise AssertionError('Fortran exited %d, expected %d '
                                 '(ERR_CFG_MATERIAL_GRID_INVALID):\n%s'
                                 % (rc_f, ERR_CFG_MATERIAL_GRID_INVALID, out_f[-1500:]))
        if COVERAGE_MSG_FRAGMENT not in out_f:
            raise AssertionError('Fortran refused (exit %d) but did not print the '
                                 'expected coverage message:\n%s' % (rc_f, out_f[-1500:]))

        rc_p, out_p = _run_python(case_dir)
        if rc_p == 0:
            raise AssertionError('python port ACCEPTED a material grid narrower than '
                                 'the mesh box (exit 0) -- build_material_grid3d should '
                                 'have refused this')
        if rc_p != ERR_CFG_MATERIAL_GRID_INVALID:
            raise AssertionError('python port exited %d, expected %d '
                                 '(ERR_CFG_MATERIAL_GRID_INVALID):\n%s'
                                 % (rc_p, ERR_CFG_MATERIAL_GRID_INVALID, out_p[-1500:]))
        if COVERAGE_MSG_FRAGMENT not in out_p:
            raise AssertionError('python port refused (exit %d) but did not print the '
                                 'expected coverage message:\n%s' % (rc_p, out_p[-1500:]))
    print('  PASS  non-covering n2mat=6 grid refused by BOTH binaries, exit %d, '
          'message matches' % ERR_CFG_MATERIAL_GRID_INVALID)


def check_covering_grid_still_runs_both():
    full = _full_grid_rows()
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'covering', full)

        rc_f, out_f = _run_fortran(case_dir)
        if rc_f is None:
            raise AssertionError('bin/eqdyna (or src/fortran/eqdyna) not found; '
                                 'build it first with ./install-eqdyna.sh')
        if rc_f != 0:
            raise AssertionError('Fortran REJECTED a covering n2mat=6 grid (exit %d):\n%s'
                                 % (rc_f, out_f[-2000:]))
        frt_f = [f for f in os.listdir(case_dir) if f.startswith('frt.txt')]
        if not frt_f:
            raise AssertionError('Fortran run wrote no frt.txt* for a covering grid')

        rc_p, out_p = _run_python(case_dir)
        if rc_p != 0:
            raise AssertionError('python port REJECTED a covering n2mat=6 grid (exit %d):\n%s'
                                 % (rc_p, out_p[-2000:]))
        if not os.path.isfile(os.path.join(case_dir, 'frt.txt0')):
            raise AssertionError('python port wrote no frt.txt0 for a covering grid:\n%s'
                                 % out_p[-2000:])
    print('  PASS  covering n2mat=6 grid (exact box match, test.tpv34\'s own shape) '
          'runs to completion on both binaries, unaffected by this change')


def main():
    print('Regression guard: n2mat==6 material-grid mesh-box coverage check '
          '(Fortran buildMaterialGrid3D <-> Python build_material_grid3d)')
    if not os.path.exists(os.path.join(ROOT, 'bin', 'eqdyna')) and \
       not os.path.exists(os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        print('FAIL test_material_grid3d_coverage: bin/eqdyna (or src/fortran/eqdyna) '
              'not found -- build it with ./install-eqdyna.sh first (this guard does '
              'not skip silently)')
        return 1
    fails = []
    for c in (check_noncovering_grid_refused_end_to_end,
              check_covering_grid_still_runs_both):
        try:
            c()
        except (AssertionError, RuntimeError) as e:
            fails.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if fails:
        print('\nFAIL test_material_grid3d_coverage (%d check(s))' % len(fails))
        return 1
    print('\nSUCCESS test_material_grid3d_coverage')
    return 0


if __name__ == '__main__':
    sys.exit(main())
