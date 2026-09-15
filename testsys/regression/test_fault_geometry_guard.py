#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the unvalidated rough-fault
geometry reader.

Background: until v5.6.0 src/readInputFiles.f90:read_fault_rough_geometry took
bFault_Rough_Geometry.txt entirely on faith. It read nnx/nnz from the header,
allocated rough_geo(3,nnx*nnz), and read exactly that many rows -- with no
check that the header described the mesh being built, that the file had that
many rows, or that the numbers were finite. Three ways to lose silently:

  * a header disagreeing with the case's fault grid: src/func_lib.f90's
    insertFaultInterface indexes rough_geo with the MESH dx/dz
    (nint((x - rough_fx_min)/dx) + 1) while rough_fx_max comes from the
    HEADER's dx, so a surface sampled at another spacing is stretched onto
    the fault and every node lands on the wrong point of it;
  * a short file: the read loop runs past the end of the data;
  * a NaN: it becomes a mesh coordinate;
  * a surface on the right grid whose per-cell fault-normal climb reaches a
    full cell: the inserted off-fault element layer tangles.

None of these produced an error message. insertFaultType=3 ("the case supplies
its own geometry") became a supported mode in v5.5.0, so hand-written and
externally-produced files are now a first-class path.

Fixed on two levels: scripts/lib.py:validateFaultRoughGeometry, called from
scripts/case.setup for every insertFaultType > 0 (unit-tested per failure mode
in testsys/unit/test_lib.py), and the Fortran checks this script exercises --
defense in depth, because these files get hand-edited and because a case set up
on one machine is routinely copied to an HPC where case.setup never runs again.

This test builds ONE tiny insertFaultType=3 case with a flat (all-zero) rough
surface -- a valid rough-geometry file whose morph is the identity -- runs it
once to establish that a GOOD file is accepted, then re-runs the same binary in
copies of that case whose geometry file has exactly one thing wrong with it.

Asserts:
  1. The clean case exits 0 and never prints the guard banner. The guard must
     not block a valid file.
  2. Every corrupted file exits non-zero AND prints the guard banner -- exit
     code and message, not one or the other (rule 3: "the script ran" and "the
     script passed" are different questions).
  3. Each failure names the specific defect, so the message is actionable.

Cheap (rule 9): a 5x3 fault node grid, term = 5 dt; the corrupted runs stop
before meshing. Builds ONLY src/eqdyna, never bin/eqdyna, since other processes
on this box may be running bin/eqdyna concurrently. Fails loudly if
mpif90/mpirun are unavailable rather than skipping.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'src')
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

GUARD_MSG = 'read_fault_rough_geometry: bFault_Rough_Geometry.txt is not usable'

# 5 x 5 fault nodes: fxmin -1000 .. fxmax 1000 at dx 500, fzmin -2000 .. 0.
# Five along EACH axis is the minimum the interior derivative check accepts --
# it needs a 4-point y''' stencil and now refuses a smaller grid outright
# rather than falling back to a weaker bound (rule 2). This case was 5 x 3.
NNX, NNZ = 5, 5

USER_PARAMS = '''#! /usr/bin/env python3
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -3.0e3, 3.0e3
par.ymin, par.ymax = -3.0e3, 3.1e3
par.zmin, par.zmax = -3.0e3, 0.0e3

par.fxmin, par.fxmax = -1.0e3, 1.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -2.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -500.0

par.dx = 500.
par.dy = par.dx
par.dz = par.dx
par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700

par.dt = 0.5*par.dx/par.vp
par.term = 5.*par.dt
par.friclaw = 1
par.tpv = 8
par.insertFaultType = 3     # the case supplies bFault_Rough_Geometry.txt

par.nucR = 300.

par.nx, par.ny, par.nz = 1, 1, 1

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin,par.fxmax,par.nfx)
par.fz  = np.linspace(par.fzmin,par.fzmax,par.nfz)

par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav = 9.8
par.fric_cohesion = 1.e6
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4]   = par.fric_cohesion
    par.on_fault_vars[iz,ix,7]   = 7378.*zcoor
    par.on_fault_vars[iz,ix,8]   = abs(0.55*par.on_fault_vars[iz,ix,7])
    if abs(xcoor-par.xsource)<=300. and abs(zcoor-par.zsource)<=300.:
        par.on_fault_vars[iz,ix,8] = 1e6+abs(1.005*0.76*par.on_fault_vars[iz,ix,7])

par.st_coor_on_fault = [[0.0, 0.0]]
par.st_coor_off_fault = [[0,1,0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
'''


def _header(lines, nnx=NNX, nnz=NNZ, dx=500.0, fxmin=-1000.0, fzmin=-2000.0):
    out = list(lines)
    out[0] = f'{nnx}\t{nnz}\t0\n'
    out[1] = f'{dx:.6f}\t{fxmin:.6f}\t{fzmin:.6f}\n'
    return out


# name -> (mutation of the good file's lines, what the Fortran must notice)
CORRUPTIONS = {
    'header_nnx_too_large':
        (lambda l: _header(l, nnx=NNX + 1),
         'not the fault grid of the mesh'),
    'header_nnz_too_small':
        (lambda l: _header(l, nnz=NNZ - 1),
         'not the fault grid of the mesh'),
    'header_dx_halved':
        (lambda l: _header(l, dx=250.0),
         'different cell size'),
    'header_wrong_corner':
        (lambda l: _header(l, fxmin=-900.0),
         'different fault corner'),
    'short_file':
        (lambda l: l[:-1],
         'short of, or has a malformed'),
    'long_file':
        (lambda l: l + [l[-1]],
         'more data rows than its header declares'),
    'nan_value':
        (lambda l: l[:5] + ['NaN\t0.0\t0.0\n'] + l[6:],
         'NaN or Inf'),
    # A surface on the RIGHT grid whose slope tangles the inserted elements:
    # wrong-but-plausible, and the one corruption here that the header and row
    # checks cannot see.
    'tangling_slope':
        (lambda l: l[:2] + [f'{i*600.0:.7e}\t1.5\t0.0\n'
                            for i in range(len(l) - 2)],
         'climbs a full fault-normal cell'),
}


def buildEqdyna():
    """Build ONLY src/eqdyna (explicit path) -- never touch bin/eqdyna."""
    env = dict(os.environ)
    env['MACHINE'] = MACHINE
    r = subprocess.run(['make'], cwd=SRC, env=env, capture_output=True,
                       text=True, timeout=600)
    if r.returncode != 0 or not os.path.exists(os.path.join(SRC, 'eqdyna')):
        raise RuntimeError(
            f'src/ build failed (exit {r.returncode}); is mpif90 '
            f'(MACHINE={MACHINE}) on PATH?\n'
            + r.stdout[-3000:] + r.stderr[-3000:])


def makeCase(tmp, name):
    import lib
    import numpy as np
    case = os.path.join(tmp, f'case_{name}')
    os.makedirs(case)
    for f in ('case.setup', 'defaultParameters.py', 'lib.py'):
        shutil.copy(os.path.join(ROOT, 'scripts', f), case)
    with open(os.path.join(case, 'user_defined_params.py'), 'w') as fh:
        fh.write(USER_PARAMS)
    # A flat surface: a perfectly valid rough-fault file whose mesh morph is
    # the identity, so the clean run is a plain planar-fault run.
    zeros = np.zeros((NNZ, NNX))
    lib.writeFaultRoughGeometry(
        zeros, zeros, zeros, 500.0, -1000.0, -2000.0,
        fname=os.path.join(case, lib.FAULT_ROUGH_GEOMETRY_FILE))
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case,
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f'case.setup failed for {name}:\n{r.stdout}\n{r.stderr}')
    return case


def corrupt(case, mutate):
    path = os.path.join(case, 'bFault_Rough_Geometry.txt')
    with open(path) as f:
        lines = f.readlines()
    with open(path, 'w') as f:
        f.writelines(mutate(lines))


def runCase(case):
    eqdyna = os.path.join(SRC, 'eqdyna')
    r = subprocess.run([MPIRUN, '-np', '1', eqdyna], cwd=case,
                       capture_output=True, text=True, timeout=300)
    return r.returncode, r.stdout + r.stderr


def main():
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        print(f'FAIL test_fault_geometry_guard: mpif90/{MPIRUN} not on PATH '
              '(this guard requires a real MPI build, it does not skip silently)')
        return 1

    fails = []
    try:
        buildEqdyna()
    except RuntimeError as e:
        print('FAIL test_fault_geometry_guard')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testFaultGeometryGuard.')
    try:
        # 1. a valid file must be accepted
        try:
            clean = makeCase(tmp, 'clean')
            rc, out = runCase(clean)
            if rc != 0:
                fails.append(f'clean: a VALID rough-geometry file was rejected '
                             f'(exit {rc}):\n{out[-2000:]}')
            if GUARD_MSG in out:
                fails.append('clean: the guard fired on a valid file')
        except RuntimeError as e:
            fails.append(str(e))

        # 2. every corruption must be caught, loudly and specifically
        for name, (mutate, expect) in CORRUPTIONS.items():
            try:
                case = makeCase(tmp, name)
            except RuntimeError as e:
                fails.append(str(e))
                continue
            corrupt(case, mutate)
            rc, out = runCase(case)
            if rc == 0:
                fails.append(f'{name}: ran to completion on a corrupted '
                             f'geometry file (exit 0) -- this is the original '
                             f'silent-failure mode')
                continue
            if GUARD_MSG not in out:
                fails.append(f'{name}: exited {rc} but never printed the guard '
                             f'message, so the cause is not discoverable:\n'
                             f'{out[-1500:]}')
                continue
            if expect not in out:
                fails.append(f'{name}: guard fired but did not name the defect '
                             f'(expected text {expect!r} in the message)')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_fault_geometry_guard')
        for f in fails:
            print(' -', f)
        return 1
    print(f'SUCCESS test_fault_geometry_guard: a valid rough-geometry file is '
          f'accepted and all {len(CORRUPTIONS)} corruption(s) '
          f'({", ".join(CORRUPTIONS)}) stop the run with a message naming the '
          f'defect')
    return 0


if __name__ == '__main__':
    sys.exit(main())
