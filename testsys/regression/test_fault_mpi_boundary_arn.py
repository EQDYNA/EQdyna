#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the fault-plane/MPI-partition
-boundary arn-doubling bug.

Background: MPI4arn (src/meshgen.f90) accumulates arn (fault-node tributary
area) from a purely 2D fault-surface quad grid that every rank touching the
fault's plane builds independently and completely. EQdyna's faults always sit
in the x-z plane (fymin==fymax==0). When an MPI partition boundary in y
coincides with that plane, every rank on the boundary already holds the
COMPLETE local arn there, so syncArnBoundary's cross-rank add-back used to
DOUBLE it (a genuine partial-sum recombination for x/z splits -- which
genuinely divide the fault -- but a duplicate for y). Every on-fault traction
term divides by arn (faulting.f90), so tractions came out at EXACTLY HALF
their correct value whenever an MPI partition plane coincided with the fault
plane (e.g. a symmetric y-domain split evenly by npy). See pathway_forward.md
item 26 for the fix's full history.

Fix: syncArnBoundary now skips arn's add-back specifically when fltxyz shows
zero nominal fault extent in that boundary's direction (the duplicate case),
while still setting fltMPI(k) and still exchanging data, since fnms/
nodalMassArr's addFaultBoundaryTerm (assembleGlobalMass.f90) depends on that
regardless of what arn does with it. A defense-in-depth hard stop,
checkFaultMPIAlignment (eqdyna3d.f90), still fires for the one residual case
this fix does not reason about: a y-boundary fault with NON-degenerate
y-extent.

This test builds ONE tiny hand-written flat-fault case (symmetric y-domain --
the hazard condition) and runs it at four decompositions that hold the mesh
geometry fixed and vary ONLY the processor grid (decomposition-invariance
shape, as prescribed):
  serial  (1,1,1) -- baseline
  xsplit  (2,1,1) -- genuinely divides the fault; unaffected either way
  zsplit  (1,1,2) -- genuinely divides the fault; unaffected either way
  ysplit  (1,2,1) -- THE hazard: an MPI boundary coincides with the fault

Asserts:
  1. All four exit 0 -- the gate must NOT block the now-fixed path.
  2. The gate's own message never appears in stdout for any of the four.
  3. The hypocenter's final on-fault (normal, strike, dip) traction, parsed
     from each run's stdout, is IDENTICAL across all four -- decomposition-
     invariance is what actually proves the fix landed, not just the gate
     staying quiet.

EXACT equality (not a tolerance) is the gate here: with the fix in place,
serial==ysplit exactly, reproduced at 8x this test's term with zero
difference at every duration tried (see pathway_forward.md item 26).

Cheap (rule 9): ~12x12x6-cell mesh, term = 5 dt, well under a second per run
once src/eqdyna is built. Builds ONLY src/eqdyna (never bin/eqdyna -- other
processes on this box may be running bin/eqdyna concurrently).
Fails loudly (rule 2) if mpif90/mpirun/the src/ build are unavailable --
never silently skips or passes.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'src', 'fortran')
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

USER_PARAMS = '''#! /usr/bin/env python3
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -3.0e3, 3.0e3
par.ymin, par.ymax = -3.0e3, 3.0e3     # SYMMETRIC y -> fault-plane-on-boundary hazard
par.zmin, par.zmax = -3.0e3, 0.0e3

par.fxmin, par.fxmax = -1.0e3, 1.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -1.0e3, 0.0e3

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

par.nucR = 300.

par.nx = __NX__
par.ny = __NY__
par.nz = __NZ__

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

CONFIGS = {
    'serial': (1, 1, 1),
    'xsplit': (2, 1, 1),
    'zsplit': (1, 1, 2),
    'ysplit': (1, 2, 1),
}

TRACT_RE = re.compile(
    r'n,s,d tract\s+\(MPa\)\s*(-?[0-9.]+E[+-][0-9]+)\s*(-?[0-9.]+E[+-][0-9]+)\s*(-?[0-9.]+E[+-][0-9]+)')

GATE_MSG = 'checkFaultMPIAlignment'


def build_eqdyna():
    """Build ONLY src/eqdyna (explicit path) -- never touch bin/eqdyna."""
    env = dict(os.environ)
    env['MACHINE'] = MACHINE
    r = subprocess.run(['make'], cwd=SRC, env=env, capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not os.path.exists(os.path.join(SRC, 'eqdyna')):
        raise RuntimeError(
            f'src/ build failed (exit {r.returncode}); is mpif90 (MACHINE={MACHINE}) on PATH?\n'
            + r.stdout[-3000:] + r.stderr[-3000:])


def make_case(tmp, name, nx, ny, nz):
    case = os.path.join(tmp, f'case_{name}')
    os.makedirs(case)
    for f in ('case.setup', 'defaultParameters.py', 'lib.py'):
        shutil.copy(os.path.join(ROOT, 'scripts', f), case)
    params = USER_PARAMS.replace('__NX__', str(nx)).replace('__NY__', str(ny)).replace('__NZ__', str(nz))
    with open(os.path.join(case, 'user_defined_params.py'), 'w') as fh:
        fh.write(params)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise RuntimeError(f'case.setup failed for {name} ({nx},{ny},{nz}):\n{r.stdout}\n{r.stderr}')
    return case


def run_case(case, n):
    eqdyna = os.path.join(SRC, 'eqdyna')
    r = subprocess.run([MPIRUN, '-np', str(n), eqdyna], cwd=case,
                       capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr


def last_traction(output):
    m = TRACT_RE.findall(output)
    if not m:
        return None
    return tuple(float(x) for x in m[-1])


def main():
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        print(f'FAIL test_fault_mpi_boundary_arn: mpif90/{MPIRUN} not on PATH '
              '(this guard requires a real MPI build, it does not skip silently)')
        return 1

    fails = []
    try:
        build_eqdyna()
    except RuntimeError as e:
        print('FAIL test_fault_mpi_boundary_arn')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testFaultMPIBoundary.')
    results = {}
    try:
        for name, (nx, ny, nz) in CONFIGS.items():
            try:
                case = make_case(tmp, name, nx, ny, nz)
            except RuntimeError as e:
                fails.append(str(e))
                continue
            rc, out = run_case(case, nx * ny * nz)
            if rc != 0:
                fails.append(f'{name} ({nx},{ny},{nz}): exited {rc} '
                              f'(gate fired on a case it should not?)\n{out[-2000:]}')
                continue
            if GATE_MSG in out:
                fails.append(f'{name}: gate message present even though exit was 0 (inconsistent)')
            tr = last_traction(out)
            if tr is None:
                fails.append(f'{name}: no hypocenter traction block found in stdout')
            results[name] = tr
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    base = results.get('serial')
    if base is None:
        fails.append('serial run produced no traction to compare against')
    else:
        for name in ('xsplit', 'zsplit', 'ysplit'):
            tr = results.get(name)
            if tr is not None and tr != base:
                fails.append(
                    f'{name} hypocenter traction {tr} != serial {base} '
                    '(decomposition-dependent result -- the arn-doubling bug is back)')

    if fails:
        print('FAIL test_fault_mpi_boundary_arn')
        for f in fails:
            print('  -', f)
        return 1
    print(f'SUCCESS test_fault_mpi_boundary_arn (serial==xsplit==zsplit==ysplit traction {base} MPa)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
