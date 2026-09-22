#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the EQDYNA_DUMP_EQUIL
equilibrium diagnostic in src/fortran/driver.f90 (dumpNodalAccel).

Background: pathway item 19(b)/39 was settled with a three-way control
showing TPV30 has NO fault-local equilibrium defect (rough geometry mean
1.949e-01 / max 1.448e+01 m/s^2; roughness zeroed 3.5e-14; rough + fault
locked 5.74e-14). The instrument behind those numbers is dumpNodalAccel:
at nt==1, after faulting and the mass divide, nodalForceArr IS the
residual acceleration of the initial state (velArr/dispArr still zero).
It writes one equilibriumDump.<rank>.txt per MPI rank -- per node:
meshCoor(1:3), residual acceleration(1:3), a fault tag (0 off-fault,
1 master, 2 slave), and the node's dof count (so 12-dof PML nodes can
be excluded in analysis) -- gated on the environment variable
EQDYNA_DUMP_EQUIL being EXACTLY '1'. Unset, or set to anything else,
must write NOTHING (fails-closed). Nothing else in testsys/ ever sets
EQDYNA_DUMP_EQUIL, so without this guard the whole code path would ship
untested: a deleted `call dumpNodalAccel`, a broken column, or a leak
of the file into ordinary runs would keep every tier green.

The diagnostic is Fortran-only by design (it instruments the Fortran
solver; there is no Python-port counterpart and none is wanted).

Fixture: the same tiny hand-written flat-fault case as
test_fault_mpi_boundary_arn.py (12x12x6 cells, dx=500, embedded here as
its own copy so the two guards stay independent), term = 2*dt since the
dump fires at nt==1. Four runs, ~1 s each once src/eqdyna is built:

  1. serial, env UNSET -> exit 0, NO equilibriumDump.* file
  2. serial, env='0'   -> exit 0, NO file (fails-closed: only exactly
     '1' enables the dump)
  3. serial, env='1'   -> exactly equilibriumDump.0000.txt with:
       - row count == EXPECTED_TOTAL_NODES (8140), the run's node
         count. Provenance: measured at landing on this frozen fixture;
         the solver itself aborts in meshGenError if the dump loop's
         bound (totalNumOfNodes) disagrees with countMeshEntities'
         independent tally, so rows==8140 pins "one row per node"
         against a number the solver cross-checks internally.
       - every row has exactly 8 whitespace-separated fields, cols 1-6
         finite floats, col 7 in {0,1,2}, col 8 a positive int
       - fault-tagged rows: tag==1 count == tag==2 count == 15
         == ((fxmax-fxmin)/dx+1)*((fzmax-fzmin)/dz+1), the analytic
         5x3 fault grid from the fixture's own geometry
       - tag==1 count also == line count of frt.txt0, the on-fault
         output written by a different code path (faulting.f90)
       - unique (x,y,z) coords == rows - 15: the only duplicated
         coordinates are the 15 split-node pairs
  4. np=2 x-split (2,1,1), env='1' -> exactly equilibriumDump.0000.txt
     and .0001.txt; per rank: 8 fields per row, tag==1 count == tag==2
     count == that rank's frt.txt<r> line count, and both > 0.

Every verdict prints the counts it compared (rows, columns, tags) --
never a bare pass/fail. Fails loudly (rule 2) if mpif90/mpirun/the
src/ build are unavailable -- never silently skips or passes. Builds
ONLY src/eqdyna (never bin/eqdyna -- other processes on this box may
be running bin/eqdyna concurrently).

Tier: regression -- this tier already drives four real mpirun runs of
the same tiny fixture (test_fault_mpi_boundary_arn.py) at well under a
second per run once src/eqdyna is built, and this guard exists because
of a past incident shape (four green gates that exercised nothing).
"""
import glob
import math
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'src', 'fortran')
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

# Fixture geometry constants the assertions derive from (keep in sync with
# USER_PARAMS below).
DX = 500.0
FXMIN, FXMAX = -1.0e3, 1.0e3
FZMIN, FZMAX = -1.0e3, 0.0e3
EXPECTED_FAULT_NODES = round((FXMAX - FXMIN) / DX + 1) * round((FZMAX - FZMIN) / DX + 1)  # 5*3
# The run's node count for this frozen fixture (see docstring for provenance).
EXPECTED_TOTAL_NODES = 8140

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
par.term = 2.*par.dt
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


def run_case(case, n, dump_env):
    """dump_env: None -> EQDYNA_DUMP_EQUIL absent; str -> set to that value."""
    env = dict(os.environ)
    env.pop('EQDYNA_DUMP_EQUIL', None)
    if dump_env is not None:
        env['EQDYNA_DUMP_EQUIL'] = dump_env
    eqdyna = os.path.join(SRC, 'eqdyna')
    r = subprocess.run([MPIRUN, '-np', str(n), eqdyna], cwd=case, env=env,
                       capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr


def dump_files(case):
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(case, 'equilibriumDump.*')))


def check_dump_content(case, rank, fails, expect_rows=None):
    """Assert structure/content of one rank's dump. Returns (rows, tag1, tag2)."""
    label = f'equilibriumDump.{rank:04d}.txt'
    path = os.path.join(case, label)
    with open(path) as fh:
        lines = fh.read().splitlines()
    rows = len(lines)
    tag1 = tag2 = 0
    coords = set()
    for ln, line in enumerate(lines, 1):
        f = line.split()
        if len(f) != 8:
            fails.append(f'{label} line {ln}: {len(f)} fields, expected 8: {line!r}')
            return rows, tag1, tag2
        try:
            vals = [float(x) for x in f[:6]]
            tag = int(f[6])
            ndof = int(f[7])
        except ValueError:
            fails.append(f'{label} line {ln}: non-numeric field: {line!r}')
            return rows, tag1, tag2
        if any(math.isnan(v) or math.isinf(v) for v in vals):
            fails.append(f'{label} line {ln}: NaN/Inf in coords/accel: {line!r}')
            return rows, tag1, tag2
        if tag not in (0, 1, 2):
            fails.append(f'{label} line {ln}: fault tag {tag} not in {{0,1,2}}')
        if ndof <= 0:
            fails.append(f'{label} line {ln}: non-positive dof count {ndof}')
        tag1 += tag == 1
        tag2 += tag == 2
        coords.add(tuple(f[:3]))
    if expect_rows is not None and rows != expect_rows:
        fails.append(f'{label}: {rows} rows != expected node count {expect_rows}')
    if tag1 == 0:
        fails.append(f'{label}: ZERO fault-tagged (tag==1) rows -- the tagger is dead')
    if tag1 != tag2:
        fails.append(f'{label}: tag==1 count {tag1} != tag==2 count {tag2} (split pairs broken)')
    if len(coords) != rows - tag2:
        fails.append(f'{label}: {len(coords)} unique coords, expected rows-{tag2}={rows - tag2} '
                     '(duplicated coordinates must be exactly the split-node pairs)')
    # frt.txt<rank>: on-fault output written by faulting.f90, an independent
    # code path from the dump's nsmp-based tagger.
    frt = os.path.join(case, f'frt.txt{rank}')
    if not os.path.exists(frt):
        fails.append(f'{label}: companion frt.txt{rank} missing, cannot cross-check fault count')
    else:
        with open(frt) as fh:
            nfrt = sum(1 for _ in fh)
        if tag1 != nfrt:
            fails.append(f'{label}: tag==1 count {tag1} != frt.txt{rank} row count {nfrt}')
    return rows, tag1, tag2


def main():
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        print(f'FAIL test_equilibrium_dump: mpif90/{MPIRUN} not on PATH '
              '(this guard requires a real MPI build, it does not skip silently)')
        return 1

    fails = []
    try:
        build_eqdyna()
    except RuntimeError as e:
        print('FAIL test_equilibrium_dump')
        print(' -', e)
        return 1

    tmp = tempfile.mkdtemp(prefix='testEquilibriumDump.')
    evidence = []
    try:
        # 1. env UNSET -> no file (fresh case dir: no stale dump can satisfy this).
        case = make_case(tmp, 'unset', 1, 1, 1)
        rc, out = run_case(case, 1, None)
        files = dump_files(case)
        if rc != 0:
            fails.append(f'env-unset run exited {rc}\n{out[-2000:]}')
        if files:
            fails.append(f'env-unset run produced dump file(s) {files} -- the gate leaks')
        evidence.append(f'env-unset: exit {rc}, {len(files)} dump files (want 0)')

        # 2. env='0' -> no file (only exactly '1' enables).
        case = make_case(tmp, 'zero', 1, 1, 1)
        rc, out = run_case(case, 1, '0')
        files = dump_files(case)
        if rc != 0:
            fails.append(f"env='0' run exited {rc}\n{out[-2000:]}")
        if files:
            fails.append(f"env='0' run produced dump file(s) {files} -- gate must require exactly '1'")
        evidence.append(f"env='0': exit {rc}, {len(files)} dump files (want 0)")

        # 3. serial, env='1' -> one file, full content assertions.
        case = make_case(tmp, 'serial_on', 1, 1, 1)
        rc, out = run_case(case, 1, '1')
        files = dump_files(case)
        if rc != 0:
            fails.append(f"serial env='1' run exited {rc}\n{out[-2000:]}")
        elif files != ['equilibriumDump.0000.txt']:
            fails.append(f"serial env='1': dump files {files}, expected exactly "
                         "['equilibriumDump.0000.txt']")
        else:
            rows, t1, t2 = check_dump_content(case, 0, fails, expect_rows=EXPECTED_TOTAL_NODES)
            if t1 != EXPECTED_FAULT_NODES:
                fails.append(f'serial: tag==1 count {t1} != analytic fault-grid count '
                             f'{EXPECTED_FAULT_NODES} (5x3 nodes from fixture geometry)')
            evidence.append(f"serial env='1': exit {rc}, rows {rows} (want {EXPECTED_TOTAL_NODES}), "
                            f'8 cols/row, tag1 {t1} tag2 {t2} (want {EXPECTED_FAULT_NODES} each)')

        # 4. np=2 x-split, env='1' -> one file per rank, per-rank structure.
        case = make_case(tmp, 'mpi_on', 2, 1, 1)
        rc, out = run_case(case, 2, '1')
        files = dump_files(case)
        if rc != 0:
            fails.append(f"np=2 env='1' run exited {rc}\n{out[-2000:]}")
        elif files != ['equilibriumDump.0000.txt', 'equilibriumDump.0001.txt']:
            fails.append(f"np=2 env='1': dump files {files}, expected exactly one per rank "
                         "['equilibriumDump.0000.txt', 'equilibriumDump.0001.txt']")
        else:
            for rank in (0, 1):
                rows, t1, t2 = check_dump_content(case, rank, fails)
                evidence.append(f'np=2 rank {rank}: rows {rows}, 8 cols/row, tag1 {t1} tag2 {t2}, '
                                f'tag1 cross-checked against frt.txt{rank} row count')
    except RuntimeError as e:
        fails.append(str(e))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for line in evidence:
        print('  .', line)
    if fails:
        print('FAIL test_equilibrium_dump')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_equilibrium_dump (env gate fails-closed; serial dump '
          f'{EXPECTED_TOTAL_NODES} rows x 8 cols, {EXPECTED_FAULT_NODES}+{EXPECTED_FAULT_NODES} '
          'fault-tagged rows cross-checked vs geometry and frt.txt; per-rank files at np=2)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
