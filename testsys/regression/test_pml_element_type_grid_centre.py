#! /usr/bin/env python3
"""
Defect A (resolution README, docs/evidence/barall-overlays-2026-10-09/
resolution/, owner decision 2026-10-11 board PR #178): an element must be
typed PML (type 2) by the centre of its GRID-LINE cell, the same unblended
coordinates setNumDof uses to give a node 12 dofs -- not by the mean of its
y-blended meshCoor.

insertFaultInterface's y-blend (func_lib.f90) stretches the whole footwall
from the y-min PML face to a dipping fault whose y-offset (`peak`) reaches
7.5 km on test.tpv10/12/13. Typing elements by the blended centre made
regular (type 1) elements in the y-min PML interface own 12-dof nodes: their
y/z force went into split x-dofs and the Method-2 initial stress there was
unbalanced. Measured at dx 500 m, nt==1 residual (EQDYNA_DUMP_EQUIL=1): up
to ~90 m/s^2 on one slanted node column before, 2.9e-2 m/s^2 (the generic
free-surface value) after. Outcome: TPV13 vs FaultMod at 8 s, slip ratio
0.61 -> 0.95, body peak-velocity ratio 4.24 -> 1.01, p90 |dt| 2.84 -> 0.53 s.

RULE 14a, both ways -- each half asserts BEHAVIOUR (rule 10a):
  python   PASS: the port's build_elements types a coarse test.tpv13 mesh
           with ZERO regular elements holding a 12-dof node.
           FAIL: the same invariant checker, fed the historical
           blended-centre typing of the SAME mesh, must find > 0 (every run).
  fortran  PASS: src/fortran/eqdyna, run at 4 ranks on the same coarse case
           with EQDYNA_DUMP_EQUIL=1, writes an nt==1 residual whose max over
           NON-fault nodes is below RESIDUAL_MAX_MS2 (an initial state in
           equilibrium up to the free-surface truncation residual).
           FAIL: `--mutant` builds a copy of src/fortran in a tempdir with
           createElement's PML test switched back to elementCenterCoor and
           runs the SAME check on it; it must go red. Not run by default
           only for cost (a second full Fortran build, ~+11 s); the mutant is
           constructed by text substitution and refuses loudly if the
           substitution site has moved, but nothing is ASSERTED on source
           text. Measured 2026-10-11 at dx 1500 m, 4 ranks, term 0.5 s:
           fixed max non-fault |a| 0.0695 m/s^2, 0 nodes > 1 m/s^2;
           mutant 131.2 m/s^2, 880 nodes > 1 m/s^2 (all 12-dof, in the
           y-min PML interface). The e2e gates on test.tpv10/12/13 (whose
           references were regenerated on this fix) also catch a revert of
           either implementation.

Cheap (rule 9): one case.setup at dx 1500 m, the port's mesh builders, an
incremental `make` of src/fortran (never bin/eqdyna, which other processes
on this box may be running) and one ~2 s 4-rank solver launch of 4 steps:
~11 s wall here, ~22 s with --mutant (measured 2026-10-11).
"""
import os
import glob
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
sys.path.insert(0, ROOT)
from testsys.common import MACHINE  # noqa: E402
SRC = os.path.join(ROOT, 'src', 'fortran')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE = 'test.tpv13'
COARSE_DX = 1500.0
COARSE_TERM = 0.5   # dt 0.114 s at dx 1500 m -> 4 steps; the dump is at nt==1
NRANKS = 4          # the case's own par.nx*ny*nz (2,2,1)
# Fixed 0.0695, mutant 131.2 (see docstring): 1 m/s^2 sits 14x above the
# fixed tree's free-surface truncation residual and 131x below the defect.
RESIDUAL_MAX_MS2 = 1.0
# The fixed createElement PML test, and the historical one the mutant restores.
FIXED_TEST = 'if (gridCenterCoor(1)>PMLb(1).or.gridCenterCoor(1)<PMLb(2) &\n' \
             '    .or.gridCenterCoor(2)>PMLb(3).or.gridCenterCoor(2)<PMLb(4) &\n' \
             '    .or.gridCenterCoor(3)<PMLb(5)) then'


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'), env.get('PATH', '')])
    return env


def _coarse_case(tmp, term=None):
    case_dir = os.path.join(tmp, 'case')
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                    case_dir, CASE], env=_env(), check=True, capture_output=True)
    p = os.path.join(case_dir, 'user_defined_params.py')
    text = open(p).read()
    text, n = re.subn(r'^par\.dx = [0-9.]+', 'par.dx = %r' % COARSE_DX, text, flags=re.M)
    if n != 1:
        raise RuntimeError('could not override par.dx in %s' % p)
    if term is not None:
        text, n = re.subn(r'^par\.term = [0-9.]+', 'par.term = %r' % term, text, flags=re.M)
        if n != 1:
            raise RuntimeError('could not override par.term in %s' % p)
    open(p, 'w').write(text)
    subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=_env(),
                   check=True, capture_output=True)
    return case_dir


def mixed_elements(conn, is_pml_elem, num_dof):
    """Number of regular elements that hold at least one 12-dof node."""
    has12 = (num_dof[conn] == 12).any(axis=1)
    return int((~is_pml_elem & has12).sum())


def check_python():
    from eqdyna import readInputFiles, meshgen
    with tempfile.TemporaryDirectory() as tmp:
        cd = _coarse_case(tmp)
        params, g = readInputFiles.build_params(cd)
        material = readInputFiles.read_bmaterial(os.path.join(cd, 'bMaterial.txt'),
                                                 g['nmat'], g['n2mat'])
    xl, yl, zl, pmlb, _ = meshgen.build_grid_lines(params)
    mc, _, nsmp = meshgen.build_node_coordinates(xl, yl, zl, params)
    conn, et, _, _ = meshgen.build_elements(xl, yl, zl, params, pmlb, nsmp, material, mc)
    nd = np.asarray(meshgen.build_equation_numbers(xl, yl, zl, params, pmlb)[0])
    if nd.shape[0] == mc.shape[0] - 1:
        nd = np.concatenate([[0], nd])
    conn = np.asarray(conn)
    n12 = int((nd == 12).sum())
    new_bad = mixed_elements(conn, np.asarray(et) == 2, nd)
    c = mc[conn].mean(axis=1)  # historical: blended-centre typing
    old_pml = ((c[:, 0] > pmlb['xmax0']) | (c[:, 0] < pmlb['xmin0']) |
               (c[:, 1] > pmlb['ymax0']) | (c[:, 1] < pmlb['ymin0']) |
               (c[:, 2] < pmlb['zmin0']))
    old_bad = mixed_elements(conn, old_pml, nd)
    print('  python: %d elements, %d 12-dof nodes; regular elements holding a '
          '12-dof node: fixed=%d, historical blended-centre=%d'
          % (len(et), n12, new_bad, old_bad))
    fails = []
    if n12 == 0:
        fails.append('no 12-dof nodes at all -- the check tested nothing')
    if new_bad != 0:
        fails.append('build_elements leaves %d regular elements holding a 12-dof node' % new_bad)
    if old_bad == 0:
        fails.append('the checker found 0 mixed elements under the historical typing too '
                     '-- it cannot come out both ways on this mesh')
    return fails


def build(src_dir):
    """Incremental `make` of src_dir (a src/fortran tree); returns the binary."""
    env = dict(os.environ)
    env['MACHINE'] = MACHINE
    r = subprocess.run(['make'], cwd=src_dir, env=env, capture_output=True,
                       text=True, timeout=600)
    exe = os.path.join(src_dir, 'eqdyna')
    if r.returncode != 0 or not os.path.exists(exe):
        raise RuntimeError('build of %s failed (exit %d)\n%s%s'
                           % (src_dir, r.returncode, r.stdout[-2000:], r.stderr[-2000:]))
    return exe


def make_mutant(tmp):
    """A copy of src/fortran whose createElement types PML by the blended
    elementCenterCoor again, built in tmp/src/fortran (the makefile calls
    ../../scripts/src_hash.py, hence the scripts symlink)."""
    dst = os.path.join(tmp, 'src', 'fortran')
    shutil.copytree(SRC, dst, ignore=shutil.ignore_patterns('*.o', '*.mod', 'eqdyna'))
    os.symlink(os.path.join(ROOT, 'scripts'), os.path.join(tmp, 'scripts'))
    p = os.path.join(dst, 'meshgen.f90')
    text = open(p).read()
    if text.count(FIXED_TEST) != 1:
        raise RuntimeError('--mutant: the fixed createElement PML test is not '
                           'where this test expects it; update FIXED_TEST')
    open(p, 'w').write(text.replace(FIXED_TEST, FIXED_TEST.replace(
        'gridCenterCoor', 'elementCenterCoor')))
    return build(dst)


def residual(exe, tmp):
    """(n nodes, n 12-dof, max non-fault |a|, n non-fault |a| > bound) of the
    nt==1 EQDYNA_DUMP_EQUIL residual of a coarse test.tpv13 run of exe."""
    cd = _coarse_case(tmp, term=COARSE_TERM)
    env = _env()
    env['EQDYNA_DUMP_EQUIL'] = '1'
    env['OMP_NUM_THREADS'] = '1'
    r = subprocess.run([MPIRUN, '-np', str(NRANKS), exe], cwd=cd, env=env,
                       capture_output=True, text=True, timeout=300)
    files = sorted(glob.glob(os.path.join(glob.escape(cd), 'equilibriumDump.*.txt')))
    if r.returncode != 0 or len(files) != NRANKS:
        raise RuntimeError('%s exited %d with %d of %d dump files\n%s'
                           % (exe, r.returncode, len(files), NRANKS,
                              (r.stdout + r.stderr)[-2000:]))
    a = np.vstack([np.loadtxt(f, ndmin=2) for f in files])
    acc = np.linalg.norm(a[:, 3:6], axis=1)
    off = a[:, 6] == 0   # fault-node tag 0: not a split-node fault node
    return (len(a), int((a[:, 7] == 12).sum()), float(acc[off].max()),
            int((acc[off] > RESIDUAL_MAX_MS2).sum()))


def check_fortran_binary(exe, label):
    with tempfile.TemporaryDirectory() as tmp:
        n, n12, amax, nbad = residual(exe, tmp)
    print('  fortran %s: %d nodes, %d 12-dof; nt==1 non-fault max|a| = %.4g m/s^2, '
          '%d nodes > %g' % (label, n, n12, amax, nbad, RESIDUAL_MAX_MS2))
    fails = []
    if n12 == 0:
        fails.append('fortran %s: no 12-dof nodes -- the check tested nothing' % label)
    if nbad:
        fails.append('fortran %s: %d non-fault nodes out of equilibrium at nt==1 '
                     '(max %.4g m/s^2 > %g)' % (label, nbad, amax, RESIDUAL_MAX_MS2))
    return fails


def check_fortran(mutant=False):
    if shutil.which('mpif90') is None or shutil.which(MPIRUN) is None:
        return ['mpif90/%s not on PATH (this guard needs a real MPI build; it '
                'does not skip silently)' % MPIRUN]
    fails = check_fortran_binary(build(SRC), 'src/fortran')
    if mutant:
        with tempfile.TemporaryDirectory() as tmp:
            bad = check_fortran_binary(make_mutant(tmp), 'MUTANT (blended-centre typing)')
        if not bad:
            fails.append('the residual check passed the blended-centre mutant -- '
                         'it cannot come out both ways')
    return fails


def main():
    fails = check_python() + check_fortran(mutant='--mutant' in sys.argv[1:])
    if fails:
        print('\nFAIL test_pml_element_type_grid_centre (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_pml_element_type_grid_centre')
    return 0


if __name__ == '__main__':
    sys.exit(main())
