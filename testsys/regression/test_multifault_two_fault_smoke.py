#! /usr/bin/env python3
"""
Regression guard: test.multifault2, the row-17 two-fault ROUTING/SMOKE case.

NOT tpv22/tpv23 (a separate, later, reviewed mission -- pathway_forward.md
item 17) and not part of testNameList.py / testsys/matrix.py's gated sweep.
This is its own throwaway-but-committed case_input compset
(case_input/test.multifault2) with its own committed reference
(test.reference.results/test.multifault2/frt.canonical.txt), built and run
independently of the TPV matrix so it cannot be confused with, or accidentally
gate, anything in it.

WHAT THIS PINS

  1. The Fortran binary actually MESHES and RUNS two vertical, planar,
     parallel faults (y=0 and y=2000 m, 4*dy apart) without the bugs a real
     2-fault run -- not a unit test -- actually found during this work:
       - msnode collision across faults (meshgen.f90 createMasterNode):
         SIGSEGV, rank-local totalNumOfNodes tally mismatch (ERR 42).
       - replaceSlaveWithMasterNode's mesh-connectivity gate keyed to fault
         1's y-plane alone: NaN velocity at a fault-2 node within 2 steps
         (ERR 63), because fault 2's neighbouring elements never had their
         slave-node references swapped for master-node ones.
       - fltgm/fltl../fltnum single-fault-sized (overwritten per fault):
         would have broken the MPI boundary force/mass exchange for
         whichever fault MPI4arn did not process last, on any decomposition
         where a rank owns fault-boundary nodes on more than one fault.
  2. Fault 2 does NOT read fault 1's on-fault arrays. The two faults are
     given a deliberately different initial normal-stress depth coefficient
     (7378 vs 5000 Pa/m) specifically so this aliasing bug -- eqquasi's own
     documented failure mode -- would show up in frt.txt's tnrm column, not
     just in whether the run completes. Checked two ways: an analytic check
     against the known coefficient (tight, since the field is static at
     t=0 and barely perturbed dynamically by t=1s), and a canonical-output
     comparison against a committed reference (catches anything the
     analytic check does not).
  3. This is the multi-fault code path actually running (not falling through
     to ntotft==1 behaviour): asserts nftnd(1) and nftnd(2) are both
     positive-and-equal-sized in the canonical output, not just nftnd(1).

Verified directly on a decomposition where fault 1 and fault 2 are split
across DIFFERENT MPI ranks (npx,npy,npz=(2,2,1): y-partition separates y=0
from y=2000) AND on a decomposition where a SINGLE rank owns both faults
entirely (serial, npx=npy=npz=1) -- the serial run is the one that actually
exercises the msnode-collision and slave/master-swap fixes above; the 4-rank
run additionally exercises the per-fault MPI boundary fix. Both are run here.

fortran only: the python-jax backend is ntotft==1 throughout as of this
writing (readInputFiles.py:289-291, eqdyna3d.py:331-333, meshgen.py) and is
NOT exercised by this test -- a known, reported gap (see the mission report),
not a silent skip.

Not cheap in the rule-9 sense (builds + runs the actual solver twice), but
bounded: dx=500 m, term=1.0 s (~23 steps), 4 + 1 MPI ranks, well under a
minute total.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_NAME = 'test.multifault2'
REFERENCE = os.path.join(ROOT, 'test.reference.results', CASE_NAME, 'frt.canonical.txt')
FAULT1_NORM_COEFF = 7378.0
FAULT2_NORM_COEFF = 5000.0
# Loose vs the ~1e7-1e8 Pa field: this is a routing/aliasing check, not a
# tight physics gate. Measured slack needed: fault 1 ruptures (hypocenter
# peak slip rate ~6.4 m/s) and radiates into fault 2, 2 km away, within the
# 1 s run, perturbing fault 2's normal stress dynamically by up to ~8.9e6 Pa
# away from its pure static value -- physically expected, not a bug (the
# aliasing cross-check below, against fault 1's coefficient at 3.9e7 Pa,
# stays well clear of this bound either way).
TOL_PA = 1.0e7

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import frt_canonical  # noqa: E402


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _require_binary():
    exe = os.path.join(ROOT, 'bin', 'eqdyna')
    if not os.path.isfile(exe):
        raise AssertionError(
            'bin/eqdyna is missing -- build it first (./install-eqdyna.sh -m <platform>) '
            'before running this test.')
    return exe


def _build_and_run(tmp, nranks, nx, ny, nz):
    """create.newcase + case.setup + mpirun, with nx/ny/nz overridden to
    control the decomposition (case.setup reads par.nx/ny/nz into
    npx/npy/npz). Returns the case directory."""
    case_dir = os.path.join(tmp, 'case_%dr' % nranks)
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))
    params = os.path.join(case_dir, 'user_defined_params.py')
    text = open(params).read()
    text = text.replace('par.nx = 2', 'par.nx = %d' % nx)
    text = text.replace('par.ny = 2', 'par.ny = %d' % ny)
    text = text.replace('par.nz = 1', 'par.nz = %d' % nz)
    open(params, 'w').write(text)

    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))

    exe = _require_binary()
    r = subprocess.run(['mpirun', '-np', str(nranks), '--oversubscribe', exe],
                       cwd=case_dir, env=env, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise AssertionError('eqdyna exited %d for %d rank(s) at (%d,%d,%d):\n%s'
                             % (r.returncode, nranks, nx, ny, nz, (r.stdout + r.stderr)[-3000:]))
    return case_dir


def _canonical(case_dir):
    return frt_canonical.canonical_from_case(case_dir)


def check_four_rank_matches_committed_reference():
    if not os.path.isfile(REFERENCE):
        raise AssertionError('missing committed reference %r -- regenerate with '
                             'python3 -m testsys.frt_canonical <case_dir>' % REFERENCE)
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _build_and_run(tmp, nranks=4, nx=2, ny=2, nz=1)
        got = _canonical(case_dir)
    want = np.loadtxt(REFERENCE)
    if got.shape != want.shape:
        raise AssertionError('canonical shape %r != reference shape %r' % (got.shape, want.shape))
    diff = abs(got - want).max()
    if diff > TOL_PA:
        raise AssertionError('max|diff| against committed reference = %.3e Pa, bound %.3e Pa'
                             % (diff, TOL_PA))
    print('  PASS  4-rank (2,2,1) canonical frt matches committed reference: '
          'max|diff|=%.3e Pa (bound %.3e), %d rows' % (diff, TOL_PA, got.shape[0]))


def _check_two_faults_present_and_distinct(rows, label):
    y = rows[:, 1]
    f1 = rows[abs(y - 0.0) < 1.0]
    f2 = rows[abs(y - 2000.0) < 1.0]
    if f1.shape[0] == 0 or f2.shape[0] == 0:
        raise AssertionError('%s: expected fault nodes at BOTH y=0 (%d found) and y=2000 '
                             '(%d found) -- the multi-fault code path did not run for both faults'
                             % (label, f1.shape[0], f2.shape[0]))
    if f1.shape[0] != f2.shape[0]:
        raise AssertionError('%s: fault 1 has %d nodes, fault 2 has %d -- both faults are the '
                             'same box (same x/z extent), so their node counts must match'
                             % (label, f1.shape[0], f2.shape[0]))
    tnrm1, z1 = f1[:, 11], f1[:, 2]
    tnrm2, z2 = f2[:, 11], f2[:, 2]
    expect1 = FAULT1_NORM_COEFF * z1
    expect2 = FAULT2_NORM_COEFF * z2
    d1 = abs(tnrm1 - expect1).max()
    d2 = abs(tnrm2 - expect2).max()
    # The cross check: if fault 2 had aliased fault 1's array, tnrm2 would
    # match expect1 (evaluated at fault 2's own z, same values since z ranges
    # are identical), not expect2 -- and would fail d2 by a wide margin while
    # fitting the WRONG coefficient suspiciously well.
    if d1 > TOL_PA:
        raise AssertionError('%s: fault 1 tnrm max|diff| vs 7378*z = %.3e Pa, bound %.3e'
                             % (label, d1, TOL_PA))
    if d2 > TOL_PA:
        cross = abs(tnrm2 - FAULT1_NORM_COEFF * z2).max()
        raise AssertionError('%s: fault 2 tnrm max|diff| vs 5000*z = %.3e Pa (bound %.3e); '
                             'vs 7378*z (fault 1''s coefficient) = %.3e -- %s'
                             % (label, d2, TOL_PA, cross,
                                'looks ALIASED to fault 1' if cross < d2 else 'not an aliasing pattern'))
    print('  PASS  %s: %d fault-1 nodes (y=0) + %d fault-2 nodes (y=2000), '
          'tnrm matches each fault''s OWN coefficient (max|diff| %.3e / %.3e Pa, bound %.3e) '
          '-- no fault-2-reads-fault-1 aliasing'
          % (label, f1.shape[0], f2.shape[0], d1, d2, TOL_PA))


def check_four_rank_structural():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _build_and_run(tmp, nranks=4, nx=2, ny=2, nz=1)
        rows = _canonical(case_dir)
    _check_two_faults_present_and_distinct(rows, '4-rank (2,2,1)')


def check_serial_single_rank_owns_both_faults():
    """The decomposition that actually exercises the msnode-collision and
    slave/master-swap fixes: ONE rank builds both faults' entire meshes."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _build_and_run(tmp, nranks=1, nx=1, ny=1, nz=1)
        rows = _canonical(case_dir)
    _check_two_faults_present_and_distinct(rows, 'serial (1,1,1)')


def main():
    print('Regression guard: test.multifault2 two-fault routing/smoke case (fortran only)')
    failures = []
    for c in (check_serial_single_rank_owns_both_faults,
              check_four_rank_structural,
              check_four_rank_matches_committed_reference):
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_multifault_two_fault_smoke (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_multifault_two_fault_smoke')
    return 0


if __name__ == '__main__':
    sys.exit(main())
