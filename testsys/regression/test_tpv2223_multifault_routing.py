#! /usr/bin/env python3
"""
Regression guard: test.tpv22, Fortran, SERIAL decomposition (one rank owns
BOTH faults) -- the one multi-fault routing check test.multifault2's own
smoke test used to carry that the everyday e2e gate for test.tpv22/
test.tpv23 does NOT exercise (item 17 section A, retiring test.multifault2
in favour of the real SCEC benchmarks as the two-fault gate cases;
pathway_forward.md item 17).

WHAT test.multifault2's SMOKE TEST CHECKED (ported forward here, this part
only -- see the dispatch report for the full list and which parts the
gate already covers on its own):

  1. Fortran actually MESHES and RUNS two faults WITHOUT the bugs a real
     2-fault run found: msnode collision across faults (meshgen.f90
     createMasterNode), the slave/master-swap gate keyed to one fault's
     y-plane, and single-fault-sized fltgm/fltl../fltnum MPI exchange
     arrays. All three are fixed in src/fortran/ (see git history); this
     guard is what keeps them fixed.
  2. This is the multi-fault code path actually running for BOTH faults
     (not falling through to ntotft==1 behaviour): both fault node counts
     are positive and distinct in y.

NOT ported forward (covered elsewhere, see the dispatch report):
  - 4-rank-split-across-faults decomposition and byte-parity against a
    committed reference: covered by the e2e gate itself
    (test.tpv22/test.tpv23 x fortran, matrix.FORTRAN_RANKS=4 -- a different
    decomposition that happens to also split the two faults across ranks,
    since par.ny=1 keeps every fault's own y-plane off the x/z partition
    boundaries that DO split -- the full node-for-node byte comparison the
    e2e gate performs is a strictly STRONGER check than multifault2's own
    coefficient-based anti-aliasing probe).
  - per-fault station file existence/routing (faultstft2_*-tagged files):
    covered by matrix.GATE_STATIONS' own 'faultstft2_050dp050.txt' entry,
    on BOTH gated backends.
  - anti-aliasing of fault #2's physics onto fault #1's: covered, far more
    strongly than multifault2's synthetic coefficient check, by this
    mission's own independent SCEC cross-validation (kaneko/SPECFEM3D,
    barall/FaultMod) at a fault-#2 station -- aliased physics could not
    independently match a different, externally-produced solution
    (testsys/parity/evidence_tpv22_23_scec_comparison.py).
  - the anti-ntotft==1 case.setup guard: test_multifault_refused.py, built
    on test.tpv8 + an appended par.ntotft=2 override, independent of any
    committed multifault2 reference -- left untouched.

THIS test is the one gap registering test.tpv22/test.tpv23 as the gate
cases left open: no gated cell ever runs Fortran with ONE rank owning both
faults (the e2e cell is always matrix.FORTRAN_RANKS[case]=4). That
decomposition is exactly what exercises createMasterNode's and
replaceSlaveWithMasterNode's cross-fault paths most directly (both faults'
nodes built and connected on the SAME rank, in the SAME call).

COST: builds test.tpv22 at its real, committed 200 m resolution (the mesh
size the bug classes above actually depend on -- a coarser synthetic mesh
would not exercise the same code paths at the same scale) but with
par.term cut to 0.5 s (~30 steps instead of ~900), serial. This is a
structural/routing check, not a physics parity gate -- it does not compare
against the 15 s frt.canonical.txt reference (a different term entirely);
it asserts the two fault node COUNTS, which are a pure function of the
mesh geometry and do not depend on how many steps ran.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import frt_canonical  # noqa: E402

MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE_NAME = 'test.tpv22'
FAULT1_Y = 0.0
FAULT2_Y = -1600.0  # test.tpv22's own stepover offset (tpv22_23_common.py)
# The known, full-200m-resolution per-fault node count -- a pure function of
# geometry (FAULT1_XRANGE/FAULT2_XRANGE/FZMIN/FZMAX/dx/dz), independent of
# par.term: this is the SAME 15251 each side of
# test.reference.results/test.tpv22/frt.canonical.txt's own 30502 rows
# (15251 fault-1 + 15251 fault-2), frozen by this same mission.
EXPECTED_NFTND_PER_FAULT = 15251
SHORT_TERM_S = 0.5  # ~30 steps at this case's dt -- routing/structure only.


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


def _build_and_run_serial(tmp):
    """create.newcase + force serial (npx=npy=npz=1) + a short par.term,
    leaving every other committed test.tpv22 parameter (dx=200m, geometry,
    friction, stress) untouched -- this is the REAL compset's own mesh, not
    a synthetic substitute."""
    case_dir = os.path.join(tmp, 'case_serial')
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial + shortened term by this test '
                       '(routing/structure check, not a physics parity gate)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n'
                       'par.term = %r\n' % SHORT_TERM_S)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))

    exe = _require_binary()
    rc, out = mpirun_capture.run_rank_files(MPIRUN, exe, case_dir, np=1, env=env,
                                            timeout=300)
    if rc != 0:
        raise AssertionError('eqdyna exited %d for serial test.tpv22 (both faults, '
                             'one rank):\n%s' % (rc, out[-3000:]))
    return case_dir


def check_serial_single_rank_builds_and_runs_both_faults():
    """ONE rank meshes, connects and steps BOTH test.tpv22 faults -- the
    decomposition that exercises createMasterNode's/replaceSlaveWithMaster
    Node's cross-fault paths most directly (test.multifault2's own smoke
    test's reason for running this decomposition at all)."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _build_and_run_serial(tmp)
        rows = frt_canonical.canonical_from_case(case_dir)
    y = rows[:, 1]
    f1 = rows[np.abs(y - FAULT1_Y) < 1.0]
    f2 = rows[np.abs(y - FAULT2_Y) < 1.0]
    if f1.shape[0] == 0 or f2.shape[0] == 0:
        raise AssertionError(
            'serial test.tpv22: expected fault nodes at BOTH y=%r (%d found) and '
            'y=%r (%d found) -- the multi-fault code path did not run for both '
            'faults on a single rank' % (FAULT1_Y, f1.shape[0], FAULT2_Y, f2.shape[0]))
    if f1.shape[0] != EXPECTED_NFTND_PER_FAULT or f2.shape[0] != EXPECTED_NFTND_PER_FAULT:
        raise AssertionError(
            'serial test.tpv22: expected %d fault nodes per fault (this case\'s '
            'known 200m-resolution geometry, independent of par.term), got '
            'fault1=%d fault2=%d -- a mesh-generation regression, not a physics '
            'one (par.term was shortened; node COUNT does not depend on it)'
            % (EXPECTED_NFTND_PER_FAULT, f1.shape[0], f2.shape[0]))
    print('  PASS  serial (1 rank owns both faults): fault1 %d nodes (y=%g), '
          'fault2 %d nodes (y=%g), both match the known 200m-resolution count'
          % (f1.shape[0], FAULT1_Y, f2.shape[0], FAULT2_Y))


def main():
    print('Regression guard: test.tpv22 Fortran serial (one rank, both faults) routing')
    failures = []
    for c in (check_serial_single_rank_builds_and_runs_both_faults,):
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_tpv2223_multifault_routing (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_tpv2223_multifault_routing')
    return 0


if __name__ == '__main__':
    sys.exit(main())
