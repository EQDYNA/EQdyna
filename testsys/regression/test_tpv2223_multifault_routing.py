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

NOT ported forward (covered elsewhere, or not covered at all -- corrected
here, victor-reyes audit, PR #76 MAJOR 3, 2026-10-02: this paragraph used to
claim two things that are not true and is rewritten to say only what is
actually tested):
  - 4-rank-split-across-faults decomposition and byte-parity against a
    committed reference: test.tpv22/test.tpv23 x fortran at
    matrix.FORTRAN_RANKS=4 IS the e2e gate's own cell for these cases, but
    both cases are RELEASE_ONLY (testsys/matrix.py, cost reason, item 17
    section B) -- that 4-rank byte-parity comparison runs under `run.py
    release` only, NOT the everyday e2e gate and NOT CI. It is not a gap
    this file needs to cover (release coverage is real), but it is NOT an
    everyday/CI guarantee either; do not read the sentence above as saying
    it is.
  - per-fault station file existence/routing (faultstft2_*-tagged files):
    covered by matrix.GATE_STATIONS' own 'faultstft2_050dp050.txt' entry --
    also only at RELEASE_ONLY's cells, same caveat as above.
  - anti-aliasing of fault #2's physics onto fault #1's against an
    INDEPENDENT, externally-produced solution:
    testsys/parity/evidence_tpv22_23_scec_comparison.py does NOT cover this.
    That script is REPORT-ONLY by its own docstring and `main()` (always
    `return 0`; it never asserts), and the barall/FaultMod cross-code
    dataset it was originally written to compare against was removed from
    this repo's scec_archive in 6e76c93 ("remove 8 barall* cross-code
    dirs"). It is evidence a human can read, not a gate. The real,
    mechanical anti-aliasing guard for this case pair is
    check_two_faults_x_ranges_are_not_aliased below (new, this commit) --
    see its own docstring.
  - the anti-ntotft==1 case.setup guard: test_multifault_refused.py, built
    on test.tpv8 + an appended par.ntotft=2 override, independent of any
    committed multifault2 reference -- left untouched.

THIS test is a gap registering test.tpv22/test.tpv23 as the gate cases
left open in the EVERYDAY sweep and CI specifically (both cases being
RELEASE_ONLY means the serial, single-rank-owns-both-faults decomposition
below never runs there otherwise): no everyday/CI cell ever runs Fortran
with ONE rank owning both faults. That decomposition is exactly what
exercises createMasterNode's and replaceSlaveWithMasterNode's cross-fault
paths most directly (both faults' nodes built and connected on the SAME
rank, in the SAME call).

COST: builds test.tpv22 at its real, committed 200 m resolution (the mesh
size the bug classes above actually depend on -- a coarser synthetic mesh
would not exercise the same code paths at the same scale) but with
par.term cut to 0.5 s (~30 steps instead of ~900), serial. This is a
structural/routing check, not a physics parity gate -- it does not compare
against the 15 s frt.canonical.txt reference (a different term entirely);
it asserts the two fault node COUNTS (a pure function of mesh geometry,
term-independent) AND, since PR #76 MAJOR 3, each fault's own along-strike
x-RANGE, which is the real per-fault-distinguishing signature: TPV22/23's
two faults overlap in x only over [-5000, 5000] (FAULT1_XRANGE=(-25000,
5000), FAULT2_XRANGE=(-5000, 25000), tpv22_23_common.py) -- if fault 2 were
aliased onto fault 1's mesh/geometry (the exact bug class commit 1ac3d60
fixed, "restore per-fault mesh extent (ift), invent nothing"), its node x
values would be bounded by fault 1's box ([-25000, 5000]) and could never
reach fault 2's true far edge (25000 m). A bare node-COUNT match (the
previous version of this test) cannot catch that: both faults have the same
dx and the same along-strike LENGTH (30 km each), so an aliased fault 2
confined to fault 1's x-range would still report the same node count,
just at the wrong x positions.
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

# Copied verbatim from case_input/test.tpv22/tpv22_23_common.py (not imported
# -- that module has import-time side effects (building par) this test does
# not want) -- the two faults' TRUE along-strike x-extents, the real
# per-fault-distinguishing signature used by
# check_two_faults_x_ranges_are_not_aliased below (MAJOR 3, PR #76).
FAULT1_XRANGE = (-25000.0, 5000.0)
FAULT2_XRANGE = (-5000.0, 25000.0)
X_EDGE_TOL_M = 250.0  # a few dx (200m) of slack for the nearest-node snap


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


_CANONICAL = []  # one solver run shared by every check below (rule 9: cheap)


def _canonical_rows():
    if not _CANONICAL:
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = _build_and_run_serial(tmp)
            _CANONICAL.append(frt_canonical.canonical_from_case(case_dir))
    return _CANONICAL[0]


def _split_faults(rows):
    y = rows[:, 1]
    f1 = rows[np.abs(y - FAULT1_Y) < 1.0]
    f2 = rows[np.abs(y - FAULT2_Y) < 1.0]
    if f1.shape[0] == 0 or f2.shape[0] == 0:
        raise AssertionError(
            'serial test.tpv22: expected fault nodes at BOTH y=%r (%d found) and '
            'y=%r (%d found) -- the multi-fault code path did not run for both '
            'faults on a single rank' % (FAULT1_Y, f1.shape[0], FAULT2_Y, f2.shape[0]))
    return f1, f2


def check_serial_single_rank_builds_and_runs_both_faults():
    """ONE rank meshes, connects and steps BOTH test.tpv22 faults -- the
    decomposition that exercises createMasterNode's/replaceSlaveWithMaster
    Node's cross-fault paths most directly (test.multifault2's own smoke
    test's reason for running this decomposition at all)."""
    f1, f2 = _split_faults(_canonical_rows())
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


def check_two_faults_x_ranges_are_not_aliased():
    """A real per-fault DISTINGUISHING value check (MAJOR 3, PR #76),
    not just a node count: fault 1 and fault 2 have different, only
    partially overlapping along-strike x-extents (FAULT1_XRANGE=(-25000,
    5000) vs FAULT2_XRANGE=(-5000, 25000), tpv22_23_common.py) -- a value
    that WOULD be identical (both confined to fault 1's box) if fault 2's
    mesh were aliased onto fault 1's, the exact bug class commit 1ac3d60
    ("restore per-fault mesh extent (ift), invent nothing") fixed. Fault 1's
    own x values must never reach fault 2's exclusive far edge (25000 m)
    and fault 2's must never be confined to fault 1's far edge (5000 m) --
    asserts both faults' x range actually reaches its OWN true edge, within
    a couple of dx of nearest-node snap."""
    f1, f2 = _split_faults(_canonical_rows())
    x1, x2 = f1[:, 0], f2[:, 0]
    x1_min, x1_max = float(x1.min()), float(x1.max())
    x2_min, x2_max = float(x2.min()), float(x2.max())
    checks = [
        ('fault 1 min x', x1_min, FAULT1_XRANGE[0]),
        ('fault 1 max x', x1_max, FAULT1_XRANGE[1]),
        ('fault 2 min x', x2_min, FAULT2_XRANGE[0]),
        ('fault 2 max x', x2_max, FAULT2_XRANGE[1]),
    ]
    bad = [(label, got, want) for label, got, want in checks
           if abs(got - want) > X_EDGE_TOL_M]
    if bad:
        raise AssertionError(
            'serial test.tpv22: fault x-range does not reach its own true edge '
            '(aliased onto the other fault''s mesh box?): %s (tol %.0fm)'
            % ('; '.join('%s=%.1f (want %.1f)' % (l, g, w) for l, g, w in bad),
               X_EDGE_TOL_M))
    # The sharpest aliasing signature: fault 1 must NEVER reach fault 2's
    # exclusive far edge (x>5000 is only reachable by fault 2's true box),
    # and fault 2 must NEVER be confined inside fault 1's far edge (x<-5000
    # untouched by fault 2's true box starts at -5000, so this is the same
    # bound -- stated explicitly, not just via the edge-match above, since an
    # edge-match alone could in principle pass by coincidence on a bug that
    # only shifts the far edge).
    if x1_max > FAULT1_XRANGE[1] + X_EDGE_TOL_M:
        raise AssertionError('serial test.tpv22: fault 1 x reaches %.1f, past its own '
                             'true far edge %.1f -- looks aliased onto fault 2''s box'
                             % (x1_max, FAULT1_XRANGE[1]))
    if x2_min < FAULT2_XRANGE[0] - X_EDGE_TOL_M:
        raise AssertionError('serial test.tpv22: fault 2 x reaches %.1f, past its own '
                             'true near edge %.1f -- looks aliased onto fault 1''s box'
                             % (x2_min, FAULT2_XRANGE[0]))
    print('  PASS  serial (1 rank owns both faults): fault1 x in [%.0f, %.0f] '
          '(true box [%.0f, %.0f]), fault2 x in [%.0f, %.0f] (true box [%.0f, %.0f]) '
          '-- no fault-2-aliased-to-fault-1 x-range collapse'
          % (x1_min, x1_max, FAULT1_XRANGE[0], FAULT1_XRANGE[1],
             x2_min, x2_max, FAULT2_XRANGE[0], FAULT2_XRANGE[1]))


def main():
    print('Regression guard: test.tpv22 Fortran serial (one rank, both faults) routing')
    failures = []
    for c in (check_serial_single_rank_builds_and_runs_both_faults,
              check_two_faults_x_ranges_are_not_aliased):
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
