#! /usr/bin/env python3
"""
Regression guard: profile_record.capture_run must always be handed the
REAL rank count of a cell's launch, never cell_cost()'s scheduling-budget
CORE cost.

THE INCIDENT (found by a paired gate sweep, 2026-09-23; diagnosed by mira):
testsys/e2e/run_e2e.py's run_one passed `ranks=cell_cost(case, backend)` to
profile_record.capture_run. cell_cost('*', 'python-jax') bills
ceil(matrix.JAX_MEASURED_CORES) == 3 (a CPU-thread scheduling cost, added by
the 2026-09-23 speed campaign) for a launch that is one serial process
writing exactly one profile.rank0.json. capture_run asserts its caller's
`ranks` equals the profile files' own reported nranks and RAISES on a
mismatch -- so every python-jax cell that PASSED PHYSICS (confirmed by a
direct compare: test.tpv8 python-jax ok=True, max|diff|=1.983643e-10) was
reported FAIL anyway, purely from this bookkeeping conflation.

FIX (testsys/e2e/run_e2e.py): a new module-level `profile_ranks(case,
backend)`, read by BOTH the profile-collection call (via
`_run_profile_capture`) and `_perf_meta`'s ranks column -- never cell_cost.

WHAT THIS FILE CHECKS, mutation-tested both ways (see each check's own
docstring):
  1. profile_ranks() returns the REAL rank count for every (case, backend)
     cell in the table -- the count matching how many profile.rank<r>.json
     files that backend's launch actually writes (fortran/jax-mpi: the
     matrix rank table; numpy/jax: 1) -- and, on cells where the two numbers
     differ, is NOT equal to cell_cost().
  2. `_run_profile_capture` -- the actual function run_one calls -- hands
     profile_ranks()'s value to capture_run, not cell_cost()'s. This is a
     BEHAVIOURAL check (it monkeypatches capture_run and inspects what it
     was actually called with), so reverting the call site back to
     `cell_cost` -- the exact regression -- turns this RED without needing
     any source-text scan.
  3. A capture_run failure is reported as ITS OWN, distinctly-labelled line
     (never indistinguishable from a physics-divergence line from
     compare.compare_cell) and still fails the cell (ok=False) -- the
     owner's mechanical-guard requirement: a missing/invalid profile SHOULD
     fail a cell, but the message must say so, not read as a physics defect.

Cheap (rule 9): imports + monkeypatched pure-Python calls, no subprocess, no
solver launch. Exits non-zero on any failure.
"""
import importlib.util
import os
import sys

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
E2E_PY = os.path.join(TESTSYS, 'e2e', 'run_e2e.py')

if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
from testsys import matrix  # noqa: E402


def _load_module(name, path):
    """A fresh module object from `path` -- never the cached sys.modules
    entry, so each check gets an independent instance it can monkeypatch
    without leaking into any other test (same discipline test_term_axis.py
    uses for this exact file)."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_profile_ranks_matches_real_launch_for_every_cell():
    """For every (case, backend) cell the table actually runs, profile_ranks
    must equal the rank count that backend's launch really uses -- fortran
    and python-jax-mpi from the matrix rank tables (real MPI launches),
    python-numpy and python-jax always 1 (serial processes, regardless of
    internal thread count)."""
    mod = _load_module('run_e2e_profile_ranks_a', E2E_PY)
    runnable, _unsupported = matrix.cells()
    assert runnable, 'matrix.cells() returned an empty table -- cannot check anything'
    for case, backend in runnable:
        got = mod.profile_ranks(case, backend)
        if backend == 'fortran':
            want = matrix.FORTRAN_RANKS[case]
        elif backend == 'python-jax-mpi':
            want = matrix.PY_MPI_RANKS[case]
        else:
            want = 1
        assert got == want, (
            'profile_ranks(%r, %r) = %r, expected %r (the real rank count '
            'that backend launch uses -- i.e. how many profile.rank<r>.json '
            'files it writes)' % (case, backend, got, want))


def check_profile_ranks_disagrees_with_cell_cost_on_python_jax():
    """The whole bug was that cell_cost('*', 'python-jax') != 1 while the
    real launch is serial (1 rank). Pin that fact down directly: if
    JAX_MEASURED_CORES ever drops so its ceil is 1, cell_cost and
    profile_ranks would coincidentally agree and this incident could recur
    silently -- so assert the precondition explicitly rather than assume
    it, per this file's own rule against a check that passes without
    exercising anything."""
    mod = _load_module('run_e2e_profile_ranks_b', E2E_PY)
    case = 'test.tpv8'
    backend = 'python-jax'
    real = mod.profile_ranks(case, backend)
    cost = mod.cell_cost(case, backend)
    assert real == 1, 'profile_ranks(%r, %r) = %r, expected 1 (a serial process)' % (case, backend, real)
    assert cost != real, (
        'cell_cost(%r, %r) = %r, same as profile_ranks() = %r -- this check '
        'assumes matrix.JAX_MEASURED_CORES (%r) ceils above 1 so the two '
        'numbers are DISTINGUISHABLE; if that constant changed, update this '
        'assumption rather than let the check pass without exercising '
        'anything' % (case, backend, cost, real, matrix.JAX_MEASURED_CORES))


def check_capture_uses_profile_ranks_not_cell_cost():
    """Calls the REAL `_run_profile_capture` (the function run_one actually
    invokes), with profile_record.capture_run monkeypatched to just record
    what `ranks` it was called with. Asserts that value is profile_ranks(),
    not cell_cost().

    MUTATION-SENSITIVE: if a future edit reverts the call inside
    `_run_profile_capture` from `profile_ranks(case, backend)` back to
    `cell_cost(case, backend)` -- the exact regression this file guards --
    the recorded `ranks` becomes 3 instead of 1 and this assertion goes RED.
    (Verified by hand this session: with that one-line revert applied
    in-session and immediately reverted back, this check FAILED as expected,
    `seen['ranks']=3 != profile_ranks()=1`; restored, it PASSES again.)"""
    mod = _load_module('run_e2e_profile_ranks_c', E2E_PY)
    case, backend = 'test.tpv8', 'python-jax'
    seen = {}

    def fake_capture_run(run_dir, *, case, backend, ranks, term, sha, tree_dirty=None):
        seen['ranks'] = ranks
        return []

    mod.profile_record.capture_run = fake_capture_run
    ok, lines = mod._run_profile_capture('/nonexistent-not-touched', case,
                                         backend, '5.0', 'deadbeefcafe', False)
    expected = mod.profile_ranks(case, backend)
    wrong = mod.cell_cost(case, backend)
    assert ok is True and lines == [], (
        '_run_profile_capture reported failure on a capture_run stub that '
        'succeeded: ok=%r lines=%r' % (ok, lines))
    assert seen.get('ranks') == expected, (
        '_run_profile_capture handed ranks=%r to capture_run, expected '
        'profile_ranks(%r, %r)=%r (cell_cost() would have been %r)'
        % (seen.get('ranks'), case, backend, expected, wrong))


def check_profile_failure_is_labelled_and_fails_the_cell():
    """capture_run raising must (a) make _run_profile_capture report
    ok=False -- a profile defect really does fail the cell, the owner's
    mechanical-guard requirement -- and (b) produce a line that names the
    profile, not one that could be mistaken for a physics-comparison
    line."""
    mod = _load_module('run_e2e_profile_ranks_d', E2E_PY)

    def raising_capture_run(run_dir, *, case, backend, ranks, term, sha, tree_dirty=None):
        raise ValueError('missing profile.rank1.json')

    mod.profile_record.capture_run = raising_capture_run
    ok, lines = mod._run_profile_capture('/nonexistent-not-touched',
                                         'test.tpv8', 'python-jax', '5.0',
                                         'deadbeefcafe', False)
    assert ok is False, 'a capture_run failure must fail the cell (ok=False)'
    assert len(lines) == 1, 'expected exactly one line describing the profile failure, got %r' % (lines,)
    line = lines[0]
    assert 'PROFILE' in line, (
        'profile-failure line %r does not name the profile -- indistinguishable '
        'from a physics-divergence line' % line)
    assert 'missing profile.rank1.json' in line, (
        'profile-failure line %r lost the underlying exception message' % line)
    assert 'max|diff|' not in line and 'diff' not in line.lower(), (
        'profile-failure line %r reads like a physics-comparison line' % line)


def main():
    checks = [check_profile_ranks_matches_real_launch_for_every_cell,
              check_profile_ranks_disagrees_with_cell_cost_on_python_jax,
              check_capture_uses_profile_ranks_not_cell_cost,
              check_profile_failure_is_labelled_and_fails_the_cell]
    failures = []
    for c in checks:
        try:
            c()
            print('  PASS  %s' % c.__name__)
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_profile_ranks_helper (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_profile_ranks_helper (%d checks)' % len(checks))
    return 0


if __name__ == '__main__':
    sys.exit(main())
