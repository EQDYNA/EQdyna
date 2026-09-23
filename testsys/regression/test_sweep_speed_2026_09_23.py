#! /usr/bin/env python3
"""
Regression guard: the 2026-09-23 sweep-speed changes (rules 2, 10).

Two independent things landed together and each gets its own check here so
one regressing does not hide behind the other going green:

  1. RELEASE_ONLY (matrix.py). test.tpv36 x python-numpy and test.tpv37 x
     python-numpy are SUPPORTED and PASSING, held out of the EVERYDAY sweep
     (run.py e2e / run_e2e.py's default selection at --term gate) for
     wall-clock cost only, and restored at --term full (run.py release,
     rule 24). The defect shape this guards against: a third cell added to
     RELEASE_ONLY with no matching update here goes RED (a silent widening
     of what "everyday" no longer covers must not pass quietly), and a
     RELEASE_ONLY cell that is also declared UNSUPPORTED goes RED (release-
     only means "supported, cost-deferred", never "does not work").

  2. Longest-first scheduling (run_e2e.py). The concurrent cell runner must
     start its most expensive cells first, using docs/perf_ledger.jsonl as
     the measured cost source, with unmeasured cells scheduled first (the
     conservative choice, never silently defaulted to "cheap"). Pinned via
     run_e2e.schedule_order directly -- no sweep is run.

Cheap (rule 9): imports + plain-data assertions + one pure-function call
against a synthetic cost dict. No solver, no MPI, no I/O beyond reading the
two modules under test. Exits non-zero on any failure.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

from testsys import matrix  # noqa: E402
import run_e2e  # noqa: E402


# --------------------------------------------------------------------------
# item 1: RELEASE_ONLY
# --------------------------------------------------------------------------
EXPECTED_RELEASE_ONLY = {
    ('test.tpv36', 'python-numpy'),
    ('test.tpv37', 'python-numpy'),
}


def check_release_only_is_exactly_the_two_flagged_cells():
    """Pins the CONTENT of the flag: a third cell added here with no update
    to EXPECTED_RELEASE_ONLY above goes RED, and so does one of the two
    being silently dropped."""
    actual = set(matrix.RELEASE_ONLY)
    assert actual == EXPECTED_RELEASE_ONLY, (
        'matrix.RELEASE_ONLY is %r, expected exactly %r -- if a cell was '
        'deliberately added or removed from the release-only flag, update '
        'EXPECTED_RELEASE_ONLY in this guard as part of that same change'
        % (actual, EXPECTED_RELEASE_ONLY))
    print('  PASS  matrix.RELEASE_ONLY == %r' % sorted(EXPECTED_RELEASE_ONLY))


def check_release_only_cells_are_never_declared_unsupported():
    """release-only means "supported, cost-deferred to the release tier" --
    never "does not work". matrix.py's own import-time consistency check
    already refuses this combination; this restates the contract at the
    level a reader of this guard will look for it."""
    overlap = set(matrix.RELEASE_ONLY) & set(matrix.UNSUPPORTED)
    assert overlap == set(), (
        '%r is both RELEASE_ONLY and UNSUPPORTED -- a cell cannot honestly '
        'claim both' % overlap)
    print('  PASS  no RELEASE_ONLY cell is declared UNSUPPORTED')


def check_release_selection_equals_everyday_selection_plus_release_only():
    """The two selections run_e2e.py's default (non-explicit) path chooses
    between: --term gate (everyday) must be exactly --term full (release)
    minus RELEASE_ONLY, and vice versa. This is the mutation-tested
    boundary: move a cell across it by hand (below) and confirm it flips."""
    release_runnable, release_unsupported = matrix.cells()
    everyday_runnable, everyday_unsupported, everyday_release_only = \
        matrix.everyday_cells()

    assert set(release_runnable) - set(everyday_runnable) == EXPECTED_RELEASE_ONLY, (
        'release selection minus everyday selection must be exactly the '
        'flagged cells')
    assert set(everyday_runnable) & EXPECTED_RELEASE_ONLY == set(), (
        'everyday selection must exclude every RELEASE_ONLY cell')
    assert set((c, b) for c, b, _ in everyday_release_only) == EXPECTED_RELEASE_ONLY
    # unsupported cells are identical on both sides -- RELEASE_ONLY only
    # ever moves a cell between "runnable" and "held back", never touches
    # the unsupported set.
    assert release_unsupported == everyday_unsupported
    print('  PASS  release == everyday + %d release-only cell(s), '
          'unsupported set unchanged' % len(EXPECTED_RELEASE_ONLY))


def check_mutation_moving_a_cell_across_the_boundary_flips_the_guard():
    """Take one real, currently-everyday cell, move it across the boundary
    by hand (a synthetic RELEASE_ONLY dict, not matrix.py's real one), and
    confirm the equality above flips red, then confirm the untouched real
    table still reports the original (unflipped) relationship -- proving
    this guard actually distinguishes the two states rather than being
    trivially true either way."""
    real_release_only = dict(matrix.RELEASE_ONLY)
    probe_cell = ('test.tpv8', 'python-jax')  # a real, everyday, non-flagged cell
    assert probe_cell not in real_release_only, (
        '%r is already RELEASE_ONLY -- pick a different probe cell for this '
        'mutation check' % (probe_cell,))

    mutated = dict(real_release_only)
    mutated[probe_cell] = 'synthetic, test-only'
    try:
        matrix.RELEASE_ONLY = mutated
        _, _, release_only_after = matrix.everyday_cells()
        moved = set((c, b) for c, b, _ in release_only_after)
        assert probe_cell in moved, (
            'moving %r into RELEASE_ONLY did not remove it from the '
            'everyday selection -- the boundary check above would not '
            'catch a real regression here' % (probe_cell,))
        assert moved != EXPECTED_RELEASE_ONLY, (
            'the mutated 3-cell release-only set was not distinguishable '
            'from the real 2-cell one')
    finally:
        matrix.RELEASE_ONLY = real_release_only

    # restored: the real table reports exactly the original relationship again
    _, _, release_only_restored = matrix.everyday_cells()
    assert set((c, b) for c, b, _ in release_only_restored) == EXPECTED_RELEASE_ONLY
    print('  PASS  moving %r across the boundary flips the selection and '
          'restoring the table restores it (mutation-tested both ways)'
          % (probe_cell,))


# --------------------------------------------------------------------------
# item 2: longest-first scheduling
# --------------------------------------------------------------------------
def check_schedule_order_is_longest_first():
    """Direct check of run_e2e.schedule_order: given measured costs, the
    returned order's cost sequence must be non-increasing."""
    costs = {
        ('a', 'fortran'): 10.0,
        ('b', 'python-numpy'): 1000.0,
        ('c', 'python-jax'): 100.0,
    }
    cells = [('a', 'fortran'), ('b', 'python-numpy'), ('c', 'python-jax')]
    ordered = run_e2e.schedule_order(cells, costs)
    assert ordered == [('b', 'python-numpy'), ('c', 'python-jax'), ('a', 'fortran')]
    seq = [costs[cb] for cb in ordered]
    assert all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1)), (
        'schedule_order returned a sequence that is not non-increasing '
        'in estimated cost: %r' % seq)
    print('  PASS  schedule_order sorts strictly longest-first: %r' % ordered)


def check_unmeasured_cells_are_scheduled_first():
    """The stated rule for a cell with no ledger entry: schedule it first
    (the conservative choice), never silently treated as cheap."""
    costs = {('known', 'fortran'): 5.0}
    cells = [('known', 'fortran'), ('unmeasured', 'python-jax')]
    ordered = run_e2e.schedule_order(cells, costs)
    assert ordered[0] == ('unmeasured', 'python-jax'), (
        'an unmeasured cell must be scheduled first, got order %r' % (ordered,))
    print('  PASS  unmeasured cell scheduled first: %r' % ordered)


def check_reversing_the_sort_direction_goes_red():
    """Mutation: a scheduler sorting ascending (shortest first) instead of
    descending must fail this guard's own non-increasing check -- this is
    what proves the check above is not vacuously true."""
    costs = {
        ('a', 'fortran'): 10.0,
        ('b', 'python-numpy'): 1000.0,
        ('c', 'python-jax'): 100.0,
    }
    cells = [('a', 'fortran'), ('b', 'python-numpy'), ('c', 'python-jax')]

    def mutated_schedule_order(cells, cost_estimates):
        # the exact bug this item exists to prevent: ascending, not
        # descending -- reverses run_e2e.schedule_order's sort direction.
        return sorted(cells, key=lambda cb: cost_estimates.get(cb, float('inf')))

    mutated_order = mutated_schedule_order(cells, costs)
    seq = [costs[cb] for cb in mutated_order]
    non_increasing = all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1))
    assert not non_increasing, (
        'the deliberately-reversed scheduler produced a non-increasing '
        'sequence by coincidence on this input -- this input no longer '
        'exercises the mutation; pick costs with a strict order')
    print('  PASS  reversing the sort direction is caught (would be RED '
          'under the real non-increasing check): %r' % mutated_order)


def check_run_e2e_actually_uses_schedule_order_for_submission():
    """Source-level tie to the real implementation (same shape as
    test_sweep_core_budget.py's own source check): the cell list handed to
    the concurrent runner must be the output of schedule_order, not the raw
    table order, or this guard's unit-level checks above would be proving
    something the sweep itself does not do."""
    path = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py')
    raw = open(path, errors='replace').read()
    src = '\n'.join(ln.split('#', 1)[0] for ln in raw.splitlines())
    assert 'cells = schedule_order(table_cells, ledger_costs)' in src, (
        'run_e2e.py no longer assigns the scheduled order to `cells` before '
        'the concurrent runner submits them -- the longest-first scheduling '
        'this guard exists for may have been removed or bypassed')
    print('  PASS  run_e2e.py submits cells in schedule_order\'s order')


def main():
    print('Regression guard: 2026-09-23 sweep-speed changes '
          '(RELEASE_ONLY + longest-first scheduling)')
    checks = [
        check_release_only_is_exactly_the_two_flagged_cells,
        check_release_only_cells_are_never_declared_unsupported,
        check_release_selection_equals_everyday_selection_plus_release_only,
        check_mutation_moving_a_cell_across_the_boundary_flips_the_guard,
        check_schedule_order_is_longest_first,
        check_unmeasured_cells_are_scheduled_first,
        check_reversing_the_sort_direction_goes_red,
        check_run_e2e_actually_uses_schedule_order_for_submission,
    ]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_sweep_speed_2026_09_23 (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_sweep_speed_2026_09_23')
    return 0


if __name__ == '__main__':
    sys.exit(main())
