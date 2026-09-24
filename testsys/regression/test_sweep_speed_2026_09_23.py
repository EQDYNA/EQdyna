#! /usr/bin/env python3
"""
Regression guard: the 2026-09-23 longest-first scheduling change (rules 2, 10).

This guard used to carry a SECOND item -- matrix.RELEASE_ONLY (test.tpv36 x
python-numpy and test.tpv37 x python-numpy, held out of the everyday sweep for
wall-clock cost) -- landed the same day as the scheduling change below. That
flag was RETIRED 2026-09-23 in the same owner decision that removed
python-numpy from the backend axis entirely (its only two occupants were both
python-numpy cells; a stale entry would now fail matrix.py's own import-time
consistency check). Its four checks are deleted with it, not left dead or
converted to check "the flag stays empty" -- matrix.everyday_cells()'s own
docstring already states that plainly, and a check whose only failure mode is
"someone re-added dead state" is not worth carrying. This file's name is left
as the 2026-09-23 changes it originally guarded; the surviving item is:

  Longest-first scheduling (run_e2e.py). The concurrent cell runner must
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

import run_e2e  # noqa: E402


# --------------------------------------------------------------------------
# longest-first scheduling
# --------------------------------------------------------------------------
def check_schedule_order_is_longest_first():
    """Direct check of run_e2e.schedule_order: given measured costs, the
    returned order's cost sequence must be non-increasing."""
    costs = {
        ('a', 'fortran'): 10.0,
        ('b', 'python-jax-mpi'): 1000.0,
        ('c', 'python-jax'): 100.0,
    }
    cells = [('a', 'fortran'), ('b', 'python-jax-mpi'), ('c', 'python-jax')]
    ordered = run_e2e.schedule_order(cells, costs)
    assert ordered == [('b', 'python-jax-mpi'), ('c', 'python-jax'), ('a', 'fortran')]
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
        ('b', 'python-jax-mpi'): 1000.0,
        ('c', 'python-jax'): 100.0,
    }
    cells = [('a', 'fortran'), ('b', 'python-jax-mpi'), ('c', 'python-jax')]

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
    print('Regression guard: 2026-09-23 longest-first scheduling change '
          '(RELEASE_ONLY checks retired with the flag itself, same date)')
    checks = [
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
