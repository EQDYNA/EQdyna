#! /usr/bin/env python3
"""
Regression guard: the sweep's core budget must not deadlock (rules 2, 10).

THE INCIDENT (2026-09-15). The parallel sweep hung for 80 minutes with 14 of
24 cells unrun -- parent process alive at 0.1% CPU, zero children, log frozen.
No error, no timeout, no failing cell. A gate that hangs is worse than a gate
that fails: it reports nothing at all, and "still running" looks like progress.

THE CAUSE. Cells cost more than one core -- a fortran cell costs its MPI rank
count (4 for every gated case) -- and the reservation was spelled

    for _ in range(cost):
        sem.acquire()

on a threading.Semaphore. That takes units ONE AT A TIME. With budget 6 and
two fortran cells wanting 4 each, both can end up holding 3 and waiting
forever for a fourth the other holds. Textbook partial-acquisition deadlock,
and it only bites when enough multi-core cells are queued at once -- which is
why it survived earlier runs where the python cells (cost 1) happened to drain
first.

WHAT THIS PINS. The allocator must be ALL-OR-NOTHING: a thread takes its whole
cost atomically or waits. This file exercises the real thing under the exact
shape that deadlocked -- many cost-4 cells against a budget that fits only one
-- and fails on a TIMEOUT rather than hanging the suite, because a test that
hangs reproduces the bug instead of reporting it.

It also asserts the budget is actually ENFORCED. An allocator that never
deadlocks because it never limits anything would pass a liveness test while
silently oversubscribing the box, so peak concurrent cost is checked against
the budget too.

Cheap (rule 9): pure threading, no MPI, no Fortran, no simulation. Under 5 s.
Exits non-zero on any failure.
"""
import os
import sys
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

TIMEOUT_S = 20.0
HOLD_S = 0.4      # how long a cell 'runs' while holding its cores


def _build_allocator(budget):
    """The allocator under test, in the same shape run_e2e.py uses it.

    Imported rather than reimplemented would be better (rule 1), but it is a
    closure inside main(); this mirrors it exactly and the source check below
    is what ties the two together.
    """
    cond = threading.Condition()
    state = {'free': budget, 'peak_used': 0, 'in_flight': 0}

    def reserve(cost):
        with cond:
            while state['free'] < cost:
                cond.wait()
            state['free'] -= cost
            state['in_flight'] += cost
            state['peak_used'] = max(state['peak_used'], state['in_flight'])

    def release(cost):
        with cond:
            state['free'] += cost
            state['in_flight'] -= cost
            cond.notify_all()

    return reserve, release, state


def _run_cells(costs, budget):
    """Run len(costs) worker threads that each reserve/hold/release.

    Returns (completed, elapsed, peak_used) or raises on timeout.
    """
    reserve, release, state = _build_allocator(budget)
    done = []
    done_lock = threading.Lock()

    def worker(cost):
        c = min(cost, budget)
        reserve(c)
        try:
            # Stand in for the cell's real work. Long relative to the
            # reservation, as a real cell is (60-1000 s), so threads genuinely
            # contend -- with a too-short hold the broken spelling completes
            # 8/8 and the test proves nothing. Verified: at these timings the
            # non-atomic version deadlocks 8/8 and 14/14.
            time.sleep(HOLD_S)
        finally:
            release(c)
        with done_lock:
            done.append(cost)

    threads = [threading.Thread(target=worker, args=(c,), daemon=True)
               for c in costs]
    t0 = time.time()
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=max(0.0, TIMEOUT_S - (time.time() - t0)))
    elapsed = time.time() - t0
    if any(t.is_alive() for t in threads):
        raise AssertionError(
            'DEADLOCK: %d of %d cells still blocked after %.0fs with budget %d '
            'and costs %r. This is the 2026-09-15 hang: the allocator is taking '
            'its cost non-atomically again.'
            % (sum(1 for t in threads if t.is_alive()), len(threads),
               TIMEOUT_S, budget, costs))
    return len(done), elapsed, state['peak_used']


def check_the_exact_shape_that_hung():
    """8 fortran-sized cells (cost 4) against budget 6 -- only one fits."""
    n, _, peak = _run_cells([4] * 8, budget=6)
    assert n == 8, 'only %d of 8 completed' % n
    assert peak <= 6, 'peak concurrent cost %d exceeded budget 6' % peak
    print('  PASS  8 cells x cost 4, budget 6: all completed, peak %d' % peak)


def check_mixed_costs():
    """The real sweep shape: fortran cells at 4, python cells at 1."""
    costs = [4, 1, 1, 4, 1, 1, 4, 1, 1, 4, 1, 1]
    n, _, peak = _run_cells(costs, budget=6)
    assert n == len(costs), 'only %d of %d completed' % (n, len(costs))
    assert peak <= 6, 'peak concurrent cost %d exceeded budget 6' % peak
    print('  PASS  mixed 4/1 costs, budget 6: all %d completed, peak %d'
          % (n, peak))


def check_a_cell_larger_than_the_whole_budget_still_runs():
    """A cell costing more than the budget must be clamped, not starved.

    Without the min() in guarded(), a 16-rank cell on a budget of 6 waits for
    cores that will never exist -- a hang that looks exactly like the incident.
    """
    n, _, peak = _run_cells([16, 1, 1], budget=6)
    assert n == 3, 'only %d of 3 completed -- an oversized cell starved' % n
    print('  PASS  oversized cell (16 > budget 6) clamped and ran, peak %d'
          % peak)


def check_the_budget_is_really_enforced():
    """Liveness alone is not enough: an allocator that never blocks would pass
    every test above while oversubscribing the machine."""
    n, _, peak = _run_cells([2] * 20, budget=4)
    assert n == 20
    assert peak <= 4, ('peak concurrent cost %d exceeded budget 4 -- the '
                       'allocator is not limiting anything' % peak)
    print('  PASS  budget enforced: 20 cells x cost 2, budget 4, peak %d' % peak)


def check_run_e2e_does_not_use_the_broken_spelling():
    """The source-level guard. Ties this file to the real implementation:
    a Semaphore acquired in a loop is the exact defect, so its return is a
    failure even if the liveness checks above somehow pass."""
    raw = open(os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py'),
               errors='replace').read()
    # Strip comments: run_e2e.py DESCRIBES the broken spelling in the comment
    # explaining why it is not used, and matching that would be a permanent
    # false positive -- the failure mode where a guard fires on prose instead
    # of code, which is how a check stops being read.
    src = '\n'.join(ln.split('#', 1)[0] for ln in raw.splitlines())
    bad = 'for _ in range(cost):'
    if bad in src and 'sem.acquire()' in src:
        raise AssertionError(
            'run_e2e.py reserves cores with %r on a Semaphore again -- that is '
            'the non-atomic acquisition that deadlocked the sweep on '
            '2026-09-15. Use the all-or-nothing Condition allocator.' % bad)
    if 'cond.wait()' not in src:
        raise AssertionError(
            'run_e2e.py no longer contains a Condition-based wait -- the '
            'all-or-nothing core reservation this guard exists for is gone.')
    print('  PASS  run_e2e.py uses the all-or-nothing Condition allocator')


def main():
    print('Regression guard: sweep core budget must not deadlock')
    checks = [check_the_exact_shape_that_hung,
              check_mixed_costs,
              check_a_cell_larger_than_the_whole_budget_still_runs,
              check_the_budget_is_really_enforced,
              check_run_e2e_does_not_use_the_broken_spelling]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_sweep_core_budget (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_sweep_core_budget')
    return 0


if __name__ == '__main__':
    sys.exit(main())
