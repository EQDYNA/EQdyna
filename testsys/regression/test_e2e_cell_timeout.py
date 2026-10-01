#! /usr/bin/env python3
"""
Regression guard: run_e2e.py's per-cell timeout (the 11h45m hang).

THE INCIDENT. `test.tpv1053d x fortran` stalled on an MPI transport issue
for 11h45m with nobody noticing until a human manually checked and killed
it. run_e2e.py's old module docstring said outright "No timeouts: this is a
shared box and dynamic-rupture runs are slow when it is busy" -- true of the
box, but it meant a genuinely stuck cell and a merely slow one were
indistinguishable from the outside: nothing would ever time either one out.

THE FIX pins three things, each checked below with the REAL functions (never
reimplemented, rule 10a -- behaviour, not source text):

  1. run_e2e.cell_deadline_s: 3x the latest docs/perf_ledger.jsonl wall_s for
     a (case, backend) that HAS been measured; run_e2e.CELL_TIMEOUT_FALLBACK_S
     for one that has not (checked against the real, committed ledger file,
     not a fabricated dict, so a never-measured case genuinely gets the
     fallback this suite ships with).
  2. run_e2e.run_one_with_deadline: a cell that outlives its deadline is
     killed -- its real subprocess, by a PID verified via
     `ps -o pid,args -p <pid>` immediately before the signal, never a second,
     separately-computed PID list -- reported FAIL(timeout) (distinct from an
     ordinary FAIL via run_e2e._verdict_label), and the sweep is NOT left
     waiting on it.
  3. BEFORE/AFTER, in one process: calling the cell's blocking subprocess
     step directly (no deadline wrapper -- the OLD code's only option) blocks
     for the full duration of a slow child with no way to cut it short; the
     SAME child, run through run_one_with_deadline with a 1s deadline, is
     killed and reported in roughly 1s regardless of how long it would have
     run otherwise. That contrast IS the regression this guard pins: remove
     the wrapper (or make it merely sleep-and-hope rather than actually
     signal the child) and `check_hung_cell_is_killed...` below goes RED
     because the sleep child is still alive at the PID check.

Cheap (rule 9): the only subprocess spawned is a plain `sleep`, never a real
Fortran/jax solver -- this is a guard on the HARNESS's timeout mechanism, not
a physics cell. Under a few seconds total. Exits non-zero on any failure.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

import run_e2e  # noqa: E402


# --------------------------------------------------------------------------
# 1. deadline computation
# --------------------------------------------------------------------------
def check_deadline_is_3x_measured_wall_clock():
    ledger_costs = {('test.fake', 'fortran'): 123.4}
    got = run_e2e.cell_deadline_s('test.fake', 'fortran', ledger_costs)
    want = run_e2e.CELL_TIMEOUT_MULTIPLIER * 123.4
    assert abs(got - want) < 1e-9, (
        'cell_deadline_s(measured=123.4) = %r, expected exactly %gx it (%r)'
        % (got, run_e2e.CELL_TIMEOUT_MULTIPLIER, want))


def check_deadline_falls_back_for_unmeasured_fabricated_entry():
    got = run_e2e.cell_deadline_s('test.fake', 'fortran', {})
    assert got == run_e2e.CELL_TIMEOUT_FALLBACK_S, (
        'cell_deadline_s with no ledger entry returned %r, expected the '
        'documented fallback CELL_TIMEOUT_FALLBACK_S=%r'
        % (got, run_e2e.CELL_TIMEOUT_FALLBACK_S))


def check_deadline_falls_back_against_the_real_committed_ledger():
    """The fallback path exercised against the ACTUAL docs/perf_ledger.jsonl
    this repo ships, not a fabricated dict -- a case name guaranteed absent
    from it must still get the documented fallback, end to end through
    load_ledger_wall_costs()."""
    ledger_costs = run_e2e.load_ledger_wall_costs()
    fake_cell = ('test.__never_run_by_any_sweep__', 'fortran')
    assert fake_cell not in ledger_costs, (
        'the fabricated cell name %r unexpectedly collides with a real '
        'docs/perf_ledger.jsonl row -- pick a different placeholder' % (fake_cell,))
    got = run_e2e.cell_deadline_s(fake_cell[0], fake_cell[1], ledger_costs)
    assert got == run_e2e.CELL_TIMEOUT_FALLBACK_S, (
        'an unmeasured cell against the real ledger got deadline %r, '
        'expected CELL_TIMEOUT_FALLBACK_S=%r' % (got, run_e2e.CELL_TIMEOUT_FALLBACK_S))


def check_deadline_multiplier_against_a_real_measured_row():
    """And the measured branch, against whatever the real ledger's FIRST row
    actually says -- proves the two paths (fabricated dict, real file) agree
    on the same function, not two different code paths."""
    ledger_costs = run_e2e.load_ledger_wall_costs()
    assert ledger_costs, 'docs/perf_ledger.jsonl produced zero rows -- cannot check the measured branch'
    (case, backend), wall = next(iter(ledger_costs.items()))
    got = run_e2e.cell_deadline_s(case, backend, ledger_costs)
    want = run_e2e.CELL_TIMEOUT_MULTIPLIER * wall
    assert abs(got - want) < 1e-6, (
        'cell_deadline_s(%r, %r) = %r, expected %gx its real ledger wall_s %r = %r'
        % (case, backend, got, run_e2e.CELL_TIMEOUT_MULTIPLIER, wall, want))


# --------------------------------------------------------------------------
# 2. verdict labelling -- FAIL(timeout) stays distinct from an ordinary FAIL
# --------------------------------------------------------------------------
def check_ordinary_failure_is_not_mislabelled_timeout():
    label = run_e2e._verdict_label(False, ['RuntimeError: solver exited 1'])
    assert label == 'FAIL', 'an ordinary failure was labelled %r, not FAIL' % label


def check_success_is_labelled_success():
    label = run_e2e._verdict_label(True, [])
    assert label == 'SUCCESS', 'a passing cell was labelled %r, not SUCCESS' % label


def check_timeout_marker_is_labelled_distinctly():
    label = run_e2e._verdict_label(False, ['FAIL(timeout): cell exceeded its 1.0s deadline -- killed pid 1'])
    assert label == 'FAIL(timeout)', 'a timed-out cell was labelled %r, not FAIL(timeout)' % label


# --------------------------------------------------------------------------
# 3. BEFORE: a direct, unwrapped call blocks for the full duration -- this is
#    exactly the old code's only option, and is what hung for 11h45m when the
#    child never exited on its own.
# --------------------------------------------------------------------------
def check_before_fix_an_unwrapped_call_blocks_for_the_full_duration():
    tmp = tempfile.gettempdir()
    env = dict(os.environ)
    t0 = time.time()
    rc = run_e2e._run(['sleep', '2'], tmp, env)
    dt = time.time() - t0
    assert rc == 0, 'sleep 2 exited %r, expected 0' % rc
    assert dt >= 1.8, (
        'a direct, unwrapped _run([\'sleep\', \'2\'], ...) returned after only '
        '%.2fs -- expected it to block for close to the full 2s, which is '
        'exactly the pre-fix behaviour this guard contrasts against '
        'check_hung_cell_is_killed_reported_and_actually_gone below' % dt)


# --------------------------------------------------------------------------
# 4. AFTER: the SAME shape of blocking child, run through
#    run_one_with_deadline with a deliberately tiny deadline, is killed,
#    reported FAIL(timeout), and is CONFIRMED gone via its own PID -- not
#    merely trusted from the returned message.
# --------------------------------------------------------------------------
def check_hung_cell_is_killed_reported_and_actually_gone():
    tmp = tempfile.gettempdir()
    env = dict(os.environ)

    def fn(cb):
        # Exactly run_e2e._run's own call shape -- the real production path
        # every cell's subprocess steps go through, not a reimplementation.
        t0 = time.time()
        rc = run_e2e._run(['sleep', '5'], tmp, env)
        return (cb[0], cb[1], rc == 0, time.time() - t0, ['sleep exited %d' % rc])

    cb = ('test.__fake_hung_cell__', 'fake-backend')
    t0 = time.time()
    result = run_e2e.run_one_with_deadline(fn, cb, 1.0)
    elapsed = time.time() - t0

    assert elapsed < 4.0, (
        'run_one_with_deadline(deadline=1.0s) took %.2fs to return -- it must '
        'report FAIL(timeout) promptly and not wait on the killed child '
        '(this is what keeps the rest of the sweep moving)' % elapsed)
    case, backend, ok, dt, lines = result
    assert (case, backend) == cb
    assert ok is False, 'a cell that outlived its deadline reported ok=%r' % (ok,)
    assert lines and lines[0].startswith('FAIL(timeout)'), (
        'timed-out cell did not report FAIL(timeout): lines=%r' % (lines,))
    assert run_e2e._verdict_label(ok, lines) == 'FAIL(timeout)'

    # The kill message names the PID it killed ("killed pid <N> (verified
    # before kill: ...)"); parse it back out and check BY PID, independent of
    # the harness's own report, that the process is actually gone.
    msg = lines[0]
    assert 'killed pid ' in msg, (
        'FAIL(timeout) message does not name a killed pid -- cannot '
        'independently verify termination: %r' % (msg,))
    pid_str = msg.split('killed pid ', 1)[1].split()[0]
    pid = int(pid_str)
    check = subprocess.run(['ps', '-o', 'pid', '-p', str(pid)],
                           capture_output=True, text=True)
    assert str(pid) not in check.stdout, (
        'pid %d (%r) is STILL ALIVE after the reported kill -- checked '
        'independently via `ps -o pid -p %d`, not the harness\'s own claim: %r'
        % (pid, msg, pid, check.stdout))

    # The registry must not leak this thread's entry once the cell is done.
    with run_e2e._ACTIVE_SUBPROCS_LOCK:
        leaked = [p for p in run_e2e._ACTIVE_SUBPROCS.values() if p.pid == pid]
    assert not leaked, 'the killed subprocess is still registered in _ACTIVE_SUBPROCS'


# --------------------------------------------------------------------------
# 5. MPI-SHAPED process tree (blocker 2 audit finding): the direct child is
#    not the only process a cell leaves behind. A single bare `sleep` child
#    (checks above) cannot distinguish a kill path that reaches the whole
#    tree from one that only reaches the direct child -- a mutation breaking
#    process-GROUP killing but keeping single-PID killing would still pass
#    those checks. This exercises `mpirun -np 2 sleep N` (if mpirun exists),
#    whose ranks get their OWN process group separate from mpirun's --
#    MEASURED on this box (ps -eo pid,ppid,pgid,sid,args on a live
#    `mpirun -np 2 sleep 15`): mpirun pgid == mpirun pid; each rank's pgid ==
#    that rank's OWN pid, not mpirun's. os.killpg on mpirun's pgid alone
#    reaches only mpirun; both ranks survive, reparent to pid 1, and keep
#    running -- exactly reproducing the orphaned-rank shape of the original
#    11h45m hang's cleanup problem. Falls back to a Python
#    double-fork/detach (a different parent, a different process group, no
#    mpirun needed) when mpirun is not on PATH, so this check still runs
#    somewhere mpirun is unavailable.
# --------------------------------------------------------------------------
def _mpi_tree_fn(cb, tmp, env):
    t0 = time.time()
    rc = run_e2e._run(['mpirun', '-np', '2', 'sleep', '20'], tmp, env)
    return (cb[0], cb[1], rc == 0, time.time() - t0, ['mpirun exited %d' % rc])


def _detached_grandchild_fn(cb, tmp, env):
    """No mpirun available: a Python parent that double-forks a detached
    grandchild (its own session, its own process group -- same orphan shape
    an MPI rank has relative to mpirun) running `sleep 20`, then itself
    sleeps so run_e2e._run has a live child to register and kill."""
    script = (
        'import os, subprocess, time, sys\n'
        'pid = os.fork()\n'
        'if pid == 0:\n'
        '    os.setsid()\n'
        '    pid2 = os.fork()\n'
        '    if pid2 == 0:\n'
        '        os.execvp("sleep", ["sleep", "20"])\n'
        '    os._exit(0)\n'
        'os.waitpid(pid, 0)\n'
        'time.sleep(20)\n'
    )
    t0 = time.time()
    rc = run_e2e._run([sys.executable, '-c', script], tmp, env)
    return (cb[0], cb[1], rc == 0, time.time() - t0, ['parent exited %d' % rc])


def check_mpi_shaped_tree_is_fully_killed_not_just_direct_child():
    tmp = tempfile.gettempdir()
    env = dict(os.environ)
    have_mpirun = shutil.which('mpirun') is not None
    fn = ((lambda cb: _mpi_tree_fn(cb, tmp, env)) if have_mpirun else
          (lambda cb: _detached_grandchild_fn(cb, tmp, env)))
    label = 'mpirun -np 2 sleep 20' if have_mpirun else 'double-forked detached grandchild'

    # Grab the registered subprocess's pid and discover ITS direct children
    # (the ranks, or the double-fork's immediate child) from a second thread
    # while the deadline wrapper is in flight, so we can independently
    # confirm by PID after the kill -- not trust run_one_with_deadline's own
    # message, which (unlike check 4 above) does not name every descendant.
    discovered = {}

    def discover():
        # Poll briefly for the subprocess to be registered and (for the
        # mpirun case) for the ranks to actually have been forked.
        deadline = time.time() + 5.0
        while time.time() < deadline:
            with run_e2e._ACTIVE_SUBPROCS_LOCK:
                procs = list(run_e2e._ACTIVE_SUBPROCS.values())
            if procs:
                p = procs[-1]
                time.sleep(0.3)  # let mpirun/fork finish spawning its child(ren)
                children = run_e2e._direct_children(p.pid)
                if children:
                    discovered['pid'] = p.pid
                    discovered['children'] = children
                    return
            time.sleep(0.1)

    discoverer = threading.Thread(target=discover, daemon=True)
    discoverer.start()

    t0 = time.time()
    result = run_e2e.run_one_with_deadline(fn, ('test.__fake_mpi_cell__', 'fake-backend'), 1.0)
    elapsed = time.time() - t0
    discoverer.join(timeout=5.0)

    assert elapsed < 4.0, (
        '%s: run_one_with_deadline(deadline=1.0s) took %.2fs to return' % (label, elapsed))
    _case, _backend, ok, _dt, lines = result
    assert ok is False and lines and lines[0].startswith('FAIL(timeout)'), (
        '%s: expected FAIL(timeout), got ok=%r lines=%r' % (label, ok, lines))

    assert 'pid' in discovered, (
        '%s: never managed to discover the registered subprocess\'s direct '
        'children before/while it was killed -- test is inconclusive, not a '
        'pass' % label)
    all_pids = [discovered['pid']] + [cpid for cpid, _ in discovered['children']]
    assert len(discovered['children']) >= 1, (
        '%s: discovered zero direct children of pid %d -- the tree this '
        'check needs (ranks, or the detached grandchild) never showed up'
        % (label, discovered['pid']))

    # Give the kill path (already invoked inside run_one_with_deadline,
    # above) its own grace window, then verify EVERY discovered pid --
    # parent and every child -- independently via ps, by PID, not by
    # trusting the harness's own report.
    deadline = time.time() + 10.0
    still_alive = list(all_pids)
    while time.time() < deadline and still_alive:
        still_alive = [pid for pid in still_alive
                       if str(pid) in subprocess.run(
                           ['ps', '-o', 'pid', '-p', str(pid)],
                           capture_output=True, text=True).stdout]
        if still_alive:
            time.sleep(0.2)
    assert not still_alive, (
        '%s: pid(s) %r STILL ALIVE after the timeout kill -- only the direct '
        'child (mpirun itself / the fork parent) was reached, not its '
        'spawned descendant(s) -- this is exactly the orphaned-process-tree '
        'regression this check exists to catch' % (label, still_alive))


def check_cell_that_finishes_within_deadline_is_unaffected():
    """A cell that finishes comfortably inside its deadline must get its OWN
    result back untouched -- the wrapper must not fire, or mislabel, a cell
    that never hung."""
    def fn(cb):
        return (cb[0], cb[1], True, 0.01, ['max|diff|=0 bound=1e-8'])

    cb = ('test.__fake_fast_cell__', 'fake-backend')
    result = run_e2e.run_one_with_deadline(fn, cb, 5.0)
    assert result == (cb[0], cb[1], True, 0.01, ['max|diff|=0 bound=1e-8']), (
        'a cell well inside its deadline was altered by the wrapper: %r' % (result,))
    assert run_e2e._verdict_label(result[2], result[4]) == 'SUCCESS'


CHECKS = (
    check_deadline_is_3x_measured_wall_clock,
    check_deadline_falls_back_for_unmeasured_fabricated_entry,
    check_deadline_falls_back_against_the_real_committed_ledger,
    check_deadline_multiplier_against_a_real_measured_row,
    check_ordinary_failure_is_not_mislabelled_timeout,
    check_success_is_labelled_success,
    check_timeout_marker_is_labelled_distinctly,
    check_before_fix_an_unwrapped_call_blocks_for_the_full_duration,
    check_hung_cell_is_killed_reported_and_actually_gone,
    check_mpi_shaped_tree_is_fully_killed_not_just_direct_child,
    check_cell_that_finishes_within_deadline_is_unaffected,
)


def main():
    print('Regression guard: run_e2e.py per-cell timeout (11h45m hang)')
    failures = []
    for c in CHECKS:
        try:
            c()
            print('  PASS  %s' % c.__name__)
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_e2e_cell_timeout (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_e2e_cell_timeout')
    return 0


if __name__ == '__main__':
    sys.exit(main())
