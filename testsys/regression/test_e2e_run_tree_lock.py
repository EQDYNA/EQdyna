#! /usr/bin/env python3
"""
Regression guard: a SECOND concurrent e2e invocation must REFUSE, not rotate
(PROJECT_RULES.md rules 21a, 8, 2, 10; pathway item 70).

THE INCIDENT THIS PINS, 2026-09-22 (docs/evidence/gate-v5.16.0/README.md).
`run_e2e.py` rotated `$REPO_ROOT/test` -> `test.prev/` at startup,
unconditionally and with no lock. wei-lin's v5.16.0 gate sweep started 21:48 in
the main checkout; a different session's e2e at 22:12:20 rotated the in-flight
sweep's live tree away. `test.tpv1053d x python-numpy` then ran 1500.1 s and
died with `FileNotFoundError .../test/test.tpv1053d.python-numpy/frt.txt0` --
`src/python/eqdyna/library_output.py:120` writes that file by absolute path at
the very END of a cell. It surfaced as `FAIL test.tpv1053d x python-numpy`,
INDISTINGUISHABLE in the summary from a parity regression. Isolated, the same
cell passes in 907 s.

WHAT IS ASSERTED, and what deliberately is not. This guard tests the LOCK, not
a sweep: nothing here builds, meshes or solves. The entry-point checks invoke
`run_e2e.py` / `run_e2e_full.py` for real, with the lock already held by this
process, and assert three things together -- non-zero exit, the refusal message
naming the holder's pid, and that `test/` was NOT renamed. The last is the one
that matters: an exit code alone would still pass if the tool refused AFTER
rotating.

NOTE ON A TRUE POSITIVE THAT LOOKS LIKE FLAKE. The entry-point checks take the
real lock on this checkout's own `test/`. If a genuine sweep is running in THIS
checkout right now, they fail -- correctly: that is the collision, and rule 21a
says the sweep belongs in its own worktree. The failure names that cause
explicitly so it is not mistaken for a broken guard.

Cheap (rule 9): flock, one short-lived child process, and two entry-point
invocations that refuse before their first gate. Under 5 s. Exits non-zero on
any failure.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from testsys import runlock  # noqa: E402

RUN_E2E = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py')
RUN_E2E_FULL = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e_full.py')


# --------------------------------------------------------------------------
# the mechanism, in a tempdir -- no repo state touched
# --------------------------------------------------------------------------
def check_a_second_acquire_refuses():
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, 'test'))
        first = runlock.acquire(tmp, 'test', argv=['first'], announce=False)
        try:
            runlock.acquire(tmp, 'test', argv=['second'], announce=False)
        except runlock.RunTreeLocked as exc:
            msg = str(exc)
        else:
            raise AssertionError(
                'a SECOND acquire on a held lock SUCCEEDED -- the lock does '
                'not lock, and two sweeps can rotate one tree again')
        finally:
            first.release()
        for want in (runlock.REFUSAL_HEADER, 'holder pid   : %d' % first.pid,
                     'holder since :', 'NOT waiting', 'NOT rotating anyway'):
            assert want in msg, ('refusal message omits %r; it read:\n%s'
                                 % (want, msg))
        assert re.search(r'holder since : \d{4}-\d\d-\d\d \d\d:\d\d:\d\d', msg), (
            'refusal names no holder START TIME in a readable form:\n%s' % msg)
        print('  PASS  second acquire refuses, naming holder pid %d and its '
              'start time' % first.pid)


def check_the_refusal_does_not_block():
    """It must refuse IMMEDIATELY. A lock that waits would park a multi-hour
    sweep behind another multi-hour sweep with nobody watching, and would hide
    the collision from whoever caused it."""
    with tempfile.TemporaryDirectory() as tmp:
        first = runlock.acquire(tmp, 'test', announce=False)
        t0 = time.time()
        try:
            runlock.acquire(tmp, 'test', announce=False)
            raise AssertionError('second acquire succeeded')
        except runlock.RunTreeLocked:
            dt = time.time() - t0
        finally:
            first.release()
        assert dt < 1.0, ('the refusal took %.1f s -- it is WAITING on the '
                          'lock instead of refusing' % dt)
        print('  PASS  refusal is non-blocking (%.3f s)' % dt)


def check_release_lets_the_next_invocation_in():
    """The counterpart to the refusal: a lock nobody can ever take again would
    pass every check above and wedge the checkout."""
    with tempfile.TemporaryDirectory() as tmp:
        runlock.acquire(tmp, 'test', announce=False).release()
        second = runlock.acquire(tmp, 'test', announce=False)
        try:
            assert second.takeover_notice is None, (
                'a cleanly-released lock was reported as a STALE takeover:\n%s'
                % second.takeover_notice)
        finally:
            second.release()
        print('  PASS  a released lock is re-acquirable, and silently (no '
              'stale-takeover notice)')


def check_a_dead_holder_is_taken_over_with_an_explicit_message():
    """SIGKILL the holder -- no atexit, no release, record left behind. The
    kernel drops the flock, so the takeover is automatic; the requirement is
    that it is never SILENT."""
    with tempfile.TemporaryDirectory() as tmp:
        code = (
            'import sys, time\n'
            'sys.path.insert(0, %r)\n'
            'from testsys import runlock\n'
            'lock = runlock.acquire(%r, "test", argv=["victim"], announce=False)\n'
            'print("HELD", lock.pid, flush=True)\n'
            'time.sleep(300)\n' % (ROOT, tmp))
        child = subprocess.Popen([sys.executable, '-c', code],
                                 stdout=subprocess.PIPE, text=True)
        try:
            line = child.stdout.readline().strip()
            assert line.startswith('HELD'), ('holder child never acquired: %r'
                                             % line)
            dead_pid = int(line.split()[1])
            child.kill()
            child.wait(timeout=10)
        finally:
            if child.poll() is None:
                child.kill()
        taken = runlock.acquire(tmp, 'test', announce=False)
        try:
            notice = taken.takeover_notice
            assert notice is not None, (
                'took over a lock whose holder was SIGKILLed mid-run and said '
                'NOTHING -- a silent takeover is exactly the fallback rule 2 '
                'forbids')
            assert runlock.TAKEOVER_HEADER in notice, notice
            assert 'pid %d' % dead_pid in notice, (
                'takeover notice does not name the dead holder pid %d:\n%s'
                % (dead_pid, notice))
            assert 'victim' in notice, (
                'takeover notice does not say what the dead holder was '
                'running:\n%s' % notice)
        finally:
            taken.release()
        print('  PASS  stale lock (holder pid %d SIGKILLed) taken over, with '
              'an explicit notice naming it' % dead_pid)


# --------------------------------------------------------------------------
# the entry points, for real, against this checkout's own trees
# --------------------------------------------------------------------------
def _refuses_without_rotating(tool, tool_argv, resource, extra_env=None,
                              banner=None):
    tree = os.path.join(ROOT, resource)
    prev = tree + '.prev'
    sentinel = 'item70_lock_guard_%d.marker' % os.getpid()
    try:
        held = runlock.acquire(ROOT, resource,
                               argv=['test_e2e_run_tree_lock.py', resource],
                               announce=False)
    except runlock.RunTreeLocked as exc:
        raise AssertionError(
            'could not take the lock on %s to run this check -- something '
            'else in THIS checkout holds it. That is the collision rule 21a '
            'forbids, not a broken guard: run your sweep in its own worktree. '
            'Holder:\n%s' % (tree, exc))
    created = not os.path.isdir(tree)
    try:
        os.makedirs(tree, exist_ok=True)
        with open(os.path.join(tree, sentinel), 'w') as fh:
            fh.write('pathway item 70 guard -- delete me if you find me\n')
        env = dict(os.environ)
        env.update(extra_env or {})
        proc = subprocess.run([sys.executable, tool] + tool_argv, cwd=ROOT,
                              env=env, capture_output=True, text=True,
                              timeout=300)
        out = proc.stdout + proc.stderr
        assert proc.returncode != 0, (
            '%s exited 0 while a second invocation held the lock on %s -- a '
            'refused run must never be green.\n%s'
            % (os.path.basename(tool), tree, out[-2000:]))
        assert runlock.REFUSAL_HEADER in out, (
            '%s exited %d but printed no refusal naming the lock holder; its '
            'output was:\n%s'
            % (os.path.basename(tool), proc.returncode, out[-2000:]))
        assert 'holder pid   : %d' % held.pid in out, (
            '%s refused without naming the holder pid %d:\n%s'
            % (os.path.basename(tool), held.pid, out[-2000:]))
        assert os.path.isfile(os.path.join(tree, sentinel)), (
            'THE DEFECT ITSELF: %s rotated %s away despite refusing -- the '
            'sentinel is gone from the live tree. Refusing AFTER rotating is '
            'the 2026-09-22 incident with an exit code attached.'
            % (os.path.basename(tool), tree))
        assert not os.path.isfile(os.path.join(prev, sentinel)), (
            'THE DEFECT ITSELF: %s moved the live %s to %s -- the sentinel '
            'turned up there.' % (os.path.basename(tool), tree, prev))
        print('  PASS  %s' % banner)
    finally:
        try:
            os.remove(os.path.join(tree, sentinel))
        except OSError:
            pass
        if created and os.path.isdir(tree) and not os.listdir(tree):
            shutil.rmtree(tree)
        held.release()


def check_run_e2e_refuses_and_does_not_rotate():
    _refuses_without_rotating(
        RUN_E2E, ['--cases', 'test.tpv8', '--backends', 'fortran'], 'test',
        banner='run_e2e.py refused a second invocation, exited non-zero, and '
               'left test/ in place')


def check_run_e2e_full_refuses_and_does_not_rotate():
    # EQDYNA_FULL_LAUNCH=yes-hours is required to get PAST run_e2e_full's own
    # opt-in refusal and reach the lock. It never reaches a solver: the lock is
    # taken before the build check, which is before the rotation, which is
    # before the first mpirun.
    _refuses_without_rotating(
        RUN_E2E_FULL, [], 'test.full',
        extra_env={'EQDYNA_FULL_LAUNCH': 'yes-hours'},
        banner='run_e2e_full.py refused a second invocation, exited non-zero, '
               'and left test.full/ in place')


# --------------------------------------------------------------------------
# shape of the sources (rule 2a: assert the shape, not a substring elsewhere)
# --------------------------------------------------------------------------
def _code_lines(path):
    raw = open(path, errors='replace').read()
    return [ln.split('#', 1)[0] for ln in raw.splitlines()]


def _first_index(lines, needle, path):
    for i, ln in enumerate(lines):
        if needle in ln:
            return i
    raise AssertionError('%s contains no %r' % (path, needle))


def check_the_lock_is_taken_before_the_rotation_in_both_tools():
    """Ordering, at the source level. A lock acquired AFTER `shutil.move` is
    not a lock -- the damage is already done by the time it is asked for."""
    for path, resource, mover in ((RUN_E2E, "'test'",
                                   'shutil.move(test_dir, prev_dir)'),
                                  (RUN_E2E_FULL, "'test.full'",
                                   'shutil.move(test_dir, prev_dir)')):
        lines = _code_lines(path)
        acq = _first_index(lines, 'runlock.acquire(REPO_ROOT, %s)' % resource, path)
        mov = _first_index(lines, mover, path)
        assert acq < mov, (
            '%s acquires the lock at line %d but rotates at line %d -- the '
            'rotation happens FIRST and the lock guards nothing'
            % (os.path.basename(path), acq + 1, mov + 1))
    print('  PASS  both tools acquire the lock BEFORE their rotation')


def check_the_rotation_itself_survived():
    """Rule 8's one preserved level of evidence. Deleting the rotation would
    make every check above pass and would silently destroy the previous run."""
    for path in (RUN_E2E, RUN_E2E_FULL):
        src = '\n'.join(_code_lines(path))
        assert 'shutil.move(test_dir, prev_dir)' in src, (
            '%s no longer rotates its run tree -- rule 8 keeps ONE level of '
            'evidence and this change removed it' % os.path.basename(path))
    print('  PASS  the rule-8 rotation is still there in both tools')


def check_the_two_tools_lock_different_trees():
    """They rotate two different trees, so one lock for both would make an
    e2e-full run block an ordinary sweep for no reason."""
    a = runlock.lock_path(ROOT, 'test')
    b = runlock.lock_path(ROOT, 'test.full')
    assert a != b, 'run_e2e and run_e2e_full would share one lock (%s)' % a
    for p in (a, b):
        assert os.path.dirname(p) == ROOT, (
            'lockfile %s is not beside the tree it guards; a lockfile INSIDE '
            'test/ would be rotated away by the rotation it serialises' % p)
    print('  PASS  separate lockfiles, both beside their tree: %s, %s'
          % (os.path.basename(a), os.path.basename(b)))


def main():
    print('Regression guard: a second concurrent e2e invocation refuses '
          'instead of rotating (item 70)')
    checks = [check_a_second_acquire_refuses,
              check_the_refusal_does_not_block,
              check_release_lets_the_next_invocation_in,
              check_a_dead_holder_is_taken_over_with_an_explicit_message,
              check_run_e2e_refuses_and_does_not_rotate,
              check_run_e2e_full_refuses_and_does_not_rotate,
              check_the_lock_is_taken_before_the_rotation_in_both_tools,
              check_the_rotation_itself_survived,
              check_the_two_tools_lock_different_trees]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_e2e_run_tree_lock (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_e2e_run_tree_lock (%d checks)' % len(checks))
    return 0


if __name__ == '__main__':
    sys.exit(main())
