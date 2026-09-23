#! /usr/bin/env python3
"""
Regression guard: a SECOND concurrent invocation of either `testsys/perf/`
tool that rebuilds a FIXED in-repo path must REFUSE, before it spends anything
(PROJECT_RULES.md rules 21a, 8, 2, 10; pathway item 74).

SIBLING OF test_e2e_run_tree_lock.py, NOT A REPLACEMENT. That guard pins the
lock mechanism and the two e2e entry points; item 70 left the same defect open
in two perf tools, and this file pins those. Nothing here re-tests the
mechanism on a single-component resource -- that is the other file's job and
it must keep passing untouched. What is new here is the NESTED resource path
(`testsys/perf/perf_case`, `src/python/eqdyna`) and the per-tool consequence
paragraph.

THE TWO COLLISIONS, and why neither shows up as an error:

  run_perf.py -- `build_perf_case()` rmtrees and rebuilds
  `testsys/perf/perf_case/test.tpv8` (and, through
  `run_tpv29_pinned_compare.py`, `perf_case_tpv29/test.tpv29`). A second run
  rebuilds that directory while the first is TIMING out of it. The first does
  not crash: it reports seconds, and the collision arrives as a perf
  regression or an improvement.

  run_jaxmpi_ab.py -- `stage()` COPIES an arm's `driver.py` and
  `MPI4NodalQuant.py` INTO `src/python/eqdyna/`, and every measurement is
  taken against whatever those two files currently hold. Two invocations swap
  each other's arm files mid-measurement and produce a number under the WRONG
  ARM LABEL, which that tool's own docstring calls worse than no number. The
  guarded resource is therefore the PACKAGE, not the `__pycache__` beneath it
  that the tool happens to delete.

WHAT IS ASSERTED. Each entry point is invoked FOR REAL with the lock already
held by this process, and four things are asserted together: non-zero exit;
a refusal naming the holder's pid; the refusal keeping the rule-2 shape (NOT
waiting, NOT falling back); and the protected directory UNCHANGED -- a
sentinel still present, and for the package, the two staged files identical by
content hash. An exit code alone would still pass if the tool refused AFTER
deleting.

AND THAT IT SPENT NOTHING. The refusal must land before the case is built,
before the vault is built, before the notes file is opened -- so the checks
also assert the tool produced none of those, and returned in seconds. This
guard NEVER runs a perf measurement and never launches mpirun.

NOTE ON A TRUE POSITIVE THAT LOOKS LIKE FLAKE, inherited from item 70's guard:
these checks take the real lock on THIS checkout's own paths. If a genuine
perf run is in flight in this checkout, they fail -- correctly; that is the
collision, and rule 21a says the run belongs in its own worktree.

Cheap (rule 9): flock in a tempdir, plus two entry-point invocations that
refuse at their first gate. Seconds. Exits non-zero on any failure.
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from testsys import runlock  # noqa: E402

RUN_PERF = os.path.join(ROOT, 'testsys', 'perf', 'run_perf.py')
RUN_JAXMPI_AB = os.path.join(ROOT, 'testsys', 'perf', 'run_jaxmpi_ab.py')

PERF_CASE_RESOURCE = os.path.join('testsys', 'perf', 'perf_case')
PKG_RESOURCE = os.path.join('src', 'python', 'eqdyna')
# The rule-2 shape every refusal must keep, whatever tool-specific consequence
# it carries. Asserted against the REAL output of each entry point, not
# against the constant the tool defines -- a constant checked against itself
# proves nothing.
REFUSAL_SHAPE = ('NOT waiting', 'NOT falling back')


def _sha(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


# --------------------------------------------------------------------------
# the nested-resource mechanism, in a tempdir -- no repo state touched
# --------------------------------------------------------------------------
def check_a_single_component_resource_is_unchanged():
    """The item-70 paths must not move. If `test`'s lockfile changed name, two
    invocations -- one at each version -- would take two DIFFERENT locks and
    believe they were serialised."""
    with tempfile.TemporaryDirectory() as tmp:
        for name in ('test', 'test.full'):
            got = runlock.lock_path(tmp, name)
            want = os.path.join(tmp, '.%s.lock' % name)
            assert got == want, ('lock_path(root, %r) moved from %r to %r -- '
                                 'item 70\'s lock is not the same lock any '
                                 'more' % (name, want, got))
        print('  PASS  single-component lockfiles unchanged (.test.lock, '
              '.test.full.lock at the root)')


def check_a_nested_resource_locks_beside_the_directory_it_guards():
    with tempfile.TemporaryDirectory() as tmp:
        for resource in (PERF_CASE_RESOURCE, PKG_RESOURCE):
            got = runlock.lock_path(tmp, resource)
            parent, name = os.path.split(resource)
            want = os.path.join(tmp, parent, '.%s.lock' % name)
            assert got == want, ('nested lock_path(%r) = %r, expected %r'
                                 % (resource, got, want))
            assert not got.startswith(os.path.join(tmp, resource) + os.sep), (
                'lockfile %s is INSIDE the directory it guards; for %s that '
                'puts it in the package run_jaxmpi_ab stages into'
                % (got, resource))
        print('  PASS  nested lockfiles sit beside their directory, not '
              'inside it')


def check_a_malformed_resource_is_refused_not_guessed():
    """A resource that escapes the repo, or names nothing, must raise. Quietly
    normalising it would put the lockfile somewhere plausible and guard
    nothing (rule 2)."""
    bad = ['', '.', '..', os.sep + 'abs', 'a' + os.sep + os.sep + 'b',
           'a' + os.sep + '..' + os.sep + 'b', 'test' + os.sep]
    with tempfile.TemporaryDirectory() as tmp:
        for resource in bad:
            try:
                got = runlock.lock_path(tmp, resource)
            except ValueError:
                continue
            raise AssertionError('lock_path accepted the malformed resource '
                                 '%r and returned %r' % (resource, got))
    print('  PASS  malformed resources (%d forms) raise instead of being '
          'normalised' % len(bad))


def check_a_second_acquire_on_a_nested_resource_refuses():
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, PKG_RESOURCE))
        first = runlock.acquire(tmp, PKG_RESOURCE, argv=['first'],
                                announce=False)
        t0 = time.time()
        try:
            runlock.acquire(tmp, PKG_RESOURCE, argv=['second'],
                            announce=False)
        except runlock.RunTreeLocked as exc:
            msg, dt = str(exc), time.time() - t0
        else:
            raise AssertionError(
                'a SECOND acquire on a held NESTED lock SUCCEEDED -- two A/B '
                'runs can stage arm files into one package again')
        finally:
            first.release()
        assert dt < 1.0, ('the refusal took %.1f s -- it is WAITING instead '
                          'of refusing' % dt)
        for want in (runlock.REFUSAL_HEADER, 'holder pid   : %d' % first.pid):
            assert want in msg, ('refusal omits %r; it read:\n%s' % (want, msg))
        print('  PASS  second acquire on a nested resource refuses in %.3f s, '
              'naming holder pid %d' % (dt, first.pid))


def check_a_custom_consequence_replaces_only_the_tail():
    """The caller supplies what ITS collision costs; the holder identification
    is not the caller's to override."""
    with tempfile.TemporaryDirectory() as tmp:
        first = runlock.acquire(tmp, 'test', argv=['first'], announce=False)
        try:
            runlock.acquire(tmp, 'test', argv=['second'], announce=False,
                            consequence=['CUSTOM TAIL LINE'])
        except runlock.RunTreeLocked as exc:
            msg = str(exc)
        else:
            raise AssertionError('second acquire succeeded')
        finally:
            first.release()
        assert 'CUSTOM TAIL LINE' in msg, msg
        assert runlock.E2E_CONSEQUENCE[0] not in msg, (
            'the custom consequence was APPENDED to the e2e text rather than '
            'replacing it; a refusal that describes two different incidents '
            'teaches the reader to skip it:\n%s' % msg)
        assert 'holder pid   : %d' % first.pid in msg, msg
        print('  PASS  a custom consequence replaces the tail and keeps the '
              'holder identification')


def check_every_guarded_resource_has_its_own_lockfile():
    paths = {r: runlock.lock_path(ROOT, r)
             for r in ('test', 'test.full', PERF_CASE_RESOURCE,
                       os.path.join('testsys', 'perf', 'perf_case_tpv29'),
                       PKG_RESOURCE)}
    assert len(set(paths.values())) == len(paths), (
        'two resources share one lockfile, so one tool would block an '
        'unrelated one: %s' % paths)
    print('  PASS  %d guarded resources, %d distinct lockfiles'
          % (len(paths), len(set(paths.values()))))


# --------------------------------------------------------------------------
# the entry points, for real, against this checkout's own paths
# --------------------------------------------------------------------------
def _hold(resource):
    try:
        return runlock.acquire(ROOT, resource,
                               argv=['test_perf_tool_locks.py', resource],
                               announce=False)
    except runlock.RunTreeLocked as exc:
        raise AssertionError(
            'could not take the lock on %s to run this check -- something '
            'else in THIS checkout holds it. That is the collision rule 21a '
            'forbids, not a broken guard: run your perf tool in its own git '
            'worktree. Holder:\n%s' % (os.path.join(ROOT, resource), exc))


def _refusal_asserts(tool, proc, out, held, elapsed):
    name = os.path.basename(tool)
    assert proc.returncode != 0, (
        '%s exited 0 while another invocation held the lock -- a refused run '
        'must never be green.\n%s' % (name, out[-2000:]))
    assert runlock.REFUSAL_HEADER in out, (
        '%s exited %d but printed no refusal naming the lock holder; its '
        'output was:\n%s' % (name, proc.returncode, out[-3000:]))
    assert 'holder pid   : %d' % held.pid in out, (
        '%s refused without naming the holder pid %d:\n%s'
        % (name, held.pid, out[-3000:]))
    for want in REFUSAL_SHAPE:
        assert want in out, (
            '%s refused but its message does not say %r -- a refusal that '
            'does not rule out waiting and falling back invites exactly '
            'those (rule 2):\n%s' % (name, want, out[-3000:]))
    assert elapsed < 120.0, (
        '%s took %.1f s to refuse -- it is doing work before Gate 0, which '
        'is the work the lock exists to prevent' % (name, elapsed))


def check_run_perf_refuses_and_does_not_rebuild_the_case():
    case_dir = os.path.join(ROOT, PERF_CASE_RESOURCE, 'test.tpv8')
    sentinel = os.path.join(case_dir, 'item74_lock_guard_%d.marker'
                            % os.getpid())
    held = _hold(PERF_CASE_RESOURCE)
    created = not os.path.isdir(case_dir)
    try:
        os.makedirs(case_dir, exist_ok=True)
        with open(sentinel, 'w') as fh:
            fh.write('pathway item 74 guard -- delete me if you find me\n')
        t0 = time.time()
        proc = subprocess.run([sys.executable, RUN_PERF], cwd=ROOT,
                              capture_output=True, text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        _refusal_asserts(RUN_PERF, proc, out, held, elapsed)
        assert os.path.isfile(sentinel), (
            'THE DEFECT ITSELF: run_perf.py rmtree\'d %s despite refusing -- '
            'the sentinel is gone. Refusing AFTER deleting is the collision '
            'with an exit code attached.' % case_dir)
        assert 'perf -- provenance' not in out and 'ms/step' not in out, (
            'run_perf.py refused but had already measured something; Gate 0 '
            'is too late:\n%s' % out[-2000:])
        print('  PASS  run_perf.py refused in %.2f s, exited %d, named the '
              'holder, and left %s untouched'
              % (elapsed, proc.returncode, os.path.basename(case_dir)))
    finally:
        try:
            os.remove(sentinel)
        except OSError:
            pass
        if created and os.path.isdir(case_dir) and not os.listdir(case_dir):
            os.rmdir(case_dir)
            parent = os.path.dirname(case_dir)
            if os.path.isdir(parent) and not os.listdir(parent):
                shutil.rmtree(parent)
        held.release()


def check_run_jaxmpi_ab_refuses_and_does_not_stage_into_the_package():
    pkg = os.path.join(ROOT, PKG_RESOURCE)
    staged = [os.path.join(pkg, n) for n in ('driver.py',
                                             'MPI4NodalQuant.py')]
    before = {p: _sha(p) for p in staged}
    cache = os.path.join(pkg, '__pycache__')
    sentinel = os.path.join(cache, 'item74_lock_guard_%d.marker' % os.getpid())
    notes = os.path.join(ROOT, 'NOTES_jaxmpi_ab_2026-09-22.md')
    notes_existed = os.path.exists(notes)
    held = _hold(PKG_RESOURCE)
    created = not os.path.isdir(cache)
    try:
        os.makedirs(cache, exist_ok=True)
        with open(sentinel, 'w') as fh:
            fh.write('pathway item 74 guard -- delete me if you find me\n')
        env = dict(os.environ)
        # The tool resolves its package from EQDYNAROOT when set. Pin it to
        # THIS checkout so the lock and the package under test are one tree;
        # an unpinned EQDYNAROOT would have this guard lock here and the tool
        # stage somewhere else, which is a real hazard but not this check's.
        env['EQDYNAROOT'] = ROOT
        t0 = time.time()
        proc = subprocess.run([sys.executable, RUN_JAXMPI_AB, '--ranks', '4'],
                              cwd=ROOT, env=env, capture_output=True,
                              text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        _refusal_asserts(RUN_JAXMPI_AB, proc, out, held, elapsed)
        for p in staged:
            assert _sha(p) == before[p], (
                'THE DEFECT ITSELF: run_jaxmpi_ab.py rewrote %s despite '
                'refusing -- an arm was staged into the package while '
                'another invocation was measuring against it, which is how a '
                'number gets the wrong arm label' % p)
        assert os.path.isfile(sentinel), (
            'run_jaxmpi_ab.py deleted %s despite refusing -- stage() ran '
            'before Gate 0' % cache)
        assert notes_existed or not os.path.exists(notes), (
            'run_jaxmpi_ab.py created %s despite refusing -- it is opening '
            'its notes file before the lock' % notes)
        print('  PASS  run_jaxmpi_ab.py refused in %.2f s, exited %d, named '
              'the holder, and left both staged files byte-identical'
              % (elapsed, proc.returncode))
    finally:
        try:
            os.remove(sentinel)
        except OSError:
            pass
        if created and os.path.isdir(cache) and not os.listdir(cache):
            os.rmdir(cache)
        held.release()


# --------------------------------------------------------------------------
# shape of the sources (rule 2a: assert the shape, not a substring elsewhere)
# --------------------------------------------------------------------------
def _code_lines(path):
    raw = open(path, errors='replace').read()
    return [ln.split('#', 1)[0] for ln in raw.splitlines()]


def _body_of(path, header):
    """The lines of one function, by indentation. Ordering must be checked
    INSIDE the function that does the damage: in run_perf.py `build_perf_case`
    is defined ABOVE `main`, so a whole-file line comparison would call a
    correct acquire late."""
    lines = _code_lines(path)
    for i, ln in enumerate(lines):
        if ln.startswith(header):
            break
    else:
        raise AssertionError('%s contains no %r' % (path, header))
    body = []
    for ln in lines[i + 1:]:
        if ln.strip() and not ln[0].isspace():
            break
        body.append(ln)
    return body


def _index(body, needle, path, header):
    for i, ln in enumerate(body):
        if needle in ln:
            return i
    raise AssertionError('%s: %s contains no %r' % (path, header, needle))


def check_each_tool_locks_before_it_destroys_anything():
    """A lock taken AFTER the rmtree or the copyfile is not a lock -- the
    damage is done by the time it is asked for."""
    for path, header, destroyers in (
            (RUN_PERF, 'def build_perf_case(',
             ['shutil.rmtree(PERF_CASE)', 'run_e2e.make_serial_case(']),
            (RUN_JAXMPI_AB, 'def main(',
             ['build_vault(work', 'stage(vault, arm)'])):
        body = _body_of(path, header)
        acq = _index(body, 'runlock.acquire(', path, header)
        for d in destroyers:
            dest = _index(body, d, path, header)
            assert acq < dest, (
                '%s: %s acquires the lock after %r -- the lock guards nothing'
                % (os.path.basename(path), header, d))
    print('  PASS  both tools acquire the lock before they delete, build or '
          'stage')


def check_the_perf_lock_follows_the_directory_actually_rebuilt():
    """run_tpv29_pinned_compare.py reassigns run_perf.PERF_CASE to
    perf_case_tpv29/ and calls build_perf_case() directly, never reaching
    main(). A lock on a hard-coded 'testsys/perf/perf_case' would guard the
    wrong directory for that caller and look like it guarded both."""
    body = _body_of(RUN_PERF, 'def build_perf_case(')
    src = '\n'.join(body)
    assert 'os.path.dirname(PERF_CASE)' in src, (
        'run_perf.build_perf_case does not derive its lock resource from '
        'PERF_CASE, so run_tpv29_pinned_compare.py rebuilds an unguarded '
        'directory:\n%s' % src)
    other = os.path.join(ROOT, 'testsys', 'perf',
                         'run_tpv29_pinned_compare.py')
    txt = open(other, errors='replace').read()
    assert 'run_perf.build_perf_case()' in txt, (
        '%s no longer goes through run_perf.build_perf_case -- it now needs '
        'a lock of its own' % os.path.basename(other))
    print('  PASS  the perf lock follows PERF_CASE, so '
          'run_tpv29_pinned_compare.py is covered by the same acquire')


def main():
    print('Regression guard: a second concurrent perf-tool invocation refuses '
          'instead of rebuilding (item 74)')
    checks = [check_a_single_component_resource_is_unchanged,
              check_a_nested_resource_locks_beside_the_directory_it_guards,
              check_a_malformed_resource_is_refused_not_guessed,
              check_a_second_acquire_on_a_nested_resource_refuses,
              check_a_custom_consequence_replaces_only_the_tail,
              check_every_guarded_resource_has_its_own_lockfile,
              check_run_perf_refuses_and_does_not_rebuild_the_case,
              check_run_jaxmpi_ab_refuses_and_does_not_stage_into_the_package,
              check_each_tool_locks_before_it_destroys_anything,
              check_the_perf_lock_follows_the_directory_actually_rebuilt]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_perf_tool_locks (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_perf_tool_locks (%d checks)' % len(checks))
    return 0


if __name__ == '__main__':
    sys.exit(main())
