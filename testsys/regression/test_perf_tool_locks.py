#! /usr/bin/env python3
"""
Regression guard: a SECOND concurrent invocation of either `testsys/perf/`
tool that rebuilds a FIXED in-repo path must REFUSE, before it spends anything
(PROJECT_RULES.md rules 21a, 21b, 8, 2, 10; pathway items 74 and 77).

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

  run_scaling.py / run_numa_scaling.py (item 77) -- `build_py_case()` and
  `build_case()` rmtree and rebuild `testsys/perf/scaling_case/<case>` and
  `testsys/perf/numa_case/<case>`. The first is called by FIVE tools
  (`run_scaling`, `run_mpi_scaling`, `run_jaxmpi_ab`, `run_setup_probe`,
  `run_shard_scaling`), so the two colliding runs need not be the same tool,
  and the lock therefore lives in the builder and not in any main().

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

BEHAVIOUR, NOT SOURCE SHAPE (2026-09-23). Three checks here used to assert
that `runlock.acquire(` appeared textually above `shutil.rmtree(d)` inside
`def build_py_case(` and friends, and that the literal
`os.path.relpath(os.path.dirname(d)` appeared there too. They caught real
defects, but by shape -- and so they also forbade the one repair that would
stop this defect class recurring a fifth time: the three near-identical
lock+rmtree+make_serial_case bodies cannot be replaced by ONE shared builder
while a guard demands `shutil.rmtree(d)` inside each tool's own function. They
are now observed instead, by running each builder with its destroyers replaced
by recorders that ask, at the instant of every write, whether the lock on the
directory being written is already held. See the BUILDERS table below.

THE HOLE A LOCK CANNOT CLOSE (item 77). `run_jaxmpi_ab.py` resolves its ROOT
from `$EQDYNAROOT`, so a stale export -- the shape `install-eqdyna.sh` leaves,
since it exports `$(pwd)` and the export survives a `cd` into a worktree --
makes it stage arm files into ANOTHER checkout's package. Gate 0 does not catch
that: the lock follows the same wrong ROOT and so locks and corrupts one
consistent wrong tree. The tool therefore REFUSES on the mismatch before any
lock, and this file pins both directions -- the refusal, and the negative
control that it does NOT fire when the roots agree or when the variable is
unset. Without that control, every other check on this tool would pass for the
wrong reason.

NOTE ON A TRUE POSITIVE THAT LOOKS LIKE FLAKE, inherited from item 70's guard:
these checks take the real lock on THIS checkout's own paths. If a genuine
perf run is in flight in this checkout, they fail -- correctly; that is the
collision, and rule 21a says the run belongs in its own worktree.

Cheap (rule 9): flock in a tempdir, plus two entry-point invocations that
refuse at their first gate. Seconds. Exits non-zero on any failure.
"""
import hashlib
import json
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
RUN_SCALING = os.path.join(ROOT, 'testsys', 'perf', 'run_scaling.py')
RUN_NUMA_SCALING = os.path.join(ROOT, 'testsys', 'perf', 'run_numa_scaling.py')
PERF_DIR = os.path.join(ROOT, 'testsys', 'perf')
PARITY_DIR = os.path.join(ROOT, 'testsys', 'parity')
PROBE_PLASTIC_TRACTION = os.path.join(PARITY_DIR, 'probe_plastic_traction.py')

PERF_CASE_RESOURCE = os.path.join('testsys', 'perf', 'perf_case')
PKG_RESOURCE = os.path.join('src', 'python', 'eqdyna')
SCALING_CASE_RESOURCE = os.path.join('testsys', 'perf', 'scaling_case')
NUMA_CASE_RESOURCE = os.path.join('testsys', 'perf', 'numa_case')
PROBE_CASE_RESOURCE = os.path.join('testsys', 'parity', 'probe_case')
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
                       PKG_RESOURCE, SCALING_CASE_RESOURCE,
                       NUMA_CASE_RESOURCE)}
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
        # main() builds its arm VAULT under $TMPDIR before it stages anything.
        # Point TMPDIR at an empty directory of our own: whether the vault was
        # built before Gate 0 is then a fact on disk, not a reading of the
        # order of two lines in main().
        tmpdir = tempfile.mkdtemp(prefix='jaxmpiab_guard.')
        env['TMPDIR'] = tmpdir
        t0 = time.time()
        proc = subprocess.run([sys.executable, RUN_JAXMPI_AB, '--ranks', '4'],
                              cwd=ROOT, env=env, capture_output=True,
                              text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        _refusal_asserts(RUN_JAXMPI_AB, proc, out, held, elapsed)
        left = sorted(os.listdir(tmpdir))
        shutil.rmtree(tmpdir, ignore_errors=True)
        assert not left, (
            'run_jaxmpi_ab.py built its arm vault under $TMPDIR despite '
            'refusing -- build_vault() ran before Gate 0, so a refused '
            'invocation has already spent a git worktree checkout: %s' % left)
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
              'the holder, left both staged files byte-identical, and built '
              'no vault ($TMPDIR still empty)' % (elapsed, proc.returncode))
    finally:
        try:
            os.remove(sentinel)
        except OSError:
            pass
        if created and os.path.isdir(cache) and not os.listdir(cache):
            os.rmdir(cache)
        held.release()


# --------------------------------------------------------------------------
# the two shared case BUILDERS (item 77), for real, against this checkout
#
# Invoked as the four other perf tools invoke them -- `rs.build_py_case(case)`
# and `numa.build_case(case)` in a fresh interpreter -- and NOT through either
# module's main(). main() first probes NUMA topology and cpu busy-ness and
# exits non-zero when the box is loaded, which would make this guard pass or
# fail on this box's tenancy rather than on the lock. Calling the builder
# directly IS the real code path for `run_mpi_scaling.py:408`,
# `run_jaxmpi_ab.py:169`, `run_setup_probe.py:142` and
# `run_shard_scaling.py:195`, which is the collision item 77 names.
# --------------------------------------------------------------------------
def _call_builder(module, func, case, env=None, cwd=ROOT):
    snippet = ('import sys; sys.path.insert(0, %r)\n'
               'import %s as m\n'
               'm.%s(%r)\n' % (PERF_DIR, module, func, case))
    t0 = time.time()
    proc = subprocess.run([sys.executable, '-c', snippet], cwd=cwd,
                          env=env or dict(os.environ), capture_output=True,
                          text=True, timeout=300)
    return proc, proc.stdout + proc.stderr, time.time() - t0


def _builder_refuses(label, tool, module, func, resource, case):
    """Hold the lock, run the builder for real, assert it refused and built
    nothing. The sentinel is the whole point: an exit code alone would still
    pass if the builder rmtree'd the case and THEN refused."""
    case_dir = os.path.join(ROOT, resource, case)
    sentinel = os.path.join(case_dir, 'item77_lock_guard_%d.marker'
                            % os.getpid())
    held = _hold(resource)
    created = not os.path.isdir(case_dir)
    try:
        os.makedirs(case_dir, exist_ok=True)
        with open(sentinel, 'w') as fh:
            fh.write('pathway item 77 guard -- delete me if you find me\n')
        proc, out, elapsed = _call_builder(module, func, case)
        _refusal_asserts(tool, proc, out, held, elapsed)
        assert os.path.isfile(sentinel), (
            'THE DEFECT ITSELF: %s.%s rmtree\'d %s despite refusing -- the '
            'sentinel is gone. Refusing AFTER deleting is the collision with '
            'an exit code attached.' % (module, func, case_dir))
        assert 'Creating case' not in out and 'case.setup' not in out, (
            '%s.%s refused but had already started building the case; Gate 0 '
            'is too late:\n%s' % (module, func, out[-2000:]))
        print('  PASS  %s refused in %.2f s, exited %d, named the holder, and '
              'left %s untouched' % (label, elapsed, proc.returncode,
                                     os.path.join(resource, case)))
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


def check_run_scaling_builder_refuses_and_does_not_rebuild_the_case():
    _builder_refuses('run_scaling.build_py_case', RUN_SCALING, 'run_scaling',
                     'build_py_case', SCALING_CASE_RESOURCE, 'test.tpv8')


def check_run_numa_scaling_builder_refuses_and_does_not_rebuild_the_case():
    _builder_refuses('run_numa_scaling.build_case', RUN_NUMA_SCALING,
                     'run_numa_scaling', 'build_case', NUMA_CASE_RESOURCE,
                     'test.tpv8')


def check_probe_plastic_traction_build_case_refuses_and_does_not_rebuild():
    """pathway item 81, the same class as items 74/77:
    `probe_plastic_traction.build_case()` rmtrees and rebuilds
    `testsys/parity/probe_case/test.drv.a6` through the e2e sweep's own
    `make_serial_case`. A second invocation deletes the first's live case
    tree mid-run, and the failure reads as a plastic-traction mismatch
    rather than as infrastructure -- exactly the false-red rule 21a exists
    to stop.

    Calls `build_case()` directly (not `main()`, which also runs the
    physics measurement) -- the same shape `_builder_refuses` uses for
    `run_scaling.build_py_case` / `run_numa_scaling.build_case`, adapted
    because this builder takes a full case_dir rather than a bare case
    name."""
    case_dir = os.path.join(ROOT, PROBE_CASE_RESOURCE, 'test.drv.a6')
    sentinel = os.path.join(case_dir, 'item81_lock_guard_%d.marker'
                            % os.getpid())
    held = _hold(PROBE_CASE_RESOURCE)
    created = not os.path.isdir(case_dir)
    try:
        os.makedirs(case_dir, exist_ok=True)
        with open(sentinel, 'w') as fh:
            fh.write('pathway item 81 guard -- delete me if you find me\n')
        mtime_before = os.path.getmtime(case_dir)
        snippet = ('import sys; sys.path.insert(0, %r)\n'
                   'import probe_plastic_traction as p\n'
                   'p.build_case(%r)\n' % (PARITY_DIR, case_dir))
        t0 = time.time()
        proc = subprocess.run([sys.executable, '-c', snippet], cwd=ROOT,
                              capture_output=True, text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        _refusal_asserts(PROBE_PLASTIC_TRACTION, proc, out, held, elapsed)
        assert os.path.isfile(sentinel), (
            'THE DEFECT ITSELF: probe_plastic_traction.build_case rmtree\'d '
            '%s despite refusing -- the sentinel is gone. Refusing AFTER '
            'deleting is the collision with an exit code attached.'
            % case_dir)
        assert os.path.getmtime(case_dir) == mtime_before, (
            'probe_plastic_traction.build_case touched %s despite refusing '
            '(mtime moved from %r)' % (case_dir, mtime_before))
        print('  PASS  probe_plastic_traction.build_case refused in %.2f s, '
              'exited %d, named the holder, and left %s untouched'
              % (elapsed, proc.returncode,
                 os.path.join(PROBE_CASE_RESOURCE, 'test.drv.a6')))
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


def check_the_builder_lock_ignores_a_foreign_eqdynaroot():
    """`run_scaling.ROOT` honours $EQDYNAROOT, but the case directory is built
    from __file__. If the lock followed ROOT it would be taken in the OTHER
    checkout -- guarding nothing here while this checkout's directory is
    rebuilt. So: hold the lock HERE, point $EQDYNAROOT at an empty tree, and
    the builder must still refuse."""
    case = 'test.tpv8'
    case_dir = os.path.join(ROOT, SCALING_CASE_RESOURCE, case)
    sentinel = os.path.join(case_dir, 'item77_root_guard_%d.marker'
                            % os.getpid())
    held = _hold(SCALING_CASE_RESOURCE)
    created = not os.path.isdir(case_dir)
    with tempfile.TemporaryDirectory() as foreign:
        try:
            os.makedirs(case_dir, exist_ok=True)
            with open(sentinel, 'w') as fh:
                fh.write('pathway item 77 guard -- delete me if you find me\n')
            env = dict(os.environ)
            env['EQDYNAROOT'] = foreign
            proc, out, elapsed = _call_builder('run_scaling', 'build_py_case',
                                               case, env=env)
            _refusal_asserts(RUN_SCALING, proc, out, held, elapsed)
            assert runlock.lock_path(ROOT, SCALING_CASE_RESOURCE) in out, (
                'run_scaling refused, but not on THIS checkout\'s lockfile -- '
                'its lock followed $EQDYNAROOT=%s:\n%s' % (foreign, out[-2000:]))
            assert os.path.isfile(sentinel), (
                'run_scaling.build_py_case rebuilt %s under a foreign '
                '$EQDYNAROOT' % case_dir)
            assert not os.listdir(foreign), (
                'run_scaling.build_py_case wrote into the foreign $EQDYNAROOT '
                '%s: %s' % (foreign, os.listdir(foreign)))
            print('  PASS  the builder lock is rooted in its own checkout, '
                  'not in $EQDYNAROOT (refused in %.2f s)' % elapsed)
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


def check_run_jaxmpi_ab_refuses_a_foreign_eqdynaroot():
    """The hole the LOCK cannot close (item 77). `stage()` writes into
    `$EQDYNAROOT/src/python/eqdyna`, and Gate 0 locks that same wrong tree --
    consistently, so the lock is no evidence of safety. The tool must refuse
    on the mismatch BEFORE the lock, and must not write into either tree.

    THE FOREIGN TREE IS A REAL SKELETON, not an empty directory (2026-09-23).
    It carries `src/python/eqdyna/` with both staged files, so an acquire
    reaching it WOULD succeed and WOULD leave `src/python/.eqdyna.lock`
    behind. That is what makes `no lockfile in the foreign tree` a load-
    bearing assertion and what replaces the retired source-order check that
    `require_root_is_this_checkout()` is written above `runlock.acquire(` in
    main(): against an empty directory the acquire would merely have crashed
    on a missing parent, and the ordering would have gone unobserved."""
    sys.path.insert(0, PERF_DIR)
    pkg = os.path.join(ROOT, PKG_RESOURCE)
    staged = [os.path.join(pkg, n) for n in ('driver.py',
                                             'MPI4NodalQuant.py')]
    before = {p: _sha(p) for p in staged}
    with tempfile.TemporaryDirectory() as foreign:
        foreign_pkg = os.path.join(foreign, PKG_RESOURCE)
        os.makedirs(foreign_pkg)
        for p in staged:
            shutil.copyfile(p, os.path.join(foreign_pkg,
                                            os.path.basename(p)))
        foreign_before = _fingerprint(foreign)
        foreign_lock = runlock.lock_path(foreign, PKG_RESOURCE)
        env = dict(os.environ)
        env['EQDYNAROOT'] = foreign
        t0 = time.time()
        proc = subprocess.run([sys.executable, RUN_JAXMPI_AB, '--ranks', '4'],
                              cwd=ROOT, env=env, capture_output=True,
                              text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        assert proc.returncode != 0, (
            'run_jaxmpi_ab.py exited 0 with $EQDYNAROOT pointing at another '
            'tree -- it would stage arm files into a checkout it does not '
            'own.\n%s' % out[-2000:])
        header = 'refusing to start: $EQDYNAROOT names a DIFFERENT checkout'
        assert header in out, (
            'run_jaxmpi_ab.py exited %d under a foreign $EQDYNAROOT but not '
            'on the mismatch -- some later accident stopped it, which is not '
            'a guarantee. Its output was:\n%s' % (proc.returncode, out[-3000:]))
        for want in (foreign, ROOT):
            assert want in out, (
                'the mismatch refusal does not name %r, so the reader cannot '
                'see which two trees disagree:\n%s' % (want, out[-3000:]))
        assert 'NOT warning and continuing' in out, (
            'the mismatch refusal does not rule out warning-and-continuing, '
            'which is the behaviour it was chosen over (rule 2):\n%s'
            % out[-3000:])
        for p in staged:
            assert _sha(p) == before[p], (
                'run_jaxmpi_ab.py rewrote %s despite refusing on the '
                '$EQDYNAROOT mismatch' % p)
        assert not os.path.exists(foreign_lock), (
            'run_jaxmpi_ab.py CREATED %s -- it took Gate 0 in the foreign '
            'checkout before checking $EQDYNAROOT, so a refused run has '
            'already written to a tree this session does not own'
            % foreign_lock)
        foreign_after = _fingerprint(foreign)
        assert foreign_after == foreign_before, (
            'run_jaxmpi_ab.py wrote into the foreign $EQDYNAROOT despite '
            'refusing:\n  before: %s\n  after : %s'
            % (foreign_before, foreign_after))
        assert elapsed < 120.0, (
            'run_jaxmpi_ab.py took %.1f s to refuse the mismatch -- the check '
            'is not the first thing main() does' % elapsed)
        nfiles = len(foreign_before['files'])
    print('  PASS  run_jaxmpi_ab.py refused a foreign $EQDYNAROOT in %.2f s, '
          'exited %d, named both trees, staged nothing, and left no %s (%d '
          'files compared, inode and mtime unchanged)'
          % (elapsed, proc.returncode, os.path.basename(foreign_lock),
             nfiles))


def check_the_mismatch_refusal_does_not_fire_when_the_roots_agree():
    """The negative control. A check that refuses unconditionally would make
    every check above pass for the wrong reason: the tool would never reach
    Gate 0 at all, and `check_run_jaxmpi_ab_refuses_and_does_not_stage...`
    would be asserting on the wrong refusal."""
    sys.path.insert(0, PERF_DIR)
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    snippet = ('import sys; sys.path.insert(0, %r)\n'
               'import run_jaxmpi_ab as ab\n'
               'ab.require_root_is_this_checkout()\n'
               'print("ROOTS AGREE")\n' % PERF_DIR)
    proc = subprocess.run([sys.executable, '-c', snippet], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=300)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0 and 'ROOTS AGREE' in out, (
        'require_root_is_this_checkout() refused with $EQDYNAROOT set to this '
        'very checkout (%s) -- it refuses unconditionally, so every other '
        'check on this tool proves nothing:\n%s' % (ROOT, out[-2000:]))
    del env['EQDYNAROOT']
    proc = subprocess.run([sys.executable, '-c', snippet], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=300)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0 and 'ROOTS AGREE' in out, (
        'require_root_is_this_checkout() refused with $EQDYNAROOT UNSET, where '
        'ROOT falls back to the tool\'s own location and no mismatch is '
        'possible:\n%s' % out[-2000:])
    print('  PASS  the $EQDYNAROOT check passes when the roots agree and when '
          'the variable is unset')


# --------------------------------------------------------------------------
# WHAT THE BUILDERS DO, OBSERVED -- not what their source says (2026-09-23)
#
# Three checks here used to read the source: they required the literal
# `shutil.rmtree(d)` / `run_e2e.make_serial_case(` to appear textually AFTER
# `runlock.acquire(` inside `def build_py_case(` and friends, and the literal
# `os.path.relpath(os.path.dirname(d)` to appear there too. They caught real
# defects -- a missing acquire, an acquire placed after the rmtree, a lock
# resource hard-coded away from the directory being destroyed -- but they
# caught them by SHAPE, and so they also forbade the one fix that would stop
# this defect class recurring: three tools carry near-identical
# lock+rmtree+make_serial_case bodies, and a SHARED builder cannot make
# `shutil.rmtree(d)` appear inside `def build_py_case(`. A guard that forbids
# the repair is not a guard.
#
# The replacement observes the same three properties by RUNNING each builder
# with its destroyers intercepted:
#
#   `_run_probe` launches the builder in a fresh interpreter with
#   `shutil.rmtree`, `os.makedirs` and `run_e2e.make_serial_case` replaced by
#   recorders that DESTROY NOTHING and, at the instant of each call, ask
#   whether the lock on the directory about to be written is ALREADY HELD --
#   by attempting a second acquire from a second open file description. flock
#   treats descriptions independently even inside one process, so a refusal
#   there means the builder itself holds it and a success means nobody does.
#
# That yields, per builder, the ordered list of destructive acts, the resource
# guarding each one, and whether the lock was held at that moment -- and the
# case directory is verified byte-for-byte unchanged afterwards, which is also
# what proves the interception took. A shared builder satisfies all of it; a
# tool that destroys before locking cannot.
# --------------------------------------------------------------------------
_PROBE_SRC = r'''
import json
import os
import shutil
import sys

ROOT = %(root)r
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))
sys.path.insert(0, %(tooldir)r)
from testsys import runlock
import run_e2e          # pre-imported, so the tool's own `import run_e2e`
                        # returns THIS module object -- the one patched below
import %(module)s as tool
%(setup)s

EVENTS = []


def lock_state(guarded_dir):
    """Is the lock on `guarded_dir` held RIGHT NOW, by anybody?"""
    resource = os.path.relpath(os.path.abspath(guarded_dir), ROOT)
    try:
        lk = runlock.acquire(ROOT, resource, argv=['probe'], announce=False)
    except runlock.RunTreeLocked:
        return resource, True
    lk.release()
    return resource, False


def note(kind, target, guarded):
    resource, held = lock_state(guarded)
    EVENTS.append(dict(kind=kind, target=os.path.abspath(target),
                       resource=resource, held=held))


def parent(p):
    return os.path.dirname(os.path.abspath(p))


# The destroyers, replaced by recorders. `os.makedirs` is called ON the
# guarded directory; the other two are called on a case directory inside it.
shutil.rmtree = lambda p, *a, **k: note('rmtree', p, parent(p))
os.makedirs = lambda p, *a, **k: note('makedirs', p, p)
run_e2e.make_serial_case = lambda n, d, *a, **k: (
    note('make_serial_case', d, parent(d)), d)[1]

returned = tool.%(call)s
_, held_after = lock_state(parent(returned))
print('PROBE ' + json.dumps(dict(events=EVENTS, held_after=held_after,
                                 returned=os.path.abspath(returned))))
'''

# Every entry point that rebuilds a fixed in-repo case directory. The second
# row is `run_tpv29_pinned_compare.py:55-56` exactly: it repoints run_perf's
# module globals at perf_case_tpv29/ and calls build_perf_case() directly,
# never reaching main(). A lock hard-coded to 'testsys/perf/perf_case' would
# guard the wrong directory for that caller and look like it guarded both.
BUILDERS = [
    dict(label='run_perf.build_perf_case',
         tooldir=PERF_DIR, module='run_perf', setup='',
         call='build_perf_case()',
         resource=PERF_CASE_RESOURCE, case='test.tpv8',
         holds_after=True,
         holds_why='run_perf then TIMES out of the directory it just built'),
    dict(label='run_perf.build_perf_case [PERF_CASE repointed, as '
               'run_tpv29_pinned_compare.py does]',
         tooldir=PERF_DIR, module='run_perf',
         setup=("tool.PERF_CASE_NAME = 'test.tpv29'\n"
                "tool.PERF_CASE = os.path.join(tool.TESTSYS, "
                "'perf_case_tpv29', 'test.tpv29')"),
         call='build_perf_case()',
         resource=os.path.join('testsys', 'perf', 'perf_case_tpv29'),
         case='test.tpv29',
         holds_after=True,
         holds_why='the pinned compare then TIMES out of that directory'),
    dict(label='run_scaling.build_py_case',
         tooldir=PERF_DIR, module='run_scaling', setup='',
         call="build_py_case('test.tpv8')",
         resource=SCALING_CASE_RESOURCE, case='test.tpv8',
         holds_after=True,
         holds_why='five tools time out of this directory after building it'),
    dict(label='run_numa_scaling.build_case',
         tooldir=PERF_DIR, module='run_numa_scaling', setup='',
         call="build_case('test.tpv8')",
         resource=NUMA_CASE_RESOURCE, case='test.tpv8',
         holds_after=True,
         holds_why='the NUMA sweep then times out of this directory'),
    dict(label='probe_plastic_traction.build_case',
         tooldir=PARITY_DIR, module='probe_plastic_traction', setup='',
         call='build_case(%r)' % os.path.join(ROOT, PROBE_CASE_RESOURCE,
                                              'test.drv.a6'),
         resource=PROBE_CASE_RESOURCE, case='test.drv.a6',
         holds_after=False,
         holds_why='its documented contract: held across the WRITE only, '
                   'released before a measurement phase that never touches '
                   'the tree destructively'),
]

_PROBED = {}


def _fingerprint(d):
    """Inode, mtime and the full file listing with sizes and content hashes.
    Printed, not just compared: a directory that was rebuilt identically and
    a directory that was never touched differ here by inode and mtime."""
    st = os.stat(d)
    files = []
    for dirpath, dirnames, filenames in os.walk(d):
        dirnames.sort()
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            files.append([os.path.relpath(p, d), os.path.getsize(p), _sha(p)])
    return dict(inode=st.st_ino, mtime=st.st_mtime, files=files)


def _run_probe(entry):
    """Run one builder with its destroyers intercepted; return what it did.

    Memoised: three checks read this and the launch is the expensive part.
    """
    if entry['label'] in _PROBED:
        return _PROBED[entry['label']]
    case_dir = os.path.join(ROOT, entry['resource'], entry['case'])
    sentinel = os.path.join(case_dir, 'lock_probe_%d.marker' % os.getpid())
    created = not os.path.isdir(case_dir)
    try:
        os.makedirs(case_dir, exist_ok=True)
        with open(sentinel, 'w') as fh:
            fh.write('behavioural lock probe -- delete me if you find me\n')
        before = _fingerprint(case_dir)
        src = _PROBE_SRC % dict(root=ROOT, tooldir=entry['tooldir'],
                                module=entry['module'], setup=entry['setup'],
                                call=entry['call'])
        t0 = time.time()
        proc = subprocess.run([sys.executable, '-c', src], cwd=ROOT,
                              capture_output=True, text=True, timeout=300)
        elapsed = time.time() - t0
        out = proc.stdout + proc.stderr
        line = [ln for ln in out.splitlines() if ln.startswith('PROBE ')]
        assert proc.returncode == 0 and line, (
            '%s could not be observed: the probe exited %d and printed no '
            'PROBE record. Its output was:\n%s'
            % (entry['label'], proc.returncode, out[-3000:]))
        data = json.loads(line[-1][len('PROBE '):])
        after = _fingerprint(case_dir)
        assert after == before, (
            'THE PROBE ITSELF DESTROYED SOMETHING, or %s writes by a route '
            'this probe does not intercept. %s changed:\n  before: %s\n  '
            'after : %s' % (entry['label'], case_dir, before, after))
        data['elapsed'] = elapsed
        data['case_dir'] = case_dir
        data['fingerprint'] = before
        _PROBED[entry['label']] = data
        return data
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


def check_every_builder_holds_the_lock_at_its_first_destructive_act():
    """Replaces the source-ordering check. A lock taken AFTER the rmtree is
    not a lock -- the damage is done by the time it is asked for -- and a
    builder with no acquire at all is the same defect with no line to point
    at. Both show up here as `held: false` at a recorded destructive act."""
    for entry in BUILDERS:
        data = _run_probe(entry)
        events = data['events']
        kinds = [e['kind'] for e in events]
        assert 'rmtree' in kinds, (
            '%s never called shutil.rmtree even though %s existed with a '
            'sentinel in it -- this probe observed nothing, so it proves '
            'nothing. Recorded: %s'
            % (entry['label'], data['case_dir'], kinds))
        unheld = [e for e in events if not e['held']]
        assert not unheld, (
            'THE DEFECT ITSELF: %s performed %s on %s with the lock on %s NOT '
            'held. Recorded sequence: %s'
            % (entry['label'], unheld[0]['kind'], unheld[0]['target'],
               unheld[0]['resource'],
               [(e['kind'], e['held']) for e in events]))
        print('  PASS  %s held %s at all %d destructive acts (%s), in %.2f s'
              % (entry['label'], events[0]['resource'], len(events),
                 ', '.join(kinds), data['elapsed']))


def check_every_builder_locks_the_directory_it_actually_rebuilds():
    """Replaces the `os.path.relpath(os.path.dirname(d)` literal. Five tools
    call `rs.build_py_case`; the lock must follow the directory that is
    actually destroyed, not a second copy of that path that can drift. Here
    the resource is read back off the WRITE TARGET the builder passed to its
    destroyers, and the held-ness of that exact resource is what the previous
    check asserted."""
    for entry in BUILDERS:
        data = _run_probe(entry)
        want = entry['resource']
        for e in data['events']:
            assert e['resource'] == want, (
                '%s %s %s, which is guarded by %r -- not by %r, the resource '
                'this guard and the tool\'s callers believe is locked. The '
                'lock and the directory have drifted apart.'
                % (entry['label'], e['kind'], e['target'], e['resource'],
                   want))
        assert os.path.dirname(data['returned']) == os.path.join(ROOT, want), (
            '%s returned %s, outside the guarded %s'
            % (entry['label'], data['returned'], want))
        print('  PASS  %s wrote only inside %s -- %d writes (%s) -- and '
              'returned %s'
              % (entry['label'], want, len(data['events']),
                 ', '.join('%s %s' % (e['kind'],
                                      os.path.relpath(e['target'], ROOT))
                           for e in data['events']),
                 os.path.relpath(data['returned'], ROOT)))


def check_each_builder_keeps_its_lock_until_its_write_is_finished():
    """A lock released mid-build is a lock that lets a second invocation in
    while the first is still writing; a lock released before the caller has
    finished MEASURING out of the tree is the item-74 collision again. Each
    builder declares which of the two it is, and the probe asks the lock."""
    for entry in BUILDERS:
        data = _run_probe(entry)
        got = data['held_after']
        assert got == entry['holds_after'], (
            '%s left the lock on %s %s when it returned; expected %s -- %s'
            % (entry['label'], entry['resource'],
               'HELD' if got else 'RELEASED',
               'HELD' if entry['holds_after'] else 'RELEASED',
               entry['holds_why']))
        print('  PASS  %s returned with the lock on %s %s -- %s'
              % (entry['label'], entry['resource'],
                 'still held' if got else 'released',
                 entry['holds_why']))


def check_the_tpv29_compare_tool_goes_through_the_guarded_builder():
    """The one call-site fact the probe cannot observe without running a perf
    measurement: that `run_tpv29_pinned_compare.py` reaches its rebuild
    THROUGH `run_perf.build_perf_case` rather than rebuilding
    perf_case_tpv29/ itself. The probe above proves that entry point is
    guarded when PERF_CASE is repointed; this proves that is the entry point
    the tool uses. If the tool ever grows its own rmtree, it needs its own
    lock and its own row in BUILDERS."""
    other = os.path.join(ROOT, 'testsys', 'perf',
                         'run_tpv29_pinned_compare.py')
    txt = open(other, errors='replace').read()
    assert 'run_perf.build_perf_case()' in txt, (
        '%s no longer goes through run_perf.build_perf_case -- it now needs '
        'a lock of its own' % os.path.basename(other))
    own = [n for n in ('shutil.rmtree', 'make_serial_case') if n in txt]
    assert not own, (
        '%s destroys or builds a case tree itself (%s) instead of delegating '
        'to the guarded builder' % (os.path.basename(other), own))
    print('  PASS  run_tpv29_pinned_compare.py (%d bytes) delegates to '
          'run_perf.build_perf_case and carries no rmtree or make_serial_case '
          'of its own' % len(txt))


def main():
    print('Regression guard: a second concurrent perf-tool invocation refuses '
          'instead of rebuilding, and no tool writes a checkout $EQDYNAROOT '
          'points at (items 74, 77, 81)')
    checks = [check_a_single_component_resource_is_unchanged,
              check_a_nested_resource_locks_beside_the_directory_it_guards,
              check_a_malformed_resource_is_refused_not_guessed,
              check_a_second_acquire_on_a_nested_resource_refuses,
              check_a_custom_consequence_replaces_only_the_tail,
              check_every_guarded_resource_has_its_own_lockfile,
              check_run_perf_refuses_and_does_not_rebuild_the_case,
              check_run_jaxmpi_ab_refuses_and_does_not_stage_into_the_package,
              check_run_scaling_builder_refuses_and_does_not_rebuild_the_case,
              check_run_numa_scaling_builder_refuses_and_does_not_rebuild_the_case,
              check_probe_plastic_traction_build_case_refuses_and_does_not_rebuild,
              check_the_builder_lock_ignores_a_foreign_eqdynaroot,
              check_run_jaxmpi_ab_refuses_a_foreign_eqdynaroot,
              check_the_mismatch_refusal_does_not_fire_when_the_roots_agree,
              check_every_builder_holds_the_lock_at_its_first_destructive_act,
              check_every_builder_locks_the_directory_it_actually_rebuilds,
              check_each_builder_keeps_its_lock_until_its_write_is_finished,
              check_the_tpv29_compare_tool_goes_through_the_guarded_builder]
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
