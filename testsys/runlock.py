#! /usr/bin/env python3
"""
One live writer per rotated run tree (PROJECT_RULES.md rule 21a, pathway
item 70).

WHAT THIS GUARDS. `testsys/e2e/run_e2e.py` and `testsys/e2e/run_e2e_full.py`
each begin by ROTATING a FIXED path under REPO_ROOT -- `test/` -> `test.prev/`
and `test.full/` -> `test.full.prev/` -- so that rule 8's one level of evidence
survives a re-run. The rotation is correct. The fixed path with no owner is
not: a second invocation in the same checkout renames the FIRST one's LIVE tree
out from under it, mid-run, and the first notices only when it next touches an
absolute path it resolved earlier.

Measured cost, 2026-09-22 (docs/evidence/gate-v5.16.0/README.md):
`test.tpv1053d x python-numpy` ran 1500.1 s and then died with
`FileNotFoundError .../test/test.tpv1053d.python-numpy/frt.txt0` --
`src/python/eqdyna/library_output.py:120` writes that file by absolute path at
the very END of a cell. The sweep printed `FAIL test.tpv1053d x python-numpy`,
character for character the string a genuine parity breach prints. Run
isolated, the same cell passed in 907 s. Four more cells were doomed the same
way and were killed by PID; ~40 minutes of a release session, plus a restarted
gate.

WHY flock AND NOT AN O_EXCL LOCKFILE. Both refuse a second holder. They differ
on what happens when the holder DIES:

  - flock is held by the open FILE DESCRIPTION and is released BY THE KERNEL
    when the holder exits, crashes, or is SIGKILLed. A dead holder leaves
    nothing to clean up and nothing to decide.
  - an O_EXCL lockfile OUTLIVES its holder. Every acquisition must then judge
    whether an existing file is live or abandoned, and the only test available
    is PID liveness -- which is wrong after PID reuse, in both directions. That
    judgement sits on the critical path of every run, and a wrong one either
    wedges the checkout until a human deletes a file, or hands out a lock that
    is genuinely held. A sweep is a multi-hour job on a shared box; a lock that
    needs a human to unstick it after a SIGKILL would be worked around, and a
    worked-around lock is not a lock.

`testsys/perf/ledger.py` already uses flock (item 53), so this is the existing
convention in this tree rather than a new one.

WHAT flock CANNOT DO is name the holder, so the holder writes an identity
record -- pid, start time, cwd, argv -- INTO the locked file, and a refuser
reads it. The record is PROVENANCE; flock is the lock. When the record is
missing or unreadable the refusal still stands and says exactly that, rather
than inventing a holder or letting the caller through (rule 2).

STALE CONTENT IS NOT A STALE LOCK. The holder TRUNCATES its record on release,
so a non-empty record at acquisition time means exactly one thing: the previous
holder died without releasing. flock has already let us in, because the kernel
dropped the dead holder's lock. Taking over is therefore correct and automatic
-- but never silent: acquire() PRINTS who died, when, and why the takeover is
safe, before it returns.

REUSE. Nothing here is e2e-specific: `resource` is any directory path relative
to REPO_ROOT. Item 74 adopted it for the two `testsys/perf/` tools of the same
class -- `run_perf.build_perf_case` (rebuilds `testsys/perf/perf_case/<name>`,
and `perf_case_tpv29/<name>` under `run_tpv29_pinned_compare.py`) and
`run_jaxmpi_ab.stage` (rewrites `src/python/eqdyna/driver.py` and
`MPI4NodalQuant.py` in place, per arm) -- and that adoption is what forced the
two API changes below.

NESTED RESOURCES (item 74). `resource` was a single directory NAME under
REPO_ROOT, which the perf resources are not: `testsys/perf/perf_case` and
`src/python/eqdyna` are both nested. It is now a relative PATH, and the
lockfile still goes BESIDE the guarded directory -- in that directory's
PARENT, named `.<basename>.lock`. For a single-component resource that is
exactly what it was (`test` -> `REPO_ROOT/.test.lock`, unchanged byte for
byte); for `src/python/eqdyna` it is `src/python/.eqdyna.lock`.

  The alternative -- leave `resource` a bare name and let each caller pass a
  deeper `repo_root` (`acquire(os.path.join(ROOT, 'src/python'), 'eqdyna')`)
  -- was rejected: it needs no change here, but it moves the path arithmetic
  into every caller under an argument named `repo_root` that is then not the
  repo root, and a caller that joins wrong puts the lockfile somewhere
  plausible and silently guards nothing. The path is validated in ONE place
  instead.

PER-TOOL CONSEQUENCE. What a collision COSTS differs by tool -- the e2e tools
lose a 1500 s cell to a rotated-away tree, `run_jaxmpi_ab` records a number
under the WRONG ARM LABEL (its own docstring: worse than no number) -- and a
refusal that describes the wrong failure teaches the reader to ignore it. So
the tail of the refusal is the caller's to supply, defaulting to the e2e text
unchanged. Any tail must still refuse in the rule-2 shape (no waiting, no
fallback); `testsys/regression/test_perf_tool_locks.py` asserts that on the
real refusals rather than trusting each caller to remember.
"""
import atexit
import errno
import fcntl
import os
import sys
import time

__all__ = ['RunTreeLocked', 'RunTreeLock', 'acquire', 'lock_path',
           'REFUSAL_HEADER', 'TAKEOVER_HEADER', 'E2E_CONSEQUENCE']

# The two strings testsys/regression/test_e2e_run_tree_lock.py asserts on.
# Named constants rather than prose buried in a format string, so the guard
# pins the CONTRACT and not an incidental wording.
REFUSAL_HEADER = 'refusing to start: another invocation already holds the lock on'
TAKEOVER_HEADER = 'NOTICE: took over a STALE lock on'


class RunTreeLocked(RuntimeError):
    """Another LIVE invocation holds the lock on this run tree."""


def lock_path(repo_root, resource):
    """The lockfile for REPO_ROOT/<resource>, in that directory's PARENT.

    BESIDE the guarded directory, never inside it: a lockfile under `test/`
    would be rotated away by the very rotation it exists to serialise, and one
    under `src/python/eqdyna/` would sit inside the package `run_jaxmpi_ab`
    stages arm files into.

    `resource` is a relative path with at least one component. A single
    component reduces to the original form exactly:
    `lock_path(root, 'test') == root + '/.test.lock'`.
    """
    parts = str(resource).split(os.sep) if resource else []
    if (not parts or os.path.isabs(resource)
            or any(p in ('', '.', '..') for p in parts)):
        raise ValueError(
            'resource must be a relative directory path under repo_root with '
            "no empty, '.' or '..' component, got %r" % (resource,))
    return os.path.join(repo_root, *(parts[:-1] + ['.%s.lock' % parts[-1]]))


def _record_text(argv):
    return ''.join('%s\t%s\n' % kv for kv in (
        ('pid', os.getpid()),
        ('started', time.strftime('%Y-%m-%d %H:%M:%S %z')),
        ('started_epoch', '%.3f' % time.time()),
        ('cwd', os.getcwd()),
        ('argv', ' '.join(sys.argv if argv is None else argv)),
    ))


def _read_record(path):
    """The holder's identity record, or None when there is not one to read.

    Never raises: a refusal must not turn into a traceback because the
    lockfile was unreadable, and it must not turn into a PASS either -- the
    caller renders None as UNKNOWN and refuses anyway.
    """
    try:
        with open(path, errors='replace') as fh:
            raw = fh.read()
    except OSError:
        return None
    rec = {}
    for line in raw.splitlines():
        if '\t' in line:
            key, val = line.split('\t', 1)
            rec[key] = val
    return rec or None


def _pid_alive(pid):
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True                      # exists, owned by another user
    except (TypeError, ValueError, OverflowError):
        return None                      # unparsable pid -- say so, do not guess
    return True


def _holder_lines(rec):
    if not rec:
        return [
            '  holder pid   : UNKNOWN',
            '  holder since : UNKNOWN',
            '                 The lockfile carries no identity record yet: the'
            ' holder took',
            '                 the lock microseconds ago and has not written'
            ' itself into it.',
            '                 The refusal stands regardless -- flock is the'
            ' lock, this',
            '                 record is only provenance.',
        ]
    return [
        '  holder pid   : %s' % rec.get('pid', 'UNKNOWN'),
        '  holder since : %s' % rec.get('started', 'UNKNOWN'),
        '  holder cwd   : %s' % rec.get('cwd', 'UNKNOWN'),
        '  holder cmd   : %s' % rec.get('argv', 'UNKNOWN'),
    ]


# The default tail: what a SECOND e2e invocation would have cost, and the
# rule-2 shape of the refusal. A caller whose collision costs something else
# passes its own tail to acquire(consequence=...); see the module docstring.
E2E_CONSEQUENCE = [
    'A second invocation in this checkout would rename that LIVE tree'
    ' out from under',
    'the holder mid-run (rule 21a, pathway item 70). On 2026-09-22 that'
    ' killed a',
    '1500.1 s cell with a FileNotFoundError, which then printed FAIL and'
    ' read as a',
    'solver regression in the sweep summary.',
    '',
    'NOT waiting. NOT rotating anyway. NOT falling back to a different'
    ' directory --',
    'each of those is the silent fallback rule 2 forbids. Run your sweep'
    ' in its own',
    'git worktree (rule 21a), or wait for the holder above to finish.']


def refusal_message(tree, path, rec, consequence=None):
    return '\n'.join(
        ['%s %s' % (REFUSAL_HEADER, tree)]
        + _holder_lines(rec)
        + ['  lockfile     : %s' % path,
           '']
        + list(E2E_CONSEQUENCE if consequence is None else consequence))


def takeover_message(tree, path, rec):
    pid = rec.get('pid', 'UNKNOWN')
    alive = _pid_alive(pid)
    if alive is None:
        why = ('the record names no parsable pid, so its liveness cannot be'
               ' stated; flock')
    elif alive:
        why = ('pid %s IS still alive but does NOT hold the flock -- it'
               ' released without' % pid)
    else:
        why = 'pid %s is no longer alive; flock' % pid
    tail = {
        None: ['                   is the authority here, and the kernel'
               ' handed us the lock.'],
        True: ['                   truncating its record, or that pid has been'
               ' reused since.',
               '                   flock is the authority here, and the kernel'
               ' handed us the lock.'],
        False: ['                   released it when that process died.'],
    }[alive]
    return '\n'.join(
        ['%s %s' % (TAKEOVER_HEADER, tree),
         '  previous holder: pid %s, started %s'
         % (pid, rec.get('started', 'UNKNOWN')),
         '  previous cmd   : %s' % rec.get('argv', 'UNKNOWN'),
         '  lockfile       : %s' % path,
         '  why it is safe : %s' % why]
        + tail
        + ['                   A LIVE holder keeps the kernel lock, so this'
           ' message cannot',
           '                   print while one exists.',
           '  That run left its evidence under %s and %s.prev (rule 8) -- read'
           ' it before' % (tree, tree),
           '  this invocation rotates them.'])


class RunTreeLock:
    """A held flock on REPO_ROOT/<resource>'s lockfile. Released at process
    exit (atexit) or explicitly; releasing twice is a no-op."""

    def __init__(self, fd, path, tree, takeover_notice=None):
        self.fd = fd
        self.path = path
        self.tree = tree
        self.takeover_notice = takeover_notice
        self.pid = os.getpid()
        self._released = False

    def release(self):
        if self._released:
            return
        self._released = True
        try:
            # Truncate FIRST, so a surviving record means "died without
            # releasing" and nothing else. See the module docstring.
            os.ftruncate(self.fd, 0)
            fcntl.flock(self.fd, fcntl.LOCK_UN)
        finally:
            os.close(self.fd)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.release()
        return False


def acquire(repo_root, resource, argv=None, announce=True, consequence=None):
    """Take the exclusive lock on REPO_ROOT/<resource>, or raise RunTreeLocked.

    Non-blocking by design: a second invocation must REFUSE, loudly and
    immediately, naming the holder. Waiting would be the wrong answer twice
    over -- it hides the collision from whoever caused it, and it parks a
    multi-hour sweep behind another multi-hour sweep with no one watching.

    `consequence` is the refusal's closing paragraph, as a list of lines:
    what THIS tool's collision destroys, and why the refusal does not
    degrade into waiting or into a second directory. Default: the e2e text.
    """
    path = lock_path(repo_root, resource)
    tree = os.path.join(repo_root, resource)
    fd = os.open(path, os.O_RDWR | os.O_CREAT | os.O_CLOEXEC, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        if exc.errno not in (errno.EACCES, errno.EAGAIN, errno.EWOULDBLOCK):
            os.close(fd)
            raise
        rec = _read_record(path)
        os.close(fd)
        raise RunTreeLocked(refusal_message(tree, path, rec, consequence))
    except BaseException:
        os.close(fd)
        raise

    notice = None
    stale = _read_record(path)
    if stale:
        notice = takeover_message(tree, path, stale)
        if announce:
            print(notice)
    os.ftruncate(fd, 0)
    os.lseek(fd, 0, os.SEEK_SET)
    os.write(fd, _record_text(argv).encode())
    os.fsync(fd)

    lock = RunTreeLock(fd, path, tree, notice)
    atexit.register(lock.release)
    return lock
