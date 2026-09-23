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

REUSE. Nothing here is e2e-specific: `resource` is any single directory name
under REPO_ROOT. The other two tools in this defect class --
`testsys/perf/run_perf.py:172` (rebuilds `testsys/perf_case/<name>`) and
`testsys/perf/run_jaxmpi_ab.py:88` (deletes
`src/python/eqdyna/__pycache__`) -- can adopt it unchanged once someone is
free to gate them; they are out of scope for item 70 as briefed.
"""
import atexit
import errno
import fcntl
import os
import sys
import time

__all__ = ['RunTreeLocked', 'RunTreeLock', 'acquire', 'lock_path',
           'REFUSAL_HEADER', 'TAKEOVER_HEADER']

# The two strings testsys/regression/test_e2e_run_tree_lock.py asserts on.
# Named constants rather than prose buried in a format string, so the guard
# pins the CONTRACT and not an incidental wording.
REFUSAL_HEADER = 'refusing to start: another invocation already holds the lock on'
TAKEOVER_HEADER = 'NOTICE: took over a STALE lock on'


class RunTreeLocked(RuntimeError):
    """Another LIVE invocation holds the lock on this run tree."""


def lock_path(repo_root, resource):
    """The lockfile for REPO_ROOT/<resource>.

    BESIDE the tree, never inside it: a lockfile under `test/` would be
    rotated away by the very rotation it exists to serialise.
    """
    if not resource or os.sep in resource or resource in ('.', '..'):
        raise ValueError('resource must be a single directory name under '
                         'repo_root, got %r' % (resource,))
    return os.path.join(repo_root, '.%s.lock' % resource)


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


def refusal_message(tree, path, rec):
    return '\n'.join(
        ['%s %s' % (REFUSAL_HEADER, tree)]
        + _holder_lines(rec)
        + ['  lockfile     : %s' % path,
           '',
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
           'git worktree (rule 21a), or wait for the holder above to finish.'])


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


def acquire(repo_root, resource, argv=None, announce=True):
    """Take the exclusive lock on REPO_ROOT/<resource>, or raise RunTreeLocked.

    Non-blocking by design: a second invocation must REFUSE, loudly and
    immediately, naming the holder. Waiting would be the wrong answer twice
    over -- it hides the collision from whoever caused it, and it parks a
    multi-hour sweep behind another multi-hour sweep with no one watching.
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
        raise RunTreeLocked(refusal_message(tree, path, rec))
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
