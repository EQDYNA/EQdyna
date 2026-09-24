#! /usr/bin/env python3
"""
Regression guard (pathway_forward.md item 109): a python e2e cell's stdout
and stderr must be KEPT in its case directory, not only inherited to the
sweep's console.

THE INCIDENT: `test.tpv36 x python-jax` exited 1 once (mira's everyday sweep
on 51eb44b, 2026-09-24 01:27-01:33) after all 464 steps, with frt.txt0 and
profile.rank0.json written. `run_standalone` used a bare `subprocess.call`
with inherited streams, so nothing reached disk and the traceback that
would have explained the exit was lost. It has not reproduced since.

THE FIX: `run_standalone` runs the solver through `_call_kept`, which tees
both streams live to the console AND to `<case_dir>/eqdyna.<backend>.stdout`
/ `.stderr`, and names both paths (plus the stderr tail) in the RuntimeError
it raises on a nonzero exit.

Behaviour, not source text (rule 10a), in two checks:
  (a) `_call_kept` on a child that writes to both streams and exits 3:
      rc is 3, each file holds exactly its stream, and the console (this
      guard runs the helper in a subprocess to capture it) saw both.
  (c) a kept-file write that fails (stdout kept on /dev/full, 1 MB child
      output) raises after the child exits cleanly instead of hanging on a
      full pipe (victor-reyes, PR #11 audit);
  (b) `run_standalone` on an EMPTY case dir, the real `python3 -m eqdyna`
      path: the solver fails reading its inputs, the RuntimeError names the
      kept stderr path, and that file exists and is non-empty.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

CHILD = ("import sys; sys.stdout.write('OUT-LINE\\n'); sys.stdout.flush(); "
         "sys.stderr.write('ERR-LINE\\n'); sys.exit(3)")

HELPER = r'''
import sys
sys.path.insert(0, {e2e!r})
import run_e2e
rc, o, e = run_e2e._call_kept([sys.executable, '-c', {child!r}], {tmp!r}, None,
                              {prefix!r})
print('RC=%d' % rc)
'''

HELPER_FULL = r'''
import builtins, os, sys
sys.path.insert(0, {e2e!r})
import run_e2e
real_open = builtins.open
def fake_open(path, mode='r', *a, **k):
    if str(path).endswith('.stdout') and 'w' in mode:
        return real_open('/dev/full', mode, *a, **k)
    return real_open(path, mode, *a, **k)
run_e2e.open = fake_open
try:
    run_e2e._call_kept([sys.executable, '-c',
                        "import sys; sys.stdout.write('x' * (1 << 20))"],
                       {tmp!r}, None, os.path.join({tmp!r}, 'full'))
    print('NO RAISE')
except RuntimeError as e:
    print('RAISED', e)
'''


def main():
    fails = []
    tmp = tempfile.mkdtemp(prefix='testCellLogsKept.')
    try:
        prefix = os.path.join(tmp, 'eqdyna.python-jax')
        r = subprocess.run(
            [sys.executable, '-c', HELPER.format(
                e2e=os.path.join(ROOT, 'testsys', 'e2e'), child=CHILD,
                tmp=tmp, prefix=prefix)],
            capture_output=True, text=True, timeout=120)
        if 'RC=3' not in r.stdout:
            fails.append('(a) child exit 3 not returned; helper stdout=%r stderr=%r'
                         % (r.stdout, r.stderr))
        try:
            kept_out = open(prefix + '.stdout').read()
            kept_err = open(prefix + '.stderr').read()
        except OSError as e:
            fails.append('(a) kept file missing: %s' % e)
            kept_out = kept_err = None
        if kept_out is not None:
            if kept_out != 'OUT-LINE\n':
                fails.append('(a) kept stdout is %r, expected exactly the child stdout' % kept_out)
            if kept_err != 'ERR-LINE\n':
                fails.append('(a) kept stderr is %r, expected exactly the child stderr' % kept_err)
        if 'OUT-LINE' not in r.stdout or 'ERR-LINE' not in r.stderr:
            fails.append('(a) console did not receive both streams: stdout=%r stderr=%r'
                         % (r.stdout, r.stderr))

        # (c) a write failure in a pump must RAISE after the child exits, not
        # hang: stdout is kept on /dev/full (every write fails with ENOSPC)
        # while the child writes 1 MB, far past a pipe buffer.
        try:
            r3 = subprocess.run(
                [sys.executable, '-c', HELPER_FULL.format(
                    e2e=os.path.join(ROOT, 'testsys', 'e2e'), tmp=tmp)],
                capture_output=True, text=True, timeout=60)
        except subprocess.TimeoutExpired:
            r3 = None
            fails.append('(c) a failing kept-file write HUNG the helper (>60 s): '
                         'a dead pump stopped draining and the child blocked')
        if r3 is not None and ('RAISED' not in r3.stdout or 'rc=0' not in r3.stdout):
            fails.append('(c) a failing kept-file write did not raise after a '
                         'clean child exit: stdout=%r stderr=%r'
                         % (r3.stdout[-300:], r3.stderr[-300:]))

        import run_e2e
        case_dir = os.path.join(tmp, 'empty_case')
        os.makedirs(case_dir)
        msg = None
        try:
            run_e2e.run_standalone(case_dir, 'python-jax')
            fails.append('(b) run_standalone on an empty case dir did NOT raise')
        except RuntimeError as e:
            msg = str(e)
        err_path = os.path.join(case_dir, 'eqdyna.python-jax.stderr')
        if msg is not None and err_path not in msg:
            fails.append('(b) RuntimeError does not name %s: %s' % (err_path, msg[:300]))
        if not os.path.isfile(err_path) or os.path.getsize(err_path) == 0:
            fails.append('(b) %s missing or empty after a failed cell' % err_path)
        err_bytes = os.path.getsize(err_path) if os.path.isfile(err_path) else -1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_e2e_python_cell_logs_kept')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_e2e_python_cell_logs_kept (tee rc=3, stdout/stderr kept exactly; '
          'failed real cell kept %d stderr bytes and named the path)' % err_bytes)
    return 0


if __name__ == '__main__':
    sys.exit(main())
