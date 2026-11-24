#! /usr/bin/env python3
"""
Regression guard: the Fortran driver refuses an unknown command-line argument
(board row 169 (a)).

eqdyna3d.f90 used to loop over command_argument_count() and act only on
--mesh-only, so a mistyped flag was ignored and a full solve ran. Now any other
argument calls abortRun(ERR_GENERIC) (exit 1, the argument named). Two checks:

  1. static: the argument loop in eqdyna3d.f90 has an abortRun branch for an
     argument that is not --mesh-only (red if the stop is removed);
  2. behavioural, both ways: with bin/eqdyna built, `eqdyna --bogus` on a real
     tpv8 case exits 1 and names the argument, and `eqdyna --mesh-only` on the
     same case still exits 0. Skipped, loudly, only when no binary is built.

Exits non-zero on any failure (rule 2).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'regression'))
SRC = os.path.join(ROOT, 'src', 'fortran', 'eqdyna3d.f90')
BIN = os.path.join(ROOT, 'bin', 'eqdyna')


def static_ok(text):
    """True if the do-iArg loop contains an abortRun for a non---mesh-only argument."""
    m = re.search(r'do iArg = 1, command_argument_count\(\)(.*?)enddo', text, re.S)
    if not m:
        return False
    body = m.group(1)
    return "'--mesh-only'" in body and re.search(r'\belse\b', body) is not None and 'abortRun' in body


def main():
    print('Regression guard: unknown command-line argument is refused')
    text = open(SRC, errors='replace').read()
    bad = 0
    # both ways on the static check: the real source passes, a copy with the branch removed fails
    stripped = re.sub(r'(do iArg = 1, command_argument_count\(\).*?)\belse\b.*?endif',
                      r'\1endif', text, count=1, flags=re.S)
    ok_real, ok_stripped = static_ok(text), static_ok(stripped)
    print('  static: real source %s, branch-removed copy %s (%d chars read)'
          % (ok_real, ok_stripped, len(text)))
    if not ok_real or ok_stripped:
        print('FAIL: static check does not distinguish the stop from its absence')
        bad += 1
    if not os.path.exists(BIN):
        print('SKIP behavioural: %s not built (run ./install-eqdyna.sh -m ubuntu)' % BIN)
        return 1 if bad else 0
    import test_mesh_screening as T
    tmp = tempfile.mkdtemp()
    try:
        d = os.path.join(tmp, 'c')
        T.make_case('test.tpv8', d)
        for arg, want_zero in (('--bogus', False), ('--mesh-only', True)):
            r = subprocess.run(['timeout', '-s', 'KILL', '600', 'mpirun', '-np', '1', BIN, arg],
                               cwd=d, env=T.env(), capture_output=True, text=True)
            out = r.stdout + r.stderr
            named = arg in out and 'unknown command-line argument' in out
            print('  %s rc=%d named=%s' % (arg, r.returncode, named))
            if want_zero and r.returncode != 0:
                print('FAIL: --mesh-only rc %d' % r.returncode)
                bad += 1
            if not want_zero and (r.returncode == 0 or not named):
                print('FAIL: %s not refused (rc %d, named %s)' % (arg, r.returncode, named))
                bad += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('PASS' if not bad else 'FAILED')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
