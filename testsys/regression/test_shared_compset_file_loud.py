#! /usr/bin/env python3
"""
Regression guard for pathway_forward.md item 44 (owner ruling 2026-09-24):
create.newcase and case.setup fail LOUDLY AND SPECIFICALLY when test.tpv37's
shared tpv36_37_common.py (a symlink into case_input/test.tpv36/) is missing,
instead of relying on test_symlink_integrity.py to notice.

Three behaviours, each run as a real subprocess on a throwaway EQDYNAROOT:
  1. dangling symlink in the compset  -> create.newcase exits nonzero naming it
     (before the fix it skipped the file silently and exited 0);
  2. symlink stored as a text file (core.symlinks=false) -> create.newcase
     exits nonzero saying so;
  3. the file absent from a built case dir -> case.setup exits nonzero naming
     tpv36_37_common and create.newcase (before: a bare ModuleNotFoundError).
No skip paths: any setup problem is a FAIL.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMMON = 'tpv36_37_common.py'


def fake_root(tmp, link_mode):
    """A minimal EQDYNAROOT: scripts/ (files only) plus case_input/test.tpv37
    whose shared file is broken in the requested way."""
    r = os.path.join(tmp, 'root_' + link_mode)
    os.makedirs(os.path.join(r, 'scripts'))
    for f in os.listdir(os.path.join(ROOT, 'scripts')):
        p = os.path.join(ROOT, 'scripts', f)
        if os.path.isfile(p):
            shutil.copy(p, os.path.join(r, 'scripts'))
    src = os.path.join(ROOT, 'case_input', 'test.tpv37')
    dst = os.path.join(r, 'case_input', 'test.tpv37')
    os.makedirs(dst)
    for f in os.listdir(src):
        if f != COMMON and os.path.isfile(os.path.join(src, f)):
            shutil.copy(os.path.join(src, f), dst)
    target = '../test.tpv36/' + COMMON
    if link_mode == 'dangling':
        os.symlink(target, os.path.join(dst, COMMON))      # no test.tpv36 dir
    elif link_mode == 'as_text':
        os.makedirs(os.path.join(r, 'case_input', 'test.tpv36'))
        shutil.copy(os.path.join(ROOT, 'case_input', 'test.tpv36', COMMON),
                    os.path.join(r, 'case_input', 'test.tpv36'))
        with open(os.path.join(dst, COMMON), 'w') as fh:
            fh.write(target)
    return r


def newcase(root, case):
    env = dict(os.environ, EQDYNAROOT=root)
    return subprocess.run([sys.executable, os.path.join(root, 'scripts', 'create.newcase'),
                           case, 'test.tpv37'], env=env, capture_output=True, text=True, timeout=60)


def main():
    fails = []
    tmp = tempfile.mkdtemp(prefix='testSharedCompset.')
    try:
        for mode, want in (('dangling', 'does not exist'), ('as_text', 'core.symlinks')):
            r = newcase(fake_root(tmp, mode), os.path.join(tmp, 'case_' + mode))
            out = r.stdout + r.stderr
            if r.returncode == 0 or COMMON not in out or want not in out:
                fails.append('create.newcase with a %s %s: exit %d, output %r'
                             % (mode, COMMON, r.returncode, out[-400:]))
        case = os.path.join(tmp, 'case_ok')
        r = newcase(ROOT, case)
        if r.returncode != 0 or not os.path.isfile(os.path.join(case, COMMON)):
            fails.append('create.newcase on the real tree failed or did not copy %s: %r'
                         % (COMMON, (r.stdout + r.stderr)[-400:]))
        else:
            os.remove(os.path.join(case, COMMON))
            s = subprocess.run([sys.executable, 'case.setup'], cwd=case,
                               capture_output=True, text=True, timeout=120)
            out = s.stdout + s.stderr
            if s.returncode == 0 or 'tpv36_37_common' not in out or 'create.newcase' not in out:
                fails.append('case.setup without %s: exit %d, output %r'
                             % (COMMON, s.returncode, out[-400:]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if fails:
        print('FAIL test_shared_compset_file_loud')
        for f in fails:
            print(' -', f)
        return 1
    print('SUCCESS test_shared_compset_file_loud: 3 of 3 broken-shared-file '
          'shapes refused by name (dangling symlink, symlink-as-text, missing '
          'from case dir)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
