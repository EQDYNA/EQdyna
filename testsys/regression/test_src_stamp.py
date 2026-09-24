#! /usr/bin/env python3
"""
Regression guard: the gates refuse a `bin/eqdyna` built from different
source (owner-approved 2026-09-24).

THE INCIDENT: a `bin/eqdyna` built on 2026-09-22, older than src/fortran,
made test_profile_env_strict fail locally for no real reason. A stale binary
could just as easily PASS on Fortran that no longer exists, and nothing
compared the binary with the tree.

THE FIX: `scripts/src_hash.py` (the ONE implementation) hashes
src/fortran/*.f90 + makefile from the working tree. The makefile embeds it
(srcStamp.inc -> eqdyna3d.f90's srcStampMod), the banner prints it, and
`testsys/run.py` (regression tier) and `testsys/e2e/run_e2e.py` (before any
fortran cell) refuse a binary whose stamp is missing or different.

Behaviour, not source text (rule 10a), both ways:
  A. check_binary on a COPY of src/fortran with a synthetic binary carrying
     the copy's stamp: fresh -> ok; a .f90 edited without rebuilding -> refused,
     naming both stamps and the fix; restored -> ok; the makefile edited ->
     refused; restored -> ok; a new .f90 -> refused; a binary with no stamp ->
     refused; no binary -> refused.
  B. `run.py regression` in a scratch tree whose bin/eqdyna is stale exits 1
     with REFUSED and never starts the regression tier.
  C. `run_e2e.py` with EQDYNA_E2E_BIN pointing at a stale synthetic binary
     exits 1 with REFUSED and launches no cell.
  D. the REAL bin/eqdyna carries this tree's stamp -- the makefile really
     embeds it (CI's unit-regression shards download the binary the build
     job just made from the same checkout). Absent binary = FAIL, not skip.
"""
import glob
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_src_hash(path=os.path.join(ROOT, 'scripts', 'src_hash.py')):
    spec = importlib.util.spec_from_file_location('src_hash_under_test', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fake_binary(path, stamp):
    with open(path, 'wb') as fh:
        fh.write(b'\x7fELF junk ')
        if stamp is not None:
            fh.write(b'EQDYNA_SRC_STAMP=' + stamp.encode())
        fh.write(b' more junk\n')


def check_a(sh, tmp, fails):
    fsrc = os.path.join(tmp, 'src', 'fortran')
    os.makedirs(fsrc)
    for p in glob.glob(os.path.join(ROOT, 'src', 'fortran', '*.f90')) + \
            [os.path.join(ROOT, 'src', 'fortran', 'makefile')]:
        shutil.copy(p, fsrc)
    binary = os.path.join(tmp, 'eqdyna')
    fake_binary(binary, sh.compute(fsrc))
    f90 = os.path.join(fsrc, 'faulting.f90')
    mk = os.path.join(fsrc, 'makefile')
    orig_f90, orig_mk = open(f90).read(), open(mk).read()

    def expect(label, want_ok, needle=None):
        ok, msg = sh.check_binary(binary, fsrc)
        if ok != want_ok or (needle and needle not in msg):
            fails.append('A %s: ok=%s (want %s), msg=%r' % (label, ok, want_ok, msg))

    expect('fresh', True)
    open(f90, 'a').write('! unbuilt edit\n')
    expect('.f90 edited, not rebuilt', False, 'was built from different source')
    expect('.f90 edited, message names the fix', False, './install-eqdyna.sh')
    open(f90, 'w').write(orig_f90)
    expect('.f90 restored', True)
    open(mk, 'a').write('# unbuilt edit\n')
    expect('makefile edited, not rebuilt', False, 'was built from different source')
    open(mk, 'w').write(orig_mk)
    expect('makefile restored', True)
    extra = os.path.join(fsrc, 'zzNew.f90')
    open(extra, 'w').write('! new file\n')
    expect('new .f90 added', False, 'was built from different source')
    os.remove(extra)
    fake_binary(binary, None)
    expect('binary carries no stamp', False, 'no source stamp')
    os.remove(binary)
    expect('no binary at all', False, 'does not exist')
    return fsrc


def check_b(tmp, fails):
    """`run.py regression` in a scratch tree with a stale bin/eqdyna."""
    t = os.path.join(tmp, 'tree_b')
    os.makedirs(os.path.join(t, 'testsys'))
    os.makedirs(os.path.join(t, 'scripts'))
    os.makedirs(os.path.join(t, 'bin'))
    shutil.copy(os.path.join(ROOT, 'testsys', 'run.py'), os.path.join(t, 'testsys'))
    shutil.copy(os.path.join(ROOT, 'scripts', 'src_hash.py'), os.path.join(t, 'scripts'))
    shutil.copytree(os.path.join(tmp, 'src'), os.path.join(t, 'src'))
    fake_binary(os.path.join(t, 'bin', 'eqdyna'), '0' * 64)
    r = subprocess.run([sys.executable, os.path.join(t, 'testsys', 'run.py'), 'regression'],
                       cwd=t, capture_output=True, text=True, timeout=60)
    out = r.stdout + r.stderr
    if r.returncode != 1 or 'REFUSED' not in out or '==== testsys: regression' in out:
        fails.append('B run.py regression with a stale bin/eqdyna: rc=%d, want 1 with '
                     'REFUSED and no regression tier started; output=%r'
                     % (r.returncode, out[-600:]))


def check_c(sh, tmp, fails):
    """run_e2e.py refuses a stale EQDYNA_E2E_BIN before any cell."""
    stale = os.path.join(tmp, 'stale_eqdyna')
    fake_binary(stale, '0' * 64)
    env = dict(os.environ, EQDYNA_E2E_BIN=stale)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py'),
                        '--cases', 'test.tpv8', '--backends', 'fortran'],
                       cwd=ROOT, env=env, capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    if r.returncode != 1 or 'REFUSED' not in out or 'was built from different source' not in out \
            or '-- cell:' in out:
        fails.append('C run_e2e.py with a stale EQDYNA_E2E_BIN: rc=%d, want 1 with '
                     'REFUSED before any cell; output tail=%r' % (r.returncode, out[-600:]))


def check_d(sh, fails):
    binary = os.path.join(ROOT, 'bin', 'eqdyna')
    ok, msg = sh.check_binary(binary)
    if not ok:
        fails.append('D the real bin/eqdyna: %s' % msg)
    return msg


def main():
    sh = load_src_hash()
    fails = []
    tmp = tempfile.mkdtemp(prefix='testSrcStamp.')
    try:
        check_a(sh, tmp, fails)
        check_b(tmp, fails)
        check_c(sh, tmp, fails)
        real = check_d(sh, fails)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if fails:
        print('FAIL test_src_stamp')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_src_stamp (A: 9 stamp outcomes incl. .f90/makefile edits '
          'refused and restores accepted; B: run.py refuses a stale binary before '
          'the regression tier; C: run_e2e refuses a stale EQDYNA_E2E_BIN before any '
          'cell; D: %s)' % real)
    return 0


if __name__ == '__main__':
    sys.exit(main())
