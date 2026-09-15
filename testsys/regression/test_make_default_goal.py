#! /usr/bin/env python3
"""
Regression guard: a bare `make` must build the binary (rules 2, 3, 10).

Guards the v5.6.0 CI-red incident. `src/makefile` used to open with

    clean:eqdyna
    	rm -r *.o

-- `clean` with no space before its prerequisite, so `clean` DEPENDED on
`eqdyna`. Two consequences, one of them load-bearing by accident:

  * `make clean` built the binary before deleting the objects (the bug that
    was reported), and
  * `clean` is the FIRST target in the file, and make builds the first target
    when invoked with no arguments -- so a bare `make` ran `clean`, which
    pulled in `eqdyna` as a prerequisite and therefore built it.

`install-eqdyna.sh:99` runs a bare `make`. Fixing `clean` removed the accident
that made that work: `make` then ran a clean that built nothing, `mv
src/eqdyna bin/` failed, and every tier that needs a binary went red -- while
the local gate stayed green, because it had always invoked `make eqdyna`
explicitly and never the bare form the install script uses.

The fix is `.DEFAULT_GOAL := eqdyna`. This test asserts the behaviour, not the
spelling: it runs a real bare `make` in a scratch copy of src/ and requires a
binary to appear, so any future reordering of targets that reintroduces the
problem fails here regardless of how it is written.

Cheap-ish (rule 9): one full compile of src/ in a temp dir, ~20 s. That is the
price of testing the actual build rather than grepping the makefile, and this
defect is invisible to a grep.
Exits non-zero on any failure (rule 2).
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, 'src')


def _copySrc(dest):
    """Copy the sources and makefile only -- never build products (rule 12)."""
    os.makedirs(dest, exist_ok=True)
    for name in sorted(os.listdir(SRC)):
        if name.endswith('.f90') or name == 'makefile':
            shutil.copy2(os.path.join(SRC, name), os.path.join(dest, name))
    return dest


def _make(cwd, args):
    env = dict(os.environ)
    env.setdefault('MACHINE', 'ubuntu')
    return subprocess.run(['make'] + args, cwd=cwd, env=env,
                          capture_output=True, text=True, timeout=900)


def main():
    print('Regression guard: a bare `make` must build the binary')
    problems = []

    with tempfile.TemporaryDirectory(prefix='testMakeDefaultGoal.') as tmp:
        work = _copySrc(os.path.join(tmp, 'src'))

        # The whole point: no target named on the command line, exactly as
        # install-eqdyna.sh invokes it.
        r = _make(work, [])
        binary = os.path.join(work, 'eqdyna')
        if r.returncode != 0:
            problems.append('bare `make` exited %d' % r.returncode)
            print((r.stdout + r.stderr)[-1500:])
        if not os.path.exists(binary):
            problems.append(
                'bare `make` completed but produced no eqdyna binary -- '
                'install-eqdyna.sh runs exactly this and its `mv src/eqdyna '
                'bin/` will fail. Check .DEFAULT_GOAL in src/makefile.')
        else:
            print('  bare `make` -> eqdyna (%d bytes)'
                  % os.path.getsize(binary))

        # And clean must NOT build: that was the original defect, and a fix
        # that restores the dependency would make this test pass for the
        # wrong reason.
        work2 = _copySrc(os.path.join(tmp, 'src2'))
        r2 = _make(work2, ['clean'])
        if r2.returncode != 0:
            problems.append('`make clean` exited %d' % r2.returncode)
        if os.path.exists(os.path.join(work2, 'eqdyna')):
            problems.append(
                '`make clean` built the binary -- clean must not depend on '
                'eqdyna (that dependency is what hid the default-goal bug)')
        else:
            print('  `make clean` -> no binary, as intended')

    if problems:
        print('\nFAIL:')
        for p in problems:
            print('  ' + p)
        return 1
    print('\nPASS: default goal builds, clean does not')
    return 0


if __name__ == '__main__':
    sys.exit(main())
