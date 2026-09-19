#! /usr/bin/env python3
"""
Regression guard for pathway_forward.md item 44 (P4).

Commit 45446d0 de-duplicated the test.tpv36 / test.tpv37 compsets by turning
case_input/test.tpv37/tpv36_37_common.py into a symlink (git mode 120000)
to case_input/test.tpv36/tpv36_37_common.py -- the repo's FIRST tracked
symlink. Item 44 records two risks it introduced, both real and previously
untested:

  1. On a checkout with `core.symlinks=false`, git materialises a tracked
     symlink as a plain text file whose CONTENT is the target path string,
     not the target file's content. `case.setup` then dies importing it,
     with a message that does not say why.
  2. test.tpv37 can no longer survive deletion or rename of
     case_input/test.tpv36/: a dangling symlink is a far less legible
     failure than a missing parameter, and nothing tested the coupling
     directly before this file.

This test enumerates tracked symlinks from git itself (`git ls-files -s`,
mode 120000) rather than a hardcoded list, so it keeps working the moment a
second symlink is added, and so that "git currently tracks zero symlinks"
is itself reported rather than silently passing.

For every tracked symlink it asserts, in order:
  a. os.path.islink(path) is True -- this is exactly what a
     core.symlinks=false checkout gets wrong, and failing here is the point:
     name the file and say plainly the checkout cannot represent symlinks.
  b. the symlink resolves (os.path.realpath) to a path that os.path.exists
     -- this is risk 2, the moment case_input/test.tpv36/ is deleted or
     renamed.
  c. the resolved target stays inside the repository root (no ../ escape
     above it).
  d. if the symlink's own name ends in .py, the resolved target is
     importable as a real Python module (not just syntax-valid) with
     scripts/ on sys.path, matching how case.setup actually loads it --
     this is what proves the case would still set up.

No skipif, no try/except that swallows: any of the above missing or wrong
is a hard FAIL, never a silent pass. Pure git + filesystem + one import;
well under 1 s.
"""
import importlib.util
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPTS = os.path.join(ROOT, 'scripts')


def tracked_symlinks():
    """[relpath, ...] for every git-tracked entry with mode 120000.

    Reads git's own index, not a hardcoded list, so this guard keeps
    working the moment a second symlink lands.
    """
    r = subprocess.run(['git', 'ls-files', '-s'], cwd=ROOT,
                       capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError('git ls-files -s failed (exit %d): %s'
                            % (r.returncode, r.stderr))
    out = []
    for line in r.stdout.splitlines():
        # format: <mode> <blob-sha> <stage>\t<path>
        meta, _, path = line.partition('\t')
        mode = meta.split()[0]
        if mode == '120000':
            out.append(path)
    return out


def check_one(relpath, fails):
    path = os.path.join(ROOT, relpath)

    if not os.path.islink(path):
        fails.append(
            f'{relpath}: git tracks this as a symlink (mode 120000) but the '
            f'working-tree entry is NOT a symlink (os.path.islink is False). '
            f'This checkout cannot represent symlinks (core.symlinks=false?) '
            f'-- git has materialised it as a plain file containing the '
            f'target path string instead of the target\'s content.')
        return  # every further check on a non-symlink is meaningless

    resolved = os.path.realpath(path)

    if not os.path.exists(resolved):
        fails.append(
            f'{relpath}: symlink target does not exist (resolved to '
            f'{resolved!r}). The file or directory it points at was '
            f'deleted or renamed.')
        return

    # Containment: resolved path must stay under ROOT.
    root_real = os.path.realpath(ROOT)
    if os.path.commonpath([root_real, resolved]) != root_real:
        fails.append(
            f'{relpath}: symlink target {resolved!r} resolves OUTSIDE the '
            f'repository root {root_real!r}.')
        return

    if relpath.endswith('.py'):
        modname = '_symlink_integrity_check_' + relpath.replace('/', '_').replace('.py', '')
        spec = importlib.util.spec_from_file_location(modname, path)
        if spec is None or spec.loader is None:
            fails.append(f'{relpath}: could not build an import spec for it')
            return
        module = importlib.util.module_from_spec(spec)
        sys.path.insert(0, SCRIPTS)
        try:
            spec.loader.exec_module(module)
        except Exception as exc:  # noqa: BLE001 -- reported as a FAIL, not swallowed
            fails.append(
                f'{relpath}: resolved target {resolved!r} is not importable '
                f'as Python ({type(exc).__name__}: {exc}) -- case.setup would '
                f'die on this import too.')
        finally:
            sys.path.remove(SCRIPTS)
            sys.modules.pop(modname, None)


def main():
    try:
        links = tracked_symlinks()
    except RuntimeError as exc:
        print('FAIL test_symlink_integrity')
        print(' -', exc)
        return 1

    if not links:
        print('FAIL test_symlink_integrity')
        print(' - git tracks ZERO symlinks (mode 120000). As of '
              'pathway_forward.md item 44 it should track exactly 1 '
              '(case_input/test.tpv37/tpv36_37_common.py). Either that '
              'symlink was replaced by a real file (fine, but then this '
              'guard has nothing left to protect and should be removed '
              'deliberately) or `git ls-files -s` / mode detection broke.')
        return 1

    fails = []
    for relpath in links:
        check_one(relpath, fails)

    if fails:
        print('FAIL test_symlink_integrity')
        for f in fails:
            print(' -', f)
        return 1

    print(f'SUCCESS test_symlink_integrity: {len(links)} tracked symlink(s) '
          f'({", ".join(links)}) are real symlinks, resolve to an existing '
          f'in-repo target, and (where .py) import cleanly')
    return 0


if __name__ == '__main__':
    sys.exit(main())
