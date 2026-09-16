#! /usr/bin/env python3
"""
Regression guard: every third-party import must be installed in CI (rules 2, 3).

Guards the v5.6.1 CI-red incident. `scripts/convertFaultGeometry` used
`scipy.interpolate` for its two interpolation paths. scipy was installed on the
development box and NOT in `.github/workflows/test.yml`, so the unit tier passed
locally and failed in CI with `ModuleNotFoundError: No module named 'scipy'`.

The general defect: a new third-party dependency can be added anywhere under
scripts/ or testsys/ and nothing connects it to CI's `pip install` line. The
local environment silently supplies it. This test is that connection.

It scans for imports and compares them against the packages CI installs, read
out of the workflow file itself so the two cannot drift. Standard-library
modules and first-party modules (anything importable from the repo) are
excluded.

An import may be absent from CI ONLY if it is lazy -- inside a function -- AND
guarded by a try/except that raises something actionable. That is the pattern
`convertFaultGeometry` and `scripts/lib.py` (imageio) both use: the feature
needing the dependency fails with a message naming it, while every other path
keeps working.

Cheap (rule 9): pure AST walk over the tree, no imports executed, well under 1 s.
Exits non-zero on any failure (rule 2).
"""
import ast
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')

# src/python is the SOLVER PACKAGE -- the code this repo actually ships and
# the code CI runs. It was never scanned: first_party() below special-cases it
# as a sys.path root (so the author knew it was there), but SCAN_DIRS stopped
# at scripts/ and testsys/, which made this guard's own docstring -- "every
# third-party import must be installed in CI" -- false for the one package
# that matters most. Harmless so far only because every import in it is numpy,
# jax or relative, and both are declared; that is luck, not coverage.
SCAN_DIRS = ('scripts', 'testsys', os.path.join('src', 'python'))

# pip name -> import name, where they differ.
PIP_TO_IMPORT = {
    'netCDF4': 'netCDF4',
    'scikit-learn': 'sklearn',
    'pillow': 'PIL',
    'pyyaml': 'yaml',
}


def ci_packages():
    """Import names CI installs, read from the workflow so it cannot drift."""
    text = open(WORKFLOW, errors='replace').read()
    names = set()
    for line in re.findall(r'pip install ([^\n]+)', text):
        for tok in line.split():
            if tok.startswith('-'):
                continue
            names.add(PIP_TO_IMPORT.get(tok, tok))
    # Always available in the CI image.
    names |= {'pip', 'setuptools'}
    return names


def first_party():
    """Modules importable from inside the repo -- not third-party.

    Walks the WHOLE tree, not just top-level directories: case directories put
    user_defined_params.py next to the script that imports it, tests import
    conftest from their own directory, and testsys modules import siblings
    (full_specs, elem_per_rank). All of those resolve at runtime via sys.path
    and none is a package CI could install.
    """
    names = set()
    for root, dirs, files in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__')]
        for entry in files:
            if entry.endswith('.py'):
                names.add(entry[:-3])
        for d in dirs:
            # A package by __init__.py, or a directory the repo puts on
            # sys.path itself (src/python/eqdyna is imported as `eqdyna` once
            # src/python is on the path -- it has no __init__.py).
            if (os.path.exists(os.path.join(root, d, '__init__.py'))
                    or os.path.relpath(root, ROOT) == os.path.join('src', 'python')):
                names.add(d)
    return names


def _is_guarded(node, tree):
    """True if this import sits inside a try/except (an actionable guard)."""
    for parent in ast.walk(tree):
        if isinstance(parent, ast.Try):
            for sub in ast.walk(parent):
                if sub is node:
                    return True
    return False


def _is_lazy(node, tree):
    """True if this import sits inside a function, not at module level."""
    for parent in ast.walk(tree):
        if isinstance(parent, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for sub in ast.walk(parent):
                if sub is node:
                    return True
    return False


def scan():
    ci = ci_packages()
    local = first_party()
    stdlib = set(getattr(sys, 'stdlib_module_names', ())) | {
        'os', 'sys', 're', 'math', 'glob', 'shutil', 'subprocess', 'tempfile',
        'types', 'argparse', 'importlib', 'json', 'textwrap', 'collections',
        'functools', 'itertools', 'pathlib', 'time', 'datetime', 'copy',
        'warnings', 'contextlib', 'hashlib', 'traceback', 'ast', 'csv',
    }
    offenders = []
    for d in SCAN_DIRS:
        for root, _, files in os.walk(os.path.join(ROOT, d)):
            # The parity/fixtures exemption that used to live here is gone with
            # the fixtures: they were frozen copies of scripts/ that had drifted
            # badly (lib.py stuck at 187 lines against the live 949), and the
            # only thing keeping them lint-clean was this exemption. Every
            # remaining .py under testsys/ is live code and is scanned.
            if '__pycache__' in root:
                continue
            for fn in sorted(files):
                path = os.path.join(root, fn)
                if not (fn.endswith('.py') or _looks_like_python(path)):
                    continue
                try:
                    tree = ast.parse(open(path, errors='replace').read())
                except SyntaxError:
                    continue
                for node in ast.walk(tree):
                    mods = []
                    if isinstance(node, ast.Import):
                        mods = [a.name.split('.')[0] for a in node.names]
                    elif isinstance(node, ast.ImportFrom) and node.level == 0:
                        mods = [(node.module or '').split('.')[0]]
                    for m in mods:
                        if not m or m in stdlib or m in local or m in ci:
                            continue
                        if _is_lazy(node, tree) and _is_guarded(node, tree):
                            continue    # documented escape hatch
                        offenders.append(
                            (os.path.relpath(path, ROOT), node.lineno, m,
                             _is_lazy(node, tree)))
    return offenders, ci


def _looks_like_python(path):
    if os.path.basename(path) in ('case.setup', 'create.newcase',
                                  'convertFaultGeometry',
                                  'generateFaultInterface',
                                  'plotRuptureDynamics'):
        return True
    return False


def main():
    print('Regression guard: third-party imports must be installed in CI')
    offenders, ci = scan()
    print('  CI installs: %s' % ', '.join(sorted(ci - {'pip', 'setuptools'})))
    if offenders:
        print('\nFAIL: %d import(s) CI does not install:' % len(offenders))
        for rel, lineno, mod, lazy in offenders:
            how = 'lazy but unguarded' if lazy else 'module-level'
            print('  %s:%d  imports %r  (%s)' % (rel, lineno, mod, how))
        print('\nEither add the package to the `pip install` line in'
              ' .github/workflows/test.yml,\n  or make the import lazy AND wrap'
              ' it in try/except that raises a message\n  naming the package and'
              ' what still works without it.')
        return 1
    print('\nPASS: no undeclared third-party imports')
    return 0


if __name__ == '__main__':
    sys.exit(main())
