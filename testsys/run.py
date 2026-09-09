#! /usr/bin/env python3
"""
Single entry point for EQdyna's tiered test system (PROJECT_RULES.md rule 3).

    python3 testsys/run.py unit          # fast pure-python unit tests (pytest, no MPI/Fortran)
    python3 testsys/run.py regression    # one guard per past incident (rule 10)
    python3 testsys/run.py e2e           # full pipeline vs test.reference.results/ (rule 7)
    python3 testsys/run.py all           # unit + regression + e2e, in that order (default)

Prints a per-test SUCCESS/FAIL line (from pytest or from each regression
script's own banner), a per-tier SUMMARY line, and exits non-zero if
anything in the requested scope failed. "The script ran" and "the script
passed" are always two different questions here (rule 3).
"""
import glob
import os
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
TIERS = ('unit', 'regression', 'e2e')


def run_unit():
    d = os.path.join(TESTSYS, 'unit')
    print('\n==== testsys: unit ====')
    return subprocess.call([sys.executable, '-m', 'pytest', '-v', d], cwd=REPO_ROOT)


def run_regression():
    """Each regression/test_*.py is a standalone script with its own
    SUCCESS/FAIL banner and sys.exit (matching the shape of the incident
    it guards) rather than pytest functions -- run each one and gate on
    its exit code."""
    print('\n==== testsys: regression ====')
    scripts = sorted(glob.glob(os.path.join(TESTSYS, 'regression', 'test_*.py')))
    if not scripts:
        print('regression: FAIL - no regression scripts found (misconfigured testsys/)')
        return 1
    overall = 0
    for script in scripts:
        name = os.path.basename(script)
        print(f'-- {name} --')
        rc = subprocess.call([sys.executable, script], cwd=REPO_ROOT)
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {name} (exit {rc})')
        overall = overall or rc
    return overall


def run_e2e():
    print('\n==== testsys: e2e ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py')],
                            cwd=REPO_ROOT)


RUNNERS = {'unit': run_unit, 'regression': run_regression, 'e2e': run_e2e}


def main(argv):
    which = argv[1] if len(argv) > 1 else 'all'
    if which not in TIERS + ('all',):
        print(f'usage: python3 testsys/run.py [{"|".join(TIERS)}|all]')
        return 2

    selected = TIERS if which == 'all' else (which,)
    results = {tier: RUNNERS[tier]() for tier in selected}

    print('\n==== testsys: SUMMARY ====')
    overall = 0
    for tier in selected:
        rc = results[tier]
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {tier} (exit {rc})')
        overall = overall or rc

    return 1 if overall else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
