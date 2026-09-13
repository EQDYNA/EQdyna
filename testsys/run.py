#! /usr/bin/env python3
"""
Single entry point for EQdyna's tiered test system (PROJECT_RULES.md rule 3).

    python3 testsys/run.py unit          # fast pure-python unit tests (pytest, no MPI/Fortran)
    python3 testsys/run.py regression    # one guard per past incident (rule 10)
    python3 testsys/run.py e2e           # full pipeline vs test.reference.results/ (rule 7)
    python3 testsys/run.py parity        # Python-port (NumPy/JAX) vs Fortran serial oracle (item 14)
    python3 testsys/run.py perf          # pinned single-core Fortran/NumPy/JAX timing, ratio-guarded
    python3 testsys/run.py all           # unit + regression + e2e, in that order (default; parity/perf are opt-in, not in "all" -- they need a Fortran build + fixtures a fresh checkout doesn't have yet)

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


def run_parity():
    print('\n==== testsys: parity ====')
    rc = subprocess.call([sys.executable, os.path.join(TESTSYS, 'parity', 'run_parity.py')],
                          cwd=REPO_ROOT)
    if rc != 0:
        return rc
    # Standalone (no-Fortran-in-the-loop) meshgen port, milestones 1+2 --
    # same fixtures, same regeneration workflow as run_parity.py above.
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'parity', 'test_standalone_meshgen.py')],
        cwd=REPO_ROOT)


def run_perf():
    print('\n==== testsys: perf ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_perf.py')],
                            cwd=REPO_ROOT)


RUNNERS = {'unit': run_unit, 'regression': run_regression, 'e2e': run_e2e,
           'parity': run_parity, 'perf': run_perf}
# 'all' stays unit+regression+e2e only (TIERS below) -- parity/perf require a
# Fortran build and generated fixtures/baseline that a fresh checkout does
# not have; they are opt-in tiers, invoked by name, not swept into 'all'.
OPTIONAL_TIERS = ('parity', 'perf')


def main(argv):
    which = argv[1] if len(argv) > 1 else 'all'
    if which not in TIERS + OPTIONAL_TIERS + ('all',):
        print(f'usage: python3 testsys/run.py [{"|".join(TIERS + OPTIONAL_TIERS)}|all]')
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
