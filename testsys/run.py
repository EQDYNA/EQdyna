#! /usr/bin/env python3
"""
Single entry point for EQdyna's tiered test system (PROJECT_RULES.md rule 3).

    python3 testsys/run.py unit          # fast pure-python unit tests (pytest, no MPI/Fortran)
    python3 testsys/run.py regression    # one guard per past incident (rule 10)
    python3 testsys/run.py e2e           # full pipeline vs test.reference.results/ (rule 7)
    python3 testsys/run.py e2e-full      # SCEC cases at spec dx/term, 16 ranks, report-only (opt-in; hours)
    python3 testsys/run.py parity        # Python-port (NumPy/JAX) vs Fortran serial oracle (item 14)
    python3 testsys/run.py accept        # standalone (no-Fortran-in-the-loop) solver vs committed 4-rank references
    python3 testsys/run.py perf          # pinned single-core Fortran/NumPy/JAX timing, ratio-guarded
    python3 testsys/run.py all           # unit + regression + e2e, in that order (default; parity/accept/perf are opt-in, not in "all" -- they need a Fortran build + fixtures a fresh checkout doesn't have yet)

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


def run_e2e_full():
    print('\n==== testsys: e2e-full ====')
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py'), '--full'],
        cwd=REPO_ROOT)


def run_parity():
    print('\n==== testsys: parity ====')
    rc = subprocess.call([sys.executable, os.path.join(TESTSYS, 'parity', 'run_parity.py')],
                          cwd=REPO_ROOT)
    if rc != 0:
        return rc
    # Standalone (no-Fortran-in-the-loop) meshgen port, milestones 1+2 --
    # same fixtures, same regeneration workflow as run_parity.py above.
    rc = subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'parity', 'test_standalone_meshgen.py')],
        cwd=REPO_ROOT)
    if rc != 0:
        return rc
    # The setup builders above are vectorized; each keeps its original
    # verbatim scalar loop as an oracle. This asserts the two are BYTE
    # identical on the same real mesh -- it catches vectorization drift,
    # which the Fortran-anchored check above tolerates inside its tolerance.
    rc = subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'parity', 'test_standalone_setup_vectorized.py')],
        cwd=REPO_ROOT)
    if rc != 0:
        return rc
    return run_functional_neutrality_check()


def run_functional_neutrality_check():
    """Cheap re-check (a byte compare, no Fortran re-run) of the
    functional-neutrality fixture make_fixtures.py generates: frt.txt0 from
    the default (pydump_noop-linked) binary must be byte-identical to
    frt.txt0 from the PYDUMP=1 (pydump-linked) binary, on the same case.
    Missing files (fixtures not regenerated since this check was added) is
    a loud FAIL, not a silent skip."""
    fixture_case = os.path.join(TESTSYS, 'parity', 'fixtures', 'test_tpv8_serial')
    default_frt = os.path.join(fixture_case, 'frt.txt0.default-build')
    pydump_frt = os.path.join(fixture_case, 'frt.txt0.pydump-build')
    for p in (default_frt, pydump_frt):
        if not os.path.isfile(p):
            print(f'FAIL functional-neutrality check: missing {p} -- '
                  f're-run testsys/parity/make_fixtures.py')
            return 1
    import filecmp
    if not filecmp.cmp(default_frt, pydump_frt, shallow=False):
        print('FAIL functional-neutrality check: frt.txt0 differs between the default and '
              'PYDUMP=1 builds')
        return 1
    print('SUCCESS functional-neutrality check (frt.txt0 byte-identical, default vs PYDUMP=1 build)')
    return 0


def run_accept():
    print('\n==== testsys: accept ====')
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'parity', 'test_standalone_acceptance.py')],
        cwd=REPO_ROOT)


def run_scaling():
    print('\n==== testsys: scaling ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_scaling.py')], cwd=REPO_ROOT)

def run_gpu():
    print('\n==== testsys: gpu ====')
    import importlib.util
    try:
        import jax
        devs = [d for d in jax.devices() if 'cuda' in str(d).lower() or 'gpu' in str(d).lower()]
    except Exception as e:
        print(f'SKIP gpu: jax not importable ({e})'); return 0
    if not devs:
        print('SKIP gpu: no CUDA device/plugin (pip install "jax[cuda12]" on a GPU box)'); return 0
    env = dict(os.environ, EQDYNA_ACCEPT_CASES='test.tpv8', EQDYNA_ACCEPT_PLATFORM='cuda')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'parity', 'test_standalone_acceptance.py')],
                           cwd=REPO_ROOT, env=env)

def run_perf():
    print('\n==== testsys: perf ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_perf.py')],
                            cwd=REPO_ROOT)


RUNNERS = {'unit': run_unit, 'regression': run_regression, 'e2e': run_e2e,
           'parity': run_parity, 'accept': run_accept, 'perf': run_perf,
           'gpu': run_gpu, 'scaling': run_scaling, 'e2e-full': run_e2e_full}
# 'all' stays unit+regression+e2e only (TIERS below) -- parity/accept/perf
# require a Fortran build and generated fixtures/baseline (accept also
# needs the committed test.reference.results/ trees) that a fresh checkout
# does not have; they are opt-in tiers, invoked by name, not swept into 'all'.
# e2e-full additionally needs EQDYNA_FULL_LAUNCH=yes-hours (see
# testsys/e2e/run_e2e_full.py) -- spec-resolution SCEC runs are hours long
# and user-scheduled, never automatic.
OPTIONAL_TIERS = ('parity', 'accept', 'perf', 'gpu', 'scaling', 'e2e-full')


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
