#! /usr/bin/env python3
"""
Regression guard for check.test.py itself (PROJECT_RULES.md rules 2, 3).

Guards a past defect: check.test.py used to swallow the one exception path
(xr.testing.assert_allclose's AssertionError) that could have turned a
mismatch into a failing exit code, and its .nc comparison was commented out
entirely -- so a real reference mismatch printed "FAIL" to the log but
`python3 check.test.py; echo $?` still reported 0 (pathway_forward.md items
3-5).

This builds a tiny run tree holding ONE case (test.tpv8) whose artifacts are
copies of the committed reference for that case -- so a clean copy must
compare exactly, and the perturbation is the only difference -- and runs the
*actual* check.test.py from the repo root against it via its --test-root
selection, asserting:
  - identical reference copies -> exit 0   (no false positive)
  - one perturbed value        -> exit != 0 (the actual regression guard)
for both compared artifacts (canonical frt, and fault.dyna.r.nc).

The sandbox no longer fakes a testNameList: check.test.py now takes
--test-root/--cases, and an unknown case name is itself a hard failure
(testsys/matrix.py has no entry for it), so the guard uses a real case name
against a fake tree rather than the reverse.

Cheap (rule 9): no build, no MPI, a few small file copies, well under 1 s.
Exits non-zero on any failure (rule 2).
"""
import os, shutil, subprocess, sys, tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from testsys import frt_canonical  # noqa: E402

CASEID = 'test.tpv8'
REF_DIR = os.path.join(ROOT, 'test.reference.results', CASEID)
FIXTURE_TXT = os.path.join(REF_DIR, 'frt.canonical.txt')
FIXTURE_NC = os.path.join(REF_DIR, 'fault.dyna.r.nc')


def _make_sandbox(tmp):
    """A run tree at <tmp>/test/test.tpv8 for check.test.py --test-root."""
    run_dir = os.path.join(tmp, 'test', CASEID)
    os.makedirs(run_dir)
    return run_dir


def _run_check_test(tmp):
    r = subprocess.run([sys.executable, 'check.test.py',
                        '--test-root', os.path.join(tmp, 'test'),
                        '--cases', CASEID],
                       cwd=ROOT, capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr


def _perturb_txt_line(src, dst):
    arr = np.loadtxt(src)
    arr[0, 4] += 10.0        # a physics column, far beyond any case bound
    frt_canonical.write_canonical(arr, dst)


def _perturb_nc(src, dst):
    import xarray as xr
    ds = xr.open_dataset(src)
    var = list(ds.data_vars)[0] if ds.data_vars else list(ds.variables)[0]
    ds = ds.copy(deep=True)
    ds[var] = ds[var] + 1000.0  # far beyond rtol/atol=1e-3
    ds.to_netcdf(dst)
    ds.close()


fails = []

if not os.path.exists(FIXTURE_TXT):
    fails.append(f'fixture missing, cannot run guard: {FIXTURE_TXT}')
if not os.path.exists(FIXTURE_NC):
    fails.append(f'fixture missing, cannot run guard: {FIXTURE_NC}')

if not fails:
    # --- a clean copy of both artifacts must pass ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        run_dir = _make_sandbox(tmp)
        shutil.copy(FIXTURE_TXT, os.path.join(run_dir, 'frt.txt0'))
        shutil.copy(FIXTURE_NC, os.path.join(run_dir, 'fault.dyna.r.nc'))
        rc, out = _run_check_test(tmp)
        if rc != 0:
            fails.append(f'false positive: an exact copy of the reference reported '
                         f'non-zero exit ({rc})\n{out}')

    # --- a perturbed frt value must fail loudly ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        run_dir = _make_sandbox(tmp)
        _perturb_txt_line(FIXTURE_TXT, os.path.join(run_dir, 'frt.txt0'))
        shutil.copy(FIXTURE_NC, os.path.join(run_dir, 'fault.dyna.r.nc'))
        rc, out = _run_check_test(tmp)
        if rc == 0:
            fails.append(f'REGRESSION: check.test.py exited 0 on a perturbed frt value\n{out}')
        if 'FAIL' not in out:
            fails.append(f'check.test.py did not print FAIL for a perturbed frt value\n{out}')

    # --- a perturbed fault.dyna.r.nc must fail loudly ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        run_dir = _make_sandbox(tmp)
        shutil.copy(FIXTURE_TXT, os.path.join(run_dir, 'frt.txt0'))
        _perturb_nc(FIXTURE_NC, os.path.join(run_dir, 'fault.dyna.r.nc'))
        rc, out = _run_check_test(tmp)
        if rc == 0:
            fails.append(f'REGRESSION: check.test.py exited 0 on a perturbed '
                         f'fault.dyna.r.nc\n{out}')
        if 'FAIL' not in out:
            fails.append(f'check.test.py did not print FAIL for a perturbed '
                         f'fault.dyna.r.nc\n{out}')

    # --- a missing artifact must fail, not narrow the comparison silently ---
    # The fortran backend is declared (matrix.ARTIFACTS) to produce frt AND
    # nc. A run that wrote only frt used to be compared on frt alone, with
    # nothing in the output saying so.
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        run_dir = _make_sandbox(tmp)
        shutil.copy(FIXTURE_TXT, os.path.join(run_dir, 'frt.txt0'))
        rc, out = _run_check_test(tmp)
        if rc == 0:
            fails.append(f'REGRESSION: check.test.py exited 0 with '
                         f'fault.dyna.r.nc missing entirely\n{out}')

if fails:
    print('FAIL test_check_gate_fails_on_mismatch')
    for m in fails:
        print('  -', m)
    sys.exit(1)
print('SUCCESS test_check_gate_fails_on_mismatch')
