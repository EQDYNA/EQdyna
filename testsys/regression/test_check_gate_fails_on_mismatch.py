#! /usr/bin/env python3
"""
Regression guard for check.test.py itself (PROJECT_RULES.md rules 2, 3).

Guards a past defect: check.test.py used to swallow the one exception path
(xr.testing.assert_allclose's AssertionError) that could have turned a
mismatch into a failing exit code, and its .nc comparison was commented out
entirely -- so a real reference mismatch printed "FAIL" to the log but
`python3 check.test.py; echo $?` still reported 0 (pathway_forward.md items
3-5).

This builds two tiny sandboxes -- one frt.txt-only, one fault.dyna.r.nc-only
-- each with a real fixture copied from test.reference.results/test.tpv8/
(so the numbers are representative, not invented), runs the *actual*
check.test.py from the repo root against each sandbox, and asserts:
  - identical ref/test copies  -> exit 0   (no false positive)
  - one perturbed value        -> exit != 0 (the actual regression guard)

Cheap (rule 9): no build, no MPI, a few small file copies, well under 1 s.
Exits non-zero on any failure (rule 2).
"""
import os, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CHECK_TEST_SRC = os.path.join(ROOT, 'check.test.py')
FIXTURE_TXT = os.path.join(ROOT, 'test.reference.results', 'test.tpv8', 'frt.txt0')
FIXTURE_NC = os.path.join(ROOT, 'test.reference.results', 'test.tpv8', 'fault.dyna.r.nc')

CASEID = 'sandboxcase'


def _make_sandbox(tmp, filename):
    """ref/test trees + testNameList.py + a copy of check.test.py, so that
    running `python3 check.test.py` with cwd=tmp picks up the sandbox's own
    (small) testNameList instead of the repo's real 5-case list."""
    ref_dir = os.path.join(tmp, 'test.reference.results', CASEID)
    test_dir = os.path.join(tmp, 'test', CASEID)
    os.makedirs(ref_dir)
    os.makedirs(test_dir)
    with open(os.path.join(tmp, 'testNameList.py'), 'w') as f:
        f.write(f"nameList = ['{CASEID}']\ncoreNumList = [1]\n")
    shutil.copy(CHECK_TEST_SRC, os.path.join(tmp, 'check.test.py'))
    return ref_dir, test_dir


def _run_check_test(tmp):
    r = subprocess.run([sys.executable, 'check.test.py'], cwd=tmp,
                        capture_output=True, text=True, timeout=30)
    return r.returncode, r.stdout + r.stderr


def _perturb_txt_line(src, dst):
    with open(src) as f:
        nums = f.read().split()
    nums[0] = str(float(nums[0]) + 10.0)  # far beyond THRESHOLD=1e-3
    with open(dst, 'w') as f:
        f.write(' '.join(nums))


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
    # --- frt.txt: identical copies must pass ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        ref_dir, test_dir = _make_sandbox(tmp, 'frt.txt0')
        shutil.copy(FIXTURE_TXT, os.path.join(ref_dir, 'frt.txt0'))
        shutil.copy(FIXTURE_TXT, os.path.join(test_dir, 'frt.txt0'))
        rc, out = _run_check_test(tmp)
        if rc != 0:
            fails.append(f'false positive: identical frt.txt reported non-zero exit ({rc})\n{out}')

    # --- frt.txt: a perturbed value must fail loudly ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        ref_dir, test_dir = _make_sandbox(tmp, 'frt.txt0')
        shutil.copy(FIXTURE_TXT, os.path.join(ref_dir, 'frt.txt0'))
        _perturb_txt_line(FIXTURE_TXT, os.path.join(test_dir, 'frt.txt0'))
        rc, out = _run_check_test(tmp)
        if rc == 0:
            fails.append(f'REGRESSION: check.test.py exited 0 on a perturbed frt.txt0\n{out}')
        if 'FAIL' not in out:
            fails.append(f'check.test.py did not print FAIL for a perturbed frt.txt0\n{out}')

    # --- fault.dyna.r.nc: identical copies must pass ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        ref_dir, test_dir = _make_sandbox(tmp, 'fault.dyna.r.nc')
        shutil.copy(FIXTURE_NC, os.path.join(ref_dir, 'fault.dyna.r.nc'))
        shutil.copy(FIXTURE_NC, os.path.join(test_dir, 'fault.dyna.r.nc'))
        rc, out = _run_check_test(tmp)
        if rc != 0:
            fails.append(f'false positive: identical fault.dyna.r.nc reported non-zero exit ({rc})\n{out}')

    # --- fault.dyna.r.nc: a perturbed value must fail loudly ---
    with tempfile.TemporaryDirectory(prefix='check_test_gate.') as tmp:
        ref_dir, test_dir = _make_sandbox(tmp, 'fault.dyna.r.nc')
        shutil.copy(FIXTURE_NC, os.path.join(ref_dir, 'fault.dyna.r.nc'))
        _perturb_nc(FIXTURE_NC, os.path.join(test_dir, 'fault.dyna.r.nc'))
        rc, out = _run_check_test(tmp)
        if rc == 0:
            fails.append(f'REGRESSION: check.test.py exited 0 on a perturbed fault.dyna.r.nc\n{out}')
        if 'FAIL' not in out:
            fails.append(f'check.test.py did not print FAIL for a perturbed fault.dyna.r.nc\n{out}')

if fails:
    print('FAIL test_check_gate_fails_on_mismatch')
    for m in fails:
        print('  -', m)
    sys.exit(1)
print('SUCCESS test_check_gate_fails_on_mismatch')
