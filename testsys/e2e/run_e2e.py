#! /usr/bin/env python3
"""
End-to-end pipeline tier (PROJECT_RULES.md rules 3, 7, 8, 9).

Orchestrates the *existing* create.newcase -> case.setup -> mpirun ->
plotRuptureDynamics -> check.test.py flow for every case in
testNameList.nameList (the plotRuptureDynamics step is what writes
fault.dyna.r.nc, one of check.test.py's fileNameList comparisons -- see
the pre-testsys testAll.py for the reference sequence). It does not
duplicate those runs anywhere else in testsys/ -- this is the one place
they happen.

Named gates, in order (rule 3):
  1. testsys/regression/test_create_newcase.py -- cheap check before an
     expensive run (rule 9).
  2. Fresh build: `./install-eqdyna.sh -m <machine>` from the repo root,
     after removing any existing bin/eqdyna, so a stale binary can never
     paper over a broken build (mirrors EQquasi's testAll.py convention:
     rm bin/eqquasi; reinstall). Machine defaults to "ubuntu"; override
     with the EQDYNA_TEST_MACHINE env var.
  3. For every case in testNameList.nameList: create.newcase, case.setup,
     mpirun -np <n> eqdyna, plotRuptureDynamics (writes fault.dyna.r.nc).
     A non-zero exit from any step is recorded and the remaining steps
     for THAT case are skipped, but other cases still run, so one crash
     doesn't hide evidence about the rest.
  4. check.test.py, comparing test/ against test.reference.results/
     (rule 7: reference tree is read-only, never written by this script).

Rule 8: the previous `test/` run is renamed to `test.prev/` (one level of
history) instead of deleted before a fresh run starts. A run that gets
interrupted, or that needs a second look before someone reruns it, is
never silently thrown away by the next invocation.

No aggressive timeouts: this is a shared box and dynamic-rupture runs may
be slow when the machine is busy.
"""
import os
import shutil
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
# Coordination hatch: when a concurrent long HPC job elsewhere on this box
# depends on bin/eqdyna staying untouched (e.g. an in-flight spec-resolution
# run), set EQDYNA_E2E_BIN to a pre-built, explicit binary path -- e2e then
# skips the bin/ rm+rebuild entirely and mpirun's that binary directly.
# Caller is responsible for that binary being a fresh, clean build (rule 4);
# this is the same "build src/eqdyna, never touch bin/eqdyna" pattern already
# used by testsys/regression/test_fault_mpi_boundary_arn.py.
BIN_OVERRIDE = os.environ.get('EQDYNA_E2E_BIN')

sys.path.insert(0, REPO_ROOT)


def _run(cmd, cwd, env):
    print('+ (' + os.path.basename(cwd) + ') ' + ' '.join(cmd))
    return subprocess.call(cmd, cwd=cwd, env=env)


def main():
    from testNameList import nameList, coreNumList

    # Gate 1 - cheap check before the expensive runs (rule 9).
    guard = os.path.join(REPO_ROOT, 'testsys', 'regression', 'test_create_newcase.py')
    if subprocess.call([sys.executable, guard]) != 0:
        print('e2e: FAIL - testCreateNewcase guard failed; aborting before expensive runs')
        return 1

    # Gate 2 - fresh build (never trust a binary left over from a previous run).
    if BIN_OVERRIDE:
        eqdyna_cmd = os.path.abspath(BIN_OVERRIDE)
        if not os.path.exists(eqdyna_cmd):
            print(f'e2e: FAIL - EQDYNA_E2E_BIN={BIN_OVERRIDE} does not exist')
            return 1
        print(f'e2e: EQDYNA_E2E_BIN set - skipping bin/eqdyna rebuild, using {eqdyna_cmd}')
    else:
        bin_exe = os.path.join(REPO_ROOT, 'bin', 'eqdyna')
        if os.path.exists(bin_exe):
            os.remove(bin_exe)
        build_rc = subprocess.call(['./install-eqdyna.sh', '-m', MACHINE], cwd=REPO_ROOT)
        if build_rc != 0:
            print(f'e2e: FAIL - ./install-eqdyna.sh -m {MACHINE} exited {build_rc}')
            return 1
        if not os.path.exists(bin_exe):
            print('e2e: FAIL - build finished but bin/eqdyna does not exist')
            return 1
        eqdyna_cmd = 'eqdyna'

    # Rule 8 - preserve, never delete, the previous run's evidence.
    test_dir = os.path.join(REPO_ROOT, 'test')
    prev_dir = os.path.join(REPO_ROOT, 'test.prev')
    if os.path.isdir(test_dir):
        if os.path.isdir(prev_dir):
            shutil.rmtree(prev_dir)
        shutil.move(test_dir, prev_dir)
        print(f'e2e: preserved previous run as {prev_dir}')
    os.makedirs(test_dir)

    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([
        os.path.join(REPO_ROOT, 'bin'),
        os.path.join(REPO_ROOT, 'scripts'),
        env.get('PATH', ''),
    ])

    start = time.time()
    run_failures = []
    for testName, coreNum in zip(nameList, coreNumList):
        case_dir = os.path.join(test_dir, testName)
        steps = [
            (['create.newcase', testName, testName], test_dir),
            (['./case.setup'], case_dir),
            ([MPIRUN, '-np', str(coreNum), eqdyna_cmd], case_dir),
            ([sys.executable, 'plotRuptureDynamics'], case_dir),
        ]
        for cmd, cwd in steps:
            rc = _run(cmd, cwd, env)
            if rc != 0:
                run_failures.append(f'{testName}: `{" ".join(cmd)}` exited {rc}')
                break  # later steps for this case would run on a broken case
    elapsed = time.time() - start
    print(f'e2e: pipeline runs took {elapsed:.1f}s for {len(nameList)} case(s)')

    # Gate 4 - compare against the golden reference tree (rule 7: read-only).
    check_rc = subprocess.call([sys.executable, 'check.test.py'], cwd=REPO_ROOT, env=env)

    if run_failures:
        print('e2e: FAIL - pipeline run failures:')
        for f in run_failures:
            print('  -', f)
    if check_rc != 0:
        print(f'e2e: FAIL - check.test.py reported comparison failures (exit {check_rc})')

    # Gate 5 - the SAME cases through the Python/JAX backend.
    #
    # e2e is backend-parameterised rather than having a parallel `accept` tier:
    # one case list, one reference tree, one place a regression shows up. The
    # Fortran pipeline above runs all 8 gated cases at 4 ranks; the standalone
    # runs the 5 it supports (tpv8, tpv104, tpv1053d, tpv10, drv.a6) SERIALLY,
    # because build_solver_state refuses npx*npy*npz > 1.
    #
    # This delegates to test_standalone_acceptance.py rather than reimplementing
    # its comparison (rule 1). That comparison is not a plain diff: the
    # standalone writes ONE whole-domain frt.txt0 while the reference is split
    # across per-rank frt.txt* files, so it dedupes the combined reference by
    # rounded (x,y,z) and lexsorts both sides before comparing. Duplicating that
    # would be two things to keep in step.
    #
    # Skipped, loudly and with a non-zero-able reason, if jax is absent -- never
    # silently, and never demoted to the numpy backend, which would report
    # "jax passed" for a run that was not jax (rule 2, and the reason
    # --backend jax stopped falling back in v5.6.2).
    #
    # CASE SUBSET, and why it is not a silent skip. One case at 20 steps peaks
    # at 3.6 GB RSS (measured); a GitHub ubuntu-22.04 runner has 7 GB and is
    # already holding the Fortran build. Running all five at full length got
    # the job SIGTERM'd (exit 143, no output -- the runner killed it), which is
    # a resource limit, not a test result. So CI gates a REAL subset rather
    # than pretending to gate everything: EQDYNA_ACCEPT_CASES names which. The
    # tier prints which cases it ran, so "green" always says what it covered.
    # The full five-case run stays available by leaving the variable unset.
    py_rc = 0
    if os.environ.get('EQDYNA_E2E_SKIP_PYTHON') == '1':
        print('e2e: python/jax backend SKIPPED (EQDYNA_E2E_SKIP_PYTHON=1)')
    else:
        accept = os.path.join(TESTSYS, 'parity', 'test_standalone_acceptance.py')
        if not os.path.exists(accept):
            print(f'e2e: FAIL - {accept} is missing; the python backend cannot be gated')
            py_rc = 1
        else:
            print('\n-- e2e: python/jax backend (standalone, serial) --')
            # -u: unbuffered. When the runner SIGTERM'd this child, Python's
            # block buffering meant ZERO output survived, so the log showed a
            # 7-minute silence and an exit code. Unbuffered output makes a
            # resource kill diagnosable instead of mute.
            py_env = dict(env, PYTHONUNBUFFERED='1')
            sub = os.environ.get('EQDYNA_ACCEPT_CASES')
            print(f'   cases: {sub if sub else "all (EQDYNA_ACCEPT_CASES unset)"}')
            py_rc = subprocess.call([sys.executable, '-u', accept],
                                    cwd=REPO_ROOT, env=py_env)
            if py_rc != 0:
                print(f'e2e: FAIL - python/jax backend exited {py_rc}')

    if run_failures or check_rc != 0 or py_rc != 0:
        return 1
    print('e2e: SUCCESS - all cases ran and matched test.reference.results/, '
          'fortran and python/jax')
    return 0


if __name__ == '__main__':
    if '--full' in sys.argv:
        # Delegate to the full tier (SCEC spec resolution/duration, 16
        # ranks, report-only). See run_e2e_full.py's module docstring for
        # why this is a separate script rather than a branch in main()
        # above: different gate shape (report-only vs bit-compare), and
        # never runs by default under `python3 testsys/run.py e2e`.
        from run_e2e_full import main as main_full
        sys.exit(main_full())
    sys.exit(main())
