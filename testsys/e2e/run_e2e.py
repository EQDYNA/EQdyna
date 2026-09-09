#! /usr/bin/env python3
"""
End-to-end pipeline tier (PROJECT_RULES.md rules 3, 7, 8, 9).

Orchestrates the *existing* create.newcase -> case.setup -> mpirun ->
check.test.py flow for every case in testNameList.nameList. It does not
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
     mpirun -np <n> eqdyna. A non-zero exit from any step is recorded and
     the remaining steps for THAT case are skipped, but other cases still
     run, so one crash doesn't hide evidence about the rest.
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

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

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
            ([MPIRUN, '-np', str(coreNum), 'eqdyna'], case_dir),
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

    if run_failures or check_rc != 0:
        return 1
    print('e2e: SUCCESS - all cases ran and matched test.reference.results/')
    return 0


if __name__ == '__main__':
    sys.exit(main())
