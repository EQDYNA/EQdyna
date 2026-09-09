#! /usr/bin/env python3
"""
Regression test for create.newcase (PROJECT_RULES.md rule 10).

Guards two past failures:
  1. IsADirectoryError when scripts/ contains a subdirectory
     (scripts/fractal_stress_diamond_square broke every case creation,
     2026-09-09).
  2. Silent PATH-shadowing / missing exec bit: the test invokes the
     repo's own scripts/create.newcase explicitly, so a wrong or
     non-executable script fails here, loudly. install-eqdyna.sh is
     checked alongside it (same guard family): it is invoked as
     `./install-eqdyna.sh` by testsys/e2e/run_e2e.py's fresh-rebuild
     step, and a lost exec bit there fails a full e2e run the same
     way a lost exec bit on create.newcase fails case creation.

Cheap targeted check (rule 9): no build, no MPI, ~1 s.
Exits non-zero on any failure (rule 2).
"""
import os, shutil, subprocess, sys, tempfile

root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ['EQDYNAROOT'] = root
script = os.path.join(root, 'scripts', 'create.newcase')
install_script = os.path.join(root, 'install-eqdyna.sh')

fails = []
if not os.access(script, os.X_OK):
    fails.append(f'{script} is not executable')
if not os.access(install_script, os.X_OK):
    fails.append(f'{install_script} is not executable')

tmp = tempfile.mkdtemp(prefix='testCreateNewcase.')
case = os.path.join(tmp, 'case')
try:
    r = subprocess.run([sys.executable, script, case, 'test.drv.a6'],
                       capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        fails.append(f'create.newcase exited {r.returncode}:\n{r.stderr}')
    else:
        for f in ['user_defined_params.py', 'case.setup', 'create.newcase',
                  'generateFaultInterface']:
            if not os.path.isfile(os.path.join(case, f)):
                fails.append(f'expected file missing in new case: {f}')
        # subdirectories of scripts/ must NOT be copied into a case
        copied_dirs = [d for d in os.listdir(case)
                       if os.path.isdir(os.path.join(case, d))]
        if copied_dirs:
            fails.append(f'directories copied into case: {copied_dirs}')
finally:
    shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print('FAIL testCreateNewcase')
    for m in fails:
        print('  -', m)
    sys.exit(1)
print('SUCCESS testCreateNewcase')
