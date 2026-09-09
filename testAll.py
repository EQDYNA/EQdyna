#! /usr/bin/env python3
"""
Thin compatibility wrapper (PROJECT_RULES.md rule 3's named gate command).

`python3 testAll.py` is documented in README.md and habitually typed by
hand; this just delegates to the e2e tier of testsys/ so that command
keeps working unchanged. The real implementation -- fresh build, rule-8
evidence preservation, per-case runs, check.test.py gate -- lives in
testsys/e2e/run_e2e.py. See also `python3 testsys/run.py [unit|regression|e2e|all]`.
"""
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    sys.exit(subprocess.call(
        [sys.executable, os.path.join(REPO_ROOT, 'testsys', 'e2e', 'run_e2e.py')],
        cwd=REPO_ROOT))
