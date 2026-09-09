"""Puts repo root and scripts/ on sys.path for the unit tier.

scripts/lib.py, scripts/defaultParameters.py, and root check.test.py all
assume they're importable from a directory that also has their sibling
modules (lib.py, testNameList.py) importable next to them -- exactly how
they're used in a real case directory or at the repo root. We replicate
that here rather than changing how the modules import.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')

for _p in (REPO_ROOT, SCRIPTS):
    if _p not in sys.path:
        sys.path.insert(0, _p)
