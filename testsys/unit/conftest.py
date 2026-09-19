"""Puts repo root, scripts/, and src/python on sys.path for the unit tier.

scripts/lib.py, scripts/defaultParameters.py, and root check.test.py all
assume they're importable from a directory that also has their sibling
modules (lib.py, testNameList.py) importable next to them -- exactly how
they're used in a real case directory or at the repo root. We replicate
that here rather than changing how the modules import.

src/python is here too so any test can `from eqdyna import ...` without
repeating this same REPO_ROOT/sys.path setup itself -- test_backend_jax.py,
test_assembleGlobalKU_scatter.py, test_eqdyna3d_numpy_affinity.py, and
test_meshgen_c_degen_wedge.py used to each carry an identical copy of this
before relying on this file for it.
"""
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')
SRC_PYTHON = os.path.join(REPO_ROOT, 'src', 'python')

for _p in (REPO_ROOT, SCRIPTS, SRC_PYTHON):
    if _p not in sys.path:
        sys.path.insert(0, _p)
