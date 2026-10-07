#! /usr/bin/env python3
"""
Board item 148: python-jax-mpi SOLVE-LEVEL coverage for a mixed-fault case,
at the REGISTERED gated multi-fault layout.

WHY THIS EXISTS. test_rank_local_mesh.py's 'tpv23-2fault' entry proves the
mixed-fault master-node-id lookup fix (MPI4NodalQuant.py:159-174) builds the
right RANK-LOCAL MESH at (2, 4) ranks on a coarse-dx stand-in case -- it never
launches mpirun, never calls jax, never time-steps. testsys/matrix.py has
opted test.tpv23 into PY_MPI_RANKS at 4 ranks since board item 145, and
run_e2e.py would run that cell under `run_e2e.py --cases test.tpv23
--backends python-jax-mpi` (item 145's own PR measured it directly:
max|diff|=7.7e-15, bound 1e-10, 374.8s wall) -- but nothing commits that
invocation as a standing check: test.tpv23 is RELEASE_ONLY (matrix.py, cost
reason, item 17 section B), so this cell never runs in the everyday sweep,
in CI, or in any committed regression script. Board item 145's own numbers
came from a one-off campaign run, not a re-runnable test. This file closes
that gap: it actually RUNS the solver (mpirun -np 4 python3 -m eqdyna ...
--backend jax --mpi, via run_e2e.run_cell -- the EXACT function the e2e
sweep itself calls, no reimplemented launch path) on the real, full-resolution
test.tpv23 case (two vertical strike-slip faults, SCEC TPV22/23) at its
registered PY_MPI_RANKS=4 / MPI4NodalQuant.DECOMP[4]=(2,2,1) split and its own
CASE_TERM_OVERRIDE=15.0s SCEC-spec term, then compares against the ONE
committed Fortran reference with the ONE shared comparison tool
(testsys/compare.py's compare_cell) -- the same artifacts (frt, nsign,
station), the same CASE_BOUND, the same station-selection convention every
other python-jax-mpi cell in the table uses. No new comparison path, no new
tolerance.

COST (why this is NOT in the everyday regression sweep). This is a real
15-second-term, 4-rank, two-fault physics solve -- board item 145 measured
374.8s wall for this exact cell. testsys/run.py's run_regression() sweep is
documented (regression_sweep_exclusions.py) to stay seconds-scale and runs on
every commit and every PR's CI shard; a 6-minute script there would break
that contract for every contributor, not just this one. This file is
EXCLUDED_FROM_SWEEP (testsys/regression_sweep_exclusions.py) for that reason
alone -- it is a real, directly runnable regression script (same SUCCESS/FAIL
banner + sys.exit shape as every other regression/test_*.py), gated by hand
and as part of the release process (paired with `run.py release`, which runs
every RELEASE_ONLY cell including this one's own (test.tpv23,
python-jax-mpi) table entry at the identical bound/term -- this file is the
MECHANICAL, re-runnable proof that that specific cell's launch path and
comparison actually land green, independent of any one-off campaign run).

Run directly:
    python3 testsys/regression/test_row148_tpv23_jaxmpi_solve.py

No mesh-only shortcut, no coarse-dx stand-in, no reduced rank count: the
whole point of this file is that board item 148 found solve-level coverage
missing at the ACTUAL gated multi-fault layout, not at some faster but
different configuration. See test_rank_local_mesh.py's 'tpv23-2fault' entry
for the (fast, mesh-only) complement this file does NOT duplicate.
"""
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTSYS)
E2E = os.path.join(TESTSYS, 'e2e')
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
sys.path.insert(0, E2E)

os.chdir(ROOT)  # run_e2e's subprocess calls assume REPO_ROOT as cwd (matches
                # how testsys/run.py invokes every regression/test_*.py)

import run_e2e  # noqa: E402
from testsys import compare, matrix  # noqa: E402

CASE = 'test.tpv23'
BACKEND = 'python-jax-mpi'


def main():
    if CASE not in matrix.PY_MPI_RANKS:
        print('FAIL row148: %s has no matrix.PY_MPI_RANKS entry -- this '
              'file exists specifically to cover its REGISTERED opt-in; a '
              'missing entry means that registration itself regressed.'
              % CASE)
        return 1
    ranks = matrix.PY_MPI_RANKS[CASE]
    term = matrix.gate_term_for(CASE)
    print('row148: running %s x %s at %d ranks (MPI4NodalQuant.DECOMP[%d]), '
          'term=%.1fs (matrix.CASE_TERM_OVERRIDE) -- this is a real physics '
          'solve and can take several minutes.' % (CASE, BACKEND, ranks, ranks, term))

    test_dir = tempfile.mkdtemp(prefix='row148_tpv23_jaxmpi_')
    try:
        env = run_e2e.base_env()
        case_dir = run_e2e.run_cell(CASE, BACKEND, test_dir, None, env, 'cpu')
        ok, lines = compare.compare_cell(CASE, BACKEND, case_dir)
    finally:
        shutil.rmtree(test_dir, ignore_errors=True)

    for line in lines:
        print(line)
    print('%s row148: %s x %s solve-level parity vs the committed Fortran '
          'reference (%d ranks, term=%.1fs)'
          % ('SUCCESS' if ok else 'FAIL', CASE, BACKEND, ranks, term))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
