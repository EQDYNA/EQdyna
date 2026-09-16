# testsys/parity/

**The parity tier no longer exists.** It was deleted on 2026-09-15 along with
the pydump Fortran instrumentation it was built on. One file remains here.

## What is left

`evidence_drv_a6_chaos.py` — a REPORT-ONLY script, never a gate, never wired
into `testsys/run.py`. It regenerates the measurements behind `test.drv.a6`'s
flip-budget bound (`testsys/matrix.py`'s `DRV_A6`), so that bound is derived
from something reproducible rather than remembered from a commit message.
It imports `testsys.compare`, `testsys.matrix` and `run_e2e` rather than
duplicating their alignment and comparison logic (rule 1).

    python3 testsys/parity/evidence_drv_a6_chaos.py

## What was removed, and why

The tier answered "at WHICH STEP does the port diverge?" by dumping Fortran's
internal per-step state and comparing it to the Python port's. That is a
debugging question, not a gating one. It did its job — the per-step dumps
root-caused a missing `rsfNucleation` branch — and then went stale: its two
test files were named for a `standalone/` package the Python restructure
removed, and its fixtures were frozen copies of `scripts/` that had drifted
badly (`lib.py` stuck at 187 lines against the live 949) and were kept
lint-clean only by an explicit carve-out in `test_ci_dependencies.py`.

Deleted with it: `run_parity.py`, `test_standalone_meshgen.py`,
`test_standalone_setup_vectorized.py`, `make_fixtures.py`, `fixtures/`, and on
the Fortran side `src/pydump.f90`, `src/pydump_noop.f90`, the `PYDUMP=1`
makefile seam and both call sites. Recover any of it from git history if a
per-step diff is ever needed again; regenerating a step dump is cheaper than
keeping a stale tier green.

There is ONE test in this repo now: the e2e sweep (`testsys/run.py e2e`),
8 cases x 3 backends, 24 cells, compared by `testsys/compare.py` against one
`frt.canonical.txt` per case.

## Open: this directory's name

`parity/` is named for a tier that is gone, and holds one file that is not a
parity test. Renaming it (`testsys/evidence/`) would touch `testsys/matrix.py`,
`testsys/compare.py`, `testsys/perf/run_perf.py`, `pathway_forward.md` and
`case_input/test.drv.a6/README.md`. Left for an owner decision rather than
folded into a release.
