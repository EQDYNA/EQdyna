# testsys/parity/

**The parity tier no longer exists.** It was deleted on 2026-09-15 along with
the pydump Fortran instrumentation it was built on. What is left are
REPORT-ONLY evidence scripts: never a gate, never wired into `testsys/run.py`,
invoked by hand.

## What is left

`evidence_drv_a6_chaos.py` — regenerates the measurements behind
`test.drv.a6`'s flip-budget bound (`testsys/matrix.py`'s `DRV_A6`), so that
bound is derived from something reproducible rather than remembered from a
commit message. It imports `testsys.compare`, `testsys.matrix` and `run_e2e`
rather than duplicating their alignment and comparison logic (rule 1).

    python3 testsys/parity/evidence_drv_a6_chaos.py

`evidence_tpv29_scec_comparison.py` — regenerates the four TPV29 cross-code
validation numbers behind pathway_forward.md items 20/28 (median rupture-time
diff, median final-slip % diff, ruptured area, Mw, at matched dx=100 m
against EQdyna's own 2015 SCEC submission in `scec_archive/tpv29/`) instead of
leaving them prose-only. Loads the archived 2015 submission plus an
already-completed local dx=100 m run (`scratch/tpv29/ny48_dx100/` by default
-- neither directory is tracked, both gitignored); does not launch a run
itself (rule 9 -- a fresh 100 m/48-rank TPV29 run is not cheap).

    python3 testsys/parity/evidence_tpv29_scec_comparison.py
    python3 testsys/parity/evidence_tpv29_scec_comparison.py --run-dir /path/to/run --no-json

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

## A kernel-level parity check DOES exist now -- it lives in regression/, not here

pathway_forward.md item 23 (Drucker-Prager viscoplastic return-mapping,
`calcElemKU.f90:127-161`, had zero isolated test) is closed by
`testsys/regression/test_drucker_prager_kernel.py`. It compiles a tiny
driver against the real, unmodified `calcElemKU.f90` (+ its 3 link
dependencies) with plain gfortran -- no MPI, no case, no mesh -- and diffs
its output against `assembleGlobalKU._drucker_prager` (the Python port)
directly, on 4 synthetic stress states (elastic, yielding, near-yield-
boundary, shear-dominated). It was NOT placed in this directory: `parity/`
is the tier that no longer exists (see above), while `regression/`'s
"standalone script, SUCCESS/FAIL banner, `sys.exit`, picked up by
`testsys/run.py regression`'s glob" shape is exactly what this needed and
wires it into the gate with zero changes to `run.py`. Reproduce directly:

    python3 testsys/regression/test_drucker_prager_kernel.py

## Open: this directory's name

`parity/` is named for a tier that is gone, and holds one file that is not a
parity test. Renaming it (`testsys/evidence/`) would touch `testsys/matrix.py`,
`testsys/compare.py`, `testsys/perf/run_perf.py`, `pathway_forward.md` and
`case_input/test.drv.a6/README.md`. Left for an owner decision rather than
folded into a release.
