# Pathway Forward

Record of tasks done and open issues for this repo. See
`PROJECT_RULES.md` for the rules each item enforces. History (what shipped,
when) lives in `pastReleaseNotes.md` and `README.md` — this file does not
duplicate it and is not append-only; items are updated or closed in place.

A blank "Last checked" means never audited and stays blank until someone
actually runs the command.

| # | Item | Rule | Surface | Re-check interval | Last checked | Command |
|---|---|---|---|---|---|---|
| 1 | ~~Untracked build artifacts~~ resolved by `.gitignore` update (8c40871); recheck each release | 12 | repo root, `src/` | every release | 2026-09-09 | `git status --porcelain` |
| 3 | ~~`check.test.py` swallows `AssertionError` ... never sets a non-zero exit~~ resolved (34691db): loop moved into `main()`, every comparison result accumulated, `sys.exit(1)` on any `FAIL`, missing file, or zero comparisons run; guarded by `testsys/regression/test_check_gate_fails_on_mismatch.py` | 2, 3 | `check.test.py` | each change to test infra | 2026-09-09 | `python3 check.test.py; echo exit=$?` |
| 4 | ~~`compare_nc_files` call commented out~~ resolved (34691db): call uncommented, `compare_nc_files` is a pure function with no import-time side effects, exercised directly in `testsys/unit/test_check_comparisons.py` and via the sandboxed `fault.dyna.r.nc` case in the regression guard above | 3 | `check.test.py` | each change to test infra | 2026-09-09 | `grep -n compare_nc_files check.test.py` |
| 5 | ~~`testAll.py` ran `rm -rf test` unconditionally~~ resolved (f847b3b): `testAll.py` is now a thin wrapper delegating to `testsys/e2e/run_e2e.py`, which renames the previous `test/` to `test.prev/` (one level of history) instead of deleting it | 8 | `testsys/e2e/run_e2e.py` | each change to test infra | 2026-09-09 | `grep -n "shutil.move(test_dir, prev_dir)" testsys/e2e/run_e2e.py` |
| 6 | ~~README Installation section instructed a recursive `chmod -R 755 install-eqdyna.sh scripts`~~ resolved (this release): `install-eqdyna.sh` itself no longer bulk-chmods `scripts/` (was doing so on every build, rule 13), so README now only needs `chmod 755 install-eqdyna.sh` | 13 | `README.md`, `install-eqdyna.sh` | each release | 2026-09-09 | `grep -n "chmod -R" README.md install-eqdyna.sh` |

| 7 | TO-DO: `netcdf_read_on_fault_eqdyna` hardcodes fault index `1` instead of loop var `ift` in all `fric(...)` assignments, so on-fault netcdf input is never applied to faults 2+ (fault-N data overwrites fault 1); `ii,jj` grid indices are also unchecked against bounds for faults 2+ (built from fault-1 extent). The restart routine below it (`netcdf_io.f90:159-179`) shows the correct `ift` pattern. Bites only `ntotft >= 2` — no current test case exercises it; fix needs a multi-fault regression case per rule 10. Does NOT explain the tpv1053d reference mismatch (that case runs `ntotft=1`; its golden reference is stale, carrying an uninitialized denormal). | 2, 10 | `src/netcdf_io.f90:76-105` | until fixed | 2026-09-09 | `grep -n "i, 1) = on_fault_vars" src/netcdf_io.f90` |

| 8 | ~~PML region-14 cross-axis bug~~ FIXED a844e92 (2026-09-13): y now vs ymax in all 4 sites; z==zmin2 gap closed. Impact was real: drv.a6 tstk moved up to 29 MPa; refs regenerated; guard test_pml_region_axes.py added | 2, 10 | `src/comdampv.f90`, `src/assembleGlobalKU.f90` | guarded | 2026-09-13 | `python3 testsys/regression/test_pml_region_axes.py` |
| 9 | BUG (latent, ntotft>=2): `output_onfault_st` opens unit 51 only when j==1 but writes unconditionally | 2 | `src/library_output.f90:29,37` | until fixed | 2026-09-12 | `grep -n "if(j==1)" src/library_output.f90` |
| 10 | BUG (latent, ntotft>=2): `MPI4arn` allocates fltl..fltu from previous fault's counts before zeroing fltnum; re-allocation aborts on fault 2 | 2 | `src/meshgen.f90:220-227` | until fixed | 2026-09-12 | `grep -n "allocate(flt" src/meshgen.f90` |
| 11 | CHECK: `mesh4num` on-fault tolerance vs `checkIsOnFault` may diverge (counting pass vs allocating pass must agree) | 3 | `src/mesh4num.f90:45-57`, `src/meshgen.f90:768-787` | until verified | 2026-09-12 | `grep -n "distToFault" src/mesh4num.f90 src/meshgen.f90` |

| 12 | ~~fric_tp_h=0 TP bug~~ FIXED (2026-09-13): thermop kernels now use per-node fric(40) (=0.02 m, TPV105-3D spec h=20mm). Verified against clean-room v5.2.0 (the SCEC-verified tag): fixed==v5.2.0, h=0 categorically different (~1 s rupture time, ~2.9 m slip). tpv1053d reference regenerated. Was: `fric_tp_h` (TP shear-zone half-width, denominator in every thermop.f90 kernel) is never assigned in src/ — runs as 0.0 (confirmed by compiled probe). Intended 0.02 is written to on_fault_vars slot 40 (fric(40)) by case.setup but never wired to the global thermop reads. All TP results from this tree used h=0. Fix is gate-moving (tpv1053d reference regen) — needs owner decision | 2 | `src/globalvar.f90` (fric_tp_h), `src/thermop.f90`, `src/netcdf_io.f90` (slot 40) | until fixed | 2026-09-13 | `grep -rn "fric_tp_h" src/` |

| 13 | ~~PML boundary inclusivity~~ RESOLVED (2026-09-13): both callers standardized to inclusive >=/<= (an on-bound point classifies with zero damping distance, harmless), and a mesh-time guard checkPMLAlignment stops loudly if any element center lies exactly on a PML bound (the cascade's remaining geometric assumption, now enforced not assumed). Bit-gate confirmed no output change. fric slot 6 retired e90ca03 (dead 'surface pore pressure' reads removed) | 2 | `src/func_lib.f90` | guarded | 2026-09-13 | `grep -n checkPMLAlignment src/meshgen.f90` |

| 14 | ~~parity reproducibility~~ RESOLVED (2026-09-13): testsys parity + perf tiers landed — pydump link-time seam (default binary byte-neutral), make_fixtures.py regenerates oracles, `run.py parity` gates NumPy+JAX vs Fortran at documented thresholds; `run.py perf` is pinned single-core and ratio-guarded vs baseline.json. Pinned verdict: JAX-CPU 0.61x Fortran wall-clock (1.6x faster), NumPy 1.8x slower. Was: `python/README-parity.md`'s parity numbers (NumPy/JAX ports vs Fortran, JAX-CPU 1.1-2.6x) were established against instrumented spike builds — `src/pydump.f90` (the oracle-dump instrumentation) is not in tracked `src/` or wired into `src/makefile`, `python/tests/data/` is missing the `pydump_*.txt` fixtures `run_parity.py` needs, and no tier under `testsys/` exercises `python/` at all. Running `run_parity.py` against the current tree fails with `FileNotFoundError: pydump_header.txt`. Not a correctness finding against the physics fixes in v5.3.6 — the port work and its parity claim are real, just not reproducible in-tree yet. Needs an owner decision: land `pydump.f90` + fixtures + a `testsys/` tier, or keep documenting the parity claim as spike-only in `python/README-parity.md`. | 2, 4, 10 | `python/`, `src/pydump.f90` (absent), `testsys/` | until decided | 2026-09-13 | `python3 python/run_parity.py tests/data/tpv8_serial 114 tests/data/tpv8_serial_frt.txt0.golden` |

**Deferred**: none.

## Goals

- **Python EQdyna (vectorized)** — long-term goal set 2026-09-12: a vectorized
  Python implementation of EQdyna, developed under `python/` at repo root
  (ALL Python solver development lives there). Phase 1: feasibility spike —
  port the hot kernel + time loop for `test.tpv8`, parity-gated against the
  Fortran `frt.txt` golden gate, with wall-clock comparison (NumPy first,
  JAX/GPU assessment after). Fortran `src/` remains the production/HPC engine;
  the Python solver targets test-scale cross-verification first, GPU later.

## Tasks done

| Date | Task | Ref |
|---|---|---|
| 2026-09-09 | GM/src-evol output every step; mpif90 on ubuntu | b1bc5e9 |
| 2026-09-09 | Project rule book seeded (14 rules); build artifacts + `scratch/` gitignored | 8c40871, 8414b5d |
| 2026-09-09 | Exec bits restored on `create.newcase`, `generateFaultInterface` (test suite was unrunnable from clean checkout — PATH-shadowed by EQquasi) | 27f8a76 |
| 2026-09-09 | v5.3.4 release audit run: build green, 4/5 reference cases pass; blocked on stale `test.tpv1053d` golden reference (uninitialized denormal captured pre-v5.3.3; HEAD == v5.3.3 byte-identical) | release not cut |
| 2026-09-09 | Root-caused tpv1053d frt.txt col 22 mismatch: old reference never wrote fric(23) (uninitialized denormal at all 2025 nodes) | 85d4d53 |
| 2026-09-09 | Fixed uninitialized `thetaPcTmp` in `NewtonRaphson` for friclaw==5 (src/faulting.f90); tpv1053d reference regenerated from the fixed binary, doubles as the rule-10 regression test | d9a50fa |
| 2026-09-09 | v5.3.4 cut: clean rebuild, full gate 5/5 SUCCESS (drv.a6, tpv8, tpv10, tpv104, tpv1053d) | v5.3.4 |
| 2026-09-09 | Built `testsys/` tiered test system (unit/regression/e2e, single entry point `testsys/run.py`); fixed check.test.py to gate on exit code and activate the dormant `.nc` comparison (items 3, 4); `testAll.py` now preserves `test/` as `test.prev/` instead of deleting it (item 5); `testCreateNewcase.py` migrated into `testsys/regression/` | 34691db, ae4162c, 17b9474, f847b3b |
| 2026-09-09 | v5.3.5 cut: fixed testsys/e2e/run_e2e.py (missing plotRuptureDynamics step), check.test.py's compare_nc_files (bit-exact `f1.identical(f2)` masquerading as a metadata check, false-failing on ordinary MPI floating-point non-determinism), scripts/lib.py (`sys.exit()` NameError in B2/B3), and install-eqdyna.sh (`chmod -R 755 scripts` rule-13 violation on every build, item 6); added install-eqdyna.sh to test_create_newcase.py's exec-bit guard; full gate 33/33 unit, 2/2 regression, 5/5 e2e cases (15/15 check.test.py comparisons) | v5.3.5 |
| 2026-09-13 | v5.3.6 cut: fixed PML region-14 cross-axis bug (item 8) and fric_tp_h=0 thermal-pressurization bug (item 12), both verified against independent oracles (static bug-pattern guard; clean-room v5.2.0 rebuild) and both reference-regenerated; landed the output-neutral refactor campaign (misc/ deletion, tabs->spaces, FRIC_SLOT_* map, func_lib.f90 extraction, MPI4arn/NodalQuant consolidation, Python post-processing dedup) verified bit-gated by clean rebuild + full e2e rerun; added vectorized python/ port (NumPy+JAX), parity claim flagged as spike-only and not yet in-tree-reproducible (new item 14); full gate green: unit SUCCESS, 3/3 regression, 5/5 e2e cases (15/15 check.test.py comparisons), ubuntu 22.04/gfortran 11.4.0/OpenMPI 4.1.1/netCDF 4.8.1 | v5.3.6 |
