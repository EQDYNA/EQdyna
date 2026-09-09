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

**Deferred**: none.

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
