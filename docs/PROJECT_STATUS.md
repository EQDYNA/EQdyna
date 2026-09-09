# EQdyna Status Board

Present-tense record of open issues and standing claims for this repo. See
`PROJECT_RULES.md` for the rules each item enforces. History (what shipped,
when) lives in `pastReleaseNotes.md` and `README.md` — this file does not
duplicate it and is not append-only; items are updated or closed in place.

A blank "Last checked" means never audited and stays blank until someone
actually runs the command.

| # | Item | Rule | Surface | Re-check interval | Last checked | Command |
|---|---|---|---|---|---|---|
| 1 | ~~Untracked build artifacts~~ resolved by `.gitignore` update (8c40871); recheck each release | 12 | repo root, `src/` | every release | 2026-09-09 | `git status --porcelain` |
| 3 | `check.test.py` swallows `AssertionError` from `xr.testing.assert_allclose` (prints, doesn't fail) and never sets a non-zero exit on a printed `FAIL` | 2, 3 | `check.test.py:30-38`, `:56-66` | each change to test infra | | `python3 check.test.py; echo exit=$?` |
| 4 | `compare_nc_files` call in `check.test.py` is commented out — `.nc` files listed in `fileNameList` are never actually diffed | 3 | `check.test.py:64` | each change to test infra | | `grep -n compare_nc_files check.test.py` |
| 5 | `testAll.py` runs `rm -rf test` unconditionally before every invocation, deleting prior run evidence even on a prior failure | 8 | `testAll.py` (root) | each change to test infra | | `grep -n "rm -rf test" testAll.py` |
| 6 | README Installation section still instructs a recursive `chmod -R 755 install-eqdyna.sh scripts` | 13 | `README.md` | each README edit | | `grep -n "chmod -R" README.md` |

| 7 | TO-DO: `netcdf_read_on_fault_eqdyna` hardcodes fault index `1` instead of loop var `ift` in all `fric(...)` assignments, so on-fault netcdf input is never applied to faults 2+ (fault-N data overwrites fault 1); `ii,jj` grid indices are also unchecked against bounds for faults 2+ (built from fault-1 extent). The restart routine below it (`netcdf_io.f90:159-179`) shows the correct `ift` pattern. Bites only `ntotft >= 2` — no current test case exercises it; fix needs a multi-fault regression case per rule 10. Does NOT explain the tpv1053d reference mismatch (that case runs `ntotft=1`; its golden reference is stale, carrying an uninitialized denormal). | 2, 10 | `src/netcdf_io.f90:76-105` | until fixed | 2026-09-09 | `grep -n "i, 1) = on_fault_vars" src/netcdf_io.f90` |

**Deferred**: none.
