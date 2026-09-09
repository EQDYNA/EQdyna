# EQdyna Status Board

Present-tense record of open issues and standing claims for this repo. See
`PROJECT_RULES.md` for the rules each item enforces. History (what shipped,
when) lives in `pastReleaseNotes.md` and `README.md` — this file does not
duplicate it and is not append-only; items are updated or closed in place.

A blank "Last checked" means never audited and stays blank until someone
actually runs the command.

| # | Item | Rule | Surface | Re-check interval | Last checked | Command |
|---|---|---|---|---|---|---|
| 1 | `bin/`, `__pycache__/`, `src/globalvar.mod`, `testNameList.pyc` untracked in working tree | 12 | repo root, `src/` | every commit | 2026-09-09 | `git status --porcelain` |
| 2 | `.gitignore` has an uncommitted edit adding `bin/`, `__pycache__/`, `*.mod`, `*.pyc`, `scratch/`; not yet landed | 12 | `.gitignore` | until committed | 2026-09-09 | `git diff .gitignore` |
| 3 | `check.test.py` swallows `AssertionError` from `xr.testing.assert_allclose` (prints, doesn't fail) and never sets a non-zero exit on a printed `FAIL` | 2, 3 | `check.test.py:30-38`, `:56-66` | each change to test infra | | `python3 check.test.py; echo exit=$?` |
| 4 | `compare_nc_files` call in `check.test.py` is commented out — `.nc` files listed in `fileNameList` are never actually diffed | 3 | `check.test.py:64` | each change to test infra | | `grep -n compare_nc_files check.test.py` |
| 5 | `testAll.py` runs `rm -rf test` unconditionally before every invocation, deleting prior run evidence even on a prior failure | 8 | `testAll.py` (root) | each change to test infra | | `grep -n "rm -rf test" testAll.py` |
| 6 | README Installation section still instructs a recursive `chmod -R 755 install-eqdyna.sh scripts` | 13 | `README.md` | each README edit | | `grep -n "chmod -R" README.md` |

**Deferred**: none.
