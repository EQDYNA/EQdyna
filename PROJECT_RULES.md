# EQdyna Project Rules

Index — read this list first; jump to a rule only when it's load-bearing.

1. Minimal changes; no unnecessary new files.
2. No silent fallbacks, swallowed errors, or placeholder data.
3. Gate every stage; pass before moving on.
4. Only fresh runs are evidence.
5. One calibrated definition of "pass" — never invent a metric.
6. Every performance number carries its provenance.
7. Reference data is read-only.
8. Never delete evidence unless the result is a confirmed pass.
9. Cheap targeted check before expensive run.
10. Every bug gets a regression test before the fix ships.
11. Docs move with the code, in the same change.
12. Untracked build artifacts never accumulate in the working tree.
13. File permission changes are reviewed individually, never bulk-applied.
14. A living status board, re-checked on a schedule.

---

## 1. Minimal changes; no unnecessary new files

Make the smallest edit that solves the problem. New files are fine when they
carry their own weight; fold content into an existing file when one already
owns that job. Never refactor or rename unrelated code in the same change,
and never leave the old version sitting alongside the new one.

**Rationale**: `misc/` already holds superseded originals (`case.setup.old`,
`makefile.old`, `install-eqdyna-ls6.sh`, `contm.f90`, `ku.f90`, `qdckd.f90`,
`qdcshl.f90`, `qdct2.f90`) next to the current names (`library.f90`,
`calcElemKU.f90`, `calcGlobalShapeFunc.f90`, `calcB.f90`). Every one of those
pairs is a rename or refactor that kept the old file instead of deleting it,
which is how a two-language, two-directory-tree repo like this one
accumulates dead weight nobody remembers to remove.

**How to apply**: when a Fortran/Python/MATLAB file is renamed or superseded,
delete the old one in the same commit — do not park it in `misc/` or anywhere
else "for reference."

---

## 2. No silent fallbacks, swallowed errors, or placeholder data

A missing input file, undefined parameter, or failed comparison fails loudly.
No substituted default, no skipped step, no result printed but not acted on.

**Rationale**: `check.test.py:30-38` calls `xr.testing.assert_allclose(f1,f2)`
inside a `try/except AssertionError`, prints the exception, and moves on —
the mismatch is visible in the log but never turns into a failing exit code.

**How to apply**: any comparison that can fail must either raise past the
caller or set a variable checked by `testAll.py`'s exit path — printing to
stdout is not a gate (see rule 3).

---

## 3. Gate every stage; pass before moving on

Named commands, named pass criteria:

- Build: `cd src && make` (via `./install-eqdyna.sh -m <machine>`) must exit 0.
- Test: `python3 testAll.py` then `python3 check.test.py` — pass means every
  printed line for every `testid` in `testNameList.nameList` reads `SUCCESS`,
  with **no** `FAIL` string, for every file in `fileNameList`.

**Rationale**: as written, `check.test.py` never converts a `FAIL` string or
a swallowed `AssertionError` (rule 2) into a non-zero process exit — a CI job
or a human skimming a green checkmark has no mechanical way to know a
comparison failed short of reading every printed line.

**How to apply**: treat "the script ran" and "the script passed" as two
different questions until `check.test.py` itself distinguishes them.

---

## 4. Only fresh runs are evidence

A number in `README.md` or `pastReleaseNotes.md` is a hypothesis until
reproduced on the hardware and git SHA it's attached to.

**Rationale**: README's news section cites "512 cores... 4 hours and 40
minutes" for 50 m TPV36 on Lonestar6 with no git SHA or compiler/netCDF
version attached — the kind of claim that gets copy-pasted into the next
release note unchanged once the code underneath has moved.

**How to apply**: before repeating a timing or verification claim in a new
release note, rerun it, or mark it explicitly as inherited/unverified.

---

## 5. One calibrated definition of "pass" — never invent a metric

`check.test.py` uses `threshold=1e-3` (absolute and relative) as the sole
numeric tolerance for both `compare_nc_files` and `compare_txt_files`.

**Rationale**: `scripts/compareTwoNc.py` is a second, standalone comparison
utility outside the test harness — if it carries its own tolerance, that
tolerance must be reconciled with `1e-3`, not treated as an independent
"looks close enough" check.

**How to apply**: any new comparison script cites `check.test.py`'s threshold
or explains in writing why it differs; it never reports "pass" against an
ad hoc number.

---

## 6. Every performance number carries its provenance

Machine (`ls6`/`ubuntu`/`macos`/`grace` per `src/makefile`), core count,
compiler (`mpif90`/`mpiifort`), NETCDF_LIB/NETCDF_INC versions, and git SHA
travel with any wall-clock number, in `README.md` or `pastReleaseNotes.md`.

**Rationale**: same TPV36 512-core figure as rule 4 — without the SHA it is
unfalsifiable against the current `src/`.

**How to apply**: `git rev-parse --short HEAD` and `module list`/compiler
`--version` output go in the commit or release note next to the number.

---

## 7. Reference data is read-only

`test.reference.results/` (`test.drv.a6`, `test.meng2023a`, `test.meng2023cb`,
`test.tpv10`, `test.tpv104`, `test.tpv1053d`, `test.tpv8`) is ground truth.
Nothing writes through it — not `testAll.py`, not a debug run, not manually.

**Rationale**: `check.test.py` sets `refRoot='test.reference.results'` and
`testRoot='test'` as two distinct trees precisely so a run's scratch output
never lands on top of the golden copy; any change that collapses that
separation (e.g. pointing `testRoot` at `refRoot` "just to check something")
destroys the only baseline the suite has.

**How to apply**: regenerating a reference result is a deliberate, reviewed
commit to `test.reference.results/`, never a side effect of running tests.

---

## 8. Never delete evidence unless the result is a confirmed pass

`testAll.py` runs `os.system('rm -rf test')` before every invocation,
unconditionally, before the new run's results have even been compared.

**Rationale**: if a run fails, the failing `test/` tree — the only record of
what actually happened — is deleted the next time anyone runs `testAll.py`,
including a `testAll.py` invoked just to "try again."

**How to apply**: do not re-run `testAll.py` over a failed `test/` directory
without first copying it aside (`cp -r test test.failed.<date>`) if the
failure hasn't been root-caused yet.

---

## 9. Cheap targeted check before expensive run

`testNameList.py` already sequences small, fast, low-core cases
(`test.drv.a6`, `test.tpv8`, `test.tpv10`, `test.tpv104`, `test.tpv1053d`,
4 cores each) ahead of any HPC-scale allocation.

**Rationale**: a TPV36-class run at 512 cores on Lonestar6 costs hours of
allocation; a mesh, friction-law, or I/O regression is almost always visible
in one of the existing 4-core cases first.

**How to apply**: `python3 testAll.py` (or the specific failing case via
`create.newcase` + `case.setup`) must pass locally before requesting a
large-core-count HPC job for the same change.

---

## 10. Every bug gets a regression test before the fix ships

**Rationale**: the test harness already has the shape for this — a case
under `case_input/`, a matching entry in `test/`, and a golden result under
`test.reference.results/`, driven by `testNameList.nameList` /
`coreNumList`.

**How to apply**: a fix to `src/*.f90` that isn't covered by an existing case
adds one: new `case_input/<case>/`, a run captured in
`test.reference.results/<case>/`, and an entry in `testNameList.py` — landed
in the same change as the fix, not a follow-up.

---

## 11. Docs move with the code, in the same change

**Rationale**: `misc/` is exactly what happens when this rule is skipped —
`ku.f90`, `qdckd.f90`, `qdct2.f90`, `qdct3.f90` are pre-rename names from the
"rename file and subroutine names for clarification" refactor noted in
`README.md`; nothing forces a reader landing on the old name back to the
current one.

**How to apply**: a file or subroutine rename updates every reference in
`README.md`, `pastReleaseNotes.md`, and the `src/makefile` dependency list in
the same commit — and deletes the old file (rule 1) rather than archiving it.

---

## 12. Untracked build artifacts never accumulate in the working tree

`bin/`, `__pycache__/`, `*.mod` (e.g. `src/globalvar.mod`), and `*.pyc` (e.g.
`testNameList.pyc`) are build/run byproducts, not source. They must be listed
in `.gitignore` and must never show up as untracked in `git status`.

**Rationale**: as of 2026-09-09, `git status` on this repo shows `bin/`,
`__pycache__/`, `src/globalvar.mod`, and `testNameList.pyc` untracked — a
recurring pattern, not a one-off, that clutters every status check and risks
an accidental `git add -A` sweeping compiled binaries into history.

**How to apply**: check with `git status --porcelain` before any commit —
none of `bin/`, `__pycache__/`, `*.mod`, `*.pyc` may appear as `??`. If a new
build/run byproduct pattern shows up, it goes into `.gitignore`, not into a
one-off manual `rm`.

---

## 13. File permission changes are reviewed individually, never bulk-applied

`chmod` is applied to the specific entry points that must be executable, not
recursively across a tracked directory.

**Rationale**: `README.md`'s own Installation section instructs
`chmod -R 755 install-eqdyna.sh scripts` — a recursive chmod over a tracked
directory. Run against an already-cloned working tree, this previously
flipped the mode bit on every file under `scripts/` (29 tracked files in one
incident), producing 29 mode-only entries in `git status`/`git diff` that
were indistinguishable from real content changes and had to be triaged by
hand before anyone could tell which files, if any, had actually changed.

**How to apply**: `chmod` only the files that are invoked directly
(`install-eqdyna.sh`, and individual scripts under `scripts/` such as
`case.setup` or `create.newcase` — not the whole directory). Before
committing, run `git diff --summary` and confirm it prints no
`mode change` lines that weren't intended.

---

## 14. A living status board, re-checked on a schedule

`pathway_forward.md` is the one file recording every open issue and
standing claim for this repo, each with a re-check interval, the date it was
last checked, and the exact command whose output was read. It is present
tense — history belongs in `pastReleaseNotes.md`, not here.

**Rationale**: `pastReleaseNotes.md` and README's news entries are
append-only and accretive by design; neither one states whether a past claim
still holds today.

**How to apply**: before citing a "this already works" or "already fixed"
claim from `README.md` or `pastReleaseNotes.md`, check `pathway_forward.md`
first, and re-run the cited command rather than trusting the recorded line.
