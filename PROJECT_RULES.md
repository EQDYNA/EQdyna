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
15. Releases follow the documented workflow, notes lead the README.
16. Test what you commit, not what is in your working tree.
17. Reviving or adding a TPV benchmark.

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

**This includes a check that degrades instead of refusing.** If a validator
cannot perform a check — a required parameter is absent, the grid is too small
for the stencil, the data does not support the test — it must FAIL, not run a
weaker version and report a pass, and not emit a warning and continue. A pass
must mean the check ran. "I could not check this" and "this is fine" are
different answers and must produce different exit codes.

**Rationale**: two separate incidents.
`check.test.py:30-38` called `xr.testing.assert_allclose(f1,f2)` inside a
`try/except AssertionError`, printed the exception, and moved on — the mismatch
was visible in the log but never turned into a failing exit code.
Then in v5.6.0, `lib.validateFaultRoughGeometry` did it twice more: it
substituted `dx` for a missing `dy` in the per-cell offset check (silently
changing the units of the quantity compared against 1.0), and it fell back to a
weaker derivative bound on grids too small for the 4-point stencil — a bound
derived from the very column under test, so a corrupted derivative set its own
tolerance. The `dy` substitution had been masking its own test coverage: when
the fallback was removed, a unit test failed and revealed that the fixture had
never passed `dy` at all, so every offset assertion in the suite had been
exercising the fallback rather than the check.

**How to apply**: any comparison that can fail must either raise past the
caller or set a variable checked by the runner's exit path (then `testAll.py`,
today `testsys/run.py`) — printing to
stdout is not a gate (see rule 3). A validator that cannot evaluate a check
appends a problem, never a warning. Reserve warnings for input that is legal
but unusual (a rough surface, a steep dip) — never for a check that did not
run. When you remove a fallback, re-run the tests that covered it: if any now
fail, they were testing the fallback.

---

## 3. Gate every stage; pass before moving on

Named commands, named pass criteria:

- Build: `cd src/fortran && make` (via `./install-eqdyna.sh -m <machine>`) must exit 0.
- Test: `python3 testsys/run.py all` — the sweep (8 cases x 3 backends).
  Pass means every
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

**How to apply**: any new comparison script cites the single calibrated
threshold (then `check.test.py`'s, today `testsys/matrix.py`'s `THRESHOLD`
and per-case `CASE_BOUND`)
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
`test.tpv10`, `test.tpv104`, `test.tpv1053d`, `test.tpv29`, `test.tpv8` — one
`frt.canonical.txt` per case) is ground truth. Nothing writes through it —
not the sweep, not a debug run, not manually.

**Rationale**: `check.test.py` sets `refRoot='test.reference.results'` and
`testRoot='test'` as two distinct trees precisely so a run's scratch output
never lands on top of the golden copy; any change that collapses that
separation (e.g. pointing `testRoot` at `refRoot` "just to check something")
destroys the only baseline the suite has.

**How to apply**: regenerating a reference result is a deliberate, reviewed
commit to `test.reference.results/`, never a side effect of running tests.

---

## 8. Never delete evidence unless the result is a confirmed pass

A failed run's output is the only record of what actually happened. Deleting
it to "try again" destroys the evidence before anyone has read it.

**Rationale**: the original offender was `testAll.py`, which ran
`os.system('rm -rf test')` unconditionally at the start of every invocation —
before the new run's results had even been compared. That file is gone, and
the behaviour is now the opposite: `testsys/e2e/run_e2e.py` ROTATES the
previous `test/` to `test.prev/` instead of deleting it, so one level of
history survives by construction. Every python cell's output also stays on
disk under `test/`; the deleted accept tier used to run them in a tempfile
directory and remove it, so a failure destroyed its own evidence.

**How to apply**: the rotation keeps ONE level. If a failure is not yet
root-caused and you are about to run the sweep twice more, copy the tree
aside first (`cp -r test test.failed.<date>`) — the second run's rotation
will otherwise overwrite `test.prev/`.

---

## 9. Cheap targeted check before expensive run

`testNameList.py` already sequences small, fast, low-core cases
(`test.drv.a6`, `test.tpv8`, `test.tpv10`, `test.tpv104`, `test.tpv1053d`,
`test.meng2023a`, `test.meng2023cb`, `test.tpv29` — 4 ranks each) ahead of any
HPC-scale allocation.

**Rationale**: a TPV36-class run at 512 cores on Lonestar6 costs hours of
allocation; a mesh, friction-law, or I/O regression is almost always visible
in one of the existing 4-core cases first.

**How to apply**: `python3 testsys/run.py unit regression` (seconds), then the
sweep or a single cell of it (`testsys/e2e/run_e2e.py --cases <case>
--backends <backend>`) must pass locally before requesting a
large-core-count HPC job for the same change. `run_e2e.py` enforces
cheap-first structurally: its Gate 1 is
`testsys/regression/test_create_newcase.py`, before any build or any case.

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

---

## 15. Releases follow the documented workflow, notes lead the README

The release workflow, in order:

1. Green gate first (rule 3): `./install-eqdyna.sh -m ubuntu` exits 0 and
   `python3 testsys/run.py all` reports 24/24 cells SUCCESS. Never tag over a
   red gate, and never tag before CI is green on the pushed commit (rule 15).
2. Bump `VERSION`.
3. Release notes: add a `* YYYYMMDD vX.Y.Z release notes` block under a
   `# News in <year>` heading at the TOP of `README.md`, ending with the
   pointer line "For past release notes, please refer to
   pastReleaseNotes.md." Move the previous release's block from `README.md`
   into `pastReleaseNotes.md` under its year heading, dropping the pointer
   line.
4. Add a Tasks-done row to `pathway_forward.md` (rule 14).
5. Commit everything above together.
6. Push the COMMIT and wait for CI to go green. Do not tag yet.
7. Tag `vX.Y.Z` and publish the GitHub Release only after CI is green on that
   commit.

   **Why the tag comes after CI, not before.** A local gate cannot model the
   runner. v5.7.0 was gated green locally through CI's own entry point, tagged,
   published -- and CI went red because a GitHub runner has 7 GB and the
   Python backend needed 10.4 GB for one case. The development box has 64 cores
   and far more memory, so the constraint was invisible to any local run. That
   was the fourth CI red in one sequence with the same shape, each time a local
   green that did not transfer: v5.6.0 *what* was committed (partial git add),
   v5.6.1 *how* the gate was invoked (tiers directly, not install-eqdyna.sh),
   v5.6.2 *what environment* it ran in (scipy undeclared), v5.7.0 *what
   resources* it had. Rule 16 closed the first three by making the local gate
   more faithful. Resources cannot be closed that way -- only by ordering.

   A released tag that points at a red commit is worse than a late tag: the
   Releases page becomes the authoritative wrong answer.
8. Push only on explicit approval from the maintainer.
9. Publish the GitHub Release for the tag
   (`gh release create vX.Y.Z --title vX.Y.Z --notes-file <notes> --latest`)
   so the Releases page always shows the current version -- after step 7.

**Rationale**: v5.3.4 (2026-09-09) was cut with its notes appended to
`pastReleaseNotes.md` instead of leading `README.md`, because the
convention existed only in the files' shape, not as a rule — the release
agent followed the wrong precedent and nothing could catch it.

README notes are USER-FACING: terse one-line bullets (the v5.3.2-era
style); the full technical detail belongs in the GitHub Release body,
the annotated tag message, and `pathway_forward.md`.

**How to apply**: at release time, `head README.md` must show the version
being released; `pastReleaseNotes.md` must contain every prior version and
not the current one.

---

## 16. Test what you commit, not what is in your working tree

After committing, the working tree must contain nothing that the tests
depended on. Before pushing a change whose gate you ran locally, confirm
`git status --porcelain` shows no modified tracked files, or re-run the gate
from a clean checkout of the commit itself (`git stash` / a fresh worktree).

**Rationale**: v5.5.0's fault-on-MPI-boundary fix was gated green locally, and
CI went red on the same commit. The commit contained the regression test and a
"FIXED" note but not `src/meshgen.f90` / `src/eqdyna3d.f90` -- a partial
`git add` left the actual fix uncommitted. The local gate passed because the
working tree had it; CI built the commit, which did not. The test was correct
and caught exactly what it was written to catch. Two commits shipped a claim
the code did not support.

**Second incident, v5.6.0**: the commit was complete, the tree was clean, and
the commit itself was re-verified in a fresh worktree — and CI still went red.
The gate had been run by invoking the tiers directly (`testsys/run.py e2e` with
`EQDYNA_E2E_BIN=src/eqdyna`, and `make eqdyna` by name), while CI runs
`./install-eqdyna.sh -m ubuntu` and then `testsys/run.py unit regression e2e-ci`.
A change to
`src/makefile` made a BARE `make` stop producing a binary, which only the
install path exercises. Testing the right commit is not enough if you invoke it
differently than CI does.

**How to apply**: `git status --porcelain` after every commit, before every
push; treat any remaining modified tracked file as a reason to stop and check
whether the gate you ran still describes the commit. Prefer `git add <paths>`
with the full list read back from the diff, or `git add -u` scoped to the
directories the change touched.

Before pushing a release, run CI's own entry point, not a convenient subset of
it — read `.github/workflows/test.yml` and reproduce the commands verbatim:

    ./install-eqdyna.sh -m ubuntu
    export EQDYNAROOT=$(pwd); export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
    python3 testsys/run.py unit regression e2e-ci

That last line is what `.github/workflows/test.yml:64` actually runs, and it
is NARROWER than `run.py all`: `e2e-ci` runs `matrix.CI_CELLS`, 10 of the 24
cells, because the rest do not fit a 7 GB runner. Run `run.py all` too — it is
the wider local gate and the one that speaks for the whole table — but do not
mistake it for a reproduction of CI. This rule was itself wrong about this
until 2026-09-16, telling the reader to reproduce CI with a command CI does
not run: the exact substitution it exists to forbid.

Shortcuts like `EQDYNA_E2E_BIN=src/eqdyna` exist to keep a gate from disturbing
a running job; they skip the build-and-install path, so a green run under them
says nothing about it. If you used one, say so, and run the real entry point
before the tag.

---

## 17. Reviving or adding a TPV benchmark

TPV29 established this recipe the expensive way. Follow it in order; each step
exists because skipping it cost something.

**1. Fetch the official spec first, and cite it by name and part.** Put it in
`scratch/specs/`. Every resolution, duration and station list comes from there.
If the spec does not state a number, it does not get invented — `test.tpv1053d`
has no spec element size, so `testsys/e2e/full_specs.py` carries an EXCLUDED
entry saying exactly that instead of a plausible guess.

**2. Decimate geometry, never interpolate.** A supplied fault surface is a
fixed sampling of ONE random realisation. Taking every n-th node keeps official
values; interpolating invents detail that is not in the benchmark. Ship the
spec resolution plus a coarser source, declare `par.faultGeometrySourceDx` and
`par.faultGeometrySourceAvailableDx`, and let
`lib.requireFaultGeometryResolution` refuse any dx that is finer than, or not a
multiple of, the source. TPV29 ships 100 m (3.5 MB) and 50 m (14 MB), both
exact decimations of the official 25 m file, so a clean checkout needs no
download.

**3. Check the CODE has that TPV's branch before trusting a run.** This is the
step that cost the most. `swtwNucleation` branched on `TPV in {201,36,37}` and
omitted 29 — the case the smoothed forced-rupture formula comes FROM — so
`case_input/test.tpv29` had been declaring `par.tpv = 36` to reach its own
physics. The misdeclaration hid a real gap: the Python port implemented only
the degenerate branch and nucleated **0 of 3321** fault nodes against a
reference of 2974, failing at 7.10e7 on both backends at the identical value.
Grep the source for the TPV number you are adding. If the case has to
impersonate another TPV to work, that is a bug in the code, not a
configuration trick.

**4. Gate COARSE, freeze ONE reference.** Gate at a dx that runs in minutes at
4 ranks (TPV29: 500 m), not at spec resolution. Freeze exactly one
`frt.canonical.txt` (`python3 -m testsys.frt_canonical <case_dir>`) — it is
decomposition-independent, so the same artifact serves Fortran at any rank
count and the serial Python backends. Then add the case to BOTH
`testNameList.py` and `testsys/matrix.py` (`CASE_BOUND` **and** `GATE`);
matrix.py fails at import if either is missing, which is deliberate.

**5. Record the spec-resolution tier without running it.** Add the
`full_specs.py` entry (dx, term, nx/ny/nz, citation) so the run is one command
away. Actually performing it is a scheduling decision — TPV29's 50 m run was
deliberately NOT done because 500/200/100 m already converged (0.076 / 0.028 /
0.006 s median rupture-time difference) and the 100 m run was the real
cross-code validation.

**6. Validate against something independent, as a SCRIPT.** TPV29 was checked
against EQdyna's own 2015 SCEC submission at matched 100 m: 0.006 s median
rupture time, 0.22% median final slip over 24 on-fault stations. Those numbers
are prose-only and their baseline (`scec_archive/`) is gitignored, so nothing
in-tree reproduces them — that is pathway item 28, and it is the one part of
the TPV29 work that was done wrong. Write the comparison as a committed script
from the start (rule 4).

**7. Run the FULL sweep, all three backends -- and ALL THREE MUST PASS.**
A new case is 3 cells, not 1.

**Supporting a TPV, or any problem, means supporting it on every backend**
(owner, 2026-09-16). A case is not added with its Python columns declared
UNSUPPORTED to be filled in later: that ships a benchmark the Fortran can run
and the port cannot, and the sweep's whole value is that the two
implementations check each other. If the port lacks a feature the case needs,
PORT THE FEATURE FIRST, then add the case.

`UNSUPPORTED` in `matrix.py` remains the honest way to record a gap that
ALREADY exists on a case already gated -- it is not a runway for new ones. The
friclaw-2 gap is the worked example: test.meng2023a and test.meng2023cb sat
declared-unsupported for months, and closing it turned out to be ~20 lines of
`time_weak`. Had this rule been in force, that gap would never have opened.

**How to apply**: `python3 testsys/e2e/run_e2e.py --cases test.tpvNN` for the
new case alone while iterating, then `python3 testsys/run.py all` before
committing the reference.
