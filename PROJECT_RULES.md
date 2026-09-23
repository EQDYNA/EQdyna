# EQdyna Project Rules

Index — read this list first; jump to a rule only when it's load-bearing.

1. Minimal changes; no unnecessary new files.
2. No silent fallbacks, swallowed errors, or placeholder data.
2a. A scripted edit to a tracked document asserts on its shape, not a substring of it.
3. Gate every stage; pass before moving on.
4. Only fresh runs are evidence.
4a. A proposed CAUSE is falsified by the outcome curve, not by the defect it predicts.
4b. An exclusion names the metric that produced it; a metric blind to a mechanism cannot exclude it.
4c. A gate downgrade names the geometry class its evidence came from, and clears no class that evidence cannot exercise.
5. One calibrated definition of "pass" — never invent a metric.
5a. A provably injective relabelling is gated at bit-identity, not at the case bound.
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
15a. A pre-tag CI check on the exact SHA is required before `git tag`, and it must be mechanical, not remembered.
15b. The local sweep and CI gate different failure classes; CI's release run re-verifies only what a local sweep structurally cannot.
16. Test what you commit, not what is in your working tree.
17. Reviving or adding a TPV benchmark.
18. A refactor that couples two previously-independent artifacts must say so.
19. A shared mutable `*_last.*` artifact is not evidence until pinned to a commit.
20. Heavy runs are launched detached and polled by artifact, never by process.
20a. A tool that writes its results file only at the end is a partial-loss hazard; prefer several smaller invocations that each land their own artifact.
20b. A heavy run's artifacts are not landed until they are committed, and a run finishing does not mean anyone is still there to commit them.
20c. A completed detached run is indistinguishable from one still in flight until something reads its artifact and updates the board — and until then its worktree is not reapable by default.
21. An agent's write surface is its own worktree; the main checkout belongs to the conductor.
21a. A gate sweep runs in its own worktree; the shared `test/` tree has no lock, and a collision reads as a solver failure.

Count, stated so a heading-shape grep does not undercount it again (that
undercount happened twice in one night, 2026-09-21/22): 21 numbered rules
(1-21) plus eleven lettered sub-rules (2a, 4a, 4b, 4c, 5a, 15a, 15b, 20a, 20b,
20c, 21a) — 32 `## ` headings
total. Verify: `grep -c '^## ' PROJECT_RULES.md` reads 32;
`grep -c '^## [0-9]*\. ' PROJECT_RULES.md` (numbered rules only, no letter
suffix) reads 21.
A count that greps only `^## [0-9]` and calls it "the rules" will silently
drop every lettered sub-rule — read this index's own list, don't re-derive
the count from heading shape alone.

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

## 2a. A scripted edit to a tracked document asserts on its shape, not a substring of it

When a script (a regex, a `sed`, an automated move-this-block edit) rewrites a
tracked document, the check that confirms the edit worked must assert on the
document's SHAPE afterward — line count, required sections present, nothing
missing — not on whether some expected substring still appears somewhere in
it. A substring check passes even when the edit deleted everything else.

**Rationale**: 2026-09-21, a scripted regex meant to move a `README.md`
release-notes block into `pastReleaseNotes.md` had a tail pattern that matched
to end-of-file. `README.md` went from 317 lines to 8, with its entire
remaining body appended into `pastReleaseNotes.md`. The editor's own
post-edit assertion passed anyway, because it tested for a substring (the
moved block's own text landing in the right file) rather than the shape of
what was left behind — a check that would have failed instantly on a line-count
or required-section assertion. `test_readme_commands.py` and
`test_stop_exit_status.py` caught it on the next gate run; that is the test
suite working as intended, not this rule's substitute.

**How to apply**: after a scripted edit to any tracked document, check the
result's shape before trusting it — line count against a sane bound, and
presence of every section the document is required to carry — in addition to,
never instead of, confirming the intended content landed. This is a norm for
how an agent verifies its own edits, not a repo-wide mechanical gate (no
generic script-edit checker exists); the tier-1 backstop for this class of
mistake is each document's own drift tests (`test_readme_commands.py`,
`test_stop_exit_status.py`, `check_pathway_tasks_done_row` and similar),
which this rule does not replace.

---

## 3. Gate every stage; pass before moving on

Named commands, named pass criteria:

- Build: `cd src/fortran && make` (via `./install-eqdyna.sh -m <machine>`) must exit 0.
- Test: `python3 testsys/run.py all` — the sweep (10 cases x 3 backends = 30
  cells, as of 2026-09-21; the table itself is `testsys/matrix.py`, which is
  authoritative over this count).
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

## 4a. A proposed CAUSE is falsified by the outcome curve, not by the defect it predicts

A causal claim — "X is why this is slow" — is settled by measuring the OUTCOME
the claim is about (the scaling curve, the per-step cost, the pass rate) with X
removed, not by measuring the defect X names. A before/after on the defect
metric (a work-spread ratio, a straggler factor, a zero-owner count) shows only
that the defect existed and was repaired; it carries no information about
whether the defect mattered. Until the outcome number moves by more than that
machine's own run-to-run variation, the mechanism stays a hypothesis and the
repair is written up as "defect fixed, consequence unmeasured."

**Rationale**: the defect metric is cheap and the outcome metric expensive,
which is exactly the pressure that makes a team stop at the cheap one. Two
numbers can both be measured correctly and be unrelated.

**Incident (2026-09-22)**: `pathway_forward.md` item 60 asserted decomposition
imbalance as the cause of the jax-MPI 32-rank plateau on evidence that was
entirely defect-side — four of 32 ranks owning zero interior elements, barrier
wait 105 ms of a 172.6 ms step. `PML_WEIGHT = 3.0`
(`src/python/eqdyna/MPI4NodalQuant.py:85`, documented in source as a guess) was
then measured at ~0.6, five times too high; recutting took zero-`Ei` ranks 4 ->
0 and predicted work spread 3.350x -> 1.004x. The defect was real and was
repaired completely. Per-step cost at 32 ranks went 134.12 -> 132.20 ms: 1.4%,
inside this box's run-to-run variation, with two mid-range points moving the
wrong way. The plateau is still unexplained, and the board stated a disproved
mechanism as fact until someone measured the curve. Second occurrence of this
shape in the repo — a fix that worked and bought nothing (see `CLAUDE.md`,
"Measure, do not infer").

**How to apply**: name the outcome metric and its run-to-run spread on that
machine BEFORE removing the suspected cause, and quote both the defect delta
and the outcome delta when reporting. One measurement per point on a shared box
supports neither an improvement nor a regression claim — say "unresolved". A
candidate cause supported only by defect-side evidence goes on the board as a
hypothesis, never into a doc as the cause.

---

## 4b. An exclusion names the metric that produced it; a metric blind to a mechanism cannot exclude it

A written exclusion — "contention is ruled out", "transport is not it", "compile
is not a factor" — carries, in the same sentence, the METRIC that produced it
and the value that metric read. An exclusion with no metric attached cannot be
re-audited: the next session inherits it as a closed door and spends its budget
everywhere else. And a metric that cannot OBSERVE a mechanism cannot exclude
it, however clean its number. Two instruments in this repo are known blind and
must be named as such wherever they appear in an exclusion:

- **`EFFECTIVE_CORES` is blind to memory stalls.** A core spinning on cache
  misses still bills one cpu-second per wall second and reads 1.00. It measures
  whether a cpu was HANDED to the process, not whether the process made
  progress.
- **A per-rank WALL-time spread under a barrier is an identity, not a
  measurement.** Every rank leaves the same barrier at the same instant, so a
  spread of 1.00x is what the loop's structure guarantees regardless of what
  the ranks did. The observable is the per-rank COMPUTE spread underneath it.

**Rationale**: an exclusion is a permission to stop looking — the cheapest
claim in a campaign to make and the most expensive to be wrong about. Rule 4a
governs a proposed CAUSE and asks for the outcome curve; this sub-rule governs
the NEGATIVE claim and asks which instrument produced it. The two fail the same
way: a number measured correctly, about the wrong quantity.

**Incident (2026-09-22)**: the 32-rank plateau mission was briefed with a
"ruled out with numbers, do not re-run" list that excluded contention on
exactly those two instruments — a per-rank wall spread of 1.00x and
`EFFECTIVE_CORES` 1.00 on 31 of 32 ranks — and `pathway_forward.md` item 60
carried the exclusion as fact ("so it is NOT contention"). Measured, contention
is the largest single term in the step: the compute spread underneath the
barrier is **2.90x** (33.73-97.94 ms) against a work spread of 1.41x, and 32
concurrent zero-communication processes reproduce the in-situ per-rank cost to
**3.7%** (94.30 vs 97.94 ms) — roughly 45 ms of a ~124 ms step, behind a door
the campaign had closed. Neither instrument could have seen it. Three missions
over two days were aimed at the mechanisms that were still open.

**How to apply**: write an exclusion as "X is excluded by <metric>, which read
<value>", never as "X is ruled out". Before inheriting one, ask whether that
metric can observe the mechanism — if it cannot, the entry is not an exclusion
and goes back on the board as open. When handing a ruled-out list to another
agent, hand the metric with each entry and ask for the contradiction rather
than the agreement. Tier: a norm for how claims are written, not mechanically
gated — no script can tell a metric from the mechanism it is blind to. The
nearest mechanical backstop is the board's own `Command` column (rule 14): an
exclusion whose command still runs is at least re-checkable.

---

## 4c. A gate downgrade names the geometry class its evidence came from, and clears no class that evidence cannot exercise

A decision to DOWNGRADE a gate — an abort turned into a warning or a NOTICE, a
refusal turned into a pass, a bound loosened — states, in the same place the
downgrade is recorded, the GEOMETRY CLASS of the evidence behind it (planar or
rough; one decomposition; one dip; one fault count). A downgrade supported by
one class covers that class only. It may not be read, later or by anyone else,
as having cleared a class its evidence is structurally unable to exercise: if
the fixture cannot produce the failure the check was guarding, its clean result
is not evidence about that failure.

This sharpens 4b for the specific case where the blind instrument is the
FIXTURE rather than the metric. 4b asks which metric produced an exclusion;
this asks which geometry the fixture had. Both fail identically — a number
measured correctly, about the wrong thing — and a downgrade is the most
expensive place for it, because the check that would have caught the uncovered
class is the thing being removed.

**Rationale**: a check's whole value is on the inputs that break it. A
downgrade decided from an input that cannot break it converts "we never tested
this class" into "this class is fine", permanently and silently, because the
guard that would have said otherwise is now a NOTICE nobody reads.

**Incident (v5.8.2, flagged retroactively 2026-09-22)**: `checkFaultMPIAlignment`
was downgraded from abort to NOTICE (`src/fortran/errorCodes.f90:85`,
`ERR_MPI_FAULT_ALIGNMENT = 51`) on the evidence of
`testsys/regression/test_dipping_fault_y_split.py` — tpv36 at (2,2,1) vs
(2,1,2), traction ratio 1.000000 at 3416 nodes. That fixture's fault is
PLANAR in grid y: it does not wander across the y rank boundary, so it cannot
exercise the mechanism the alignment check exists for on a ROUGH fault, where
`fltxyz`'s y-extent reads (0,0) while the physical surface crosses y=0 (pathway
item 66, `src/fortran/meshgen.f90:426`). The downgrade was recorded as if the
check had been cleared; what it had actually been cleared for was one geometry
class, and the class it was blind to is a wrong-answer path (a
mis-classification doubles `arn` and halves every on-fault traction on those
nodes) that still ships behind a NOTICE today.

**How to apply**: write a downgrade as "downgraded on <fixture/run>, which
exercised <geometry class>", and name in the same sentence the classes NOT
exercised. If a named class is one the check was written for, the downgrade
does not extend to it — keep the abort for that class, or open a board row
(rule 14) saying it is unguarded and why. When you INHERIT a downgrade, read
the cited fixture's geometry before relying on it, the same way 4b asks you to
read the cited metric. Tier: reviewable, not mechanically gated — nothing in
this repo can classify a fixture's geometry by itself. The nearest mechanical
backstop is that a downgrade must cite a test file by path, so a reviewer has
something concrete to open; a downgrade citing no fixture at all is a finding
on its own.

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

## 5a. A provably injective relabelling is gated at bit-identity, not at the case bound

When a change touches only an INDEX array — a gather/scatter map, a
relabelling of node or element IDs into a new numbering — and the relabelling
is provably injective (each old index maps to exactly one new slot, nothing
merges or drops), the correct pass criterion is exact equality on the
compared output, not the case's numeric `CASE_BOUND`. A relabelling changes
where a value is stored, not what it is; if the output differs by even one
ULP under a pure re-indexing, the map is wrong, and a tolerance-based bound
would let that bug pass as long as it stayed inside the case's ordinary
roundoff budget.

**Rationale (2026-09-22)**: `MPI4NodalQuant`'s rank-local carry
(`pathway_forward.md` item 60/43/48/61, commit `01a1040`) renumbers every
index array into the rank-local `touched`/`my_eq` sets the restriction had
already computed — a pure relabelling, no arithmetic on the physics. It was
gated at bit-identity rather than at `test.tpv8`'s bound, and held: four rank
files (`frt.txt0`-`frt.txt3`) compared byte-for-byte and md5-equal against
master's global-extent solver on the same case and placement, even though the
four ranks now carry four DIFFERENT carry sizes (16.14/17.57/17.53/16.18 MB)
where at global extent all four were byte-identical. Gating that at the
1.0e-08 case bound instead would still have passed and would have hidden an
off-by-one in the relabelling as long as it landed inside roundoff.

**How to apply**: before gating a relabelling-only change, state the
invariant that makes it injective (an explicit index-space partition, a
`KeyError`-on-unclassified check, or equivalent), and compare the affected
artifact byte-for-byte (`md5sum` or equivalent) rather than through
`testsys/compare.py`'s bound. Cite both: the exact-equality result AND the
case's own bound run unchanged, so a reader can tell a relabelling proof from
an ordinary physics-preserving change.

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

`test.reference.results/` is ground truth — one `frt.canonical.txt` per case.
Nothing writes through it — not the sweep, not a debug run, not manually.
As of 2026-09-21 it holds 11 case directories, and the enumeration below is a
reader's aid, not the authority: `ls test.reference.results/` is.

- The 10 GATED cases, every one of them covered by this rule: `test.drv.a6`,
  `test.meng2023a`, `test.meng2023cb`, `test.tpv10`, `test.tpv104`,
  `test.tpv1053d`, `test.tpv29`, `test.tpv8`, and — added when they were
  gated in v5.9.0/v5.10.0 and read-only on exactly the same terms as the
  other eight — `test.tpv36` and `test.tpv37`.
- `test.tpv30` is an eleventh committed reference belonging to a case that is
  DELIBERATELY NOT REGISTERED in `testNameList.py`/`matrix.py` (pathway item
  19(b); the equilibrium question that first blocked it was closed 2026-09-21
  — measured as the fault frictionally failing at t=0, not a force-balance
  defect — but the case stays unregistered while the python-port divergence
  and the rule-17 steps-4/7 registration decision stand). Being ungated does
  not make it writable: it is the frozen artifact
  that investigation's evidence is measured against, so it is read-only for
  the same reason as the rest. Do not regenerate it to close the divergence.

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
`test.meng2023a`, `test.meng2023cb`, `test.tpv29`, `test.tpv36`,
`test.tpv37` — 10 cases, 4 ranks each, in that order) ahead of any HPC-scale
allocation. Read the order from `testNameList.nameList`, not from this list.

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
   `python3 testsys/run.py all` reports every cell in the table SUCCESS —
   31/31 as of 2026-09-21 (10 cases x 3 backends, plus one opt-in
   `python-jax-mpi` cell). The number is
   `len(matrix.CASE_BOUND) * 3 + len(matrix.PY_MPI_RANKS)`, not a constant in
   this rule: a releaser who accepts a green count smaller than the current
   table has accepted a partial sweep. Never tag over a
   red gate, and never tag before CI is green on the pushed commit (rule 15).

   **Corrected 2026-09-21**: item 43 (jax MPI parallelism) added a fourth
   backend axis, `python-jax-mpi`, opt-in per case via `matrix.PY_MPI_RANKS`
   (currently one case, `test.tpv8` at 4 ranks; the other 9 are DECLARED
   UNSUPPORTED for that mode with a recorded reason). The table therefore now
   declares 40 cells, of which 31 run — this is a correction of the FORMULA
   for what the table declares, because the table grew; it does not relax
   what a green gate must cover. This is not a mandate to opt every case into
   the new mode: see CLAUDE.md's "There is ONE test" and rule 17 step 7's own
   note on the distinction between a backend implementation and an optional
   execution mode of one.
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
7. Tag `vX.Y.Z`, push the tag, and publish the GitHub Release as ONE
   uninterrupted action, only after CI is green on the commit (formerly two
   separate steps, 7 and 9 — merged 2026-09-16, see rationale below). Run all
   three with no pause in between, and do not stop to watch CI between the
   tag push and `gh release create`:

       git tag -a vX.Y.Z -m "<release summary>" <sha>
       git push origin vX.Y.Z
       gh release create vX.Y.Z --notes-file <notes-file> --latest --verify-tag

   `--verify-tag` makes `gh release create` use the tag just pushed instead of
   minting its own from the branch tip, which would be a lightweight tag and
   fail `check_tag_is_annotated`
   (`testsys/regression/test_release_complete.py:106-113`).

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

   **Why the tag push and the GitHub Release cannot be two separate steps
   (2026-09-16, formerly step 9).** Tags are pushed refs, so pushing one fires
   its own `push`-triggered workflow run, independent of the branch push that
   already went green. v5.8.2's release commit `bccfb845` has two: run
   `35121095979` (created 16:18:56Z, the branch push) concluded SUCCESS; run
   `35122388271` (created 16:30:56Z, twelve minutes later, the TAG push,
   identical commit) concluded FAILURE. Its only failing check was
   `test_release_complete.py`'s `check_network_side`
   (`testsys/regression/test_release_complete.py:132-139`), which SKIPS while
   no tag exists for `VERSION` and fires the instant one does — the guard
   exists specifically because `gh release create` (old step 9) was skipped on
   both v5.8.0 and v5.8.1, so it correctly refuses to pass a tag with no
   Release behind it. The old ordering (tag at step 7, `gh release create` as
   a separate, manually-remembered step 9) guaranteed a window between the two
   where exactly that condition holds, and the tag push itself schedules a CI
   run that lands inside it. Collapsing the two into one chained command
   removes the window: the guard's network check runs late enough in CI
   (after checkout, build, unit and regression tiers) that the immediately
   following `gh release create` has already landed by the time it executes.
   (Three older release commits — `335e21d`, `238f1ac`, `9950ae1` — each show
   two push-triggered runs one second apart, both FAILURE; that is a separate,
   unexplained duplicate-push artifact, not this mechanism, and was not
   investigated further.)
8. Push only on explicit approval from the maintainer.

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

## 15a. A pre-tag CI check on the exact SHA is required before `git tag`, and it must be mechanical, not remembered

Before running `git tag`, a COMPLETED, successful CI run must already exist
for the exact commit SHA about to be tagged — not the branch tip in general,
not "the parent commit was green", not an inherited belief that CI is green.
This tightens rule 15 step 7 (tag only after CI is green on the commit); it
does not weaken or replace it. `iris-vermeulen`'s `--pre-tag <sha>` guard is
this rule's enforcement: a release is not gated until that guard reports a
completed, successful run for the SHA being tagged, run before `git tag`, not
inferred from a run that happens to complete afterward.

**A commit that touches ONLY files in `.github/workflows/test.yml`'s
`paths-ignore` list cannot trigger its own CI run, by construction, and must
be recognized as such rather than tagged anyway.** Such a commit has to be
tagged either by pointing at an earlier, already-green commit (fine only if
that earlier commit is what actually gets tagged), or by accepting the nearest
ancestor's completed green run as the evidence — the guard's
`--ack-paths-ignored-parent` flag is the sanctioned form of this: it verifies
that ancestor's run and prints the evidence SHA into the record, so the
acceptance is written down rather than remembered — or by accepting that the
only CI evidence for it will be the tag-push run itself, which happens after
the tag already exists — never by asserting, from memory or by pattern-match
to the normal case, that "a release commit always also touches a non-ignored
file so it still triggers CI regardless."

**Rationale**: 2026-09-21, tag `v5.13.1` was pushed at `dfee14d` when the
ONLY CI run that had ever executed against that exact SHA was the run the tag
push itself triggered — there was no completed, successful pre-tag run for
`dfee14d` at the time `git tag` ran. This is a violation of rule 15 step 7 as
written, not a near-miss softened by its mitigating facts, which are real but
do not make it compliant: the parent commit `e9c3fa2` was green 7/7
(`35658746273`), `git diff e9c3fa2 dfee14d` is a single line in
`pathway_forward.md`, and the tag-push run `35667535066` did conclude 7/7
success — so the breach closed without ever producing a red tag, but it was
still a tag pushed ahead of its own evidence. The structural cause outlives
this one mistake: `pathway_forward.md` and `PROJECT_RULES.md` are in
`test.yml`'s `paths-ignore`, so a Tasks-done-row-only commit — exactly the
shape rule 15 step 4 asks for — can never have its own pre-tag CI run by
push alone. `test.yml`'s own comment claiming "a release commit always also
touches non-ignored files (e.g. VERSION), so it still triggers CI regardless"
is true for a normal release commit and false for a follow-up row commit, and
nothing before this rule mechanically noticed the difference.

**How to apply**: run the pre-tag guard against the exact SHA before every
`git tag`. If it reports no completed run for that SHA, stop — do not tag on
the strength of a parent commit's green run, and do not manufacture a trigger
by touching an unrelated non-ignored file just to get CI to run. If the SHA's
own diff is paths-ignore-only, either retarget the tag at the last commit that
did get a real pre-tag run; or re-run the guard with
`--ack-paths-ignored-parent`, which accepts the nearest ancestor commit's
completed green run as the evidence and prints that ancestor's SHA — keep the
printed line with the release record (amended 2026-09-21: the guard has
provided this third path since it landed and this rule omitted it — a rule
that hides a mechanism its own enforcement offers sends the next releaser
down a path the tool already solved; verified working the same day, see
pathway items 50/51); or accept and record — as this rule requires,
not as an afterthought — that the tag's only CI evidence is the tag-push run
itself.

---

## 15b. The local sweep and CI gate different failure classes; CI's release run re-verifies only what a local sweep structurally cannot

Rule 15 step 1 and rule 15a together make every release pay for two nearly
identical physics sweeps, run strictly in series. This sub-rule divides the
two gates by the failure class each one can actually see. It does not touch
15a's core: nothing red is ever pushed, the tag goes last, and it goes only on
a COMPLETED, successful CI conclusion for the exact SHA being tagged.

**The division.**

- **CI is the MERGE gate.** It owns every failure class this development box
  cannot represent: a clean checkout of the COMMIT (rule 16's v5.5.0 partial
  `git add`), the real entry point `./install-eqdyna.sh -m ubuntu` on a bare
  machine (v5.6.1), the declared apt/pip dependency set resolving from nothing
  (v5.6.2, scipy), the 7 GB runner memory ceiling (v5.7.0), ubuntu-22.04's own
  gfortran/mpich/libnetcdf/numpy/jax versions, and the network-side release
  guards that need `GH_TOKEN`. It runs the full `matrix.CI_CELLS` — 25 of the
  30 e2e cells, `testsys/regression/test_ci_workflow_coverage.py` authoritative
  over that count — on every push that changes anything the solver or the
  harness reads.
- **The local sweep is the RELEASE gate for physics.** `python3 testsys/run.py
  all` is a strict SUPERSET in cells (31 of the 40 declared, against CI's 25),
  it runs on the tree being released, and it is the only gate that covers the
  6 cells CI cannot fit on a 7 GB runner.

**What CI's release run may take as established, and the two conditions that
must BOTH hold before it may.** CI's run on a release commit may reduce to the
release profile below only when:

  (a) `git diff --name-only <ancestor>..<sha>` lands entirely inside the
      release-metadata whitelist — `VERSION`, `README.md`,
      `pastReleaseNotes.md`, `pathway_forward.md`, `PROJECT_RULES.md` — so no
      file the build, the solver, the harness, a compset or a reference reads
      has changed since `<ancestor>`; **and**
  (b) `<ancestor>` has a COMPLETED, successful **full-matrix** CI run of its
      own — not a release-profile run, not a parent-of-parent, not an
      inherited belief.

Both are machine-decidable and both are decided by the workflow, never by an
operator at release time. If either fails — one line in `src/`, `testsys/`,
`case_input/`, `test.reference.results/` or `.github/` — the full matrix runs.
There is no third outcome and no judgment call: an operator-chosen reduced CI
run is exactly the silently-skipped step rule 2 forbids.

**The release profile** is `build` (clean checkout plus the real install
entry point), `unit-regression` (`fetch-depth: 0` and `GH_TOKEN` — rule 15's
and 15a's own guards live in that job, including `check_pretag_ci.py`'s
negative test and `test_release_complete.py`'s network-side check), and ONE
smoke cell per backend (`test.tpv8` x fortran, python-numpy, python-jax). Not
zero e2e cells: the smoke cells are what proves the runner's toolchain still
executes a case end to end from a clean checkout, which `build` alone does not
show. A reduced run is still a run — 15a's pre-tag guard is still executed
against the exact SHA and must still report completed and successful before
`git tag`.

**Rationale**: the two gates were built for different reasons and grew into
each other. Measured 2026-09-22: `run.py all` is 2904.6 s on a quiet box
(9413.8 s contended, per the board's tenancy note) and CI is ~54 min, of which
25 cells are cells the local sweep has just run on a tree whose only difference
from an already-merge-gated parent is release metadata. CI's own wall clock is
floored by one cell that cannot be split — `test.tpv29` x `python-numpy`,
~950 s of one sequential simulation — so parallelising further buys nothing.
Roughly half of each ~2 h release is spent re-deciding a question already
settled on the same bytes.

**Residual risk — stated here rather than left implicit.** Under the release
profile the other 24 CI cells' physics is NOT re-verified on the runner at the
release SHA. It is verified at `<ancestor>`. The exposure is TIME, not tree:
everything that can change between `<ancestor>`'s green run and the tag lives
outside this repo and is unpinned — GitHub's `ubuntu-22.04` image (gfortran,
mpich, libnetcdf, rolled on GitHub's own cadence) and whatever
`pip install numpy scipy netCDF4 matplotlib xarray jax` resolves to that day.
Three specific failure classes move from per-release coverage to per-merge
coverage, and a release cut days after its last code commit no longer sees any
of them at its own SHA:

  (i) a cell's peak RSS crossing the 7 GB runner ceiling — the v5.7.0 failure
      mode exactly. It is now caught at merge, monitored by
      `matrix.MEASURED_PEAK_RSS_GB`, which is measured on THIS box and is an
      ESTIMATE of the runner's footprint, not a measurement of it;
  (ii) a toolchain-version change that moves a cell's `max|diff|` past its
      `CASE_BOUND` without any commit in this repo;
  (iii) a dependency declaration gap that only the full python jobs exercise.

Bounding the time exposure needs a number this rule does not have. **PROPOSED,
owner's call, not a measurement**: the release profile is unavailable when
`<ancestor>`'s full-matrix green run is more than 7 days old, and the full
matrix runs instead. Until the owner sets that number, treat it as 0 — full
matrix — rather than as unbounded.

Adopt this or do not, but do not adopt it and then describe a release as
fully CI-verified: under the profile it is verified on the release SHA for
build, checkout, unit, regression and three smoke cells, and on `<ancestor>`
for everything else.

**How to apply**: the releaser runs the local sweep as rule 15 step 1 already
requires, commits, pushes, and READS which profile CI chose from the run's own
job list rather than assuming. The Tasks-done row (rule 15 step 4) records
three things: `<ancestor>`'s SHA, the date and run id of its full-matrix green
run, and which profile the release SHA's own run used. If the full matrix ran,
nothing was deferred and this rule did nothing that release.

**Tier: NOT mechanical today, and the workflow change must NOT land before it
is.** Condition (a) belongs in `.github/workflows/test.yml` as a job-selection
step plus a `testsys/regression/` guard proving the profile's job set is
exactly {build, unit-regression, three tpv8 smoke cells} and that any
non-whitelist path in the diff forces the full matrix; condition (b) is an
extension of `check_pretag_ci.py`, which already asks "is there a completed
successful run for this SHA" and must additionally ask "was that ancestor's
run the FULL matrix". The evidence command that settles whether a given
release SHA is eligible, and the one a reviewer runs against the claim:

    git diff --name-only <ancestor>..<sha> | \
      grep -vE '^(VERSION|README\.md|pastReleaseNotes\.md|pathway_forward\.md|PROJECT_RULES\.md)$'

Expected output: empty. Non-empty means the full matrix, and the profile was
never available. The workflow edit itself is routed to whoever conducts the
release, the guards to `iris-vermeulen`; this rule's author writes neither.

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
`EQDYNA_E2E_BIN=src/eqdyna`, and `make eqdyna` by name), while CI at the time
ran `./install-eqdyna.sh -m ubuntu` and then
`testsys/run.py unit regression e2e-ci` (CI's shape as of v5.6.0; today's
seven-job layout is below).
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
it — read `.github/workflows/test.yml` and reproduce the commands verbatim.

**CI is not one invocation, and no single command reproduces it.** As of
2026-09-21 it is SEVEN parallel jobs (`build`, `unit-regression`,
`e2e-ci-fortran-a`, `e2e-ci-fortran-b`, `e2e-ci-python-cheap`,
`e2e-ci-python-meng`, `e2e-ci-python-tpv29`), and the e2e jobs call
`run_e2e.py --ci` directly — CI never runs `run.py e2e-ci`, and never runs
`run.py unit regression e2e-ci` as one line. What it actually runs, with the
line of `test.yml` each command sits on:

    ./install-eqdyna.sh -m ubuntu                                           # :65
    export EQDYNAROOT=$(pwd); export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
    python3 testsys/run.py unit regression                                  # :149
    python3 testsys/e2e/run_e2e.py --ci --backends fortran --cases test.tpv29,test.tpv1053d,test.tpv104,test.tpv8,test.tpv37        # :198
    python3 testsys/e2e/run_e2e.py --ci --backends fortran --cases test.drv.a6,test.meng2023a,test.meng2023cb,test.tpv10,test.tpv36 # :242
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv8 --backends python-numpy,python-jax                                        # :287
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv10,test.tpv104,test.tpv1053d --backends python-jax                          # :288
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv36,test.tpv37 --backends python-numpy,python-jax                            # :289
    python3 testsys/e2e/run_e2e.py --ci --cases test.meng2023a,test.meng2023cb --backends python-numpy,python-jax                   # :327
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv29 --backends python-jax                                                    # :328
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv29 --backends python-numpy                                                  # :368

Those line numbers move; re-read the workflow rather than trusting them, and
note the python jobs also export thread-pinning env vars in the same step.
(Command list and line numbers re-verified against `.github/workflows/test.yml`
2026-09-22 — every line number above had drifted by 14-20 lines since the
2026-09-21 correction; the commands themselves were unchanged. Rule 15b, if
adopted, changes WHICH of these run on a release commit — it does not change
this list, which is what CI runs on any commit that touches code.)

The union of the `--ci` invocations is `matrix.CI_CELLS`, **25 of the 30
cells** as of 2026-09-21 — the remaining 5 do not fit a 7 GB runner. That
union is proved mechanically by
`testsys/regression/test_ci_workflow_coverage.py`, which is authoritative over
any count written in prose, this rule's included. `run.py e2e-ci` is a local
convenience that runs the same cell list; it is not what CI invokes, so a
green `run.py e2e-ci` is a close model of CI, not a reproduction of it. Run
`run.py all` too — it is the wider local gate and the one that speaks for the
whole table — but do not mistake it for CI either.

This rule was wrong about its own subject twice. On 2026-09-16 it told the
reader to reproduce CI with a command CI does not run — the exact substitution
it exists to forbid. The 2026-09-16 correction then left behind a claim that
`test.yml:64` runs `run.py unit regression e2e-ci`; :64-65 is the BUILD step,
and that single line appears nowhere in the workflow. Corrected 2026-09-21
against the file itself. A rule whose whole purpose is quoting CI verbatim has
to be re-read against CI whenever CI's job layout changes.

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

**Scope note added 2026-09-21 (item 43)**: this step is about the three
backend IMPLEMENTATIONS -- fortran, python-numpy, python-jax. It does not
require every case to opt into a new OPTIONAL execution mode of one of those
backends (e.g. `python-jax-mpi`, per-case opt-in via `matrix.PY_MPI_RANKS`).
A case fully supported per this step may still be declared UNSUPPORTED for
such a mode; expanding that opt-in is a suite-cost decision for the owner,
not a requirement this step imposes on new landings.

**How to apply**: `python3 testsys/e2e/run_e2e.py --cases test.tpvNN` for the
new case alone while iterating, then `python3 testsys/run.py all` before
committing the reference.

---

## 18. A refactor that couples two previously-independent artifacts must say so

If a change makes one artifact's correctness depend on the filesystem
location or content of another that it did not depend on before — a symlink,
a shared config file, a directory that must exist relative to another — the
same commit states the new dependency in writing and names what breaks if the
target moves, is renamed, or the platform can't represent the mechanism used
(e.g. a checkout with `core.symlinks=false`).

**This rule's symlink half is enforced**, by
`testsys/regression/test_symlink_integrity.py`: it enumerates every
git-tracked symlink from `git ls-files -s` (mode 120000, never a hardcoded
list, so a second symlink is picked up automatically and "zero tracked
symlinks" is itself a hard FAIL rather than a quiet pass), and for each
asserts the working-tree entry really is a symlink, the target resolves to an
existing in-repo path, and — where the target is a `.py` file — it actually
imports. **Its non-symlink half is not enforced by any guard.** A shared
config file two artifacts now both read, or a directory-level coupling with
no symlink in it, is caught by review only, the same as any other
refactor-introduced coupling; do not read "rule 18" on a commit as proof the
non-symlink case was checked.

**Rationale**: `45446d0` (2026-09-19) de-duplicated the `test.tpv36` /
`test.tpv37` compsets by turning `case_input/test.tpv37/tpv36_37_common.py`
into a symlink (this repo's first tracked one) to
`case_input/test.tpv36/tpv36_37_common.py`. Two risks followed and were named
in `pathway_forward.md` item 44 at merge time, not discovered later: a
`core.symlinks=false` checkout materialises the symlink as a text stub
containing the target path string, so `case.setup` dies importing it with no
useful message; and `test.tpv37` can no longer survive deletion or rename of
`case_input/test.tpv36/`, a directory-level coupling that did not exist
before. Naming both at merge time is what let the second one become a guard
instead of a future mystery failure.

**How to apply**: when a refactor introduces this kind of dependency, name it
in the commit (or the board row it lands under) the way item 44 was named —
not as an afterthought landing note, but as its own statement of what now
breaks and how. If the coupling is a tracked symlink,
`test_symlink_integrity.py` already covers it; run
`git ls-files -s | awk '$1==120000'` yourself first if you are unsure whether
you just added one. If it is a shared config or a directory coupling with no
symlink, say so in the commit and route the "should this be guarded
mechanically" question to a reviewer — this rule does not claim a check
exists for that case, because none does.

---

## 19. A shared mutable `*_last.*` artifact is not evidence until pinned to a commit

A tool-written file whose name says "last" (`testsys/perf/scaling_last.json`,
`testsys/perf/scatter_bandwidth_last.json`, `testsys/perf/tpv29_pinned_last.json`,
or any future `*_last.*`) holds the most recent run **by design** and carries
no guarantee of being the run anyone cited. Being tracked in git does not make
it durable evidence — only one specific commit of it is.

- A `pathway_forward.md` row may not cite a `*_last.*` path as evidence
  without a commit SHA pinning the content it means. A bare path is a
  citation to whatever ran most recently, not to what the row describes.
- Durable evidence goes to a dated, immutable name, never the `_last` file
  itself: `docs/perf_snapshots/<tool>_<date>_<what>.{json,log}`. That is what
  a row cites when it needs to survive the next run of the tool.
- Before committing a tool-written shared artifact, inspect the committed
  version first: `git show HEAD:<path>` and compare a property of the content
  — not the mtime, not "looks similar" — against what you are about to write.
  If another artifact or board row depends on the committed content, confirm
  you are not about to replace it.
- Two agents or sessions must not run a tool that writes a committed shared
  artifact concurrently. The write has no lock: whichever commits second
  silently discards the other's run, and nothing raises an error either time.
  A conductor dispatching parallel missions that both touch the same
  `*_last.*` path must serialize them or redirect their output to different
  paths — git will not detect the collision, let alone merge it.

**Rationale**: `testsys/perf/scaling_last.json` was last committed by
`88f5227` (2026-09-19) from a python-only run — the string `fortran` appears
in it zero times. Item 33's board row cites this exact path for its
per-config skip records. On 2026-09-21 a Fortran re-baseline sweep
(`run_scaling.py`) was about to overwrite it with a run containing `fortran`
five times, which would have silently destroyed the evidence item 33's row
points at while the row kept reading as maintained — caught only by
inspecting `git show HEAD:testsys/perf/scaling_last.json` before committing,
not by any mechanical check. At the same moment a second agent had the same
file dirty in a separate worktree; whichever of the two commits landed second
would have overwritten the other's results with no warning either way. A
`*_last` filename is not a bug in the tool — it correctly means "the last
run" — the defect is citing a deliberately-transient file as durable evidence
and running such a tool from two places at once.

**Incident (2026-09-21)**: the day's run was preserved under
`docs/perf_snapshots/scaling_2026-09-21_fortran_rebaseline_tpv104.json` (plus
its `.log`), and `scaling_last.json` was restored to its committed content
with `git checkout --` so item 33's citation kept resolving to the data it
describes (`3a46ab8`).

**Enforceable vs procedural.** Only the first bullet is mechanical: a
regression test can scan `pathway_forward.md` for a `*_last.*` citation and
fail any that has no commit SHA next to it. The other three bullets —
inspecting `HEAD:<path>` before committing, redirecting durable results to
`docs/perf_snapshots/`, and not running two writers of the same path
concurrently — rest on discipline. Nothing in this repo today detects a
skipped `git show HEAD:<path>` check or two agents racing on the same file;
reading this rule as enforcing them would be false, the same standard this
rule book applies to any other unchecked claim.

**How to apply**: before citing a `*_last.*` path in `pathway_forward.md`,
add the commit SHA that produced the content you mean, in the row's own text
(e.g. `` `testsys/perf/scaling_last.json` (committed at `88f5227`) ``); before
running a tool that writes one, `git show HEAD:<path>` and diff a real
property of the content, not the mtime, against what your run is about to
produce; before dispatching a second agent, check whether its mission writes
the same shared path and serialize or redirect if so.

---

## 20. Heavy runs are launched detached and polled by artifact, never by process

A heavy run — anything expected to exceed a few minutes: a sweep, a scaling
point, a solver run, a native or XLA compile — is launched DETACHED, with
stdout+stderr redirected to a FILE, and its progress is judged by its
ARTIFACT, never by its process. The mechanism:

    setsid nohup <cmd> > <logfile> 2>&1 < /dev/null &

The one-command tell, run BEFORE walking away:

    ls -l /proc/<pid>/fd/1

If that shows `pipe:[...]` instead of a real path, the run's output is going
to a pipe that dies with the launching turn — nothing will ever read it, and
the run's entire record dies with it. Verify fd 1 points at a path first; a
correctly detached run reads `-> /.../scaling_A.log` or similar.

Liveness is the artifact advancing — the log file's mtime moving, new rows
appearing, a NOTES checkpoint landing — not CPU%, not the process existing,
and not a foreground wait on the command.

## 20a. A tool that writes its results file only at the very end of all its work is a partial-loss hazard

Prefer several smaller invocations that each land their own artifact over one
long invocation that lands one — a killed end-writer loses everything, a
killed N-th invocation loses only the N-th. The 2026-09-22 scaling deliverable
was split into two runs for exactly this reason. (This sub-rule previously
had no `## ` heading of its own — a bolded lead sentence only — which is
exactly the shape gap that let a heading-based grep miscount the rule set
twice in one night; see 20b below and the index's count note.)

**Rationale**: an agent's turn — and with it every pipe and foreground child
it holds — can end at any time. A rule that assumes the launcher outlives the
run has now lost work FOUR times in one campaign.

**Incident (four, cumulative)**: (1) round 6's `NOTES_*` checkpoints, lost
with their worktree. (2) mira's untracked files, lost the same way.
(3) 2026-09-21: item 32's dx=250 serial run, killed at ~40 min when its
agent's turn ended — "external stop, not a solver error". (4) 2026-09-21: a
32-rank scaling point ran a FULL HOUR with `mpirun`'s stdout on
`pipe:[80141375]` that nobody was left reading — no file, no NOTES
checkpoint, no ledger row. An hour of a 32-way run, gone.

**How to apply**: before launching anything over a few minutes, wrap it in
the `setsid nohup ... > <logfile> ... &` form above, check
`ls -l /proc/<pid>/fd/1` resolves to a real path, then poll the logfile and
the run's output artifacts on a schedule. When designing or invoking a tool
for a long campaign, split it so each stage commits its own artifact (20a).
Tier: the fd-1 check is mechanical per launch; the rule as a whole is a norm
for how sessions launch work — no repo-wide gate can see a foreground pipe
after the fact, which is exactly why the check happens at launch time.

---

## 20b. A heavy run's artifacts are not landed until they are committed, and finishing does not mean anyone is still there to commit them

Rule 20 covers launching a heavy run detached so it survives the launching
turn. Rule 20a covers writing its artifact incrementally so a kill loses only
the latest stage. Neither says whose job it is to commit the artifact once
the run has actually finished — and a run finishing green is not the same
event as a person being present to `git add`/`git commit` it. Until that
commit happens, the evidence exists only in one working tree, owned by
whichever session happens to still be around, which may be nobody.

**Rationale**: a detached run's whole point is to survive its launching
turn — but that same property means it can also *finish* after its launching
turn has ended, with no session left that knows to commit what it produced.
Rule 19 already treats an uncommitted, unpinned perf-ledger row as not
evidence; this rule states the general case, for any heavy-run artifact, not
only the ledger.

**Incident (2026-09-22)**: a full `testsys/run.py all` sweep (3839.6s, 31 of
40 cells) finished green and sat with two uncommitted artifacts — an
appended `docs/perf_ledger.jsonl` and a new perf snapshot — because the
agent that launched it had already exited before it finished. Caught before
the next heavy run could have overwritten or raced it; the next version of
this gap loses the full sweep's evidence outright, the same shape of loss as
20a's incident (4), one stage later in the lifecycle.

**How to apply**: when launching a heavy run detached, name who checks on it
and commits its artifacts, and treat "the run finished" as the trigger to
commit immediately — not as a state that persists safely until someone next
looks. Do not treat a finished background job as self-sufficient just
because it produced its file. Tier: a norm, not mechanically gated — no
repo-wide check can see evidence sitting uncommitted in a worktree nobody is
in; the nearest mechanical backstop is rule 19's own guard (unpinned ledger
row rejected), which this rule generalises from the ledger to any artifact.

---

## 20c. A completed detached run is indistinguishable from one still in flight until something reads its artifact and updates the board — and until then its worktree is not reapable by default

Rule 20 gets a heavy run launched so it survives its launching turn. Rule 20a
gets its artifact written incrementally so a crash costs one increment. Rule
20b gets a finished run's artifacts committed rather than left dirty in a
tree nobody is watching. None of the three says who notices that the run has
FINISHED, reads what it produced, and updates the board to say so — a board
row can keep asserting "in flight" indefinitely about a process that exited
hours or days ago, and nothing in the repo will contradict it until someone
independently checks.

**Rationale**: rule 20a's signal is the artifact ADVANCING; the matching
signal at the other end — the artifact having reached a real terminal state
(an `Exit status: 0` in a time.log, a snapshot file whose last row stopped
changing) — is evidence of completion, and no rule currently assigns anyone
to read it. The gap sits one stage past both 20a (produce the artifact) and
20b (commit it): a run can finish, produce every file it was designed to
produce, and still have its actual DELIVERABLE — the comparison, the
analysis, the number the run existed to answer — never run, because the
agent that would have run it was gone before it could. A worktree holding
that unread result is then a single reap away from losing it, with nothing
that flags it as different from any other stale worktree.

**Incident (2026-09-22, one stage later than 20a's four and 20b's own)**:
item 32's dx=250 mesh-refinement experiment (worktree
`agent-a062e5e0475453665`) finished both its solver runs —
`dx250_serial/time.log` reads `Exit status: 0`, written 2026-09-22 00:04 —
while the board's P1 item-32 row kept asserting "dx=250 run in flight
2026-09-21, result pending" for roughly 21 hours, because the agent that
launched it was rate-limited before running the comparison. The run's entire
deliverable — a 30-second `flips.py` comparison of artifacts already on
disk — went unrun the whole time, and the worktree holding the only copy of
those artifacts was not flagged as anything other than an ordinary reapable
worktree. It was recovered only because an unrelated worktree-hygiene pass
happened to open it and check.

**How to apply**: when a heavy run is launched detached (rule 20), name —
beside the launch, not only in a person's head — what "finished" looks like
(an exit code in a log, a row count settling) and what the very next command
is once that state is observed; a rule-20/20a-compliant launch is not done
until that next command has actually run and the board reflects its answer.
Treat a board row's "in flight" as a claim to verify against the process
table and the artifact's own mtime before trusting it — the same discipline
rule 4 already asks of any inherited conclusion, applied to the board's own
words. A worktree holding the only copy of a run's raw artifacts behind a
result nobody has read yet is not eligible for default reaping; it must be
checked for a finished-but-unread run first, not assumed safe because
nothing flagged it.

**Tier**: a norm, not mechanically gated — no repo-wide check can watch a
board claim go stale in real time or know that a given worktree is "the only
copy" of something. The nearest mechanical backstop is the same one 20 and
20a already lean on — `ls -l /proc/<pid>/fd/1`, the artifact's own mtime and
exit code — checked by whoever next touches that worktree; a board row that
names its own evidence command (rule 14) is at least falsifiable by
re-running it, which is what caught this incident.

---

## 21. An agent's write surface is its own worktree; the main checkout belongs to the conductor

An agent session working in a linked worktree (`.claude/worktrees/<name>/`)
commits THERE and nowhere else. The main checkout — `/home/utig5/dliu/EQdyna`,
the tree whose `--git-dir` is the repository's own `.git` — is the conductor's:
no agent commits in it, merges in it, moves its HEAD, or leaves modified
tracked files in it. Wanting to "just land this one small thing" on master is
exactly the case this refuses; the conductor merges, agents branch.

Rules 20b and 20c govern a worktree's ARTIFACTS — getting them committed, and
not reaping a tree whose result nobody has read. This rule governs the tree
ITSELF: which checkout a session is allowed to write to at all. Rule 19's
fourth bullet is the nearest existing statement (two sessions must not write
one shared path concurrently) but it is scoped to a shared `*_last.*` artifact;
the checkout is the shared object this rule is about.

**Rationale**: worktrees exist so N sessions can hold N different HEADs over
one object store. A commit landed in the main checkout breaks that for every
session at once, and breaks it SILENTLY: another session's next read of a
source file simply returns different bytes, mid-build or mid-run, with no
error raised anywhere and nothing in the repo recording that it happened. The
victim's symptom — a gate that fails against source it never saw, a diff it
cannot explain — points at its own change, not at the tree moving underneath
it, so the cost is paid in debugging the wrong thing.

**Incident (2026-09-22, RELAYED by the conductor to this rule's author, not
witnessed by the author and not independently re-derived)**: two agents
committed directly into the main checkout `/home/utig5/dliu/EQdyna` instead of
their own worktrees, and a third had HEAD move underneath it mid-run as a
result. Twice in one day, same shape — which is the frequency that makes this a
rule rather than a note. The relayed provenance is stated here deliberately:
the mechanism is not in doubt, but no commit SHAs or session names were
recorded with it, so the incident cannot be re-audited from the repo (rule 4).

**How to apply**: before your first commit of a session, run

    git rev-parse --git-dir --git-common-dir

In a linked worktree the two differ — `.git/worktrees/<name>` against `.git`.
If they are EQUAL you are standing in the conductor's checkout: do not commit,
say so, and ask for a worktree. The same check settles it after a `cd` you did
not expect to change trees.

**Tier today: NOT mechanical — this rule is hortatory until the hook below
lands, and it is the weaker artifact for it.** The enforcement is a small one
and is specified here so it can be built rather than argued about: a
`pre-commit` hook that refuses when `--git-dir` equals `--git-common-dir` and
an agent-session env marker is set. Linked worktrees share
`$GIT_COMMON_DIR/hooks`, so ONE installed hook covers every present and future
worktree; it needs `git config core.hooksPath` pointed at a tracked directory
(e.g. `testsys/hooks/`) by `install-eqdyna.sh`, because git does not install
hooks on clone, plus a `testsys/regression/` guard asserting the hook is
tracked, executable, and actually refuses on the equal-dir case. Until that
exists, every clone is unguarded and this rule is enforced by reading it.
Tracked as pathway item 68; the build is routed to `iris-vermeulen`, not done
by this rule's author.

---

## 21a. A gate sweep runs in its own worktree; the shared `test/` tree has no lock, and a collision reads as a solver failure

Any run of the e2e sweep — `python3 testsys/run.py all`, or `run.py` with the
`e2e`/`e2e-ci` tier, or `testsys/e2e/run_e2e.py` directly — runs in a linked
worktree of its own, never in the main checkout, for as long as more than one
session can touch this repository. A worktree has its own `REPO_ROOT`, and
therefore its own `test/`; that separation is the only thing that makes the
rotation below safe.

Rule 21 assigns the main checkout to the conductor. That is an ownership
statement about COMMITS, and it does not make the main checkout a safe place to
RUN a sweep: the conductor's own gate collides with another session's e2e run
exactly as an agent's would. Ownership is not exclusivity.

**The mechanism**, written out so the next reader does not re-derive it.
`testsys/e2e/run_e2e.py:587-594` rotates the run tree at startup,
unconditionally and with no lock of any kind:

    test_dir = os.path.join(REPO_ROOT, 'test')
    prev_dir = os.path.join(REPO_ROOT, 'test.prev')
    if os.path.isdir(test_dir):
        if os.path.isdir(prev_dir):
            shutil.rmtree(prev_dir)
        shutil.move(test_dir, prev_dir)

The rotation itself is correct — it is rule 8's one preserved level of
evidence. The defect is that it acts on a FIXED path with no owner. A second
invocation in the same checkout renames the FIRST invocation's live working
tree out from under it, mid-run; the first invocation notices only when it next
touches an absolute path it resolved earlier. The Python backend writes
`frt.txt0` by absolute path at the very END of a cell
(`src/python/eqdyna/library_output.py:120`, reached from `eqdyna3d.py:522`), so
a cell can burn its entire runtime and then die on its last write.

**This generalises past `run_e2e.py`.** Any tool that rotates, deletes or
rebuilds a FIXED path under `REPO_ROOT` has the same shape, and every such tool
is unsafe to run twice concurrently in one checkout:
`testsys/e2e/run_e2e_full.py:129-133` (`test.full` → `test.full.prev`, the
identical pair), `testsys/perf/run_perf.py:172` (deletes
`testsys/perf_case/<name>`), `testsys/perf/run_jaxmpi_ab.py:88` (deletes
`src/python/eqdyna/__pycache__`). A tool that works inside
`tempfile.mkdtemp()` — as every `testsys/regression/` test does — is not in
this class and needs no worktree. This is the filesystem counterpart of rule
19's fourth bullet, which says the same thing about a shared committed
`*_last.*` artifact.

**The failure mode is a FALSE RED, and that is worse than it sounds.** Nothing
is silently corrupted: a cell whose directory has vanished fails loudly and
immediately. The damage is that the loud failure is INDISTINGUISHABLE AT A
GLANCE from a real solver failure — it surfaces in the sweep summary as
`FAIL <case> x <backend>`, the same string a genuine parity breach produces. A
conductor who reads that line and reverts a good change, re-gates a landing
that was fine, or holds a release has been misled by infrastructure. The
evidence that tells the two apart — a `FileNotFoundError` on a path under
`test/`, rather than a `max|diff|` over the case bound — is buried in the
cell's own log, not in the summary anyone actually reads.

**Incident (2026-09-22; the conductor who wrote this rule is its violator)**:
wei-lin launched the v5.16.0 gate sweep (`python3 testsys/run.py all`) in the
MAIN CHECKOUT at 21:48, knowing the isolation convention and not following it.
At 22:12:20 a different Claude session started its own e2e in the same checkout
and rotated the in-flight sweep's `test/` to `test.prev/` — that session's own
artifact dates it exactly
(`docs/perf_snapshots/e2e_cells_2026-09-22_221220_2462614.json`, landed on
master in `88a4012`). `test.tpv1053d x python-numpy` then ran **1500.1 s** and
died with `FileNotFoundError: [Errno 2] No such file or directory:
'/home/utig5/dliu/EQdyna/test/test.tpv1053d.python-numpy/frt.txt0'`. Four
further python-numpy cells (`test.tpv29`, `test.tpv36`, `test.tpv37`,
`test.drv.a6`) were already running into paths that no longer resolved, doomed
for the same reason, and were killed by PID. 25 cells had reported SUCCESS
before the collision. Cost: roughly 40 minutes of a release session, plus a
restarted release gate. The rule exists because it was done, not because it was
foreseen.

**How to apply**: `git worktree add` a tree for the sweep and launch it from
there — `git rev-parse --git-dir --git-common-dir` must differ (rule 21's own
check settles this too), and the sweep's own `REPO_ROOT` must be that worktree,
not `/home/utig5/dliu/EQdyna`. Before acting on ANY `FAIL` line, open that
cell's log and confirm the failure is numeric: a `FileNotFoundError`, or a
missing directory under `test/`, is a collision and not a result — nothing may
be reverted, re-gated or held on it. If you must diagnose a collision after the
fact, `test.prev/` holds the rotated tree and the other session's timestamped
snapshot in `docs/perf_snapshots/` names who rotated it and when.

**Tier: NOT mechanical, and this rule says so about itself.** A written rule
that depends on every session remembering is strictly weaker than a lock, and
both sessions in the incident above were able to read this book. The mechanical
fix — a lockfile on `$REPO_ROOT/test`, or a PID/SHA-stamped run directory so
two invocations cannot name the same tree — is NOT DONE; it is tracked as
pathway item 70 and routed to `iris-vermeulen`. Until it lands, the collision
is forbidden but not prevented, and this rule is enforced by being read.
