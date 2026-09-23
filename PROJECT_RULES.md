# EQdyna Project Rules

Index — read this list first; jump to a rule only when it's load-bearing.

1. Minimal changes; no unnecessary new files.
2. No silent fallbacks, swallowed errors, or placeholder data.
2a. A scripted edit to a tracked document asserts on its shape, not a substring of it.
3. Gate every stage; pass before moving on.
3a. A tier red on master is a stop-everything condition, regardless of cause.
3b. A red is found by polling the run list, not by filtering to the SHA you already expect.
3c. A change is gated on a case that EXECUTES the path it changed.
4. Only fresh runs are evidence.
4a. A proposed CAUSE is falsified by the outcome curve, not by the defect it predicts.
4b. An exclusion names the metric that produced it; a metric blind to a mechanism cannot exclude it.
4c. A gate downgrade names the geometry class its evidence came from, and clears no class that evidence cannot exercise.
4d. A BLOCKING precondition decays; re-measure it before citing it.
4e. A multi-rank performance ratio names its placement policy, and a mechanism claim needs a packed-vs-spread control before it is cited.
5. One calibrated definition of "pass" — never invent a metric.
5a. A provably injective relabelling is gated at bit-identity, not at the case bound.
6. Every performance number carries its provenance.
6a. `EFFECTIVE_CORES` admits a measurement; it does not make two measurements comparable.
7. Reference data is read-only.
8. Never delete evidence unless the result is a confirmed pass.
9. Cheap targeted check before expensive run.
10. Every bug gets a regression test before the fix ships.
10a. A regression guard asserts on BEHAVIOUR, not on the source text of the fix.
11. Docs move with the code, in the same change.
12. Untracked build artifacts never accumulate in the working tree.
13. File permission changes are reviewed individually, never bulk-applied.
14. A living status board, re-checked on a schedule.
14a. A board row's evidence command must be capable of both outcomes.
15. Releases follow the documented workflow, notes lead the README.
15a. A pre-tag CI check on the exact SHA is required before `git tag`, and it must be mechanical, not remembered.
15b. The local sweep and CI gate different failure classes; CI's release run re-verifies only what a local sweep structurally cannot.
15c. The tag push and `gh release create` are one action; the release-completeness guard is the check that a releaser split them.
15d. Step 5's single release commit does not cover the two files 21c owns; those split out.
15e. A release commit is gated on itself before it is pushed, not carried on step 1's earlier green.
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
21b. No session writes the main checkout — conductors branch too, and its HEAD moves only by fast-forward sync.
21c. `PROJECT_RULES.md` and `pathway_forward.md` have exactly one writer per session.
21d. A dispatch carries its isolation and its scope in writing, or it is not issued.
22. A scope restriction is itself a rule, and it can conflict with another rule.
23. Fortran is the reference implementation; the port follows its NUMERICS, not its file layout.
24. A release tag requires a committed full-term local sweep at the exact SHA, not only green CI.

Count, stated so a heading-shape grep does not undercount it again (that
undercount happened twice in one night, 2026-09-21/22): 24 numbered rules
(1-24) plus twenty-five lettered sub-rules (2a, 3a, 3b, 3c, 4a, 4b, 4c, 4d,
4e, 5a, 6a, 10a, 14a, 15a, 15b, 15c, 15d, 15e, 20a, 20b, 20c, 21a, 21b, 21c,
21d) — 49 `## ` headings total. Verify: `grep -c '^## ' PROJECT_RULES.md`
reads 49; `grep -c '^## [0-9]*\. ' PROJECT_RULES.md` (numbered rules only, no
letter suffix) reads 24.
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

## 3a. A tier red on master is a stop-everything condition, regardless of cause

Rule 3 names the commands and the pass criterion. This sub-rule states what a
FAILURE of that criterion means once it shows up at master HEAD: work stops
until the red clears — whatever the cause. "It's just the known
release-ordering thing, ignore it" is not a standing exemption; nothing in
this project marks a class of red as pre-approved to skip past.

**Rationale**: the difference between "red because of a code bug" and "red
for some other, already-understood reason" exists only in the head of
whoever is currently context-loaded on the incident. The next reader — a
different agent, a different day, reading the tier's exit code and nothing
else — sees one signal: red. Treating one class of red as safe to wait out
is exactly how a real break later gets read as "the known one" and gets
waved through.

**Incident (2026-09-23)**: see rule 15c. `test_release_complete.py` went red
at master HEAD for a benign, already-diagnosed reason for the second time in
this project's history (v5.8.2, 2026-09-16, and v5.16.1). Both times the
cause turned out to be nothing — a release action mid-split, not a defect —
but the rule this incident argues for is for the time it will not be.

**How to apply**: on any red tier at master HEAD, stop and diagnose before
starting or continuing anything unrelated — do not let "I already know why
that one's red" become a reason to proceed with other work while it stays
red.

**Tier: norm, not a gate.** Nothing scripts "an agent read a red tier and
kept going anyway" — the guard only reports the red itself
(`testsys/run.py unit regression`'s exit code). Making this mechanical would
need a recorded acknowledgment step, e.g. a required `--ack-red <reason>`
argument logged before any subsequent command that assumes green runs
anyway; nothing like that exists today.

---

## 3b. A red is found by polling the run list, not by filtering to the SHA you already expect

Rule 3a says a red tier at master HEAD stops everything, regardless of
cause. This sub-rule states how that red gets SEEN in the first place: by
reading `gh run list` for the branch as a whole and scanning every row back
to the last known-green one, not by querying only the one SHA already being
gated. A `--commit <SHA>` (or equivalent) query answers "is my commit green";
it silently excludes every other red sitting on the same branch history,
which is exactly what a stop-everything rule needs to see to fire at all.

**Rationale**: a rule that says "stop on any red" is only as good as the
habit of looking for one. Filtering a run list down to the commit you came to
check is the natural thing to do when gating your own change, and it is
precisely the query shape that cannot see a red anywhere else on the branch.

**Incident (2026-09-23)**: CI run 35833383579 on `73ac3a5` (item 81's landing)
FAILED on `unit-regression`, and it went unread for two hours because `gh run
list` was being filtered down to the SHA being gated and a red on the same
branch slid past outside that filter. Rule 3a says that stops everything; it
stopped nothing, because nobody looked. The red itself turned out to be a
flake in `test_stop_exit_status.py` — a two-rank MPI stdout-capture race, not
a solver defect (see pathway item 89) — but that does not soften the finding:
the rule is about a red being SEEN, not about what it turns out to mean once
it is.

**How to apply**: before relying on rule 3a's "nothing is red" precondition,
run `gh run list` unfiltered or filtered only to the branch, and read every
row back to the last known-green one, before separately checking the run you
came to check. Do this every time, not only when something already feels
wrong.

**Tier: norm, not a gate**, same as 3a — nothing scripts "the operator polled
the list before filtering."

---

## 3c. A change is gated on a case that EXECUTES the path it changed

The run that clears a change must be one in which the changed code actually
executes, and the report must name the OBSERVED artifact that moved — the line
in the output file, the header text, the count that differs. A green run on a
case where the edit is arithmetically or structurally a no-op measures the rest
of the system and says nothing about the change. If no gated case triggers the
path, say so in those words and gate on the smallest case that does, rather
than quoting the green you have.

**Rationale**: rule 3 requires a named command and a named pass criterion, and
a no-op case satisfies both while proving nothing. This is the papercuts file's
recurring shape — *a green result that tested nothing* — in the one place it is
hardest to see, because the run is genuinely green, the parity bound genuinely
holds, and the only thing missing is that the changed lines never ran.

**Incident (2026-09-23)**: the item-85 mission proved a `dsin(dip)` down-dip
correction in `output_onfault_st` on `test.tpv8` — a dip-90 case, where
`dsin(dip)` is 1 and the correction is exactly a no-op. The gate was green and
carried no information about the fix. The item-88 mission, the same defect one
subroutine down in the same file, had to be told explicitly to prove itself on
a case that triggers the path; it then gated on `test.tpv8`'s **11 off-fault
body-station files, two of them at 12 km depth**, and quoted the header line
the change produces (`# location = -3.0 km off fault, 12.0 km along strike,
12.0 km depth`) plus a 7-field data line against a 7-name legend. That is the
difference between a gate and a formality, and it cost a dispatch instruction
to get.

**How to apply**: in the same sentence as the gate result, name the case, the
reason that case reaches the changed branch, and the artifact text or number
that differs from before. "Parity unchanged" is necessary and never sufficient
— a metadata-only change must show the metadata.

**Tier**: not mechanical in general — no script knows which branch a change was
meant to reach. It is checkable per change by a reviewer, and the nearest
mechanical backstop is rule 10/10a: a guard that goes RED when the fix is
reverted proves a triggering case exists, because the guard is running on one.

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

## 4d. A BLOCKING precondition decays; re-measure it before citing it

A precondition that BLOCKS work — "this needs an idle box", "that API is
unavailable", "the runner cannot hold this case" — is a measurement with a
timestamp, not a property of the world. Before citing one to defer a task a
second time, re-run the measurement that established it and write the new date
beside the answer. A blocker inherited across sessions is a hypothesis under
rule 4 like any other, and it is the most expensive kind, because it ends the
investigation before it starts and leaves nothing to reproduce.

**Rationale**: rule 4 already says an inherited conclusion is a hypothesis and
warns hardest about "impossible / inherent / already fixed". This sub-rule
names the case where that warning is systematically not applied: nobody
re-checks a blocker, because a blocker is not a result — it reads as a fact
about the environment, and the row that carries it stops looking like work.
The asymmetry is what makes it worth a rule: re-measuring costs one command,
and not re-measuring costs the whole waiting period.

**Incident (2026-09-23)**: item 33 sat blocked for two days on "REQUIRES AN
IDLE BOX", declared UNSATISFIABLE on 2026-09-21 when **0 of 64 cpus** were
under the strict `--busy-ceiling 0.2`, measured twice, and re-scoped to an
owner decision (raise the ceiling, or obtain a reserved allocation). The
precondition was measured again at 06:28 on 2026-09-23 before anything else
ran: **53 of 64 cpus under 10% busy**, exactly 7 pegged by foreign single-core
jobs. Three full sweeps then completed at the strict ceiling, never overridden,
and the measurement the row had waited for since 2026-09-16 landed at
`b978637`. No owner decision was needed and none of the re-scoping options was
taken. The cost of the check was one command; the cost of not making it was two
days, and it also left three superseded numbers in circulation (see item 33).

**How to apply**: a board row whose `blocked on` cell names an environmental
condition carries the date that condition was last MEASURED and the command
that measures it — same contract as rule 14a's evidence command, applied to
the blocker rather than to the claim. When you pick up such a row, re-run that
command first and record the result whichever way it comes out. If the probe
has no committed entry point, say so in the row and treat exposing it as part
of the work.

**Tier**: partly mechanical. The DATE on a blocking precondition is checkable —
a `blocked on` cell citing a measurement older than the row's re-check interval
is a finding a script can raise. Whether the condition still holds is not
checkable without running the probe, which is exactly why the probe must exist
as a command rather than as a paragraph.

---

## 4e. A multi-rank performance ratio names its placement policy, and a mechanism claim needs a packed-vs-spread control before it is cited

A jax-vs-Fortran (or any backend-vs-backend) ratio measured across multiple
MPI ranks is confounded by WHERE those ranks land on the box's NUMA nodes,
not only by rank count. A scaling table that reports ratios at increasing
rank counts without saying which placement policy produced each row is
silently reporting a placement artifact as a rank-count trend. Before citing
a mechanism — "the port's halo cost dominates above N ranks", "decomposition
imbalance explains the plateau" — for a ratio measured across ranks, run a
packed and a spread placement at the SAME rank count and report both; a
mechanism claim resting on one placement is a hypothesis under rule 4a, not a
finding.

**Rationale**: this is rule 4a's discipline (a proposed cause is falsified by
the outcome curve, not by the defect it predicts) applied to a confound that
is specific to this box and easy to miss, because it looks like an
implementation detail of the measurement tool rather than a variable of the
experiment.

**Incident (2026-09-23)**: `docs/SESSION_LOG_2026-09-23_autopilot.md` section
EE read a 1D-slab jax-MPI-vs-Fortran table (ranks 1/2/4/8/16, ratios
0.83x/1.07x/1.07x/1.56x/1.58x) and attributed the turn at 8 ranks to the
slab's halo surface no longer shrinking with rank count — a real per-rank
number, visible in the run's own exchange-time breakdown. The table was taken
entirely under `run_mpi_scaling.py`'s default `least_loaded_cpus`, which
sorts candidates `(busy, node, cpu)` (`testsys/perf/run_mpi_scaling.py:103`)
and therefore packs onto the lowest-numbered NUMA nodes on a quiet box. A
same-SHA (`50c1277`) spread placement (2 ranks per node across all 8 nodes)
at 16 ranks measured jax-MPI at **0.83x** of Fortran; two packed points at 16
ranks (3 nodes) measured **1.68x** and **1.71x** — three runs, same code,
minutes to hours apart. The halo-surface mechanism was cited, and a
3D-decomposition design (`56d2401`) registered a prediction built on it,
before the packed-vs-spread control existed.

**How to apply**: when a scaling row or table crosses more than one rank
count, state the placement policy (packed/spread, or the CPU list) each row
was taken under. Before writing a mechanism sentence for a rank-count trend,
run the same rank count under both a packed and a spread placement; if the
ratio moves by more than the mechanism's predicted effect, the ratio is about
placement, not the mechanism, and the row says so.

**Tier**: not mechanical for the control itself — nothing can check that an
agent actually ran the spread arm before writing a causal sentence. A
retroactive backstop is checkable: `sorted(Counter(cpu // 8 for cpu in
r['cpus']).values())` over a ledger row's `cpus` field names its placement
after the fact, so a table whose every row groups the same way (never
spread) is a script-detectable gap even though the missing control itself is
not.

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

## 6a. `EFFECTIVE_CORES` admits a measurement; it does not make two measurements comparable

`EFFECTIVE_CORES >= 0.99`, computed from a process's own `getrusage`, is an
ADMISSION filter: it says the process was handed the cpus it asked for. It is
not a comparability proof, and a perf claim may not rest on it as one. On a
shared box the only defensible claim is a PAIRED ratio — both arms in the same
repetition, on the same cpu set, arm order shuffled between repetitions — and
an ABSOLUTE ms/step or wall-clock number is quotable only with the box's state
recorded beside it (what else was running, at what rank count, and the
`--exclude-cpus` set in force).

**Rationale**: rusage bills one cpu-second per wall second to a core spinning
on cache misses exactly as it does to one making progress. It cannot see
memory bandwidth, an SMT sibling, or a frequency drop. Rule 4b already names
this instrument blind when it is used to EXCLUDE a mechanism; this sub-rule
says the identical blindness applies when it is used to admit a number that is
then quoted absolutely. Rule 6's provenance fields (machine, cores, compiler,
SHA) do not capture contention either, which is why `pathway_forward.md`'s
header carries a tenancy warning about its own ledger rows.

**Incident (2026-09-23)**: the identical `eqdyna3d.build_solver_state` on
`test.tpv104` measured **19.03 s at eff 1.00** pinned to cpu 0 while a 16-rank
MPI job was measuring elsewhere on the box, and **4.42 s at eff 1.00** on a
quiet box — 4.3x apart at the same "effective cores", both admitted by the
filter. The same night, base-arm per-step in `testsys/perf/run_jaxmpi_ab.py`
ranged **84.30-110.74 ms** at 16 ranks on ONE cpu set and **32.43-57.40 ms** at
32 ranks, every point admitted; only the within-repetition paired ratios
carried meaning (16 ranks: mean 1.034, sd 0.042 over 5 pairs; 32 ranks: mean
1.011, sd 0.013 over 3 pairs). Most of the spread was this team's own
concurrent sessions.

**How to apply**: quote a ratio paired inside its own repetition and cpu set,
never a ratio assembled from two repetitions. Keep the `>= 0.99` filter — it is
doing real work; it correctly rejected two 32-rank base points at 0.83 and 0.86
the same night — but report it as "the process got its cpus", never as "these
two numbers are comparable". When an absolute number must be quoted, write the
box state in the same sentence, and say plainly that absolute numbers from
different repetitions on this box are not comparable.

**Tier**: not mechanical. No script can know what else a shared box was doing
at measurement time; the nearest backstop is the perf ledger's own provenance
fields (rule 6, rule 19's pinning) plus this rule being cited by the row that
quotes the number.

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

**Two references per case, when the case has two terms (added 2026-09-23,
owner decision relayed by the conductor — codified here, not yet built: see
rule 24).** A case whose full term differs from the 5 s everyday term (as of
2026-09-23: `test.tpv29` 20 s, `test.tpv36`/`test.tpv37` 6 s; `test.tpv30` 20 s
if it is ever gated) carries TWO committed references under
`test.reference.results/<case>/`, not one — the existing full-term
`frt.canonical.txt`, compared against by the release sweep (rule 24), and a
second, 5 s reference, compared against by the everyday local gate (rule 9).
Each is its own reviewed commit under this rule, generated from the SAME code
at its own term, never to make a cell pass. The gate-term reference is a
separate artifact with its own commit; it does not replace the full-term
reference, and the full-term reference is never regenerated at the shorter
term to save time.

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

## 10a. A regression guard asserts on BEHAVIOUR, not on the source text of the fix

The guard rule 10 requires asserts what the code DOES — the refusal, the exit
code, the written value, the second caller that must lose the race. It does not
assert that a particular literal appears inside a particular function body. A
guard written against source text pins one implementation of the fix, and the
defect it protects against is then only reachable by the implementations that
guard permits: any correct alternative — most often the DEDUPLICATION that
would end the defect class — turns it RED and reads as a regression.

**Rationale**: rule 10 asks for a test "asserted tightly enough to catch a
return". Source-text matching is tight in the wrong dimension: it is
simultaneously too strict (it forbids refactors that preserve every behaviour)
and too weak (a file can contain the blessed substring and still take the lock
after the delete, if the substring moved). Rule 2a already draws this line for
scripted edits to tracked DOCUMENTS — assert on shape, not on a substring —
and the same argument applies to a test reading a source file.

**Incident (2026-09-23)**: `testsys/regression/test_perf_tool_locks.py:644-695`
guards the four-times-patched "rmtree a fixed in-repo path without a lock"
class (pathway items 70, 74, 77, 81) by reading four named function bodies and
requiring the literal strings `shutil.rmtree(d)` and `run_e2e.make_serial_case(`
inside them, plus the literal expression `os.path.relpath(os.path.dirname(d)`
in two of them. The single change that would prevent a FIFTH occurrence —
collapsing the four near-identical builders into one shared builder — turns the
guard RED with `contains no 'shutil.rmtree(d)'`. Proved by mutation, not
argued. The guard against the recurrence forbids the end of the recurrence.

**How to apply**: express the guard as an invocation. For a lock: run the
builder twice concurrently and assert the second REFUSES; for an ordering
requirement: remove the lock and assert the guard goes red. Where reading
source really is the only access — a Fortran branch that cannot be run from a
unit tier, say — assert the PROPERTY (this branch declares a column count equal
to the number of fields it writes, as
`test_station_header_column_count.py` does) rather than the spelling, and say
in the docstring which refactors the assertion deliberately forbids.

**Tier**: not mechanically enforced — no check can tell a source-text
assertion that pins an implementation from one that has no behavioural route.
A cheap review signal exists: a guard that opens a source file and compares
strings is the shape to look at, and each one should carry a docstring line
saying why behaviour was not reachable.

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

## 14a. A board row's evidence command must be capable of both outcomes

A row's `Command` column must be able to read as UNFIXED and, separately, as
FIXED — the row states, or makes obvious, what each output looks like. A
command that names a file, a header line, or a location the eventual fix will
not touch reads the same before the fix and after it. A "Last checked" date
next to such a command is worse than a blank one: a blank date says plainly
that nobody has looked; a date beside a command that cannot move says the row
is maintained when nothing was actually re-verified.

**Rationale (2026-09-23)**: two rows in this file carried exactly that shape.
Item 71's command, `grep -c insertFaultType
testsys/regression/test_rough_fault_normal_consistency.py`, read "0 = still
uncovered" — but the coverage landed correctly as a SIBLING file
(`testsys/regression/test_fractal_fault_geometry_derivatives.py`, `2590c65`),
deliberately, so a failure names which path broke. The named command reads 0
forever, before the fix and after it, and it cost one full dispatch before
anyone noticed the row had been closed by a file it never looked at. Item
67's command told a reader to identify the station file "by its `# location
= on fault` header line" — no such line is ever written: `stLocStamp` is
computed at `src/fortran/library_output.f90:55` and never reaches a `write`
statement (see item 85). The command could not have been run as written by
anyone, ever, and yet the row carried a recent "Last checked" date.

**How to apply**: when writing or updating a row's `Command` column, state or
demonstrate both readings before trusting either — what it prints against the
tree as it stands, and what it would print if the fix already existed
somewhere in the tree. A command whose two runs would print the same thing is
not evidence for that row; rewrite it against the artifact the fix will
actually change. Do this at the moment the row is WRITTEN, not only when it
is later re-checked — both incidents here were wrong from birth, not decayed
into wrongness.

**Tier**: not mechanical — no script can tell, from a command's text alone,
whether it distinguishes the two states it is asked to distinguish. The
nearest backstop is this rule itself, read by whoever next re-checks the row:
a reviewer who runs the command and cannot say what a different outcome would
have looked like has found a rule-14a violation, not a clean row.

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

**Scope limit, measured 2026-09-23 — this gate reads ONE workflow's
conclusion, which is not the same claim as "CI is green for this SHA."**
v5.16.0 passed this rule honestly: both "Automatic Testing of EQdyna" runs on
`894cdc1` succeeded. On that same SHA the "Publish EQdyna Docker image"
workflow was RED (run 35819232914) — `.github/workflows/publish.yml` used
`actions/checkout@v2` at the default fetch-depth 1 and `Dockerfile:33` is
`COPY . /opt/eqdyna`, so a shallow `.git` travelled into the image and
`test_history_table.py` and `test_pretag_ci_negative.py` failed inside it
reading history that was not there. The guards were RIGHT and the environment
was wrong; the fix was `fetch-depth: 0` (`a99902f`) and deliberately NOT
teaching the guards to skip when their evidence is missing, which is rule 2.
The consequence the gate did not catch: the gate step runs before the push
step, so `ghcr.io/eqdyna/eqdyna:v5.16.0` was never published and `:latest`
still points at v5.15.0.

Whether this rule should require EVERY workflow for the SHA rather than the
testing workflow is a methodology change at the release boundary and is NOT
decided here — it is the owner's, recorded as pathway item 78 with the
recommendation that it should. Until it is decided, this rule requires what it
has always required, and a releaser who wants the wider claim reads
`gh run list --commit <SHA>` themselves and says so in the release record.

---

## 15b. The local sweep and CI gate different failure classes; CI never re-verifies physics, at release or otherwise

Rule 15 step 1 and rule 15a together used to make every release pay for two
nearly identical physics sweeps, run strictly in series, and this sub-rule
used to divide them by ancestor and by a staleness window. **Rewritten
2026-09-23** (owner decision, relayed by the conductor): the division is no
longer "when may CI's release run reduce to a smoke profile" — CI does not run
a physics sweep on any commit, release or otherwise, so there is nothing left
to reduce. It does not touch 15a's core: nothing red is ever pushed, the tag
goes last, and it goes only on a COMPLETED, successful CI conclusion for the
exact SHA being tagged.

**The division, as it now stands.**

- **CI is the MERGE gate, and its e2e coverage is a PORTABILITY check, not
  physics coverage.** It owns every failure class this development box cannot
  represent: a clean checkout of the COMMIT (rule 16's v5.5.0 partial
  `git add`), the real entry point `./install-eqdyna.sh -m ubuntu` on a bare
  machine (v5.6.1), the declared apt/pip dependency set resolving from nothing
  (v5.6.2, scipy), the 7 GB runner memory ceiling (v5.7.0), ubuntu-22.04's own
  gfortran/mpich/libnetcdf/numpy/jax versions, and the network-side release
  guards that need `GH_TOKEN`. Its only e2e cells are `matrix.CI_CELLS`,
  redefined to exactly one smoke set — `test.tpv8` x {fortran on the runner's
  mpich, python-numpy, python-jax} at the 5 s everyday term — proving the
  toolchain still builds and executes one case end to end from a clean
  checkout on a different MPI. `testsys/regression/test_ci_workflow_coverage.py`
  is retained and stays authoritative over what `CI_CELLS` actually is; this
  rule does not restate its count.
- **The local sweep is the ONLY physics gate this project has**, at either of
  two terms: the 5 s everyday term for ordinary gating (rule 9), and the full
  term at release, committed as evidence before the tag (rule 24). Nothing in
  CI substitutes for either.

**What used to be an owner-set staleness bound is no longer a bound to set.**
This rule previously proposed that a release commit could inherit an
ancestor's full-matrix CI green so long as that ancestor's run was no more
than some number of days old, "until the owner sets that number, treat it as
0." The owner has not set a number, and under the division above there is no
ancestor full-matrix CI run to go stale in the first place — CI does not run
one, on any commit, release or merge. The bound stays owner-held and stays
treated as 0: a release requires the full-term local sweep at the exact
release SHA (rule 24), with no ancestor-based or time-based exemption. This
rule does not choose a number in the owner's place.

**Rationale**: measured 2026-09-23, CI's e2e sweep duplicated roughly 25 of
the table's 30 cells on every push at a median 55.4 minutes (n=27) — cells the
local sweep had already run, on the same tree, for the same change — and it
was never a full-matrix run regardless of how fresh the ancestor was: the 7 GB
runner ceiling has excluded 5-9 cells from it throughout (rule 16). The
previous version of this rule managed that duplication with a staleness
window; the redesign removes the duplication instead.

**How to apply**: read CI's job list as portability evidence only — a green
smoke run says the toolchain builds and runs `test.tpv8` end to end from a
clean checkout on this runner's MPI; it says nothing about any other case,
term, or backend combination. Physics evidence comes only from the local
sweep, at whichever term rule 9 (everyday) or rule 24 (release) requires for
the work at hand. Do not describe a release as "CI-verified" for physics —
under this division it never is, by design.

**Tier: NOT mechanical today, and the workflow change must NOT be assumed
landed before it is verified.** As this rule is written, the `.github/workflows/test.yml`
edit removing the e2e sweep and redefining `matrix.CI_CELLS` to the tpv8 smoke
set is described by the owner as being written, not confirmed built — do not
cite this rule as evidence the edit exists; check `test.yml` and
`testsys/matrix.py` directly. Once it lands, `test_ci_workflow_coverage.py`
already proves CI runs exactly `CI_CELLS`, whatever it is redefined to; rule
24's own guard is the mechanical backstop for the local side.

---

## 15c. The tag push and `gh release create` are one action; the release-completeness guard is the check that a releaser split them

Step 7 already runs the tag push and the GitHub Release creation as one
uninterrupted command sequence with no pause in between. This sub-rule names
the check that proves a releaser actually did that, and the reason a failure
to is worse than "the page appeared a little late."

After `git push origin vX.Y.Z`, the very next command is `gh release create
... --verify-tag` — no status check, no editor, no wait for CI, nothing else
runs in between. This is checkable with the guard already in tree, not
aspirational: run `python3 testsys/run.py unit regression` (or
`test_release_complete.py` on its own) right after the tag push. If its sole
failure is `check_network_side`'s `no GitHub Release for v<X.Y.Z>`
(`testsys/regression/test_release_complete.py:237-266`), that is not a new
finding to route elsewhere — it IS the guard reporting, in real time, that
step 7 has already been split. The fix is to finish `gh release create`
immediately, not to open a pathway row and move on to something else.

**Consequence, stated because it is what makes this worth obeying**: for
however long that gap lasts, the regression tier — the gate rule 3 (and 3a)
say is a stop-everything signal — is RED on master for a reason that has
nothing to do with any code under test. A red that turns out to be "nothing,
wait for the release script" trains whoever reads it to discount the NEXT
red, which might not be nothing. That is the failure the tier exists to
prevent, and it is a worse cost than the missing Release page itself.

**Incident (2026-09-23)**: `v5.16.1` reproduced the exact pattern step 7 was
merged 2026-09-16 (from v5.8.2) to prevent. Verified fresh, not relayed:
`git for-each-ref refs/tags/v5.16.1 --format='%(taggerdate:iso-strict)'` reads
`2026-09-23T01:54:22-05:00` (`2026-09-23T06:54:22Z`); `gh release view
v5.16.1 --json publishedAt` reads `2026-09-23T07:16:39Z` — a 22-minute gap
between the tag existing and the Release being published, during which
`check_network_side` was red at master HEAD for this reason and no other.

**How to apply**: never let anything — including watching CI, including
writing the release notes body, including this very board — sit between the
tag push and `gh release create` in step 7's sequence. If something must
intervene, treat the resulting red exactly per rule 3a: stop, finish the
Release immediately, then resume.

---

## 15d. Step 5's single release commit does not cover the two files 21c owns; those split out

Step 5 says "commit everything above together" — `VERSION`, the README notes,
and (step 4) the `pathway_forward.md` Tasks-done row, as one commit. That is
no longer obeyable as written: `pathway_forward.md` is one of the two files
21c's `pre-commit` hook refuses to see staged alongside anything else. A
release therefore lands in at least two commits, not one: the Tasks-done row
by itself, landed first; then `VERSION` + the README notes + release notes,
together, on top of it. The release commit's own message can cite the board
commit's SHA rather than asserting a row that, at the moment step 5 runs,
does not yet exist on the branch.

**Rationale**: rule 15 predates 21c's mechanical separation. Followed
literally today, step 5 is refused by a check this same book made mechanical
(item 80/21c, 2026-09-23) — a release rule a releaser cannot obey without
tripping a different mechanical rule gets "fixed" by disabling the hook
instead of splitting the commit, which is exactly what 21c exists to prevent.

**Incident (2026-09-23)**: v5.16.2's Tasks-done row was written by a separate
dispatch for this reason alone — the releaser named the extra dispatch and
the resulting two-commit release as a cost paid deliberately rather than
routed around the hook, and asked this book to say so rather than leave it to
be rediscovered at the next release.

**How to apply**: at release time, expect two commits: (1) the board-only
commit adding the Tasks-done row (rule 14/21c), pushed first; (2) `VERSION` +
README + release notes, committed and pushed together on top of it, only
after (1) is on the branch. `check_pathway_tasks_done_row`
(`testsys/regression/test_release_complete.py:164-177`) only greps the
working-tree text for a dated row naming the version — it needs the row
reachable from the tag's history, not in the same commit as `VERSION`.

**Tier**: mechanical as a consequence of 21c's own hook
(`testsys/hooks/pre-commit`, `test_precommit_board_separation_guard.py`); this
sub-rule adds no new check of its own, it names what the existing one already
implies for the release sequence.

---

## 15e. A release commit is gated on itself before it is pushed, not carried on step 1's earlier green

Step 1 gates the tree BEFORE the `VERSION` bump (step 2). Steps 2-5 then edit
`VERSION`, the runtime banner rule 11 requires move with it, the release
notes and the board. Nothing between step 1 and step 6 (push) re-runs the
tier on the tree those edits produced — so step 1's green is evidence about a
tree that no longer exists by the time the release commit is pushed, and the
first thing to actually gate the committed diff is CI, after the push.

Before step 6, run `python3 testsys/run.py unit regression` again, fresh, on
the exact release commit — a second run of the same tier rule 16 already
requires be run on the committed tree, made a second run here because the
tree changed under step 1's green after step 1 finished.

**Rationale**: a gate that only ever runs on the pre-edit tree cannot see a
defect step 2-4's own edits introduce, however cheap that defect would have
been to catch locally.

**Incident (2026-09-23)**: v5.16.2's release commit `7144a98` bumped
`VERSION` 5.16.1 -> 5.16.2 without moving the runtime banner
(`src/fortran/eqdyna3d.f90:17`) — a rule-11 violation step 1's earlier green
could not have caught, because `VERSION` had not moved yet when step 1 ran.
`test_version_banner.py`, already in the regression tier, went red at master
HEAD — rule 3a's exact condition — and CI run `35838379925` on `7144a98`
failed at `unit-regression`, correctly holding the tag. Caught independently,
before gating the tag, by running the tier at master HEAD; fixed in `2964325`
(banner -> 5.16.2, tier green, `test.tpv8 x fortran` oracle re-run since the
fix touched `src/fortran`).

**How to apply**: after steps 2-4 land and before step 5's commit is pushed,
run the regression tier fresh on that exact tree and require it green, in
addition to — not instead of — CI's own post-push check (rule 15a). A red
result here is rule 3a's stop-everything condition, applied one step
earlier than CI would have applied it.

**Tier: mechanical.** The guard this incident needed already existed
(`test_version_banner.py`); the gap was procedural — WHEN it ran, not
whether it existed — so this rule adds no new check, only the step that
calls the existing one at the right point in the sequence.

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

**CI is not one invocation, and no single command reproduces it.** Corrected
2026-09-23 against `.github/workflows/test.yml` itself (the prior text named
`e2e-ci-python-cheap` and six jobs total; that job was renamed and split
2026-09-23 once its tpv36/tpv37 numpy cells were found serialising to 42 of
its 59 min — see `CLAUDE.md`). As verified fresh at this rewrite, CI is NINE
parallel jobs (`build`, `unit-regression`, `e2e-ci-fortran-a`,
`e2e-ci-fortran-b`, `e2e-ci-python-jax-tpv8`, `e2e-ci-python-tpv36-numpy`,
`e2e-ci-python-tpv37-numpy`, `e2e-ci-python-meng`, `e2e-ci-python-tpv29`), and
the e2e jobs call `run_e2e.py --ci` directly — CI never runs `run.py e2e-ci`,
and never runs `run.py unit regression e2e-ci` as one line. What it actually
runs, with the line of `test.yml` each command sits on:

    ./install-eqdyna.sh -m ubuntu                                           # :75
    export EQDYNAROOT=$(pwd)                                                # :250
    export PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH                   # :251
    python3 testsys/run.py unit regression                                  # :253
    python3 testsys/e2e/run_e2e.py --ci --backends fortran --cases test.tpv29,test.tpv1053d,test.tpv104,test.tpv8,test.tpv37        # :302
    python3 testsys/e2e/run_e2e.py --ci --backends fortran --cases test.drv.a6,test.meng2023a,test.meng2023cb,test.tpv10,test.tpv36 # :346
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv8 --backends python-numpy,python-jax                                        # :401
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv10,test.tpv104,test.tpv1053d --backends python-jax                          # :402
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv36,test.tpv37 --backends python-jax                                         # :403
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv36 --backends python-numpy                                                  # :443
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv37 --backends python-numpy                                                  # :479
    python3 testsys/e2e/run_e2e.py --ci --cases test.meng2023a,test.meng2023cb --backends python-numpy,python-jax                   # :517
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv29 --backends python-jax                                                    # :518
    python3 testsys/e2e/run_e2e.py --ci --cases test.tpv29 --backends python-numpy                                                  # :558

Those line numbers move; re-read the workflow rather than trusting them, and
note the python jobs also export thread-pinning env vars in the same step.
This is CI's job layout as it stands at this rewrite — the redesign rule 15b
and rule 24 describe (CI stops running the e2e sweep entirely and shrinks to
a `test.tpv8` smoke set) is an owner-approved decision **not yet built as this
sentence is written**; when it lands, this whole job list and the count below
change and this rule needs a matching rewrite, not a patch.

The union of the `--ci` invocations above is `matrix.CI_CELLS`, **25 of the 30
cells**, verified fresh at this rewrite — the remaining 5 do not fit a 7 GB
runner. That union is proved mechanically by
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

**Tier: MECHANICAL FOR COMMITS since v5.16.1 (`5b7a278`, 2026-09-22), and
hortatory for everything else.** `testsys/hooks/pre-commit` refuses when
`git rev-parse --git-dir` equals `--git-common-dir` (both resolved with
`pwd -P` first), `install-eqdyna.sh:133` sets
`git config core.hooksPath testsys/hooks`, and
`testsys/regression/test_precommit_main_checkout_guard.py` guards both halves.
Verified by the conductor 2026-09-23, not inherited: the hook REFUSED with cwd
at `/home/utig5/dliu/EQdyna` ("They are EQUAL, which means this tree is the
shared main checkout"), and `git -c core.hooksPath=testsys/hooks commit
--allow-empty` inside a linked worktree succeeded, exit 0 — the load-bearing
case, since a hook that blocks legitimate worktree commits is a hook that gets
turned off.

**Two limits, stated so nobody reads the hook as more than it is.** (i) It is
INERT in any clone where `install-eqdyna.sh` has not been run, because git
installs no hooks on clone and `core.hooksPath` is per-clone config, not
tracked state. (ii) It gates COMMITS ONLY. Nothing stops a session from
editing tracked files, running a sweep (rule 21a), or leaving the main
checkout dirty — `pre-commit` never fires for any of those, and rules 21/21b
remain enforced by being read for every write that is not a commit.

Retained for the record, and superseded by what landed. The enforcement was
specified here so it could be built rather than argued about: a
`pre-commit` hook that refuses when `--git-dir` equals `--git-common-dir` and
an agent-session env marker is set. Linked worktrees share
`$GIT_COMMON_DIR/hooks`, so ONE installed hook covers every present and future
worktree; it needs `git config core.hooksPath` pointed at a tracked directory
(e.g. `testsys/hooks/`) by `install-eqdyna.sh`, because git does not install
hooks on clone, plus a `testsys/regression/` guard asserting the hook is
tracked, executable, and actually refuses on the equal-dir case. That is what
was built, minus the env marker, which 21b's closing note dropped before the
build started. Tracked as pathway item 68, CLOSED 2026-09-23.

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
`testsys/e2e/run_e2e.py:615` (`:587-594` before the lock landed; re-read the
file rather than trusting either number) rotates the run tree at startup —
since v5.16.1 behind the lock the tier note below describes, and before that
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
`testsys/e2e/run_e2e_full.py:148` (`test.full` → `test.full.prev`, the
identical pair — locked since v5.16.1, acquire at `:121`),
`testsys/perf/run_perf.py:172` (deletes `testsys/perf/perf_case/<name>` —
the path was written `testsys/perf_case/` here and on the board until
2026-09-23; `run_perf.TESTSYS` is a misnomer for `testsys/perf`),
`testsys/perf/run_jaxmpi_ab.py:88` (deletes
`src/python/eqdyna/__pycache__`). **The last two have been MECHANICAL since
v5.16.1** — the paragraph below gives the detail; this sentence read "still
UNGUARDED" for a day after that landing and contradicted its own rule, which is
what a tier claim maintained in two places does.
`testsys/runlock.acquire(REPO_ROOT, <relative path>)` is reusable for any tool
still in this class. A tool that works inside
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

**Tier: MECHANICAL for `test/` and `test.full/` since v5.16.1 (`5b7a278`,
2026-09-22); hortatory for every other fixed path named above.**
`testsys/runlock.py` takes an exclusive `flock` on the run tree;
`testsys/e2e/run_e2e.py:570` acquires it as Gate 0, BEFORE the build, and
`testsys/e2e/run_e2e_full.py:121` does the same for `test.full`;
`testsys/regression/test_e2e_run_tree_lock.py` guards it. Verified by the
conductor 2026-09-23 against the real entry point, not inherited: with the
lock held, the invocation returned **exit 1 in 0.2 s with no build attempted**,
printing the holder's pid and "NOT waiting. NOT rotating anyway. NOT falling
back to a different directory" — which is the rule-2 behaviour this needed
(refuse loudly; do not degrade into a second tree).

**The two alternatives this rule originally offered are NOT interchangeable,
and the second one is now rejected.** A lockfile was the correct choice; a
PID/SHA-stamped run directory would have let two invocations proceed under
different names and thereby defeated rule 8's single preserved level of
evidence — `test.prev/` only means something when `test/` is one fixed tree
with one owner. Do not "fix" the remaining unguarded tools by stamping their
paths.

**The two `testsys/perf/` tools above are MECHANICAL too since v5.16.1
(`8a4269a`, 2026-09-23).** The lockfile sits beside the guarded directory
rather than in one central place — `testsys/perf/.perf_case.lock` and
`src/python/.eqdyna.lock` — because `runlock.acquire` now takes a relative
PATH, not a bare name. `run_perf.py` acquires inside `build_perf_case()`, not
in `main()`, deliberately: `run_tpv29_pinned_compare.py` reassigns
`run_perf.PERF_CASE` and calls that function directly without ever reaching
`main()`. `run_jaxmpi_ab.py` acquires on `src/python/eqdyna`, the PACKAGE and
not the `__pycache__` this rule originally named, because `stage()` copies
arm-specific `driver.py`/`MPI4NodalQuant.py` into it — a collision there
yields a number under the WRONG ARM LABEL rather than a crash, which is the
worse failure. Guard: `testsys/regression/test_perf_tool_locks.py` — 10 checks
at that landing, 17 since `f4718e5` extended it over the two builders below.

**Those three holes are CLOSED too, at `f4718e5` (2026-09-23, pathway item 77)
— after the v5.16.1 tag, so on master and unreleased.**
`run_scaling.build_py_case` (`run_scaling.py:428-447` → `testsys/perf/scaling_case/<case>`)
and `run_numa_scaling.build_case` (`:223-245` → `testsys/perf/numa_case/<case>`)
acquire inside the BUILDER, not `main()`: `build_py_case` has five direct
callers (`run_mpi_scaling.py:408`, `run_jaxmpi_ab.py:169`,
`run_setup_probe.py:142`, `run_shard_scaling.py:195`), so a lock placed in any
one `main()` leaves four entry points unguarded — and that also closes the
third hole, since `run_jaxmpi_ab.py` reaches the case tree through that same
builder. Two details worth carrying: `run_scaling.py` derives `REPO_ROOT` from
its own file rather than from `ROOT`, because `ROOT` honours `$EQDYNAROOT`
while the case dir is built from `TESTSYS` — a lock rooted at `ROOT` would lock
one checkout while rebuilding another; and each tool memoises its lock, because
`flock` is per open file description and a second call in one process would
otherwise refuse against its own pid.

**The `$EQDYNAROOT` hole is closed by REFUSAL, not by a warning.** Run from a
worktree whose `EQDYNAROOT` still names the main checkout, `run_jaxmpi_ab.py`
used to stage arm-specific files into THAT tree's package (rule 21b) and, since
the lock followed `ROOT`, to lock and corrupt the same wrong tree consistently —
a lock on the wrong tree is not a smaller bug than no lock. `main()` now refuses
a mismatch as its FIRST statement, before argparse and before Gate 0, so no
lockfile lands in the foreign tree ahead of the refusal. Refusing is right here
rather than warning because every one of the 11 `EQDYNAROOT` setters in this
repo points at its own tree and nothing under `testsys/perf/` sets it at all:
no supported workflow produces a mismatch, so a mismatch is rule 2's
"fail loudly", not a case to degrade through.

**One hole in this class remains**, carried as pathway item 81:
`testsys/parity/probe_plastic_traction.py:77-82` rmtrees and rebuilds
`testsys/parity/probe_case/test.drv.a6` through the same
`run_e2e.make_serial_case`, with no lock. One caller, and a collision costs a
diagnostic probe rather than a reported measurement — smaller, not different.

What is still enforced by reading: running a sweep in the main checkout at all
(rule 21b forbids it; the lock does not know which checkout it is in — two
sessions in ONE checkout now serialise instead of colliding, which is a
different and lesser guarantee). Pathway item 70 is closed for the e2e trees,
item 74 for the two perf tools above and item 77 for the two shared case
builders; the remainder is item 81.

---

## 21b. No session writes the main checkout — conductors branch too, and its HEAD moves only by fast-forward sync

The main checkout `/home/utig5/dliu/EQdyna` — the tree whose `--git-dir` equals
`--git-common-dir` — is written by NOBODY. No session, agent or conductor,
commits there, merges there, runs a gate there (rule 21a), or leaves modified
tracked files there. Every session works in a linked worktree, and a conductor
is not an exception: the conductor merges, but it merges on a branch in a
worktree of its own and pushes, the same as everyone else.

The main checkout's HEAD moves by exactly ONE mechanism, and that mechanism
creates no commit:

    git fetch origin && git merge --ff-only origin/master     # or: git pull --ff-only

on a clean tree. If it refuses — because the tree is dirty, or because the
history has diverged and a real merge would be needed — the work belongs on a
branch in a worktree and gets pushed from there; it is never resolved in place.

**This AMENDS rule 21, it does not replace it.** Rule 21's agent half stands
verbatim. What 21b removes is the conductor's exemption, and the reason is
scope, not blame: rule 21 assigns the checkout to "the conductor", singular,
and says nothing about two. On 2026-09-22 two conductor sessions were live in
this repository at once (the second stood down about 22:30) and, under 21 as
written, both held an equally good claim on the same tree — the rule could not
adjudicate, because it never contemplated the case. Read strictly, 21 also did
NOT forbid a conductor's own commits in the conductor's own checkout: wei-lin's
direct `commit:` entries `87af56e` and `bbf62c8` were compliant with the rule
as it stood. They are reclassified by this sub-rule PROSPECTIVELY. A rule the
repo cannot adjudicate does not get applied backwards to whoever is asked about
it first.

**Why "nobody writes it" and not "one named owner".** The test has to be
decidable from the repository, not from seniority or from who asked first.
`git rev-parse --git-dir --git-common-dir` answers "is this tree the main
checkout" in every clone, today, with no state to maintain. NOTHING in a git
repository answers "am I the conductor" — a named-owner rule needs a lease
artifact (a claim file, a branch, a lockfile) that a crashed or rate-limited
session leaves stale, that a second conductor starting in parallel has no
reason to read, and that would itself want a root file rule 1 does not have a
slot for. Between an ownership question that needs an oracle and a prohibition
that needs one existing command, the prohibition is the enforceable one. The
cost is one extra `git worktree add` per conductor session; the thing bought is
that "may I write here" has the same answer for every session, always.

**Incident (2026-09-22, the same night, demonstrated in BOTH halves)**: two
conductor sessions shared `/home/utig5/dliu/EQdyna`.

- **Commits.** `88a4012` is a direct `commit:` entry in master's reflog, made
  by the second session in the main checkout — evidence of the commit half,
  recorded here because rule 21's own incident was relayed without SHAs and is
  not re-auditable (rule 4). This one is.
- **Runs.** That same session's e2e run at 22:12:20 rotated the first session's
  in-flight release gate `test/` to `test.prev/`
  (`docs/perf_snapshots/e2e_cells_2026-09-22_221220_2462614.json`), killing a
  1500.1 s cell with a `FileNotFoundError` and costing ~40 minutes plus a
  restarted release gate. That is rule 21a's incident, and 21a closes only that
  half: it forbids sweeping in a shared checkout, and says nothing about who may
  commit there.

Two sessions, one tree, two distinct loss modes in under an hour. That is the
frequency that makes this a sub-rule rather than a note.

**How to apply**: before your first write of ANY session — conductor sessions
included — run `git rev-parse --git-dir --git-common-dir`. Equal means you are
standing in the main checkout: do not commit, do not merge, do not launch a
sweep (rule 21a), `git worktree add` and move. To bring the main checkout
forward after a merge landed on origin, use the `--ff-only` form above and let
it refuse rather than converting a refusal into a merge commit.

**Note for rule 21's hook (pathway item 68) — this sub-rule settles what the
hook keys on.** Rule 21 sketched a `pre-commit` hook that refuses when
`--git-dir` equals `--git-common-dir` **and an agent-session env marker is
set**. Drop the marker. Under 21b the hook refuses on EQUAL DIRS ALONE, with no
role test at all, because no role is permitted to commit there — and that is
the only version that can be built, since the equal-dirs check cannot
distinguish a conductor from an agent and nothing in the repo sets a trustworthy
role marker. A marker-keyed hook would have exempted precisely the writer this
sub-rule forbids, and a hook that fires on the person it was meant to protect is
a hook that gets disabled. The permitted HEAD move is safe under an
unconditional refusal by construction: a true fast-forward creates no commit and
so never invokes `pre-commit`. Build routed to `iris-vermeulen` under item 68;
this rule's author writes no hook.

**Tier: MECHANICAL FOR COMMITS since v5.16.1 — item 68's hook landed
(`5b7a278`) and refuses on equal dirs alone, exactly as this sub-rule
specified.** See rule 21's tier note for the verification and for the two
limits that survive: the hook is inert in a clone where `install-eqdyna.sh`
has not been run, and it fires only on `commit`. A sweep launched in the main
checkout (rule 21a), an edit left dirty there, or a file written there is
still caught by nothing. For those, the one command this rule needs exists in
every clone (`git rev-parse --git-dir --git-common-dir`) and still has to be
run by someone who chooses to.

---

## 21c. `PROJECT_RULES.md` and `pathway_forward.md` have exactly one writer per session

These two files are written by ONE session — the one that owns the rule book
and the board — and by nobody else. A mission agent that finds a rule wrong, a
tier stale, or a board row unrunnable REPORTS it, quoting the text it would
change and naming the row it wants; it does not edit either file, however
correct its proposed wording. This is the document-level counterpart of rule
21: 21 says which TREE a session may write, this says which two FILES inside
that tree are not a mission's to write even in its own worktree.

**Rationale**: a rule book takes concurrent edits worse than code does. Two
sessions adding "the same" rule produce two phrasings of one rule at two
numbers, and a renumber breaks every commit message and board row that cites
the old one — the duplication this book's own index exists to prevent. The
board is worse still: both writers are right about different rows, and the
conductor cannot tell a fresh row from a stale one by reading the file.

**Incident (2026-09-23)**: both subagents dispatched in one night's autopilot
edited both files unasked — item 70's mission rewrote rule 21a's tier
paragraph and its own board row, item 68's mission edited the board. The
conductor stripped both files back to `894cdc1` before landing (integration
commit `02e620d`, "drop agent edits to PROJECT_RULES.md/pathway_forward.md"),
so nothing reached master; the cost was paid in the integration pass and in
two closures going unrecorded until the owning session ran. Both briefs said
the rule book was the conductor's to flip; neither FORBADE the edit, and two
independent agents did it anyway on the same night. That is the frequency that
makes this a sub-rule rather than a note. The stripped proposals remain
readable at `git diff 894cdc1 fix/item70-e2e-test-lock -- PROJECT_RULES.md
pathway_forward.md` and the same against `iris/item68-precommit-hook`.

**How to apply**: dispatch briefs say "do not edit `PROJECT_RULES.md` or
`pathway_forward.md`; report the rule or row you want changed and quote the
replacement text", not "these are mine to flip" — the second is a statement
about the conductor and was read as silence about the agent. A mission's
proposed text is welcome; it travels in the report, and the owning session
lands it.

**Tier: MECHANICAL FOR SEPARATION since v5.16.1 (`765f22e`, 2026-09-23);
hortatory for AUTHORSHIP, and the heading above still promises authorship.**
Read the two apart or this tier will be cited for more than it holds: a hook
sees staged PATHS and never sessions, so what is gated is that a rule-book or
board change lands as its OWN commit — one a conductor can drop with a single
`git revert` — while WHO wrote it remains a matter of the dispatch brief and
of this rule being read. `testsys/hooks/pre-commit` (installed via
`core.hooksPath` by `install-eqdyna.sh`, rule 21) refuses any commit whose
staged set includes `PROJECT_RULES.md` or `pathway_forward.md` together with
any other path: `pre-commit: REFUSED -- this commit MIXES the rule book /
board with other files`, exit 1. Guard:
`testsys/regression/test_precommit_board_separation_guard.py`;
`test_precommit_main_checkout_guard.py` is unmodified and green. Verified
fresh by the conductor in a linked worktree 2026-09-23 — mixed commit exit 1,
the same commit board-only exit 0 — not inherited from the building mission.

**Three limits, measured rather than assumed (2026-09-23), so nobody reads
this tier as wider than it is:**

1. Git does not run `pre-commit` for an automatic merge commit, so the
   board/code merge path is unaffected — and must STAY unaffected: a change
   that made the hook fire on merges would refuse every integration of a
   branch that touched the board.
2. `git commit --amend` is blind to it — `git diff --cached` compares against
   the commit being replaced — so an amend that folds other paths into a board
   commit passes the hook. That hole is MECHANICAL over PUSHED history since
   `76c6bfb` (2026-09-23, pathway item 80, post-v5.16.1):
   `testsys/check_board_separation.py <range>` is the first step of
   `test.yml`'s `unit-regression` job, the only job already at
   `fetch-depth: 0`. The range is `before..sha` on push and base..head on
   pull_request, falling back to `origin/master..HEAD`, and NARROWING to
   `HEAD^..HEAD` when that fallback is empty — the shape a force-push leaves,
   where "0 commits, SUCCESS" is indistinguishable from a suppressed check. It
   refuses, rather than passes, on an unhandled event, an absent
   `origin/master`, or a root-commit HEAD. Guard:
   `testsys/regression/test_ci_board_separation_step.py`. What is still
   enforced by reading is only an amend nobody has pushed yet. **Do not widen
   the range to this repository's own history without reading the number
   first**: 54 of 738 non-merge commits mix the board with other paths
   (measured 2026-09-23 at `d6d846f`), all of them from before the hook
   existed.
3. The hook is inert until `install-eqdyna.sh` has run in that clone, exactly
   as rule 21's own mechanical claim is.

Pathway item 75 is closed by this landing.

---

## 21d. A dispatch carries its isolation and its scope in writing, or it is not issued

A sub-agent is launched only from a brief that already states, in the prompt
itself, the worktree it may write and the files it may touch there. A dispatch
whose prompt is empty, a placeholder, or a bare "continue" is not
underspecified — it is UNISSUED. Write the brief and dispatch again; never send
it and steer afterwards.

**Rationale**: rules 21, 21b and 21c decide what a session may write, and every
one of them reaches an agent through its brief and through nothing else. An
agent that receives no brief inherits no surface: it starts in whatever tree the
launcher is standing in — under 21b, the tree nobody may write — with full
tools. The window between "dispatched" and "first tool use" is the whole of the
protection.

**Incident (2026-09-23)**: a conductor dispatched a general-purpose agent with
the literal prompt `placeholder`, having reached for a continue-an-agent
facility that does not exist. The cost was ZERO — 0 tool uses, main checkout
verified clean afterwards — and it is recorded BECAUSE it was free. An
ambiguous prompt in place of an inert one puts a full-tool agent in the main
checkout with no scope, and that is the same slip with a bill attached. Scoped
deliberately to the act of dispatching: this rule claims nothing about what a
briefed agent then does, which is 21, 21b and 21c's ground and is unchanged.

**Tier: HORTATORY, and nothing in this repository can change that.** A dispatch
leaves no artifact here — no commit, no file, no reflog entry — so no check in
`testsys/` can see one, let alone refuse it. What WOULD make it mechanical is a
dispatch wrapper that refuses a prompt naming no worktree, and that wrapper
lives in the agent harness, not in this repo. Until one exists this is enforced
by the person typing, and it is written down so the next occurrence is
recognised as the second and not the first.

---

## 22. A scope restriction is itself a rule, and it can conflict with another rule

A brief, dispatch, or session-level restriction — "do not edit anything under
`src/`", "read-only in this tree" — carries the same force as any other rule
in this book, and it can collide with one: a restriction written to keep a
change small can make an ORTHOGONAL rule (most often rule 11, docs move with
the code) impossible to satisfy without breaking the restriction. When that
happens, halting and reporting the conflict is correct behaviour, not a
failure to just pick a side — guessing which rule to break is worse than
asking, because the guess is invisible to whoever wrote the restriction.

**Rationale**: a scope restriction is drafted for the common case (keep a
release from turning into a refactor) without checking it against every rule
the scoped work might also owe. The conflict surfaces only when the two
collide in a real change, and by then the agent standing inside the
restriction has no way to know whether the restriction or the rule is the one
that should bend — that is the restriction's AUTHOR's call, not the agent's,
and it stays that way even when the fix is one line.

**Incident (2026-09-23)**: the v5.16.2 release brief forbade any edit under
`src/`, to keep the release from turning into a refactor. Fixing rule 15e's
incident (rule 11: the runtime banner moves with `VERSION` in the same
commit) REQUIRES an edit to `src/fortran/eqdyna3d.f90:17` — a one-line
version-string change, not a refactor, but still inside the forbidden tree.
The release engineer halted and reported rather than guessing which rule to
break, costing about forty minutes; that was the correct call, and the
defect was in the restriction — it excluded a whole tree wholesale where it
meant to exclude a refactor — not in the halt.

This generalises what rule 21d already asks for one half of: 21d requires a
dispatch to state its scope in writing; this rule says the scope, once
stated, is checked against the rest of the book rather than assumed
compatible with it, and a collision is reported rather than silently
resolved by whichever side the agent happens to guess.

**How to apply**: when writing a scope restriction, name the carve-outs the
rules already require — the version literals rule 11 names, the regression
test rule 10 requires alongside a fix — rather than excluding a whole tree
wholesale. When RECEIVING a restriction that conflicts with a rule you also
owe, stop and report the conflict rather than resolving it by guessing; the
fix belongs to whoever wrote the restriction.

**Tier: norm, not mechanically gated.** No script can read a prose scope
restriction and diff it against every rule in this book. The nearest
backstop is rule 21d's existing requirement that a dispatch state its scope
in writing, which at least gives a reviewer text to check by hand against the
rule the scoped work also owes.

---

## 23. Fortran is the reference implementation; the port follows its NUMERICS, not its file layout

Owner, 2026-09-23: *"Fortran is gold. Python and Jax follow faithfully."*
Amended by the owner the same day, and the amendment is part of the rule:
*"If there is perf gain to change structure for python and Jax, allowed."*

**1. The numerics are authoritative, and the relationship is not symmetric.**
Where `src/fortran/` and `src/python/eqdyna/` disagree on a NUMBER — a
tolerance-exceeding parity failure, a different branch taken, a different
value written — the port is presumed wrong and the burden is on whoever
claims otherwise. That burden is discharged with evidence (a spec citation, an
independent code, a measurement), never with a plausible argument about which
form looks more correct. Until it is discharged, the change goes into the
port.

**2. This rule says who defers to whom, NOT who is right against the
benchmark.** "Gold" is about authority between two implementations of the same
solver; it is not a finding about physics. A Fortran line the port translated
faithfully can still be the defect, and an open, unattributed divergence stays
open: the TPV30 python-port divergence (rule 7's eleventh reference,
`pathway_forward.md` item 19(b)) is live as this rule is written and is NOT
decided by it. When Fortran is shown wrong against a spec or an independent
code, Fortran is fixed and the port follows — the port is never quietly bent
to reproduce a Fortran bug, and this rule is never cited to close a physics
question by decree.

**3. Structure is NOT bound one-to-one.** The port may fuse, split or reshape
modules relative to Fortran where that buys MEASURED performance or is forced
by a backend constraint (vectorisation shape, jax tracing, a fused loop).
`assembleGlobalKU.py` is the standing precedent and a good one: it covers four
Fortran files (`assembleGlobalKU`, `calcElemKU`, `calcHourglassResist`,
`calcElemMass`) because they are one fused loop in the port, and that is
recorded where a reader meets the correspondence. A departure is not a gap.

**4. Three states per Fortran module, no fourth, and the forbidden thing is
the SILENT gap.** Same contract `testsys/matrix.py` already uses for
SUPPORTED / DECLARED UNSUPPORTED. Every `src/fortran/*.f90` is exactly one of:

- **(a) Counterpart** — a same-named `.py` in `src/python/eqdyna/`, meant to be
  read beside it.
- **(b) Documented departure** — a `.py` of another name or shape covers it,
  and the record names *which Fortran files it covers*, *what shape it takes
  instead*, and *what motivated it* (the measurement, per 5, or the backend
  constraint).
- **(c) Declared absence** — not ported, with a recorded reason, exactly as an
  UNSUPPORTED matrix cell records a reason.

There is no fourth state and no skip. A Fortran file with neither a same-named
`.py` nor an entry saying why is the violation this rule exists to name. The
record lives in ONE place — the correspondence table commissioned 2026-09-23
under `docs/`, and `CLAUDE.md`'s correspondence paragraph until that table
lands; two places would drift, and a reader who finds the shorter one would
read a departure as a gap.

**5. "Measured" means measured.** A departure justified by performance cites
the measurement: the artifact, the number, and the box state rules 6 and 6a
require, paired inside its own repetition. A structural change made for
EXPECTED performance that was never demonstrated is an undocumented departure
wearing a justification — record it under its real reason (readability, a
backend constraint, an inherited port shape) rather than as a perf win.

**6. A new Fortran module obliges a port decision when it lands.** The commit
that adds a `.f90` lands one of the three states for it in the same change,
the way rule 11 makes docs move with the code. "Decide later" is state four.

**Supersedes a reading of `CLAUDE.md`**: "When you change physics, change it in
BOTH or say plainly which one you changed and why" reads as two PEER
implementations. The obligation in that sentence stands — change both, or say
which you changed and why — but the relationship does not: on a numeric
disagreement the port is presumed wrong, per 1. The next reader should not have
to adjudicate that.

**Rationale**: the sweep's whole value is that the two implementations check
each other (rule 17 step 7), and a check between two peers has no tiebreaker.
Naming Fortran as the reference gives every parity failure a default direction
and stops the recurring argument about which side to change. Binding LAYOUT as
well would buy nothing and cost the port its only real freedom — the port runs
on numpy, jax and jax-MPI, and a loop shape that is right in Fortran can be
several times slower in any of them.

**Incident (the divergence, 2026-09-20/21)**: `driver.f90:30` DIVIDES by the
mass; the port used `force * inv_mass` for friclaw 4/5. A reciprocal-multiply
is the obvious, arithmetically-equivalent-looking simplification, and it was a
real divergence worth **38 rupture-arrival flips on `test.drv.a6`**. The port
moved to match Fortran. Under a peer reading, that is a debate about which form
is better; under this rule it is not a debate.

**Incident (the silent gap, 2026-09-23)**: of the 26 files in `src/fortran/`,
14 have no same-named `.py`. Three of them (`calcElemKU`, `calcElemMass`,
`calcHourglassResist`) are state (b), recorded as the `assembleGlobalKU.py`
fusion. The other **eleven** — `calcB`, `calcGlobalShapeFunc`,
`calcLocalShapeFunc`, `calcQAttenuationCoeff`, `computePMLDampingVector`,
`countMeshEntities`, `checkInputConsistency`, `errorCodes`, `library`,
`library_degeneration`, `netcdf_io` — carry no entry of any kind, so nothing in
the repo says whether each is a deliberate departure, an intentional
non-port, or a port gap nobody has noticed. Verified fresh at `928c4d4` by
listing both directories, not inherited from the audit that is classifying them.
That is 11 of 26 Fortran modules in state four, which is the state this rule
abolishes.

**Incident (the unmeasured justification, 2026-09-22)**: `PML_WEIGHT = 3.0`
(`src/python/eqdyna/MPI4NodalQuant.py:85`) was a documented GUESS carried as if
it were a cost model; measured, the ratio is ~0.6, five times lower. Recutting
the decomposition took predicted work spread 3.350x -> 1.004x and zero-`Ei`
ranks 4 -> 0 — and per-step cost at 32 ranks moved 134.12 -> 132.20 ms, 1.4%,
inside this box's run-to-run variation. A structural decision in the port,
reasoned rather than measured, bought nothing. That is why clause 5 asks for
the number and rule 4a asks for the outcome curve and not the defect.

**How to apply**: on a parity failure, change the port unless you can show
Fortran is wrong, and say in the report which of the two you changed. Before
committing a structural departure in `src/python/eqdyna/`, write its state-(b)
entry in the correspondence record in the SAME change, with the measurement or
the constraint spelled out. Before adding a `.f90`, decide its state. When you
find a Fortran module in no state at all, do not guess which it is — record it
as unclassified and route it, the same as any other open item.

**Tier: partly mechanical.** The existence half is scriptable today —
enumerate `src/fortran/*.f90` basenames, subtract the same-named
`src/python/eqdyna/*.py`, and every residual must appear in the correspondence
record; a residual that does not is a hard FAIL, and the check must also fail
when the record names a Fortran file that no longer exists. No such guard
exists as this rule lands, so the eleven above are found by reading, not by
the tier. The judgment half — whether a departure's cited measurement actually
supports its claim, and whether a numeric disagreement was resolved in the
right direction — is reviewable only.

---

## 24. A release tag requires a committed full-term local sweep at the exact SHA, not only green CI

Rule 15b now divides the gates so CI never runs a physics sweep, on any
commit, release or merge — it is a portability check only. That leaves the
local sweep as the only physics gate this project has, and this rule states
what a release tag requires of it, on top of rule 15a's existing CI check
(which stays, unchanged).

- The EVERYDAY local gate uses the 5 s term for every case (rule 9, rule 3) —
  this rule does not change ordinary, non-release gating.
- At RELEASE time, the FULL committed term is required: every runnable cell of
  `testsys/matrix.py`'s table, each at that case's own full term (`test.tpv29`
  20 s; `test.tpv36`/`test.tpv37` 6 s; `test.tpv30` 20 s if it is ever gated;
  5 s everywhere else — rule 7's two-references provision), run locally on the
  exact release tree and committed as evidence, at
  `docs/evidence/sweep-<shortsha>/summary.json`, before `git tag` runs.
- A release tag requires BOTH of the following for the exact SHA being
  tagged, or for an ancestor whose diff from that SHA lands entirely inside
  `docs/evidence/`, the perf ledger, or `pathway_forward.md`: (i) a committed
  full-term local sweep as above; and (ii) a completed, successful CI run per
  rule 15a. Neither substitutes for the other, and CI's smoke cells (rule 15b)
  are not physics evidence toward (i).

**Rationale**: CI's e2e sweep duplicated roughly 25 of the table's cells on
every release at a median 55.4 minutes per push (n=27), verifying nothing the
local sweep had not already verified on the same tree — and it was never a
full-matrix run in the first place, since the 7 GB runner ceiling has excluded
5-9 cells from it throughout (rule 16). Removing that duplication (rule 15b)
leaves a gap where CI's own sweep used to stand in as part of the release
gate; this rule is what fills it, so that "the sweep ran somewhere" is never
allowed to mean "CI ran a subset of it."

**Incident (2026-09-23)**: owner decision, relayed by the conductor: CI is
redefined to a `test.tpv8` portability smoke and the e2e sweep is removed from
`test.yml`. Without a rule naming a replacement physics gate at release, the
tag guard (rule 15a) would still pass on CI-green alone, and a release could
ship having run full-term physics nowhere at all. This rule closes that gap by
naming the local sweep, committed as evidence, as the non-optional replacement.

**How to apply**: before `git tag`, run the full-term local sweep on the exact
release SHA, write and commit `docs/evidence/sweep-<shortsha>/summary.json`,
then run rule 15a's pre-tag guard. Record both the evidence commit's SHA and
the CI run id in the Tasks-done row (rule 15 step 4). A tag with CI green and
no committed full-term evidence is not a compliant release under this rule,
and neither is full-term evidence committed on a SHA the tag does not point
at (or an ancestor differing by more than the whitelisted paths above).

**Tier: mechanical once `testsys/regression/check_pretag_ci.py` is extended to
refuse a tag lacking this evidence — NOT YET LANDED as this rule is written.**
The owner describes the guard as being written, not built; do not cite this
rule as proof the refusal exists. Until it lands, this rule is enforced the
way rule 15a's own guard was enforced before it existed: by whoever runs the
release checklist by hand, checking `docs/evidence/sweep-<shortsha>/` exists
and matches the SHA before typing `git tag`.
