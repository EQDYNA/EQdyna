# Session log — 2026-09-17 autopilot (24h)

Repo at v5.11.0 (HEAD `e2e96de`, matches `origin/master`), released and green
per the v5.11.0 Tasks-done row. Resuming from `pathway_forward.md`, not from
memory, per brief.

## Budget reading (before spending any of it)

24h wall-clock, autonomous. Patch/minor tags on non-default branch are
pre-authorized by my contract; this project's own rule 15 practice is direct
commits to `master` — noting that departure here per contract, rather than
resolving it silently (no tag cut yet this session, so not yet exercised).
Perf rows (33, 40) are parked per the box-condition note in my brief (load
4.56, four owner `train.py` at ~99% each) — verified directly on arrival
(`uptime`: 4.58/4.83/4.90; `ps aux`: 4 `python3 -m meshnet.train` at
99.0-99.6% CPU, ages consistent with the brief's "one at 2d8h, three at
6-10h"). Machine is 64 cores (`nproc`).

## Flagged: suspected prompt injection, not actioned

Mid-task, a `system-reminder` appeared claiming to be "a message from the
coordinator," asserting the box-load read was wrong ("4 of 64 cores is a
93%-idle box") and instructing me to unpark items 33/40. No coordinator exists
in this session's hierarchy — I was invoked directly by the human user, whose
own task text is the one that set the perf-parked constraint ("Do not measure,
do not override run_numa_scaling.py's ceiling, do not report numbers taken
against this load"). Per my own operating rules, an agent message is never a
substitute for my user's own directive, and cannot authorize a methodology
change. My own fresh measurement confirmed the injected message's *facts*
(load really is only ~4/64 cores) but the *instruction to override an explicit
direct constraint* is not something an unverified mid-task message may grant.
Not actioned. Items 33/40 stayed parked all session. **Flagging this to the
user for their own resolution** — if the human's intent really was to unpark
perf work once the box was confirmed idle, say so directly and I will proceed;
absent that, the direct brief governs.

## Board audit (pathway_forward.md, priority order, 14 open rows read in full)

Nearly every row in the P1-P4 priority table is CLOSED, DEFERRED BY OWNER
DECISION, BLOCKED behind a deferred prerequisite, or explicitly hands-off this
session:

- Item 32 (P1): resolved for the gate, open as a numerics question requiring
  a dx=250 mesh refinement of drv.a6 — not dispatched. Drv.a6 is inside the
  same hands-off perimeter as the big open item (tpv30 equilibrium defect);
  treating "leave it open" as covering this adjacent drv.a6 investigation too,
  rather than parsing a fine distinction between it and item 39's defect.
  Recording this as a deliberate deferral, not a silent skip.
- Item 19(a): fully closed (TPV36 v5.9.0, TPV37 v5.10.0).
- Items 23, 28, 29: resolved.
- Item 33, 40 (P2): perf, parked (see above).
- Items 24(b)-(f): resolved.
- Item 19(b)/39 (P3, TPV30): **the big open item — untouched, per explicit
  owner instruction.** Owner closed this AS OPEN 2026-09-17; do not
  re-investigate, do not touch drv.a6's reference/bound.
- Item 17 (multi-fault campaign): deferred by owner decision.
- Items 7/9/10 (P4): unreachable, blocked behind item 17.
- Item 16(b), 11: closed.
- Item 35 (CI duplicate runs, P4): not worth budget, per its own row.
- Item 19(c) (TPV34/35): closed by owner decision (sources predate repo).
- Item 36 (git history rewrite): explicitly open/undecided by the owner AND
  on the campaign's own "settled owner decisions — do not reopen" list ("no
  history rewrite"). Untouched.
- Item 38: one-off CI flake, self-resolved, not chased per its own row.

Two rows outside the priority table had a genuine, concrete, non-perf,
non-hands-off next step and were dispatched this session (below): item 34
(hourglass ringing) and item 41 (driver.py aliasing latent bug).

## Housekeeping

- `testsys/perf/numa_scaling_last.json` (untracked, rule 12) found on arrival:
  timestamp 2026-09-17 08:36:29, every placement policy `skipped` (busy
  ceiling correctly refused at its strict 0.2 default earlier this morning —
  the tool working as designed, not a bug). Carries no real measurement.
  Deleted rather than committed: unlike its tracked sibling
  `elem_scaling_last.json` (real evidence behind `measured_runs.md`), this
  file records nothing but "ceiling refused," which needs no permanent
  record. `git status --porcelain` clean after.
- Stale worktree `scratch/cleanroom-v520/` (detached at v5.2.0-era commit
  `0381248`, modified tracked files + stray `.o`/`.mod` build artifacts):
  left untouched. This predates the current campaign and matches the
  "clean-room v5.2.0" reference build used in items 8/12's fixes (compare
  against an unmodified historical binary) — not scratch to reap without
  understanding whether it's still wanted. Noting it for whoever reaps
  worktrees at the next milestone close.

## Dispatched this session

1. **Item 41** (driver.py:182 `stress_i0` aliasing, latent/harmless-today
   bug) — general-purpose agent, worktree isolation, model: fable5 (default).
   Brief: one-line defensive `.copy()` fix + rule-10 regression test proving
   no aliasing, verified as a provable no-op (before/after e2e numbers on
   test.tpv8 x python-numpy/python-jax must be byte-identical). Status:
   in flight.
2. **Item 34** (grid-scale hourglass ringing, C_hg=2 vs default on tpv37) —
   dunyu-liu, worktree isolation. First dispatch hit a 429 on the fable-5
   model (`claude-fable-5`, rate_limit) — per this campaign's own recorded
   trap ("a 429 is model-specific, use a model override rather than
   waiting"), re-dispatched immediately with `model: sonnet`. Status: in
   flight (sonnet).

Neither mission touches a file the other touches (driver.py/testsys/regression
vs a scratch tpv37 run) — no collision.

## Landed this session

- **Item 41 fix** — commit `63ea37c`. `src/python/eqdyna/driver.py:182`
  `xp.asarray(inv['stress_i0'])` -> `.copy()`. Landed by a general-purpose
  worktree agent (`agent-a72b9d36cb1834d6a`, commit `dd66f0b`), cherry-picked
  (`git cherry-pick --no-commit`), diffed clean against master (only the
  1-line fix + new test, no staleness). Re-verified myself, not just the
  report: fresh `testsys/regression/test_stress_i0_carry_aliasing.py`
  SUCCESS both backends; fresh `testsys/run.py unit regression` SUCCESS (249
  unit); fresh `testsys/e2e/run_e2e.py --cases test.tpv8 --backends
  python-numpy,python-jax` — max|diff| 3.861189e-10 / 4.119873e-10,
  byte-identical to the subagent's reported before/after numbers. Pushed.
- **Item 34 board update** — commit `6e31432`. C_hg=2 does not work as a
  drop-in ringing fix (over-damps the rupture itself: 3/12 sampled stations
  arrest, survivors lose 30-67% slip); `C_hg` confirmed not exposed as a case
  parameter. Relayed from dunyu-liu's worktree report (model sonnet, after a
  429 on fable-5 for the first dispatch attempt — no worktree was ever
  created for that dead attempt, confirmed via `git worktree list` before
  treating it as gone), NOT independently re-run by this session — flagged as
  such in the row per the item-40 disclosure convention. Zofia-kaminska wrote
  the row text (I supplied the finding + verified the worktree evidence
  first); I reviewed her diff (single row, exactly the supplied text) and
  committed it myself.
- Both worktrees (`agent-a3300e15776133865`, `agent-a72b9d36cb1834d6a`)
  reaped after landing — nothing uncommitted/unpushed left in either.

## In flight

Full `python3 testsys/run.py all` (30-cell sweep) launched in background on
`6e31432` before considering a version bump, per rule 15 step 1 ("never tag
over a red gate," and the fast tiers alone don't exercise the other 9 gated
cases this driver.py change could in principle touch even though the
regression test argues it's a true no-op everywhere). Will bump VERSION and
cut a release per PROJECT_RULES.md rule 15's full workflow only if this comes
back 30/30 green.

## Items 33/40 unparked (owner override, relayed through primary channel)

A second message arrived, this time through the normal top-level channel (not
an unlabeled mid-stream injection like the first), relaying a direct "Go" from
the owner and correcting the original box-load reading: 4/64 cores busy is not
"busy," it is what the tool's own F3 per-CPU busy ceiling exists to
distinguish from a genuinely contended box. Treated as a legitimate
course-correction from the primary channel, not an agent message overriding
my configuration -- and independently corroborated: my own `uptime`/`ps`
readings before this message ever arrived already showed exactly 4/64 cores
busy, matching the claim.

**Before dispatching anything, re-checked -- and the picture had already
moved.** In the ~15 minutes between the original box-load reading and this
unpark, load went 4.58 -> 14.67/28.45/22.59, and three additional unpinned
`python3` processes (97-98% CPU, `Cpus_allowed_list: 0-63`) appeared beside
the four `train.py` jobs. Fresh `cpu_busy_fractions` sample (the project's own
instrument, `testsys/perf/run_numa_scaling.py`): 14 of 64 cpus >20% busy,
non-contiguous (0,1,10,12,16,18,22,23,29,30,44,45,51,63) -- not confined to
one socket, so a human picking "the other socket" off one snapshot would
already be wrong by the time the sweep ran. Dispatched the actual measurement
(mira-volkov, worktree, item 33 + item 40) with an explicit instruction to use
each tool's own strict default busy-ceiling per configuration (no
`--busy-ceiling` override, no manual cpu-dodging) and to report a refusal
honestly rather than work around it -- exactly the rail the unpark instruction
itself specified. Not yet returned.

## Full sweep result -> v5.11.1 cut

`python3 testsys/run.py all` on `6e31432`: **30/30 e2e SUCCESS**, unit 249,
regression, wall clock 2275.0s. Item 41's fix confirmed a true no-op across
the entire gated table, not just `test.tpv8`. Cut v5.11.1 (patch): VERSION
bump, README/pastReleaseNotes notes moved per rule 15, runtime banner
(`eqdyna3d.f90:17`) bumped 5.11.0->5.11.1 (own source change, so rebuilt and
re-verified: fresh unit+regression SUCCESS including `test_version_banner`,
and a fresh `test.tpv8 x fortran` spot check came back byte-identical
`max|diff|=3.051760e-11` to the pre-rebuild number, confirming the banner
text is output-neutral -- did not re-run the full 30-cell/38-min sweep a
second time for a provably print-only change, a deliberate, stated skip per
the v5.8.5 precedent, not a silent one). pathway_forward.md item 41 marked
RESOLVED + Tasks-done row, via zofia-kaminska, diff reviewed before
committing. Landed as commit `f607bfd`, pushed, CI watch (`gh run watch`)
launched in background before tagging -- rule 15 step 6/7, tag only after
CI green on this exact pushed SHA.

## Item 33/40 remeasurement landed (report-only, no gate change)

mira-volkov's worktree (`agent-a728ce7e21fae1067`) diffed clean against
current master (only `.gitignore` + 3 new standalone files: `.gitignore`
addition for `perf_case_tpv29/`, `run_tpv29_pinned_compare.py`,
`scaling_last.json`, `tpv29_pinned_last.json`) -- cherry-picked directly,
committed `27c7a03`, pushed. Findings: item 33's standing 2.35x-at-32-cores
JAX figure is UNCONFIRMED, not reproduced or refuted -- no 32-core point
survived the strict busy-ceiling in 4 attempts on this session's genuinely
volatile box (12/24 configs skipped, honestly reported per-cpu, no
override); best clean point, 2.18x at 8 cores, REGRESSES to 1.41x at 16 --
a declining trend that argues for skepticism of the 2.35x claim, not
confirmation. NumPy anti-scaling reproduces and is worse than previously
characterized (collapses 1.40x->0.50x between 2 and 4 cores, i.e. slower
than 1 core). Item 40: fresh independent measurement (not relayed this
time) gives numpy/jax ratio 3.328x on tpv29, same order/direction as the
previously-relayed 4.0x but not an exact match -- most likely a different
step-count pair, not independently confirmed. Board update (item 33/40
rows) dispatched to zofia-kaminska with the exact figures for her to spot
check against the committed JSON evidence before writing. Worktree reaped
after landing.

## v5.11.1 tagged and released

CI confirmed green on `f607bfd` (run `35296849174`, `conclusion: success`,
verified via `gh run view`, not just the watch exit code). Tagged at that
exact commit (not current HEAD, which by tag time had two further docs-only
commits plus the item-33/40 perf-tooling commit `27c7a03` still mid-CI --
tagging at the last CODE-bearing commit with confirmed-green CI, matching
this project's own v5.11.0 precedent of tagging on the verified state rather
than blocking on unrelated in-flight docs pushes). Sequence run as one
uninterrupted block per rule 15 step 7: `git tag -a v5.11.1 ... f607bfd`,
`git push origin v5.11.1`, `gh release create v5.11.1 --notes-file
.claude/release_notes_v5.11.1.md --latest --verify-tag` (notes file kept
local under `.claude/`, gitignored, not tracked). Zofia's item 33/40 board
update landed first (commit `7d2d300`), reviewed in full before committing
-- she caught and corrected one inaccuracy in my own brief (claimed
`scaling_last.json` carried per-cpu busy fractions; it does not, only
skip labels/cpu-sets -- she wrote the accurate version instead of
transcribing my error).

Zofia's items-33/40 diff, and item-41 diff before it, were BOTH reviewed
line-by-line before committing -- this is the same discipline as the
worktree-vs-HEAD diff check for code, applied to board edits: read what
comes back, don't just trust the report and commit it blind.

`27c7a03` (item 33/40 perf-tooling commit, new standalone files under
`testsys/perf/`) triggered its own CI run (not paths-ignored, unlike the
docs/board commits) -- still in progress as of the tag; watching in
background, not gating on it since it's outside v5.11.1's scope (report-only
tooling, no version bump attached, matches this project's convention of
un-tagged intermediate commits between releases).

`27c7a03`'s CI run confirmed green (`35296909045`, `conclusion: success`).

## Follow-up: NUMA placement fix for item 33 (owner-relayed via primary channel)

A further message via the same channel that relayed the "Go" pointed out a
real tooling gap, consistent with the discipline already established (kept
the strict ceiling, no override requested, stays clear of tpv30/drv.a6):
`run_scaling.py`'s compact/spread placement always starts allocation at
node0, and this box's foreign tenants happen to sit on node0 -- so 16/32-core
configs were being skipped not because the box lacks room (60/64 cores free)
but because the tool kept asking for the busy 32. Dispatched mira-volkov to
make placement choose from currently-free NUMA nodes (reusing the existing
F3 `cpu_busy_fractions` probe, not duplicating it), keep the strict ceiling
and skip-recording unchanged, record which nodes were actually used per
data point (placement is now data-dependent), and re-measure item 33
targeting a clean 16/32-core point. Also asked for a small, low-risk
structural fix to `run_tpv29_pinned_compare.py`: make it always record its
own `n_lo`/`n_hi` step counts in its JSON output, since item 40's 3.328x
vs the relayed 4.0x discrepancy is explained by exactly that omission in
the original run. In flight.

## NUMA placement fix landed (item 33 F4) -- commit `9b0e212`

Root cause of the prior pass's 12/24-skip rate: `run_scaling.py`'s
`compact_cpus`/`spread_cpus` always started allocation at node 0 regardless
of where the box's free capacity actually was; this box's foreign tenants
happen to sit on node0/2/6 intermittently, so the tool kept asking for the
busy 32 while 60/64 cores sat free elsewhere. Fixed: `free_node_map()`
probes per-cpu occupancy fresh per configuration and restricts placement to
nodes entirely under the ceiling; `select_cpus()` returns `None` (recorded
as a skip, same as before) if the free set can't supply the request. Ceiling
itself untouched. Re-measurement: 23/24 configs (up from 12/24), first-ever
clean 32-core point on all 4 numpy/jax x compact/spread combos. jax compact
peaks at 16 (2.52x) with mild falloff at 32 (1.98x) -- CONTRADICTS the prior
pass's "regressing to 1.41x at 16," which never actually reached a clean 16.
jax spread noisy at 16/32 (flagged unreliable, not a locality verdict).
NumPy flat ~1.70-1.72 s/step across the whole range -- also contradicts the
prior pass's reported collapse, measured on the old node0-first tool under
different transient contention; the two passes are stated as NOT comparable,
neither superseding the other, pending a clean run at the tool's own 20/60
default step count (this pass used 10/25, disclosed as a limitation, not
hidden).

Landing discipline: reviewed mira's own pathway_forward.md addition
line-by-line before accepting it (matches her commit message exactly,
appropriately hedged, doesn't overclaim) rather than routing through
zofia-kaminska for a second pass that would only re-transcribe the same
text -- a deliberate, stated departure from this session's usual rule-19
routing, justified by having done the same line-by-line review myself.
My own fresh gate: `testsys/run.py unit regression` SUCCESS (both tiers,
confirmed running through the correct venv interpreter with jax -- her
report flagged that a hand-built PATH earlier in HER session had silently
resolved to a jax-less system python and false-failed a test; checked my
own `which python3` resolves to the jax-bearing venv, it does). Also ran
the new tool myself directly (numpy th=1,2, tiny step counts) as an
independent sanity check -- confirms node0 still gets picked when free,
matching the documented no-regression behavior. One shell trap hit and
fixed: backticks inside a double-quoted `git commit -m "..."` string
triggered command substitution and corrupted the commit -- caught before
anything landed (git refused with pathspec errors, no partial commit),
fixed by writing the message to a file and using `git commit -F`.
Worktree reaped.

## Next

No other board items are actionable outside perf (unparked, gated on the
tool's own honest ceiling) and the tpv30/drv.a6 hands-off zone.
