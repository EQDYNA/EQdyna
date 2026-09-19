# Session log -- 2026-09-19 autopilot (wei-lin)

Resumed mid-flight. Predecessor had vanished from the session registry with
round 3's worktree still held and a `python3 testsys/run.py all` already
running unbounded inside it. State at start: `origin/master` at `1f84521`,
tree clean, 23 commits past `v5.11.1`.

## Round 3 landed -- `7ebe558`

Inherited worktree `agent-a84efd4fb6364632e` was based on current master
(no staleness), with round 3's doc-consolidation edits already made but
uncommitted (3 session logs, 172 insertions / 351 deletions). Did NOT kill
or redo the in-flight sweep -- it was genuine live work (PID confirmed
actively burning CPU, not hung). Let it finish: **30/30 e2e SUCCESS, unit +
regression SUCCESS, wall clock 2171.3s.** Verified the diff was docs-only
(no src/testsys/case_input touched, no `CASE_BOUND`/reference/refusal-logic
change). Committed in the worktree, cherry-picked onto master (fast-forward,
clean), pushed, force-released the worktree's own lock (held by this
session's own resumed PID `2872539`, not a foreign agent -- recorded here,
not a quiet cleanup) and reaped it.

## Board rows (via zofia-kaminska, isolated worktrees, sequential to avoid
two writers on `pathway_forward.md`)

- **Item 42 (new)** -- `98f1020`: local full-sweep gate reliability under
  shared-box contention. Synthesized two data points, stated honestly as
  contradictory-looking but reconcilable: 2026-09-18's session needed 3
  attempts (root-caused to ITS OWN concurrent heavy jobs stacking, not the
  box or cell); this session's inherited unbounded run completed clean
  30/30 in 2171.3s under the box's steady ~25 load average with no
  self-contention. Verdict: the gate reliably reaches 30/30 given (a) no
  timeout under ~40 min and (b) no second heavy job co-scheduled -- a
  procedural fix, not evidence the sweep is unachievable. P3, flagged as
  environmental/procedural rather than silently ranked.
- **Item 36** -- `69b86e8`: corrected from a stale "Left explicitly OPEN"
  to DECIDED (no history rewrite), citing the on-record decision already in
  `docs/SESSION_LOG_2026-09-17_v5.9.0-followon.md:48` that had never been
  reflected in the board row itself. The specific reasoning relayed to this
  session (fork/cloner counts, scec_archive offsite-copy framing) is flagged
  as relayed-not-primary-source, same disclosure pattern as item 19(a).
- **v5.12.0 Tasks-done row** -- folded into the release commit directly
  (not landed as her separate commit) so rule 15 step 5's "commit
  everything together" held for real.

## Release: v5.12.0

Stated grant: tagging on `master`, EQdyna's default branch. Autonomous
mode's base grant is non-default-branch only; this repo's own history shows
~20 prior patch/minor tags cut the same way across this campaign, and this
session's own queue explicitly directed the release on master via rule 15.
Treated as an already-standing widened grant for this repo, stated here per
the disclosure requirement so the maintainer can correct it. Patch/minor
only -- no major bump, no force-tag, nothing outside rule 15's own workflow.

Bumped 5.11.1 -> 5.12.0 (minor: accumulated real fixes -- items 7/9/10,
item 35 -- plus item 33 tooling and the 3-round refactor, not just one
patch-sized change). Gate reused this session's own fresh 30/30 sweep
rather than re-running, since every intervening commit was docs-only
(verified: `git diff --stat a3c18b9 HEAD -- . ...` empty against every
non-doc path). `./install-eqdyna.sh -m ubuntu` re-run fresh, exit 0.

**First release commit (`c782c2e`) failed CI**: `test_version_banner.py`
caught a missed banner bump (`eqdyna3d.f90:17` still said 5.11.1) -- a real
regression guard doing its job, not tightened or bypassed. Fixed same
session (`737358f`): rebuilt, fresh unit+regression SUCCESS, and a
`test.tpv8 x fortran` spot check came back byte-identical
(`3.051760e-11`, unchanged) confirming the banner text is output-neutral --
did not re-run the full 30-cell sweep for a provably print-only change
(v5.8.5/v5.11.1 precedent). CI green on `737358f`
(run `35431026081`). Tagged `v5.12.0` at `737358f`, pushed, `gh release
create --verify-tag` -- annotated tag confirmed via remote peel
(`git ls-remote --tags origin` shows both the ref and its `^{}` peeled
commit). `test_release_complete.py` SUCCESS, all 6 checks.

**Docker publish then failed** (tag push's own `publish.yml` run,
`35433598116`): flagged mid-task by the coordinator, independently
corroborated before acting (own `gh run view --log-failed`, own read of
`Dockerfile`/`publish.yml`) rather than trusted on report alone. Real
interaction, not a regression in either component: item 41's guard
(`test_stress_i0_carry_aliasing.py`) exercises both backends, and the image
deliberately ships jax-free (`Dockerfile:25-27`, sound, unchanged) --
neither had met the other until this tag. Fixed (`9d1977d`): added
`pip3 install jax` to both of `publish.yml`'s CI-only gate overlays,
matching `test.yml`'s own jax-install pattern; did NOT make the guard skip
without jax (rule 2 -- `test_stop_exit_status.py` already burned this
project once on a guard that skipped instead of failing). Verified via
`gh workflow run publish.yml` (`workflow_dispatch`, build+gate, no push) --
SUCCESS -- rather than burning a second tag to test it, per
`Docker.guide.md`'s own stated purpose for that trigger. Corrected
`README.md`'s v5.12.0 block and the already-published GitHub Release body
in place (`gh release edit`), same shape as the v5.8.3/v5.8.4 precedent in
`pastReleaseNotes.md`. **`v5.12.0`'s tag was NOT moved** -- it still points
at `737358f`, whose main `test.yml` CI is green; only the container never
shipped, and the fix rides the next tag. Both follow-up commits (`9d1977d`,
`c726fa7`) confirmed green on `test.yml` before closing this out.

## Final state

`origin/master` at `c726fa7`. Tags: `v5.12.0` (annotated, `737358f`,
Release published, Docker image correction on record) is the newest.
Worktrees: only `scratch/cleanroom-v520` remains (pre-existing, unrelated
to this session, left alone -- not mine to reap without knowing its
owner/purpose).

## Board check (queue's named rows 17/19/24/29/33/34/36/38)

All either closed this session (36, 42-new), already owner-deferred/hands-off
(17, 19(b)/24(a) tpv30 defect, 19(c) TPV34/35, 29 stashes), externally
blocked (38, TPV38 awaiting SCEC TSurf), or open but not actionable without
resources this session doesn't have (33 -- needs a genuinely idle box for
the matched-cores ms/step figure; this box carried ~25 load the entire
session). Nothing else in the queue. Stopping rather than manufacturing
work.

---

# Continuation, same date (wei-lin, new invocation)

Resumed at `origin/master` == local `b47e2fe`, tree clean. Board rows
17/19/24/29/33/34/36/38 re-confirmed unactionable (deferred/closed/hands-off),
consistent with the prior session's close-out above. Standing work per owner
instruction: serial refactor rounds. Dispatched round 4 (`kai-fischer`,
worktree `agent-ad4246f41491e3c61`, isolated) scoped to case_input
tpv36/tpv37 near-duplication, `src/python/eqdyna/` docstrings, and non-perf
`testsys/` helper duplication -- `testsys/perf/` explicitly off-limits, live
item-33 work in progress there (`haruto-nakamura`, worktree
`agent-a98251453a27e55fe`, detached HEAD at `737358f`).

## Item 33/43 landed -- `88f5227`, `07452a6`

haruto reported item 33 complete (relayed via coordinator). Corroborated
independently before acting (nothing taken on the report alone):

- **Metric-bias fix** (`testsys/perf/run_scaling.py`): confirmed
  `737358f` is an ancestor of current master and the file was byte-identical
  at that base vs current master (no staleness) before copying the diff
  over. The bug is real: Fortran's branch reported `wall/n_hi` for a single
  run while jax's branch already differenced two runs to cancel fixed cost
  (interpreter/mesh/MPI/netCDF startup) -- a one-way bias in jax's favour,
  worst at the short step counts this tool uses (measured: old metric read
  1168 ms/step at np=1, new by-difference metric reads 931.65). Fix makes
  Fortran difference two runs the same way. Landed as `88f5227` after a
  fresh `python3 testsys/run.py unit regression` (SUCCESS both tiers -- the
  applicable gate for a perf-tooling-only change, not the full sweep, since
  no solver/`CASE_BOUND`/reference file is touched).
- **Item 43 board row**: zofia-kaminska had already authored it in worktree
  `agent-a9c2f1b6c53ec6308` (base `b47e2fe`, confirmed non-stale for
  `pathway_forward.md` against current master before copying). Verified the
  row's numbers against haruto's report line-for-line before landing as
  `07452a6`: **Fortran overtakes jax between 2 and 4 cores, reaches 14.39x
  speedup at 16 cores against jax's 2.76x plateau** (`test.tpv104`,
  `737358f`/v5.12.0, compact placement, strict busy ceiling). Scatter
  inventory corrected in the same row: 3 sites, all in
  `src/python/eqdyna/backend.py` (line 55 `scatter_add` 21/38 HLO scatters,
  line 63 `setat` 16/38, line 73 `addat` 1/38) -- the earlier "41 sites
  across 5 files" matched docstring prose, not code. Zero of 38 scatters get
  XLA's implicit partition annotation at any core count; fusions partition
  only to round(sqrt(N)). Amdahl on the measured partition counts caps any
  decomposition scheme at 6.7-10x, below Fortran's measured 14.39x --
  recorded as the reason to fix the scatter-add kernel before attempting
  `shard_map`, with the row explicitly marking `shard_map` itself as
  UNVERIFIED (the Amdahl estimate is not a demonstration).

Both worktrees' only diffs were exactly what got landed (`git status
--porcelain -uall` checked before removal, matched); force-removed and
recorded here by name: `agent-a98251453a27e55fe` (haruto) and
`agent-a9c2f1b6c53ec6308` (zofia).

## Round 4 refactor -- still in flight

`kai-fischer`'s worktree (`agent-ad4246f41491e3c61`) delegated its mandatory
full 30-cell sweep to a background job (confirmed genuinely running via PID,
not hung) and will re-notify on completion. Interim diff on disk, NOT yet
reviewed for merge: 7 files, +21/-277, including a new `tpv36_37_common.py`
duplicated **in both** `case_input/test.tpv36/` and `test.tpv37/` -- flagged
for hard scrutiny at landing time (same filename in two places is the shape
of a second duplication under a new name, not a real shared module, unless
the two copies are byte-identical AND a single location is what everything
imports). Not merged; nothing decided yet.

---

# Continuation, same date (wei-lin, third invocation)

Resumed with `origin/master` == local at `572f640`, tree clean (`git status
--porcelain` 0), round 4 gated-and-committed at `45446d0` on branch
`worktree-agent-ad4246f41491e3c61`, rebased onto `572f640`, with its 6
tpv36/37 e2e cells in flight.

## Round 4 landed -- `45446d0` (fast-forward, pushed)

The "duplicated in both dirs" flag above **resolved as a non-issue, by
inspection not assumption**: `case_input/test.tpv37/tpv36_37_common.py` is a
SYMLINK (`git ls-tree` shows mode `120000`), not a second copy. There is one
real 136-line file, in `test.tpv36/`. `create.newcase`'s `shutil.copy`
resolves it, so a created case directory still gets a real 6466-byte file and
stays self-contained.

Gate, all four axes:

- **Axis 4 (stale base)** -- `git merge-base --is-ancestor 572f640 45446d0`
  YES; `git diff --stat 572f640 45446d0` is exactly 9 files, +158/-277, with
  no reverted lines and no dtype/flag change.
- **Axis 1 (no signal removal)** -- read the whole `testsys/unit/` half of the
  diff line by line. It moves `REPO_ROOT`/`sys.path` boilerplate out of four
  test files into `conftest.py`. **Zero assertions removed**, zero tolerances
  touched, one `sys.path.insert` deleted from inside a test body that
  `conftest` now covers. Required the post-edit unit COUNT, not just exit 0,
  because a silently uncollected file reads as a pass.
- **Axis 3 (own re-run, not the report)** -- `python3 testsys/run.py unit
  regression` run by me in the candidate worktree: SUCCESS both tiers.
- **Axis 2 (gate exercises the new path)** -- the new path is the generated
  case inputs, so the e2e cells for both cases on all three backends are the
  gate. 6/6 SUCCESS, exit 0, 1391.1s: `test.tpv36` fortran `0.000000e+00` /
  numpy `7.894064e-09` / jax `5.867293e-09`; `test.tpv37` fortran
  `0.000000e+00` / numpy `6.929140e-09` / jax `5.256410e-09`; bound `1.0e-06`,
  3477 fault nodes each. **Both fortran cells are exactly 0.000000e+00 --
  bit-identical, not merely in-bound** -- which is the real confirmation that
  byte-identical generated inputs produce a bit-identical run. Every python
  figure is unmoved to the last digit from the 2026-09-17 TPV37 gating run
  already on the board.

No tag cut for this landing, matching this campaign's own cadence (the two
prior landings `88f5227`/`07452a6` carry none either; tags come at rule-15
milestones, not per merge).

Worktree `agent-ad4246f41491e3c61` checked before reaping -- `git status
--porcelain -uall` 0, `git log master..HEAD` 0 commits, i.e. nothing unlanded
-- then unlocked, force-removed and its branch deleted, recorded here by name.

## Two accepted risks, flagged at merge time, not discovered after

Neither vetoed the merge; both are on the board as their own row (routed to
`zofia-kaminska`, rule 19):

1. `45446d0` introduces the repo's **first tracked symlink**. On a
   `core.symlinks=false` checkout it materialises as a text stub and
   `case.setup` dies on the import. Ubuntu CI and the owner's boxes are
   unaffected, which is exactly why nothing catches it.
2. `test.tpv37` can no longer survive deletion or rename of
   `case_input/test.tpv36/`, and nothing tests that coupling. The 6 gated
   cells catch it on any sweep, which is why it is a note and not a veto, but
   a dangling symlink is a far less legible failure than a missing parameter.

## Dispatched next

- `zofia-kaminska` (worktree) -- the round-4 landing row plus the new row for
  the two risks above, with the literal 6-cell evidence block and today's
  date. Step 0 stale-base re-sync of `pathway_forward.md` required in her
  brief.
- `kai-fischer` (worktree) -- refactor **round 5**, scoped to
  `src/python/eqdyna/`: docstrings that contradict the code (this project has
  already lost a board item to one -- backend.py's "41 scatter sites across 5
  files" against an actual 3) and local duplication. Hard constraint is ZERO
  numerical change, gated A/B in her own worktree: `test.tpv8` on all three
  backends must reproduce each max|diff| to every digit (jax `4.119873e-10`,
  fortran `0.000000e+00`), and the unit test count must be unchanged.
  `backend.py`'s three scatter sites are code-frozen (open item 43, owner
  decision); `compare.py`/`matrix.py`/`testNameList.py`/references/`perf/` are
  out of scope. Round 6 does not start until round 5 lands or is rejected.

## Board rows landed -- `98cc4fd`

`zofia-kaminska`'s two rows (round-4 landing in Tasks done; **new item 44,
P4**, for the two accepted risks). Her Step 0 re-sync showed no drift under
her. She flagged, correctly without acting on it, that no rule in
`PROJECT_RULES.md` covers a refactor introducing a new filesystem dependency
between previously-independent artifacts, and proposed a rule 18.

## Enforcement first, then the rule -- `8d5261e`, then `4849bfc`

Deliberate order, and the reason is on the record here: this project's own
claim is that its rules are **enforced by `testsys/regression/`**, not merely
written down, so commissioning rule 18 before anything could check it would
have produced exactly the kind of rule `zofia-kaminska`'s own audits flag as
unenforceable-as-written.

- **`8d5261e`** (`iris-vermeulen`) --
  `testsys/regression/test_symlink_integrity.py`, 157 lines. Enumerates every
  tracked symlink from `git ls-files -s` (mode 120000, never a hardcoded
  list) and asserts each is a real symlink, resolves to an existing path,
  stays inside the repo, and -- where `.py` -- actually imports. **Zero
  tracked symlinks is itself a hard FAIL**, so the guard cannot pass quietly
  if its own mechanism breaks. No `skipif`, no swallowing `except`.
  Re-verified by me, not taken on her report: `unit regression` SUCCESS both
  tiers, and I reproduced the **bite** myself by moving
  `case_input/test.tpv36/tpv36_37_common.py` aside -- exit 1, naming the file
  and the resolved missing target -- then restored it and got exit 0 and a
  clean `git status` back. A guard nobody has watched fail is a guard nobody
  knows works.
- **`4849bfc`** (`zofia-kaminska`) -- rule 18 written AROUND that guard, plus
  item 44 updated in place (risk 2 GUARDED; risk 1 only partly, since nothing
  exercises an actual `core.symlinks=false` clone and the owner call between
  a `case.setup` guard and documenting the requirement is still open), plus
  `CLAUDE.md`'s now-stale "17 rules". The rule names the guard for its
  symlink half and states plainly that its non-symlink half (a shared config,
  a directory coupling with no symlink) is enforced by REVIEW ONLY -- an
  honest enforcement gap beats an overstated one. She left the rule counts in
  `pastReleaseNotes.md` and the dated session logs alone: those are history
  as-of-that-date, not drift.

## Round 5 halted on MY error, corrected, re-dispatched

`kai-fischer` stopped before making a single edit, exactly as briefed, because
his measured baseline disagreed with the one I gave him. **He was right and my
brief was wrong.** I told him `test.tpv8 x fortran` should read
`0.000000e+00`; it reads `3.051760e-11`, reproducibly, to six figures across
two independent invocations, well inside its `1.0e-08` bound. I had carried
the figure over from the tpv36/tpv37 fortran cells, which ARE bit-exact for an
unrelated reason (byte-identical generated inputs), and assumed it
generalised. This log already records the true value: the v5.12.0 banner-fix
entry above cites a `test.tpv8 x fortran` spot check at exactly
`3.051760e-11`. I had the number in front of me and used a different one.

Recorded loudly rather than quietly fixed, because the failure it nearly
caused is the expensive kind: had he trusted the brief over his measurement,
the "correct" move would have been to hunt a nonexistent regression, or worse,
to treat a drifting number as normal later. The halt cost ~20 minutes; the
alternative costs a session. His full baseline, now the recorded one:
`test.tpv8` fortran `3.051760e-11`, python-numpy `3.861189e-10`, python-jax
`4.119873e-10`; unit 249; regression now 22 scripts (the symlink guard).

Re-dispatched with the corrected table and a single cheap pre-edit
confirmation (the fortran cell alone) rather than a full three-backend
re-baseline, since nothing landed since his run touches a solver file. One
wasted dispatch on my side too: a fork launched with a placeholder prompt
instead of a directive -- it did nothing and returned nothing, noted here so
its notification is not mistaken later for a real mission.

## Housekeeping -- leftover agent branches

33 `worktree-agent-*` branches had accumulated from earlier sessions whose
worktrees were already reaped. Checked before deleting, not assumed scratch:
`git cherry master <branch>` on each. 30 were fully landed; the 3 that showed
unlanded patches were each verified individually -- an experiment branch whose
net diff against its own base is EMPTY (experiment + its revert), the items
7/9/10 fix branch whose only remaining delta against master is one stale
version-banner line (the fix itself landed as `4d5e58e`), and the v5.12.0
Tasks-done row branch whose row is in master (folded into the release commit,
as this log already records). All 33 deleted; only `kai-fischer`'s live round-5
branch remains. Worktrees on disk: round 5's, and the pre-existing
`scratch/cleanroom-v520` (not mine, left alone).
