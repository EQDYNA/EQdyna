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
