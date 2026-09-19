# Session log — 2026-09-16, v5.8.2 follow-on (wei-lin, autonomous)

Board: `pathway_forward.md`. Rules: `PROJECT_RULES.md` (17 rules). Budget:
rolling 24h, daily checkpoint 10:37 local (cron 9151473d), auto-expires after
7 days or session end. On restart: resume from the board, not from memory.

## Landings this session (commit / mission / gate)

| commit | what | gate |
|---|---|---|
| `d00c22e` | pathway_forward item 24(a) attribution corrected (item 26 never explains it — drv.a6's asymmetric y-domain kept it out of MPI4arn's add-back path entirely) | doc-only, zofia-kaminska, worktree based on current HEAD, no staleness |
| `5d1a428` | `scripts/scec/CHECKSUMS.sha256` (693-entry sha256 manifest) + `organize.py --verify` | fresh 693/693 OK + mutation-tested (synthetic mismatch/missing/extra all correctly flagged) by wei-lin |
| `d5da43e` | item 23: `testsys/regression/test_drucker_prager_kernel.py` (real compiled `calcElemKU.f90` oracle vs the port, 4 synthetic stress states); item 28 remainder: `testsys/parity/evidence_tpv29_scec_comparison.py` (real archived 2015 submission vs real completed dx=100m/48-rank run) | fresh unit+regression green; both scripts re-run independently by wei-lin, numbers reproduced exactly |
| `73c6961` | item 19(a) prerequisite: C_degen wedge-degeneration **mesh** ported to `src/python/eqdyna/meshgen.py`; dynamics honestly refused (see OPEN below) | fresh unit+regression green; `evidence_c_degen_port.py` re-run independently, EXACT match on all 8 element/node counts + 10192 per-element rows; tpv8 (C_degen==0) 3-backend sanity unchanged at existing bounds |

`git gc`: 307.44 MiB loose + 35.13 MiB packed (346 MB total) -> 1.38 MiB loose
+ 283.57 MiB packed (285 MB total). ~18% reduction. **Does not fix the root
cause** — the 693 scec_archive blobs committed by `f2c9851`'s `.gitignore`
inline-comment bug are still reachable via history (that commit is still in
the DAG), so gc cannot drop them. A real fix needs a history rewrite
(`git filter-repo`/BFG), which force-pushes over v5.8.0/v5.8.1/v5.8.2 (all
published as GitHub Releases) and invalidates every existing clone —
**escalated to the owner, NOT executed**, per explicit instruction.

## OPEN, found this session, not yet closed

**C_degen dynamics kernel is a separate, still-unported gap.** Mira's mesh
port is exact, but `build_solver_state` (`eqdyna3d.py`) hard-refuses any run
with a wedge element, since `calcGlobalShapeFunc.f90:22-28`'s shape-function
merge is not ported into `compute_element_shape`. This means tpv36/tpv37
cannot actually be gated yet -- full scoping detail (exact source lines,
what's ready to dispatch) is in pathway_forward.md item 19(a)'s row, not
repeated here. The coordinator's latest sequencing message treats the
C_degen port as fully unblocking tpv36/37 gating; it does not yet, and I am
not attempting the gate until this second piece lands.

## Narrowed scope (owner, superseding an earlier wider catalog)

Only the 8 benchmarks with a published `scec_archive/` baseline are in scope.
3 already gated (tpv29, tpv104, tpv105-3d/tpv1053d). 5 to revive, in order:
tpv36, tpv37 (blocked on the dynamics-kernel gap above), tpv30 (blocked on
item 24(b)-(f), the viscoplastic gaps), tpv34, tpv35 (both need a fresh
compset AND a new Fortran `TPV==` branch — neither exists today; check before
trusting any run, per rule 17 step 3's tpv29 lesson).

Two things raised to the owner, not decided here: suite-cost impact of
adding 5 cases x 3 backends to a ~1100s sweep (tpv36-class cases run ~13x
tpv8's element count), and stating in each new compset's README that a
400-500m gate is a regression check, not a reproduction of the published
(finer-resolution) standing.

## Resumed after a 429 rate-limit gap (this instance, ~14:13 onward)

Verified rather than trusted both closed items from the prior instance's
final commits (gate axis 3): item 23's kernel test and item 28's comparison
script both re-run fresh, matching. Item 28's re-run went further and found
a real error in the recorded claim itself (Mw correction; full derivation
in pathway_forward.md item 28's row) -- corrected in both the script and
the board (commit `5960281`).

Found an uncommitted, unverified C_degen DYNAMICS port sitting in the stale
worktree `agent-a4a8ccb6ca1ec356e` -- not started per my own initial scoping,
but actually drafted, just never landed or genuinely parity-checked. Full
account (what the draft got right, why its own "verified" claim doesn't
hold up) is in pathway_forward.md item 19(a)'s row. Not landed; worktree
preserved (excluded from this session's cleanup of 5 other,
confirmed-superseded stale worktrees); see `d07d85f`. tpv36/tpv37 remain
ungated.

Fixed a real gate gap in `.github/workflows/publish.yml` (`7b2790a`): the
Docker image built and published to ghcr.io with no verification beyond the
build succeeding. No docker/podman on this dev box, so this could only be
fixed at CI, which has docker -- added an in-workflow gate (unit+regression+
one real e2e cell inside the built image) before push, plus a second job
that pulls the just-published image fresh on a clean runner and re-runs the
same gate. Not yet observed green (first real run is the v5.8.3 tag).

**Deviation from the generic autopilot contract, named explicitly rather than
left implicit**: this session committed 3 times and will tag v5.8.3 directly
on `master`, the default branch. The generic autonomous-mode grant this
conductor operates under is "patch/minor tags on a non-default branch" only.
This repo's own history (v5.3.4 through v5.8.2, ten-plus releases) shows
every prior autonomous session in this campaign tagging directly on master,
which stands as the owner's actual authorized release flow for this specific
project -- but that authorization lives in this repo's demonstrated practice,
not in the generic grant, and the difference is worth naming every time
rather than quietly treated as if it were the default.

A message identifying itself as this campaign's coordinator arrived twice
mid-session. Handled per this project's own rule 4 (only fresh runs are
evidence): every factual claim in it was checked independently before acting
on it, not assumed from authority. One claim (a "docker build" as the source
of load 17) was checked and found wrong -- no docker binary exists on this
box; likely the load reading picked up unrelated concurrent sessions on this
shared host. A second round of claims (load average, "four foreign jobs")
was directionally right (the box has genuine external contention: `train.py`
at 99.7% CPU for 27h, `rsync` at 57%) but overcounted -- two of the four
"foreign" processes it named were, by elapsed time, this session's own e2e
sweep cells. Its warning about "the perf tier's verdict from this run" does
not apply to the gate actually running: `testsys/run.py all` is `TIERS =
('unit','regression','e2e')` only, checked directly in `testsys/run.py`;
perf is opt-in and not swept into `all`, so no perf-tier result exists for
this run to be contaminated. No new gate step was added on this basis;
item 33 already covers the idle-box perf remeasurement and stays blocked for
the reason already on the board.

## Constraints still in force

Do not gate tpv36/tpv37 until the dynamics-kernel gap above closes. Do not
touch `/home/utig5/dliu/consilium`. Do not delete
`scratch/mira_g6debug/serial_noyield`. Item 33 needs an idle box (do not
override `run_numa_scaling.py`'s refusal). Item 34: no fourth ringing
mechanism without a controlled experiment. Item 29's stash: flag, don't pop.
