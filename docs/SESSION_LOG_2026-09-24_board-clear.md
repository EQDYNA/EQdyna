# Session log 2026-09-24 -- "the 7, clear the board" (wei-lin)

Budget as read: close board rows 91, 120, 66, 100, 92, 109, 34 on commands that ran, or restate each so its own
command decides it. Owner rulings relayed mid-session added 94, 56, 76, 87, 44, 112, 73, 49, 63, root cleanup,
docs rule, README rewrite, docs site. Grant: merge + minor/patch tags on master of this repo per rule 25
(squash PRs for src/testsys/.github; docs/board direct). No tag, publish or Pages flip this session without the owner.

## Milestone 1 -- master red fixed; 66 and 100 decided; board carries the queue (22:40)

- **Master red (rule 3a).** CI run 36082084976 on b1c6d95 failed `test.tpv8 x python-jax` on profile_schema's 5%
  sum check (0.803 s of 15.398 s); physics green. The brief said CI was green, which was stale. Root cause, measured:
  (1) the 'write stations' phase (#23) had no bucket; (2) `Profile.phase` started its clock after the entry sync,
  whose first call imports jax (3.45 s here). PR #27 -> f186860: unaccounted 1.1245 s / 35.79 s (3.1%) ->
  0.0008 s / 24.91 s. Guard mutation-checked. The owner session landed #28 (8425260, `SUM_FLOOR_S=2.0`) in parallel.
  Contradiction recorded: the floor admits 13% on the 15 s cell although the attributed gap is 0.003%, and keep/drop
  is owner-held (row 122). Lesson compounded as rule 5b.
- **Row 66 closed.** tpv29, 32 ranks, (4,2,4) vs (8,1,4): max|diff| 9.708830e-15 over 3321 nodes (bound 1e-10); both
  runs vs reference <= 1.24e-14; the y split sat on the fault (7224 vs 3698 raw rows). The row's own command was
  unrunnable (run_e2e fixes ranks at `FORTRAN_RANKS`), so the case build was replicated by hand.
- **Row 100 recorded.** Serial tpv104 peak RSS: fortran 946632 kB, numpy 2398380 kB, jax 3848884 kB
  (was 0.95 / 4.16 / 10.41 GB). Load 33-36/64, so the wall clock is not comparable.
- **Board** f8332d1, **rules** 012150c (26, 1a, 5b), both by zofia-kaminska. Rule 15 "notes lead the README" vs
  rule 26 "3-line pointer" is reconcilable: `test_release_complete.py:178` accepts a pointer line.
- **Collision observed.** The owner session edited the main checkout mid-session (profile_schema floor), against
  rule 21b, and was invisible to `git worktree list` as an agent. Nothing reverted; it landed as #28.
- **In flight:** mira-volkov (94), anya-petrov (125 README), victor-reyes (#29 audit). Branches ready, serial
  behind #29: wei/row44-common-loud, wei/row56-no-device-auto, wei/row87-push-filter (live-checked: its push started
  no run), iris/perf-round2 (91/76/92-probe/112), wei/user-docs-style (lands with 125).

## Milestone 2 -- root, README, row 94, row 44 landed (2026-09-25 10:40)

- **#29 -> 482734d** root cleanup: 10 NOTES_*.md -> docs/notes/, `test_root_allowlist.py`. Audit HOLD fixed first:
  the move had re-pointed `run_jaxmpi_ab.py --notes` at a nonexistent `docs/notes/` path. The tool opens it in append
  mode, so it would silently have started a second ledger. The guard now also fails on a missing allowlisted entry.
  The board-only reference update (242a962) is NOT landed: its cherry-pick conflicts, and 1 stale board ref remains
  (`pathway_forward.md:84` NOTES_tpv3435_spec.md). Routed to the next board pass.
- **#31 -> fee1315** README 382 -> 127 lines, detail moved to docs/user/, `test_user_docs_style.py` (15 negative +
  2 positive controls). Audit HOLD: 4 dropped or wrong facts, restored; the guard missed backticked file:line and bare
  SHAs; 7 case READMEs carried "Provenance: SHA 5e76e8e" (now "shortly after v5.4.0", describe v5.4.0-3). Stranger-clone
  gate on 4f3f520: clone -> install rc=0 -> unit regression SUCCESS -> tpv8 quickstart rc=0.
- **#30 -> 7d8b498** row 94 depth clamp (mira-volkov). Audit HOLD with 6 findings. The fix agent died on the rate limit
  mid-edit, so the conductor salvaged the uncommitted Fortran WIP to `mira/row94-fixes-wip` (d847ccd), and a second
  dispatch finished it. Re-audit MERGE 6/6. Conductor oracle on cbe3457: 8/8 cells
  (tpv8/10/36/37 x fortran/jax), references untouched. The fortran/jax body-file counts are 14/15 on tpv8:
  station 11 sits on the y-partition boundary, a pre-existing Fortran gap and rule-23 divergence that needs a NEW row.
  tpv36/37 stations 20-22 are real outside-mesh drops. New exit code 49 `ERR_MESH_GRID_TOO_LARGE`.
- **#32 -> 5dcf6fc** row 44 loud shared-compset failure.
- **Violation, self-reported:** #32 was opened while #30 was open (rule 25 serial). It was held as a draft on the
  owner's catch and re-synced after #30.
- **Dispatch miss:** row 126 went to anya-petrov, who declined it as outside her role (publication staging); it was
  re-routed to haruto-nakamura. The cost was one cold start.
- In flight: #33 row 56 (re-audit), haruto (126). Queued: wei/perf-round2-land (91/76/92-probe/112), row 87.

## Milestone 3 -- 56, perf round 2, CI trigger + docs site landed; README executed (2026-09-25 13:00)

- **#33 -> 4faf0dc** row 56. The audit caught a CRITICAL: run_e2e's GPU sweep set only JAX_PLATFORMS=cuda, which the new
  cpu default overwrote, so `run.py gpu` would have gone green on the CPU. Fixed by passing --device explicitly. The GPU
  choice is now spelled `cuda` (the owner's words), which reverses my earlier choice to keep `gpu`.
- **#34 -> 7607fd9** rows 91/76/92-probe/112 (iris). Audit HOLD: a backfill of a shard-shaped snapshot filed it as
  tool=run_scaling, the exact mislabel 91a claimed to close. Fixed, and the check was mutation-verified. All 659 historical
  ledger rows still validate; an append without `contention` is refused. Limitation: run_e2e rows sample contention after
  the sweep, box-wide.
- **#35 -> 9ea1d07** rows 87 + 126, batched per the owner. Live evidence for 87: the batch branch's push started no push
  run. Docs site: MkDocs --strict 0.40 s; parameters are generated from defaultParameters.py and checked for freshness.
  Row 126 was first dispatched to anya-petrov, who declined it as outside her role; re-routed to haruto-nakamura.
  **Owner step pending:** Settings -> Pages -> Source = GitHub Actions.
- **Board 7d0f4b7 / rules e1c9234** (zofia): 94/56/44/76/91/112/123/124/125 closed; new rows 127 (station-11
  y-partition gap), 128 (README executed), 129 (main-checkout sweep at 08:42, left untouched for the owner).
- **Row 128, PR #36 open.** The executing gate (iris) found 3 README defects that existence checks could never see:
  `create.newcase ~/runs/tpv8` failed in a fresh HOME; `pip install jax` was prose, not a step; and, in a clean user env,
  `pip install jax` pulled numpy 2.2.6 under apt-built netCDF4/cftime (`numpy.dtype size changed`), so the README now
  uses a rootless virtualenv. The full gate passes on 17cfe3c: 19 lines, exit 0, 77.8 s.
- **Lessons:** my session's python3 is a venv (/home/utig5/dliu/gns/gns/venv_cotopaxi), so no earlier run of mine exercised
  the system python a new user gets. Twice a shared branch checked out in two worktrees left one tree stale under a
  moved ref: 'dirty' was an artifact, and a gate run on the stale tree was invalid.
- Worktrees: 20 reaped after a per-tree check (no uncommitted or unpushed work). Kept: 2 pre-existing owner-held, 2 live.
