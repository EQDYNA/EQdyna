# Session log — 2026-10-06 — TPV34/TPV35 validation (rule 17 step 6)

## Budget read (stated before spending any of it)

Owner's word, verbatim: "validate!" (PR #95, e358920, board row 19(c)). Read
as authorizing exactly the validation scope laid out by the monitor's
2026-10-06 read: matched-resolution, full-term runs of TPV35 (100 m, 18 s)
and TPV34 (50 m, 20 s — measure-first, stop and report the LS6 option if it
doesn't fit this box), each compared against its own SCEC submission by a
committed script, plus confirming the TPV35 near/far-side sign mapping and
backfilling the missing cycle_ledger rows + this session log. Not authorized:
a release, a version bump, or scope past that (e.g. fixing the TPV34 Medium
audit finding or the material-grid-coverage assertion — those stay open
follow-ups per row 19(c)).

## Orientation (before dispatch)

- Board: `pathway_forward.md` row 19(c) — both TPV34 (PR #93, 809a023) and
  TPV35 (PR #91, 05b4c1e) landed gated at 500 m/`GATE_TERM_S`; rule 17 step 6
  explicitly open for both (TPV35: not done at all; TPV34: done but at
  mismatched resolution, so its front-lag/slip-gap numbers can't be separated
  from resolution).
- Disk: this box's network home (`utig5.ig.utexas.edu:/bigpool/home`) has
  23 T avail of 172 T; local `/` has 666 G avail of 900 G. The task brief's
  "370 GB avail / meshnet holds 360 GB" does not match either measured mount
  on this box as of this session (`/home/staff` is a symlink to the same
  172 T pool, not a distinct constrained filesystem) — flagged as a
  discrepancy, not reconciled; the actual constraint for TPV34 at 50 m is
  CVM-H grid-extraction output size/time (next section), measured directly
  rather than inferred from the stale figure.
- `meshnet_opt.train`/`meshnet.train` processes (dliu, multiple PIDs, running
  since Oct 02-04 + two fresh rollout jobs started 07:14/07:21 today) confirmed
  live via `ps aux`; left untouched per the owner's order.
- TPV35 is the lighter lift: `tpv35Tools.faultGridForCase` already accepts
  `dx=100` (the shipped spec resolution) with no new data needed —
  `requireFaultGeometryResolution`'s `availableDx` is an error-message hint,
  not a restriction; the real check is "dx is an integer multiple of, and not
  finer than, the 100 m source." Only `par.dx`/`par.term` need to change from
  the gate's 500 m/5 s.
- TPV34 is the heavy lift: `n2mat==6` accepts **only** `par.dx=500` today
  (`case_input/test.tpv34/README.md`) because the shipped CVM-H grid was
  extracted once at 500 m (120×64×52 = 399,360 rows). A 50 m grid needs a
  fresh `vx_lite` extraction, ~1000× more rows (3 axes × 10) — this is the
  real size/time question, to be measured before committing to a full run,
  per the owner's "measure first" instruction.

## Missions dispatched

- TPV35 100 m/18 s run + evidence script: worktree `.claude/worktrees/agent-ae27ab786de91a9e9`,
  branch `worktree-agent-ae27ab786de91a9e9`. Script merged as PR #99 (`283a555`).
- TPV34 100 m/20.02 s run + extended evidence script: worktree
  `.claude/worktrees/agent-ab970015c048c7ae4`, branch `worktree-agent-ab970015c048c7ae4`.
  Script merged as PR #100 (`ef42c11`).
- Owner rulings on row 19(c) carried through fast-lane board PRs #94-#98
  (bookkeeping, "validate!", "100 m suffices", "go").
- Conductor (direct, not dispatched): found and fixed the TPV35 cplot overlay
  defect (tricontour over the ruptured-only point set triangulates the convex
  hull of a non-convex patch, bridging unruptured bays and drawing spurious
  diagonal contours through them) by masking triangles whose longest edge
  exceeds 2*dz, in `testsys/parity/evidence_tpv35_scec_comparison.py`'s
  `make_overlays`. Re-ran both scripts fresh against the existing run dirs
  (not re-simulating) before merging either PR — both reproduced the briefed
  numbers to the same significant figures.

## Numbers

**TPV35** (100 m/17.87 s fortran, 32 ranks, vs `eqdyna3d-v4.1.3-100m-2017`,
matched resolution): rupture-time median dt -0.003 s, max |dt| 0.006 s;
initial stress ratio 1.0000; hypocentre h-slip 0.616 vs 0.617 m; off-fault
peak h-vel ratio median 1.003 over 41 named stations (outlier 0.571 is inside
the PML past x=20 km); ruptured fraction 0.266 both sides. Near/far-side
sign mapping confirmed against `TPV35_Description_v05.pdf` Part 2 text.

**TPV34** (100 m/20.02 s fortran, 16 ranks, vs `eqdyna3d-v3.2.3-50m-2016`,
NOT matched resolution — owner-accepted): rupture-time median |dt| 0.006 s,
p90 0.022 s over 45451 shared nodes; ruptured fraction 0.9868 (run, 100 m)
vs 0.9934 (2016, 50 m); final h-slip ratio median 0.996; peak h-slip-rate
ratio median 0.946; initial stress ratio median 1.0000; off-fault peak h-vel
ratio median 0.986 over 56 stations. CVM-H 15.1.1 (this grid) vs the
submission's 15.1.0.

Rule 17 step 6 is now met for both cases (board row 19(c) updated,
`pathway_forward.md`). Both scripts are report-only and never wired into
`testsys/run.py` or any gate.

## Open items carried forward (unchanged by this session unless noted)

- TPV34 Medium: nothing enforces the `n2mat==6` material grid box covers the
  mesh box (owner/Zofia follow-up row, untouched here).
- TPV35 advisory: near/far-side sign mapping inferred from filenames —
  addressed in this session (see Numbers) — CLOSED.
- Leftover from a prior, unrelated session: local branch
  `worktree-agent-a19871c37d5317fc0` carries 4 commits
  ("evidence: conductor re-verification of test.tpv34 (fresh run, this box)")
  not on master and with no open PR (`git cherry master` shows all 4 as `+`).
  Not touched this session — flagged for the next conductor to triage
  (keep as evidence + board row, or confirm superseded and delete).
