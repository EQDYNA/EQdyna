# Conductor handoff, 2026-11-24 (stopped at the 100-call cap, 21:10 CST)

Landed this run (master `8c2fc7f`): #294 board rows 142/169 (`8ad30cb`), #295 ledger rows #289-#293 plus
correspondence line 44 (`386f998`), #296 route five launchers through `mpi_child_env` (`2b9f8df`),
#298 HPC scaling suite dx snap, parity-first job script, vista `mpirun`, venv install (`428b030`, CI run
45534684623 success), #299 route the remaining launchers, `KNOWN_BYPASS` now empty (`8c2fc7f`, CI run
45536266479 success). The invoker was told #298 merged; it rebuilds and resubmits on Vista itself.

## Open work, one line each

- **PR #284 (row 167 step 3), draft, HELD, head `3628339`, branch `wedge/step3-tpv36-37-refs`, worktree
  `.claude/worktrees/agent-a5ff1f811f6cc183e`.** Audit done (no BLOCKER), four-cell re-run done 15:35 CST.
  Owed: fine-grid pre/post overlays. Pre tree = `.claude/worktrees/cond-pre` (`d720df6`), post tree =
  `.claude/worktrees/cond-pr` (`3628339`). dx 250 m done (scratch, volatile):
  `/tmp/claude-16759/-home-utig5-dliu-EQdyna/e618789a-38a5-4235-9ca3-c8a3a0fe9ba6/scratchpad/fine/overlays/`
  (`tpv36_dx250_cplot.png`, `tpv37_dx250_cplot.png`, `tpv3?_dx250_ts_own.png`; runs in `.../fine/runs/`).
  Result at dx 250 m: TPV36 rupture time max |dt| 0.72 s, median 0, 99.8% extent overlap (13868 vs 13888
  ruptured nodes), the difference confined to shallow nodes; station series differ at 0-3 km depth
  (relative 0.96-0.997 at strike 4/8/12 km, dp000) and converge with depth (0.003 at dp120/dp180).
  TPV37 at dx 250 m: identical (max |dt| 0.00 s, station series <= 2e-8 relative), i.e. no wedge effect at
  this dx; check this against #284's 100 changed frt rows per case at the gate dx before merging.
  `scec_compare.py --plot ts-fault` refuses most stations at dx 249.999 m (not on the node grid); the
  same-name station diff is `.../fine/stdiff.py`. Merge plan per the owner brief: full-lane CI, close step 3
  only, row 167 stays open. The PR body still carries the hold marker "NOT merged: pending owner
  acceptance": the owner brief (merge on CI once overlays, audit, re-run exist) and that marker conflict;
  resolve it by asking the owner, do not merge past it silently.
- **Fine-grid queue**: `queue2.sh` pid 233718 (script, log, `queue2.pid` in `.../scratchpad/fine/`) is
  running `post_tpv36_dx200` (started 20:55), then `pre_tpv36_dx200`, `post_tpv37_dx200`, `pre_tpv37_dx200`
  (about 20-25 min each). Do not start a second queue. When it ends, repeat the overlay step at dx 200 m,
  then reap worktrees `cond-e`, `cond-pr`, `cond-pre` and branch `board/cond-2026-11-24e`.
- **Board PR to write** (row 142, one board-only PR): #295 pointer, the Vista failure (solver exit 33, dx
  not snapped, fixed by #298), and the note that sbatch `-o a.eqdyna.log%j` lands in the suite out-dir
  (`#SBATCH -o` does not expand `$SCRATCH`; not fixed in code).
- **Docs PR to write** (fast lane): `docs/user/getting-started.md:40` still says `--submit` refuses on vista.
- **Ledger rows** for #294, #295, #296, #298, #299 are not yet appended to `docs/cycle_ledger.jsonl`.
- **Housekeeping left**: worktrees `cond-vista` (branch `hpc-suite-snap-dx-pr`, merged via #298) and
  `cond-handoff`; remote branch `worktree-agent-a509d8577e86f8671` still on origin (delete failed, re-check);
  `fix/fault-edge-master-node` unexplained, not mine.
- **Asked once, unanswered**: the invoking session never cited a board commit for its Vista instruction.

BLOCKED(owner) and untouched: rows 149, 155, 153; PRs #156, #187, #212; worktrees `board153-ckpt2`,
`agent-ad9168f90bcfe0823`, `tpv13-surface`, `agent-ab1e5b7c96c5121cb`; branch `board/batch-2026-11-23`
(invoking session). No agent is live. No wake armed.
