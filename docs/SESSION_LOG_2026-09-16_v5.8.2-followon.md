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
port is exact, but `build_solver_state` (`eqdyna3d.py`) now HARD-REFUSES any
run whose mesh contains an elemTypeArr 11/12 (wedge) element, because
`calcGlobalShapeFunc.f90:22-28`'s shape-function merge (`globalShapeFunc(i,3)
+= (i,4); (i,4)=0`, same for 7/8, for elemTypeArr 11/12) is not ported into
`assembleGlobalMass.py`'s `compute_element_shape`. **This means tpv36/tpv37
CANNOT actually be gated yet** — any attempt to run them dynamically on the
Python backends will hit this refusal immediately, since both cases are
100% wedge-degenerate (C_degen=15). The coordinator's latest sequencing
message treats the C_degen port as fully unblocking tpv36/37 gating; it does
not yet, and I am not attempting the gate until this second piece lands.
Scoped and ready to dispatch: port `calcGlobalShapeFunc.f90:22-28`'s
elemTypeArr 11/12 branch into `compute_element_shape`
(`src/python/eqdyna/assembleGlobalMass.py`), verify against a real dynamic
run (not just mesh counts this time — the previous verification couldn't
reach this code at all, since the refusal fires first).

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

## Constraints still in force

Do not gate tpv36/tpv37 until the dynamics-kernel gap above closes. Do not
touch `/home/utig5/dliu/consilium`. Do not delete
`scratch/mira_g6debug/serial_noyield`. Item 33 needs an idle box (do not
override `run_numa_scaling.py`'s refusal). Item 34: no fourth ringing
mechanism without a controlled experiment. Item 29's stash: flag, don't pop.
