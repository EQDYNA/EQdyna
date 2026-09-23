# A/B evidence for `663d758` (jax-MPI per-step barrier moved under the profile)

The raw transcripts behind the numbers quoted in `663d758`'s commit message
and in `docs/NOTES_jaxmpi_ab_2026-09-22.md`. Preserved out of an agent
worktree before it was reaped; they existed nowhere in history.

**FILES** (raw bytes, md5 of the raw log, all round-trip verified)

| file | raw | md5 |
|---|---|---|
| `ab_run.log.gz`  | 6,896  | `4233fb783fa7a6478b900027998093dd` |
| `ab_run2.log.gz` | 13,825 | `5bdde73a7aac9140652a578a5c9d6ae4` |
| `ab_run3.log.gz` | 1,902  | `e4006271270754e28131b16cfaf95013` |
| `parity.log.gz`  | 7,423  | `3429a5808472893dbaed957c93e654a3` |

`ab_run{,2,3}.log` are the three repetitions of the four-arm A/B
(base / b / a / ab) at 4, 8 and 12 ranks on `test.tpv104`, produced by
`testsys/perf/run_jaxmpi_ab.py` and aggregated by
`testsys/perf/agg_jaxmpi_ab.py`. `parity.log` is the
`run_e2e.py --cases test.tpv8 --backends python-jax-mpi` run on the merged
tree.

**WHAT THEY SUPPORT**

    ranks    arm b (landed)     arm a (declined)    noise floor
      4      1.005 +- 0.008     1.050 +- 0.083         5.2%
      8      1.032 +- 0.003     1.003 +- 0.019         1.3%
     12      1.030 +- 0.022     0.969 +- 0.013         1.7%

All 36 points at `EFFECTIVE_CORES >= 0.99` per rank, verified AFTER launch
from the ranks' own stdout. Nothing was rejected at 4/8/12.

**THE BOUND ON THE CLAIM.** 16 and 32 ranks are NOT measured. The 16-rank
base point was rejected at min `EFFECTIVE_CORES` 0.25 — cpus 16-19 each read
`/proc/stat` busy 0.00 immediately before launch and then delivered a quarter
core, the same pathology already on record for cpus 0/1. The wave argument
that motivated arm (a) was made at 32 ranks and is untested above 12. A
sibling session started an e2e sweep at 21:50 which took the box from loadavg
~11 to ~35; the 4/8/12 data all completed by ~21:47, before that began.

**WHAT IS NOT HERE.** Arm (a), the non-blocking halo exchange, is not landed.
It survives unmerged on branch `mira/jaxmpi-step-2026-09-22` (`0ab3fa7`) and
on `mira/jaxmpi-merged-2026-09-22` (`99a3f18`). Neither branch may be merged
wholesale: both predate `193ffd8`/`d7fe2c2` and would revert
`docs/evidence/gate-44afd4a`, `scripts/figures/scec_compare.py`, the v31
evidence, the dx=100 figures and rules 4c/21, and drop 33 perf-ledger rows.

**VERIFY**

    gunzip -c ab_run2.log.gz | grep VERDICT
    gunzip -c parity.log.gz  | tail -5
