# Gate evidence for `d4f625e` (item 61: per-rank XLA compilation-cache directory)

**WHAT:** the full-suite run that gated `d4f625e` before it landed on master,
plus the per-cell record from the same run.

**WHY IT IS HERE:** same reason as `gate-01a1040/` — the run lived only in a
session scratchpad, and a gate that cannot be re-read is a gate taken on trust.

**FILES**
- `run_all_sweep_2026-09-22.log.gz` — `python3 testsys/run.py all`, gzipped
  (2,291,505 bytes raw, md5 `cc2623fff96d57e03824cff956967852`; round-trip
  verified before committing). 31 of 40 cells ran, 31 passed, 0 failed, 9
  declared-unsupported, unit/regression/e2e all SUCCESS, 9413.8 s wall.
- `e2e_cells_2026-09-22_172026_3930207.json` (in `docs/perf_snapshots/`) — the
  per-cell snapshot from this run; its rows are in `docs/perf_ledger.jsonl`.
  Three smaller snapshots from the same session (`_142438_`, `_142530_`,
  `_142555_`) are the single-cell runs described below.

**THE CELL THIS GATE EXISTS FOR:** `test.tpv8 x python-jax-mpi`,
`max|diff| = 1.220703e-10` against bound `1.0e-08` — identical to every printed
digit before and after the change. The point of the landing is not the number;
it is that the number was previously obtainable only with the compilation cache
disabled BY HAND, and is now obtainable with the cache on.

**THE TWO SINGLE-CELL RUNS, which are the part a full sweep does not show:**
the mpi cell was run twice on its own with `EQDYNA_JAX_CACHE_DIR` unset (cache
ON, no workaround) — SUCCESS in 33.2 s and 28.8 s. The second ran concurrently
with a separate serial `test.tpv8 x python-jax` process, which is the condition
the wedge needs, and did not wedge.

**WHAT THIS RUN DOES NOT SHOW.** The hang itself was not reproduced; it is
intermittent and no committed single-command repro exists. The deterministic
coverage of the new path is the three unit tests in
`testsys/unit/test_mpi_decomposition.py`, not this sweep.

**TENANCY.** The box carried several other agents' jobs throughout (load average
~100 on 64 cores), which is why the wall clock is 9413.8 s against
`gate-01a1040`'s 2904.6 s for the same 31 cells. Parity verdicts are
load-independent; the wall clocks in this run's ledger rows are NOT
clean-tenancy numbers and must not be compared against quiet-box baselines.

**VERIFY**
    gunzip -c run_all_sweep_2026-09-22.log.gz | md5sum   # cc2623fff96d57e03824cff956967852
    gunzip -c run_all_sweep_2026-09-22.log.gz | tail -20 # the SUMMARY block
