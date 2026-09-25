# Perf round 2 -- testsys residuals (2026-09-24)

Branch: iris/perf-round2. Scope: testsys/perf, testsys/e2e/run_e2e.py,
testsys/regression -- no src/ edits, no pathway_forward.md /
PROJECT_RULES.md edits (per brief).

## Item 1 (row 91a, ledger half)
- `ledger.rows_from_scaling_snapshot` takes `tool=` (default 'run_scaling');
  reads `r.get('engine', 'python-jax')` / `r.get('policy', r.get('mode'))`
  so a run_shard_scaling row (no 'engine'/'policy' keys) converts without a
  second near-identical function.
- Shard's `ranks` declared via the EXISTING 'threads' parallelism value
  (engine defaults to 'python-jax', already mapped to 'threads' in
  `_SCALING_PARALLELISM_BY_ENGINE`) -- no new PARALLELISM member invented.
- `run_shard_scaling.py` now calls `ledger.box_tenancy` +
  `ledger.append_rows(ledger.rows_from_scaling_snapshot(..., tool=
  'run_shard_scaling'))` at the end of `main()`.
- Test: `testsys/regression/test_perf_item91_guards.py`, extended
  `check_91a_snapshot_naming_and_ledger` with a synthetic shard-shaped
  meta -> asserts tool/backend/parallelism/mode on the resulting row.

## Item 2 (row 91, wall-clock guard)
- `testsys/perf/run_mpi_scaling.py`: `per_step_jax_mpi`'s
  `ps = (hi[0]-lo[0])/(n_hi-n_lo)` now raises RuntimeError on `ps <= 0`,
  same pattern as the existing `ps_solve <= 0` guard two lines below it.
- Test: new `check_91_wall_clock_guard` in test_perf_item91_guards.py
  (registered in CHECKS). Drives the real function with a faked
  `jax_mpi_once` returning a DECREASING wall clock; also proves a genuine
  increasing case still returns.

## Item 3 (row 91d) -- NOT a fix, reported only
`grep -n "add_argument('--exclude-cpus'" testsys/perf/*.py`:
```
testsys/perf/run_mpi_scaling.py:521:    ap.add_argument('--exclude-cpus', default='',
testsys/perf/run_setup_probe.py:136:    ap.add_argument('--exclude-cpus', default='0,1,16,17,18,19')
testsys/perf/run_jaxmpi_ab.py:220:    ap.add_argument('--exclude-cpus', default='0,1')
```
Three different defaults across three tools. Owner ruling needed; not
changed here.

## Item 4 (row 76, contention field)
- `ledger.py`: new `contention` field, described precisely in the module
  docstring (what's sampled: run_numa_scaling.cpu_busy_fractions, a
  two-read /proc/stat idle-time delta, sample_s=0.3s; scope: the SPECIFIC
  cpus a measurement used, not the whole box). Universal (unlike
  `placement`/`tree_dirty`): required on append for EVERY tool's new rows,
  tolerated absent on read (historical rows in docs/perf_ledger.jsonl
  untouched, still load -- verified: `validate(row, appending=False)` on a
  contention-less row passes).
- New helpers: `_busy_dict` (normalises either busy-check shape in live
  use), `contention_from_check` (built from a tool's OWN pre-flight busy
  check, never resampled -- reuses row 92's sampler per the brief),
  `contention_from_tenancy` (for run_e2e's unpinned cells: the whole box
  IS "the cpus used").
- `box_tenancy` now also returns `cpus` (the ids its reading covers), so
  `contention_from_tenancy` needs no second sample.
- Wired into `rows_from_scaling_snapshot`, `rows_from_mpi_scaling_snapshot`
  (both jax and fortran sub-rows), `rows_from_e2e_results`, and
  `run_jaxmpi_ab.py`'s own hand-built rows (it does not go through the
  shared converter -- found by reading it, not assumed).
- Test: new `testsys/regression/test_perf_row76_contention.py` (7 checks):
  required-on-append / tolerated-on-read, shape enforcement (empty cpus,
  busy>total, total!=len(cpus), non-dict), both busy-check shapes, the
  box_tenancy/contention_from_tenancy pairing, and the run_jaxmpi_ab wiring.
- NOTE: `_rows_from_snapshot_file`'s `backfill`/`reissue` CLI, which seeds
  ledger rows from OLD committed snapshot files, was NOT specifically
  re-verified against a pre-2026-09-24 snapshot that lacks a 'busy' key on
  its rows -- `contention_from_check` would raise a plain `KeyError` (not
  a curated ValueError) on such a snapshot. Flagged, not fixed: no such
  snapshot is backfilled in this session, and the brief scoped this item to
  live capture, not the backfill CLI.

## Item 5 (row 92, busy probe)
- New `testsys/perf/busy_probe.py`: `free_node_map` moved here verbatim
  (docstring kept), plus a `main()` printing the live per-cpu busy map and
  a per-NUMA-node free/total summary.
- `run_scaling.py` now does `import busy_probe; free_node_map =
  busy_probe.free_node_map` instead of defining the function a second time.
  `run_shard_scaling.py`'s `room()` calls `rs.free_node_map` unchanged (same
  object via the module attribute).
- `grep -c "def main" testsys/perf/busy_probe.py` -> 1.
- Test: new `testsys/regression/test_perf_row92_busy_probe.py` (5 checks):
  one main(), `run_scaling.free_node_map is busy_probe.free_node_map`
  (identity, not just equal output), busy-ceiling filtering behaviour,
  override bypass, and the unreadable-utilisation raise.

## Item 6 (row 112)
- (i) `testsys/e2e/run_e2e.py` `cell_cost`/`profile_ranks`: the fortran and
  python-jax-mpi branches now raise ValueError naming the case/table on a
  non-positive configured rank count, instead of `max(1, ...)` silently
  flooring it. The python-jax branch's `max(1, math.ceil(JAX_MEASURED_CORES))`
  (a measured-cores ceiling, not a configured count) is left alone, as
  instructed, and covered by a "still returns a plain int" smoke check.
  Test: new `testsys/regression/test_perf_row112_e2e_ranks.py` (5 checks).
  MUTATION CHECK: reverted `cell_cost`'s fortran/python-jax-mpi branches to
  the old `max(1, ...)` in place, reran the test -- 3 of 5 checks went RED
  (both "raises" checks + the jax-mpi raise), then restored the fix and
  reran green. See gate output below for the restored-green run.
- (ii) `ledger.tree_dirty_once`: memoized per (paths, root) for the life of
  the process. Call sites switched from `profile_record.ledger.tree_dirty()`
  to `profile_record.ledger.tree_dirty_once()` at run_perf.py (3 sites),
  run_scaling.py (2 sites), run_mpi_scaling.py (1 site) -- the exact lines
  the brief named. `per_step_jax_mpi` is also called by run_jaxmpi_ab.py;
  since the memoization lives in the shared `ledger` module (one dict per
  process), that caller is fixed for free without touching it.
  Test: `test_perf_row76_contention.py`'s `check_tree_dirty_once_memoizes`
  (drives the real memoization: the raw tree_dirty is invoked once across
  3 calls) and `check_capture_sites_use_tree_dirty_once` (source-level,
  documented as a rule-10a carve-out: within ONE call there is no
  behavioural difference between the raw and memoized function -- only a
  SECOND call in the same process tells them apart, and driving each
  tool's full main() twice per process is not a fast-tier fit).
- (iii) `testsys/regression/test_profile_env_strict.py`: added a
  blank-only value (`'   '`) and a 33-char value (one over the Fortran
  side's 32-char `envval` buffer, `src/fortran/eqdyna3d.f90:27`) to both
  the Python bogus-value list and the Fortran `for bad in (...)` list.
  Both were ALREADY correctly rejected by the existing code on both sides
  (Python: falls through to the bogus-value raise; Fortran: the 33-char
  value hits the `status=-1` truncation branch, the blank-only value hits
  the `envLength/=len_trim(envval)` "carries blanks" branch) -- this is
  pure coverage, no defect found, confirmed by directly running
  `profile_emit.enabled()` against both values outside the test harness
  (see report) and by the regression script itself passing end-to-end
  against the real built binary.

## Fixture updates forced by the new required `contention` field
Four PRE-EXISTING regression files build synthetic ledger rows/snapshot
metas by hand and needed a `contention`/`busy` fixture value added (not a
weakening -- these are exactly the rows that should now carry the field):
`test_perf_ledger.py` (`good_row`, an mpi-scaling row, a scaling row, two
e2e tenancy dicts), `test_perf_mpi_placement.py` (`base_row`, a snapshot
row), `test_perf_parallelism_discriminator.py` (`good_row`, 3 scaling rows,
1 mpi row, 1 e2e tenancy dict). All four re-verified green after the fix.

## Real capture demo (then reverted)
`python3 testsys/perf/run_scaling.py --case test.tpv8 --skip-fortran
--backends numpy --py-threads 1 --policies compact --n-lo 5 --n-hi 10
--repeats 1 --i-know-the-box-is-busy` (box was ~54-56/64 cpus busy from
concurrent sessions; `--i-know-the-box-is-busy` used only to avoid an
indefinite SKIP-and-retry loop on a shared box, per box discipline "no perf
sweep beyond the one smallest capture"). Appended ledger row (then
reverted, see report) carried `"contention":{"busy":1,"cpus":[0],"total":1}`
-- cpu 0 measured 100% busy at capture, correctly distinct from the
whole-box `"tenancy_busy":54,"tenancy_total":64`. Reverted with `git
checkout -- docs/perf_ledger.jsonl docs/run_profiles.jsonl
testsys/perf/scaling_last.json` plus deleting the new untracked snapshot
file; `git status --porcelain -- docs/` empty afterward.

## Box discipline
OMP_NUM_THREADS/OPENBLAS_NUM_THREADS/MKL_NUM_THREADS=2 exported for every
command in this session. No perf sweep run beyond the one smallest capture
above. Box load ~37/64 at session start (`uptime`), rose to ~54-56/64 by
the time of the capture demo (other users' training jobs and other
worktrees' e2e sweeps, confirmed via `ps aux`) -- unrelated to this work.

## Gate
`python3 testsys/run.py unit regression`: SUCCESS unit (512 passed), SUCCESS
regression (68 scripts, 0 FAIL) -- full log has zero `^FAIL` lines. Four
pre-existing regression files initially went red from the new required
`contention` field (fixture rows lacking it) and one (`test_ci_shard_
coverage`) from the 3 new files being unassigned to a CI shard; all fixed
and reverified individually, then the full gate rerun green end to end.
