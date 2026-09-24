# Item 91 checkpoint notes -- six testsys/perf/ defects

Branch: iris/item91-perf-tools. Worktree-isolated (confirmed
`git rev-parse --git-dir --git-common-dir` prints two different paths).

## Sub-items fixed

- **(a)** `testsys/perf/run_shard_scaling.py:70` -- `%Y-%m-%d` -> `%Y-%m-%d_%H%M%S`.
  Wired ledger append (`ledger.rows_from_scaling_snapshot` + `append_rows`),
  reusing the existing generic reader (no edit to `ledger.py`, out of scope):
  rows now carry `engine='python-jax'`/`policy=<mode>` so the shared
  converter reads them. Also fixed the same silent-baseline pattern in the
  final fortran/jax summary table (`el[min(el)]` -> explicit `el_base_n`).
- **(b)** three more per-step-by-difference sites, enumerated by reading
  every `(t_hi - t_lo) / float(n_hi - n_lo)`-shaped line in the six in-scope
  files:
    - `testsys/perf/run_numa_scaling.py` `per_step_and_fixed` (~line 293,
      shared by this file's own `per_step` AND `run_shard_scaling.per_step`,
      which calls it directly -- one fix, two call sites).
    - `testsys/perf/probe_scatter_bandwidth.py` `per_iter` (~line 175).
    - `testsys/perf/run_perf.py` `steady_state_per_step` (~line 145).
  All three now raise `RuntimeError` on `<= 0`, matching
  `run_mpi_scaling.py:415-419`'s reference (read, not edited).
- **(c)** `testsys/perf/run_numa_scaling.py` `time_one` (~line 244-284):
  bare `taskset -c` -> `numactl --physcpubind=... --membind=...`, membind
  nodes derived from a new pure `_nodes_of(node_map, cpus)` helper. Signature
  gained `node_map`; confirmed (grep) it has no external callers besides this
  file's own `per_step`/`main`.
- **(e)** `run_jaxmpi_ab.py` (~line 228) and `run_setup_probe.py` (~line 154):
  git-sha subprocess result now passed through a new pure `require_git_sha`
  in each file, which raises `SystemExit` on `returncode != 0` or empty
  stdout instead of returning `''`.
- **(f)** `run_perf.py` main() (~line 217): the affinity WARNING-then-proceed
  replaced with a new pure `require_affinity_pinned`, called from main(),
  which raises `SystemExit` on mismatch.
- **(g)** grepped `speedup` across `testsys/perf/*.py`; found the silent-
  first-non-skipped-baseline pattern in `run_shard_scaling.py` (per-mode
  `base.setdefault`) and `run_numa_scaling.py` (`if base is None: base = ps`).
  Factored each into a pure function (`record_point` / `record_baseline`)
  that returns `(speedup, baseline_n_or_label, note)`; `note` is non-None
  exactly when the baseline point is not the intended first one, and is
  printed AND carried onto every row/result (`baseline_n` / `baseline_label`)
  so the snapshot never silently swaps baselines. `run_scaling.py` has the
  identical pattern (`base = None` at two sites) but is OUT OF SCOPE (not in
  the files-I-may-edit list) -- left unfixed, reported below.
  Same identical-pattern check for `run_jaxmpi_ab.py`'s speedup (per-arm,
  not per-n) showed it already prints an explicit "NO USABLE BASELINE"
  verdict when the base arm is rejected -- not a defect, left untouched
  (sub-item e only for that file per the mission).

(d) is explicitly out of scope (owner ruling needed).

## Guard

`testsys/regression/test_perf_item91_guards.py` -- one file, one check
function per sub-item (91a/b/c/e/f/g), all behavioural: real functions
driven with inputs that produce the bad case (monkeypatched subprocess/time
boundaries, no MPI, no box-state dependence, no real git/numactl/jax).
Registered in `testsys/ci_shard.py` shard "1".

RED-without-fix verified for every sub-item by temporarily reverting the
exact fix (sed/python one-liner), re-running the guard, observing the
matching check(s) FAIL, then `git checkout -- <file>` to restore. Evidence
lines are in the final report.

## Not fixed / deferred

- `run_scaling.py`'s two `base = None` speedup-baseline sites (same 91g
  shape as run_shard_scaling/run_numa_scaling): file not in the SCOPE list
  for this mission, left untouched. Flagged for owner/next session.

