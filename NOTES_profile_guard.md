# NOTES_profile_guard.md -- checkpoint trail

Mission: second half of item 2 (always-on per-rank profile). Emitter is
mira-volkov's (`origin/mira/profile-emitter-2026-09-23` @ c5b4d7d, rebased
onto master here). This mission: the mechanical GUARD, the COLLECTION
(append-only ledger + tree-dirty field), and the ZERO-COST GATE on top of it.
Working on `origin/iris/profile-guard-2026-09-23`.

## Checkpoint 1 -- setup
- `git fetch && git reset --hard origin/mira/profile-emitter-2026-09-23 &&
  git rebase origin/master` -> clean, HEAD d063223 (master was e93e8f7,
  touched only docs per the brief).
- Read docs/run_profile.md, src/python/eqdyna/profile_emit.py,
  src/fortran/library_output.f90 `output_profile` (#10, :407-471,
  EQDYNA_PROFILE=0 check at :428), and the four untested testsys files
  (profile_schema.py, profile_record.py, profile_query.py,
  profile_overhead.py) in full.
- Built src/fortran fresh (`./install-eqdyna.sh -m ubuntu`) -- no prior
  bin/eqdyna in this worktree.

## Checkpoint 2 -- affinity finding (owner asked to check before trusting
any pinned python number)
- `run_e2e.py` itself has ZERO occurrences of `sched_get/setaffinity`,
  `numactl`, or `affinity` (grepped) -- it pins nothing and checks nothing;
  it relies entirely on OS-level fork/exec inheritance from however it was
  itself launched.
- Reproduced NOTES_profile_emitter.md's own Checkpoint 7 finding (a
  python-numpy subprocess showing `Cpus_allowed_list 0-255` despite a
  numactl-pinned parent) using a bare `python3` on this shell's PATH.
  ROOT CAUSE FOUND: this shell's PATH resolves bare `python3` to
  `/home/utig5/dliu/gns/gns/venv_cotopaxi/bin/python3` (an unrelated
  project's venv), ahead of `/usr/bin/python3`. Isolated tests (plain
  subprocess, ThreadPoolExecutor worker thread, a real numpy/torch import
  AND a real BLAS/torch op) all showed CORRECT 16-cpu-mask inheritance
  through that venv, so the venv itself does not widen affinity.
  Re-running the EXACT SAME real `run_e2e.py --cases test.tpv8 --backends
  python-numpy` invocation with `/usr/bin/python3` explicit and a strict
  PPID-verified child lookup (not a bare `pgrep -f`) showed CORRECT
  behaviour end to end, twice (cpu 21, then cpu 17, both inside a
  `numactl --physcpubind=16-31` pin). Conclusion: the "0-255" reading (both
  mine and mira's) is a MEASUREMENT ARTIFACT of `pgrep -f "python3 -u -m
  eqdyna"` matching an unrelated, unpinned `eqdyna` process from another
  concurrent session on this shared, busy box (load 30-34 all session) --
  not a defect in run_e2e.py or `_narrow_numpy_affinity`
  (src/python/eqdyna/eqdyna3d.py:103-164, which only ever narrows a mask,
  never widens one; verified by reading it in full). Logged as two entries
  in ~/code/papercuts.md (pgrep-without-PPID-check; stray venv on PATH).
- SEPARATE, real, already-known gap (found by a peer session, papercuts.md,
  2026-09-23): `numactl --cpunodebind` wrapping the WHOLE `run_e2e.py`
  process does NOT confine its `mpirun`-launched cells (fortran,
  python-jax-mpi) -- OpenMPI's own default per-rank binding remaps them,
  observed as cpus_allowed [0-7]/[8-15]/[16-23]/[24-31] under a 16-31 outer
  pin. Confirmed in my own 4-backend evidence run below (fortran/jax-mpi
  fixture cpus_allowed span 0-31, not 16-31) -- expected, not a new defect.
  My OWN new tool (`testsys/perf/profile_overhead.py`'s MPI arms) does NOT
  have this problem: it pins each rank explicitly via
  `OMPI_COMM_WORLD_LOCAL_RANK` + `--bind-to none` (the same technique
  `run_scaling.fortran_cmd` already uses and documents), not an outer
  numactl wrapping the whole mpirun.

## Checkpoint 3 -- mechanical guard (owner requirement 1)
- `testsys/regression/test_profile_guard.py` (new): 9 checks, GREEN on real
  emitter output from all 4 backends (committed fixtures under
  `testsys/regression/fixtures/profile_guard/<backend>/`, copied verbatim
  from one real `run_e2e.py --cases test.tpv8 --backends
  fortran,python-numpy,python-jax,python-jax-mpi --term gate` sweep on this
  box, this session), RED on: a missing rank file, a missing field, a
  malformed field type, and a bucket inflated 511x (both the
  float-precision remainder check and the wide 5% tolerance check tested
  independently). All 9 PASS. Verbatim output in the final report.
- Wired the SAME guard into the real e2e path: `testsys/e2e/run_e2e.py`'s
  `run_one()` now calls `profile_record.capture_run(...)` for every cell
  that PASSES physics comparison -- which validates (via
  `profile_schema.validate_run_dir` inside `capture_run`) AND appends in one
  call. No `_or_warn` wrapper: a profile defect can turn a physics-passing
  cell into a failed one, matching `profile_record.capture_run`'s own
  documented contract.

## Checkpoint 4 -- collection (owner requirement 2)
- `docs/run_profiles.jsonl` created (did not exist before this session).
  Verified end-to-end with a real `python-jax` x `test.tpv8` e2e cell: one
  621-byte row, `cpus_allowed`/`numa_nodes` correctly reflecting the real
  pin, `tree_dirty: true` (this session's own edits -- expected).
- `tree_dirty` field added to BOTH records:
  - `testsys/profile_record.py`: unconditionally REQUIRED (zero legacy rows
    exist yet), computed via a NEW `ledger.tree_dirty()` (one
    implementation, reused, not duplicated).
  - `testsys/perf/ledger.py`: scoped to `tool='run_e2e'` (the tool the
    incident happened in), required on append for that tool only, tolerated
    absent on read for every other tool and every legacy row -- the exact
    same contract this file already uses for `platform`/`parallelism`/
    `placement`. `rows_from_e2e_results` and `run_e2e.py`'s `_perf_meta` both
    updated. Verified live: after these edits, the (deliberately dirty,
    mid-session) tree produced a `docs/perf_ledger.jsonl` row with
    `tree_dirty: true` -- the exact scenario the field exists to flag.
  - `testsys/perf/profile_query.py`: `--exclude-dirty` added to all three
    subcommands.
- Wired into `testsys/perf/run_perf.py` (one of the three named perf tools):
  captures a profile row for each of fortran/numpy/jax's single-`nsteps`
  run (before `steady_state_per_step`'s n_lo/n_hi runs overwrite the same
  case directory), term='perf-pinned-single-core'. Verified live: 3 new
  rows appeared in `docs/run_profiles.jsonl` mid-run, and `run_perf.py`
  itself still printed `SUCCESS perf` (ratio 0.581 vs baseline 0.626,
  degradation 0.929x, well under the 1.5x gate) -- the wiring did not
  disturb the tier's own verdict.
- NOT wired this session: `run_mpi_scaling.py`, `run_scaling.py`. Both are
  large (684/662 lines), have their OWN case-build/lock machinery I have
  not fully read, and I judged a rushed, unverified edit to either riskier
  than leaving them as an explicit gap. Recommendation for the next
  session: `run_scaling.py`'s `time_one_py`/`run_fortran` already isolate
  one case dir per (n, policy) -- call `profile_record.capture_run` right
  after the LARGER (n_hi) run in `per_step_py`/`per_step_fortran`, before
  `shutil.rmtree` (fortran path) deletes the directory, with
  `term='perf-scaling-probe'`. Same shape for `run_mpi_scaling.py`.

## Checkpoint 5 -- zero-cost gate (owner requirement 3)
- `testsys/perf/profile_overhead.py` extended: `--backend` now accepts
  `fortran` and `python-jax-mpi` (`--ranks`, default 4) alongside the
  existing `python-numpy`/`python-jax`. MPI arms difference EACH RANK'S OWN
  solve time (CLAUDE.md's explicit convention), read from an
  EQDYNA_PROFILE-INDEPENDENT source on both arms of the switch: fortran's
  `compTime<rank>` file (gated on `writeCompTime`, hardcoded 1, not on
  EQDYNA_PROFILE) and jax-mpi's own unconditional `driver.run_mpi rank R:
  T s,` stdout line (driver.py:670). This is a deliberate departure from
  `run_scaling.py`'s OWN Fortran scaling gate (which differences the outer
  mpirun wall clock) -- that tool answers a different question
  (rank-count scaling); this one answers "did turning profiling on cost
  this run anything", which is exactly the question CLAUDE.md's note is
  about.
- BUG FOUND running this for real, first time ever executed: `timed_run`
  passed the matrix.py-style backend label (`python-numpy`/`python-jax`)
  straight into `eqdyna3d.run_case`, which expects the internal name
  (`numpy`/`jax`) -- raised `ValueError` from `_resolve_solver` on the very
  first real invocation. Fixed (one-line engine-name map, mirroring
  `run_e2e.run_standalone`'s own mapping) -- this is test tooling
  (testsys/perf/), not src/, so fixed directly per this mission's mandate.
- `testsys/run.py`: new opt-in tier `profile-overhead`, gated on
  `EQDYNA_PROFILE_OVERHEAD_CPUS` being set (refuses, does not default, if
  unset), runs all 4 backends in sequence, never concurrently.
- Real numbers: see final report (this file is the trail, not the ledger).

## CPU discipline log
- Every heavy run this session pinned to `numactl --physcpubind=16-31
  --membind=2,3` (nodes 2-3), one at a time, tracked by PID, never
  `pkill -f`. Box load stayed 30-34 throughout (a genuinely busy shared
  box, matching the brief).
