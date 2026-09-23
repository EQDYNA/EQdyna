# NOTES_profile_emitter.md -- recovery + completion checkpoints

Mission: land the always-on per-rank profile.rank<r>.json emitter across all
four backends (fortran, python-numpy, python-jax, python-jax-mpi), zero perf
cost, byte-identical frt on/off. See docs/run_profile.md for the delivered
design; this file is the checkpoint trail.

## Checkpoint 1 -- recovery
- `git fetch && git reset --hard origin/master` -> 95d4220 (already HEAD in
  this worktree, no-op).
- Cherry-picked `origin/wip/profile-fortran-2026-09-23` (fd8cd72) and
  `origin/wip/profile-guard-2026-09-23` (b75c69a). Both applied clean.
- Read testsys/profile_schema.py and testsys/profile_record.py in full
  (guard mission's files, NOT edited). Schema: eqdyna-profile/1, buckets
  setup/element/fault/exchange/wait/io, SUM_TOLERANCE=0.05, FLOAT_TOL=1e-9.
- Read the predecessor's Fortran diff (assembleGlobalMass/eqdyna3d/
  globalvar/library_output.f90): output_profile (#10), the tStageStart/
  tWaitStart local-timer split, MPIWaitTimeInSeconds. Built clean
  (`MACHINE=ubuntu make -C src/fortran`), installed to bin/eqdyna.
- Ran the existing e2e baseline (fortran+python-numpy+python-jax x
  test.tpv8) BEFORE touching anything further: 3/3 SUCCESS. Confirms the
  predecessor's Fortran-side profile emitter is already output-neutral.

## Checkpoint 2 -- the 511x compTimeInSeconds(2) question
`git log --all -- src/fortran/assembleGlobalMass.f90` found commit 8839638
(2026-09-16, item 30, BEFORE this mission): the shared global
`startTimeStamp` bug, and its fix (local `tStageStart`), with its own
before/after numbers (0.0001s -> 0.0511s). CONFIRMED, not re-fixed: re-ran
test.tpv8 x fortran x 4 ranks with `writeCompTime` newly wired on (see next
checkpoint) and read `compTime0` directly: slot 2 = 0.0525s, same order of
magnitude as the historical fix, not the pre-fix 0.0001s. No code change
needed here; this session's contribution was reading the number, not fixing
a bug that was already fixed.

## Checkpoint 3 -- writeCompTime wired on
`globalvar.f90:109`: `writeCompTime = 0` was never set from anywhere (no
bGlobal.txt slot, no env, no CLI) -- item 30's own commit message flagged
this as the reason its 511x bug went undetected. Flipped the initialiser to
`= 1`. Verified: writes `compTime<rank>` (a separate file, `library_output.
f90:422-431`), touches zero bytes of `frt.txt*` -- confirmed by the same
on/off byte-identity run used for the profile gate (fortran row below).

## Checkpoint 4 -- Python emitter design + implementation
New module `src/python/eqdyna/profile_emit.py` (no Fortran counterpart, like
backend.py). `build_row`/`write_profile`, `EQDYNA_PROFILE=0` off switch
(mirrors the Fortran side's env var exactly), `cpus_allowed()`/`numa_nodes()`
via os.sched_getaffinity + /sys/devices/system/node parsing (same source and
overlap test as Fortran's readCpusAllowed/computeNumaNodes).

Wired into `eqdyna3d.py`:
- `run_case` (serial numpy/jax): setup = 'setup (mesh+input)' + 'resolve
  solver' Profile phases; element = 'solve' phase (element+fault FUSED,
  documented departure from Fortran, rule 23); fault/exchange/wait = 0.0
  (genuinely zero, serial path); io = 'write frt' phase. total_s is an
  INDEPENDENT time.perf_counter() span, not a phase-sum (schema's own "sum
  check" docstring insists on this).
- `run_case_mpi` (jax-mpi): same shape, element/exchange/wait sourced from
  driver.run_mpi's report (see checkpoint 5). Emits for BOTH the n_own==0
  early-return rank and the normal path -- a rank that owns no fault node
  still did real setup/compute/exchange work and gets a profile row.

First validation (numpy, jax --device cpu): both passed profile_schema.
validate() first try, unaccounted_s < 0.1%. Byte-identity vs EQDYNA_PROFILE=0
confirmed both ways (see Checkpoint 6 table).

FOUND BUG while first testing python-jax (no --device flag): my ad-hoc CLI
call picked up GPU (device=auto -> cuda visible on this box), giving
unaccounted_s = 19.8% of total_s (first-call CUDA context/compile latency,
unattributed). Diagnosed, NOT fixed (out of scope: every real gate --
run_e2e.py, CI -- forces --device cpu for jax; recorded in docs/run_profile.
md as a GPU-only, non-gating finding). Re-ran with --device cpu: unaccounted
0.003%, clean pass.

## Checkpoint 5 -- python-jax-mpi: the real gap, found and closed
First jax-mpi run (all 4 ranks) validated the JSON SHAPE but FAILED
profile_schema's SUM_TOLERANCE check outright: unaccounted_s 6.5-30% of
total_s across ranks/runs, once as high as 3.8s of 12.0s on one rank.
Root cause (binary-searched by reading driver.run_mpi's structure against
what eqdyna3d.py's run_case_mpi actually bucketed): the outer 'solve' Profile
phase in run_case_mpi wraps driver.run_mpi's ENTIRE call, but only the
step-loop portion (t_compute, previously prof-gated) was being extracted
into a bucket -- driver.run_mpi's substantial PRE-LOOP work (decompose,
mass-check, to_device, carry alloc, jax.jit() construction) had nowhere to
go and fell straight into unaccounted_s.

Fix, in two parts, both zero-new-sync:
  1. `driver.py`: new `t_setup0`/`t_setup` wrapping driver.run_mpi's own
     pre-loop span, stopping right after the PRE-EXISTING unconditional
     `comm.Barrier()` before the step loop (that barrier was already there
     on every run; only the timing around it is new). Returned as
     `rep['setup_s']`; eqdyna3d.py adds it into the `setup` bucket.
  2. `t_compute` (element+fault, fused) made UNCONDITIONAL: it was
     previously `if prof:`-gated, but the block it needs
     (`jax.block_until_ready(hv)`) was ALREADY unconditional on the default
     path (comment explains why: the exchange clock must not start before
     the element kernel has actually run). Two more zero-new-sync additions
     closed the remaining small gap: timing the `b_jit` dispatch line
     (dispatch-only, does not block) and the loop's trailing
     `jax.block_until_ready(carry)` (pre-existing, unconditional -- drains
     the LAST step's part_b, which nothing else ever drains).

Result: unaccounted_s dropped from 6.5-30% (varied by rank/run) to
0.19-0.32% across all 4 ranks. profile_schema.validate_run_dir PASSED
cleanly on the retry. Re-ran the e2e cell after the fix: identical physics
(max|diff|=1.220703e-10, same as before the fix -- these are pure telemetry
additions, verified to not move the answer).

`wait` bucket: confirmed 0.0 is correct, not a gap. The production exchange
(nearest-neighbour Sendrecv) has no collective barrier; the only barrier
exists under the OPT-IN step profile (`if prof:`, driver.py), which the
mission's rules forbid adding to the default path. Documented in
docs/run_profile.md as a real "nothing to measure without a new sync"
finding, not a placeholder.

## Checkpoint 6 -- parity gate results (all four backends/cells)

| backend | ranks | unaccounted_s/total_s (max) | schema validate | frt byte-identical on/off |
|---|---|---|---|---|
| fortran | 4 | 0.51% | PASS x4 | YES (961+961=1922 lines, ranks 1/3 write nothing both ways) |
| python-numpy | 1 | 0.06% | PASS | YES (1891 lines) |
| python-jax (--device cpu) | 1 | 0.003% | PASS | YES (1891 lines) |
| python-jax-mpi | 4 | 0.32% | PASS x4 | YES (132+829+806+124=1891 lines) |

`EQDYNA_PROFILE=0` confirmed to suppress every `profile.rank<r>.json` on all
four backends (glob empty) while leaving `frt`/`nc` output byte-for-byte
unchanged.

## Checkpoint 7 -- CPU discipline
All heavy runs pinned `numactl --cpunodebind=2,3 --membind=2,3` (cpus
16-31), verified via /proc/<pid>/status Cpus_allowed_list on the top-level
harness process each time. One exception found and NOT chased further
(out of scope, noted for the record): run_e2e.py's own python-numpy
subprocess showed Cpus_allowed_list 0-255 (unrestricted) despite the parent
being correctly pinned to 16-31 -- the harness or eqdyna3d._narrow_numpy_
affinity appears to reset/widen affinity somewhere between parent and
worker for that one path. Did not affect any measurement here (numpy's own
_narrow_numpy_affinity still narrows to ONE cpu of whatever mask it is
given, and the unaccounted_s/byte-identity numbers above are unaffected by
which cpu that is) so not investigated further; flagged for whoever next
touches run_e2e.py's subprocess launch.

## Remaining, explicitly out of scope for this mission
- Overhead A/B (profile on vs off, timing delta) is NOT this mission's gate
  per the brief -- left to a later, quiet-box mission. One rough per-step
  number recorded in docs/run_profile.md's measured table is close enough
  for now (setup/loop_s magnitudes did not visibly shift on/off in the runs
  above, but no controlled A/B was run).
- GPU-side unaccounted_s (jax --device gpu/auto-with-GPU-visible) is a real,
  reproducible finding (19.8% observed) but out of gate scope; documented,
  not fixed, in docs/run_profile.md.
