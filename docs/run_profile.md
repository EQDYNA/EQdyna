# run_profile.md -- the always-on per-rank runtime profile

Schema `eqdyna-profile/1` (`testsys/profile_schema.py`, owned by the parallel
profile-guard mission -- read, not edited, here). One file per rank:
`<run_dir>/profile.rank<r>.json`. Six buckets, always present:
`setup`, `element`, `fault`, `exchange`, `wait`, `io`, plus `loop_s`,
`total_s`, `unaccounted_s`, `cpus_allowed`, `numa_nodes`, `sampling_every`.

Default ON, every backend, every run. `EQDYNA_PROFILE=0` turns it off; that
switch exists ONLY for the overhead A/B a later, quiet-box mission runs --
it is not a normal operating mode.

**Rule 23 (PROJECT_RULES.md): Fortran is the reference and DEFINES the
buckets.** The port departs from it in exactly one place, and (as of
2026-09-23) only on the jax backends: `fault` folded into `element` on
python-jax and python-jax-mpi. python-numpy splits `fault` out along
Fortran's own boundary and matches it. Measured and documented below, not
guessed at.

## Fortran -- `src/fortran/`

Emitter: `output_profile` (#10), `src/fortran/library_output.f90:408-600`.
Off switch: `EQDYNA_PROFILE=0`, checked at `library_output.f90:428`. Called
once per rank at end of run, `src/fortran/eqdyna3d.f90:112-115`.

| bucket | source | file:line |
|---|---|---|
| `setup` | `compTimeInSeconds(1)+(2)` | `eqdyna3d.f90:72,76` (mesh/input read through `assembleGlobalMass`) |
| `element` | `compTimeInSeconds(3)+(4)+(5)` | `driver.f90:211` (velDispUpdate), `assembleGlobalKU.f90:70`, `calcHourglassResist.f90:99` |
| `fault` | `compTimeInSeconds(6)` | `faulting.f90:28` |
| `exchange` | `mpiCommLoop - mpiWaitLoop` | `eqdyna3d.f90:84-96` snapshots `MPICommTimeInSeconds`/`MPIWaitTimeInSeconds` at the setup/loop boundary; the loop-phase deltas are disjointed by subtracting the nested wait sub-timer |
| `wait` | `mpiWaitLoop` | `assembleGlobalMass.f90:79,216-218,221` (`tWaitStart`/`MPIWaitTimeInSeconds`, a sub-timer NESTED inside `MPI4NodalQuant`'s own `MPICommTimeInSeconds` span, subtracted back out for `exchange` above) |
| `io` | `compTimeInSeconds(8)` | `eqdyna3d.f90:98-105` (`output_onfault_st`/`output_offfault_st`/`output_frt`/plastic/finalSurfDisp) |
| `loop_s` | `loopS` | `eqdyna3d.f90:87-89` (wraps the `driver` call only) |
| `total_s` | `compTimeInSeconds(9)` | `eqdyna3d.f90:106` (`simuStartTime` to end) -- an INDEPENDENT timer, not a bucket sum |

### The 511x `compTimeInSeconds(2)` bug: CONFIRMED FIXED, not by this
mission -- by commit `8839638` (2026-09-16, item 30), before this mission
started. `globalvar.f90`'s single shared `startTimeStamp` let
`MPI4NodalQuant` (called twice inside `assembleGlobalMass`) silently reset
the outer timer `eqdyna3d.f90:74` had started, so slot 2 measured only the
tail of the LAST `MPI4NodalQuant` call. The fix replaced every nested
writer's use of the shared global with a LOCAL `tStageStart`
(`assembleGlobalMass.f90:72`, comment at `:68-71` explains why). Commit
8839638's own measured numbers: `compTimeInSeconds(2)` 0.0001 s before,
0.0511 s after (511x). Re-confirmed live in this mission on a fresh
test.tpv8 x 4-rank run (rank 0, `writeCompTime` newly wired on --
`compTime0`): slot 2 = **0.0525 s**, the same order of magnitude as the
2026-09-16 fix, not the 0.0001 s pre-fix value. No further fix was needed;
this mission only re-ran the measurement.

### `writeCompTime` wired on
`globalvar.f90:109` initialised `writeCompTime = 0` and NOTHING ever read or
set it away from 0 (no `bGlobal.txt` slot, no env, no CLI) -- the legacy
`compTime<rank>` dump (`output_timeanalysis`, `library_output.f90:422-431`,
gated by `eqdyna3d.f90:117`) was unreachable without editing source and
rebuilding, which is why the 511x bug survived undetected in it. Flipped the
initialiser to `writeCompTime = 1` (`globalvar.f90:109`, this mission,
2026-09-23). Verified output-neutral: `compTime<rank>` is its own file and
touches no `frt.txt*` byte; the profile on/off byte-identity gate below does
not even exercise this flag (it is independent of `EQDYNA_PROFILE`).

## Python -- `src/python/eqdyna/`

New module, no Fortran counterpart (like `backend.py`): `profile_emit.py`.
Builds and writes the row; `enabled()` is the `EQDYNA_PROFILE=0` off switch
(same env var, same contract, both languages). Wired into `eqdyna3d.py`'s
`run_case` (serial numpy/jax) and `run_case_mpi` (jax-mpi).

### BUCKET DEPARTURE FROM FORTRAN (rule 23, documented) -- jax only, since 2026-09-23
Every Python backend fuses velDispUpdate + both element kernels + faulting
into ONE step function (`driver.py`'s module docstring: `make_step` for the
serial path, `a_body`/`b_body` for MPI). Splitting `element` from `fault`
back apart the way Fortran does needs a NEW timer INSIDE that fused step.

For **jax** that timer would be a NEW `block_until_ready` between the
element kernels and faulting -- forbidden by the "no new sync" rule on the
default path, and a timer placed anyway would measure the ONE-TIME trace,
not the per-step cost (wrong, not imprecise). So **`fault` stays 0.0,
folded into `element`, on python-jax and python-jax-mpi**. Read `element` as
"element + fault, fused" for those two backends only, not as a
like-for-like column against Fortran's `element`.

For **numpy** that reasoning does not hold: numpy executes the step body
eagerly and synchronously (no queue, no device to drain), so a
`time.perf_counter()` pair around the `FLT.faulting` call inside
`driver.make_step_parts`'s `part_b` costs no new sync and changes no
arithmetic -- landed 2026-09-23, guarded by `fault_timer is not None`
(created only `if not B.is_jax(xp)`) so the jax path builds neither the
timer nor the branch's timed variant of the call. **python-numpy now
matches Fortran's bucket boundary exactly.**

### python-numpy (serial, nranks=1) -- matches Fortran's `fault` boundary
`eqdyna3d.py`'s `run_case`, `src/python/eqdyna/eqdyna3d.py:519-580` (grep
`_profile_emit.write_profile` in that function for the exact call site).
`driver.run` returns `fault_s` (the accumulated `FLT.faulting` time across
all steps) alongside its usual dict; `run_case` splits it back out of
`Profile['solve']`.

| bucket | source |
|---|---|
| `setup` | `Profile['setup (mesh+input)'] + Profile['resolve solver']` |
| `element` | `Profile['solve'] - out['fault_s']` |
| `fault` | `out['fault_s']` -- `driver.f90:28`'s boundary, measured via `perf_counter()` around `FLT.faulting` in `make_step_parts`'s `part_b` |
| `exchange` | `0.0` -- genuinely zero: this path is serial by construction (`build_solver_state` refuses `npx/npy/npz>1`), not a folded/unmeasured cost |
| `wait` | `0.0` -- same reason as `exchange` |
| `io` | `Profile['write frt'] + Profile['write stations']` -- Fortran's `compTimeInSeconds(8)` spans the whole output stage, stations included (`eqdyna3d.f90:175-182`). The phase-to-bucket map is `eqdyna3d._SERIAL_PHASE_BUCKET`; `serial_buckets` refuses a phase it does not map |
| `total_s` | independent `time.perf_counter()` span wrapping the whole of `run_case`, NOT a sum of the `Profile` phases. Each `Profile.phase` starts its clock BEFORE its entry device sync, so the first phase (`setup`) carries the one-time jax import |

### python-jax (serial, nranks=1) -- `fault` folded into `element`
Same `run_case` call site as numpy above; same table except `element` is
`Profile['solve']` unsplit and `fault` is `0.0`, for the jax reason given
above (`driver.run`'s `fault_timer` is `None` on this backend, so
`out['fault_s']` is always `0.0` and the subtraction is a no-op).

`Profile` (the existing `--profile` wall-clock helper, `eqdyna3d.py`'s
`Profile` class) is now read unconditionally, not only under `--profile`; it
was already synchronising jax devices at each phase boundary (`_sync()`), so
reading it always-on adds no new device sync of its own.

### python-jax-mpi (`driver.run_mpi`)
`eqdyna3d.py`'s `run_case_mpi`, `src/python/eqdyna/eqdyna3d.py:440-517`.

| bucket | source | file:line |
|---|---|---|
| `setup` | `Profile['setup (mesh+input)'] + rep['setup_s']` | `driver.py:357` (`t_setup0`) to `:530` (`t_setup`), wraps `decompose`/mass-check/`to_device`/carry-alloc/`jax.jit()` CONSTRUCTION -- all synchronous host-side calls, no sync added |
| `element` | `rep['compute_s']` (`t_compute`) | `driver.py:551-556` (unconditional since this landing; was `if prof:`-gated), `:606-610` (b_jit dispatch), `:614-618` (final drain) -- see below |
| `fault` | `0.0` (folded into `element`, same reason as serial) | |
| `exchange` | `rep['mpi_s']` (`t_mpi`) | `driver.py:564/599` -- was ALREADY unconditional before this landing; unchanged |
| `wait` | `0.0` -- see "what could not be measured" below | |
| `io` | `Profile['write frt']`, or `0.0` for a rank owning zero fault nodes | `eqdyna3d.py`'s `_emit` helper inside `run_case_mpi` |
| `total_s` | independent `time.perf_counter()` span wrapping the whole of `run_case_mpi` | |

**`t_compute` made unconditional (this landing).** It was previously
accumulated only under the opt-in `MQ.step_profile()` flag
(`EQDYNA_MPI_STEP_PROFILE` or similar -- see `MPI4NodalQuant.step_profile`).
The mission's "no new sync" rule is satisfied because the block it times
(`jax.block_until_ready(hv)`, `driver.py:555`) was ALREADY unconditional on
the production path before this change -- it exists so the MPI exchange
clock does not start before the element kernel has actually finished
(comment immediately above it). Timing around an already-mandatory sync
costs two `perf_counter()` calls, not a new synchronisation point. Two more
additions of the same kind, both needed to close a real accounting gap
(below): the `b_jit` dispatch line (`driver.py:600`, now `:609`) is timed
even though it does not block (dispatch-only, still zero new sync), and the
loop's trailing `jax.block_until_ready(carry)` (`driver.py:601`, now
`:616-618`) -- pre-existing and unconditional, only the timing around it is
new -- is folded into `element` because it drains the LAST step's `part_b`,
which no other measurement point ever drains.

**What could not be measured without a new sync: `wait`.** The production
exchange (`MPI4NodalQuant.exchange`, nearest-neighbour `Sendrecv`) has no
collective barrier -- only the OPT-IN step profile adds one
(`driver.py:582-585`, `comm.Barrier()` gated `if prof:`). There is nothing to
time without ADDING that barrier, which the mission explicitly forbids on
the default path. This is reported as `wait=0.0`, and it is a REAL
statement (no barrier exists to wait at), not a placeholder: any per-step
blocking on a slow neighbour happens inside `Sendrecv` itself and is already
inside `exchange`, mirroring how Fortran's own `wait` is a sub-timer nested
inside (and subtracted back out of) `exchange` -- the same nesting
relationship, minus the sub-timer, because there is no separate call here to
subtract one out of.

## Measured on test.tpv8 (canonical case), 2026-09-23, this box, NUMA nodes 2-3

| backend | ranks | max `|unaccounted_s / total_s|` | schema validate |
|---|---|---|---|
| fortran | 4 | 0.507% (rank 1) | PASS (all 4 ranks) |
| python-numpy | 1 | 0.063% | PASS |
| python-jax (`--device cpu`) | 1 | 0.003% | PASS |
| python-jax-mpi | 4 | 0.32% (rank 0) | PASS (all 4 ranks) |

All four comfortably inside `profile_schema.SUM_TOLERANCE` (5%). The
python-jax-mpi number needed the `setup_s`/`t_compute` fixes above: before
them, unconditional `t_compute` did not yet exist and the pre-loop
`decompose`/`to_device`/`jit`-construction cost had nowhere to go, measuring
6.5-30% unaccounted (varied by run) -- diagnosed and fixed in this same
session, not carried forward as a known gap.

**A GPU-side finding, out of scope for the gate but worth recording:**
`python-jax` run with `--device gpu` (the only way to reach a GPU since
board row 56 removed `--device auto`) showed `unaccounted_s` as high as 19.8% of `total_s` -- first-call
CUDA context/compile latency not attributed to any bucket. Not investigated
further because every gate that matters here (`run_e2e.py`, CI) forces
`--device cpu` for jax, per `CLAUDE.md`'s "jax-GPU is nondeterministic run to
run" rule; a GPU-profiled run is not a supported cell of this table.

## Parity gate: profile on vs profile off, byte-identical `frt` output

Verified this session, `test.tpv8`, this box:

| backend | ranks | frt files compared | lines | byte-identical |
|---|---|---|---|---|
| fortran | 4 | `frt.txt0`, `frt.txt2` (ranks 1/3 own 0 fault nodes, write nothing, both settings) | 961+961=1922 | YES |
| python-numpy | 1 | `frt.txt0` | 1891 | YES |
| python-jax (`--device cpu`) | 1 | `frt.txt0` | 1891 | YES |
| python-jax-mpi | 4 | `frt.txt0..3` | 132+829+806+124=1891 | YES |

`EQDYNA_PROFILE=0` was also confirmed to suppress every `profile.rank<r>.json`
file on all four backends (glob returns nothing) while leaving `frt`/`nc`
output untouched -- profile adds files, changes nothing else.
