# NOTES_jax_speed.md -- python-jax critical-path speed investigation, 2026-09-24

Base SHA: `0e68d86` (origin/master, "session log UU: TPV30 gated (#5);
always-on profile landed (#6)"). No `src/`, `testsys/`, or `.github/` file was
touched by this session -- verdict below is "no safe lever found today", not
a landed change, so `git diff 0e68d86 --stat` for every path except this file
is empty.

## Verdict

No code or config change landed. Every candidate lever named in the brief was
measured on this box, uncontended, pinned via `numactl --physcpubind`, on the
three target cells (`test.tpv36`, `test.tpv37`, `test.drv.a6`, all
`python-jax`, `--device cpu`, all built via `testsys/e2e/run_e2e.make_serial_case`
so the case-build path matches the sweep's own exactly). Two of three came
back already-optimal; the third (thread count vs contention) came back a
falsified hypothesis, not a fix -- the same shape CLAUDE.md's "Measure, do
not infer" section warns about (item 60: a plausible mechanism that does not
move the outcome curve). Nothing here is landed as a wire-in because nothing
here is a win.

## Method

`jax.lax.fori_loop` in `backend.time_loop` traces the trip count as a runtime
argument, not a compile-time constant (docstring: "one compiled executable
serves every step count"), so the compile/setup-vs-per-step split from
CLAUDE.md's difference method does not need two FULL-length runs -- two SHORT
runs (nsteps=20, 40) at the same case hit the identical compiled executable
and isolate per-step cost cleanly and far more cheaply than running the full
gate step count twice. Verified this holds (see attribution table: tpv36 and
tpv37 share the same mesh/`E`~0.95M and land on the same ms/step, as expected
of a difference method that only depends on the compiled body).

Every run: `numactl --physcpubind=<cpus> --membind=4,5 python3 -m eqdyna
<case_dir> <nsteps> --backend jax` (`JAX_PLATFORMS=cpu`,
`PYTHONPATH=src/python`), from a scratch case dir built once per case, never
inside a completed/golden case tree. Never two heavy runs at once except the
one deliberate contention A/B below, which is the whole point of that A/B.

## Attribution table (16 cpus, nodes 4-5, uncontended -- this box's "clean" number)

Fixed cost = `setup` bucket (mesh+input+resolve, `profile.rank0.json`) + the
compiled body's one-time trace/warm-cache-load cost (isolated by differencing
element-bucket time at nsteps=20 vs 40) + `io`. Per-step = the same
difference, divided by 20.

| cell | E (elements) | setup_s | compiled-body fixed_s | io_s | ms/step | full-run nsteps | clean estimate | **in-sweep measured (954434c)** |
|---|---|---|---|---|---|---|---|---|
| test.tpv36 x python-jax | ~0.95M | 5.93 | 3.11 | 0.17 | 302 | 464 | 149.3 s | **286.0 s** |
| test.tpv37 x python-jax | ~0.95M | 5.93 (assumed = tpv36, same mesh) | ~3.0 (assumed) | 0.17 | 278.6 | 464 | ~139 s | **276.7 s** |
| test.drv.a6 x python-jax | ~2M | not separately isolated (see caveat) | included below | -- | 733 | 120 | ~92 s + fixed | **217.5 s** |

tpv36 buckets are the median of 3 repeated (n=20,n=40) samples (element_s:
9.15/15.19 median; raw samples 9.208/15.176, 9.121/15.214, 9.209/15.214 --
run-to-run spread ~1%, well inside noise). tpv37 and drv.a6 got one (n=20,
n=40) sample each (raw `driver.run` stdout, not the json bucket, for speed);
their `ms/step` numbers are from that one pair and are the ones used for the
per-step column, so drv.a6's "clean estimate" above is per-step-only
(120*0.733=88.0s) without a separately-measured fixed cost for that case --
call it ~92-96 s all-in, not a mixed-precision claim.

**The gap between "clean estimate" and "in-sweep measured" is 92-137 s per
cell, 1.6-2.3x.** That gap is the real finding, and it does NOT close on any
lever tested below.

## Levers tried

### 1. XLA persistent compilation cache -- ALREADY LANDED, already warm, not a lever today
`backend.enable_compilation_cache()` is called unconditionally from
`backend.run_time_loop` (serial jax path, no `subdir`, cache at
`~/.cache/eqdyna-jax`) and already holds `jit_body-*`/`jit__lambda_-*` entries
(242 entries total in the shared dir on this box, mixed with the MPI path's
per-rank `jit_a_body`/`jit_b_body` subdirectories). The fixed compiled-body
cost measured above (2.2-3.1 s) is in the same range backend.py's own
docstring gives for a WARM run (2.8-6.5 s CPU), not a cold compile (4.7-14 s)
-- so the cache is doing its job. Nothing to change.

### 2. Thread count / intra-op parallelism vs the sweep's core budget -- MEASURED, HYPOTHESIS FALSIFIED
The brief's hypothesis: `run_e2e` bills a python-jax cell
`ceil(matrix.JAX_MEASURED_CORES)` = 3 cores in its scheduling ledger
(`testsys/e2e/run_e2e.py:275`), but `run_standalone` launches the cell with
**no** cpu-affinity restriction at all (unlike `python-numpy`, which
self-narrows via `eqdyna3d._narrow_numpy_affinity`) -- so a jax cell inherits
whatever affinity mask its parent had, sees all 16 cpus of nodes 4-5 on this
box, and XLA's CPU thread pool auto-sizes to that. Several such cells running
concurrently would then oversubscribe by ~5x (each spawning ~16 threads for a
16-cpu pool shared 3-5 ways), which could plausibly cost far more than the
scheduler's 3-core bill assumes.

Solo scaling (tpv36, uncontended, `--physcpubind`):

| cpus pinned | ms/step |
|---|---|
| 1 | 798 |
| 3 | 502-523 |
| 16 (no restriction) | 279-336 |

So jax genuinely uses the extra cores (798->502->~300, real speedup, not
flat) -- restricting a cell to its billed 3 cores would by itself make that
cell 1.5-2.7x SLOWER solo, so pinning-to-billed-cores is not a free win on
its own.

The actual question is contention, not solo scaling. A/B, 3 processes
launched simultaneously on the SAME 16-cpu pool (`tpv36`, nsteps=30):

| mode | batch wall (3 procs) | per-proc ms/step |
|---|---|---|
| shared (all 3 unrestricted, inherit all 16 cpus, ~48 threads on 16 cpus) | 30.63 s | 611-633 |
| pinned (3 disjoint slices, 5/5/6 cpus, no oversubscription) | 29.63 s | 538-641 |

3% batch-wall difference, per-proc spread of both modes overlapping (pinned's
own 538-641 spread is wider than the shared-vs-pinned gap). **Thread
oversubscription inside a fixed 16-cpu pool is not the cost** -- both regimes
converge to ~16/3 cores' worth of throughput per process, whether the OS
scheduler is timesharing unrestricted threads or the cpuset is drawn ahead of
time. Landing an affinity-pinning wrapper in `testsys/e2e/run_e2e.py` (the
`run_standalone` call site) would cost real engineering and buy nothing
measurable -- dropped, per rule 4a (falsify against the outcome, not the
mechanism).

### 3. Where the 92-137 s gap actually is -- not closed, flagged not fixed
Clean-estimate-vs-in-sweep gap does not match anything reproduced by a 3-way
same-pool A/B. The brief's own numbers say these three cells were measured
"in-sweep, box load 30-57" -- i.e. system load average 30-57 at measurement
time, which is far more concurrent demand than 3 same-size jax cells sharing
16 cpus (my A/B above). Two explanations fit and neither is a jax-code or
jax-config defect:
  - the everyday sweep's own `default_jobs_budget` tenancy heuristic may be
    packing more concurrent cells (jax + numpy + fortran-MPI ranks together)
    onto the box than its `cell_cost` billing accounts for, OR
  - load 30-57 includes processes outside this sweep entirely (another
    session's job on a shared box), which no change inside EQdyna can fix.
Distinguishing these needs a repeat of the ACTUAL 23-cell sweep with a
per-process cpu-time/wall-time trace at the moment tpv36/37/drv.a6 run, which
is a scheduling-tenancy investigation, not a jax-speed investigation -- out
of this session's scope per the brief's own instruction not to touch
`testsys/` scheduling without a measured reason, and I have none yet. Not
guessed at further.

### 4. Levers not tried (already known-optimal from static reading, not re-verified here)
  - index dtype: `backend.to_device` already narrows every non-float array to
    int32 (`_PROMOTED_FLOAT` keeps the float loop-invariants as jit
    arguments, not HLO constants -- both already landed, per the module's own
    docstring and comments, not new work).
  - scatter-add lowering: `backend.scatter_add` uses `.at[idx].add(val)`,
    required for bit-parity with numpy's `np.add.at` block order per its own
    docstring -- not a place to touch without breaking the bit-identity
    contract this whole port exists to hold.
  - host<->device transfers: `run_time_loop` wraps the ENTIRE step loop in
    one `jax.jit(body)` call via `lax.fori_loop`; there is no per-step Python
    dispatch or per-step host round-trip to remove -- the "per step" cost
    measured above is 100% inside one compiled XLA executable already.

## Bit-identity

No code changed, so there is nothing to re-check for parity. `git status
--short` at the end of this session (this file only, plus itself) confirms
it; the scratch case dirs and profile jsons this investigation built live
under the session scratchpad, not under `tests/data/` or any golden dir.

## Expected everyday-sweep effect

**None measured, because nothing landed.** The sweep's wall stays governed by
tpv36's ~286 s cell (or whichever of the three is slowest that run) exactly
as at `954434c`. If the owner wants the 92-137 s in-sweep-vs-clean gap
closed, the next investigation is a tenancy trace of the ACTUAL sweep (not a
synthetic 3-way A/B), scoped separately from this jax-speed mission.
