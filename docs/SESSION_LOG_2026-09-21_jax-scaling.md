# Session log — 2026-09-21 — jax 16-core scaling (pathway item 43)

Conductor: wei-lin. Start: 11:01. Base: `51b0fa1` = `v5.13.0`, tree clean,
`origin/master` level (verified, not assumed: `git status --porcelain | wc -l`
= 0; `git log --oneline -1 origin/master` = `51b0fa1`).

ONE objective for 24h: jax backend scales ~linearly to 16 cores (bar 12-14x),
against Fortran's measured 14.39x. Closes item 43. Not a board sweep.

## Budget reading (stated before spending any)

- **Wall clock is the scarce resource, and it is mostly Mira's, not mine.** The
  work is one deep mission (shard_map over a 16-device CPU mesh) already in
  flight. My spend goes to the merge boundary and to keeping the loop alive,
  not to re-deriving her measurements.
- **Cores are the second scarce resource and they are contended.** `nproc` 64,
  but load average at start was 19.45 / 42.03 / 45.54 from a foreign tenant
  (user `kehua`, ~21 single-core jobs — the same contamination items 33/42/43
  all record). Every scaling number this campaign produces must be read against
  a recorded contention sample; a number without one is not evidence.
- **Therefore heavy jobs SERIALIZE.** Item 42: three timeout near-misses across
  two sessions, all from self-inflicted stacking. Mira's measurement passes and
  my gate sweep (`run.py all`, 30 cells, 36-90 min) may never overlap. I hold
  my sweep until she reports a landing candidate.
- **Tag authority as granted, restated for correction before the first tag:**
  patch and minor on `master` of this repo. No major, no publish, no
  force-update of an existing tag. Campaign sign-off is the owner's.

## Agents

- `mira-volkov` `ad64043fb3e42e07a` — ADOPTED, not respawned. Worktree
  `/home/utig5/dliu/wt-jaxshard`, branch `jax-shardmap`, created 10:59
  (2 min before my first probe), base `8919405` = `v5.13.0`.
- **Base staleness assessed, not waved past (gate axis 4):** her branch is 2
  commits behind `51b0fa1`, and both are DOCS ONLY — `ec4d464`
  (`docs/SESSION_LOG_2026-09-19_autopilot.md`, +88) and `51b0fa1`
  (`pathway_forward.md`, +4). Zero source overlap, so her code work is
  unaffected. Consequence recorded here because it is the one way this bites:
  a wholesale copy of `pathway_forward.md` out of her worktree would revert
  those 4 lines. Board rows are zofia-kaminska's (rule 19) and no mission
  writes that file, so the exposure is closed by ownership, not by luck.
- **Channel limitation, stated plainly:** this deployment has no
  `SendMessage`/`ListAgents` tool — only `Read`/`Edit`/`Write`/`Bash`/`Agent`.
  I cannot query her for status and cannot enumerate live agents. Adoption is
  therefore by observation of artifacts she owns (worktree mtimes, git index,
  commits, `NOTES_*`, her processes), which is the liveness definition I use
  anyway. This is the one thing about this run I would ask the owner to fix:
  without a message channel I cannot relay a mid-flight correction to her.

## Gate shape for this objective (decided up front)

A `shard_map` rewrite lands in `src/python/eqdyna/backend.py`, whose 3 scatter
sites are reached from 8 call sites in `assembleGlobalKU.py` and 1 in
`faulting.py`. That is every python cell in the table, both backends. So:

- Merge gate = full `python3 testsys/run.py all` (30 cells), unbounded or
  >=90 min, nothing else heavy on the box. Not a subset.
- Parity first: every cell inside its own `testsys/matrix.py` bound. No
  reference regenerated (rule 7), no bound loosened.
- jax-CPU only. No gate may assume GPU bit-identity (8.9e-08 run-to-run).
- Per-step cost by DIFFERENCE over two step counts, both backends — the bias
  `88f5227` removed (it read 91.25 vs a true 64.74 ms/step at 16 ranks).
- Every verdict carries a property of the compared content (cell count, row
  count, a grep count). "Green" with no such property is a vacuous gate; three
  of those passed in the last campaign.

## Stopping rule (owner's, recorded so it is not softened later)

Stop on a named mechanism WITH a measurement: "plateaus at Nx because <thing>
costs <measured amount> per step." The withdrawn 6x kill criterion does not
apply. The Amdahl 6.7-10x cap in item 43 was computed off AUTO-partition counts
and does not bind explicit sharding — it is not grounds to decline the attempt.

## Timeline

- **11:01** — Orient. Git state verified clean and level. Item 43 read in full.
  Mira confirmed live (worktree 2 min old, git index advancing) — not
  respawned. Contention baseline recorded above.

- **11:03-11:11** — Mira liveness confirmed by artifacts, since I have no
  message channel. 16 polls at 30 s. Dirty-file count went 0 -> 2 -> 3 at
  11:10 (`src/python/eqdyna/backend.py`, `driver.py`,
  `assembleGlobalKU.py`). Her running processes show the shard_map line of
  attack: a `psum`-cost-vs-message-size microbenchmark across
  `--xla_force_host_platform_device_count` 4/8/16 x {base, eigen1, gloo}, and
  an `EQDYNA_JAX_DEVICES=2` run of `test.tpv8`. So she is executing the
  ORIGINAL brief correctly and is ~10 min into it.

- **11:12 — OWNER REDIRECT arrives** (relayed by coordinator): stop treating
  jax's built-in parallelism as the path; do it the explicit MPI way, one
  Python process per rank, reusing the Fortran decomposition (npx/npy/npz plus
  `MPI4NodalQuant`'s boundary lists) that already delivers 14.39x. shard_map
  demoted to a comparison data point.

- **11:12 — BLOCKER, escalated, not worked around: I cannot deliver the
  redirect.** This deployment exposes only `Read`/`Edit`/`Write`/`Bash`/`Agent`
  — no `SendMessage`, no `ListAgents`. I can observe Mira's worktree but cannot
  speak to her. Spawning a second mira-volkov is the one thing the coordinator
  explicitly forbade ("she is yours, do not go around her"), and it would also
  put two agents on `backend.py` — the collision this role exists to prevent.
  So the redirect must be relayed through the same channel the previous two
  arrived on. Recorded as a channel defect, not a Mira defect: she is doing the
  right work against the brief she holds.

- **11:13-11:20 — groundwork on the new primary path, chosen to collide with
  nothing she holds** (no repo file touched, no heavy job co-scheduled). Three
  results, each a number with what produced it:

  1. **Prerequisite was MISSING and is now in place.** `mpi4py` and `mpi4jax`
     were both absent from every interpreter on this box, so the redirect's
     primary path was unrunnable as stated. Open MPI 4.1.1 and `mpif90` do
     exist. Installed `mpi4py 4.1.2` (cp310 manylinux wheel) into the project's
     actual jax interpreter — which is `/home/utig5/dliu/gns/gns/venv_cotopaxi`
     (`python3` resolves there; `sys.base_prefix` is `/usr`). **Environment
     change, declared:** additive only, offline from a cached wheel
     (`--no-deps --no-index --find-links`), so nothing else in that venv moved.
     Reverse with `python3 -m pip uninstall mpi4py`. An isolated
     `--system-site-packages` venv was tried first and rejected: it inherits
     `/usr`, where jax does not live. `mpi4jax 0.9.1.post2` is sdist-only
     (needs a source build against jax's C API, and it version-pins jax) —
     cached at `<scratchpad>/mpitest/`, NOT installed. That is Mira's call.
  2. **Halo cost: 0.116 ms for a 160 kB ring `Sendrecv`, at 16 ranks, all 16
     measuring identically.** Against Fortran's measured 64.74 ms/step at 16
     ranks that is **0.18% of a step**; ten exchanges per step would be 1.8%.
     The exchange cannot be what bounds this path. This is the redirect's
     central bet, and it holds at the smoke level.
  3. **Memory risk retired at this scale: 0.50 GB peak RSS per rank, 16 ranks
     ~8 GB**, against 597 GB available. Even the worst known figure
     (drv.a6 x jax at 9.57 GB) x 16 = ~153 GB fits.

  **And one risk CONFIRMED, which is the finding that matters:** each rank
  spawns **8 threads** — 128 threads on 16 cores — *despite*
  `OMP_NUM_THREADS=1` AND `--xla_cpu_multi_thread_eigen=false`. The redirect
  expected process-per-rank to remove the oversubscription problem; measured,
  it does not. jaxlib holds its own intra-op pool that neither flag reaches.
  This is now the first thing the MPI path must solve, and it is the same
  hazard as check 1 of the original brief, relocated rather than removed.

  Gate honesty on the above, per the vacuous-gate rule: 15 of 16 rank lines
  were captured in stdout, not 16 (rank 11's line was lost in interleaving).
  All 16 ranks provably ran regardless — every printed rank reported
  `size=16`, and a 16-rank ring `Sendrecv` cannot complete with a rank
  missing. Stated rather than rounded to "16/16".

  Verify: `mpirun -np 16 --bind-to core python3 <scratchpad>/rankcheck.py`.

- **11:25-11:40 — RETRACTION of my own 11:20 finding, and the measurement that
  replaces it.** I reported 8 threads/rank as a confirmed oversubscription
  hazard; the coordinator reasonably elevated it to the campaign's central
  finding and told Mira to fix it first. **It does not survive being measured,
  and the instruction should be withdrawn.** I reported a thread COUNT and
  called it contention, which is the exact error CLAUDE.md's "Measure, do not
  infer" section is about — four claims in this repo failed this way and two
  nearly became fixes. This would have been the fifth, and it was mine.

  What the count is worth, measured per rank from `/proc/self/task/*/stat`
  over a 12 s window (`<scratchpad>/threadprobe.py`): threads=8, but peak
  thread states are **`{R: 2, S: 7}`** — never more than 2 runnable — and
  **`EFFECTIVE_CORES` (cpu_used/wall) = 0.96-1.00 on all 16 ranks** with
  `cpus_allowed=1`. Eight threads consume one core's worth. The pool is
  overwhelmingly parked. Under core binding the 8 threads cost nothing.

  Scope of the retraction, stated so it is not over-read: this was measured
  BOUND (`cpus_allowed=1`). Whether those threads become runnable and
  contend when UNBOUND is unmeasured — the `--bind-to none` arm produced no
  parseable output (see the papercut below). Since the recommended
  configuration is bound, the hazard is moot for the chosen path. It is not
  "disproven in general"; it is "costs nothing in the configuration we will
  use."

- **11:40 — my read on `--bind-to core --map-by core`, which the coordinator
  asked for: NECESSARY BUT NOT SUFFICIENT, and blind binding is actively worse
  than none here.** Two measurements.

  1. **What mpirun picks vs what is busy.** `--bind-to core --map-by core
     --report-bindings` selects socket 0, cores 0-15, always. Per-cpu
     utilization over a 4 s window: **0 of 64 cpus are below 5% utilized** —
     there is no idle cpu on this box, only less-loaded ones — 29 of 64 are
     over 50%, 12 of 64 over 90%. Within mpirun's own pick of 0-15: cpu 5 and
     cpu 8 at **1.00**, cpu 11 at 0.81 (an earlier sample had cpus 0, 6, 10,
     14 at 1.00). So binding lands ranks on top of the tenant.
  2. **What that costs, measured, not inferred.** In the 2-rank `--bind-to
     core` run, rank 0 landed on a saturated cpu and read
     **`EFFECTIVE_CORES=0.51`, ms/step 18.31 against rank 1's 9.10 — 2.01x
     slower on identical work** (`elems=3818584` both). In a halo-synchronised
     solver the slowest rank sets the step, so one unlucky rank in sixteen
     halves the whole measurement. That is why blind binding is worse than
     none: it converts a fluctuating average into a pinned straggler.

  **What does work, measured:** `mpirun --cpu-set <16 measured-least-loaded
  cpus> --bind-to cpu-list:ordered` → 16/16 ranks, `cpus_allowed=1` each,
  `EFFECTIVE_CORES` 0.96-1.00 across the board. The 0.51 straggler is gone.
  The cpu set must be re-measured immediately before each run: no cpu here is
  ever idle, the least-loaded set MOVES (item 33 already records the busy set
  shifting mid-session), and a stale pin list silently reintroduces case 2.
  The repo already has this mechanism — `run_scaling.py`'s free-node-aware
  `numactl` pinning — and it should be reused, not reimplemented.

- **11:40 — the finding that actually bounds this campaign, and it is about the
  BOX, not about jax.** With all 16 ranks correctly pinned to least-loaded
  cpus and all 16 at `EFFECTIVE_CORES` ~1.00 — i.e. nobody starved of cpu
  time — per-rank cost still spread **7.94 to 21.25 ms/step, a 2.68x spread**,
  on identical per-rank work. Ranks are not short of CPU; they are short of
  **memory bandwidth**, which is shared with the tenant and cannot be pinned.

  Mean 12.62 ms, max 21.25, so **max/mean = 1.68x lost to straggler alone**.
  For a perfectly parallel 16-rank job with a zero-cost collective, that turns
  an ideal 16x into **~9.5x** before any serial fraction or any collective is
  charged. If that transfers to the solver, **the 12-14x bar is not measurable
  on this box at this tenancy, however good the implementation is.**

  **The caveat that keeps this honest, and it is a big one:** my probe is
  `(a*c+0.5).sum()` over a 30.5 MB array — a pure streaming kernel, maximally
  bandwidth-bound, so 1.68x is an UPPER bound on straggler severity, not the
  solver's figure. Fortran's own measured **14.39x at 16 ranks on this same
  box under load 22-34** (item 43) is direct evidence the real solver is much
  less bandwidth-sensitive than this probe. So: not "the bar is unreachable" —
  rather **"no 16-rank number on this box is interpretable without its
  per-rank spread printed beside it."** That is now a harness requirement, not
  a caveat.

- **Papercut, and a vacuous-gate near-miss of my own:** X11's `Invalid
  MIT-MAGIC-COOKIE-1 key` is emitted with NO trailing newline, so it
  concatenates onto the first line of real stdout. My `grep -E "^rank"`
  therefore matched nothing and I was one step from reporting the probe as
  broken when it had run perfectly. Both arms of the bind-mode comparison were
  lost this way. Use `tr '\r' '\n' | grep -o "rank .*"`, never a `^` anchor,
  on any mpirun stdout from this box.

- **11:45 — the recommended mechanism VERIFIED in code, not assumed.** I told
  the coordinator to reuse `run_scaling.py`'s pinning rather than reimplement
  it; before letting that stand I read it. `free_node_map()`
  (`testsys/perf/run_scaling.py:190`) is already **per-cpu, not whole-node** —
  fixed 2026-09-18 for exactly this box's behaviour, with the reasoning in its
  own docstring ("~20-28 foreign single-core jobs with no cpu affinity of
  their own, so the Linux scheduler smears them across all 8 nodes and no node
  is ever seen fully idle"). It re-probes fresh on every call, and a per-cpu
  `require_idle` re-check runs immediately before the mpirun pin. It also
  already has `--repeats`. So the recommendation stands on read code. This
  matters because per-NODE selection would NOT work here: cpu 44 and cpu 46
  are both on NUMA node 5 and measured 7.94 vs 21.22 ms/step — 2.7x apart
  within one node. Topology (8 nodes x 8 cpus, 2 sockets, no SMT) is not the
  discriminator; the tenant's per-core occupancy is.

- **11:50 — and the constraint that decides what is measurable in the next 24
  h: the harness's own strict ceiling cannot select a single cpu today.** Two
  independent 8 s samples: **0 of 64 cpus under the `--busy-ceiling 0.2`
  default, both times.** At 0.35 -> 11 then 5 cpus. At 0.45 -> 32 then 28. So
  a strict-ceiling 16-rank point is IMPOSSIBLE at today's tenancy: it will
  SKIP, which is exactly what item 33 kept recording ("no 32-core point
  survived the busy check in any of the 4 attempts"). A valid 16-rank number
  today needs the ceiling at ~0.45, and that choice must be printed with the
  number, not buried.

  Set stability, measured because the re-probe design depends on it: between
  consecutive 8 s samples the under-0.45 set went 6 cpus out, 2 in, 26 stable
  — ~80% stable over 8 s but drifting, so `free_node_map`'s re-probe-before-pin
  is necessary. **It is not sufficient, and this is the subtle part:** cpus 46
  and 47 both sat UNDER 0.45 in the pre-check and still produced the 21.2
  ms/step stragglers. A pre-sample does not predict a rank's throughput,
  because the binding constraint is shared memory bandwidth over the run, not
  cpu occupancy at one instant. **Therefore the per-rank spread must be
  printed from the run itself. No pre-check can substitute for it.**

## Standing requirements for any scaling number this campaign produces

Derived from the measurements above, not from preference:

1. **Rank count and per-rank element count printed** beside every number. A
   multi-process harness that silently runs on one rank produces plausible
   figures (coordinator's concern; it is the right one).
2. **Per-rank ms/step printed for every rank, plus max/mean.** A 16-rank mean
   with a 2.68x hidden spread is not a scaling measurement.
3. **`EFFECTIVE_CORES` = cpu_used/wall per rank.** This is what separates "the
   rank was starved of cpu" (0.51, fixable by pinning) from "the rank had a
   full core and was still slow" (1.00, bandwidth, not fixable here).
4. **The busy ceiling actually used, printed.** 0.2 selects nothing today; a
   number taken at 0.45 or under `--i-know-the-box-is-busy` is a different
   claim and must say so.
5. **Per-step cost by DIFFERENCE over two step counts**, both backends
   (`88f5227`'s fix — it read 91.25 vs a true 64.74 ms/step at 16 ranks).
6. **Thread count per rank printed**, per the coordinator — kept even though
   the hazard retracted above, because the retraction is conditional on
   binding and the print is what would catch it regressing.

## 12:05 — PREREQUISITE MET: Fortran re-baselined on today's box. The denominator holds.

Promoted from queued follow-up to prerequisite by the coordinator, on the
reasoning that my bandwidth probe implied a ~9.5x straggler ceiling while
Fortran had measured 14.39x on this same box — both cannot describe one
machine unless the real solver is far less bandwidth-bound than a 30.5 MB
streaming sum. That was my own stated caveat, and this settles it.

**Fortran source is bit-identical between `737358f` (`v5.12.0`, the 14.39x
baseline) and HEAD** — `git diff --stat 737358f..HEAD -- src/fortran/` is
empty — and `bin/eqdyna` is still the 2026-09-19 15:00 build. So this
measurement changes exactly one variable: today's tenancy.

Run from the clean main checkout (`git status --porcelain | wc -l` = 0),
per-step by difference at n_lo/n_hi 40/160, compact placement, ceiling 0.45:

| np | ms/step today | speedup today | 09-19 (item 43) | delta |
|---|---|---|---|---|
| 1 | 937.41 | 1.00x | 931.65 | +0.6% |
| 2 | 467.08 | 2.01x | 466.64 | +0.1% |
| 4 | SKIPPED (cpu 1 at 87%) | — | 237.62 | — |
| 8 | SKIPPED (cpu 1 at 87%) | — | 117.07 | — |
| 16 | **65.89** | **14.23x** | 64.74 / 14.39x | **+1.8% / -1.1%** |

**VERDICT: Fortran reads ~14x today, not ~9x.** Per the coordinator's own
decision tree, that means: the ~9.5x straggler ceiling **does not bind the
real solver**, near-linear remains a live target, the bar stays at 12-14x,
and item 43's 14.39x denominator is confirmed rather than replaced. My
probe's upper-bound caveat was the correct reading of it.

**The tension is now resolved with numbers on both sides, not papered over.**
A pure-streaming kernel loses 1.68x to stragglers on this box; the FEM solver
loses ~1% — and it did so while spread across **three NUMA nodes** (np=16
took cpus 0,1,3,4,5,6,8,10,11,12,13,15,17,18,19,20 = nodes 0,1,2) at whole-box
load 26-30 with 15 other tenants named in the tool's own roster. Element-local
FEM work has cache reuse a 30.5 MB streaming sum cannot; the solver is
therefore not bandwidth-limited at 16 ranks, measured, not assumed.

**Gate honesty:** 2 of 5 points were SKIPPED, correctly — the tool refused
`np=4`/`np=8` because cpu 1 read **87%** against the 45% ceiling, printing
"REFUSING TO MEASURE ... Taken while one of them is busy the number measures
that process, not us." That is the gate working, and it is the third
independent confirmation that this box cannot deliver a full strict-ceiling
sweep. The three surviving points are the ones that matter: the 1-core anchor
and the 16-core target, both reproducing to within 1.8%. `--repeats 2` was
requested; the tool emitted one line per config.

**A reframing the owner should rule on (not a blocker, no action taken).**
The bar "12-14x at 16 cores" is a *ratio* target, and because jax's 1-core
kernel is already faster than Fortran's (611.00 vs 937.41 ms/step), the two
readings of the goal are not the same thing:

- **Ratio parity** (12-14x) requires jax at ~43.6 ms/step — which would be
  **1.5x FASTER than Fortran's absolute 16-core time** of 65.89.
- **Absolute parity** with Fortran at 16 cores needs only 65.89 ms/step, i.e.
  **9.3x** — below the stated bar, and arguably what the science cares about.

So 9.3x is the point where jax stops being slower than the production solver,
and 12-14x is the point where it scales as well. The owner's wording ("scale
largely linearly, like the Fortran/MPI solver does") asks for the ratio, so
the bar is unchanged and stands at 12-14x; recording the 9.3x crossover
because it is a real milestone that will be passed on the way, and a run that
reaches 9.3x and stops has still delivered a jax backend that beats Fortran
at 16 cores.

Reproduce: `python3 testsys/perf/run_scaling.py --case test.tpv104
--fortran-ranks 1,2,4,8,16 --py-threads '' --policies compact --n-lo 40
--n-hi 160 --busy-ceiling 0.45 --repeats 2`

## 12:15 — GATE GAP: the 30-cell sweep cannot exercise an MPI-parallel python backend

Found by reading `testsys/e2e/run_e2e.py`, not at merge time. The two backend
paths are not symmetric:

- **Fortran, line 170:** `[MPIRUN, '-np', str(matrix.FORTRAN_RANKS[case_name]),
  eqdyna_cmd]` — MPI-aware, rank count per case from `matrix.FORTRAN_RANKS`.
- **Python, line 147:** `subprocess.call([sys.executable, '-u', '-m', 'eqdyna',
  ...])` — a **single process**. No `mpirun`, no rank count, no equivalent of
  `FORTRAN_RANKS`.

So if Mira's work makes jax run under `mpirun -np 16`, `python3 testsys/run.py
all` will still run all 20 python cells single-process. That means **the full
30-cell sweep, on its own, is exactly the "passes by falling through to the OLD
code path" gate my own axis 2 forbids.** It would come back 30/30 green while
testing none of the new code.

Consequences, both directions, because the news is not all bad:

- **The sweep remains VALID and REQUIRED as a no-regression gate.** It exercises
  the single-process path, which must stay bit-identical. That is a real and
  necessary check — it is just not sufficient.
- **A second gate is needed and does not exist:** at least one cell that
  actually runs a python backend under `mpirun -np N`, N>1, compared against
  **the same committed `frt.canonical.txt`** at that case's existing bound.
- **The reference design already supports this, so no reference moves and rule 7
  is not in play.** `frt.canonical.txt` is deduped on rounded (x,y,z) then
  lexsorted precisely so "a result is a statement about the PHYSICS and not
  about the decomposition" — which is what lets a 4-rank Fortran run, a serial
  Fortran run and a serial Python run all compare to one artifact. A 16-rank
  Python run is the same class of thing. The missing piece is purely a harness
  hook, not a new ground truth.

**Two landing preconditions for Mira, to be relayed:**

1. **Rank-1 must fall back to today's exact code path**, so the 30-cell sweep
   stays meaningful as the no-regression gate. If `-np 1` takes a new branch,
   the sweep stops being a control.
2. **Her landing must come with the harness hook** (an `EQDYNA_PY_RANKS` /
   `--python-ranks` path in `run_e2e.py`, mirroring `FORTRAN_RANKS`) **plus one
   gated cell that uses it**, or it cannot be merged under axis 2. I am not
   asking her to design the test pyramid — but the invocation contract is hers,
   and nobody can write the hook until that contract exists.

`iris-vermeulen` owns the gate cell itself and will be dispatched once Mira's
invocation contract is fixed. Dispatching her now would have her guess at an
interface that does not exist yet, which is how a gate ends up testing the
harness instead of the solver.

**Cost disclosed:** Mira was mid-pass (short single-process jax probes, her own
worktree, `cwd` verified). The coordinator authorised pausing her for this and
I have no channel to deliver the pause, so her probes ran concurrently with my
sweep. My numbers are protected — `run_scaling.py` busy-checks per cpu and
routed around her — but any jax figure she took between 11:24 and 11:35 is
contaminated and should be discarded, not averaged. Stating it rather than
hoping the pinning covered both sides.

