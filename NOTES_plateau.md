# What bounds `python-jax-mpi` at ~5.3x on 32 ranks

Base `12c5be0`. Case `test.tpv104`. All runs `EQDYNA_JAX_CACHE_DIR=off`
(the shared persistent XLA cache wedges multi-rank runs; session log
2026-09-22). Every number here comes from a committed artifact under
`scratch/plateau/` or `docs/perf_snapshots/`; nothing is transcribed from a
console I did not keep.

## Step 0

`git fetch origin && git reset --hard origin/master` -> `git log --oneline -1`
= `12c5be0`. Confirmed before anything else ran.

## The framing arithmetic, re-derived rather than accepted

Fit `T(N) = T1/N + C` on the committed 1->32 table
(`mpi_scaling_2026-09-22_012512_tpv104_1to32_cacheoff.json`), `T1 = 756.92`:

| N | T measured | T1/N | C |
|---|---|---|---|
| 2 | 444.87 | 378.46 | **66.41** |
| 4 | 315.31 | 189.23 | **126.08** |
| 8 | 211.29 | 94.62 | **116.67** |
| 16 | 150.04 | 47.31 | **102.73** |
| 32 | 143.27 | 23.65 | **119.62** |
| 32 (w=0.6) | 132.20 | 23.65 | **108.55** |

Reproduces the brief's values to the decimal. So the target is a ~103-126 ms
per-step cost that DOES NOT SHRINK with rank count, against a scalable part of
~24 ms at 32 ranks.

## HYPOTHESIS H1 (pre-measurement, recorded so it can fail honestly)

**The rank-independent term is replicated global-sized nodal work, not
orchestration and not transport.**

Textual basis, `driver.py:382-389` inside `run_mpi`:

```python
N = S['N']; NEQ = S['NEQ']; nftnd_l = int(finv_l['nftnd'])
carry = (z(NEQ + 1), z((N, 3)), z((N, 3)), z(NEQ + 1),
         xp.asarray(inv_l['stress_i0']).copy(), z((inv_l['Ep'], 15)), ...)
```

`v1`, `velArr`, `dispArr` and `force` take their shapes from the **global**
`S['N']` / `S['NEQ']`, on every rank, at every rank count. Only `stress_i`
(`Ei`) and `s_p` (`Ep`) shrink. `mass` is global too (`driver.py:367`, never
restricted). So at least these are per-step work every rank does regardless of
N:

- `driver.py:137` — `force = setat(force, slice(None), 0.0)`: zero all NEQ+1.
- `driver.py:177` — `force[1:] / mass[1:]`: a global-length divide by a
  global-length mass.
- the a_jit/b_jit boundary — the full carry is materialised per step, twice,
  where the serial path's fused `fori_loop` never materialises it at all.

`MPI4NodalQuant.py:36-52` states this trade explicitly as a cost ("each rank
holding the full-length nodal arrays", "NEQ+1 doubles = 30.5 MB on
test.tpv104, times nranks"), and `backend.py:330-337` already names exactly
this mechanism as the **measured ceiling of the shard_map route**: "the NODAL
stages (velDispUpdate, faulting, the mass divide) stay REPLICATED ... It also
caps the speedup at Amdahl's law on the nodal fraction." H1 is that the MPI
route inherited the same ceiling, partially: `decompose()` DOES restrict the
node-index arrays (`MPI4NodalQuant.py:183-190`), so velDispUpdate's gathers
shrink, but the ARRAYS they write into, the force zero and the mass divide do
not.

**H1 is a plausible mechanism and therefore worth nothing until measured.**
Five such claims have already collapsed in this campaign. What follows is the
measurement.

## Experiment 1 — split the step (in flight)

`testsys/perf/probe_mpi_step_split.py`, report-only, no MPI at all.
`MPI4NodalQuant.decompose` needs no communication (every rank holds the full
mesh), so ONE process can build rank r of N's exact local `inv`/`finv`/`carry`
and run the exact same `a_jit`/`b_jit`. Configurations measured, one artifact
each (rule 20a): `0:1`, `4:8`, `31:32`, `0:32`.

Discriminator, stated BEFORE the numbers arrive:

- per-rank compute at 31:32 ~= 140 ms -> C is INSIDE the computation, MPI is
  not implicated at all, and the next question is which region.
- per-rank compute at 31:32 ~= 24 ms -> C is in the orchestration/host layer
  despite exchange measuring 1-11 ms.

Sub-block attribution (`head` = velDispUpdate + force zero, `elem` =
assembleGlobalKU + calcHourglassResist, `div` = the mass divide) is timed at
the same shapes so a rank-independent cost lands on a code REGION. Those
sub-jits materialise state the real step keeps fused, so their sum is an upper
bound on the step, never an identity.

Donation is checked by asking the BUFFERS (`Array.is_deleted()` on the carry
entries the call consumed), not by hoping a warning was emitted. Bytes/step
come from `nbytes` per carry entry plus XLA's own `memory_analysis()`.

## Results — experiment 1, uncontended single-process compute

Artifacts: `scratch/plateau/split_r{0,8,16,24,31}_n32.json`,
`split_r4_n8.json`, `split_r0_n1.json`. 12 timed reps each, median,
`taskset -c <one cpu>`, `EQDYNA_JAX_CACHE_DIR=off`.

**The probe is faithful.** At 1 rank it measures `ab = 753.00 ms/step` against
the committed `756.92` for the real `run_mpi` at `-np 1` — 0.5%. So a
single-process reproduction of a rank's local work IS the right instrument.

| rank of N | Ei | Ep | a+b | part_a | part_b | head | elem | div |
|---|---|---|---|---|---|---|---|---|
| 0 of 1 | 516096 | 218904 | **753.00** | 740.70 | 6.16 | 85.12 | 515.22 | 7.40 |
| 4 of 8 | 76608 | 23331 | **108.47** | 103.47 | 4.88 | 37.49 | 61.20 | 7.34 |
| 0 of 32 | 0 | 12217 | **36.85** | 32.80 | 4.10 | 24.01 | 10.97 | 7.41 |
| 8 of 32 | 19029 | 5874 | **50.05** | 45.45 | 4.34 | 27.36 | 12.13 | 7.37 |
| 16 of 32 | 19026 | 5875 | **51.62** | 46.97 | 4.50 | 35.98 | 20.87 | 15.65 |
| 24 of 32 | 19749 | 5634 | **50.11** | 45.85 | 4.08 | 28.79 | 12.25 | 7.23 |
| 31 of 32 | 0 | 12216 | **36.72** | 32.68 | 4.11 | 26.95 | 11.15 | 7.54 |

### C IS NOT INSIDE THE COMPUTATION. First hard finding.

| N | measured T | max-over-ranks local compute | T − compute |
|---|---|---|---|
| 1 | 756.92 | 753.00 | **3.92** |
| 8 | 211.29 | 108.47 | **102.82** |
| 32 | 143.27 | 51.62 | **91.65** |

The residual is ~4 ms at one rank and ~92-103 ms at 8 and 32 — it IS C, and it
appears the moment `nranks > 1`. The whole jitted step at 32 ranks, both
halves, on the SLOWEST rank, is 51.62 ms of a 143.27 ms step.

### Imbalance is closed a third way, independently

Uncontended compute across ranks of 32 spans **36.72 to 51.62 ms, a 1.41x
spread**, on the COMMITTED `PML_WEIGHT = 3.0` — the configuration whose
predicted work spread is 3.350x and which leaves `Ei = 0` on four ranks. A
1.41x compute spread converts to at most ~15 ms of barrier wait. The session
log's 4.65-105.09 ms of measured barrier wait is therefore not imbalance, and
this is why recutting to `w = 0.6` moved the 32-rank point 1.4%.

**`Ei = 0` costs a rank almost nothing**: a zero-interior rank is the FASTEST
(36.7 ms), because it still owns 12217 PML elements and PML is cheaper than
`PML_WEIGHT` claimed. `Ei > 0` was never a goal, confirming the session log's
own arithmetic from the other direction.

### What inside compute IS rank-independent — only a third of C

`head` (velDispUpdate + the force zero, `driver.py:40-83`, `:137`) goes
85.12 → 24-36 ms from 1 rank to 32: a 2.4-3.5x drop for a 32x work cut. `div`
(the mass divide, `driver.py:177`) is 7.40 ms at 1 rank and 7.41 at 32 —
**flat to two decimals**, as H1 predicted, because `mass` is global
(`driver.py:367`) and so is `force`. `part_b` is 6.16 → 4.10 ms.

So ~28-31 ms of rank-independent cost does sit inside compute and H1 named it
correctly — but C is ~92-103 ms. **H1 accounts for about a third of C and is
therefore not the answer.** Recorded as a partial hit, not a hit.

## Results — experiment 2, bytes/step and donation

**Donation is EFFECTIVE**, checked by asking the buffers rather than by
trusting the absence of a warning: `Array.is_deleted()` is `True` for `v1`,
`velArr`, `dispArr` and `force` after both `a_jit` and `b_jit`, at every rank
count measured. XLA confirms it independently — `b_jit`'s `memory_analysis()`
reports `alias = out = arg = 99.21 MB` with `temp` of 0.25 MB at 32 ranks. The
`~120 MB/step` copy `driver.py:410-417` warns about is NOT being paid.

Bytes per step per rank, from XLA's own `memory_analysis()`:

| | 1 rank | 32 ranks |
|---|---|---|
| `a_jit` arg / temp / out | 808.0 / 610.5 / 151.0 MB | 114.3 / 52.3 / 99.2 MB |
| `b_jit` arg / temp / out / alias | 151.0 / 2.9 / 151.0 / 151.0 MB | 99.2 / 0.25 / 99.2 / 99.2 MB |
| carry total | 150.97 MB | 99.21 MB |
| **carry GLOBAL-sized part** | **97.75 MB** | **97.75 MB** |

The last row is the structural fact: `v1 + velArr + dispArr + force` =
**97.75 MB, identical at 1 rank and at 32**, because their shapes come from
`S['N']`/`S['NEQ']` (`driver.py:382-389`). At 32 ranks they are **98.5% of the
entire carry** and 99.21 of the 114.3 MB `a_jit` reads. Per-rank per-step
traffic floors at ~100 MB in and ~100 MB out however many ranks there are, so
the AGGREGATE over 32 ranks (~6.4 GB/step, against ~2.2 GB/step at 1 rank)
GROWS with N instead of staying fixed.

## Results — experiment 3, bandwidth saturation, no MPI at all

N concurrent INDEPENDENT copies of the probe, each doing rank r of 32's local
work, zero communication, one per cpu on the 32 least-loaded cpus (chosen set
and worst busy fraction in `scratch/plateau/cpusel.txt`).
Artifacts `scratch/plateau/conc{8,32}_r*_n32.json`.

| concurrency | ab_ms min | median | **MAX** |
|---|---|---|---|
| 1 (uncontended) | 36.72 | 50.05 | 51.62 |
| 8 | 38.80 | 70.60 | 73.58 |
| 32 | 41.76 | 66.60 | **94.30** |

Paired, same rank, same work, 1 process vs 32:

| rank | uncontended | 32 concurrent | ratio |
|---|---|---|---|
| 0 | 36.85 | 43.99 | 1.19x |
| 8 | 50.05 | 54.70 | 1.09x |
| 16 | 51.62 | 72.18 | 1.40x |
| 24 | 50.11 | 62.97 | 1.26x |
| 31 | 36.72 | 49.87 | 1.36x |

**Contention is real and it is not the whole story.** It inflates a rank's
compute 1.09-1.40x and lifts max-over-ranks compute from 51.62 to 94.30 ms. It
also more than doubles the SPREAD (41.76-94.30 = 2.26x) on work that varies by
at most 1.41x — and in a barrier-synchronised solver the max sets the step, so
a bandwidth-induced spread is indistinguishable from a straggler while having
nothing to do with the decomposition. But 94.30 ms still leaves ~49 ms of the
143.27 ms step unattributed, and saturation is essentially complete by 8
processes (73.58 max at 8 vs 94.30 at 32) — which is the right SHAPE for a
curve that is flat from 4 ranks up rather than degrading progressively.

## In-situ 32-rank attribution, real MPI (`EQDYNA_MPI_STEP_PROFILE=1`)

Artifacts `docs/perf_snapshots/plateau_2026-09-22/attrib_n{40,160}.out` (all 32
ranks in each), snapshots `mpi_scaling_2026-09-22_06{0234,0401}_tpv104_prof{8,32}.json`,
2 ledger rows at sha `aa57928`. Direct `mpirun`, because
`run_mpi_scaling.py` captures rank stdout and prints it only on failure, so
the per-rank profile lines never reach a file through that path.

**Instrumentation does not perturb the number**: profiled 8 ranks reads
**212.43** ms/step against the committed unprofiled **211.29**, and profiled 32
ranks reads **123.07** (harness) / **124.61** (direct `mpirun`) in the same
session. Today's box is busier than the committed table's (load 17.8 vs 10.6),
which is why 32 ranks reads 123-125 here and 143.27 there.

### The table, by difference over 40 -> 160 steps (XLA compile removed)

Per-rank buckets sum to the step by construction, and they do:

| rank | compute | barrier | mpi | halo d2h | host | TOTAL |
|---|---|---|---|---|---|---|
| **20 (max compute — sets the step)** | **97.94** | 2.56 | 5.78 | 0.22 | 18.11 | 124.60 |
| 23 | 97.42 | 1.43 | 7.00 | 0.21 | 18.56 | 124.61 |
| 9 | 97.24 | 4.80 | 3.31 | 0.29 | 18.96 | 124.62 |
| 1 | 41.46 | 70.32 | 1.66 | 0.24 | 10.93 | 124.62 |
| 0 | 40.11 | 76.22 | 1.01 | 0.27 | 7.01 | 124.61 |
| 30 (min compute) | 33.73 | 75.11 | 9.16 | 0.27 | 6.34 | 124.62 |

Across all 32 ranks: compute 33.73-97.94, barrier 1.43-76.22, mpi 1.01-9.16,
halo d2h 0.21-0.32, host 5.98-18.96.

**On the rank that sets the step, 78.6% of the step is jitted compute**, 14.5%
is host residual, 4.6% MPI, 2.1% barrier, 0.2% the halo copy. The barrier is
the FAST ranks waiting for the slow one — it is 76 ms on the rank with the
least compute and 2.6 ms on the rank with the most, which is what a barrier
does and not a cost of its own.

## THE MECHANISM, NAMED

**The port's per-rank per-step memory traffic does not shrink with rank count,
so from 4 ranks up the step is memory-system-bound on replicated GLOBAL-sized
nodal state instead of compute-bound on the element work that was divided.
The rank-independent term C is that traffic, plus the 32-way contention it
causes; it is inside `compute`, and it is bytes, not instructions.**

Where it lives:

| file:line | what |
|---|---|
| `driver.py:382-389` | the carry is allocated at the GLOBAL `S['N']`/`S['NEQ']` on every rank. `v1+velArr+dispArr+force` = **97.75 MB, identical at 1 rank and at 32**, and **98.5% of the whole carry at 32 ranks**. |
| `driver.py:367` | `mass` is global, never restricted. |
| `driver.py:177` | the mass divide: **7.40 ms at 1 rank, 7.23-7.54 ms at 32** — flat. 91.6 MB of traffic at 12.4 GB/s. |
| `driver.py:137` | the force zero, all NEQ+1, every step, every rank. |
| `MPI4NodalQuant.py:183-190` | the decomposition restricts the element and node INDEX arrays; it does not restrict the ARRAYS they address. |
| `MPI4NodalQuant.py:36-52` | states the trade, and calls the full-length nodal arrays a MEMORY cost. Measured, it is a BANDWIDTH cost, and it is the plateau. |

Five measurements, each of which the naming would fail without:

1. **C is in compute, not orchestration.** Rank 20: compute 97.94 of a 124.60
   ms step; mpi 5.78; d2h 0.22; barrier 2.56; host 18.11. Transport is
   confirmed not the bound.
2. **The same rank-local work run ALONE takes 51.62 ms** (uncontended probe,
   max over ranks) against 97.94 in situ — a **1.90x inflation** that no
   property of the rank in isolation explains.
3. **32 concurrent INDEPENDENT processes with zero communication reproduce
   it**: max 94.30 vs in-situ 97.94, **3.7% apart**. So the inflation is not
   MPI, not the barrier, and not the exchange.
4. **The rank-independent floor is proportional to the GLOBAL problem size,
   not to per-rank work.** `div` is 1.94 ns per global equation on tpv104
   (7.40 ms / 3 818 583) and 1.86 ns on tpv8 (2.649 ms / 1 425 873) — 4% apart
   — and is rank-count-independent (7.40 at 1 rank, 7.23-7.54 at 32). This was
   pre-registered as the falsifiable prediction and it held.
5. **Donation is effective, so it is not the cause.** `is_deleted()` True on
   all four global arrays after both jits; `b_jit` `alias = arg = out =
   99.21 MB` with `temp` 0.25 MB. The ~120 MB/step copy `driver.py:410-417`
   warns about is not being paid.

Why Fortran does not have this: `meshgen.f90` gives each rank its OWN mesh with
rank-local numbering, so its nodal arrays SHRINK with rank count and its
working set fits cache progressively better — which is exactly the mildly
superlinear 37.33x it measures. The port keeps the proven serial mesh on every
rank and restricts it, which is why its indices cannot be wrong and why its
bytes cannot shrink. That is the trade, and at 32 ranks the trade is the bound.

### What I could NOT decompose, stated as a limit

I did not separate DRAM-bandwidth saturation from shared-last-level-cache
capacity, and I did not read hardware counters. What I measured is that
CONCURRENCY ALONE, with zero communication, reproduces the in-situ per-rank
compute. A fixed workload (rank 16 of 32) degrades only **1.37x** from
concurrency 1 to 32 (50.09 -> 68.80, and flat from K=2 to K=16), so uniform
bandwidth saturation is not a sufficient account of the **1.90x** on the worst
rank: part of the excess is POSITION-dependent (which cpu/NUMA node a rank
lands on), not uniform. Naming that split needs `perf` counters and is the
next measurement, not a conclusion I am entitled to here.

## CONTRADICTIONS of the framing I was given

1. **"Contention/placement: ALREADY RULED OUT WITH NUMBERS — per-rank spread
   1.00x at 32 ranks, EFFECTIVE_CORES 1.00 on 31 of 32." This ruling-out does
   not hold, and contention is the largest single term I measured.**
   - Per-rank *wall* spread is 1.00x because every rank sits on the same
     barrier. In a barrier-synchronised loop that is an IDENTITY, not evidence
     about compute. The compute spread underneath is **2.90x** (33.73-97.94).
   - `EFFECTIVE_CORES` cannot see this: a core stalled on cache misses still
     bills one cpu-second per wall second and reports 1.00. The metric is
     blind to the mechanism it was used to exclude.
   This is the one door the campaign closed that should be reopened.
2. **The dichotomy "either per-step work every rank does regardless of N, or
   per-step fixed overhead in the host/dispatch path" is a false one.** It is
   the former, but the work is BYTES, not instructions: the identical
   instruction stream costs 51.62 ms alone and 97.94 ms with 31 siblings.
   Neither branch describes a cost that is not a property of the rank in
   isolation.
3. **`T(N) = T1/N + C` mis-attributes part of C.** At 32 ranks the max-over-
   ranks *uncontended* compute is 51.62 ms, not the 23.65 ms `T1/N` predicts —
   because ~28-31 ms of it is the rank-independent replicated nodal work (H1).
   So C splits roughly into ~30 ms of replicated compute plus ~45 ms of
   contention plus ~18 ms of host residual, and the fit folds all three into
   one constant.
4. **"Barrier wait was a SYMPTOM" — of what is now specific.** It is a symptom
   of a MEMORY-induced compute spread (2.90x), not of an element-count spread
   (1.41x). That is why the `w=0.6` recut moved the 32-rank point 1.4%: it
   corrected a 1.41x spread that was never the problem.
5. **`Ei = 0` is not even a mild handicap.** A zero-interior rank is the
   FASTEST measured (36.72-36.85 ms uncontended) because it still owns ~12 200
   PML elements and PML is cheaper than `PML_WEIGHT` claims. The session log's
   "`Ei>0` is a canary, not a goal" is confirmed from the timing side.
6. **Agreeing with the framing:** transport really is not the bound (mpi
   1.01-9.16 ms, halo d2h 0.21-0.32 ms of a 124.6 ms step), XLA compile really
   is not a cliff (4.40-4.64 s), and imbalance really is closed.

## Parity

`driver.py` and `MPI4NodalQuant.py` are on a gated cell's path, so the gate was
re-run rather than argued about:

```
test.tpv8 x python-jax-mpi  SUCCESS  max|diff| = 1.220703e-10  bound 1.0e-08
test.tpv8 x python-jax      SUCCESS  max|diff| = 4.119873e-10  bound 1.0e-08
                            1891 fault nodes compared, 2/2 cells
```

Default behaviour is unchanged by construction: with `EQDYNA_MPI_STEP_PROFILE`
unset both additions are branches that do not execute, and the profile adds no
`block_until_ready`. No reference was regenerated; no bound was touched.

## Operational record

- All heavy runs launched detached writing to FILES; `/proc/1829034/fd/1`
  verified as a real path, not `pipe:[...]` (rule 20).
- Each invocation lands its own artifact; 8 and 32 ranks were separate
  invocations so a crash at 32 could not discard the 8-rank point (rule 20a).
- `EQDYNA_JAX_CACHE_DIR=off` on every run (the shared persistent XLA cache
  wedges multi-rank runs). No run wedged; no log mtime froze.
- Item 42 honoured: the attribution run WAITS on the sweep's own completion
  marker before taking a cpu rather than being launched alongside it.
- `--exclude-cpus 0,1`, ceiling 0.45; every selected cpu read busy 0.00
  (`cpusel*.txt`). Box load 12-18 during the runs, one foreign tenant.
- `scratch/` is gitignored, so every artifact was copied to
  `docs/perf_snapshots/plateau_2026-09-22/` to survive this worktree.

