# Measured run speeds

Real timings from actual runs, kept because per-step speed is how we size jobs
and judge scaling. Every row carries its provenance (rule 6). Rates are
**simulated seconds per wall minute** — the quantity that scales with rank
count, unlike total wall time which also moves with `term`.

Read the rate straight off a run's own log: it prints `TimeElapsed (s)` each
output step, so

```
rate = (last TimeElapsed) / (wall minutes since case.setup finished)
projected total = wall * term / TimeElapsed
```

## How to get elements per rank (read this before quoting one)

Per-step cost is driven by elements **per rank**, so every performance claim
needs that number. Three numbers in the tree look like it and are not:

| source | what it actually is |
|---|---|
| `library.f90:memory_estimate` at HEAD: `1.54*totalNumOfElements/1e6*npx*npy*npz GB` | a **global** memory estimate — `totalNumOfElements` is rank 0's local count, multiplied straight back up by the rank count. Near-invariant with rank count at fixed `dx` *by construction*; that is the line working as designed, not a bug. Cannot be inverted to a per-rank count. |
| the working tree's rewritten `memory_estimate`, `Cells per rank (rank 0)` | a true per-rank count, but **rank 0's**. Rank 0 is `(mex,mey,mez)=(0,0,0)`, which gets the *smallest* slice whenever a dimension does not divide evenly (`meshgen.f90:476-480` hands the extra node line to the high-id ranks). Per-step cost is gated by the **slowest** rank, i.e. the max. |
| `scripts/case.setup:estimate_HPC_resource` | wrong twice: counts only the `±nuni_y` cells in y (dropping the geometric coarsening zone and the PML), then multiplies a global count by the rank count. |

Worked example, TPV29 `dx = 500`, `4,1,4`: rank 0 has **53,176**; the max rank
has **57,960** (+9.0%); the true global is **876,024**, not
`53,176 x 16 = 850,816` (−2.9%). Both numbers that have been quoted for this
case — 53,176 and 54,751 — come from this one configuration: the first is
rank 0's share, the second is `global/nranks`.

**Use `testsys/perf/elem_per_rank.py`.** It is a line-by-line replica of
`meshgen.f90:getLocalOneDimCoorArrAndSize` plus `countMeshEntities.f90` and
reports max / min / mean per rank plus the global count, with no run needed:

```
python3 testsys/perf/elem_per_rank.py --case-dir <set-up case> [--np NPX NPY NPZ]
```

Validated two ways: it reproduces the numbers the 1024-rank 50 m bundle was
sized against (`scratch/tpv29/hpc50m/README.md`: node lines 985 x 247 x 493,
119,095,488 global, 119,164 max/rank, +2.5% imbalance), and it matches the
per-rank `totalNumOfElements` that every rank writes into `compTime<me>` at 1,
4, 5, 16, 25 and 48 ranks — `pred_per_rank_match` is `true` on all 15 rows of
the sweep below, i.e. the predicted multiset of per-rank counts equals the
measured one exactly, rank for rank.

The runtime ground truth is that `compTime<me>` file
(`library_output.f90:output_timeanalysis`); it needs `writeCompTime = 1`
(`globalvar.f90:109`), which is `0` in the default build. It also carries the
per-stage `MPI_WTIME` breakdown used below.

## TPV29 (rough fault, elastic, 70x60x35 km domain, term = 20 s)

Element counts recomputed with `elem_per_rank.py` — the previous "44 M" for
100 m was wrong (it is 28.1 M; the 500 m and 200 m figures were right), and the
"~0.9 M per rank at 100 m / 48 ranks" note was wrong for the same reason
(0.59 M).

| dx (m) | ranks | decomposition | global elements | max elem/rank | rate (sim-s / wall-min) | wall (total) |
|---|---|---|---|---|---|---|
| 500 | 4 | 2,2,1 | 0.876 M | 219,006 | — | 2.5 min |
| 200 | 16 | 4,2,2 | 6.48 M | 410,625 | — | 11.6 min |
| 100 | 16 | 4,1,4 | 28.1 M | 1,780,920 | 0.173 | 1.93 h (projected from 21% at 23.8 min) |
| 100 | 48 | 8,1,6 | 28.1 M | 593,640 | 0.573 | **34.9 min (measured)** |
| 50 | 48 | 8,1,6 | 119.1 M | 2,481,156 | 0.529 | **5.04 h (projected from 914 of 4800 steps)** |
| 50 | 1024 | 16,8,8 | 119.1 M | 119,164 | — | ~0.7 h (estimate, LS6, unrun) |

The 100 m / 48-rank row is now a real run: `scratch/tpv29/ny48_dx100`,
2400 steps between 17:23:18 and 17:58:13 on 2026-09-14 = **0.873 s/step**, i.e.
**1.47 µs per element per timestep** at 593,640 elem/rank. That is an
independent case (rough fault, a different compset, timed with nothing but
`date` and the step counter) landing inside the 1467–1881 ns/elem band the
`test.tpv8` 48-rank rows below give — a useful confirmation that the 48-rank
penalty is a property of the machine, not of `test.tpv8`.

The 50 m / 48-rank row is `scratch/tpv29/spec50m`, started 19:36:48 and stopped
by hand at 20:34:23 on 2026-09-14 after 914 of 4800 steps (the run was ended to
free the node, not because it failed): 3455 s / 914 steps = **3.780 s/step** at
2,481,156 elem/rank = **1524 ns per element per timestep**. Memory was 3.82 GB
per rank, 183 GB across the 48 ranks. The last ~1 min of the 57 overlapped an
e2e run starting; at 1/57 of the wall that is below the resolution of this
measurement.

### The cleanest test of the linear claim: hold the rank count fixed

The two TPV29 rows above are the same case, same binary, same node, same 48-rank
decomposition (8,1,6) — the *only* thing that changes is resolution, and hence
elements per rank. That isolates the variable the claim is actually about, which
the `test.tpv8` sweep below could not do (there, rank count and per-rank size
move together with the machine's bandwidth behaviour).

| dx | mean elem/rank | s/step | ns/elem/step |
|---|---|---|---|
| 100 m | 584,918 | 0.873 | 1492 |
| 50 m | 2,481,156 | 3.780 | 1524 |

**4.242× the elements per rank costs 4.330× the time — +2.1% off proportional**,
and the per-element rate moves by 2.1% across a 4.2× range in per-rank size.

This is the strongest support the claim has: *at fixed rank count*, wall time per
step really is a linear function of elements per rank, to about 2%. It does not
contradict the refutation below, which is about a different variable — there, the
rank count changes and the fitted slope moves 1.48× with it. Stated precisely:

> The slope is constant in elements-per-rank (±2%), but the slope itself is a
> function of how many ranks share the node (±48%).

So "cells per np" predicts cost well when you are choosing a resolution for a
fixed job size, and predicts it badly when you are choosing how many ranks to
run on one node.

Note the 500 m and 200 m rows predate the fault-on-MPI-boundary fix and used
`ny=2` with a symmetric y-domain, so their *results* are not comparable — their
*timings* are still valid as cost measurements.

## Is speed just a linear function of cells per rank?

**Verdict: qualified — true of the compute kernel, false of the wall clock.**
Per-step cost is *not* a function of elements per rank alone. How the elements
arrive matters, and the deviation reaches **+51%** on this node.

### The sweep

`test.tpv8`, `dx ∈ {500, 250, 125}` crossed with `ranks ∈ {1, 4, 5, 16, 25, 48}`
so that resolution and rank count are decorrelated — the same elements/rank is
reached from several `(dx, ranks)` pairs. `npy = 1` in every decomposition, so
no MPI partition boundary can land on the fault plane (the known defect class,
avoided by construction rather than relied on as fixed). Range covered:
**4,896 → 5,858,146 elements/rank, a factor of 1,197**.

Per-step cost is taken from the run's own `MPI_WTIME` clocks in `compTime<me>`:
`LOOP = (comp(9) − comp(1) − comp(8)) / nstep` (whole program, less setup/mesh,
less output), split into `kernel` (comp 3..6 — velDispUpdate, assembleGlobalKU,
calcHourglassResist, faulting; accumulated only inside `driver.f90`'s loop) and
`halo` (MPICommTime — the sendrecv plus the six `mpi_barrier`s per step). Wall
clock is carried alongside: it agrees with LOOP to within 1.3% on every row
at or above 125,091 elem/rank and to 4.3% at 58,752, and diverges below that
(+11% at 25,410, +48% at 4,896) exactly as it should — what wall adds is the
mesh build, netCDF read and `mpirun` startup that LOOP excludes, and those are
fixed costs being divided by a small step time.

`nstep` is recorded explicitly on every row and confirmed from the run's own
output: the last `TimeElapsed (s)` the run printed equals `term` to the seven
digits the format carries, and `term/dt = nstep` exactly, on all 15 rows. (Do
not use the *number* of `TimeElapsed` lines as the step count —
`faulting.f90:showSourceDynamics` fires once per fault node pair matching the
hypocentre, so the line count comes out at 1x, 2x or 4x `nstep` depending on
`dx`. The last *value* is the reliable one.)

| dx (m) | ranks | decomp | nstep | elem/rank (max) | imbal | ms/step LOOP | ms/step kernel | ms/step halo | ms/step wall | ns/elem/step |
|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 48 | 8x1x6 | 200 | 4,896 | +0.0% | 9.187 | 6.762 | 3.333 | 13.556 | 1877 |
| 500 | 16 | 4x1x4 | 200 | 14,688 | +0.0% | 19.023 | 18.219 | 0.987 | 22.204 | 1295 |
| 250 | 48 | 8x1x6 | 200 | 25,410 | +2.3% | 37.268 | 33.972 | 7.129 | 41.556 | 1467 |
| 500 | 4 | 2x1x2 | 200 | 58,752 | +0.0% | 75.293 | 73.279 | 1.680 | 78.522 | 1282 |
| 250 | 16 | 4x1x4 | 200 | 74,536 | +0.0% | 97.683 | 94.805 | 3.996 | 101.620 | 1311 |
| 125 | 48 | 8x1x6 | 200 | 125,091 | +2.5% | 235.277 | 219.703 | 75.717 | 240.985 | 1881 |
| 500 | 1 | 1x1x1 | 200 | 235,008 | +0.0% | 292.817 | 289.791 | 0.000 | 296.687 | 1246 |
| 125 | 25 | 5x1x5 | 200 | 242,385 | +3.4% | 323.668 | 312.447 | 23.975 | 329.100 | 1335 |
| 250 | 5 | 5x1x1 | 200 | 243,936 | +2.3% | 306.649 | 300.598 | 10.401 | 310.662 | 1257 |
| 250 | 4 | 2x1x2 | 200 | 298,144 | +0.0% | 383.556 | 373.860 | 12.407 | 387.756 | 1286 |
| 125 | 16 | 4x1x4 | 200 | 375,273 | +2.5% | 494.389 | 478.136 | 26.845 | 500.698 | 1317 |
| 125 | 5 | 5x1x1 | 200 | 1,182,545 | +0.9% | 1542.576 | 1498.021 | 45.740 | 1550.779 | 1305 |
| 250 | 1 | 1x1x1 | 200 | 1,192,576 | +0.0% | 1538.485 | 1522.818 | 0.000 | 1546.982 | 1290 |
| 125 | 4 | 2x1x2 | 120 | 1,473,633 | +0.6% | 1909.059 | 1866.886 | 53.259 | 1924.014 | 1296 |
| 125 | 1 | 1x1x1 | 80 | 5,858,146 | +0.0% | 7565.092 | 7447.337 | 0.000 | 7644.894 | 1291 |

### The fit

```
Global OLS       time/step = 7.80 ms + 1290.0 ± 3.7 ns × (elem/rank),  R² = 0.99990
Forced through 0 time/step =            1292.3 ns × (elem/rank),        R² = 0.99988
Relative-error   time/step = 2.40 ms + 1298.3 ns × (elem/rank)
   (1/y² weighted)          rms relative residual 8.9%, max 29.9%
Power law        time/step = 2.11e-6 × (elem/rank)^0.9644 ± 0.0169
```

**R² = 0.9999 here is worthless and must not be quoted.** Over a 1,197× range an
unweighted fit is determined entirely by the three largest points; it reports a
near-perfect line while the smallest configuration sits **−54%** off it. The
OLS residuals make this plain: −53.7%, −40.6%, −11.0%, −8.9% … +28.1%. The
honest statements are the relative-error fit (rms residual 8.9%, worst 30%) and
the power-law exponent **p = 0.964 ± 0.017, which is 2.1σ below 1**.

The intercept is real but small: **+2.4 ms/step** on the relative-error fit,
which is **26% of the fastest configuration's step time** (9.19 ms at 48 ranks
× 4,896 elem) and under 0.2% of the slowest. So fixed per-step overhead does
bite, but only below ~25k elements/rank — and it is *not* a constant, because
it grows with rank count (see below).

### Per-rank-count sub-fits — the load-bearing result

| ranks | n points | slope (ns/elem/step) | R² |
|---|---|---|---|
| 1 | 3 | 1292.8 | 1.00000 |
| 4 | 3 | 1296.6 | 1.00000 |
| 5 | 2 | 1316.8 | — |
| 16 | 3 | 1318.6 | 1.00000 |
| 48 | 3 | 1914.6 | 0.99788 |

**The slope is not a constant: it spans 1.48× across rank counts on one node.**
Under the hypothesis it would be one number. 1 → 16 ranks costs +2%; 48 ranks
costs +48%.

### Matched sets — the hypothesis stated as a falsifiable prediction

Same elements per rank, reached at different rank counts *and* different
resolutions. The hypothesis says every member must cost the same per step.

| elem/rank | dx (m) | ranks | ms/step | ns/elem/step | vs 1 rank |
|---|---|---|---|---|---|
| 235,008 | 500 | 1 | 292.8 | 1246 | — |
| 243,936 (+3.8%) | 250 | 5 | 306.6 | 1257 | **+0.9%** |
| 242,385 (+3.1%) | 125 | 25 | 323.7 | 1335 | **+7.2%** |
| 1,192,576 | 250 | 1 | 1538.5 | 1290 | — |
| 1,182,545 (−0.8%) | 125 | 5 | 1542.6 | 1305 | **+1.1%** |

A 4× change in resolution with a 5× change in rank count, at fixed
elements/rank, costs **1%**. That part of the hypothesis is solidly true. The
25-rank member costs 7% more, and the 48-rank rows (no matched partner, but
1467–1881 ns/elem against 1246–1296 for 1–5 ranks) cost up to 51% more.

### Mechanism — why it deviates

The `kernel`/`halo` split separates the two causes cleanly:

| ranks | kernel ns/elem/step | halo ns/elem/step | halo as % of step |
|---|---|---|---|
| 1 | 1233, 1277, 1271 | 0 | 0.0% |
| 4 | 1247, 1254, 1267 | 29, 42, 36 | 2.2–3.2% |
| 5 | 1232, 1267 | 43, 39 | 3.0–3.4% |
| 16 | 1240, 1272, 1274 | 67, 54, 72 | 4.1–5.4% |
| 25 | 1289 | 99 | 7.4% |
| 48 | 1381, 1337, 1756 | 681, 281, 605 | 19.1–36.3% |

1. **The compute kernel is proportional to elements per rank, to ±4.6%.**
   Across 1–25 ranks and a 400× range in per-rank size, `kernel` costs
   **1232–1289 ns/element/step**. This is the true content of the hypothesis
   and it holds well.
2. **Halo exchange is a rank-count term, not an element term.** It is
   identically zero at 1 rank and climbs monotonically with rank count —
   ~35 ns/elem at 4–5 ranks, ~65 at 16, ~99 at 25, 281–681 at 48. At 48 ranks
   it is a fifth to a third of the entire step. `MPI4NodalQuant` issues six
   `mpi_barrier`s per timestep, so this term also absorbs every source of
   inter-rank skew, and skew grows with the number of ranks being
   synchronised.
3. **At 48 ranks the kernel itself slows down**, to 1337–1756 ns/elem
   (+6% to +39% over the 1–25-rank band), with the worst case at the largest
   per-rank size (125,091 elem). Same instruction stream, same elements per
   rank, more cores pulling on the same memory — this is memory-bandwidth and
   shared-L3 saturation, not element count. The box is 2×32-core EPYC 7532
   with 16 MB of L3 per 4-core CCX, so at 48 ranks four ranks share each L3
   slice and eight share each memory controller.

### Verdict

- **Refuted as stated.** "Wall time per timestep is proportional to elements
  per rank, independent of how those elements arrive" is false: at fixed
  elements/rank, going from 1 to 48 ranks on one node costs up to **+51%** per
  element, and the fitted slope varies **1.48×** with rank count.
- **Qualified, and useful, in this form:** per-step cost is
  `≈ 1.25 µs × (elements per rank)` for the compute kernel, good to ±5% from
  1 to 25 ranks over a 400× range in per-rank size — plus a halo term that
  grows with rank count and is negligible (<5%) up to 16 ranks, plus a
  bandwidth penalty that appears when the node is filled.
- **Practical rule for sizing jobs on this box:** use 1.3 µs/element/step up to
  ~16 ranks per node; use 1.9 µs/element/step at 48 ranks per node; and do not
  size below ~25k elements/rank, where the fixed per-step overhead (2.4 ms) is
  a visible fraction of the step.
- **The intercept is real but secondary.** 2.4 ms/step of fixed cost only
  matters below ~25k elements/rank. The dominant deviation is the rank-count
  dependence of the *slope*, not the intercept.

### Caveats

- One node, one case (`test.tpv8`), one compiler. The rank-count penalty is an
  on-node effect (shared L3 and memory controllers); it says nothing about
  scaling across nodes, where per-rank bandwidth is restored and the halo goes
  over a network instead of shared memory.
- The box was shared throughout: two GPU training jobs held the background load
  at 6–11 for the whole sweep (recorded per row), so the 48-rank figures —
  which use 48 of 64 cores and are the most sensitive to any competition for
  bandwidth and to barrier skew — are an **upper bound** on the penalty a
  dedicated node would show. The low-rank rows are insensitive to this and are
  internally consistent to ±4%.
- Two rows (`250/4`, `125/4`) had another eqdyna job appear before they
  finished (`contended = 4`). Both show `loop_spread = +0.0%` across their
  ranks and per-element costs (1286, 1296 ns) inside the 1–5-rank band, so
  neither is contaminated; they are kept and flagged.
- `nstep` differs on the two largest configurations (120 and 80 instead of 200)
  to keep the sweep affordable. This is controlled: `125/4` at nstep = 120
  gives 1296 ns/elem against 1282 and 1286 for the two 4-rank rows at
  nstep = 200, so per-step cost is step-count independent as expected.
- `compTimeInSeconds(2)` (mass assembly) is corrupted and is **not** subtracted
  from LOOP: `assembleGlobalMass → MPI4NodalQuant` resets the shared global
  `startTimeStamp` (`assembleGlobalMass.f90:67`), so `eqdyna3d.f90:69` times
  only the tail of the call. The real mass-assembly cost therefore sits inside
  LOOP; it is a one-off, O(elements), worth roughly one timestep, so it biases
  per-step cost by ~1/nstep (≤1%) and, being element-proportional, adds no
  spurious intercept.

### Reproducing

```
# build with writeCompTime=1 (one-character change to globalvar.f90:109)
cp src/*.f90 src/makefile <scratch>/src && sed -i 's/writeCompTime = 0/writeCompTime = 1/' <scratch>/src/globalvar.f90
cd <scratch>/src && MACHINE=ubuntu make eqdyna

python3 testsys/perf/run_elem_scaling.py --bin <scratch>/src/eqdyna \
    --work <scratch>/work --out testsys/perf/elem_scaling_last.json \
    --nstep 200 --load-ceiling 14 --max-others 0 \
    --grid 500:48,500:16,250:48,500:4,250:16,125:48,500:1,250:5x1x1,\
125:5x1x5,250:4,125:16,250:1,125:5x1x1,125:4,125:1 \
    --nstep-override 125:4:120,125:1:80
python3 testsys/perf/analyze_elem_scaling.py      # table, fits, elem_scaling.png
```

Raw rows (every per-rank element count and timing slot):
`testsys/perf/elem_scaling_last.json`. Plot: `testsys/perf/elem_scaling.png`.
Case trees and run logs: `scratch/perf_elemscaling/work/`.

For the older controlled strong-scaling sweep with pinned ranks, use
`python3 testsys/run.py scaling` (harness in `run_scaling.py`), which runs one
configuration at a time and writes `scaling_last.json`.

## Provenance (rule 6)

- **Elements-per-rank sweep and its fits**: `test.tpv8`, EQdyna at `12c8274`
  with `src/` dirty (the fault-MPI-boundary fix in the working tree) plus
  `writeCompTime = 1`; host `cotopaxi`, 2 × AMD EPYC 7532 32-core (64 cores,
  8 NUMA nodes, 16 MB L3 per 4-core CCX), Linux 6.8.0-107, Ubuntu 22.04;
  `mpif90` GNU Fortran 11.4.0, `FFLAGS = -fopenmp -ffree-line-length-none -O3`
  (`src/makefile`, `MACHINE=ubuntu`); Open MPI 4.1.1, netCDF 4.8.1;
  `OMP_NUM_THREADS=1`; ranks pinned 1:1 to cores `0..n-1` with
  `mpirun --bind-to core --map-by core`; 2026-09-14 18:13–19:21 CDT; shared
  machine, background load 6–11 (recorded per row).
- **TPV29 rows**: same host/toolchain, EQdyna at `4944efe` + the fault-MPI fix,
  2026-09-14, shared machine (load 8–23 during those runs). The 100 m/48-rank
  row was measured at background load ~7 from the GPU jobs plus nothing else.
