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

It is validated two ways: it reproduces the numbers the 1024-rank 50 m bundle
was sized against (`scratch/tpv29/hpc50m/README.md`: node lines
985 x 247 x 493, 119,095,488 global, 119,164 max/rank, +2.5% imbalance), and it
matches the per-rank `totalNumOfElements` that every rank writes into
`compTime<me>` at 1, 4, 5, 16, 25 and 48 ranks in the sweep below
(`pred_per_rank_match` is `true` on every row).

The runtime ground truth is that `compTime<me>` file
(`library_output.f90:output_timeanalysis`); it needs `writeCompTime = 1`
(`globalvar.f90:109`), which is `0` in the default build.

## TPV29 (rough fault, elastic, 70x60x35 km domain, term = 20 s)

Element counts recomputed with `elem_per_rank.py` — the previous "44 M" for
100 m was wrong (it is 28.1 M; the 500 m and 200 m figures were right).

| dx (m) | ranks | decomposition | global elements | max elem/rank | rate (sim-s / wall-min) | wall (total) |
|---|---|---|---|---|---|---|
| 500 | 4 | 2,2,1 | 0.876 M | 219,006 | — | 2.5 min |
| 200 | 16 | 4,2,2 | 6.48 M | 410,625 | — | 11.6 min |
| 100 | 16 | 4,1,4 | 28.1 M | 1,780,920 | 0.173 | 1.93 h (projected from 21% at 23.8 min) |
| 100 | 48 | 8,1,6 | 28.1 M | 593,640 | (pending) | (pending) |
| 50 | 1024 | 16,8,8 | 119.1 M | 119,164 | — | ~0.7 h (estimate, LS6, unrun) |

Note the 500 m and 200 m rows predate the fault-on-MPI-boundary fix and used
`ny=2` with a symmetric y-domain, so their *results* are not comparable — their
*timings* are still valid as cost measurements.

## Is speed just a linear function of cells per rank?

(section filled in by the sweep below — see `elem_scaling_last.json`)

Provenance for the TPV29 rows: host cotopaxi (64 cores), Ubuntu 22.04,
gfortran 11.4.0, Open MPI 4.1.1, netCDF 4.8.1, EQdyna at 4944efe + the
fault-MPI fix, 2026-09-14, shared machine (load 8-23 during these runs).
