# Wedge-mesh free-surface trace on test.tpv36 / test.tpv37 (diagnostic, 2026-11-12)

Evidence only. No source, case, test or reference change. Base: master 29c9a1d. Solver: the row149-force-budget
instrumented build (516da2c plus its dump patch; `git diff 516da2c 29c9a1d -- src scripts case_input` is empty, so master's solver
plus a dump). Scratch cases are `create.newcase` copies with dx, term and partition edited in the scratch copy only (`set_case.py`).

## Question

On the 60-degree TPV12 wedge mesh (row 149 wedge test, PR #269) the free-surface trace node `faultst000dp000` carried
+30 / +85 / +111 MPa of tension at 500 / 200 / 150 m, pinned at the cohesion, slip 1.57 m against 7.83 m at 150 m.
Do the registered 15-degree wedge cases test.tpv36 / test.tpv37 show tension, pinning or slip loss at their trace?

## Verdict: NO. No tension, no pinning, no slip loss at the trace on either case, at the gate dx or at 250 / 200 m.

| case | dx (m) | trace n-stress, max over run (MPa) | trace n-stress, min (MPa) | trace arrival (s) | trace slip at 13 s (m) |
|---|---|---|---|---|---|
| tpv36 wedge | 500 (gate) | +0.019 | -3.47 | 9.059 | 0.0892 |
| tpv36 wedge | 250 | +0.0023 | -5.40 | 7.371 | 0.1272 |
| tpv36 wedge | 200 | +0.0015 | -6.26 | 7.363 | 0.1331 |
| tpv36 archive (own 50 m, n-stress corrected) | 50 | +0.0002 | -16.96 | 7.366 | 0.1248 |

(Full table, all stations: `trace_summary.csv`. Positive n-stress is tension; the trace starts at exactly 0 so the small positive
maxima are the first steps, not tension. Run to 13 s for tpv36, 9 s for tpv37.)

- **tpv36**: the trace ruptures at 7.37 s (50 m archive 7.366 s) and slips 0.127 m at 250 m and 0.133 m at 200 m against 0.125 m in the 50 m
  archive; the 500 m grid is late (9.06 s) and short (0.089 m), an ordinary coarse-grid lag. Slip converges, it does not collapse.
  The compressive trace value grows with refinement (-3.5 / -5.4 / -6.3 MPa, then -17.0 at 50 m): a point value that is grid-dependent
  in the SAME direction as the 50 m submission, bounded by compression, with no pinning. One row down (`dp010`, 1 km down dip)
  slip at 13 s is 0.52 / 0.77 / 0.84 m at 500 / 250 / 200 m against 1.085 m at 50 m (-52 / -29 / -23 %), monotone convergence.
- **tpv37**: the trace never ruptures (cohesion 15 MPa at the surface; 50 m archive also never ruptures it). Trace n-stress stays between
  -0.47 and +0.02 MPa at 500 / 250 / 200 m (archive 50 m: -0.43 / 0.00), compressive, grid-independent to ~0.03 MPa.
- **Control that the binary and runs are sound**: the 500 m tpv36 and tpv37 runs reproduce the committed gate reference stations
  (`dp180`, `dp010`, `120dp030`: h-slip, v-slip, n-stress, 463 samples to 5 s) with max difference 0.000e+00 in all 12 rows
  (`control_vs_gate_reference.csv`). Exact zero is expected of the same deterministic CPU Fortran code at the same decomposition
  (2,1,2) and term-independent physics; provenance: different directories and files (run dirs vs `test.reference.results`).

## Why it does not carry over (what differs from the 60-degree TPV12 case)

Measured differences, cause not isolated (not a treatment mission; no run varied one parameter at a time):

| | TPV12 wedge test | tpv36/37 |
|---|---|---|
| dip | 60 | 15 |
| wedge cell | dy = dx/2, dz = 0.87 dx | dy = 0.97 dx, dz = 0.26 dx (nearly flat cell) |
| trace cohesion | 0.2 MPa (pin level) | 4 MPa (tpv36), 15 MPa (tpv37); strength cannot fall to 0.2 |
| trace initial stress | non-zero shear / normal ramp, rupture arrives early | exactly zero normal and shear at the trace; rupture arrives at 7.4 s, loads it dynamically |
| trace slip | 7.8 m (insertion) | 0.12 m (the trace barely slips even at 50 m) |

The row-149 trace tension was a TPV12-regime statement: a tens-of-MPa tensile point value at a node held at 0.2 MPa strength.
Neither tpv36 nor tpv37 reaches that regime in any run here.

## What the gate does and does not see

The gate runs at the 5 s term (`GATE_TERM_S`); the tpv36 trace is first reached at 7.37 s and tpv37's never. The committed station set
(`testsys/matrix.py` GATE_STATIONS: `dp180`, `dp010`, `120dp030` plus two off-fault) has no `dp000` and nothing within 1 km of the
surface at a time the rupture has arrived, so the trace is outside the gate in space and time. This is not a defect found: this run shows
the trace is benign through 13 s (tpv36) / 9 s (tpv37) at three grids. It is a statement that the gate would not catch one.
Not measured: tpv36/37 at finer than 200 m, term beyond 13 s (tpv36) / 9 s (tpv37), the force budget of the trace node.

## Not done / limits

- No independent oracle at the trace: the shared dataset has no TPV36/37 submission from another code; the only baseline is EQdyna's own 50 m
  archive (`scec_archive`, read-only), so this is MEASURED against own data, not VERIFIED. `scec_compare.py` overlays were not produced for that reason
  (an own-archive overlay of tpv36/37 exists under `case_input/test.tpv36/cmp_tpv36_dd.png`); figures here are `fig_tpv36_trace.png`, `fig_tpv37_trace.png`.
- Grids: 500 m at 4 ranks (2,1,2); 250 / 200 m at 16 ranks (4,1,4), one thread per rank, `flock`-serialised, load measured before each
  (13-30 on 64 cores, a foreign 8-rank job present). Walls in `runs.csv`.
- Rebuild: `queue.sh cases_coarse.txt`, `queue.sh cases_fine.txt` (env `WT_ROOT`, `WT_TREE`, `WT_LOCK`), then
  `wedge_trace_ana.py RUNS ARCHIVE REFERENCE OUT 13.0`. Scratch runs are deleted after this lands.
