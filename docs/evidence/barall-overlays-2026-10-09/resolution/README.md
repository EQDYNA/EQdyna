# Resolution study: is the FaultMod divergence mesh resolution or a defect? (2026-10-09)

The parent overlays (`../README.md`) found that none of tpv12/13/26/27/31/32/33 met the declared
MATCH criteria. They were measured at the gate dx and the 5 s gate term. This directory re-runs
every case at the gate dx, half of it and a quarter of it, each run to the **spec term**. It
then classifies each case as **RESOLUTION** (converges onto FaultMod as dx shrinks) or
**DEFECT** (does not converge, or diverges against the spec for a reason dx cannot remove).

The benchmarks have no analytic answer. FaultMod (Barall 2009, 2013, 2014, 2015, the SCEC CVWS
submissions) is the comparison code, not the truth.

## Runs

- **How each run was made.** Every run is EQdyna fortran on 4 ranks, a scratch run of the
  committed case with only `par.dx` changed.
  - They are not gated cells.
  - No committed reference was read, edited or regenerated (rule 7).
  - dx-derived quantities (`par.dt = 0.5*dx/vp`, fault grids, `nfx`/`nfz`) follow `par.dx`
    because the case files compute them from it.
- **Geometry and spec constraints.** All seven cases have planar faults, and none of their case
  files calls `lib.requireFaultGeometryResolution`, so there is no geometry floor. None of the
  chosen dx is finer than the case's own spec tier in `testsys/e2e/full_specs.py`.

| case | spec term | dx (m) | wall (s, 4 ranks, shared box) |
|---|---|---|---|
| tpv12 | 8 s | 500 / 250 / 125 | 80 / 561 / 3629 |
| tpv13 | 8 s | 500 / 250 / 125 | 60 / 594 / 3764 |
| tpv26 | 13 s | 500 / 250 / 125 | 129 / 968 / 7431 |
| tpv27 | 13 s | 500 / 250 / 125 | 138 / 985 / 7426 |
| tpv31 | 15 s | 500 / 250 / 125 | 119 / 846 / 6131 |
| tpv32 | 15 s | 500 / 250 / 125 | 124 / 847 / 5951 |
| tpv33 | 13 s | 400 / 200 / 100 | 90 / 627 / 4036 |

FaultMod's runs used dx 100 m (tpv12/13/26/27), 50 m (tpv31), 25 m (tpv32) and 12.5 m (tpv33).

## Criteria

These are the parent's declared MATCH criteria, unchanged:

1. At least 95% of the nodes FaultMod ruptures in-window are also ruptured in ours.
2. On exactly-shared nodes, median |dt| ≤ 0.2 s and p90 ≤ 0.5 s.
3. No station is ruptured by one code only.

Every number is computed on two windows: the 5 s gate term and the spec term.

- **Stations scored.** Stations are limited to those on the gate-dx grid, so all three dx of a
  case are scored on the same station set.
- **Arrival (one change from the parent's `numbers` mode).** Station arrival is measured on the
  vector rate:
  - on the fault, the first time hypot(h-rate, v-rate) > 1 mm/s;
  - off the fault, the first time max(|h|,|v|,|n| velocity) > 1 cm/s.

  The parent's h-rate-only arrival is blind on the dip-slip TPV12/13 fault, and there it
  produced false "one-code-only" stations. `numbers` mode is unchanged.

The full table, with every window, dx and column, is `table.md`. Per-station values are in
`metrics.json`; a missing arrival is `null`.

## Per-case trend (spec-term window)

| case | dx (m) | ruptured | median \|dt\| (s) | p90 (s) | mean (s) | one-only st | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|---|---|---|
| tpv12 | 500 | 93.7% | 0.399 | 0.606 | +0.556 | 0 | 0.92 | 0.95 | no |
| | 250 | 94.3% | 0.139 | 0.204 | +0.205 | 0 | 0.97 | 0.99 | no |
| | 125 | **94.3%** | 0.014 | 0.035 | +0.016 | 0 | 0.99 | 1.02 | no (ruptured) |
| tpv13 | 500 | 91.2% | 0.342 | 2.838 | +0.774 | 0 | 0.61 | 4.24 | no |
| | 250 | 93.5% | 0.127 | 0.483 | +0.320 | 0 | 0.70 | 3.95 | no |
| | 125 | 95.4% | 0.027 | 0.381 | +0.076 | 0 | **0.79** | **3.58** | yes (rupture time only) |
| tpv26 | 500 | 99.8% | 0.226 | 1.541 | −0.428 | 0 | 1.01 | 1.24 | no |
| | 250 | 99.9% | 0.148 | 1.548 | −0.462 | 0 | 1.01 | 1.53 | no |
| | 125 | 100% | 0.066 | 0.163 | −0.100 | 0 | 1.00 | 1.04 | **yes** |
| tpv27 | 500 | 99.8% | 0.644 | 1.896 | −0.738 | 2 | 0.98 | 1.26 | no |
| | 250 | 100% | 0.410 | 0.689 | −0.459 | 0 | 1.00 | 1.15 | no |
| | 125 | 100% | 0.184 | 0.360 | −0.245 | 0 | 1.00 | 1.07 | **yes** |
| tpv31 | 500 | 100% | 0.308 | 0.510 | +0.329 | 0 | 0.98 | 1.04 | no |
| | 250 | 100% | 0.058 | 0.141 | +0.062 | 0 | 0.99 | 1.01 | **yes** |
| | 125 | 100% | 0.018 | 0.048 | −0.005 | 0 | 1.00 | 1.03 | **yes** |
| tpv32 | 500 | 100% | 0.328 | 0.502 | +0.345 | 0 | 0.97 | 0.94 | no |
| | 250 | 100% | 0.077 | 0.156 | +0.082 | 0 | 0.98 | 1.03 | **yes** |
| | 125 | 100% | 0.018 | 0.049 | +0.012 | 0 | 0.99 | 1.03 | **yes** |
| tpv33 | 400 | 92.2% | 0.997 | 1.499 | +1.054 | 0 | 1.00 | 0.62 | no |
| | 200 | 96.0% | 0.170 | 0.346 | +0.161 | 0 | 0.99 | 0.81 | **yes** |
| | 100 | 98.0% | 0.052 | 0.114 | −0.002 | 0 | 1.00 | 0.93 | **yes** |

"slip ratio" and "body pv ratio" are station medians of ours / FaultMod. The 5 s-window rows
follow the same trend; see `table.md`. The one exception is tpv27's two 5 s one-only stations,
covered under Defect C.

## Verdicts

| case | verdict | one line |
|---|---|---|
| tpv12 | **DEFECT** (case setup) | Timing is RESOLUTION: median 0.40 → 0.14 → 0.014 s. Ruptured share is capped at 94.3% by Defect B. |
| tpv13 | **DEFECT** | Rupture timing converges. Fault slip (0.79) and surface velocity (3.6x) do not; see Defect A, plus B and C. |
| tpv26 | RESOLUTION | The p90 of 1.5 s was a spurious early supershear front at x > 0. It disappears at 125 m. |
| tpv27 | RESOLUTION | Timing 0.64 → 0.41 → 0.18 s. The "+55% slip" was a 5 s-window artifact; final slip ratio is 1.00. Carries Defect C. |
| tpv31 | RESOLUTION | MATCH at 250 m. |
| tpv32 | RESOLUTION | MATCH at 250 m. |
| tpv33 | RESOLUTION | The 1 s lag is gone at 200 m. Median body pv ratio goes 0.62 → 0.81 → 0.93; the stations within 0.8 km of the fault are still converging (per-station ranges 0.1–0.3 → 0.2–0.6 → 0.5–0.8). |

### TPV26: the early front at x > 0 is resolution

The parent table's tpv26 p90 (1.5 s) came from the long x > 0 side of the fault (the hypocentre
is at x = −5 km). Ours arrives early there.

| dx (m) | median dt, x 10–20 km |
|---|---|
| 500 | −1.34 s |
| 250 | −1.36 s |
| 125 | −0.13 s |

The near-flat 500 → 250 step looked like a defect, but the change is a threshold, not gradual.

- The S ratio is about 0.6, so a supershear transition is admissible.
- The static process zone at 10 km depth is about 0.8 km. That is about 3 elements at 250 m,
  and fewer once it contracts dynamically.
- Under-resolving it moves the transition earlier. At 125 m the transition no longer happens
  inside the fault, and the case meets MATCH.
- The forced-nucleation formula checks against TPV26_27_Description_v13 Part 5:
  `src/fortran/faulting.f90:405-436`.

## Defects (localized, not fixed)

### A. TPV13: a Method-2 initial-stress imbalance at the model's lateral boundary launches a P wave at t = 0

Spec TPV12_13_Description_v6, Part 6 (p.18): Method 2 is total stress plus explicit gravity, and
"If your code allows the mesh boundaries to move, then you must also apply external tractions."

`test.tpv13` uses Method 2 (`par.C_elastic = 0`, `case_input/test.tpv13/user_defined_params.py:118`)
inside a PML-bounded box. No external boundary traction exists anywhere in the code.

Relevant code:

- **TPV13 stress builder.** `src/fortran/meshgen.f90:1453-1480` (`setPlasticStress`). It writes
  the PML-slot copy into `+15` (s(16..21)) at `:1475-1478`.
- **PML element force.** `src/fortran/assembleGlobalKU.f90:254-282` (`calcPMLElemKU`). The split
  stresses s(1..15) drive f(1..9), and s0 = s(16..21) drives f(10..12).
- **Gravity body force.** `src/fortran/assembleGlobalKU.f90:21`.
- **Fault initial traction under Method 2.** `src/fortran/faulting.f90:129-137, 189-190`
  (`*C_elastic`). It comes from the body stress, not from `on_fault_vars`.

Evidence, all at dx 500 m and the spec 8 s term unless stated otherwise:

- **Not plasticity.** With yield removed (`par.coheplas = 1.0e12`), TPV13 still diverges from
  the Method-1 TPV12.

  | quantity | TPV13, yield removed | TPV12 | FaultMod TPV13 |
  |---|---|---|---|
  | normal stress, fault station dp039 | −185 MPa at 6 s | −41 MPa | −43 MPa |
  | slip at dp039 by 8 s | 12.7 m | 15.6 m | 17.4 m |
  | surface peak velocity, body-030st120dp000 | 13–20 m/s | 2.4–3.5 m/s | — |

  The real TPV13 run reaches −182 MPa and then +25 MPa (tensile) at dp039. Its surface peak
  velocity is 10–12 m/s, against FaultMod's 2–3.5 m/s.
- **It comes from the y-min boundary.**
  - The departure from TPV12 first appears at 3.33 s at the footwall station farthest from the
    fault. It then moves toward the fault and across it.
  - Moving `par.ymin` from −20 km to −35 km delays the departure by 2.65 s.
  - That is 15 km at the P speed of 5.716 km/s (2.62 s). It is a P wave leaving the y-min
    boundary at t ≈ 0.
- **Not resolution.** Across 500 → 250 → 125 m the slip ratio is 0.61 → 0.70 → 0.79 and the
  body pv ratio is 4.24 → 3.95 → 3.58. Neither approaches 1 at the rate the timing does
  (median 0.34 → 0.13 → 0.027 s).
- **TPV27 (also Method 2, elastic in the bulk).** It shows no fault-station departure from TPV26
  on normal stress. Its body stations depart by about 5% after 4–5.7 s.
- **TPV30.** The gated case also uses Method 2 and is likely exposed too. This was **not**
  measured here.

### B. TPV12/13: border fault nodes are pinned, against the spec

Spec TPV12_13_Description_v6, p.3: "A node which lies exactly on the border of the 30000 m ×
15000 m rectangle is considered to be inside the rectangle, and so should be permitted to
rupture."

Both cases set `mu_s = 1000` on |x| = 15 km and on the 15 km down-dip row:

- `case_input/test.tpv12/user_defined_params.py:207-208`
- `case_input/test.tpv13/user_defined_params.py:262-263`

At 125 m, every shared node that FaultMod ruptures and ours does not is a border node.

| case | border nodes among the misses | total shared nodes |
|---|---|---|
| tpv12 | 107 of 107 | 1879 |
| tpv13 | 87 of 87 | 1877 |

This caps tpv12's ruptured share at 94.3% at every dx, below the 95% criterion. Timing
otherwise converges to 0.014 s.

### C. Method-2 cases: the free surface sinks from step 1 (TPV13, TPV27)

Every body station (all at the free surface) in the two `C_elastic = 0` cases moves from the
first step.

| dx (m) | TPV13 \|v\| | TPV27 \|v\| | TPV12 and TPV26 (Method 1) |
|---|---|---|---|
| 500 | 1.0e-2 m/s | 1.08e-2 m/s | exactly 0 |
| 250 | 1.1e-2 m/s | 1.08e-2 m/s | exactly 0 |
| 125 | 1.1e-2 m/s | 1.08e-2 m/s | exactly 0 |

- **The signal.** The vertical velocity takes the same values per step at 500 m and 125 m:
  −7.48e-3 m/s at step 1, then −1.08e-2 and −9.75e-3 at steps 2 and 3. The vertical
  displacement then drifts at a steady −7e-3 m/s.
- **Effect on the scores.** This is what makes tpv27's two 3 km off-fault stations at x = 15 km
  read "one code only" in the 5 s window. Our "arrival" is at the first step (0.03–0.125 s),
  while FaultMod's is at 5.4–5.5 s.
- **Candidate cause (not proven).**
  - `setPlasticStress` is called at element-centre depth + 7.3215 m
    (`src/fortran/meshgen.f90:181`).
  - At the surface node, a stress offset of ρg·7.3215 m against a lumped mass ρ·dx/2 gives a
    per-step velocity increment of g·7.3215/vp ≈ 1.2e-2 m/s. That is independent of dx, which
    is the measured signature.

### Minor, not a verdict driver: TPV31 layer-jump nodes

The comment at `case_input/test.tpv31/user_defined_params.py:28-36` says no grid point falls
exactly on a stress-table jump depth (2400 / 5000 / 10000 m). At dx 250 and 125 m, nodes sit
exactly at 5000 and 10000 m. `np.interp` on the duplicated x-values gives the deeper value
there, while FaultMod takes the two-sided average. TPV31 still meets MATCH at 250 m.

### Checked and clean

- **t = 0 fault stresses.** Initial shear and normal stress at every fault station match
  FaultMod at depth for all seven cases. Surface-node differences are the dx/2 depth
  convention.
- **TPV26/27 forced nucleation.** Checked against the spec, as above.

## Figures

The figures are drawn by `scripts/figures/scec_compare.py` only. Each one overlays all three dx
plus FaultMod, and carries its term and run description on the figure.

- `<case>_cplot_res.png`: rupture-time contours, each model on its own grid.
- `<case>_ts_fault_res.png`: fault stations, h, v and n components.
- `<case>_ts_body_res.png`: body stations, h, v and n components.

## Rebuild

1. **Runs.** For each `(case, dx)` above, run from the repo root with
   `EQDYNAROOT=$PWD PATH=$PWD/bin:$PWD/scripts:$PATH`:

   ```
   scripts/create.newcase RUNS/<case>_dx<dx> test.<case>
   # edit RUNS/<case>_dx<dx>/user_defined_params.py: change only the number on the `par.dx = ...` line
   (cd RUNS/<case>_dx<dx> && ./case.setup && mpirun -np 4 $EQDYNAROOT/bin/eqdyna)
   python3 -m testsys.frt_canonical RUNS/<case>_dx<dx>
   ```

   The case's own `par.term` is the spec term, and the 5 s gate override is not applied here.
   Directory names must end in `_dx<number>`; anything else under RUNS is ignored.

2. **Defect-A diagnostics.** Append these lines to the dx 500 copy of `test.tpv13`:
   - yield removed: `par.coheplas = 1.0e12`;
   - wider footwall: additionally `par.ymin = -35.0e3`.

   Name the run directories `tpv13_dx500_nocoh` and `tpv13_dx500_nocoh_ymin35`.

3. **Metrics and figures.** FaultMod data is read, read-only, from `$SCEC_CVWS` (default
   `~/shared_dataset`).

   ```
   python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS docs/evidence/barall-overlays-2026-10-09/resolution
   ```

   Add `--no-figures` for the tables only, or list cases to restrict the run.
