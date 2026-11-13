# Surface strength on the wedge-mesh trace: one-variable tests on TPV12 (60 deg) and test.tpv36 (15 deg) (diagnostic, 2026-11-12)

Evidence only. Row 149 stays P4 paused. No source, case, test or reference change; scratch cases are `create.newcase` copies
edited in the scratch copy only. Base: master 7f7984c. Solver: the row149-force-budget instrumented build (master's solver plus a dump).

## Question (owner decision, board 297b9c9)

"Raise TPV12's surface strength to 4 MPa and see whether the tension and slip loss go away; or lower TPV36's to 0.2 MPa and see
whether they appear." Both are run.

## The one change per test (state before running)

- **Test 1, TPV12 wedge 60 deg.** TPV12 sets cohesion UNIFORM over the fault: `par.fric_cohesion = 0.2e6`, copied into every fault
  node at `on_fault_vars[iz,ix,4]` (`case_input/test.tpv12/user_defined_params.py:137`, `:223`; the row-149 wedge variant
  `docs/evidence/row149-wedge-mesh-2026-11-11/user_defined_params.py:148`). The surface value is the only value, so the change is
  that one constant: `par.fric_cohesion     = 0.2e6` -> `4.0e6` in the wedge variant (and, for a matched control, in the insertion
  case at `:137`). The TPV36 shape (4 MPa tapering to 0 at 8 km) was NOT used: it would also remove the 0.2 MPa below 8 km.
- **Test 2, test.tpv36 wedge 15 deg.** TPV36 sets cohesion `cohesionSlope*(8e3 - downDipDistance)` for down-dip distance <= 8 km,
  0 below (`case_input/test.tpv36/tpv36_37_common.py:113-114`), `cohesionSlope=0.0005e6` (`case_input/test.tpv36/user_defined_params.py:7`),
  i.e. 4 MPa at the trace. The change is that slope: `0.0005e6` -> `0.000025e6`, 0.2 MPa at the trace, same linear shape and 8 km
  extent. dx, term and partition are set per run exactly as in `docs/evidence/wedge-trace-tpv36-37-2026-11-12/set_case.py`.

No other line differs from the baseline runs (each run's `set_case` prints every changed line; checked into `runs.csv` via the run dir).

## Kill criteria (written and committed before any run; not adjusted afterwards)

Strength **is the cause** iff ALL of:
1. Test 1 removes the tension: TPV12 wedge `faultst000dp000` n-stress at 4 s is compressive, or within a few MPa (<= 5 MPa) of the
   insertion run at the same dx, at 200 and 150 m;
2. Test 1 removes the slip loss: wedge dip slip at 4 s at 1.5 km down dip (`faultst000dp015`, 150 m grid) within 5 % of the insertion
   run (both at 4 MPa; the 0.2 MPa insertion run is quoted beside it);
3. Test 2 produces tension or pinning at the tpv36 trace: `faultst000dp000` n-stress max over the run > +1 MPa after arrival, or
   |shear| held at the strength (ratio >= 0.99) while the node is in tension, at 250 and 200 m.

Strength is **not the cause** iff test 1 leaves tension (> +5 MPa above the insertion value) AND slip loss (>= 5 % at 1.5 km) and test 2
shows neither. Any other combination is reported as **mixed**, with the numbers.

## Measured quantities

Trace n-stress (positive = tension) at the evaluation time (TPV12 4 s, tpv36 13 s) and its maximum after arrival; pinning ratio
|shear| / (c - mu(slip) min(sigma_n, 0)) with the slip-weakening mu of the case; slip at the trace and down dip (TPV12: 1.5 and 3 km,
`dp015`, `dp030`; tpv36: 1 km, `dp010`); trace arrival (1 mm/s); per-side stress-divergence / hourglass parts of the trace-node
traction from the instrumented dump. Baselines re-derived from raw runs by the same script before any comparison.

## Amendment 1 (2026-11-12, before any Test 1 result at 200/150 m; the kill criteria above are NOT changed)

The pre-registered Test 1 change (`par.fric_cohesion` 0.2e6 -> 4.0e6, uniform over the fault) was run first at 500 m on the wedge
mesh (4 ranks, 120 s, rc=0) and is **degenerate**: every fault node has slip 0 at every station, `frt.txt*` rupture time = 99999
everywhere, n-stress = 0 at the trace. TPV12's nucleation patch (10.5-13.5 km down dip, `mu_s` 0.54) is only just over strength with
0.2 MPa cohesion, so a uniform 4 MPa removes the nucleation too; criteria 1 and 2 would then pass vacuously (0 MPa of tension, 0 m of
slip loss on a fault that never moves), a green that tests nothing. The 200/150 m cells of that design were not run. That run is kept as
`w12u4_dx500`.

Test 1 is therefore run as the surface-only reading of the owner's wording ("raise TPV12's surface strength to 4 MPa"): cohesion
`0.2 + 3.8 * max(0, 1 - downDipDistance/8 km)` MPa, i.e. 4 MPa at the trace falling linearly to the **unchanged** 0.2 MPa at 8 km down
dip and below (the TPV36 taper shape with a 0.2 MPa floor, so nothing below 8 km differs from the baseline, nucleation patch
included). One line changes (`par.on_fault_vars[iz,ix,4] = ...`, wedge variant `:234`, insertion case `:223`), in the scratch copy only.
The criteria read the same: trace n-stress at 4 s, slip at `dp015`, both against the insertion run at 4 MPa-at-trace, plus a guard that
a run with no rupture is reported DEGENERATE, not PASS. Test 2 is unchanged (its taper already extends only to 8 km).

Other execution changes, none of which touch a criterion: 16 ranks (4,2,2) for the 200/150 m TPV12 cells, validated by a 500 m
16-rank control of the baseline wedge case (`w12c02_dx500_p422`) against the committed 4-rank run; the 4 MPa tpv36 baseline is rerun
only at 500 m (the 250 m and 200 m values come from the committed wedge-trace evidence); the insertion 0.2 MPa rerun at 150 m is the
last queue entry (committed row-149 numbers are quoted otherwise).

## Verdict (2026-11-12): strength is NOT the cause (criteria applied as written)

Baselines re-derived from the raw runs by `ss_ana.py` reproduce the committed numbers to every printed digit (w12 dx500/200/150 dp000 n-stress
+30.0603 / +84.9702 / +110.5110 MPa; insertion 0.2 MPa dx150 -2.2977 MPa; tpv36 4 MPa dx500 n max +0.019248 MPa), and the 16-rank (4,2,2) 500 m
wedge control equals the committed 4-rank (2,2,1) run at dp000/dp015/dp075 to the printed 6 digits, so the rank change is not a variable.
The headline numbers below were also re-read from the raw `faultst*.txt` station files with an independent awk pass.

| quantity (4 s, TPV12 wedge; 13 s, tpv36) | 500 m | 250 m | 200 m | 150 m |
|---|---|---|---|---|
| TPV12 wedge c=0.2 baseline, dp000 n-stress (MPa) | +30.06 | | +84.97 | +110.51 |
| **TPV12 wedge c=4 taper (test 1)**, dp000 n-stress | +28.31 | | +76.79 | +98.82 |
| TPV12 insertion c=4 taper, dp000 n-stress | -3.88 | | +2.06 | +3.82 |
| C1 excess, wedge minus insertion (limit 5 MPa) | +32.2 | | +74.7 FAIL | +95.0 FAIL |
| wedge c=4 shear/strength at trace, tensile-and-pinned fraction after arrival | 0.83 | | 0.95 | 0.95 |
| dp000 dip slip: wedge c=4 / insertion c=4 (m) | 1.59 / 5.21 | | 1.39 / 6.29 | 1.32 / 6.48 |
| dp015 (1.5 km) dip slip: wedge c=4 / insertion c=4 / insertion c=0.2 (m) | 3.80 / 5.74 / n.a. | | no dp015 station | 5.18 / 7.10 / 8.38 |
| C2 slip loss at dp015, wedge vs insertion c=4 (limit 5 %) | 33.7 % | | | 27.0 % FAIL |
| dp030 (3 km) dip slip: wedge c=4 / insertion c=4 (m) | 6.05 / 6.55 | | 6.97 / 7.87 | 7.15 / 8.06 |
| tpv36 c=4 baseline, dp000 n max / slip at 13 s (m) | +0.019 / 0.089 | +0.002 / 0.127 | +0.001 / 0.133 | |
| **tpv36 c=0.2 (test 2)**, dp000 n max after arrival (MPa) | -0.147 | -0.116 | -0.107 | |
| tpv36 c=0.2, dp000 slip at 13 s (m), arrival (s) | 0.481, 6.08 | 0.437, 6.09 | 0.423, 6.10 | |
| tpv36 c=0.2, tensile-and-pinned fraction (C3: > 1 MPa tension or pinned in tension) | 0 | 0 | 0 | |

- Criterion 1 fails: raising the trace strength from 0.2 to 4 MPa lowers the wedge trace tension by only 6 to 11 % (30.06 -> 28.31, 84.97 -> 76.79, 110.51 -> 98.82 MPa),
  and the node stays pinned (shear/strength 1.000, tensile and at strength for 95 % of the post-arrival record at 200 and 150 m), 75 to 95 MPa above the insertion run.
- Criterion 2 fails: at 150 m the wedge dp015 slip is 27 % short of the insertion run at the same 4 MPa (38 % short of the 0.2 MPa insertion run); at the trace itself
  the loss is 80 %.
- Criterion 3 fails (no tension, no pinning): lowering tpv36 to 0.2 MPa at the trace leaves the trace compressive after arrival at 500, 250 and 200 m (n max
  -0.15 / -0.12 / -0.11 MPa), never pinned in tension. It does change the rupture (trace arrival 9.06 -> 6.08 s at 500 m, trace slip 0.089 -> 0.481 m,
  1 km slip 0.52 -> 1.17 m at 500 m), so the strength was a live variable there; it just does not produce the TPV12 signature.
- By the pre-registered rule: Test 1 leaves tension (> 5 MPa above insertion) AND slip loss (>= 5 % at 1.5 km), Test 2 shows neither -> **not the cause**.
  Mesh/geometry (60 deg wedge elements at the trace) stays the open suspect; nothing here points at a treatment.

Scope and limits: Test 1 is the amended taper (see Amendment 1), not the uniform-4 MPa variant, which never nucleates (`w12u4_dx500`). The tpv36 4 MPa baseline
is rerun at 500 m only (250/200 m values are the committed wedge-trace numbers). One mesh family per test, term as in the baselines. Diagnostic only: row 149 stays P4
paused; no source, case, test or reference changed. Figures: `fig_test1_tpv12.png`, `fig_test2_tpv36.png`; tables `test1_tpv12.csv`, `test2_tpv36.csv` (every station, every run);
run wall times and ranks in `runs.csv`.
