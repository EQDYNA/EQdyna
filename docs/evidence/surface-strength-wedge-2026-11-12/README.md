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

## Verdict

PENDING (runs not started when this file was committed).
