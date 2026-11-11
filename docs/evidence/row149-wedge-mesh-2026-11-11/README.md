# Row 149: does the trace-node normal stress still rise with refinement on a wedge mesh? (scratch test, 2026-11-11)

Evidence only. No source change, no case converted, no reference moved. Base: master 595e70e.
Solver binary: the row149-force-budget instrumented build (516da2c plus its instrumentation patch); the tracked
difference between 516da2c and 595e70e under src, scripts and case_input is empty, so it is master's solver plus the dump.

## Question

On the insertion mesh (`par.insertFaultType = 1`) the normal stress at the free-surface trace station
`faultst000dp000` of TPV12/13 rises as the grid is refined. Does it still rise when the same 60-degree,
surface-breaking fault is meshed by wedge degeneration (`par.C_degen = 60`, `par.insertFaultType = 0`)?
A measurement, not a treatment.

## Kill criterion (written into the run notes before the first wedge run, not adjusted afterwards)

> The issue is "still there" on the wedge mesh if faultst000dp000 n-stress at 4 s RISES with refinement by more
> than ~1 MPa between the two finest grids (insertion reference: +1.69 MPa for TPV12 200 -> 150 m).

## Verdict

**Still there by the criterion, and it is not the same defect: on the wedge mesh the trace value is far worse.**

| mesh path | dp000 at 4 s, 200 m | 150 m | rise 200 to 150 m |
|---|---|---|---|
| insertion | -3.9865 MPa | -2.2977 MPa | +1.6888 MPa |
| wedge | +84.9702 MPa | +110.5110 MPa | **+25.5408 MPa** |

(FaultMod at 100 m: -3.0050 MPa.) The wedge trace node carries 30 / 85 / 111 MPa of TENSION at 500 / 200 / 150 m,
roughly proportional to 1/dx (value times dx = 15.0 / 17.0 / 16.6 GPa m), while it slips only 1.92 / 1.65 / 1.57 m
against 6.53 / 7.64 / 7.83 m on the insertion mesh. The node one row down is compressive on both meshes.
So the wedge mesh does not remove the rise; it replaces a 1-2 MPa drift by a grid-divergent tensile point value at a
nearly pinned trace node. Away from the trace the two meshes solve the same problem (below), so the dp000 values
of the two paths are **not comparable as the same phenomenon**, and this test does not tell whether the
insertion-mesh rise is a property of the insertion mesh.

## Step 0: code read (paths relative to the repo root)

1. **Wedge path at 60 degrees.** Generic in angle, nothing tied to 15 degrees or to a buried fault.
   - `C_degen` is an integer (`src/fortran/globalvar.f90:107`, read at `readInputFiles.f90:19`); accepted values are
     0 or anything above 3 (`readInputFiles.f90:31`, `scripts/case.setup:44-62`). 60 passes; a non-integer dip cannot.
   - The fault plane is `z = -y tan(dip)` (`library_degeneration.f90:8-15`); `wedge()` (`:17-71`) collapses every brick whose
     centre lies on that plane inside the fault box into a footwall wedge (type 11) and a hanging-wall wedge (type 12);
     `meshgen.f90:164-178` tags them. Wedges are integrated as collapsed one-point hexahedra by the same element
     routine and the same hourglass control (`assembleGlobalKU.f90:26`).
   - The plane passes through grid nodes only if `dz/dy = tan(dip)`. test.tpv36 sets `dy = dx cos(dip)`,
     `dz = dx sin(dip)` (`case_input/test.tpv36/tpv36_37_common.py:35-36`), so along-dip node spacing is dx.
     The uniform y belt must span the fault (`meshgen.f90:655-678`, `nuni_y_plus = round(W/dx)+5`), and the fault
     box needs `fymin, fymax = 0, W cos(dip)` (`scripts/lib.py:106-131`).
   - The on-fault test is distance to the plane below dx/100 inside the fault box (`meshgen.f90:1235-1242`); the
     trace row (z = 0) is included, and test.tpv36 is itself surface-breaking.
2. **Vertical spacing.** `mod4dip` sets `dz = dx sin(dip)` for both paths (`scripts/defaultParameters.py:374-383`), so dz and
   the fault-node positions are the SAME as on the insertion mesh; only dy differs (dx/2 instead of dx).
   The time step had to change: `mod4dip` gives `dt = 0.5 dz/vp`, which exceeds the stability limit of a
   dx by dx/2 by 0.866 dx cell, so the variant uses `dt = 0.5 min(dy, dz)/vp` (0.58 of the insertion step).
3. **What carries over.** Initial stress, friction and the nucleation patch are per fault node from down-dip
   distance (`case_input/test.tpv12/user_defined_params.py`), unchanged because the node grid is unchanged; `C_nuclea = 0`.
   Stations are matched by (x, z) (`meshgen.f90:1380-1381`); fault normal, strike and dip vectors are the same
   (`meshgen.f90:1392-1414`). `C_normal_stress_caps = 0` on both paths (`defaultParameters.py:262`).
4. **TPV13 is blocked on the wedge path without a source change, so it was not run.** `meshgen.f90:201` sets the
   initial element stress (`C_elastic = 0`) for the current element only, after `wedge()` has split the brick in two:
   the footwall wedge would start stress-free. The shipped wedge cases are elastic and never reach this.
   Estimated cost of lifting it: one call plus its mirror in the port plus a test, hours not days. Not done here.

## Runs (fortran, 4 ranks, term 4.5 s, one thread per rank recorded in each run's `run.meta`)

Matched quantity: dx = along-strike node spacing = along-dip node spacing, and dz, so that fault nodes, stations
and on-fault inputs sit at identical positions on both meshes. The only geometric difference is dy.

| run | mesh path | ranks | dx | dy | dz | along-dip | dt (s) | dp000 at 4 s (MPa) | top-500 m mean (MPa) | wall (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| t12_dx500_diag | insertion | 2,2,1 | 500 | 500 | 433.013 | 500 | 0.03788 | -9.2340 | -15.681 | prior evidence |
| t12_dx200 | insertion | 2,2,1 | 200 | 200 | 173.205 | 200 | 0.01515 | -3.9865 | -14.888 | 476 |
| t12_dx150 | insertion | 2,2,1 | 150 | 150 | 129.904 | 150 | 0.01136 | -2.2977 | -14.620 | 1108 |
| w12_dx500 | wedge | 2,2,1 | 500 | 250 | 433.013 | 500 | 0.02187 | +30.0603 | +5.058 | 101 |
| w12_dx500_p212 | wedge | 2,1,2 | 500 | 250 | 433.013 | 500 | 0.02187 | +30.0603 | +5.058 | 88 |
| w12_dx200 | wedge | 2,2,1 | 200 | 100 | 173.205 | 200 | 0.00875 | +84.9702 | +2.874 | 1728 |
| w12_dx150 | wedge | 2,2,1 | 150 | 75 | 129.904 | 150 | 0.00656 | +110.5110 | +2.064 | 4185 |

The (2,1,2) control (the decomposition the shipped wedge cases use) reproduces (2,2,1) to every printed digit, so
the y rank boundary crossing the fault is not the cause. The top-500 m mean (method of the row149-force-budget
evidence) is dominated by the trace node on the wedge mesh and is tensile; with two nodes in 500 m at dx = 500 it
is a two-point average. Full numbers: `wedge_vs_insertion.csv`; profile down dip: `trace_profile.csv`.

## Same-problem check (done before reading dp000)

Rupture-time field, all fault nodes matched on (x, z), wedge against insertion at the same dx
(`rupture_time_field.csv`): median absolute difference 0.282 / 0.108 / 0.079 s at 500 / 200 / 150 m, maximum
0.638 / 0.364 / 0.195 s; the wedge mesh ruptures slightly earlier and the two converge with refinement.

Stations at 150 m (`same_problem.csv`; at 200 m only dp030 and dp120 sit on a fault node, on both meshes):

| station | rupture time (s) wedge / insertion | dip slip at 4 s (m) wedge / insertion | n-stress at 4 s (MPa) wedge / insertion |
|---|---|---|---|
| st000dp015 | 2.001 / 2.034 | 5.94 / 8.38 | -30.58 / -29.53 |
| st000dp030 | 1.745 / 1.773 | 7.95 / 9.19 | -37.23 / -36.33 |
| st000dp075 | 0.925 / 0.955 | 10.31 / 10.27 | -57.84 / -57.96 |
| st045dp075 | 1.411 / 1.489 | 9.61 / 9.50 | -57.47 / -57.56 |
| st000dp120 | 0.013 / 0.023 | 8.65 / 8.49 | -90.98 / -90.62 |
| st120dp000 | 3.162 / 3.193 | 0.31 / 1.00 | +22.14 / -0.41 |

At 7.5 km down dip and deeper the two meshes agree to about 2 percent in slip and 0.4 MPa in normal stress: the
same problem. Toward the surface the wedge mesh loses slip (14 percent at 3 km, 29 percent at 1.5 km, 80 percent
at the trace) and every trace station is tensile. That is gross disagreement confined to the top of the fault,
which is why the dp000 comparison is declared not like for like.

## Stress divergence against hourglass at the trace node (x = 0, 4 s, MPa, per side)

| run | stress div. hanging wall | stress div. footwall | hourglass hanging wall | hourglass footwall |
|---|---|---|---|---|
| insertion 200 m | -6.650 | -8.618 | +2.885 | +4.410 |
| insertion 150 m | -5.479 | -7.140 | +3.303 | +4.720 |
| wedge 500 m | +31.188 | +29.677 | -0.165 | -0.578 |
| wedge 200 m | +85.136 | +88.079 | -0.064 | -3.211 |
| wedge 150 m | +110.595 | +115.385 | -0.046 | -4.911 |

On the insertion mesh the trace value is a small difference of a compressive stress-divergence part and a tensile
hourglass part. On the wedge mesh the hourglass part is small and of the other sign; the tension is in the
stress divergence of the elements beside the trace. The sum of the parts reproduces the station value to the
printed digits on every run (`trace_node_from_parts_MPa`).

## Limits

- TPV12 only; TPV13 blocked as stated. One friction law, one dip, fortran only.
- The wedge variant changes dy and dt as well as the element shape; both follow from the wedge constraint and are
  not separable from it on this mesh generator.
- The 1/dx scaling of the wedge trace value is an observation over three grids, not a diagnosed mechanism; no cause
  was investigated and no fix is proposed here.
- The station rupture time uses a 1 mm/s slip-rate threshold on the station series.

## Files and regeneration

- `user_defined_params.py`: the scratch variant; every changed line is tagged `WEDGE`. Produced from
  `case_input/test.tpv12/user_defined_params.py` by `make_params.py`; `set_case.py` sets dx and the partition per run.
- `run_case.sh`, `queue.sh`, `launch.sh`, `cases_coarse.txt`, `cases_fine.txt`: how the runs were made. They read
  `ROW149_WEDGE_ROOT` (scratch folder holding `params/` and `scripts/`), `ROW149_TREE` (a built tree carrying the
  instrumentation patch of the row149-force-budget evidence) and `ROW149_LOCK`.
- `wedge_ana.py` with arguments insertion_root, insertion_500_root, wedge_root, out_dir regenerates the four CSV
  files and `fig_dp000_vs_dx.png` from the raw run directories, and prints both verdict lines.
