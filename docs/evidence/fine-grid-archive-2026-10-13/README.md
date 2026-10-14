# Fine-grid acceptance runs against own published EQdyna: TPV30, 34, 35, 36, 37, 104, 105-3D (2026-10-13)

Owner acceptance criterion: every released TPV gets a fine-grid run that matches either the owner's own published submission or Barall/FaultMod on the same grid, judged visually on cplot and time series. This directory covers the seven cases whose only baseline is an own published EQdyna submission in `scec_archive/` (gitignored, read-only). The baseline for each is the `ARCHIVE` map in `docs/evidence/barall-overlays-2026-10-09/make_overlays.py`: the finest-dx submission, using its `-nstress-corrected` copy where one exists. Template: `docs/evidence/fine-grid-overlays-2026-10-11/`.

## Verdicts

| TPV | baseline (own published EQdyna) | dx | spec term | run? | median \|dt\| (s) | p90 (s) | verdict |
|---|---|---|---|---|---|---|---|
| 104 | EQdyna 4.1 (2016), 100 m, n-stress corrected | 100 m | 12 s | yes, 5161 s | 0.004 | 0.008 | **MATCH**; zip built (after PR #198) |
| 35 | EQdyna 4.1.3 (2017), 100 m | 100 m | 18 s | yes, 6924 s | 0.003 | 0.008 | **MATCH** on cplot and ts-fault; ts-body **partial** (see below); zip built |
| 105-3D | EQdyna 5.1.0 (2020), 100 m | 100 m | 15 s | yes, 16587 s | 0.005 | 0.011 | **MATCH**; zip built |
| 30 | EQdyna 3.1 (2015), 25 m | 25 m | 20 s | no, est. 183 h | — | — | not run |
| 34 | EQdyna 4.0 (2016), 25 m | 25 m | 20 s | no, est. 247 h | — | — | not run |
| 36 | EQdyna 5.3.3 (2024), 50 m, n-stress corrected | 50 m | 22 s | no, est. 91 h | — | — | not run |
| 37 | EQdyna 5.3.3 (2024), 50 m, n-stress corrected | 50 m | 22 s | no, est. 91 h | — | — | not run |

Numbers are for the spec-term window, on exactly shared fault nodes. The full tables, which also give the 5 s window, are `<case>/table_*.md`; per-station values are in `<case>/metrics_*.json`.

### TPV104 (rate-state, slip law, velocity-strengthening border)

- **cplot: MATCH, line-for-line.** 100.0% ruptured (62861/62872 in the 12 s window); latest arrival 11.92 s against 11.93 s. Every contour from 2 s to 10 s overlays, including the 16.5-17 km bottom arrest line and the lateral-edge fragments. 13 nodes are ruptured in ours only, and max |dt| is 2.96 s, at the arrest edge.
- **ts-fault: MATCH.** h-slip, slip rate, shear stress and v-components overlay at all 4 stations shown (9 comparable). Final slip ratio is 1.00 at all 9 stations, and arrival agrees to within 0.01 s. The one visible difference is late normal stress, after about 8 s, at faultst-120dp030 and faultst120dp120: the published trace drifts by up to 0.2 MPa and ours stays flatter.
- **ts-body: MATCH.** All 6 components overlay at the 4 stations shown (6 comparable). Peak velocity runs about 7% low in ours at the sharp v-velocity pulse near 7.5 s (median peak ratio 0.93). Arrival is within 0.025 s.
- **Zip: built** (`scratch/submissions/tpv104_eqdyna_100m.zip`, 21 files, author field `<FILL IN>`) after PR #198. The packager now emits exactly the TPV103/104 upload spec's 9 on-fault columns, with psi taken from our column 9. That column agrees with the own published psi to within 0.6% at 3 stations and 4 times. The two thermal-pressurization columns are dropped, and the drop is recorded in PROVENANCE.md. The zip carries all 13 of our on-fault stations; the spec and the own published set list 9.

### TPV35 (data-driven heterogeneous stress, two-sided 1D velocity structure)

- **cplot: MATCH, line-for-line.** 100.0% ruptured (16647/16647 in the 18 s window); latest arrival 10.74 s in both. The 2-10 s contours overlay, including the isolated patch at x 4-9 km. The published 6 s contour adds a few short dashed fragments near x -17 and -10 km at 10 km depth that ours does not draw; they sit at a contour-level crossing (max |dt| 3.19 s, 10 nodes ruptured in ours only).
- **ts-fault: MATCH, indistinguishable at line width.** 7 of 7 stations are comparable (dp081, 8.1 km, which is on the 100 m grid) and 4 are shown. That includes the 2 mm v-slip and the 0.5 MPa n-stress transients. Final slip ratio is 1.00, and the arrival-time median difference is 0.00 s.
- **ts-body: PARTIAL match** (`tpv35_ts_body_eqdyna413.png`). The figure needs PR #199, which pairs a named receiver (`4064_donna`, ...) to our coordinate-named `body*` file. The pairing goes through the official list `case_input/test.tpv35/tpv35_station_locations.txt`, rounded to our 100 m names. 42 of 42 published receivers pair, and 41 are common with our run; 4072_redh has no output file in ours, and that is not investigated. Over the 18 s window, the arrival \|dt\| median is 0.003 s and the max 0.12 s. The peak-velocity ratio (ours / published) has median 0.99 and range 0.43-1.06.
  - All 6 components overlay at the 4 stations shown. The visible gap is the sharpest h- and v-velocity peaks, where ours runs 15-25% low; at body080st046, for example, h-velocity peaks at 0.11 against 0.13 m/s.
  - 8 of 41 stations are more than 10% low on peak velocity: body059st213 0.43, body-108st-021 0.77, body062st-022 0.79, body041st187 0.80, body055st059 0.83, body080st046 0.83, body051st113 0.84, body026st118 0.87.
  - **body059st213 does not match.** It sits at x 21.3 km, near the strike end. Its h-displacement before 9 s is about half the published one, and the late static offset has the opposite sign. Only its h-component and that of body062st-022 were viewed; the other 6 low stations were not viewed. No tuning was attempted.
- **Zip:** `scratch/submissions/tpv35_eqdyna_100m.zip` (50 files; author field `<FILL IN>`).

### TPV105-3D (rate-state with thermal pressurization)

- **cplot: MATCH, line-for-line.** 100.0% ruptured (68059/68071 in the 15 s window); latest arrival 14.99 s against 14.93 s. The 2-14 s contours overlay. The only differences are short dashed 10 s fragments along the 17 km arrest line, near x -15 and +15 km, and near x 19.5 km at 12 km depth. Max \|dt\| is 3.02 s, and 1 node is ruptured in ours only.
- **ts-fault: MATCH.** 13 of 13 stations are comparable and 4 are shown. All 7 quantities overlay at line width, including the normal-stress drop to 28 MPa and its recovery. Final slip ratio is 1.00 and the arrival median difference is +0.01 s.
- **ts-body: MATCH.** 6 of 6 stations are comparable and 4 are shown; all 6 components overlay. The median peak-velocity ratio is 1.01.
- **Zip:** `scratch/submissions/tpv1053d_eqdyna_100m.zip` (21 files: 13 on-fault, 6 off-fault, cplot, PROVENANCE.md, the same station set as the own published one; author field `<FILL IN>`). This needs the PR that lets the packager keep TPV105-3D's 11 spec columns, temperature and pore pressure included; before it, the packager refused them as having "no spec wording". Commit `b95a986`, verified against the run banner stamp `1d351422917f`.

## Runs

- **Setup.** Each run is a scratch run of the committed case made with `scripts/create.newcase`, then three in-place edits to `user_defined_params.py`: `par.dx` 500 to 100, `par.term` 5 to the spec term, and `par.nx` 2 to 4. That gives 8 ranks: (4,2,1) for tpv104 and (4,1,2) for tpv35, each the gate layout doubled in x. dt and the fault grid are computed from `par.dx` downstream of the edited line, so they follow the new dx.
- **TPV105-3D setup.** This case has no `par.dx` line of its own; it inherits every parameter from `scripts/defaultParameters.py`. So the run directory's own copy of that file had `dx` set 500 to 100, which makes dt (0.00833 s), the fault grid and the initial-stress loop all follow. `par.term` was set to 15 and `par.nx` to 4, giving (4,2,1) on 8 ranks and 29.8 M cells. Wall time was 16587 s (4.6 h against the 2.2 h estimate), at box load 50-75 on 64 cores from other users' jobs.
- **Binary and box.** Fortran built from master `6881b50` in a fresh worktree (TPV105-3D: master `b95a986`, same source stamp `1d351422917f`). `OMP_NUM_THREADS=1`. The box was shared, at load 27-32 on 64 cores. The two runs went one at a time.
- **Scope.** These are not gated cells. No committed reference was read, written or regenerated (rule 7).
- **Figures and numbers.** `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS <out> <case> --baseline <scec_archive dir> "own published EQdyna <ver> (<year>)" <tag>`, which calls `scripts/figures/scec_compare.py`.
- **Script change in this PR.** `SPEC_TERM` gains the seven archive cases' spec terms. For archive cases, `GATE_DX` is `None`, so a single-dx fine run is scored on its own node grid; previously it would be filtered to a 500 m grid it never shares. FaultMod cases are unchanged.

## Wall-time estimate (cases not run)

No gate run keeps per-step timings, only `wall_s` in `docs/evidence/sweep-*/summary.json`. So the estimate comes from the gate run's wall time (sweep `d79c8cc`, 4 ranks, 5 s, dx 500 m):

    est = gate_wall x (spec_term / 5) x (500 / dx)^3 x 0.5   (8 ranks against 4)

- **Why the exponent is 3, not 4.** The mesh coarsens in the fault-normal direction (`enlarging_ratio`), and the measured spec-term walls of tpv26 at 250 to 125 m and tpv33 at 200 to 100 m grew 7.7x and 6.4x.
- **Check on earlier runs.** The model predicts 8062 s for the measured tpv26 run at 100 m on 8 ranks, which took 7242 s.
- **Check on today's runs.** Predicted against measured: tpv104 1.42 h against 1.43 h; tpv35 1.59 h against 1.92 h.

| TPV | gate wall (s) | dx | spec term | est. at 8 ranks |
|---|---|---|---|---|
| 105-3D | 41.7 | 100 m | 15 s | 2.2 h |
| 36 | 149.1 | 50 m | 22 s | 91 h |
| 37 | 149.1 | 50 m | 22 s | 91 h |
| 30 | 41.2 | 25 m | 20 s | 183 h |
| 34 | 55.6 | 25 m | 20 s | 247 h |

Additional blockers beyond cost:

- **TPV30 at 25 m** needs a 25 m rough-fault geometry file. `case_input/test.tpv30/` ships only the 100 m one, and `lib.requireFaultGeometryResolution` refuses any finer dx.
- **TPV34 at 25 m** needs a fresh CVM-H extraction; the shipped grid admits only 500 m (`testsys/e2e/full_specs.py`, the EXCLUDED entry).

Every case not run is over 8 h and needs the owner's OK.

## Files

- `tpv104/tpv104_cplot_eqdyna41.png`, `tpv104_ts_fault_eqdyna41.png`, `tpv104_ts_body_eqdyna41.png`, `table_eqdyna41.md`, `metrics_eqdyna41.json`
- `tpv35/tpv35_cplot_eqdyna413.png`, `tpv35_ts_fault_eqdyna413.png`, `tpv35_ts_body_eqdyna413.png`, `table_eqdyna413.md`, `metrics_eqdyna413.json`; the table and metrics now carry the body peak-velocity ratio
- `tpv1053d/tpv1053d_cplot_eqdyna510.png`, `tpv1053d_ts_fault_eqdyna510.png`, `tpv1053d_ts_body_eqdyna510.png`, `table_eqdyna510.md`, `metrics_eqdyna510.json`
