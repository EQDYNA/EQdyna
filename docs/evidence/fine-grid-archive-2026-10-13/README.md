# Fine-grid acceptance runs against own published EQdyna: TPV30, 34, 35, 36, 37, 104, 105-3D (2026-10-13)

Owner acceptance criterion: every released TPV gets a fine-grid run that matches either the owner's own published submission or Barall/FaultMod on the same grid, judged visually on cplot and time series. This directory covers the seven cases whose only baseline is an own published EQdyna submission in `scec_archive/` (gitignored, read-only). The baseline for each is the `ARCHIVE` map in `docs/evidence/barall-overlays-2026-10-09/make_overlays.py`: the finest-dx submission, using its `-nstress-corrected` copy where one exists. Template: `docs/evidence/fine-grid-overlays-2026-10-11/`.

## Verdicts

| TPV | baseline (own published EQdyna) | dx | spec term | run? | median \|dt\| (s) | p90 (s) | verdict |
|---|---|---|---|---|---|---|---|
| 104 | EQdyna 4.1 (2016), 100 m, n-stress corrected | 100 m | 12 s | yes, 5161 s | 0.004 | 0.008 | **MATCH**; zip blocked, see below |
| 35 | EQdyna 4.1.3 (2017), 100 m | 100 m | 18 s | yes, 6924 s | 0.003 | 0.008 | **MATCH**; zip built |
| 105-3D | EQdyna 5.1.0 (2020), 100 m | 100 m | 15 s | no, est. 2.2 h | — | — | not run (over 2 h) |
| 30 | EQdyna 3.1 (2015), 25 m | 25 m | 20 s | no, est. 183 h | — | — | not run |
| 34 | EQdyna 4.0 (2016), 25 m | 25 m | 20 s | no, est. 247 h | — | — | not run |
| 36 | EQdyna 5.3.3 (2024), 50 m, n-stress corrected | 50 m | 22 s | no, est. 91 h | — | — | not run |
| 37 | EQdyna 5.3.3 (2024), 50 m, n-stress corrected | 50 m | 22 s | no, est. 91 h | — | — | not run |

Numbers are for the spec-term window, on exactly shared fault nodes. The full tables, which also give the 5 s window, are `<case>/table_*.md`; per-station values are in `<case>/metrics_*.json`.

### TPV104 (rate-state, slip law, velocity-strengthening border)

- **cplot: MATCH, line-for-line.** 100.0% ruptured (62861/62872 in the 12 s window); latest arrival 11.92 s against 11.93 s. Every contour from 2 s to 10 s overlays, including the 16.5-17 km bottom arrest line and the lateral-edge fragments. 13 nodes are ruptured in ours only, and max |dt| is 2.96 s, at the arrest edge.
- **ts-fault: MATCH.** h-slip, slip rate, shear stress and v-components overlay at all 4 stations shown (9 comparable). Final slip ratio is 1.00 at all 9 stations, and arrival agrees to within 0.01 s. The one visible difference is late normal stress, after about 8 s, at faultst-120dp030 and faultst120dp120: the published trace drifts by up to 0.2 MPa and ours stays flatter.
- **ts-body: MATCH.** All 6 components overlay at the 4 stations shown (6 comparable). Peak velocity runs about 7% low in ours at the sharp v-velocity pulse near 7.5 s (median peak ratio 0.93). Arrival is within 0.025 s.
- **Zip: NOT built.** `scripts/scec_submission.py` refuses the run with "field-name tokens with no spec wording: ['psi', 'temperature', 'pressure']". The refusal is correct as written. The TPV103/104 upload spec (`scratch/specs/uploadTPV103.txt`) defines 9 on-fault columns ending in `psi` ("State variable psi (dimensionless)"), but our friclaw-4 station files write 11, adding temperature and pore pressure, which the spec does not ask for. Packaging therefore needs a packager change: a `psi` label, and dropping the two thermal-pressurization columns that are off-spec. That is a `scripts/` change and is not made in this evidence PR.

### TPV35 (data-driven heterogeneous stress, two-sided 1D velocity structure)

- **cplot: MATCH, line-for-line.** 100.0% ruptured (16647/16647 in the 18 s window); latest arrival 10.74 s in both. The 2-10 s contours overlay, including the isolated patch at x 4-9 km. The published 6 s contour adds a few short dashed fragments near x -17 and -10 km at 10 km depth that ours does not draw; they sit at a contour-level crossing (max |dt| 3.19 s, 10 nodes ruptured in ours only).
- **ts-fault: MATCH, indistinguishable at line width.** 7 of 7 stations are comparable (dp081, 8.1 km, which is on the 100 m grid) and 4 are shown. That includes the 2 mm v-slip and the 0.5 MPa n-stress transients. Final slip ratio is 1.00, and the arrival-time median difference is 0.00 s.
- **ts-body: no figure.** The published set uses the spec's named surface receivers (`4064_donna`, ...), while our run writes coordinate-named `body*` files. `scec_compare.py` finds 41 stations in ours and 0 in the baseline, so nothing is paired and the tool exits 1. This is a naming gap in the comparison tool, not a physics gap.
- **Zip:** `scratch/submissions/tpv35_eqdyna_100m.zip` (50 files; author field `<FILL IN>`).

## Runs

- **Setup.** Each run is a scratch run of the committed case made with `scripts/create.newcase`, then three in-place edits to `user_defined_params.py`: `par.dx` 500 to 100, `par.term` 5 to the spec term, and `par.nx` 2 to 4. That gives 8 ranks: (4,2,1) for tpv104 and (4,1,2) for tpv35, each the gate layout doubled in x. dt and the fault grid are computed from `par.dx` downstream of the edited line, so they follow the new dx.
- **Binary and box.** Fortran built from master `6881b50` in a fresh worktree. `OMP_NUM_THREADS=1`. The box was shared, at load 27-32 on 64 cores. The two runs went one at a time.
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

Every case except 105-3D is over 8 h and needs the owner's OK.

## Files

- `tpv104/tpv104_cplot_eqdyna41.png`, `tpv104_ts_fault_eqdyna41.png`, `tpv104_ts_body_eqdyna41.png`, `table_eqdyna41.md`, `metrics_eqdyna41.json`
- `tpv35/tpv35_cplot_eqdyna413.png`, `tpv35_ts_fault_eqdyna413.png`, `table_eqdyna413.md`, `metrics_eqdyna413.json`
