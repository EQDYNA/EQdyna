# TPV32: fine-grid (100 m) acceptance run, fortran, 8 ranks, spec term 15 s, against FaultMod

## Baseline

No own published EQdyna submission exists for TPV32. `scec_archive/` has no tpv32 entry, and a scratch copy of `scripts/scec/fetch_all_dliu.py` (run 2026-10-12, so the tracked `enumeration.json` stayed untouched) found none of the owner's CVWS account ids (including the one the 2014 EQdyna 3.1 runs were filed under) among the 5 users on the TPV32 listing. The listing does hold EQdyna entries, but they are another author's submissions, not the owner's own. The baseline is therefore the finest FaultMod set: `~/shared_dataset/scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014` (FaultMod 25 m, the only FaultMod set, read-only).

## Run

This is a scratch run of `case_input/test.tpv32` with `par.dx` changed from 500 m to 100 m and `par.nx, par.ny, par.nz` changed from 2,1,2 to 4,1,2 (8 ranks, with ny=1 kept). Everything else is unchanged, including the spec term 15 s. It used the fortran backend, built from master `b600274`, on a shared box (load ~42-55/64), running concurrently with the other two cases. Wall time was 5799 s, against a pre-run estimate of ~1.8 h. The estimate scaled the measured tpv26 100 m run (9.5e-7 s per cell-step at 8 ranks) by `case.setup`'s cell count (3.6M) and the step count (1890). This is not a gated cell: no committed reference was read or written, and no run data is in git.

## Numbers (shared nodes only, no interpolation)

| window | dx (m) | ruptured | median \|dt\| | p90 | mean | max | ours-only nodes | one-code-only stations | fault-st dt median | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 100.0% (30578/30578) | 0.013 | 0.032 | +0.013 | 0.05 | 0 | 0  | +0.01 | 0.99 | 1.00 | yes |
| 15 s | 100 | 100.0% (44850/44850) | 0.012 | 0.035 | +0.004 | 0.06 | 0 | 0  | +0.00 | 0.99 | 1.02 | yes |

Station ratios over all shared stations, ours / FaultMod, on [0, 15 s], as median [min-max]: final slip 0.995 [0.980-1.005], peak slip rate 0.917 [0.709-1.017], body peak velocity 1.015 [0.781-1.174]. They come from `../tpv31/station_ratios.py`, stored in `station_ratios_faultmod.json`. Final slip and peak slip rate use |(h, v)|. Body peak velocity is max |v| over three components; a station is skipped when the baseline peak is < 1e-3.

## Verdict (visual, every PNG inspected)

**cplot: MATCH, line-for-line.** Every contour overlays, including the near-surface 1-2 s and the lateral 7 s edges. Max |dt| is 0.06 s over 44850 nodes, with a 100% overlap.

**ts-fault: MATCH.** h-slip, slip rate and stresses overlay at all four sampled stations. One small difference is at surface station faultst060dp000: final v-slip is -0.055 m against FaultMod's -0.058 m (about 5%), and its peak h-slip rate is 2.7 against 3.0 m/s.

**ts-body: MATCH.** Displacements and v/n velocities overlay to line width. h-velocity shows small high-frequency differences after 10 s at the dp000 stations.

**Overall: MATCHES FaultMod 25 m.**

## Reproduce

`python3 -m testsys.frt_canonical RUNS/tpv32_dx100`, then `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS <out> tpv32 --baseline <FaultMod dir> FaultMod faultmod`. That driver is the PR #190 version, and every figure comes from `scripts/figures/scec_compare.py`. Files: `tpv32_cplot_faultmod.png`, `tpv32_ts_fault_faultmod.png` (4 sampled stations), `tpv32_ts_body_faultmod.png` (4 sampled stations), `metrics_faultmod.json`, `table_faultmod.md`.
