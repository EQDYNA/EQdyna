# TPV31: fine-grid (100 m) acceptance run, fortran, 8 ranks, spec term 15 s, against FaultMod

## Baseline

No own published EQdyna submission exists for TPV31. `scec_archive/` has no tpv31 entry, and a scratch copy of `scripts/scec/fetch_all_dliu.py` (run 2026-10-12, so the tracked `enumeration.json` stayed untouched) found none of the owner's CVWS account ids (including the one the 2014 EQdyna 3.1 runs were filed under) among the 8 users on the TPV31 listing. The listing does hold EQdyna entries, but they are another author's submissions, not the owner's own. The baseline is therefore the finest FaultMod set: `~/shared_dataset/scec_cvws.tpv3132/raw/tpv31/barall.2-faultmod-25m-2014` (FaultMod 25 m, finer than the 50 m FaultMod set, read-only).

## Run

This is a scratch run of `case_input/test.tpv31` with `par.dx` changed from 500 m to 100 m and `par.nx, par.ny, par.nz` changed from 2,1,2 to 4,1,2 (8 ranks, with ny=1 kept). Everything else is unchanged, including the spec term 15 s. It used the fortran backend, built from master `b600274`, on a shared box (load ~42-55/64), running concurrently with the other two cases. Wall time was 5967 s, against a pre-run estimate of ~1.9 h. The estimate scaled the measured tpv26 100 m run (9.5e-7 s per cell-step at 8 ranks) by `case.setup`'s cell count (3.6M) and the step count (1950). This is not a gated cell: no committed reference was read or written, and no run data is in git.

## Numbers (shared nodes only, no interpolation)

| window | dx (m) | ruptured | median \|dt\| | p90 | mean | max | ours-only nodes | one-code-only stations | fault-st dt median | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 100.0% (30804/30804) | 0.013 | 0.031 | +0.010 | 0.05 | 0 | 0  | +0.02 | 1.00 | 0.99 | yes |
| 15 s | 100 | 100.0% (44850/44850) | 0.013 | 0.037 | +0.001 | 0.07 | 0 | 0  | +0.01 | 1.00 | 1.03 | yes |

Station ratios over all shared stations, ours / FaultMod, on [0, 15 s], as median [min-max]: final slip 0.998 [0.981-1.003], peak slip rate 0.956 [0.764-1.049], body peak velocity 1.024 [0.924-1.171]. They come from `../tpv31/station_ratios.py`, stored in `station_ratios_faultmod.json`. Final slip and peak slip rate use |(h, v)|. Body peak velocity is max |v| over three components; a station is skipped when the baseline peak is < 1e-3.

## Verdict (visual, every PNG inspected)

**cplot: MATCH, line-for-line.** All contours from 1 s to 7 s overlay, including the free-surface breakout and the bottom-edge arrival. The only visible offset is about 0.1-0.2 km on the outermost contour at x = ±14-15 km near 1-4 km depth. Max |dt| is 0.07 s over 44850 nodes, with a 100% rupture-extent overlap.

**ts-fault: MATCH.** h-slip, h-slip rate, h-shear and n-normal overlay at all four sampled stations. So do the surface v-slip and v-slip rate at faultst060dp000. Peak slip rate runs slightly low at 100 m (median 0.96, down to 0.76 at one station), as expected against a 25 m grid.

**ts-body: MATCH.** Displacements overlay to line width. The velocities agree except for small high-frequency spikes around 4-5 s at the dp000 stations, where peak ratios range 0.92-1.17.

**Overall: MATCHES FaultMod 25 m.**

## Reproduce

`python3 -m testsys.frt_canonical RUNS/tpv31_dx100`, then `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS <out> tpv31 --baseline <FaultMod dir> FaultMod faultmod`. That driver is the PR #190 version, and every figure comes from `scripts/figures/scec_compare.py`. Files: `tpv31_cplot_faultmod.png`, `tpv31_ts_fault_faultmod.png` (4 sampled stations), `tpv31_ts_body_faultmod.png` (4 sampled stations), `metrics_faultmod.json`, `table_faultmod.md`.
