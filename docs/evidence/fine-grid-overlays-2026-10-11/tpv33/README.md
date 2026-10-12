# TPV33: fine-grid (100 m) acceptance run, fortran, 6 ranks, spec term 13 s, against FaultMod

## Baseline

No own published EQdyna submission exists for TPV33. `scec_archive/` has no tpv33 entry, and a scratch copy of `scripts/scec/fetch_all_dliu.py` (run 2026-10-12, so the tracked `enumeration.json` stayed untouched) found none of the owner's CVWS account ids (including the one the 2014 EQdyna 3.1 runs were filed under) among the 6 users on the TPV33 listing. The listing holds no EQdyna entry. The baseline is therefore the finest FaultMod set: `~/shared_dataset/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015` (FaultMod 12.5 m, an 8x finer grid, read-only).

## Run

This is a scratch run of `case_input/test.tpv33` with `par.dx` changed from 400 m to 100 m and `par.nx, par.ny, par.nz` changed from 2,1,2 to 3,1,2 (6 ranks, with ny=1 kept). Everything else is unchanged, including the spec term 13 s. It used the fortran backend, built from master `b600274`, on a shared box (load ~42-55/64), running concurrently with the other two cases. Wall time was 2609 s, against a pre-run estimate of ~1.0 h. The estimate scaled the measured tpv26 100 m run (9.5e-7 s per cell-step at 8 ranks) by `case.setup`'s cell count (1.76M) and the step count (1560). This is not a gated cell: no committed reference was read or written, and no run data is in git.

## Numbers (shared nodes only, no interpolation)

| window | dx (m) | ruptured | median \|dt\| | p90 | mean | max | ours-only nodes | one-code-only stations | fault-st dt median | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 99.2% (10472/10553) | 0.049 | 0.110 | +0.005 | 0.86 | 2 | 0  | -0.01 | 1.01 | 0.96 | yes |
| 13 s | 100 | 98.0% (11481/11713) | 0.052 | 0.114 | -0.002 | 2.66 | 9 | 0  | -0.01 | 1.00 | 0.93 | yes |

Station ratios over all shared stations, ours / FaultMod, on [0, 13 s], as median [min-max]: final slip 1.000 [0.908-1.012], peak slip rate 0.707 [0.464-0.923], body peak velocity 0.962 [0.492-1.155]. They come from `../tpv31/station_ratios.py`, stored in `station_ratios_faultmod.json`. Final slip and peak slip rate use |(h, v)|. Body peak velocity is max |v| over three components; a station is skipped when the baseline peak is < 1e-3.

## Verdict (visual, every PNG inspected)

**cplot: MATCH, close but not line-for-line.** The 1-3 s contours overlay. On the +x side, our 4 s contour leads FaultMod's by about 0.2-0.3 km. The 2-3 s contours lag up-dip by about 0.2 km (2 s top at 1.9 against 1.65 km depth). The free-surface (1 km) and bottom (9 km) arrest lines and the late 5 s lobe overlay. 98.0% of the baseline's ruptured shared nodes ruptured in ours. The 232 that did not sit on the arrest edges, and max |dt| is 2.66 s.

**ts-fault: MATCH on slip and stress, LOW on peak slip rate.** Final slip and stress drops match (slip ratio median 1.00). Peak slip rate is median 0.71 of FaultMod's, as low as 0.46 (e.g. about 2.3 against 3.8-4.1 m/s at faultst-100dp020 and faultst020dp080). This is the expected under-resolution of the slip-rate peak at 100 m against 12.5 m: the stress drop and slip are unchanged. faultst020dp080's final h-shear is about 0.6 MPa lower than FaultMod's.

**ts-body: MATCH through about 8 s.** The waveforms overlay. After about 9 s, the static h/v/n displacements at the far stations (body-100st-120dp000, body100st040dp000) drift apart by up to about 0.02 m. This is consistent with this case's smaller domain (x -26..18 km, y ±22 km, z 20 km) and absorbing-boundary arrivals, but it has not been measured.

**Overall: MATCHES FaultMod 12.5 m in rupture time, slip and waveforms.** The visible differences are peak slip rate (100 m against 12.5 m) and the late-time static offsets at far body stations.

## Reproduce

`python3 -m testsys.frt_canonical RUNS/tpv33_dx100`, then `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS <out> tpv33 --baseline <FaultMod dir> FaultMod faultmod`. That driver is the PR #190 version, and every figure comes from `scripts/figures/scec_compare.py`. Files: `tpv33_cplot_faultmod.png`, `tpv33_ts_fault_faultmod.png` (4 sampled stations), `tpv33_ts_body_faultmod.png` (4 sampled stations), `metrics_faultmod.json`, `table_faultmod.md`.
