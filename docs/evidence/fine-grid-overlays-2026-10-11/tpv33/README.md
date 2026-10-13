# TPV33: fine-grid (100 m) acceptance run, fortran, 6 ranks, spec term 13 s, against FaultMod

## Baseline

No own published EQdyna submission exists for TPV33. `scec_archive/` has no tpv33 entry, and a scratch copy of `scripts/scec/fetch_all_dliu.py` (run 2026-10-12, so the tracked `enumeration.json` stayed untouched) found none of the owner's CVWS account ids (including the one the 2014 EQdyna 3.1 runs were filed under) among the 6 users on the TPV33 listing. The listing holds no EQdyna entry. The baseline is therefore the finest FaultMod set: `~/shared_dataset/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015` (FaultMod 12.5 m, an 8x finer grid, read-only).

## Run

This is a scratch run of `case_input/test.tpv33` with `par.dx` changed from 400 m to 100 m and `par.nx, par.ny, par.nz` changed from 2,1,2 to 3,1,2 (6 ranks, with ny=1 kept). Everything else is unchanged, including the spec term 13 s. It used the fortran backend, built from master `b600274`, on a shared box (load ~42-55/64), running concurrently with the other two cases. Wall time was 2609 s, against a pre-run estimate of ~1.0 h. The estimate scaled the measured tpv26 100 m run (9.5e-7 s per cell-step at 8 ranks) by `case.setup`'s cell count (1.76M) and the step count (1560). This is not a gated cell: no committed reference was read or written, and no run data is in git.

## Numbers (shared nodes only, no interpolation)

| window | dx (m) | ruptured | median \|dt\| | p90 | mean | max | ours-only nodes | one-code-only stations | fault-st dt median | slip ratio | body pv ratio | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 99.2% (10472/10553) | 0.049 | 0.110 | +0.005 | 0.86 | 2 | 0  | -0.01 | 1.01 | 0.96 | not yet accepted |
| 13 s | 100 | 98.0% (11481/11713) | 0.052 | 0.114 | -0.002 | 2.66 | 9 | 0  | -0.01 | 1.00 | 0.93 | not yet accepted |

Station ratios over all shared stations, ours / FaultMod, on [0, 13 s], as median [min-max]: final slip 1.000 [0.908-1.012], peak slip rate 0.707 [0.464-0.923], body peak velocity 0.962 [0.492-1.155]. They come from `../tpv31/station_ratios.py`, stored in `station_ratios_faultmod.json`. Final slip and peak slip rate use |(h, v)|. Body peak velocity is max |v| over three components; a station is skipped when the baseline peak is < 1e-3.

## Verdict (visual, every PNG inspected)

**not yet accepted: contour offset up to ~0.7 km, peak slip rate 0.71×; finer run pending.** The owner's acceptance criterion is a fine-grid cplot and time series that visibly match "really well". TPV33 at 100 m does not meet it yet; TPV31 and TPV32 do (line for line). The gaps against FaultMod 12.5 m:

- cplot: the 3 s and 4 s contours sit 0.2-0.4 km inside FaultMod's, and the 5 s contour on +x lags by about 0.7 km. 98.0% of the baseline's ruptured shared nodes ruptured in ours (max |dt| 2.66 s, on the arrest edges).
- ts-fault: peak slip rate is median 0.71 of FaultMod's (as low as 0.46), and v-slip at faultst-100dp020 is about 25% high. Final slip matches (ratio median 1.00).
- ts-body: waveforms overlay through about 8 s; after about 9 s the static displacements at the far stations (body-100st-120dp000, body100st040dp000) drift apart by up to about 0.02 m (not yet explained).

Whether these gaps are resolution (100 m against 12.5 m) or something in the setup is open: a finer (50 m) run and a spec check of the setup are pending.

## Reproduce

`python3 -m testsys.frt_canonical RUNS/tpv33_dx100`, then `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution RUNS <out> tpv33 --baseline <FaultMod dir> FaultMod faultmod`. That driver is the PR #190 version, and every figure comes from `scripts/figures/scec_compare.py`. Files: `tpv33_cplot_faultmod.png`, `tpv33_ts_fault_faultmod.png` (4 sampled stations), `tpv33_ts_body_faultmod.png` (4 sampled stations), `metrics_faultmod.json`, `table_faultmod.md`.

## 50 m run (2026-10-12): match

A 50 m run, `RUNS/tpv33_dx50` (12 ranks, 13 s spec term, fortran, commit `c84a21b`, wall 5279 s), was built with the same driver, `make_overlays.py resolution RUNS <out> tpv33`. Against FaultMod 12.5 m:

| window | dx (m) | ruptured | median \|dt\| | p90 | slip ratio | body pv ratio | verdict |
|---|---|---|---|---|---|---|---|
| 5 s | 50 | 99.9% (42125/42174) | 0.017 | 0.058 | 1.01 | 1.03 | match |
| 13 s | 50 | 99.3% (46501/46850) | 0.019 | 0.065 | 1.01 | 1.00 | match |

**Verdict (visual, every PNG inspected): match at 50 m.**
- The 1-4 s contours lie on FaultMod's line for line. The arrest front on the right side halves its 100 m offset and now sits within ~0.2 km.
- The fault series (slip, shear and normal stress) overlay FaultMod. The peak slip-rate spikes stay below FaultMod's 12.5 m spikes, as expected at 4x the grid spacing.
- The body series overlay through ~10 s. After that FaultMod's static displacement drifts away from both EQdyna runs alike, so the cause is not resolution.

The median rupture-time difference falls from 52 ms to 19 ms when dx halves, so EQdyna converges toward FaultMod. TPV33 is accepted at 50 m, the same resolution TPV27 needed. Files: `tpv33_{cplot,ts_fault,ts_body}_faultmod_50m.png`, `metrics_faultmod_50m.json`, `table_faultmod_50m.md`.
