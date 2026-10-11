# FaultMod (SCEC CVWS) overlays for the new TPVs, 2026-10-09

These are cross-code overlays of EQdyna fortran against the FaultMod submissions to the SCEC
Code Verification Web Server (Barall 2009, 2013, 2014, 2015). They cover test.tpv12, 13, 26,
27, 31, 32 and 33, which are the cases gating the next release (PR #158). They also cover the
TPV24/25 branch-fault feature run, which is not gated yet.

The benchmarks have no analytic answer, so agreement here means agreement between codes, not
accuracy. Every figure comes from the one overlay tool, `scripts/figures/scec_compare.py`, using
the reader, plot and CLI layers added in PR #159. That PR must be on the tree before anything
here will rebuild.

## Runs, term, dx

| case | run | dx | term |
|---|---|---|---|
| tpv12, 13, 26, 27, 31, 32 | fortran, 4 ranks, gate config | 500 m | **5 s gate term** (spec term is longer) |
| tpv33 | fortran, 4 ranks, gate config | 400 m | **5 s gate term** |
| tpv24 | fortran, branch-feature case on `board/row153-ckpt2` | 250 m | 12 s |
| tpv25 | fortran, branch-feature case on `board/row153-ckpt2` | 500 m | 12 s |

- **Gated cases.** Each gated case's `frt.canonical.txt` is byte-identical to the committed
  reference in `test.reference.results/test.<case>/`. No reference was regenerated.
- **Refreshed 2026-10-10.** The tpv12, tpv13 and tpv27 figures and numbers were rebuilt
  against the references re-frozen by PR #165 (tpv12 border unpin) and PR #167 (Method-2
  +7.3215 m depth offset removed, tpv13/27/30). See "Refresh 2026-10-10" below.
- **Term.** Each figure states its own term on the figure.
- **Comparison window.** FaultMod's arrivals are windowed to our term. Only nodes that FaultMod
  ruptures before term − 0.25 s are compared.

Rebuild everything with:

```
export SCEC_CVWS=~/shared_dataset        # read-only store; the default
python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py figures RESULTS OUTDIR
python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py numbers RESULTS
```

`RESULTS/<case>/` holds the run's station files and its `frt.canonical.txt`. The per-figure
commands are listed at the end of this file.

## Declared MATCH criteria (fixed before reading the table)

On exactly-shared node coordinates, with nothing interpolated, a case is **MATCH** only if all
of the following hold:

1. At least 95% of the nodes FaultMod ruptures in-window are also ruptured in ours.
2. The median |dt| is at most 0.2 s.
3. The p90 |dt| is at most 0.5 s.
4. No comparable station is ruptured by one code only within the term.

Anything else is **DIVERGES**.

How the station numbers are taken:

- **Which stations.** Only stations that sit exactly on our node grid are used.
- **Arrival.** On the fault, arrival is the first time |h-slip-rate| > 1 mm/s. Off the fault, it
  is the first time |h-vel| > 1 cm/s.
- **Slip.** Slip is read at our term.

## Verdicts

"ruptured" is the share of FaultMod's in-window nodes that ours also ruptures. "mean" is the
mean of (ours − FaultMod). Times are in seconds.

| case | verdict | ruptured | median \|dt\| | p90 | mean | where it diverges |
|---|---|---|---|---|---|---|
| tpv12 | **DIVERGES** | 92.0% (1709/1857) | 0.391 | 0.567 | +0.449 | ours is late; bottom edge. See note 1. |
| tpv13 | **DIVERGES** | 80.0% (1425/1781) | 0.320 | 0.519 | +0.352 | See note 2. |
| tpv26 | **DIVERGES** (core matches) | 100% (1383) | 0.069 | 0.36 | −0.09 | ours leads at the front. See note 3. |
| tpv27 | **DIVERGES** | 100% (1199) | 0.231 | 0.593 | −0.251 | See note 4. |
| tpv31 | **DIVERGES** | 97.4% (1196/1228) | 0.335 | 0.54 | +0.35 | See note 5. |
| tpv32 | **DIVERGES** | 97.9% (1186/1212) | 0.345 | 0.53 | +0.37 | See note 6. |
| tpv33 | **DIVERGES** | 82.5% (539/653) | 1.00 | 1.39 | +1.01 | See note 7. |
| tpv24 main | **MATCH** | 100% (1650) | 0.065 | 0.30 | −0.12 | See note 8. |
| tpv25 main | **DIVERGES** past the junction | 100% (1650) | 0.074 | 0.69 | −0.17 | See note 9. |
| tpv24 branch | **DIVERGES** (geometry) | 95.7% (2161/2257), nearest node | 0.127 | 0.23 | +0.12 | See note 10. |
| tpv25 branch | **DIVERGES** (geometry) | 4/44 | 0.09 | — | — | both codes largely arrest. See note 11. |

Notes:

1. **tpv12.**
   - Rupture: max |dt| is 3.43 s.
   - Of the 148 nodes FaultMod ruptures in-window and ours does not, 114 are at depth ≥ 14 km
     (the bottom edge, FaultMod arrivals from 2.2 s) and the rest at |x| ≥ 12 km.
   - Fault stations: arrivals are 0.27–0.44 s late, and slip is within −23% to +19%.
2. **tpv13.**
   - faultst120dp000: FaultMod ruptures at 3.16 s with −2.19 m of slip. Ours does not rupture by
     5 s.
   - faultst120dp075: slip is 0.09 m in ours vs 0.50 m in FaultMod.
   - body-030st000dp000: ours reaches 0.10 m/s, while FaultMod is about 0.
   - Defect A (the Method-2 lateral-boundary stress imbalance, `resolution/README.md`) is
     still open, and the TPV13 fault border is still pinned (PR #165 unpinned TPV12 only).
3. **tpv26.**
   - 158 nodes are ruptured in ours that FaultMod ruptures only after 5.25 s.
   - faultst100dp100: ours arrives at 4.67 s with 0.62 m of slip, vs 5.54 s in FaultMod.
   - Body stations: arrivals are within 0.06 s, and peak velocity is within 8–18%.
4. **tpv27.**
   - 307 nodes rupture early in ours.
   - st100dp100 is 1.12 s early.
   - Slip at st±050/−150dp100 is 1.65 m vs 1.07 m (+54%).
   - Body peak velocity is 0.137 vs 0.094 m/s (+46%).
5. **tpv31.**
   - st060dp000: FaultMod ruptures at 4.64 s, ours does not.
   - Slip is 10–30% low.
   - Body arrivals are 0.2–0.3 s late.
6. **tpv32.**
   - st060dp000 is not ruptured in ours.
   - st000dp000: slip is 0.22 m vs 0.75 m.
7. **tpv33.**
   - The st000/020 stations are not ruptured in ours.
   - Body peak velocities are 3–6x lower.
8. **tpv24 main.**
   - For x < 0: median 0.019 s, p90 0.072 s.
   - For x ≥ 0 (past the junction): median 0.24 s, and ours is 0.25 s early.
   - Fault stations: arrivals are within 0.34 s, and slip is within 10%.
   - Body peak velocity is within 11%.
9. **tpv25 main.**
   - For x < 0: median 0.032 s.
   - For x ≥ 0: median 0.31 s, p90 1.0 s, and ours is 0.43 s early.
   - faultst090dp100: 5.51 s vs 6.38 s.
10. **tpv24 branch.** Times are measured along the branch from each model's own junction. See
    "Branch geometry" below.
11. **tpv25 branch.** By 11.75 s, FaultMod ruptures 44 of the 589 in-domain nodes and ours
    ruptures 7.

**Summary for the release gate.** None of tpv12/13/26/27/31/32/33 meets MATCH.

- **Plastic cases.** TPV13 and TPV27 diverge.
- **Lagging cases.** TPV12, 31, 32 and 33 lag by 0.3–1.0 s.
- **TPV26** matches in the core at the same 500 m dx.
- **Is the lag a resolution effect?** That is not established by these runs. One data point:
  TPV13's own PR re-measured at dx 250 m, and median |dt| fell from 0.32 to 0.12 s. That points
  to resolution, but the other cases have not been re-run finer.

## Refresh 2026-10-10 (after PR #165 and PR #167)

Fresh fortran runs, 4 ranks, gate dx, 5 s gate term, built from `fe36128`. Each run's
`frt.canonical.txt` is byte-identical to its committed reference:

| case | rows | sha256 |
|---|---|---|
| tpv12 | 1891 | `ae2764ff0181e056d265eceec4dd75cdaa4de7f00fdd945b035eed0a847dc38b` |
| tpv13 | 1891 | `3c31b2514a6369736f685ac8e75d30780d0eedadff875184eefa6703d2f12dc7` |
| tpv27 | 3321 | `bce624bec7c0b381b056508be5c92ceb3765fa72527a2c5523e68cfc8fb6e097` |
| tpv30 | 3321 | `e024a8a3f0b98c9f5a1a2737c45a75dce1498530477661436bd95602f6528625` |

Old (2026-10-09 figures) vs new, gate dx, 5 s window. The gate verdict uses the declared MATCH
criteria above. The class is the resolution study's (`resolution/README.md`). It was not re-run
here, so it is carried over, and amended only where a fix changes what it rests on.

| case | ruptured old → new | median \|dt\| old → new | p90 old → new | gate verdict | class |
|---|---|---|---|---|---|
| tpv12 | 90.5% (1681) → 92.0% (1709) of 1857 | 0.39 → 0.391 | 0.57 → 0.567 | DIVERGES → DIVERGES | DEFECT B → fixed by PR #165 (+28 nodes at gate dx). Timing is RESOLUTION (study). Whether 95% is now reached at finer dx is not re-measured. |
| tpv13 | 80.0% (1425/1781) → unchanged | 0.32 → 0.320 | 0.52 → 0.519 | DIVERGES → DIVERGES | DEFECT → DEFECT. Defect A is open, and the border is still pinned. Defect C removal leaves all four numbers unchanged to 3 decimals. |
| tpv27 | 100% (1199) → 100% (1199) | 0.235 → 0.231 | 0.59 → 0.593 | DIVERGES → DIVERGES | RESOLUTION → RESOLUTION. The "carries Defect C" caveat is cleared by PR #167. |
| tpv30 | — | — | — | not overlaid | No FaultMod (or any CVWS) TPV30 submission is in `~/shared_dataset`, so there is no overlay. None was fetched. |

The TPV13/TPV27 stability matches PR #167's own measurement: the deep rupture-time change was
0.0000 s on both.

## Not overlaid, and why

- **TPV30.** The shared-dataset store holds no CVWS submission for TPV30.

- **Body stations that FaultMod puts elsewhere (TPV12/13).** Eight body stations are not
  overlaid. They sit at different points in the two codes: ours are at -0.3 km / 0.4 km depth,
  FaultMod's at −500 m / 300 m. Figures never interpolate.
- **TPV12 stations at st000.** The 8 comparable body stations are all at st000, and neither code
  shows motion there by 5 s. They are plotted but carry no information.
- **TPV24/25 branch stations.** There is no `fault2` ts-fault or ts-body figure. No branch
  station falls on a shared on-grid node: our along-branch spacing is 288.7 m. Only the branch
  cplot is drawn, and those numbers use the nearest FaultMod node (at most 111 m away).
- **Late-time behaviour.** Anything after the 5 s gate term is not covered for the gated cases.

## Branch geometry (TPV24/25 work-in-progress case, `board/row153-ckpt2`)

`case_input/test.tpv24/tpv24_25_common.py` hardcodes `BRANCH_FXMIN = 1000.0`, while the wedge
line starts at `fxmin − dx`. The result is that the measured junction moves with dx:

| dx | junction |
|---|---|
| 1000 m | x = 0 |
| 500 m | x = 500 m |
| 250 m | x = 750 m |

The spec junction is x = 0.

The branch is also decimated to a 10 km x-extent. Its along-branch length is 10.68 km at dx 250
and 10.97 km at dx 500, against the spec's 12 km.

The branch cplots carry this warning on the figure. The fault-2 numbers above are therefore not
a statement about the branch physics.

## Per-figure commands

Run from the repo root. `$RESULTS` and `$SCEC_CVWS` are as above.

- `tpv12_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv12 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv12/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv12_cplot.png
  ```

- `tpv12_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv12 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv12/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv12_ts_body.png
  ```

- `tpv12_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv12 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv12/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv12_ts_fault.png
  ```

- `tpv13_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv13 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv13/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv13_cplot.png
  ```

- `tpv13_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv13 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv13/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv13_ts_body.png
  ```

- `tpv13_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv13 $SCEC_CVWS/scec_cvws.tpv1213/raw/tpv13/barall-faultmod-100m-2009 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv13_ts_fault.png
  ```

- `tpv24_fault1_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --fault 1 --models $RESULTS/tpv24 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv24/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 250 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv24_fault1_cplot.png
  ```

- `tpv24_fault1_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --fault 1 --models $RESULTS/tpv24 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv24/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 250 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv24_fault1_ts_body.png
  ```

- `tpv24_fault1_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --fault 1 --models $RESULTS/tpv24 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv24/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 250 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv24_fault1_ts_fault.png
  ```

- `tpv24_fault2_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --fault 2 --models $RESULTS/tpv24 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv24/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 250 m, run to 12 s (branch-feature run, not a gated case yet); NOTE this run puts the branch junction at x = 750 m, not the spec x = 0 (WIP case: branch fxmin fixed at 1000 m while the wedge line starts at fxmin - dx), and the branch is decimated to a 10 km x-extent; along-strike is measured from this run's junction' --out tpv24_fault2_cplot.png
  ```

- `tpv25_fault1_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --fault 1 --models $RESULTS/tpv25 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv25/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv25_fault1_cplot.png
  ```

- `tpv25_fault1_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --fault 1 --models $RESULTS/tpv25 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv25/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv25_fault1_ts_body.png
  ```

- `tpv25_fault1_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --fault 1 --models $RESULTS/tpv25 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv25/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 12 s (branch-feature run, not a gated case yet)' --out tpv25_fault1_ts_fault.png
  ```

- `tpv25_fault2_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --fault 2 --models $RESULTS/tpv25 $SCEC_CVWS/scec_cvws.tpv2425/raw/tpv25/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 12 s (branch-feature run, not a gated case yet); NOTE this run puts the branch junction at x = 500 m, not the spec x = 0 (WIP case: branch fxmin fixed at 1000 m while the wedge line starts at fxmin - dx), and the branch is decimated to a 10 km x-extent; along-strike is measured from this run's junction' --out tpv25_fault2_cplot.png
  ```

- `tpv26_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv26 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv26/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv26_cplot.png
  ```

- `tpv26_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv26 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv26/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv26_ts_body.png
  ```

- `tpv26_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv26 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv26/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv26_ts_fault.png
  ```

- `tpv27_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv27 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv27/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv27_cplot.png
  ```

- `tpv27_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv27 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv27/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv27_ts_body.png
  ```

- `tpv27_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv27 $SCEC_CVWS/scec_cvws.tpv2627/raw/tpv27/barall-faultmod-100m-2013 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv27_ts_fault.png
  ```

- `tpv31_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv31 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv31/barall-faultmod-50m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv31_cplot.png
  ```

- `tpv31_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv31 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv31/barall-faultmod-50m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv31_ts_body.png
  ```

- `tpv31_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv31 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv31/barall-faultmod-50m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv31_ts_fault.png
  ```

- `tpv32_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv32 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv32_cplot.png
  ```

- `tpv32_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv32 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv32_ts_body.png
  ```

- `tpv32_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv32 $SCEC_CVWS/scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014 --dpi 200 --note 'EQdyna fortran, gate config dx 500 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv32_ts_fault.png
  ```

- `tpv33_cplot.png`

  ```
  python3 scripts/figures/scec_compare.py --plot cplot --models $RESULTS/tpv33 $SCEC_CVWS/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015 --dpi 200 --note 'EQdyna fortran, gate config dx 400 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv33_cplot.png
  ```

- `tpv33_ts_body.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-body --component h v n --models $RESULTS/tpv33 $SCEC_CVWS/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015 --dpi 200 --note 'EQdyna fortran, gate config dx 400 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv33_ts_body.png
  ```

- `tpv33_ts_fault.png`

  ```
  python3 scripts/figures/scec_compare.py --plot ts-fault --component h v n --models $RESULTS/tpv33 $SCEC_CVWS/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015 --dpi 200 --note 'EQdyna fortran, gate config dx 400 m, run to 5 s (5 s gate term; spec term is longer)' --out tpv33_ts_fault.png
  ```

