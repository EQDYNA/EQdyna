# EQdyna-archive overlays, 2026-10-10

Owner decision 2026-10-10 (board PR #175): for every TPV with the owner's archived SCEC
submission in `scec_archive/`, the overlay baseline is that submission (**EQdyna-archive**)
instead of FaultMod. TPVs without one keep FaultMod (`../barall-overlays-2026-10-09/`).

This is a **cross-version consistency check** of EQdyna against its own earlier published
results, not an independent validation: agreement certifies nothing against the benchmark.

## Runs

- EQdyna fortran, 4 ranks, gate config (dx 500 m), **5 s gate term**, master `adf0c1c`.
  The archive submissions are run at 25-100 m to the full spec term; only their arrivals
  before term - 0.25 s are compared.
- Each run passed its e2e cell (8/8, 279 s wall for all eight) and its `frt.canonical.txt`
  is **byte-identical** to `test.reference.results/test.<case>/frt.canonical.txt`
  (line counts 1891-4005). No reference was regenerated.
- In scope: every registered case with an archived submission. `tpv22`/`tpv23` have none
  (only third-party submissions), so they stay on FaultMod and are not here.

Baseline per case (`ARCHIVE` in `make_overlays.py`): the finest-resolution submission in the
benchmark's main listing, its `-nstress-corrected` copy where one exists.

## Verdicts

Same declared MATCH criteria as `../barall-overlays-2026-10-09/README.md` (>= 95% of baseline
in-window nodes ruptured, median |dt| <= 0.2 s, p90 <= 0.5 s, no station ruptured by one
code only), on exactly-shared node coordinates, nothing interpolated, station arrivals on the
vector rate. "mean" is ours - baseline; positive means ours is late.

| case | baseline (EQdyna-archive) | verdict | ruptured | median \|dt\| (s) | p90 (s) | mean (s) | one-code-only stations |
|---|---|---|---|---|---|---|---|
| tpv30 | v3.1 25 m (2015) | **DIVERGES** | 100.0% (1298/1298) | 0.221 | 0.537 | -0.226 | 1 (faultst-150dp050) |
| tpv29 | v3.1 50 m (2015) | **MATCH** | 100.0% (1449/1449) | 0.070 | 0.254 | -0.061 | 0 |
| tpv34 | v4.0 25 m (2016) | **DIVERGES** | 96.4% (1158/1201) | 0.322 | 0.473 | +0.324 | 3 |
| tpv35 | v4.1.3 100 m (2017) | **DIVERGES** | 92.5% (396/428) | 0.341 | 0.474 | +0.289 | — (no comparable station) |
| tpv36 | v5.3.3 50 m (2024), nstress-corrected | **DIVERGES** | 93.8% (1183/1261) | 0.138 | 0.439 | +0.131 | 0 |
| tpv37 | v5.3.3 50 m (2024), nstress-corrected | **DIVERGES** | 93.8% (1181/1259) | 0.138 | 0.438 | +0.130 | 0 |
| tpv104 | v4.1 100 m (2016), nstress-corrected | **DIVERGES** | 68.4% (733/1072) | 0.391 | 0.861 | +0.459 | 2 |
| tpv1053d | v5.1.0 100 m (2020) | **DIVERGES** | 72.7% (893/1229) | 0.449 | 0.830 | +0.472 | 2 |

Full table (station dt, slip and peak-velocity ratios): `table.md`; per-station numbers:
`metrics.json`.

Reading it:

- **tpv29 matches** its own 2015 50 m submission at the gate dx.
- **tpv30** (first overlay it has ever had): every in-window archive node is ruptured; ours
  leads by 0.23 s on average with the front outrunning the 25 m archive along the -x edge.
- **tpv104 / tpv1053d** (rate-state, friclaw 4/5): ours is uniformly late (+0.46 / +0.47 s)
  and ruptures ~70% of the archive's in-window area -- the concentric contours in
  `tpv104_cplot.png` keep their shape and lag. Consistent with an under-resolved
  rate-state process zone at 500 m vs the archive's 100 m; this run does not test that cause.
- **tpv34, tpv35, tpv36, tpv37**: late by 0.13-0.32 s on average. tpv36/37 fail on the
  ruptured fraction only (93.8%), tpv35 on fraction and median, and tpv34 on median and three
  one-code-only stations.

No resolution study was run here; "DIVERGES" at the gate dx says nothing about convergence.

## Figures

`<case>_cplot.png`, `<case>_ts_fault.png`, `<case>_ts_body.png` for each case, except
**tpv35: cplot only**. Its archive's on-fault stations are at 8.1 km down-dip (`dp081`), off
the 500 m grid, and its off-fault receivers are named sites, not a `body*` grid, so no station
is exactly shared. They are not snapped or interpolated.

## Rebuild

```
python3 testsys/e2e/run_e2e.py --cases test.tpv30,test.tpv29,test.tpv34,test.tpv35,test.tpv36,test.tpv37,test.tpv104,test.tpv1053d --backends fortran --jobs 8
python3 -m testsys.frt_canonical test/test.<case>     # then cmp against the committed reference
# RESULTS/<case-without-test.> -> test/test.<case>
python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py figures RESULTS OUTDIR tpv30 tpv29 tpv34 tpv35 tpv36 tpv37 tpv104 tpv1053d
python3 docs/evidence/barall-overlays-2026-10-09/make_overlays.py table   RESULTS OUTDIR tpv30 tpv29 tpv34 tpv35 tpv36 tpv37 tpv104 tpv1053d
```

`scec_archive/` is gitignored and read-only (`$SCEC_ARCHIVE`, default: the main checkout's).
Figures come only from `scripts/figures/scec_compare.py`, which needed no change: every archive
layout here was already supported.
