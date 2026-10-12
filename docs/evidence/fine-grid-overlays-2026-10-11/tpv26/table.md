### tpv26 -- fine-grid (100 m) acceptance run, fortran, 8 ranks, spec term 13 s

Run: a scratch run of the committed case with `par.dx` changed from 500 to 100. The rank layout `par.nx, par.ny, par.nz` was also changed, from 2,1,2 to 4,1,2 (8 ranks, with ny=1 kept so the fault stays off partition boundaries). This changes the decomposition, not the physics: frt is canonicalised, so the result does not depend on the decomposition. The recipe is otherwise the same as `docs/evidence/barall-overlays-2026-10-09/resolution/README.md`. Binary built from master `3720ccf`, which includes PR #180. Wall time was 7242 s on a shared box (load ~40/64), with tpv27 running concurrently. Baseline: `scec_cvws.tpv2627/raw/tpv26/barall-faultmod-100m-2013` (Barall/FaultMod on a 100 m grid, the same resolution on both sides). No committed reference was read or written, and this is not a gated cell.

| window | dx (m) | median \|dt\| | p90 | mean | ruptured % | slip ratio | body pv ratio | MATCH (numeric) |
|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 0.027 | 0.064 | -0.034 | 100.0% (34535/34535) | 1.01 | 1.00 | yes |
| 13 s | 100 | 0.047 | 0.118 | -0.075 | 100.0% (68460/68465) | 1.00 | 1.03 | yes |

At 125 m (prior evidence) the 13 s figures were median 0.066, p90 0.163 and mean -0.100. Every number moved toward FaultMod at 100 m. Per-station values are in `metrics.json`.

Where the residual sits (13 s window, shared nodes):

| down-dip band | mean dt | median \|dt\| |
|---|---|---|
| 0-2 km | -0.047 | 0.046 |
| 2-15 km | -0.058 | 0.041 |
| 15-20 km | -0.205 | 0.125 |

All 215 nodes with |dt| > 1 s lie on the 17.0-17.2 km bottom arrest edge.

## Verdict against the owner's criterion

**cplot: MATCH, visually very good.** The 2 s and 4 s contours overlap exactly. Our later contours (8-12 s) lead FaultMod slightly on the +x side, by about 0.3-0.8 km. Our 12 s contour also has a small tongue near (x ~ 19 km, 5 km depth) at the lateral fault edge. That difference is under 1 s; none of the |dt| > 1 s nodes are there. The bottom arrest line is ragged in both codes.

**ts-fault: MATCH.** h-slip, h-slip-rate, h-shear and n-normal overlap at all four sampled stations. There are minor differences on the small v-component: v-slip ends about 2 mm (~6%) high at faultst150dp150 and low at faultst100dp100, and the peak h-slip-rate at faultst150dp150 is about 5.3 m/s against 6.4 m/s.

**ts-fault, surface stations (`tpv26_ts_fault_surface.png`): MATCH. The TPV12/13 surface signature is not present.** At faultst-050dp000 and faultst150dp000, n-normal stress stays flat near 0 in both codes, with no early spike. v-shear tracks FaultMod, apart from a late (>11 s) offset of about 0.1 MPa at faultst-050dp000.

**ts-body: MATCH, visually excellent.** All six components at the four shown stations overlap to within line width.

**Overall: MATCHES "really well".** Only small differences remain: a late lateral lead and the bottom arrest edge.

## Files

- `tpv26_cplot.png` -- rupture-time contour overlay.
- `tpv26_ts_fault.png` -- on-fault time series (h/v/n), 4 of 12 shared stations, sampled by the tool.
- `tpv26_ts_fault_surface.png` -- on-fault time series at the two surface (dp000) and two 5 km (dp050) stations, named by `--stations`.
- `tpv26_ts_body.png` -- off-fault time series (h/v/n), 4 of 6 shared stations.
- `metrics.json` -- per-window and per-station numbers.

The cplot/ts-fault/ts-body figures and the numbers were produced by `docs/evidence/barall-overlays-2026-10-09/make_overlays.py resolution`, unmodified, pointed at a one-dx scratch RUNS directory. Their captions say "4 ranks" because that text is hardcoded in the script; the run used 8. The surface figure is a direct `scripts/figures/scec_compare.py --plot ts-fault --stations ...` call. All comparisons go through that one tool.
