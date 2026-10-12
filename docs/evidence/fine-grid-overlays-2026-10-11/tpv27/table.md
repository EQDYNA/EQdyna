### tpv27 -- fine-grid (100 m) acceptance run, fortran, 8 ranks, spec term 13 s

Run: a scratch run of the committed case with `par.dx` changed from 500 to 100. The rank layout `par.nx, par.ny, par.nz` was also changed, from 2,1,2 to 4,1,2 (8 ranks, with ny=1 kept). This changes the decomposition, not the physics. The recipe is otherwise the same as `docs/evidence/barall-overlays-2026-10-09/resolution/README.md`. Binary built from master `3720ccf`. Wall time was 7391 s on a shared box (load ~40/64), with tpv26 running concurrently. Baseline: `scec_cvws.tpv2627/raw/tpv27/barall-faultmod-100m-2013` (Barall/FaultMod on a 100 m grid). No committed reference was read or written, and this is not a gated cell.

| window | dx (m) | median \|dt\| | p90 | mean | ruptured % | slip ratio | body pv ratio | MATCH (numeric) |
|---|---|---|---|---|---|---|---|---|
| 5 s | 100 | 0.089 | 0.157 | -0.090 | 100.0% (29808/29808) | 1.04 | 1.06 | yes |
| 13 s | 100 | 0.134 | 0.292 | -0.184 | 100.0% (62174/62174) | 1.00 | 1.05 | yes |

At 125 m (prior evidence) the 13 s figures were median 0.184, p90 0.360 and mean -0.245. The trend is toward FaultMod, but slow. EQdyna stays systematically early. Per-station values are in `metrics.json`.

Where the residual sits (13 s window, shared nodes):

| down-dip band | mean dt | median \|dt\| |
|---|---|---|
| 0-2 km | **-1.098** | **0.994** |
| 2-15 km | -0.133 | 0.122 |
| 15-20 km | -0.317 | 0.236 |

1100 nodes have |dt| > 1 s, mostly in the shallow band.

## Verdict against the owner's criterion

**cplot: NOT a "really well" match.** The 2 s contour overlaps. From 4 s on, EQdyna's front leads FaultMod's everywhere: by about 0.3-0.6 km laterally on the 6/8/10 s contours, and much more toward the free surface. Our 10 s contour reaches about 0.8 km depth; FaultMod's stops near 1.5 km. Neither code ruptures the top ~1 km within 13 s. In the shallow 0-2 km band, our rupture arrives about 1.1 s early on average. The plasticity case is where the two codes differ; the elastic twin tpv26 matches (see `../tpv26/table.md`). That points at the shallow off-fault plastic response (spec cohesion 1.36 MPa, Tv = 0.03 s, low confining stress near the surface), not at the fault solver. No cause was diagnosed here, and no code was changed.

**ts-fault: MATCH with visible minor offsets.** h-slip, h-slip-rate, h-shear and n-normal overlap at the four sampled stations. Our arrivals are slightly early. v-slip ends about 3 mm (~12%) high at faultst150dp150 and about 2 mm low at faultst100dp100.

**ts-fault, surface stations (`tpv27_ts_fault_surface.png`): the TPV12/13 surface signature is not present.** n-normal stress at faultst-050dp000 and faultst150dp000 stays near 0 in both codes, with no early spike. Neither code slips at these stations by 13 s, so the two agree. h-shear tracks. v-shear diverges late (>11 s) at faultst150dp000: about 0.25 MPa in ours against 0.13 MPa in FaultMod. This is small, and the timing fits the early shallow front.

**ts-body: MATCH, visually very good.** All six components at the four shown stations overlap. Our phases arrive slightly early, consistent with the rupture lead.

**Overall: NOT a clean "match really well" by the owner's bar.** Slip, stations and body motion are good. The rupture-time field is systematically early, by about 1 s in the shallow 2 km. It converges with dx, but slowly (the mean went from -0.245 at 125 m to -0.184 at 100 m). Recommend holding tpv27 sign-off until the shallow lead is explained: either a resolution effect, which a 50 m point would show, or a plasticity defect.

## Files

- `tpv27_cplot.png` -- rupture-time contour overlay.
- `tpv27_ts_fault.png` -- on-fault time series (h/v/n), 4 of 12 shared stations, sampled by the tool.
- `tpv27_ts_fault_surface.png` -- on-fault time series at the two surface (dp000) and two 5 km (dp050) stations, named by `--stations`.
- `tpv27_ts_body.png` -- off-fault time series (h/v/n), 4 of 6 shared stations.
- `metrics.json` -- per-window and per-station numbers.

The figures and numbers were produced the same way as tpv26: `make_overlays.py resolution`, unmodified. Its captions say "4 ranks" because that text is hardcoded; the run used 8. The surface figure is a direct `scripts/figures/scec_compare.py` call.
