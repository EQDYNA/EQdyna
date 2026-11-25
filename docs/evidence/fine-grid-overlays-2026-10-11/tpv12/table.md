### tpv12 -- fine-grid (100 m) acceptance run, fortran, 4 ranks, spec term 8 s

Run: scratch case (`par.dx` changed 500 -> 100 only, same recipe as
`docs/evidence/barall-overlays-2026-10-09/resolution/README.md`). Baseline:
`scec_cvws.tpv1213/raw/tpv12/barall-faultmod-100m-2009` (Barall/FaultMod, 100 m
grid -- same resolution on both sides). No committed reference read or
written; this is not a gated cell.

| dx (m) | median \|dt\| | p90 | ruptured % | slip ratio | body pv ratio | MATCH |
|---|---|---|---|---|---|---|
| 100 | 0.019 | 0.040 | 100.0% (45005/45005) | 1.01 | 1.04 | yes (rupture time, slip, body motion) |

(8 s spec-term window; the 5 s window is effectively identical -- see
`metrics.json` in this directory.)

## Verdict against the owner's criterion

**cplot: MATCH, visually excellent.** Rupture-time contours overlap EQdyna vs
FaultMod almost exactly everywhere, including the bottom edge and the lateral
(|x|=15 km) border that previously failed (Defect B, fixed PR #165) and the
PML-typing fix (PR #180, landed same day as this run, not previously measured
at this dx).

**ts-body: MATCH, visually excellent.** All six components (h/v/n
displacement and velocity) at all four shown stations overlap the baseline to
within line width.

**ts-fault: PARTIAL.** h-slip, h-slip-rate, h-shear-stress, v-slip and
v-slip-rate all overlap the baseline closely. **v-shear-stress and n-normal
stress at the three shown surface (dp000) stations do not match**: ours shows
an early transient spike (n-normal stress up to ~13-19 MPa at
faultst000dp000/faultst045dp000, ~4 MPa at faultst120dp000) that FaultMod does
not show (FaultMod stays near 0, briefly negative). This is the same
signature reported for TPV13 at 100 m (on-fault n-normal/v-shear divergence
at surface stations, ~+10 vs ~2 MPa) that a parallel job is diagnosing on
`fix/tpv13-surface-hourglass` (not touched here). **TPV12 shows the same
divergence**, same component (n-normal stress, v-shear stress), same station
depth (dp000, the free surface).

**Overall: NOT a clean "match really well" by the owner's bar.** Rupture
timing, slip and off-fault motion are excellent; the fault-surface normal/
shear stress channel is not, and shares a cause with the open TPV13 surface
investigation. Recommend holding TPV12 submission until that fix lands and
closes this gap too (or an owner call that cplot + slip/timing suffice for
sign-off regardless).

## Files

- `tpv12_cplot.png` -- rupture-time contour overlay.
- `tpv12_ts_fault.png` -- on-fault station time series (h/v/n), 3 of 3
  exactly-shared stations shown.
- `tpv12_ts_body.png` -- off-fault station time series (h/v/n), 4 of 8
  exactly-shared stations shown.

Produced by `docs/evidence/barall-overlays-2026-10-09/make_overlays.py
resolution` (unmodified) pointed at a one-dx scratch RUNS directory; figures
are `scripts/figures/scec_compare.py` only (rule: comparisons go through that
tool, no per-TPV script).
