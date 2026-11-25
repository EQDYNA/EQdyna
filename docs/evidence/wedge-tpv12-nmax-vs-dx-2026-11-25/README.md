# TPV12 trace-node transient normal-stress spike versus dx, wedge and insertion meshes (row 167 research) -- 2026-11-25

Evidence only; same-script re-run of the row 167 wedge case (`docs/evidence/row149-wedge-mesh-2026-11-11/run_case.sh`), master src
(`a954df6`-equivalent binaries of the `test/row169-bk` build), 4 ranks, threads 1. Wedge: `par.C_degen = 60`, `insertFaultType = 0`, term 4.5 s.
Control: `test.tpv12` as shipped, `par.dip = 60`, `insertFaultType = 1` (fault inserted, no wedge degeneration), term 8 s, cut to t <= 4.5 s for the table.
Column 8 of `faultst000dp000` = n-stress (+ = tension). Script: `curve.py`; outputs `curve_wedge.txt`, `curve_insertion.txt`. MEASURED (same script, single runs,
no independent oracle); not audited.

| mesh | dx (m) | n_max (MPa) | at t (s) | sr_max (m/s) | abs(tau) at peak (MPa) |
|---|---|---|---|---|---|
| wedge | 500 | +5.0989 | 2.580 | 8.922 | 0.200 |
| wedge | 300 | +7.7259 | 2.493 | 10.473 | 0.200 |
| wedge | 250 | +8.7191 | 2.460 | 10.987 | 0.200 |
| wedge | 200 | +9.8515 | 2.414 | 11.634 | 0.200 |
| insertion | 500 | +9.3827 | 2.879 | 8.662 | 0.200 |
| insertion | 250 | +13.7161 | 2.538 | 10.223 | 0.200 |

Wedge 500 m and 200 m reproduce the 2026-11-24 evidence (+5.0989, +9.8515) to all printed figures: same binaries, same inputs, deterministic CPU run; the 500 m file here
is a fresh run of this session.

## Findings
1. Wedge: n_max = 5.094 + 5.197 ln(500/dx) MPa, max residual 0.023 MPa over 4 points. Out-of-sample: fitted on 500 + 200 only, the log law predicts 7.748 (300) and
   8.694 (250) against measured 7.726 and 8.719. Power laws fit worse (p = 0.5: 0.140 MPa, p = 1: 0.276 MPa max residual). Peak slip rate also grows ~ log (8.92 -> 11.63 m/s).
   There is no converged value on this curve: the trace-node peak keeps growing at about 5.2 MPa per e-fold of refinement.
2. RULED OUT: the spike is a wedge-mesh wiring/degeneration artifact. The insertion-mesh control (no wedge cells, no C_degen) shows the same spike and it is LARGER
   (+9.38 at 500 m vs wedge +5.10; +13.72 at 250 m vs +8.72), also growing with refinement (+4.33 MPa per ln 2 on two points, vs +3.6 for the wedge).
   This closes the question the 2026-11-24 note left open ("cause NOT explained") only in the negative: the real fix (PR #281) is not what leaves it.
3. Shear is pinned at the 0.200 MPa cohesion at every peak, in every run: the node is fully slip-weakened when n peaks, so the peak is a normal-stress transient under
   zero frictional resistance, during the slip-rate rise.
4. n_max is affine in sr_max within the wedge series (slope ~1.7 MPa per m/s) but NOT across meshes (insertion 500 m: sr_max 8.66, n_max 9.38; wedge 500 m: 8.92, 5.10), so
   sr_max does not by itself determine n_max. "n/sr_max is not constant" stands; affine within one mesh is a fit, not a mechanism.
5. TPV36 (data from `wedge-trace-audit-dipsweep-2026-11-13/dip_sweep.csv`, not re-run): trace slip at 13 s is 0.10-0.14 m at 30/45/60 deg, against about 7-8 m at the TPV12 trace
   (README of that dir). The trace there barely slips, which is consistent with a slip-rate-driven transient being absent, but this is an observation, not a controlled test (TPV36 differs
   in friction, stress, nucleation, term).

## Verdict
Still unexplained physically. Measured: the trace-node n-stress peak of TPV12 at 60 deg diverges logarithmically with dx on the wedge mesh and is present, larger, on the plain insertion
mesh, so it is not a wedge-wiring defect. A log divergence is the signature of a field singularity at the free-surface / dipping-fault intersection (not shown here: that needs an
analytic or independent-code comparison, e.g. a fine-grid overlay). No gate depends on this node's peak; a gate that did could not converge.

## Run notes (corrections to the predecessor's handoff)
The wedge run that had finished when the handoff was written was dx 300, not dx 250 (dx 250 ran later: rc 0, 771 s). dx 400 is refused by the mesh builder (rc 33, `dx400_refused.meta`);
dx 300 is accepted (rc 0, 434 s). Raw trace files: `wedge_dx*.txt`, `ins_dx*.txt`.
