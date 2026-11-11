# Row 149 force budget at the fault trace (TPV12 / TPV13, 60 deg dip, free surface)
 
## Verdict
 
COMPLETE as a diagnostic (no longer partial); nothing under `src/`, `testsys/` or any reference is changed.
 
- The trace normal stress of the dipping cases is the small difference of two large parts on BOTH sides of the fault:
  element stress divergence (compressive, shrinking with dx) and the stiffness hourglass resistance (tensile, the only
  part that grows with dx). Damping, gravity and inertia are below 0.25 MPa. Neither side alone carries the drift.
- Outcome test (rule 4a): replacing the stiffness hourglass by the viscous form (one line, scratch build) cuts the
  TPV12 200 -> 150 m rise of dp000 at 4 s from +1.69 to +0.65 MPa and the run stays stable. The hourglass stiffness is
  therefore a real driver of the outcome curve; the budget share (1/5 hourglass, 4/5 stress) is NOT a causal share.
- A fixed-length mean converges where the point value does not: top 500 m mean at 4 s agrees to 1.8 % (TPV12, 200 /
  150 m) and 3.4 % (TPV13, 250 / 200 / 150 m).
- The vertical cases do not drift because the fault-normal parts are antisymmetric across a vertical fault and cancel
  identically in the traction (TPV8: every normal part is zero; dip = 90 stand-in: per-side parts 5 to 20 times smaller
  and of opposite sign). A vertical case therefore cannot test the n-stress curve at all; it only shows the per-side parts.
 
## Runs
 
Base SHA 516da2c plus `diag_budget_instrumentation.patch`; Fortran, 4 ranks, par.term 4.5, one run per cell
(run-to-run spread not measured). dp000 n-stress near 4 s in MPa; morning = the uninstrumented-physics build of the
row149-trace-traction evidence.
 
| case | dx (m) | wall (s) | dp000, this build | dp000, morning | station rows (this / morning) |
|---|---|---|---|---|---|
| TPV12 | 200 | 476 | -3.9865 | -3.9865 | 297 / 297 |
| TPV12 | 150 | 1123 | -2.2977 | -2.2977 | 396 / 396 |
| TPV13 | 200 | 485 | +1.9555 | +1.9555 | 297 / 297 |
| TPV13 | 150 | 1145 | +5.4215 | +5.4215 | 396 / 396 |
| TPV8 (registered vertical strike-slip) | 200 | 164 | 0.0000 | not run | 257 |
| TPV8 | 150 | 464 | 0.0000 | not run | 343 |
| TPV12 with par.dip = 90 (STAND-IN, not a registered case) | 200 | 393 | 0.0000 | not run | 257 |
| TPV12 with par.dip = 90 (STAND-IN) | 150 | 928 | 0.0000 | not run | 343 |
| TPV12, viscous hourglass spike (`spike_hg2_viscous.patch`, clean build) | 200 | 383 | -7.9028 | n/a | 297 |
| TPV12, viscous hourglass spike | 150 | 901 | -7.2573 | n/a | 396 |
 
The spike run logs print `C_hg 2` (4 ranks each), the baseline `C_hg 1`.
 
## Closure (closure.csv, all eight budget runs)
 
- sum of element parts vs assembled nodal force: <= 1.4e-15 relative
- m dv/dt vs post-faulting nodal force: <= 5.9e-15 relative
- station metric rebuilt from the parts vs the value the code uses: <= 9.1e-13 MPa
 
Metric as coded (faulting.f90:93-95, 145-155): per direction e,
T_e = (mS mM ((vM - vS).e + [normal only] (uM - uS).e / dt) / dt + mS FM.e - mM FS.e) / ((mS + mM) arn) + T0 C_elastic.
Master is the y+ side (hanging wall when dipping). Parts are per side, F_part.e / arn, signed so they add to T_e.
 
## Diagnosis: fault-normal parts at the trace node, near 4 s, MPa (at4s.csv)
 
| part | TPV12 200 m | TPV12 150 m | TPV13 200 m | TPV13 150 m |
|---|---|---|---|---|
| stress divergence, hanging wall | -6.651 | -5.479 | -5.837 | -4.226 |
| stress divergence, footwall | -8.619 | -7.141 | -7.628 | -5.050 |
| hourglass (displacement), hanging wall | +2.878 | +3.300 | +8.719 | +10.320 |
| hourglass (displacement), footwall | +4.407 | +4.717 | +8.644 | +9.788 |
| damping (mass + stiffness + hourglass) | < 0.01 | < 0.01 | < 0.01 | < 0.01 |
| gravity, hanging wall / footwall | 0 / 0 | 0 / 0 | -0.721 / +0.723 | -0.541 / +0.542 |
| inertial residual | +-0.221 | +-0.122 | +-0.213 | +-0.137 |
| fault traction (station value) | -3.987 | -2.298 | +1.955 | +5.422 |
 
Rise 200 -> 150 m: TPV12 +1.69 = +1.33 stress + 0.37 hourglass; TPV13 +3.47 = +2.09 stress + 1.37 hourglass.
By side: the stress part rises more on the footwall (+1.48 / +2.58) than on the hanging wall (+1.17 / +1.61); the
hourglass part rises more on the hanging wall (+0.42 / +1.60) than on the footwall (+0.31 / +1.14). The gravity parts
of TPV13 are equal and opposite and cancel in the traction.
 
Second node row (its depth moves with dx, so rows are not at the same place): hourglass is of opposite sign to the
trace node, TPV12 -1.415 -> -1.284 (HW), -2.295 -> -1.986 (FW); TPV13 -1.242 -> -1.182 (HW), -4.236 -> -3.366 (FW).
Total n-traction of row 2: TPV12 -15.21 -> -13.01, TPV13 -15.56 -> -12.33.
 
Top-row element stress (one point, elem_stress.csv), 200 -> 150 m: TPV12 sig_nn -11.04 -> -9.34 (HW), -6.90 -> -5.58
(FW); sig_zn -4.39 -> -3.86, -1.71 -> -1.55. TPV13 sig_nn -11.74 -> -8.95 (HW), -6.14 -> -4.37 (FW); sig_zn -5.90 ->
-4.72, -1.51 -> -0.68. The trace-node stress part equals sig_nn - sig_zn (HW) and sig_nn + sig_zn (FW) of the top
element to three figures in both cases: the element carries a free-surface traction it cannot release.
 
Drucker-Prager return factor (TPV13, dumped as rjust): 1.00000 in the top two element rows on both sides at 4 s at
both grids, i.e. no trace element is yielding at 4 s; accumulated plastic strain there is 0.6e-4 to 3.7e-4 from
earlier yielding. Plasticity doubles to triples the hourglass part (8.7 to 10.3 against 2.9 to 4.7 MPa).
 
Vertical controls. TPV8: every fault-normal part is 0.000 on both sides and the normal-traction change is zero to
1e-27 MPa; its dp000 normal stress can test nothing about row 149. What it can show is the slip-direction budget at
a vertical trace: stress -0.093 -> -0.114, hourglass -0.041 -> -0.038 MPa per side, two orders below the dipping
cases. Stand-in (TPV12 with dip = 90, dip-slip): per-side stress +-0.638 -> +-0.396, hourglass -+0.043 -> -+0.164,
equal and opposite on the two sides, sum zero (max 4e-13 MPa). The dipping geometry breaks that mirror symmetry and
adds a slip-induced normal-stress change that the vertical geometry does not have.
 
All exponents and rises quoted here come from two grid sizes (200 and 150 m); they carry no error bar.
 
## Mechanism claims against the dp000-vs-dx outcome curve
 
| claim | status | evidence |
|---|---|---|
| plastic return misses the hourglass stiffness | FALSIFIED (morning) | rise +3.47 -> +3.21 MPa in TPV13 |
| the stiffness hourglass drives the rise | SURVIVES, partly | TPV12 rise +1.69 -> +0.65 MPa with the viscous form; +0.65 remains |
| the budget share is the causal share (hourglass only 1/5) | FALSIFIED | removing the hourglass stiffness removed 62 % of the rise, not 21 % |
| non-smooth corner field of the 60 deg wedge | UNTESTED against the curve | consistent with element-stress exponents 0.34 to 0.74; may be the remaining +0.65 |
| one-point element cannot meet free surface and fault together | UNTESTED against the curve | identity sig_nn -/+ sig_zn holds to three figures |
 
Corner eigen-analysis (`corner_eigen.py`, Poisson ratio 0.25), now with the singular-value scan run: for a sliding
interface at 60 deg with friction 0 and 0.1 there is no root with Re < 1 (scan of Re 0.05 to 0.98, Im 0 to 1.5 on a
0.005 grid, smallest normalized singular value 3.9e-4, no grid minimum below 2e-3), so no power singularity. The scan
confirms the root near 1.605 (stress correction like r^0.61). It does NOT confirm the secant root 1.376+0.290i quoted
earlier for friction 0; treat the r^0.38 figure as unverified. Static, linear analysis of a dynamic sliding state.
 
## Side finding (TPV13 initial trace traction): CONFIRMED from the dump
 
TPV12 (`C_elastic = 1`) starts with stress-free elements and a per-node initial traction that is exactly 0 at the
trace (`faulting.f90:147`, case parameters). TPV13 (`C_elastic = 0`) carries its initial stress in the elements, set
at the element-centre depth (`meshgen.f90:201`). The dump at the first step shows top-row sig_nn = -0.739 (200 m) and
-0.554 MPa (150 m) on both sides, and three times that in the second row: -3.69 kPa per metre of dx, the same as the
station start values -0.7379 / -0.5536. The initial sig_zn of those elements (-0.721 / -0.541) is balanced by the
gravity part at the trace node (-0.721 / -0.541).
 
## Two measured tests
 
Fixed-length mean normal traction at 4 s (means.csv; piecewise linear through the nodal values at x = 0, initial
traction included), MPa:
 
| case | length | 250 m grid | 200 m grid | 150 m grid | spread |
|---|---|---|---|---|---|
| TPV12 | top 250 m | not run | -10.833 | -10.389 | 4.2 % |
| TPV12 | top 500 m | not run | -14.888 | -14.620 | 1.8 % |
| TPV12 | top 1000 m | not run | -19.429 | -19.321 | 0.6 % |
| TPV13 | top 250 m | -9.931 | -8.737 | -8.075 | 21 % |
| TPV13 | top 500 m | -15.376 | -15.071 | -14.869 | 3.4 % |
| TPV13 | top 1000 m | -18.801 | -18.690 | -18.514 | 1.5 % |
| TPV12 viscous spike | top 500 m | not run | -15.459 | -15.572 | 0.7 % |
 
Option A criterion (written first): top 500 m within 5 % across grids for both cases. Result: PASS, 1.8 % and 3.4 %.
Limits: TPV12 has two grids only; the TPV12 500 m grid gives -15.681 (7 % from the 150 m grid) but has a single cell
inside the averaging length. The mean itself moves 0.6 to 0.95 MPa (4 to 6 %) between the two hourglass forms.
 
Option B criterion (written first): TPV12 rise of dp000 at 4 s, 200 -> 150 m, must drop from +1.69 to +1.4 MPa or
less and the run must stay stable. Result: PASS, -7.9028 -> -7.2573, rise +0.65 MPa, both runs completed (rc 0).
Limits: the spike changes the hourglass form in EVERY element, not only at the trace; it moves the level by -3.9 and
-5.0 MPa, away from the FaultMod 100 m value (-3.005); TPV13, rupture times and any later time were not examined.
 
## Proposal (nothing here is implemented)
 
Not proposed, by owner decision: a tension cap, adopting a treatment in this change, moving a reference.
 
**A. Report a fixed-length mean beside the point value. RECOMMENDED.**
- Mechanism: the resultant over a fixed down-dip length is insensitive to how the top element splits its force
  between stress and hourglass; the point value is not.
- Code: no solver change. One post-processing script reading `fault.dyna.r.nc` (or the dp stations) from either
  backend; Fortran and python-jax need no edit. User documentation states that dp000 n-stress on a dipping
  surface-breaking fault is a mesh-dependent point value.
- References moved: none. Cost: hours. Test and result: above, PASS.
 
**B. Replace the stiffness hourglass in the elements that touch the trace.**
- Mechanism: the displacement-proportional hourglass force at the trace node is tensile and grows under refinement.
- Code, Fortran: `calcHourglassResist.f90` (per-element form), `globalvar.f90` (flag array), `meshgen.f90` (set the
  flag), `assembleGlobalMass.f90`. python-jax: `assembleGlobalKU.py` (fused loop), `assembleGlobalMass.py`, `meshgen.py`.
- References moved: every surface-breaking case is at risk (rule 7). Cost: days plus a stability study, because the
  viscous form does not resist a static hourglass mode.
- Test and result: global one-line spike PASSES its kill criterion (+1.69 -> +0.65). Next falsification, not run: the
  same change restricted to trace elements must keep the rise at or below +1.0 MPa in TPV12 AND cut the TPV13 rise
  (+3.47) by half, with rupture times inside the case bound; otherwise dead.
 
**C. Integrate the trace elements fully (or selective-reduced) instead of one point plus hourglass.**
- Mechanism: only a stress that varies inside the element can meet the free surface and the fault at 60 deg at once.
- Code, Fortran: `calcElemKU.f90`, `assembleGlobalKU.f90`, `calcGlobalShapeFunc.f90`, `calcHourglassResist.f90`,
  `assembleGlobalMass.f90`, per-integration-point stress and plastic state. python-jax: `assembleGlobalKU.py`,
  `assembleGlobalMass.py`, `meshgen.py`, `driver.py`.
- References moved: all surface-breaking cases. Cost: weeks; a new method inside the solver. Not measured.
- Falsification: TPV12 rise below 0.5 MPa with the dip = 90 control unchanged; above 1.0 MPa the option is dead.
 
**D. Change nothing; document.** One doc sentence; makes no convergence claim.
 
Recommendation: A, with the documentation sentence of D. It is the only option that passed its test on both cases
and moves no reference. B is now ranked above C: its cheap spike survived, so the trace-only spike is the next
measurement worth making before any solver change is considered.
 
## Figures
 
`fig_spike_and_means.png` (point value, spike, fixed-length means; panels 2 and 3 share one scale),
`fig_parts_vs_dx.png`, `fig_surface_violation.png`, `fig_timeseries_normal.png`, `fig_timeseries_slipdir.png`.
Known limits of the figures: time-series rows share a scale within a case, not across cases; the TPV8 normal panel
is zero by symmetry (axis offset 1e-12 MPa); the TPV8 panel of `fig_surface_violation.png` is flat because the
plotted components vanish for vertical strike-slip; legends overlap some curves in the first panels.
 
## Not run
 
- TPV12 at 250 m and 100 m; TPV13 with the viscous hourglass; a trace-elements-only hourglass spike
- force parts of the spike runs (the dumps exist; only dp000 and the means were reduced)
- rupture-time (frt) comparison of the spike against the baseline; any time other than 4 s
- any python-jax run; any full-integration element
 
## Reproduce
 
    # fresh tree at 516da2c
    patch -p1 < diag_budget_instrumentation.patch      # or: python3 apply_budget_edits.py <tree>
    ./install-eqdyna.sh -m ubuntu                      # spike: also patch -p1 < spike_hg2_viscous.patch, clean build
    # per case: create at dx, par.term = 4.5, 4 ranks, into <runs_root>/<name> (names in the RUNS lists)
    python3 budget_ana.py <runs_root> <morning_root> <out_dir>
    python3 budget_extra.py <runs_root> <out_dir> <fig_dir> <morning_root>
    python3 budget_figs.py <out_dir> <fig_dir>
    python3 corner_eigen.py 0.25
