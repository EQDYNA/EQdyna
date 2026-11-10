# Row 149 force budget at the fault trace (TPV12, 60 deg dip, free surface)

## Verdict

PARTIAL. The elastic case (TPV12) and a vertical stand-in were measured at 200 and 150 m;
TPV13, the registered vertical control (TPV8) and the outcome-curve falsification spike were
NOT run (the work stopped at its tool budget). Read every mechanism statement below as
"consistent with the measurements", not as "falsified against the dp000-vs-dx curve".

What is measured: the budget closes to round-off, the instrumented build reproduces the
station value exactly, and the 200 -> 150 m rise of the TPV12 trace normal stress
(+1.69 MPa) is +1.33 MPa of element stress divergence and +0.36 MPa of hourglass
resistance. Damping, gravity and inertia are below 0.25 MPa and shrinking.

What is inferred: the stress part behaves like a bounded, non-smooth corner field
(element-stress exponents 0.34 to 0.74 in dx, where a smooth field gives 1 and a log
singularity gives 0); the hourglass part at the trace node is the only part that GROWS
under refinement. A point value at the trace therefore converges slowly at best, while
fixed-length averages (morning evidence, TPV13) converge.

## Runs

Base SHA 516da2c plus `diag_budget_instrumentation.patch`; Fortran, 4 ranks, par.term 4.5.

| case | dx (m) | wall (s) | dp000 n-stress near 4 s, this build (MPa) | uninstrumented, morning (MPa) | rows compared |
|---|---|---|---|---|---|
| TPV12 | 200 | 476 | -3.9865 | -3.9865 | 297 / 297 |
| TPV12 | 150 | 1123 | -2.2977 | -2.2977 | 396 / 396 |
| TPV12 with par.dip = 90 (STAND-IN, not a registered case) | 200 | 393 | 0.0000 | not run | 257 |
| TPV12 with par.dip = 90 (STAND-IN, not a registered case) | 150 | 928 | 0.0000 | not run | 343 |

The stand-in changes one parameter of TPV12 (dip 60 -> 90) and nothing else; its initial
stress and nucleation are TPV12's and it is not a SCEC benchmark.

## Closure (closure.csv)

- sum of element parts vs assembled nodal force: <= 1.4e-15 relative
- m dv/dt vs post-faulting nodal force: <= 3.1e-15 relative
- station metric rebuilt from the parts vs the code's value: <= 9.1e-13 MPa

Metric as coded (faulting.f90:93-95, 145-155): per direction e,
T_e = (mS mM ((vM - vS).e + [normal only] (uM - uS).e / dt) / dt + mS FM.e - mM FS.e) / ((mS + mM) arn) + T0 C_elastic.
Master is the y+ side (hanging wall when dipping). Parts are reported per side as
F_part.e / arn with the sign that makes them add to T_e.

## Diagnosis (TPV12, trace node x = 0, near 4 s, MPa; at4s.csv, elem_stress.csv)

| part, fault-normal | 200 m | 150 m |
|---|---|---|
| stress divergence, hanging wall | -6.65 | -5.48 |
| stress divergence, footwall | -8.62 | -7.14 |
| hourglass (displacement), hanging wall | +2.88 | +3.30 |
| hourglass (displacement), footwall | +4.41 | +4.72 |
| damping (stiffness + hourglass) | < 0.01 | < 0.01 |
| gravity | 0 | 0 |
| inertial residual | +-0.22 | +-0.12 |
| station value | -3.99 | -2.30 |

Second fault row: Tn -15.21 -> -13.01, hourglass -1.42 -> -1.28 (hanging wall) and
-2.30 -> -1.99 (footwall), i.e. opposite in sign to the trace node: a couple.

Top-row element stress, 200 -> 150 m: sig_nn -11.04 -> -9.34 (HW), -6.90 -> -5.58 (FW);
sig_zn -4.39 -> -3.86 (HW), -1.71 -> -1.55 (FW); sig_zz -3.57 -> -3.03 (FW).
sig_zn and sig_zz must vanish at a free surface; the one-point element carries them at
its centroid, and the trace-node stress part equals the top element's sig_nn -/+ sig_zn
(checked to three figures). Two-point exponents in dx: sig_zn 0.44 / 0.34, sig_zz 0.57,
sig_zd 0.59 / 0.65, sig_nn 0.58 / 0.74. Two grid sizes give one exponent each with no
error bar; treat them as indicative.

Corner analysis (corner_eigen.py, Poisson ratio 0.25, written before the 150 m data):
for a sliding interface meeting a free surface at 60 / 120 deg the smallest eigenvalue
is 1 with next roots near 1.38 and 1.61 (complex), so stress is bounded but its
correction terms go as r^0.38 to r^0.61; at 90 deg only integer roots exist. This matches
the measured exponents and predicts that the vertical geometry is smooth. It is a
static, linear analysis of a dynamic sliding state and is not a proof.

Vertical stand-in: the normal-traction change is zero by antisymmetry (max 4.1e-13 MPa
over the run), so it CANNOT test the n-stress outcome curve. It can only show the
per-side parts: stress divergence +-0.64 -> +-0.40, hourglass 0.043 -> 0.164, top-element
sig_zn 0.22 -> 0.12, sig_zz 0.64 -> 0.37 (exponents about 2), against 2.9 to 4.7 hourglass and
1.6 to 4.4 sig_zn in the dipping case. The closest REGISTERED control is TPV8 (vertical,
elastic, surface-breaking, same material); it was not run.

Figures: `fig_timeseries_normal.png`, `fig_timeseries_slipdir.png`, `fig_parts_vs_dx.png`,
`fig_surface_violation.png`. The TPV8 and TPV13 panels are empty because those runs do not exist.

## Not run

- TPV13 200 / 150 m instrumented (station targets +1.9555 / +5.4215); no plastic return factor was dumped
- TPV8 200 / 150 m (registered vertical control)
- 250 m, 100 m
- any spike that changes the solver and re-measures dp000 vs dx (so no rule-4a falsification of a cause)
- fixed-length means for TPV12 (the morning's are TPV13 only)
- any python-jax run

## Reproduce

    # in a fresh tree checked out at 516da2c:
    patch -p1 < diag_budget_instrumentation.patch     # or: python3 apply_budget_edits.py <tree>
    ./install-eqdyna.sh -m ubuntu
    # per case: create the case at dx, par.term = 4.5, run with 4 ranks into <runs_root>/<name>
    python3 budget_ana.py <runs_root> <morning_root> <out_dir>
    python3 budget_figs.py <out_dir> <fig_dir>

Run directory names expected by budget_ana.py are in its RUNS list.


## Side finding: why the TPV12 trace starts at exactly 0.0 and TPV13 does not

Traced. TPV12 is Method 1 (`C_elastic = 1`): the elements start stress-free (the dump shows every element stress
exactly 0 at step 1) and the initial traction is the per-node table value added in `src/fortran/faulting.f90:147`
(`+ nsdInitTractionVector(1)*C_elastic`). That value is set in `case_input/test.tpv12/user_defined_params.py:229`
(`sigEff = -SIGN_SHALLOW_PA_PER_M*downDipDistance`, stored at `:245`), which is exactly 0 at down-dip distance 0.
TPV13 is Method 2 (`case_input/test.tpv13/user_defined_params.py:118`, `C_elastic = 0`): the `:147` term is multiplied
by zero and the initial traction is whatever the element stresses produce. Those are set per element at the
element-centre depth by `src/fortran/meshgen.f90:201` (`setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)), elemCount)`,
TPV13 branch at `:1495`), so the trace node starts at the top-row centroid stress, which is proportional to dx:
-1.8410 / -0.7379 / -0.5536 MPa at 500 / 200 / 150 m is -3.68 kPa per metre of dx at all three. Not confirmed from a
TPV13 budget dump (TPV13 was not run here).

## Proposal (nothing here is implemented; no source, test or reference file is changed)

Ruled out by the owner and not proposed: a tension cap, adopting any treatment in this change, moving any reference.

**Option A. Report a fixed-length average beside the point value (post-processing only). RECOMMENDED.**
- Mechanism: the point value at the trace samples a corner field that varies like r^0.4 to r^0.6; the resultant over a
  fixed down-dip length does not. Measured this morning for TPV13 (250 / 200 / 150 m): top 500 m -15.38 / -15.07 /
  -14.87 MPa (3 %), top 1000 m -18.80 / -18.69 / -18.51; top 250 m -9.93 / -8.74 / -8.08 is NOT converged.
- Code: none in the solver. One analysis script that integrates the piecewise-linear nodal n-stress of the
  `fault.dyna.r.nc` profile (or of the dp000, dp005, ... stations) over the top L metres, used by the cross-code
  overlay; Fortran and python-jax outputs are read by the same script, so both backends are covered without edits.
- References moved: none. Cost: hours. It does not make dp000 itself converge and it is not the SCEC point station;
  the user documentation must say that dp000 n-stress on a dipping surface-breaking fault is not a converged quantity.
- Falsification test (written before any result): for L = 500 m the averaged n-stress at 4 s must agree within 5 %
  across 250 / 200 / 150 m for BOTH TPV12 and TPV13. TPV13 passes on this morning's numbers; TPV12 is not measured. A
  spread above 5 % in TPV12 kills the option as a convergence statement.

**Option B. Change the hourglass control in the elements that touch the trace (scratch spike first).**
- Mechanism: the hourglass part is the only part of the trace-node traction that does not shrink with dx (TPV12 +2.88
  to +3.30 and +4.41 to +4.72 MPa per side from 200 to 150 m; TPV13 trace-node total +6.66 / +8.68 / +10.06 at 250 /
  200 / 150 m this morning). It is the Kosloff-Frazier stiffness form acting on TOTAL displacement
  (`calcHourglassResist.f90`, `C_hg = 1` compiled in `globalvar.f90:105`), with opposite sign on node row 2.
- Code, Fortran: `calcHourglassResist.f90` (per-element switch to the existing viscous branch, or a scaled stiffness,
  for elements with a node on the trace), `globalvar.f90` (flag array), `meshgen.f90` (set the flag),
  `assembleGlobalMass.f90` (`calcSSPhi4Hrgls`). python-jax: `assembleGlobalKU.py` (the fused element loop holds the
  hourglass force), `assembleGlobalMass.py`, `meshgen.py`.
- References moved: every surface-breaking case has trace elements, so frt, nc and station references of all of them
  are at risk and must be shown inside their bounds or regenerated as their own reviewed change (rule 7).
- Cost: days, plus a stability check (the viscous form does not resist a static hourglass mode).
- Evidence against, stated plainly: this morning's plastic-hourglass spike moved the TPV13 level by about -5.5 MPa but
  left the 200 to 150 m rise unchanged (+3.21 against +3.47). In TPV12 the hourglass part carries only one fifth of
  the rise.
- Falsification test (written before any result): scratch spike with the viscous form (`C_hg = 2`, one line) on
  TPV12 at 200 and 150 m. The 200 to 150 m rise of dp000 at 4 s is +1.69 MPa now; the stress part alone accounts for
  +1.33. If the rise does not drop to +1.4 MPa or less, the hourglass stiffness is not what drives the outcome curve
  and the option is dead. If the run goes unstable the option is dead in that form.

**Option C. Integrate the trace elements fully (2x2x2, or selective-reduced) instead of one point plus hourglass.**
- Mechanism: with one point the element stress is constant, and the trace-node force is then (element normal stress)
  -/+ (the element's own free-surface traction sig_zn), measured to three digits below. Only a stress that varies
  inside the element can meet the free surface and the fault at a 60 degree corner at once.
- Code, Fortran: `calcElemKU.f90`, `assembleGlobalKU.f90`, `calcGlobalShapeFunc.f90`, `calcHourglassResist.f90` (skip
  those elements), `assembleGlobalMass.f90`, and the per-element stress and plastic state (`stressArr`,
  `stressCompIndexArr`, `pstrain`) become per-integration-point for the flagged elements. python-jax:
  `assembleGlobalKU.py`, `assembleGlobalMass.py`, `meshgen.py`, `driver.py`.
- References moved: all surface-breaking cases. Cost: weeks; this is a new method inside the solver, not a patch.
- Falsification test (written before any result): TPV12 dp000 at 4 s must differ by less than 0.5 MPa between 200 and
  150 m (now 1.69), and the dip = 90 control must stay at zero normal-stress change. Because the corner field itself
  is non-smooth, a full-integration element may still converge only like dx^0.4; if the 200 to 150 m difference stays
  above 1.0 MPa the option is dead.

**Option D. Change nothing; document.** Record in the user documentation that the trace-node n-stress of a dipping
surface-breaking fault is a mesh-dependent point value, and keep it out of every convergence claim. Cost: one doc
change. Falsification: none needed; it makes no convergence claim.

Recommendation: A (with D's documentation sentence), because it is the only option already supported by a
measurement against the outcome and it moves no reference. B's one-line spike is the cheapest next measurement and
ranks B against C; it was not run here.
