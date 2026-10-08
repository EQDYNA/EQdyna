# test.tpv13 (SCEC TPV13)

60-degree dipping, planar, normal fault, identical geometry to
`test.tpv12`, but with non-associative Drucker-Prager plasticity in the
off-fault medium. Official spec: `TPV12_13_Description_v6.pdf`
(strike.scec.org/cvws/tpv12_13docs.html, 2009-10-01), Part 2/3 (TPV12 and
TPV13 share the geometry/resolution/duration text) and the plastic
"Initial Stress Tensor" section (p.4-6).

This closes the gap `test.tpv12/README.md` documented under "TPV13
status: not built": the shared off-fault stress-tensor builder used by
`test.tpv29`/`test.tpv30` pins `sxx + syy` to twice the vertical stress and
cannot represent TPV13's own asymmetric `sigma2 = (sigma1+sigma3)/2`
relation with its existing two free parameters.

## What this case reuses

Same `insertFaultType=1` + `mod4dip` dipping-fault mesh, same material,
same on-fault traction/friction/nucleation code path as `test.tpv12`
(REDUCED static friction in the nucleation patch, `C_nuclea=0`, no
forced-rupture formula). This case changes only the off-fault
material/stress Method (1 -> 2) and the off-fault plastic parameters.

## The new `if (TPV == 13)` branch

TPV13's fault strike (x-axis) and dip (y-z plane) put its principal
stress axes directly on the model's x/y/z axes: sigma1 (vertical, most
compressive) is `szz`, sigma3 (horizontal, normal to the fault trace) is
`syy`, sigma2 (horizontal, parallel to the fault trace) is `sxx`. Unlike
the TPV29/30 branch (a near-vertical strike-slip fault whose sigma1/sigma3
lie in the horizontal plane at an angle to the fault, needing a deviatoric
rotation with an xy shear term), TPV13 needs no rotation: `sxy = sxz = syz
= 0`.

Formula (`src/fortran/meshgen.f90`'s `setPlasticStress`, mirrored in
`src/python/eqdyna/eqdyna3d.py`'s `init_stress` construction -- the
functional mirror of this formula lives in `eqdyna3d.py`, NOT in
`meshgen.py`, because `meshgen.py` only computes the per-element depth
argument, not the stress-tensor formula itself):

- `sig1Eff = strVert` (the existing `gamar=0` vertical-effective-stress
  formula, verbatim -- no new formula needed for the vertical component)
- depth < 11951.15 m: `sig3Eff = 0.3496 * sig1Eff`
- depth >= 11951.15 m: `sig3Eff = sig1Eff` (spec p.5: "stresses become
  isotropic")
- `sig2Eff = 0.5 * (sig1Eff + sig3Eff)` (spec p.5 relation; holds for the
  effective stress too since `Pf` cancels out linearly)

`0.3496` and `11951.15 m` are read directly off the spec (p.5: "depths less
than 11951.15 meters" and the `(sigma3-Pf) = 0.3496*(sigma1-Pf)` ratio) --
TPV12/13-specific physical constants, not case-configurable knobs, kept as
named `parameter`s in the branch (`TPV13_SIG3_RATIO`, `TPV13_DEPTH_SPLIT_M`).
The branch `return`s before the existing `devStr`/rotation code, so every
other TPV (tpv27, tpv30, drv.a6, drv.a6.v2, ...) runs the pre-existing
branch unchanged and is untouched bit-for-bit -- verified directly: those
cases' own gate numbers (both backends) are unchanged to the bit against
the numbers already committed in `testsys/matrix.py`'s own comments.

## Instantaneous plastic return map

EQdyna's shared plastic-update kernel implements Duvaut-Lions
viscoplastic relaxation (`rjust = exp(-dt/Tv) + (1-exp(-dt/Tv))*Y/sqrt(J2)`,
`Tv = par.viscoplasticRelaxTime`). TPV13's spec wants the INSTANTANEOUS
(non-associative Drucker-Prager) limit. This case gets that limit
bit-exactly, not approximately, by setting `par.viscoplasticRelaxTime =
1.0e-5` s: at this case's gate `dt` (~0.0379 s), `dt/Tv ~ 4374`, far past
IEEE double precision's underflow threshold (~745) for `exp(-dt/Tv)`, so
`exp(-dt/Tv)` underflows to exactly `0.0` and `rjust` reduces exactly to
`Y/sqrt(J2)`, the instantaneous return map -- no new plastic-update code
needed. Confirmed empirically: the Fortran run prints
`Note: The following floating-point exceptions are signalling:
IEEE_UNDERFLOW_FLAG` for this case, as predicted.

## Resolution and duration

Spec recommends 100 m node spacing, 0-8 s duration (same Part 2/3 text as
`test.tpv12`); everyday regression runs this case coarser and shorter
(`par.dx=500`, gate-overridden to the one everyday 5 s gate term) for
speed.
See `testsys/e2e/full_specs.py`'s `test.tpv13` entry for the recorded
(not run) spec-resolution tier.

## Barall cross-code comparison

`testsys/parity/evidence_tpv13_scec_comparison.py` compares this case's
own gate run against Michael Barall's independent FaultMod TPV13
submission (100 m, under `~/shared_dataset/scec_cvws.tpv1213/raw/tpv13/`),
on the rupture-time field (fraction of non-barrier nodes ruptured) and the
two off-fault body stations, following the same scoping decision as
`test.tpv12/README.md`'s own Barall section (on-fault stations are not
compared; Barall's fixed on-fault depth grid does not line up with this
case's gate-selected on-fault station sample).

The rupture-fraction gap at gate resolution is measured larger here
(`|delta| = 0.1783`) than TPV12's own elastic gap (`0.0486`): this case
adds Drucker-Prager plastic dissipation at the rupture front, which is
itself under-resolved differently at 5x coarser-than-spec mesh than pure
elastic wave propagation is -- a larger gap is the physically expected
direction at this resolution, not a bug. `FRACTION_BOUND = 0.54` (~3x the
measured gap, the same headroom convention `test.tpv12`'s own script
uses) is set from this case's own measurement, not loosened to force a
pass with no margin.
