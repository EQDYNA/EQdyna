# test.tpv31

SCEC TPV31: a single, vertical, planar, right-lateral strike-slip fault
(30 km along strike x 15 km deep) in a linear elastic half-space, with a
DISCONTINUOUS 1D velocity structure (3 jumps at 2400/5000/10000 m depth),
a (mu/mu0)-scaled depth-independent-shape initial-stress tensor, and a
static (time-independent) additive nucleation shear-stress patch -- no
forced-rupture-time formula. Spec: TPV31_32_Description_v03
(strike.scec.org/cvws); `scratch/specs/` holds the fetched PDF.

## Gap confirmation (row 151, 2026-10-06): NOT BLOCKED

Three questions, three existing mechanisms, zero src/ changes:

1. **Material.** TPV31's table is piecewise-LINEAR with 3 discontinuities
   (verbatim "2400-"/"2400+" duplicate rows, spec Part 2 p.4). EQdyna's
   `n2mat==4` (1D layered, piecewise-constant-per-element; already live via
   test.meng2023a) reproduces this exactly: this case precomputes one row
   per `par.dz` mesh layer via `np.interp` on the spec's verbatim table,
   evaluated at each layer's CENTER depth (never exactly at a jump depth).
2. **Initial stress.** The spec's own Part 4 "Method 2" (store only the
   on-fault traction, no off-fault stress tensor, no gravity body force
   needed for a purely elastic case) is EXACTLY the `on_fault_vars[...,7]`/
   `[...,8]` manual-fault-stress mechanism test.tpv8/test.tpv26/test.tpv29
   already use -- normal=-sigma33, strike shear=sigma13 (fault is the z=0
   plane, spec Part 2 p.6-7).
3. **Nucleation.** The spec's nucleation (Part 2 p.7) is an additive,
   TIME-INDEPENDENT shear-stress bump with a cosine taper -- not a
   forced-rupture-time formula. `faulting.f90:swtwNucleation` is a
   confirmed no-op for any unmatched `par.tpv` (read in full, 2026-10-06),
   so this case leaves `par.tpv=31` unmatched and `par.C_nuclea=0`, adding
   the bump directly into `on_fault_vars[...,8]` -- the same static-patch
   style as test.tpv8/test.tpv10.

## Geometry and hypocenter (spec Part 1, p.3)

Fault 30 km (strike) x 15 km (depth), vertical, reaching the surface,
planar -- no rough-fault geometry file needed (geometry decimation: N/A).
Hypocenter at the CENTER of the fault, 7.5 km deep:
`par.xsource,ysource,zsource = 0.0, 0.0, -7.5e3`.

## Material (spec Part 2, p.4)

| depth (m) | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---|---|
| 0 - 2400- | 4050 | 2250 | 2580 |
| 2400+ - 5000- | ramp to 5200/3050/2620 | | |
| 5000+ - 10000- | 5750 | 3450 | 2720 |
| 10000+ - 15000 | 6500 | 3800 | 3000 |

## Initial stress tensor and nucleation (spec Part 2, p.6-7)

sigma22=0; sigma11=sigma33=(-60.00 MPa)(mu/mu0); sigma13=(30.00 MPa)(mu/mu0);
sigma23=sigma12=0. mu0=32.03812032 GPa. Nucleation: tau_nuke(r)=(4.95 MPa)
(mu/mu0) for r<=1400 m, cosine-tapered to 0 by r=2000 m, mu evaluated at the
LOCAL node depth (matches the spec's own worked check: 34.95 MPa total shear
at the hypocenter, depth 7500 m).

## Friction (spec Part 2, p.8)

mu_s=0.580, mu_d=0.450, D0=0.18 m (shared with TPV32). Frictional cohesion
C0 = 0.000425 MPa/m * (2400 m - depth) for depth<=2400 m, else 0
(1.02 MPa at the surface, 0 by 2400 m).

## Resolution tier

Gated at `par.dx=500` m as a coarse regression gate. The spec's own 50 m
resolution request (Part 2 p.9, mandatory for TPV31) is recorded, not run,
in `testsys/e2e/full_specs.py['test.tpv31']`.

## Stations

The 30 on-fault / 18 off-fault stations in `st_coor_on_fault`/
`st_coor_off_fault` are the spec's own (Part 5/6, p.12-21), used verbatim.
