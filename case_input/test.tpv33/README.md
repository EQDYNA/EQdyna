# test.tpv33

SCEC TPV33: a single, vertical, planar, right-lateral strike-slip fault
(16 km along strike x 10 km deep) in a linear elastic half-space, centered
in a 1.6 km thick FAULT-PARALLEL LOW-VELOCITY ZONE (fault zone guided
waves), a constant-amplitude initial-stress tensor with a static additive
nucleation shear-stress patch -- no forced-rupture-time formula. Spec:
TPV33_Description_v04 (strike.scec.org/cvws); `scratch/specs/` holds the
fetched PDF.

## Gap confirmation (row 152, 2026-10-07): NOT BLOCKED

The question posed: does the EXISTING `n2mat==6` 3D material-grid mechanism
(landed for SCEC TPV34, later hardened by a grid-covers-mesh assertion)
already cover a thin FAULT-PARALLEL low-velocity zone, or does TPV33 need a
genuinely new material-grid mode?

**Confirmed: reuse, zero src/ changes.** Evidence, by direct code read:

- `src/fortran/meshgen.f90`, subroutine `setElementMaterial`: the `n2mat==6`
  branch does a per-axis NEAREST-CELL lookup against a self-describing 3D
  grid (`matGridOrigin`/`matGridSpacing`/`matGridCount`, one index per axis,
  independently). Nothing in this lookup requires variation along all three
  axes -- a grid that is uniform-spaced along one axis and carries only two
  (identical-valued) edge samples along the other two reproduces a field
  that varies along ONE axis only, exactly.
- `src/fortran/readInputFiles.f90`, subroutine `buildMaterialGrid3D`: derives
  `matGridOrigin`/`matGridSpacing`/`matGridCount` per axis purely from the
  rows present (min value, minimal positive offset, count), and requires
  `nmat == nx*ny*nz` (complete grid, no missing/duplicate cells) plus the
  grid-covers-mesh check that each axis's nearest-neighbour reach contains
  the mesh box from `bModelGeometry.txt`. None of this assumes 3-axis
  variation; a 2x2xN grid (two x edges, two z edges, N uniform y samples)
  satisfies every one of these checks the same way TPV34's dense grid does.
- `src/python/eqdyna/checkInputConsistency.py`, function
  `build_material_grid3d`, mirrors the Fortran logic and error codes exactly
  (same invalid-grid error raised on both sides), so the same grid design is
  valid on both gated backends.

TPV33's velocity structure (spec Part 2, p.4) depends ONLY on the spec's
fault-normal coordinate `z`, which maps directly (no sign flip) onto
EQdyna's `y` axis -- the same axis mapping TPV34's own
`extract_cvmh_grid.py` documents (`x_eq=x, y_eq=z, z_eq=-y`) -- and is
CONSTANT along the spec's `x` (EQdyna `x`) and depth `y` (EQdyna `z`). This
case's material grid (`user_defined_params.py`) is therefore a `2 x N x 2`
grid: two edge samples along EQdyna x and z (so every element's x/z lookup
always resolves to one of two identical-valued edge rows, giving true
x/z-independence), and `N = (ymax-ymin)/dx` uniform element-centre samples
along EQdyna y reproducing the spec's three velocity bands at this gate's
`par.dx`. This is the same nearest-cell discretization principle as TPV34's
dense 3-axis CVM-H grid, just exploiting that two of TPV33's three axes
carry no material variation at all -- not a new architecture.

## Geometry and hypocenter (spec Part 1, p.3)

Fault 16 km (strike, -12000 m to 4000 m) x 10 km (depth), vertical, reaching
the surface, planar -- no rough-fault geometry file needed (geometry
decimation: N/A). Hypocenter 6 km from the left edge, 6 km deep:
`par.xsource,ysource,zsource = -6.0e3, 0.0, -6.0e3`.

## Material (spec Part 2, p.4)

| fault-normal distance | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---|---|
| \|z\| > 800 m (host rock) | 6000 / 5626 | 3464 / 3248 | 2670 |
| \|z\| <= 800 m (low-velocity zone, 1.6 km wide) | 3750 | 2165 | 2670 |

(two host-rock values: 6000/3464 m/s on the `z>800` side, 5626/3248 m/s on
the `z<-800` side -- both reproduced verbatim via the `_vp_of_y`/`_vs_of_y`
band functions.)

## Initial stress tensor and nucleation (spec Part 2, p.5)

`sigma0 = 60 MPa` constant (NOT mu/mu0-scaled, unlike TPV31/34/35 -- TPV33
has no material-dependent stress scaling). `tau0 = 30 MPa*(1-Rtau) +
tau_nuke`, with `Rx`/`Ry`/`Rtau` the spec's own taper (breakpoints -9800/
1100 m along strike, 2300/8000 m depth) and `tau_nuke(r) = 3.150 MPa` for
`r<=550` m, cosine-tapered to 0 by `r=800` m. No gravity (spec Part 1, p.1:
"There is no gravity in the model."). Worked check (spec p.6): total
initial shear at the hypocenter = 33.15 MPa, yield = 33.00 MPa.

## Friction (spec Part 3, p.8)

mu_s=0.550, mu_d=0.450, D0=0.18 m. No frictional cohesion anywhere (C0=0,
unlike TPV31/34).

## Resolution tier

Gated at `par.dx=400` m as a coarse regression gate. The spec's own
12.5-25 m fault-plane resolution request (mandatory, Part 2 p.7 -- "TPV33
requires twice the resolution of any earlier benchmark"), 50 m through the
low-velocity zone, 100 m outside it, is recorded, not run, in
`testsys/e2e/full_specs.py['test.tpv33']`.

## Stations

28 on-fault stations (Part 4, p.9-10, 7 strikes x 4 depths, used verbatim).
Off-fault: a VERIFIED SUBSET of the spec's 85 stations (Part 5, p.15-19) --
the four earth-surface dense transects (36 stations) plus the 12 "other
stations surrounding the fault" (also surface), 48 of 85 total, used
verbatim. The spec's depth-6-km transects (offsets including an asymmetric
+-0.1 km pair at the fault trace) are NOT transcribed -- a documented gap,
left open the same way test.tpv26/27/31/32 leave independent cross-code
validation open; none of the off-fault station picks selected for this
case's gate need them.

## Validation

No independent `evidence_tpv33_*.py` script is shipped, the same documented
choice already made for test.tpv26/test.tpv27/test.tpv31/test.tpv32 (see
`testsys/parity/README.md`). This case's comparison is the cross-backend
fortran/python-jax agreement gated by this case's own frt and station
comparison bounds.
