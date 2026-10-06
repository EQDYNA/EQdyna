# test.tpv34

SCEC TPV34: **Imperial Fault, Model 1** (https://strike.scec.org/cvws/tpv34docs.html,
`TPV34_Description_v10.pdf`, 2016-01-31). A vertical right-lateral planar
fault, 30 km along strike (x in [-15, 15] km) by 15 km deep, reaching the
free surface, in the **3D CVM-H velocity structure of the Imperial Valley**,
linear slip-weakening friction. The initial stresses scale with the local
shear modulus (`tau0 = 30 MPa * mu/mu0`, `sigma0 = 60 MPa * mu/mu0`,
`mu0 = 32.03812032 GPa`) and nucleation is an **added shear stress at t = 0**
(4.95 MPa * mu/mu0 inside r = 1.4 km of (x = 0, depth 7.5 km), cosine taper
to 2 km), so the fault fails spontaneously: `par.C_nuclea = 0`, and
`swtwNucleation` needs -- and has -- no TPV34 branch: the nucleation here
is a stress field, not a forced-rupture-time formula.

## The velocity structure: extracted from CVM-H, shipped as data

SCEC supplies no data files for TPV34; the spec (Part 2) gives the recipe
instead: convert each point to UTM zone 11 NAD27, query CVM-H with `vx_lite
-s -z dep`, take Vp/Vs/rho from columns 17-19, refuse non-positive values,
clamp to Vp 2984 / Vs 1400 m/s / rho 2220.34 kg/m3 below those velocities,
and use CVM-H throughout the whole model domain. `extract_cvmh_grid.py` is
that recipe, run ONCE, on the uniform grid of **element centres** of this
case's box at the gate spacing:

| file | content | used for |
|---|---|---|
| `extract_cvmh_grid.py` | the spec Part 2 recipe, committed | regenerating the grid (needs vx_lite + the CVM-H data) |
| `tpv34_cvmh_grid_500m.txt.gz` | 120 x 64 x 52 = 399,360 rows `x y z vp vs rho` on the element centres of the 60 x 32 x 26 km box at 500 m; 31,029 rows (7.8 %) at the spec's clamp floor; Vs 1400..4499 m/s; provenance in its `#` header | `par.mat` (the `n2mat == 6` 3D material grid) and, through `tpv34Tools.faultState`, every fault node's `tau0`/`sigma0` |

Provenance, stated rather than hidden:

* **Model:** CVM-H data files as distributed by UCVM (SCECcode/cvmh main,
  `model/config`), kept once read-only in `~/shared_dataset/scec_cvmh.15.1.1/`
  with per-file md5 in its `MANIFEST.json`. **The spec names CVM-H 15.1.0;
  UCVM ships 15.1.1.** The hypocentre point query (Vp 5676.89, Vs 3436.40,
  rho 2651.02 -> mu/mu0 = 0.977) reproduces the mu/mu0 implied by EQdyna's
  2016 SCEC submission at that point (`faultst000dp075`: 34.004 / 58.628 MPa
  -> 0.977), which was made with 15.1.0, so the two releases agree there to
  three digits. Nothing more is claimed about their difference elsewhere.
* **Code:** `vx_lite` built from SCECcode/cvmh main (reports Version 11.9.0).
* **Grid, not interpolation.** The solver's `n2mat == 6` branch
  (`meshgen.f90 setElementMaterial`, `meshgen.py build_elements`) gives each
  element the grid cell **nearest its centre**, clamped to the grid. Every
  uniform-belt element centre IS a grid point, so those elements carry their
  own CVM-H sample exactly; the geometrically stretched and PML elements
  outside the belt read the nearest sample. Nothing is interpolated.
  Consequently the case accepts **only `par.dx = 500 m`**, the spacing the
  grid was sampled at (`lib.requireFaultGeometryResolution`, `availableDx =
  [500]`): a finer or coarser mesh needs a fresh extraction with the script,
  not a resampling of this file.

## Fault stresses from the same samples

`tpv34Tools.faultState` takes mu at each fault node as the **mean of the
adjacent element centres on both sides of the fault** (the spec's own
suggested method; 8 cells for an interior node, 4 at the free surface),
then `tau0`, `sigma0`, the nucleation increment and the cohesion
`C0 = 0.000425 MPa/m * (2400 m - depth)` above 2400 m. At 500 m the
hypocentre node averages to mu/mu0 = 0.991 (tau0 34.65 MPa, sigma0
59.48 MPa) where the point value is 0.977: the medium varies over the
+-250 m the average spans. mu/mu0 ranges 0.136 (clamp-floor sediments at the
surface, sigma0 8.2 MPa) to 1.596 on this fault. 25 nodes start above
yield (tau0 > mu_s * sigma0, cohesion aside), the r <= 1.4 km disc.

## Solver generalisation this case introduced: `par.n2mat == 6`

`bMaterial.txt` rows `[x, y, z, vp, vs, rho]` (EQdyna frame, metres).
`readInputFiles.f90 buildMaterialGrid3D` / `checkInputConsistency.
build_material_grid3d` derive the grid from the rows themselves (origin =
min, spacing = smallest positive offset, count from the extent) and refuse,
with `ERR_CFG_MATERIAL_GRID_INVALID` (17), an incomplete block, a duplicate
or off-grid row, or a non-positive property. The existing `n2mat` 3/4/5
paths are untouched.

## Gate configuration (what the committed reference is)

`par.dx = 500` m, `par.term` = the sweep's gate term (5 s), `dt =
0.5*dx/vmaxPML`, 4 ranks `(nx, ny, nz) = (2, 1, 2)`, `friclaw = 1`, `tpv =
34`. Side/bottom fault borders carry `mu_s = 1000` (spec: slip is zero on
the borders), as tpv8/tpv29/tpv35 do. Stations: the spec's 35 on-fault
(x in {-12, -6, 0, 6, 12} km x 7 depths; the 2.4 km row is rounded to the
nearest fault plane, 2.5 km -> `dp025`, because `setOnFaultStation` needs an
exact node) and 56 off-fault (`body<Z>st<X>dp<D>`, spec coordinates).

The gate is a coarse REGRESSION check, not a spec-accuracy claim: the spec
asks for 25-50 m node spacing (recorded as EXCLUDED, no single number, in
`testsys/e2e/full_specs.py`) and 20 s; at 500 m the clamp-floor Vs of
1400 m/s is under 3 nodes per wavelength at 1 Hz and the cohesive zone is
under one element.

## Independent validation

`testsys/parity/evidence_tpv34_scec_comparison.py` compares a gate-config
run with EQdyna's 2016 SCEC submissions (`scec_archive/tpv34/`, 50 m, 20 s):
what is comparable at 500 m / 5 s (initial stresses at the 35 on-fault
stations, early rupture-front arrival), labelled as such; a report, not a
gate.
