# test.tpv26

SCEC TPV26: a single, vertical, planar, right-lateral strike-slip fault
(40 km along strike x 20 km deep) in a linear elastic half-space, with a
depth-dependent (b11, b33, b13) initial-stress tensor, frictional cohesion
raised in the top 5 km to suppress a free-surface instability, and smoothed
forced-rupture nucleation. Spec: TPV26_27_Description_v13 (strike.scec.org/
cvws/tpv26_27docs.html, Jan 9 2014); `scratch/specs/` holds the fetched PDF.

## Gap confirmation (row 150, 2026-10-06)

A sibling mission (SCEC TPV12/TPV13) found those cases need SCEC "Method
2" nucleation (full volumetric stress tensor + explicit gravity body force
everywhere, nucleation by locally-reduced friction) which the code cannot
represent (the element-assembly routine hard-wires gravity off for an elastic
run; the stress-initialization routine cannot hold an asymmetric tensor).
Direct spec read confirms TPV26/27 are NOT blocked by that gap: nucleation is
the ORDINARY smoothed forced-rupture formula (Part 5 p.12, byte-identical in
symbols/constants to TPV22/23/29/30/36/37/201's), and the initial stress
tensor is the SCEC "Method 1" (stress-change) kind already implemented for
TPV29/30's own depth-taper machinery -- only the on-fault-resolved
shear/effective-normal stress is needed for this case's elastic solve (spec
Part 1 bullet 2), not an explicit off-fault stress tensor or gravity body
force.

## Geometry and hypocenter

Fault 40 km (strike) x 20 km (depth), vertical, reaching the surface.
Hypocenter 15 km from the left edge, 10 km deep -- the SAME hypocenter
EQdyna already uses for test.tpv29/test.tpv30 (`par.xsource,ysource,zsource =
-5.0e3, 0.0, -10.0e3`). Planar, so no rough-fault geometry file is needed
(unlike tpv29/30) -- the fault grid is the `par.fx`/`par.fz` linspace, same
style as test.tpv8.

## Initial stress tensor (spec Part 3, p.6-7)

Effective stress is parameterized by a lithostatic vertical component and a
depth-tapered horizontal deviatoric component rotated by an angle, resolved
onto the fault by this compset's own `on_fault_vars` loop (same code shape as
test.tpv30, different numeric coefficients):

| coefficient | value |
|---|---|
| b11 | 0.926793 |
| b33 | 1.073206 |
| b13 | -0.169029 |
| Omega taper | 1.0 (depth <= 15 km) -> 0.0 (depth >= 20 km), linear between |
| gravity g | 9.8 m/s^2 exactly |

## Friction (spec Part 3, p.8)

mu_s=0.18, mu_d=0.12, D0=0.30 m. Frictional cohesion C0 = 4.00 MPa at the
surface, tapering linearly to 0.40 MPa by 5 km depth, constant below --
`par.on_fault_vars[...,4]`.

## Resolution tier

Gated at `par.dx=500` m as a coarse regression gate. The spec's
own 100 m / 50 m resolution request (Part 3 p.9) is recorded, not run, in
`testsys/e2e/full_specs.py['test.tpv26']` -- no finer `on_fault_vars` grid is
shipped in this directory yet; this planar fault needs no geometry-file
download the way tpv29/30's rough surface does, only a finer loop at setup
time, so running the full tier is a scheduling decision, not blocked on data.

## Stations

The 12 on-fault / 6 off-fault stations and filenames in `st_coor_on_fault`/
`st_coor_off_fault` are the spec's own (Part 8/9, p.21-26), used verbatim.
