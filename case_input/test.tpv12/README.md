# test.tpv12 (SCEC TPV12)

60-degree dipping, planar, normal fault in a uniform linear-elastic
half-space. Official spec: `TPV12_13_Description_v6.pdf`
(strike.scec.org/cvws/tpv12_13docs.html, 2009-10-01).

This is the elastic half of the TPV12/13 pair (no gravity, per-node
fault tractions set directly). TPV13, the Drucker-Prager-plastic sibling
(explicit gravity, a full off-fault stress tensor), is not built in this
repository yet -- see "TPV13 status" below.

## What this case reuses

Geometry, the dipping-fault mesh construction, and the nucleation
mechanism are the same, already-validated machinery `test.tpv10` uses
(same dip, same fault footprint, same material: rho 2700, Vs 3300,
Vp 5716). This case changes only its own input parameters and spec-specific
initial-stress/friction/station values.

## Initial stress and friction

Depth-dependent normal/shear traction on the fault, resolved by down-dip
distance (spec Part 2, p.5-6):

- down-dip distance < 13800 m: `(sigma_n - Pf) = 7390.01 Pa/m * ddist`,
  `tau = 0.549847 * (sigma_n - Pf)`
- down-dip distance >= 13800 m: `(sigma_n - Pf) = 14427.98 Pa/m * ddist`,
  `tau = 0` (stresses become isotropic)

Friction: `mu_s = 0.70` outside the nucleation patch, `mu_s = 0.54` inside
it (a 3 km x 3 km square centered 12 km down-dip, along-strike centered);
`mu_d = 0.10`, `d0 = 0.50 m`, frictional cohesion `c0 = 0.2 MPa` everywhere.
Nucleation is the spec's own mechanism: a lower static friction coefficient
in the patch, applied to the same stress field used everywhere else --
no artificial forced-rupture formula is used.

## Resolution and duration

Spec recommends 100 m node spacing; everyday regression runs this case
coarser and shorter for speed (see `user_defined_params.py`'s own
`par.dx`/`par.term`). Published SCEC run duration: 0-8 s after nucleation.

## TPV13 status: not built

TPV13 adds non-associative Drucker-Prager plasticity off the fault and
needs an explicit-gravity, full-stress-tensor initial condition the same
way `test.tpv30` builds its own off-fault stress field.

That shared off-fault stress-tensor builder always produces a vertical
principal stress symmetric between the two horizontal components (its
`sxx + syy` is pinned to twice the vertical stress, for any choice of its
own two free parameters). TPV12/13's own spec instead asks, in its
shallower stress regime, for one horizontal component to equal the AVERAGE
of the vertical and the other horizontal component -- an asymmetric
relation the existing builder cannot represent with its current two free
parameters. Representing it correctly would need a change to that shared
builder, which this case's elastic half does not touch and does not need.
That is why TPV13 is not built here.
