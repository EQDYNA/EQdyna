# test.tpv32

SCEC TPV32: identical to test.tpv31 (same fault geometry, stress tensor,
nucleation, friction, stations -- see that case's README for the shared
derivation and the "not blocked" gap confirmation) except for a
CONTINUOUS, piecewise-linear 1D velocity structure (spec Part 2, p.5; no
discontinuities). Spec: TPV31_32_Description_v03 (strike.scec.org/cvws);
`scratch/specs/` holds the fetched PDF.

## Material (spec Part 2, p.5)

| depth (m) | Vp (m/s) | Vs (m/s) | rho (kg/m^3) |
|---|---|---|---|
| 0 | 2200 | 1050 | 2200 |
| 500 | 3000 | 1400 | 2450 |
| 1000 | 3600 | 1950 | 2550 |
| 1600 | 4400 | 2500 | 2600 |
| 2400 | 4800 | 2800 | 2600 |
| 3600 | 5250 | 3100 | 2620 |
| 5000 | 5500 | 3250 | 2650 |
| 9000 | 5750 | 3450 | 2720 |
| 11000 | 6100 | 3600 | 2750 |
| 15000 | 6300 | 3700 | 2900 |

Reused via the same `n2mat==4` mechanism as test.tpv31: one row per
`par.dz` mesh layer, `np.interp`-evaluated at the layer's center depth on
the verbatim table above (no duplicate/jump rows this time -- truly
continuous).

## Resolution tier

Gated at `par.dx=500` m as a coarse regression gate. The spec allows TPV32
node spacing in [25, 50] m (Part 2 p.9) -- the spec's own preliminary tests
needed 25 m near the low-velocity surface layer for acceptable results, so
the range is wider than TPV31's mandatory 50 m. This case records 50 m (the
coarser end of that range) as its full-resolution tier in
`testsys/e2e/full_specs.py['test.tpv32']` -- a scheduling choice, not a
claim that 50 m is the spec's preferred value.
