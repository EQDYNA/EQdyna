# test.tpv35

SCEC TPV35: the **Parkfield 2004 Mw 6.0 validation benchmark**
(https://strike.scec.org/cvws/tpv35docs.html, `TPV35_Description_v05.pdf`).
A vertical right-lateral planar fault, 40 km along strike (x in [-30, 10] km)
by 15.5 km deep, reaching the free surface, in a linear elastic medium with a
**different 1D velocity profile on each side of the fault**, linear
slip-weakening friction. Unlike the TPV8/TPV29-style benchmarks, the yield
stress `mu_s(x,z)` and the initial shear stress `tau0(x,z)` are **supplied as
data** on a 100 m grid (Ma, Custodio, Archuleta & Liu, 2008, JGR 113,
B02301), and nucleation is built into that data: around (x = 0, depth 8.1 km)
`mu_s * 60 MPa < tau0`, so the fault fails at t = 0 without any artificial
nucleation (`par.C_nuclea = 0`).

## Official data, shipped verbatim

All four files of `tpv35_data_files.zip` (sha256
`c8d4c28df326ad24a8222cb832ec8b6d38053dafdd5e288b8f2529663852c2f3`, fetched
2026-10-05 from `tpv35docs.html -> download/`) ship in this directory
unmodified, and `tpv35Tools.py` reads them at `case.setup` time:

| file | content | used for |
|---|---|---|
| `tpv35_input_data.txt` | 401 x 156 nodes at 100 m, `nx ny x y mu_s tau0(MPa)` | `on_fault_vars` slots 1 (mu_s) and 8 (tau0, Pa) |
| `tpv35_velocity_structure_near_side.txt` | 8 layers + halfspace (thick Vp Vs rho), SCEC z < 0 | `par.mat` rows with side -1 |
| `tpv35_velocity_structure_far_side.txt` | 9 layers + halfspace, SCEC z > 0 | `par.mat` rows with side +1 |
| `tpv35_station_locations.txt` | 43 surface stations `name code z x` | `par.st_coor_off_fault` |

Frame: `x_eq = x_tpv`, `y_eq = z_tpv` (fault-normal; near side is y < 0),
`z_eq = -y_tpv` -- the same proper rotation `test.tpv29` uses.

**Decimated, never interpolated**. The
case's fault grid must be an exact integer-stride subsample of the official
grid, so `par.dx` must be a multiple of 100 m that also tiles 40 km x
15.5 km: that is **100 m (spec) or 500 m (gate)** and nothing else;
`lib.requireFaultGeometryResolution` / `tpv35Tools.faultGridForCase` refuse
anything else with the reason. The spec allows interpolation for off-grid
nodes; this case never has any. The spec's optional 50 m resolution would
need a finer official grid, which SCEC does not supply.

## Two-sided material: `par.n2mat == 5`

This case introduced the `nmat > 1, n2mat == 5` material table, columns
`[bottom depth (m), Vp, Vs, rho, side]` with `side = -1` for element centres
at `y < fymin` of fault 1 and `+1` for `y > fymin`. Within a side the rule
is the existing 1D one (ascending bottoms, first match). It is accepted only
for a single vertical planar fault (`ntotft == 1`, `C_degen == 0`); the
solver refuses anything else (`ERR_CFG_MATERIAL_TABLE_INVALID`). The
existing `n2mat == 4` path is untouched.

## Gate configuration (what the committed reference is)

`par.dx = 500` m, `par.term` = the sweep's gate term (5 s), `dt =
0.5*dx/7300`, 4 ranks `(nx, ny, nz) = (2, 1, 2)`, `friclaw = 1`, `tpv = 35`.
At 500 m the nucleation patch holds 7 fault nodes (x in {-500, 0, 500},
depth in {7.5, 8.0, 8.5} km). Side/bottom fault borders carry `mu_s = 1000`
(spec: "slip is zero on the borders"), as tpv8/tpv29 do.

The gate is a coarse REGRESSION check, not a spec-accuracy
claim: at 500 m the minimum Vs (1100 m/s) is resolved by ~3 nodes per
wavelength at 1 Hz and the slip-weakening cohesive zone is under one
element. The spec-resolution tier (100 m, 18 s) is recorded in
`testsys/e2e/full_specs.py` and is a scheduling decision.

## Validation against recordings

TPV35 is a validation benchmark: SCEC links NGA-West2 and Ma et al. (2008)
recordings for comparison. No such comparison exists here yet; per the
project's own pattern (TPV29) it belongs in a committed
evidence script after gating, and is NOT a gate.
