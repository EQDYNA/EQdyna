# test.restart2

The mode = 2 (EQsimu restart) end-to-end gate. A 20 km x 4 km
planar vertical strike-slip fault (81 x 17 nodes, dx = 250 m), rate-and-state
aging law (friclaw 3), no artificial nucleation, restarted from a real
EQquasi fault.r.nc. Gated on fortran (4 ranks) and python-jax (serial);
python-jax-mpi refuses mode 2 and is declared UNSUPPORTED in
testsys/matrix.py.

## Fixture provenance

fault.r.nc (67 KB, float32, rank-2, dims nid_dip = 17 x nid_strike = 81) is
a crop of the EQquasi run bp1001.qd.250, snapshot Q3/fault.r.nc (161 x 321
nodes, dx = 250 m, BP5 rate-and-state parameters), dip rows 144-160 (z = -4 km
to the free surface), strike columns 120-200 (x = -10 km to +10 km). No value
is altered; the 12 variables the reader needs are copied as stored and the
three slip variables it does not read are dropped. Regenerate with
testsys/crop_restart_fixture.py (arguments: source file, destination file).

The window sits above the source's locked band (dip rows 100-140, slip rates
1e-28 to 1e-34 m/s), where the aging law's L/V + (theta - L/V) exp(-V dt/L)
cancels to theta = 0 and the solver stops with a NaN (fric.f90,
rate_state_ageing_law).

## What it exercises

Fric slots 8, 49, 7, 47, 20, 23 (stresses, slip rate, state, theta_pc) and
31-36 (node velocities, which seed v1 in init_vel and the python carry) come
from the file. The per-node values in user_defined_params.py are what a
mode = 1 twin of the case starts from. Nothing ruptures in the 5 s gate term
(slip rates 1e-6 to 6e-6 m/s); the gate compares the frt columns (slip,
tractions, peak slip rate, velocities), the on/off-fault station series and
fault.dyna.r.nc.

Source provenance: sha256 of the EQquasi source file (read-only, outside the repo) is
2ccbd164edc58bf3d5e004afc73480c792aad8dff42856777b1ca1f461dd5ceb; the fixture is its
unaltered slice (`testsys/crop_restart_fixture.py`). Negative control (conductor, 2026-11-05):
zeroing vxm/vxs/vym/vys/vzm/vzs in the fixture turns the fortran cell red
(station worst e=3.3e3 vs bound 1e-10), so the gate is sensitive to the velocity seeding.
