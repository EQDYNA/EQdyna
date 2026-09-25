# Row 120: python-jax-mpi station output

## Head SHA at start
9ea1d0719f54023598401201e39c9ee853326af0 (origin/master, includes row 94 / PR #30, 7d8b498)

## Measurement: where tpv8's 4-rank (2,2,1) split lands, and station 11

Built a scratch serial case for test.tpv8 (case.setup, par.nx=ny=nz=1) and
inspected `meshgen.build_grid_lines` + `meshgen.local_line` directly:

```
global nx,ny,nz 97 52 49
xline[0],[-1] -26473.00 .. 26473.00
yline[0],[-1] -14772.33 .. 17292.21
DECOMP[4] = (2, 2, 1)
mex=0 local nx=49 off=0  first=-26473.00 last=0.0
mex=1 local nx=49 off=48 first=0.0       last=26473.00
mey=0 local ny=26 off=0  first=-14772.33 last=500.0
mey=1 local ny=27 off=25 first=500.0     last=17292.21
```

`x=0` and `y=+500 m` are EXACT shared rank-boundary plane coordinates (both
neighbours' local line ends land on the same double, bit for bit -- a plain
slice of the same global array, `meshgen.local_line`).

`par.st_coor_off_fault` station 11 (1-indexed) = `[0, 0.5, -0.3]` km ->
y=+500 m, exactly the mey=0/mey=1 boundary.

## Ownership decision: station 11 (and the analogous x=0 collisions)

**Decision: bug-compatible with Fortran (rule 23), not fixed and not
"rescued" on the python side.**

`src/fortran/meshgen.f90`'s `setSurfaceStation` has edge branches ONLY for
`ix` (Part2 `ix==1`, Part3 `ix==nx`); every branch, including those two,
requires `iy>1 .and. iy<nodeXyzIndex(5)` (strictly interior in y). There is
no `iy==1`/`iy==ny` branch anywhere in the subroutine. A rank-local mesh
therefore drops any off-fault station whose y sits exactly on ITS OWN
y-extent's edge -- true of BOTH ranks sharing a y-boundary plane, since the
shared node is the local edge for each of them. Station 11 (y=+500 m) is
exactly such a node at tpv8's (2,2,1) split: dropped by every rank, on both
the Fortran binary and this port, once the port matches per rank instead of
per global domain.

The pre-existing serial-Python docstring
(`eqdyna3d.report_dropped_stations`, row 94 audit) already named this
exact gap before this mission touched anything: "the Fortran side's tpv8
station 11 is dropped by a pre-existing y-partition-boundary gap in
setSurfaceStation... the serial Python port carries the identical un-gated
iy==0/iy==ny-1 case in principle, just unmasked at different domain sizes
than that MPI split." This mission's job is to make python-jax-mpi
REPRODUCE that Fortran behaviour, not repair it -- board row 94's own
successor row (per the mission brief) owns fixing it, and only in a way
that changes BOTH implementations together with evidence, not silently on
one side.

**Implementation choice that makes this automatic, no special-casing
needed**: `meshgen.build_station_matching` (the existing serial port of
`setSurfaceStation`/`createMasterNode`) already has NO y-edge branch either
-- it was written faithfully off the Fortran already. Calling it with
THIS RANK'S LOCAL `xline`/`yline` (instead of the global lines) reproduces
Fortran's per-rank behaviour, warts included, with zero new branch logic.
Same reasoning explains the OTHER collision this measurement found: `x=0`
is also an exact rank boundary (mex=0/mex=1), and `ix` DOES have edge
branches in both Fortran and the port, so a station or fault node at x=0
CAN be matched by both x-neighbours and both write the same filename --
the identical "last rank to call wins" shape `testsys/matrix.py`'s
`GATE_STATIONS` docstring already documents for test.drv.a6/tpv104/
tpv1053d's Fortran runs. Not a new hazard: Fortran's own 4-rank tpv8
reference (the one already committed under
`test.reference.results/test.tpv8/stations/`) was produced under the exact
same race. This port reproduces the selection rule, not a specific winner.

## Files touched (by line range at this checkpoint)

- `src/python/eqdyna/meshgen.py`
  - `build_station_matching` (~1511-1626): added optional `zline_global`
    kwarg (defaults to `zline`, a no-op for the serial caller and for any
    npz==1 run) so the off-fault depth-band snap can use the GLOBAL z line
    (matching Fortran's own `zGridFull` reconstruction) while the x/y loop
    still walks THIS RANK'S local lines for node numbering.
- `src/python/eqdyna/eqdyna3d.py`
  - `build_solver_state` (~389-465): station matching now runs on BOTH
    paths (serial: global lines; MPI: this rank's local lines, called
    AFTER `part.slice_lines`), replacing the `part is not None` branch that
    used to leave `st_on_idx`/`st_off_idx` empty. `st_off_actual_m` (row 94
    header stamp) computed uniformly for both paths once `meshCoor` exists.
    `report_dropped_stations` stays serial-only (a per-rank call would
    misreport a station matched by ANOTHER rank as globally dropped;
    documented as a scope-cut, not silently omitted).
  - `run_case_mpi` (~681-760): removed the stderr "writes none" warning;
    added a `write stations` profile phase calling
    `library_output.write_onfault_stations`/`write_offfault_stations`
    (S's per-rank `st_on_idx`/`st_off_idx`), BEFORE the `n_own==0` early
    return -- an off-fault station can sit in a box that never touches the
    fault at all (test.tpv8's mey=1 boxes), so gating station writes behind
    frt's fault-ownership test would silently drop them.
- `src/python/eqdyna/driver.py`
  - `run_mpi` return dict (~door): now returns `on_st_hist`/`off_st_hist`
    (previously computed into the carry tuple and discarded on return).
    Docstring comments near `n_on_st_l`/`n_off_st_l` and near the return
    updated from "does not port station output" to the per-rank reality.
- `testsys/matrix.py`
  - `ARTIFACTS['python-jax-mpi']`: `('frt',)` -> `('frt', 'nsign',
    'station')`.
  - Removed `STATION_ARTIFACT_UNSUPPORTED_REASON` (dict, its
    `coverage_report` sub_gaps block, and its two consistency asserts).

## Body-file counts, tpv8, at the gate term

(measured below once the e2e run completes -- see the Gate section of the
final report for the actual numbers: fortran / python-jax-mpi / python-jax)

## Per-step cost delta (rule 4e placement)

(measured once the mpi e2e run completes; the new work is a python-side
`meshgen.build_station_matching` scalar loop over ~this rank's node count,
run ONCE in `build_solver_state`'s setup phase, and a per-step `xp.stack`/
`B.setat` already present in `make_step_parts` for ANY nonzero
`st_on_idx`/`st_off_idx` width -- this landing is what makes that width
nonzero on the MPI path for the first time, not new code in the per-step
loop itself.)
