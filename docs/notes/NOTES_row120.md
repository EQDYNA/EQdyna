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

Measured from the actual `run_e2e.py --cases test.tpv8 --backends
fortran,python-jax,python-jax-mpi --jobs 1` run (2026-09-25):

| backend         | body*.txt | faultst*.txt | total |
|-----------------|-----------|---------------|-------|
| fortran (4 rank)|        14 |             8 |    22 |
| python-jax-mpi  |        14 |             8 |    22 |
| python-jax (serial) |    15 |             8 |    23 |

python-jax-mpi's body count matches Fortran's EXACTLY (14, not the serial
port's 15) -- `body005st000dp003.txt` (station 11) exists ONLY under
`test/test.tpv8.python-jax/`, absent from both `test/test.tpv8/` (fortran)
and `test/test.tpv8.python-jax-mpi/` -- direct confirmation that the MPI
port reproduces Fortran's y-partition-boundary drop, not the serial port's
(correct, for a single global domain) full coverage.

## Per-step cost delta (rule 4e placement)

Measured directly (not inferred), same case dir
(`test/test.tpv8.python-jax-mpi`, tpv8's real mesh, 114 steps, 4 ranks,
`EQDYNA_MPI_SYNC=halo`, `JAX_PLATFORMS=cpu`), same harness script called
before and after, back to back, on the shared box (load ~45/64):

  - BEFORE (`build_solver_state` monkeypatched to force
    `st_on_idx`/`st_off_idx` empty on the MPI path, i.e. this landing
    reverted in place): rank ms/step = 69.9036 / 69.9166 / 69.9065 / 69.9069
    (mean 69.908).
  - AFTER (this landing, unmodified): rank ms/step = 74.7412 / 74.7721 /
    74.7493 / 74.7438 (mean 74.752).
  - Delta: **+4.84 ms/step, +6.9%**, uniform across all 4 ranks regardless
    of which side of the fault they carry stations for (ranks 0/2 carry
    only the 3 on-fault matches, ranks 1/3 only the 2 off-fault matches) --
    consistent with CPU jax's per-op DISPATCH overhead dominating over the
    (tiny, 2-5 station) array width, not a cost that scales with station
    count.
  - PLACEMENT: the added work (the on-fault block's `fric[xh, SLOT]`
    gathers + `xp.stack` + `B.setat`, and the off-fault block's analogous
    `dispArr`/`velArr` gathers) sits INSIDE `make_step_parts`'s `part_a`/
    `part_b`, which is jitted together with the element kernel, hourglass
    resistance, faulting and the mass divide -- the SAME jit `driver.run_mpi`
    already builds. It is not a new collective, not a new host round-trip,
    and not a new bucket: like `fault` already does on jax (profile_emit.py's
    documented fold), this cost lands in the `element` bucket
    (`compute_s`), not `io`/`exchange`/`wait`/`setup`. No new per-step
    collective was added (requirement 2): `on_st_hist`/`off_st_hist` are
    written into the RANK-LOCAL carry by the EXISTING `make_step_parts`
    (shared with the serial `run`), and returned once, at the end of the
    loop, by `run_mpi`'s own return dict -- nothing new crosses an MPI
    boundary per step.
  - Baseline/after scripts: `row120_perf_before.py` / `row120_perf_after.py`
    (scratch, not committed -- see this note for the exact monkeypatch and
    invocation, reproducible from any checkout of this branch).

## Gate (final, this checkpoint)

`python3 testsys/run.py unit regression`: SUCCESS unit (exit 0), SUCCESS
regression (exit 0) -- includes the new
`test_row120_mpi_station_output.py` (shard 2), 6 checks incl. 2 mutations,
0 failures.

`python3 testsys/e2e/run_e2e.py --cases test.tpv8 --backends
fortran,python-jax,python-jax-mpi --jobs 1`:
```
test.tpv8        python-jax    SUCCESS      27.1  max|diff|=1.983643e-10 ...
test.tpv8        python-jax-mpi SUCCESS      17.8  max|diff|=1.220703e-10 ...
test.tpv8        fortran       SUCCESS      13.0  max|diff|=3.051760e-11 ...
ran       : 3 of 33 cells in the 11 case x 3 backend table (3 passed, 0 failed)
```
python-jax-mpi's station line:
```
station: 3 on-fault + 2 off-fault file(s), worst e=1.1176e-10 bound=1.0e-07 (on v-slip-rate at faultst000dp120.txt)
```
(all 13 columns individually `ok`, e in [0, 1.12e-10], every S_q either a
real physical scale or FLOOR-clamped at 1e-6 -- see the full per-column
table in the run log.) `git status test.reference.results/` is clean (no
reference regenerated, rule 7).

## Head SHA at this checkpoint

55dd4b0 (mira/row120-jaxmpi-stations)
