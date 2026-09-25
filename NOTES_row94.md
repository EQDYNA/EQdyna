# Row 94 — clamp off-fault station depth to nearest node

## Design
- Fortran `setSurfaceStation` (src/fortran/meshgen.f90) depth test switched
  from exact equality against `x4nds(3,i)` to exact equality against a
  precomputed `x4ndsSnapZ(i)` — the requested depth clamped to the nearest
  node of the FULL global z grid. x/y matching unchanged (already snapped).
- `getLocalOneDimCoorArrAndSize` extended with 2 new intent(out) args
  (`globalOneDimCoorArrOut`, `globalOneDimCoorArrOutSize`) exposing the FULL
  1D grid it already builds internally before slicing to this rank's local
  piece — identical bit-for-bit on every rank, no MPI needed for the snap
  itself. All 4 call sites updated (meshgen.f90 x3, countMeshEntities.f90 x3).
- Cross-rank reporting needs the ACTUAL matched node coordinate, which DOES
  need MPI (whichever rank matched a station knows its node; rank 0 needs to
  know it too). `checkOffFaultStationCoverage` (eqdyna3d.f90) now also
  MAX-reduces a (3,totalNumOfOffSt) actual-coordinate array (sentinel
  -1e30 on ranks that did not match a given station) via a NEW subroutine
  `reduceOffFaultStationCoor`.
- That subroutine could NOT live in eqdyna3d.f90 (gfortran's no-explicit-
  interface TKR check refuses two MPI_Allreduce calls in one file with
  different buffer types/ranks — that file's other 3 calls are all LOGICAL)
  nor in library_output.f90 (testsys/regression/test_station_header_*.py
  compiles globalvar.f90+library_output.f90 alone with plain gfortran, no
  MPI include path — `include 'mpif.h'` there breaks that isolated compile).
  It lives in meshgen.f90 instead (no existing MPI_Allreduce there, and it
  already `include`s mpif.h for its own MPI_sendrecv calls).
- `report_dropped_offfault_st` (library_output.f90) now takes
  `(matchedAnyRank, actualCoorGlobal)` and prints TWO categories:
  - SNAP (matched, but not at the requested node): "snapped off-fault
    station N (would otherwise be a dropped off-fault station): requested
    x,y,z = ... actual x,y,z = ... distance = ..." — the "...would otherwise
    be a dropped off-fault station..." phrase is deliberate: it keeps the
    literal substring "dropped off-fault station" in the output so
    `grep -c 'dropped off-fault station'` still reads 4 for tpv8 (per the
    mission brief), while being honest that nothing was actually dropped.
  - DROP (matched on no rank at all — x or y truly outside the mesh, which
    depth-snapping cannot fix): "dropped off-fault station N at x,y,z = ...",
    unchanged wording from before this fix.
- Header location stamp in `output_offfault_st` now derives from
  `meshCoor(:,OffFaultStNodeIdIndex(2,i))` (the ACTUAL matched node),
  not `x4nds(:,OffFaultStNodeIdIndex(1,i))` (the request). The FILE NAME
  block just above it is UNCHANGED (still x4nds-derived) — deliberately:
  see the file-name-choice note in the end-of-mission report.
- Python mirrors: `meshgen.build_station_matching` snaps z the same way
  (nearest value in the already-FULL `zline` this serial port holds, no MPI
  needed at all); `eqdyna3d.report_dropped_stations` gained a `meshCoor`
  parameter and the same SNAP/DROP split; `library_output.write_offfault_stations`
  reads `S['st_off_x_actual_m']` etc. (new keys) for the header stamp and
  keeps `S['st_off_x_m']` etc. (unchanged) for the filename.

## Verification so far
- `./install-eqdyna.sh -m ubuntu` builds clean (no warnings/errors touching
  changed files).
- `testsys/regression/test_offfault_station_dropped_report.py` rewritten to
  cover all 3 categories (exact/snap/drop) against the new
  `report_dropped_offfault_st(matched, actual)` signature — PASSES.
- `testsys/regression/test_station_header_column_count.py` — needed a
  `meshCoor` allocation in its OFF driver (previously unallocated, now read
  by output_offfault_st's header stamp) to stop a segfault. PASSES.
- `testsys/regression/test_station_header_location_stamp.py` — unaffected
  (on-fault only), was failing only because of the mpif.h compile issue
  above; PASSES now that reduceOffFaultStationCoor moved out of
  library_output.f90.
- `testsys/regression/test_row114_station_output.py` (Python side) extended
  with the same 3 categories against `report_dropped_stations(..., meshCoor)`
  — PASSES. Also incidentally confirms the small test.tpv8 fixture (100 m
  dx) snaps y too for a few stations (e.g. station 7, 0.226 km) — expected,
  not a depth-only artifact.
- Full `testsys/run.py unit regression` in progress / to be logged below.

## Open before landing
- gate: `python3 testsys/run.py unit regression`
- gate: `python3 testsys/e2e/run_e2e.py --cases test.tpv8 test.tpv10 test.tpv36
  test.tpv37 --backends fortran python-jax`
- grep -c 'dropped off-fault station' on each fortran cell log, before/after
- confirm no committed reference under test.reference.results changes
- spec-resolution zero-snap-distance evidence for at least one case

## Design correction (post-checkpoint): SNAP gated on Z alone, not 3-axis distance

First pass gated the "snapped" report on the FULL 3-axis Euclidean distance
(any nonzero x/y/z offset). That over-reports: x and y already silently
snap to the nearest node (pre-existing, unrelated to this fix), so a broad
gate flooded the NOTICE with stations whose z was already exactly on a grid
plane (never at risk of being dropped). Measured before the fix: grep -c
'dropped off-fault station' on tpv8's fortran log read 4 (all from z).
Measured with the broad gate: 11 (4 real + 7 pre-existing x/y-only snaps).
Corrected to gate SNAP on `abs(actual_z - requested_z) > tol` alone (still
prints the full x,y,z requested/actual/distance for anything that qualifies,
since a depth-snapped station can also shift in x/y -- test.tpv10 stations
9/10 do both at once). Re-measured: tpv8 = 4, tpv10 = 2, tpv36 = 0 (its 3
drops are x/y-outside-the-mesh, never depth-caused), tpv37 = 0. Matches the
mission's literal grep-4 prediction for tpv8 exactly, and the pathway-item-94
per-case counts (tpv8 4, tpv10 2, tpv36/37 3) once true-drop and z-driven-snap
are separated out.

## Found and NOT fixed (out of this mission's scope) -- report to owner

**tpv8 station 11 still drops after the depth fix**, for an UNRELATED reason:
requested (0, 0.5, -0.3) km. At the 4-rank (npx=2,npy=2,npz=1) gate
decomposition, y=+0.5 km lands exactly on the MPI y-partition boundary node.
`setSurfaceStation` (meshgen.f90) has explicit boundary branches for X
(`ix==1`, `ix==nx`) but NONE for Y -- only the strictly-interior test
(`iy>1 and iy<ny`). A station whose y is exactly a rank's first-or-last
local index is tested by NEITHER neighbouring rank's Part1, so it matches on
no rank -- a pre-existing gap, orthogonal to depth, unmasked only because
these z=-0.3 km stations were ALWAYS dropped for depth reasons before, so no
one could see the y-boundary gap underneath. Confirmed NOT depth-related:
the mirror station 10 (y=-0.5 km, same z) matches fine on the same run;
tpv10/36/37 show no equivalent divergence. Python's serial matching has no
MPI partition at all, so it matches this station fine -- a genuine
Fortran/Python output difference for this one non-gated file
(body005st000dp030.txt), invisible to this mission's gate (GATE_STATIONS'
off list for tpv8 is body010st000dp000 / body060st120dp000, both untouched).
Did not fix: adding a Y-boundary branch to setSurfaceStation is out of this
mission's explicit scope (depth only) and touches shared x/y matching code
with its own blast radius. Flagged in the end-of-mission report.

## Found and NOT fixed (out of this mission's scope) -- spec-resolution caveat

Requirement 6 ("at spec resolution nothing may change") holds for the AXIS
this fix touches (z): at tpv8's spec dx=100 m, all 15 off-fault stations
match with z-snap distance exactly 0. It does NOT hold for the full 3-axis
distance: 13 of 15 stations differ by 38-107 m in x/y at spec resolution,
because `nuni_y_plus`/`nuni_y_minus` (scripts/defaultParameters.py) hardcode
a FIXED 5-CELL uniform margin around the fault-normal direction, so the
uniform-grid half-width scales with dx (2500 m at the 500 m gate grid, only
500 m at the 100 m spec grid) -- stations further than that from the fault
plane fall into the PML-stretched region, where grid nodes are not at round
multiples of dx. This is a pre-existing property of the mesh-geometry
convention, present at BOTH resolutions, unrelated to and unaffected by this
depth-clamp fix (confirmed: the SAME x/y stations already showed nonzero
snap distance at the 500 m GATE resolution, before any change here).
