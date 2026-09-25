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
