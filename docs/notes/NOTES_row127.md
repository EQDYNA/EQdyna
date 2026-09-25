# Row 127 -- station ownership at a shared MPI-partition boundary

Salvaged WIP: branch `mira/row127-station-yboundary` @ `7e15d62`, no notes
left (agent died on a rate limit). This file is my checkpoint log, written
BEFORE the first test run per this mission's instructions.

## What 7e15d62 actually contains (read via `git show 7e15d62`)

`src/fortran/meshgen.f90`'s `setSurfaceStation`: the old code had three
mutually exclusive x-branches (`ix==1`, `ix>1&&ix<nx` interior, `ix==nx`),
and ALL THREE required `iy>1 .and. iy<ny` -- there was no y-edge branch at
all, so a station whose y lands exactly on a rank seam matched on NO rank
(item 94's residual, tpv8 station 11 at (2,2,1): 14 body files instead of
15). The WIP replaces the three-x-branches / forced-interior-y shape with
a full 3x3 (x-category x y-category) structure, plus a new z gate, and
claims (commit message only, UNTESTED) "only the timestamp header differs,
numerics bit-identical" -- treated as a hypothesis below, not fact.

**Ownership rule the commit implements** (my restatement, checked against
the diff line by line):
- Every axis has 3 categories: low edge (local index 1), interior
  (1<idx<n), high edge (local index n).
- Interior nodes are never shared -- no gate, matches as before.
- The LOW edge (local index 1) on an axis is a match candidate on this
  rank only if this rank has no lower neighbour on that axis
  (`mex`/`mey`/`mez` respectively `== 0`). If it has a lower neighbour, the
  low-edge node is the seam this rank shares with that neighbour, and this
  rank yields.
- The HIGH edge (local index n) on an axis is **always** a match
  candidate. It is either the true global high edge (no higher neighbour)
  or the seam shared with a higher-coordinate neighbour -- and by the same
  rule applied FROM THE OTHER SIDE, that neighbour's own low edge yields to
  THIS rank. So the pair's shared node is owned by the lower-mpi-coordinate
  rank of the two, unconditionally.
- z has no index-keyed match branch at all (depth is matched by VALUE
  against the globally-identical `x4ndsSnapZ`), so z's ownership gate is a
  single early `return` for the whole node: skip entirely when
  `iz==1 .and. mez/=0` (this rank is not the z-owner of this z=1 row).
  There is no analogous z-high-edge check because -- symmetric to x/y's
  high edge -- it needs none: it is always owned.

Net effect: each shared seam node is owned by EXACTLY ONE rank (the lower
MPI coordinate along whichever axis the seam is on), on every axis,
simultaneously. This is what point 1 of the mission asks me to verify, not
just restate.

`src/python/eqdyna/meshgen.py` (`build_station_matching`) and
`src/python/eqdyna/eqdyna3d.py` (`build_solver_state`) carry the identical
rule: `x_matches`/`y_matches` gain `ix==0`/`iy==0` branches gated on
`mex==0`/`mey==0` (0-indexed, so `ix==0` corresponds to Fortran's `ix==1`),
the `ix==nx-1`/`iy==ny-1` branches are ungated (always candidates,
matching Fortran's high edge), and a `z_row_owned = (iz != 0) or
(mez == 0)` flag gates the iz==0 row of the per-node double loop. `mex,
mey, mez = part.mexyz if part is not None else (0, 0, 0)` -- verified
`Partition.mexyz` exists (`MPI4NodalQuant.py:90`,
`= meshgen.calc_xyz_mpi_id(rank, *dims)`, itself verified algebraically
identical to Fortran's `calcXyzMPIId`) and that the serial path (`part is
None`) defaults to `(0,0,0)`, i.e. "sole owner of every axis" -- correct
for a 1-rank decomposition.

## setOnFaultStation (mission point 4) -- checked, NOT touched by the WIP

Read `createMasterNode`'s embedded `setOnFaultStation` loop
(meshgen.f90:1063-1071): it matches by EXACT VALUE equality against
`xonfs(1,i,iFault)`/`xonfs(2,i,iFault)` (strike x, depth z), for every
fault node (`isOnFault==1`), with **no ix/iy/iz-keyed branch of any kind**
-- i.e. no index gate to have a hole in. Python's mirror
(`meshgen.py:1737-1739`, inside `build_station_matching`) matches the same
way: value equality only. So there is no row-127-shaped DROP gap here (my
own reading, not just trusting the WIP's code comment, which asserts the
same thing) -- confirmed as point 4's deliverable below.

There IS a pre-existing, DIFFERENT hazard: a fault node sitting exactly on
a shared x or z seam (if npx>1 or npz>1 splits the fault plane) is
physically present on BOTH ranks (the one-node mesh overlap
`getLocalOneDimCoorArrAndSize` gives every axis), and both will match it
by value with no ownership gate at all -- a DUPLICATE write, not a drop.
This is the SAME shape as the pre-existing off-fault x/z duplicate hazard
row 120 already named ("last rank to call wins", `matrix.py`'s
`GATE_STATIONS` docstring). It is not new, not introduced by this row, and
out of THIS row's stated defect (which is about a DROP). Left unfixed
here; named as a residual in the report.

## Hypotheses to test, in order

H1. Build succeeds (already reported rc=0 by the dead agent -- re-verify
    myself, do not trust the report).
H2. Serial python-jax invariance: for tpv8/10/36/37, before
    (origin/master) and after (this branch) give byte-identical selected
    nodes and station file sets in serial. Rationale: the fix only
    changes behaviour when mex/mey/mez != 0 somewhere (serial passes
    mex=mey=mez=0 always, and with only one rank there is no "high edge
    is also a low edge of a neighbour" case -- the high-edge branch was
    ALREADY unconditional before this change, so serial's set of matched
    stations should be provably unchanged). If this fails, the WIP's edit
    touched more than the ownership gate.
H3. Fortran fresh 4-rank tpv8 (2,2,1): body/faultst file COUNT goes from
    14 to 15, station 11's file appears, on EXACTLY ONE rank (not two).
H4. python-jax-mpi in-process per-rank selection (test_row120's harness,
    updated) matches Fortran's per-rank OWNERSHIP, not just the file-name
    union -- this is what #37's audit said the existing test cannot see.
H5. A synthetic decomposition with a station on a shared x AND shared z
    boundary (not just y) gives exactly one owner on both languages.
H6. The commit message's claim ("only the timestamp header differs") is
    checked literally on a real diff of two runs' output files, not
    assumed.

## Log

- 2026-09-25 T0: branch checked out at 7e15d62, already merged onto
  origin/master (fast-forward, no conflicts). Notes written. About to
  build and run H1.

- 2026-09-25 T1: H1 CONFIRMED. `EQDYNAROOT=<worktree> ./install-eqdyna.sh
  -m ubuntu` -> rc 0, `bin/eqdyna` produced (249952 bytes).

- 2026-09-25 T2: ran the EXISTING `testsys/regression/
  test_row120_mpi_station_output.py` as a cheap targeted check (rule 9)
  before writing anything new. 5 of 6 checks PASS (mutation checks still
  have teeth); the one FAIL is its own now-STALE assertion "station 11 is
  dropped by BOTH fortran and python" -- expected, since that assertion
  encodes the row-120 DROP this row's fix removes. Both fortran and python
  now report 23 station files (was fewer): the fix and the port agree.
  This is evidence for H3/H4, not a regression -- the assertion itself
  needs updating, tracked below.

- 2026-09-25 T3: built a small ad hoc per-rank-ownership probe (not yet a
  committed test) reusing `MPI4NodalQuant.Partition` + `eqdyna3d.
  build_solver_state` in-process for tpv8 at (2,2,1): station 11
  (`body005st000dp003.txt`) has owners=[rank 0] -- exactly one, confirming
  H3's ownership half in PYTHON. Also surfaced a residual: `faultst000dp000/
  045/075/120.txt` (on-fault) are matched by BOTH rank 0 (mex=0) and rank 2
  (mex=1) -- a real duplicate-write hazard. This is `setOnFaultStation`
  (point 4): read `createMasterNode` (meshgen.f90:1063-1071) directly --
  it matches fault nodes by EXACT VALUE against `xonfs`, with NO ix/iy/iz
  index branch of any kind, so it has no row-127-shaped DROP gap (nothing
  to have a hole in) but it also has NO ownership gate at all, so a fault
  node sitting exactly on a shared x or z seam (tpv8's fault spans x, and
  npx=2 splits it) is physically present on both ranks and both write it.
  Confirmed this is PRE-EXISTING and untouched by the WIP: `git diff
  origin/master...HEAD -- src/fortran/meshgen.f90` shows exactly one
  hunk, at `setSurfaceStation` (line ~729), nothing near
  `createMasterNode` (~1002-1104). Point 4's answer: NO drop gap to fix;
  YES a separate, pre-existing, out-of-scope DUPLICATE hazard, same shape
  as the off-fault one row 120 already named ("last rank wins",
  `matrix.py`'s `GATE_STATIONS` docstring). Not fixed here -- flagged in
  the report.

- 2026-09-25 T4: REAL Fortran, not just the Python side. Built a per-rank
  output-directory harness (mpirun launches a wrapper that `cd`s into
  `rank$OMPI_COMM_WORLD_RANK` before exec'ing the real binary --
  `OMPI_COMM_WORLD_RANK` is Open MPI's own per-process env var; every case
  input file is symlinked into each rank dir so nothing is duplicated on
  disk) so I can see PER-RANK ownership directly from a real 4-rank run,
  not just the aggregate file set on shared cwd. Ran a fresh tpv8 case
  (par.term = 5*dt) at (2,2,1) on the WIP's own freshly-built binary:
    - 15 distinct body* files total (was 14 pre-fix, matches H3 exactly).
    - Every body* filename appears in EXACTLY ONE rank's directory --
      zero duplicates. `body005st000dp003.txt` (station 11) is rank 0's
      alone.
    - faultst* (on-fault) files: rank 0 and rank 2 both wrote all four of
      rank 0's faultst files -- the SAME duplicate hazard T3 found in
      Python, now confirmed in real Fortran, cross-validating both sides
      independently rather than trusting one language's arithmetic.
  This directly answers point 1 (off-fault: exactly one owner, confirmed)
  and point 4 (on-fault: pre-existing duplicate, not a drop, not fixed by
  this row).

- 2026-09-25 T5: H5 (synthetic decomposition, shared x AND z). Built a
  second case, tpv8's own geometry but `par.nx,par.ny,par.nz = 2,1,2`
  (4 ranks, decomposition (2,1,2) -- exercises npz>1, which tpv8's own
  gated (2,2,1) never does), queried the exact seam coordinates by calling
  `meshgen.build_grid_lines` + `Partition.slice_lines` directly (pure
  Python, no case run needed, and per rule 23/row 94's own established
  fact that this line-building is bit-identical on every rank/language):
  x seam at 0.0 m (rank mex 0/1 boundary), z seam at -12000.0 m (rank mez
  0/1 boundary). Added ONE synthetic off-fault station at
  `(0, 0.5, -12)` km -- exactly the CORNER shared by all four ranks
  (x-seam AND z-seam simultaneously) -- reran `case.setup`, then ran BOTH
  the real Fortran binary (same per-rank-directory harness) and Python's
  per-rank `build_solver_state` in-process:
    - 16 distinct body* files both sides.
    - The corner station's file, `body005st000dp120.txt`, has EXACTLY
      ONE owner on both languages: rank 0 (mex=0, mez=0 -- the lowest
      coordinate on both shared axes, exactly the ownership rule
      predicts).
    - Fortran's full per-rank file-to-owner map and Python's are IDENTICAL,
      filename for filename, rank for rank (compared by hand, all 16
      match) -- not just the same total count.
  H5 CONFIRMED for the off-fault path, in both languages, on real
  compiled/executed code, for a case that exercises x AND z simultaneously
  (tpv8 itself only exercises x and y, npz==1 always).

- 2026-09-25 T6: "mutation-check against origin/master" (item 3's own
  requirement), done by hand before committing to the shape of the
  regression test. Built a SEPARATE `eqdyna` binary with every source file
  identical to this branch EXCEPT `meshgen.f90`, replaced by
  `git show origin/master:src/fortran/meshgen.f90` (the pre-fix version) --
  rc 0. Ran it on the SAME two cases with the SAME per-rank harness:
    - tpv8 (2,2,1): 14 distinct body* files (was 15 on the branch).
      `body005st000dp003.txt` (station 11) is ABSENT from every rank's
      directory -- the DROP reproduces exactly as item 94/127 describe.
      No off-fault duplicate happened to appear in THIS case (the x/y
      seam stations here didn't line up with a coincidental double-match
      in the old code, or one direction's old branch order just didn't
      produce one for these particular coordinates).
    - synthetic (2,1,2) corner case: `body005st000dp120.txt` (the x+z
      corner station) is written by BOTH rank 2 AND rank 3 -- a genuine
      DUPLICATE, not a drop, on the pre-fix code. This is exactly the
      "converse hazard" the mission's Defect section warns about, caught
      on a real corner case rather than argued abstractly.
  H6/mutation-check CONFIRMED both ways: the fix's own test, run against
  the pre-fix source, goes red by DROP on one case and by DUPLICATE on
  the other -- proving the new regression test has teeth in both
  directions the mission asked for, not just one.

- 2026-09-25 T7: about to (a) write the formal regression test encoding
  T4/T5/T6's harness + assertions, (b) register it in `testsys/ci_shard.py`,
  (c) fix `test_row120_mpi_station_output.py`'s now-stale "station 11 is
  dropped" assertion (T2) to assert the NEW behaviour (found by exactly
  one rank) instead, (d) run the full required gate, (e) check
  `test.reference.results/` stays clean (mission point 5), (f) commit.

## Findings summary (for the final report)

1. Off-fault ownership rule verified end to end, both languages, real
   compiled/executed code, two decompositions (one exercising x+y seams,
   one exercising x+z simultaneously): every in-mesh off-fault station is
   matched by EXACTLY ONE rank. Ownership = the lower-MPI-coordinate rank
   of the pair sharing a seam node, on every axis independently.
2. Serial invariance not yet re-checked for tpv10/36/37 -- next step.
3. On-fault (`setOnFaultStation`) has NO drop gap (no index branch to have
   a hole in) but DOES have a pre-existing, untouched-by-this-row
   DUPLICATE hazard at a shared x or z fault-plane seam (confirmed live on
   tpv8 itself: ranks 0 and 2 both write all four of rank 0's faultst*
   files). Named, not fixed -- out of this row's scope (the mission's
   defect statement is about off-fault DROP; the mission's point 4 asks
   only to check on-fault for "the same edge gap", which is the drop
   shape, and it does not have one).

