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
