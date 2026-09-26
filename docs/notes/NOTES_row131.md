# Rows 131/132 -- on-fault station ownership duplicate + Fortran axis-too-thin guard

Branch `mira/row131-onfault-owner` @ origin/master 18c5740 (row 127, PR #38).
Written BEFORE the first build/test run, per this mission's discipline.

## Row 131 -- the defect, read directly (not trusted from the mission prose)

`createMasterNode`'s embedded `setOnFaultStation` loop (meshgen.f90, was
~1072-1080 pre-edit) matches a fault node by EXACT VALUE against
`xonfs(1,i,iFault)`/`xonfs(2,i,iFault)` (strike x, depth z) -- confirmed by
direct read, no ix/iy/iz-keyed branch of any kind, same finding
docs/notes/NOTES_row127.md (T3/T4) already recorded as a NAMED, unfixed
residual. `getLocalOneDimCoorArrAndSize` overlaps neighbouring ranks by
exactly one node on every axis (x, y, and z when npz>1) -- established fact
from row 127, re-checked here by reading the subroutine again, not assumed.
So a fault node sitting exactly on a shared x seam (npx>1, tpv8's own case)
or z seam (npz>1) is isOnFault==1 on BOTH ranks and both execute the
value-match loop, both increment `numOfOnFaultStCount`, both write the same
`faultst*` file (last writer on a shared cwd wins; row 127's per-rank-cwd
harness shows this directly as a real double-write, not just an inference).

## The fix

Same ownership rule as row 127: the lower-MPI-coordinate rank of a seam pair
owns the shared node. Applied ONLY to the station-selection sub-block
(the `setOnFaultStation` loop / `anonfs` append), NOT to anything else in
`createMasterNode` -- `nftnd0`, `nsmp`, `msnode`, `fltgm`, `un`/`us`/`ud`,
equation numbers, area calc, all stay unconditional on every rank, because
each rank needs its OWN full split-node pair for the solve regardless of
who writes the station file. This is the mission's own explicit caution
("only which rank WRITES the station file changes") and I have checked by
reading the diff that nothing else in the subroutine is touched.

Fortran (`src/fortran/meshgen.f90`, `createMasterNode`): compute
`mex,mey,mez` via `calcXyzMPIId` once at subroutine entry; gate the
`setOnFaultStation` loop behind
`isOnFaultStationOwner = .not. ((nodeXyzIndex(1)==1 .and. mex/=0) .or. (nodeXyzIndex(3)==1 .and. mez/=0))`.
No y-gate: the mission scopes this to "a shared x OR z seam" and the fault
plane's own y is fixed (checkIsOnFault tests `nodeCoor(2)==0.0d0` exactly),
so y never varies across a fault node the way x/z (strike/depth) do.

Python (`src/python/eqdyna/meshgen.py`, `build_station_matching`): identical
gate on the 0-indexed `ix`/`iz` loop variables, right where `anonfs.append`
happens, using the SAME `mex`/`mez` already threaded into this function for
the off-fault (row 127) gates. `fault_seq` (Fortran's `nftnd0`) still
increments unconditionally, matching the Fortran change exactly.

## Row 132(1) -- Fortran axis-too-thin guard

`setSurfaceStation`'s x/y branches are `elseif`-chained: `1<ix<nx`,
`ix==1`, `ix==nx`, mutually exclusive. If a rank's local nx==1 (its axis
holds fewer than 2 nodes -- possible when npx/npy/npz is large relative to
the grid), `ix==1` is tested first and always wins, so `ix==nx` (which per
row 127's rule is ALWAYS an ownership candidate) is unreachable for that
rank -- a silent zero-owner drop, not a crash. Python's `check_partition_1d`
already refuses this shape loudly. Fix: guard in
`getLocalOneDimCoorArrAndSize`, right where `localOneDimCoorArrSize` is
first computed for THIS rank's axis -- `localOneDimCoorArrSize < 2` calls
`abortRun(ERR_MPI_AXIS_TOO_THIN, ...)`. New code 50 (51-59 MPI/domain-
decomposition range, decade-0 slot, per errorCodes.f90's own convention).
Every rank computes its own local size, so whichever rank(s) are too thin
abort and MPI_Abort tears down the whole job -- equivalent in effect to a
global check without needing one.

## Row 132(2) -- mutation-check tautology in test_row127_station_ownership.py

The DROP-shape mutation sub-check (~374-378 in the row-127 commit) never
calls `_ownership_violations` -- it computes `expected_names - observed_names`
by hand, so it is testing Python set arithmetic, not the guard function. Fix:
extend `_ownership_violations` with an optional `expected` set (any name in
`expected` that owns zero ranks is reported, as an empty-owner-list entry)
and rewrite the drop sub-check to call it.

## Hypotheses to test, in order

H1. Build succeeds (fresh, not trusting row 127's old build).
H2. Fortran, tpv8 (2,2,1), real per-rank harness: faultst* file count and
    per-rank owners -- BEFORE this edit (already recorded once in
    NOTES_row127.md T4: ranks 0 and 2 both write all four of rank 0's
    faultst files) vs AFTER (expect: each faultst* file exactly one owner,
    same 4-file count, no drop).
H3. Python, same case, same decomposition, in-process per-rank selection:
    same before/after comparison, cross-checked against Fortran's per-rank
    map (not just the union).
H4. frt output (canonical rupture-time field) is UNCHANGED end to end --
    the physics (nftnd0/nsmp/msnode/un/us/ud) must not move at all. Checked
    two ways: (a) the full e2e gate's frt comparison against the committed
    reference, (b) a direct byte diff of every `faultst*` file's NUMERIC
    content (header date line excluded) between this branch and master, on
    a SERIAL (1-rank) run where the ownership gate is a no-op (mex=mey=mez=0
    always) -- proving the gate changes nothing when there is only one rank
    -- and between two MULTI-rank runs' surviving (single-owner) files,
    proving the retained file's own numbers are identical to what the
    corresponding rank produced before (only the DUPLICATE write is
    removed, not the content of the real one).
H5. row132(1): a synthetic decomposition that starves one rank's local axis
    below 2 nodes triggers `ERR_MPI_AXIS_TOO_THIN` (Fortran) and
    `check_partition_1d`'s ValueError (Python), both loudly, neither
    silently continuing.
H6. row132(2): the DROP-shape mutation sub-check in
    test_row127_station_ownership.py, rewritten to call
    `_ownership_violations` (extended with an `expected` set so a 0-owner
    name is observable through the same function, not a hand-rolled
    set-difference), still reports the drop -- and the pre-existing
    duplicate-shape check still reports the duplicate -- so both directions
    keep teeth through the SAME function.
H7. Full gate (unit+regression, e2e fortran+python-jax on 4 cases,
    python-jax-mpi on tpv8) passes; `git status test.reference.results/`
    stays clean (no committed reference should move -- this is a station
    OWNERSHIP change, not a physics change).

## Log

- T0: build (H1) CONFIRMED. `install-eqdyna.sh -m ubuntu` -> rc 0, `bin/eqdyna`
  produced.

- T1: Fortran fix applied to `createMasterNode` (meshgen.f90): `mex,mey,mez`
  computed once via `calcXyzMPIId`; the `setOnFaultStation` loop (only) gated
  behind `isOnFaultStationOwner = .not. ((ix==1 .and. mex/=0) .or. (iz==1
  .and. mez/=0))`. Nothing else in the subroutine touched (checked by re-
  reading the diff: nftnd0/nsmp/msnode/fltgm/un/us/ud/equation numbers/area
  calc all still unconditional).

- T2: Python mirror applied to `build_station_matching` (meshgen.py): same
  gate on 0-indexed ix/iz, right at `anonfs.append`; `fault_seq` still
  increments unconditionally (matches Fortran's nftnd0 change exactly).

- T3: H2/H3 (per-rank ownership, both languages) CONFIRMED via an ad hoc
  in-process probe (test_row127_station_ownership's own build_eqdyna/
  make_case/run_fortran_per_rank/python_per_rank/_ownership_violations,
  imported, not reimplemented): tpv8 (2,2,1) -- zero on-fault violations in
  BOTH languages; 8 distinct faultst* files total, split rank0={dp000,dp045,
  dp075,dp120}, rank2={dp045(120?),...} -- exactly one owner per file, no
  drop (still 8 files, same as before the fix, just no longer duplicated).
  Synthetic corner (2,1,2) -- also zero violations both languages, including
  the 4-way corner case (faultst000dp120.txt, shared by all 4 ranks under
  the OLD code -- see T5 below).

- T4: row 132(1) (Fortran axis-too-thin guard). Read setSurfaceStation's
  elseif chain again to confirm the exact failure mode (ix==1 always wins
  over ix==nx when nx==1). Searched for a REAL case (tpv8, coarsened dx/
  nPML/dis4uniF/dis4uniB via case.setup overrides) that could exercise this
  within the shared box's 8-rank cap -- measured floors of 17 (x), 25 (y,
  unaffected by dis4uniF/B or nPML overrides -- case.setup evidently does
  not wire those through to what build_grid_lines reads; nPML is a
  globalvar.f90 HARD DEFAULT of 6 regardless), 9 (z), all >8 -- no real
  case run fits the box. Used a standalone probe instead (real production
  `getLocalOneDimCoorArrAndSize`, hand-picked but self-consistent inputs:
  dx=1, fault span=1, xmin=-2, xmax=2, rat=1.025, nPML=1, npx=8 -> global
  line = 7 nodes / 8 ranks -> ranks 0 AND 1 both size 1; rank 1 is NOT an
  edge (0 and 7 are). Added `ERR_MPI_AXIS_TOO_THIN` (errorCodes.f90; first
  tried code 50, WRONG -- the test_stop_exit_status.py RANGES table has
  gaps at the x0 slot of every decade (40, 50, 60...) that are NOT part of
  the labelled range and so never render in docs/user/troubleshooting.md;
  the decade-0 slot is for when a decade FILLS UP, not a default "put new
  codes here" location -- corrected to 53, next free in the 51-59 MPI/
  decomposition range which had room. PAPERCUT candidate: the errorCodes.f90
  header comment ("decade 0 is left free... so a related code can be added
  next to its neighbours") reads as "always free to use" on a first pass;
  it means "spare capacity for when the labelled range is full", and the
  render-table silently drops anything in the gap instead of erroring --
  a code registered outside every RANGES tuple should fail loudly, not
  vanish from the docs). Guard added in getLocalOneDimCoorArrAndSize right
  where localOneDimCoorArrSize is first computed. `test_stop_exit_status.py
  --update` regenerated docs/user/troubleshooting.md; new committed test
  `test_row132_axis_too_thin.py` (registered in testsys/ci_shard.py) -- 5/5
  PASS: Fortran probe exits 53 with a FATAL block naming the reason;
  Python's `check_partition_1d(7, 8)` raises citing "rank 1 holds 1
  node(s)" (the middle rank) explicitly, not only rank 0.

- T5: row 132(2). Extended `_ownership_violations` with an optional
  `expected` set (a 0-owner name in `expected` is reported as `name: []`).
  Rewrote the DROP mutation sub-check to call it instead of a hand-rolled
  set difference. Rewrote the on-fault check from "hazard reproduces
  IDENTICALLY in both languages" (which, after the row 131 fix, PASSES
  VACUOUSLY -- both sides now agree by both being EMPTY, not by both
  reproducing the hazard) to "zero on-fault ownership violations" in each
  language separately, plus an explicit per-rank-map equality check.
  test_row127_station_ownership.py: 25/25 PASS.

- T6: mutation-check against origin/master's CODE (not just a hand
  argument) for the row 131 on-fault fix, same method row127's own T6d
  used: swapped this branch's meshgen.f90 for `git show
  origin/master:src/fortran/meshgen.f90`, rebuilt (all OTHER files
  unchanged), ran the SAME test module's own build_eqdyna/make_case/
  run_fortran_per_rank/python_per_rank/_ownership_violations (imported,
  not reimplemented) plus a `git archive origin/master -- src/python`
  checkout for the Python side. RESULT: on unfixed code, BOTH languages
  reproduce the duplicate identically --
    tpv8 (2,2,1): fortran/python on-fault violations both
      {faultst000dp045/075/120/000.txt: [0, 2]} (ranks 0 and 2 both write
      all 4 of rank 0's files -- matches NOTES_row127.md's T4 exactly).
    synthetic corner (2,1,2): fortran/python on-fault violations both
      {faultst000dp120.txt: [0,1,2,3], faultst000dp045/075/000.txt: [1,3]}
      -- the corner station (on BOTH the x and z seam at once) is
      duplicated across ALL FOUR ranks on unfixed code.
  Restored this branch's meshgen.f90 immediately after capturing the
  master binary (`diff` confirmed byte-identical restore before
  rebuilding); binaries/checkouts used only in /tmp, cleaned up after.
  CONFIRMED: this row's fix closes a real, reproduced-on-demand defect in
  BOTH languages, and the updated test (T5) would have caught it (its
  "zero on-fault ownership violations" assertions are exactly what FAILS
  on this master-meshgen binary/checkout).

- T7: H4 (frt/physics unchanged), direct evidence beyond the e2e gate.
  (a) SERIAL invariance: built a real tpv8 case forced to (1,1,1) (explicit
  `par.nx,par.ny,par.nz=1,1,1` -- tpv8's own compset default is (2,2,1), not
  serial, checked by reading the failure when this was omitted: "MPI ranks
  (npx*npy*npz): 4" then an MPI_Sendrecv rank-out-of-range abort). Ran this
  branch's `bin/eqdyna` and a freshly rebuilt origin/master-meshgen.f90
  binary (same swap-build-restore procedure as T6) via `mpirun -np 1`, both
  rc=0. Same 23-file set (matches row127's own T6c serial count); for every
  file, `diff <(grep -v '# date' branch/f) <(grep -v '# date' master/f)` is
  EMPTY -- 0 of 23 differ numerically. Confirms the ownership gate is a
  true no-op in serial (mex=mey=mez=0 always -> the new `.or.` condition is
  always False), exactly as the gate's own logic predicts.
  (b) MULTI-RANK retained-file content: real tpv8 (2,2,1), 4-rank run,
  per-rank-cwd harness, branch binary THEN (same case dir, rank subdirs
  cleared of real files, symlinks kept) origin/master-meshgen binary.
  Branch rank 0 and master rank 0 write the IDENTICAL 4 faultst* filenames;
  every one is numerically IDENTICAL (date line excluded) between the two
  binaries -- 0 of 4 differ. This is the direct proof the fix removes only
  the DUPLICATE WRITE (on rank 2, which master's rank 2 still produces and
  branch's does not) and does not alter what the surviving, single, correct
  writer puts in the file.
  Both master binaries and temp case dirs removed after use (/tmp only,
  nothing left in the worktree).

- T8: registered `test_row132_axis_too_thin.py` in `testsys/ci_shard.py`
  (`ci_shard.py verify`: 81 on-disk, 81 assigned, 0 missing/stale).

- T9: FULL REQUIRED GATE, fresh:
    - `python3 testsys/run.py unit regression` -> SUCCESS both tiers (exit 0).
    - `python3 testsys/e2e/run_e2e.py --cases test.tpv8,test.tpv10,test.tpv36,
      test.tpv37 --backends fortran,python-jax --jobs 1` -> 8/8 cells SUCCESS
      (wall clock 878.3s); every frt/nc/station comparison against the
      COMMITTED reference passed at its case's own bound (e.g. tpv8 fortran
      max|diff|=3.05e-11 vs bound 1e-08; tpv36/tpv37 fortran EXACT 0.0e+00).
    - `python3 testsys/e2e/run_e2e.py --cases test.tpv8 --backends
      python-jax-mpi --jobs 1` -> 1/1 cell SUCCESS (14.9s; max|diff|=1.22e-10
      vs bound 1e-08; station worst e=1.12e-10 vs bound 1e-07).
    - `git status test.reference.results/` -> clean, nothing to stop for
      (row 131/132 are ownership/refusal changes, no reference moved).
  Reverted the gate's own scratch (`docs/perf_ledger.jsonl`,
  `docs/run_profiles.jsonl`, two new `docs/perf_snapshots/*.json`) before
  committing, per this mission's instruction.

## Findings summary (for the final report)

1. Row 131 (on-fault station duplicate): FIXED in both languages, same
   ownership rule as row 127, gating ONLY the station-selection sub-block.
   Mutation-checked against origin/master (both languages): reproduces
   identically on unfixed code (tpv8 (2,2,1): ranks 0/2 both write all 4 of
   rank 0's faultst files; synthetic (2,1,2) corner: 4-way duplicate).
   frt/physics proven unchanged: serial invariance (0/23 files differ,
   date header excluded) and multi-rank retained-file content identity
   (rank 0's own 4 faultst files are byte-identical branch vs master, date
   header excluded).
2. Row 132(1) (Fortran axis-too-thin guard): FIXED via new
   `ERR_MPI_AXIS_TOO_THIN=53` and a guard in `getLocalOneDimCoorArrAndSize`.
   No real case in this repo fits the shared box's 8-rank cap while
   triggering the shape (measured floors: tpv8 x=17,y=25,z=9, all >8) --
   tested via a standalone probe that calls the real production routine
   with a hand-derived, real, self-consistent decomposition (7 nodes / 8
   ranks; rank 1, a MIDDLE rank, goes thin). Python's `check_partition_1d`
   already refused this loudly; now Fortran does too, both exercised in
   `test_row132_axis_too_thin.py`.
3. Row 132(2) (mutation-check tautology): fixed -- `_ownership_violations`
   gained an `expected` parameter so the DROP shape is observable through
   the SAME function real callers use, not a hand-rolled set difference.
4. Full gate green: unit+regression, 8/8 e2e cells (fortran+python-jax x
   4 cases), 1/1 python-jax-mpi cell, `test.reference.results/` untouched.

## Correction (PR #41 audit, 2026-09-25)

The "No y-gate" reasoning above was wrong. y is constant ON the fault, but
what causes the hazard is which RANKS hold the node. When an npy boundary
lands on the fault plane (checkFaultMPIAlignment's DUPLICATE case, and
test_fault_mpi_boundary_arn's symmetric-y (1,2,1) case), both ranks build
every fault node, and with an x/z-only gate both wrote every `faultst*`
file.

Reproduced before the fix: `faultst000dp000.txt` written by ranks 0 and 1,
in both languages. The gate now also has `(nodeXyzIndex(2)==1 .and. mey/=0)`
in Fortran and `(iy == 0 and mey != 0)` in Python. After the fix, rank 0
alone writes it, and Fortran == Python.

`test_row127_station_ownership.py` now runs that y-seam case against the
serial run's file set. Because it passes `expected=`, a station written by
zero ranks fails too. The mutation check (drop the y term in both languages)
turns both "exactly one owner" checks red.

Scope of "FIXED" in item 1 below: x, y and z seams.
