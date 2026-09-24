# Row 114 -- station output port (faultst*/body*)

## Status: implementation drafted, NOT YET RUN against Fortran.

## Plan (locked)
- Serial backends (numpy/jax), ONE code path in driver.py's make_step_parts:
  - off-fault station history recorded in part_a, right after velDispUpdate
    (matches driver.f90:20-21 storeOffFaultStData's position).
  - on-fault station history recorded in part_b, right after FLT.faulting()
    (matches faulting.f90:24 storeOnFaultStationQuantSCEC's position --
    AFTER solveSWTW/solveRSF have updated fric for this step).
  - Both histories are new carry entries (on_st_hist, off_st_hist), sized
    (n_on, ncols, nsteps) / (n_off, 7, nsteps). Column formulas derived
    directly from library_output.f90's write statements -- see
    library_output.py's module docstring for the full derivation (in
    particular: friclaw>=3's creep-rate add-back at record time, since
    fric[SLIP_STRIKE]/[SLIPRATE_STRIKE] etc. are written BEFORE solveRSF
    adds the background slip rate).
- eqdyna3d.build_solver_state: reads bStations.txt unconditionally, runs
  meshgen.build_station_matching on the GLOBAL grid lines for the SERIAL
  path only (part is None). MPI path leaves st_on_idx/st_off_idx EMPTY.
- eqdyna3d.run_case: writes faultst*/body* via library_output's new
  write_onfault_stations/write_offfault_stations, guarded the same way
  frt_path is under backend.timing_only().
- eqdyna3d.run_case_mpi: prints a loud (rank 0, stderr) warning when
  bStations.txt names any station and writes NO station files -- chosen
  over a hard raise because test.tpv8 (the one case opted into
  python-jax-mpi) has stations in its DEFAULT bStations.txt and a raise
  would break that already-gated cell.

## nStressOutSign (board item 22a collision)
Column 8 (on-fault n-stress) uses `S['nStressOutSign']`, hardcoded -1.0
(pre-22a Fortran sign) in eqdyna3d.build_solver_state, with a TODO to read
it from readInputFiles.read_bglobal's g['nStressOutSign'] once that PR lands
on origin/master. Poll `git log --oneline origin/master -5` for a commit
whose subject contains "22a" before running the final parity gate.

## Verification -- TODO, in order
1. Rebuild bin/eqdyna in this worktree.
2. python3 testsys/run.py unit regression (must stay green).
3. New regression test (behavioural, rule 10a) asserting station filenames
   + column counts on a small case.
4. Per-case fortran vs python-jax sweep (testNameList minus drv.a6), write
   docs/perf_snapshots/station_spread_<date>_fortran_vs_jax.json.
5. tpv8/tpv29 python-numpy vs fortran.
6. tpv8 jax ms/step before/after (two step counts, differenced).
7. drv.a6 once, report only (chaotic case).

## Open questions / risks noted so far
- Header text (date line aside) is best-effort Fortran-format-descriptor
  reproduction; NOT yet round-tripped against a real Fortran header byte
  for byte-identity of the free-text lines (only DATA columns are the
  contractual parity target per the mission brief).
