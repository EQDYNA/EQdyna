! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
!
!=======================================================================
! errorCodes -- the single registry of EQdyna's fatal exit codes, and the
! single way to exit on one.
!
! WHY THIS MODULE EXISTS
!
! Before it, refusing a run was unreliable in three separate ways:
!
!   1. `stop` with no code, and `stop 'some message'`, both exit with
!      status ZERO under gfortran (verified directly: a program whose only
!      statement is `stop 'refused'` prints "STOP refused" and exits 0).
!      Thirteen sites used one of those two forms, including the two
!      alignment gates whose entire purpose is to refuse a mesh that would
!      produce a wrong answer. A batch script, testsys/run.py, or CI
!      checking $? scored those refusals as passes.
!
!   2. Codes above 255 wrap. A process exit status is 8 bits, so the old
!      `stop 1001` was reported by the shell as 233 and `stop 2002` as
!      210 -- the numbers in the source and the numbers an operator saw
!      had no visible relationship. Every code here is <= 125, so the
!      code in the source is the code the shell reports.
!
!   3. A rank-local `stop` in an MPI job terminates one rank and leaves
!      the other N-1 blocked in their next collective, so the failure
!      presented as a HANG rather than an error -- the worst failure mode
!      for a long queue job. abortRun uses MPI_Abort, which tears down the
!      whole communicator.
!
! CODE RANGES  (1-125; 126, 127 and 128+n are reserved by the shell)
!
!    1- 9   reserved (1 = generic/unclassified)
!   11-19   configuration and parameter consistency
!   21-29   input files
!   31-39   fault geometry
!   41-49   mesh generation and element quality
!   51-59   MPI and domain decomposition
!   61-69   numerics and runtime state
!   71-79   external libraries
!
! Decade 0 is left free in every range so a related code can be added
! next to its neighbours instead of at the end of the list.
!
! ADDING A CODE
!   1. Add the parameter below, in its range, with a one-line comment.
!      The comment IS the documentation -- README.md's "Exit codes" table is
!      generated from these lines, and testsys/regression/test_stop_exit_status.py
!      fails if the two disagree, so the table cannot drift.
!   2. Call abortRun(ERR_..., 'what was wrong and what to do about it').
! testsys/regression/test_stop_exit_status.py fails the build if a bare
! `stop` or a `stop 'string'` reappears in src/.
!=======================================================================
MODULE errorCodes

    implicit none

    ! --- generic ------------------------------------------------------
    integer, parameter :: ERR_GENERIC                = 1   ! unclassified fatal error

    ! --- 11-19 configuration and parameter consistency ----------------
    integer, parameter :: ERR_CFG_Q_NEEDS_ELASTIC    = 11  ! C_Q=1 requires C_elastic=1
    integer, parameter :: ERR_CFG_Q_NEEDS_UNIFORM    = 12  ! C_Q=1 requires rat=1.0 (uniform elements)
    integer, parameter :: ERR_CFG_PLASTIC_OUTPUT     = 13  ! output_plastic=1 requires C_elastic=0
    integer, parameter :: ERR_CFG_PROFILE_ENV_INVALID = 14 ! EQDYNA_PROFILE is set to a value other than unset, "", "1" or "0"
    integer, parameter :: ERR_CFG_NSTRESS_SIGN_INVALID = 15 ! bGlobal.txt's station n-stress sign is neither +1 nor -1

    ! --- 21-29 input files --------------------------------------------
    integer, parameter :: ERR_INPUT_FILE_MISSING     = 21  ! a required FE_*.txt / data file is absent
    integer, parameter :: ERR_INPUT_FILE_STALE       = 22  ! an input file is present but predates this binary's input format

    ! --- 31-39 fault geometry -----------------------------------------
    integer, parameter :: ERR_GEOM_ROUGH_INVALID     = 31  ! bFault_Rough_Geometry.txt does not match this mesh

    ! --- 41-49 mesh generation and element quality --------------------
    integer, parameter :: ERR_MESH_STRESS_ARR_SMALL  = 41  ! sizeOfStressDofIndexArr exceeds 5*sizeOfEqNumIndexArr
    integer, parameter :: ERR_MESH_COUNT_MISMATCH    = 42  ! meshgen's node/element/equation tallies disagree
    integer, parameter :: ERR_MESH_EQNUM_MISMATCH    = 43  ! eqNumIndexArrLocTag /= sizeOfEqNumIndexArr
    integer, parameter :: ERR_MESH_FAULT_MISMATCH    = 44  ! nftnd0 /= nftnd (meshgen vs countMeshEntities)
    integer, parameter :: ERR_MESH_MULTIFAULT_MSNODE = 45  ! master-node construction cannot handle ntotft>1
    integer, parameter :: ERR_MESH_BAD_WEDGE         = 46  ! degenerate wedge built with unequal node ids
    integer, parameter :: ERR_MESH_BAD_JACOBIAN      = 47  ! non-positive Jacobian determinant (inverted element)
    integer, parameter :: ERR_MESH_MATERIAL_UNSET    = 48  ! an element has no material property assigned
    integer, parameter :: ERR_MESH_GRID_TOO_LARGE    = 49  ! a 1D grid line (or its per-rank slice) exceeds the fixed-size 10000 buffer

    ! --- 51-59 MPI and domain decomposition ---------------------------
    integer, parameter :: ERR_MPI_FAULT_ALIGNMENT    = 51  ! a rank boundary in y coincides with the fault plane (not raised as of v5.8.2; downgraded to a NOTICE, see syncArnBoundary)
    integer, parameter :: ERR_MPI_BAD_NEIGHBOR       = 52  ! point-to-point exchange with a rank outside 0..npx*npy*npz-1

    ! --- 61-69 numerics and runtime state -----------------------------
    integer, parameter :: ERR_NUM_PML_ALIGNMENT      = 61  ! element centre lies exactly on a PML bound
    integer, parameter :: ERR_NUM_PML_DAMPING        = 62  ! negative PML damping vector component
    integer, parameter :: ERR_NUM_VELOCITY_NAN       = 63  ! NaN velocity during time stepping
    integer, parameter :: ERR_NUM_NEGATIVE_DEPTH     = 64  ! negative depth passed to a B-function

    ! --- 71-79 external libraries -------------------------------------
    integer, parameter :: ERR_NETCDF                 = 71  ! a NetCDF call returned an error

CONTAINS

!-----------------------------------------------------------------------
! abortRun -- print a structured fatal-error block and terminate the whole
! job with `code` as the process exit status.
!
! Every rank that reaches this prints, because a rank-local failure is
! often the interesting one and suppressing non-master output would hide
! it. MPI_Abort tears down the communicator so the surviving ranks fail
! fast instead of blocking forever in their next collective.
!
! `code` is passed to MPI_Abort as the error code. How faithfully that number
! reaches the caller depends on the launcher: mpirun (Open MPI, hydra) reports
! it directly, but srun reports the MAXIMUM status across tasks, so a straggler
! killed during MPI_Abort teardown yields 137/143 and masks the code, and
! sbatch reports the wrapper script's status unless it ends with `exit $?`.
! So: a NON-ZERO status is the reliable signal and the printed FATAL block is
! authoritative; the specific number is advisory under srun/ibrun (ls6, grace).
! If MPI is not initialised (or already
! finalised) the routine falls back to `call exit(code)`, which both
! gfortran and ifort provide -- a plain `stop code` cannot be used here
! because the standard requires a constant stop-code, not a variable.
!-----------------------------------------------------------------------
    subroutine abortRun(code, reason)

        implicit none
        include 'mpif.h'

        integer, intent(in)          :: code
        character(len=*), intent(in) :: reason

        integer :: rank, ierr
        logical :: mpiStarted, mpiDone, mpiReady

        rank       = -1
        mpiStarted = .false.
        mpiDone    = .false.
        ierr       = 0
        ! MPI_Initialized stays .true. AFTER MPI_Finalize, so it alone is not a
        ! licence to call MPI: MPI_Abort on a finalised communicator is
        ! undefined behaviour. Nothing between MPI_Finalize and the normal exit
        ! calls abortRun today, but an error check added to an output routine
        ! would land exactly there.
        call MPI_Initialized(mpiStarted, ierr)
        if (ierr /= 0) mpiStarted = .false.
        if (mpiStarted) then
            call MPI_Finalized(mpiDone, ierr)
            if (ierr /= 0) mpiDone = .true.   ! cannot confirm; do not touch MPI
        endif
        mpiReady = mpiStarted .and. (.not. mpiDone)
        if (mpiReady) then
            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            if (ierr /= 0) rank = -1
        endif

        write(*,*)
        write(*,*) '==================== EQdyna: FATAL ===================='
        if (rank >= 0) then
            write(*,*) ' rank      : ', rank
        endif
        write(*,*) ' exit code : ', code
        write(*,*) ' reason    : ', trim(reason)
        write(*,*) ' See the "Exit codes" table in README.md for this code.'
        write(*,*) '======================================================='
        write(*,*)
        flush(6)

        if (mpiReady) then
            call MPI_Abort(MPI_COMM_WORLD, code, ierr)
        endif
        ! Reached only when MPI is unavailable, or if MPI_Abort returns.
        call exit(code)

    end subroutine abortRun

subroutine requireValidNeighbor(neighbor, site, ixyz, ib, extra)
! Refuse a point-to-point exchange with a rank that does not exist, and say
! EVERYTHING needed to place it: the call site, the direction, the side, this
! rank's 3-D coordinates, the decomposition, and the offending neighbour id.
!
! Without this, an out-of-range neighbour surfaces only as the MPI library's
! own message --
!     *** An error occurred in MPI_Sendrecv
!     *** MPI_ERR_RANK: invalid rank
! -- which names neither the call site nor the value, and leaves the reader to
! guess among every exchange in the code. That cost real time on test.tpv36
! (2026-09-16): two decompositions tried, a serial control run, and a read of
! both sendrecv sites, to learn something the MPI library already knew.
    use globalvar, only : me, npx, npy, npz
    implicit none
    integer (kind = 4), intent(in) :: neighbor, ixyz, ib
    character(len=*), intent(in) :: site, extra
    integer (kind = 4) :: mex, mey, mez, nranks
    character(len=600) :: reason

    nranks = npx*npy*npz
    if (neighbor >= 0 .and. neighbor < nranks) return

    mex = int(me/(npy*npz))
    mey = int((me - mex*npy*npz)/npz)
    mez = int(me - mex*npy*npz - mey*npz)

    write(reason, '(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,A)') &
        'invalid MPI neighbour in ', trim(site), ': neighbour=', neighbor, &
        ' outside 0..', nranks-1, '.  me=', me, ' (mex,mey,mez)=(', mex,   &
        ',', mey, ',', mez, ')  npx,npy,npz=', npx, ',', npy, ',', npz,    &
        '  ixyz=', ixyz, ' ib=', ib, '.  ', trim(extra)
    call abortRun(ERR_MPI_BAD_NEIGHBOR, trim(reason))
end subroutine requireValidNeighbor

END MODULE errorCodes
