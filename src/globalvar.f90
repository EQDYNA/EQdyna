! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
MODULE globalvar

    implicit none
    integer, parameter :: dp = selected_real_kind(15,307) ! precision of double precision

    !=====================================================================
    ! fric() first-index slot map. fric(:,i,ift) holds all per-fault-node-
    ! pair scalar quantities for node pair i on fault ift, in one flat
    ! record; the slot constants below name every index that is read or
    ! written by name (i.e. with a literal first index) somewhere in src/.
    ! Evidence for each slot: src/netcdf_io.f90 (read table + inline
    ! comments), src/library_output.f90:output_frt (written columns),
    ! src/faulting.f90 (solveSWTW/solveRSF/NewtonRaphson usage),
    ! scripts/defaultParameters.py:162-216 (on_fault_vars, Python mirror
    ! of the same slot numbering).
    ! Values are unchanged from the literals used throughout src/; this
    ! block is compile-time aliasing only.
    !=====================================================================
    integer, parameter ::                        &
        FRIC_SLOT_SW_FS             =  1, &  ! slip-weakening static friction coefficient
        FRIC_SLOT_SW_FD             =  2, &  ! slip-weakening dynamic friction coefficient
        FRIC_SLOT_SW_D0             =  3, &  ! slip-weakening critical slip distance D0, m
        FRIC_SLOT_COHESION          =  4, &  ! fault cohesion, Pa
        FRIC_SLOT_TW_T0             =  5, &  ! time-weakening rupture-time offset t0, s
        FRIC_SLOT_NORM_STRESS_ADD   =  6, &  ! additive normal-stress term (non-restart path); no writer in current src/, always 0
        FRIC_SLOT_INIT_NORM         =  7, &  ! initial effective normal stress, Pa
        FRIC_SLOT_INIT_STRIKE_SHEAR =  8, &  ! initial strike-direction shear stress, Pa
        FRIC_SLOT_RSF_A             =  9, &  ! RSF direct-effect parameter a
        FRIC_SLOT_RSF_B             = 10, &  ! RSF evolution-effect parameter b
        FRIC_SLOT_RSF_DC            = 11, &  ! RSF characteristic slip distance Dc, m
        FRIC_SLOT_RSF_V0            = 12, &  ! RSF reference slip rate, m/s
        FRIC_SLOT_RSF_R0            = 13, &  ! RSF reference friction coefficient
        FRIC_SLOT_RSF_FW            = 14, &  ! RSF flash-heating fully-weakened friction coefficient
        FRIC_SLOT_RSF_VW            = 15, &  ! RSF flash-heating weakening velocity, m/s
        FRIC_SLOT_TP_A_HY           = 16, &  ! thermal pressurization hydraulic diffusivity
        FRIC_SLOT_TP_A_TH           = 17, &  ! thermal pressurization thermal diffusivity
        FRIC_SLOT_TP_ROUC           = 18, &  ! thermal pressurization rho*c, heat capacity per volume
        FRIC_SLOT_TP_LAMBDA         = 19, &  ! thermal pressurization Lambda, pore-pressure/temperature coupling
        FRIC_SLOT_STATE             = 20, &  ! RSF state variable theta
        FRIC_SLOT_THETA_PC          = 23, &  ! normal-stress-evolution state variable theta_pc (Shi & Day, 2013)
        FRIC_SLOT_THETA_PC_DOT      = 24, &  ! time derivative of FRIC_SLOT_THETA_PC
        FRIC_SLOT_VINI_N            = 25, &  ! background/creep slip-rate offset, normal component
        FRIC_SLOT_VINI_X            = 26, &  ! background/creep slip-rate offset, strike (x) component
        FRIC_SLOT_VINI_Z            = 27, &  ! background/creep slip-rate offset, dip (z) component
        FRIC_SLOT_VEL_MASTER_X      = 31, &  ! master-node velocity, x component, m/s
        FRIC_SLOT_VEL_MASTER_Y      = 32, &  ! master-node velocity, y component, m/s
        FRIC_SLOT_VEL_MASTER_Z      = 33, &  ! master-node velocity, z component, m/s
        FRIC_SLOT_VEL_SLAVE_X       = 34, &  ! slave-node velocity, x component, m/s
        FRIC_SLOT_VEL_SLAVE_Y       = 35, &  ! slave-node velocity, y component, m/s
        FRIC_SLOT_VEL_SLAVE_Z       = 36, &  ! slave-node velocity, z component, m/s
        FRIC_SLOT_TP_H              = 40, &  ! thermal-pressurization shear-zone half-width, m
        FRIC_SLOT_TP_TINI           = 41, &  ! initial fault temperature, K
        FRIC_SLOT_TP_PINI           = 42, &  ! initial fault pore pressure, Pa
        FRIC_SLOT_CREEP_VMIN        = 46, &  ! creeping/initial slip-rate lower bound, m/s
        FRIC_SLOT_PEAK_SLIPRATE     = 47, &  ! current/peak (trial) slip-rate magnitude, m/s
        FRIC_SLOT_SHEAR_MAG         = 48, &  ! in-plane (strike+dip) shear-traction magnitude, Pa
        FRIC_SLOT_INIT_DIP_SHEAR    = 49, &  ! initial dip-direction shear stress, Pa
        FRIC_SLOT_TP_NORM_TP        = 51, &  ! TP-derived pore-pressure contribution to normal traction, Pa
        FRIC_SLOT_TP_TEMP           = 52, &  ! current fault temperature, K
        FRIC_SLOT_SLIP_STRIKE       = 71, &  ! cumulative slip, strike component, m
        FRIC_SLOT_SLIP_DIP          = 72, &  ! cumulative slip, dip component, m
        FRIC_SLOT_SLIP_NORM         = 73, &  ! cumulative slip, normal component, m
        FRIC_SLOT_SLIPRATE_STRIKE   = 74, &  ! slip rate, strike component, m/s
        FRIC_SLOT_SLIPRATE_DIP      = 75, &  ! slip rate, dip component, m/s
        FRIC_SLOT_SLIPRATE_MAX      = 76, &  ! running maximum slip-rate magnitude, m/s
        FRIC_SLOT_CUM_SLIP          = 77, &  ! cumulative slip magnitude, m
        FRIC_SLOT_TRACT_NORM        = 78, &  ! current effective normal traction, Pa
        FRIC_SLOT_TRACT_STRIKE      = 79, &  ! current strike-direction shear traction, Pa
        FRIC_SLOT_TRACT_DIP         = 80, &  ! current dip-direction shear traction, Pa
        FRIC_SLOT_NUC_DTAU0         = 81     ! nucleation stress-perturbation amplitude, Pa

    !=====================================================================
    ! Character / string identifiers
    !=====================================================================
    character(len=30) :: mm                    ! MPI-rank suffix appended to per-process output file names
    character(len=30) :: sttmp, dptmp, bodytmp ! scratch string buffers
    character(len=30) :: projectname = 'San-Ti'
    character(len=30) :: author      = 'Sophon'
    character(len=90) :: stLocStamp            ! on-fault-station output path/label

    !=====================================================================
    ! Model control flags and case identifiers
    !=====================================================================
    logical           :: fltMPI(6)              ! per-direction flags: does this MPI rank touch a fault boundary
    integer (kind = 4) :: mode                   ! run mode (fresh start vs. restart from prior cycle)
    integer (kind = 4) :: TPV = -1               ! SCEC TPV benchmark ID, or case tag
    integer (kind = 4) :: friclaw                ! friction law selector: 1/2 slip/time-weakening, 3/4/5 RSF variants
    integer (kind = 4) :: insertFaultType        ! >0: non-planar/rough fault geometry inserted
    integer (kind = 4) :: nucfault               ! fault index on which artificial nucleation is applied
    integer (kind = 4) :: ntotft                 ! total number of faults
    integer (kind = 4) :: C_elastic              ! 1 = elastic version; 0 = plastic version
    integer (kind = 4) :: C_nuclea                ! 1 = allow artificial nucleation; 0 = disabled
    integer (kind = 4) :: C_Q  = 0                ! only with C_elastic==1: 1 = allow Q attenuation; 0 = do not
    integer (kind = 4) :: C_hg = 1                ! hourglass control: 1 = KF78, 2 = viscous HG
    integer (kind = 4) :: C_dc = 0                ! double-couple source: 1 = yes, 0 = no
    integer (kind = 4) :: C_degen                  ! degenerate-element flag: 0 = brick, 1 = wedge, 2 = tetra
    integer (kind = 4) :: output_plastic          ! 1 = write plastic-strain output
    integer (kind = 4) :: writeCompTime = 0       ! 1 = write per-stage wall-clock timing
    integer (kind = 4) :: outputGroundMotion       ! 1 = write ground-motion output
    integer (kind = 4) :: outputFinalSurfDisp = 0 ! 1 = write final surface displacement

    !=====================================================================
    ! Mesh & geometry
    !=====================================================================
    integer (kind = 4) :: nsd=3, ndof=3, ned=3, nesd=3            ! spatial/DOF dimensionality
    integer (kind = 4) :: nen=8, nee=24, nrowsh=4, nrowb=6, nrowc=6, nstr=6 ! element shape/stress bookkeeping (8-node brick)
    integer (kind = 4) :: totalNumOfNodes, totalNumOfElements, totalNumOfEquations
    integer (kind = 4) :: sizeOfEqNumIndexArr, sizeOfStressDofIndexArr
    integer (kind = 4) :: npx, npy, npz                            ! processor grid dimensions
    integer (kind = 4) :: nnx, nnz                                 ! fault-plane node counts along strike/dip
    integer (kind = 4) :: nmat, n2mat                              ! number of material blocks / material properties per block
    integer (kind = 4) :: nftmx, nonmx                             ! max fault-node-pair count / max on-fault station count (array sizing)
    integer (kind = 4) :: nt                                       ! current time step index
    integer (kind = 4) :: nstep                                    ! total number of time steps
    integer (kind = 4) :: dis4uniF, dis4uniB                        ! distance (in cells) to uniform-mesh region, front/back
    integer (kind = 4) :: surface_nnode = 0                        ! number of free-surface output nodes
    integer (kind = 4) :: fltnum(6) = 0                             ! fault-node counts per MPI-boundary direction
    real (kind = dp) :: dx, dy, dz                ! grid cell sizes, m
    real (kind = dp) :: rat                       ! geometrical enlarging ratio of cell size outside the uniform-mesh region
    real (kind = dp) :: xmin,  xmax,  ymin,  ymax,  zmin,  zmax   ! domain boundaries, m
    real (kind = dp) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1  ! inner (uniform-mesh) domain boundaries, m
    real (kind = dp) :: dxtmp                     ! scratch cell-size value
    real (kind = dp) :: gamar, roumax             ! rough-fault geometry generation parameters
    real (kind = dp) :: rough_fx_min, rough_fx_max, rough_fz_min ! left/right/top boundaries of the inserted rough-fault interface, m
    real (kind = dp) :: str1ToFaultAngle          ! angle between max compressive stress and fault strike, degrees
    real (kind = dp) :: devStrToStrVertRatio      ! deviatoric-to-vertical stress ratio

    !=====================================================================
    ! MPI / parallel decomposition
    !=====================================================================
    integer (kind = 4) :: me                       ! this MPI rank's ID
    integer (kind = 4) :: totalNumOfMPIProcs       ! total number of MPI ranks
    integer (kind = 4) :: masterProcsId = 0        ! rank ID of the master process

    !=====================================================================
    ! Fault & friction (module-level scalars mirrored per-node in fric())
    !=====================================================================
    real (kind = dp) :: slipRateThres              ! slip-rate threshold used to flag rupture-time arrival, m/s
    real (kind = dp) :: nucR, nucRuptVel, nucdtau0, nucT ! nucleation patch radius (m), forced rupture velocity (m/s), stress-perturbation amplitude (Pa), duration (s)
    real (kind = dp) :: xsource, ysource, zsource  ! hypocenter/nucleation-patch coordinates, m
    real (kind = dp) :: fdip, fstrike              ! fault dip / strike angles, degrees
    real (kind = dp) :: max_norm = -40.0d6         ! upper normal-stress cap enforced on non-planar/elastic faults, Pa
    real (kind = dp) :: min_norm = -10.0d6         ! lower normal-stress cap enforced on non-planar/elastic faults, Pa
    ! Scalar defaults mirrored per-node into fric(); see fric() slot map above.
    ! Currently unused as scalars in src/ (kept: not part of this refactor's deletion scope).
    real (kind = dp) :: fric_sw_fs, fric_sw_fd, fric_sw_D0
    real (kind = dp) :: fric_rsf_a, fric_rsf_deltaa0, fric_rsf_b, fric_rsf_Dc
    real (kind = dp) :: fric_rsf_r0, fric_rsf_v0, fric_rsf_vinix, fric_rsf_viniz
    real (kind = dp) :: fric_rsf_fw, fric_rsf_vw, fric_rsf_deltavw0
    real (kind = dp) :: fric_tp_a_th, fric_tp_pini, fric_tp_Tini
    real (kind = dp) :: fric_ww, fric_w, fric_ini_sliprate

    !=====================================================================
    ! Plasticity / off-fault material response
    !=====================================================================
    real (kind = dp) :: critd0, cohes, brangle, bulk, coheplas, tv, ccosphi, sinphi

    !=====================================================================
    ! PML & damping
    !=====================================================================
    integer (kind = 4) :: nPML = 6                 ! thickness of the PML layer, counted in nodes
    real (kind = dp) :: R = 0.01d0                  ! theoretical PML reflection coefficient.
        ! Other options for pairs of (nPML/R) are:
        !   nPML/R=6/0.01;10/0.001;20/0.0001. Collino & Tsogka (2001)
    real (kind = dp) :: PMLb(8), vmaxPML
    real (kind = dp) :: rdampm=0.0d0, rdampk        ! Rayleigh mass-/stiffness-proportional damping coefficients
    real (kind = dp) :: kapa_hg = 0.1d0             ! viscous-hourglass coefficient; typical range 0.05-0.15 (Goudreau & Hallquist, 1982)

    !=====================================================================
    ! Materials
    !=====================================================================
    real (kind = dp) :: grav = 9.8d0                ! gravitational acceleration, m/s^2
    real (kind = dp) :: rhow                        ! fluid (pore water) density, kg/m^3
    real (kind = dp) :: w                           ! scratch/general-purpose scalar
    real (kind = dp) :: critt0=0.2d0                ! critical time for nucleation to occur, s
    real (kind = dp) :: srcrad0=2500.d0             ! nucleation source radius, m
    real (kind = dp) :: vrupt0=1500.d0              ! forced rupture velocity for artificial nucleation, m/s

    !=====================================================================
    ! Time stepping / bookkeeping
    !=====================================================================
    real (kind = dp) :: dt                          ! time-step size, s
    real (kind = dp) :: totalSimuTime               ! total simulated physical time, s
    real (kind = dp) :: timeElapsed = 0.0d0         ! elapsed simulated time, s
    real (kind = dp) :: tol = 1.0d-5                ! floating-point comparison tolerance

    !=====================================================================
    ! Outputs
    !=====================================================================
    integer (kind = 4) :: totalNumOfOffSt, numOfOnFaultStCount, numOfOffFaultStCount ! on-/off-fault station counts
    integer (kind = 4) :: numcount(9)               ! per-category output record counters

    !=====================================================================
    ! Timers / performance
    !=====================================================================
    integer (kind = 8) :: dateTimeStamp(8)          ! wall-clock date/time components (DATE_AND_TIME output)
    real (kind = dp) :: compTimeInSeconds(9) = 0.0d0 ! per-stage wall-clock cost, s
    real (kind = dp) :: startTimeStamp               ! MPI_WTIME() value at the start of the current stage
    real (kind = dp) :: simuStartTime                ! MPI_WTIME() value at simulation start
    real (kind = dp) :: MPICommTimeInSeconds = 0.0d0 ! cumulative MPI-communication wall-clock time, s
    real (kind = dp) :: totmemcost, memcost = 0.0d0  ! memory-cost bookkeeping
    real (kind = dp) :: pi = 4*atan(1.0_dp)

    !=====================================================================
    ! Allocatable arrays
    !=====================================================================
    real (kind = dp), allocatable, dimension(:) :: nodalForceArr, nodalMassArr, &
        v1,      stressArr,  miuonf, vponf,  eleporep,       &
        pstrain, eledet,     fnms,   fxmin,  fxmax,  fymin,  fymax,  &
        fzmin,   fzmax
    real (kind = dp), allocatable, dimension(:,:) :: meshCoor, &
        dispArr,      velArr,      mat,    localShapeFunc,    fnft,   arn, &
        slp4fri, elemass,ss,     plane1, plane2, &
        Tatnode,patnode,OffFaultStGramSCEC,   material,       rough_geo,      &
        x4nds
    real (kind = dp), allocatable, dimension(:,:,:) :: fric,    &
        un,     us,     ud,     onFaultQuantHistSCECForm, eleshp, phi,    &
        fltxyz, xonfs
    real (kind = dp), allocatable, dimension(:,:,:,:):: onFaultTPHist

    integer (kind = 4), allocatable, dimension(:) :: nftnd,     &
        eqNumIndexArr,    stressCompIndexArr,    elemTypeArr,     eqNumStartIndexLoc,  numOfDofPerNodeArr,   surfaceNodeIdArr,&
        nonfs,  n4yn,   fltl,   fltr,   fltf,   fltb,   fltd,   &
        fltu,   fltgm
    integer (kind = 4), allocatable, dimension(:,:) :: nodeElemIdRelation,     &
        anonfs, idhist, OffFaultStNodeIdIndex
    integer (kind = 4), allocatable, dimension(:,:,:) :: nsmp

    integer (kind = 4) :: np = 1000000              ! legacy default array-sizing hint (see individual allocate() calls for actual sizes)

end MODULE globalvar
