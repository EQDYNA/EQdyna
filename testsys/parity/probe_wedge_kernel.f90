! NOT part of src/fortran/ (that tree is read-only, the reference). This is a
! NEW, standalone probe program (mira-volkov parity evidence for the C_degen
! wedge-degeneration DYNAMICS port) that links against the SAME compiled .o
! object files eqdyna3d.f90 uses and calls the REAL, UNMODIFIED, production
! subroutines directly:
!
!   calcLocalShapeFunc  (calcLocalShapeFunc.f90)   -- sets w/localShapeFunc
!   calcGlobalShapeFunc (calcGlobalShapeFunc.f90)  -- det/eleshp/xs, incl. the
!                                                      elemTypeArr==11/12 merge
!   contm               (assembleGlobalMass.f90)   -- per-node lumped mass
!   calcSSPhi4Hrgls      (assembleGlobalMass.f90)   -- ss/phi hourglass tensors
!                                                      (calls vlm, func_lib.f90)
!   calcElemMass        (calcElemMass.f90)         -- gravity/body-force term
!   calcElemKU          (calcElemKU.f90)           -- interior stiffness/force
!   calcHourglassResist (calcHourglassResist.f90)  -- KF78 hourglass force
!
! No mesh, no case, no MPI decomposition: this program builds ONE synthetic
! element's module state directly (elemTypeArr, mat, nodeElemIdRelation,
! numOfDofPerNodeArr, eqNumStartIndexLoc, eqNumIndexArr as an IDENTITY map so
! nodalForceArr(3*(i-1)+j) lands exactly on local node i / dof j -- see
! evidence_wedge_kernel.py for why that makes the comparison direct), reads
! one or more (elemType, xl, mate, vl, dl, priorStress, dt, rdampk) cases from
! stdin (list-directed, matching test_drucker_prager_kernel.py's proven
! stdin protocol -- SAME INPUT BYTES on both sides, not a re-typed literal),
! and writes every intermediate/final array list-directed to stdout.
!
! C_elastic is fixed to 1 (elastic) for every case: the Drucker-Prager
! plastic branch is already isolated and verified against Fortran by
! testsys/regression/test_drucker_prager_kernel.py, and eqdyna3d.py refuses
! C_degen>3 (the only source of wedge elements) combined with C_elastic==0 --
! so a wedge element in this port's actual runtime scope NEVER takes that
! branch, and re-testing it here would not exercise anything this port can
! reach. C_Q is fixed to 0 (the non-attenuating stress-update branch;
! C_Q==1 is unrelated to C_degen and untouched by this port).
program probe_wedge_kernel
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    integer (kind = 4) :: iMPIerr, ncases, ic, i, j, etype
    real (kind = dp) :: xl(3,8), mate(5), vl(24), dl(24), stress(12)
    real (kind = dp) :: al(3,8), elresf(24), elmass(24)
    real (kind = dp) :: det, xs(3,3), globalShapeFunc(4,8)
    real (kind = dp) :: constk, porep, pstrmag
    real (kind = dp) :: s_dt, s_rdampk
    logical :: lcubic

    call MPI_Init(iMPIerr)
    call mpi_comm_rank(MPI_COMM_WORLD, me, iMPIerr)
    call mpi_comm_size(MPI_COMM_WORLD, totalNumOfMPIProcs, iMPIerr)

    allocate(localShapeFunc(nrowsh, nen))
    call calcLocalShapeFunc   ! sets module w / localShapeFunc (real subroutine call)

    C_elastic = 1
    C_Q = 0
    totalNumOfElements = 1
    allocate(elemTypeArr(1), mat(1,5), ss(6,1), phi(nen,4,1))
    allocate(nodeElemIdRelation(nen,1), numOfDofPerNodeArr(nen))
    allocate(eqNumStartIndexLoc(nen), eqNumIndexArr(nee))
    allocate(nodalForceArr(nee), velArr(3,nen), dispArr(3,nen))

    do i = 1, nen
        nodeElemIdRelation(i,1) = i
        numOfDofPerNodeArr(i) = 3
        eqNumStartIndexLoc(i) = 3*(i-1)
    enddo
    do i = 1, nee
        eqNumIndexArr(i) = i   ! identity map: (node i, local dof j) -> global eq 3*(i-1)+j
    enddo

    read(*,*) ncases
    do ic = 1, ncases
        read(*,*) etype
        read(*,*) xl
        read(*,*) mate
        read(*,*) vl
        read(*,*) dl
        read(*,*) stress(1:6)
        read(*,*) s_dt, s_rdampk

        dt = s_dt
        rdampk = s_rdampk
        elemTypeArr(1) = etype
        mat(1,1:5) = mate

        lcubic = .true.   ! calcGlobalShapeFunc's own dummy arg is unused dead code (grepped)
        call calcGlobalShapeFunc(xl, det, globalShapeFunc, 1, xs, lcubic)
        call contm(globalShapeFunc, det, elmass, mat(1,3))
        call calcSSPhi4Hrgls(1, xl, xs, globalShapeFunc)

        ! rdampm==0 (always, in this codebase) and C_elastic==1 make `al`
        ! exactly the assembleGlobalKU.f90:19-20 value for a wedge run: zero.
        al = 0.0d0
        elresf = 0.0d0
        call calcElemMass(elmass, al, elresf)

        stress(7:12) = 0.0d0   ! only C_Q==1's memory-variable slots; unused here
        constk = -det           ! assembleGlobalKU.f90:28 passes -eledet(nel)
        porep = 0.0d0
        pstrmag = 0.0d0
        call calcElemKU(globalShapeFunc(1:3,1:8), mate, vl, dl, stress, elresf, &
                        constk, porep, pstrmag, xl)

        do i = 1, nen
            do j = 1, 3
                velArr(j,i) = vl(3*(i-1)+j)
                dispArr(j,i) = dl(3*(i-1)+j)
            enddo
        enddo
        nodalForceArr = 0.0d0
        call calcHourglassResist

        write(*,'(A)') 'BEGIN_CASE'
        write(*,*) det, globalShapeFunc, xs, elmass, ss(1:6,1), phi(1:nen,1:4,1), &
                    elresf, nodalForceArr
        write(*,'(A)') 'END_CASE'
    enddo

    call MPI_Finalize(iMPIerr)
end program probe_wedge_kernel
