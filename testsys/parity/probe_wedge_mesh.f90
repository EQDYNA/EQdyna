! NOT part of src/fortran/ (that tree is read-only, the reference). This is
! a NEW, separate probe program (mira-volkov parity evidence, item 3 of
! evidence_c_degen_port.py) that links against the SAME compiled .o object
! files eqdyna3d.f90 uses, calling the identical subroutine sequence
! eqdyna3d.f90's `program EQdyna` calls up through (and including)
! `call meshgen`, then dumps globalvar's elemTypeArr/mat/nodeElemIdRelation/
! nftnd to a text file for comparison against the Python port. No existing
! Fortran source file is modified.
!
! `probeAllocInit` below is a byte-for-byte copy of eqdyna3d.f90's
! `allocInit` subroutine body -- duplicated here (not linked from
! eqdyna3d.o) only because eqdyna3d.o's own `program EQdyna` entry point
! would collide with this file's `program ProbeWedgeMesh` at link time.
program ProbeWedgeMesh
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    integer (kind = 4) :: i, iMPIerr, unit1
    integer (kind = 4) :: c1, c2, c11, c12, c13, cother

    call MPI_Init(iMPIerr)
    call mpi_comm_rank(MPI_COMM_WORLD, me, iMPIerr)
    call mpi_comm_size(MPI_COMM_WORLD, totalNumOfMPIProcs, iMPIerr)

    call readglobal
    call readmodelgeometry

    allocate(fxmin(ntotft), fxmax(ntotft), fymin(ntotft), fymax(ntotft), fzmin(ntotft), fzmax(ntotft), material(nmat, n2mat))
    allocate(nonfs(ntotft))
    allocate(fltxyz(2, 4, ntotft))

    call readfaultgeometry
    call readmaterial
    call readstations1

    allocate(OffFaultStNodeIdIndex(2, totalNumOfOffSt), xonfs(2, maxval(nonfs), ntotft), x4nds(3, totalNumOfOffSt))

    call readstations2
    if (insertFaultType > 0) call read_fault_rough_geometry
    call checkInputConsistency

    allocate(nftnd(ntotft), localShapeFunc(nrowsh, nen))

    call calcLocalShapeFunc
    call countMeshEntities
    call probeAllocInit
    call meshgen

    c1 = 0; c2 = 0; c11 = 0; c12 = 0; c13 = 0; cother = 0
    do i = 1, totalNumOfElements
        if (elemTypeArr(i) == 1) then
            c1 = c1 + 1
        else if (elemTypeArr(i) == 2) then
            c2 = c2 + 1
        else if (elemTypeArr(i) == 11) then
            c11 = c11 + 1
        else if (elemTypeArr(i) == 12) then
            c12 = c12 + 1
        else if (elemTypeArr(i) == 13) then
            c13 = c13 + 1
        else
            cother = cother + 1
        endif
    enddo

    open(newunit=unit1, file='probe_wedge_mesh_out.txt', status='replace')
    write(unit1,*) 'totalNumOfElements ', totalNumOfElements
    write(unit1,*) 'totalNumOfNodes ', totalNumOfNodes
    write(unit1,*) 'count_type1 ', c1
    write(unit1,*) 'count_type2 ', c2
    write(unit1,*) 'count_type11 ', c11
    write(unit1,*) 'count_type12 ', c12
    write(unit1,*) 'count_type13 ', c13
    write(unit1,*) 'count_other ', cother
    write(unit1,*) 'nftnd1 ', nftnd(1)
    write(unit1,*) 'BEGIN_ELEMENTS'
    do i = 1, totalNumOfElements
        if (elemTypeArr(i) == 11 .or. elemTypeArr(i) == 12 .or. elemTypeArr(i) == 13) then
            write(unit1,'(I8,1X,I4,1X,8I10,1X,5ES24.15)') i, elemTypeArr(i), &
                nodeElemIdRelation(1:8,i), mat(i,1:5)
        endif
    enddo
    write(unit1,*) 'END_ELEMENTS'
    close(unit1)

    call MPI_Finalize(iMPIerr)
end program ProbeWedgeMesh

subroutine probeAllocInit
    use globalvar
    use errorCodes
    implicit none

    write(mm,'(i6)') me
    mm = trim(adjustl(mm))

    allocate(eqNumIndexArr(sizeOfEqNumIndexArr), eqNumStartIndexLoc(totalNumOfNodes), &
            numOfDofPerNodeArr(totalNumOfNodes), meshCoor(ndof,totalNumOfNodes), &
            fnms(totalNumOfNodes), surfaceNodeIdArr(totalNumOfNodes), &
            nodeElemIdRelation(nen,totalNumOfElements), mat(totalNumOfElements,5), &
            elemTypeArr(totalNumOfElements), eleporep(totalNumOfElements), &
            pstrain(totalNumOfElements), eledet(totalNumOfElements), &
            elemass(nee,totalNumOfElements), eleshp(nrowsh-1,nen,totalNumOfElements), &
            ss(6,totalNumOfElements), phi(nen,4,totalNumOfElements))

    meshCoor = 0.0d0
    fnms     = 0.0d0
    eqNumIndexArr = 0
    eqNumStartIndexLoc = 0
    numOfDofPerNodeArr = 0
    surfaceNodeIdArr   = 0
    nodeElemIdRelation = 0
    elemTypeArr = 0
    mat      = 0.0d0
    eleporep = 0.0d0
    pstrain  = 0.0d0
    eledet   = 0.0d0
    elemass  = 0.0d0
    eleshp   = 0.0d0
    ss       = 0.0d0
    phi      = 0.0d0

    nftmx = maxval(nftnd)
    if(nftmx<=0) nftmx=1
    nonmx = sum(nonfs)

    allocate(onFaultTPHist(2,nftmx,nstep,ntotft))
    onFaultTPHist = 0.0d0

    allocate(nsmp(2,nftmx,ntotft), fnft(nftmx,ntotft), un(3,nftmx,ntotft), &
                us(3,nftmx,ntotft), ud(3,nftmx,ntotft), fric(100,nftmx,ntotft), &
                arn(nftmx,ntotft),  anonfs(3,nonmx), fltgm(nftmx), &
                Tatnode(nftmx,ntotft), patnode(nftmx,ntotft))
    fltgm   = 0
    nsmp    = 0
    fnft    = 99999.0d0
    fric    = 0.0d0
    un      = 0.0d0
    us      = 1000.0d0
    ud      = 0.0d0
    arn     = 0.0d0
    anonfs  = 0
    Tatnode = 0.0d0
    patnode = 0.0d0

    allocate(stressCompIndexArr(totalNumOfElements))
    allocate(stressArr(5*sizeOfEqNumIndexArr))
    stressCompIndexArr = 0
    stressArr = 0.0d0
end subroutine probeAllocInit
