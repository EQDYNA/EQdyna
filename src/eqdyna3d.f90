! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
program EQdyna
    use globalvar
    implicit none
    include 'mpif.h'
        
    integer (kind = 4) :: i, iMPIerr

    call MPI_Init(iMPIerr)
    call mpi_comm_rank(MPI_COMM_WORLD,me,iMPIerr)
    call mpi_comm_size(MPI_COMM_WORLD,totalNumOfMPIProcs,iMPIerr)

    if (me == masterProcsId) then 
        write(*,*) '====================================================================='
        write(*,*) '==================   Welcome to EQdyna 5.3.3  ======================='
        write(*,*) '===== Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>           ====' 
        write(*,*) '====    & Dunyu Liu <dliu@ig.utexas.edu> under MIT License.      ===='
        write(*,*) '============== https://github.com/EQDYNA/EQdyna.git   ==============='
        write(*,*) '=                                                                   ='
        write(*,*) '=   EQdyna is a parallel finite element software to simulate        ='
        write(*,*) '=    earthquake spontaneous dynamic rupture, seismic wave           ='
        write(*,*) '=    propagation and high frequency deterministic ground motions.   ='
        write(*,*) '=                                                                   ='
        write(*,*) '=          Model and system parameters can be adjusted in           ='
        write(*,*) '=            User_defined_params.py                                 ='
        write(*,*) '====================================================================='    
    endif 
    
    simuStartTime = MPI_WTIME()
    startTimeStamp = MPI_WTIME()    
    
    call readglobal
    call readmodelgeometry

    allocate(fxmin(ntotft),fxmax(ntotft),fymin(ntotft),fymax(ntotft),fzmin(ntotft),fzmax(ntotft),material(nmat,n2mat))
    allocate(nonfs(ntotft))
    allocate(fltxyz(2,4,ntotft))

    call readfaultgeometry
    call readmaterial
    call readstations1

    allocate(OffFaultStNodeIdIndex(2,totalNumOfOffSt), xonfs(2,maxval(nonfs),ntotft), x4nds(3,totalNumOfOffSt))

    call readstations2
    if (insertFaultType > 0) call read_fault_rough_geometry
    call warning
    
    allocate(nftnd(ntotft),localShapeFunc(nrowsh,nen))
    
    call calcLocalShapeFunc
    call mesh4num
    call allocInit
    call memory_estimate
    call meshgen
    call checkMeshMaterial
    !call checkArrSize  ! disabled for productive runs
    call netcdf_read_on_fault_eqdyna
    if (mode==2) call netcdf_read_on_fault_eqdyna_restart
    if (outputGroundMotion == 1 .or. outputFinalSurfDisp == 1) call find_surfaceNodeIdArr
    
    call allocInitAfterMeshGen
    compTimeInSeconds(1) = MPI_WTIME() - startTimeStamp

    startTimeStamp = MPI_WTIME()
    call assembleGlobalMass
    compTimeInSeconds(2) = MPI_WTIME() - startTimeStamp

    call init_vel ! Initiate on-fault node velocities
    
    call driver
    
    startTimeStamp = MPI_WTIME()
    call output_onfault_st
    call output_offfault_st
    call output_frt
    if (output_plastic == 1) call output_plastic_strain  
    if (outputFinalSurfDisp == 1) call output_finalSurfDisp

    compTimeInSeconds(8) = MPI_WTIME() - startTimeStamp 
    compTimeInSeconds(9) = MPI_WTIME() - simuStartTime 
   
    if (writeCompTime == 1) call output_timeanalysis
    
    call MPI_Finalize(iMPIerr)
    stop
    
end program EQdyna

subroutine allocInit
    use globalvar 
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

    nftmx = maxval(nftnd) !max fault nodel num for all faults, used for arrays.
    if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
    nonmx = sum(nonfs)    !max possible on-fault stations number

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
end subroutine allocInit

subroutine allocInitAfterMeshGen
    use globalvar 
    implicit none 
    integer (kind = 4) :: iSt, iDof, dispOrVel, rowCount, nodeId
    if(numOfOnFaultStCount<=0) numOfOnFaultStCount=1 
    allocate(onFaultQuantHistSCECForm(12,nstep,numOfOnFaultStCount))
    onFaultQuantHistSCECForm = 0.0d0

    allocate(nodalForceArr(totalNumOfEquations), v1(totalNumOfEquations), &
            nodalMassArr(totalNumOfEquations), &
            velArr(ndof,totalNumOfNodes), dispArr(ndof,totalNumOfNodes))
    nodalForceArr = 0.0d0
    nodalMassArr  = 0.0d0
    v1      = 0.0d0
    velArr  = 0.0d0
    dispArr = 0.0d0

    if(numOfOffFaultStCount>0) then 
        allocate(idhist(3,numOfOffFaultStCount*ndof*2), &
                OffFaultStGramSCEC(numOfOffFaultStCount*ndof*2+1,nstep))
        idhist = 0
        OffFaultStGramSCEC = 0.0d0
        rowCount = 0
        do iSt = 1, numOfOffFaultStCount
            do iDof = 1,ndof
                do dispOrVel = 1, 2
                    rowCount = rowCount + 1
                    nodeId = OffFaultStNodeIdIndex(2,iSt)
                    idhist(1,rowCount) =  nodeId
                    !if(idhist(1,l)<=0) idhist(1,l)=1  !avoid zero that cannot be in array below
                    idhist(2,rowCount) = iDof  
                    idhist(3,rowCount) = dispOrVel    
                enddo
            enddo
        enddo            
    endif
end subroutine allocInitAfterMeshGen

subroutine init_vel
    ! initiate the 1d velocity array v1. 
    ! if mode==2, non-zero values for fric(31-36,i,ift) loaded from 
    !    the restart file.
    ! if mode==1, fric(31-36,i) will be zeros.
    use globalvar
    implicit none
    integer (kind = 4) :: i, ift, tmp

    do ift = 1, ntotft
        do i = 1,nftnd(ift)
            tmp = eqNumStartIndexLoc(nsmp(1,i,ift))! slave nodeid i 
            v1(eqNumIndexArr(tmp+1)) = fric(34,i,ift) ! vxs
            v1(eqNumIndexArr(tmp+2)) = fric(35,i,ift) ! vys
            v1(eqNumIndexArr(tmp+3)) = fric(36,i,ift) ! vzs
            tmp = eqNumStartIndexLoc(nsmp(2,i,ift))! master nodeid i
            v1(eqNumIndexArr(tmp+1)) = fric(31,i,ift) ! vxm
            v1(eqNumIndexArr(tmp+2)) = fric(32,i,ift) ! vym
            v1(eqNumIndexArr(tmp+3)) = fric(33,i,ift) ! vzm
        enddo
    enddo
end subroutine init_vel

subroutine checkMeshMaterial
    use globalvar
    implicit none
    integer (kind = 4) :: i
    do i = 1, totalNumOfElements
        if (mat(i,1) == 0.0d3 .or. mat(i,2) == 0.0d3 .or. mat(i,3) == 0.0d3) then
            write(*,*) 'Element ', i, ' is not assigned material property. Exiting ... ...'
            stop 
        endif 
    enddo 
end subroutine checkMeshMaterial

subroutine checkArrSize
    use globalvar
    implicit none

    write(*,*) 'EqNumIndexArr size is ', sizeOfEqNumIndexArr
    write(*,*) 'Size of stress dof index array is ', sizeOfStressDofIndexArr
    write(*,*) 'Is', 5*sizeOfEqNumIndexArr, 'too large compared to ', sizeOfStressDofIndexArr
end subroutine checkArrSize
