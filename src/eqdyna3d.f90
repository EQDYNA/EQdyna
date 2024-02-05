! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
program EQdyna
    use globalvar
    implicit none
    include 'mpif.h'
        
    integer (kind = 4) :: i, j, k, l, itmp, alloc_err, ierr

    call MPI_Init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,me,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)

    if (me == master) then 
        write(*,*) '====================================================================='
        write(*,*) '==================   Welcome to EQdyna 5.3.2  ======================='
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
    
    timebegin=MPI_WTIME()
    time1=MPI_WTIME()    

    write(mm,'(i6)') me
    mm=trim(adjustl(mm))
    
    call readglobal
    call readmodelgeometry

    allocate(fxmin(ntotft),fxmax(ntotft),fymin(ntotft),fymax(ntotft),fzmin(ntotft),fzmax(ntotft),material(nmat,n2mat))
    allocate(nonfs(ntotft))
    allocate(fltxyz(2,4,ntotft))

    call readfaultgeometry
    call readmaterial
    call readstations1

    itmp = maxval(nonfs)
    allocate(an4nds(2,n4nds), xonfs(2,itmp,ntotft), x4nds(3,n4nds))

    call readstations2
    if (insertFaultType > 0) call read_fault_rough_geometry
    call warning
    
    nplpts=0    !initialize number of time history plot
    if (nhplt>0) then
        nplpts=int(nstep/nhplt)+2
    endif

    allocate(nftnd(ntotft),shl(nrowsh,nen))
    
    call qdcshl

    call mesh4num
    call allocInit
    call memory_estimate
    call meshgen
    call checkMeshMaterial

    call netcdf_read_on_fault_eqdyna
    if (mode==2) call netcdf_read_on_fault_eqdyna_restart
    
    if (outputGroundMotion == 1) call find_surface_node_id
    
    if (C_degen == 1) then 
        do i = 1, numnp
           if (x(2,i)>dx/2.0d0) then 
                x(1,i) = x(1,i) - dx/2.0d0
           endif
           if (x(2,i)<-dx/2.0d0) then 
                x(1,i) = x(1,i) + dx/2.0d0
           endif 
        enddo 
    endif
    
    if(n4onf<=0) n4onf=1 
    allocate(fltsta(12,nplpts-1,n4onf),stat=alloc_err)
    fltsta = 0.0d0

    allocate(nodalForceArr(totalNumOfEquations),v1(totalNumOfEquations),d1(totalNumOfEquations), nodalMassArr(totalNumOfEquations), v(ndof,numnp),d(ndof,numnp),stat=alloc_err)

    nodalForceArr    = 0.0d0
    nodalMassArr    = 0.0d0
    v1      = 0.0d0
    d1      = 0.0d0
    v       = 0.0d0
    d       = 0.0d0

    allocate(frichis(2,nftmx,nplpts,ntotft))
    frichis = 0.0d0

    if(n4out>0) then 
        ndout=n4out*ndof*noid!3 components of 2 quantities: v and d
        !   write(*,*) 'ndout= ',ndout    
        allocate(idhist(3,ndout),dout(ndout+1,nplpts),stat=alloc_err)
        if(alloc_err /=0) then
            write(*,*) 'me= ',me,'insufficient space to allocate array idhist or dout'
        endif

        idhist=0
        dout=0.0d0
        l=0
        do i=1,n4out
            do j=1,ndof
                do k=1,noid
                    l = l + 1
                    idhist(1,l) = an4nds(2,i) !node number (>1, <NUMNP)
                    if(idhist(1,l)<=0) idhist(1,l)=1  !avoid zero that cannot be in array below
                    idhist(2,l) = j    !degree of freedom number (<=NDOF)
                    idhist(3,l) = k    !kinematic quantity specifier 
                    !(disp, vel, or acc)
                enddo
            enddo
        enddo            
    endif
    
    time2=MPI_WTIME()        
    timeused(1)=time2-time1

    time1 = MPI_WTIME()
    call qdct2
    timeused(2) = timeused(2) + MPI_WTIME() - time1

    call init_vel ! Initiate on-fault node velocities.
    
    call driver
    
    if (me == master) then 
        write(*,*) '=                                                                   ='
        write(*,*) '=                      Writing out results                          ='
        write(*,*) '====================================================================='    
    endif
    
    time1=MPI_WTIME()
    
    call output_onfault_st
    call output_offfault_st
    call output_frt

    time2=MPI_WTIME()
    timeused(8)=time2-time1 
    timeover=MPI_WTIME()
    timeused(9)=timeover-timebegin 
    
    if (timeinfo == 1) call output_timeanalysis
    if (output_plastic == 1) call output_plastic_strain
    
    call MPI_Finalize(ierr)
    stop
    
end program EQdyna

subroutine allocInit
    use globalvar 
    implicit none 
    allocate(equationNumIndexArr(maxm),locateEqNumStartIndex(numnp),dof1(numnp),x(ndof,numnp), fnms(numnp), surface_node_id(numnp)) 

    allocate(ien(nen,totalNumOfElements), mat(totalNumOfElements,5), et(totalNumOfElements), eleporep(totalNumOfElements), pstrain(totalNumOfElements), &
                eledet(totalNumOfElements), elemass(nee,totalNumOfElements), eleshp(nrowsh-1,nen,totalNumOfElements), &
                ss(6,totalNumOfElements), phi(nen,4,totalNumOfElements))

    eleporep = 0.0d0
    pstrain  = 0.0d0
    eledet   = 0.0d0
    elemass  = 0.0d0

    nftmx=maxval(nftnd) !max fault nodel num for all faults, used for arrays.
    if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
    nonmx=sum(nonfs)    !max possible on-fault stations number

    allocate(nsmp(2,nftmx,ntotft), fnft(nftmx,ntotft), un(3,nftmx,ntotft),&
                us(3,nftmx,ntotft), ud(3,nftmx,ntotft), fric(100,nftmx,ntotft),&
                arn(nftmx,ntotft), r4nuc(nftmx,ntotft), anonfs(3,nonmx),&
                arn4m(nftmx,ntotft), state(nftmx,ntotft), fltgm(nftmx),&
                Tatnode(nftmx,ntotft), patnode(nftmx,ntotft))

    fltgm   = 0  
    nsmp    = 0    
    fnft    = 99999.0d0 
    fric    = 0.0d0
    un      = 0.0d0
    us      = 1000.0d0
    ud      = 0.0d0
    arn     = 0.0d0
    arn4m   = 0.0d0
    r4nuc   = 0.0d0
    anonfs  = 0
    state   = 0.0d0
    Tatnode = 0.0d0 
    patnode = 0.0d0
        
    allocate(stressCompIndexArr(totalNumOfElements))
    allocate(s1(5*maxm))
    s1      = 0.0d0
end subroutine allocInit

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