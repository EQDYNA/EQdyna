!/* Copyright (C) 2006-2023, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
subroutine meshgen
    ! Create regular node grids and hexahedral elements for this MPI process.
    ! Create split nodes.
    ! Distort the regular mesh for complex fault geometries.
    ! Set material propertiesa, initial stress, and other element-wise properties. 
    
    use globalvar
    implicit none
    include 'mpif.h'
    ! incremental variables 
    integer (kind = 4) :: nnode, nelement, neq0, ntag, ntags
    integer (kind = 4) :: nxt,nyt,nzt,nx,ny,nz,ix,iy,iz, &
                       edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,edgezn, &
                       nxuni,nyuni,nzuni,ift, &
                       n1,n2,n3,n4,m1,m2,m3,m4,&
                       mex,mey,mez,itmp1,&
                       nodeDofNum, msnode, nodeXyzIndex(10)
    integer (kind = 4), dimension(ntotft) :: nftnd0,ixfi,izfi,ifs,ifd
    integer (kind = 4), allocatable :: fltrc(:,:,:,:)
    ! Temporary real variables
    real (kind = dp) :: nodeCoor(10), elementCenterCoor(3), &
                       a,b,area,aa1,bb1,cc1,dd1,p1,q1, ycoort, pfx = 0.0d0, pfz = 0.0d0
    real (kind = dp), allocatable :: xlinet(:), ylinet(:), zlinet(:), &
                                     xline(:),  yline(:) , zline(:)

    write(mm,'(i6)') me
    mm = trim(adjustl(mm))
    !dy=dx
    !dz=dx

    call calcXyzMPIId(mex, mey, mez)
    
    ! get xline, yline, zline
    call getSize1DCoorX(nxt, nxuni, edgex1)
        allocate(xlinet(nxt))
    call get1DCoorX(mex, nxt, nxuni, edgex1, xlinet, nx)
        allocate(xline(nx))
    call get1DCoorXLocal(mex, nxt, nx, xline, xlinet)
    
    call getSize1DCoorY(nyt, nyuni, edgey1)
        allocate(ylinet(nyt))
    call get1DCoorY(mey, nyt, nyuni, edgey1, ylinet, ny)
        allocate(yline(ny))
    call get1DCoorYLocal(mey, nyt, ny, yline, ylinet)
    
    call getSize1DCoorZ(nzt, nzuni, edgezn)
        allocate(zlinet(nzt))
    call get1DCoorZ(mez, nzt, nzuni, edgezn, zlinet, nz)
        allocate(zline(nz))
    call get1DCoorZLocal(mez, nzt, nz, zline, zlinet)

    ! Move adjacent plane1 and plane2 to create hexahedral meshes
    allocate(plane1(ny+ntotft,nz),plane2(ny+ntotft,nz),fltrc(2,nxuni,nzuni,ntotft))
    plane1 = 0
    plane2 = 0
    
    numcount    = 0
    numcount(1) = nx
    numcount(2) = ny
    numcount(3) = nz
    
    ! Initialize scalars

    msnode   = nx*ny*nz
    n4onf    = 0
    n4out    = 0
    an4nds   = 0
    ! Initialize arrays
    nftnd0   = 0
    x        = 0.0d0
    ien      = 0
    id1      = 0
    ids      = 0
    locid    = 0
    dof1     = 0
    et       = 0
    
    ixfi = 0
    izfi = 0
    ! increments
    nnode    = 0
    nelement = 0
    neq0     = 0
    ntag     = 0
    ntags    = 0
    
    allocate(n4yn(n4nds))
    n4yn = 0
    
    ! Loop over x,y,z grids to create nodes and elements
    ! Create node
    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                call initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
                call createNode(nodeCoor, xline(ix), yline(iy), zline(iz), nnode, nodeXyzIndex)
                if (insertFaultType > 0) then 
                    call insertFaultInterface(nodeCoor, ycoort, pfx, pfz)
                    x(2,nnode) = ycoort
                endif 
                
                call setNumDof(nodeCoor, nodeDofNum)

                locid(nnode) = ntag
                dof1(nnode)  = nodeDofNum
                
                call setEquationNumber(nodeXyzIndex, nodeCoor, ntag, neq0, nodeDofNum)
                call setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nnode)
                call createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nnode, msnode, nftnd0, neq0, ntag, &
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
                
                ! Create element
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                    call createElement(nelement, ntags, iy, iz, elementCenterCoor)
                    call setElementMaterial(nelement, elementCenterCoor)
                    
                    if (C_degen == 1) then 
                        call wedge(elementCenterCoor(1), elementCenterCoor(2), elementCenterCoor(3), nelement, ntags, iy, iz, nftnd0(1))
                    elseif (C_degen == 2) then 
                        call tetra(elementCenterCoor(1), elementCenterCoor(2), elementCenterCoor(3), nelement, ntags, iy, iz, nftnd0(1))
                    endif         
                    
                    call replaceSlaveWithMasterNode(nodeCoor, nelement, nftnd0) 
                    if (C_elastic==0 .and. TPV==2802) call plastic_set_mat_stress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0, nelement)          
                endif!if element
            enddo!iy
        enddo!iz
        plane1 = plane2
    enddo!ix
    
    maxs=ntags
    
    call meshGenError(nx, ny, nz, nnode, msnode, nelement, neq0, ntag, nftnd0)
    
    ! compute on-fault area associated with each fault node pair and distance from source
    do ift=1,ntotft
        if(nftnd0(ift)>0) then
        !...distance from the source point
        do i=1,ifd(ift)
            do j=1,ifs(ift)
                i1 = fltrc(1,j,i,ift)
                j1 = fltrc(2,j,i,ift)
                r4nuc(j1,ift) = sqrt((x(1,i1)-xsource)**2 + (x(2,i1)-ysource)**2 &
                + (x(3,i1)-zsource)**2)
            enddo
        enddo
        !...element areas and distribute evenly to its four nodes
        do i=2,ifd(ift)
            do j=2,ifs(ift)
            !...4 nodes of quadrilateral
                n1 = fltrc(1,j,i,ift) !use nodal number in nnode
                n2 = fltrc(1,j-1,i,ift)
                n3 = fltrc(1,j-1,i-1,ift)
                n4 = fltrc(1,j,i-1,ift)
                m1 = fltrc(2,j,i,ift) !use nodal number in nftnd0
                m2 = fltrc(2,j-1,i,ift)
                m3 = fltrc(2,j-1,i-1,ift)
                m4 = fltrc(2,j,i-1,ift)
                !...calculate area of the quadrilateral
                !......if fault is not in coor axes plane
                aa1=sqrt((x(1,n2)-x(1,n1))**2 + (x(2,n2)-x(2,n1))**2 &
                + (x(3,n2)-x(3,n1))**2)
                bb1=sqrt((x(1,n3)-x(1,n2))**2 + (x(2,n3)-x(2,n2))**2 &
                + (x(3,n3)-x(3,n2))**2)
                cc1=sqrt((x(1,n4)-x(1,n3))**2 + (x(2,n4)-x(2,n3))**2 &
                + (x(3,n4)-x(3,n3))**2)
                dd1=sqrt((x(1,n1)-x(1,n4))**2 + (x(2,n1)-x(2,n4))**2 &
                + (x(3,n1)-x(3,n4))**2)
                p1=sqrt((x(1,n4)-x(1,n2))**2 + (x(2,n4)-x(2,n2))**2 &
                + (x(3,n4)-x(3,n2))**2)
                q1=sqrt((x(1,n3)-x(1,n1))**2 + (x(2,n3)-x(2,n1))**2 &
                + (x(3,n3)-x(3,n1))**2)
                area=0.25 * sqrt(4*p1*p1*q1*q1 - &
               (bb1*bb1 + dd1*dd1 - aa1*aa1 -cc1*cc1)**2) 
                !...distribute above area to 4 nodes evenly
                area = 0.25 * area
                arn(m1,ift) = arn(m1,ift) + area
                arn(m2,ift) = arn(m2,ift) + area
                arn(m3,ift) = arn(m3,ift) + area
                arn(m4,ift) = arn(m4,ift) + area
                enddo
            enddo
            !...save area for moment (rate) calculation before MPI. B.D. 8/11/10
            arn4m = arn
        endif!end if nftnd0=0/ not
        !
        time1 = MPI_WTIME()
        
        call MPI4arn(nx, ny, nz, mex, mey, mez, nftnd0(ift), ift)
       
    enddo !ift=i:nft 

    time2 = MPI_WTIME()
    btime = btime + (time2 - time1) 
    
end subroutine meshgen
!==================================================================================================
!**************************************************************************************************
!==================================================================================================

subroutine setElementMaterial(nelement, elementCenterCoor)
! Subroutine velocityStructure will asign Vp, Vs and rho 
!   based on input from bMaterial.txt, which is created by 
!   case input file user_defined_param.py.

    use globalvar
    implicit none
    integer (kind = 4) :: nelement, i
    real (kind = dp) :: elementCenterCoor(3)
    
    if (nmat==1 .and. n2mat == 3) then
    ! homogenous material
        mat(nelement,1)  = material(1,1)
        mat(nelement,2)  = material(1,2)
        mat(nelement,3)  = material(1,3)
    elseif (nmat>1 .and. n2mat == 4) then 
        ! 1D velocity structure
        ! material = material(i,j), i=1,nmat, and j=1,4
        ! for j
        !   1: bottom depth of a layer (should be positive in m).
        !   2: Vp, m/s
        !   3: Vs, m/s
        !   4: rho, kg/m3
        if (abs(elementCenterCoor(3)) < material(1,1)) then
            mat(nelement,1)  = material(1,2)
            mat(nelement,2)  = material(1,3)
            mat(nelement,3)  = material(1,4)
        else
            do i = 2, nmat
                if (abs(elementCenterCoor(3)) < material(i,1) &
                    .and. abs(elementCenterCoor(3)) >= material(i-1,1)) then
                    
                    mat(nelement,1)  = material(i,2)
                    mat(nelement,2)  = material(i,3)
                    mat(nelement,3)  = material(i,4)
                endif 
            enddo
        endif
    endif 
    
    ! calculate lambda and mu from Vp, Vs and rho.
    ! mu = Vs**2*rho
    mat(nelement,5)  = mat(nelement,2)**2*mat(nelement,3)
    ! lambda = Vp**2*rho-2*mu
    mat(nelement,4)  = mat(nelement,1)**2*mat(nelement,3)-2.0d0*mat(nelement,5)

end subroutine setElementMaterial

subroutine MPI4arn(nx, ny, nz, mex, mey, mez, totalNumFaultNode, iFault)
! Add up arn from neighbor MPI blocks.
    use globalvar
    implicit none 
    include 'mpif.h'
    
    integer (kind = 4) :: bndl,bndr,bndf,bndb,bndd,bndu, nx, ny, nz, itmp1, mex, mey, mez
    integer (kind = 4) :: istatus(MPI_STATUS_SIZE), ierr, totalNumFaultNode, iFault, i
    real(kind = dp),allocatable,dimension(:) :: btmp, btmp1
    ! Initialize fltMPI(6) to .false.
    fltMPI=.false.
    
    if(fltnum(1) /= 0) allocate(fltl(fltnum(1)))
    if(fltnum(2) /= 0) allocate(fltr(fltnum(2)))
    if(fltnum(3) /= 0) allocate(fltf(fltnum(3)))
    if(fltnum(4) /= 0) allocate(fltb(fltnum(4)))
    if(fltnum(5) /= 0) allocate(fltd(fltnum(5)))
    if(fltnum(6) /= 0) allocate(fltu(fltnum(6)))
    
    fltnum = 0
    do i = 1, totalNumFaultNode
        if(mod(fltgm(i),10)==1 ) then 
            fltnum(1) = fltnum(1) + 1
            fltl(fltnum(1)) = i
        endif
        if(mod(fltgm(i),10)==2 ) then 
            fltnum(2) = fltnum(2) + 1
            fltr(fltnum(2)) = i
        endif
        if(mod(fltgm(i),100)-mod(fltgm(i),10)==10 ) then
            fltnum(3) = fltnum(3) + 1
            fltf(fltnum(3)) = i
        endif
        if(mod(fltgm(i),100)-mod(fltgm(i),10)==20 ) then
            fltnum(4) = fltnum(4) + 1
            fltb(fltnum(4)) = i
        endif
        if(fltgm(i)-mod(fltgm(i),100)==100 ) then
            fltnum(5) = fltnum(5) + 1
            fltd(fltnum(5)) = i
        endif
        if(fltgm(i)-mod(fltgm(i),100)==200 ) then
            fltnum(6) = fltnum(6) + 1
            fltu(fltnum(6)) = i
        endif
    enddo
    if (npx > 1) then
        bndl=1  !left boundary
        bndr=nx ! right boundary
        if (mex == master) then
            bndl=0
        elseif (mex == npx-1) then
            bndr=0
        endif
 
        if (bndl/=0) then
            call mpi_sendrecv(fltnum(1), 1, MPI_INTEGER, me-npy*npz, 1000+me, &
                itmp1, 1, MPI_INTEGER, me-npy*npz, 1000+me-npy*npz, &
                MPI_COMM_WORLD, istatus, ierr)
            !write(*,*) 'me= ',me, 'itmp1=',itmp1                
            if(fltnum(1)>0 .and. itmp1>0 ) then  
                if(fltnum(1) /= itmp1) stop 'error in fltnum(1)'
                fltMPI(1)=.true.
                allocate(btmp(fltnum(1)),btmp1(fltnum(1)))
                btmp = 0.
                btmp1 = 0.
                do i = 1, fltnum(1)
                    btmp(i)=arn(fltl(i),iFault)
                enddo
                call mpi_sendrecv(btmp, fltnum(1), MPI_DOUBLE_PRECISION, me-npy*npz, 1000+me, &
                    btmp1, fltnum(1), MPI_DOUBLE_PRECISION, me-npy*npz, 1000+me-npy*npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                do i = 1, fltnum(1)
                    arn(fltl(i),iFault) = arn(fltl(i),iFault) + btmp1(i)
                enddo
                deallocate(btmp,btmp1)       
            endif
        endif !if bhdl/=0

        if (bndr/=0) then
            call mpi_sendrecv(fltnum(2), 1, MPI_INTEGER, me+npy*npz, 1000+me, &
                itmp1, 1, MPI_INTEGER, me+npy*npz, 1000+me+npy*npz, &
                MPI_COMM_WORLD, istatus, ierr)
            !write(*,*) 'me= ',me, 'itmp1=',itmp1                
            if(fltnum(2)>0 .and. itmp1>0 ) then  
                if(fltnum(2) /= itmp1) stop 'error in fltnum(2)'
                fltMPI(2)=.true.
                allocate(btmp(fltnum(2)),btmp1(fltnum(2)))
                btmp  = 0.
                btmp1 = 0.
                do i = 1, fltnum(2)
                    btmp(i)=arn(fltr(i),iFault)
                enddo
                call mpi_sendrecv(btmp, fltnum(2), MPI_DOUBLE_PRECISION, me+npy*npz, 1000+me, &
                    btmp1, fltnum(2), MPI_DOUBLE_PRECISION, me+npy*npz, 1000+me+npy*npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                do i = 1, fltnum(2)
                    arn(fltr(i),iFault) = arn(fltr(i),iFault) + btmp1(i)
                enddo
                deallocate(btmp,btmp1)       
            endif
        endif !bndr/=0
    endif !npx>1
!*****************************************************************************************
    if (npy > 1) then   !MPI:no fault message passing along y due to the x-z vertical fault.
        bndf=1  !front boundary
        bndb=ny ! back boundary
        if (mey == master) then
            bndf=0
        elseif (mey == npy-1) then
            bndb=0
        endif

        if (bndf/=0) then
            call mpi_sendrecv(fltnum(3), 1, MPI_INTEGER, me-npz, 2000+me, &
                itmp1, 1, MPI_INTEGER, me-npz, 2000+me-npz, &
                MPI_COMM_WORLD, istatus, ierr)
            !write(*,*) 'me= ',me, 'itmp1=',itmp1                
            if(fltnum(3)>0 .and. itmp1>0) then
                if(fltnum(3) /= itmp1) stop 'error in fltnum(3)'
                fltMPI(3)=.true.
                allocate(btmp(fltnum(3)),btmp1(fltnum(3)))
                btmp = 0.
                btmp1 = 0.
                do i = 1, fltnum(3)
                    btmp(i)=arn(fltf(i),iFault)
                enddo
                call mpi_sendrecv(btmp, fltnum(3), MPI_DOUBLE_PRECISION, me-npz, 2000+me, &
                    btmp1, fltnum(3), MPI_DOUBLE_PRECISION, me-npz, 2000+me-npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                do i = 1, fltnum(3)
                    arn(fltf(i),iFault) = arn(fltf(i),iFault) + btmp1(i)
                enddo
                deallocate(btmp,btmp1)       
            endif
        endif !bhdf/=0

        if (bndb/=0) then
            call mpi_sendrecv(fltnum(4), 1, MPI_INTEGER, me+npz, 2000+me, &
                itmp1, 1, MPI_INTEGER, me+npz, 2000+me+npz, &
                MPI_COMM_WORLD, istatus, ierr)
            !write(*,*) 'me= ',me, 'itmp1=',itmp1                
            if(fltnum(4)>0 .and. itmp1>0) then  
                if(fltnum(4) /= itmp1) stop 'error in fltnum(4)'
                fltMPI(4)=.true.
                allocate(btmp(fltnum(4)),btmp1(fltnum(4)))
                btmp  = 0.
                btmp1 = 0.
                do i = 1, fltnum(4)
                    btmp(i)=arn(fltb(i),iFault)
                enddo
                call mpi_sendrecv(btmp, fltnum(4), MPI_DOUBLE_PRECISION, me+npz, 2000+me, &
                    btmp1, fltnum(4), MPI_DOUBLE_PRECISION, me+npz, 2000+me+npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                do i = 1, fltnum(4)
                    arn(fltb(i),iFault) = arn(fltb(i),iFault) + btmp1(i)
                enddo
                deallocate(btmp,btmp1)       
            endif
        endif !bndb/=0
     endif !npy>1
!*****************************************************************************************
    if (npz > 1) then
        bndd=1  !lower(down) boundary
        bndu=nz ! upper boundary
        if (mez == master) then
            bndd=0
        elseif (mez == npz-1) then
            bndu=0
        endif
        if (bndd/=0) then
            call mpi_sendrecv(fltnum(5), 1, MPI_INTEGER, me-1, 3000+me, &
                itmp1, 1, MPI_INTEGER, me-1, 3000+me-1, &
                MPI_COMM_WORLD, istatus, ierr)
            !write(*,*) 'me= ',me, 'itmp1=',itmp1                
            if(fltnum(5)>0 .and. itmp1>0) then
                if(fltnum(5) /= itmp1) stop 'error in fltnum(5)'
                fltMPI(5)=.true.
                allocate(btmp(fltnum(5)),btmp1(fltnum(5)))
                btmp = 0.
                btmp1 = 0.
                do i = 1, fltnum(5)
                    btmp(i)=arn(fltd(i),iFault)
                enddo
                call mpi_sendrecv(btmp, fltnum(5), MPI_DOUBLE_PRECISION, me-1, 3000+me, &
                    btmp1, fltnum(5), MPI_DOUBLE_PRECISION, me-1, 3000+me-1, &
                    MPI_COMM_WORLD, istatus, ierr)
                do i = 1, fltnum(5)
                    arn(fltd(i),iFault) = arn(fltd(i),iFault) + btmp1(i)
                enddo
                deallocate(btmp,btmp1)       
            endif
        endif !bhdd/=0
        
        if (bndu/=0) then
            call mpi_sendrecv(fltnum(6), 1, MPI_INTEGER, me+1, 3000+me, &
                itmp1, 1, MPI_INTEGER, me+1, 3000+me+1, &
                MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(6)>0 .and. itmp1>0) then  
                    if(fltnum(6) /= itmp1) stop 'error in fltnum(6)'
                    fltMPI(6)=.true.
                    allocate(btmp(fltnum(6)),btmp1(fltnum(6)))
                    btmp  = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(6)
                        btmp(i)=arn(fltu(i),iFault)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(6), MPI_DOUBLE_PRECISION, me+1, 3000+me, &
                        btmp1, fltnum(6), MPI_DOUBLE_PRECISION, me+1, 3000+me+1, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(6)
                    arn(fltu(i),iFault) = arn(fltu(i),iFault) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
        endif !bndu/=0
    endif !npz>1
end subroutine MPI4arn 

subroutine meshGenError(nx, ny, nz, nnode, msnode, nelement, neq0, ntag, nftnd0)
! Check consistency between mesh4 and meshgen
    use globalvar
    implicit none
    integer (kind = 4) :: nx, ny, nz, nnode, msnode, nelement, neq0, nftnd0(ntotft), ntag
    integer (kind = 4) :: i
    if (maxs>=(5*maxm)) then
        write(*,*) '5*maxm',maxm,'is not enough for maxs',maxs
        stop 2002
    endif
    if(nnode/=nx*ny*nz.or.msnode/=numnp.or.nelement/=numel.or.neq0/=neq) then
        write(*,*) 'Inconsistency in node/element/equation/between meshgen and mesh4num: stop!',me
        write(*,*) 'nnode&numnp=',nnode,numnp
        write(*,*) 'nelement,numel=',nelement,numel
        write(*,*) 'neq0,neq=',neq0,neq
        write(*,*) 'nnode,nx,ny,nz',nnode,nx,ny,nz
        write(*,*) 'msnode,numnp',msnode,numnp
        stop 2003
    endif
    if(ntag/=maxm) then
        write(*,*) 'Inconsistency in ntag and maxm: stop!',me
        write(*,*) ntag,maxm
        stop 2004
    endif
    do i=1,ntotft
        if(nftnd0(i)/=nftnd(i)) then
            write(*,*) 'Inconsistency in fault between meshgen and mesh4num: stop!',me,i
            stop 2005
        endif
    enddo
end subroutine meshGenError

subroutine calcXyzMPIId(mex, mey, mez)
    use globalvar 
    implicit none
    integer (kind = 4) :: mex, mey, mez
     
    mex=int(me/(npy*npz))
    mey=int((me-mex*npy*npz)/npz)
    mez=int(me-mex*npy*npz-mey*npz)
end subroutine calcXyzMPIId

subroutine getSize1DCoorX(nxt, nxuni, edgex1)
    use globalvar
    implicit none
    integer (kind = 4) :: nxuni, ix, edgex1, nxt
    real (kind = dp) :: xstep, xcoor
    
    
    nxuni = (fltxyz(2,1,1) - fltxyz(1,1,1)) / dx + 1
    xstep = dx
    xcoor = fltxyz(1,1,1)
    do ix = 1, np
        xstep = xstep * rat
        xcoor = xcoor - xstep
        if(xcoor <= xmin) exit
    enddo
    edgex1 = ix + nPML
    xstep = dx
    xcoor = fltxyz(2,1,1)
    do ix = 1, np
        xstep = xstep * rat
        xcoor = xcoor + xstep
        if(xcoor >= xmax) exit
    enddo
    nxt = nxuni + edgex1 + ix + nPML
    
    end subroutine getSize1DCoorX

subroutine get1DCoorX(mex, nxt, nxuni, edgex1, xlinet, nx)
    use globalvar
    implicit none
    integer (kind = 4) :: nxuni, nxt, edgex1, nx, j1, rlp, rr, mex, ix
    real (kind = dp) :: xlinet(nxt), xstep
    
    xlinet(edgex1+1) = fltxyz(1,1,1)
    xstep = dx
    do ix = edgex1, 1, -1
        xstep = xstep * rat
        xlinet(ix) = xlinet(ix+1) - xstep
    enddo
    xmin = xlinet(1)
    do ix = edgex1+2,edgex1+nxuni
        xlinet(ix) = xlinet(ix-1) + dx
    enddo
    xstep = dx
    do ix = edgex1+nxuni+1,nxt
        xstep = xstep * rat
        xlinet(ix) = xlinet(ix-1) + xstep
    enddo
    xmax = xlinet(nxt)

    j1 = nxt + npx - 1
    rlp = j1/npx
    rr = j1 - rlp*npx
    if(mex<(npx-rr)) then
        nx = rlp
    else
        nx = rlp + 1    !evenly distributed to last rr
    endif

end subroutine get1DCoorX

subroutine get1DCoorXLocal(mex, nxt, nx, xline, xlinet)
    use globalvar
    implicit none
    integer (kind = 4) :: mex, j1, rlp, rr, nxt, nx, k1, ix
    real (kind = dp) :: xline(nx), xlinet(nxt)
    
    j1  = nxt + npx - 1
    rlp = j1/npx
    rr  = j1 - rlp*npx
    if(mex<=npx-rr) then
        do ix=1,nx
          xline(ix) = xlinet(rlp*mex-mex+ix)
        enddo
    else
        do ix=1,nx
            k1 = mex - (npx - rr)
            xline(ix) = xlinet(rlp*mex-mex+k1+ix)
        enddo
    endif    
end subroutine get1DCoorXLocal

subroutine getSize1DCoorY(nyt, nyuni, edgey1)
    use globalvar
    implicit none
    integer (kind = 4) :: nyuni, iy, edgey1, nyt
    real (kind = dp) :: ystep, ycoor
    nyuni = dis4uniF + dis4uniB + 1
    ystep = dy
    !ycoor = -dy*(dis4uniF+fltxyz(2,1,1)/dx+1)
    ycoor = -dy*(dis4uniF)
    do iy = 1, np
        ystep = ystep * rat
        ycoor = ycoor - ystep
        if(ycoor <= ymin) exit
    enddo
    edgey1 = iy + nPML
    ystep = dy
    ycoor = dy*(dis4uniB)
    do iy = 1, np
        ystep = ystep * rat
        ycoor = ycoor + ystep
        if(ycoor >= ymax) exit
    enddo
    nyt = nyuni + edgey1 + iy + nPML
    
end subroutine getSize1DCoorY

subroutine get1DCoorY(mey, nyt, nyuni, edgey1, ylinet, ny)
    use globalvar
    implicit none
    integer (kind = 4) :: nyuni, nyt, edgey1, ny, j1, rlp, rr, mey, iy
    real (kind = dp) :: ylinet(nyt), ystep

    ylinet(edgey1+1) = -dy*(dis4uniF)
    ystep = dy
    do iy = edgey1, 1, -1
        ystep = ystep * rat    
        ylinet(iy) = ylinet(iy+1) - ystep
    enddo
    ymin = ylinet(1)
    do iy = edgey1+2,edgey1+nyuni
        ylinet(iy) = ylinet(iy-1) + dy
    enddo
    ystep = dy
    do iy = edgey1+nyuni+1,nyt
        ystep = ystep * rat
        ylinet(iy) = ylinet(iy-1) + ystep
    enddo
    ymax = ylinet(nyt)

    j1 = nyt + npy - 1
    rlp = j1/npy
    rr = j1 - rlp*npy
    if(mey<(npy-rr)) then
        ny = rlp
    else
        ny = rlp + 1    
    endif
end subroutine get1DCoorY

subroutine get1DCoorYLocal(mey, nyt, ny, yline, ylinet)
    use globalvar
    implicit none
    integer (kind = 4) :: mey, j1, rlp, rr, nyt, ny, k1, iy
    real (kind = dp) :: yline(ny), ylinet(nyt)
    
    j1  = nyt + npy - 1
    rlp = j1/npy
    rr  = j1 - rlp*npy
    if(mey<=npy-rr) then
        do iy=1,ny
            yline(iy) = ylinet(rlp*mey-mey+iy)
        enddo
    else
        do iy=1,ny
            k1 = mey - (npy - rr)
            yline(iy) = ylinet(rlp*mey-mey+k1+iy)
        enddo
    endif
end subroutine get1DCoorYLocal

subroutine getSize1DCoorZ(nzt, nzuni, edgezn)
    use globalvar
    implicit none
    integer (kind = 4) :: nzt, nzuni, edgezn, iz
    real (kind = dp) :: zstep, zcoor
    
    zstep = dz
    zcoor = fltxyz(1,3,1)
    do iz=1,np
        zstep = zstep * rat
        zcoor = zcoor - zstep
        if(zcoor <= zmin) exit
    enddo
    edgezn = iz + nPML
    nzuni = (fltxyz(2,3,1)-fltxyz(1,3,1))/dz + 1 
    nzt = edgezn + nzuni
end subroutine getSize1DCoorZ
    
subroutine get1DCoorZ(mez, nzt, nzuni, edgezn, zlinet, nz)
    use globalvar
    implicit none
    integer (kind = 4) :: nzuni, nzt, edgezn, nz, j1, rlp, rr, mez, iz
    real (kind = dp) :: zlinet(nzt), zstep

    zlinet(nzt) = zmax
    do iz = nzt-1,nzt-nzuni+1,-1
        zlinet(iz) = zlinet(iz+1) - dz
    enddo
    zstep = dz
    do iz = nzt-nzuni,1,-1
        zstep = zstep * rat
        zlinet(iz) = zlinet(iz+1) -zstep
    enddo
    zmin = zlinet(1) 
    
    j1 = nzt + npz - 1
    rlp = j1/npz
    rr = j1 - rlp*npz
    if(mez<(npz-rr)) then
        nz = rlp
    else
        nz = rlp + 1   
    endif
end subroutine get1DCoorZ

subroutine get1DCoorZLocal(mez, nzt, nz, zline, zlinet)
    use globalvar
    implicit none
    integer (kind = 4) :: mez, j1, rlp, rr, nzt, nz, k1, iz
    real (kind = dp) :: zline(nz), zlinet(nzt)
    
    j1  = nzt + npz - 1
    rlp = j1/npz
    rr  = j1 - rlp*npz
    if(mez<=npz-rr) then
        do iz=1,nz
            zline(iz) = zlinet(rlp*mez-mez+iz)
        enddo
    else
        do iz=1,nz
            k1 = mez - (npz - rr)
            zline(iz) = zlinet(rlp*mez-mez+k1+iz)
        enddo
    endif
end subroutine get1DCoorZLocal

subroutine setNumDof(nodeCoor, nodeDofNum)
    use globalvar
    implicit none
    integer (kind = 4) :: nodeDofNum
    real (kind = dp) :: nodeCoor(10)
    nodeDofNum = ndof !Default
    if (nodeCoor(1)>PMLb(1) .or. nodeCoor(1)<PMLb(2) .or. nodeCoor(2)>PMLb(3) &
        .or. nodeCoor(2)<PMLb(4) .or. nodeCoor(3)<PMLb(5)) then
        nodeDofNum = 12 ! Modify if inside PML
    endif   
end subroutine setNumDof

subroutine setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nnode)
    use globalvar
    implicit none
    integer (kind = 4) :: nodeXyzIndex(10), ix, iy, nnode, i
    real (kind = dp) :: nodeCoor(10), xline(nodeXyzIndex(4)), yline(nodeXyzIndex(5))
    
    ix = nodeXyzIndex(1)
    iy = nodeXyzIndex(2)
    !Part1. Stations inside the region.
    if(ix>1.and.ix<nodeXyzIndex(4) .and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,n4nds
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(abs(nodeCoor(1)-x4nds(1,i))<tol .or.&
                    (x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
                    (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) .or. &
                    (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
                    (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            n4out = n4out + 1
                            an4nds(1,n4out) = i
                            an4nds(2,n4out) = nnode
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif
    !...identify output nodes (off-fault)
    !Part2. Stations along ix==1
    !Big Bug!!!iz==nz is no valid for 3D MPI.
    if(ix==1.and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,n4nds
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(abs(nodeCoor(1)-x4nds(1,i))<tol .or. &
                    (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
                    (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            n4out = n4out + 1
                            an4nds(1,n4out) = i
                            an4nds(2,n4out) = nnode
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif
    !...identify output nodes (off-fault)
    !Part3. Stations along ix==nx
    if(ix==nodeXyzIndex(4) .and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,n4nds
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
                    (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            n4out = n4out + 1
                            an4nds(1,n4out) = i
                            an4nds(2,n4out) = nnode
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif    
end subroutine setSurfaceStation

subroutine setEquationNumber(nodeXyzIndex, nodeCoor, ntag, neq0, nodeDofNum)
    use globalvar 
    implicit none
    integer (kind = 4) :: iDof, nodeDofNum, ntag, neq0, nodeXyzIndex(10)
    real (kind = dp) :: nodeCoor(10)
    
    do iDof = 1,nodeDofNum
        if(abs(nodeCoor(1)-xmin)<tol.or.abs(nodeCoor(1)-xmax)<tol.or.abs(nodeCoor(2)-ymin)<tol &
        .or.abs(nodeCoor(2)-ymax)<tol.or.abs(nodeCoor(3)-zmin)<tol) then
            ntag      = ntag+1
            id1(ntag) = -1 
            ! Dof = -1 for fixed boundary nodes, no eq #. 
        else
            neq0      = neq0 + 1
            ntag      = ntag+1
            id1(ntag) = neq0
            
            !Count # of DOF on MPI boundaries
            if (nodeXyzIndex(1)==1) then !Left
                numcount(4)=numcount(4)+1
            endif
            if (nodeXyzIndex(1)==nodeXyzIndex(4)) then !Right
                numcount(5)=numcount(5)+1
            endif
            if (nodeXyzIndex(2)==1) then !Front
                numcount(6)=numcount(6)+1
            endif
            if (nodeXyzIndex(2)==nodeXyzIndex(5)) then !Back
                numcount(7)=numcount(7)+1
            endif
            if (nodeXyzIndex(3)==1) then !Down
                numcount(8)=numcount(8)+1
            endif
            if (nodeXyzIndex(3)==nodeXyzIndex(6)) then !Up
                numcount(9)=numcount(9)+1
            endif
        endif
    enddo
    
end subroutine setEquationNumber


subroutine createElement(nelement, ntags, iy, iz, elementCenterCoor)
    use globalvar
    implicit none
    integer (kind = 4) :: nelement, ntags, iy, iz, i, j 
    real (kind = dp) :: elementCenterCoor(3)
    
    elementCenterCoor = 0.0d0 
    
    nelement        = nelement + 1
    et(nelement)    = 1 
    ien(1,nelement) = plane1(iy-1,iz-1)
    ien(2,nelement) = plane2(iy-1,iz-1)
    ien(3,nelement) = plane2(iy,iz-1)
    ien(4,nelement) = plane1(iy,iz-1)
    ien(5,nelement) = plane1(iy-1,iz)
    ien(6,nelement) = plane2(iy-1,iz)
    ien(7,nelement) = plane2(iy,iz)
    ien(8,nelement) = plane1(iy,iz)
    ids(nelement)   = ntags
    
    do i=1,nen
        do j=1,3
            elementCenterCoor(j)=elementCenterCoor(j)+x(j,ien(i,nelement))
        enddo
    enddo
    elementCenterCoor = elementCenterCoor/8.0d0
    
    if (elementCenterCoor(1)>PMLb(1).or.elementCenterCoor(1)<PMLb(2) &
    .or.elementCenterCoor(2)>PMLb(3).or.elementCenterCoor(2)<PMLb(4) &
    .or.elementCenterCoor(3)<PMLb(5)) then
        et(nelement) = 2
        ntags        = ntags+15+6
    else
        ntags        = ntags+12
    endif
    
    ! assign ien(nelement), ids(nelement), et(nelement)
    ! return elementCenterCoor, 
    ! update nelement, ntags
end subroutine createElement

 subroutine replaceSlaveWithMasterNode(nodeCoor, nelement, nftnd0)
    use globalvar 
    implicit none
    integer (kind = 4) :: nelement, iFault, iFaultNodePair, nftnd0(ntotft), k
    real (kind = dp) :: nodeCoor(10)
    ! The default grids only contain slave nodes. 
    ! This subroutine will replace slave nodes with corresponding master nodes.
    
    ! Currently only work for vertical fault on xz plane. 
    if(et(nelement) == 1 .and. &
        (nodeCoor(1)>(fltxyz(1,1,1)-tol) .and. nodeCoor(1)<(fltxyz(2,1,1)+dx+tol) .and. &
         nodeCoor(3)>(fltxyz(1,3,1)-tol) .and. nodeCoor(2)>0.0d0 .and. abs(nodeCoor(2)-dy)<tol)) then
        do iFault = 1, ntotft
            do iFaultNodePair = 1, nftnd0(iFault)
                do k = 1,nen
                    if(ien(k,nelement)==nsmp(1,iFaultNodePair,iFault)) then
                        ien(k,nelement) = nsmp(2,iFaultNodePair,iFault)  !use master node for the node!
                    endif
                enddo
            enddo
        enddo
    endif      
    
end subroutine replaceSlaveWithMasterNode

subroutine checkIsOnFault(nodeCoor, iFault, isOnFault)
    use globalvar
    implicit none
    integer (kind = 4) :: isOnFault, iFault
    real (kind = dp) :: nodeCoor(10)
    isOnFault = 0
    if(nodeCoor(1)>=(fltxyz(1,1,iFault)-tol).and.nodeCoor(1)<=(fltxyz(2,1,iFault)+tol).and. &
        nodeCoor(2)>=(fltxyz(1,2,iFault)-tol).and.nodeCoor(2)<=(fltxyz(2,2,iFault)+tol).and. &
        nodeCoor(3)>=(fltxyz(1,3,iFault)-tol) .and. nodeCoor(3)<=(fltxyz(2,3,iFault)+tol)) then
        isOnFault = 1
    endif 
end subroutine checkIsOnFault

subroutine createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nnode, msnode, nftnd0, neq0, ntag,&
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
use globalvar 
implicit none
integer (kind = 4) :: iFault, iFaultNodePair, isOnFault, nnode, msnode, nftnd0(ntotft), neq0, i, nxuni, nzuni, ntag
integer (kind = 4) :: fltrc(2,nxuni,nzuni,ntotft), ixfi(ntotft), izfi(ntotft), ifs(ntotft), ifd(ntotft), nodeXyzIndex(10)
real (kind = dp) :: nodeCoor(10), ycoort, pfx, pfz

do iFault = 1, ntotft

    isOnFault = 0 
    call checkIsOnFault(nodeCoor, iFault, isOnFault)

    if (isOnFault == 1) then
        nftnd0(iFault)                = nftnd0(iFault) + 1 ! # of split-node pairs + 1
        nsmp(1,nftnd0(iFault),iFault) = nnode              ! set Slave node nodeID to nsmp  
        msnode                        = nodeXyzIndex(4)*nodeXyzIndex(5)*nodeXyzIndex(6) + nftnd0(iFault) ! create Master node at the end of regular grids
        if (iFault>1) stop 'msnode cannot handle iFault>1'
        
        locid(msnode)                 = ntag
        dof1(msnode)                  = 3         
        nsmp(2,nftnd0(iFault),iFault) = msnode !set Master node nodeID to nsmp
        plane2(nodeXyzIndex(5)+iFault,nodeXyzIndex(3)) = msnode
        
        x(1,msnode) = nodeCoor(1)
        x(2,msnode) = nodeCoor(2)
        x(3,msnode) = nodeCoor(3)
        if (insertFaultType > 0) then 
            x(2,msnode) = ycoort
        endif
        
        !set Equation Numbers for the newly created Master Node.
        do i = 1, ndof
            neq0      = neq0 + 1
            ntag      = ntag + 1
            id1(ntag) = neq0
        enddo
        
        ! Count split-node pair # for MPI
        ! fltgm, fltnum are global vars.
        if(nodeXyzIndex(1) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 1
            fltnum(1) = fltnum(1) + 1
        endif
        if(nodeXyzIndex(1) == nodeXyzIndex(4)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 2
            fltnum(2) = fltnum(2) + 1
        endif
        if(nodeXyzIndex(2) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 10
            fltnum(3) = fltnum(3) + 1
        endif
        if(nodeXyzIndex(2) == nodeXyzIndex(5)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 20
            fltnum(4) = fltnum(4) + 1
        endif
        if(nodeXyzIndex(3) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 100
            fltnum(5) = fltnum(5) + 1
        endif
        if(nodeXyzIndex(3) == nodeXyzIndex(6)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 200
            fltnum(6) = fltnum(6) + 1
        endif             

        ! setOnFaultStation
        do i = 1, nonfs(iFault)
            if(abs(nodeCoor(1)-xonfs(1,i,iFault))<tol .and. &
                abs(nodeCoor(3)-xonfs(2,i,iFault))<tol) then
                n4onf = n4onf + 1
                anonfs(1,n4onf) = nftnd0(iFault)
                anonfs(2,n4onf) = i
                anonfs(3,n4onf) = iFault
                exit
            endif
        enddo  

        ! set unit vectors to split-node pair    
        un(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(2,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(3,nftnd0(iFault),iFault) = dcos(fltxyz(2,4,iFault))        
        us(1,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))
        us(2,nftnd0(iFault),iFault) = -dcos(fltxyz(1,4,iFault))
        us(3,nftnd0(iFault),iFault) = 0.0d0
        ud(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(2,nftnd0(iFault),iFault) = dsin(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(3,nftnd0(iFault),iFault) = dsin(fltxyz(2,4,iFault))
        
        if (insertFaultType >0) then
            un(1,nftnd0(iFault),iFault) = -pfx/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(2,nftnd0(iFault),iFault) = 1.0d0/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(3,nftnd0(iFault),iFault) = -pfz/(pfx**2 + 1.0d0 + pfz**2)**0.5    
            us(1,nftnd0(iFault),iFault) = 1.0d0/(1.0d0 + pfx**2)**0.5
            us(2,nftnd0(iFault),iFault) = pfx/(1.0d0 + pfx**2)**0.5
            us(3,nftnd0(iFault),iFault) = 0.0d0
            ud(1,nftnd0(iFault),iFault) = us(2,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault) &
                - us(3,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault)
            ud(2,nftnd0(iFault),iFault) = us(3,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault) &
                - us(1,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault)
            ud(3,nftnd0(iFault),iFault) = us(1,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault) &
                - us(2,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault)
        endif                             
        
        !...prepare for area calculation
        if(ixfi(iFault)==0) ixfi(iFault)=nodeXyzIndex(1)
        if(izfi(iFault)==0) izfi(iFault)=nodeXyzIndex(3)
        ifs(iFault)=nodeXyzIndex(1)-ixfi(iFault)+1
        ifd(iFault)=nodeXyzIndex(3)-izfi(iFault)+1
        fltrc(1,ifs(iFault),ifd(iFault),iFault) = msnode    !master node
        fltrc(2,ifs(iFault),ifd(iFault),iFault) = nftnd0(iFault) !fault node num in sequence
    endif 
enddo 

end subroutine createMasterNode

subroutine createNode(nodeCoor, xcoor, ycoor, zcoor, nnode, nodeXyzIndex)
    use globalvar
    implicit none
    integer (kind = 4) :: nnode, nodeXyzIndex(10), iy, iz
    real (kind = dp) :: nodeCoor(10), xcoor, ycoor, zcoor
    iy = nodeXyzIndex(2)
    iz = nodeXyzIndex(3)
    nodeCoor(1)   = xcoor
    nodeCoor(2)   = ycoor
    nodeCoor(3)   = zcoor
    nnode         = nnode + 1        
    plane2(iy,iz) = nnode
    x(1,nnode)    = nodeCoor(1)
    x(2,nnode)    = nodeCoor(2)
    x(3,nnode)    = nodeCoor(3) 
end subroutine createNode

subroutine insertFaultInterface(nodeCoor, ycoort, pfx, pfz)
    ! This subroutine is to modify ycoor if a rough_fault interface is inserted.
    
    use globalvar
    implicit none
    real (kind = dp) :: nodeCoor(10), peak, ycoort, pfx, pfz
    real (kind = dp) :: fx1, fx2, fz1
    integer (kind = 4) :: ixx, izz
    
    fx1 = rough_fx_min
    fx2 = rough_fx_max
    fz1 = rough_fz_min
    if ((nodeCoor(1) < fx2 + tol) .and. (nodeCoor(1) > fx1 - tol) .and. (nodeCoor(3) > fz1 - tol)) then 
        ixx = (nodeCoor(1) - fx1)/dx + 1
        izz = (nodeCoor(3) - fz1)/dx + 1
    elseif ((nodeCoor(1) < fx1 - tol) .and. (nodeCoor(3) > fz1 - tol) ) then
        ixx = 1
        izz = (nodeCoor(3) - fz1)/dx + 1
    elseif ((nodeCoor(1) > fx2 + tol) .and. (nodeCoor(3) > fz1 - tol)) then 
        ixx = nnx
        izz = (nodeCoor(3) - fz1)/dx + 1
    elseif ((nodeCoor(1) < fx2 + tol) .and. (nodeCoor(1) > fx1 - tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = (nodeCoor(1) - fx1)/dx + 1
        izz = 1
    elseif ((nodeCoor(1) < fx1 - tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = 1
        izz = 1 
    elseif ((nodeCoor(1) > fx2 + tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = nnx
        izz = 1
    endif 
    
    peak = rough_geo(1,nnz*(ixx-1)+izz)
    pfx = rough_geo(2,nnz*(ixx-1)+izz)
    pfz = rough_geo(3,nnz*(ixx-1)+izz)    
    
    if (nodeCoor(2) > -tol) then
        ycoort = nodeCoor(2)*(ymax - peak)/ymax + peak
    elseif (nodeCoor(2) < -tol) then 
        ycoort = nodeCoor(2)*(peak - ymin)/(-ymin) + peak 
    endif 
    
end subroutine insertFaultInterface

subroutine initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
    use globalvar 
    implicit none 
    integer (kind = 4) :: ix, iy, iz, nx, ny, nz, nodeXyzIndex(10)
    nodeXyzIndex(1) = ix
    nodeXyzIndex(2) = iy
    nodeXyzIndex(3) = iz
    nodeXyzIndex(4) = nx
    nodeXyzIndex(5) = ny
    nodeXyzIndex(6) = nz
end subroutine initializeNodeXyzIndex