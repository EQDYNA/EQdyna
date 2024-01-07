SUBROUTINE qdct2
use globalvar
implicit none
include 'mpif.h'
!...program to calcuate lumped mass matrix for 8-node hexahedral
!    (brick),elastic continue element and assemble.
logical :: lcubic
integer(kind=4) :: nel,i,j,k,k1,m,itmp,itmp1,itmp2
real(kind = dp) :: det,constm,vol,ce,co
real(kind = dp),dimension(nee) :: eleffm
real(kind = dp),dimension(nesd,nen) :: xl
real(kind = dp),dimension(3,3) :: xs
real(kind = dp),dimension(nrowb,nee) :: b
real(kind = dp),dimension(nrowsh,nen) :: shg
!...MPI
integer  ierr, rlp, rr, jj
integer istatus(MPI_STATUS_SIZE)
real (kind = dp), allocatable, dimension(:) :: btmp, btmp1, btmp2, btmp3
!DL variables for PML layer
integer (kind=4):: non,itag,eqn
!*******3D MPI partitioning*************
integer(kind=4)::mexyz(3), npxyz(3), numxyz(3), ix,iy,iz,nodenumtemp,ntagMPI,izz
integer(kind=4)::bnd(2),bndf,bndb,bndd,bndu,rrr,jjj

do nel=1,numel
    !...initialize and localize
    eleffm = 0.0d0
    do i=1,nen
        k = ien(i,nel)
        do j=1,nesd
        xl(j,i) = x(j,k)
        enddo
    enddo
   
    !...if degenerate
    lcubic = .true.
    itmp = 0
    outloop: do i=2,nen
                do j=1,i-1
                    if(ien(j,nel)==ien(i,nel)) itmp=itmp+1
                        if(itmp>=2) then
                            lcubic = .false.
                        exit outloop
                    endif
                enddo
            enddo outloop
    !...compute global shape function and derivatives at 1 Gaussian
    !    point.
    call qdcshg(xl, det, shg, nel, xs, lcubic)
    call contm(shg, det, eleffm, mat(nel,3))
    call assembleElementMassDetShg(nel, eleffm, det, shg)

    call calcSSPhi4Hrgls(nel, xl, xs, shg)
enddo  !nel

call MPI4NodalQuant(alhs, 3)
call MPI4NodalQuant(fnms, 1)

end SUBROUTINE qdct2

subroutine MPI4NodalQuant(quantArray, numDof)
    ! handle MPI communication for Nodal quantities - nodal force and nodal mass.
    use globalvar
    implicit none
    include 'mpif.h'
    integer (kind = 4) ::  ierr, istatus(MPI_STATUS_SIZE), i, ixyz, numDof, rrr, &
        ix,iy,iz, nodenumtemp, ntagMPI, dest, sendtag, source, recvtag, ib, iSign, &
        bnd(2), mexyz(3), npxyz(3), numxyz(3), abc(3)
    real (kind = dp) :: quantArray(neq) 
    real (kind = dp), allocatable, dimension(:) :: btmp, btmp1
    
    mexyz(1)=int(me/(npy*npz))
    mexyz(2)=int((me-mexyz(1)*npy*npz)/npz)
    mexyz(3)=int(me-mexyz(1)*npy*npz-mexyz(2)*npz)
    
    npxyz(1) = npx 
    npxyz(2) = npy
    npxyz(3) = npz
    
    do i = 1, 3
        numxyz(i) = numcount(i) ! number of nodes along each direction inside the MPI process.
    enddo
    
    abc(1) = numxyz(2)*numxyz(3)
    abc(2) = numxyz(1)*numxyz(3)
    abc(3) = numxyz(1)*numxyz(2)

    ! loop over MPI interfaces along x, y, z directions
    do ixyz = 1, 3
        if (npxyz(ixyz)>1)then
            bnd(1) = 1            !- boundary
            bnd(2) = numxyz(ixyz) !+ boundary
            if (mexyz(ixyz) == master) then
                bnd(1) = 0
            elseif (mexyz(ixyz) == npxyz(ixyz)-1) then
                bnd(2) = 0
            endif
            !send dest: me-npy*npz: envelop:20000+me 
            !recv source: me-npy*npz: envelop:20000+me-npy*npz 
            !left boundary of me  ---  right boundary of me-npy*npz
            !right boundary of me  ---  left boundary of me+npy*npz
            !front boundary of me  ---  back boundary of me-npz
            !back boundary of me  ---  front boundary of me+npz
            !lower boundary of me  ---  upper boundary of me-1
            !upper boundary of me  ---  lower boundary of me+1
            !me :     x+ ==me+npy*npz
            !        x- ==me-npy*npz
            !        y+ ==me+npz
            !        y- ==me-npz    
            !        z+ ==me+1    
            !        z- ==me-1
            do ib = 1, 2
                if (ib == 1) iSign=1
                if (ib == 2) iSign=-1
                if (bnd(ib) /= 0)then 
                    if (ixyz == 1) then
                        dest    = me - npxyz(2)*npxyz(3)*iSign
                        sendtag = 200000*numDof + me
                        source  = me - npxyz(2)*npxyz(3)*iSign
                        recvtag = 200000*numDof + me-npxyz(2)*npxyz(3)*iSign
                    elseif (ixyz == 2) then 
                        dest    = me - npxyz(3)*iSign
                        sendtag = 210000*numDof + me
                        source  = me - npxyz(3)*iSign
                        recvtag = 210000*numDof + me - npxyz(3)*iSign
                    elseif (ixyz == 3) then 
                        dest    = me - iSign
                        sendtag = 220000*numDof + me
                        source  = me - iSign
                        recvtag = 220000*numDof + me - iSign
                    endif 
                    
                    if (numDof == 1) then 
                        rrr = abc(ixyz) + fltnum(2*(ixyz-1)+ib)
                    elseif (numDof == 3) then
                        rrr = numcount(3+2*(ixyz-1)+ib) + fltnum(2*(ixyz-1)+ib)*3
                    endif 
                    allocate(btmp(rrr),btmp1(rrr))
                    
                    ntagMPI = 0
                    if (ixyz == 1) then 
                        do iz=1,numxyz(3)
                            do iy=1,numxyz(2)
                                nodenumtemp=(bnd(ib)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                                !                         nodeID, numDof, fetch, ...
                                call processNodalQuantArr(nodenumtemp, numDof, 1, btmp, rrr, quantArray, ntagMPI) 
                            enddo
                        enddo
                    elseif (ixyz == 2) then 
                        do ix=1,numxyz(1)
                            do iz=1,numxyz(3)
                                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bnd(ib)
                                call processNodalQuantArr(nodenumtemp, numDof, 1, btmp, rrr, quantArray, ntagMPI)   
                            enddo
                        enddo
                    elseif (ixyz == 3) then 
                        do ix=1,numxyz(1)
                            do iy=1,numxyz(2)
                                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bnd(ib)-1)*numxyz(2)+iy
                                call processNodalQuantArr(nodenumtemp, numDof, 1, btmp, rrr, quantArray, ntagMPI)   
                            enddo
                        enddo
                    endif 
                    
                    ! !Check
                    ! if (numcount(3+2*(ixyz-1)+ib)/=ntagMPI) then 
                        ! stop 'rr&ntagMPI-qdct2-bnd(1)'
                        ! write(*,*) 'rr=',numcount(3+2*(ixyz-1)+ib),'ntagMPI=',ntagMPI
                    ! endif
        !
                    if (fltMPI(2*(ixyz-1)+ib)) then
                        do ix = 1, fltnum(2*(ixyz-1)+ib)
                            if (ixyz == 1 .and. ib==1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                            elseif (ixyz == 1 .and. ib==2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                            elseif (ixyz == 2 .and. ib == 1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltf(ix)
                            elseif (ixyz == 2 .and. ib == 2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltb(ix)
                            elseif (ixyz == 3 .and. ib == 1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltd(ix)
                            elseif (ixyz == 3 .and. ib == 2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltu(ix)
                            endif 
                            call processNodalQuantArr(nodenumtemp, numDof, 1, btmp, rrr, quantArray, ntagMPI)
                          
                        enddo
                    endif
                
                    call mpi_sendrecv(btmp,  rrr, MPI_DOUBLE_PRECISION, dest, sendtag, &
                                      btmp1, rrr, MPI_DOUBLE_PRECISION, source, recvtag, &
                                      MPI_COMM_WORLD, istatus, ierr)
                    
                    ntagMPI=0
                    if (ixyz == 1) then 
                        do iz=1,numxyz(3)
                            do iy=1,numxyz(2)
                                nodenumtemp=(bnd(ib)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                                !                                        add
                                call processNodalQuantArr(nodenumtemp, numDof, 2, btmp1, rrr, quantArray, ntagMPI)  
                            enddo
                        enddo
                    elseif (ixyz == 2) then 
                        do ix=1,numxyz(1)
                            do iz=1,numxyz(3)
                                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bnd(ib)
                                call processNodalQuantArr(nodenumtemp, numDof, 2, btmp1, rrr, quantArray, ntagMPI)   
                            enddo
                        enddo
                    elseif (ixyz == 3) then 
                        do ix=1,numxyz(1)
                            do iy=1,numxyz(2)
                                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bnd(ib)-1)*numxyz(2)+iy
                                call processNodalQuantArr(nodenumtemp, numDof, 2, btmp1, rrr, quantArray, ntagMPI)    
                            enddo
                        enddo
                    endif 
                    
                    if (fltMPI(2*(ixyz-1)+ib)) then
                        do ix=1,fltnum(2*(ixyz-1)+ib)
                            if (ixyz == 1 .and. ib==1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                            elseif (ixyz == 1 .and. ib==2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                            elseif (ixyz == 2 .and. ib == 1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltf(ix)
                            elseif (ixyz == 2 .and. ib == 2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltb(ix)
                            elseif (ixyz == 3 .and. ib == 1) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltd(ix)
                            elseif (ixyz == 3 .and. ib == 2) then 
                                nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltu(ix)
                            endif 
                            call processNodalQuantArr(nodenumtemp, numDof, 2, btmp1, rrr, quantArray, ntagMPI) 
                        enddo
                    endif
                    deallocate(btmp, btmp1)
                endif 
            enddo 
        endif
        call mpi_barrier(MPI_COMM_WORLD, ierr)
    enddo 
end subroutine MPI4NodalQuant

subroutine processNodalQuantArr(nodeID, numDof, operation, resArr, resArrSize, quantArray, ntagMPI)
    use globalvar
    implicit none
    integer (kind = 4) :: nodeID, numDof, iDof, ntagMPI, resArrSize, operation
    real(kind = dp) :: resArr(resArrSize), quantArray(neq)
    !character (len = 20) :: operation
    
    ! operation = 'fetch'; 'add'
    
    if (numDof == 1) then 
        ntagMPI = ntagMPI + 1
        if (operation == 1) then 
            resArr(ntagMPI) = quantArray(nodeID)
        elseif (operation == 2) then 
            quantArray(nodeID) = quantArray(nodeID) + resArr(ntagMPI)
        endif 
    elseif (numDof == 3) then 
        if (operation == 1) then 
            do iDof = 1, dof1(nodeID)
                if (id1(locid(nodeID)+iDof)>0) then
                    ntagMPI = ntagMPI + 1
                    resArr(ntagMPI) = quantArray(id1(locid(nodeID)+iDof))
                endif
            enddo
        elseif (operation == 2) then 
            do iDof = 1, dof1(nodeID)
                if (id1(locid(nodeID)+iDof)>0) then
                    ntagMPI = ntagMPI + 1
                    quantArray(id1(locid(nodeID)+iDof)) = quantArray(id1(locid(nodeID)+iDof)) + &
                        resArr(ntagMPI)
                endif
            enddo   
        endif
    endif 
end subroutine processNodalQuantArr

subroutine assembleElementMassDetShg(elemID, elementMass, det, shg)
    use globalvar 
    implicit none 
    integer (kind = 4) :: i, j, eqn, nodeID, ixyz, elemID
    real (kind = dp) :: elementMass(nee), det, shg(nrowsh, nen)
    
    do i = 1, nen 
        nodeID = ien(i,elemID)
        if (dof1(nodeID)==12) then
            do ixyz = 1, 3
                do j = 3*(ixyz-1)+1, 3*(ixyz-1)+3
                    eqn = id1(locid(nodeID)+j)
                    if (eqn > 0) then
                        alhs(eqn) = alhs(eqn) + elementMass(3*(i-1)+ixyz)
                    endif
                enddo
                eqn = id1(locid(nodeID)+ixyz+9)
                if (eqn>0) then
                    alhs(eqn) = alhs(eqn) + elementMass(3*(i-1)+ixyz)
                endif
            enddo                       
        elseif (dof1(nodeID) == ndof) then
            do j = 1, ndof
                eqn = id1(locid(nodeID)+j)
                alhs(eqn) = alhs(eqn) + elementMass((i-1)*3+j)
            enddo
        endif
        
        fnms(nodeID) = fnms(nodeID) + elementMass((i-1)*ned+1)
    enddo
    
    do i = 1, nee
        elemass(i,elemID) = elementMass(i)
    enddo

    eledet(elemID) = det
    
    do i = 1, nen
        do j=1, nrowsh-1
            eleshp(j,i,elemID) = shg(j,i)
        enddo
    enddo
    
end subroutine assembleElementMassDetShg

subroutine calcSSPhi4Hrgls(elemID, xl, xs, shg)
    use globalvar
    implicit none
    integer (kind = 4) :: elemID, i, j, k
    integer (kind = 4), dimension(8,4) :: ha = reshape((/ &
            1,1,-1,-1,-1,-1,1,1, 1,-1,-1,1,-1,1,1,-1, &
            1,-1,1,-1,1,-1,1,-1, -1,1,-1,1,1,-1,1,-1/), &
            (/8,4/))
    real (kind = dp) :: xl(nesd,nen), shg(nrowsh, nen), xs(3,3), &
        vol, ce, co
    
    ! calc SS matrix for hourglass control 
    call vlm(xl,vol)
    ce = mat(elemID,5)*(3*mat(elemID,4)+2*mat(elemID,5))/(mat(elemID,4)+mat(elemID,5))
    ce = 16.d0* ce / 15.d0    !close to plane-strain
    !E=miu*(3*lam+2*miu)/(lam+miu)
    co = ce * vol / 48.d0
    ss(1,elemID) = co*(xs(1,1)**2+xs(2,1)**2+xs(3,1)**2)
    ss(2,elemID) = co*(xs(1,1)*xs(1,2)+xs(2,1)*xs(2,2)+xs(3,1)*xs(3,2))
    ss(3,elemID) = co*(xs(1,1)*xs(1,3)+xs(2,1)*xs(2,3)+xs(3,1)*xs(3,3))
    ss(4,elemID) = co*(xs(1,2)**2+xs(2,2)**2+xs(3,2)**2)
    ss(5,elemID) = co*(xs(1,2)*xs(1,3)+xs(2,2)*xs(2,3)+xs(3,2)*xs(3,3))
    ss(6,elemID) = co*(xs(1,3)**2+xs(2,3)**2+xs(3,3)**2)
    
    ! calc 4 sets of Phi for hourglass control
    do i = 1, 4    
        ! calc phi prime
        do j = 1, nen
            vol = 0.0d0    !temp variable
            do k = 1, nen
                vol = vol + ha(k,i)*(xl(1,k)*shg(1,j) + &
                    xl(2,k)*shg(2,j) + xl(3,k)*shg(3,j))
            enddo
            phi(j,i,elemID) = ha(j,i) - vol
        enddo
        
        vol = 0.0d0
        do j = 1, nen
            vol = vol + phi(j,i,elemID)**2
        enddo
        vol = sqrt(vol/8.0d0)
        ! normalize to get phi from phi prime
        do j = 1, nen
            phi(j,i,elemID) = phi(j,i,elemID) / vol
        enddo
    enddo
end subroutine calcSSPhi4Hrgls

subroutine contm(shg,det,elmass,constm)
    use globalvar
    implicit none     
    ! to calc lumped mass for an element
    integer (kind=4) :: j,n,k
    real (kind = dp) :: det,dsum,totmas,constm,temp1,temp2
    real (kind = dp),dimension(nee) :: elmass
    real (kind = dp),dimension(nen) :: work
    real (kind = dp),dimension(nrowsh,nen) :: shg
    !...initialize
    dsum   = 0.0d0
    totmas = 0.0d0
    work   = 0.0d0
    !...calculate
    totmas = constm*w*det
    do j=1,nen
        temp2 = totmas*shg(nrowsh,j)**2
        dsum = dsum + temp2
        work(j) = work(j) + temp2
    enddo
    !...scale diagonal to conserve total mass
    temp1 = totmas/dsum
    !...store terms in a column
    do j=1,nen
        temp2 = temp1*work(j)
        n = (j - 1)*ned
        do k=1,ned
            elmass(n + k) = temp2
        enddo
    enddo
end subroutine contm