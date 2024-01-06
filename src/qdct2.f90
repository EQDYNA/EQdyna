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
integer(kind=4),dimension(8,4) :: ha = reshape((/ &
    1,1,-1,-1,-1,-1,1,1, 1,-1,-1,1,-1,1,1,-1, &
    1,-1,1,-1,1,-1,1,-1, -1,1,-1,1,1,-1,1,-1/), &
    (/8,4/))
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
    m = 1    !material number
    !...if element degenerate. B.D. 11/27/08
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
    call qdcshg(xl,det,shg,nel,xs,lcubic)
    !...form mass matrix for left-hand-side
    constm = (1.0d0+rdampm*0.5d0*dt)*mat(nel,3)
    if(constm /= 0.0d0) then
        call contm(shg,det,eleffm,constm)
    endif
    !...assemble into left-hand-side
    
    ! *.* D.L. Jan/23/2015. Including rayleigh damping into the left-hand-side when using SUBROUTINE contm
    do i=1,nen
            ! do j=1,ned
                ! k=lm(j,i,nel)
                ! k1=j+(i-1)*ned
                ! if(k>0) then
                    ! alhs(k) = alhs(k) + eleffm(k1)
                ! endif
            ! enddo
        ! *.* Update alhs according to dof1. D.L. Jan/23/2015    
        non=ien(i,nel)
        if (dof1(non)==12) then
            do j=1,3
                itag=locid(non)+j
                eqn=id1(itag)
                if (eqn>0) then
                    alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+1)
                endif
            enddo
            do j=4,6
                itag=locid(non)+j
                eqn=id1(itag)
                if (eqn>0) then
                    alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+2)
                endif
            enddo    
            do j=7,9
                itag=locid(non)+j
                eqn=id1(itag)
                if (eqn>0) then
                    alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+3)
                endif
            enddo    
            j=10
            itag=locid(non)+j
            eqn=id1(itag)
            if (eqn>0) then
                alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+1)
            endif
            j=11
            itag=locid(non)+j
            eqn=id1(itag)
            if (eqn>0) then
                alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+2)
            endif    
            j=12
            itag=locid(non)+j
            eqn=id1(itag)
            if (eqn>0) then
                alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+3)
            endif                        
        elseif (dof1(non)==3) then
            do j=1,3
                itag=locid(non)+j
                eqn=id1(itag)
                alhs(eqn)=alhs(eqn)+eleffm((i-1)*3+j)
            enddo
        endif
        ! *.* D.L.
    enddo

    eledet(nel) = det
    !...derivatives of global shape function
    do i=1,nen
        do j=1,nrowsh-1
            eleshp(j,i,nel) = shg(j,i)
        enddo
    enddo
    !...mass matrix for right-hand-side to use.
    !    Note: different "constm" from above.
    constm = mat(nel,3)
    eleffm = 0.0d0    !must initialize again here
    call contm(shg,det,eleffm,constm)
    do i=1,nee
        elemass(i,nel) = eleffm(i)
    enddo
    !...assemble above mass to nodal mass for "faulting".
    !    Note: ned degree of an element node have a same
    !    mass, so only one used here.
    do i=1,nen
        k = ien(i,nel)
        k1= 1 + (i-1) * ned
        fnms(k) = fnms(k) + eleffm(k1)
    enddo
    !...SS matrix for "hrglss"
    call vlm(xl,vol)
    ce=mat(nel,5)*(3*mat(nel,4)+2*mat(nel,5))/(mat(nel,4)+mat(nel,5))
    ce = 16.d0* ce / 15.d0    !close to plane-strain
    !E=miu*(3*lam+2*miu)/(lam+miu)
    co = ce * vol / 48.d0
    ss(1,nel)=co*(xs(1,1)**2+xs(2,1)**2+xs(3,1)**2)
    ss(2,nel)=co*(xs(1,1)*xs(1,2)+xs(2,1)*xs(2,2)+xs(3,1)*xs(3,2))
    ss(3,nel)=co*(xs(1,1)*xs(1,3)+xs(2,1)*xs(2,3)+xs(3,1)*xs(3,3))
    ss(4,nel)=co*(xs(1,2)**2+xs(2,2)**2+xs(3,2)**2)
    ss(5,nel)=co*(xs(1,2)*xs(1,3)+xs(2,2)*xs(2,3)+xs(3,2)*xs(3,3))
    ss(6,nel)=co*(xs(1,3)**2+xs(2,3)**2+xs(3,3)**2)
    !...Phi for "hrglss"
    ! Note: there are 4 sets of Phi for one element in 3D.
    do i=1,4    !4 sets of phi
        !...phi prime
        do j=1,nen
            vol = 0.0d0    !temp variable
            do k=1,nen
                vol = vol + ha(k,i)*(xl(1,k)*shg(1,j) + xl(2,k)  &
                *shg(2,j) + xl(3,k)*shg(3,j))
            enddo
            phi(j,i,nel) = ha(j,i) - vol
        enddo
        !...sum
        vol = 0.0d0
        do j=1,nen
            vol = vol + phi(j,i,nel)**2
        enddo
        vol = sqrt(vol/8.0d0)
        !...normalize to get phi from phi prime
        do j=1,nen
            phi(j,i,nel) = phi(j,i,nel) / vol
        enddo
        !
    enddo
    !
enddo  !nel

call MPI4NodalQuant(alhs)

!prepare for MPI partitioning
mexyz(1)=int(me/(npy*npz))
mexyz(2)=int((me-mexyz(1)*npy*npz)/npz)
mexyz(3)=int(me-mexyz(1)*npy*npz-mexyz(2)*npz)

npxyz(1) = npx 
npxyz(2) = npy
npxyz(3) = npz

numxyz(1)=numcount(1)
numxyz(2)=numcount(2)
numxyz(3)=numcount(3)
    
!===============================3DMPI===============================!
!=====================Partitioning along x axis=====================!
if (npxyz(1)>1) then
    rr=numxyz(2)*numxyz(3)+fltnum(1)
    jj=numxyz(2)*numxyz(3)+fltnum(2)
       bnd(1)=1  !-x boundary
       bnd(2)=numxyz(1) !+x boundary 
       if (mexyz(1) == master) then
         bnd(1)=0
       elseif (mexyz(1) == npxyz(1)-1) then
         bnd(2)=0
       endif    
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

    if (bnd(1)/=0) then
        ntagMPI=0
        do iz=1,numxyz(3)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(bnd(1)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(1)) then
            do ix=1,fltnum(1)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-npxyz(2)*npxyz(3), 300000+me, &
                btmp1, rr, MPI_DOUBLE_PRECISION, me-npxyz(2)*npxyz(3), 300000+me-npxyz(2)*npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do iz=1,numxyz(3)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(bnd(1)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)                
            enddo
        enddo
        if (fltMPI(1)) then
            do ix=1,fltnum(1)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
            enddo
        endif
    endif !bnd(1)/=0
    !write(*,*) 'me',me,'Finished BNDL-fnms'
    if (bnd(2)/=0)then
        ntagMPI=0
        do iz=1,numxyz(3)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(bnd(2)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(2)) then
            do ix=1,fltnum(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+npxyz(2)*npxyz(3), 300000+me, &
                btmp3, jj, MPI_DOUBLE_PRECISION, me+npxyz(2)*npxyz(3), 300000+me+npxyz(2)*npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do iz=1,numxyz(3)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(bnd(2)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)                
            enddo
        enddo
        if (fltMPI(2)) then
            do ix=1,fltnum(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
            enddo
        endif        
    endif !bnd(2)/=0
    !write(*,*) 'me',me,'Finished BNDR-fnms'    
    deallocate(btmp, btmp1, btmp2, btmp3)
endif !npxyz(1)>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along y axis=====================!
if (npxyz(2)>1) then
    rr=numxyz(1)*numxyz(3)+fltnum(3)
    jj=numxyz(1)*numxyz(3)+fltnum(4)
       bndf=1  !-y boundary
       bndb=numxyz(2) !+y boundary 
       if (mexyz(2) == master) then
         bndf=0
       elseif (mexyz(2) == npxyz(2)-1) then
         bndb=0
       endif    
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

    if (bndf/=0) then
        ntagMPI=0
        do ix=1,numxyz(1)
            do iz=1,numxyz(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndf
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(3)) then
            do iy=1,fltnum(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltf(iy)
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-npxyz(3), 310000+me, &
                btmp1, rr, MPI_DOUBLE_PRECISION, me-npxyz(3), 310000+me-npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do ix=1,numxyz(1)
            do iz=1,numxyz(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndf
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)                
            enddo
        enddo
        if (fltMPI(3)) then
            do iy=1,fltnum(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltf(iy)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
            enddo
        endif
    endif !bndf/=0
    !write(*,*) 'me',me,'Finished BNDF-fnms'
    if (bndb/=0)then
        ntagMPI=0
        do ix=1,numxyz(1)
            do iz=1,numxyz(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndb
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(4)) then
            do iy=1,fltnum(4)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltb(iy)
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+npxyz(3), 310000+me, &
                btmp3, jj, MPI_DOUBLE_PRECISION, me+npxyz(3), 310000+me+npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do ix=1,numxyz(1)
            do iz=1,numxyz(3)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndb
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)                
            enddo
        enddo
        if (fltMPI(4)) then
            do iy=1,fltnum(4)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltb(iy)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
            enddo
        endif        
    endif !bndb/=0
    !write(*,*) 'me',me,'Finished BNDB-fnms'
    deallocate(btmp, btmp1, btmp2, btmp3)
endif !npxyz(2)>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along z axis=====================!
if (npxyz(3)>1) then
    rr=numxyz(1)*numxyz(2)+fltnum(5)
    jj=numxyz(1)*numxyz(2)+fltnum(6)
       bndd=1  !-z boundary
       bndu=numxyz(3) !+z boundary 
       if (mexyz(3) == master) then
         bndd=0
       elseif (mexyz(3) == npxyz(3)-1) then
         bndu=0
       endif    
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

    if (bndd/=0) then
        ntagMPI=0
        do ix=1,numxyz(1)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndd-1)*numxyz(2)+iy
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(5)) then
            do iz=1,fltnum(5)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltd(iz)
                btmp(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-1, 320000+me, &
                btmp1, rr, MPI_DOUBLE_PRECISION, me-1, 320000+me-1,&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do ix=1,numxyz(1)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndd-1)*numxyz(2)+iy
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)                
            enddo
        enddo
        if (fltMPI(5)) then
            do iz=1,fltnum(5)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltd(iz)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
            enddo
        endif
    endif !bndd/=0
    !write(*,*) 'me',me,'Finished BNDD-fnms'
    if (bndu/=0)then
        ntagMPI=0
        do ix=1,numxyz(1)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndu-1)*numxyz(2)+iy
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        enddo
        if (fltMPI(6)) then
            do iz=1,fltnum(6)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltu(iz)
                btmp2(ntagMPI)=fnms(nodenumtemp)
            enddo
        endif
         call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+1, 320000+me, &
                btmp3, jj, MPI_DOUBLE_PRECISION, me+1, 320000+me+1,&
                MPI_COMM_WORLD, istatus, ierr)
        ntagMPI=0
        do ix=1,numxyz(1)
            do iy=1,numxyz(2)
                ntagMPI=ntagMPI+1
                nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndu-1)*numxyz(2)+iy
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)                
            enddo
        enddo
        if (fltMPI(6)) then
            do iz=1,fltnum(6)
                ntagMPI=ntagMPI+1
                nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltu(iz)
                fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
            enddo
        endif        
    endif !bndu/=0
    !write(*,*) 'me',me,'Finished BNDU-fnms'    
    deallocate(btmp, btmp1, btmp2, btmp3)
endif !npxyz(3)>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
end SUBROUTINE qdct2


subroutine MPI4NodalQuant(quantArray)
    ! handle MPI communication for Nodal quantities - nodal force and nodal mass.
    use globalvar
    implicit none
    include 'mpif.h'
    integer (kind = 4) ::  ierr, rlp, rr, jj, istatus(MPI_STATUS_SIZE), i
    real (kind = dp) :: quantArray(neq) 
    real (kind = dp), allocatable, dimension(:) :: btmp, btmp1, btmp2, btmp3
    integer (kind = 4)::ix,iy,iz,nodenumtemp,ntagMPI,izz, dest, sendtag, source, recvtag, ib, iSign
    integer (kind = 4)::bnd(2),bndf,bndb,bndd,bndu,rrr,jjj, mexyz(3), npxyz(3), numxyz(3)
    
    !prepare for MPI partitioning
    mexyz(1)=int(me/(npy*npz))
    mexyz(2)=int((me-mexyz(1)*npy*npz)/npz)
    mexyz(3)=int(me-mexyz(1)*npy*npz-mexyz(2)*npz)
    
    npxyz(1) = npx 
    npxyz(2) = npy
    npxyz(3) = npz
    
    do i = 1, 3
        numxyz(i) = numcount(i)
    enddo

    !write(*,*) 'me',me,'begin MPI partioning'
    !*************************************************************
    !************Partioning along x axis**************************
    !Loop sequence x->z->y
    if (npxyz(1)>1)then
        rr     = numcount(3+2*1-1)
        jj     = numcount(3+2*1)
        bnd(1) = 1!-x boundary
        bnd(2) = numxyz(1)!+x boundary
        if (mexyz(1) == master) then
            bnd(1) = 0
        elseif (mexyz(1) == npxyz(1)-1) then
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
                !ib = 1
                dest    = me - npxyz(2)*npxyz(3)*iSign
                sendtag = 200000 + me
                source  = me - npxyz(2)*npxyz(3)*iSign
                recvtag = 200000 + me-npxyz(2)*npxyz(3)*iSign
                
                !Check; if rr==0, bnd(1) should ==0.
                !if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
                rrr = numcount(3+2*(1-1)+ib) + fltnum(2*(1-1)+ib)*3
                allocate(btmp(rrr),btmp1(rrr))
                ntagMPI = 0
                do iz=1,numxyz(3)
                    do iy=1,numxyz(2)
                        nodenumtemp=(bnd(ib)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                        do izz=1,dof1(nodenumtemp)
                            if (id1(locid(nodenumtemp)+izz)>0) then
                                ntagMPI=ntagMPI+1
                                btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                            endif
                        enddo    
                    enddo
                enddo
                !Check
                if (numcount(3+2*(1-1)+ib)/=ntagMPI) then 
                    stop 'rr&ntagMPI-qdct2-bnd(1)'
                    write(*,*) 'rr=',numcount(3+2*(1-1)+ib),'ntagMPI=',ntagMPI
                endif
    !
                if (fltMPI(2*(1-1)+ib)) then
                    do ix = 1, fltnum(2*(1-1)+ib)
                        if (ib==1) then 
                            nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                        elseif (ib==2) then 
                            nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                        endif 
                        do izz = 1, 3
                            ntagMPI = ntagMPI + 1
                            btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        enddo
                    enddo
                endif
            
                call mpi_sendrecv(btmp,  rrr, MPI_DOUBLE_PRECISION, dest, sendtag, &
                                  btmp1, rrr, MPI_DOUBLE_PRECISION, source, recvtag, &
                                  MPI_COMM_WORLD, istatus, ierr)
                
                ntagMPI=0
                do iz=1,numxyz(3)
                    do iy=1,numxyz(2)
                        nodenumtemp=(bnd(ib)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                        do izz=1,dof1(nodenumtemp)
                            if (id1(locid(nodenumtemp)+izz)>0) then
                                ntagMPI=ntagMPI+1
                                quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                    btmp1(ntagMPI)
                            endif
                        enddo    
                    enddo
                enddo
                  if (fltMPI(2*(1-1)+ib)) then
                      do ix=1,fltnum(2*(1-1)+ib)
                        if (ib==1) then 
                            nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltl(ix)
                        elseif (ib==2) then 
                            nodenumtemp = numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                        endif 
                        do izz=1,3
                            ntagMPI = ntagMPI + 1
                            quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
                        enddo
                    enddo
                  endif
                deallocate(btmp, btmp1)
            endif 
        enddo 
        !write(*,*) 'me',me,'Finished BNDL-quantArray'        
        !
        ! if (bnd(2)/=0)then 
            ! ib = 1
            ! dest    = me + npxyz(2)*npxyz(3)
            ! sendtag = 200000 + me
            ! source  = me + npxyz(2)*npxyz(3)
            ! recvtag = 200000 + me + npxyz(2)*npxyz(3)
            ! !Check
            ! if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'            
            ! jjj=jj+fltnum(2*1)*3
                
            ! allocate(btmp2(jjj),btmp3(jjj))
            ! ntagMPI=0
            ! do iz=1,numxyz(3)
                ! do iy=1,numxyz(2)
                    ! nodenumtemp=(bnd(2)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                    ! do izz=1,dof1(nodenumtemp)
                        ! if (id1(locid(nodenumtemp)+izz)>0) then
                            ! ntagMPI=ntagMPI+1
                            ! btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        ! endif
                    ! enddo    
                ! enddo
            ! enddo
            ! !Check
            ! if (jj/=ntagMPI) then 
                ! stop 'jj&ntagMPI-qdct2-bnd(2)'
                ! write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
            ! endif
! !
              ! if (fltMPI(2*1)) then
                  ! do ix=1,fltnum(2*1)
                    ! nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                    ! do izz=1,3
                    ! ntagMPI=ntagMPI+1
                    ! btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                      ! enddo
                ! enddo
              ! endif
             ! call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, dest, sendtag, &
                               ! btmp3, jjj, MPI_DOUBLE_PRECISION, source, recvtag,&
                               ! MPI_COMM_WORLD, istatus, ierr)
! !Update
            ! ntagMPI=0
            ! do iz=1,numxyz(3)
                ! do iy=1,numxyz(2)
                    ! nodenumtemp=(bnd(2)-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+iy
                    ! do izz=1,dof1(nodenumtemp)
                        ! if (id1(locid(nodenumtemp)+izz)>0) then
                            ! ntagMPI=ntagMPI+1
                            ! quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                ! btmp3(ntagMPI)
                        ! endif
                    ! enddo    
                ! enddo
            ! enddo
              ! if (fltMPI(2*1)) then
                  ! do ix=1,fltnum(2*1)
                    ! nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltr(ix)
                    ! do izz=1,3
                    ! ntagMPI=ntagMPI+1
                    ! quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
                      ! enddo
                ! enddo
              ! endif
            ! deallocate(btmp2, btmp3)
        ! endif !bnd(2)/=0
    endif!if npxyz(1)>1
    call mpi_barrier(MPI_COMM_WORLD, ierr)
!****************************************************************************
!******************Paritioning along y axis**********************************
    if (npxyz(2)>1)then
        rr=numcount(3+3)
        jj=numcount(3+4)
        bndf=1!-y boundary
        bndb=numxyz(2)!+y boundary
        if (mexyz(2)==master) then
            bndf=0
        elseif (mexyz(2)==npxyz(2)-1) then
            bndb=0
        endif
        
        if (bndf/=0)then 
            !Check
            if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
            rrr=rr+fltnum(3)*3
            allocate(btmp(rrr),btmp1(rrr))
            ntagMPI=0
            do ix=1,numxyz(1)
                do iz=1,numxyz(3)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndf
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then 
                            ntagMPI=ntagMPI+1
                            btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        endif
                    enddo    
                enddo
            enddo
            !Check
            if (rr/=ntagMPI) then 
                stop 'rr&ntagMPI-qdct2-bndf'
                write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
            endif
!
              if (fltMPI(3)) then
                  do iy=1,fltnum(3)!using the free index, here iy.
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltf(iy)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                      enddo
                enddo
              endif
             call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-npxyz(3), 210000+me, &
                btmp1, rrr, MPI_DOUBLE_PRECISION, me-npxyz(3), 210000+me-npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
!Update
            ntagMPI=0
            do ix=1,numxyz(1)
                do iz=1,numxyz(3)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndf
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then 
                            ntagMPI=ntagMPI+1
                            quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                btmp1(ntagMPI)
                        endif
                    enddo    
                enddo
            enddo
              if (fltMPI(3)) then
                  do iy=1,fltnum(3)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltf(iy)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
                      enddo
                enddo
              endif
            deallocate(btmp, btmp1)            
        endif!bndf/=0
        !write(*,*) 'me',me,'Finished BNDF-quantArray'    
        !
        if (bndb/=0)then 
            !Check
            if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'
            jjj=jj+fltnum(4)*3
            allocate(btmp2(jjj),btmp3(jjj))
            ntagMPI=0
            do ix=1,numxyz(1)
                do iz=1,numxyz(3)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndb
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then
                            ntagMPI=ntagMPI+1
                            btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        endif
                    enddo    
                enddo
            enddo
            !Check
            if (jj/=ntagMPI) then 
                stop 'jj&ntagMPI-qdct2-bnd(1)'
                write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
            endif
!
              if (fltMPI(4)) then
                  do iy=1,fltnum(4)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltb(iy)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                      enddo
                enddo
              endif
             call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+npxyz(3), 210000+me, &
                btmp3, jjj, MPI_DOUBLE_PRECISION, me+npxyz(3), 210000+me+npxyz(3),&
                MPI_COMM_WORLD, istatus, ierr)
!Update
            ntagMPI=0
            do ix=1,numxyz(1)
                do iz=1,numxyz(3)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(iz-1)*numxyz(2)+bndb
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then
                            ntagMPI=ntagMPI+1
                            quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                btmp3(ntagMPI)
                        endif
                    enddo    
                enddo
            enddo
              if (fltMPI(4)) then
                  do iy=1,fltnum(4)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltb(iy)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
                      enddo
                enddo
              endif
            deallocate(btmp2, btmp3)            
        endif!bndb/=0
        !write(*,*) 'me',me,'Finished BNDB-quantArray'
        !
      endif!if npxyz(2)>1
    call mpi_barrier(MPI_COMM_WORLD, ierr)
!****************************************************************************
!******************Paritioning along z axis**********************************
    if (npxyz(3)>1)then
        rr=numcount(3+5)
        jj=numcount(3+6)
        bndd=1!-z boundary
        bndu=numxyz(3)!+z boundary
        if (mexyz(3)==master) then
            bndd=0
        elseif (mexyz(3)==npxyz(3)-1) then
            bndu=0
        endif
        
        if (bndd/=0)then 
            !Check
            if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
            rrr=rr+fltnum(5)*3
            allocate(btmp(rrr),btmp1(rrr))
            ntagMPI=0
            do ix=1,numxyz(1)
                do iy=1,numxyz(2)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndd-1)*numxyz(2)+iy
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then
                            ntagMPI=ntagMPI+1
                            btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        endif
                    enddo    
                enddo
            enddo
            !Check
            if (rr/=ntagMPI) then 
                stop 'rr&ntagMPI-qdct2-bndb-quantArray'
                write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
            endif
!
              if (fltMPI(5)) then
                  do iz=1,fltnum(5)!using the free index, here iz.
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltd(iz)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    btmp(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                      enddo
                enddo
              endif
             call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-1, 220000+me, &
                btmp1, rrr, MPI_DOUBLE_PRECISION, me-1, 220000+me-1,&
                MPI_COMM_WORLD, istatus, ierr)
!Update
            ntagMPI=0
            do ix=1,numxyz(1)
                do iy=1,numxyz(2)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndd-1)*numxyz(2)+iy
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then 
                            ntagMPI=ntagMPI+1
                            quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                btmp1(ntagMPI)
                        endif
                    enddo    
                enddo
            enddo
              if (fltMPI(5)) then
                  do iz=1,fltnum(5)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltd(iz)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
                      enddo
                enddo
              endif
            deallocate(btmp, btmp1)            
        endif!bndd/=0
        !write(*,*) 'me',me,'Finished BNDD-quantArray'
        !
        if (bndu/=0)then 
            !Check
            if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'
            jjj=jj+fltnum(6)*3
            allocate(btmp2(jjj),btmp3(jjj))
            ntagMPI=0
            do ix=1,numxyz(1)
                do iy=1,numxyz(2)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndu-1)*numxyz(2)+iy
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then
                            ntagMPI=ntagMPI+1
                            btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                        endif
                    enddo    
                enddo
            enddo
            !Check
            if (jj/=ntagMPI) then 
                stop 'jj&ntagMPI-qdct2-bndu'
                write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
            endif
!
              if (fltMPI(6)) then
                  do iz=1,fltnum(6)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltu(iz)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    btmp2(ntagMPI)=quantArray(id1(locid(nodenumtemp)+izz))
                      enddo
                enddo
              endif
             call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+1, 220000+me, &
                btmp3, jjj, MPI_DOUBLE_PRECISION, me+1, 220000+me+1,&
                MPI_COMM_WORLD, istatus, ierr)
!Update
            ntagMPI=0
            do ix=1,numxyz(1)
                do iy=1,numxyz(2)
                    nodenumtemp=(ix-1)*numxyz(2)*numxyz(3)+(bndu-1)*numxyz(2)+iy
                    do izz=1,dof1(nodenumtemp)
                        if (id1(locid(nodenumtemp)+izz)>0) then
                            ntagMPI=ntagMPI+1
                            quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+&
                                btmp3(ntagMPI)
                        endif
                    enddo    
                enddo
            enddo
              if (fltMPI(6)) then
                  do iz=1,fltnum(6)
                    nodenumtemp=numxyz(1)*numxyz(2)*numxyz(3)+fltu(iz)
                    do izz=1,3
                    ntagMPI=ntagMPI+1
                    quantArray(id1(locid(nodenumtemp)+izz))=quantArray(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
                      enddo
                enddo
              endif
            deallocate(btmp2, btmp3)            
        endif!bndu/=0
        !write(*,*) 'me',me,'Finished BNDU-quantArray'
      endif!if npxyz(3)>1
    call mpi_barrier(MPI_COMM_WORLD, ierr)
end subroutine MPI4NodalQuant