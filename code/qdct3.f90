SUBROUTINE qdct3(numel,numnp,neq,mat,ien,d,v,rdampk,rdampm,rho,ccosphi,sinphi,&
  mushr,elestress,eleporep,elemass,eleshp,eledet,pstrain,c,lm,brhs, me, master, nprocs)
  use globalvar
  implicit none
  include 'mpif.h'
  !
  !### program to calculate element contributions to residual force
  !        for the 8-node hexahedral, elastic continuum element 
  !        and assemble into the global right-hand-side vector.
  !		Note: only for 1-Gaussian point case now! B.D. 8/20/05
  !
  !  Explicitly use central difference method: al() always zero if no
  !	Rayleigh dampling (rdampm = 0), no global a() is needed.
  !	B.D. 7/21/05
  !
  logical :: formkd,zerovl,formma,zeroal
  integer (kind=4) :: nel,m,i,j,ntemp,k,k1,numel,numnp,neq
  real (kind=8) :: det,constk,pstrinc
  real (kind=8),dimension(nee) :: elresf,eleffm
  real (kind=8),dimension(ned,nen) :: dl,vl,al
  !...element arrays
  integer (kind=4),dimension(numel) :: mat	    
  integer (kind=4),dimension(nen,numel) :: ien
  integer (kind=4),dimension(ned,nen,numel) :: lm
  real (kind=8),dimension(numel) :: eledet,eleporep,pstrain
  real (kind=8),dimension(nee,numel) :: elemass
  real (kind=8),dimension(nstr,numel) :: elestress
  real (kind=8),dimension(nrowsh-1,nen,numel) :: eleshp
  !......nodes' arrays
  real (kind=8),dimension(ndof,numnp) :: d,v
  real (kind=8),dimension(neq) :: brhs
  !...material properties
  real (kind=8):: grav=-9.8  !gravity acceleration
  real (kind=8),dimension(numat) :: rho,rdampm,rdampk,ccosphi,sinphi,mushr
  real (kind=8),dimension(nrowc,nrowc,numat) :: c
  integer me, master, nprocs, rlp, rr, ierr,jj

!do we need to do this? numel now is smaller num, local. B.D. 4/15/09
!  rlp = numel/nprocs
!  rr  = numel-rlp*nprocs
!    if (me ==nprocs-1) then
!      jj=(me+1)*rlp+rr
!    else
!     jj=(me+1)*rlp
!    endif

!  if (numel /= rlp*nprocs) then
!        write(*,*) 'The number of processors', nprocs, 'is not suited to the total elements', numel
!        stop
!  endif

  !
  !*** loop over elements ***
  !
!$omp parallel do default(shared) private(nel,formma,formkd,m,ntemp,dl,vl,al,&
!$omp	j,i,zeroal,zerovl,elresf,det,eleffm,constk,k,k1)
!!$omp do private(formkd,zerodl,formma,zeroal,nel,j,i,k,m,ntemp,k1,dl,vl,al,det,constk,elresf,eleffm)
  do nel=1,numel
!  do nel=me*rlp+1,jj
    !
    formma = .false.
    formkd = .false.
    m = mat(nel)
    !...localize dl,vl,al
    do j=1,nen
      ntemp = ien(j,nel)
      do i=1,ned
	dl(i,j) = d(i,ntemp)
	vl(i,j) = v(i,ntemp)
	al(i,j) = 0.0
      enddo
    enddo
    !...compute effective dl accounting for Rayleigh damping
    do j=1,nen
      do i=1,ned
        !...for rate formulation, damping done in qdckd.f90. B.D. 1/5/12
        !dl(i,j) = dl(i,j) + rdampk(m)*vl(i,j)
        al(i,j) = al(i,j) + rdampm(m)*vl(i,j)
        if(i==3) then  !for inelastic off-fault, gravity included
           al(i,j) = al(i,j) - grav
        endif
      enddo
    enddo
    !...determine if element makes inertial contribution
    zeroal = .true.
    outer1: do j=1,nen	!Giving names to control constructs
      inner1: do i=1,ned
             !  k=id(i,j)
	       if(al(i,j) /= 0.0) then
	         zeroal = .false.
	         exit outer1
	       endif
             enddo inner1
           enddo outer1
    if ( (.not.zeroal) .and. (rho(m) /= 0.0) ) then
      formma = .true.
    endif
    !...determine if element makes stiffness contribution
    zerovl = .true.
    outer2: do j=1,nen
      inner2: do i=1,ned
		if(vl(i,j) /= 0.0) then
		  zerovl = .false.
		  exit outer2
		endif
	      enddo inner2
	    enddo outer2
  !  if (.not.zerodl) then
      formkd = .true.
  !  endif
    !...gravity effect and determine formma again
    !    at present, no gravity. zerog = .true. B.D. 3/26/05
    !	  zerog = .true.
    !     call ztest(grav,nesd,zerog)
    !  if ((.not.zerog) .and. (lfbody.ne.0) .and. (rho(m).ne.zero)	&
    ! &   .and. (imass.ne.2)) then
    !     formma = .true.
    !     do 400 i=1,ned
    !     temp = grav(i)*g1(lfbody)
    !         do 300 j=1,nen
    !           al(i,j) = al(i,j) - temp
    !300    continue
    !400  continue
    !  endif
    !*** if either, start time-consuming computing ***
    if (formma .or. formkd) then
      !...initialize element right hand vector
      elresf = 0.0
      !...directly assign shg and det to avoid recalculation to speed up.
      ! 	only for 1-point Gaussian case. 	B.D. 3/25/05
      !	   No shg needed in this new formation. B.D. 7/3/05
      !do i=1,nrowsh
      !	do j=1,nen
      !	  shg(i,j) = eleshap(i,j,nel)
      !	enddo
      !enddo
      det = eledet(nel)
      !*** form inertial and/or body force ***
      if (formma) then
        !......assign eleffm from previous stored. B.D. 3/26/05
	do i=1,nee
	  eleffm(i)=elemass(i,nel)
	enddo
        !......call to compute elresf	
        call contma(eleffm,al,elresf)
      endif
      !*** form internal force: most time-consuming part ***
      if (formkd) then
        !...... form internal force
        constk = - det
	call qdckd(eleshp(1,1,nel),c(1,1,m),vl,elestress(1,nel),elresf,constk,&
             zerovl,rdampk(m),ccosphi(m),sinphi(m),eleporep(nel),mushr(m),pstrinc)
        pstrain(nel) = pstrain(nel) + pstrinc
      endif
      !*** only either formma or formkd, assemble ***
      do i=1,nen
        do j=1,ned
          k=lm(j,i,nel)
	  k1=j+(i-1)*ned
	  if(k > 0) then
	    brhs(k) = brhs(k) + elresf(k1)
	  endif
        enddo
      enddo
    endif
  !
  enddo
!!$omp end do nowait
!$omp end parallel do
  !
  !*** form surface force ***
  !     note: assembly of surface loads is performed inside qdcsuf
  !	at present, no surface force applied! B.D. 7/2/05
  !if ( (nsurf.gt.0) .and. (lfsurf.gt.0)) then
  !   call qdcsuf(ielno,ien,x,xl,iside,mat,th,press,shear,elresf, &
  !             brhs,lm,g1(lfsurf),nsurf,nen,nsd,nesd,ned,nee,iopt)
  !endif
  !
end SUBROUTINE qdct3
