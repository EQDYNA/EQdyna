SUBROUTINE qdct3(numel,numnp,neq,mat,ien,d,v,rdampk,rdampm,ccosphi,sinphi,&
mushr,eleporep,elemass,eleshp,eledet,pstrain,c,brhs,me,master,nprocs, &! Delete lm elestress
maxm,id1,locid,dof1,et,v1,d1,PMLb,x,maxs,ids,s1,nt)!Adding maxm,id1,loci,v1,d1,PMLb,maxs,ids
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
real(kind=8),dimension(numel,5) :: mat	    
integer (kind=4),dimension(nen,numel) :: ien
real (kind=8),dimension(numel) :: eledet,eleporep,pstrain
real (kind=8),dimension(nee,numel) :: elemass
real (kind=8),dimension(nrowsh-1,nen,numel) :: eleshp
!......nodes' arrays
real (kind=8),dimension(ndof,numnp) :: d,v
real (kind=8),dimension(neq) :: brhs
!...material properties
real (kind=8):: grav=-9.8  !gravity acceleration
real (kind=8),dimension(numat) :: rdampm,rdampk,ccosphi,sinphi,mushr
real (kind=8),dimension(nrowc,nrowc,numat) :: c
integer me, master, nprocs, rlp, rr, ierr,jj
!*.* variables for PML. D.L. Jan/23/2015
integer (kind=4):: maxm,maxs,non,itag,eqn,label,nt
integer (kind=4),dimension(maxm)::id1
integer (kind=4),dimension(numel)::ids
integer (kind=4),dimension(numnp)::locid,dof1
integer (kind=4),dimension(numel)::et
real (kind=8),dimension(neq)::v1,d1
real (kind=8),dimension(maxs)::s1
real (kind=8),dimension(8)::PMLb
real (kind=8),dimension(96)::evPML,efPML
real(kind=8),dimension(3,8)::ex
real (kind=8),dimension(nsd,numnp) :: x
real(kind=8),dimension(12)::es
real(kind=8),dimension(21)::esPML
real(kind=8)::xc(3),matelement(5)
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
	if (nel.lt.0) then 
		write(*,*) 'me=',me,'nel=',nel,'numel=',numel
			stop 3
		endif
	!  do nel=me*rlp+1,jj
	!
	formma = .false.
	formkd = .false.
	m = 1
	!...localize dl,vl,al
	do j=1,nen
		ntemp = ien(j,nel)
		do i=1,ned
			dl(i,j) = d(i,ntemp)
			vl(i,j) = v(i,ntemp)
			al(i,j) = 0.0
		enddo
	enddo
	!...compute effective dl accounting for Rayleigh damping!
	do j=1,nen
		do i=1,ned
		!...for rate formulation, damping done in qdckd.f90. B.D. 1/5/12
		!dl(i,j) = dl(i,j) + rdampk(m)*vl(i,j)
			al(i,j) = al(i,j) + rdampm(m)*vl(i,j)
			if(i==3.and.C_elastic==0) then  !for inelastic off-fault, gravity included
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
	if ( (.not.zeroal) .and. (mat(nel,3) /= 0.0) ) then
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
	if (.not.zerovl) then
		formkd = .true.
	endif
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
		if (formkd) then
			if (et(nel)==1) then
				!*** form internal force: most time-consuming part ***
				!...... form internal force
				constk = - det
				do i=1,12!DL using es intead of elestress(1,nel)
					es(i)=s1(ids(nel)+i)
				enddo
				do i=1,8
					do j=1,3
						ex(j,i)=x(j,ien(i,nel))
					enddo
				enddo
				do i=1,5 
					matelement(i)=mat(nel,i)
				enddo
				call qdckd(eleshp(1,1,nel),matelement,vl,dl,es,elresf,constk,&
						zerovl,rdampk(m),ccosphi(m),sinphi(m),eleporep(nel),pstrinc,ex,PMLb)
				pstrain(nel) = pstrain(nel) + pstrinc
				do i=1,12!DL update the stress
					s1(ids(nel)+i)=es(i)
				enddo				
				do i=1,nen					
					non=ien(i,nel)
					do j=1,ned
						itag=locid(non)+j
						k=id1(itag)					
						if(k > 0) then
							brhs(k) = brhs(k) + elresf((i-1)*ned+j)
						endif
					enddo
				enddo
			elseif (et(nel)==2) then ! PML element. 
				efPML=0.0
				evPML=0.0
				esPML=0.0 
				do i=1,8
					do j=1,3
					efPML((i-1)*12+9+j)=elresf((i-1)*3+j)
					enddo
				enddo
				do i=1,8
					do j=1,3
						ex(j,i)=x(j,ien(i,nel))
					enddo
				enddo
				do i=1,21
					esPML(i)=s1(ids(nel)+i)
				enddo
				do i=1,5 
					matelement(i)=mat(nel,i)
				enddo		
				call PMLwhg(vl,efPML,evPML,esPML,ex,PMLb,matelement,eleshp(1,1,nel),det,nt,rdampk(m),nel,me)
				do i=1,8
					non=ien(i,nel)
					if (dof1(non)==12) then
						do j=1,12
							itag=locid(non)+j
							eqn=id1(itag)
							if (eqn.gt.0) then
								brhs(eqn)=brhs(eqn)+efPML((i-1)*12+j)
							endif
						enddo
					elseif (dof1(non)==3) then
						itag=locid(non)+1
						eqn=id1(itag)
						brhs(eqn)=brhs(eqn)+efPML((i-1)*12+1)+efPML((i-1)*12+2)+efPML((i-1)*12+3)+efPML((i-1)*12+10)
						itag=locid(non)+2
						eqn=id1(itag)
						brhs(eqn)=brhs(eqn)+efPML((i-1)*12+4)+efPML((i-1)*12+5)+efPML((i-1)*12+6)+efPML((i-1)*12+11)
						itag=locid(non)+3
						eqn=id1(itag)
						brhs(eqn)=brhs(eqn)+efPML((i-1)*12+7)+efPML((i-1)*12+8)+efPML((i-1)*12+9)+efPML((i-1)*12+12)				
					endif
				enddo
				!update the stress components 
				do i=1,21
					s1(ids(nel)+i)=esPML(i)
				enddo
			endif!et(nel)=1/2
			
			! if (me==10.and.nel==206232) then 
				! write(*,*) 'NEL1-3',vl(1,1),vl(2,1),vl(3,1)
				! write(*,*) 'NEL4-6',vl(1,2),vl(2,2),vl(3,2)
				! write(*,*) 'NEL7-9',vl(1,3),vl(2,3),vl(3,3)
				! write(*,*) 'NEL10-12',vl(1,4),vl(2,4),vl(3,4)
				! write(*,*) 'NEL13-15',vl(1,5),vl(2,5),vl(3,5)
				! write(*,*) 'NEL16-18',vl(1,6),vl(2,6),vl(3,6)
				! write(*,*) 'NEL19-21',vl(1,7),vl(2,7),vl(3,7)
				! write(*,*) 'NEL22-24',vl(1,8),vl(2,8),vl(3,8)
			! endif
		endif!formkd
	endif!either formma or formkd
enddo!nel loop
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
