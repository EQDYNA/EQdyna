SUBROUTINE qdct3(numel,numnp,neq,mat,ien,d,v,eleporep,elemass,eleshp,eledet,pstrain,brhs,me, &
				maxm,id1,locid,dof1,et,v1,d1,PMLb,x,maxs,ids,s1,nt)
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
logical :: formkd,zerovl,formma,zeroal
integer(kind=4)::me,nel,m,i,j,ntemp,k,k1,numel,numnp,neq,maxm,&
	maxs,non,itag,eqn,label,nt,ien(nen,numel),id1(maxm),&
	ids(numel),locid(numnp),dof1(numnp),et(numel)
real(kind=8)::grav=-9.8,det,constk,pstrinc,xc(3),matelement(5),esPML(21),es(12),ex(3,8),evPML(96),efPML(96),PMLb(8),&
	x(nsd,numnp),s1(maxs),v1(neq),d1(neq),brhs(neq),d(ndof,numnp),v(ndof,numnp),eleshp(nrowsh-1,nen,numel),&
	elemass(nee,numel),eledet(numel),eleporep(numel),pstrain(numel),mat(numel,5),elresf(nee),eleffm(nee),&
	dl(ned,nen),vl(ned,nen),al(ned,nen)
	
!$omp parallel do default(shared) private(nel,formma,formkd,m,ntemp,dl,vl,al,&
!$omp	j,i,zeroal,zerovl,elresf,det,eleffm,constk,k,k1)
do nel=1,numel
	if (nel.lt.0) then 
		write(*,*) 'me=',me,'nel=',nel,'numel=',numel
			stop 3
		endif
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
			al(i,j) = al(i,j) + rdampm*vl(i,j)
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
		elresf = 0.0
		det = eledet(nel)
		if (formma) then
			do i=1,nee
				eleffm(i)=elemass(i,nel)
			enddo
			call contma(eleffm,al,elresf)
		endif
		if (formkd) then
			if (et(nel)==1) then
				constk=-det
				do i=1,12
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
						zerovl,eleporep(nel),pstrinc,ex,PMLb)
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
				call PMLwhg(vl,efPML,evPML,esPML,ex,PMLb,matelement,eleshp(1,1,nel),det,nt,nel,me)
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
!$omp end parallel do
end SUBROUTINE qdct3
