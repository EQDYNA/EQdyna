SUBROUTINE qdct3
use globalvar
implicit none
include 'mpif.h'

logical :: formkd,zerovl,formma,zeroal
integer(kind=4)::nel,m,i,j,ntemp,k,k1, non,itag,eqn,label
real(kind=8)::det,constk,pstrinc,xc(3),matelement(5),esPML(21),es(12),ex(3,8),efPML(96),&
	elresf(nee), eleffm(nee), dl(ned,nen), vl(ned,nen), al(ned,nen)
	

do nel=1,numel

	formma = .false.
	formkd = .false.
	
	m = 1
	!...localize dl,vl,al
	do j=1,nen
		ntemp = ien(j,nel)
		do i=1,ned
			dl(i,j) = d(i,ntemp)
			vl(i,j) = v(i,ntemp)
			al(i,j) = 0.0d0
		enddo
	enddo
	!...compute effective dl accounting for Rayleigh damping!
	do j=1,nen
		do i=1,ned
		!...for rate formulation, damping done in qdckd.f90. B.D. 1/5/12
		!dl(i,j) = dl(i,j) + rdampk(m)*vl(i,j)
			al(i,j) = al(i,j) + rdampm*vl(i,j)
			if(i==3.and.C_elastic==0) then  !for inelastic off-fault, gravity included
				al(i,j) = al(i,j) + grav
			endif
		enddo
	enddo
	!...determine if element makes inertial contribution
	zeroal = .true.
	outer1: do j=1,nen	!Giving names to control constructs
	inner1: do i=1,ned
	!  k=id(i,j)
		if(al(i,j) /= 0.0d0) then
			zeroal = .false.
			exit outer1
		endif
		enddo inner1
		enddo outer1
	if ( (.not.zeroal) .and. (mat(nel,3) /= 0.0d0) ) then
		formma = .true.
	endif
	!...determine if element makes stiffness contribution
	zerovl = .true.
	outer2: do j=1,nen
	inner2: do i=1,ned
				if(vl(i,j) /= 0.0d0) then
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
	!if (formma .or. formkd) then
		elresf = 0.0d0 
		det = eledet(nel)
		!if (formma) then
			do i=1,nee
				eleffm(i)=elemass(i,nel)
			enddo
			call contma(eleffm,al,elresf)
		!endif
		!if (formkd) then
			if (et(nel)==1.or.et(nel)>10) then
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
						zerovl,eleporep(nel),pstrinc,ex)
				! if (nel==584309.and.me==31) then 
					! write(*,*) es(1),es(2),es(3),es(4),es(5),es(6)
					! write(*,*) es(1+6),es(2+6),es(3+6),es(4+6),es(5+6),es(6+6)
					! write(*,*) elresf(1),elresf(2),elresf(3)
					! write(*,*) elresf(4),elresf(5),elresf(6)
					! write(*,*) elresf(7),elresf(8),elresf(9)
				! endif
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
				efPML=0.0d0
				esPML=0.0d0 
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
				call PMLwhg(vl,efPML,esPML,ex,matelement,eleshp(1,1,nel),det,nel)
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
		!endif!formkd
	!endif!either formma or formkd
enddo!nel loop
end SUBROUTINE qdct3
