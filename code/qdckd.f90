SUBROUTINE qdckd(shg,cc,vl,stress,elresf,constk,zerovl,bdamp,&
ccosphi,sinphi,porep,mushr,pstrmag,ex,PMLb)
use globalvar
implicit none
!
!### program to form internal force ("-k*d") for a continuum element
!        with "nen" nodes
!	Avoid calling subroutines for matrix operations to speed up.
!	B.D. 3/24/05
!	Further speed up by taking into account zero in B and D 
!		for 3-D, isotropic elasticity.
!	B.D. 8/20/05.
!	Rewrite again: strain-displacement matrix "b" was stored before
!		time loops. B.D. 7/2/05
!     Revise for stress rate-velocity formulation for elastoplastic
!		calculations. B.D. 1/5/12
!
logical :: zerovl
integer (kind=4) :: i,j,j1,j2,j3
real (kind=8) :: constk,bdamp,temp,ccosphi,sinphi,porep,mushr
real (kind=8),dimension(nee) :: elresf,work,vl
real (kind=8),dimension(nstr) :: strainrate,stressrate,stress,strtemp
real (kind=8),dimension(nrowsh-1,nen) :: shg
real (kind=8),dimension(nrowb,nee) :: bb	!correspond to b
real (kind=8),dimension(nrowc,nrowc) :: cc	!correspond to c
!...plasticity vlrables. B.D. 1/5/12
real (kind=8) :: strmea,taomax, yield, rjust,pstrmea,pstrmag,xc(3),ex(3,8),PMLb(6)
real (kind=8),dimension(nstr) :: strdev,pstrinc
!
pstrmag = 0.0
stressrate = 0.0
!...calcuate b from shg
call qdcb(shg,bb)
if(.not.zerovl) then !only nonzero velocity, update stress. B.D. 1/5/12
	!...calculate strainrate
	! Take into account zero in bb. B.D. 8/20/05
	strainrate = 0.0	!initialize
	!Important: OpenMP directives should not be used in this routine as
	!  they have been used in the parent roution qdct3.f90. Adding these
	!  directives did not cause slowdown problem with SUN compiler, but 
	!  results in significant reduction in comupting speed with Intel 
	!  compiler!! B.D. 8/7/10
	!!$omp parallel do default(shared) privlte(i,j1,j2,j3)
	do i=1,nen
		j1 = ned * (i - 1) + 1
		j2 = ned * (i - 1) + 2
		j3 = ned * (i - 1) + 3      
		strainrate(1) = strainrate(1) + bb(1,j1) * vl(j1)
		strainrate(2) = strainrate(2) + bb(2,j2) * vl(j2)
		strainrate(3) = strainrate(3) + bb(3,j3) * vl(j3)
		strainrate(4) = strainrate(4) + bb(4,j2) * vl(j2) + bb(4,j3) * vl(j3)
		strainrate(5) = strainrate(5) + bb(5,j1) * vl(j1) + bb(5,j3) * vl(j3)
		strainrate(6) = strainrate(6) + bb(6,j1) * vl(j1) + bb(6,j2) * vl(j2)
	enddo
	!!$omp end parallel do
	!...calculate stressrate
	! Take into account zero in cc. B.D. 8/20/05
	do i=1,3
		do j=1,3
			stressrate(i) = stressrate(i) + cc(i,j)*strainrate(j)
		enddo
	enddo
	do i=4,6
		stressrate(i) = cc(i,i) * strainrate(i)
	enddo
	!...calculate stress from stress rate & previous stress. B.D. 1/5/12
	do i=1,6
		stress(i) = stress(i) + stressrate(i) * dt
		strdev(i) = stress(i)
	enddo
	!...Drucker-Prager plasticity in shear. B.D. 1/5/12
	! refer to the benchmark problem description.
	strmea = (stress(1)+stress(2)+stress(3))/3.
	do i=1,3
		strdev(i) = stress(i) - strmea
	enddo
	taomax = 0.5*(strdev(1)**2+strdev(2)**2+strdev(3)**2) &
		+ strdev(4)**2+strdev(5)**2+strdev(6)**2  !second invlriant of deviator
	taomax = sqrt(taomax)  !kind of max shear
	yield = ccosphi - sinphi * (strmea + porep)  !yield stress
	if(yield<0) yield = 0  !non-negative
	if(taomax > yield) then  !yielding, stress adjust in deviator domain
		!rjust = yield/taomax  !adjust ratio
		!implement viscoplasticity now. B.D. 6/2/12
		rjust=yield/taomax + (1-yield/taomax)*exp(-dt/tv)
		do i=1,6  !adjust stress
			stress(i) =strdev(i) * rjust
			!calculate plastic strain increment components
			pstrinc(i) = (strdev(i) - stress(i))/mushr
			!back to stress domain by ading mean, which does not change
			if(i<=3) then
				stress(i) = stress(i) + strmea
			endif
		enddo
		!calculate plastic strain increment scalar (amplitude)
		pstrmea = (pstrinc(1)+pstrinc(2)+pstrinc(3))/3.
		do i=1,6
			pstrinc(i)=pstrinc(i) - pstrmea
		enddo
		pstrmag = 0.5*(pstrinc(1)**2+pstrinc(2)**2+pstrinc(3)**2) &
				+pstrinc(4)**2+pstrinc(5)**2+pstrinc(6)**2
		pstrmag = sqrt(pstrmag)
	endif
	!
endif !endif not zerovl (update stress)
!...calcuate element internal force
! Take into account zero in bbT. B.D. 8/20/05
! Now, damping force and internal force are separated in
!   rate formulation. But can put together for total
!   internal force here. B.D. 1/5/12
temp = constk * w
do i=1,nstr
	strtemp(i) = temp * (stress(i) + bdamp * stressrate(i))
	!strtemp(i) = temp * stress(i)
enddo
!!$omp parallel do default(shared) privlte(i,j1,j2,j3)
do i=1,nen
	j1 = ned * (i - 1) + 1
	j2 = ned * (i - 1) + 2
	j3 = ned * (i - 1) + 3      
	work(j1) = bb(1,j1)*strtemp(1) + bb(5,j1)*strtemp(5) &
			+ bb(6,j1)*strtemp(6)
	work(j2) = bb(2,j2)*strtemp(2) + bb(4,j2)*strtemp(4) &
			+ bb(6,j2)*strtemp(6)
	work(j3) = bb(3,j3)*strtemp(3) + bb(4,j3)*strtemp(4) &
			+ bb(5,j3)*strtemp(5)
enddo		
!!$omp end parallel do
!!$omp parallel do default(shared) privlte(i)
!...add into already element forces
do i=1,nee
	elresf(i) = elresf(i) + work(i)
enddo
	!xc=0.0
	!do i=1,3
	!	do j=1,8
	!		xc(i)=xc(i)+ex(i,j)
	!	enddo
	!enddo
	!xc=xc/8
	!if(xc(1)<51.and.xc(1)>49.and.xc(2)<51.and.xc(2)>49.and.(xc(3)<(PMLb(5)+200))) then
		!write(*,*) 'qdckd 1'
		!write(*,*) 'vel1',vl(1),vl(2),vl(3)
		!write(*,*) 'vel2',vl(4),vl(5),vl(6)
		!write(*,*) 'vel3',vl(7),vl(8),vl(9)
		!write(*,*) 'vel4',vl(10),vl(11),vl(12)
		!write(*,*) 'vel5',vl(13),vl(14),vl(15)
		!write(*,*) 'vel6',vl(16),vl(17),vl(18)
		!write(*,*) 'vel7',vl(19),vl(20),vl(21)
		!write(*,*) 'vel8',vl(22),vl(23),vl(24)
		!write(*,*) 'str',stress(1),stress(2),stress(3),stress(4),stress(5),stress(6)
	!endif
!!$omp end parallel do
end SUBROUTINE qdckd