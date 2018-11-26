SUBROUTINE faulting(ift,nftnd,numnp,neq,lstr,lstr1,fnms,brhs,d,v,x,maxm,id1,locid,dof1,& ! Delete id. Add maxm,id1,loci
					n4onf,mushr,momnt,momntrat,maxslprat,fltsta, &
					nsmp,fnft,fltslp,un,us,ud,fric,arn,r4nuc,arn4m,slp4fri,anonfs,itmp,&
					me, master,nprocs)
use globalvar
implicit none
!
!### program to implement faulting boundary conditions.
!	Call this routine after "form_rhs" to correct right-hand-side force
!	vector by taking into fault boundary. Basic coding was done in Feb,
!	2005. Now rewrite. B.D. 7/3/05
! ...totally rewrite this routine in November of 2006 by implementing
!     Day et al. (2005) formulation which does not require to treat
!     slipping, healing seperately. Much nicer! B.D. 11/23/06
! ...extend to 3D case. B.D. 1/26/07
!  
logical :: lstr,lstr1
integer (kind=4) :: ift,nftnd,numnp,neq,i,i1,j,k,n,isn,imn,n4onf,itmp
real (kind=8) ::slipn,slips,slipd,slip,slipraten,sliprates,sliprated,&
sliprate,xmu,mmast,mslav,mtotl,fnfault,fsfault,fdfault,tnrm,tstk, &
tdip,taox,taoy,taoz,ttao,taoc,ftix,ftiy,ftiz,trupt,tr,&
tmp1,tmp2,tmp3,tmp4,momnt,momntrat,maxslprat,tnrm0
integer (kind=4),dimension(3,itmp) :: anonfs
integer (kind=4),dimension(2,nftnd) :: nsmp
real (kind=8),dimension(nftnd) :: fnft,arn,r4nuc,arn4m,slp4fri
real (kind=8),dimension(3,nftnd) :: un,us,ud,fltslp
real (kind=8),dimension(6,nftnd) :: fric
real (kind=8),dimension(10,nplpts-1,n4onf) :: fltsta
real (kind=8),dimension(6,2,3)::fvd=0.0
real (kind=8),dimension(neq) :: brhs
real (kind=8),dimension(ndof,numnp) :: d,v,x
!  integer (kind=4),dimension(ndof,numnp) :: id
real (kind=8),dimension(numnp) :: fnms
real (kind=8),dimension(numat) :: mushr
integer master, me, nprocs
integer (kind=4):: maxm
integer (kind=4),dimension(maxm)::id1
integer (kind=4),dimension(numnp)::locid,dof1  
!
momnt=0  !for moment and rate calculation
momntrat=0 
maxslprat=0  !for max slip rate
!...do not use OpenMP for better output for ExGM 100runs. B.D. 8/12/10
!*** loop over slave nodes ***
!!$omp parallel do default(shared) private(i,j,k,fnfault,fsfault,fdfault,isn,imn,fvd, &
!!$omp	tmp1,tmp2,tmp3,tmp4,slipn,slips,slipd,slip,slipraten,sliprates,sliprated,sliprate, & 
!!$omp	mslav,mmast,mtotl,tnrm,tstk,tdip,ttao,taoc,taox,taoy,taoz,ftix,ftiy,ftiz,&
!!$omp	xmu,trupt)
!...the above definition of private is very important. xmu was not defined as private
!	ealier and resulted in problems: it can be imaged that it should if different
!	OpenMP threads mess up xmu! B.D. 10/31/09
do i=1,nftnd	!just fault nodes
	isn = nsmp(1,i)
	imn = nsmp(2,i)
	!...get nodal force,velocity, and diplacement in x,y,z.
	!   B.D. 1/26/07
	!...aslo add Rayleigh stiffness damping before using d.
	!   assume stifness coefficient is first material: rdampk(1).
	!   B.D. 11/26/06
	!...it seems the damping should not be used here.
	!   B.D. 1/28/07
	do j=1,2  !1-slave, 2-master
		do k=1,3  !1-x comp, 2-y comp, 3-z comp
		!          fvd(k,j,1) = brhs(id(k,nsmp(j,i)))  !1-force
			fvd(k,j,1) = brhs(id1(locid(nsmp(j,i))+k))  !1-force !DL 
			fvd(k,j,2) = v(k,nsmp(j,i)) !2-vel
			fvd(k,j,3) = d(k,nsmp(j,i)) !3-di,iftsp
			!fvd(k,j,3) = d(k,nsmp(j,i,ift)) + rdampk(1)*fvd(k,j,2) !3-di,iftsp
		enddo
	enddo
	!...resolve x,y,z components onto normal, strike and dip components.
	!   B.D. 1/26/07
	do j=1,3    !1-force,2-vel,3-disp
		do k=1,2  !1-slave,2-master
			fvd(4,k,j) = fvd(1,k,j)*un(1,i) + fvd(2,k,j)*un(2,i) + fvd(3,k,j)*un(3,i)  !4-norm
			fvd(5,k,j) = fvd(1,k,j)*us(1,i) + fvd(2,k,j)*us(2,i) + fvd(3,k,j)*us(3,i)  !5-strike
			fvd(6,k,j) = fvd(1,k,j)*ud(1,i) + fvd(2,k,j)*ud(2,i) + fvd(3,k,j)*ud(3,i)  !6-dip
		enddo
	enddo
	!
	!...no opening and no penetrating constraints for SCEC TPV12/13. B.D. 10/22/09
	!   should apply before trial traction.
	!   use the average for normal v and d. B.D. 10/22/09
	!Actually, my implementation of Day et al. (2005) here may already take care
	!of this and an extra constraint caused incorrect slip (rate) at surface fault
	!station in TPV13 (dip-slip). So should not apply again here. B.D. 11/2/09
	!tmp1 = fvd(4,1,2)
	!tmp2 = fvd(4,2,2)
	!tmp3 = fvd(4,1,3)
	!tmp4 = fvd(4,2,3)
	!fvd(4,1,2) = 0.5*(tmp1+tmp2)
	!fvd(4,2,2) = fvd(4,1,2)
	!fvd(4,1,3) = 0.5*(tmp3+tmp4)
	!fvd(4,2,3) = fvd(4,1,3)
	!...slip and slip rate should be calculated from norm, strike, and dip so that
	!    the above no opening or penetrating constraints included. B.D. 10/22/09
	slipn = fvd(4,2,3) - fvd(4,1,3)
	slips = fvd(5,2,3) - fvd(5,1,3)
	slipd = fvd(6,2,3) - fvd(6,1,3)
	slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
	fltslp(1,i) = slips  !save for final slip output
	fltslp(2,i) = slipd
	fltslp(3,i) = slipn  !normal should be zero, but still keep to ensure
	slipraten = fvd(4,2,2) - fvd(4,1,2)
	sliprates = fvd(5,2,2) - fvd(5,1,2)
	sliprated = fvd(6,2,2) - fvd(6,1,2)
	sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
	!...path-itegrated slip for slip-weakening. B.D. 8/12/10
	slp4fri(i) = slp4fri(i) + sliprate * dt
	!...calculate moment rate and moment if needed. B.D. 8/11/10
	!  also, max slip rate for early termination.
	!...for homogeneous material, i.e., only myshr(1). B.D. 1/3/12
	! or for heterogeneous case, but mushr(1) for rupture fault.
	if(time > term-dt) then  !only at the end, do this for LVFZ3D Plastic.
		momnt = momnt + mushr(1) * arn4m(i) * slip
		momntrat = momntrat + mushr(1) * arn4m(i) *sliprate
		if(sliprate>maxslprat) maxslprat=sliprate
	endif
	!
	!...nodal mass. Mass of each element may not be distributed among its 
	! nodes evenly. Instead, distribution is related to element shape. 
	!   Note: nodal mass should not be directly obtained from left-hand-side
	! diagnoal mass matrix, because that's the effective mass, which takes 
	! damping coefficient into accout. Instead, I computed nodal mass from 
	! element mass and assembled in "qdct2.f90".B.D.7/3/05
	mslav = fnms(isn)		
	mmast = fnms(imn)
	mtotl = mslav + mmast
	!
	!...trial traction to enforce continuity. B.D. 11/23/06
	!...divided by the associated area to get traction from force for EQdyna3d v2.1.2.
	!   initial stress, not initial force (f*fault) used here. B.D. 2/28/08
	!...no fault initial stress in elastoplastic rheology. B.D. 1/8/12
	mtotl = mtotl * arn(i)
	tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
		+ mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl          
	tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
		- mmast*fvd(5,1,1)) / mtotl
	tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
		- mmast*fvd(6,1,1)) / mtotl
	ttao = sqrt(tstk*tstk + tdip*tdip) !total shear magnitude	    
	!
	!...friction law to determine friction coefficient
	!   slip-weakening only so far. B.D. 1/26/07
	!... based on choices, call corresponding friction laws.
	! B.D. 10/8/08
	if(friclaw == 1) then
		call slip_weak(slp4fri(i),fric(1,i),xmu)
	elseif(friclaw == 2) then
		trupt =  time - fnft(i)
		call time_weak(trupt,fric(1,i),xmu)
	endif
	!......for nucleation zone of the nucleation fault,which initiates rupture,
	!	rupture propagates at a fixed speed to drop "xmu". B.D. 8/31/06
	if(ift == nucfault .and. xmu > fric(2,i)) then	
		!only nucleation fault and before finishing dropping, do...
		if(r4nuc(i) <= srcrad0) then !only within nucleation zone, do...
			tr = r4nuc(i) / vrupt0
			if(tr <= time) then !only ready or already fail, do...
				trupt = time - tr
			call time_weak(trupt,fric(1,i),xmu)
			endif
		endif
	endif
	!
	!...adjust tstk,tdip and tnrm based on jump conditions on fault.
	!   before calculate taoc, first adjust tnrm if needed. 
	!   after this, they are true (corrected) values. B.D. 11/23/06
	!   cohesion is added here. B.D. 2/28/08
	!...for SCEC TPV10/11, no opening is allowed. B.D. 11/24/08
	!if(tnrm > 0) tnrm = 0   !norm must be <= 0, otherwise no adjust
	!taoc = cohes - xmu * tnrm
	!...for ExGM 100 runs, no opening allowed means following.
	!  B.D. 8/12/10
	if((tnrm+fric(6,i))>0) then
		tnrm0 = 0.0
	else
		tnrm0 = tnrm+fric(6,i)
	endif
	taoc = fric(4,i) - xmu *tnrm0
	!taoc = cohes - xmu * tnrm0
	!if(tnrm > 0) tnrm = 0   !norm must be <= 0, otherwise no adjust
	!taoc = fistr(5,i) - xmu * tnrm
	if(ttao > taoc) then
		tstk = tstk * taoc / ttao
		tdip = tdip * taoc / ttao
		if(fnft(i)>600) then	!fnft should be initialized by >10000
			if(sliprate >= 0.001) then	!first time to reach 1mm/s
				fnft(i) = time	!rupture time for the node
			endif
		endif
	endif
	!
	!...add the above fault boundary force and initial force to elastic
	!	force of the split nodes. 
	!   first resolve normal, strike and dip back to x-,y-,z-. 
	!   then subtract them from slave, add to master as the above calculation
	!   based on this convention. see Day et al. (2005). B.D. 11/23/06
	!...due to traction, not force used in friction law above, need area to 
	!   convert traction to force for v2.1.2. B.D. 2/28/08
	taox = (tnrm*un(1,i) + tstk*us(1,i) + tdip*ud(1,i))*arn(i)
	taoy = (tnrm*un(2,i) + tstk*us(2,i) + tdip*ud(2,i))*arn(i)
	taoz = (tnrm*un(3,i) + tstk*us(3,i) + tdip*ud(3,i))*arn(i)
	!DL
	brhs(id1(locid(isn)+1)) = brhs(id1(locid(isn)+1)) + taox !brhs(id1(loci(1,imn)+1))
	brhs(id1(locid(isn)+2)) = brhs(id1(locid(isn)+2)) + taoy
	brhs(id1(locid(isn)+3)) = brhs(id1(locid(isn)+3)) + taoz
	brhs(id1(locid(imn)+1)) = brhs(id1(locid(imn)+1)) - taox
	brhs(id1(locid(imn)+2)) = brhs(id1(locid(imn)+2)) - taoy
	brhs(id1(locid(imn)+3)) = brhs(id1(locid(imn)+3)) - taoz
	!DL
	!brhs(id(1,isn)) = brhs(id(1,isn)) + taox
	!brhs(id(2,isn)) = brhs(id(2,isn)) + taoy
	!brhs(id(3,isn)) = brhs(id(3,isn)) + taoz
	!brhs(id(1,imn)) = brhs(id(1,imn)) - taox
	!brhs(id(2,imn)) = brhs(id(2,imn)) - taoy
	!brhs(id(3,imn)) = brhs(id(3,imn)) - taoz
	!
	!...Store fault forces and slip/slipvel for fault nodes 
	!		at set time interval.
	! note: forces will be transferred to stress later
	! B.D. 8/21/05
	!...now, they are directly traction (stress) in version 2.1.2.
	!   and can be written out here. B.D. 2/28/08
	!...Store only, no write out. B.D. 10/25/09
	if(n4onf>0 .and. lstr) then	
		do j=1,n4onf
			if(anonfs(1,j)==i.and.anonfs(3,j)==ift) then !only selected stations. B.D. 10/25/09    
				fltsta(1,locplt-1,j) = time
				fltsta(2,locplt-1,j) = sliprates
				fltsta(3,locplt-1,j) = sliprated
				fltsta(4,locplt-1,j) = slipraten
				fltsta(5,locplt-1,j) = slips
				fltsta(6,locplt-1,j) = slipd
				fltsta(7,locplt-1,j) = slipn
				fltsta(8,locplt-1,j) = tstk
				fltsta(9,locplt-1,j) = tdip
				fltsta(10,locplt-1,j) = tnrm+fric(6,i)
			endif
		enddo 
	endif
	!...slip rate output. B.D. 8/11/10
	!    if(lstr1) then
	!      myrec = myrec + 1
	!      tmp1 = sqrt(x(2,isn)*x(2,isn)+x(3,isn)*x(3,isn))
	!      write(ioutrat,rec=myrec) real( x(1,isn)),real(tmp1),&
	!           real(sliprates),real(-sliprated)
	!      call flush(ioutrat)
	!    endif
	!  if(lstr) then
	!      write(ioutsl,'(1x,3f10.1,12f10.3)') (x(j,isn),j=1,3),(((fvd(j,k,n),j=1,3), &
	!      		k=1,2),n=2,3)
	!flush_ for IBM systems
	!call flush_(ioutsl)
	!  call flush(ioutsl)
	!   if(i==1) write(ioutst,*) 'time=', time
	!   write(ioutst,'(1x,3f10.1,4e18.7e4)') (x(j,isn),j=1,3),tstk,tdip,taoc,tnrm      		
	!call flush_(ioutst)
	!call flush(ioutst)
	!  endif
	!    
	!   enddo	!ending i1
enddo	!ending i
!!$omp end parallel do
end SUBROUTINE faulting	 