!/* Copyright (C) 2006-2020, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./

subroutine faulting

	use globalvar
	implicit none

	character(len=30) :: foutmov
	integer (kind = 4) :: i,i1,j,k,n,isn,imn,itmp,ifout,ift
	real (kind = dp) :: slipn,slips,slipd,slip,slipraten,sliprates,sliprated,&
						sliprate,xmu,mmast,mslav,mtotl,fnfault,fsfault,fdfault,tnrm,tstk, &
						tdip,taox,taoy,taoz,ttao,taoc,ftix,ftiy,ftiz,trupt,tr,&
						tmp1,tmp2,tmp3,tmp4,tnrm0,rcc,fa,fb
	real (kind = dp):: fvd(6,2,3)
	real (kind = dp) :: rr,R0,T,F,G,dtao0,dtao=0.0d0 !RSF
	real (kind = dp) :: statetmp, v_trial, T_coeff!RSF
	integer (kind=4) :: iv,ivmax  !RSF
	real (kind = dp) :: tstk0, tdip0, tstk1, tdip1, ttao1, taoc_old, taoc_new !RSF
	real (kind = dp) :: dxmudv, rsfeq, drsfeqdv, vtmp !RSF
	real (kind = dp) :: accn,accs,accd, accx, accy, accz, Rx, Ry, Rz, mr

	!===================================================================!
	do ift = 1, ntotft

		do i=1,nftnd(ift)	!just fault nodes
		!-------------------------------------------------------------------!
			!RSF nucleate by imposing a horizontal shear traction perturbation
			!2016.08.28
			if((C_nuclea==1).and.(friclaw == 3 .or. friclaw == 4 .or. friclaw == 5).and.(ift == nucfault)) then
			   R0 = 1500.0d0
			   dtao0 = 50.0d6
			   T = 1.0d0
			   F = 0.0d0
			   rr=sqrt((x(1,nsmp(1,i,ift))-xsource)**2+(x(3,nsmp(1,i,ift))-zsource)**2)    
			   if (rr<R0)    F=dexp(rr**2/(rr**2-R0**2))
			   G = 1.0d0
			   if (time<=T)  G=dexp((time-T)**2/(time*(time-2*T)))
			   dtao=dtao0*F*G
			endif
		!-------------------------------------------------------------------!	
			fnfault = fric(7,i,ift) !initial forces on the fault node
			fsfault = fric(8,i,ift)+dtao !norm, strike, dip components directly
			fdfault = 0.0d0
			isn = nsmp(1,i,ift)
			imn = nsmp(2,i,ift)
			do j=1,2  !1-slave, 2-master
				do k=1,3  !1-x comp, 2-y comp, 3-z comp
				!          fvd(k,j,1) = brhs(id(k,nsmp(j,i)))  !1-force
					fvd(k,j,1) = brhs(id1(locid(nsmp(j,i,ift))+k))  !1-force !DL 
					fvd(k,j,2) = v(k,nsmp(j,i,ift)) !2-vel
					fvd(k,j,3) = d(k,nsmp(j,i,ift)) !3-di,iftsp
					!fvd(k,j,3) = d(k,nsmp(j,i,ift)) + rdampk(1)*fvd(k,j,2) !3-di,iftsp
				enddo
			enddo
			!...resolve x,y,z components onto normal, strike and dip components.
			!   B.D. 1/26/07
			do j=1,3    !1-force,2-vel,3-disp
				do k=1,2  !1-slave,2-master
					fvd(4,k,j) = fvd(1,k,j)*un(1,i,ift) + fvd(2,k,j)*un(2,i,ift) + fvd(3,k,j)*un(3,i,ift)  !4-norm
					fvd(5,k,j) = fvd(1,k,j)*us(1,i,ift) + fvd(2,k,j)*us(2,i,ift) + fvd(3,k,j)*us(3,i,ift)  !5-strike
					fvd(6,k,j) = fvd(1,k,j)*ud(1,i,ift) + fvd(2,k,j)*ud(2,i,ift) + fvd(3,k,j)*ud(3,i,ift)  !6-dip
				enddo
			enddo
			slipn = fvd(4,2,3) - fvd(4,1,3)
			slips = fvd(5,2,3) - fvd(5,1,3)
			slipd = fvd(6,2,3) - fvd(6,1,3)
			slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
			fric(71,i,ift) = slips  
			fric(72,i,ift) = slipd
			fric(73,i,ift) = slipn  
			slipraten = fvd(4,2,2) - fvd(4,1,2)
			sliprates = fvd(5,2,2) - fvd(5,1,2)
			sliprated = fvd(6,2,2) - fvd(6,1,2)
			fric(74,i,ift) = sliprates  
			fric(75,i,ift) = sliprated
			sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
			if (sliprate>fric(76,i,ift)) then 
				fric(76,i,ift)=sliprate
			endif
			fric(77,i,ift) = fric(77,i,ift) + sliprate * dt
		
			mslav = fnms(isn)		
			mmast = fnms(imn)
			mtotl = mslav + mmast
			mtotl = mtotl * arn(i,ift)
			

			if (C_elastic==0) then!Plastic 	
				tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
					+ mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl          
				tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
					- mmast*fvd(5,1,1)) / mtotl
				tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
					- mmast*fvd(6,1,1)) / mtotl
			else!Elastic
				tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
					+ mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl + fnfault         
				tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
					- mmast*fvd(5,1,1)) / mtotl + fsfault
				tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
					- mmast*fvd(6,1,1)) / mtotl + fdfault
			endif		
			ttao = sqrt(tstk*tstk + tdip*tdip)   
		
			if (friclaw==1.or.friclaw==2)then!Differ 1&2 and 3&4	
				if(friclaw == 1) then
					call slip_weak(fric(77,i,ift),fric(1,i,ift),xmu)
				elseif(friclaw == 2) then
					trupt =  time - fnft(i,ift)
					call time_weak(trupt,fric(1,i,ift),xmu)
				endif

				if (C_Nuclea==1) then	
					if(r4nuc(i,ift)<=srcrad0) then !only within nucleation zone, do...
						tr=(r4nuc(i,ift)+0.081*srcrad0*(1./(1-(r4nuc(i,ift)/srcrad0)*(r4nuc(i,ift)/srcrad0))-1))/(0.7*3464.)
					else
						tr=1.0e9 
					endif
					if(time<tr) then 
						fb=0.0
					elseif ((time<(tr+critt0)).and.(time>=tr)) then 
						fb=(time-tr)/critt0
					else 
						fb=1.0
					endif
					tmp1=fric(1,i,ift)+(fric(2,i,ift)-fric(1,i,ift))*fb
					tmp2=xmu
					xmu=min(tmp1,tmp2)  !minimum friction used. B.D. 2/16/13	
				endif

				if((tnrm+fric(6,i,ift))>0) then
					tnrm0 = 0.0
				else
					tnrm0 = tnrm+fric(6,i,ift)
				endif
				taoc = fric(4,i,ift) - xmu *tnrm0

				if(ttao > taoc) then
					tstk = tstk * taoc / ttao
					tdip = tdip * taoc / ttao
					if(fnft(i,ift)>600) then	!fnft should be initialized by >10000
						if(sliprate >= 0.001) then	!first time to reach 1mm/s
							fnft(i,ift) = time	!rupture time for the node
						endif
					endif
				endif

				taox = (tnrm*un(1,i,ift) + tstk*us(1,i,ift) + tdip*ud(1,i,ift))*arn(i,ift)
				taoy = (tnrm*un(2,i,ift) + tstk*us(2,i,ift) + tdip*ud(2,i,ift))*arn(i,ift)
				taoz = (tnrm*un(3,i,ift) + tstk*us(3,i,ift) + tdip*ud(3,i,ift))*arn(i,ift)

				if (C_elastic==0) then!Plastic	
					brhs(id1(locid(isn)+1)) = brhs(id1(locid(isn)+1)) + taox !brhs(id1(loci(1,imn)+1))
					brhs(id1(locid(isn)+2)) = brhs(id1(locid(isn)+2)) + taoy
					brhs(id1(locid(isn)+3)) = brhs(id1(locid(isn)+3)) + taoz
					brhs(id1(locid(imn)+1)) = brhs(id1(locid(imn)+1)) - taox
					brhs(id1(locid(imn)+2)) = brhs(id1(locid(imn)+2)) - taoy
					brhs(id1(locid(imn)+3)) = brhs(id1(locid(imn)+3)) - taoz
				else!Elastic
					ftix = (fnfault*un(1,i,ift) + fsfault*us(1,i,ift) + fdfault*ud(1,i,ift))*arn(i,ift)
					ftiy = (fnfault*un(2,i,ift) + fsfault*us(2,i,ift) + fdfault*ud(2,i,ift))*arn(i,ift)
					ftiz = (fnfault*un(3,i,ift) + fsfault*us(3,i,ift) + fdfault*ud(3,i,ift))*arn(i,ift)  
					brhs(id1(locid(isn)+1)) = brhs(id1(locid(isn)+1)) + taox - ftix
					brhs(id1(locid(isn)+2)) = brhs(id1(locid(isn)+2)) + taoy - ftiy
					brhs(id1(locid(isn)+3)) = brhs(id1(locid(isn)+3)) + taoz - ftiz
					brhs(id1(locid(imn)+1)) = brhs(id1(locid(imn)+1)) - taox + ftix
					brhs(id1(locid(imn)+2)) = brhs(id1(locid(imn)+2)) - taoy + ftiy
					brhs(id1(locid(imn)+3)) = brhs(id1(locid(imn)+3)) - taoz + ftiz
				endif	
			elseif (friclaw>=3)then
			
				tnrm = tnrm + fric(51,i,ift) ! consider termopressurization.
					
				slipn = slipn + fric(16,i,ift) * time 
				slips = slips + fric(17,i,ift) * time
				slipd = slipd + fric(18,i,ift) * time
				slip = sqrt(slips**2 + slipd**2) !slip mag

				slipraten =  slipraten + fric(16,i,ift) 
				sliprates =  sliprates + fric(17,i,ift)
				sliprated =  sliprated + fric(18,i,ift)
				sliprate = sqrt(sliprates**2+sliprated**2)
				if(fnft(i,ift)>600.0d0) then	!fnft should be initialized by >10000
					if(sliprate >= 0.001d0) then	!first time to reach 1mm/s
						fnft(i,ift) = time	!rupture time for the node
					endif
				endif
				v_trial = sliprate
				mr =   mmast * mslav / (mmast+mslav) !reduced mass   
				T_coeff = arn(i,ift)* dt / mr
				statetmp = fric(21,i,ift)  !RSF: a temporary variable to store the currently value of state variable. B.L. 1/8/16
				if (abs(x(1,isn)-xsource)<1.0d0 .and. abs(x(3,isn)-zsource)<1.0d0) then 
					write(*,*) 'faulting,i,tnrm,tstk,v_trial,p,T,x,z',tnrm/1.0d6,tstk/1.0d6,v_trial,fric(51,i,ift),fric(52,i,ift),x(1,isn),x(3,isn)
				endif
				if(friclaw == 3) then
					call rate_state_ageing_law(v_trial,fric(21,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
				elseif (friclaw == 4 .or. friclaw==5) then
					call rate_state_slip_law(v_trial,fric(21,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
				endif 

				taoc_old = fric(4,i,ift) - xmu * MIN(tnrm, 0.0d0)
				tstk0=tstk
				tdip0=tdip
				tstk1 = tstk0 - taoc_old*0.5d0 * (sliprates / sliprate) + fric(17,i,ift)/T_coeff
				tdip1 = tdip0 - taoc_old*0.5d0 * (sliprated / sliprate) + fric(18,i,ift)/T_coeff
				  
				ttao1 = sqrt(tstk1*tstk1 + tdip1*tdip1)
			  
				ivmax = 30  !RSF: maximum 30 loops for iteration, once a criterion is met, jump out of this loop. B.L. 1/8/16
			  
				do iv = 1,ivmax
					fric(21,i,ift) = statetmp
					if(friclaw == 3) then
						call rate_state_ageing_law(v_trial,fric(21,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
					else
						call rate_state_slip_law(v_trial,fric(21,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
					endif 
			  
					taoc_new = fric(4,i,ift) - xmu * MIN(tnrm, 0.0d0)
			  
					rsfeq = v_trial + T_coeff * (taoc_new*0.5d0 - ttao1)
					drsfeqdv = 1.0d0 + T_coeff * (-dxmudv * MIN(tnrm,0.0d0))*0.5d0  
			 
					if(abs(rsfeq/drsfeqdv) < 1.d-14 * abs(v_trial) .and. abs(rsfeq) < 1.d-6 * abs(v_trial)) exit 
						vtmp = v_trial - rsfeq / drsfeqdv
					if(vtmp <= 0.0d0) then
						v_trial = v_trial/2.0d0
					else
						v_trial = vtmp
					endif  
				enddo !iv

				if(v_trial < fric(19,i,ift)) then
					v_trial = fric(19,i,ift)
					taoc_new = ttao1 * 2.0d0
				endif

				tstk = taoc_old*0.5d0 * (sliprates / sliprate) + taoc_new*0.5d0 * (tstk1 / ttao1) 
				tdip = taoc_old*0.5d0 * (sliprated / sliprate) + taoc_new*0.5d0 * (tdip1 / ttao1) 
				fric(78,i,ift)=tnrm 
				fric(79,i,ift)=tstk

				fric(49,i,ift) = v_trial
				fric(50,i,ift) = (tstk**2 + tdip**2)**0.5 
				frichis(1,i,nt,ift) = fric(49,i,ift)
				frichis(2,i,nt,ift) = fric(50,i,ift)
				accn = -slipraten/dt - slipn/dt/dt
				accs = (v_trial * (tstk1 / ttao1) - sliprates)/dt
				accd = (v_trial * (tdip1 / ttao1) - sliprated)/dt
			 
				accx = accn*un(1,i,ift) + accs*us(1,i,ift) + accd*ud(1,i,ift)
				accy = accn*un(2,i,ift) + accs*us(2,i,ift) + accd*ud(2,i,ift)
				accz = accn*un(3,i,ift) + accs*us(3,i,ift) + accd*ud(3,i,ift)
				if (C_elastic==0) then
					stop 777
				elseif (C_elastic==1) then
					Rx = brhs(id1(locid(isn)+1)) + brhs(id1(locid(imn)+1))
					Ry = brhs(id1(locid(isn)+2)) + brhs(id1(locid(imn)+2))
					Rz = brhs(id1(locid(isn)+3)) + brhs(id1(locid(imn)+3))
					brhs(id1(locid(isn)+1)) = (-accx + Rx/mmast) * mr
					brhs(id1(locid(isn)+2)) = (-accy + Ry/mmast) * mr
					brhs(id1(locid(isn)+3)) = (-accz + Rz/mmast) * mr
					brhs(id1(locid(imn)+1)) = (accx + Rx/mslav) * mr
					brhs(id1(locid(imn)+2)) = (accy + Ry/mslav) * mr
					brhs(id1(locid(imn)+3)) = (accz + Rz/mslav) * mr
					endif
				endif

				if(n4onf>0.and.lstr) then	
					do j=1,n4onf
						if(anonfs(1,j)==i.and.anonfs(3,j)==ift) then !only selected stations. B.D. 10/25/09    
							fltsta(1,locplt-1,j)  = time
							fltsta(2,locplt-1,j)  = sliprates
							fltsta(3,locplt-1,j)  = sliprated
							fltsta(4,locplt-1,j)  = fric(21,i,ift)
							fltsta(5,locplt-1,j)  = slips
							fltsta(6,locplt-1,j)  = slipd
							fltsta(7,locplt-1,j)  = slipn
							fltsta(8,locplt-1,j)  = tstk
							fltsta(9,locplt-1,j)  = tdip
							fltsta(10,locplt-1,j) = tnrm
							fltsta(11,locplt-1,j) = fric(51,i,ift) + fric_tp_pini
							fltsta(12,locplt-1,j) = fric(52,i,ift) 
						endif
					enddo 
				endif   
			if (x(1,isn)==xsource.and.x(2,isn)==ysource.and.x(3,isn)==zsource)then
				write(*,*)'S1:slip,ft',slips,fnft(i,ift),r4nuc(i,ift),tr
				write(*,*)'source,taoc,ttao',(taoc_old+taoc_new)/2,ttao
				!write(*,*)'source,tnrm,tstk,tdip',tnrm,tstk,tdip
				!write(*,*)'source,brhs isn',brhs(id1(locid(isn)+1)),brhs(id1(locid(isn)+2)),brhs(id1(locid(isn)+3))	
			endif
		enddo	!ending i
	enddo !ift
	!-------------------------------------------------------------------!
	!-------------Late Sep.2015/ D.Liu----------------------------------!
	!-----------Writing out results on fault for evert nstep------------!
	!if(mod(nt,315)==1.and.nt<5000) then 
	!	write(mm,'(i6)') me
	!	mm = trim(adjustl(mm))
	!	foutmov='fslipout_'//mm
	!	open(9002+me,file=foutmov,form='formatted',status='unknown',position='append')
	!		write(9002+me,'(1x,4f10.3)') ((fltslp(j,ifout),j=1,3),fltslr(1,ifout),ifout=1,nftnd)
	!endif
	!----nftnd for each me for plotting---------------------------------!
	!if (nt==1) then
	!	write(mm,'(i6)') me	
	!	mm = trim(adjustl(mm))			
	!	foutmov='fnode.txt'//mm
	!	open(unit=9800,file=foutmov,form='formatted',status='unknown')
	!		write(9800,'(2I7)') me,nftnd 
	!	close(9800)			
	!endif 	
	!-------------------------------------------------------------------!	
end subroutine faulting	 
