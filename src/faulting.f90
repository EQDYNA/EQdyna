!/* Copyright (C) 2006-2023, Earthquake Modeling Lab @ Texas A&M University. 
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
    real (kind = dp) :: fvd(6,2,3)
    real (kind = dp) :: dtau0,dtau
    real (kind = dp) :: statetmp, v_trial, T_coeff!RSF
    integer (kind=4) :: iv,ivmax  !RSF
    real (kind = dp) :: tstk0, tdip0, tstk1, tdip1, ttao1, taoc_old, taoc_new !RSF
    real (kind = dp) :: dxmudv, rsfeq, drsfeqdv, vtmp, theta_pc_tmp !RSF
    real (kind = dp) :: accn,accs,accd, accx, accy, accz, Rx, Ry, Rz, mr, theta_pc, theta_pc_dot
    real (kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4), nsdTractionVector(4)
    !===================================================================!
    do ift = 1, ntotft

        do i=1,nftnd(ift)    !just fault nodes
        !-------------------------------------------------------------------!
            if (TPV == 104 .or. TPV == 105) then
                if (C_nuclea == 1 .and.ift == nucfault) then
                    call nucleation(dtau, xmu, x(1,nsmp(1,i,ift)), x(2,nsmp(1,i,ift)), & 
                                x(3,nsmp(1,i,ift)), fric(5,i,ift), fric(1,i,ift), &
                                fric(2,i,ift))
                endif
            endif 
            
            call getNsdSlipSliprateTraction(ift, i, nsdSlipVector, nsdSliprateVector, nsdTractionVector, dtau)
            slipn = nsdSlipVector(1)
            slips = nsdSlipVector(2)
            slipd = nsdSlipVector(3)
            slip = nsdSlipVector(4)
            
            slipraten = nsdSliprateVector(1)
            sliprates = nsdSliprateVector(2)
            sliprated = nsdSliprateVector(3)
            sliprate = nsdSliprateVector(4)
            
            tnrm = nsdTractionVector(1)
            tstk = nsdTractionVector(2)
            tdip = nsdTractionVector(3)
            ttao = nsdTractionVector(4)
          
            fnfault = fric(7,i,ift) !initial forces on the fault node
            fsfault = fric(8,i,ift) + dtau !norm, strike, dip components directly
            fdfault = fric(49,i,ift)
            isn = nsmp(1,i,ift)
            imn = nsmp(2,i,ift)
            mslav = fnms(isn)        
            mmast = fnms(imn)
            mtotl = mslav + mmast
            mtotl = mtotl * arn(i,ift)
            
        ! !-------------------------------------------------------------------!    
            ! fnfault = fric(7,i,ift) !initial forces on the fault node
            ! fsfault = fric(8,i,ift) + dtau !norm, strike, dip components directly
            ! fdfault = fric(49,i,ift)
            ! isn = nsmp(1,i,ift)
            ! imn = nsmp(2,i,ift)
            ! do j=1,2  !1-slave, 2-master
                ! do k=1,3  !1-x comp, 2-y comp, 3-z comp
                ! !          fvd(k,j,1) = brhs(id(k,nsmp(j,i)))  !1-force
                    ! fvd(k,j,1) = brhs(id1(locid(nsmp(j,i,ift))+k))  !1-force !DL 
                    ! fvd(k,j,2) = v(k,nsmp(j,i,ift)) !2-vel
                    ! fvd(k,j,3) = d(k,nsmp(j,i,ift)) !3-di,iftsp
                    ! !fvd(k,j,3) = d(k,nsmp(j,i,ift)) + rdampk(1)*fvd(k,j,2) !3-di,iftsp
                ! enddo
            ! enddo
            ! !...resolve x,y,z components onto normal, strike and dip components.
            ! !   B.D. 1/26/07
            ! do j=1,3    !1-force,2-vel,3-disp
                ! do k=1,2  !1-slave,2-master
                    ! fvd(4,k,j) = fvd(1,k,j)*un(1,i,ift) + fvd(2,k,j)*un(2,i,ift) + fvd(3,k,j)*un(3,i,ift)  !4-norm
                    ! fvd(5,k,j) = fvd(1,k,j)*us(1,i,ift) + fvd(2,k,j)*us(2,i,ift) + fvd(3,k,j)*us(3,i,ift)  !5-strike
                    ! fvd(6,k,j) = fvd(1,k,j)*ud(1,i,ift) + fvd(2,k,j)*ud(2,i,ift) + fvd(3,k,j)*ud(3,i,ift)  !6-dip
                ! enddo
            ! enddo
            ! slipn = fvd(4,2,3) - fvd(4,1,3)
            ! slips = fvd(5,2,3) - fvd(5,1,3)
            ! slipd = fvd(6,2,3) - fvd(6,1,3)
            ! slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
            ! fric(71,i,ift) = slips  
            ! fric(72,i,ift) = slipd
            ! fric(73,i,ift) = slipn  
            ! slipraten = fvd(4,2,2) - fvd(4,1,2)
            ! sliprates = fvd(5,2,2) - fvd(5,1,2)
            ! sliprated = fvd(6,2,2) - fvd(6,1,2)
            ! fric(74,i,ift) = sliprates  
            ! fric(75,i,ift) = sliprated
            ! sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
            ! if (sliprate>fric(76,i,ift)) then 
                ! fric(76,i,ift)=sliprate
            ! endif
            ! fric(77,i,ift) = fric(77,i,ift) + sliprate * dt
        
            ! mslav = fnms(isn)        
            ! mmast = fnms(imn)
            ! mtotl = mslav + mmast
            ! mtotl = mtotl * arn(i,ift)
            

            ! if (C_elastic==0) then!Plastic     
                ! tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
                    ! + mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl          
                ! tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
                    ! - mmast*fvd(5,1,1)) / mtotl
                ! tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
                    ! - mmast*fvd(6,1,1)) / mtotl
            ! else!Elastic
                ! tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
                    ! + mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl + fnfault         
                ! tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
                    ! - mmast*fvd(5,1,1)) / mtotl + fsfault
                ! tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
                    ! - mmast*fvd(6,1,1)) / mtotl + fdfault
            ! endif        
            ! ttao = sqrt(tstk*tstk + tdip*tdip)   
        
            if (friclaw==1 .or. friclaw==2)then!Differ 1&2 and 3&4    
                if(friclaw == 1) then
                    call slip_weak(fric(77,i,ift),fric(1,i,ift),xmu)
                elseif(friclaw == 2) then
                    trupt =  time - fnft(i,ift)
                    call time_weak(trupt,fric(1,i,ift),xmu)
                endif
                
                ! Artificial nucleation 
                if (TPV == 201 .or. TPV == 202) then 
                    if (C_nuclea == 1 .and.ift == nucfault) then
                        call nucleation(dtau, xmu, x(1,nsmp(1,i,ift)), x(2,nsmp(1,i,ift)), & 
                                    x(3,nsmp(1,i,ift)), fric(5,i,ift), fric(1,i,ift), &
                                    fric(2,i,ift))
                    endif
                endif 

                if((tnrm+fric(6,i,ift))>0) then
                    tnrm0 = 0.0d0
                else
                    tnrm0 = tnrm+fric(6,i,ift)
                endif
                taoc = fric(4,i,ift) - xmu *tnrm0

                if(ttao > taoc) then
                    tstk = tstk * taoc / ttao
                    tdip = tdip * taoc / ttao
                    if(fnft(i,ift)>600) then    !fnft should be initialized by >10000
                        if(sliprate >= 0.001d0 .and. mode==1) then    !first time to reach 1mm/s
                            fnft(i,ift) = time    !rupture time for the node
                        elseif (sliprate >=0.05d0 .and. mode==2) then
                            fnft(i,ift) = time
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
                
                ! modify normal stress to effective normal stress.
                if (friclaw == 5) then 
                    ! for friclaw==5, where thermopressurization is applied, 
                    ! needs to add presure from fric(51).
                    tnrm = tnrm + fric(51,i,ift) ! consider termopressurization.
                else
                    ! otherwise, update normal stress with assigned pore pressure fric(6).
                    tnrm = tnrm + fric(6,i,ift)
                endif 
                
                ! If non-planar fault geometry and elastic material, enforce normal stress caps.
                if (rough_fault == 1 .and. C_elastic == 1) then
                    !tnrm = min(min_norm, tnrm) ! Maintain a minimum normal stress level.
                    max_norm      = -40.0d6
                    min_norm      = -10.0d6
                
                    if (tnrm>=min_norm) then 
                        tnrm = min_norm
                    elseif (tnrm<=max_norm) then
                        tnrm = max_norm
                    endif
                endif 
                
                ! avoid positive effective normal stress.
                if (tnrm > 0.0d0) tnrm = 0.0d0
                
                ! Add the background slip rate on top. 
                slipn = slipn + fric(25,i,ift) * time 
                slips = slips + fric(26,i,ift) * time
                slipd = slipd + fric(27,i,ift) * time
                slip = sqrt(slips**2 + slipd**2) !slip mag
                slipraten =  slipraten + fric(25,i,ift) 
                sliprates =  sliprates + fric(26,i,ift)
                sliprated =  sliprated + fric(27,i,ift)
                sliprate = sqrt(sliprates**2+sliprated**2)
                
                    
                if(fnft(i,ift)>600.0d0) then    !fnft should be initialized by >10000
                    if(sliprate >= 0.001d0 .and. mode==1) then    !first time to reach 1mm/s
                        fnft(i,ift) = time    !rupture time for the node
                    elseif (sliprate>=0.05d0 .and. mode==2) then
                        fnft(i,ift) = time
                    endif
                endif
                v_trial = sliprate
                
                ! retrieve the state variable for normal stress theta_pc_tmp from fric(23).
                ! this accounts for normal stress change. 
                theta_pc_tmp = fric(23,i,ift)
                ! get updated trial state variable for normal stress [fric(23)] and its rate [fric(24)].
                call rate_state_normal_stress(v_trial, fric(23,i,ift), theta_pc_dot, tnrm, fric(1,i,ift))    
                fric(24,i,ift) = theta_pc_dot
                
                
                mr       = mmast * mslav / (mmast+mslav) !reduced mass   
                T_coeff  = arn(i,ift)* dt / mr
                
                ! retrieve the RSF state variable and assign it to a tempraroy statetmp. 
                statetmp = fric(20,i,ift)  
                ! get updated trial RSF state variable [fric(20)],
                !   and trial friction coefficient, xmu,
                !   and trial derivative d(xmu)/dt, dxmudv,
                !   for friclaw=3,4,5.
                if(friclaw == 3) then
                    call rate_state_ageing_law(v_trial,fric(20,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                elseif (friclaw == 4 .or. friclaw==5) then
                    call rate_state_slip_law(v_trial,fric(20,i,ift),fric(1,i,ift),xmu,dxmudv) !RSF
                endif            
                
                
                ! compute trial traction.
                ! for cases with large fluctuations of effective normal stress, 
                !   use the state variable for effective normal stress, theta_pc_tmp, 
                !   rather than tnrm, when friclaw==5/rough_fault==1.
                ! [NOTE]: friclaw=5 doesn't support normal stress evolution yet. See TPV1053D.
                if (friclaw==5) then 
                    taoc_old = fric(4,i,ift) - xmu * tnrm
                else
                    taoc_old = xmu * theta_pc_tmp
                endif
                
                tstk0    = tstk
                tdip0    = tdip
                ! get shear tractions, tstk1 and tdip1, and total shear traction, ttao1, updated.
                tstk1    = tstk0 - taoc_old*0.5d0 * (sliprates / sliprate) + fric(26,i,ift)/T_coeff
                tdip1    = tdip0 - taoc_old*0.5d0 * (sliprated / sliprate) + fric(27,i,ift)/T_coeff
                ttao1    = sqrt(tstk1*tstk1 + tdip1*tdip1)
              
                ! Netwon solver for slip rate, v_trial, for the next time step.
                ivmax    = 20  ! Maximum iterations.
                do iv = 1,ivmax
                    ! in each iteration, reupdate the new state variable [fric(20)] given the new 
                    !   slip rate, v_trial.
                    fric(20,i,ift)  = statetmp
                    if(friclaw == 3) then
                        call rate_state_ageing_law(v_trial,fric(20,i,ift),fric(1,i,ift),xmu,dxmudv)
                    else
                        call rate_state_slip_law(v_trial,fric(20,i,ift),fric(1,i,ift),xmu,dxmudv)
                    endif 
                    
                    ! [NOTE]: the code doesn't support normal stress evolution under thermo pressurization.
                    if (friclaw < 5) then
                        fric(23,i,ift)  = theta_pc_tmp 
                        call rate_state_normal_stress(v_trial, fric(23,i,ift), theta_pc_dot, tnrm, fric(1,i,ift))    
                        taoc_new        = xmu*theta_pc_tmp
                        rsfeq           = v_trial + T_coeff * (taoc_new*0.5d0 - ttao1)
                        drsfeqdv        = 1.0d0 + T_coeff * (dxmudv * theta_pc_tmp)*0.5d0  
                    else
                        taoc_new        = fric(4,i,ift) - xmu * MIN(tnrm, 0.0d0)
                        rsfeq           = v_trial + T_coeff * (taoc_new*0.5d0 - ttao1)
                        drsfeqdv        = 1.0d0 + T_coeff * (-dxmudv * MIN(tnrm,0.0d0))*0.5d0  
                    endif
                    
                    ! exiting criteria:
                    !   1. relative residual, rsfeq/drsfeqdv, is smaller than 1e-14*v_trial
                    !   2. residual, rsfeq, is smaller than 1e-6*v_trial
                    if(abs(rsfeq/drsfeqdv) < 1.d-14 * abs(v_trial) .and. abs(rsfeq) < 1.d-6 * abs(v_trial)) exit 
                    !if(abs(rsfeq) < 1.d-5 * abs(v_trial)) exit 
                        vtmp = v_trial - rsfeq / drsfeqdv
                    
                    ! additional constraints for solving trial slip rate, v_trial
                    !   if vtmp smaller than zero, reset it to half of v_trial in the last try. 
                    if(vtmp <= 0.0d0) then
                        v_trial = v_trial/2.0d0
                    else
                        v_trial = vtmp
                    endif  
                    
                enddo !iv
                
                ! If cannot find a solution for v_trial, manually set it to a small value, typically the creeping rate.
                ! Also reset taoc_new to 2 X ttao1.
                ! Without this, TPV1053D blew up at the surface station (-4.2,0)
                if(v_trial < fric(46,i,ift)) then
                    v_trial  = fric(46,i,ift)
                    taoc_new = ttao1*2.0d0
                endif
                
                tstk = taoc_old*0.5d0 * (sliprates / sliprate) + taoc_new*0.5d0 * (tstk1 / ttao1) 
                tdip = taoc_old*0.5d0 * (sliprated / sliprate) + taoc_new*0.5d0 * (tdip1 / ttao1) 
                
                ! store tnrm, tstk, tdip ... 
                ! [effective normal stress, shear_strike, and shear_dip]
                fric(78,i,ift) = tnrm 
                fric(79,i,ift) = tstk
                fric(80,i,ift) = tdip
                
                ! store final slip rate and final total traction ...
                fric(47,i,ift) = v_trial
                fric(48,i,ift) = (tstk**2 + tdip**2)**0.5 
                
                frichis(1,i,nt,ift) = fric(47,i,ift)
                frichis(2,i,nt,ift) = fric(48,i,ift)
                
                ! 3 components of relative acceleration bewteen m-s nodes in the fault plane coordinate sys. 
                accn = -slipraten/dt - slipn/dt/dt
                accs = (v_trial * (tstk1 / ttao1) - sliprates)/dt
                accd = (v_trial * (tdip1 / ttao1) - sliprated)/dt
                
                ! 3 components of relative acceleration bewteen m-s nodes in the FEM xyz coordinate sys. 
                accx = accn*un(1,i,ift) + accs*us(1,i,ift) + accd*ud(1,i,ift)
                accy = accn*un(2,i,ift) + accs*us(2,i,ift) + accd*ud(2,i,ift)
                accz = accn*un(3,i,ift) + accs*us(3,i,ift) + accd*ud(3,i,ift)
                
                ! determine total forces acting on the node pair ...
                Rx = brhs(id1(locid(isn)+1)) + brhs(id1(locid(imn)+1))
                Ry = brhs(id1(locid(isn)+2)) + brhs(id1(locid(imn)+2))
                Rz = brhs(id1(locid(isn)+3)) + brhs(id1(locid(imn)+3))
                
                ! calculate xyz components of nodal forces that can generate 
                !  the above calculated accelerations for the m-s node pair. 
                brhs(id1(locid(isn)+1)) = (-accx + Rx/mmast) * mr ! asx
                brhs(id1(locid(isn)+2)) = (-accy + Ry/mmast) * mr ! asy
                brhs(id1(locid(isn)+3)) = (-accz + Rz/mmast) * mr ! asz
                brhs(id1(locid(imn)+1)) = (accx  + Rx/mslav) * mr ! amx
                brhs(id1(locid(imn)+2)) = (accy  + Ry/mslav) * mr ! amy
                brhs(id1(locid(imn)+3)) = (accz  + Rz/mslav) * mr ! amz
                
                ! store normal velocities for master-slave node pair ...
                ! v(k,nsmp(j,i,ift)) - k:xyz, j:slave1,master2
                fric(31,i,ift) = v(1,nsmp(2,i,ift)) + (accx+Rx/mslav)*dt  !vmx 
                fric(32,i,ift) = v(2,nsmp(2,i,ift)) + (accy+Ry/mslav)*dt  !vmy 
                fric(33,i,ift) = v(3,nsmp(2,i,ift)) + (accz+Rz/mslav)*dt  !vmz 
                fric(34,i,ift) = v(1,nsmp(1,i,ift)) + (-accx+Rx/mmast)*dt !vsx 
                fric(35,i,ift) = v(2,nsmp(1,i,ift)) + (-accy+Ry/mmast)*dt !vsy 
                fric(36,i,ift) = v(3,nsmp(1,i,ift)) + (-accz+Rz/mmast)*dt !vsz 
            endif

            if(n4onf>0.and.lstr) then    
                do j=1,n4onf
                    if(anonfs(1,j)==i.and.anonfs(3,j)==ift) then !only selected stations. B.D. 10/25/09    
                        fltsta(1,locplt-1,j)  = time
                        fltsta(2,locplt-1,j)  = sliprates
                        fltsta(3,locplt-1,j)  = sliprated
                        fltsta(4,locplt-1,j)  = fric(20,i,ift)
                        fltsta(5,locplt-1,j)  = slips
                        fltsta(6,locplt-1,j)  = slipd
                        fltsta(7,locplt-1,j)  = slipn
                        fltsta(8,locplt-1,j)  = tstk
                        fltsta(9,locplt-1,j)  = tdip
                        fltsta(10,locplt-1,j) = tnrm
                        fltsta(11,locplt-1,j) = fric(51,i,ift) + fric(42,i,ift) ! + fric_tp_pini
                        fltsta(12,locplt-1,j) = fric(52,i,ift) 
                    endif
                enddo 
            endif   
            
        enddo    !ending i
    enddo !ift
    !-------------------------------------------------------------------!
    !-------------Late Sep.2015/ D.Liu----------------------------------!
    !-----------Writing out results on fault for evert nstep------------!
    !if(mod(nt,315)==1.and.nt<5000) then 
    !    write(mm,'(i6)') me
    !    mm = trim(adjustl(mm))
    !    foutmov='fslipout_'//mm
    !    open(9002+me,file=foutmov,form='formatted',status='unknown',position='append')
    !        write(9002+me,'(1x,4f10.3)') ((fltslp(j,ifout),j=1,3),fltslr(1,ifout),ifout=1,nftnd)
    !endif
    !----nftnd for each me for plotting---------------------------------!
    !if (nt==1) then
    !    write(mm,'(i6)') me    
    !    mm = trim(adjustl(mm))            
    !    foutmov='fnode.txt'//mm
    !    open(unit=9800,file=foutmov,form='formatted',status='unknown')
    !        write(9800,'(2I7)') me,nftnd 
    !    close(9800)            
    !endif     
    !-------------------------------------------------------------------!    
end subroutine faulting     

! Subroutine rate_state_normal_stress calculates the effect of normal stress change
! from RSF. The formulation follows Shi and Day (2013), eq B8. {"Frictional sliding experiments with variable normal stress show that the shear strength responds gradually to abrupt changes of normal stress (e.g., Prakash and Clifton, 1993; Prakash, 1998)."}

! theta_pc_dot = - V/L_pc*[theta_pc - abs(tnrm)]

! Input: slip_rate, L_pc, theta_pc, tnrm. 
! Output: theta_pc, the state variable which is used to calculate shear stress in eq B2
! B2: abs(traction) = friction * theta_pc.
subroutine rate_state_normal_stress(V2, theta_pc, theta_pc_dot, tnrm, fricsgl)
    use globalvar
    implicit none
    real (kind = dp) :: V2, theta_pc, theta_pc_dot, tnrm, L
    real (kind = dp),dimension(100) :: fricsgl
    
    L  = fricsgl(11) ! Use Dc in RSF as L_pc
    
    theta_pc_dot = - V2/L*(theta_pc - abs(tnrm))
    ! the following eq is to update theta_pc with theta_pc_doc.
    ! this is now consistent with EQquasi.
    theta_pc = theta_pc + theta_pc_dot*dt
    
    ! the following eq, which is not used, directly writes out the analytic solution
    ! of the OED.
    ! theta_pc = abs(tnrm) + (theta_pc - abs(tnrm))*dexp(-V2*dt/L)
    
end subroutine rate_state_normal_stress

subroutine nucleation(dtau, xmu, xx, yy, zz, twt0, fs, fd)
    ! Subroutine nucleation handles the artificial nucleation for 
    !   various friction laws.
    ! It will return friction coefficient xmu or dtau as a function of time 
    !   and fault node locations.
    
    ! nucR, dtau0, nucRuptVel are loaded from input file bGlobal.txt, which is 
    !   generated from user_defined_param.py.
    use globalvar
    implicit none
    real(kind = dp) :: T, F, G, rr, dtau, xmu, xx, yy, zz
    real(kind = dp) :: tr, tc, tmp1, tmp2, twt0, fs, fd
    
    dtau = 0.0d0 
    
    if (TPV == 105 .or. TPV == 104) then
        T  = 1.0d0
        F  = 0.0d0
        G  = 1.0d0
        rr = sqrt((xx-xsource)**2 + (yy-ysource)**2 + (zz-zsource)**2)
        
        if (rr < nucR) then 
            F=dexp(rr**2/(rr**2-nucR**2))
        endif 

        if (time<=T)  then 
            G=dexp((time-T)**2/(time*(time-2*T)))
        endif 
    
        dtau = nucdtau0*F*G
    
    elseif (TPV == 201 .or. TPV == 202) then
        rr = sqrt((xx-xsource)**2 + (yy-ysource)**2 + (zz-zsource)**2)
        if(rr <= nucR) then 
            if (TPV == 201) then 
                tr = (rr+0.081d0*nucR*(1.0d0/(1.0d0-(rr/nucR)**2)-1.0d0))/(0.7d0*3464.d0)
            elseif (TPV == 202) then 
                tr = rr/nucRuptVel
            endif 
        else
            tr = 1.0d9 
        endif
        
        if(time<tr) then 
            tc = 0.0d0
        elseif ((time<(tr+twt0)).and.(time>=tr)) then 
            tc = (time-tr)/twt0
        else 
            tc = 1.0d0
        endif
        
        tmp1 = fs+(fd-fs)*tc 
        tmp2 = xmu
        xmu  = min(tmp1,tmp2)  
    else
        write(*,*) "Artificial nucleation mode is not supported yet"
        write(*,*) "Exiting ... ..."
        stop
    endif 
end subroutine nucleation

subroutine getNsdSlipSliprateTraction(iFault, iFaultNodePair, nsdSlipVector, nsdSliprateVector, nsdTractionVector, dtau)
! get slip, sliprate, and traction vectors on one pair of fault split-nodes
    use globalvar
    implicit none
    integer(kind = 4) :: iFault, iFaultNodePair, iSlaveNodeID, iMasterNodeID
    integer(kind = 4) :: j, k
    real(kind = dp) :: initNormal, initStrikeShear, initDipShear
    real(kind = dp) :: xyzNodalQuant(3,2,3), nsdNodalQuant(3,2,3)
    real(kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4), nsdTractionVector(4)
    real(kind = dp) :: massSlave, massMaster, totalMass
    
    real(kind = dp) :: dtau
    
    initNormal      = fric(7, iFaultNodePair, iFault)
    initStrikeShear = fric(8, iFaultNodePair, iFault)+dtau
    initDipShear    = fric(49, iFaultNodePair, iFault)
    
    iSlaveNodeID  = nsmp(1, iFaultNodePair, iFault)
    iMasterNodeID = nsmp(2, iFaultNodePair, iFault)
    massSlave     = fnms(iSlaveNodeID)        
    massMaster    = fnms(iMasterNodeID)
    totalMass     = (massSlave + massMaster)*arn(iFaultNodePair, iFault)
    
    do j=1,2  ! slave, master
        do k=1,3  ! x,y,z
            xyzNodalQuant(k,j,1) = brhs(id1(locid(nsmp(j,iFaultNodePair,iFault))+k))  !1-force !DL 
            xyzNodalQuant(k,j,2) = v(k,nsmp(j,iFaultNodePair,iFault)) !2-vel
            xyzNodalQuant(k,j,3) = d(k,nsmp(j,iFaultNodePair,iFault)) !3-di,iftsp
        enddo
    enddo
    
    do j=1,3    !1-force,2-vel,3-disp
        do k=1,2  !1-slave,2-master
            nsdNodalQuant(1,k,j) = xyzNodalQuant(1,k,j)*un(1,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(2,k,j)*un(2,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(3,k,j)*un(3,iFaultNodePair, iFault)  !norm
            nsdNodalQuant(2,k,j) = xyzNodalQuant(1,k,j)*us(1,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(2,k,j)*us(2,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(3,k,j)*us(3,iFaultNodePair, iFault)  !strike
            nsdNodalQuant(3,k,j) = xyzNodalQuant(1,k,j)*ud(1,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(2,k,j)*ud(2,iFaultNodePair, iFault) &
                                    + xyzNodalQuant(3,k,j)*ud(3,iFaultNodePair, iFault)  !dip
        enddo
    enddo
    
    do j=1,3 !n,s,d
        nsdSlipVector(j) = nsdNodalQuant(j,2,3) - nsdNodalQuant(j,1,3)
    enddo
    nsdSlipVector(4) = sqrt(nsdSlipVector(1)**2 + nsdSlipVector(2)**2 + nsdSlipVector(3)**2)
    
    do j=1,3 !n,s,d
        nsdSliprateVector(j) = nsdNodalQuant(j,2,2) - nsdNodalQuant(j,1,2)
    enddo
    nsdSliprateVector(4) = sqrt(nsdSliprateVector(1)**2 + nsdSliprateVector(2)**2 + nsdSliprateVector(3)**2)
    
    ! keep records
    fric(71, iFaultNodePair, iFault) = nsdSlipVector(2) !s
    fric(72, iFaultNodePair, iFault) = nsdSlipVector(3) !d
    fric(73, iFaultNodePair, iFault) = nsdSlipVector(1) !n
    fric(74, iFaultNodePair, iFault) = nsdSliprateVector(2) !s
    fric(75, iFaultNodePair, iFault) = nsdSliprateVector(3) !d
    if (nsdSliprateVector(4)>fric(76, iFaultNodePair, iFault)) fric(76, iFaultNodePair, iFault) = nsdSliprateVector(4) !mag
    fric(77, iFaultNodePair, iFault) = fric(77, iFaultNodePair, iFault) + nsdSliprateVector(4)*dt ! cummulated slip
    
    ! n
    nsdTractionVector(1) = (massSlave*massMaster*((nsdNodalQuant(1,2,2)-nsdNodalQuant(1,1,2))+(nsdNodalQuant(1,2,3)-nsdNodalQuant(1,1,3))/dt)/dt &
                            + massSlave*nsdNodalQuant(1,2,1) - massMaster*nsdNodalQuant(1,1,1))/totalMass &
                            + initNormal*C_elastic         
    ! s
    nsdTractionVector(2) = (massSlave*massMaster*(nsdNodalQuant(2,2,2)-nsdNodalQuant(2,1,2))/dt &
                            + massSlave*nsdNodalQuant(2,2,1) - massMaster*nsdNodalQuant(2,1,1))/totalMass &
                            + initStrikeShear*C_elastic
    ! d
    nsdTractionVector(3) = (massSlave*massMaster*(nsdNodalQuant(3,2,2)-nsdNodalQuant(3,1,2))/dt &
                            + massSlave*nsdNodalQuant(3,2,1) - massMaster*nsdNodalQuant(3,1,1)) /totalMass &
                            + initDipShear*C_elastic
    ! shear traction magnitude
    nsdTractionVector(4) = sqrt(nsdTractionVector(1)**2+nsdTractionVector(2)**2)
    
     !-------------------------------------------------------------------!    
            ! fnfault = fric(7,i,ift) !initial forces on the fault node
            ! fsfault = fric(8,i,ift) + dtau !norm, strike, dip components directly
            ! fdfault = fric(49,i,ift)
            ! isn = nsmp(1,i,ift)
            ! imn = nsmp(2,i,ift)
            ! do j=1,2  !1-slave, 2-master
                ! do k=1,3  !1-x comp, 2-y comp, 3-z comp
                ! !          fvd(k,j,1) = brhs(id(k,nsmp(j,i)))  !1-force
                    ! fvd(k,j,1) = brhs(id1(locid(nsmp(j,i,ift))+k))  !1-force !DL 
                    ! fvd(k,j,2) = v(k,nsmp(j,i,ift)) !2-vel
                    ! fvd(k,j,3) = d(k,nsmp(j,i,ift)) !3-di,iftsp
                    ! !fvd(k,j,3) = d(k,nsmp(j,i,ift)) + rdampk(1)*fvd(k,j,2) !3-di,iftsp
                ! enddo
            ! enddo
            ! !...resolve x,y,z components onto normal, strike and dip components.
            ! !   B.D. 1/26/07
            ! do j=1,3    !1-force,2-vel,3-disp
                ! do k=1,2  !1-slave,2-master
                    ! fvd(4,k,j) = fvd(1,k,j)*un(1,i,ift) + fvd(2,k,j)*un(2,i,ift) + fvd(3,k,j)*un(3,i,ift)  !4-norm
                    ! fvd(5,k,j) = fvd(1,k,j)*us(1,i,ift) + fvd(2,k,j)*us(2,i,ift) + fvd(3,k,j)*us(3,i,ift)  !5-strike
                    ! fvd(6,k,j) = fvd(1,k,j)*ud(1,i,ift) + fvd(2,k,j)*ud(2,i,ift) + fvd(3,k,j)*ud(3,i,ift)  !6-dip
                ! enddo
            ! enddo
            ! slipn = fvd(4,2,3) - fvd(4,1,3)
            ! slips = fvd(5,2,3) - fvd(5,1,3)
            ! slipd = fvd(6,2,3) - fvd(6,1,3)
            ! slip = sqrt(slipn**2 + slips**2 + slipd**2) !slip mag
            ! fric(71,i,ift) = slips  
            ! fric(72,i,ift) = slipd
            ! fric(73,i,ift) = slipn  
            ! slipraten = fvd(4,2,2) - fvd(4,1,2)
            ! sliprates = fvd(5,2,2) - fvd(5,1,2)
            ! sliprated = fvd(6,2,2) - fvd(6,1,2)
            ! fric(74,i,ift) = sliprates  
            ! fric(75,i,ift) = sliprated
            ! sliprate = sqrt(slipraten**2+sliprates**2+sliprated**2)
            ! if (sliprate>fric(76,i,ift)) then 
                ! fric(76,i,ift)=sliprate
            ! endif
            ! fric(77,i,ift) = fric(77,i,ift) + sliprate * dt
        
            ! mslav = fnms(isn)        
            ! mmast = fnms(imn)
            ! mtotl = mslav + mmast
            ! mtotl = mtotl * arn(i,ift)
            

            ! if (C_elastic==0) then!Plastic     
                ! tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
                    ! + mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl          
                ! tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
                    ! - mmast*fvd(5,1,1)) / mtotl
                ! tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
                    ! - mmast*fvd(6,1,1)) / mtotl
            ! else!Elastic
                ! tnrm = (mslav*mmast*((fvd(4,2,2)-fvd(4,1,2))+(fvd(4,2,3)-fvd(4,1,3))/dt)/dt &
                    ! + mslav*fvd(4,2,1) - mmast*fvd(4,1,1)) / mtotl + fnfault         
                ! tstk = (mslav*mmast*(fvd(5,2,2)-fvd(5,1,2))/dt + mslav*fvd(5,2,1) &
                    ! - mmast*fvd(5,1,1)) / mtotl + fsfault
                ! tdip = (mslav*mmast*(fvd(6,2,2)-fvd(6,1,2))/dt + mslav*fvd(6,2,1) &
                    ! - mmast*fvd(6,1,1)) / mtotl + fdfault
            ! endif        
            ! ttao = sqrt(tstk*tstk + tdip*tdip)   
end subroutine getNsdSlipSliprateTraction
