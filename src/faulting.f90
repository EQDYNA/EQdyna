! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine faulting
    use globalvar
    implicit none
    include 'mpif.h'
    
    integer (kind = 4) :: i, j, ift
    real (kind = dp) :: dtau=0.0d0
    real (kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4), nsdTractionVector(4)
    real (kind = dp) :: nsdInitTractionVector(3)
    
    time1 = MPI_WTIME()
    do ift = 1, ntotft
        do i = 1, nftnd(ift)   
            call getNsdSlipSliprateTraction(ift, i, nsdSlipVector, nsdSliprateVector, nsdTractionVector, nsdInitTractionVector, dtau)
            if (friclaw<=2) call solveSWTW(ift, i, friclaw, nsdTractionVector, nsdInitTractionVector)
            if (friclaw>=3) call solveRSF(ift, i, friclaw, nsdSlipVector, nsdSliprateVector, nsdTractionVector)
            call showSourceDynamics(ift,i)
            call storeOnFaultStationQuantSCEC(ift, i, nsdSlipVector, nsdSliprateVector, nsdTractionVector)
            call storeRuptureTime(ift, i, nsdSliprateVector)
        enddo 
    enddo  
    timeused(6) = timeused(6) + MPI_WTIME() - time1 
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

subroutine getNsdSlipSliprateTraction(iFault, iFaultNodePair, nsdSlipVector, nsdSliprateVector, nsdTractionVector, nsdInitTractionVector, dtau)
! get slip, sliprate, and traction vectors on one pair of fault split-nodes
    use globalvar
    implicit none
    integer(kind = 4) :: iFault, iFaultNodePair
    integer(kind = 4) :: j, k
    real(kind = dp) :: xyzNodalQuant(3,2,3), nsdNodalQuant(3,2,3)
    real(kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4), nsdTractionVector(4)
    real(kind = dp) :: nsdInitTractionVector(3)
    real(kind = dp) :: massSlave, massMaster, totalMass
    
    real(kind = dp) :: dtau
    
    nsdInitTractionVector(1) = fric(7,iFaultNodePair,iFault) !normal
    nsdInitTractionVector(2) = fric(8,iFaultNodePair,iFault)+dtau !strike
    nsdInitTractionVector(3) = fric(49,iFaultNodePair,iFault) !dip
    
    massSlave     = fnms(nsmp(1,iFaultNodePair,iFault))        
    massMaster    = fnms(nsmp(2,iFaultNodePair,iFault))
    totalMass     = (massSlave + massMaster)*arn(iFaultNodePair,iFault)
    
    do j = 1, 2  ! slave, master
        do k = 1, 3  ! x,y,z
            xyzNodalQuant(k,j,1) = brhs(id1(locid(nsmp(j,iFaultNodePair,iFault))+k))  !1-force !DL 
            xyzNodalQuant(k,j,2) = v(k,nsmp(j,iFaultNodePair,iFault)) !2-vel
            xyzNodalQuant(k,j,3) = d(k,nsmp(j,iFaultNodePair,iFault)) !3-di,iftsp
        enddo
    enddo
    
    do j = 1, 3    !1-force,2-vel,3-disp
        do k = 1, 2  !1-slave,2-master
            nsdNodalQuant(1,k,j) = xyzNodalQuant(1,k,j)*un(1,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(2,k,j)*un(2,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(3,k,j)*un(3,iFaultNodePair,iFault)  !norm
            nsdNodalQuant(2,k,j) = xyzNodalQuant(1,k,j)*us(1,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(2,k,j)*us(2,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(3,k,j)*us(3,iFaultNodePair,iFault)  !strike
            nsdNodalQuant(3,k,j) = xyzNodalQuant(1,k,j)*ud(1,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(2,k,j)*ud(2,iFaultNodePair,iFault) &
                                    + xyzNodalQuant(3,k,j)*ud(3,iFaultNodePair,iFault)  !dip
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
    fric(71,iFaultNodePair,iFault) = nsdSlipVector(2) !s
    fric(72,iFaultNodePair,iFault) = nsdSlipVector(3) !d
    fric(73,iFaultNodePair,iFault) = nsdSlipVector(1) !n
    fric(74,iFaultNodePair,iFault) = nsdSliprateVector(2) !s
    fric(75,iFaultNodePair,iFault) = nsdSliprateVector(3) !d
    if (nsdSliprateVector(4)>fric(76,iFaultNodePair,iFault)) fric(76,iFaultNodePair,iFault) = nsdSliprateVector(4) !mag
    fric(77,iFaultNodePair,iFault) = fric(77,iFaultNodePair,iFault) + nsdSliprateVector(4)*dt ! cummulated slip
    
    ! n
    nsdTractionVector(1) = (massSlave*massMaster*((nsdNodalQuant(1,2,2)-nsdNodalQuant(1,1,2))+(nsdNodalQuant(1,2,3)-nsdNodalQuant(1,1,3))/dt)/dt &
                            + massSlave*nsdNodalQuant(1,2,1) - massMaster*nsdNodalQuant(1,1,1))/totalMass &
                            + nsdInitTractionVector(1)*C_elastic         
    ! s
    nsdTractionVector(2) = (massSlave*massMaster*(nsdNodalQuant(2,2,2)-nsdNodalQuant(2,1,2))/dt &
                            + massSlave*nsdNodalQuant(2,2,1) - massMaster*nsdNodalQuant(2,1,1))/totalMass &
                            + nsdInitTractionVector(2)*C_elastic
    ! d
    nsdTractionVector(3) = (massSlave*massMaster*(nsdNodalQuant(3,2,2)-nsdNodalQuant(3,1,2))/dt &
                            + massSlave*nsdNodalQuant(3,2,1) - massMaster*nsdNodalQuant(3,1,1)) /totalMass &
                            + nsdInitTractionVector(3)*C_elastic
                            
    if (friclaw>=3 .and. C_nuclea==1 .and. iFault==nucfault) call rsfNucleation(iFault, iFaultNodePair, nsdTractionVector, nsdSliprateVector)
    
    ! shear traction magnitude
    nsdTractionVector(4) = sqrt(nsdTractionVector(2)**2+nsdTractionVector(3)**2)
    
end subroutine getNsdSlipSliprateTraction

subroutine solveSWTW(iFault, iFaultNodePair, iFrictionLaw, nsdTractionVector, nsdInitTractionVector)
                
    use globalvar
    implicit none
    
    integer(kind = 4) :: iFault, iFaultNodePair, iFrictionLaw
    integer(kind = 4) :: j
    real(kind = dp) :: fricCoeff, trupt, dtau, effectiveNormalStress, trialShearTraction
    real(kind = dp) :: nsdTractionVector(4), nsdInitTractionVector(3)
    real(kind = dp) :: xyzTractionVector(3), xyzInitTractionVector(3)
    
    
    if (iFrictionLaw == 1) then
        call slip_weak(fric(77,iFaultNodePair,iFault),fric(1,iFaultNodePair,iFault),fricCoeff)
    elseif(iFrictionLaw == 2) then
        trupt =  time - fnft(iFaultNodePair,iFault)
        call time_weak(trupt,fric(1,iFaultNodePair,iFault),fricCoeff)
    endif
    
    if (C_nuclea==1 .and. iFault==nucFault) call swtwNucleation(iFault, iFaultNodePair, fricCoeff)

    if((nsdTractionVector(1)+fric(6,iFaultNodePair,iFault))>0) then
        effectiveNormalStress = 0.0d0
    else
        effectiveNormalStress = nsdTractionVector(1)+fric(6,iFaultNodePair,iFault)
    endif
    trialShearTraction = fric(4,iFaultNodePair,iFault) - fricCoeff*effectiveNormalStress

    if(nsdTractionVector(4) > trialShearTraction) then
        ! adjust strike shear traction 
        nsdTractionVector(2) = nsdTractionVector(2)*trialShearTraction/nsdTractionVector(4)
        ! adjust dip shear traction
        nsdTractionVector(3) = nsdTractionVector(3)*trialShearTraction/nsdTractionVector(4)
    endif
    
    do j=1,3 !x,y,z
        xyzTractionVector(j) = (nsdTractionVector(1)*un(j,iFaultNodePair,iFault) &
                                + nsdTractionVector(2)*us(j,iFaultNodePair,iFault) &
                                + nsdTractionVector(3)*ud(j,iFaultNodePair,iFault)) * arn(iFaultNodePair,iFault)
        xyzInitTractionVector(j) = (nsdInitTractionVector(1)*un(j,iFaultNodePair,iFault) &
                                + nsdInitTractionVector(2)*us(j,iFaultNodePair,iFault) &
                                + nsdInitTractionVector(3)*ud(j,iFaultNodePair,iFault)) * arn(iFaultNodePair,iFault)
    enddo
    
    do j=1,3 !x,y,z
        brhs(id1(locid(nsmp(1,iFaultNodePair,iFault))+j)) = brhs(id1(locid(nsmp(1,iFaultNodePair,iFault))+j)) + xyzTractionVector(j) - xyzInitTractionVector(j)*C_elastic
        brhs(id1(locid(nsmp(2,iFaultNodePair,iFault))+j)) = brhs(id1(locid(nsmp(2,iFaultNodePair,iFault))+j)) - xyzTractionVector(j) + xyzInitTractionVector(j)*C_elastic
    enddo
    
    do j=1,3 !n,s,d
        fric(78+j-1,iFaultNodePair,iFault) = nsdTractionVector(j)
    enddo
end subroutine solveSWTW

subroutine solveRSF(iFault, iFaultNodePair, iFrictionLaw, nsdSlipVector, nsdSliprateVector, nsdTractionVector)
    use globalvar
    implicit none
    integer(kind = 4) :: iFault, iFaultNodePair, iFrictionLaw
    integer(kind = 4) :: j, iv, ivmax
    real(kind = dp) :: effectiveNormalStress, trialSliprate, v_trial, taoc_new, taoc_old, rsfeq, drsfeqdv, xmu, dxmudv, vtmp, massMaster, massSlave, T_coeff, mr, theta_pc_dot, theta_pc_tmp, statetmp
    real(kind = dp) :: nsdTractionVector(4), trialTractVec(4)
    real(kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4)
    real(kind = dp) :: nsdAccVec(3), xyzAccVec(3), xyzR(3)
    
    
    ! adjust normal stress 
    if (iFrictionLaw==5) then
        nsdTractionVector(1) = nsdTractionVector(1) + fric(51,iFaultNodePair,iFault)
    else
        nsdTractionVector(1) = nsdTractionVector(1) + fric(6,iFaultNodePair,iFault)
    endif 

    ! If non-planar fault geometry and elastic material, enforce normal stress caps.
    if (insertFaultType >0 .and. C_elastic == 1) then
        !tnrm = min(min_norm, tnrm) ! Maintain a minimum normal stress level.
        max_norm      = -40.0d6
        min_norm      = -10.0d6
    
        if (nsdTractionVector(1)>=min_norm) then 
            nsdTractionVector(1) = min_norm
        elseif (nsdTractionVector(1)<=max_norm) then
            nsdTractionVector(1) = max_norm
        endif
    endif 

    ! avoid positive effective normal stress.
    if (nsdTractionVector(1) > 0.0d0) nsdTractionVector(1) = 0.0d0

    !-----------------
    ! Add the background slip rate on top. 
    do j=1,3 !n,s,d
        nsdSlipVector(j) = nsdSlipVector(j) + fric(25+j-1,iFaultNodePair,iFault)*time
        nsdSliprateVector(j) = nsdSliprateVector(j) + fric(25+j-1,iFaultNodePair,iFault)
    enddo
    nsdSlipVector(4) = sqrt(nsdSlipVector(2)**2+nsdSlipVector(3)**2)
    nsdSliprateVector(4) = sqrt(nsdSliprateVector(2)**2+nsdSliprateVector(3)**2)
        
    if(fnft(iFaultNodePair,iFault)>600.0d0) then    !fnft should be initialized by >10000
        if(nsdSliprateVector(4) >= 0.001d0 .and. mode==1) then    !first time to reach 1mm/s
            fnft(iFaultNodePair,iFault) = time    !rupture time for the node
        elseif (nsdSliprateVector(4)>=0.05d0 .and. mode==2) then
            fnft(iFaultNodePair,iFault) = time
        endif
    endif
    
    ! Given tractions, state variables, find the sliprate for next time step. 
    v_trial = nsdSliprateVector(4)
    
    ! retrieve the state variable for normal stress thetaPcTmp from fric(23).
    ! this accounts for normal stress change. 
    theta_pc_tmp = fric(23,iFaultNodePair,iFault)
    ! get updated trial state variable for normal stress [fric(23)] and its rate [fric(24)].
    call rate_state_normal_stress(v_trial, fric(23,iFaultNodePair,iFault), theta_pc_dot, nsdTractionVector(1), fric(1,iFaultNodePair,iFault))    
    fric(24,iFaultNodePair,iFault) = theta_pc_dot
    ! retrieve the RSF state variable and assign it to a tempraroy statetmp. 
    statetmp = fric(20,iFaultNodePair,iFault) 
    
    ! get updated trial RSF state variable [fric(20)],
    !   and trial friction coefficient, xmu,
    !   and trial derivative d(xmu)/dt, dxmudv,
    !   for friclaw=3,4,5.
    if(friclaw == 3) then
        call rate_state_ageing_law(v_trial,fric(20,iFaultNodePair,iFault),fric(1,iFaultNodePair,iFault),xmu,dxmudv) !RSF
    elseif (friclaw == 4 .or. friclaw==5) then
        call rate_state_slip_law(v_trial,fric(20,iFaultNodePair,iFault),fric(1,iFaultNodePair,iFault),xmu,dxmudv) !RSF
    endif            
    ! compute trial traction.
    ! for cases with large fluctuations of effective normal stress, 
    !   use the state variable for effective normal stress, thetaPcTmp, 
    !   rather than tnrm, when friclaw==5/insertFaultType>0.
    ! [NOTE]: friclaw=5 doesn't support normal stress evolution yet. See TPV1053D.
    if (friclaw==5) then 
        taoc_old = fric(4,iFaultNodePair,iFault) - xmu * nsdTractionVector(1)
    else
        taoc_old = xmu * theta_pc_tmp
    endif
    massSlave  = fnms(nsmp(1,iFaultNodePair,iFault))
    massMaster = fnms(nsmp(2,iFaultNodePair,iFault))
    mr = massMaster*massSlave/(massMaster+massSlave) !reduced mass
    T_coeff = arn(iFaultNodePair,iFault)*dt/mr
   
    ! get shear tractions, tstk1 and tdip1, and total shear traction, ttao1, updated.
    do j=2,3 ! s,d/ tstk1, tdip1
        trialTractVec(j) = nsdTractionVector(j) - taoc_old*0.5d0*(nsdSliprateVector(j)/nsdSliprateVector(4)) + fric(25+j-1,iFaultNodePair,iFault)/T_coeff
    enddo
    trialTractVec(4) = sqrt(trialTractVec(2)**2 + trialTractVec(3)**2) !ttao1
    
    ! solve sliprate with Newton-Raphson method
    call NewtonRaphson(iFault, iFaultNodePair, v_trial, taoc_new, statetmp, theta_pc_tmp, theta_pc_dot, nsdTractionVector, trialTractVec, T_coeff)
    fric(20,iFaultNodePair,iFault) = statetmp
    fric(23,iFaultNodePair,iFault) = theta_pc_tmp
    do j=2,3 !s,d
        nsdTractionVector(j) = taoc_old*0.5d0*(nsdSliprateVector(j)/nsdSliprateVector(4)) + taoc_new*0.5d0*(trialTractVec(j)/trialTractVec(4))
    enddo
        
    ! store tnrm, tstk, tdip ... 
    ! [effective normal stress, shear_strike, and shear_dip]
    do j=1,3 !n,s,d
        fric(78+j-1,iFaultNodePair,iFault) = nsdTractionVector(j)
    enddo
    ! store final slip rate and final total traction ...
    fric(47,iFaultNodePair,iFault) = v_trial
    fric(48,iFaultNodePair,iFault) = sqrt(nsdTractionVector(2)**2 + nsdTractionVector(3)**2) 
    
    frichis(1,iFaultNodePair,nt,iFault) = fric(47,iFaultNodePair,iFault)
    frichis(2,iFaultNodePair,nt,iFault) = fric(48,iFaultNodePair,iFault)
    
    ! 3 components of relative acceleration bewteen m-s nodes in the fault plane coordinate sys. 
    nsdAccVec(1) = -nsdSliprateVector(1)/dt - nsdSlipVector(1)/dt/dt
    nsdAccVec(2) = (v_trial * (trialTractVec(2)/trialTractVec(4)) - nsdSliprateVector(2))/dt
    nsdAccVec(3) = (v_trial * (trialTractVec(3)/trialTractVec(4)) - nsdSliprateVector(3))/dt
    
    ! 3 components of relative acceleration bewteen m-s nodes in the FEM xyz coordinate sys. 
    do j=1,3 !x,y,z
        xyzAccVec(j) = nsdAccVec(1)*un(j,iFaultNodePair,iFault) + nsdAccVec(2)*us(j,iFaultNodePair,iFault) + nsdAccVec(3)*ud(j,iFaultNodePair,iFault)
        ! determine total forces acting on the node pair ...
        xyzR(j) = brhs(id1(locid(nsmp(1,iFaultNodePair,iFault))+j)) + brhs(id1(locid(nsmp(2,iFaultNodePair,iFault))+j))
    enddo 
    
    do j=1,3 !x,y,z
        ! calculate xyz components of nodal forces that can generate 
        !  the above calculated accelerations for the m-s node pair. 
        brhs(id1(locid(nsmp(1,iFaultNodePair,iFault))+j)) = (-xyzAccVec(j) + xyzR(j)/massMaster)*mr ! Acc.Slave.x/y/z
        brhs(id1(locid(nsmp(2,iFaultNodePair,iFault))+j)) = (xyzAccVec(j)  + xyzR(j)/massSlave)*mr ! Acc.Master.x/y/z
        ! store normal velocities for master-slave node pair ...
        ! v(k,nsmp(j,iFaultNodePair,iFault)) - k:xyz, j:slave1,master2
        fric(31+j-1,iFaultNodePair,iFault) = v(j,nsmp(2,iFaultNodePair,iFault)) + (xyzAccVec(j)+xyzR(j)/massSlave)*dt !Velocity.Master.x/y/z
        fric(34+j-1,iFaultNodePair,iFault) = v(j,nsmp(1,iFaultNodePair,iFault)) + (-xyzAccVec(j)+xyzR(j)/massMaster)*dt !Velocity.Slave.x/y/z
    enddo
    
end subroutine solveRSF

subroutine storeRuptureTime(iFault, iFaultNodePair, nsdSliprateVector)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault, iFaultNodePair
    real (kind = dp) :: nsdSliprateVector(4)
    
    if(fnft(iFaultNodePair,iFault)>5000.0d0) then    !fnft is initialized to be 99999
        if(nsdSliprateVector(4) >= 0.001d0 .and. mode==1) then    !first time to reach 1mm/s
            fnft(iFaultNodePair,iFault) = time    !rupture time for the node
        elseif (nsdSliprateVector(4) >=0.05d0 .and. mode==2) then
            fnft(iFaultNodePair,iFault) = time
        endif
    endif  

end subroutine storeRuptureTime

subroutine showSourceDynamics(iFault,iFaultNodePair)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault, iFaultNodePair
    
    if (abs(x(1,nsmp(1,iFaultNodePair,iFault))-xsource)<tol .and. &
        abs(x(3,nsmp(1,iFaultNodePair,iFault))-zsource)<tol) then
        write(*,*) 'Hypocenter dynamics'
        write(*,'(X,A,3E15.7)') 'X, Y, Z      (m)   ', x(1,nsmp(1,iFaultNodePair,iFault)), x(2,nsmp(1,iFaultNodePair,iFault)), x(3,nsmp(1,iFaultNodePair,iFault))
        write(*,'(X,A,E15.7)') 'Time         (s)   ', time
        write(*,'(X,A,3E15.7)') 'n,s,d tract  (MPa) ', fric(78,iFaultNodePair,iFault)/1.d6, fric(79,iFaultNodePair,iFault)/1.d6, fric(80,iFaultNodePair,iFault)/1.d6
        write(*,'(X,A,E15.7)') 'state_normal (MPa) ', fric(23,iFaultNodePair,iFault)/1.d6
        write(*,'(X,A,3E15.7)') 'n,s,d slip   (m)   ', fric(73,iFaultNodePair,iFault), fric(71,iFaultNodePair,iFault), fric(72,iFaultNodePair,iFault)
        write(*,'(X,A,3E15.7)') 's,d, peak sr (m/s) ', fric(74,iFaultNodePair,iFault), fric(75,iFaultNodePair,iFault), fric(76,iFaultNodePair,iFault)
        write(*,'(X,A,E15.7)') 'cummul slip  (m)   ', fric(77,iFaultNodePair,iFault)
        write(*,'(X,A,3E15.7)') 'sw_fs, sw_fd, sw_D0', fric(1,iFaultNodePair,iFault), fric(2,iFaultNodePair,iFault), fric(3,iFaultNodePair,iFault)
        write(*,'(X,A,2E15.7)') 'rsf_a, rsf_b       ', fric(9,iFaultNodePair,iFault), fric(10,iFaultNodePair,iFault)
        write(*,'(X,A,E15.7)') 'rsf_state          ', fric(20,iFaultNodePair,iFault)
        write(*,'(X,A,E15.7)') 'Nuc add tau0 (MPa) ', fric(81,iFaultNodePair,iFault)/1.d6
    endif 
    
end subroutine showSourceDynamics

subroutine rsfNucleation(iFault, iFaultNodePair, nsdTractionVector, nsdSliprateVector)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault, iFaultNodePair, j
    real (kind = dp) :: dtau, radius, T, F, G, nsdTractionVector(4), nsdSliprateVector(4), ttao, backSliprate
    
    dtau = 0.0d0 
    radius = sqrt((x(1,nsmp(1,iFaultNodePair,iFault))-xsource)**2 + &
                  (x(2,nsmp(1,iFaultNodePair,iFault))-ysource)**2 + &
                  (x(3,nsmp(1,iFaultNodePair,iFault))-zsource)**2)
    
    T    = 1.0d0
    F    = 0.0d0
    G    = 1.0d0
    if (radius < nucR) then 
        F = dexp(radius**2/(radius**2-nucR**2))
    endif 
    if (time<=T)  then 
        G = dexp((time-T)**2/(time*(time-2.0d0*T)))
    endif
    
    if (TPV == 105 .or. TPV == 104) then
        dtau = nucdtau0*F*G
    elseif (TPV == 2802) then
        !tpv2802 is drv.a6, plastic.
        if (nt == 1) then  
            fric(81,iFaultNodePair,iFault) = nucdtau0 !nsdTractionVector(2)*perturb
            ttao = sqrt(nsdTractionVector(2)**2 + nsdTractionVector(3)**2)
            backSliprate = sqrt((nsdSliprateVector(2)+fric(26,iFaultNodePair,iFault))**2 + &
                            (nsdSliprateVector(3)+fric(27,iFaultNodePair,iFault))**2)
            fric(20,iFaultNodePair,iFault) = fric(9,iFaultNodePair,iFault)*dlog(2.0d0*fric(12,iFaultNodePair,iFault)/backSliprate &
                                        *dsinh(ttao/abs(nsdTractionVector(1))/fric(9,iFaultNodePair,iFault)))
            fric(23,iFaultNodePair,iFault) = abs(nsdTractionVector(1))
        endif 
        dtau = fric(81,iFaultNodePair,iFault)*F*G
    endif
    
    ! if (radius < 1.d3) then 
        ! write(*,*) '========='
        ! write(*,*) x(1,nsmp(1,iFaultNodePair,iFault)), x(2,nsmp(1,iFaultNodePair,iFault)), x(3,nsmp(1,iFaultNodePair,iFault))
        ! write(*,*) radius
        ! write(*,*) xsource, ysource, zsource
        ! write(*,*) F, G, dtau
    ! endif 
    
    nsdTractionVector(2) = nsdTractionVector(2) + dtau
    
end subroutine rsfNucleation

subroutine swtwNucleation(iFault, iFaultNodePair, fricCoeff)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault, iFaultNodePair
    real (kind = dp) :: radius, tr, tc, fricCoeff
    
    radius = sqrt((x(1,nsmp(1,iFaultNodePair,iFault))-xsource)**2 + &
                  (x(2,nsmp(1,iFaultNodePair,iFault))-ysource)**2 + &
                  (x(3,nsmp(1,iFaultNodePair,iFault))-zsource)**2)
    
    tr = 1.0d9 
    if(radius <= nucR) then 
        if (TPV == 201) tr = (radius+0.081d0*nucR*(1.0d0/(1.0d0-(radius/nucR)**2)-1.0d0))/(0.7d0*3464.d0)
        if (TPV == 202) tr = radius/nucRuptVel
    endif
    
    tc = 1.0d0 
    if(time<tr) then 
        tc = 0.0d0
    elseif ((time<(tr+fric(5,iFaultNodePair,iFault))).and.(time>=tr)) then 
        tc = (time-tr)/fric(5,iFaultNodePair,iFault)
    endif

    fricCoeff = min(fric(1,iFaultNodePair,iFault)+(fric(2,iFaultNodePair,iFault)-fric(1,iFaultNodePair,iFault))*tc, fricCoeff)
    ! tmp1 = fs+(fd-fs)*tc 
    ! tmp2 = xmu
    ! xmu  = min(tmp1,tmp2)  
end subroutine

subroutine showNewton(iFault,iFaultNodePair,iv,v_trial,nsdTractionVector)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault,iFaultNodePair,iv
    real (kind = dp) :: v_trial, nsdTractionVector(4)
    if (x(1,nsmp(1,iFaultNodePair,iFault))==-22500.d0 .and. &
        x(3,nsmp(1,iFaultNodePair,iFault))==-21500.d0) then
        write(*,*) 'Newton: iteration step iv = ', iv
        write(*,*) 'Newton: trial slip rate is ', v_trial
        write(*,*) 'Newton: nsdTractionVector', nsdTractionVector(1)/1.0d6, nsdTractionVector(2)/1.0d6, nsdTractionVector(3)/1.0d6, nsdTractionVector(4)/1.0d6 
        write(*,*) 'Newton: state variable is ', fric(20,iFaultNodePair,iFault)
    endif
end subroutine showNewton

subroutine NewtonRaphson(iFault, iFaultNodePair, v_trial, taoc_new, state0, thetaPc0, thetaPcDot, nsdTractionVector, trialTractVec, T_coeff)
    use globalvar
    implicit none 
    integer (kind = 4) :: iFault, iFaultNodePair, ivmax=20, iv
    real (kind = dp) :: v_trial, newSliprate, state0, stateTmp, thetaPc0, thetaPcTmp, thetaPcDot, xmu, dxmudv, nsdTractionVector(4), trialTractVec(4), T_coeff, rsfeq, drsfeqdv, taoc_new
    ! Netwon solver for slip rate, v_trial, for the next time step.
    
    do iv = 1,ivmax
        ! in each iteration, reupdate the new state variable [fric(20)] given the new 
        !   slip rate, v_trial.
        !fric(20,iFaultNodePair,iFault)  = statetmp
        stateTmp = state0
        if(friclaw == 3) then
            call rate_state_ageing_law(v_trial,stateTmp,fric(1,iFaultNodePair,iFault),xmu,dxmudv)
        else
            call rate_state_slip_law(v_trial,stateTmp,fric(1,iFaultNodePair,iFault),xmu,dxmudv)
        endif 
        
        ! [NOTE]: the code doesn't support normal stress evolution under thermo pressurization.
        if (friclaw < 5) then
            !fric(23,iFaultNodePair,iFault)  = thetaPcTmp 
            thetaPcTmp = thetaPc0
            call rate_state_normal_stress(v_trial, thetaPcTmp, thetaPcDot, nsdTractionVector(1), fric(1,iFaultNodePair,iFault))    
            taoc_new        = xmu*thetaPcTmp
            rsfeq           = v_trial + T_coeff * (taoc_new*0.5d0 - trialTractVec(4))
            drsfeqdv        = 1.0d0 + T_coeff * (dxmudv * thetaPcTmp)*0.5d0  
        else
            taoc_new        = fric(4,iFaultNodePair,iFault) - xmu * MIN(nsdTractionVector(1), 0.0d0)
            rsfeq           = v_trial + T_coeff * (taoc_new*0.5d0 - trialTractVec(4))
            drsfeqdv        = 1.0d0 + T_coeff * (-dxmudv * MIN(nsdTractionVector(1),0.0d0))*0.5d0  
        endif
        
        ! exiting criteria:
        !   1. relative residual, rsfeq/drsfeqdv, is smaller than 1e-14*v_trial
        !   2. residual, rsfeq, is smaller than 1e-6*v_trial
        if(abs(rsfeq/drsfeqdv) < 1.d-14 * abs(v_trial) .and. abs(rsfeq) < 1.d-6 * abs(v_trial)) exit 
        !if(abs(rsfeq) < 1.d-5 * abs(v_trial)) exit 
            newSliprate = v_trial - rsfeq / drsfeqdv
        
        ! additional constraints for solving trial slip rate, v_trial
        !   if vtmp smaller than zero, reset it to half of v_trial in the last try. 
        if(newSliprate <= 0.0d0) then
            v_trial = v_trial/2.0d0
        else
            v_trial = newSliprate
        endif  
        !call showNewton(iFault,iFaultNodePair,iv,v_trial,nsdTractionVector)
    enddo !iv

    ! If cannot find a solution for v_trial, manually set it to a small value, typically the creeping rate.
    ! Also reset taoc_new to 2 X ttao1.
    ! Without this, TPV1053D blew up at the surface station (-4.2,0)
    if(TPV==105 .and. v_trial < fric(46,iFaultNodePair,iFault)) v_trial = fric(46,iFaultNodePair,iFault) ! necessary for tpv1053d
    
    state0   = stateTmp
    thetaPc0 = thetaPcTmp
    
end subroutine 

subroutine storeOnFaultStationQuantSCEC(iFault, iFaultNodePair, nsdSlipVector, nsdSliprateVector, nsdTractionVector)
    use globalvar
    implicit none
    integer (kind = 4) :: iFault, iFaultNodePair, j
    real (kind = dp) :: nsdSlipVector(4), nsdSliprateVector(4), nsdTractionVector(4)
    if(n4onf>0 .and. lstr) then    
        do j = 1, n4onf
            if(anonfs(1,j)==iFaultNodePair .and. anonfs(3,j)==iFault) then 
                fltsta(1,locplt-1,j)  = time
                fltsta(2,locplt-1,j)  = nsdSliprateVector(2)
                fltsta(3,locplt-1,j)  = nsdSliprateVector(3)
                fltsta(4,locplt-1,j)  = fric(20,iFaultNodePair,iFault)
                fltsta(5,locplt-1,j)  = nsdSlipVector(2)
                fltsta(6,locplt-1,j)  = nsdSlipVector(3)
                fltsta(7,locplt-1,j)  = nsdSlipVector(1)
                fltsta(8,locplt-1,j)  = nsdTractionVector(2) !tstk
                fltsta(9,locplt-1,j)  = nsdTractionVector(3) !tdip
                fltsta(10,locplt-1,j) = nsdTractionVector(1) !tnrm
                fltsta(11,locplt-1,j) = fric(51,iFaultNodePair,iFault) + fric(42,iFaultNodePair,iFault) ! + fric_tp_pini
                fltsta(12,locplt-1,j) = fric(52,iFaultNodePair,iFault) 
            endif
        enddo 
    endif   
    
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
end subroutine storeOnFaultStationQuantSCEC
