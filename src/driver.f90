! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine driver

    use globalvar
    implicit none
    include 'mpif.h'

    do nt = 1, nstep

        time = time + dt
        
        if (mod(nt,100) == 1 .and. me == master) then
            write(*,*) '=                                                                   ='
            write(*,*) '=     Current time in dynamic rupture                               ='
            write(*,'(X,A,40X,f7.3,4X,A)') '=',  time  , 's'
        endif
        
        call velDispUpdate
        call offFaultStationSCEC
        
        nodalForceArr=0.0d0
        
        call assembleGlobalKU
        call hrglss   
        call MPI4NodalQuant(nodalForceArr, 3)
        if (friclaw == 5) call thermop
        call faulting
        nodalForceArr(1:totalNumOfEquations) = nodalForceArr(1:totalNumOfEquations)/nodalMassArr(1:totalNumOfEquations)
        if ((mod(nt,10) == 1) .and. (outputGroundMotion == 1)) call output_gm
    enddo 

end subroutine driver

subroutine init_vel
    ! initiate the 1d velocity array v1. 
    ! if mode==2, non-zero values for fric(31-36,i,ift) loaded from 
    !    the restart file.
    ! if mode==1, fric(31-36,i) will be zeros.
    use globalvar
    implicit none
    integer (kind = 4) :: i, ift, tmp
    
    do ift = 1, ntotft
        do i = 1,nftnd(ift)
            tmp = eqNumStartIndexLoc(nsmp(1,i,ift))! slave nodeid i 
            v1(eqNumIndexArr(tmp+1)) = fric(34,i,ift) ! vxs
            v1(eqNumIndexArr(tmp+2)) = fric(35,i,ift) ! vys
            v1(eqNumIndexArr(tmp+3)) = fric(36,i,ift) ! vzs
            tmp = eqNumStartIndexLoc(nsmp(2,i,ift))! master nodeid i
            v1(eqNumIndexArr(tmp+1)) = fric(31,i,ift) ! vxm
            v1(eqNumIndexArr(tmp+2)) = fric(32,i,ift) ! vym
            v1(eqNumIndexArr(tmp+3)) = fric(33,i,ift) ! vzm
        enddo
    enddo
end subroutine init_vel

subroutine doubleCouplePointSource
    use globalvar 
    implicit none
    integer (kind = 4) :: i
    
    if (C_dc==1)then
        do i=1,totalNumOfNodes
        !x positive. Adding a point force in y+ direction.
        if (meshCoor(1,i)==100.0.and.meshCoor(2,i)==0.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'right',i,me
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+1.0e14/2
            endif
        endif
        !x negative. Adding a point force in y- direction.
        if (meshCoor(1,i)==-100.0.and.meshCoor(2,i)==0.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'left',i,me    
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-1.0e14/2
            endif
        endif
        !y positive. Adding a point force in x+ direction.
        if (meshCoor(1,i)==0.0.and.meshCoor(2,i)==100.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'up',i,me    
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+1.0e14/2
            endif
        endif
        !y negative. Adding a point force in x- direction.
        if (meshCoor(1,i)==0.0.and.meshCoor(2,i)==-100.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'down',i,me    
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))-&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))-1.0e14/2
            endif
        endif
        enddo!Enddo double-couple point source.    
    endif!ldc(logical double couple)    
end subroutine doubleCouplePointSource

subroutine velDispUpdate
    use globalvar
    implicit none 
    include 'mpif.h'
    
    integer (kind = 4) :: i, j, eqNumTmp
    real (kind = dp) :: dampv(9)
   
    time1 = MPI_WTIME()
    do i = 1, totalNumOfNodes
        if (numOfDofPerNodeArr(i)==3) then
            do j = 1, 3
                !eqNumStartIndexLoc(i) = eqNumSt
                eqNumTmp = eqNumIndexArr(eqNumStartIndexLoc(i)+j)
                v1(eqNumTmp) = v1(eqNumTmp) + nodalForceArr(eqNumTmp)*dt
                velArr(j,i) = v1(eqNumTmp)
                dispArr(j,i) = dispArr(j,i) + v1(eqNumTmp)*dt
            enddo
        elseif (numOfDofPerNodeArr(i)==12) then
            eqNumTmp=eqNumIndexArr(eqNumStartIndexLoc(i)+1)
            if (eqNumTmp>0) then
                call comdampv(meshCoor(1,i), meshCoor(2,i), meshCoor(3,i), dampv)
            endif
            do j = 1, 9
                eqNumTmp = eqNumIndexArr(eqNumStartIndexLoc(i)+j)
                if (eqNumTmp>0) then
                    v1(eqNumTmp) = (nodalForceArr(eqNumTmp)+v1(eqNumTmp)*(1.0d0/dt-dampv(j)/2.0d0))/(1.0d0/dt+dampv(j)/2.0d0)
                endif
            enddo
            do j = 10, 12
                eqNumTmp = eqNumIndexArr(eqNumStartIndexLoc(i)+j)
                if (eqNumTmp>0) then
                    v1(eqNumTmp) = v1(eqNumTmp) + nodalForceArr(eqNumTmp)*dt
                endif                
            enddo
            ! Update final velocity.
            eqNumTmp = eqNumIndexArr(eqNumStartIndexLoc(i)+1)
            if (eqNumTmp>0) then
                velArr(1,i) = v1(eqNumIndexArr(eqNumStartIndexLoc(i)+1)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+2)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+3)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+10))
                velArr(2,i) = v1(eqNumIndexArr(eqNumStartIndexLoc(i)+4)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+5)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+6)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+11))
                velArr(3,i) = v1(eqNumIndexArr(eqNumStartIndexLoc(i)+7)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+8)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+9)) + &
                        v1(eqNumIndexArr(eqNumStartIndexLoc(i)+12))
                dispArr(1,i) = dispArr(1,i) + velArr(1,i)*dt
                dispArr(2,i) = dispArr(2,i) + velArr(2,i)*dt
                dispArr(3,i) = dispArr(3,i) + velArr(3,i)*dt
            elseif (eqNumTmp == -1) then
                velArr(1:3,i) = 0.0d0
                dispArr(1:3,i) = 0.0d0
            endif    
        endif
        if ((velArr(1,i)/=velArr(1,i)).or.velArr(2,i)/=velArr(2,i).or.velArr(3,i)/=velArr(3,i)) then 
            write(*,*) 'Velocity NaN at point ', meshCoor(1,i), meshCoor(2,i), meshCoor(3,i), ' at time step ', nt
            stop
        endif
    enddo
    timeused(3) = timeused(3) + MPI_WTIME() - time1
end subroutine velDispUpdate

subroutine offFaultStationSCEC
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j, l, k, k1
    
        if (ndout>0) then
            dout(1,nt) = time
            do i = 1, ndout 
                j = idhist(1,i)
                if (j<=0) j=1  !avoid zero that cannot be used below
                    k = idhist(2,i)
                    l = idhist(3,i)
                if (l ==1) then
                    dout(i+1,nt) = dispArr(k,j)
                elseif(l==2) then
                    dout(i+1,nt) = velArr(k,j)
                elseif(l==3) then
                    k1 = eqNumIndexArr(eqNumStartIndexLoc(j)+k)
                    dout(i+1,nt) = nodalForceArr(k1)
                endif
            enddo
        endif
end subroutine offFaultStationSCEC
 
