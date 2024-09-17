! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine driver

    use globalvar
    implicit none
    include 'mpif.h'

    do nt = 1, nstep

        timeElapsed = timeElapsed + dt
        
        if (mod(nt,100) == 1 .and. me == masterProcsId) then
            write(*,*) '=                                                                   ='
            write(*,*) '=     Current time in dynamic rupture                               ='
            write(*,'(X,A,40X,f7.3,4X,A)') '=',  timeElapsed  , 's'
        endif
        
        call velDispUpdate
        call storeOffFaultStData
        
        nodalForceArr = 0.0d0
        
        call assembleGlobalKU
        call hrglss   
        call MPI4NodalQuant(nodalForceArr, 3)
        if (friclaw == 5) call thermop
        call faulting
        nodalForceArr(1:totalNumOfEquations) = nodalForceArr(1:totalNumOfEquations)/nodalMassArr(1:totalNumOfEquations)
        if ((mod(nt,10) == 1) .and. (outputGroundMotion == 1)) call output_gm
    enddo 

end subroutine driver

subroutine doubleCouplePointSource
    use globalvar 
    implicit none
    integer (kind = 4) :: i
    
    if (C_dc==1)then
        do i=1,totalNumOfNodes
        !x positive. Adding a point force in y+ direction.
        if (meshCoor(1,i)==100.0.and.meshCoor(2,i)==0.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'right',i,me
            if (timeElapsed<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+&
                (1.0e14/0.2/2*timeElapsed-1.0e14/2/2/3.1415*sin(2*3.1415*timeElapsed/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+1.0e14/2
            endif
        endif
        !x negative. Adding a point force in y- direction.
        if (meshCoor(1,i)==-100.0.and.meshCoor(2,i)==0.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'left',i,me    
            if (timeElapsed<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-&
                (1.0e14/0.2/2*timeElapsed-1.0e14/2/2/3.1415*sin(2*3.1415*timeElapsed/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-1.0e14/2
            endif
        endif
        !y positive. Adding a point force in x+ direction.
        if (meshCoor(1,i)==0.0.and.meshCoor(2,i)==100.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'up',i,me    
            if (timeElapsed<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+&
                (1.0e14/0.2/2*timeElapsed-1.0e14/2/2/3.1415*sin(2*3.1415*timeElapsed/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+1.0e14/2
            endif
        endif
        !y negative. Adding a point force in x- direction.
        if (meshCoor(1,i)==0.0.and.meshCoor(2,i)==-100.0.and.meshCoor(3,i)==-2000.)then
            !write(*,*) 'down',i,me    
            if (timeElapsed<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))-&
                (1.0e14/0.2/2*timeElapsed-1.0e14/2/2/3.1415*sin(2*3.1415*timeElapsed/0.2))
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
   
    startTimeStamp = MPI_WTIME()
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
            write(*,*) 'brhs at this node are', nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i))+1)
            write(*,*) 'alhs at this node are', nodalMassArr(eqNumIndexArr(eqNumStartIndexLoc(i))+1:eqNumIndexArr(eqNumStartIndexLoc(i))+3)
            stop
        endif
    enddo
    compTimeInSeconds(3) = compTimeInSeconds(3) + MPI_WTIME() - startTimeStamp
end subroutine velDispUpdate

subroutine storeOffFaultStData
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, nodeId, quantType, k, eqNum
    
        if (numOfOffFaultStCount*ndof*2>0) then
            OffFaultStGramSCEC(1,nt) = timeElapsed
            do i = 1, numOfOffFaultStCount*ndof*2 
                nodeId = idhist(1,i)
                !if (j<=0) j=1  !avoid zero that cannot be used below
                    k = idhist(2,i)
                    quantType = idhist(3,i)
                if (quantType ==1) then
                    OffFaultStGramSCEC(i+1,nt) = dispArr(k,nodeId)
                elseif(quantType == 2) then
                    OffFaultStGramSCEC(i+1,nt) = velArr(k,nodeId)
                elseif(quantType == 3) then
                    eqNum = eqNumIndexArr(eqNumStartIndexLoc(nodeId)+k)
                    OffFaultStGramSCEC(i+1,nt) = nodalForceArr(eqNum)
                endif
            enddo
        endif
end subroutine storeOffFaultStData
 
