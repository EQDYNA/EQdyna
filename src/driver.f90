! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine driver

    use globalvar
    implicit none
    include 'mpif.h'

    integer (kind = 4) :: i

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
        
        call ku
        call hrglss   
        call MPI4NodalQuant(nodalForceArr, 3)
        if (friclaw == 5) call thermop
        call faulting
        nodalForceArr(1:totalNumOfEquations) = nodalForceArr(1:totalNumOfEquations)/nodalMassArr(1:totalNumOfEquations) ! timeused(7), depreciated
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
        if (x(1,i)==100.0.and.x(2,i)==0.0.and.x(3,i)==-2000.)then
            !write(*,*) 'right',i,me
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))+1.0e14/2
            endif
        endif
        !x negative. Adding a point force in y- direction.
        if (x(1,i)==-100.0.and.x(2,i)==0.0.and.x(3,i)==-2000.)then
            !write(*,*) 'left',i,me    
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+2))-1.0e14/2
            endif
        endif
        !y positive. Adding a point force in x+ direction.
        if (x(1,i)==0.0.and.x(2,i)==100.0.and.x(3,i)==-2000.)then
            !write(*,*) 'up',i,me    
            if (time<0.2)then
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+&
                (1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
            else
                nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))=nodalForceArr(eqNumIndexArr(eqNumStartIndexLoc(i)+1))+1.0e14/2
            endif
        endif
        !y negative. Adding a point force in x- direction.
        if (x(1,i)==0.0.and.x(2,i)==-100.0.and.x(3,i)==-2000.)then
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
    
    integer (kind = 4) :: i, j, eqNumTmp, itag
    real (kind = dp) :: dampv(9)
   
    time1   = MPI_WTIME()
    do i=1,totalNumOfNodes
        if (numOfDofPerNodeArr(i)==3) then
            do j=1,3
                itag=eqNumStartIndexLoc(i)+j
                eqNumTmp=eqNumIndexArr(itag)
                v1(eqNumTmp)=v1(eqNumTmp)+nodalForceArr(eqNumTmp)*dt
                d1(eqNumTmp)=d1(eqNumTmp)+v1(eqNumTmp)*dt
                v(j,i)=v1(eqNumTmp)
                d(j,i)=d1(eqNumTmp)
            enddo
        elseif (numOfDofPerNodeArr(i)==12) then
            itag=eqNumStartIndexLoc(i)+1
            eqNumTmp=eqNumIndexArr(itag)
            if (eqNumTmp>0) then
                call comdampv(x(1,i),x(2,i),x(3,i),dampv)
            endif
            do j=1,9
                itag=eqNumStartIndexLoc(i)+j
                eqNumTmp=eqNumIndexArr(itag)
                if (eqNumTmp>0) then
                    v1(eqNumTmp)=(nodalForceArr(eqNumTmp)+v1(eqNumTmp)*(1/dt-dampv(j)/2))/(1/dt+dampv(j)/2)
                endif
            enddo
            do j=10,12
                itag=eqNumStartIndexLoc(i)+j
                eqNumTmp=eqNumIndexArr(itag)
                if (eqNumTmp>0) then
                    v1(eqNumTmp)=v1(eqNumTmp)+nodalForceArr(eqNumTmp)*dt
                endif                
            enddo
            ! Update final velocity.
            itag=eqNumStartIndexLoc(i)
            eqNumTmp=eqNumIndexArr(itag+1)
            if (eqNumTmp>0) then
                v(1,i)=v1(eqNumIndexArr(itag+1))+v1(eqNumIndexArr(itag+2))+v1(eqNumIndexArr(itag+3))+v1(eqNumIndexArr(itag+10))
                v(2,i)=v1(eqNumIndexArr(itag+4))+v1(eqNumIndexArr(itag+5))+v1(eqNumIndexArr(itag+6))+v1(eqNumIndexArr(itag+11))
                v(3,i)=v1(eqNumIndexArr(itag+7))+v1(eqNumIndexArr(itag+8))+v1(eqNumIndexArr(itag+9))+v1(eqNumIndexArr(itag+12))
                d(1,i)=d(1,i)+v(1,i)*dt
                d(2,i)=d(2,i)+v(2,i)*dt
                d(3,i)=d(3,i)+v(3,i)*dt
            elseif (eqNumTmp==-1) then
                v(1,i)=0.0d0
                v(2,i)=0.0d0
                v(3,i)=0.0d0
                d(1,i)=0.0d0
                d(2,i)=0.0d0
                d(3,i)=0.0d0                
            endif    
        endif
        if ((v(1,i)/=v(1,i)).or.v(2,i)/=v(2,i).or.v(3,i)/=v(3,i)) then 
            write(*,*) x(1,i),x(2,i),x(3,i), 'nt=',nt
            stop 'NAN'
        endif
    enddo
    timeused(3) = timeused(3) + MPI_WTIME() - time1
end subroutine velDispUpdate

subroutine offFaultStationSCEC
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j, l, k, k1
    
    if (mod(nt,nhplt) == 0) then    
        lstr    = .true.    
        locplt  = locplt+ 1    !when nt=1, locplt=2 due to 1 in eqdy3d.f90
    else
        lstr    = .false.
    endif
    if (lstr) then
        if((ndout>0).and.(locplt>1)) then
            dout(1,locplt)=time
            do i=1,ndout 
                j=idhist(1,i)
                if(j<=0) j=1  !avoid zero that cannot be used below
                    k=idhist(2,i)
                    l=idhist(3,i)
                if(l==1) then
                    dout(i+1,locplt)=d(k,j)
                elseif(l==2) then
                    dout(i+1,locplt)=v(k,j)
                elseif(l==3) then
                    k1=eqNumIndexArr(eqNumStartIndexLoc(j)+k)
                    dout(i+1,locplt)=nodalForceArr(k1)
                endif
            enddo
        endif
    endif
end subroutine offFaultStationSCEC
 
