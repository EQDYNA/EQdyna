! Temporary instrumentation for the Python-port feasibility spike (not part
! of the production solver). Dumps the static/geometric state EQdyna has
! built by the time the explicit time loop starts, so the Python port can
! load real Fortran-computed mesh/mass/shape-function data instead of
! re-deriving mesh4num/meshgen/assembleGlobalMass in Python for this spike.
subroutine pydump_state
    use globalvar
    implicit none
    integer (kind = 4) :: i, j, k, u

    u = 90001
    open(unit=u, file='pydump_header.txt', status='unknown')
    write(u,*) totalNumOfNodes
    write(u,*) totalNumOfElements
    write(u,*) totalNumOfEquations
    write(u,*) nen
    write(u,*) ned
    write(u,*) nftnd(1)
    write(u,*) ntotft
    write(u,*) nstep
    write(u,*) dt
    write(u,*) w
    write(u,*) rdampk
    write(u,*) rdampm
    write(u,*) kapa_hg
    write(u,*) R
    write(u,*) nPML
    write(u,*) vmaxPML
    write(u,*) (PMLb(i), i=1,8)
    write(u,*) grav
    write(u,*) C_elastic
    write(u,*) roumax
    write(u,*) rhow
    write(u,*) gamar
    write(u,*) slipRateThres
    write(u,*) xsource, ysource, zsource
    write(u,*) nucR, nucRuptVel, nucdtau0, nucT
    write(u,*) TPV
    write(u,*) C_nuclea
    write(u,*) nucfault
    write(u,*) friclaw
    write(u,*) timeElapsed
    close(u)

    open(unit=u, file='pydump_meshCoor.txt', status='unknown')
    do i = 1, totalNumOfNodes
        write(u,'(3(1x,e24.16e3))') meshCoor(1,i), meshCoor(2,i), meshCoor(3,i)
    enddo
    close(u)

    open(unit=u, file='pydump_conn.txt', status='unknown')
    do i = 1, totalNumOfElements
        write(u,'(8(1x,i8),1x,i3,5(1x,e24.16e3))') (nodeElemIdRelation(j,i), j=1,nen), elemTypeArr(i), &
            (mat(i,j), j=1,5)
    enddo
    close(u)

    open(unit=u, file='pydump_elemgeo.txt', status='unknown')
    do i = 1, totalNumOfElements
        write(u,'(1x,e24.16e3, 24(1x,e24.16e3), 6(1x,e24.16e3), 32(1x,e24.16e3))') eledet(i), &
            ((eleshp(j,k,i), j=1,3), k=1,8), (ss(j,i), j=1,6), ((phi(j,k,i), j=1,8), k=1,4)
    enddo
    close(u)

    open(unit=u, file='pydump_nodeinfo.txt', status='unknown')
    do i = 1, totalNumOfNodes
        write(u,'(1x,i2,1x,i10,12(1x,i10))') numOfDofPerNodeArr(i), eqNumStartIndexLoc(i), &
            (eqNumIndexArr(eqNumStartIndexLoc(i)+j), j=1,numOfDofPerNodeArr(i))
    enddo
    close(u)

    open(unit=u, file='pydump_nodalmass.txt', status='unknown')
    do i = 1, totalNumOfEquations
        write(u,'(1x,e24.16e3)') nodalMassArr(i)
    enddo
    close(u)

    open(unit=u, file='pydump_fnms.txt', status='unknown')
    do i = 1, totalNumOfNodes
        write(u,'(1x,e24.16e3)') fnms(i)
    enddo
    close(u)

    open(unit=u, file='pydump_fault.txt', status='unknown')
    do i = 1, nftnd(1)
        write(u,'(2(1x,i10),3(1x,e24.16e3),3(1x,e24.16e3),3(1x,e24.16e3),1x,e24.16e3,100(1x,e24.16e3))') &
            nsmp(1,i,1), nsmp(2,i,1), un(1,i,1),un(2,i,1),un(3,i,1), us(1,i,1),us(2,i,1),us(3,i,1), &
            ud(1,i,1),ud(2,i,1),ud(3,i,1), arn(i,1), (fric(j,i,1), j=1,100)
    enddo
    close(u)

    open(unit=u, file='pydump_v1.txt', status='unknown')
    do i = 1, totalNumOfEquations
        write(u,'(1x,e24.16e3)') v1(i)
    enddo
    close(u)

    ! Milestone 5 (standalone station matching): anonfs(1:3,1:numOfOnFaultStCount)
    ! maps on-fault station index i -> (fault-node-seq, xonfs column j, iFault);
    ! OffFaultStNodeIdIndex(1:2,1:numOfOffFaultStCount) maps off-fault station
    ! index -> (x4nds column i, matched nodeCount).
    open(unit=u, file='pydump_stations.txt', status='unknown')
    write(u,*) numOfOnFaultStCount
    write(u,*) numOfOffFaultStCount
    do i = 1, numOfOnFaultStCount
        write(u,'(3(1x,i10))') anonfs(1,i), anonfs(2,i), anonfs(3,i)
    enddo
    do i = 1, numOfOffFaultStCount
        write(u,'(2(1x,i10))') OffFaultStNodeIdIndex(1,i), OffFaultStNodeIdIndex(2,i)
    enddo
    close(u)
end subroutine pydump_state

! Per-step checkpoint (nt<=5) for first-divergence diagnosis: end-of-step
! nodal acceleration (nodalForceArr, already divided by mass), full velArr/
! dispArr, and the fault-node friction-state columns the port needs to
! reproduce (71-73 slip, 74-77 sliprate/cumslip, 47 final sliprate, 78-80
! traction, 20 RSF state (unused, friclaw=1), 23 normal-stress state).
subroutine pydump_step(step)
    use globalvar
    implicit none
    integer (kind = 4) :: step, i, j
    character(len=8) :: stmp
    write(stmp,'(i0)') step

    open(unit=90101, file='pydump_step'//trim(stmp)//'_accel.txt', status='unknown')
    do i = 1, totalNumOfEquations
        write(90101,'(1x,e24.16e3)') nodalForceArr(i)
    enddo
    close(90101)

    open(unit=90102, file='pydump_step'//trim(stmp)//'_veldisp.txt', status='unknown')
    do i = 1, totalNumOfNodes
        write(90102,'(6(1x,e24.16e3))') velArr(1,i),velArr(2,i),velArr(3,i), &
                                         dispArr(1,i),dispArr(2,i),dispArr(3,i)
    enddo
    close(90102)

    open(unit=90103, file='pydump_step'//trim(stmp)//'_fault.txt', status='unknown')
    do i = 1, nftnd(1)
        write(90103,'(11(1x,e24.16e3))') fric(71,i,1),fric(72,i,1),fric(73,i,1), &
            fric(74,i,1),fric(75,i,1),fric(76,i,1),fric(77,i,1),fric(47,i,1), &
            fric(78,i,1),fric(79,i,1),fric(80,i,1)
    enddo
    close(90103)
end subroutine pydump_step
