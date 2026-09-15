! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine meshgen
    ! Create regular node grids and hexahedral elements for this MPI process.
    ! Create split nodes.
    ! Distort the regular mesh for complex fault geometries.
    ! Set material propertiesa, initial stress, and other element-wise properties. 
    
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    ! incremental variables 
    integer (kind = 4) :: nodeCount=0, elemCount=0, equationNumCount=0, eqNumIndexArrLocTag=0, stressDofCount=0
    integer (kind = 4) :: nxt,nyt,nzt,nx,ny,nz,ix,iy,iz, &
                       edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,edgezn, &
                       nxuni,nyuni,nzuni,ift, &
                       n1,n2,n3,n4,m1,m2,m3,m4,&
                       mex,mey,mez,itmp1,&
                       numOfDofPerNodeTmp, msnode, nodeXyzIndex(10),isOnFt
    integer (kind = 4), dimension(ntotft) :: nftnd0,ixfi,izfi,ifs,ifd
    integer (kind = 4), allocatable :: fltrc(:,:,:,:)
    ! Temporary real variables
    real (kind = dp) :: nodeCoor(10), elementCenterCoor(3), &
                       a,b,area,aa1,bb1,cc1,dd1,p1,q1, ycoort, pfx = 0.0d0, pfz = 0.0d0
    real (kind = dp) :: xline(10000), yline(10000), zline(10000), modelBoundCoor(3,2)

    call calcXyzMPIId(mex, mey, mez)
    call getLocalOneDimCoorArrAndSize(nxt, nxuni, edgex1, mex, nx, xline, modelBoundCoor, 1)
    call getLocalOneDimCoorArrAndSize(nyt, nyuni, edgey1, mey, ny, yline, modelBoundCoor, 2)
    call getLocalOneDimCoorArrAndSize(nzt, nzuni, edgezn, mez, nz, zline, modelBoundCoor, 3)
    xmin = modelBoundCoor(1,1)
    xmax = modelBoundCoor(1,2)
    ymin = modelBoundCoor(2,1)
    ymax = modelBoundCoor(2,2)
    zmin = modelBoundCoor(3,1)
    zmax = modelBoundCoor(3,2)
    ! Move adjacent plane1 and plane2 to create hexahedral meshes
    allocate(plane1(ny+ntotft,nz),plane2(ny+ntotft,nz),fltrc(2,nxuni,nzuni,ntotft))
    plane1 = 0
    plane2 = 0
    
    numcount    = 0
    numcount(1) = nx
    numcount(2) = ny
    numcount(3) = nz
    
    ! Initialize scalars

    msnode   = nx*ny*nz
    numOfOnFaultStCount    = 0
    numOfOffFaultStCount    = 0
    OffFaultStNodeIdIndex   = 0
    ! Initialize arrays
    nftnd0   = 0
    
    ixfi = 0
    izfi = 0
    
    allocate(n4yn(totalNumOfOffSt))
    n4yn = 0
    
    ! Loop over x,y,z grids to create nodes and elements
    ! Create node
    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                call initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
                call createNode(nodeCoor, xline(ix), yline(iy), zline(iz), nodeCount, nodeXyzIndex)
                if (insertFaultType > 0) then 
                    call insertFaultInterface(nodeCoor, ycoort, pfx, pfz)
                    meshCoor(2,nodeCount) = ycoort
                endif 
                
                call setNumDof(nodeCoor, numOfDofPerNodeTmp)

                eqNumStartIndexLoc(nodeCount) = eqNumIndexArrLocTag
                numOfDofPerNodeArr(nodeCount) = numOfDofPerNodeTmp
                
                call setEquationNumber(nodeXyzIndex, nodeCoor, eqNumIndexArrLocTag, equationNumCount, numOfDofPerNodeTmp)
                call setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nodeCount)
                call createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nodeCount, msnode, nftnd0, equationNumCount, eqNumIndexArrLocTag, &
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
                
                ! Create element
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                    call createElement(elemCount, stressDofCount, iy, iz, elementCenterCoor)
                    call checkPMLAlignment(elementCenterCoor)
                    call setElementMaterial(elemCount, elementCenterCoor)
                    
                    if (C_degen > 3.0d0) then 
                        call wedge(elementCenterCoor(1), elementCenterCoor(2), elementCenterCoor(3), elemCount, stressDofCount, iy, iz, nftnd0(1))
                        isOnFt=0
                        call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(1,elemCount)), 1, isOnFt)
                        if (isOnFt==1 .and. elemTypeArr(elemCount)==1) elemTypeArr(elemCount) = 13

                        isOnFt=0
                        call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(2,elemCount)), 1, isOnFt)
                        if (isOnFt==1 .and. elemTypeArr(elemCount)==1) elemTypeArr(elemCount) = 13

                    endif         
                    
                    call replaceSlaveWithMasterNode(nodeCoor, elemCount, nftnd0) 
                    if (C_elastic == 0) call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0, elemCount)          
                 endif!if element
            enddo!iy
        enddo!iz
        plane1 = plane2
    enddo!ix
    
    sizeOfStressDofIndexArr=stressDofCount
    
    call meshGenError(nx, ny, nz, nodeCount, msnode, elemCount, equationNumCount, eqNumIndexArrLocTag, nftnd0)
    
    ! compute on-fault area associated with each fault node pair and distance from source
    do ift=1,ntotft
        if(nftnd0(ift)>0) then
        !...element areas and distribute evenly to its four nodes
        do i=2,ifd(ift)
            do j=2,ifs(ift)
            !...4 nodes of quadrilateral
                n1 = fltrc(1,j,i,ift) !use nodal number in nodeCount
                n2 = fltrc(1,j-1,i,ift)
                n3 = fltrc(1,j-1,i-1,ift)
                n4 = fltrc(1,j,i-1,ift)
                m1 = fltrc(2,j,i,ift) !use nodal number in nftnd0
                m2 = fltrc(2,j-1,i,ift)
                m3 = fltrc(2,j-1,i-1,ift)
                m4 = fltrc(2,j,i-1,ift)
                !...calculate area of the quadrilateral
                !......if fault is not in coor axes plane
                aa1=sqrt((meshCoor(1,n2)-meshCoor(1,n1))**2 + (meshCoor(2,n2)-meshCoor(2,n1))**2 &
                + (meshCoor(3,n2)-meshCoor(3,n1))**2)
                bb1=sqrt((meshCoor(1,n3)-meshCoor(1,n2))**2 + (meshCoor(2,n3)-meshCoor(2,n2))**2 &
                + (meshCoor(3,n3)-meshCoor(3,n2))**2)
                cc1=sqrt((meshCoor(1,n4)-meshCoor(1,n3))**2 + (meshCoor(2,n4)-meshCoor(2,n3))**2 &
                + (meshCoor(3,n4)-meshCoor(3,n3))**2)
                dd1=sqrt((meshCoor(1,n1)-meshCoor(1,n4))**2 + (meshCoor(2,n1)-meshCoor(2,n4))**2 &
                + (meshCoor(3,n1)-meshCoor(3,n4))**2)
                p1=sqrt((meshCoor(1,n4)-meshCoor(1,n2))**2 + (meshCoor(2,n4)-meshCoor(2,n2))**2 &
                + (meshCoor(3,n4)-meshCoor(3,n2))**2)
                q1=sqrt((meshCoor(1,n3)-meshCoor(1,n1))**2 + (meshCoor(2,n3)-meshCoor(2,n1))**2 &
                + (meshCoor(3,n3)-meshCoor(3,n1))**2)
                area = 0.25d0 * sqrt(4*p1*p1*q1*q1 - &
               (bb1*bb1 + dd1*dd1 - aa1*aa1 -cc1*cc1)**2) 
                !...distribute above area to 4 nodes evenly
                area = 0.25d0 * area
                arn(m1,ift) = arn(m1,ift) + area
                arn(m2,ift) = arn(m2,ift) + area
                arn(m3,ift) = arn(m3,ift) + area
                arn(m4,ift) = arn(m4,ift) + area
                enddo
            enddo
        endif

        call MPI4arn(nx, ny, nz, mex, mey, mez, nftnd0(ift), ift)
    enddo  
end subroutine meshgen
!==================================================================================================
!**************************************************************************************************
!==================================================================================================

subroutine setElementMaterial(elemCount, elementCenterCoor)
! Subroutine velocityStructure will asign Vp, Vs and rho 
!   based on input from bMaterial.txt, which is created by 
!   case input file user_defined_param.py.

    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: elemCount, i
    real (kind = dp) :: elementCenterCoor(3), vptmp, vstmp, rhotmp
    
    if (nmat == 1 .and. n2mat == 3) then
    ! homogenous material
        mat(elemCount,1)  = material(1,1)
        mat(elemCount,2)  = material(1,2)
        mat(elemCount,3)  = material(1,3)
    elseif (nmat > 1 .and. n2mat == 4) then 
        ! 1D velocity structure
        ! material = material(i,j), i=1,nmat, and j=1,4
        ! for j
        !   1: bottom depth of a layer (should be positive in m).
        !   2: Vp, m/s
        !   3: Vs, m/s
        !   4: rho, kg/m3
        if (abs(elementCenterCoor(3)) < material(1,1)) then
            mat(elemCount,1)  = material(1,2)
            mat(elemCount,2)  = material(1,3)
            mat(elemCount,3)  = material(1,4)
        else
            do i = 2, nmat
                if (abs(elementCenterCoor(3)) < material(i,1) &
                    .and. abs(elementCenterCoor(3)) >= material(i-1,1)) then
                    
                    mat(elemCount,1)  = material(i,2)
                    mat(elemCount,2)  = material(i,3)
                    mat(elemCount,3)  = material(i,4)
                endif 
            enddo
        endif
    endif 
    
    ! calculate lambda and mu from Vp, Vs and rho.
    ! mu = Vs**2*rho
    mat(elemCount,5)  = mat(elemCount,2)**2*mat(elemCount,3)
    ! lambda = Vp**2*rho-2*mu
    mat(elemCount,4)  = mat(elemCount,1)**2*mat(elemCount,3)-2.0d0*mat(elemCount,5)

end subroutine setElementMaterial

subroutine MPI4arn(nx, ny, nz, mex, mey, mez, totalNumFaultNode, iFault)
! Add up arn from neighbor MPI blocks.
    use globalvar
    use errorCodes
    implicit none 
    include 'mpif.h'
    
    integer (kind = 4) :: bndl,bndr,bndf,bndb,bndd,bndu, nx, ny, nz, mex, mey, mez
    integer (kind = 4) :: totalNumFaultNode, iFault, i
    ! Initialize fltMPI(6) to .false.
    fltMPI=.false.
    
    if(fltnum(1) /= 0) allocate(fltl(fltnum(1)))
    if(fltnum(2) /= 0) allocate(fltr(fltnum(2)))
    if(fltnum(3) /= 0) allocate(fltf(fltnum(3)))
    if(fltnum(4) /= 0) allocate(fltb(fltnum(4)))
    if(fltnum(5) /= 0) allocate(fltd(fltnum(5)))
    if(fltnum(6) /= 0) allocate(fltu(fltnum(6)))
    
    fltnum = 0
    do i = 1, totalNumFaultNode
        if(mod(fltgm(i),10)==1 ) then 
            fltnum(1) = fltnum(1) + 1
            fltl(fltnum(1)) = i
        endif
        if(mod(fltgm(i),10)==2 ) then 
            fltnum(2) = fltnum(2) + 1
            fltr(fltnum(2)) = i
        endif
        if(mod(fltgm(i),100)-mod(fltgm(i),10)==10 ) then
            fltnum(3) = fltnum(3) + 1
            fltf(fltnum(3)) = i
        endif
        if(mod(fltgm(i),100)-mod(fltgm(i),10)==20 ) then
            fltnum(4) = fltnum(4) + 1
            fltb(fltnum(4)) = i
        endif
        if(fltgm(i)-mod(fltgm(i),100)==100 ) then
            fltnum(5) = fltnum(5) + 1
            fltd(fltnum(5)) = i
        endif
        if(fltgm(i)-mod(fltgm(i),100)==200 ) then
            fltnum(6) = fltnum(6) + 1
            fltu(fltnum(6)) = i
        endif
    enddo
    ! FIX (pathway_forward.md, "fault-plane-on-MPI-boundary halving"): arn is
    ! accumulated per fault node from the LOCAL ifs x ifd quad grid (lines
    ! ~115-152 above), independent of y -- so whenever this rank's boundary
    ! in a direction the fault has ZERO nominal extent in (per fltxyz; always
    ! y for EQdyna's x-z-plane faults) carries fault nodes, this rank has
    ! already built the FULL local fault-node arn there (both ranks rebuild
    ! the whole ifs x ifd grid identically), so summing across that boundary
    ! DUPLICATES rather than DIVIDES the area. Directions the fault has real
    ! extent in (x, z -- or any future non-degenerate direction) are genuine
    ! divisions: each rank only owns a slice of the fault there, and summing
    ! is correct (verified: scratch/mira_faultmpi_guard, x-split and z-split
    ! reproduce the serial arn exactly; y-split alone doubled it before this
    ! fix, matched after).
    !
    ! A prior attempt skipped the ENTIRE `call syncArnBoundary` for such
    ! directions and was reverted: fltMPI(k) is set at the top of
    ! syncArnBoundary and is ALSO the gate `addFaultBoundaryTerm`
    ! (assembleGlobalMass.f90) uses to decide whether fnms/nodalMassArr need
    ! a cross-rank correction on that same boundary -- skipping the whole
    ! call silently skipped that correction too and broke a rough-fault
    ! 2,2,1 run (zero mass / NaN velocity). This fix instead still calls
    ! syncArnBoundary (fltMPI(k) still gets set, the exchange still happens,
    ! fnms/nodalMassArr's correction is untouched) and only skips arn's OWN
    ! add-back when the direction is degenerate (see syncArnBoundary below).
    ! Direct evidence fnms is NOT affected by this bug (unlike arn): dumped
    ! fnms at the same physical fault node is bit-identical between serial
    ! and a fault-splitting decomposition (scratch/mira_faultmpi_guard) --
    ! fnms's volume-element accumulation is genuinely owned exclusively by
    ! one side's rank per split-node copy, so its cross-rank correction on
    ! this boundary is a no-op there in the cases checked; left unchanged.
    if (npx > 1) then
        bndl=1  !left boundary
        bndr=nx ! right boundary
        if (mex == masterProcsId) then
            bndl=0
        elseif (mex == npx-1) then
            bndr=0
        endif

        if (bndl/=0) then
            if(fltnum(1)>0 ) call syncArnBoundary(fltl, fltnum(1), 1, me-npy*npz, 1000, iFault, 1)
        endif !if bhdl/=0

        if (bndr/=0) then
            if(fltnum(2)>0 ) call syncArnBoundary(fltr, fltnum(2), 2, me+npy*npz, 1000, iFault, 1)
        endif !bndr/=0
    endif !npx>1
!*****************************************************************************************
    if (npy > 1) then   !MPI:no fault message passing along y due to the x-z vertical fault.
        bndf=1  !front boundary
        bndb=ny ! back boundary
        if (mey == masterProcsId) then
            bndf=0
        elseif (mey == npy-1) then
            bndb=0
        endif

        if (bndf/=0) then
            if(fltnum(3)>0) call syncArnBoundary(fltf, fltnum(3), 3, me-npz, 2000, iFault, 2)
        endif !bhdf/=0

        if (bndb/=0) then
            if(fltnum(4)>0) call syncArnBoundary(fltb, fltnum(4), 4, me+npz, 2000, iFault, 2)
        endif !bndb/=0
     endif !npy>1
!*****************************************************************************************
    if (npz > 1) then
        bndd=1  !lower(down) boundary
        bndu=nz ! upper boundary
        if (mez == masterProcsId) then
            bndd=0
        elseif (mez == npz-1) then
            bndu=0
        endif
        if (bndd/=0) then
            if(fltnum(5)>0) call syncArnBoundary(fltd, fltnum(5), 5, me-1, 3000, iFault, 3)
        endif !bhdd/=0

        if (bndu/=0) then
            if(fltnum(6)>0) call syncArnBoundary(fltu, fltnum(6), 6, me+1, 3000, iFault, 3)
        endif !bndu/=0
    endif !npz>1
contains
    subroutine syncArnBoundary(idxArr, n, k, neighbor, tagBase, ift, dimId)
    ! Send this rank's arn values for the given boundary's fault nodes to
    ! neighbor, receive neighbor's values for the same nodes, and accumulate
    ! -- UNLESS the fault has zero nominal extent in the physical dimension
    ! (dimId: 1=x, 2=y, 3=z) this boundary lies along, in which case this
    ! rank's local arn is already the fault's full local contribution and
    ! the neighbor's value is a duplicate, not a partial sum (see the NOTE
    ! above MPI4arn's npx/npy/npz blocks). fltMPI(k) is still set and the
    ! exchange still happens either way: addFaultBoundaryTerm
    ! (assembleGlobalMass.f90) depends on fltMPI(k)/the send-recv having run,
    ! independent of what this subroutine does with arn.
        integer (kind = 4), intent(in) :: idxArr(:), n, k, neighbor, tagBase, ift, dimId
        integer (kind = 4) :: j, jMPIstatus(MPI_STATUS_SIZE), jMPIerr
        real (kind = dp), allocatable, dimension(:) :: sendBuf, recvBuf

        fltMPI(k)=.true.
        allocate(sendBuf(n),recvBuf(n))
        sendBuf = 0.
        recvBuf = 0.
        do j = 1, n
            sendBuf(j)=arn(idxArr(j),ift)
        enddo
        call mpi_sendrecv(sendBuf, n, MPI_DOUBLE_PRECISION, neighbor, tagBase+me, &
            recvBuf, n, MPI_DOUBLE_PRECISION, neighbor, tagBase+neighbor, &
            MPI_COMM_WORLD, jMPIstatus, jMPIerr)
        if (fltxyz(2,dimId,ift) /= fltxyz(1,dimId,ift)) then
            do j = 1, n
                arn(idxArr(j),ift) = arn(idxArr(j),ift) + recvBuf(j)
            enddo
        endif
        deallocate(sendBuf,recvBuf)
    end subroutine syncArnBoundary
end subroutine MPI4arn

subroutine meshGenError(nx, ny, nz, nodeCount, msnode, elemCount, equationNumCount, eqNumIndexArrLocTag, nftnd0)
! Check consistency between mesh4 and meshgen
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: nx, ny, nz, nodeCount, msnode, elemCount, equationNumCount, nftnd0(ntotft), eqNumIndexArrLocTag
    integer (kind = 4) :: i
    if (sizeOfStressDofIndexArr>=(5*sizeOfEqNumIndexArr)) then
        write(*,*) '5*sizeOfEqNumIndexArr',sizeOfEqNumIndexArr,'is not enough for sizeOfStressDofIndexArr',sizeOfStressDofIndexArr
        call abortRun(ERR_MESH_STRESS_ARR_SMALL, &
            'stressArr is too small for this mesh; 5*sizeOfEqNumIndexArr does not cover sizeOfStressDofIndexArr.')
    endif
    if(nodeCount/=nx*ny*nz.or.msnode/=totalNumOfNodes.or.elemCount/=totalNumOfElements.or.equationNumCount/=totalNumOfEquations) then
        write(*,*) 'Inconsistency in node/element/equation/between meshgen and countMeshEntities: stop!',me
        write(*,*) 'nodeCount&totalNumOfNodes=',nodeCount,totalNumOfNodes
        write(*,*) 'elemCount,totalNumOfElements=',elemCount,totalNumOfElements
        write(*,*) 'equationNumCount,totalNumOfEquations=',equationNumCount,totalNumOfEquations
        write(*,*) 'nodeCount,nx,ny,nz',nodeCount,nx,ny,nz
        write(*,*) 'msnode,totalNumOfNodes',msnode,totalNumOfNodes
        call abortRun(ERR_MESH_COUNT_MISMATCH, &
            'meshgen node/element/equation tallies disagree with countMeshEntities (counts printed above).')
    endif
    if(eqNumIndexArrLocTag/=sizeOfEqNumIndexArr) then
        write(*,*) 'Inconsistency in eqNumIndexArrLocTag and sizeOfEqNumIndexArr: stop!',me
        write(*,*) eqNumIndexArrLocTag,sizeOfEqNumIndexArr
        call abortRun(ERR_MESH_EQNUM_MISMATCH, &
            'eqNumIndexArrLocTag /= sizeOfEqNumIndexArr (values printed above).')
    endif
    do i=1,ntotft
        if(nftnd0(i)/=nftnd(i)) then
            write(*,*) 'Inconsistency in fault between meshgen and countMeshEntities:',me,i
            call abortRun(ERR_MESH_FAULT_MISMATCH, &
                'nftnd0 from meshgen disagrees with nftnd from countMeshEntities (rank and fault printed above).')
        endif
    enddo
end subroutine meshGenError

subroutine calcXyzMPIId(mex, mey, mez)
    use globalvar 
    use errorCodes
    implicit none
    integer (kind = 4) :: mex, mey, mez
     
    mex=int(me/(npy*npz))
    mey=int((me-mex*npy*npz)/npz)
    mez=int(me-mex*npy*npz-mey*npz)
end subroutine calcXyzMPIId

subroutine getLocalOneDimCoorArrAndSize(globalOneDimCoorArrSize, numOfNodesWithUniformGridsize, &
    frontEdgeNodeId, MPIXyzId, localOneDimCoorArrSize, localOneDimCoorArr, modelBoundCoor, dimId)
    use globalvar
    use errorCodes
    implicit none 
    integer (kind = 4) :: globalOneDimCoorArrSize, numOfNodesWithUniformGridsize, dimId
    integer (kind = 4) :: frontEdgeNodeId, localOneDimCoorArrSize, MPIXyzId
    integer (kind = 4) :: numOfNodesPerMPI, residualNumOfNodes
    integer (kind = 4) :: i, numOfMPIXyz
    real (kind = dp) :: gridSize, frontEdgeCoor, backEdgeCoor, &
            minCoor, maxCoor, coorTmp, gridSizeTmp, localOneDimCoorArr(10000), &
            modelBoundCoor(3,2)
    real (kind = dp), allocatable :: globalOneDimCoorArr(:)

    if (dimId == 1) then 
        numOfNodesWithUniformGridsize = nint((fltxyz(2,1,1) - fltxyz(1,1,1))/dx) + 1
        gridSize      = dx
        frontEdgeCoor = fltxyz(1,1,1)
        backEdgeCoor  = fltxyz(2,1,1)
        minCoor       = xmin
        maxCoor       = xmax
        numOfMPIXyz   = npx
    elseif (dimId == 2) then 
        numOfNodesWithUniformGridsize = dis4uniF + dis4uniB + 1
        gridSize      = dy
        frontEdgeCoor = -dis4uniF*dy
        backEdgeCoor  = dis4uniB*dy
        minCoor       = ymin
        maxCoor       = ymax
        numOfMPIXyz   = npy
    elseif (dimId == 3) then 
        numOfNodesWithUniformGridsize = nint((fltxyz(2,3,1) - fltxyz(1,3,1))/dz) + 1
        gridSize      = dz
        frontEdgeCoor = fltxyz(1,3,1)
        backEdgeCoor  = fltxyz(2,3,1)
        minCoor       = zmin
        maxCoor       = zmax
        numOfMPIXyz   = npz 
    endif 

    coorTmp = frontEdgeCoor
    gridSizeTmp = gridSize
    do i = 1, np
        gridSizeTmp = gridSizeTmp * rat
        coorTmp = coorTmp - gridSizeTmp
        if (coorTmp <= minCoor) exit
    enddo 
    frontEdgeNodeId = i + nPML
  
    coorTmp = backEdgeCoor
    gridSizeTmp = gridSize
    do i = 1, np
        gridSizeTmp = gridSizeTmp * rat
        coorTmp = coorTmp + gridSizeTmp
        if (coorTmp >= maxCoor) exit
    enddo
    if (dimId == 3) i = -nPML 
    globalOneDimCoorArrSize = numOfNodesWithUniformGridsize + frontEdgeNodeId + i + nPML
    allocate(globalOneDimCoorArr(globalOneDimCoorArrSize))

    numOfNodesPerMPI = int((globalOneDimCoorArrSize+numOfMPIXyz-1)/numOfMPIXyz)
    residualNumOfNodes = (globalOneDimCoorArrSize+numOfMPIXyz-1) - numOfNodesPerMPI * numOfMPIXyz

    if (MPIXyzId<(numOfMPIXyz-residualNumOfNodes)) then 
        localOneDimCoorArrSize = numOfNodesPerMPI
    else 
        localOneDimCoorArrSize = numOfNodesPerMPI + 1
    endif 
    if (localOneDimCoorArrSize > 10000) write(*,*) 'localOneDimCoorArrSize should be < 10000'

    globalOneDimCoorArr(frontEdgeNodeId+1) = frontEdgeCoor
    gridSizeTmp = gridSize
    do i = frontEdgeNodeId, 1, -1
        gridSizeTmp = gridSizeTmp * rat 
        globalOneDimCoorArr(i) = globalOneDimCoorArr(i+1) - gridSizeTmp
    enddo 
    do i = frontEdgeNodeId+2, frontEdgeNodeId + numOfNodesWithUniformGridsize
        globalOneDimCoorArr(i) = globalOneDimCoorArr(i-1) + gridSize
    enddo 
    if (dimId < 3) then 
        gridSizeTmp = gridSize
        do i = frontEdgeNodeId+numOfNodesWithUniformGridsize+1, globalOneDimCoorArrSize
            gridSizeTmp = gridSizeTmp * rat
            globalOneDimCoorArr(i) = globalOneDimCoorArr(i-1) + gridSizeTmp
        enddo 
    endif 
    if (dimId == 1) then 
        modelBoundCoor(1,1) = globalOneDimCoorArr(1)
        modelBoundCoor(1,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)
        PMLb(1) = globalOneDimCoorArr(globalOneDimCoorArrSize-nPML)
        PMLb(2) = globalOneDimCoorArr(nPML+1)
        PMLb(6) = globalOneDimCoorArr(globalOneDimCoorArrSize) - globalOneDimCoorArr(globalOneDimCoorArrSize-1)
    elseif (dimId == 2) then 
        modelBoundCoor(2,1) = globalOneDimCoorArr(1)
        modelBoundCoor(2,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)     
        PMLb(3) = globalOneDimCoorArr(globalOneDimCoorArrSize-nPML)
        PMLb(4) = globalOneDimCoorArr(nPML+1) 
        PMLb(7) = globalOneDimCoorArr(globalOneDimCoorArrSize) - globalOneDimCoorArr(globalOneDimCoorArrSize-1)
    elseif (dimId == 3) then
        modelBoundCoor(3,1) = globalOneDimCoorArr(1) 
        modelBoundCoor(3,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)    
        PMLb(5) = globalOneDimCoorArr(nPML+1)
        PMLb(8) = globalOneDimCoorArr(2) - globalOneDimCoorArr(1)
    endif 

    if (MPIXyzId <= (numOfMPIXyz - residualNumOfNodes)) then 
        do i = 1, localOneDimCoorArrSize
            localOneDimCoorArr(i) = globalOneDimCoorArr((numOfNodesPerMPI-1)*MPIXyzId+i)
        enddo
    else
        do i = 1, localOneDimCoorArrSize
            localOneDimCoorArr(i) = globalOneDimCoorArr((numOfNodesPerMPI-1)*MPIXyzId+i+(MPIXyzId-numOfMPIXyz+residualNumOfNodes))
        enddo
    endif
end subroutine getLocalOneDimCoorArrAndSize

subroutine setNumDof(nodeCoor, numOfDofPerNodeTmp)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: numOfDofPerNodeTmp
    real (kind = dp) :: nodeCoor(10)
    numOfDofPerNodeTmp = ndof !Default
    if (nodeCoor(1)>PMLb(1) .or. nodeCoor(1)<PMLb(2) .or. nodeCoor(2)>PMLb(3) &
        .or. nodeCoor(2)<PMLb(4) .or. nodeCoor(3)<PMLb(5)) then
        numOfDofPerNodeTmp = 12 ! Modify if inside PML
    endif   
end subroutine setNumDof

subroutine setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nodeCount)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: nodeXyzIndex(10), ix, iy, nodeCount, i
    real (kind = dp) :: nodeCoor(10), xline(nodeXyzIndex(4)), yline(nodeXyzIndex(5))
    
    ix = nodeXyzIndex(1)
    iy = nodeXyzIndex(2)
    !Part1. Stations inside the region.
    if(ix>1.and.ix<nodeXyzIndex(4) .and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,totalNumOfOffSt
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(abs(nodeCoor(1)-x4nds(1,i))<tol .or.&
                    (x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
                    (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) .or. &
                    (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
                    (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            numOfOffFaultStCount = numOfOffFaultStCount + 1
                            OffFaultStNodeIdIndex(1,numOfOffFaultStCount) = i
                            OffFaultStNodeIdIndex(2,numOfOffFaultStCount) = nodeCount
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif
    !...identify output nodes (off-fault)
    !Part2. Stations along ix==1
    !Big Bug!!!iz==nz is no valid for 3D MPI.
    if(ix==1.and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,totalNumOfOffSt
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(abs(nodeCoor(1)-x4nds(1,i))<tol .or. &
                    (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
                    (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            numOfOffFaultStCount = numOfOffFaultStCount + 1
                            OffFaultStNodeIdIndex(1,numOfOffFaultStCount) = i
                            OffFaultStNodeIdIndex(2,numOfOffFaultStCount) = nodeCount
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif
    !...identify output nodes (off-fault)
    !Part3. Stations along ix==nx
    if(ix==nodeXyzIndex(4) .and. iy>1.and.iy<nodeXyzIndex(5)) then  !at surface only
        do i=1,totalNumOfOffSt
            if(n4yn(i)==0) then
                if (abs(nodeCoor(3)-x4nds(3,i))<tol) then
                    if(x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
                    (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) then
                        if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                        (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
                        (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                        (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                        (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) then
                            n4yn(i) = 1
                            numOfOffFaultStCount = numOfOffFaultStCount + 1
                            OffFaultStNodeIdIndex(1,numOfOffFaultStCount) = i
                            OffFaultStNodeIdIndex(2,numOfOffFaultStCount) = nodeCount
                            exit     !if node found, jump out the loop
                        endif
                    endif
                endif
            endif
        enddo
    endif    
end subroutine setSurfaceStation

subroutine setEquationNumber(nodeXyzIndex, nodeCoor, eqNumIndexArrLocTag, equationNumCount, numOfDofPerNodeTmp)
    use globalvar 
    use errorCodes
    implicit none
    integer (kind = 4) :: iDof, numOfDofPerNodeTmp, eqNumIndexArrLocTag, equationNumCount, nodeXyzIndex(10)
    real (kind = dp) :: nodeCoor(10)
    
    do iDof = 1,numOfDofPerNodeTmp
        if(abs(nodeCoor(1)-xmin)<tol.or.abs(nodeCoor(1)-xmax)<tol.or.abs(nodeCoor(2)-ymin)<tol &
        .or.abs(nodeCoor(2)-ymax)<tol.or.abs(nodeCoor(3)-zmin)<tol) then
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag+1
            eqNumIndexArr(eqNumIndexArrLocTag) = -1 
            ! Dof = -1 for fixed boundary nodes, no eq #. 
        else
            equationNumCount      = equationNumCount + 1
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag+1
            eqNumIndexArr(eqNumIndexArrLocTag) = equationNumCount
            
            !Count # of DOF on MPI boundaries
            if (nodeXyzIndex(1)==1) then !Left
                numcount(4)=numcount(4)+1
            endif
            if (nodeXyzIndex(1)==nodeXyzIndex(4)) then !Right
                numcount(5)=numcount(5)+1
            endif
            if (nodeXyzIndex(2)==1) then !Front
                numcount(6)=numcount(6)+1
            endif
            if (nodeXyzIndex(2)==nodeXyzIndex(5)) then !Back
                numcount(7)=numcount(7)+1
            endif
            if (nodeXyzIndex(3)==1) then !Down
                numcount(8)=numcount(8)+1
            endif
            if (nodeXyzIndex(3)==nodeXyzIndex(6)) then !Up
                numcount(9)=numcount(9)+1
            endif
        endif
    enddo
    
end subroutine setEquationNumber


subroutine createElement(elemCount, stressDofCount, iy, iz, elementCenterCoor)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: elemCount, stressDofCount, iy, iz, i, j 
    real (kind = dp) :: elementCenterCoor(3)
    
    elementCenterCoor = 0.0d0 
    
    elemCount        = elemCount + 1
    elemTypeArr(elemCount)    = 1 
    nodeElemIdRelation(1,elemCount) = plane1(iy-1,iz-1)
    nodeElemIdRelation(2,elemCount) = plane2(iy-1,iz-1)
    nodeElemIdRelation(3,elemCount) = plane2(iy,iz-1)
    nodeElemIdRelation(4,elemCount) = plane1(iy,iz-1)
    nodeElemIdRelation(5,elemCount) = plane1(iy-1,iz)
    nodeElemIdRelation(6,elemCount) = plane2(iy-1,iz)
    nodeElemIdRelation(7,elemCount) = plane2(iy,iz)
    nodeElemIdRelation(8,elemCount) = plane1(iy,iz)
    stressCompIndexArr(elemCount)   = stressDofCount
    
    do i=1,nen
        do j=1,3
            elementCenterCoor(j) = elementCenterCoor(j) + meshCoor(j,nodeElemIdRelation(i,elemCount))
        enddo
    enddo
    elementCenterCoor = elementCenterCoor/8.0d0
    
    if (elementCenterCoor(1)>PMLb(1).or.elementCenterCoor(1)<PMLb(2) &
    .or.elementCenterCoor(2)>PMLb(3).or.elementCenterCoor(2)<PMLb(4) &
    .or.elementCenterCoor(3)<PMLb(5)) then
        elemTypeArr(elemCount) = 2
        stressDofCount        = stressDofCount+15+6
    else
        stressDofCount        = stressDofCount+12
    endif
    
    ! assign nodeElemIdRelation(elemCount), ids(elemCount), et(elemCount)
    ! return elementCenterCoor, 
    ! update elemCount, stressDofCount
end subroutine createElement

 subroutine replaceSlaveWithMasterNode(nodeCoor, elemCount, nftnd0)
    use globalvar 
    use errorCodes
    implicit none
    integer (kind = 4) :: elemCount, iFault, iFaultNodePair, nftnd0(ntotft), k
    real (kind = dp) :: nodeCoor(10)
    ! The default grids only contain slave nodes. 
    ! This subroutine will replace slave nodes with corresponding master nodes.
    
    if ((elemTypeArr(elemCount) == 1 .and. &
        (nodeCoor(1)>(fltxyz(1,1,1)-tol) .and. nodeCoor(1)<(fltxyz(2,1,1)+dx+tol) .and. &
         nodeCoor(3)>(fltxyz(1,3,1)-tol) .and. nodeCoor(2)>0.0d0 .and. abs(nodeCoor(2)-dy)<tol)) &
         .or. elemTypeArr(elemCount)==12 .or. elemTypeArr(elemCount)==13 ) then
        do iFault = 1, ntotft
            do iFaultNodePair = 1, nftnd0(iFault)
                do k = 1,nen
                    if(nodeElemIdRelation(k,elemCount)==nsmp(1,iFaultNodePair,iFault)) then
                        nodeElemIdRelation(k,elemCount) = nsmp(2,iFaultNodePair,iFault)  !use master node for the node!
                    endif
                enddo
            enddo
        enddo
    endif      
    
end subroutine replaceSlaveWithMasterNode

subroutine checkIsOnFault(nodeCoor, iFault, isOnFault)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: isOnFault, iFault
    real (kind = dp) :: nodeCoor(10), distToFault
    isOnFault = 0
 
    if(nodeCoor(1)>=(fltxyz(1,1,iFault)-tol).and.nodeCoor(1)<=(fltxyz(2,1,iFault)+tol).and. &
        nodeCoor(2)>=(fltxyz(1,2,iFault)-tol).and.nodeCoor(2)<=(fltxyz(2,2,iFault)+tol).and. &
        nodeCoor(3)>=(fltxyz(1,3,iFault)-tol) .and. nodeCoor(3)<=(fltxyz(2,3,iFault)+tol)) then
        if (C_degen==0.0d0 .and. nodeCoor(2)==0.0d0) then 
            isOnFault = 1
        elseif (C_degen>3.) then
            if (fltxyz(1,2,iFault)>=fltxyz(2,2,iFault)) write(*,*) 'ymax should be > ymin. Wrong geo, exit'
            distToFault = abs(nodeCoor(3)+nodeCoor(2)*dtan(C_degen/180.d0*pi)) 
            distToFault = distToFault/(1.d0+dtan(C_degen/180.d0*pi)**2)**0.5
            if (distToFault < dx/100.d0) isOnFault = 1   
    endif 
    endif
end subroutine checkIsOnFault

subroutine createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nodeCount, msnode, nftnd0, equationNumCount, eqNumIndexArrLocTag,&
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
use globalvar 
use errorCodes
implicit none
integer (kind = 4) :: iFault, iFaultNodePair, isOnFault, nodeCount, msnode, nftnd0(ntotft), equationNumCount, i, nxuni, nzuni, eqNumIndexArrLocTag
integer (kind = 4) :: fltrc(2,nxuni,nzuni,ntotft), ixfi(ntotft), izfi(ntotft), ifs(ntotft), ifd(ntotft), nodeXyzIndex(10)
real (kind = dp) :: nodeCoor(10), ycoort, pfx, pfz

do iFault = 1, ntotft

    isOnFault = 0 
    call checkIsOnFault(nodeCoor, iFault, isOnFault)

    if (isOnFault == 1) then
        nftnd0(iFault)                = nftnd0(iFault) + 1 ! # of split-node pairs + 1
        nsmp(1,nftnd0(iFault),iFault) = nodeCount              ! set Slave node nodeID to nsmp  
        msnode                        = nodeXyzIndex(4)*nodeXyzIndex(5)*nodeXyzIndex(6) + nftnd0(iFault) ! create Master node at the end of regular grids
        if (iFault>1) then
            write(*,*) 'meshgen: got iFault =', iFault
            call abortRun(ERR_MESH_MULTIFAULT_MSNODE, &
                'master-node construction cannot handle more than one fault (ntotft>1).')
        endif
        
        eqNumStartIndexLoc(msnode) = eqNumIndexArrLocTag
        numOfDofPerNodeArr(msnode) = 3         
        nsmp(2,nftnd0(iFault),iFault) = msnode !set Master node nodeID to nsmp
        plane2(nodeXyzIndex(5)+iFault,nodeXyzIndex(3)) = msnode
        
        meshCoor(1,msnode) = nodeCoor(1)
        meshCoor(2,msnode) = nodeCoor(2)
        meshCoor(3,msnode) = nodeCoor(3)
        if (insertFaultType > 0) then 
            meshCoor(2,msnode) = ycoort
        endif
        
        !set Equation Numbers for the newly created Master Node.
        do i = 1, ndof
            equationNumCount      = equationNumCount + 1
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag + 1
            eqNumIndexArr(eqNumIndexArrLocTag) = equationNumCount
        enddo
        
        ! Count split-node pair # for MPI
        ! fltgm, fltnum are global vars.
        if(nodeXyzIndex(1) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 1
            fltnum(1) = fltnum(1) + 1
        endif
        if(nodeXyzIndex(1) == nodeXyzIndex(4)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 2
            fltnum(2) = fltnum(2) + 1
        endif
        if(nodeXyzIndex(2) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 10
            fltnum(3) = fltnum(3) + 1
        endif
        if(nodeXyzIndex(2) == nodeXyzIndex(5)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 20
            fltnum(4) = fltnum(4) + 1
        endif
        if(nodeXyzIndex(3) == 1) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 100
            fltnum(5) = fltnum(5) + 1
        endif
        if(nodeXyzIndex(3) == nodeXyzIndex(6)) then
            fltgm(nftnd0(iFault)) = fltgm(nftnd0(iFault)) + 200
            fltnum(6) = fltnum(6) + 1
        endif             

        ! setOnFaultStation
        do i = 1, nonfs(iFault)
            if(abs(nodeCoor(1)-xonfs(1,i,iFault))<tol .and. &
                abs(nodeCoor(3)-xonfs(2,i,iFault))<tol) then
                numOfOnFaultStCount = numOfOnFaultStCount + 1
                anonfs(1,numOfOnFaultStCount) = nftnd0(iFault)
                anonfs(2,numOfOnFaultStCount) = i
                anonfs(3,numOfOnFaultStCount) = iFault
                exit
            endif
        enddo  

        ! set unit vectors to split-node pair    
        un(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(2,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(3,nftnd0(iFault),iFault) = dcos(fltxyz(2,4,iFault))        
        us(1,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))
        us(2,nftnd0(iFault),iFault) = -dcos(fltxyz(1,4,iFault))
        us(3,nftnd0(iFault),iFault) = 0.0d0
        ud(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(2,nftnd0(iFault),iFault) = dsin(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(3,nftnd0(iFault),iFault) = dsin(fltxyz(2,4,iFault))
        
        if (insertFaultType >0) then
            un(1,nftnd0(iFault),iFault) = -pfx/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(2,nftnd0(iFault),iFault) = 1.0d0/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(3,nftnd0(iFault),iFault) = -pfz/(pfx**2 + 1.0d0 + pfz**2)**0.5    
            us(1,nftnd0(iFault),iFault) = 1.0d0/(1.0d0 + pfx**2)**0.5
            us(2,nftnd0(iFault),iFault) = pfx/(1.0d0 + pfx**2)**0.5
            us(3,nftnd0(iFault),iFault) = 0.0d0
            ud(1,nftnd0(iFault),iFault) = us(2,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault) &
                - us(3,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault)
            ud(2,nftnd0(iFault),iFault) = us(3,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault) &
                - us(1,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault)
            ud(3,nftnd0(iFault),iFault) = us(1,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault) &
                - us(2,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault)
        endif                             
        
        !...prepare for area calculation
        if(ixfi(iFault)==0) ixfi(iFault)=nodeXyzIndex(1)
        if(izfi(iFault)==0) izfi(iFault)=nodeXyzIndex(3)
        ifs(iFault)=nodeXyzIndex(1)-ixfi(iFault)+1
        ifd(iFault)=nodeXyzIndex(3)-izfi(iFault)+1
        fltrc(1,ifs(iFault),ifd(iFault),iFault) = msnode    !master node
        fltrc(2,ifs(iFault),ifd(iFault),iFault) = nftnd0(iFault) !fault node num in sequence
    endif 
enddo 

end subroutine createMasterNode

subroutine createNode(nodeCoor, xcoor, ycoor, zcoor, nodeCount, nodeXyzIndex)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: nodeCount, nodeXyzIndex(10), iy, iz
    real (kind = dp) :: nodeCoor(10), xcoor, ycoor, zcoor
    iy = nodeXyzIndex(2)
    iz = nodeXyzIndex(3)
    nodeCoor(1)   = xcoor
    nodeCoor(2)   = ycoor
    nodeCoor(3)   = zcoor
    nodeCount         = nodeCount + 1        
    plane2(iy,iz) = nodeCount
    meshCoor(1,nodeCount) = nodeCoor(1)
    meshCoor(2,nodeCount) = nodeCoor(2)
    meshCoor(3,nodeCount) = nodeCoor(3) 
end subroutine createNode

! insertFaultInterface moved to func_lib.f90.

subroutine initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
    use globalvar 
    use errorCodes
    implicit none 
    integer (kind = 4) :: ix, iy, iz, nx, ny, nz, nodeXyzIndex(10)
    nodeXyzIndex(1) = ix
    nodeXyzIndex(2) = iy
    nodeXyzIndex(3) = iz
    nodeXyzIndex(4) = nx
    nodeXyzIndex(5) = ny
    nodeXyzIndex(6) = nz
end subroutine initializeNodeXyzIndex

subroutine setPlasticStress(depth, elemCount)
    use globalvar
    use errorCodes
    implicit none
    
    real(kind = dp) :: depth, vstmp, vptmp, routmp, strVert, devStr
    integer(kind = 4) :: elemCount, etTag
    
    etTag = 0
    if (elemTypeArr(elemCount)==2) etTag = 1 ! adjustment for PML elements
    
    eleporep(elemCount) = 0.0d0  !rhow*tmp2*gama  !pore pressure>0
    strVert            = -(roumax- rhow*(gamar+1.0d0))*depth*grav ! should be negative   
    devStr             = abs(strVert)*devStrToStrVertRatio ! positive
    
    stressArr(stressCompIndexArr(elemCount)+3+15*etTag) = strVert
    stressArr(stressCompIndexArr(elemCount)+1+15*etTag) = strVert - devStr*dcos(2.0d0*str1ToFaultAngle)
    stressArr(stressCompIndexArr(elemCount)+2+15*etTag) = strVert + devStr*dcos(2.0d0*str1ToFaultAngle)
    stressArr(stressCompIndexArr(elemCount)+6+15*etTag) = devStr*dsin(2.0d0*str1ToFaultAngle)
    if (stressArr(stressCompIndexArr(elemCount)+2+15*etTag) >= 0.0d0) write(*,*) 'WARNING: positive Sigma3 ... ...'
    
end subroutine setPlasticStress
