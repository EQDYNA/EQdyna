! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine countMeshEntities
    use globalvar
    implicit none
    include 'mpif.h'

    integer(kind = 4) :: nodeCount=0, elementCount=0, equationNumCount=0, &
            nxt, nyt, nzt, nx, ny, nz, ix, iy, iz, &
        edgex1,edgey1, iDof,edgezn, eqNumIndexArrSizeCount=0,numOfDof, nxuni,nyuni,nzuni,ift,mex,mey,mez,isOnFt
    real (kind = dp) :: xcoor, ycoor, zcoor, xline(10000), yline(10000), zline(10000), modelBoundCoor(3,2), nodeCoor(10)

    call calcXyzMPIId(mex, mey, mez)
    call getLocalOneDimCoorArrAndSize(nxt, nxuni, edgex1, mex, nx, xline, modelBoundCoor, 1)
    call getLocalOneDimCoorArrAndSize(nyt, nyuni, edgey1, mey, ny, yline, modelBoundCoor, 2)
    call getLocalOneDimCoorArrAndSize(nzt, nzuni, edgezn, mez, nz, zline, modelBoundCoor, 3)

    nftnd = 0

    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                xcoor = xline(ix)
                ycoor = yline(iy)
                zcoor = zline(iz)    

                nodeCount = nodeCount + 1
                nodeCoor(1) = xcoor
                nodeCoor(2) = ycoor
                nodeCoor(3) = zcoor
                call setNumDof(nodeCoor, numOfDof)

                do iDof=1,numOfDof
                    if(abs(xcoor-modelBoundCoor(1,1))<tol.or.abs(xcoor-modelBoundCoor(1,2))<tol.or.abs(ycoor-modelBoundCoor(2,1))<tol &
                        .or.abs(ycoor-modelBoundCoor(2,2))<tol.or.abs(zcoor-modelBoundCoor(3,1))<tol) then
                        ! -1 for fixed boundary nodes; no equation number needed.
                        eqNumIndexArrSizeCount=eqNumIndexArrSizeCount+1
                    else
                        equationNumCount = equationNumCount + 1
                        eqNumIndexArrSizeCount=eqNumIndexArrSizeCount+1                    
                    endif
                enddo

                do ift=1,ntotft
                    isOnFt = 0
                    call checkIsOnFault(nodeCoor, ift, isOnFt)
                    if(isOnFt==1) then
                        nftnd(ift) = nftnd(ift) + 1
                        nodeCount = nodeCount + 1
                        !...establish equation numbers for this master node
                        do iDof=1,ndof
                            equationNumCount = equationNumCount + 1
                            eqNumIndexArrSizeCount=eqNumIndexArrSizeCount+1!DL
                        enddo
                        exit !can only be on 1 fault, thus if isOnFt==1, exit do loop
                    endif  !if isOnFt
                enddo  !do ift
                
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                
                    elementCount = elementCount + 1            
                    
                    if (C_degen > 3.d0) then 
                        call wedge4num(xcoor-dx/2.0d0, ycoor-dy/2.0d0, zcoor-dz/2.0d0, elementCount)
                    endif                 
                endif 

            enddo   
        enddo   
    enddo       

    sizeOfEqNumIndexArr  = eqNumIndexArrSizeCount
    totalNumOfNodes = nodeCount
    totalNumOfElements = elementCount
    totalNumOfEquations = equationNumCount
end subroutine countMeshEntities
