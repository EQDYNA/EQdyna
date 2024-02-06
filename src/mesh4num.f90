! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine mesh4num
    use globalvar
    implicit none
    include 'mpif.h'

    logical,dimension(ntotft) :: ynft
    integer(kind = 4) :: nodeCount=0, elementCount=0, equationNumCount=0, &
            nxt, nyt, nzt, nx, ny, nz, ix, iy, iz, &
        edgex1,edgey1, iDof,edgezn, eqNumIndexArrSizeCount=0,numOfDof, nxuni,nyuni,nzuni,ift,mex,mey,mez
    real (kind = dp) :: xcoor, ycoor, zcoor, xline(10000), yline(10000), zline(10000), modelBoundCoor(3,2)

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
                numOfDof=ndof

                if (xcoor>PMLb(1).or.xcoor<PMLb(2).or.ycoor>PMLb(3).or.ycoor<PMLb(4) &
                    .or.zcoor<PMLb(5)) then
                    numOfDof = 12
                endif            

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
                    ynft(ift) = .false.
                    if(xcoor>=(fltxyz(1,1,ift)-tol).and.xcoor<=(fltxyz(2,1,ift)+tol).and. &
                        ycoor>=(fltxyz(1,2,ift)-tol).and.ycoor<=(fltxyz(2,2,ift)+tol).and. &
                        zcoor>=(fltxyz(1,3,ift)-tol) .and. zcoor<=(fltxyz(2,3,ift)+tol)) then
                        if(ift==1) then
                            ynft(ift) = .true.
                        elseif(ift==2.and.abs(xcoor*dtan(brangle)-abs(ycoor))<tol) then
                            ynft(ift)=.true.
                        endif
                        if(ynft(ift)) then
                            nftnd(ift) = nftnd(ift) + 1
                            nodeCount = nodeCount + 1
                            !...establish equation numbers for this master node
                            do iDof=1,ndof
                                equationNumCount = equationNumCount + 1
                                eqNumIndexArrSizeCount=eqNumIndexArrSizeCount+1!DL
                            enddo
                            exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
                        endif  !if ynft(ift)
                    endif  !if flt range
                enddo  !do ift 
                
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                
                    elementCount = elementCount + 1            
                    
                    if (C_degen == 1) then 
                        call wedge4num(xcoor-dx/2.0d0, ycoor-dx/2.0d0, zcoor-dx/2.0d0, elementCount)
                    elseif (C_degen == 2) then 
                        call tetra4num(xcoor-dx/2.0d0, ycoor-dx/2.0d0, zcoor-dx/2.0d0, elementCount)
                    endif                 
                endif  !if element
            enddo    !iy
        enddo    !iz
    enddo        !ix
    sizeOfEqNumIndexArr  = eqNumIndexArrSizeCount
    totalNumOfNodes = nodeCount
    totalNumOfElements = elementCount
    totalNumOfEquations = equationNumCount
end subroutine mesh4num
