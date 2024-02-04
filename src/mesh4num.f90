! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine mesh4num
    use globalvar
    implicit none
    include 'mpif.h'
    
    integer(kind = 4)::nnode,nelement,nxt,nyt,nzt,nx,ny,nz,ix,iy,iz,&
        edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,edgezn
    logical,dimension(ntotft)::ynft
    integer(kind=4)::ntag,zz,nsx,nsy
    integer(kind=4)::nxuni,nyuni,nzuni,ift,n1,n2,n3,n4,m1,m2,m3,m4,rlp,rr,mex,mey,mez,msnode
    real (kind = dp)::xcoor,ycoor,zcoor,xstep,ystep,zstep
    real (kind = dp) :: xline(10000), yline(10000), zline(10000), modelBoundCoor(3,2)

    call calcXyzMPIId(mex, mey, mez)
    call getLocalOneDimCoorArrAndSize(nxt, nxuni, edgex1, mex, nx, xline, modelBoundCoor, 1)
    call getLocalOneDimCoorArrAndSize(nyt, nyuni, edgey1, mey, ny, yline, modelBoundCoor, 2)
    call getLocalOneDimCoorArrAndSize(nzt, nzuni, edgezn, mez, nz, zline, modelBoundCoor, 3)

    nnode = 0
    nelement = 0
    neq = 0
    nftnd = 0
    ntag=0

    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                xcoor = xline(ix)
                ycoor = yline(iy)
                zcoor = zline(iz)    

                nnode = nnode + 1
                zz=ndof

                if (xcoor>PMLb(1).or.xcoor<PMLb(2).or.ycoor>PMLb(3).or.ycoor<PMLb(4) &
                    .or.zcoor<PMLb(5)) then
                    zz = 12
                endif            

                do i1=1,zz
                    if(abs(xcoor-modelBoundCoor(1,1))<tol.or.abs(xcoor-modelBoundCoor(1,2))<tol.or.abs(ycoor-modelBoundCoor(2,1))<tol &
                        .or.abs(ycoor-modelBoundCoor(2,2))<tol.or.abs(zcoor-modelBoundCoor(3,1))<tol) then
                        i = -1  !-1 for fixed boundary nodes
                        ntag=ntag+1
                    else
                        neq = neq + 1
                        ntag=ntag+1                    
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
                            nnode = nnode + 1
                            !...establish equation numbers for this master node
                            do i1=1,ndof
                                neq = neq + 1
                                ntag=ntag+1!DL
                            enddo
                            exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
                        endif  !if ynft(ift)
                    endif  !if flt range
                enddo  !do ift 
                
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                
                    nelement = nelement + 1            
                    
                    if (C_degen == 1) then 
                        call wedge4num(xcoor-dx/2.0d0, ycoor-dx/2.0d0, zcoor-dx/2.0d0, nelement)
                    elseif (C_degen == 2) then 
                        call tetra4num(xcoor-dx/2.0d0, ycoor-dx/2.0d0, zcoor-dx/2.0d0, nelement)
                    endif                 
                endif  !if element
            enddo    !iy
        enddo    !iz
    enddo        !ix
    maxm  = ntag
    numnp = nnode
    numel = nelement
end subroutine mesh4num
