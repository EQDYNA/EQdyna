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
    real (kind = dp)::xcoor,ycoor,zcoor,xstep,ystep,zstep,mdx(3)
    real (kind = dp),allocatable,dimension(:)::xlinet,ylinet,zlinet,xline,yline,zline

    mex=int(me/(npy*npz))
    mey=int((me-mex*npy*npz)/npz)
    mez=int(me-mex*npy*npz-mey*npz)
    !...determine num of nodes along x
    nxuni = (fltxyz(2,1,1) - fltxyz(1,1,1)) / dx + 1
    xstep = dx
    xcoor = fltxyz(1,1,1)
    do ix = 1, np
        xstep = xstep * rat
        xcoor = xcoor - xstep
        if(xcoor <= xmin) exit
    enddo
    edgex1 = ix + nPML
    xstep = dx
    xcoor = fltxyz(2,1,1)
    do ix = 1, np
        xstep = xstep * rat
        xcoor = xcoor + xstep
        if(xcoor >= xmax) exit
    enddo
    nxt = nxuni + edgex1 + ix + nPML
    allocate(xlinet(nxt))
    !
    !...predetermine x-coor
    xlinet(edgex1+1) = fltxyz(1,1,1)
    xstep = dx
    do ix = edgex1, 1, -1
        xstep = xstep * rat
        xlinet(ix) = xlinet(ix+1) - xstep
    enddo
    xmin1=xlinet(1)
    do ix = edgex1+2,edgex1+nxuni
        xlinet(ix) = xlinet(ix-1) + dx
    enddo
    xstep = dx
    do ix = edgex1+nxuni+1,nxt
        xstep = xstep * rat
        xlinet(ix) = xlinet(ix-1) + xstep
    enddo
    xmax1=xlinet(nxt)

    j1 = nxt + npx - 1
    rlp = j1/npx
    rr = j1 - rlp*npx
    if(mex<(npx-rr)) then
        nx = rlp
    else
        nx = rlp + 1    !evenly distributed to last rr
    endif
    allocate(xline(nx))
    if(mex<=npx-rr) then
        do ix=1,nx
            xline(ix) = xlinet(rlp*mex-mex+ix)
        enddo
    else
        do ix=1,nx
            k1 = mex - (npx - rr)
            xline(ix) = xlinet(rlp*mex-mex+k1+ix)
        enddo
    endif

    nyuni = dis4uniF + dis4uniB + 1
    ystep = dy
    !ycoor = -dy*(dis4uniF+fltxyz(2,1,1)/dx+1)
    ycoor = -dy*(dis4uniF)
    do iy = 1, np
        ystep = ystep * rat
        ycoor = ycoor - ystep
        if(ycoor <= ymin) exit
    enddo
    edgey1 = iy + nPML
    ystep = dy
    ycoor = dy*(dis4uniB)
    do iy = 1, np
        ystep = ystep * rat
        ycoor = ycoor + ystep
        if(ycoor >= ymax) exit
    enddo
    nyt = nyuni + edgey1 + iy + nPML
    !...pre-determine y-coor
    allocate(ylinet(nyt))
    !...predetermine y-coor
    ylinet(edgey1+1) = -dy*(dis4uniF)
    ystep = dy
    do iy = edgey1, 1, -1
        ystep = ystep * rat
        ylinet(iy) = ylinet(iy+1) - ystep
    enddo
    ymin1=ylinet(1)
    do iy = edgey1+2,edgey1+nyuni
        ylinet(iy) = ylinet(iy-1) + dy
    enddo
    ystep = dy
    do iy = edgey1+nyuni+1,nyt
        ystep = ystep * rat
        ylinet(iy) = ylinet(iy-1) + ystep
    enddo
    ymax1=ylinet(nyt)

      j1 = nyt + npy - 1
      rlp = j1/npy
      rr = j1 - rlp*npy
      if(mey<(npy-rr)) then
        ny = rlp
      else
        ny = rlp + 1    !evenly distributed to last rr
      endif
      allocate(yline(ny))
      if(mey<=npy-rr) then
        do iy=1,ny
          yline(iy) = ylinet(rlp*mey-mey+iy)
        enddo
      else
        do iy=1,ny
          k1 = mey - (npy - rr)
          yline(iy) = ylinet(rlp*mey-mey+k1+iy)
        enddo
      endif
    !
    !...determine num of nodes along z
    zstep = dz
    zcoor = fltxyz(1,3,1)
    do iz=1,np
        zstep = zstep * rat
        zcoor = zcoor - zstep
        if(zcoor <= zmin) exit
    enddo
    edgezn = iz + nPML
    nzuni = (fltxyz(2,3,1)-fltxyz(1,3,1))/dz + 1 
    nzt = edgezn + nzuni
    !...predetermine z-coor

    allocate(zlinet(nzt))

    zlinet(nzt) = zmax
    do iz = nzt-1,nzt-nzuni+1,-1
        zlinet(iz) = zlinet(iz+1) - dz
    enddo
    zstep = dz
    do iz = nzt-nzuni,1,-1
        zstep = zstep * rat
        zlinet(iz) = zlinet(iz+1) -zstep
    enddo
    zmin1=zlinet(1)

    j1 = nzt + npz - 1
    rlp = j1/npz
    rr = j1 - rlp*npz

    if(mez<(npz-rr)) then
        nz = rlp
    else
        nz = rlp + 1    !evenly distributed to last rr
    endif

    allocate(zline(nz))

    if(mez<=npz-rr) then
        do iz=1,nz
          zline(iz) = zlinet(rlp*mez-mez+iz)
        enddo
    else
        do iz=1,nz
          k1 = mez - (npz - rr)
          zline(iz) = zlinet(rlp*mez-mez+k1+iz)
        enddo
    endif
      
    PMLb(1)=xlinet(nxt-nPML)!PMLxmax
    PMLb(2)=xlinet(nPML+1)!PMLxmin
    PMLb(3)=ylinet(nyt-nPML)!PMLymax
    PMLb(4)=ylinet(nPML+1)!PMLymin
    PMLb(5)=zlinet(nPML+1)!PMLzmin
    mdx(1)=xlinet(nxt)-xlinet(nxt-1)
    mdx(2)=ylinet(nyt)-ylinet(nyt-1)
    mdx(3)=zlinet(2)-zlinet(1)
    PMLb(6)=mdx(1)
    PMLb(7)=mdx(2)
    PMLb(8)=mdx(3)

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
                    if(abs(xcoor-xmin1)<tol.or.abs(xcoor-xmax1)<tol.or.abs(ycoor-ymin1)<tol &
                        .or.abs(ycoor-ymax1)<tol.or.abs(zcoor-zmin1)<tol) then
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
