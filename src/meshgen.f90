!/* Copyright (C) 2006-2023, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
subroutine meshgen

    use globalvar
    implicit none
    include 'mpif.h'

    integer(kind=4) :: nnode, nelement, neq0, &
                       nxt,nyt,nzt,nx,ny,nz,ix,iy,iz, &
                       edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,edgezn, &
                       nxuni,nyuni,nzuni,ift, &
                       n1,n2,n3,n4,m1,m2,m3,m4,alloc_err,&
                       istatus(MPI_STATUS_SIZE), &
                       bndl,bndr,bndf,bndb,bndd,bndu, &
                       mex,mey,mez,rlp,rr,ierr,jj,itmp1,&
                       ntag,ntags,zz,ivp1,ivp2,ixe,iye,ize, &
                       itemp,iye1,nxe,nye,nze,&
                       nsx,nsy,nfx,nfz,msnode
    integer(kind=4),dimension(ntotft) :: nftnd0,ixfi,izfi,ifs,ifd
    integer(kind=4),allocatable :: fltrc(:,:,:,:)

    real(kind = dp) :: tol,xcoor,ycoor,zcoor, &
                       xstep,ystep,zstep,tmp1,tmp2,tmp3,&
                       a,b,area,aa1,bb1,cc1,dd1,p1,q1
    real(kind = dp),allocatable,dimension(:) :: xlinet,ylinet,zlinet, &
                       xline,yline,zline,btmp,btmp1,btmp2,btmp3
    real(kind = dp),allocatable::vpstruct(:,:)
    !Heteogeneous stress setup
    real(kind = dp)::initialinput(3,91001),xc(3)
    real (kind = dp)::omega
    logical,dimension(ntotft) :: ynft
    character (len=30)  ::meshfile
    character (len=100) ::fname
    real (kind = dp) :: ycoort, pfx = 0.0d0, pfz = 0.0d0
    !-------------------------------------------------------------------!
    !------------Early Oct.2015/ D.Liu----------------------------------!
    !------------Input of initial stress for TSN------------------------!
    ! open(3001,file='0InputFEM.txt',form='formatted',status='old')
    ! do i=1,91001
    ! read(3001,*) (initialinput(j,i),j=1,3)
    ! enddo
    ! close(3001)
    !-------------------------------------------------------------------!
    allocate(n4yn(n4nds))
    write(mm,'(i6)') me
    mm = trim(adjustl(mm))

    n4yn=0
    dy=dx
    dz=dx
    tol=dx/100.0d0

    !Prepare for 3D MPI partitioning
    !Search Tag: 3DMPI
    !2/5/2016/L.Bin
    !Modified by D.Liu
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
    xmin = xlinet(1)
    do ix = edgex1+2,edgex1+nxuni
        xlinet(ix) = xlinet(ix-1) + dx
    enddo
    xstep = dx
    do ix = edgex1+nxuni+1,nxt
        xstep = xstep * rat
        xlinet(ix) = xlinet(ix-1) + xstep
    enddo
    xmax = xlinet(nxt)
    !
    !***********MPI***************  
    !...MPI partitioning based on ix and xlinet. B.D. 4/18/09
    !...need overlap between adjacent procs. B.D. 5/10/09
    !...more evenly distributed. B.D. 11/1/09
    !...due to overlap, total line num should be: nxt+nprocs-1
    !    B.D. 1/19/10
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
    nyt = nyuni + edgey1 + iy + nPML!3DMPI
    !...pre-determine y-coor
    allocate(ylinet(nyt))
    !...predetermine y-coor
    ylinet(edgey1+1) = -dy*(dis4uniF)
    ystep = dy
    do iy = edgey1, 1, -1
        ystep = ystep * rat    
        ylinet(iy) = ylinet(iy+1) - ystep
    enddo
    ymin = ylinet(1)
    do iy = edgey1+2,edgey1+nyuni
        ylinet(iy) = ylinet(iy-1) + dy
    enddo
    ystep = dy
    do iy = edgey1+nyuni+1,nyt
        ystep = ystep * rat
        ylinet(iy) = ylinet(iy-1) + ystep
    enddo
    ymax = ylinet(nyt)
    !***********MPI***************  
    !...MPI partitioning based on iy and ylinet.
    !...need overlap between adjacent procs.
    !...more evenly distributed. 
    !...due to overlap, total line num should be: nyt+npy-1
    !Search Tag: 3DMPI
    !2/5/2016/L.Bin
    !Modified by D.Liu
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
    nzuni = (fltxyz(2,3,1)-fltxyz(1,3,1))/dx + 1 
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
    zmin = zlinet(1) !these needed for fixed boundaryies
    !***********MPI***************  
    !...MPI partitioning based on iz and zlinet.
    !...need overlap between adjacent procs.
    !...more evenly distributed. 
    !...due to overlap, total line num should be: nzt+npz-1
    !Search Tag: 3DMPI
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

    !...prepare for digitizing
    allocate(plane1(ny+ntotft,nz),plane2(ny+ntotft,nz),fltrc(2,nxuni,nzuni,ntotft))
    nnode=0
    !3DMPI
    msnode=nx*ny*nz
    numcount=0
    numcount(1)=nx
    numcount(2)=ny
    numcount(3)=nz
    nelement = 0
    n4onf = 0
    n4out = 0
    plane1 = 0
    plane2 = 0
    neq0 = 0
    nftnd0 = 0
    an4nds = 0
    x = 0.0
    ien = 0
    ixfi = 0
    izfi = 0
    id1= 0
    ids= 0
    locid=0
    dof1=0
    et=0
    ntag = 0
    ntags = 0
    !  
    !...digitize along constant x plane (normal to fault strike for MPI)
    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                xcoor = xline(ix)
                ycoor = yline(iy)
                zcoor = zline(iz)
                
                !...create nodes
                nnode = nnode + 1        
                plane2(iy,iz) = nnode
                x(1,nnode) = xcoor
                x(2,nnode) = ycoor
                x(3,nnode) = zcoor
                if (rough_fault == 1) then 
                    call insert_rough_fault(xcoor, ycoor, zcoor, ycoort, pfx, pfz)
                    x(2,nnode) = ycoort
                endif 
                
                zz=ndof
                if (xcoor>PMLb(1).or.xcoor<PMLb(2).or.ycoor>PMLb(3).or.ycoor<PMLb(4) &
                .or.zcoor<PMLb(5)) then
                    zz = 12
                endif        
                locid(nnode)=ntag
                dof1(nnode)=zz
                !...establish equation numbers for this node
                do i1=1,zz
                !...elastoplastic off-fault, stress assigned in entire model,
                ! need to fix model boundaries (except free surface).
                    if(abs(xcoor-xmin)<tol.or.abs(xcoor-xmax)<tol.or.abs(ycoor-ymin)<tol &
                    .or.abs(ycoor-ymax)<tol.or.abs(zcoor-zmin)<tol) then
                        ntag=ntag+1
                        id1(ntag)=-1!-1 for fixed boundary nodes
                    else
                        neq0 = neq0 + 1
                        ntag=ntag+1
                        id1(ntag)=neq0
                        !Count numbers of d.o.f on interfaces between MPIs
                        if (ix==1) then !Left
                            numcount(4)=numcount(4)+1
                        endif
                        if (ix==nx) then !Right
                            numcount(5)=numcount(5)+1
                        endif
                        if (iy==1) then !Front
                            numcount(6)=numcount(6)+1
                        endif
                        if (iy==ny) then !Back
                            numcount(7)=numcount(7)+1
                        endif
                        if (iz==1) then !Down
                            numcount(8)=numcount(8)+1
                        endif
                        if (iz==nz) then !Up
                            numcount(9)=numcount(9)+1
                        endif
                    endif
                enddo ! enddo i1=1,ndof
                !...identify output nodes (off-fault)
                !Part1. Stations inside the region.
                if(ix>1.and.ix<nx .and. iy>1.and.iy<ny) then  !at surface only
                    do i=1,n4nds
                        if(n4yn(i)==0) then
                            if (abs(zcoor-x4nds(3,i))<tol) then
                                if(abs(xcoor-x4nds(1,i))<tol .or.&
                                (x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<xcoor.and. &
                                (xcoor-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) .or. &
                                (x4nds(1,i)>xcoor.and.x4nds(1,i)<xline(ix+1).and. &
                                (x4nds(1,i)-xcoor)<(xline(ix+1)-x4nds(1,i)))) then
                                    if(abs(ycoor-x4nds(2,i))<tol .or. &
                                    (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<ycoor.and. &
                                    (ycoor-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                                    (x4nds(2,i)>ycoor.and.x4nds(2,i)<yline(iy+1).and. &
                                    (x4nds(2,i)-ycoor)<(yline(iy+1)-x4nds(2,i)))) then
                                        n4yn(i) = 1
                                        n4out = n4out + 1
                                        an4nds(1,n4out) = i
                                        an4nds(2,n4out) = nnode
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
                if(ix==1.and. iy>1.and.iy<ny) then  !at surface only
                    do i=1,n4nds
                        if(n4yn(i)==0) then
                            if (abs(zcoor-x4nds(3,i))<tol) then
                                if(abs(xcoor-x4nds(1,i))<tol .or. &
                                (x4nds(1,i)>xcoor.and.x4nds(1,i)<xline(ix+1).and. &
                                (x4nds(1,i)-xcoor)<(xline(ix+1)-x4nds(1,i)))) then
                                    if(abs(ycoor-x4nds(2,i))<tol .or. &
                                    (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<ycoor.and. &
                                    (ycoor-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                                    (x4nds(2,i)>ycoor.and.x4nds(2,i)<yline(iy+1).and. &
                                    (x4nds(2,i)-ycoor)<(yline(iy+1)-x4nds(2,i)))) then
                                        n4yn(i) = 1
                                        n4out = n4out + 1
                                        an4nds(1,n4out) = i
                                        an4nds(2,n4out) = nnode
                                        exit     !if node found, jump out the loop
                                    endif
                                endif
                            endif
                        endif
                    enddo
                endif
                !...identify output nodes (off-fault)
                !Part3. Stations along ix==nx
                if(ix==nx .and. iy>1.and.iy<ny) then  !at surface only
                    do i=1,n4nds
                        if(n4yn(i)==0) then
                            if (abs(zcoor-x4nds(3,i))<tol) then
                                if(x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<xcoor.and. &
                                (xcoor-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) then
                                    if(abs(ycoor-x4nds(2,i))<tol .or. &
                                    (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<ycoor.and. &
                                    (ycoor-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
                                    (x4nds(2,i)>ycoor.and.x4nds(2,i)<yline(iy+1).and. &
                                    (x4nds(2,i)-ycoor)<(yline(iy+1)-x4nds(2,i)))) then
                                        n4yn(i) = 1
                                        n4out = n4out + 1
                                        an4nds(1,n4out) = i
                                        an4nds(2,n4out) = nnode
                                        exit     !if node found, jump out the loop
                                    endif
                                endif
                            endif
                        endif
                    enddo
                endif    

                do ift=1,ntotft 
                    ynft(ift) = .false.
                    if(xcoor>=(fltxyz(1,1,ift)-tol).and.xcoor<=(fltxyz(2,1,ift)+tol).and. &
                    ycoor>=(fltxyz(1,2,ift)-tol).and.ycoor<=(fltxyz(2,2,ift)+tol).and. &
                    zcoor>=(fltxyz(1,3,ift)-tol) .and. zcoor<=(fltxyz(2,3,ift)+tol)) then
                        if(ift==1) then
                            ynft(ift) = .true.
                        endif
                        if(ynft(ift)) then
                            nftnd0(ift) = nftnd0(ift) + 1
                            if(nftnd(ift)==0) then
                            write(*,*) 'inconsistency between mesh4num and meshgen for nftnd0',me
                            write(*,*) 'current node x,y,z:',xcoor,ycoor,zcoor
                            stop 2001
                            endif
                            nsmp(1,nftnd0(ift),ift) = nnode !slave node
                            msnode=nx*ny*nz+nftnd0(ift)                                 
                            locid(msnode)=ntag
                            dof1(msnode)=3         
                            nsmp(2,nftnd0(ift),ift) = msnode !master node
                            plane2(ny+ift,iz) = msnode
                            x(1,msnode) = xcoor
                            x(2,msnode) = ycoor
                            x(3,msnode) = zcoor
                            
                            if (rough_fault == 1) then 
                                x(2,msnode) = ycoort
                            endif
                            
                            !3DMPI
                            if(ix == 1) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 1
                                fltnum(1) = fltnum(1) + 1
                            endif
                            if(ix == nx) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 2
                                fltnum(2) = fltnum(2) + 1
                            endif
                            if(iy == 1) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 10
                                fltnum(3) = fltnum(3) + 1
                            endif
                            if(iy == ny) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 20
                                fltnum(4) = fltnum(4) + 1
                            endif
                            if(iz == 1) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 100
                                fltnum(5) = fltnum(5) + 1
                            endif
                            if(iz == nz) then
                                fltgm(nftnd0(ift)) = fltgm(nftnd0(ift)) + 200
                                fltnum(6) = fltnum(6) + 1
                            endif             

                            if(ift==1) then
                                tmp1 = xcoor
                            elseif(ift==2) then
                                tmp1 = sqrt(xcoor*xcoor+ycoor*ycoor)
                            endif
                            do i1=1,nonfs(ift)
                                if(abs(tmp1-xonfs(1,i1,ift))<tol .and. &
                                abs(zcoor-xonfs(2,i1,ift))<tol) then
                                    n4onf = n4onf + 1
                                    anonfs(1,n4onf) = nftnd0(ift)
                                    anonfs(2,n4onf) = i1
                                    anonfs(3,n4onf) = ift
                                    exit
                                elseif(abs(zcoor-xonfs(2,i1,ift))<tol .and. &
                                abs(tmp1-xonfs(1,i1,ift))<dx/2) then
                                    n4onf = n4onf + 1
                                    anonfs(1,n4onf) = nftnd0(ift)
                                    anonfs(2,n4onf) = i1
                                    anonfs(3,n4onf) = ift
                                    tmp2=xonfs(1,i1,ift)*dcos(brangle)
                                    tmp3=-xonfs(1,i1,ift)*dsin(brangle)
                                    x(1,nnode)=tmp2 !move mesh to stations!
                                    x(1,nnode-1)=tmp2
                                    x(2,nnode)=tmp3
                                    x(2,nnode-1)=tmp3
                                    exit          
                                endif
                            enddo  
                            !...establish equation numbers for this master node
                            do i1=1,ndof
                                neq0 = neq0 + 1
                                ntag=ntag+1
                                id1(ntag)=neq0!assume no explicit boundary nodes
                            enddo
                            !...assign unit vectors to the split node pair    
                            un(1,nftnd0(ift),ift) = dcos(fltxyz(1,4,ift))*dsin(fltxyz(2,4,ift))
                            un(2,nftnd0(ift),ift) = -dsin(fltxyz(1,4,ift))*dsin(fltxyz(2,4,ift))
                            un(3,nftnd0(ift),ift) = dcos(fltxyz(2,4,ift))        
                            us(1,nftnd0(ift),ift) = -dsin(fltxyz(1,4,ift))
                            us(2,nftnd0(ift),ift) = -dcos(fltxyz(1,4,ift))
                            us(3,nftnd0(ift),ift) = 0.0d0
                            ud(1,nftnd0(ift),ift) = dcos(fltxyz(1,4,ift))*dcos(fltxyz(2,4,ift))
                            ud(2,nftnd0(ift),ift) = dsin(fltxyz(1,4,ift))*dcos(fltxyz(2,4,ift))
                            ud(3,nftnd0(ift),ift) = dsin(fltxyz(2,4,ift))
                            
                            if (rough_fault == 1) then
                                un(1,nftnd0(ift),ift) = -pfx/(pfx**2 + 1.0d0 + pfz**2)**0.5
                                un(2,nftnd0(ift),ift) = 1.0d0/(pfx**2 + 1.0d0 + pfz**2)**0.5
                                un(3,nftnd0(ift),ift) = -pfz/(pfx**2 + 1.0d0 + pfz**2)**0.5    
                                us(1,nftnd0(ift),ift) = 1.0d0/(1.0d0 + pfx**2)**0.5
                                us(2,nftnd0(ift),ift) = pfx/(1.0d0 + pfx**2)**0.5
                                us(3,nftnd0(ift),ift) = 0.0d0
                                ud(1,nftnd0(ift),ift) = us(2,nftnd0(ift),ift)*un(3,nftnd0(ift),ift) &
                                    - us(3,nftnd0(ift),ift)*un(2,nftnd0(ift),ift)
                                ud(2,nftnd0(ift),ift) = us(3,nftnd0(ift),ift)*un(1,nftnd0(ift),ift) &
                                    - us(1,nftnd0(ift),ift)*un(3,nftnd0(ift),ift)
                                ud(3,nftnd0(ift),ift) = us(1,nftnd0(ift),ift)*un(2,nftnd0(ift),ift) &
                                    - us(2,nftnd0(ift),ift)*un(1,nftnd0(ift),ift)
                            endif                             
                            
                            !...prepare for area calculation
                            if(ixfi(ift)==0) ixfi(ift)=ix
                            if(izfi(ift)==0) izfi(ift)=iz
                            ifs(ift)=ix-ixfi(ift)+1
                            ifd(ift)=iz-izfi(ift)+1
                            fltrc(1,ifs(ift),ifd(ift),ift) = msnode    !master node
                            fltrc(2,ifs(ift),ifd(ift),ift) = nftnd0(ift) !fault node num in sequence

                            ! fric(1,nftnd0(ift),ift) = fric_sw_fs
                            ! fric(2,nftnd0(ift),ift) = fric_sw_fd        
                            ! fric(3,nftnd0(ift),ift) = fric_sw_D0 !D0
                            ! fric(4,nftnd0(ift),ift) = 0.0          !cohesion
                            ! fric(5,nftnd0(ift),ift) = 0.03      !Viscoplastic relaxation time
                            
                            !if (C_elastic == 0) then 
                            !    fric(6,nftnd0(ift),ift) = 0.0d0!rhow*grav*abs(zcoor)*gama        !pore pressure
                            !    !if (TPV==2802) fric(6,nftnd0(ift),ift) = rhow*grav*min(abs(zcoor),5.0d3)
                            !elseif (C_elastic == 1) then 
                            !    fric(6,nftnd0(ift),ift) = 0.0d0
                            !endif 
                            
                            !if (C_elastic==1) then
                                ! fric(7,nftnd0(ift),ift) = -120.0d6  
                                ! fric(8,nftnd0(ift),ift) = 40.0d6 
                                ! fric(50,nftnd0(ift),ift) = fric(8,nftnd0(ift),ift)
                                ! if (TPV==104) then 
                                    ! fric(7,nftnd0(ift),ift) = -120.0d6  
                                    ! fric(8,nftnd0(ift),ift) = 40.0d6 
                                    ! fric(50,nftnd0(ift),ift) = fric(8,nftnd0(ift),ift)
                                ! endif
                                !if (TPV==105) then 
                                !    fric(7,nftnd0(ift),ift) = -max(min(grav*1670.0d0*abs(zcoor),45.0d6),grav*1670.0d0*dx/2.0d0) 
                                !    fric(8,nftnd0(ift),ift) = -fric(7,nftnd0(ift),ift)*0.41d0
                                !    fric(50,nftnd0(ift),ift) = fric(8,nftnd0(ift),ift)
                                !endif
                            !endif

                            if (friclaw >= 3)then 
                                !if (TPV==105) then
                                !    call fb1(xcoor,fric_ww,fric_w, tmp1)
                                !    call fb2(-zcoor,fric_ww,fric_w, tmp2)
                                !endif
                                ! if (TPV == 104) then 
                                    ! if (abs(xcoor)<=15e3) then 
                                        ! tmp1=1.0
                                    ! elseif ((abs(xcoor)<18e3).and.(abs(xcoor)>15e3)) then 
                                        ! tmp1=0.5*(1+dtanh(3e3/(abs(xcoor)-18e3)+3e3/(abs(xcoor)-15e3)))
                                    ! else
                                        ! tmp1=0.0
                                    ! endif
                                    ! if (abs(zcoor--7.5e3)<=7.5e3) then 
                                        ! tmp2=1.0
                                    ! elseif ((abs(zcoor--7.5e3)<10.5e3).and.(abs(zcoor--7.5e3)>7.5e3)) then 
                                        ! tmp2=0.5*(1+dtanh(3e3/(abs(zcoor--7.5e3)-10.5e3)+3e3/(abs(zcoor--7.5e3)-7.5e3)))
                                    ! else
                                        ! tmp2=0.0
                                    ! endif                
                                ! endif
                                if (TPV == 2800 .or. TPV == 2801 .or. TPV == 2802) then 
                                    ! if (abs(xcoor)<=15.0d3) then 
                                        ! tmp1=1.0d0
                                    ! elseif ((abs(xcoor)<18.0d3).and.(abs(xcoor)>15.0d3)) then 
                                        ! tmp1=0.5d0*(1.0d0+dtanh(3.0d3/(abs(xcoor)-18.0d3)+3.0d3/(abs(xcoor)-15.0d3)))
                                    ! else
                                        ! tmp1=0.0d0
                                    ! endif
                                    ! if (abs(zcoor--7.5d3)<=7.5d3) then 
                                        ! tmp2=1.0d0
                                    ! elseif ((abs(zcoor--7.5d3)<10.5d3).and.(abs(zcoor--7.5d3)>7.5d3)) then 
                                        ! tmp2=0.5d0*(1+dtanh(3.0d3/(abs(zcoor--7.5d3)-10.5d3)+3.0d3/(abs(zcoor--7.5d3)-7.5d3)))
                                    ! else
                                        ! tmp2=0.0d0
                                    ! endif                                
                                    if (abs(xcoor)<=20.0d3) then 
                                        tmp1 = 1.0d0
                                    elseif ((abs(xcoor)>20.0d3).and.(abs(xcoor)<22.0d3)) then 
                                        tmp1 = 1.0d0 - (abs(xcoor)-20.0d3)/2.0d3
                                    elseif (abs(xcoor)>=22.0d3) then 
                                        tmp1 = 0.0d0
                                    endif
                                    if (abs(zcoor)>=2.0d3 .and. abs(zcoor)<=14.0d3) then 
                                        tmp2 = 1.0d0
                                    elseif (abs(zcoor)<2.0d3) then 
                                        tmp2 = 1.0d0 - (2.0d3 - abs(zcoor))/2.0d3
                                    elseif ((abs(zcoor)>14.0d3) .and. (abs(zcoor)<15.0d3)) then
                                        tmp2 = 1.0d0 - (abs(zcoor)-14.0d3)/1.0d3
                                    elseif (abs(zcoor)>=15.0d3) then 
                                        tmp2 = 0.0d0
                                    endif                
                                    tmp3 = fric_rsf_deltaa0
                                    if (abs(zcoor)>15.0d3) then 
                                        tmp3 = 2.0d0*fric_rsf_deltaa0
                                    endif 
                                endif                                
                                ! fric(9,nftnd0(ift),ift) = fric_rsf_a + tmp3*(1.0d0-tmp1*tmp2)
                                ! fric(10,nftnd0(ift),ift) = fric_rsf_b!b 
                                ! fric(11,nftnd0(ift),ift) = fric_rsf_Dc!RSF critical distance.
                                ! fric(12,nftnd0(ift),ift) = fric_rsf_v0!RSF:V0
                                ! fric(13,nftnd0(ift),ift) = fric_rsf_r0!RSF:miu0
                                ! fric(14,nftnd0(ift),ift) = fric_rsf_fw !RSF: fw for strong rate weakenging
                                ! !call fb1(xcoor,fric_ww,fric_w, tmp1)
                                ! !call fb2(-zcoor,fric_ww,fric_w, tmp2)
                                ! fric(15,nftnd0(ift),ift) = fric_rsf_vw +fric_rsf_deltavw0*(1.0d0-tmp1*tmp2) !RSF: Vw for strong rate weakening
                                ! call fb1(xcoor,fric_ww,fric_w, tmp1)
                                ! call fb3(-zcoor,fric_ww,fric_w, tmp2)
                                ! fric(20,nftnd0(ift),ift) = fric_tp_a_hy + fric_tp_deltaa_hy0*(1.0d0-tmp1*tmp2)
                                ! fric(16,nftnd0(ift),ift)=0.0d0 !RSF: initial normal slip rate 
                                ! fric(17,nftnd0(ift),ift)=fric_rsf_vinix!RSF:s 
                                ! fric(18,nftnd0(ift),ift)=fric_rsf_viniz!RSF:d
                                ! fric_ini_sliprate = sqrt(fric_rsf_vinix**2 + fric_rsf_viniz**2)
                                ! fric(19,nftnd0(ift),ift)=fric_ini_sliprate!RSF:mag    
                                ! fric(49,nftnd0(ift),ift)=fric_ini_sliprate!RSF:mag    
                                ! if(friclaw == 3) then 
                                    ! !Theta=d0/v0*dexp(a*dlog(2*dsinh()))
                                    ! fric(21,nftnd0(ift),ift) = fric(11,nftnd0(ift),ift)/fric(12,nftnd0(ift),ift)*dexp((fric(9,nftnd0(ift),ift)*dlog(2.0d0*dsinh &
                                                                ! (sqrt(fric(8,nftnd0(ift),ift)**2+0.0d0**2)/abs(fric(7,nftnd0(ift),ift))/fric(9,nftnd0(ift),ift))) &
                                                                ! -fric(13,nftnd0(ift),ift)- fric(9,nftnd0(ift),ift)*dlog(sqrt((fric(16,nftnd0(ift),ift))**2+(fric(17,nftnd0(ift),ift))**2+(fric(18,nftnd0(ift),ift))**2)/fric(12,nftnd0(ift),ift)))/fric(10,nftnd0(ift),ift))
                                ! elseif(friclaw == 4 .or. friclaw == 5) then
                                    ! !Theta0=a*log(2V0/Vini*sinh(TAOini/a/SigmaNini))
                                    ! fric(21,nftnd0(ift),ift)=fric(9,nftnd0(ift),ift)*dlog(2.0d0*fric(12,nftnd0(ift),ift)/sqrt((fric(16,nftnd0(ift),ift))**2+(fric(17,nftnd0(ift),ift))**2+(fric(18,nftnd0(ift),ift))**2) &
                                                            ! *dsinh(sqrt(fric(8,nftnd0(ift),ift)**2+0.0d0**2)/abs(fric(7,nftnd0(ift),ift))/fric(9,nftnd0(ift),ift)))
                                ! endif                                
                                ! if (xcoor==xsource.and.zcoor==zsource) then
                                    ! write(*,*) 'fric,9,a',fric(9,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,10,b',fric(10,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,11,Dc',fric(11,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,12,v0',fric(12,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,13,r0',fric(13,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,14,fw',fric(14,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,15,rsf_vw',fric(15,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,20,a_hy',fric(20,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,17,inisliprate',fric(17,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,7,norm',fric(7,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,8,tstk',fric(8,nftnd0(ift),ift)
                                    ! write(*,*) 'fric,50,tstk',fric(50,nftnd0(ift),ift)
                                ! endif
                            endif
    !-----------------------------END ZONE IV---------------------------!                        
                            !special values below.                        
                            if(abs(xcoor-fltxyz(1,1,ift))<tol .or. abs(xcoor-fltxyz(2,1,ift))<tol &
                            .or. abs(zcoor-fltxyz(1,3,ift))<tol) then !-x,+x,-z for 1, +x,-z for 2
                                fric(1,nftnd0(ift),ift) = 10000.    !fault edge, pinned
                            endif
                            !if(abs(zcoor)<tol) fric(6,nftnd0(ift),ift) = 0.5 * dz  !surface pore pressure
                            !if(ift==2)  then !branch, fault 2, only 12 km strike rupturable
                            !  tmp1=sqrt(xcoor**2+ycoor**2)
                            !  if(tmp1>12000) fric(1,nftnd0(ift),ift) = 10000.
                            !endif               
                            exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
                        endif  !if ynft(ift)
                    endif  !if flt range
                enddo  !do ift 
                !...create elements
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                    nelement = nelement + 1
                    if(nelement>numel) then
                        write(*,*) 'more elements in meshgen than in mesh4num'
                        write(*,*) 'x,y,z',xcoor,ycoor,zcoor
                        stop
                    endif
                    
                    ! creating hexahedral elements, which et=1.
                    et(nelement)    = 1 
                    ien(1,nelement) = plane1(iy-1,iz-1)
                    ien(2,nelement) = plane2(iy-1,iz-1)
                    ien(3,nelement) = plane2(iy,iz-1)
                    ien(4,nelement) = plane1(iy,iz-1)
                    ien(5,nelement) = plane1(iy-1,iz)
                    ien(6,nelement) = plane2(iy-1,iz)
                    ien(7,nelement) = plane2(iy,iz)
                    ien(8,nelement) = plane1(iy,iz)

                    ! using element center coords to determine PML elements, et = 2.
                    xc=0.0d0 
                    ids(nelement)=ntags
                    do i=1,nen
                        do j=1,3
                            xc(j)=xc(j)+x(j,ien(i,nelement))
                        enddo
                    enddo
                    xc=xc/8.0d0
                    
                    if (xc(1)>PMLb(1).or.xc(1)<PMLb(2) &
                        .or.xc(2)>PMLb(3).or.xc(2)<PMLb(4) &
                        .or.xc(3)<PMLb(5)) then
                        
                        et(nelement) = 2
                        ntags        = ntags+15+6
                    else
                        ntags        = ntags+12
                    endif
                    
                    ! velocityStructure will assign Vp, Vs and rho
                    !   to the elem id nelement given its location xc.

                    call velocityStructure(nelement, xc)
                    
                    if (C_degen == 1) then 
                        call wedge(xc(1), xc(2), xc(3), nelement, ntags, iy, iz, nftnd0(1))
                    elseif (C_degen == 2) then 
                        call tetra(xc(1), xc(2), xc(3), nelement, ntags, iy, iz, nftnd0(1))
                    endif         
                    
                    if(et(nelement) == 1 .and. (xcoor>(fltxyz(1,1,1)-tol).and.xcoor<(fltxyz(2,1,1)+dx+tol).and. &
                        zcoor>(fltxyz(1,3,1)-tol).and.ycoor>0.and.abs(ycoor-dy)<tol)) then
                        do i=1,nftnd0(1)
                            do k=1,nen
                                if(ien(k,nelement)==nsmp(1,i,1)) then
                                    ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
                                endif
                            enddo
                        enddo
                    endif                

                    if (C_elastic==0 .and. TPV==2800) then                
                        tmp1 = -0.5d0*(zline(iz)+zline(iz-1))  !z<0, thus tmp1>0
                        tmp2 = tmp1 * grav
                        if(et(nelement)==1)then
                            eleporep(nelement)= rhow*tmp2  !pore pressure>0
                            s1(ids(nelement)+3)=-mat(nelement,3)*tmp2  !vertical, comp<0
                            s1(ids(nelement)+1)=-mat(nelement,3)*tmp2 
                            s1(ids(nelement)+2)=-mat(nelement,3)*tmp2 
                            s1(ids(nelement)+6)= (mat(nelement,3)-rhow)*tmp2 /3.0d0
                        elseif(et(nelement)==2)then
                            eleporep(nelement)= rhow*tmp2  !pore pressure>0
                            s1(ids(nelement)+3+15)=-mat(nelement,3)*tmp2  !vertical, comp<0
                            s1(ids(nelement)+1+15)=-mat(nelement,3)*tmp2
                            s1(ids(nelement)+2+15)=-mat(nelement,3)*tmp2 
                            s1(ids(nelement)+6+15)=(mat(nelement,3)-rhow)*tmp2 /3.0d0                        
                        endif                        
                    endif
                    if (C_elastic==0 .and. TPV==2801) then                
                        tmp1 = -0.5d0*(zline(iz)+zline(iz-1))  !z<0, thus tmp1>0
                        tmp2 = tmp1 * grav
                        if(et(nelement)==1)then
                            eleporep(nelement)= rhow*tmp2  !pore pressure>0
                            s1(ids(nelement)+3)=-mat(nelement,3)*tmp2  !vertical, comp<0
                            if (-s1(ids(nelement)+3)>100.0d6) s1(ids(nelement)+3) = -100.0d6
                            s1(ids(nelement)+1)=s1(ids(nelement)+3)
                            s1(ids(nelement)+2)=s1(ids(nelement)+3) 
                            s1(ids(nelement)+6)=(mat(nelement,3) -rhow)*0.4d0*tmp2
                        elseif(et(nelement)==2)then
                            eleporep(nelement)= rhow*tmp2  !pore pressure>0
                            s1(ids(nelement)+3+15)=-mat(nelement,3)*tmp2  !vertical, comp<0
                            if (-s1(ids(nelement)+3+15)>100.0d6) s1(ids(nelement)+3+15) = -100.0d6
                            s1(ids(nelement)+1+15)=s1(ids(nelement)+3+15)
                            s1(ids(nelement)+2+15)=s1(ids(nelement)+3+15) 
                            s1(ids(nelement)+6+15)=(mat(nelement,3) -rhow)*0.4d0*tmp2                        
                        endif                        
                    endif    
                    if (C_elastic==0 .and. TPV==2802) then 
                        tmp1 = -0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0  !z<0, thus tmp1>0
                        call plastic_set_mat_stress(tmp1, nelement)
                    endif 
    !--------------------------END ZONE III-----------------------------!                
                endif!if element
            enddo!iy
        enddo!iz
        plane1 = plane2
    enddo!ix
    maxs=ntags
    !-------------------------------------------------------------------!
    !-----------------Checking between mesh4 and meshgen----------------!
    !------------Updated by Oct.5.2015/ D. Liu--------------------------!
    if (maxs>=(5*maxm)) then
        write(*,*) '5*maxm',maxm,'is not enough for maxs',maxs
        stop 2002
    endif
    if(nnode/=nx*ny*nz.or.msnode/=numnp.or.nelement/=numel.or.neq0/=neq) then
        write(*,*) 'Inconsistency in node/element/equation/between meshgen and mesh4num: stop!',me
        write(*,*) 'nnode&numnp=',nnode,numnp
        write(*,*) 'nelement,numel=',nelement,numel
        write(*,*) 'neq0,neq=',neq0,neq
        write(*,*) 'nnode,nx,ny,nz',nnode,nx,ny,nz
        write(*,*) 'msnode,numnp',msnode,numnp
        stop 2003
    endif
    if(ntag/=maxm) then
        write(*,*) 'Inconsistency in ntag and maxm: stop!',me
        write(*,*) ntag,maxm
        stop 2004
    endif
    do i=1,ntotft
        if(nftnd0(i)/=nftnd(i)) then
            write(*,*) 'Inconsistency in fault between meshgen and mesh4num: stop!',me,i
            stop 2005
        endif
    enddo
    !-------------------------------------------------------------------!
    !...calculate area associated with each fault node pair and distance from source. 
    !B.D. 4/19/09
    !but only needed for nftnd>0. B.D. 10/16/09
    !for multiple faults. B.D. 1/7/12
    do ift=1,ntotft
        if(nftnd0(ift)>0) then
        !...distance from the source point
        do i=1,ifd(ift)
            do j=1,ifs(ift)
                i1 = fltrc(1,j,i,ift)
                j1 = fltrc(2,j,i,ift)
                r4nuc(j1,ift) = sqrt((x(1,i1)-xsource)**2 + (x(2,i1)-ysource)**2 &
                + (x(3,i1)-zsource)**2)
            enddo
        enddo
        !...element areas and distribute evenly to its four nodes
        do i=2,ifd(ift)
            do j=2,ifs(ift)
            !...4 nodes of quadrilateral
                n1 = fltrc(1,j,i,ift) !use nodal number in nnode
                n2 = fltrc(1,j-1,i,ift)
                n3 = fltrc(1,j-1,i-1,ift)
                n4 = fltrc(1,j,i-1,ift)
                m1 = fltrc(2,j,i,ift) !use nodal number in nftnd0
                m2 = fltrc(2,j-1,i,ift)
                m3 = fltrc(2,j-1,i-1,ift)
                m4 = fltrc(2,j,i-1,ift)
                !...calculate area of the quadrilateral
                !......if fault is not in coor axes plane
                aa1=sqrt((x(1,n2)-x(1,n1))**2 + (x(2,n2)-x(2,n1))**2 &
                + (x(3,n2)-x(3,n1))**2)
                bb1=sqrt((x(1,n3)-x(1,n2))**2 + (x(2,n3)-x(2,n2))**2 &
                + (x(3,n3)-x(3,n2))**2)
                cc1=sqrt((x(1,n4)-x(1,n3))**2 + (x(2,n4)-x(2,n3))**2 &
                + (x(3,n4)-x(3,n3))**2)
                dd1=sqrt((x(1,n1)-x(1,n4))**2 + (x(2,n1)-x(2,n4))**2 &
                + (x(3,n1)-x(3,n4))**2)
                p1=sqrt((x(1,n4)-x(1,n2))**2 + (x(2,n4)-x(2,n2))**2 &
                + (x(3,n4)-x(3,n2))**2)
                q1=sqrt((x(1,n3)-x(1,n1))**2 + (x(2,n3)-x(2,n1))**2 &
                + (x(3,n3)-x(3,n1))**2)
                area=0.25 * sqrt(4*p1*p1*q1*q1 - &
               (bb1*bb1 + dd1*dd1 - aa1*aa1 -cc1*cc1)**2) 
                !...distribute above area to 4 nodes evenly
                area = 0.25 * area
                arn(m1,ift) = arn(m1,ift) + area
                arn(m2,ift) = arn(m2,ift) + area
                arn(m3,ift) = arn(m3,ift) + area
                arn(m4,ift) = arn(m4,ift) + area
                enddo
            enddo
            !...save area for moment (rate) calculation before MPI. B.D. 8/11/10
            arn4m = arn
        endif!end if nftnd0=0/ not
        !
        time1 = MPI_WTIME()
        !3DMPI: Prepare for MPI partitioning
        fltMPI=.false.

        if(fltnum(1) /= 0) allocate(fltl(fltnum(1)))
        if(fltnum(2) /= 0) allocate(fltr(fltnum(2)))
        if(fltnum(3) /= 0) allocate(fltf(fltnum(3)))
        if(fltnum(4) /= 0) allocate(fltb(fltnum(4)))
        if(fltnum(5) /= 0) allocate(fltd(fltnum(5)))
        if(fltnum(6) /= 0) allocate(fltu(fltnum(6)))
        fltnum = 0
        do i = 1, nftnd0(ift)
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
        if (npx > 1) then
            bndl=1  !left boundary
            bndr=nx ! right boundary
            if (mex == master) then
                bndl=0
            elseif (mex == npx-1) then
                bndr=0
            endif
     
            if (bndl/=0) then
                call mpi_sendrecv(fltnum(1), 1, MPI_INTEGER, me-npy*npz, 1000+me, &
                    itmp1, 1, MPI_INTEGER, me-npy*npz, 1000+me-npy*npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(1)>0 .and. itmp1>0 ) then  
                    if(fltnum(1) /= itmp1) stop 'error in fltnum(1)'
                    fltMPI(1)=.true.
                    allocate(btmp(fltnum(1)),btmp1(fltnum(1)))
                    btmp = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(1)
                        btmp(i)=arn(fltl(i),ift)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(1), MPI_DOUBLE_PRECISION, me-npy*npz, 1000+me, &
                        btmp1, fltnum(1), MPI_DOUBLE_PRECISION, me-npy*npz, 1000+me-npy*npz, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(1)
                        arn(fltl(i),ift) = arn(fltl(i),ift) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
            endif !if bhdl/=0

            if (bndr/=0) then
                call mpi_sendrecv(fltnum(2), 1, MPI_INTEGER, me+npy*npz, 1000+me, &
                    itmp1, 1, MPI_INTEGER, me+npy*npz, 1000+me+npy*npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(2)>0 .and. itmp1>0 ) then  
                    if(fltnum(2) /= itmp1) stop 'error in fltnum(2)'
                    fltMPI(2)=.true.
                    allocate(btmp(fltnum(2)),btmp1(fltnum(2)))
                    btmp  = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(2)
                        btmp(i)=arn(fltr(i),ift)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(2), MPI_DOUBLE_PRECISION, me+npy*npz, 1000+me, &
                        btmp1, fltnum(2), MPI_DOUBLE_PRECISION, me+npy*npz, 1000+me+npy*npz, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(2)
                        arn(fltr(i),ift) = arn(fltr(i),ift) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
            endif !bndr/=0
        endif !npx>1
    !*****************************************************************************************
        if (npy > 1) then   !MPI:no fault message passing along y due to the x-z vertical fault.
            bndf=1  !front boundary
            bndb=ny ! back boundary
            if (mey == master) then
                bndf=0
            elseif (mey == npy-1) then
                bndb=0
            endif

            if (bndf/=0) then
                call mpi_sendrecv(fltnum(3), 1, MPI_INTEGER, me-npz, 2000+me, &
                    itmp1, 1, MPI_INTEGER, me-npz, 2000+me-npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(3)>0 .and. itmp1>0) then
                    if(fltnum(3) /= itmp1) stop 'error in fltnum(3)'
                    fltMPI(3)=.true.
                    allocate(btmp(fltnum(3)),btmp1(fltnum(3)))
                    btmp = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(3)
                        btmp(i)=arn(fltf(i),ift)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(3), MPI_DOUBLE_PRECISION, me-npz, 2000+me, &
                        btmp1, fltnum(3), MPI_DOUBLE_PRECISION, me-npz, 2000+me-npz, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(3)
                        arn(fltf(i),ift) = arn(fltf(i),ift) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
            endif !bhdf/=0

            if (bndb/=0) then
                call mpi_sendrecv(fltnum(4), 1, MPI_INTEGER, me+npz, 2000+me, &
                    itmp1, 1, MPI_INTEGER, me+npz, 2000+me+npz, &
                    MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(4)>0 .and. itmp1>0) then  
                    if(fltnum(4) /= itmp1) stop 'error in fltnum(4)'
                    fltMPI(4)=.true.
                    allocate(btmp(fltnum(4)),btmp1(fltnum(4)))
                    btmp  = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(4)
                        btmp(i)=arn(fltb(i),ift)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(4), MPI_DOUBLE_PRECISION, me+npz, 2000+me, &
                        btmp1, fltnum(4), MPI_DOUBLE_PRECISION, me+npz, 2000+me+npz, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(4)
                        arn(fltb(i),ift) = arn(fltb(i),ift) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
            endif !bndb/=0
         endif !npy>1
    !*****************************************************************************************
        if (npz > 1) then
            bndd=1  !lower(down) boundary
            bndu=nz ! upper boundary
            if (mez == master) then
                bndd=0
            elseif (mez == npz-1) then
                bndu=0
            endif
            if (bndd/=0) then
                call mpi_sendrecv(fltnum(5), 1, MPI_INTEGER, me-1, 3000+me, &
                    itmp1, 1, MPI_INTEGER, me-1, 3000+me-1, &
                    MPI_COMM_WORLD, istatus, ierr)
                !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                if(fltnum(5)>0 .and. itmp1>0) then
                    if(fltnum(5) /= itmp1) stop 'error in fltnum(5)'
                    fltMPI(5)=.true.
                    allocate(btmp(fltnum(5)),btmp1(fltnum(5)))
                    btmp = 0.
                    btmp1 = 0.
                    do i = 1, fltnum(5)
                        btmp(i)=arn(fltd(i),ift)
                    enddo
                    call mpi_sendrecv(btmp, fltnum(5), MPI_DOUBLE_PRECISION, me-1, 3000+me, &
                        btmp1, fltnum(5), MPI_DOUBLE_PRECISION, me-1, 3000+me-1, &
                        MPI_COMM_WORLD, istatus, ierr)
                    do i = 1, fltnum(5)
                        arn(fltd(i),ift) = arn(fltd(i),ift) + btmp1(i)
                    enddo
                    deallocate(btmp,btmp1)       
                endif
            endif !bhdd/=0
            
            if (bndu/=0) then
                call mpi_sendrecv(fltnum(6), 1, MPI_INTEGER, me+1, 3000+me, &
                    itmp1, 1, MPI_INTEGER, me+1, 3000+me+1, &
                    MPI_COMM_WORLD, istatus, ierr)
                    !write(*,*) 'me= ',me, 'itmp1=',itmp1                
                    if(fltnum(6)>0 .and. itmp1>0) then  
                        if(fltnum(6) /= itmp1) stop 'error in fltnum(6)'
                        fltMPI(6)=.true.
                        allocate(btmp(fltnum(6)),btmp1(fltnum(6)))
                        btmp  = 0.
                        btmp1 = 0.
                        do i = 1, fltnum(6)
                            btmp(i)=arn(fltu(i),ift)
                        enddo
                        call mpi_sendrecv(btmp, fltnum(6), MPI_DOUBLE_PRECISION, me+1, 3000+me, &
                            btmp1, fltnum(6), MPI_DOUBLE_PRECISION, me+1, 3000+me+1, &
                            MPI_COMM_WORLD, istatus, ierr)
                        do i = 1, fltnum(6)
                        arn(fltu(i),ift) = arn(fltu(i),ift) + btmp1(i)
                        enddo
                        deallocate(btmp,btmp1)       
                    endif
            endif !bndu/=0
        endif !npz>1
    enddo !ift=i:nft 
    !
    time2 = MPI_WTIME()
    btime = btime + (time2 - time1) 

    deallocate(xlinet,xline,yline,zline,fltrc)
    
end subroutine meshgen

subroutine velocityStructure(nelement, xc)
! Subroutine velocityStructure will asign Vp, Vs and rho 
!   based on input from bMaterial.txt, which is created by 
!   case input file user_defined_param.py.

    use globalvar
    implicit none
    integer (kind = 4) :: nelement, i
    real (kind = dp) :: xc(3)
    
    if (nmat==1 .and. n2mat == 3) then
    ! homogenous material
        mat(nelement,1)  = material(1,1)
        mat(nelement,2)  = material(1,2)
        mat(nelement,3)  = material(1,3)
    elseif (nmat>1 .and. n2mat == 4) then 
        ! 1D velocity structure
        ! material = material(i,j), i=1,nmat, and j=1,4
        ! for j
        !   1: bottom depth of a layer (should be positive in m).
        !   2: Vp, m/s
        !   3: Vs, m/s
        !   4: rho, kg/m3
        if (abs(xc(3)) < material(1,1)) then
            mat(nelement,1)  = material(1,2)
            mat(nelement,2)  = material(1,3)
            mat(nelement,3)  = material(1,4)
        else
            do i = 2, nmat
                if (abs(xc(3)) < material(i,1) &
                    .and. abs(xc(3)) >= material(i-1,1)) then
                    
                    mat(nelement,1)  = material(i,2)
                    mat(nelement,2)  = material(i,3)
                    mat(nelement,3)  = material(i,4)
                endif 
            enddo
        endif
    endif 
    
    ! calculate lambda and mu from Vp, Vs and rho.
    ! mu = Vs**2*rho
    mat(nelement,5)  = mat(nelement,2)**2*mat(nelement,3)
    ! lambda = Vp**2*rho-2*mu
    mat(nelement,4)  = mat(nelement,1)**2*mat(nelement,3)-2.0d0*mat(nelement,5)

end subroutine velocityStructure