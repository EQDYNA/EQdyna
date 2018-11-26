subroutine mesh4num(fltxyz,dx,xmin,xmax,ymin,ymax,zmin,zmax,&
nnode,nelement,neq,PMLb,maxm,nftnd,master,me,nprocs)!Adding PMLb and maxm
!... program to generate mesh for 3D models.
! This model is a strike-slip fault in the half space, for SCEC TPV16 & 17.
! B.D. 1/1/12
use globalvar
implicit none
include 'mpif.h'
!...node # and related
integer (kind=4)::nnode,nelement,nxt,nx,ny,nz,ix,iy,iz, &
	edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,np=10000,ncnt,nfrt,nbck,edgezn
!...fault definition: 
integer (kind=4),dimension(ntotft) :: nftnd
real (kind=8),dimension(2,4,ntotft) :: fltxyz
real (kind=8) :: x1,x2,x3=0,y1,y2,y3=0,z1,z2,z3=0 !3 points for fault plane (inc origin)
!...entire model size
real (kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax
!...grid size
real (kind=8) :: dx,dy,dz	!dx is given in parcon.f90, dy/dz will be from dx
real (kind=8) :: rat=1.025 !enlarge ratio for buffers to use
!...uniform element num from fault y-ccor
integer (kind=4) :: dis4uniF=20,dis4uniB=20
!...mesh results
integer (kind=4)::neq
!...working variables
logical,dimension(ntotft) :: ynft
integer (kind=4) :: nxuni,nyuni,nzuni,maxtotal,ift,n1,n2,n3,n4,m1,m2,m3,m4,alloc_err
integer (kind=4) :: me, master, nprocs, rlp, rr, ierr,jj,itmp1,itmp2
integer istatus(MPI_STATUS_SIZE)
real (kind=8) :: tol,xcoor,ycoor,zcoor,xstep,ystep,zstep, &
		tmp1,tmp2,tmp3,a,b,fltybtm,linea,lineb,linea1,lineb1,ycoor1,&
		xmin1,xmax1,ymin1,ymax1,zmin1
real (kind=8),allocatable,dimension(:) :: xlinet,xline,yline,zline
!DL
integer(kind=4)::maxm,ntag,zz,nl=12
real(kind=8)::PMLb(8),mdx(3)
!
!dy = dx * abs(dtan(brangle))
dy=dx
dz = dx
tol = dx/100

!...determine num of nodes along x
nxuni = (fltxyz(2,1,1) - fltxyz(1,1,1)) / dx + 1
xstep = dx
xcoor = fltxyz(1,1,1)
do ix = 1, np
	xstep = xstep * rat
	xcoor = xcoor - xstep
	if(xcoor < xmin) exit
enddo
edgex1 = ix + 1
xstep = dx
xcoor = fltxyz(2,1,1)
do ix = 1, np
	xstep = xstep * rat
	xcoor = xcoor + xstep
	if(xcoor > xmax) exit
enddo
nxt = nxuni + edgex1 + ix + 1
allocate(xlinet(nxt))
!
!...predetermine x-coor
xlinet(edgex1+1) = fltxyz(1,1,1)
xstep = dx
do ix = edgex1, 1, -1
	xlinet(ix) = xlinet(ix+1) - xstep
	xstep = xstep * rat
enddo
xmin1=xlinet(1)
do ix = edgex1+2,edgex1+nxuni+1
	xlinet(ix) = xlinet(ix-1) + dx
enddo
xstep = dx
do ix = edgex1+nxuni+2,nxt
	xlinet(ix) = xlinet(ix-1) + xstep
	xstep = xstep * rat
enddo
xmax1=xlinet(nxt)
!  
!
!***********MPI***************  
!...MPI partitioning based on ix and xlinet. B.D. 4/18/09
!...need overlap between adjacent procs. B.D. 5/10/09
!...more evenly distributed. B.D. 11/1/09
!...due to overlap, total line num should be: nxt+nprocs-1
!	B.D. 1/19/10
j1 = nxt + nprocs - 1
rlp = j1/nprocs
rr = j1 - rlp*nprocs
if(me<(nprocs-rr)) then
	nx = rlp
else
	nx = rlp + 1	!evenly distributed to last rr
endif
allocate(xline(nx))
if(me<=nprocs-rr) then
	do ix=1,nx
		xline(ix) = xlinet(rlp*me-me+ix)
	enddo
else
	do ix=1,nx
		k1 = me - (nprocs - rr)
		xline(ix) = xlinet(rlp*me-me+k1+ix)
	enddo
endif
!  if(me==nprocs-1) then
!    nx=rlp + rr
!  else
!    nx=rlp + 1
!  endif
!  allocate(xline(nx))
!  do ix=1,nx
!    xline(ix) = xlinet(me*rlp+ix)
!  enddo
!write(*,*) 'me=',me,'xline(1),xline(nx)',xline(1),xline(nx),'in meshgen'
!
!...determine num of nodes along y
!nyuni = dis4uniF + fltxyz(2,1,1)/dx + dis4uniB + 1
nyuni = dis4uniF + dis4uniB + 1
ystep = dy
!ycoor = -dy*(dis4uniF+fltxyz(2,1,1)/dx+1)
ycoor = -dy*(dis4uniF+1)
do iy = 1, np
	ystep = ystep * rat
	ycoor = ycoor - ystep
	if(ycoor < ymin) exit
enddo
edgey1 = iy + 1
ystep = dy
ycoor = dy*(dis4uniB+1)
do iy = 1, np
	ystep = ystep * rat
	ycoor = ycoor + ystep
	if(ycoor > ymax) exit
enddo
ny = nyuni + edgey1 + iy + 1
!...pre-determine y-coor
allocate(yline(ny))
!...predetermine y-coor
!yline(edgey1+1) = -dy*(dis4uniF+fltxyz(2,1,1)/dx+1)
yline(edgey1+1) = -dy*(dis4uniF+1)
ystep = dy
do iy = edgey1, 1, -1
	yline(iy) = yline(iy+1) - ystep
	ystep = ystep * rat
enddo
ymin1=yline(1)
do iy = edgey1+2,edgey1+nyuni+1
	yline(iy) = yline(iy-1) + dy
enddo
ystep = dy
do iy = edgey1+nyuni+2,ny
	yline(iy) = yline(iy-1) + ystep
	ystep = ystep * rat
enddo
ymax1=yline(ny)
!
!...determine num of nodes along z
zstep = dz
zcoor = fltxyz(1,3,1)
do iz=1,np
	zstep = zstep * rat
	zcoor = zcoor - zstep
	if(zcoor < zmin) exit
enddo
edgezn = iz + 1
nzuni = (fltxyz(2,3,1)-fltxyz(1,3,1))/dx + 1 
nz = edgezn + nzuni
!...predetermine z-coor
allocate(zline(nz))
zline(nz) = zmax
do iz = nz-1,nz-nzuni,-1
	zline(iz) = zline(iz+1) - dz
enddo
zstep = dz
do iz = nz-nzuni-1,1,-1
	zline(iz) = zline(iz+1) -zstep
	zstep = zstep * rat
enddo
zmin1=zline(1)
!DL  
PMLb(1)=xlinet(nxt-nl)!PMLxmax
PMLb(2)=xlinet(nl+1)!PMLxmin
PMLb(3)=yline(ny-nl)!PMLymax
PMLb(4)=yline(nl+1)!PMLymin
PMLb(5)=zline(nl+1)!PMLzmin
mdx(1)=xlinet(nxt)-xlinet(nxt-1)
mdx(2)=yline(ny)-yline(ny-1)
mdx(3)=zline(2)-zline(1)
PMLb(6)=mdx(1)
PMLb(7)=mdx(2)
PMLb(8)=mdx(3)
write(*,*) 'PML info',me
write(*,*) 'xmax',PMLb(1),'xmin',PMLb(2),'ymax',PMLb(3),'ymin',PMLb(4),'zmin',PMLb(5)
write(*,*) 'Max dx',mdx(1),'Max dy',mdx(2),'Max dz',mdx(3)
write(*,*) 'Max',PMLb(6),PMLb(7),PMLb(8)
	if (PMLb(5).gt.fltxyz(1,3,1)) then
		stop 111
	elseif (PMLb(1).lt.fltxyz(2,1,1)) then
		stop 112
	endif	
!...prepare for digitizing
nnode = 0
nelement = 0
neq = 0
nftnd = 0
!DL 
ntag=0
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
			zz=ndof
			!DL. Change dof of this node to 12.
			if (xcoor>PMLb(1).or.xcoor<PMLb(2).or.ycoor>PMLb(3).or.ycoor<PMLb(4) &
				.or.zcoor<PMLb(5)) then
				zz = 12
			endif	
			!if ((abs(zcoor-PMLb(5))<tol).and.(xcoor<PMLb(1)+tol).and.(xcoor>PMLb(2)-tol).and.(ycoor<PMLb(3)+tol).and.(ycoor>PMLb(4)-tol)) then! z face
			!	zz=15
			!elseif ((abs(xcoor-PMLb(1))<tol).and.(ycoor<PMLb(3)+tol).and.(ycoor>PMLb(4)-tol).and.(zcoor>PMLb(5)-tol)) then!12 wall
			!	zz=15
			!elseif ((abs(xcoor-PMLb(2))<tol).and.(ycoor<PMLb(3)+tol).and.(ycoor>PMLb(4)-tol).and.(zcoor>PMLb(5)-tol)) then!34 wall
			!	zz=15	
			!elseif ((abs(ycoor-PMLb(3))<tol).and.(xcoor<PMLb(1)+tol).and.(xcoor>PMLb(2)-tol).and.(zcoor>PMLb(5)-tol)) then!41 wall
			!	zz=15		
			!elseif ((abs(ycoor-PMLb(4))<tol).and.(xcoor<PMLb(1)+tol).and.(xcoor>PMLb(2)-tol).and.(zcoor>PMLb(5)-tol)) then!23 wall
			!	zz=15					
			!endif
			!...establish equation numbers for this node
			do i1=1,zz
				!...elastoplastic off-fault, stress assigned in entire model,
				! need to fix model boundaries (except free surface).
				if(abs(xcoor-xmin1)<tol.or.abs(xcoor-xmax1)<tol.or.abs(ycoor-ymin1)<tol &
					.or.abs(ycoor-ymax1)<tol.or.abs(zcoor-zmin1)<tol) then
					i = -1  !-1 for fixed boundary nodes
					ntag=ntag+1
				else
					neq = neq + 1
					ntag=ntag+1					
				endif
			enddo
			!...fault plane: split nodes
			!...for multiple faults. B.D. 1/6/12
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
			!...create elements
			if(ix>=2 .and. iy>=2 .and. iz>=2) then
			nelement = nelement + 1
			!...when the current node is one element above the branch fault in y-coor,
			! one hexahedron degenerates into two wedges. B.D. 1/7/12
			!          if(xcoor>=(fltxyz(1,1,2)-tol).and.xcoor<=(fltxyz(2,1,2)+tol).and. &
			!            ycoor>(fltxyz(1,2,2)-tol).and.ycoor<=(fltxyz(2,2,2)+dy+tol).and. &
			!            zcoor>=(fltxyz(1,3,2)-tol).and.zcoor<=(fltxyz(2,3,2)+tol)) then
			!            !above is x,y,z ranges for possible degeneation. B.D. 1/7/12
			!            if(abs(ycoor+xcoor*dtan(brangle)-dy)<tol) then !degenerate
			!              nelement = nelement + 1 !one more wedge element
			!            endif
			!          endif
			endif  !if element
		enddo	!iy
	enddo	!iz
enddo		!ix
maxm=ntag
end subroutine mesh4num