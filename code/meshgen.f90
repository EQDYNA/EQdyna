subroutine meshgen(fltxyz,xmin,xmax,ymin,ymax,zmin,zmax,&
nnode1,nelement1,ien,mat,s1,eleporep,rho,vs,xnode,neq1,& !Delete id !change elestress as s1
id1,maxm,locid,dof1,et,PMLb,maxs,ids, & !Adding id1,maxm,loci,PMLb,maxs,ids
nr4mpi,er4mpi,n4out,an4nds,nftnd1,n4onf,xonfs,nonfs,&
nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,anonfs,&
arn4m,master,me,nprocs)
!... program to generate mesh for 3D models.
! This model is a strike-slip fault in the half space, for SCEC TPV16 & 17.
! B.D. 1/1/12
use globalvar
implicit none
include 'mpif.h'
!...node # and related
integer (kind=4)::nnode,nelement,nnode1,nelement1,nxt,nx,ny,nz,ix,iy,iz, &
edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,np=10000,ncnt,nfrt,nbck,edgezn
!...fault definition: 
integer (kind=4) :: nftmx,nonmx
integer (kind=4),dimension(ntotft) :: nftnd,nftnd1,ixfi,izfi
real (kind=8),dimension(2,4,ntotft) :: fltxyz
real (kind=8) :: x1,x2,x3=0,y1,y2,y3=0,z1,z2,z3=0 !3 points for fault plane (inc origin)
integer (kind=4),dimension(3,nonmx) :: anonfs
integer (kind=4),dimension(2,nftmx,ntotft) :: nsmp
real (kind=8),dimension(nftmx,ntotft) :: fnft,arn,r4nuc,arn4m,slp4fri
real (kind=8),dimension(3,nftmx,ntotft) :: un,us,ud
real (kind=8),dimension(8,nftmx,ntotft) :: fric
!...entire model size
real (kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax
!...grid size
real (kind=8) :: dy,dz	!dx is given in parcon.f90, dy/dz will be from dx
!...uniform element num from fault y-ccor
!For Double-couple point source. Ma and Liu (2006) 
!integer (kind=4) :: dis4uniF=50,dis4uniB=50
!For TPV8 
integer (kind=4) :: dis4uniF=100,dis4uniB=100
!...output on-fault stations. B.D. 10/25/09
integer (kind=4) :: n4onf
integer,dimension(ntotft)::nonfs
real (kind=8),dimension(2,8,ntotft) :: xonfs
!...output node coors. 6 Off-fault stations in this example. B.D. 10/25/11
integer (kind=4)::n4nds=6,n4out
integer (kind=4),dimension(6)::n4yn=0
integer (kind=4),dimension(2,6)::an4nds
real (kind=8),dimension(2,6)::x4nds=reshape((/0.0,-3.0,&!8.0,6.0
												0.0,-2.0,&
												0.0,-1.0,&
												0.0,1.0,&
												12.0,-3.0,&
												12.0,3.0/),(/2,6/))
!...mesh results
integer (kind=4)::neq,neq1
integer (kind=4),dimension(2,2)::nr4mpi,er4mpi  !node, equation range for MPI
integer (kind=4),dimension(8,nelement1)::ien 
!integer (kind=4),dimension(ndof,nnode1)::id
integer (kind=4),dimension(nelement1)::mat
real (kind=8),dimension(nelement1)::eleporep
!real (kind=8),dimension(nstr,nelement1)::elestress
real (kind=8),dimension(ndof,nnode1)::xnode
!...for initial stress and pore pressure of elements: inelastic. B.D. 1/8/12
real (kind=8)::grav=9.8,rhow=1000.,b22=0.44327040,b33=0.50911055, &
	b23=-0.15487679
real (kind=8),dimension(numat)::rho,vs!...working variables
logical,dimension(ntotft) :: ynft
integer (kind=4) :: nxuni,nyuni,nzuni,maxtotal,ift,n1,n2,n3,n4,m1,m2,m3,m4,alloc_err
integer (kind=4) :: me, master, nprocs, rlp, rr, ierr,jj,itmp1,itmp2
integer istatus(MPI_STATUS_SIZE)
integer (kind=4),dimension(ntotft) :: ifs,ifd
real (kind=8), allocatable, dimension(:) :: btmp, btmp1,btmp2,btmp3
real (kind=8) :: tol,xcoor,ycoor,zcoor,xstep,ystep,zstep, &
	tmp1,tmp2,tmp3,a,b,fltybtm,linea,lineb,linea1,lineb1,ycoor1
real (kind=8) :: area,a1,b1,c1,d1,p1,q1
integer (kind=4),allocatable,dimension(:,:) :: plane1,plane2
integer (kind=4),allocatable,dimension(:,:,:,:) :: fltrc
real (kind=8),allocatable,dimension(:) :: xlinet,xline,yline,zline
character (len=6) :: mm
character (len=30) :: meshfile

!*.* Variables for PML. D.L. Jan/23/2015
integer (kind=4)::maxm,maxs,label,ntag,ntags,zz
integer (kind=4),dimension(maxm)::id1
integer (kind=4),dimension(nnode1)::locid,dof1
integer (kind=4),dimension(nelement1)::et
real(kind=8),dimension(8)::PMLb
integer(kind=4),dimension(nelement1)::ids
real(kind=8),dimension(4*maxm)::s1
real(kind=8)::xc(3)
!*.* D.L. 
!...because eos and vanerne with linux do not allow 'write(mm,*) me',
! I have to use a way for conversion from num to char as below for
! what I want. Currently, only for MPI processes <= 300. For more, need
! to extend further.  B.D. 8/13/10
!...actually, use formatted writing, solve the problem
  write(mm,'(i6)') me
!  if(me<10) then
!    mm = char(me+48)
!  elseif(me<20) then
!    mm = char(49)//char(me-10+48)
!  ...
!  elseif(me<290) then
!    mm = char(50)//char(56)//char(me-180+48)
!  elseif(me<300) then
!    mm = char(50)//char(57)//char(me-190+48)
!  endif
mm = trim(adjustl(mm))

!on-fault station coordinates are given here 
!(along strike, z).
xonfs=reshape((/0.0,0.0,&
				4.5,0.0,&
				12.0,0.0,&
				0.0,-4.5,&
				0.0,-7.5,&
				4.5,-7.5,&
				12.0,-7.5,&
				0.0,-12.0/),(/2,8,1/))
xonfs = xonfs * 1000  !convert from km to m

!allocate fault arrays if nonzero fault node. B.D. 10/16/09
!only allocate here for nonzero case. B.D. 10/17/09
!  itmp1 = max(nftnd1(1),nftnd1(2))
!  itmp2 = nonfs(1) + nonfs(2)
!  if(itmp1>0) then
!    allocate(nsmp(2,itmp1,ntotft),fnft(itmp1,ntotft),&
!      un(3,itmp1,ntotft),us(3,itmp1,ntotft),ud(3,itmp1,ntotft),&
!      fric(6,itmp1,ntotft),arn(itmp1,ntotft),r4nuc(itmp1,ntotft),&
!      anonfs(3,itmp2),arn4m(itmp1,ntotft),slp4fri(itmp1,ntotft))
!    nsmp = 0	!initialize here, already in eqdyna3d.f90
!    fnft = 1000. !as failure time initialization: should larger than actual!
!    fric = 0.0
!    un = 0.0
!    us = 1000.0
!    ud = 0.0
!    arn = 0.0
!    arn4m = 0.0
!    r4nuc = 0.0
!    anonfs = 0
!    slp4fri = 0.0
!  endif
 
!convert to m from km for output nodes' coor
x4nds = x4nds * 1000

!dy = dx * abs(dtan(brangle))
dy = dx
dz = dx
tol = dx/100

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
ny = nyuni + edgey1 + iy + nPML
!...pre-determine y-coor
allocate(yline(ny))
!...predetermine y-coor
!yline(edgey1+1) = -dy*(dis4uniF+fltxyz(2,1,1)/dx+1)
yline(edgey1+1) = -dy*(dis4uniF)
ystep = dy
do iy = edgey1, 1, -1
    ystep = ystep * rat    
	yline(iy) = yline(iy+1) - ystep
enddo
ymin = yline(1)
do iy = edgey1+2,edgey1+nyuni
    yline(iy) = yline(iy-1) + dy
enddo
ystep = dy
do iy = edgey1+nyuni+1,ny
    ystep = ystep * rat
	yline(iy) = yline(iy-1) + ystep
enddo
ymax = yline(ny)
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
nz = edgezn + nzuni
!...predetermine z-coor
allocate(zline(nz))
zline(nz) = zmax
do iz = nz-1,nz-nzuni+1,-1
    zline(iz) = zline(iz+1) - dz
enddo
zstep = dz
do iz = nz-nzuni,1,-1
	zstep = zstep * rat
	zline(iz) = zline(iz+1) -zstep
enddo
	zmin = zline(1) !these needed for fixed boundaryies

!...prepare for digitizing
allocate(plane1(ny+ntotft,nz),plane2(ny+ntotft,nz),fltrc(2,nxuni,nzuni,ntotft))
nnode = 0
nelement = 0
n4out = 0
n4onf = 0
plane1 = 0
plane2 = 0
neq = 0
nftnd = 0
an4nds = 0
xnode = 0.0
ien = 0
ixfi = 0
izfi = 0
!DL  
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
			!*.* Comments D.L. Jan/23/2015
			if (ntag>maxm) then
				write(*,*) 'inconsistancy between ntag and maxm'
				stop
			endif
			!*.* D.L.
			xcoor = xline(ix)
			ycoor = yline(iy)
			zcoor = zline(iz)
			!...create nodes
			nnode = nnode + 1		
			plane2(iy,iz) = nnode
			xnode(1,nnode) = xcoor
			xnode(2,nnode) = ycoor
			xnode(3,nnode) = zcoor
			!write(*,*) nnode,xnode(1,nnode),xnode(2,nnode),xnode(3,nnode)
			!*.* Find the PML boundary D.L. Jan/23/2015
			zz=ndof
			if (xcoor>PMLb(1).or.xcoor<PMLb(2).or.ycoor>PMLb(3).or.ycoor<PMLb(4) &
			.or.zcoor<PMLb(5)) then
				zz = 12
			endif		
			locid(nnode)=ntag
			dof1(nnode)=zz
			!		
			!...establish equation numbers for this node
			do i1=1,zz
			!...elastoplastic off-fault, stress assigned in entire model,
			! need to fix model boundaries (except free surface).
				if(abs(xcoor-xmin)<tol.or.abs(xcoor-xmax)<tol.or.abs(ycoor-ymin)<tol &
				.or.abs(ycoor-ymax)<tol.or.abs(zcoor-zmin)<tol) then
					ntag=ntag+1
					id1(ntag)=-1!-1 for fixed boundary nodes
				else
					neq = neq + 1
					ntag=ntag+1
					id1(ntag)=neq
				endif
				!...find two node ranges at edges of the proc for MPI.
				! only do here as faults should not be at edges of iy or iz.
				! B.D. 4/19/09
				! node range for fnms() in faulting, equation range for brhs(), 
				! alhs(). B.D. 5/2/09
				if(ix==1) then
					if(iy==1 .and. iz==1 .and. i1==1) then
						nr4mpi(1,1) = nnode
					elseif(iy==2 .and. iz==2 .and. i1==1) then
						er4mpi(1,1) = neq
					elseif(iy==ny .and. iz==nz .and. i1==zz) then
						nr4mpi(2,1) = nnode
					elseif(iy==ny-1 .and. iz==nz .and. i1==zz) then
						er4mpi(2,1) = neq
					endif
				elseif(ix==nx) then
					if(iy==1 .and. iz==1 .and. i1==1) then
						nr4mpi(1,2) = nnode
					elseif(iy==2 .and. iz==2 .and. i1==1) then
						er4mpi(1,2) = neq
					elseif(iy==ny .and. iz==nz .and. i1==zz) then
						nr4mpi(2,2) = nnode
					elseif(iy==ny-1 .and. iz==nz .and. i1==zz) then
						er4mpi(2,2) = neq
					endif
				endif
			enddo ! enddo i1=1,ndof
			!...identify output nodes (off-fault)
			if(ix>1.and.ix<nx .and. iy>1.and.iy<ny .and. iz==nz) then  !at surface only
				do i=1,n4nds
					if(n4yn(i)==0) then
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
								if(xcoor/=x4nds(1,i)) then
								!...don't reset xcoor,ycoor for degeneration reason.
								! but do reset xnode( ) to conform station coor. B.D. 1/9/12
								!  xcoor=x4nds(1,i)
								!  xnode(1,nnode)=xcoor
									xnode(1,nnode) = x4nds(1,i)
								endif
								if(ycoor/=x4nds(2,i)) then
								!  ycoor=x4nds(2,i)
								!  xnode(2,nnode)=ycoor
									xnode(2,nnode) = x4nds(2,i)
								endif
								exit 	!if node found, jump out the loop
							endif
						endif
					endif
				enddo
			endif
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
						if(nftnd1(ift)==0) then
						write(*,*) 'inconsistency between mesh4num and meshgen for nftnd',me
						write(*,*) 'current node x,y,z:',xcoor,ycoor,zcoor
						stop
						endif
						nsmp(1,nftnd(ift),ift) = nnode !slave node
						nnode = nnode + 1
						!D.L.			 
						locid(nnode)=ntag
						dof1(nnode)=3
						!			 
						nsmp(2,nftnd(ift),ift) = nnode !master node
						plane2(ny+ift,iz) = nnode
						xnode(1,nnode) = xcoor
						xnode(2,nnode) = ycoor
						xnode(3,nnode) = zcoor
						!write(*,*) nnode,xnode(1,nnode),xnode(2,nnode),xnode(3,nnode)
						!...identify output on-fault stations. B.D. 10/25/09
						!...fault-dependent for multiple faults. B.D. 1/8/12
						if(ift==1) then
							tmp1 = xcoor
						elseif(ift==2) then
							tmp1 = sqrt(xcoor*xcoor+ycoor*ycoor)
						endif
						do i1=1,nonfs(ift)
							if(abs(tmp1-xonfs(1,i1,ift))<tol .and. &
							abs(zcoor-xonfs(2,i1,ift))<tol) then
								n4onf = n4onf + 1
								anonfs(1,n4onf) = nftnd(ift)
								anonfs(2,n4onf) = i1
								anonfs(3,n4onf) = ift
								exit
							elseif(abs(zcoor-xonfs(2,i1,ift))<tol .and. &
							abs(tmp1-xonfs(1,i1,ift))<dx/2) then
								n4onf = n4onf + 1
								anonfs(1,n4onf) = nftnd(ift)
								anonfs(2,n4onf) = i1
								anonfs(3,n4onf) = ift
								tmp2=xonfs(1,i1,ift)*dcos(brangle)
								tmp3=-xonfs(1,i1,ift)*dsin(brangle)
								xnode(1,nnode)=tmp2 !move mesh to stations!
								xnode(1,nnode-1)=tmp2
								xnode(2,nnode)=tmp3
								xnode(2,nnode-1)=tmp3
								exit          
							endif
						enddo  
						!...establish equation numbers for this master node
						do i1=1,ndof
							neq = neq + 1
							ntag=ntag+1
							id1(ntag)=neq!assume no explicit boundary nodes
						enddo
						!...assign unit vectors to the split node pair	
						un(1,nftnd(ift),ift) = dcos(fltxyz(1,4,ift))*dsin(fltxyz(2,4,ift))
						un(2,nftnd(ift),ift) = -dsin(fltxyz(1,4,ift))*dsin(fltxyz(2,4,ift))
						un(3,nftnd(ift),ift) = dcos(fltxyz(2,4,ift))		
						us(1,nftnd(ift),ift) = -dsin(fltxyz(1,4,ift))
						us(2,nftnd(ift),ift) = -dcos(fltxyz(1,4,ift))
						us(3,nftnd(ift),ift) = 0.0
						ud(1,nftnd(ift),ift) = dcos(fltxyz(1,4,ift))*dcos(fltxyz(2,4,ift))
						ud(2,nftnd(ift),ift) = dsin(fltxyz(1,4,ift))*dcos(fltxyz(2,4,ift))
						ud(3,nftnd(ift),ift) = dsin(fltxyz(2,4,ift))
						!...prepare for area calculation
						if(ixfi(ift)==0) ixfi(ift)=ix
						if(izfi(ift)==0) izfi(ift)=iz
						ifs(ift)=ix-ixfi(ift)+1
						ifd(ift)=iz-izfi(ift)+1
						fltrc(1,ifs(ift),ifd(ift),ift) = nnode	!master node
						fltrc(2,ifs(ift),ifd(ift),ift) = nftnd(ift) !fault node num in sequence
						!...assign friction parameters for SCEC TPV19. B.D. 1/8/12
						!...now for Ma and Andrews (2010) model. B.D. 6/1/12
						fric(1,nftnd(ift),ift) = 0.760!10000.!mus
						fric(2,nftnd(ift),ift) = 0.448   !mud
						fric(3,nftnd(ift),ift) = 0.5    !D0
						fric(4,nftnd(ift),ift) = 1.0e6  !cohesion
						fric(5,nftnd(ift),ift) = 0.08  !Not used in TPV8 !T for forced rupture,i.e.,nucleation
						fric(6,nftnd(ift),ift) = 0.0!rhow*grav*(-zcoor) !pore pressure	
!------------------------------ZONE I-------------------------------!
!------------Assign initial stresses for elastic version------------!
!-----------------------------EQ V3.2.1-----------------------------!	
!-----------------------Sep.19.2015/ D.Liu--------------------------!
!-------------------------------TPV 8-------------------------------!						
						if (C_elastic==1) then
							fric(7,nftnd(ift),ift) = 7378.0*zcoor! Initial normal stress for ElasticVersion
							fric(8,nftnd(ift),ift) = abs(0.55*7378.0*zcoor)!Initial shear stress for ElasticVersion
							if (abs(xcoor)<=1500.and.abs(zcoor--12e3)<=1500.)then
								fric(8,nftnd(ift),ift) =1.0e6+1.005*0.760*abs(7378.0*zcoor)
							endif
						endif
!-----------------------------END ZONE I----------------------------!						
						!special values below.						
						if(abs(xcoor-fltxyz(1,1,1))<tol .or. abs(xcoor-fltxyz(2,1,ift))<tol &
						.or. abs(zcoor-fltxyz(1,3,ift))<tol) then !-x,+x,-z for 1, +x,-z for 2
							fric(1,nftnd(ift),ift) = 10000.	!fault edge, pinned
						endif
						!if(abs(zcoor)<tol) fric(6,nftnd(ift),ift) = 0.5 * dz  !surface pore pressure
						!if(ift==2)  then !branch, fault 2, only 12 km strike rupturable
						!  tmp1=sqrt(xcoor**2+ycoor**2)
						!  if(tmp1>12000) fric(1,nftnd(ift),ift) = 10000.
						!endif               
						exit !can only be on 1 fault, thus if ynft(ift), exit do loop       
					endif  !if ynft(ift)
				endif  !if flt range
			enddo  !do ift 
			!...create elements
			if(ix>=2 .and. iy>=2 .and. iz>=2) then
				nelement = nelement + 1
				if(nelement>nelement1) then
					write(*,*) 'more elements in meshgen than in mesh4num'
					write(*,*) 'x,y,z',xcoor,ycoor,zcoor
					pause
				endif
				mat(nelement) = 1  !homogeneous half space
				et(nelement) = 1 !brick element. D.L. Jan/23/2015
				ien(1,nelement) = plane1(iy-1,iz-1)
				ien(2,nelement) = plane2(iy-1,iz-1)
				ien(3,nelement) = plane2(iy,iz-1)
				ien(4,nelement) = plane1(iy,iz-1)
				ien(5,nelement) = plane1(iy-1,iz)
				ien(6,nelement) = plane2(iy-1,iz)
				ien(7,nelement) = plane2(iy,iz)
				ien(8,nelement) = plane1(iy,iz)
				!*.* Compute et for PML element. D.L. Jan/23/2015
				xc=0.0
				label=0
				ids(nelement)=ntags
				do i=1,8
					do j=1,3
					!label=label+dof1(ien(i,nelement)) !label=24,corner:
						xc(j)=xc(j)+xnode(j,ien(i,nelement))
					enddo
				enddo
				xc=xc/8
				!if (label>24) then
				if (xc(1)>PMLb(1).or.xc(1)<PMLb(2).or.xc(2)>PMLb(3).or.xc(2)<PMLb(4) &
					.or.xc(3)<PMLb(5)) then
					et(nelement) = 2
					ntags=ntags+15+6
				else
					ntags=ntags+12
				endif
!------------------------------ZONE II------------------------------!
!Assign material types according to coordinates of center of the element!
!-----------------------------EQ V3.2.1-----------------------------!	
!-----------------------Sep.19.2015/ D.Liu--------------------------!
			!For Double-couple point source. Ma and Liu (2006) 				
				! if (xc(3)<-1000.)then
					! mat(nelement)=2
				! endif
!----------------------------END ZONE II----------------------------!
				!special treatment for using master nodes above the main fault.
				! B.D. 1/7/12
				if((xcoor>(fltxyz(1,1,1)-tol).and.xcoor<(fltxyz(2,1,1)+dx+tol).and. &
				zcoor>(fltxyz(1,3,1)-tol).and.ycoor>0.and.abs(ycoor-dy)<tol)) then
					do i=1,nftnd(1)
						do k=1,8
							if(ien(k,nelement)==nsmp(1,i,1)) then
								ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
							endif
						enddo
					enddo
				endif
				!...assign initial stress and pore pressure to element for inelastic
				! off-fault rheology. for SCEC TPV18-21. B.D. 1/8/12
				!...now change to Ma & Andrews (2010) models. B.D. 6/1/12
!---------------------------ZONE III--------------------------------!
!------------Assign initial stresses for plastic version------------!
!-----------------Use 1D s1 to store initial stress-----------------!
!-----------------------Jan.23.2015/ D.Liu--------------------------!
!-----------------------------EQ V3.2.1-----------------------------!	
!-----------------------Sep.19.2015/ D.Liu--------------------------!	
				if (C_elastic==0) then				
					tmp1 = -0.5*(zline(iz)+zline(iz-1))  !z<0, thus tmp1>0
					tmp2 = tmp1 * grav
					if(et(nelement)==1)then
						eleporep(nelement) = rhow*tmp2  !pore pressure>0
						s1(ids(nelement)+3)= -rho(mat(nelement))*tmp2  !vertical, comp<0
						s1(ids(nelement)+1)=s1(ids(nelement)+3)
						s1(ids(nelement)+2)=s1(ids(nelement)+3)
						s1(ids(nelement)+6)=-0.4 * (s1(ids(nelement)+3)+eleporep(nelement))
					elseif(et(nelement)==2)then
						eleporep(nelement) = rhow*tmp2  !pore pressure>0
						s1(ids(nelement)+3+15)= -rho(mat(nelement))*tmp2  !vertical, comp<0
						s1(ids(nelement)+1+15)=s1(ids(nelement)+3+15)
						s1(ids(nelement)+2+15)=s1(ids(nelement)+3+15)
						s1(ids(nelement)+6+15)=-0.4 * (s1(ids(nelement)+3+15)+eleporep(nelement))				
					endif
				endif
!--------------------------END ZONE III-----------------------------!				
				!...when the current node is one element above the branch fault in y-coor,
				! one hexahedron degenerates into two wedges. B.D. 1/7/12
				!          if(xcoor>=(fltxyz(1,1,2)-tol).and.xcoor<=(fltxyz(2,1,2)+tol).and. &
				!            ycoor>(fltxyz(1,2,2)+dy-tol).and.ycoor<=(fltxyz(2,2,2)+dy+tol).and. &
				!            zcoor>=(fltxyz(1,3,2)-tol).and.zcoor<=(fltxyz(2,3,2)+tol)) then
				!            !above is x,y,z ranges for possible degeneation. B.D. 1/7/12
				!            if(abs(ycoor+xcoor*dtan(brangle)-dy)<tol) then !degenerate
				!              ien(1,nelement) = plane1(iy,iz-1)
				!              ien(2,nelement) = plane1(iy-1,iz-1)
				!              ien(3,nelement) = plane2(iy-1,iz-1)
				!              ien(4,nelement) = ien(3,nelement)
				!              ien(5,nelement) = plane1(iy,iz)
				!              ien(6,nelement) = plane1(iy-1,iz)
				!              ien(7,nelement) = plane2(iy-1,iz)
				!              ien(8,nelement) = ien(7,nelement) 
				!              nelement = nelement + 1 !one more wedge element
				!              mat(nelement) = 1
				!              ien(1,nelement) = plane2(iy-1,iz-1)
				!              ien(2,nelement) = plane2(iy,iz-1)
				!              ien(3,nelement) = plane1(iy,iz-1)
				!              ien(4,nelement) = ien(3,nelement)
				!              ien(5,nelement) = plane2(iy-1,iz)
				!              ien(6,nelement) = plane2(iy,iz)
				!              ien(7,nelement) = plane1(iy,iz)
				!              ien(8,nelement) = ien(7,nelement)
				!              !assign initial stress & pore pressure for this element. B.D. 1/8/12
				!              !same depth as nelement-1; only for possible nonzero components.
				!              eleporep(nelement) = eleporep(nelement-1)
				!              elestress(1,nelement) = elestress(1,nelement-1)
				!              elestress(2,nelement) = elestress(2,nelement-1)
				!              elestress(3,nelement) = elestress(3,nelement-1)
				!              elestress(6,nelement) = elestress(6,nelement-1)
				!            endif
				!          endif
				!          !special treatment for using master nodes above the branch fault (in y-coor).
				!          ! both degenerated immediately above branch & regular one more ebove. B.D. 1/13/12
				!          if(xcoor>(fltxyz(1,1,2)+tol).and.xcoor<=(fltxyz(2,1,2)+dx+tol).and. &
				!            zcoor>(fltxyz(1,3,2)-tol).and.ycoor>(-xcoor*dtan(brangle)+tol).and.&
				!            ycoor<-xcoor*dtan(brangle)+3*dy) then
				!            do i=1,nftnd(2)
				!              do k=1,8
				!                if(ien(k,nelement)==nsmp(1,i,2)) then
				!                  ien(k,nelement) = nsmp(2,i,2)  !use master node for the node!
				!                endif
				!              enddo
				!            enddo
				!          	endif
			endif!if element
		enddo!iy
    enddo!iz
    plane1 = plane2
enddo!ix
maxs=ntags
if (maxs>=(4*maxm)) then
	write(*,*) '4*maxm',maxm,'is not enough for maxs',maxs
	stop
endif
!   write(*,*) 'nnd,nele,nftnd,neq',nnode,nelement,nftnd,neq
!   write(*,*) 'nx,ny,nz,dx,dy,dz',nx,ny,nz,dx,dy,dz
   
  !...verify consistence with mesh4num.f90. B.D. 4/30/09
  ! for multiple faults now. 1/7/12
if(nnode/=nnode1 .or. nelement/=nelement1 .or. neq/=neq1) then
	write(*,*) 'Inconsistency in node/element/equation between meshgen and mesh4num: stop!',me
	write(*,*) nnode,nnode1,nelement,nelement1,neq,neq1
	stop
endif
if(ntag/=maxm) then
	write(*,*) 'Inconsistency in ntag and maxm: stop!',me
	write(*,*) ntag,maxm
	stop
endif
do i=1,ntotft
	if(nftnd(i)/=nftnd1(i)) then
		write(*,*) 'Inconsistency in fault between meshgen and mesh4num: stop!',me,i
		stop
	endif
enddo

!...calculate area associated with each fault node pair and distance from source. 
!B.D. 4/19/09
!but only needed for nftnd>0. B.D. 10/16/09
!for multiple faults. B.D. 1/7/12
do ift=1,ntotft
	if(nftnd(ift)>0) then
	!...distance from the source point
	do i=1,ifd(ift)
		do j=1,ifs(ift)
			i1 = fltrc(1,j,i,ift)
			j1 = fltrc(2,j,i,ift)
			r4nuc(j1,ift) = sqrt((xnode(1,i1)-xsource)**2 + (xnode(2,i1)-ysource)**2 &
			+ (xnode(3,i1)-zsource)**2)
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
			m1 = fltrc(2,j,i,ift) !use nodal number in nftnd
			m2 = fltrc(2,j-1,i,ift)
			m3 = fltrc(2,j-1,i-1,ift)
			m4 = fltrc(2,j,i-1,ift)
			!...calculate area of the quadrilateral
			!......if fault is not in coor axes plane
			a1=sqrt((xnode(1,n2)-xnode(1,n1))**2 + (xnode(2,n2)-xnode(2,n1))**2 &
			+ (xnode(3,n2)-xnode(3,n1))**2)
			b1=sqrt((xnode(1,n3)-xnode(1,n2))**2 + (xnode(2,n3)-xnode(2,n2))**2 &
			+ (xnode(3,n3)-xnode(3,n2))**2)
			c1=sqrt((xnode(1,n4)-xnode(1,n3))**2 + (xnode(2,n4)-xnode(2,n3))**2 &
			+ (xnode(3,n4)-xnode(3,n3))**2)
			d1=sqrt((xnode(1,n1)-xnode(1,n4))**2 + (xnode(2,n1)-xnode(2,n4))**2 &
			+ (xnode(3,n1)-xnode(3,n4))**2)
			p1=sqrt((xnode(1,n4)-xnode(1,n2))**2 + (xnode(2,n4)-xnode(2,n2))**2 &
			+ (xnode(3,n4)-xnode(3,n2))**2)
			q1=sqrt((xnode(1,n3)-xnode(1,n1))**2 + (xnode(2,n3)-xnode(2,n1))**2 &
			+ (xnode(3,n3)-xnode(3,n1))**2)
			area=0.25 * sqrt(4*p1*p1*q1*q1 - &
      	   (b1*b1 + d1*d1 - a1*a1 -c1*c1)**2) 
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
	endif
	!
	time1 = MPI_WTIME()
	!*******************MPI***************************
	!...need to deal with edge nodes of two MPI regions for area.
	!   make use of fltrc(2,i,ifd) for i=1 and ifs as edges. B.D. 4/20/09
	! mpi_send arn() for fltrc(2,1,1:ifd) to me-1, for fltrc(2,ifs,1:ifd) to me+1;
	! sum current me's arn() of fltrc(2,ifs,1:ifd) and received fltrc(2,1,1:ifd) from 
	! me+1, and sum me's arn() of fltrc(2,1,1:ifd) and received fltrc(2,ifs,1:ifd) from
	! me-1.
	!...need to deal with special case: either me or me-1/me+1 has zero fault node.
	! B.D. 10/17/09 
	if (nprocs > 1) then
		if (me == master) then
			call mpi_sendrecv(nftnd(ift), 1, MPI_INTEGER, 1, 1000+me, &
                        itmp1, 1, MPI_INTEGER, 1, 1000+me+1, &
        		MPI_COMM_WORLD, istatus, ierr)
			if(nftnd(ift)>0 .and. itmp1>0) then  
				allocate(btmp(ifd(ift)),btmp1(ifd(ift)))
				do i = 1, ifd(ift)
					btmp(i)=arn(fltrc(2,ifs(ift),i,ift),ift)
				enddo
				call mpi_sendrecv(btmp, ifd(ift), MPI_DOUBLE_PRECISION, 1, 1000+me, &
                        btmp1, ifd(ift), MPI_DOUBLE_PRECISION, 1, 1000+me+1, &
                        MPI_COMM_WORLD, istatus, ierr)
				do i=1, ifd(ift)
					arn(fltrc(2,ifs(ift),i,ift),ift) = arn(fltrc(2,ifs(ift),i,ift),ift) + btmp1(i)
				enddo
				deallocate(btmp,btmp1)       
			endif
		elseif (me == nprocs-1) then
			call mpi_sendrecv(nftnd(ift), 1, MPI_INTEGER, nprocs-2, 1000+me, &
                        itmp1, 1, MPI_INTEGER, nprocs-2, 1000+me-1, &
        		MPI_COMM_WORLD, istatus, ierr)
			if(nftnd(ift)>0 .and. itmp1>0) then
				allocate(btmp(ifd(ift)),btmp1(ifd(ift)))
				do i = 1, ifd(ift)
					btmp(i)=arn(fltrc(2,1,i,ift),ift)
				enddo
				call mpi_sendrecv(btmp, ifd(ift), MPI_DOUBLE_PRECISION, nprocs-2, 1000+me, &
                        btmp1, ifd(ift), MPI_DOUBLE_PRECISION, nprocs-2, 1000+me-1, &
                        MPI_COMM_WORLD, istatus, ierr)
				do i=1, ifd(ift)
					arn(fltrc(2,1,i,ift),ift) = arn(fltrc(2,1,i,ift),ift) + btmp1(i)
				enddo
				deallocate(btmp,btmp1)       
			endif
		elseif (me > 0 .and. me < nprocs-1) then
			call mpi_sendrecv(nftnd(ift), 1, MPI_INTEGER, me-1, 1000+me, &
                       itmp1, 1, MPI_INTEGER, me-1, 1000+me-1,&
                       MPI_COMM_WORLD, istatus, ierr)
			call mpi_sendrecv(nftnd(ift), 1, MPI_INTEGER, me+1, 1000+me, &
                       itmp2, 1, MPI_INTEGER, me+1, 1000+me+1, &
                       MPI_COMM_WORLD, istatus, ierr)
			if(nftnd(ift)>0 .and. itmp1>0) then
				allocate(btmp(ifd(ift)), btmp1(ifd(ift)))
				do i = 1,ifd(ift)
					btmp(i)=arn(fltrc(2,1,i,ift),ift)
				enddo
				call mpi_sendrecv(btmp, ifd(ift), MPI_DOUBLE_PRECISION, me-1, 1000+me, &
                         btmp1, ifd(ift), MPI_DOUBLE_PRECISION, me-1, 1000+me-1,&
                         MPI_COMM_WORLD, istatus, ierr)
				do i=1,ifd(ift)
					arn(fltrc(2,1,i,ift),ift) = arn(fltrc(2,1,i,ift),ift) + btmp1(i)
				enddo
				deallocate(btmp,btmp1)       
			endif
			if(nftnd(ift)>0 .and. itmp2>0) then
				allocate(btmp(ifd(ift)),btmp1(ifd(ift)))
				do i = 1,ifd(ift)
					btmp(i)=arn(fltrc(2,ifs(ift),i,ift),ift)
				enddo
				call mpi_sendrecv(btmp, ifd(ift), MPI_DOUBLE_PRECISION, me+1, 1000+me, &
                         btmp1, ifd(ift), MPI_DOUBLE_PRECISION, me+1, 1000+me+1, &
                         MPI_COMM_WORLD, istatus, ierr)
				do i=1,ifd(ift)
					arn(fltrc(2,ifs(ift),i,ift),ift) = arn(fltrc(2,ifs(ift),i,ift),ift) + btmp1(i)
				enddo
				deallocate(btmp,btmp1)       
			endif
		endif
	endif
!
enddo !do ift 
!
time2 = MPI_WTIME()
btime = btime + (time2 - time1) 
!  if(me==master) then
!  write(*,*) 'me= ',me,'in meshgen, after MPI for arn()'
!  endif  
  
!  write(*,*) 'Entire model region has:'
!  write(*,*) nnode,' nodes; ',nelement,' elements.'
!  write(*,*) 'outnodes all or not: '
!  write(*,'( 12i3)') (n4yn(i),i=1,n4nds)

!***
!...write out mesh data for matlab to plot
!***
!  open(11,file='g1km/vert.txt',status='unknown')
!  write(11,'(1x,3f13.2)') ((xnode(j,i),j=1,3),i=1,nnode)
!  close(11)
!  open(12,file='g1km/fac.txt',status='unknown')
!  write(12,'(1x,4i10)') (((fac(i,j,k),j=1,4),i=1,6),k=1,nelement)
!  close(12)

!...write out mesh information, including output node coor for check. B.D. 10/27/07
meshfile = 'meshinfo.txt'//mm 
open(33,file=meshfile,status='unknown')
write(33,*) 'me= ',me
write(33,*) 'x-range for this MPI:',xline(1),xline(nx)
write(33,'( a45,5i10)') 'total in this process:node,ele,neq,nftnd ',nnode,nelement,neq,(nftnd(i),i=1,ntotft)
write(33,'( a60)') 'off-fault stations # in order and x, y, z- xoor'
write(33,'( 2i10,2f12.2)') ((an4nds(j,i),j=1,2),(x4nds(j,an4nds(1,i)), &
  	j=1,2),i=1,n4out)
write(33,'( a65,i5)') 'on-fault stations # in order and strike and down-dip distance',n4onf
write(33,'( 3i10,2f12.2)') ((anonfs(j,i),j=1,3),(xonfs(j,anonfs(2,i),anonfs(3,i)),j=1,2), &
	i=1,n4onf)
!write(33,*) 'branch fault'
!write(33,'( 10f8.1)') (arn(i,2),i=1,nftnd(2))
close(33)

!if(me==73) then
!open(31,file='area73.txt',status='unknown')
!write(31,*) 'main fault'
!write(31,'( 10f8.1)') (arn(i,1),i=1,nftnd(1))
!write(31,*) 'branch fault'
!write(31,'( 10f8.1)') (arn(i,2),i=1,nftnd(2))
!close(31)
!endif
!...deallocate working arrays to free memory afte mesh. B.D. 10/28/09
! actually doen't matter. they will free up after exit!
deallocate(xlinet,xline,yline,zline,plane1,plane2,fltrc)
end subroutine meshgen
