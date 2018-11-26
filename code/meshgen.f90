subroutine meshgen(nnode1,nelement1,ien,mat,s1,eleporep,xnode,neq1,id1,maxm,locid,&
				dof1,et,PMLb,maxs,ids,n4out,nftnd1,n4onf,nsmp,un,us,ud,fric,arn,&
				r4nuc,anonfs,arn4m,me,nsurnd1,surnode,surid,miuonf,vponf)
use globalvar
implicit none
include 'mpif.h'
integer(kind=4)::nnode,nelement,nnode1,nelement1,nxt,nyt,nzt,nx,ny,nz,ix,iy,iz, &
edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,np=10000,edgezn
integer(kind=4)::neq,neq1,n4onf,n4out,n4yn(n4nds),anonfs(3,nonmx),&
	nsmp(2,nftmx,ntotft),ien(8,nelement1),id1(maxm),locid(nnode1),dof1(nnode1),&
	et(nelement1),ids(nelement1)
integer(kind=4)::nxuni,nyuni,nzuni,ift,n1,n2,n3,n4,m1,m2,m3,m4,alloc_err,&
	istatus(MPI_STATUS_SIZE),bndl,bndr,bndf,bndb,bndd,bndu,me,mex,mey,mez,rlp,rr,ierr,jj,itmp1,&
	maxm,maxs,label,ntag,ntags,zz,ivp1,ivp2,ixe,iye,ize,itemp,iye1,nxe,nye,nze,nsurnd1,nsurnd,&
	surid(nsurnd1),nsx,nsy,nfx,nfz,msnode
integer(kind=4),dimension(ntotft)::nftnd,nftnd1,ixfi,izfi,ifs,ifd
integer(kind=4),allocatable::plane1(:,:),plane2(:,:),fltrc(:,:,:,:)

real(kind=8)::mat(nelement1,5),eleporep(nelement1),xnode(ndof,nnode1),&
	fric(20,nftmx,ntotft),PMLb(8),xc(3),s1(4*maxm),surnode(nsurnd1,2),miuonf(nftmx),vponf(nftmx),Dampx,Dampz
real(kind=8),dimension(nftmx,ntotft)::fnft,arn,r4nuc,arn4m,slp4fri
real(kind=8),dimension(3,nftmx,ntotft)::un,us,ud
real(kind=8)::tol,xcoor,ycoor,zcoor,xstep,ystep,zstep,tmp1,tmp2,tmp3,a,b,area,a1,b1,c1,d1,p1,q1
real(kind=8),allocatable,dimension(:)::xlinet,ylinet,zlinet,xline,yline,zline,btmp,btmp1,btmp2,btmp3
real(kind=8),allocatable::vpstruct(:,:)
!Heteogeneous stress setup
real(kind=8)::initialinput(3,91001)
!...for initial stress and pore pressure of elements: inelastic. B.D. 1/8/12
real (kind=8)::grav=9.8,omega
logical,dimension(ntotft) :: ynft
character(len=30)::meshfile,mm
character(len=100)::fname
!-------------------------------------------------------------------!
!------------Early Oct.2015/ D.Liu----------------------------------!
!------------Input of initial stress for TSN------------------------!
! open(3001,file='0InputFEM.txt',form='formatted',status='old')
! do i=1,91001
! read(3001,*) (initialinput(j,i),j=1,3)
! enddo
! close(3001)
!-------------------------------------------------------------------!
write(mm,'(i6)') me
mm = trim(adjustl(mm))
!===================================================================!
!==============================ZONE A===============================!
!on-fault station coordinates (along strike, z).
xonfs=reshape((/0.0,-3.0,&
				0.0,-7.5,&
				0.0,-12.0,&
				9.0,-7.5,&
				12.0,-3.0,&
				12.0,-12.0,&
				-9.0,-7.5,&
				-12.0,-3.0,&
				-12.0,-12.0/),(/2,9,1/))
!off-fault station coordinates(along strike,normal,z(negative)).				
x4nds=reshape((/0.0,9.0,0.0,&
				0.0,-9.0,0.0,&
				12.0,6.0,0.0,&
				12.0,-6.0,0.0,&
				-12.0,6.0,0.0,&
				-12.0,-6.0,0.0/),(/3,6/))
!==========================END ZONE A===============================!
!===================================================================!
!==============================ZONE B===============================!
!Sep.21.2015. D.Liu/ Loading vp structures--------------------------!
	! write(*,*) 'nxe,me',nxe,me
	! write(*,*) 'limits',xline(1),xline(nx)
	! write(*,*) (xline(1)+250-(xlinet(1)+250))/500+1,(xline(nx)+250-(xlinet(1)+250))/500+1
	!nxe=nx-1
	!nye=ny-1 
	!nze=nz-1!!Probematic!!!Should change later.
	! fname="/scratch/dunyuliu/256coreVP#9/vpstruct"//mm
	! allocate(vpstruct(nxe*5250,10))
	! open(3002,file=fname,form='formatted',status='old')
	! do ivp1=1,nxe*105*50
		! read(3002,*) (vpstruct(ivp1,ivp2),ivp2=1,10)
	! enddo
	! close(3002)
	! write(*,*) 'Finish loading vp',me
!==========================END ZONE B===============================!
!===================================================================!
xonfs=xonfs*1000.  !convert from km to m
x4nds=x4nds*1000.
n4yn=0
!dy = dx * abs(dtan(brangle))
dy=dx
dz=dx
tol=dx/100

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
!	B.D. 1/19/10
j1 = nxt + npx - 1
rlp = j1/npx
rr = j1 - rlp*npx
if(mex<(npx-rr)) then
    nx = rlp
else
    nx = rlp + 1	!evenly distributed to last rr
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
    ny = rlp + 1	!evenly distributed to last rr
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
    nz = rlp + 1	!evenly distributed to last rr
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
n4out = 0
n4onf = 0
plane1 = 0
plane2 = 0
neq = 0
nftnd = 0
nsurnd=0
an4nds = 0
xnode = 0.0
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
nsx=(surxmax-surxmin)/dx+1
nsy=(surymax-surymin)/dx+1
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
			xnode(1,nnode) = xcoor
			xnode(2,nnode) = ycoor
			xnode(3,nnode) = zcoor
			if (abs(xcoor-3300.)<tol.and.abs(ycoor-5700.)<tol.and.abs(zcoor--6000.)<tol)then 	
				write(*,*) 'TestNodeCOrner',me,mex,mey,mez,nnode
			endif
			if (abs(xcoor-3300.)<tol.and.abs(ycoor-400.)<tol.and.abs(zcoor--11500.)<tol)then 	
				write(*,*) 'TestNodePMML',me,mex,mey,mez,nnode
			endif			
			if (abs(xcoor-3200.)<tol.and.abs(ycoor-5600.)<tol.and.abs(zcoor--5900.)<tol)then 	
				write(*,*) 'TestNodeIn',me,mex,mey,mez,nnode
			endif			
			if (zcoor==0.0)then
				do i=1,nsx
					do j=1,nsy
						if (xcoor==(surxmin+(i-1)*dx).and.ycoor==(surymin+(j-1)*dx)) then 
							nsurnd=nsurnd+1
							surid(nsurnd)=nnode							
							surnode(nsurnd,1)=xcoor 
							surnode(nsurnd,2)=ycoor 							
						endif
					enddo
				enddo
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
					neq = neq + 1
					ntag=ntag+1
					id1(ntag)=neq
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
									exit 	!if node found, jump out the loop
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
									exit 	!if node found, jump out the loop
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
									exit 	!if node found, jump out the loop
								endif
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
						stop 2001
						endif
						nsmp(1,nftnd(ift),ift) = nnode !slave node
               			msnode=nx*ny*nz+nftnd(ift)								 
						locid(msnode)=ntag
						dof1(msnode)=3		 
						nsmp(2,nftnd(ift),ift) = msnode !master node
						plane2(ny+ift,iz) = msnode
						xnode(1,msnode) = xcoor
						xnode(2,msnode) = ycoor
						xnode(3,msnode) = zcoor
						!write(*,*) nnode,xnode(1,nnode),xnode(2,nnode),xnode(3,nnode)
               			!3DMPI
						if(ix == 1) then
                 			fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 1
                  			fltnum(1) = fltnum(1) + 1
               			endif
               			if(ix == nx) then
               			   	fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 2
             			 	fltnum(2) = fltnum(2) + 1
               			endif
               			if(iy == 1) then
                  			fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 10
                  			fltnum(3) = fltnum(3) + 1
               			endif
               			if(iy == ny) then
                  			fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 20
                  			fltnum(4) = fltnum(4) + 1
               			endif
               			if(iz == 1) then
                  			fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 100
                  			fltnum(5) = fltnum(5) + 1
           				endif
               			if(iz == nz) then
                  			fltgm(nftnd(ift)) = fltgm(nftnd(ift)) + 200
                  			fltnum(6) = fltnum(6) + 1
               			endif             
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
						fltrc(1,ifs(ift),ifd(ift),ift) = msnode	!master node
						fltrc(2,ifs(ift),ifd(ift),ift) = nftnd(ift) !fault node num in sequence
!------------------------------ZONE I-------------------------------!
!------------Assign initial stresses for elastic version------------!
!----------------------EQ V4.1-2016/10/03 D.Liu---------------------!
!-------------------------------TPV 27v2-------------------------------!
						fric(1,nftnd(ift),ift) = 0.18
						fric(2,nftnd(ift),ift) = 0.12
						! if (xcoor<-20e3.and.xcoor>-40e3.and.zcoor<-3e3.and.zcoor>-8e3)then
							! fric(2,nftnd(ift),ift) = 0.5   !mud
						! endif						
						fric(3,nftnd(ift),ift) = 0.30    !D0
						fric(4,nftnd(ift),ift) = 0.0  !cohesion
						fric(5,nftnd(ift),ift) = 0.03  !Viscoplastic relaxation time
						fric(6,nftnd(ift),ift) = 0.0!pore pressure						
						if (C_elastic==1) then
							! fric(7,nftnd(ift),ift) = 7378.0*zcoor! Initial normal stress for ElasticVersion
							! if (abs(zcoor)<tol) fric(7,nftnd(ift),ift)=-7378.0*dx/4
							! fric(8,nftnd(ift),ift) = 0.55*abs(fric(7,nftnd(ift),ift))!Initial shear stress for ElasticVersion
							! if (abs(xcoor)<=1500.and.abs(zcoor--12e3)<=1500.)then
								! fric(8,nftnd(ift),ift) =1.0e6+1.005*0.760*abs(7378.0*zcoor)
							! endif
							fric(7,nftnd(ift),ift)=-120e6
							fric(8,nftnd(ift),ift)=40e6
						endif
!-------------------------------Tianjin ref-------------------------------!						
						! if (C_elastic==1) then
							! fric(7,nftnd(ift),ift) = 9.8*1700*zcoor! Initial normal stress for ElasticVersion
							! if (zcoor==0.0) then 
								! fric(7,nftnd(ift),ift) = -9.8*1700*dx/4
							! endif
							! if (abs(xcoor)<60e3) then
								! Dampx=1.0
							! else
								! Dampx=1+((abs(xcoor)-60e3)/30e3)**4
							! endif
							! Dampz=1+(abs(zcoor)/30e3)**4
							! fric(8,nftnd(ift),ift) = (0.6+0.05)*abs(fric(7,nftnd(ift),ift))/&
								! Dampz/Dampx!Initial shear stress for ElasticVersion
						! endif
!------------------------------Tianjin Hete-------------------------!
!-----------------------------Oct.4.2015/ D.Liu---------------------!
						! if (C_elastic==1) then
							! nfx=(xcoor--90e3)/dx+1
							! nfz=(zcoor--20e3)/dx+1
							! fric(8,nftnd(ift),ift)=initialinput(1,(nfz-1)*901+nfx)
							! fric(7,nftnd(ift),ift)=-fric(8,nftnd(ift),ift)&
							! /initialinput(2,(nfz-1)*901+nfx)						
							! fric(1,nftnd(ift),ift)=initialinput(3,(nfz-1)*901+nfx)
						! endif						
!-----------------------------END ZONE I----------------------------!
!-----------------------------ZONE IV RSF---------------------------!
!---2016.08.28. Bin
!TPV104 2016.10.05 D.Liu
						if (friclaw==3.or.friclaw==4)then 
							if (abs(xcoor)<=15e3) then 
								tmp1=1.0
							elseif ((abs(xcoor)<18e3).and.(abs(xcoor)>15e3)) then 
								tmp1=0.5*(1+dtanh(3e3/(abs(xcoor)-18e3)+3e3/(abs(xcoor)-15e3)))
							else
								tmp1=0.0
							endif
							if (abs(zcoor--7.5e3)<=7.5e3) then 
								tmp2=1.0
							elseif ((abs(zcoor--7.5e3)<10.5e3).and.(abs(zcoor--7.5e3)>7.5e3)) then 
								tmp2=0.5*(1+dtanh(3e3/(abs(zcoor--7.5e3)-10.5e3)+3e3/(abs(zcoor--7.5e3)-7.5e3)))
							else
								tmp2=0.0
							endif							
							fric(9,nftnd(ift),ift)=0.01+0.01*(1-tmp1*tmp2)!a
							if (xcoor==17e3.and.zcoor==-17e3) then 
								write(*,*) 'tmp1',fric(9,nftnd(ift),ift),tmp1,tmp2 
							endif
							fric(10,nftnd(ift),ift)=0.014!b 
							fric(11,nftnd(ift),ift)=0.4d0!RSF critical distance.
							fric(12,nftnd(ift),ift)=1.0d-6!RSF:V0
							fric(13,nftnd(ift),ift)=0.6d0!RSF:miu0
							if(friclaw==4) then 
								fric(14,nftnd(ift),ift)  = 0.2d0 !RSF: fw for strong rate weakenging
								fric(15,nftnd(ift),ift)  = 0.1d0+0.9*(1-tmp1*tmp2) !RSF: Vw for strong rate weakening
							endif	
							fric(16,nftnd(ift),ift)=0.0d0 !RSF: initial normal slip rate 
							fric(17,nftnd(ift),ift)=1.0e-16!RSF:s 
							fric(18,nftnd(ift),ift)=0.0d0!RSF:d
							fric(19,nftnd(ift),ift)=1.0e-16!RSF:mag							
						endif
!-----------------------------END ZONE IV---------------------------!						
						!special values below.						
						if(abs(xcoor-fltxyz(1,1,ift))<tol .or. abs(xcoor-fltxyz(2,1,ift))<tol &
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
                if((xcoor==-5e3).and.(ycoor==0).and.(zcoor==-10e3))then
					write(*,*) 'me,nele',me,nelement
				endif
				nelement = nelement + 1
				if(nelement>nelement1) then
					write(*,*) 'more elements in meshgen than in mesh4num'
					write(*,*) 'x,y,z',xcoor,ycoor,zcoor
					pause
				endif
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
						xc(j)=xc(j)+xnode(j,ien(i,nelement))
					enddo
				enddo
				xc=xc/8
				if (xc(1)>PMLb(1).or.xc(1)<PMLb(2).or.xc(2)>PMLb(3).or.xc(2)<PMLb(4) &
					.or.xc(3)<PMLb(5)) then
					et(nelement) = 2
					ntags=ntags+15+6
				else
					ntags=ntags+12
				endif
				! if (abs(xc(1)-3250)<tol.and.abs(xc(2)-350)<tol.and.abs(xc(3)--11450)<tol)then
					! write(*,*) 'TestElem',me,nelement
				! endif
!------------------------------ZONE II------------------------------!
!Assign material types according to coordinates of center of the element!
!---------------------EQ V4.1-2016/10/02 D.Liu----------------------!
			!!For Double-couple point source. Ma and Liu (2006) 	
				! mat(nelement,1)=2800.
				! mat(nelement,2)=1500.
				! mat(nelement,3)=2600.
				! if (xc(3)<-1000.)then
					! mat(nelement,1)=6000.
					! mat(nelement,2)=3464.
					! mat(nelement,3)=2700.
				! endif
				! mat(nelement,5)=mat(nelement,2)**2*mat(nelement,3)!miu=vs**2*rho
				! mat(nelement,4)=mat(nelement,1)**2*mat(nelement,3)-2*mat(nelement,5)!lam=vp**2*rho-2*miu
!TPV8
!Feb.19.2016				
				mat(nelement,1)=6000.
				mat(nelement,2)=3464.
				mat(nelement,3)=2670.				
				mat(nelement,5)=mat(nelement,2)**2*mat(nelement,3)!miu=vs**2*rho
				mat(nelement,4)=mat(nelement,1)**2*mat(nelement,3)-2*mat(nelement,5)!lam=vp**2*rho-2*miu				
!For Tianjin 
!Sep.21.2015
				! if (xc(3)>-10e3) then
					! ixe=ix-1
					! iye=iy-1
					! ize=(xc(3)-(-10e3+dx/2))/dx+1
					! iye1=(iye-1)/10
					! itemp=iye-iye1*10
					! mat(nelement,1) =vpstruct((ixe-1)*50*105+(ize-1)*105+iye1+1,itemp)!Vp
				! elseif (xc(3)<-10e3.and.xc(3)>-20e3) then
					! mat(nelement,1) =6.42-abs(xc(3)--20e3)/10e3*(6.42-6.2)
				! elseif (xc(3)<-20e3) then	
					! mat(nelement,1) =6.42
				! endif					
				! mat(nelement,2)=0.7858-1.2344*mat(nelement,1)+0.7949*mat(nelement,1)**2&
					! -0.1238*mat(nelement,1)**3+0.0064*mat(nelement,1)**4
				! mat(nelement,3)=1.6612*mat(nelement,1)-0.4721*mat(nelement,1)**2+0.0671*&
					! mat(nelement,1)**3-0.0043*mat(nelement,1)**4+0.000106*mat(nelement,1)**5
				! if (mat(nelement,2)<0.5) then 
					! mat(nelement,2)=0.5
				! endif
				! mat(nelement,1)=mat(nelement,1)*1.0e3
				! mat(nelement,2)=mat(nelement,2)*1.0e3
				! mat(nelement,3)=mat(nelement,3)*1.0e3					
				! !mat(nelement,2)=mat(nelement,1)/sqrt(3.)!For Possion solid.
				! mat(nelement,5)=mat(nelement,2)**2*mat(nelement,3)!miu=vs**2*rho
				! mat(nelement,4)=mat(nelement,1)**2*mat(nelement,3)-2*mat(nelement,5)!lam=vp**2*rho-2*miu
!----------------------------END ZONE II----------------------------!
				!special treatment for using master nodes above the main fault.
				! B.D. 1/7/12
				if((xcoor>(fltxyz(1,1,1)-tol).and.xcoor<(fltxyz(2,1,1)+dx+tol).and. &
				zcoor>(fltxyz(1,3,1)-tol).and.ycoor>0.and.abs(ycoor-dy)<tol)) then
					do i=1,nftnd(1)
						do k=1,8
							if(ien(k,nelement)==nsmp(1,i,1)) then
								ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
								miuonf(i)=mat(nelement,5)
								vponf(i)=mat(nelement,1)
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
!---------------------EQ V4.1-2016/10/02 D.Liu----------------------!	
				if (C_elastic==0) then				
					tmp1 = -0.5*(zline(iz)+zline(iz-1))  !z<0, thus tmp1>0
					tmp2 = tmp1 * grav
					if (tmp1<=15e3) then
						omega=1
					elseif (tmp1>=15e3.and.tmp1<=20e3)then 
						omega=(20e3-tmp1)/5e3
					else
						omega=0
					endif
					
					if(et(nelement)==1)then
						eleporep(nelement)=rhow*tmp2  !pore pressure>0
						s1(ids(nelement)+3)=-mat(nelement,3)*tmp2  !vertical, comp<0
						s1(ids(nelement)+1)=omega*(b11*(s1(ids(nelement)+3)+eleporep(nelement))-eleporep(nelement))+&
							(1-omega)*s1(ids(nelement)+3)
						s1(ids(nelement)+2)=omega*(b33*(s1(ids(nelement)+3)+eleporep(nelement))-eleporep(nelement))+&
							(1-omega)*s1(ids(nelement)+3)
						s1(ids(nelement)+6)=omega*b13*(s1(ids(nelement)+3)+eleporep(nelement))
					elseif(et(nelement)==2)then
						eleporep(nelement)=rhow*tmp2  !pore pressure>0
						s1(ids(nelement)+3+15)=-mat(nelement,3)*tmp2  !vertical, comp<0
						s1(ids(nelement)+1+15)=omega*(b11*(s1(ids(nelement)+3+15)+eleporep(nelement))-eleporep(nelement))+&
							(1-omega)*s1(ids(nelement)+3+15)
						s1(ids(nelement)+2+15)=omega*(b33*(s1(ids(nelement)+3+15)+eleporep(nelement))-eleporep(nelement))+&
							(1-omega)*s1(ids(nelement)+3+15)
						s1(ids(nelement)+6+15)=omega*b13*(s1(ids(nelement)+3+15)+eleporep(nelement))							
					endif
					if (me==30.and.nnode==604363) then 
						write(*,*) s1(ids(nelement)+1),s1(ids(nelement)+2),s1(ids(nelement)+3),s1(ids(nelement)+4),s1(ids(nelement)+5),s1(ids(nelement)+6)
						write(*,*) tmp1,tmp2,grav,omega,zline(iz),zline(iz-1),et(nelement)
					endif
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
if (maxs>=(4*maxm)) then
	write(*,*) '4*maxm',maxm,'is not enough for maxs',maxs
	stop 2002
endif
if(nnode/=nx*ny*nz.or.msnode/=nnode1.or.nelement/=nelement1.or.neq/=neq1.or.nsurnd/=nsurnd1) then
	write(*,*) 'Inconsistency in node/element/equation/nsurnd between meshgen and mesh4num: stop!',me
	write(*,*) 'nnode&nnode1=',nnode,nnode1
	write(*,*) 'nelement,nelement1=',nelement,nelement1
	write(*,*) 'neq,neq1=',neq,neq1
	write(*,*) 'nsurnd,nsurnd1=',nsurnd,nsurnd1
	write(*,*) 'nnode,nx,ny,nz',nnode,nx,ny,nz
	write(*,*) 'msnode,nnode1',msnode,nnode1
	stop 2003
endif
if(ntag/=maxm) then
	write(*,*) 'Inconsistency in ntag and maxm: stop!',me
	write(*,*) ntag,maxm
	stop 2004
endif
do i=1,ntotft
	if(nftnd(i)/=nftnd1(i)) then
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
	endif!end if nftnd=0/ not
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
 	do i = 1, nftnd(ift)
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

! !...write out mesh information, including output node coor for check. B.D. 10/27/07
! meshfile = 'meshinfo.txt'//mm 
! open(33,file=meshfile,status='unknown')
! write(33,*) 'me=',me
! write(33,*) 'mex=',mex,'mey=',mey,'mez=',mez
! write(33,*) 'nx=',nx,'ny=',ny,'nz=',nz
! write(33,*) 'x-range for this MPI:',xline(1),xline(nx)
! write(33,*) 'y-range for this MPI:',yline(1),yline(ny)
! write(33,*) 'z-range for this MPI:',zline(1),zline(nz)
! write(33,*) 'numcount1-3',numcount(1),numcount(2),numcount(3)
! write(33,*) 'numcount4-6',numcount(4),numcount(5),numcount(6)
! write(33,*) 'numcount7-9',numcount(7),numcount(8),numcount(9)
! write(33,*) 'fltnum1-3',fltnum(1),fltnum(2),fltnum(3)
! write(33,*) 'fltnum4-6',fltnum(4),fltnum(5),fltnum(6)
! write(33,'( a45,5i10)') 'total in this process:node,ele,neq,nftnd ',nnode,nelement,neq,(nftnd(i),i=1,ntotft)
! write(33,'( a60)') 'off-fault stations # in order and x, y, z- xoor'
! write(33,'( 2i10,2f12.2)') ((an4nds(j,i),j=1,2),(x4nds(j,an4nds(1,i)), &
  	! j=1,2),i=1,n4out)
! write(33,'( a65,i5)') 'on-fault stations # in order and strike and down-dip distance',n4onf
! write(33,'( 3i10,2f12.2)') ((anonfs(j,i),j=1,3),(xonfs(j,anonfs(2,i),anonfs(3,i)),j=1,2), &
	! i=1,n4onf)
! !write(33,*) 'branch fault'
! !write(33,'( 10f8.1)') (arn(i,2),i=1,nftnd(2))
! close(33)

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
