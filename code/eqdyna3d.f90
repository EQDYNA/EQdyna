!========================================================================!
! Complete partitioning for better performance of parallel EQdyna.
! This is collaborative work by Drs. Benchun Duan and Xingfu Wu 
!     at Texas A&M Univ. This is version 2.3.2 of EQdyna 3D.
! April, 2009.
!========================================================================!
! Hybrid MPI/OpenMP Implementation of EQDY3D
! Author: Xingfu Wu, CSE Department, Texas A&M University
! Feb. 5. 2009
!========================================================================!
!			*** EQDY3D ***				 !
! A program to simulate earthquke dynamic problems in 3D using Finite 	 !
! Element Method. This is the main program unit and global driver.	 !
!	Author: Benchun Duan (B.D.)					 !
!	Version: This is version 1.0. Start from 8/16/05 B.D.		 !
!									 !
!========================================================================!
			
!========================================================================!
!		*** Version 4.0 featuring 3D MPI partitioning***
!		*** Feb.19.2016 ***
!		*** Written by B. Luo; Implemented by D. Liu ***
!		*** Search Tag: 3DMPI ***
			
!========================================================================!
PROGRAM eqdy3d
use globalvar
implicit none
include 'mpif.h'
!*** declaration ***
real (kind=8) :: timebegin, timeover
character (len=7) outpath
character (len=30) :: outgl,outft,outsl,outft1, outft2, outft3, outoff,outgm,&
	outrat,outdp
character (len=4), dimension(20) :: title	!one line of input file
integer (kind=4) :: n,i,j,k,l,itemp1,itemp2,itemp3,alloc_err
real (kind=8) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9, &
	temp10,temp11,temp12,temp13,temp14,temp15
!...important integers  		
integer (kind=4) :: numnp,numel,neq,n4out,ndout=0,n4onf,nftmx,nonmx
integer (kind=4),allocatable,dimension(:) :: nftnd
!...model and fault dimension
real (kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax
real (kind=8),allocatable,dimension(:,:,:) :: fltxyz !4: x,y,z,strike/dip
!......time histories
integer (kind=4),dimension(2,187)::an4nds
integer (kind=4),allocatable,dimension(:,:) :: idhist
real (kind=4),allocatable,dimension(:,:) :: dout
!......nodes' arrays
!integer (kind=4),allocatable,dimension(:,:) :: id
real (kind=8),allocatable,dimension(:,:) :: x
real (kind=8),allocatable,dimension(:) :: brhs
real (kind=8),allocatable,dimension(:,:) :: d,v
!...element arrays
real (kind=8),allocatable,dimension(:,:) :: mat	    
integer (kind=4),allocatable,dimension(:,:) :: ien
integer (kind=4),allocatable,dimension(:,:,:) :: lm
real (kind=8),allocatable,dimension(:) :: eleporep,pstrain
real (kind=8),allocatable,dimension(:,:) :: shl!Delete elestress
!...fault arrays
integer (kind=4),allocatable,dimension(:,:) :: anonfs
integer (kind=4),allocatable,dimension(:,:,:) :: nsmp
real (kind=8),allocatable,dimension(:,:) :: fnft,arn,r4nuc,arn4m,slp4fri
real (kind=8),allocatable,dimension(:,:,:) :: fric,un,us,ud,fltsta,fltslp
!...material properties
real (kind=8),allocatable,dimension(:) :: emod,pois,mushr,rho,rdampm,rdampk,&
ccosphi,sinphi,vp,vs
real (kind=8),allocatable,dimension(:,:,:) :: c
!...
integer (kind=4),dimension(2,2)::nr4mpi,er4mpi  !node,equation range for MPI
!...moment, moment rate, max slip rate at the end. B.D. 8/11/10
real (kind=8)::momntall,momntratall,maxslpratall,magn
!...on-fault station info. B.D. 8/11/10
!  integer (kind=4),dimension(2) :: nonfs=(/8,6/)
!  real (kind=8),dimension(2,8,2) :: xonfs  !8 is max of nonfs
integer (kind=4),dimension(1) :: nonfs=(/8/)
real (kind=8),dimension(2,8,1) :: xonfs
character (len=90) :: loca
integer time_array(8)
!...working variables. B.D. 1/15/12
integer (kind=4) :: itmp1,itmp2
real (kind=8) :: tmp1,tmp2
logical:: sts=.false.,dpyn=.false.
!*.* Variables for new features in PML. D.L. Jan/23/15
integer(kind=4)::maxm,maxs,nsurnd,nnodPML
real(kind=8),dimension(8)::PMLb
real (kind=8),allocatable,dimension(:) :: v1,d1,s1 ! 1D array for velocity,displacement and stress 
integer (kind=4),allocatable,dimension(:) :: id1,ids,et,locid,dof1 !et: element type.
!
real(kind=8),allocatable::surnode(:,:),miuonf(:),vponf(:)!
integer(kind=4),allocatable::surid(:),nPMLid(:),nPMLvp(:)
integer(kind=4)::isur
character(len=30)::fsur
!*.* D.L.
!$     integer  omp_get_max_threads
!$     external omp_get_max_threads
!!$     integer  omp_get_num_procs
!!$     external omp_get_num_procs
integer master, me, nprocs
integer istatus(MPI_STATUS_SIZE)
real (kind=8),allocatable,dimension(:) :: aa
!integer (kind=4),allocatable,dimension(:) :: kk,jj,ll

!new variables for MPI. B.D. 4/15/09	    
integer (kind=4) :: rlp,rr,ierr,jj1,i1,i2,i3
character (len=6) :: mm
!-------------------------------------------------------------------!
!-------------Global control of the EQ V3.2.1 system----------------!
!-----------------------Sep.19.2015/ D.Liu--------------------------!
if (C_elastic==0.and.C_Q==1) then
	write(*,*) 'Q model can only work with elastic code'
	stop 1001 
endif
if (C_Q==1.and.rat>1) then
	write(*,*) 'Q model can only work with uniform element size'
	write(*,*) 'rat should be 1.0'
	stop 1002
endif
!-----------End of Global control of the EQ V3.2.1 system-----------!
!-------------------------------------------------------------------!
! MPI initialize
call MPI_Init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,me,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
master = 0
!
!...because eos and vanerne with linux do not allow 'write(mm,*) me',
! I have to use a way for conversion from num to char as below for
! what I want. Currently, only for MPI processes <= 300. For more, need
! to extend further.  B.D. 8/13/10
!...actually, using formatted write to solve the problem.
!  write(mm,*) me
write(mm,'(i6)') me
!  if(me<10) then
!    mm = char(me+48)
!  elseif(me<20) then
!    mm = char(49)//char(me-10+48)
!  ...
!  elseif(me<300) then
!    mm = char(50)//char(57)//char(me-190+48)
!  endif
mm = trim(adjustl(mm))

!*** record starting time ***
timebegin = MPI_WTIME()
time1 = MPI_WTIME()		
!
!*** set up I/O ***
!
!......directly give file name
outgl  = 'fem.txt'//mm
!  outft = 'fltst.txt'//mm
outft1 = 'frt.txt'//mm
outft2 = 'frt2.txt'//mm
outdp = 'srfdp.txt'//mm
!  outoff = 'offst.txt'//mm
!  outsl = 'fnlslp.txt'//mm
!  outrat = 'ftrat.txt'//mm
! outft2 = 'fst.txt'//mm
!  outft3 = 'fsl.txt'//mm
!  outgm = 'gmt.txt'

!
!*** INPUT PHASE ***
!
!...parameter control file is incorporated into the code.
!   no input file is needed. B.D. 10/28/09
allocate(emod(numat),pois(numat),mushr(numat),rho(numat),vp(numat),vs(numat),&
		rdampm(numat),rdampk(numat),ccosphi(numat),sinphi(numat),&
		nftnd(ntotft),fltxyz(2,4,ntotft))
call parcon(xmin,xmax,ymin,ymax,zmin,zmax,fltxyz,&
			emod,pois,mushr,rho,vp,vs,rdampm,rdampk,ccosphi,sinphi)

!...time hostories output steps
nplpts = 0	!initialize number of time history plot
if (nhplt > 0) then
	nplpts = int(nstep/nhplt) + 2
endif
!
!...material property c and shl 
allocate(c(nrowc,nrowc,numat),shl(nrowsh,nen))
call qdct1(shl,emod,pois,c)

!*** PRE-MESHING to GET NUMBERS for ARRAYS' DEFINITTION
call mesh4num(fltxyz,xmin,xmax,ymin,ymax,zmin,zmax, &
				numnp,numel,neq,PMLb,maxm,nftnd,master,me,nprocs,nsurnd,nnodPML)!Adding PMLb and maxm
!  write(*,*) 'me=', me, '# of fault node',(nftnd(i),i=1,ntotft),'# of elements',numel
write(*,*) 'me=',me,'maxm=',maxm
!*** ALLOCATE ARRAYS ***
!...nodes' arrays
allocate(id1(maxm),locid(numnp),dof1(numnp),x(ndof,numnp),stat=alloc_err)!Adding id1(maxm),loci(2,numnp) !Delete id(ndof,numnp),
allocate(surid(nsurnd),surnode(nsurnd,2),stat=alloc_err)
allocate(nPMLvp(nnodPML),nPMLid(nnodPML),stat=alloc_err)
surid=0
surnode=0.0 
nPMLvp=0
nPMLid=0
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate arrays id or x'
endif
!...element's arrays
allocate(ien(8,numel),mat(numel,5),et(numel),&! Delete lm(ned,nen,numel),elestress(nstr,numel),
		eleporep(numel),pstrain(numel),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate arrays ien, mat, or lm'
endif
!elestress = 0.0  !initialize
eleporep = 0.0
pstrain = 0.0

allocate(fltgm(nftmx))  !MPIxyz
fltgm = 0  !MPIxyz
!after allocate, initialize here before call meshgen
nsmp = 0    !initialize here
fnft = 1000. !as failure time initialization: should larger than actual!
miuonf=0.0!...fault arrays
nftmx = maxval(nftnd) !max fault nodel num for all faults, used for arrays.
if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
nonmx = sum(nonfs)    !max possible on-fault stations number
allocate(nsmp(2,nftmx,ntotft),fnft(nftmx,ntotft),miuonf(nftmx),vponf(nftmx),un(3,nftmx,ntotft),&
			us(3,nftmx,ntotft),ud(3,nftmx,ntotft),fric(8,nftmx,ntotft),&
			arn(nftmx,ntotft),r4nuc(nftmx,ntotft),anonfs(3,nonmx),&
			arn4m(nftmx,ntotft),slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft))
vponf=0.0
fric = 0.0
un = 0.0
us = 1000.0
ud = 0.0
arn = 0.0
arn4m = 0.0
r4nuc = 0.0
anonfs = 0
slp4fri = 0.0
fltslp = 0.0
!
!*** MESH GENERATION **** 
write(*,*) 'before meshgen me=',me    
allocate(ids(numel))
allocate(s1(4*maxm))
s1=0.0
call meshgen(fltxyz,xmin,xmax,ymin,ymax,zmin,zmax, &
			numnp,numel,ien,mat,s1,eleporep, &!change elestress as s1
			x,neq, & !Delete id
			id1,maxm,locid,dof1,et,PMLb,maxs,ids, &!Adding id1,maxm,loci,PMLb,maxs,ids
			nr4mpi,er4mpi,n4out,an4nds, & 
			nftnd,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,arn4m,master,me,nprocs,nsurnd,surnode,surid,miuonf,vponf,nnodPML,nPMLid,nPMLvp)			
write(*,*) 'after meshgen me=',me 
!-------------------------------------------------------------------!
!---------------------OutPut surface nodes coordinates--------------!
!----------------------Oct.1.2015/ D.Liu----------------------------!
! if (nsurnd>0) then
	! fsur='surnode.txt'//mm 
	! open(1001,file=fsur,form='formatted',status='unknown')
	! write(1001,'(1x,2f10.1)') (surnode(isur,1),surnode(isur,2),isur=1,nsurnd)
	! close(1001)
! endif
!-------------------------------------------------------------------!
!...allocate on-fault station array. B.D. 10/25/09
!  if(n4onf>0) then
if(n4onf<=0) n4onf=1 !cannot be zero in fortran, min 1
allocate(fltsta(10,nplpts-1,n4onf),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate array fltsta'
endif
fltsta = 0.0
!  endif  

!...form assembly mapping. 
!call formlm(numel,numnp,lm,id,ien)  

!...... allocate brhs,v,d for this proc
allocate(brhs(neq),v1(neq),d1(neq),v(ndof,numnp),d(ndof,numnp),stat=alloc_err)!Adding v1,d1
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate arrays brhs,v,or d'
endif
brhs = 0.0	!initialize (initial conditions if need)
v1 = 0.0!DL
d1 = 0.0!DL
v = 0.0
d = 0.0
locplt = 1
!...initialize nodal output time-history data
!...not use global arrays, thus if n4nout=0, set =1 for universal call.
! B.D. 1/16/12
!if(n4out > 0) then
if(n4out==0) n4out=1
ndout = n4out * ndof * noid	!3 components of 2 quantities: v and d
!   write(*,*) 'ndout= ',ndout    
allocate(idhist(3,ndout),dout(ndout+1,nplpts),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate array idhist or dout'
endif
idhist = 0
dout = 0	!initialize for safety
i1 = 0
do i=1,n4out
	do j=1,ndof
		do k=1,noid
			i1 = i1 + 1
			idhist(1,i1) = an4nds(2,i) !node number (>1, <NUMNP)
			if(idhist(1,i1)<=0) idhist(1,i1)=1  !avoid zero that cannot be in array below
				idhist(2,i1) = j	!degree of freedom number (<=NDOF)
				idhist(3,i1) = k	!kinematic quantity specifier 
				!(disp, vel, or acc)
		enddo
	enddo
enddo
!...store initial kinematic data for time histories
dout(1,locplt) = time
do i=1,ndout
	itemp1 = idhist(1,i)	!node number
	itemp2 = idhist(2,i)	!i dof
	itemp3 = idhist(3,i)	!code of value
	if(itemp3 == 1) then
		dout(i+1,locplt) = d(itemp2,itemp1)
	elseif(itemp3 == 2) then
		dout(i+1,locplt) = v(itemp2,itemp1)
	elseif(itemp3 == 3) then
		k = id1(locid(itemp1)+itemp2) !id(itemp2,itemp1) !*.* D.L. Jan/23/15
		dout(i+1,locplt) = brhs(k)
	endif
enddo
!  endif      				

!......timing CPU elapsed for initialize and mesh generation
time2 = MPI_WTIME()		
timeused(1) = time2 - time1 	!time in second

! wait for all MPI processes to get initial data
!call mpi_barrier(MPI_COMM_WORLD, ierr)
!......open files for writing
! if (me == master) then
!  open(unit=ioutst,file=outft2,status='unknown')	!stress
!  open(unit=ioutsl,file=outft3,status='unknown')	!slip velocity
!  open(unit=ioutgm,file=outgm,status='unknown')		!ground motion
! else
!  open(unit=ioutst,file=outft2,position='append')	!stress
!  open(unit=ioutsl,file=outft3,position='append')	!slip velocity
!  open(unit=ioutgm,file=outgm,position='append')		!ground motion
! endif

!...for fault slip rate time histories output. B.D. 8/11/10
!  if(nftnd>0) then
!    open(unit=ioutrat,file=outrat,form='unformatted', &
!      access='direct',recl=4,status='replace')
!   open(unit=ioutst,file=outft2,status='replace')
!  endif

!
!*** SOLUTION PHASE ***
!
if (iexec == 1) then	!otherwise data check only
	!write(*,*) 'before driver me=', me,'ndout=',ndout
	call driver(numel,numnp,neq,nftnd,ndout,dout,idhist, &
				x,brhs,d,v,mat,ien,eleporep,pstrain, & ! Delete id elestress
				id1,maxm,locid,dof1,et,v1,d1,PMLb,maxs,ids,s1, & !Adding id1(:),maxm,loci(:,:),v1,d1,PMLb
				shl, & ! Delete lm
				emod,pois,mushr,rdampm,rdampk,c,ccosphi,sinphi,&
				er4mpi,nr4mpi,n4onf,momntall,momntratall,maxslpratall,&
				nftmx,nonmx,nsmp,fnft,fltslp,un,us,ud,fric,arn,r4nuc,&
				anonfs,arn4m,slp4fri,fltsta,master,me,nprocs,nsurnd,surid,miuonf,&
				nnodPML,nPMLid,nPMLvp)
endif
write(*,*) 'after driver me=',me    
!*** record ending time ***
timeover = MPI_WTIME()
! if (me == master) then
!  close(ioutst)
!  close(ioutsl)
!  close(ioutgm)
! endif
!*** total execution time ***
timeused(7) = timeover - timebegin 
!  write(*,*) 'total seconds used:', timeused(7)
! output the number of threads used
!  close(ioutrat)
!
!*** OUTPUT PHASE ***
!
!...summary text file. B.D. 8/11/10
! calculate magnitude from moment first.
magn =( 2./3.)*log10(momntall*1.e7) - 10.7
open(51,file='event_summary.txt',status='unknown')
write(51,'( i10)') 321
write(51,'( f10.1)') dx
write(51,'( i10)') 151
write(51,'( f10.1)') dx
!write(51,'( i10)') int(nstep/nhplt1)+1
write(51,'( i10)') int(locplt/nhplt1)
write(51,'( f10.5)') dt*2
write(51,'( e15.6)') momntall
write(51,'( f10.3)') magn
close(51)
!...on-fault stations. B.D. 8/11/10
if(n4onf > 0) then
	do i=1,n4onf
		j=anonfs(3,i)
		if(j==1)  then  !main fault stations
			if(xonfs(1,anonfs(2,i),j)== 0 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst000dp000.txt',status='unknown')
				loca = '# location=on main fault, 0 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)==4.5e3 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst045dp000.txt',status='unknown')
				loca = '# location = on main fault, 4.5 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== 12e3 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst120dp000.txt',status='unknown')
				loca = '# location = on mian fault, 12.0 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== 0 .and. xonfs(2,anonfs(2,i),j)==-4.5e3) then
				open(51,file='faultst000dp045.txt',status='unknown')
				loca = '# location = on mian fault, 0 km along strike, 4.5 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)==0 .and. xonfs(2,anonfs(2,i),j)==-7.5e3) then
				open(51,file='faultst000dp075.txt',status='unknown')
				loca = '# location = on main fault, 0.0 km along strike, 7.5 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== 4.5e3 .and. xonfs(2,anonfs(2,i),j)==-7.5e3) then
				open(51,file='faultst045dp075.txt',status='unknown')
				loca = '# location = on main fault, 4.5 km along strike, 7.5 km down-dip'				
			elseif(xonfs(1,anonfs(2,i),j)== 12e3 .and. xonfs(2,anonfs(2,i),j)==-7.5e3) then
				open(51,file='faultst120dp075.txt',status='unknown')
				loca = '# location = on main fault, 12.0 km along strike, 7.5 km down-dip'				
			elseif(xonfs(1,anonfs(2,i),j)== 0 .and. xonfs(2,anonfs(2,i),j)==-12e3) then
				open(51,file='faultst000dp120.txt',status='unknown')				
				loca = '# location = on mian fault, 0 km along strike, 12 km down-dip'
			endif
		endif
		write(51,*) '#For SCEC TPV8'
		write(51,*) '#Author=D.Liu'
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
			'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
			':',time_array(7)
		write(51,*) '#code = EQdyna3d'
		write(51,*) '#code_version=4.0'
		write(51,*) '#element_size =',dx
		write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
		write(51,'( a19,i6)') '# num_time_steps =', locplt-1
		!write(51,*) loca !Disable the loca to avoid writting in two lines.
		write(51,*) '# Time series in 8 columns in format e15.7'
		write(51,*) '# Column #1 = Time (s)'
		write(51,*) '# Column #2 = horizontal slip (m)'
		write(51,*) '# Column #3 = horizontal slip rate (m/s)'
		write(51,*) '# Column #4 = horizontal shear stress (MPa)'
		write(51,*) '# Column #5 = down-dip slip (m)'
		write(51,*) '# Column #6 = down-dip slip rate (m/s)'
		write(51,*) '# Column #7 = down-dip shear stress (MPa)'
		!write(51,*) '# Column #8 = normal stress (MPa)'
		write(51,*) '# The line below lists the names of the data fields:'
		write(51,*) 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress'
		write(51,*) '#'
		do j=1,locplt-1
			write(51,'( f12.5,6e18.7e4)') fltsta(1,j,i),fltsta(5,j,i),fltsta(2,j,i),fltsta(8,j,i)/1.e6,&
					-fltsta(6,j,i),-fltsta(3,j,i),-fltsta(9,j,i)/1.e6
		enddo
		close(51)
	enddo
endif 
!...output off-fault stations. B.D. 8/11/10
if(n4out > 0) then
	do i=1,n4out
		! if(an4nds(1,i)==1) then 
			! open(51,file='body-1BeiJing.txt',status='unknown')
		! elseif(an4nds(1,i)==2) then
			! open(51,file='body-2TianJin.txt',status='unknown')
		! elseif(an4nds(1,i)==3) then
			! open(51,file='body-3LangFang.txt',status='unknown')
		! elseif(an4nds(1,i)==4) then
			! open(51,file='body-4NanKai.txt',status='unknown')
		! elseif(an4nds(1,i)==5) then
			! open(51,file='body-5TJRS.txt',status='unknown')
		! elseif(an4nds(1,i)==6) then
			! open(51,file='body-6YongQing.txt',status='unknown')
		! elseif(an4nds(1,i)==7) then
			! open(51,file='body-7TangShan.txt',status='unknown')
		! elseif(an4nds(1,i)==8) then
			! open(51,file='body-8NingHe.txt',status='unknown')
		! elseif(an4nds(1,i)==9) then
			! open(51,file='body-9BaoDi.txt',status='unknown')
		! elseif(an4nds(1,i)==10) then
			! open(51,file='body-10XiangHe.txt',status='unknown')
		! elseif(an4nds(1,i)==11) then 
			! open(51,file='body-11JingHai.txt',status='unknown')
		! elseif(an4nds(1,i)==12) then
			! open(51,file='body-12DaCheng.txt',status='unknown')
		! elseif(an4nds(1,i)==13) then
			! open(51,file='body-13QingXian.txt',status='unknown')
		! elseif(an4nds(1,i)==14) then
			! open(51,file='body-14WenAn.txt',status='unknown')
		! elseif(an4nds(1,i)==15) then
			! open(51,file='body-15WuQing.txt',status='unknown')
		! elseif(an4nds(1,i)==16) then
			! open(51,file='body-16RenQiu.txt',status='unknown')
		! elseif(an4nds(1,i)==17) then
			! open(51,file='body-17XiongXian.txt',status='unknown')	
		! elseif(an4nds(1,i)==18) then
			! open(51,file='body-18BinHai.txt',status='unknown')
		! elseif(an4nds(1,i)==19) then!19~54-A1
			! open(51,file='body-A1.txt',status='unknown')
		! elseif(an4nds(1,i)==20) then!A2
			! open(51,file='body-A2.txt',status='unknown')
		! elseif(an4nds(1,i)==21) then!A3
			! open(51,file='body-A3.txt',status='unknown')
		! elseif(an4nds(1,i)==22) then!A4
			! open(51,file='body-A4.txt',status='unknown')
		! elseif(an4nds(1,i)==23) then!A5
			! open(51,file='body-A5.txt',status='unknown')
		! elseif(an4nds(1,i)==24) then!A6
			! open(51,file='body-A6.txt',status='unknown')	
		! elseif(an4nds(1,i)==25) then!A7
			! open(51,file='body-A7.txt',status='unknown')
		! elseif(an4nds(1,i)==26) then!A8
			! open(51,file='body-A8.txt',status='unknown')
		! elseif(an4nds(1,i)==27) then!A9
			! open(51,file='body-A9.txt',status='unknown')
		! elseif(an4nds(1,i)==28) then!A10
			! open(51,file='body-A10.txt',status='unknown')
		! elseif(an4nds(1,i)==29) then!A11
			! open(51,file='body-A11.txt',status='unknown')
		! elseif(an4nds(1,i)==30) then!A12
			! open(51,file='body-A12.txt',status='unknown')
		! elseif(an4nds(1,i)==31) then!A13
			! open(51,file='body-A13.txt',status='unknown')
		! elseif(an4nds(1,i)==32) then!A14
			! open(51,file='body-A14.txt',status='unknown')
		! elseif(an4nds(1,i)==33) then!A15
			! open(51,file='body-A15.txt',status='unknown')
		! elseif(an4nds(1,i)==34) then!A16
			! open(51,file='body-A16.txt',status='unknown')
		! elseif(an4nds(1,i)==35) then!A17
			! open(51,file='body-A17.txt',status='unknown')
		! elseif(an4nds(1,i)==36) then!A18
			! open(51,file='body-A18.txt',status='unknown')
		! elseif(an4nds(1,i)==37) then!A19
			! open(51,file='body-A19.txt',status='unknown')
		! elseif(an4nds(1,i)==38) then!A20
			! open(51,file='body-A20.txt',status='unknown')
		! elseif(an4nds(1,i)==39) then!A21
			! open(51,file='body-A21.txt',status='unknown')
		! elseif(an4nds(1,i)==40) then!A22
			! open(51,file='body-A22.txt',status='unknown')
		! elseif(an4nds(1,i)==41) then!A23
			! open(51,file='body-A23.txt',status='unknown')
		! elseif(an4nds(1,i)==42) then!A24
			! open(51,file='body-A24.txt',status='unknown')
		! elseif(an4nds(1,i)==43) then!A25
			! open(51,file='body-A25.txt',status='unknown')
		! elseif(an4nds(1,i)==44) then!A26
			! open(51,file='body-A26.txt',status='unknown')
		! elseif(an4nds(1,i)==45) then!A27
			! open(51,file='body-A27.txt',status='unknown')
		! elseif(an4nds(1,i)==46) then!A28
			! open(51,file='body-A28.txt',status='unknown')
		! elseif(an4nds(1,i)==47) then!A29
			! open(51,file='body-A29.txt',status='unknown')
		! elseif(an4nds(1,i)==48) then!A30
			! open(51,file='body-A30.txt',status='unknown')
		! elseif(an4nds(1,i)==49) then!A31
			! open(51,file='body-A31.txt',status='unknown')
		! elseif(an4nds(1,i)==50) then!A32
			! open(51,file='body-A32.txt',status='unknown')
		! elseif(an4nds(1,i)==51) then!A33
			! open(51,file='body-A33.txt',status='unknown')
		! elseif(an4nds(1,i)==52) then!A34
			! open(51,file='body-A34.txt',status='unknown')
		! elseif(an4nds(1,i)==53) then!A35
			! open(51,file='body-A35.txt',status='unknown')
		! elseif(an4nds(1,i)==54) then!A36
			! open(51,file='body-A36.txt',status='unknown')
! !B		
		! elseif(an4nds(1,i)==55) then!55~80-B1
			! open(51,file='body-B1.txt',status='unknown')
		! elseif(an4nds(1,i)==56) then!B2
			! open(51,file='body-B2.txt',status='unknown')
		! elseif(an4nds(1,i)==57) then!B3
			! open(51,file='body-B3.txt',status='unknown')
		! elseif(an4nds(1,i)==58) then!B4
			! open(51,file='body-B4.txt',status='unknown')
		! elseif(an4nds(1,i)==59) then!B5
			! open(51,file='body-B5.txt',status='unknown')
		! elseif(an4nds(1,i)==60) then!B6
			! open(51,file='body-B6.txt',status='unknown')	
		! elseif(an4nds(1,i)==61) then!B7
			! open(51,file='body-B7.txt',status='unknown')
		! elseif(an4nds(1,i)==62) then!B8
			! open(51,file='body-B8.txt',status='unknown')
		! elseif(an4nds(1,i)==63) then!B9
			! open(51,file='body-B9.txt',status='unknown')
		! elseif(an4nds(1,i)==64) then!B10
			! open(51,file='body-B10.txt',status='unknown')
		! elseif(an4nds(1,i)==65) then!B11
			! open(51,file='body-B11.txt',status='unknown')
		! elseif(an4nds(1,i)==66) then!B12
			! open(51,file='body-B12.txt',status='unknown')
		! elseif(an4nds(1,i)==67) then!B13
			! open(51,file='body-B13.txt',status='unknown')
		! elseif(an4nds(1,i)==68) then!B14
			! open(51,file='body-B14.txt',status='unknown')
		! elseif(an4nds(1,i)==69) then!B15
			! open(51,file='body-B15.txt',status='unknown')
		! elseif(an4nds(1,i)==70) then!B16
			! open(51,file='body-B16.txt',status='unknown')
		! elseif(an4nds(1,i)==71) then!B17
			! open(51,file='body-B17.txt',status='unknown')
		! elseif(an4nds(1,i)==72) then!B18
			! open(51,file='body-B18.txt',status='unknown')
		! elseif(an4nds(1,i)==73) then!B19
			! open(51,file='body-B19.txt',status='unknown')
		! elseif(an4nds(1,i)==74) then!B20
			! open(51,file='body-B20.txt',status='unknown')
		! elseif(an4nds(1,i)==75) then!B21
			! open(51,file='body-B21.txt',status='unknown')
		! elseif(an4nds(1,i)==76) then!B22
			! open(51,file='body-B22.txt',status='unknown')
		! elseif(an4nds(1,i)==77) then!B23
			! open(51,file='body-B23.txt',status='unknown')
		! elseif(an4nds(1,i)==78) then!B24
			! open(51,file='body-B24.txt',status='unknown')
		! elseif(an4nds(1,i)==79) then!B25
			! open(51,file='body-B25.txt',status='unknown')
		! elseif(an4nds(1,i)==80) then!B26
			! open(51,file='body-B26.txt',status='unknown')
		! elseif(an4nds(1,i)==81) then!B27
			! open(51,file='body-B27.txt',status='unknown')
		! elseif(an4nds(1,i)==82) then!B28
			! open(51,file='body-B28.txt',status='unknown')
		! elseif(an4nds(1,i)==83) then!B29
			! open(51,file='body-B29.txt',status='unknown')
		! elseif(an4nds(1,i)==84) then!B30
			! open(51,file='body-B30.txt',status='unknown')
		! elseif(an4nds(1,i)==85) then!B31
			! open(51,file='body-B31.txt',status='unknown')
		! elseif(an4nds(1,i)==86) then!B32
			! open(51,file='body-B32.txt',status='unknown')
		! elseif(an4nds(1,i)==87) then!B33
			! open(51,file='body-B33.txt',status='unknown')
		! elseif(an4nds(1,i)==88) then!B34
			! open(51,file='body-B34.txt',status='unknown')
		! elseif(an4nds(1,i)==89) then!B35
			! open(51,file='body-B35.txt',status='unknown')
		! elseif(an4nds(1,i)==90) then!B36
			! open(51,file='body-B36.txt',status='unknown')	
! !C
		! elseif(an4nds(1,i)==91) then!91~126-C1
			! open(51,file='body-C1.txt',status='unknown')
		! elseif(an4nds(1,i)==92) then!C2
			! open(51,file='body-C2.txt',status='unknown')
		! elseif(an4nds(1,i)==93) then!C3
			! open(51,file='body-C3.txt',status='unknown')
		! elseif(an4nds(1,i)==94) then!C4
			! open(51,file='body-C4.txt',status='unknown')
		! elseif(an4nds(1,i)==95) then!C5
			! open(51,file='body-C5.txt',status='unknown')
		! elseif(an4nds(1,i)==96) then!C6
			! open(51,file='body-C6.txt',status='unknown')	
		! elseif(an4nds(1,i)==97) then!C7
			! open(51,file='body-C7.txt',status='unknown')
		! elseif(an4nds(1,i)==98) then!C8
			! open(51,file='body-C8.txt',status='unknown')
		! elseif(an4nds(1,i)==99) then!C9
			! open(51,file='body-C9.txt',status='unknown')
		! elseif(an4nds(1,i)==100) then!C10
			! open(51,file='body-C10.txt',status='unknown')
		! elseif(an4nds(1,i)==101) then!C11
			! open(51,file='body-C11.txt',status='unknown')
		! elseif(an4nds(1,i)==102) then!C12
			! open(51,file='body-C12.txt',status='unknown')
		! elseif(an4nds(1,i)==103) then!C13
			! open(51,file='body-C13.txt',status='unknown')
		! elseif(an4nds(1,i)==104) then!C14
			! open(51,file='body-C14.txt',status='unknown')
		! elseif(an4nds(1,i)==105) then!C15
			! open(51,file='body-C15.txt',status='unknown')
		! elseif(an4nds(1,i)==106) then!C16
			! open(51,file='body-C16.txt',status='unknown')
		! elseif(an4nds(1,i)==107) then!C17
			! open(51,file='body-C17.txt',status='unknown')
		! elseif(an4nds(1,i)==108) then!C18
			! open(51,file='body-C18.txt',status='unknown')
		! elseif(an4nds(1,i)==109) then!C19
			! open(51,file='body-C19.txt',status='unknown')
		! elseif(an4nds(1,i)==110) then!C20
			! open(51,file='body-C20.txt',status='unknown')
		! elseif(an4nds(1,i)==111) then!C21
			! open(51,file='body-C21.txt',status='unknown')
		! elseif(an4nds(1,i)==112) then!C22
			! open(51,file='body-C22.txt',status='unknown')
		! elseif(an4nds(1,i)==113) then!C23
			! open(51,file='body-C23.txt',status='unknown')
		! elseif(an4nds(1,i)==114) then!C24
			! open(51,file='body-C24.txt',status='unknown')
		! elseif(an4nds(1,i)==115) then!C25
			! open(51,file='body-C25.txt',status='unknown')
		! elseif(an4nds(1,i)==116) then!C26
			! open(51,file='body-C26.txt',status='unknown')
		! elseif(an4nds(1,i)==117) then!C27
			! open(51,file='body-C27.txt',status='unknown')
		! elseif(an4nds(1,i)==118) then!C28
			! open(51,file='body-C28.txt',status='unknown')
		! elseif(an4nds(1,i)==119) then!C29
			! open(51,file='body-C29.txt',status='unknown')
		! elseif(an4nds(1,i)==120) then!C30
			! open(51,file='body-C30.txt',status='unknown')
		! elseif(an4nds(1,i)==121) then!C31
			! open(51,file='body-C31.txt',status='unknown')
		! elseif(an4nds(1,i)==122) then!C32
			! open(51,file='body-C32.txt',status='unknown')
		! elseif(an4nds(1,i)==123) then!C33
			! open(51,file='body-C33.txt',status='unknown')
		! elseif(an4nds(1,i)==124) then!C34
			! open(51,file='body-C34.txt',status='unknown')
		! elseif(an4nds(1,i)==125) then!C35
			! open(51,file='body-C35.txt',status='unknown')
		! elseif(an4nds(1,i)==126) then!C36
			! open(51,file='body-C36.txt',status='unknown')
! !D
		! elseif(an4nds(1,i)==127) then!126~161-D1
			! open(51,file='body-D1.txt',status='unknown')
		! elseif(an4nds(1,i)==128) then!D2
			! open(51,file='body-D2.txt',status='unknown')
		! elseif(an4nds(1,i)==129) then!D3
			! open(51,file='body-D3.txt',status='unknown')
		! elseif(an4nds(1,i)==130) then!D4
			! open(51,file='body-D4.txt',status='unknown')
		! elseif(an4nds(1,i)==131) then!D5
			! open(51,file='body-D5.txt',status='unknown')
		! elseif(an4nds(1,i)==132) then!D6
			! open(51,file='body-D6.txt',status='unknown')	
		! elseif(an4nds(1,i)==133) then!D7
			! open(51,file='body-D7.txt',status='unknown')
		! elseif(an4nds(1,i)==134) then!D8
			! open(51,file='body-D8.txt',status='unknown')
		! elseif(an4nds(1,i)==135) then!D9
			! open(51,file='body-D9.txt',status='unknown')
		! elseif(an4nds(1,i)==136) then!D10
			! open(51,file='body-D10.txt',status='unknown')
		! elseif(an4nds(1,i)==137) then!D11
			! open(51,file='body-D11.txt',status='unknown')
		! elseif(an4nds(1,i)==138) then!D12
			! open(51,file='body-D12.txt',status='unknown')
		! elseif(an4nds(1,i)==139) then!D13
			! open(51,file='body-D13.txt',status='unknown')
		! elseif(an4nds(1,i)==140) then!D14
			! open(51,file='body-D14.txt',status='unknown')
		! elseif(an4nds(1,i)==141) then!D15
			! open(51,file='body-D15.txt',status='unknown')
		! elseif(an4nds(1,i)==142) then!D16
			! open(51,file='body-D16.txt',status='unknown')
		! elseif(an4nds(1,i)==143) then!D17
			! open(51,file='body-D17.txt',status='unknown')
		! elseif(an4nds(1,i)==144) then!D18
			! open(51,file='body-D18.txt',status='unknown')
		! elseif(an4nds(1,i)==145) then!D19
			! open(51,file='body-D19.txt',status='unknown')
		! elseif(an4nds(1,i)==146) then!D20
			! open(51,file='body-D20.txt',status='unknown')
		! elseif(an4nds(1,i)==147) then!D21
			! open(51,file='body-D21.txt',status='unknown')
		! elseif(an4nds(1,i)==148) then!D22
			! open(51,file='body-D22.txt',status='unknown')
		! elseif(an4nds(1,i)==149) then!D23
			! open(51,file='body-D23.txt',status='unknown')
		! elseif(an4nds(1,i)==150) then!D24
			! open(51,file='body-D24.txt',status='unknown')
		! elseif(an4nds(1,i)==151) then!D25
			! open(51,file='body-D25.txt',status='unknown')
		! elseif(an4nds(1,i)==152) then!D26
			! open(51,file='body-D26.txt',status='unknown')
		! elseif(an4nds(1,i)==153) then!D27
			! open(51,file='body-D27.txt',status='unknown')
		! elseif(an4nds(1,i)==154) then!D28
			! open(51,file='body-D28.txt',status='unknown')
		! elseif(an4nds(1,i)==155) then!D29
			! open(51,file='body-D29.txt',status='unknown')
		! elseif(an4nds(1,i)==156) then!D30
			! open(51,file='body-D30.txt',status='unknown')
		! elseif(an4nds(1,i)==157) then!D31
			! open(51,file='body-D31.txt',status='unknown')
		! elseif(an4nds(1,i)==158) then!D32
			! open(51,file='body-D32.txt',status='unknown')
		! elseif(an4nds(1,i)==159) then!D33
			! open(51,file='body-D33.txt',status='unknown')
		! elseif(an4nds(1,i)==160) then!D34
			! open(51,file='body-D34.txt',status='unknown')
		! elseif(an4nds(1,i)==161) then!D35
			! open(51,file='body-D35.txt',status='unknown')
		! elseif(an4nds(1,i)==162) then!D36
			! open(51,file='body-D36.txt',status='unknown')
! !E
		! elseif(an4nds(1,i)==163) then!126~161-E1
			! open(51,file='body-E1.txt',status='unknown')
		! elseif(an4nds(1,i)==164) then!E2
			! open(51,file='body-E2.txt',status='unknown')
		! elseif(an4nds(1,i)==165) then!E3
			! open(51,file='body-E3.txt',status='unknown')
		! elseif(an4nds(1,i)==166) then!E4
			! open(51,file='body-E4.txt',status='unknown')
		! elseif(an4nds(1,i)==167) then!E5
			! open(51,file='body-E5.txt',status='unknown')
		! elseif(an4nds(1,i)==168) then!E6
			! open(51,file='body-E6.txt',status='unknown')	
		! elseif(an4nds(1,i)==169) then!E7
			! open(51,file='body-E7.txt',status='unknown')
		! elseif(an4nds(1,i)==170) then!E8
			! open(51,file='body-E8.txt',status='unknown')
		! elseif(an4nds(1,i)==171) then!E9
			! open(51,file='body-E9.txt',status='unknown')
		! elseif(an4nds(1,i)==172) then!E10
			! open(51,file='body-E10.txt',status='unknown')
		! elseif(an4nds(1,i)==173) then!E11
			! open(51,file='body-E11.txt',status='unknown')
		! elseif(an4nds(1,i)==174) then!E12
			! open(51,file='body-E12.txt',status='unknown')
		! elseif(an4nds(1,i)==175) then!E13
			! open(51,file='body-E13.txt',status='unknown')
		! elseif(an4nds(1,i)==176) then!E14
			! open(51,file='body-E14.txt',status='unknown')
		! elseif(an4nds(1,i)==177) then!E15
			! open(51,file='body-E15.txt',status='unknown')
		! elseif(an4nds(1,i)==178) then!E16
			! open(51,file='body-E16.txt',status='unknown')
		! elseif(an4nds(1,i)==179) then!E17
			! open(51,file='body-E17.txt',status='unknown')
		! elseif(an4nds(1,i)==180) then!E18
			! open(51,file='body-E18.txt',status='unknown')
		! elseif(an4nds(1,i)==181) then!E19
			! open(51,file='body-E19.txt',status='unknown')
		! elseif(an4nds(1,i)==182) then!E20
			! open(51,file='body-E20.txt',status='unknown')
		! elseif(an4nds(1,i)==183) then!E21
			! open(51,file='body-E21.txt',status='unknown')
		! elseif(an4nds(1,i)==184) then!E22
			! open(51,file='body-E22.txt',status='unknown')
		! elseif(an4nds(1,i)==185) then!E23
			! open(51,file='body-E23.txt',status='unknown')
		! elseif(an4nds(1,i)==186) then!E24
			! open(51,file='body-E24.txt',status='unknown')
		! elseif(an4nds(1,i)==187) then!E25
			! open(51,file='body-E25.txt',status='unknown')
		! endif
		if(an4nds(1,i)==1) then 
			open(51,file='body-030st000dp000.txt',status='unknown')
			loca = '#location=-3.0km off main fault(near side),0.0km along strike,0 km depth'
		elseif(an4nds(1,i)==2) then
			open(51,file='body-020st000dp000.txt',status='unknown')
			loca = '# location=-2.0km off main fault(near side),0.0km along strike,0km depth'
		elseif(an4nds(1,i)==3) then
			open(51,file='body-010st000dp000.txt',status='unknown')
			loca = '#location=-1.0km off main fault(near side),0.0km along strike,0km depth'
		elseif(an4nds(1,i)==4) then
			open(51,file='body010st000dp000.txt',status='unknown')
			loca = '#location=1.0km off main fault(near side),0.0km along strike,0km depth'
		elseif(an4nds(1,i)==5) then
			open(51,file='body-030st120dp000.txt',status='unknown')
			loca = '#location=-3.0km off main fault(near side),12.0km along strike,0km depth'
		elseif(an4nds(1,i)==6) then
			open(51,file='body030st120dp000.txt',status='unknown')
			loca = '#location=3.0km off main fault(near side),12.0km along strike,0km depth'
		endif		
		write(51,*) '#TPV8'
		write(51,*) '#Author=D.Liu'
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
				'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
				':',time_array(7)
		write(51,*) '#code=EQdyna3d'
		write(51,*) '#code_version = 4.0'
		write(51,*) '# element_size =',dx
		write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
		write(51,'( a19,i6)') '# num_time_steps =',locplt
		!write(51,*) loca
		write(51,*) '# Time series in 7 columns in format e15.7'
		write(51,*) '# Column #1 = Time (s)'
		write(51,*) '# Column #2 = horizontal displacement (m)'
		write(51,*) '# Column #3 = horizontal velocity (m/s)'
		write(51,*) '# Column #4 = vertical displacement (m)'
		write(51,*) '# Column #5 = vertical velocity (m/s)'
		write(51,*) '# Column #6 = normal displacement (m)'
		write(51,*) '# Column #7 = normal velocity (m/s)'
		write(51,*) '#'
		write(51,*) '# The line below lists the names of the data fields:'
		write(51,*) 't h-disp h-vel v-disp v-vel n-disp n-vel'
		do j=1,locplt
			write(51,'( f12.5,6e18.7e4)') dout(1,j),dout((i-1)*6+2,j), &
			dout((i-1)*6+3,j),-dout((i-1)*6+6,j),-dout((i-1)*6+7,j), &
			dout((i-1)*6+4,j),dout((i-1)*6+5,j)
		enddo
		close(51)
	enddo
endif
!...write out surface displacement (final) and plastic strain
! B.D. 5/29/12
!dpyn = .false.
!do i=1,numel
!	if(abs(x(3,ien(5,i)))<1) then  !at surface
!		if(x(1,ien(5,i))>fltxyz(1,1,1)-3000 .and. x(1,ien(5,i))<fltxyz(2,1,1)+3000 &
!			.and. x(2,ien(5,i))>=-10000 .and. x(2,ien(5,i))<=10000) then
!			if(dpyn==.false.) then
!				dpyn=.true.
!				open(unit=ioutdp,file=outdp,status='unknown')
!			endif
!			write(ioutdp,'(1x,3f10.1,3f10.3,e10.2)') (x(j,ien(5,i)),j=1,3),&
!					(d(j,ien(5,i)),j=1,3),pstrain(i)
!		endif
!	endif
!enddo
!if(dpyn==.true.) close(ioutdp)
!...write out rupture time for multiple faults.
! also, final slip strike, dip, normal components. B.D. 6/2/12 
if(nftnd(1) > 0) then
	open(unit=ioutrt1,file=outft1,status='unknown')	!rupture time
	write(ioutrt1,'(1x,3f10.1,f10.3,3f10.3,f20.3)') ((x(j,nsmp(1,i,1)),j=1,3),fnft(i,1),&
		(fltslp(j,i,1),j=1,3),miuonf(i),i=1,nftnd(1))
	! flush_ for IBM systems
	!call flush_(ioutrt)
	!call flush(ioutrt)
	close(ioutrt1)
endif
! no 2nd fault, comment here!
!  if(nftnd(2) > 0) then
!    open(unit=ioutrt2,file=outft2,status='unknown')	!rupture time
!    write(ioutrt2,'(1x,4f15.5)') ((x(j,nsmp(1,i,2)),j=1,3),fnft(i,2),i=1,nftnd(2))
! flush_ for IBM systems
!call flush_(ioutrt)
!call flush(ioutrt)
!    close(ioutrt2)
!  endif
!...write out final slip. B.D. 8/11/10
!  if(nftnd > 0) then
!    open(ioutsl,file=outsl, status='unknown')
!    write(ioutsl,'(1x,5f15.5)') ((x(j,nsmp(1,i)),j=1,3),(ftslp(j,i),j=1,2),i=1,nftnd)
!    close(ioutsl)
!  endif
!...general information output
!......execution times, etc. 
!if (me == master) then
!open(unit=ioutgl,file=outgl,status='unknown')		!general output
!else
!  open(unit=ioutgl,file=outgl,position='append')	!general output
!endif
! if (me == master) then
!write(ioutgl,2000) timeused
!write(ioutgl,*) 'mpi_send + mpi_recv time:',btime,' seconds'
!write(ioutgl,*) 'The model has total elements:',numel
!write(ioutgl,*) 'The model has total nodes:',numnp
!$  write(ioutgl,*) 'Number of OpenMP threads used:', omp_get_max_threads()
!!$  write(ioutgl,*) 'Number of processors:', omp_get_num_procs()
!write(ioutgl,*) 'Number of MPI processes used:', nprocs 
!......output general nodal time histories
!write(ioutgl,*)'Proc#',me,': nodes in n4nds and in numnp first'
!write(ioutgl,'(1x,2i10)') ((an4nds(j,i),j=1,2),i=1,n4out)
!write(ioutgl,*) 'nodal time histories: time first, then node by node'
!write(ioutgl,*) 'Simulation is terminated at ', time, ' s of 15 s.'
!write(ioutgl,*) 'Moment rate and max slip vel are:', momntratall, maxslpratall
!close(ioutgl)
!    
!input format  
1000 format(20a4)
!output format
2000 format(1x,' e x e c u t i o n   t i m i n g   i n f o r m a t i o n' &
	 ///5x, &
		 '              Input   p h a s e        = ',1pe15.6//5x, &
		 '              qdct2: form_lhs          = ',1pe15.6//5x, &
		 '              update: d and v          = ',1pe15.6//5x, &
		 '              qdct3: form_rhs          = ',1pe15.6 //5x, &
		 '              hourglass subroutine     = ',1pe15.6 //5x, &
		 '              faulting subroutines     = ',1pe15.6 //5x, &
	 '   Total execution time in seconds     = ',1pe15.6 //)
!
! MPI finalize
!call mpi_barrier(MPI_COMM_WORLD, ierr)
call MPI_Finalize(ierr)
stop
end PROGRAM eqdy3d
