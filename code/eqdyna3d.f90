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
real (kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax,dx
real (kind=8),allocatable,dimension(:,:,:) :: fltxyz !4: x,y,z,strike/dip
!......time histories
integer (kind=4),dimension(2,78)::an4nds
integer (kind=4),allocatable,dimension(:,:) :: idhist
real (kind=4),allocatable,dimension(:,:) :: dout
!......nodes' arrays
!integer (kind=4),allocatable,dimension(:,:) :: id
real (kind=8),allocatable,dimension(:,:) :: x
real (kind=8),allocatable,dimension(:) :: brhs
real (kind=8),allocatable,dimension(:,:) :: d,v
!...element arrays
integer (kind=4),allocatable,dimension(:) :: mat	    
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
ccosphi,sinphi,vs
real (kind=8),allocatable,dimension(:,:,:) :: c
!...
integer (kind=4),dimension(2,2)::nr4mpi,er4mpi  !node,equation range for MPI
!...moment, moment rate, max slip rate at the end. B.D. 8/11/10
real (kind=8)::momntall,momntratall,maxslpratall,magn
!...on-fault station info. B.D. 8/11/10
!  integer (kind=4),dimension(2) :: nonfs=(/8,6/)
!  real (kind=8),dimension(2,8,2) :: xonfs  !8 is max of nonfs
integer (kind=4),dimension(1) :: nonfs=(/6/)
real (kind=8),dimension(2,6,1) :: xonfs
character (len=90) :: loca
integer time_array(8)
!...working variables. B.D. 1/15/12
integer (kind=4) :: itmp1,itmp2
real (kind=8) :: tmp1,tmp2
logical:: sts=.false.,dpyn=.false.
!*.* Variables for new features in PML. D.L. Jan/23/15
integer(kind=4)::maxm,maxs
real(kind=8),dimension(6)::PMLb!Adding PMLb(6)
real (kind=8),allocatable,dimension(:) :: v1,d1,s1 ! 1D array for velocity,displacement and stress 
integer (kind=4),allocatable,dimension(:) :: id1,ids,et,locid,dof1 !et: element type.
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
allocate(emod(numat),pois(numat),mushr(numat),rho(numat),vs(numat),&
		rdampm(numat),rdampk(numat),ccosphi(numat),sinphi(numat),&
		nftnd(ntotft),fltxyz(2,4,ntotft))
call parcon(dx,xmin,xmax,ymin,ymax,zmin,zmax,fltxyz,&
			emod,pois,mushr,rho,vs,rdampm,rdampk,ccosphi,sinphi)

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
call mesh4num(fltxyz,dx,xmin,xmax,ymin,ymax,zmin,zmax, &
				numnp,numel,neq,PMLb,maxm,nftnd,master,me,nprocs)!Adding PMLb and maxm
!  write(*,*) 'me=', me, '# of fault node',(nftnd(i),i=1,ntotft),'# of elements',numel
write(*,*) 'me=',me,'maxm=',maxm
!*** ALLOCATE ARRAYS ***
!...nodes' arrays
allocate(id1(maxm),locid(numnp),dof1(numnp),x(ndof,numnp),stat=alloc_err)!Adding id1(maxm),loci(2,numnp) !Delete id(ndof,numnp),
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate arrays id or x'
endif
!...element's arrays
allocate(ien(8,numel),mat(numel),et(numel),&! Delete lm(ned,nen,numel),elestress(nstr,numel),
		eleporep(numel),pstrain(numel),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'insufficient space to allocate arrays ien, mat, or lm'
endif
!elestress = 0.0  !initialize
eleporep = 0.0
pstrain = 0.0
!...fault arrays
nftmx = maxval(nftnd) !max fault nodel num for all faults, used for arrays.
if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
nonmx = sum(nonfs)    !max possible on-fault stations number
allocate(nsmp(2,nftmx,ntotft),fnft(nftmx,ntotft),un(3,nftmx,ntotft),&
			us(3,nftmx,ntotft),ud(3,nftmx,ntotft),fric(6,nftmx,ntotft),&
			arn(nftmx,ntotft),r4nuc(nftmx,ntotft),anonfs(3,nonmx),&
			arn4m(nftmx,ntotft),slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft))
!after allocate, initialize here before call meshgen
nsmp = 0    !initialize here
fnft = 1000. !as failure time initialization: should larger than actual!
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
!  write(*,*) 'before meshgen me=',me    
allocate(ids(numel))
allocate(s1(2*maxm))
s1=0.0
call meshgen(fltxyz,dx,xmin,xmax,ymin,ymax,zmin,zmax, &
			numnp,numel,ien,mat,s1,eleporep,rho,vs, &!change elestress as s1
			x,neq, & !Delete id
			id1,maxm,locid,dof1,et,PMLb,maxs,ids, &!Adding id1,maxm,loci,PMLb,maxs,ids
			nr4mpi,er4mpi,n4out,an4nds, & 
			nftnd,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,arn4m,master,me,nprocs)			
write(*,*) 'after meshgen me=',me  
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
call mpi_barrier(MPI_COMM_WORLD, ierr)
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
  write(*,*) 'before driver me=', me,'ndout=',ndout
	call driver(numel,numnp,neq,nftnd,ndout,dout,idhist, &
				x,brhs,d,v,mat,ien,eleporep,pstrain, & ! Delete id elestress
				id1,maxm,locid,dof1,et,v1,d1,PMLb,maxs,ids,s1, & !Adding id1(:),maxm,loci(:,:),v1,d1,PMLb
				shl, & ! Delete lm
				emod,pois,mushr,rho,rdampm,rdampk,c,ccosphi,sinphi,&
				er4mpi,nr4mpi,n4onf,momntall,momntratall,maxslpratall,&
				nftmx,nonmx,nsmp,fnft,fltslp,un,us,ud,fric,arn,r4nuc,&
				anonfs,arn4m,slp4fri,fltsta,master,me,nprocs)
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
write(51,'( f10.1)') 100.0
write(51,'( i10)') 151
write(51,'( f10.1)') 100.0
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
			if(xonfs(1,anonfs(2,i),j)== -6000 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst-060dp000.txt',status='unknown')
				loca = '# location = on mian fault, -6.0 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)==0 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst000dp000.txt',status='unknown')
				loca = '# location = on main fault, 0.0 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== 12000 .and. xonfs(2,anonfs(2,i),j)==0) then
				open(51,file='faultst120dp000.txt',status='unknown')
				loca = '# location = on mian fault, 12.0 km along strike, 0 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== -6000 .and. xonfs(2,anonfs(2,i),j)==-7500) then
				open(51,file='faultst-060dp075.txt',status='unknown')
				loca = '# location = on mian fault, -6.0 km along strike, 7.5 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)==0 .and. xonfs(2,anonfs(2,i),j)==-7500) then
				open(51,file='faultst000dp075.txt',status='unknown')
				loca = '# location = on main fault, 0.0 km along strike, 7.5 km down-dip'
			elseif(xonfs(1,anonfs(2,i),j)== 12000 .and. xonfs(2,anonfs(2,i),j)==-7500) then
				open(51,file='faultst120dp075.txt',status='unknown')
				loca = '# location = on mian fault, 12.0 km along strike, 7.5 km down-dip'
			endif
			!     elseif(j==2) then  !branch fault stations
			!      if(xonfs(1,anonfs(2,i),j)==2000 .and. xonfs(2,anonfs(2,i),j)==0) then
			!        open(51,file='branchst020dp000.txt',status='unknown')
			!        loca = '# location = on branch fault, 2.0 km along strike, 0 km down-dip'
			!      elseif(xonfs(1,anonfs(2,i),j)==5000 .and. xonfs(2,anonfs(2,i),j)==0) then
			!        open(51,file='branchst050dp000.txt',status='unknown')
			!        loca = '# location = on branch fault, 5.0 km along strike, 0 km down-dip'
			!      elseif(xonfs(1,anonfs(2,i),j)==9000 .and. xonfs(2,anonfs(2,i),j)==0) then
			!        open(51,file='branchst090dp000.txt',status='unknown')
			!        loca = '# location = on branch fault, 9.0 km along strike, 0 km down-dip'
			!      elseif(xonfs(1,anonfs(2,i),j)==2000 .and. xonfs(2,anonfs(2,i),j)==-7500) then
			!        open(51,file='branchst020dp075.txt',status='unknown')
			!        loca = '# location = on branch fault, 2.0 km along strike, 7.5 km down-dip'
			!      elseif(xonfs(1,anonfs(2,i),j)==5000 .and. xonfs(2,anonfs(2,i),j)==-7500) then
			!        open(51,file='branchst050dp075.txt',status='unknown')
			!        loca = '# location = on branch fault, 5.0 km along strike, 7.5 km down-dip'
			!      elseif(xonfs(1,anonfs(2,i),j)==9000 .and. xonfs(2,anonfs(2,i),j)==-7500) then
			!        open(51,file='branchst090dp075.txt',status='unknown')
			!        loca = '# location = on branch fault, 9.0 km along strike, 7.5 km down-dip'
			!      endif
		endif
		write(51,*) '# For LVFZ project, revised from SCEC TPV19'
		write(51,*) '# author = Benchun Duan'
		!    write(51,*) '# date = 8-20-10'
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
			'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
			':',time_array(7)
		write(51,*) '# code = EQdyna'
		write(51,*) '# code_version = 3.0'
		!write(51,*) '# element_size = 50 m x 50 m on main fault, 57.7 m x 50 m on branch'
		write(51,*) '# element_size = 100 m x 100 m on fault'
		write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
		write(51,'( a19,i6)') '# num_time_steps =', locplt-1
		write(51,*) loca
		write(51,*) '# Time series in 8 columns in format e15.7'
		write(51,*) '# Column #1 = Time (s)'
		write(51,*) '# Column #2 = horizontal slip (m)'
		write(51,*) '# Column #3 = horizontal slip rate (m/s)'
		write(51,*) '# Column #4 = horizontal shear stress (MPa)'
		write(51,*) '# Column #5 = down-dip slip (m)'
		write(51,*) '# Column #6 = down-dip slip rate (m/s)'
		write(51,*) '# Column #7 = down-dip shear stress (MPa)'
		write(51,*) '# Column #8 = normal stress (MPa)'
		write(51,*) '# The line below lists the names of the data fields:'
		write(51,*) 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress'
		write(51,*) '#'
		do j=1,locplt-1
			write(51,'( f12.5,7e18.7e4)') fltsta(1,j,i),fltsta(5,j,i),fltsta(2,j,i),fltsta(8,j,i)/1.e6,&
					-fltsta(6,j,i),-fltsta(3,j,i),-fltsta(9,j,i)/1.e6,fltsta(10,j,i)/1.e6
		enddo
		close(51)
	enddo
endif 
!...output off-fault stations. B.D. 8/11/10
if(n4out > 0) then
	do i=1,n4out
		if(an4nds(1,i)==1) then 
			open(51,file='body-060st-120dp000.txt',status='unknown')
			loca = '# location = -6.0 km off main fault (near side) , -12.0 km along strike, 0 km depth'
		elseif(an4nds(1,i)==2) then
			open(51,file='body-060st000dp000.txt',status='unknown')
			loca = '# location = -6.0 km off main fault (near side), 0.0 km along strike, 0 km depth'
		elseif(an4nds(1,i)==3) then
			open(51,file='body-060st100dp000.txt',status='unknown')
			loca = '# location = -6.0 km off main fault (near side), 10.0 km along strike, 0 km depth'
		elseif(an4nds(1,i)==4) then
			open(51,file='body-100st-120dp000.txt',status='unknown')
			loca = '# location = -10.0 km off main fault (near side), -12.0 km along strike, 0 km depth'
		elseif(an4nds(1,i)==5) then
			open(51,file='body-100st000dp000.txt',status='unknown')
			loca = '# location = -10.0 km off main fault (near side), 0.0 km along strike, 0 km depth'
		elseif(an4nds(1,i)==6) then
			open(51,file='body-100st100dp000.txt',status='unknown')
			loca = '# location = -10.0 km off main fault (near side), 10.0 km along strike, 0 km depth'
		endif
		write(51,*) '# for LVFZ project, revised from SCEC TPV19'
		write(51,*) '# author = Benchun Duan'
		!   write(51,*) '# date = 8-20-10'
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
				'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
				':',time_array(7)
		write(51,*) '# code = EQdyna'
		write(51,*) '# code_version = 3.0'
		write(51,*) '# element_size = 100 m x 100 m on fault'
		write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
		write(51,'( a19,i6)') '# num_time_steps =',locplt
		write(51,*) loca
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
dpyn = .false.
do i=1,numel
	if(abs(x(3,ien(5,i)))<1) then  !at surface
		if(x(1,ien(5,i))>fltxyz(1,1,1)-3000 .and. x(1,ien(5,i))<fltxyz(2,1,1)+3000 &
			.and. x(2,ien(5,i))>=-10000 .and. x(2,ien(5,i))<=10000) then
			if(dpyn==.false.) then
				dpyn=.true.
				open(unit=ioutdp,file=outdp,status='unknown')
			endif
			write(ioutdp,'(1x,3f10.1,3f10.3,e10.2)') (x(j,ien(5,i)),j=1,3),&
					(d(j,ien(5,i)),j=1,3),pstrain(i)
		endif
	endif
enddo
if(dpyn==.true.) close(ioutdp)
!...write out rupture time for multiple faults.
! also, final slip strike, dip, normal components. B.D. 6/2/12 
if(nftnd(1) > 0) then
	open(unit=ioutrt1,file=outft1,status='unknown')	!rupture time
	write(ioutrt1,'(1x,3f10.1,f10.3,3f10.3)') ((x(j,nsmp(1,i,1)),j=1,3),fnft(i,1),&
		(fltslp(j,i,1),j=1,3),i=1,nftnd(1))
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
open(unit=ioutgl,file=outgl,status='unknown')		!general output
!else
!  open(unit=ioutgl,file=outgl,position='append')	!general output
!endif
! if (me == master) then
write(ioutgl,2000) timeused
write(ioutgl,*) 'mpi_send + mpi_recv time:',btime,' seconds'
write(ioutgl,*) 'The model has total elements:',numel
write(ioutgl,*) 'The model has total nodes:',numnp
!$  write(ioutgl,*) 'Number of OpenMP threads used:', omp_get_max_threads()
!!$  write(ioutgl,*) 'Number of processors:', omp_get_num_procs()
write(ioutgl,*) 'Number of MPI processes used:', nprocs 
!......output general nodal time histories
write(ioutgl,*)'Proc#',me,': nodes in n4nds and in numnp first'
write(ioutgl,'(1x,2i10)') ((an4nds(j,i),j=1,2),i=1,n4out)
write(ioutgl,*) 'nodal time histories: time first, then node by node'
write(ioutgl,*) 'Simulation is terminated at ', time, ' s of 15 s.'
write(ioutgl,*) 'Moment rate and max slip vel are:', momntratall, maxslpratall
close(ioutgl)
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
call mpi_barrier(MPI_COMM_WORLD, ierr)
call MPI_Finalize(ierr)
stop
end PROGRAM eqdy3d