!===================================================================!
!============Introduction to & Major versions of EQdyna3d===========!
!----------------------------Version 1.0----------------------------!
! EQdyna3d is a 3D finite element code to simulate spontaneous 
!	dynamic earthquake ruptures. It has been under development by 
!	Dr. Benchun Duan (B.Duan) since 2005/08/16.
!-------------------------------------------------------------------!
! Hybrid 1D MPI/OpenMP Implementation of EQdyna3d on 2009/02 by
! Dr. Xingfu Wu , CSE Department, Texas A&M University.
!----------------------------Version 2.3.2--------------------------! 
! Partitioning was implemented for better performance of the parallel
! 	EQdyna3d by Drs. Benchun Duan and Xingfu Wu on 2009/04
!----------------------------Version 3.1----------------------------!
! Drucker-Prager Plastic yielding was implemented by B.Duan.
!----------------------------Version 3.2.1--------------------------!
! Perfectly Matched Layer absorbing boundary conditions &
!	Coarsed-grain Q modeling were implemented by Dunyu Liu (D.Liu)
!	on 2015/01/23
! Also, elastic and plastic versions of the code was integrated.
!----------------------------Version 4.0----------------------------! 
! 3D MPI partitioning has been achieved by Bin Luo (B.Luo) on .
! Later the schemes were transferred to Ver 4.0 by D.Liu on 2016/02/19
!----------------------------Version 4.1----------------------------!
! The system has been reshaped for easier use and maintenance by D. Liu on 2016/10/2 
! It is testified against TPV8.
!===================================================================!
PROGRAM eqdy3d
use globalvar
implicit none
include 'mpif.h'
!===================================================================!
character(len=30)::outgl,outft1,outft2,mm,sttmp,dptmp,bodytmp,fsur
character(len=90)::loca	
!===================================================================!
integer(kind=4)::n,i,i1,j,k,l,itemp1,itemp2,itemp3,alloc_err,ierr,me,istatus(MPI_STATUS_SIZE),time_array(8) 		
integer(kind=4)::numnp,numel,neq,n4out,ndout=0,n4onf,maxm,maxs,nsurnd
integer(kind=4),allocatable,dimension(:)::nftnd,surid,id1,ids,et,locid,dof1
integer(kind=4),allocatable,dimension(:,:)::ien,anonfs
integer(kind=4),allocatable,dimension(:,:,:)::nsmp
!===================================================================!
real(kind=8)::timebegin,timeover,PMLb(8)
real(kind=8),allocatable,dimension(:,:)::x,d,v,mat,shl,&
	fnft,arn,r4nuc,arn4m,slp4fri,surnode,state
real(kind=8),allocatable,dimension(:)::brhs,v1,d1,s1,miuonf,vponf,eleporep,pstrain 
real(kind=8),allocatable,dimension(:,:,:)::fric,un,us,ud,fltsta,fltslp 
!==============================PHASE1===============================!
!-------------Global control of the EQdyna3d_v4.0-------------------!
if (C_elastic==0.and.C_Q==1) then
	write(*,*) 'Q model can only work with elastic code'
	stop 1001 
endif
if (C_Q==1.and.rat>1) then
	write(*,*) 'Q model can only work with uniform element size'
	write(*,*) 'rat should be 1.0'
	stop 1002
endif
nstep=idnint(term/dt)
rdampk=rdampk*dt
fltxyz(1,1,1)=fxmin
fltxyz(2,1,1)=fxmax
fltxyz(1,2,1)=fymin
fltxyz(2,2,1)=fymax
fltxyz(1,3,1)=fzmin
fltxyz(2,3,1)=fzmax
fltxyz(1,4,1)=fstrike*pi/180.
fltxyz(2,4,1)=fdip*pi/180.
ccosphi=coheplas*dcos(atan(bulk))
sinphi=dsin(atan(bulk))
!==============================PHASE2===============================!
!-------------------------------------------------------------------!
!INITIATION of 3D MPI.
call MPI_Init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,me,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
timebegin=MPI_WTIME()
time1=MPI_WTIME()	

write(mm,'(i6)') me
mm=trim(adjustl(mm))

nplpts=0	!initialize number of time history plot
if (nhplt>0) then
	nplpts=int(nstep/nhplt)+2
endif
!*** record starting time ***
	
!......directly give file name
outgl='fem.txt'//mm
outft1='frt.txt'//mm
outft2='frt2.txt'//mm

allocate(nftnd(ntotft),shl(nrowsh,nen))
call qdcshl(shl)
!PRE-MESHING to GET NUMBERS for ARRAYS' DEFINITTION and DIMENSIONS
call mesh4num(numnp,numel,neq,PMLb,maxm,nftnd,me,nsurnd)
write(*,*) 'me=',me,'maxm=',maxm
!ALLOCATE NODAL ARRAYS
allocate(id1(maxm),locid(numnp),dof1(numnp),x(ndof,numnp),surid(nsurnd),surnode(nsurnd,2),stat=alloc_err) 
surid=0
surnode=0.0 
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'Insufficient memory to allocate nodal arrays'
endif
!ALLOCATE ELEMENTRAL ARRAYS
allocate(ien(8,numel),mat(numel,5),et(numel),eleporep(numel),pstrain(numel),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'Insufficient memory to allocate elementral arrays'
endif
eleporep = 0.0
pstrain = 0.0
!ALLOCATE FAULTING ARRAYS
nftmx=maxval(nftnd) !max fault nodel num for all faults, used for arrays.
if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
nonmx=sum(nonfs)    !max possible on-fault stations number
allocate(nsmp(2,nftmx,ntotft),fnft(nftmx,ntotft),miuonf(nftmx),vponf(nftmx),un(3,nftmx,ntotft),&
			us(3,nftmx,ntotft),ud(3,nftmx,ntotft),fric(20,nftmx,ntotft),&
			arn(nftmx,ntotft),r4nuc(nftmx,ntotft),anonfs(3,nonmx),&
			arn4m(nftmx,ntotft),slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft),state(nftmx,ntotft))
!ALLOCATE 3DMPI ARRAYS			
allocate(fltgm(nftmx))  !MPIxyz
fltgm=0  !MPIxyz
!INITIATION before calling meshgen
nsmp = 0    !initialize here
fnft = 1000. !as failure time initialization: should larger than actual!
miuonf=0.0
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
state=0.0
!==============================PHASE3===============================!
!--------------------------MESH GENERATION--------------------------!
write(*,*) 'before meshgen me=',me    
allocate(ids(numel))
allocate(s1(4*maxm))
s1=0.0
call meshgen(numnp,numel,ien,mat,s1,eleporep,x,neq,id1,maxm,locid,dof1,et,PMLb,maxs,ids, &
			n4out,nftnd,n4onf,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,arn4m,me,nsurnd,surnode,surid,miuonf,vponf)			
write(*,*) 'after meshgen me=',me 

if(n4onf<=0) n4onf=1 
allocate(fltsta(10,nplpts-1,n4onf),stat=alloc_err)
if(alloc_err /=0) then
	write(*,*) 'me= ',me,'Insufficient memory to allocate fltsta'
endif
fltsta = 0.0

!ALLOCATE BRHS V D
allocate(brhs(neq),v1(neq),d1(neq),v(ndof,numnp),d(ndof,numnp),stat=alloc_err)
if(alloc_err/=0) then
	write(*,*) 'me= ',me,'Insufficient memory to allocate arrays brhs,v,or d'
endif
brhs=0.0
v1=0.0
d1=0.0
v=0.0
d=0.0
!INITIALIZE nodal output time-history data
if(n4out>0) then 
	ndout=n4out*ndof*noid!3 components of 2 quantities: v and d
	!   write(*,*) 'ndout= ',ndout    
	allocate(idhist(3,ndout),dout(ndout+1,nplpts),stat=alloc_err)
	if(alloc_err /=0) then
		write(*,*) 'me= ',me,'insufficient space to allocate array idhist or dout'
	endif

	idhist=0
	dout=0.0
	i1=0
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
 
!INITIALIZE kinematic data for these time histories
! dout(1,locplt) = time
! do i=1,ndout
	! itemp1=idhist(1,i)	!node number
	! itemp2=idhist(2,i)	!i dof
	! itemp3=idhist(3,i)	!code of value
	! if(itemp3==1) then
		! dout(i+1,locplt)=d(itemp2,itemp1)
	! elseif(itemp3==2) then
		! dout(i+1,locplt)=v(itemp2,itemp1)
	! elseif(itemp3==3) then
		! k=id1(locid(itemp1)+itemp2)
		! dout(i+1,locplt)=brhs(k)
	! endif
! enddo    				
endif
!TIME timeused(1) used before SOLVER
time2=MPI_WTIME()		
timeused(1)=time2-time1
!==============================PHASE4===============================!
!------------------------------SOLVER-------------------------------!
if (iexec==1) then	!otherwise data check only
	!write(*,*) 'before driver me=', me,'ndout=',ndout
	call driver(numel,numnp,neq,nftnd,ndout,x,brhs,d,v,mat,ien,eleporep,pstrain, &
			id1,maxm,locid,dof1,et,v1,d1,PMLb,maxs,ids,s1,shl,n4onf,nsmp,fnft,fltslp,un, &
			us,ud,fric,arn,r4nuc,anonfs,arn4m,slp4fri,fltsta,me,nsurnd,surid,miuonf,state)
endif
write(*,*) 'after driver me=',me    
!write(*,*) 'Total seconds consumed = ',timeused(7)
!==============================PHASE5===============================!
time1=MPI_WTIME()
!------------------------------OUTPUT-------------------------------!
!---------------------OutPut surface nodes coordinates--------------!
!----------------------Oct.1.2015/ D.Liu----------------------------!
! if (nsurnd>0) then
	! fsur='surnode.txt'//mm 
	! open(1001,file=fsur,form='formatted',status='unknown')
	! write(1001,'(1x,2f10.1)') (surnode(i1,1),surnode(i1,2),i1=1,nsurnd)
	! close(1001)
! endif
!-------------------------------------------------------------------!
!----------------------------EVENT SUMMARY--------------------------!
open(51,file='event_summary.txt',status='unknown')
write(51,'( i10)') 321
write(51,'( f10.1)') dx
write(51,'( i10)') 151
write(51,'( f10.1)') dx
write(51,'( i10)') int(locplt/nhplt1)
write(51,'( f10.5)') dt*2
close(51)
!-------------------------------------------------------------------!
!--------------------------ON-FAULT STATIONS------------------------!
if(n4onf>0) then
	do i=1,n4onf
		j=anonfs(3,i)
		if(j==1)  then  !main fault stations
			sttmp = '      '
			dptmp = '      '
			write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)/100.d0) 
			write(dptmp,'(i4.3)') int(abs(xonfs(2,anonfs(2,i),j))/100.d0) 
			open(51,file='faultst'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

			sttmp = '      '
			dptmp = '      '
			write(sttmp,'(f5.1)') xonfs(1,anonfs(2,i),j)/1000.d0 
			write(dptmp,'(f5.1)') abs(xonfs(2,anonfs(2,i),j)/1000.d0) 
			loca = '# location = on fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km down-dip'		
		endif
		write(51,*) '# ',projectname
		write(51,*) '# Author=',author
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
			'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
			':',time_array(7)
		write(51,*) '# code = EQdyna3d'
		write(51,*) '# code_version = 4.1'
		write(51,*) '# element_size =',dx
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
		write(51,*) '# Column #8 = normal stress (MPa)'
		write(51,*) '# Column #9 = state variable psi (dimensionless)'
		write(51,*) '# The line below lists the names of the data fields:'
		write(51,*) 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress psi'
		write(51,*) '#'
		do j=1,locplt-1
			write(51,'( f12.5,8e18.7e4)') fltsta(1,j,i),fltsta(5,j,i),fltsta(2,j,i),fltsta(8,j,i)/1.e6,&
					-fltsta(6,j,i),-fltsta(3,j,i),-fltsta(9,j,i)/1.e6,fltsta(10,j,i)/1.e6,fltsta(4,j,i)
		enddo
		close(51)
	enddo
endif 
!-------------------------------------------------------------------!
!-------------------------OFF-FAULT STATIONS------------------------!
if(n4out>0) then
	do i=1,n4out
		bodytmp = '      '
		sttmp = '      '
		dptmp = '      '
		write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i))/100.d0) 
		write(sttmp,'(i4.3)') int(x4nds(1,an4nds(1,i))/100.d0) 
		write(dptmp,'(i4.3)') int(abs(x4nds(3,an4nds(1,i)))/100.d0) 
		write(*,*) 'xcoor',bodytmp,sttmp,dptmp,me
		open(51,file='body'//trim(adjustl(bodytmp))//'st'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

		bodytmp = '      '
		sttmp = '      '
		dptmp = '      '
		write(bodytmp,'(f5.1)') x4nds(2,an4nds(1,i))/1000. 
		write(sttmp,'(f5.1)') x4nds(1,an4nds(1,i))/1000. 
		write(dptmp,'(f5.1)') abs(x4nds(3,an4nds(1,i)))/1000. 
		loca = '# location = '//trim(adjustl(bodytmp))//' km off fault, '//trim(adjustl(sttmp))//' km along strike'//trim(adjustl(dptmp))//' km depth'
		write(51,*) '# ',projectname
		write(51,*) '# Author=',author
		call date_and_time(values=time_array)
		write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
				'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
				':',time_array(7)
		write(51,*) '# code = EQdyna3d'
		write(51,*) '# code_version = 4.1'
		write(51,*) '# element_size =',dx
		write(51,'( a14,f8.4,a3)') '# time_step=', dt, ' s'
		write(51,'( a19,i6)') '# num_time_steps=',locplt
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
!-------------------------------------------------------------------!
!------------------------RUPTURE TIME CONTOURS----------------------!
if(nftnd(1) > 0) then
	open(unit=ioutrt1,file=outft1,status='unknown')	!rupture time
	write(ioutrt1,'(1x,3f10.1,f10.3,3f10.3,f20.3)') ((x(j,nsmp(1,i,1)),j=1,3),fnft(i,1),&
		(fltslp(j,i,1),j=1,3),miuonf(i),i=1,nftnd(1))
	close(ioutrt1)
endif
time2=MPI_WTIME()
timeused(8)=time2-time1 
timeover=MPI_WTIME()
timeused(9)=timeover-timebegin 
!-------------------------------------------------------------------!
!----------------------------TIME EXPENSE---------------------------!
open(unit=ioutrt1,file='timeinfo'//mm,status='unknown')	!rupture time
write(ioutrt1,'(1x,10e18.7e4,2i10)') (timeused(i),i=1,9),btime,numel,neq
close(ioutrt1)
!-------------------------------------------------------------------!
call MPI_Finalize(ierr)
stop
end PROGRAM eqdy3d
