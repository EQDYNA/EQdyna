SUBROUTINE driver(numel,numnp,neq,nftnd,ndout,dout,idhist,x,brhs,d,v, &
mat,ien,eleporep,pstrain,id1,maxm,locid,dof1,et,v1,d1,PMLb,maxs,ids,s1, &
shl,e,pois,mushr,rdampm,rdampk,c,ccosphi,sinphi,er4mpi,nr4mpi,n4onf, &
momntall,momntratall,maxslpratall,itmp1,itmp2,nsmp,fnft,fltslp,&
un,us,ud,fric,arn,r4nuc,anonfs,arn4m,slp4fri,fltsta,master,me,nprocs,&
nsurnd,surid,miuonf,nnodPML,nPMLid,nPMLvp)
use globalvar
implicit none
include 'mpif.h'
!### solution driver program
! As shown above, only fault arrays are transferred.
! Modify to explictly use central difference method as DYNA3D, 
!	rather than alpha-method. B.D. 7/21/05
logical :: lstr,lstr1
integer (kind=4) :: ntstep=0,ifault,nsq,n,i,j,k,k1,l,nel,m
integer (kind=4) ::numel,numnp,neq,ndout,n4onf
!...fault arrays, dynamic defined in case zero fault node. B.D. 10/16/09
integer (kind=4),dimension(ntotft) :: nftnd
!real (kind=8),dimension(ndof,numnp):: acc
!......time histories
integer (kind=4),dimension(3,ndout) :: idhist
real (kind=4),dimension(ndout+1,nplpts) :: dout
!......nodes' arrays
!integer (kind=4),dimension(ndof,numnp) :: id
real (kind=8),dimension(nsd,numnp) :: x
real (kind=8),dimension(neq) :: alhs,brhs
real (kind=8),dimension(ndof,numnp) :: d,v
!...element arrays
real (kind=8),dimension(numel,5) :: mat	    
integer (kind=4),dimension(nen,numel) :: ien
!integer (kind=4),dimension(ned,nen,numel) :: lm
real (kind=8),dimension(nrowsh,nen) :: shl
!real (kind=8),dimension(nstr,numel) :: elestress
real (kind=8),dimension(numel) :: eledet,eleporep,pstrain
real (kind=8),dimension(nee,numel) :: elemass
real (kind=8),dimension(nrowsh-1,nen,numel) :: eleshp
!...fault arrays
integer (kind=4),dimension(3,itmp2) :: anonfs
integer (kind=4),dimension(2,itmp1,ntotft) :: nsmp
real (kind=8),dimension(itmp1,ntotft) :: fnft,arn,r4nuc,arn4m,slp4fri
real (kind=8),dimension(3,itmp1,ntotft) :: un,us,ud,fltslp
real (kind=8),dimension(8,itmp1,ntotft) :: fric
real(kind=8)::miuonf(itmp1)
real (kind=8),dimension(10,nplpts-1,n4onf) :: fltsta
!...material properties
real (kind=8),dimension(numat) :: e,pois,mushr,rdampm,rdampk,&
ccosphi,sinphi
real (kind=8),dimension(nrowc,nrowc,numat) :: c
!... equations,node mass,weight,solutions,coordinate,shape
real (kind=8),dimension(numnp) :: fnms
!...hourglass control arrays
real (kind=8),dimension(6,numel) :: ss
real (kind=8),dimension(nen,4,numel) :: phi
!...OpenMP and MPI
!$ 	integer  omp_get_max_threads
!$     external omp_get_max_threads
!!$     integer  omp_get_num_procs
!!$     external omp_get_num_procs
integer master, me, nprocs, ierr, rlp, rr, jj
real (kind=8), allocatable, dimension(:) :: btmp, btmp1, btmp2, btmp3
!real (kind=8) :: btime
integer (kind=4),dimension(nprocs) :: kk	    
integer istatus(MPI_STATUS_SIZE)
!...
integer (kind=4),dimension(2,2)::er4mpi,nr4mpi  !equation,node range for MPI
!... moment, moment rate, max slip rat
real (kind=8) :: momnt,momntrat,maxslprat,momntall,momntratall,maxslpratall
!...working variables.
integer (kind=4) :: itmp1,itmp2
!*.* Variables for PML layer. D.L. Jan/23/2015
integer (kind=4):: maxm,maxs,itag,eqn,nsurnd
integer (kind=4),dimension(maxm)::id1
integer (kind=4),dimension(numel)::ids
integer (kind=4),dimension(numnp)::locid,dof1
integer (kind=4),dimension(numel)::et
real (kind=8),dimension(neq)::v1,d1
real (kind=8),dimension(maxs)::s1
real(kind=8),dimension(8)::PMLb
real(kind=8),dimension(9)::dampv
integer(kind=4)::surid(nsurnd)
real(kind=8)::maxsur(nsurnd,6)
character(len=6)::mm
character(len=30)::foutmov
!Vp for PML nodes. Oct.5.2015/D. Liu
integer(kind=4)::nnodPML,ncount
integer(kind=4)::nPMLid(nnodPML),nPMLvp(nnodPML)
real(kind=8)::vptemp
!--3DMPI------------------------------------------------------------!
integer(kind=4)::izz,ntagMPI,ix,iy,iz,nx,ny,nz,mex,mey,mez,nodenumtemp
integer(kind=4)::bndl,bndr,bndf,bndb,bndd,bndu,rrr,jjj
!-------------------------------------------------------------------!
!*.* D.L.
!*** write out a reminder for code user ***
!
if (me == master) then
!$  write(*,*) 'Number of OpenMP threads used:', omp_get_max_threads()
!!$  write(*,*) 'Number of processors available:', omp_get_num_procs()
	write(*,*) 'Number of MPI processes used:', nprocs
	write(*,*) 'Please be patient! Program is running on Process...', me
endif
!if (me == master) then
!  write(*,*) 'Please be patient! Program is running on Process...', me
!endif  
!
!... form effective mass matrix
!
time1 = MPI_WTIME()
alhs = 0.0
call qdct2(numel,numnp,neq,shl,ien,x,mat,rdampm,alhs,& !delete lm
	eledet,elemass,eleshp,fnms,ss,phi,er4mpi,nr4mpi,master,me,nprocs,maxm,id1,locid,dof1)!Adding maxm & id1 & loci	
time2 = MPI_WTIME()
timeused(2) = timeused(2) + (time2 - time1)
write(*,*) 'me=',me,'Finished qdct2','time-qdct2',timeused(2)
!
!...allocate/deallicate outside time loop to improve performance.
! found by TAMU supercomputing facility colleagues. B.D. 1/12/12
!  jj = er4mpi(2,1)-er4mpi(1,1)
!  rr = er4mpi(2,2)-er4mpi(1,2)
!  allocate(btmp(rr))
!  allocate(btmp1(rr))
!  allocate(btmp2(jj))
!  allocate(btmp3(jj))
!*** time step loop ***
!
!btime = 0.0  !timing MPI
maxsur=0.0
do n=1,nstep
	!
	time = time + dt	 	! dynamic problem timing
	ntstep = ntstep + 1 
	!...print on screen for monitoring
	if(mod(n,nhshw)==0) then
		if (me == master) then
			write(*,'( a6,f9.3)') 'time=',time
		endif
	endif
	momnt=0.0  !initialize in case no fault node for the MPI
	momntrat=0.0
	maxslprat=0.0
	!
	!*** central difference method to update ***
	!  brhs stores accelarations now, only need to update
	!	 velocity and displacement. Time step dt assumes to 
	!	 be constant now. See DYNA3D theoretical manual 21.3.
	!	 B.D. 7/21/05.
	!  To account for initial vel and acc, update here rather 
	!	 than at the end of time step. B.D. 7/21/05
	!
	time1 = MPI_WTIME()
	!*.* Update velocity and displacement from brhs. D.L. Jan/23/2015
	ncount=0
	do i=1,numnp
		if (dof1(i)==3) then
			do j=1,3
				itag=locid(i)+j
				eqn=id1(itag)
				v1(eqn)=v1(eqn)+brhs(eqn)*dt
				d1(eqn)=d1(eqn)+v1(eqn)*dt
				v(j,i)=v1(eqn)
				d(j,i)=d1(eqn)
			enddo
	!		if (x(1,i)==-24e3.and.x(3,i)==0.0.and.x(2,i)>-90e3.and.x(2,i)<30e3)then
	!			write(mm,'(i6)') me
	!			mm = trim(adjustl(mm))
	!			foutmov='0ProfileA_'//mm
	!			open(5000+me,file=foutmov,form='formatted',position='append')
	!			write(5000+me,'(1x,3e18.7e4)') v(1,i),v(2,i),v(3,i)
	!			if (n==1) then 
	!				write(mm,'(i6)') me
	!				mm = trim(adjustl(mm))
	!				foutmov='0ProfileAcoord_'//mm
	!				open(6000+me,file=foutmov,form='formatted',position='append')
	!					write(6000+me,'(1x,3e18.7e4)') x(1,i),x(2,i),x(3,i)
	!			endif
	!		elseif(x(1,i)==0.and.x(3,i)==0.0.and.x(2,i)>-90e3.and.x(2,i)<30e3)then
	!			write(mm,'(i6)') me
	!			mm = trim(adjustl(mm))
	!			foutmov='0ProfileB_'//mm
	!			open(5000+me,file=foutmov,form='formatted',position='append')
	!			write(5000+me,'(1x,3e18.7e4)') v(1,i),v(2,i),v(3,i)	
	!			if (n==1) then 
	!				write(mm,'(i6)') me
	!				mm = trim(adjustl(mm))
	!				foutmov='0ProfileBcoord_'//mm
	!				open(6000+me,file=foutmov,form='formatted',position='append')
	!					write(6000+me,'(1x,3e18.7e4)') x(1,i),x(2,i),x(3,i)
	!			endif				
	!		elseif(x(1,i)==18e3.and.x(3,i)==0.0.and.x(2,i)>-90e3.and.x(2,i)<30e3)then
	!			write(mm,'(i6)') me
	!			mm = trim(adjustl(mm))
	!			foutmov='0ProfileC_'//mm
	!			open(5000+me,file=foutmov,form='formatted',position='append')
	!			write(5000+me,'(1x,3e18.7e4)') v(1,i),v(2,i),v(3,i)
	!			if (n==1) then 
	!				write(mm,'(i6)') me
	!				mm = trim(adjustl(mm))
	!				foutmov='0ProfileCcoord_'//mm
	!				open(6000+me,file=foutmov,form='formatted',position='append')
	!					write(6000+me,'(1x,3e18.7e4)') x(1,i),x(2,i),x(3,i)
	!			endif				
	!		elseif(x(1,i)==60e3.and.x(3,i)==0.0.and.x(2,i)>-90e3.and.x(2,i)<30e3)then
	!			write(mm,'(i6)') me
	!			mm = trim(adjustl(mm))
	!			foutmov='0ProfileD_'//mm
	!			open(5000+me,file=foutmov,form='formatted',position='append')
	!			write(5000+me,'(1x,3e18.7e4)') v(1,i),v(2,i),v(3,i)
	!			if (n==1) then 
	!				write(mm,'(i6)') me
	!				mm = trim(adjustl(mm))
	!				foutmov='0ProfileDcoord_'//mm
	!				open(6000+me,file=foutmov,form='formatted',position='append')
	!					write(6000+me,'(1x,3e18.7e4)') x(1,i),x(2,i),x(3,i)
	!			endif				
	!		endif
		elseif (dof1(i)==12) then
			itag=locid(i)+1
			eqn=id1(itag)
			if (eqn>0) then
!-------------------------------------------------------------------!
!----Oct.6.2015/ D.Liu--Use Vp in element for nodes in PML----------!			
				ncount=ncount+1 
				!Only check inconsistancy at the first time step.
				if (n==1.and.nPMLid(ncount)/=i) then 
					write(*,*) 'nPMLid(ncount)=',nPMLid(ncount)
					write(*,*) 'i=',i
					write(*,*) 'CoordnPML',x(1,nPMLid(ncount)),x(2,nPMLid(ncount)),x(3,nPMLid(ncount))
					write(*,*) 'Coordi',x(1,i),x(2,i),x(3,i)
					stop 3001
				endif
				!vptemp=mat(nPMLvp(ncount),1)
				call comdampv(x(1,i),x(2,i),x(3,i),PMLb,dampv)
			endif
			do j=1,9
				itag=locid(i)+j
				eqn=id1(itag)
				if (eqn>0) then
					v1(eqn)=(brhs(eqn)+v1(eqn)*(1/dt-dampv(j)/2))/(1/dt+dampv(j)/2)
				endif
			enddo
			do j=10,12
				itag=locid(i)+j
				eqn=id1(itag)
				if (eqn>0) then
					v1(eqn)=v1(eqn)+brhs(eqn)*dt
				endif				
			enddo
			! Update final velocity.
			itag=locid(i)
			eqn=id1(itag+1)
			if (eqn>0) then
				v(1,i)=v1(id1(itag+1))+v1(id1(itag+2))+v1(id1(itag+3))+v1(id1(itag+10))
				v(2,i)=v1(id1(itag+4))+v1(id1(itag+5))+v1(id1(itag+6))+v1(id1(itag+11))
				v(3,i)=v1(id1(itag+7))+v1(id1(itag+8))+v1(id1(itag+9))+v1(id1(itag+12))
				d(1,i)=d(1,i)+v(1,i)*dt
				d(2,i)=d(2,i)+v(2,i)*dt
				d(3,i)=d(3,i)+v(3,i)*dt
			elseif (eqn==-1) then
				v(1,i)=0.0
				v(2,i)=0.0
				v(3,i)=0.0
				d(1,i)=0.0
				d(2,i)=0.0
				d(3,i)=0.0				
			endif	
		endif
		if ((v(1,i)/=v(1,i)).or.v(2,i)/=v(2,i).or.v(3,i)/=v(3,i)) then 
			write(*,*) x(1,i),x(2,i),x(3,i),me,mex,mey,mez
			stop 'NAN'
		endif
	enddo
!-------------------------------------------------------------------!
!-------------------Surface Outputs--Sep.24.2015--------------------!
!----------------------------D. Liu --------------------------------!
	if (nsurnd>0) then
		do i=1,nsurnd
			if (maxsur(i,1)<abs(v(1,surid(i)))) then
				maxsur(i,1)=abs(v(1,surid(i)))
			endif
			if (maxsur(i,2)<abs(v(2,surid(i)))) then
				maxsur(i,2)=abs(v(2,surid(i)))
			endif		
			if (maxsur(i,3)<abs(v(3,surid(i)))) then
				maxsur(i,3)=abs(v(3,surid(i)))
			endif	
			if (maxsur(i,4)<abs(brhs(id1(locid(surid(i))+1)))) then
				maxsur(i,4)=abs(brhs(id1(locid(surid(i))+1)))
			endif	
			if (maxsur(i,5)<abs(brhs(id1(locid(surid(i))+2)))) then
				maxsur(i,5)=abs(brhs(id1(locid(surid(i))+2)))
			endif
			if (maxsur(i,6)<abs(brhs(id1(locid(surid(i))+3)))) then
				maxsur(i,6)=abs(brhs(id1(locid(surid(i))+3)))
			endif			
		enddo
		if(n==nstep) then
			write(mm,'(i6)') me
			mm = trim(adjustl(mm))
			foutmov='fmaxsur_'//mm
			open(unit=2001+me,file=foutmov,form='formatted',status='unknown')	
				write(2001+me,'(1x,6e18.7e4)') ((maxsur(i,j),j=1,6),i=1,nsurnd)
		endif			
		! if(mod(n,ninterval)==0) then
			! write(mm,'(i6)') me
			! mm = trim(adjustl(mm))
			! foutmov='fsurout_'//mm
			! open(unit=2001+me,file=foutmov,form='formatted',status='unknown',position='append')	
				! write(2001+me,'(1x,6f10.3)') ((d(j,surid(i)),j=1,3),(v(j,surid(i)),j=1,3),i=1,nsurnd)
		! endif	
	endif	
!-------------------------------------------------------------------!	
	
	!*.* D.L.
	time2 = MPI_WTIME()
	timeused(3) = timeused(3) + (time2 - time1)
	
	!*** store desired results at set time intervals ***
	if (mod(n,nhplt) == 0) then	
		lstr = .true.	
		locplt = locplt + 1	!when n=1, locplt=2 due to 1 in eqdy3d.f90
	else
		lstr = .false.
	endif
	if(mod(n,nhplt1) == 0) then
		lstr1 = .true.
	else
		lstr1 = .false.
	endif
	if(lstr) then
	! if (me == master) then
	!  write(ioutsl,'( a6,f12.5)') 'time=',time
	! flush_ for IBM systems
	!call flush_(ioutsl)
	! call flush(ioutsl)
	! write(ioutst,'( a6,f12.5)') 'time=',time
	!call flush_(ioutst)
	!  call flush(ioutst)
	! endif
	endif
	if (lstr) then
	!...nodal output
		if((ndout > 0) .and. (locplt > 1)) then
			dout(1,locplt) = time
			do i=1,ndout
				j = idhist(1,i)
				if(j<=0) j=1  !avoid zero that cannot be used below
					k = idhist(2,i)
					l = idhist(3,i)
				if(l == 1) then
					dout(i+1,locplt) = d(k,j)
				elseif(l == 2) then
					dout(i+1,locplt) = v(k,j)
				elseif(l == 3) then
					k1 = id1(locid(j)+k)!id(k,j) !*.* D.L. Jan/23/15
					dout(i+1,locplt) = brhs(k1)
				endif
			enddo
		endif
		!	
	endif
	!...ground motion nearby the fault output. B.D. 9/17/08
	!if(lstr) then
	!  write(ioutgm,*) 'time=',time
	!  do l=1,numnp
	!    if(abs(x(3,l))<10) then
	!    if(x(1,l)>=-250000 .and. x(1,l)<=250000 .and. &
	!       x(2,l)>=-100000 .and. x(2,l)<=100000) then
	!       if(n==1) then	!coordinate for the first time step
	!         write(ioutgm,'( 2f9.1,9f7.2)') (x(j,l),j=1,2),(acc(j,l),v(j,l),d(j,l),j=1,ndof)
	!       else	!just a, v and d
	!         write(ioutgm,'( 9f7.2)') (acc(j,l),v(j,l),d(j,l),j=1,ndof)
	!       endif
	!     endif
	!     endif
	!   enddo            
	!   call flush(ioutgm)
	! endif
	!
	!*** initialize for right hand force to use ***
	!
	brhs = 0.0	!initialize it every interation
	!
	!*** some steps may be needed in future.
	!
	!...evaluate load-time functions at time n+1.
	!   at present, no load. B.D. 7/2/05
	!if(nlvect > 0) then
	!  call lfac
	!endif
	!...overwrite predictors to account for kinematic boundary.
	!   at present, no boundary, no execution of this routine.
	!if(dt /= 0) then
	!  call compbc
	!endif
	!... evaluate load-time functions at time n+1+alpha.
	!	no load so far. B.D. 7/2/05
	!if (nltftn > 0) then  !actually, so far nltftn=0
	!  call lfac
	!endif
	!... form nodal contribution to residual force vector
	!if (nlvect > 0) then  !so far, nltftn=0
	!  call load
	!endif
	!
	!*** form element contribution to residual force vector ***
	! This is most time-consuming routine in this analysis!
	!
	!   !$omp parallel default(shared) shared(v,d,brhs)
	!   !$omp sections 
	time1 = MPI_WTIME()
	!   !$omp section
	call qdct3(numel,numnp,neq,mat,ien,d,v,rdampk,rdampm,ccosphi,sinphi,&
				mushr,eleporep,elemass,eleshp,eledet,pstrain,c,brhs,& ! Delete lm
				me,master,nprocs,maxm,id1,locid,dof1,et,v1,d1,PMLb,x,maxs,ids,s1,n)!Adding maxm & id1 & loci,v1,d1,PMLb,x,maxs,ids,s1
	time2 = MPI_WTIME()
	timeused(4) = timeused(4) + (time2 - time1)
	!
	!*** add hourglass resistence to residual force vector ***
	!  Note: this should be before faulting routine.
	!  B.D. 7/22/05
	!
	time1 = MPI_WTIME()
	! !$omp section
	call hrglss(numel,numnp,neq,ien,d,v,rdampk,mat,ss,phi,brhs,me,master,nprocs,&
		maxm,id1,locid,dof1,et,eledet) ! Delete lm !Add maxm,id1,loci
	time2 = MPI_WTIME()
	timeused(5) = timeused(5) + (time2 - time1)
	! !$omp end sections
	! !$omp end parallel
	time1 = MPI_WTIME()
	
	!****************MPI**********************
	!...Massage Passing for general nprocs
	!Basic facts: me need to send from er4mpi(1,1) to er4mpi(2,1) of 
	!brhs() to its left me-1, and to recv me-1's brhs() from er4mpi(1,2) 
	!to er4mpi(2,2), add these correspondences in each. me also need 
	!send from er4mpi(1,2) to er4mpi(2,2) of brhs() to its right me+1, 
	!and recv me+1's brhs() from er4mpi(1,1) to er4mpi(2,1), add these 
	!correspondences in each.
!**********************************************************************
!***********3D partioning comments*********************!
	!prepare for MPI partitioning
  	mex=int(me/(npy*npz))
  	mey=int((me-mex*npy*npz)/npz)
  	mez=int(me-mex*npy*npz-mey*npz)
	nx=numcount(1)
	ny=numcount(2)
	nz=numcount(3)
!Loop sequence x->z->y
	if (npx>1)then
		rr=numcount(3+1)
		jj=numcount(3+2)
		bndl=1!-x boundary
		bndr=nx!+x boundary
		if (mex==master) then
			bndl=0
		elseif (mex==npx-1) then
			bndr=0
		endif
		
		if (bndl/=0)then 
			!Check
			if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
			rrr=rr+fltnum(1)*3
			allocate(btmp(rrr),btmp1(rrr))
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndl-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-driver-bndl'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(1)) then
          		do ix=1,fltnum(1)
            		nodenumtemp=nx*ny*nz+fltl(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-npy*npz, 100000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-npy*npz, 100000+me-npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndl-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp1(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(1)) then
          		do ix=1,fltnum(1)
            		nodenumtemp=nx*ny*nz+fltl(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)			
		endif !bndl/=0
		!
		if (bndr/=0)then 
			!Check
			if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'
			jjj=jj+fltnum(2)*3
			allocate(btmp2(jjj),btmp3(jjj))
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndr-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-driver-bndr'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(2)) then
          		do ix=1,fltnum(2)
            		nodenumtemp=nx*ny*nz+fltr(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+npy*npz, 100000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+npy*npz, 100000+me+npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndr-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp3(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(2)) then
          		do ix=1,fltnum(2)
            		nodenumtemp=nx*ny*nz+fltr(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)				
		endif !bndr/=0
		!
  	endif!if npx>1
	call mpi_barrier(MPI_COMM_WORLD, ierr)
!****************************************************************************
!******************Paritioning along y axis**********************************
	if (npy>1)then
		rr=numcount(3+3)
		jj=numcount(3+4)
		bndf=1!-y boundary
		bndb=ny!+y boundary
		if (mey==master) then
			bndf=0
		elseif (mey==npy-1) then
			bndb=0
		endif
		
		if (bndf/=0)then 
			!Check
			if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
			rrr=rr+fltnum(3)*3
			allocate(btmp(rrr),btmp1(rrr))
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndf
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-driver-bndf'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(3)) then
          		do iy=1,fltnum(3)!using the free index, here iy.
            		nodenumtemp=nx*ny*nz+fltf(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-npz, 110000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-npz, 110000+me-npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndf
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp1(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(3)) then
          		do iy=1,fltnum(3)
            		nodenumtemp=nx*ny*nz+fltf(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)				
		endif!bndf/=0
		!
		if (bndb/=0)then 
			!Check
			if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'
			jjj=jj+fltnum(4)*3
			allocate(btmp2(jjj),btmp3(jjj))
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndb
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-driver-bndb'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(4)) then
          		do iy=1,fltnum(4)
            		nodenumtemp=nx*ny*nz+fltb(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+npz, 110000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+npz, 110000+me+npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndb
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp3(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(4)) then
          		do iy=1,fltnum(4)
            		nodenumtemp=nx*ny*nz+fltb(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)				
		endif!bndb/=0
		!
  	endif!if npy>1
	call mpi_barrier(MPI_COMM_WORLD, ierr)
!****************************************************************************
!******************Paritioning along z axis**********************************
	if (npz>1)then
		rr=numcount(3+5)
		jj=numcount(3+6)
		bndd=1!-z boundary
		bndu=nz!+z boundary
		if (mez==master) then
			bndd=0
		elseif (mez==npz-1) then
			bndu=0
		endif
		
		if (bndd/=0)then 
			!Check
			if (rr.eq.0) stop 'inconsistency in rr MPI Phase 1'
			rrr=rr+fltnum(5)*3
			allocate(btmp(rrr),btmp1(rrr))
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndd-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-driver-bndd'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(5)) then
          		do iz=1,fltnum(5)!using the free index, here iz.
            		nodenumtemp=nx*ny*nz+fltd(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-1, 120000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-1, 120000+me-1,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndd-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp1(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(5)) then
          		do iz=1,fltnum(5)
            		nodenumtemp=nx*ny*nz+fltd(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)				
		endif!bndd/=0
		!
		if (bndu/=0)then 
			!Check
			if (jj.eq.0) stop 'inconsistency in jj MPI Phase 1'
			jjj=jj+fltnum(6)*3
			allocate(btmp2(jjj),btmp3(jjj))
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndu-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-driver-bndu'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(6)) then
          		do iz=1,fltnum(6)
            		nodenumtemp=nx*ny*nz+fltu(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=brhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+1, 120000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+1, 120000+me+1,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndu-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+&
								btmp3(ntagMPI)
						endif
					enddo	
				enddo
			enddo
      		if (fltMPI(6)) then
          		do iz=1,fltnum(6)
            		nodenumtemp=nx*ny*nz+fltu(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			brhs(id1(locid(nodenumtemp)+izz))=brhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)				
		endif!bndu/=0
		!
  	endif!if npz>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
!	jj = er4mpi(2,1)-er4mpi(1,1)
!	rr = er4mpi(2,2)-er4mpi(1,2)
!	if (nprocs > 1) then
!		if (me == master .or. me == nprocs-1) then 
!			if (me == master) then
!				allocate(btmp(rr))
!				allocate(btmp1(rr))
!				do i = 1, rr 
!					btmp(i)=brhs(er4mpi(1,2)+i-1)
!				enddo
!				call mpi_sendrecv(btmp,  rr, MPI_DOUBLE_PRECISION, 1, 100000+me, &
!					btmp1, rr, MPI_DOUBLE_PRECISION, 1, 100000+me+1, &
!					MPI_COMM_WORLD, istatus, ierr)
!				do i=1, rr
!					brhs(er4mpi(1,2)+i-1) = brhs(er4mpi(1,2)+i-1) + btmp1(i)
!				enddo
!			endif
!			if (me == nprocs-1) then
!				allocate(btmp(jj))
!				allocate(btmp1(jj))
!				do i = 1, jj 
!					btmp(i)=brhs(er4mpi(1,1)+i-1)
!				enddo
!				call mpi_sendrecv(btmp, jj, MPI_DOUBLE_PRECISION, nprocs-2, 100000+me, &
!						btmp1, jj, MPI_DOUBLE_PRECISION, nprocs-2, 100000+me-1, &
!						MPI_COMM_WORLD, istatus, ierr)
!				do i=1, jj
!					brhs(er4mpi(1,1)+i-1) = brhs(er4mpi(1,1)+i-1) + btmp1(i)
!				enddo
!			endif
!		elseif (me > 0 .and. me < nprocs-1) then
!			allocate(btmp(rr))
!			allocate(btmp1(rr))
!			allocate(btmp2(jj))
!			allocate(btmp3(jj))
!			do i = 1,rr
!				btmp(i)=brhs(er4mpi(1,2)+i-1)
!			enddo
!			do i = 1,jj
!				btmp2(i)=brhs(er4mpi(1,1)+i-1)
!			enddo
!			call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me-1, 100000+me, &
!					btmp3, jj, MPI_DOUBLE_PRECISION, me-1, 100000+me-1,&
!					MPI_COMM_WORLD, istatus, ierr)
!			do i=1, jj
!				brhs(er4mpi(1,1)+i-1) = brhs(er4mpi(1,1)+i-1) + btmp3(i)
!			enddo
!			call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me+1, 100000+me, &
!					btmp1, rr, MPI_DOUBLE_PRECISION, me+1, 100000+me+1, &
!					MPI_COMM_WORLD, istatus, ierr)
!			do i=1, rr
!				brhs(er4mpi(1,2)+i-1) = brhs(er4mpi(1,2)+i-1) + btmp1(i)
!			enddo
!			deallocate(btmp2)
!			deallocate(btmp3)
!		endif
!		deallocate(btmp)
!		deallocate(btmp1)
!	endif
	time2 = MPI_WTIME()
	btime = btime + (time2 - time1)
	!
	!*** call faulting subrotuine to revise residual force ***
	!  This is main revision on general dynamic code to
	!  study earthquke rupture problems. The main purpose
	!  is to revise right-hand-side vector "brhs" by 
	!  faulting boundary. B.D. 7/3/05
	!
	!call mpi_barrier(MPI_COMM_WORLD, ierr)
	do i=1,ntotft
		if (nftnd(i)>0) then !only nonzero fault node, does faulting. B.D. 10/16/09
			time1 = MPI_WTIME()
			call faulting(i,nftnd(i),numnp,neq,lstr,lstr1,fnms,brhs,d,v,x,maxm,id1,locid,dof1,n4onf,&
					momnt,momntrat,maxslprat,fltsta,nsmp(1,1,i),fnft(1,i),fltslp(1,1,i),&
					un(1,1,i),us(1,1,i),ud(1,1,i),fric(1,1,i),arn(1,i),r4nuc(1,i),arn4m(1,i),&
					slp4fri(1,i),anonfs,itmp2, &
					me, master,nprocs,n,miuonf)
			time2 = MPI_WTIME()
			timeused(6) = timeused(6) + (time2 - time1) 
		endif
	enddo
!Implementation of Double-couple point source.
!Sep.12.2015/D.L.
if (C_dc==1)then
	do i=1,numnp
	!x positive. Adding a point force in y+ direction.
	if (x(1,i)==100.0.and.x(2,i)==0.0.and.x(3,i)==-2000.)then
		!write(*,*) 'right',i,me
		if (time<0.2)then
			brhs(id1(locid(i)+2))=brhs(id1(locid(i)+2))+&
			(1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
		else
			brhs(id1(locid(i)+2))=brhs(id1(locid(i)+2))+1.0e14/2
		endif
	endif
	!x negative. Adding a point force in y- direction.
	if (x(1,i)==-100.0.and.x(2,i)==0.0.and.x(3,i)==-2000.)then
		!write(*,*) 'left',i,me	
		if (time<0.2)then
			brhs(id1(locid(i)+2))=brhs(id1(locid(i)+2))-&
			(1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
		else
			brhs(id1(locid(i)+2))=brhs(id1(locid(i)+2))-1.0e14/2
		endif
	endif
	!y positive. Adding a point force in x+ direction.
	if (x(1,i)==0.0.and.x(2,i)==100.0.and.x(3,i)==-2000.)then
		!write(*,*) 'up',i,me	
		if (time<0.2)then
			brhs(id1(locid(i)+1))=brhs(id1(locid(i)+1))+&
			(1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
		else
			brhs(id1(locid(i)+1))=brhs(id1(locid(i)+1))+1.0e14/2
		endif
	endif
	!y negative. Adding a point force in x- direction.
	if (x(1,i)==0.0.and.x(2,i)==-100.0.and.x(3,i)==-2000.)then
		!write(*,*) 'down',i,me	
		if (time<0.2)then
			brhs(id1(locid(i)+1))=brhs(id1(locid(i)+1))-&
			(1.0e14/0.2/2*time-1.0e14/2/2/3.1415*sin(2*3.1415*time/0.2))
		else
			brhs(id1(locid(i)+1))=brhs(id1(locid(i)+1))-1.0e14/2
		endif
	endif
	enddo!Enddo double-couple point source.	
endif!ldc(logical double couple)	
	!*** solve uncoupled equation for explicit analysis ***
	!  After, brhs stores acceleration a. B.D. 7/4/05	
	!$omp parallel do default(shared) private(i)
	do i = 1,neq
		brhs(i) = brhs(i)/alhs(i)
	enddo
	!$omp end parallel do
	!...sum moment/moment rate, find max slip rate, judge term or not. B.D. 8/11/10
	! no early termination, cal moment at the end only. B.D. 5/29/12
	if(time > (nstep-1)*dt) then
		time1 = MPI_WTIME()
		call MPI_REDUCE(momnt,momntall,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
				MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(momntrat,momntratall,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
				MPI_COMM_WORLD,ierr)
		call MPI_REDUCE(maxslprat,maxslpratall,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
				MPI_COMM_WORLD,ierr)
		call MPI_BCAST(momntall,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(momntratall,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(maxslpratall,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
		!call mpi_barrier(MPI_COMM_WORLD, ierr)
		time2 = MPI_WTIME()
		btime = btime + (time2 - time1)
	!    if(momntratall<1.e13 .and. maxslpratall<0.001) then
	!      exit
	!    endif
	endif
!
enddo 	!end time step loop n
!   deallocate(btmp)
!   deallocate(btmp1)
!   deallocate(btmp2)
!   deallocate(btmp3)
if (me == master) then
	write(*,*) 'mpi_send + mpi_recv time:', btime
endif
end SUBROUTINE driver
