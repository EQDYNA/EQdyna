SUBROUTINE driver(numel,numnp,neq,nftnd,ndout,x,brhs,d,v,mat,ien,eleporep,pstrain,&
				id1,maxm,locid,dof1,et,v1,d1,PMLb,maxs,ids,s1,shl,n4onf,nsmp,fnft,fltslp,un,&
				us,ud,fric,arn,r4nuc,anonfs,arn4m,slp4fri,fltsta,me,nsurnd,surid,miuonf)
use globalvar
implicit none
include 'mpif.h'
logical::lstr
character(len=30)::foutmov,mm
integer(kind=4)::ntstep=0,n,i,j,k,k1,l,m,numel,numnp,neq,ndout,n4onf,&
	nftnd(ntotft),ien(nen,numel),anonfs(3,nonmx),nsmp(2,nftmx,ntotft),&
	id1(maxm),ids(numel),locid(numnp),dof1(numnp),et(numel),surid(nsurnd)
integer(kind=4)::me,ierr,rr,jj,izz,ntagMPI,ix,iy,iz,nx,ny,nz,mex,mey,mez,nodenumtemp,&
	bndl,bndr,bndf,bndb,bndd,bndu,rrr,jjj,istatus(MPI_STATUS_SIZE),maxm,maxs,itag,eqn,nsurnd
real(kind=8)::PMLb(8),dampv(9)
real(kind=8)::x(nsd,numnp),alhs(neq),brhs(neq),d(ndof,numnp),v(ndof,numnp),mat(numel,5),&
	shl(nrowsh,nen),eledet(numel),eleporep(numel),pstrain(numel),elemass(nee,numel),&
	eleshp(nrowsh-1,nen,numel),fric(8,nftmx,ntotft),miuonf(nftmx),fltsta(10,nplpts-1,n4onf),&
	fnms(numnp),ss(6,numel),phi(nen,4,numel),v1(neq),d1(neq),s1(maxs),maxsur(nsurnd,6)
real(kind=8),dimension(nftmx,ntotft)::fnft,arn,r4nuc,arn4m,slp4fri
real(kind=8),dimension(3,nftmx,ntotft)::un,us,ud,fltslp
real(kind=8),allocatable,dimension(:)::btmp,btmp1,btmp2,btmp3
!-------------------------------------------------------------------!
if (me == master) then
!!$  write(*,*) 'Number of OpenMP threads used:', omp_get_max_threads()
!!$  write(*,*) 'Number of processors available:', omp_get_num_procs()
	write(*,*) 'Number of MPI processes used:', nprocs
	write(*,*) 'Please be patient! Program is running on Process...', me
endif
time1=MPI_WTIME()
alhs=0.0
call qdct2(numel,numnp,neq,shl,ien,x,mat,alhs,eledet,elemass,eleshp,fnms,ss,phi,me,maxm,id1,locid,dof1)
time2 = MPI_WTIME()
timeused(2)=timeused(2)+(time2-time1)
!...allocate/deallicate outside time loop to improve performance.
! found by TAMU supercomputing facility colleagues. B.D. 1/12/12
maxsur=0.0
do n=1,nstep
	time=time+dt
	ntstep=ntstep+1 
	!...print on screen for monitoring
	if(mod(n,nhshw)==0) then
		if (me == master) then
			write(*,'( a6,f9.3)') 'time=',time
		endif
	endif
	time1=MPI_WTIME()
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
		elseif (dof1(i)==12) then
			itag=locid(i)+1
			eqn=id1(itag)
			if (eqn>0) then
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
	time2=MPI_WTIME()
	timeused(3)=timeused(3)+(time2-time1)
	 
	!*** store desired results at set time intervals ***
	if (mod(n,nhplt) == 0) then	
		lstr=.true.	
		locplt=locplt+ 1	!when n=1, locplt=2 due to 1 in eqdy3d.f90
	else
		lstr=.false.
	endif
	if (lstr) then
		if((ndout>0).and.(locplt>1)) then
			dout(1,locplt)=time
			do i=1,ndout
				j=idhist(1,i)
				if(j<=0) j=1  !avoid zero that cannot be used below
					k=idhist(2,i)
					l=idhist(3,i)
				if(l==1) then
					dout(i+1,locplt)=d(k,j)
				elseif(l==2) then
					dout(i+1,locplt)=v(k,j)
				elseif(l==3) then
					k1=id1(locid(j)+k)
					dout(i+1,locplt)=brhs(k1)
				endif
			enddo
		endif
	endif

	brhs=0.0!initialize it every interation
	time1=MPI_WTIME()
	call qdct3(numel,numnp,neq,mat,ien,d,v,eleporep,elemass,eleshp,eledet,pstrain,brhs,& 
				me,maxm,id1,locid,dof1,et,v1,d1,PMLb,x,maxs,ids,s1,n)
	! if (me==30) then 
		! write(*,*) '11slave',brhs(id1(604363)+1),brhs(id1(604363)+2),brhs(id1(604363)+3)
		! write(*,*) 'master',brhs(id1(1276251)+1),brhs(id1(1276251)+2),brhs(id1(1276251)+3)
	! endif				
	time2 = MPI_WTIME()
	timeused(4)=timeused(4)+(time2-time1)
	time1 = MPI_WTIME()
	call hrglss(numel,numnp,neq,ien,d,v,mat,ss,phi,brhs,me,maxm,id1,locid,dof1,et,eledet)
	time2 = MPI_WTIME()
	timeused(5)=timeused(5)+(time2-time1)
	time1 = MPI_WTIME()
!===============================3DMPI===============================!
!=====================Partitioning along x axis=====================!
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
	!call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along y axis=====================!
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
	!call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along z axis=====================!
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
	!call mpi_barrier(MPI_COMM_WORLD, ierr)
	time2 = MPI_WTIME()
	btime=btime+(time2-time1)
	!  Faulting TSN. B.Duan 2005/07/03
	! if (me==30) then 
		! write(*,*) 'slave',brhs(id1(604363)+1),brhs(id1(604363)+2),brhs(id1(604363)+3)
		! write(*,*) 'master',brhs(id1(1276251)+1),brhs(id1(1276251)+2),brhs(id1(1276251)+3)
	! endif
	do i=1,ntotft
		if (nftnd(i)>0) then !only nonzero fault node, does faulting. B.D. 10/16/09
			time1=MPI_WTIME()
			call faulting(i,nftnd(i),numnp,neq,lstr,fnms,brhs,d,v,x,maxm,id1,locid,dof1,n4onf,&
					fltsta,nsmp(1,1,i),fnft(1,i),fltslp(1,1,i),&
					un(1,1,i),us(1,1,i),ud(1,1,i),fric(1,1,i),arn(1,i),r4nuc(1,i),arn4m(1,i),&
					slp4fri(1,i),anonfs,nonmx,me,n,miuonf)
			time2=MPI_WTIME()
			timeused(6)=timeused(6)+(time2-time1) 
		endif
	enddo
	!Implementation of Double-couple point source.Sep.12.2015/D.L.
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
	!$omp parallel do default(shared) private(i)
	do i=1,neq
		brhs(i)=brhs(i)/alhs(i)
	enddo
	!$omp end parallel do
enddo 	!end time step loop n
if (me==master) then
	write(*,*) 'mpi_send + mpi_recv time:',btime
endif
end SUBROUTINE driver
