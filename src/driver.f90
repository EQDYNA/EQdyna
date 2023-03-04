!/* Copyright (C) 2006-2020, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
subroutine driver

use globalvar
implicit none
include 'mpif.h'

integer (kind = 4) ::ntstep=0,i,j,k,k1,l,m,ierr,rr,jj,izz,ntagMPI,ix,iy,iz,nx,ny,nz,mex,mey,mez,nodenumtemp,&
					bndl,bndr,bndf,bndb,bndd,bndu,rrr,jjj,istatus(MPI_STATUS_SIZE),itag,eqn
real (kind = dp) :: dampv(9)
real (kind = dp), allocatable,dimension(:) :: btmp,btmp1,btmp2,btmp3

! Tatnode = fric_tp_Tini
! patnode = fric_tp_pini

time1       = MPI_WTIME()
call qdct2

time2       = MPI_WTIME()
timeused(2) = timeused(2)+(time2-time1)

! Initiate on-fault node velocities.
call init_vel

do nt=1,nstep

	time    = time + dt
	ntstep  = ntstep + 1 
	
	if (mod(nt,100) == 1 .and. me == master) then
		write(*,*) '=                                                                   ='
		write(*,*) '=     Current time in dynamic rupture                               ='
		write(*,'(X,A,40X,f7.3,4X,A)') '=',  time  , 's'
	endif
	
	time1   = MPI_WTIME()
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
				call comdampv(x(1,i),x(2,i),x(3,i),dampv)
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
				v(1,i)=0.0d0
				v(2,i)=0.0d0
				v(3,i)=0.0d0
				d(1,i)=0.0d0
				d(2,i)=0.0d0
				d(3,i)=0.0d0				
			endif	
		endif
		if ((v(1,i)/=v(1,i)).or.v(2,i)/=v(2,i).or.v(3,i)/=v(3,i)) then 
			write(*,*) x(1,i),x(2,i),x(3,i),me,mex,mey,mez, 'nt=',nt
			stop 'NAN'
		endif
	enddo
!-------------------------------------------------------------------!	
	time2       = MPI_WTIME()
	timeused(3) = timeused(3)+(time2-time1)
	 
	!*** store desired results at set time intervals ***
	if (mod(nt,nhplt) == 0) then	
		lstr    = .true.	
		locplt  = locplt+ 1	!when nt=1, locplt=2 due to 1 in eqdy3d.f90
	else
		lstr    = .false.
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

	brhs=0.0d0
	
	time1=MPI_WTIME()
	
	call qdct3
				
	time2 = MPI_WTIME()
	timeused(4)=timeused(4)+(time2-time1)
	
	time1 = MPI_WTIME()
	
	call hrglss
	
	time2 = MPI_WTIME()
	timeused(5)=timeused(5)+(time2-time1)

!===============================3DMPI===============================!
!=====================Partitioning along x axis=====================!
 	time1 = MPI_WTIME() 	
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
	time1=MPI_WTIME()
	
	if (friclaw == 5) then 
		call thermop
	endif	
	!write(*,*) 'after thermop, me', me
	call faulting
	!write(*,*) 'after faulting, me', me	
		
	time2=MPI_WTIME()
	timeused(6)=timeused(6)+(time2-time1) 
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
	time1=MPI_WTIME()
	!!$omp parallel do default(shared) private(i)
	do i=1,neq
		brhs(i)=brhs(i)/alhs(i)
	enddo
	!!$omp end parallel do
	time2=MPI_WTIME()
	timeused(7)=timeused(7)+(time2-time1) 	
	
	if ((mod(nt,10) == 1) .and. (outputGroundMotion == 1)) call output_gm
	
enddo 	!end time step loop nt
if (me==master) then
	write(*,*) 'mpi_send + mpi_recv time:',btime
endif
end SUBROUTINE driver

subroutine init_vel
	! initiate the 1d velocity array v1. 
	! if mode==2, non-zero values for fric(31-36,i,ift) loaded from 
	!    the restart file.
	! if mode==1, fric(31-36,i) will be zeros.
    use globalvar
    implicit none
    integer (kind = 4) :: i, ift, tmp
    
    do ift = 1, ntotft
        do i = 1,nftnd(ift)
			tmp = locid(nsmp(1,i,ift))! slave nodeid i 
			v1(tmp+1) = fric(34,i,ift) ! vxs
			v1(tmp+2) = fric(35,i,ift) ! vys
			v1(tmp+3) = fric(36,i,ift) ! vzs
			tmp = locid(nsmp(2,i,ift))! master nodeid i
			v1(tmp+1) = fric(31,i,ift) ! vxm
			v1(tmp+2) = fric(32,i,ift) ! vym
			v1(tmp+3) = fric(33,i,ift) ! vzm
		enddo
	enddo
end subroutine init_vel