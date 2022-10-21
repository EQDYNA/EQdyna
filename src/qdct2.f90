SUBROUTINE qdct2
use globalvar
implicit none
include 'mpif.h'
!...program to calcuate lumped mass matrix for 8-node hexahedral
!	(brick),elastic continue element and assemble.
logical :: lcubic
integer(kind=4) :: nel,i,j,k,k1,m,itmp,itmp1,itmp2
real(kind=8) :: det,constm,vol,ce,co
real(kind=8),dimension(nee) :: eleffm
real(kind=8),dimension(nesd,nen) :: xl
real(kind=8),dimension(3,3) :: xs
real(kind=8),dimension(nrowb,nee) :: b
real(kind=8),dimension(nrowsh,nen) :: shg
integer(kind=4),dimension(8,4) :: ha = reshape((/ &
	1,1,-1,-1,-1,-1,1,1, 1,-1,-1,1,-1,1,1,-1, &
	1,-1,1,-1,1,-1,1,-1, -1,1,-1,1,1,-1,1,-1/), &
	(/8,4/))
!...MPI
integer  ierr, rlp, rr, jj
integer istatus(MPI_STATUS_SIZE)
real (kind=8), allocatable, dimension(:) :: btmp, btmp1, btmp2, btmp3
!DL variables for PML layer
integer (kind=4):: non,itag,eqn
!*******3D MPI partitioning*************
integer(kind=4)::mex,mey,mez,ix,iy,iz,nodenumtemp,ntagMPI,izz,nx,ny,nz
integer(kind=4)::bndl,bndr,bndf,bndb,bndd,bndu,rrr,jjj

do nel=1,numel
	!...initialize and localize
	eleffm = 0.0
	do i=1,nen
		k = ien(i,nel)
		do j=1,nesd
		xl(j,i) = x(j,k)
		enddo
	enddo
	m = 1	!material number
	!...if element degenerate. B.D. 11/27/08
	lcubic = .true.
	itmp = 0
	outloop: do i=2,nen
				do j=1,i-1
					if(ien(j,nel)==ien(i,nel)) itmp=itmp+1
						if(itmp>=2) then
							lcubic = .false.
						exit outloop
					endif
				enddo
			enddo outloop
	!...compute global shape function and derivatives at 1 Gaussian
	!	point.
	call qdcshg(xl,det,shg,nel,xs,lcubic)
	!...form mass matrix for left-hand-side
	constm = (1.0d0+rdampm*0.5*dt)*mat(nel,3)
	if(constm /= 0.0d0) then
		call contm(shg,det,eleffm,constm)
	endif
	!...assemble into left-hand-side
	
	! *.* D.L. Jan/23/2015. Including rayleigh damping into the left-hand-side when using SUBROUTINE contm
	do i=1,nen
			! do j=1,ned
				! k=lm(j,i,nel)
				! k1=j+(i-1)*ned
				! if(k>0) then
					! alhs(k) = alhs(k) + eleffm(k1)
				! endif
			! enddo
		! *.* Update alhs according to dof1. D.L. Jan/23/2015	
		non=ien(i,nel)
		if (dof1(non)==12) then
			do j=1,3
				itag=locid(non)+j
				eqn=id1(itag)
				if (eqn>0) then
					alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+1)
				endif
			enddo
			do j=4,6
				itag=locid(non)+j
				eqn=id1(itag)
				if (eqn>0) then
					alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+2)
				endif
			enddo	
			do j=7,9
				itag=locid(non)+j
				eqn=id1(itag)
				if (eqn>0) then
					alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+3)
				endif
			enddo	
			j=10
			itag=locid(non)+j
			eqn=id1(itag)
			if (eqn>0) then
				alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+1)
			endif
			j=11
			itag=locid(non)+j
			eqn=id1(itag)
			if (eqn>0) then
				alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+2)
			endif	
			j=12
			itag=locid(non)+j
			eqn=id1(itag)
			if (eqn>0) then
				alhs(eqn)=alhs(eqn)+eleffm(3*(i-1)+3)
			endif						
		elseif (dof1(non)==3) then
			do j=1,3
				itag=locid(non)+j
				eqn=id1(itag)
				alhs(eqn)=alhs(eqn)+eleffm((i-1)*3+j)
			enddo
		endif
		! *.* D.L.
	enddo

	eledet(nel) = det
	!...derivatives of global shape function
	do i=1,nen
		do j=1,nrowsh-1
			eleshp(j,i,nel) = shg(j,i)
		enddo
	enddo
	!...mass matrix for right-hand-side to use.
	!	Note: different "constm" from above.
	constm = mat(nel,3)
	eleffm = 0.0d0	!must initialize again here
	call contm(shg,det,eleffm,constm)
	do i=1,nee
		elemass(i,nel) = eleffm(i)
	enddo
	!...assemble above mass to nodal mass for "faulting".
	!	Note: ned degree of an element node have a same
	!	mass, so only one used here.
	do i=1,nen
		k = ien(i,nel)
		k1= 1 + (i-1) * ned
		fnms(k) = fnms(k) + eleffm(k1)
	enddo
	!...SS matrix for "hrglss"
	call vlm(xl,vol)
	ce=mat(nel,5)*(3*mat(nel,4)+2*mat(nel,5))/(mat(nel,4)+mat(nel,5))
	ce = 16.d0* ce / 15.d0	!close to plane-strain
	!E=miu*(3*lam+2*miu)/(lam+miu)
	co = ce * vol / 48.d0
	ss(1,nel)=co*(xs(1,1)**2+xs(2,1)**2+xs(3,1)**2)
	ss(2,nel)=co*(xs(1,1)*xs(1,2)+xs(2,1)*xs(2,2)+xs(3,1)*xs(3,2))
	ss(3,nel)=co*(xs(1,1)*xs(1,3)+xs(2,1)*xs(2,3)+xs(3,1)*xs(3,3))
	ss(4,nel)=co*(xs(1,2)**2+xs(2,2)**2+xs(3,2)**2)
	ss(5,nel)=co*(xs(1,2)*xs(1,3)+xs(2,2)*xs(2,3)+xs(3,2)*xs(3,3))
	ss(6,nel)=co*(xs(1,3)**2+xs(2,3)**2+xs(3,3)**2)
	!...Phi for "hrglss"
	! Note: there are 4 sets of Phi for one element in 3D.
	do i=1,4	!4 sets of phi
		!...phi prime
		do j=1,nen
			vol = 0.0d0	!temp variable
			do k=1,nen
				vol = vol + ha(k,i)*(xl(1,k)*shg(1,j) + xl(2,k)  &
				*shg(2,j) + xl(3,k)*shg(3,j))
			enddo
			phi(j,i,nel) = ha(j,i) - vol
		enddo
		!...sum
		vol = 0.0d0
		do j=1,nen
			vol = vol + phi(j,i,nel)**2
		enddo
		vol = sqrt(vol/8.0d0)
		!...normalize to get phi from phi prime
		do j=1,nen
			phi(j,i,nel) = phi(j,i,nel) / vol
		enddo
		!
	enddo
	!
enddo  !nel

!*******************MPI**********************
!...assemble alhs() for common nodes between adjacent MPI procs.
!   4/19/09

!***********************************************!
!********3D MPI partioning**********************!
!prepare for MPI partitioning
	mex=int(me/(npy*npz))
	mey=int((me-mex*npy*npz)/npz)
	mez=int(me-mex*npy*npz-mey*npz)
	nx=numcount(1)
	ny=numcount(2)
	nz=numcount(3)

do i=1,neq
	if(alhs(i)==0.0) then
		write(*,*) 'i=',i,'alhs=0'
		write(*,*) 'me,mex,mey,mez',me,mex,mey,mez
		stop
	endif
enddo

!write(*,*) 'me',me,'begin MPI partioning'
!*************************************************************
!************Partioning along x axis**************************
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
			!Check; if rr==0, bndl should ==0.
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
							btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-qdct2-bndl'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(1)) then
          		do ix=1,fltnum(1)
            		nodenumtemp=nx*ny*nz+fltl(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
			!send dest: me-npy*npz: envelop:20000+me 
			!recv source: me-npy*npz: envelop:20000+me-npy*npz 
			!left boundary of me  ---  right boundary of me-npy*npz
			!right boundary of me  ---  left boundary of me+npy*npz
			!front boundary of me  ---  back boundary of me-npz
			!back boundary of me  ---  front boundary of me+npz
			!lower boundary of me  ---  upper boundary of me-1
			!upper boundary of me  ---  lower boundary of me+1
			!me : 	x+ ==me+npy*npz
			!		x- ==me-npy*npz
			!		y+ ==me+npz
			!		y- ==me-npz	
			!		z+ ==me+1	
			!		z- ==me-1			
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-npy*npz, 200000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-npy*npz, 200000+me-npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndl-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)
		endif !bndl/=0
		!write(*,*) 'me',me,'Finished BNDL-alhs'		
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
							btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-qdct2-bndr'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(2)) then
          		do ix=1,fltnum(2)
            		nodenumtemp=nx*ny*nz+fltr(ix)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+npy*npz, 200000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+npy*npz, 200000+me+npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do iz=1,nz
				do iy=1,ny
					nodenumtemp=(bndr-1)*ny*nz+(iz-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)
		endif !bndr/=0
		!write(*,*) 'me',me,'Finished BNDR-alhs'	
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
							btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-qdct2-bndf'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(3)) then
          		do iy=1,fltnum(3)!using the free index, here iy.
            		nodenumtemp=nx*ny*nz+fltf(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-npz, 210000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-npz, 210000+me-npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndf
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)			
		endif!bndf/=0
		!write(*,*) 'me',me,'Finished BNDF-alhs'	
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
							btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-qdct2-bndl'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(4)) then
          		do iy=1,fltnum(4)
            		nodenumtemp=nx*ny*nz+fltb(iy)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+npz, 210000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+npz, 210000+me+npz,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iz=1,nz
					nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndb
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)			
		endif!bndb/=0
		!write(*,*) 'me',me,'Finished BNDB-alhs'
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
							btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (rr/=ntagMPI) then 
				stop 'rr&ntagMPI-qdct2-bndb-alhs'
				write(*,*) 'rr=',rr,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(5)) then
          		do iz=1,fltnum(5)!using the free index, here iz.
            		nodenumtemp=nx*ny*nz+fltd(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp, rrr, MPI_DOUBLE_PRECISION, me-1, 220000+me, &
				btmp1, rrr, MPI_DOUBLE_PRECISION, me-1, 220000+me-1,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndd-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then 
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp1(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp, btmp1)			
		endif!bndd/=0
		!write(*,*) 'me',me,'Finished BNDD-alhs'
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
							btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
						endif
					enddo	
				enddo
			enddo
			!Check
			if (jj/=ntagMPI) then 
				stop 'jj&ntagMPI-qdct2-bndu'
				write(*,*) 'jj=',jj,'ntagMPI=',ntagMPI
			endif
!
      		if (fltMPI(6)) then
          		do iz=1,fltnum(6)
            		nodenumtemp=nx*ny*nz+fltu(iz)
            		do izz=1,3
					ntagMPI=ntagMPI+1
	    			btmp2(ntagMPI)=alhs(id1(locid(nodenumtemp)+izz))
          			enddo
        		enddo
      		endif
     		call mpi_sendrecv(btmp2, jjj, MPI_DOUBLE_PRECISION, me+1, 220000+me, &
				btmp3, jjj, MPI_DOUBLE_PRECISION, me+1, 220000+me+1,&
				MPI_COMM_WORLD, istatus, ierr)
!Update
			ntagMPI=0
			do ix=1,nx
				do iy=1,ny
					nodenumtemp=(ix-1)*ny*nz+(bndu-1)*ny+iy
					do izz=1,dof1(nodenumtemp)
						if (id1(locid(nodenumtemp)+izz)>0) then
							ntagMPI=ntagMPI+1
							alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+&
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
	    			alhs(id1(locid(nodenumtemp)+izz))=alhs(id1(locid(nodenumtemp)+izz))+btmp3(ntagMPI)
          			enddo
        		enddo
      		endif
			deallocate(btmp2, btmp3)			
		endif!bndu/=0
		!write(*,*) 'me',me,'Finished BNDU-alhs'
  	endif!if npz>1
	call mpi_barrier(MPI_COMM_WORLD, ierr)
!===============================3DMPI===============================!
!=====================Partitioning along x axis=====================!
if (npx>1) then
	rr=ny*nz+fltnum(1)
	jj=ny*nz+fltnum(2)
   	bndl=1  !-x boundary
   	bndr=nx !+x boundary 
   	if (mex == master) then
     	bndl=0
   	elseif (mex == npx-1) then
     	bndr=0
   	endif	
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

	if (bndl/=0) then
		ntagMPI=0
		do iz=1,nz
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(bndl-1)*ny*nz+(iz-1)*ny+iy
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(1)) then
			do ix=1,fltnum(1)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltl(ix)
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-npy*npz, 300000+me, &
				btmp1, rr, MPI_DOUBLE_PRECISION, me-npy*npz, 300000+me-npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do iz=1,nz
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(bndl-1)*ny*nz+(iz-1)*ny+iy
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)				
			enddo
		enddo
		if (fltMPI(1)) then
			do ix=1,fltnum(1)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltl(ix)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
			enddo
		endif
	endif !bndl/=0
	!write(*,*) 'me',me,'Finished BNDL-fnms'
	if (bndr/=0)then
		ntagMPI=0
		do iz=1,nz
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(bndr-1)*ny*nz+(iz-1)*ny+iy
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(2)) then
			do ix=1,fltnum(2)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltr(ix)
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+npy*npz, 300000+me, &
				btmp3, jj, MPI_DOUBLE_PRECISION, me+npy*npz, 300000+me+npy*npz,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do iz=1,nz
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(bndr-1)*ny*nz+(iz-1)*ny+iy
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)				
			enddo
		enddo
		if (fltMPI(2)) then
			do ix=1,fltnum(2)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltr(ix)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
			enddo
		endif		
	endif !bndr/=0
	!write(*,*) 'me',me,'Finished BNDR-fnms'	
	deallocate(btmp, btmp1, btmp2, btmp3)
endif !npx>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along y axis=====================!
if (npy>1) then
	rr=nx*nz+fltnum(3)
	jj=nx*nz+fltnum(4)
   	bndf=1  !-y boundary
   	bndb=ny !+y boundary 
   	if (mey == master) then
     	bndf=0
   	elseif (mey == npy-1) then
     	bndb=0
   	endif	
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

	if (bndf/=0) then
		ntagMPI=0
		do ix=1,nx
			do iz=1,nz
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndf
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(3)) then
			do iy=1,fltnum(3)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltf(iy)
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-npz, 310000+me, &
				btmp1, rr, MPI_DOUBLE_PRECISION, me-npz, 310000+me-npz,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do ix=1,nx
			do iz=1,nz
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndf
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)				
			enddo
		enddo
		if (fltMPI(3)) then
			do iy=1,fltnum(3)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltf(iy)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
			enddo
		endif
	endif !bndf/=0
	!write(*,*) 'me',me,'Finished BNDF-fnms'
	if (bndb/=0)then
		ntagMPI=0
		do ix=1,nx
			do iz=1,nz
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndb
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(4)) then
			do iy=1,fltnum(4)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltb(iy)
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+npz, 310000+me, &
				btmp3, jj, MPI_DOUBLE_PRECISION, me+npz, 310000+me+npz,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do ix=1,nx
			do iz=1,nz
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(iz-1)*ny+bndb
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)				
			enddo
		enddo
		if (fltMPI(4)) then
			do iy=1,fltnum(4)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltb(iy)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
			enddo
		endif		
	endif !bndb/=0
	!write(*,*) 'me',me,'Finished BNDB-fnms'
	deallocate(btmp, btmp1, btmp2, btmp3)
endif !npy>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
!=====================Partitioning along z axis=====================!
if (npz>1) then
	rr=nx*ny+fltnum(5)
	jj=nx*ny+fltnum(6)
   	bndd=1  !-z boundary
   	bndu=nz !+z boundary 
   	if (mez == master) then
     	bndd=0
   	elseif (mez == npz-1) then
     	bndu=0
   	endif	
    allocate(btmp(rr), btmp1(rr), btmp2(jj), btmp3(jj))

	if (bndd/=0) then
		ntagMPI=0
		do ix=1,nx
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(bndd-1)*ny+iy
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(5)) then
			do iz=1,fltnum(5)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltd(iz)
				btmp(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me-1, 320000+me, &
				btmp1, rr, MPI_DOUBLE_PRECISION, me-1, 320000+me-1,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do ix=1,nx
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(bndd-1)*ny+iy
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)				
			enddo
		enddo
		if (fltMPI(5)) then
			do iz=1,fltnum(5)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltd(iz)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp1(ntagMPI)
			enddo
		endif
	endif !bndd/=0
	!write(*,*) 'me',me,'Finished BNDD-fnms'
	if (bndu/=0)then
		ntagMPI=0
		do ix=1,nx
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(bndu-1)*ny+iy
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		enddo
		if (fltMPI(6)) then
			do iz=1,fltnum(6)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltu(iz)
				btmp2(ntagMPI)=fnms(nodenumtemp)
			enddo
		endif
     	call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me+1, 320000+me, &
				btmp3, jj, MPI_DOUBLE_PRECISION, me+1, 320000+me+1,&
				MPI_COMM_WORLD, istatus, ierr)
		ntagMPI=0
		do ix=1,nx
			do iy=1,ny
				ntagMPI=ntagMPI+1
				nodenumtemp=(ix-1)*ny*nz+(bndu-1)*ny+iy
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)				
			enddo
		enddo
		if (fltMPI(6)) then
			do iz=1,fltnum(6)
				ntagMPI=ntagMPI+1
				nodenumtemp=nx*ny*nz+fltu(iz)
				fnms(nodenumtemp)=fnms(nodenumtemp)+btmp3(ntagMPI)
			enddo
		endif		
	endif !bndu/=0
	!write(*,*) 'me',me,'Finished BNDU-fnms'	
	deallocate(btmp, btmp1, btmp2, btmp3)
endif !npz>1
call mpi_barrier(MPI_COMM_WORLD, ierr)
end SUBROUTINE qdct2
