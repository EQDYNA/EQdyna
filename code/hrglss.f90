SUBROUTINE hrglss(numel,numnp,neq,ien,d,v,rdampk,mat,lm,ss,phi,brhs, me,master,nprocs)
  use globalvar
  implicit none
  include 'mpif.h'
  !
  !### program to calculate hourglass resistence and to add to
  !	residual force vector for the 8-node hexahedral
  !	element and assemble into the global right-hand-side 
  !	vector. (one-gaussian quadrature).
  !  Notice about an important assumption: infinitesimal 
  !	deformation so that coordinates x() are not updated, 
  !	so that lots of variables can be clculated once and 
  !	save for time step loops to save time.
  !  B.D. 8/21/05
  !
  logical :: lcubic,zerodl 
  integer (kind=4) :: nel,i,j,k,m,is,itmp,numel,numnp,neq
  real (kind=8),dimension(ned) :: phid
  real (kind=8),dimension(ned,nen) :: dl,vl,fhr
  integer (kind=4),dimension(nen,numel) :: ien
  real (kind=8),dimension(ndof,numnp) :: d,v
  real (kind=8),dimension(neq) :: brhs
  real (kind=8),dimension(numat) :: rdampk
  integer (kind=4),dimension(numel) :: mat	    
  integer (kind=4),dimension(ned,nen,numel) :: lm
  !...hourglass control arrays
  real (kind=8),dimension(6,numel) :: ss
  real (kind=8),dimension(nen,4,numel) :: phi
  integer me, master, nprocs, rlp, rr, ierr,jj
 ! real (kind=8),dimension(nprocs) :: btmp
  !real (kind=8) :: btmp
  !integer (kind=4),dimension(nprocs) :: kk	    
  !integer istatus(MPI_STATUS_SIZE)

!  rlp = numel/nprocs
!  rr = numel-nprocs*rlp
!    if (me ==nprocs-1) then
!      jj=(me+1)*rlp+rr
!    else
!     jj=(me+1)*rlp
!    endif

 !
  !*** loop over elements ***
  !
!$omp parallel default(shared) private(nel,i,j,k,m,is,lcubic,itmp,dl,vl,&
!$omp	zerodl,fhr,phid)
!$omp do schedule (static)
  do nel=1,numel
!  do nel=me*rlp+1,jj
    m = mat(nel)
    !...if element degenerate. B.D. 11/20/08
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
    !...only hourglass control for hexahedron element
    if(lcubic) then    
      !...localize dl and account for Rayleigh damping
      do i=1,nen
        k=ien(i,nel)
        do j=1,ned
          dl(j,i) = d(j,k)
    	  vl(j,i) = v(j,k)
	  dl(j,i) = dl(j,i) + rdampk(m) * vl(j,i)
        enddo
      enddo
      !...determine if hourglass calculation is really needed
      zerodl = .true.
      outer: do j=1,nen
        inner: do i=1,ned
!		if(abs(dl(i,j)) > 1.e-12) then
                if(dl(i,j) /= 0.0) then
		  zerodl = .false.
		  exit outer
		endif
	    enddo inner
	  enddo outer
	  if(.not.zerodl) then   
      !
      !...4 sets of hourglass operators
     ! do is=1,4 ! unroll the loop four times
        fhr = 0.0	!initialize 
        !... calculate sum(phi*dl)
        do i=1,ned
          phid(i) = 0.0
          do j=1,nen
            phid(i) = phid(i) + phi(j,1,nel) * dl(i,j)
          enddo
        enddo
        !... calculate hourglass resistence
        do i=1,nen
          fhr(1,i) = phi(i,1,nel) * (ss(1,nel)*phid(1)  &
	    + ss(2,nel)*phid(2) + ss(3,nel)*phid(3))
	  fhr(2,i) = phi(i,1,nel) * (ss(2,nel)*phid(1)  &
	    + ss(4,nel)*phid(2) + ss(5,nel)*phid(3))
	  fhr(3,i) = phi(i,1,nel) * (ss(3,nel)*phid(1)  &
	    + ss(5,nel)*phid(2) + ss(6,nel)*phid(3))
        enddo
        !... assemble to global right-hand-side force vector
        do i=1,nen
          do j=1,ned
            k=lm(j,i,nel)
  	    if(k > 0) then
	      brhs(k) = brhs(k) - fhr(j,i)
            endif
          enddo
        enddo
        !
        fhr = 0.0	!initialize 
        !... calculate sum(phi*dl)
        do i=1,ned
          phid(i) = 0.0
          do j=1,nen
            phid(i) = phid(i) + phi(j,2,nel) * dl(i,j)
          enddo
        enddo
        !... calculate hourglass resistence
        do i=1,nen
          fhr(1,i) = phi(i,2,nel) * (ss(1,nel)*phid(1)  &
  	    + ss(2,nel)*phid(2) + ss(3,nel)*phid(3))
	  fhr(2,i) = phi(i,2,nel) * (ss(2,nel)*phid(1)  &
	    + ss(4,nel)*phid(2) + ss(5,nel)*phid(3))
	  fhr(3,i) = phi(i,2,nel) * (ss(3,nel)*phid(1)  &
	    + ss(5,nel)*phid(2) + ss(6,nel)*phid(3))
        enddo
        !... assemble to global right-hand-side force vector
        do i=1,nen
          do j=1,ned
            k=lm(j,i,nel)
      	    if(k > 0) then
	      brhs(k) = brhs(k) - fhr(j,i)
            endif
          enddo
        enddo
        !
        fhr = 0.0	!initialize 
        !... calculate sum(phi*dl)
        do i=1,ned
          phid(i) = 0.0
          do j=1,nen
            phid(i) = phid(i) + phi(j,3,nel) * dl(i,j)
          enddo
        enddo
        !... calculate hourglass resistence
        do i=1,nen
          fhr(1,i) = phi(i,3,nel) * (ss(1,nel)*phid(1)  &
	    + ss(2,nel)*phid(2) + ss(3,nel)*phid(3))
	  fhr(2,i) = phi(i,3,nel) * (ss(2,nel)*phid(1)  &
	    + ss(4,nel)*phid(2) + ss(5,nel)*phid(3))
	  fhr(3,i) = phi(i,3,nel) * (ss(3,nel)*phid(1)  &
	    + ss(5,nel)*phid(2) + ss(6,nel)*phid(3))
        enddo
        !... assemble to global right-hand-side force vector
        do i=1,nen
          do j=1,ned
            k=lm(j,i,nel)
  	    if(k > 0) then
	      brhs(k) = brhs(k) - fhr(j,i)
            endif
          enddo
        enddo
        !
        fhr = 0.0	!initialize 
        !... calculate sum(phi*dl)
        do i=1,ned
          phid(i) = 0.0
          do j=1,nen
            phid(i) = phid(i) + phi(j,4,nel) * dl(i,j)
          enddo
        enddo
        !... calculate hourglass resistence
        do i=1,nen
          fhr(1,i) = phi(i,4,nel) * (ss(1,nel)*phid(1)  &
    	    + ss(2,nel)*phid(2) + ss(3,nel)*phid(3))
	  fhr(2,i) = phi(i,4,nel) * (ss(2,nel)*phid(1)  &
	    + ss(4,nel)*phid(2) + ss(5,nel)*phid(3))
	  fhr(3,i) = phi(i,4,nel) * (ss(3,nel)*phid(1)  &
	    + ss(5,nel)*phid(2) + ss(6,nel)*phid(3))
        enddo
        !... assemble to global right-hand-side force vector
        do i=1,nen
          do j=1,ned
            k=lm(j,i,nel)
  	    if(k > 0) then
	      brhs(k) = brhs(k) - fhr(j,i)
!              if (nprocs > 1) then
!	      call mpi_allgather(k, 1, MPI_INTEGER, &
!                 kk, 1, MPI_INTEGER, &
!                 MPI_COMM_WORLD, ierr)
!	      do m=1,nprocs-1
!		do is=m+1,nprocs
!		  if (kk(m) == kk(is)) then
!			if (me == m-1) then
!			  call mpi_send(brhs(k), 1, MPI_DOUBLE_PRECISION, &
!		                is-1, 1000, MPI_COMM_WORLD, ierr)
!			  call mpi_recv(btmp, 1, MPI_DOUBLE_PRECISION, &
!               			is-1, 1001, MPI_COMM_WORLD, istatus, ierr)
!			  brhs(k) = btmp
!			else
!			  if (me == is-1) then
!			    call mpi_recv(btmp, 1, MPI_DOUBLE_PRECISION, &
!                		m-1, 1000, MPI_COMM_WORLD, istatus, ierr)
!!			    brhs(k) = brhs(k)+btmp
!			    call mpi_send(brhs(k), 1, MPI_DOUBLE_PRECISION, &
!                 		m-1, 1001, MPI_COMM_WORLD, ierr)
!			  endif
!			endif
!		  endif
!		enddo
!	     enddo
!	      call mpi_allgather(brhs(k), 1, MPI_DOUBLE_PRECISION, &
!                 btmp, 1, MPI_DOUBLE_PRECISION, &
!                 MPI_COMM_WORLD, ierr)
! 	      do m=1,nprocs
!		brhs(kk(m))=btmp(m)
!	      enddo 
!             endif
            endif
          enddo
        enddo
        !
           endif	!if(.not.zerodl)
      !
      endif	!if(lcubic)
      !
!    enddo	!is unroll the loop
  !
  enddo
!$omp end do nowait
!$omp end parallel
  !
end SUBROUTINE hrglss
       
  
  
