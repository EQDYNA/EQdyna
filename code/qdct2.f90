SUBROUTINE qdct2(numel,numnp,neq,shl,ien,x,mat,rdampm,rho,e,lm,alhs, &
	eledet,elemass,eleshp,fnms,ss,phi,er4mpi,nr4mpi,master,me,nprocs)
  use globalvar
  implicit none
  include 'mpif.h'
  !
  !...program to calcuate lumped mass matrix for 8-node hexahedral
  !	(brick),elastic continue element and assemble.
  !
  logical :: lcubic
  integer(kind=4) :: nel,i,j,k,k1,m,itmp,itmp1,itmp2,&
  		numel,numnp,neq
  real(kind=8) :: det,constm,vol,ce,co
  real(kind=8),dimension(nee) :: eleffm
  real(kind=8),dimension(nesd,nen) :: xl
  real(kind=8),dimension(3,3) :: xs
  real(kind=8),dimension(nrowb,nee) :: b
  real(kind=8),dimension(nrowsh,nen) :: shg
  real(kind=8),dimension(8,4) :: ha = reshape((/ &
    1,1,-1,-1,-1,-1,1,1, 1,-1,-1,1,-1,1,1,-1, &
    1,-1,1,-1,1,-1,1,-1, -1,1,-1,1,1,-1,1,-1/), &
    (/8,4/))
   !...material properties
  real (kind=8),dimension(numat) :: e,rho,rdampm
  !...element arrays
  integer (kind=4),dimension(numel) :: mat	    
  integer (kind=4),dimension(nen,numel) :: ien
  integer (kind=4),dimension(ned,nen,numel) :: lm
  real (kind=8),dimension(nrowsh,nen) :: shl
  real (kind=8),dimension(numel) :: eledet
  real (kind=8),dimension(nee,numel) :: elemass
  real (kind=8),dimension(nrowsh-1,nen,numel) :: eleshp
  !... equations,node mass,weight,solutions,coordinate,shape
  real (kind=8),dimension(numnp) :: fnms
  real (kind=8),dimension(nsd,numnp) :: x
  real (kind=8),dimension(neq) :: alhs
  !...hourglass control arrays
  real (kind=8),dimension(6,numel) :: ss
  real (kind=8),dimension(nen,4,numel) :: phi
  !...MPI
  integer  master, me, nprocs, ierr, rlp, rr, jj
  integer (kind=4),dimension(2,2)::er4mpi,nr4mpi !equation, node range for MPI
  integer istatus(MPI_STATUS_SIZE)
  real (kind=8), allocatable, dimension(:) :: btmp, btmp1, btmp2, btmp3
  !
!!$omp parallel do private(nel,i,j,k,k1,m, xl)
!!$omp parallel do default(shared) private(nel,i,j,k,k1,m, xl,eleffm,lcubic,&
!!$omp	itmp,det,shl,shg,xs,constm,vol,ce,co)
  do nel=1,numel
    !...initialize and localize
    eleffm = 0.0
    do i=1,nen
      k = ien(i,nel)
      do j=1,nesd
        xl(j,i) = x(j,k)
      enddo
    enddo
    m = mat(nel)	!material number
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
    call qdcshg(xl,det,shl,shg,nel,xs,lcubic)
    !...form mass matrix for left-hand-side
    constm = (1.0+rdampm(m)*0.5*dt)*rho(m)
    if(constm /= 0.0) then
      call contm(shg,det,eleffm,constm)
    endif
    !...assemble into left-hand-side
    do i=1,nen
      do j=1,ned
        k=lm(j,i,nel)
        k1=j+(i-1)*ned
        if(k>0) then
          alhs(k) = alhs(k) + eleffm(k1)
        endif
      enddo
    enddo
    !...following to compute and store global variables for
    !	"qdct3","faulting",and "hrglss" to use: assumption here
    !	is that nodal coordinates do not need to update (valid
    !	for infinitesmal deformation.
    !...determinant
    eledet(nel) = det
    !...derivatives of global shape function
    do i=1,nen
      do j=1,nrowsh-1
        eleshp(j,i,nel) = shg(j,i)
      enddo
    enddo
    !...mass matrix for right-hand-side to use.
    !	Note: different "constm" from above.
    constm = rho(m)
    eleffm = 0.0	!must initialize again here
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
    ce = 16.* e(m) / 15.	!close to plane-strain
    co = ce * vol / 48.
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
        vol = 0.0	!temp variable
	do k=1,nen
	  vol = vol + ha(k,i)*(xl(1,k)*shg(1,j) + xl(2,k)  &
	    *shg(2,j) + xl(3,k)*shg(3,j))
	enddo
	phi(j,i,nel) = ha(j,i) - vol
      enddo
      !...sum
      vol = 0.0
      do j=1,nen
        vol = vol + phi(j,i,nel)**2
      enddo
      vol = sqrt(vol/8.0)
      !...normalize to get phi from phi prime
      do j=1,nen
        phi(j,i,nel) = phi(j,i,nel) / vol
      enddo
      !
    enddo
    !
  enddo  !nel
!!$omp end parallel do
!$omp parallel do default(shared) private(i)
    do i=1,neq
      if(alhs(i)==0.0) then
   	 write(*,*) 'i=',i,'alhs=0'
   	 stop
      endif
    enddo
!$omp end parallel do
  
 !*******************MPI**********************
  !...assemble alhs() for common nodes between adjacent MPI procs.
  !   4/19/09
   jj = er4mpi(2,1)-er4mpi(1,1)
   rr = er4mpi(2,2)-er4mpi(1,2)
  if (nprocs > 1) then
   if (me == master .or. me == nprocs-1) then
     if (me == master) then
      allocate(btmp(rr))
      allocate(btmp1(rr))
      do i = 1, rr
          btmp(i)=alhs(er4mpi(1,2)+i-1)
      enddo
      call mpi_sendrecv(btmp,  rr, MPI_DOUBLE_PRECISION, 1, 200000+me, &
                        btmp1, rr, MPI_DOUBLE_PRECISION, 1, 200000+me+1, &
                        MPI_COMM_WORLD, istatus, ierr)
      do i=1, rr
           alhs(er4mpi(1,2)+i-1) = alhs(er4mpi(1,2)+i-1) + btmp1(i)
      enddo
     endif
     if (me == nprocs-1) then
      allocate(btmp(jj))
      allocate(btmp1(jj))
      do i = 1, jj
          btmp(i)=alhs(er4mpi(1,1)+i-1)
      enddo
      call mpi_sendrecv(btmp, jj, MPI_DOUBLE_PRECISION, nprocs-2, 200000+me, &
                        btmp1, jj, MPI_DOUBLE_PRECISION, nprocs-2, 200000+me-1, &
                        MPI_COMM_WORLD, istatus, ierr)
      do i=1, jj
           alhs(er4mpi(1,1)+i-1) = alhs(er4mpi(1,1)+i-1) + btmp1(i)
      enddo
     endif
   elseif (me > 0 .and. me < nprocs-1) then
     allocate(btmp(rr))
     allocate(btmp1(rr))
     allocate(btmp2(jj))
     allocate(btmp3(jj))
     do i = 1,rr
          btmp(i)=alhs(er4mpi(1,2)+i-1)
     enddo
     do i = 1,jj
          btmp2(i)=alhs(er4mpi(1,1)+i-1)
     enddo
     call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me-1, 200000+me, &
                            btmp3, jj, MPI_DOUBLE_PRECISION, me-1, 200000+me-1,&
                        MPI_COMM_WORLD, istatus, ierr)
     do i=1, jj
           alhs(er4mpi(1,1)+i-1) = alhs(er4mpi(1,1)+i-1) + btmp3(i)
     enddo
     call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me+1, 200000+me, &
                            btmp1, rr, MPI_DOUBLE_PRECISION, me+1, 200000+me+1, &
                MPI_COMM_WORLD, istatus, ierr)
     do i=1, rr
           alhs(er4mpi(1,2)+i-1) = alhs(er4mpi(1,2)+i-1) + btmp1(i)
      enddo
   deallocate(btmp2)
   deallocate(btmp3)
   endif
   deallocate(btmp)
   deallocate(btmp1)
  !...assemble fnms() for common nodes between adjacent MPI procs.
  !   4/19/09
   jj = nr4mpi(2,1)-nr4mpi(1,1)
   rr = nr4mpi(2,2)-nr4mpi(1,2)
   if (me == master .or. me == nprocs-1) then
     if (me == master) then
      allocate(btmp(rr))
      allocate(btmp1(rr))
      do i = 1, rr
          btmp(i)=fnms(nr4mpi(1,2)+i-1)
      enddo
      call mpi_sendrecv(btmp,  rr, MPI_DOUBLE_PRECISION, 1, 200000+me, &
                        btmp1, rr, MPI_DOUBLE_PRECISION, 1, 200000+me+1, &
                        MPI_COMM_WORLD, istatus, ierr)
      do i=1, rr
           fnms(nr4mpi(1,2)+i-1) = fnms(nr4mpi(1,2)+i-1) + btmp1(i)
      enddo
     endif
     if (me == nprocs-1) then
      allocate(btmp(jj))
      allocate(btmp1(jj))
      do i = 1, jj
          btmp(i)=fnms(nr4mpi(1,1)+i-1)
      enddo
      call mpi_sendrecv(btmp, jj, MPI_DOUBLE_PRECISION, nprocs-2, 200000+me, &
                        btmp1, jj, MPI_DOUBLE_PRECISION, nprocs-2, 200000+me-1, &
                        MPI_COMM_WORLD, istatus, ierr)
      do i=1, jj
           fnms(nr4mpi(1,1)+i-1) = fnms(nr4mpi(1,1)+i-1) + btmp1(i)
      enddo
     endif
   elseif (me > 0 .and. me < nprocs-1) then
     allocate(btmp(rr))
     allocate(btmp1(rr))
     allocate(btmp2(jj))
     allocate(btmp3(jj))
     do i = 1,rr
          btmp(i)=fnms(nr4mpi(1,2)+i-1)
     enddo
     do i = 1,jj
          btmp2(i)=fnms(nr4mpi(1,1)+i-1)
     enddo
     call mpi_sendrecv(btmp2, jj, MPI_DOUBLE_PRECISION, me-1, 200000+me, &
                            btmp3, jj, MPI_DOUBLE_PRECISION, me-1, 200000+me-1,&
                        MPI_COMM_WORLD, istatus, ierr)
     do i=1, jj
           fnms(nr4mpi(1,1)+i-1) = fnms(nr4mpi(1,1)+i-1) + btmp3(i)
     enddo
     call mpi_sendrecv(btmp, rr, MPI_DOUBLE_PRECISION, me+1, 200000+me, &
                            btmp1, rr, MPI_DOUBLE_PRECISION, me+1, 200000+me+1, &
                MPI_COMM_WORLD, istatus, ierr)
     do i=1, rr
           fnms(nr4mpi(1,2)+i-1) = fnms(nr4mpi(1,2)+i-1) + btmp1(i)
      enddo
   deallocate(btmp2)
   deallocate(btmp3)
   endif
   deallocate(btmp)
   deallocate(btmp1)
  endif 
  !
end SUBROUTINE qdct2
