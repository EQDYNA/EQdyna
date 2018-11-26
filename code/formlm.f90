SUBROUTINE formlm(numel,numnp,lm,id,ien)
  use globalvar
  implicit none
  integer (kind=4) :: numel,numnp,i,j,k,node
  integer (kind=4),dimension(ndof,numnp) :: id
  integer (kind=4),dimension(nen,numel) :: ien
  integer (kind=4),dimension(ned,nen,numel) :: lm
  !
  !### program to form lm array
  !
!$omp parallel do default(shared) private(k,j,node,i)
  do k=1,numel
    do j=1,nen
      node=ien(j,k)
      do i=1,ned
        lm(i,j,k) = id(i,node)
      enddo
    enddo
  enddo
!$omp end parallel do
end SUBROUTINE formlm
