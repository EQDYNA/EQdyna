subroutine qdcshl(shl)
  use globalvar
  implicit none
  !
  !...the local shape function subroutine for element group of
  !	8-node trilinear hexahedral (brick) element: 
  !	only one-point Gaussian rule is considered here. 
  ! 	B.D. 8/19/05
  !
  integer(kind=4) :: i,j
  real(kind=8) :: cst=1.0/8.0	!constant used intensively
  real(kind=8),dimension(3,8) :: acoor = reshape((/-1,-1,-1, &
    		1,-1,-1, 1,1,-1, -1,1,-1, -1,-1,1, 1,-1,1, &
		1,1,1, -1,1,1/),(/3,8/))
  real (kind=8),dimension(nrowsh,nen) :: shl
  !
  w=8.0	!weight for 1-point Gaussian rule
  !
  do i=1,8
    shl(4,i) = cst	!the 4th is shape function itself
    do j=1,3
      shl(j,i) = cst * acoor(j,i)	!derivatives
    enddo
  enddo
  !
end subroutine qdcshl
