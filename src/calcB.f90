subroutine calcB(globalShapeFunc,b)
  use globalvar
  implicit none
  ! calculate the strain-displacement matrix B.
  integer (kind=4) :: i,j,k,l
  real (kind = dp) :: globalShapeFunc(nrowsh-1,nen), b(nrowb,nee)
  !...initialize so that zero elements are taken care of
  b = 0.0d0
  !...loop over element nodes
  do i = 1, nen
    !...index
    j = 3 * (i - 1) + 1
    k = 3 * (i - 1) + 2
    l = 3 * (i - 1) + 3
    !...assign to one node: only non-zero 
    b(1,j) = globalShapeFunc(1,i)
    b(2,k) = globalShapeFunc(2,i)
    b(3,l) = globalShapeFunc(3,i)
    b(4,k) = globalShapeFunc(3,i)
    b(4,l) = globalShapeFunc(2,i)
    b(5,j) = globalShapeFunc(3,i)
    b(5,l) = globalShapeFunc(1,i)
    b(6,j) = globalShapeFunc(2,i)
    b(6,k) = globalShapeFunc(1,i)
  enddo
end subroutine calcB
