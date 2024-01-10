subroutine qdcb(shg,b)
  use globalvar
  implicit none
  !
  !### program to set up the strain-displacement matrix "b" for
  !	3-D continuum elements.
  !
  integer (kind=4) :: i,j,k,l
  real (kind = dp) :: shg(nrowsh-1,nen), b(nrowb,nee)
  !
  !...initialize so that zero elements are taken care of
  b = 0.0d0
  !...loop over element nodes
  do i = 1, nen
    !...index
    j = 3 * (i - 1) + 1
    k = 3 * (i - 1) + 2
    l = 3 * (i - 1) + 3
    !...assign to one node: only non-zero 
    b(1,j) = shg(1,i)
    b(2,k) = shg(2,i)
    b(3,l) = shg(3,i)
    b(4,k) = shg(3,i)
    b(4,l) = shg(2,i)
    b(5,j) = shg(3,i)
    b(5,l) = shg(1,i)
    b(6,j) = shg(2,i)
    b(6,k) = shg(1,i)
  enddo
end subroutine qdcb
