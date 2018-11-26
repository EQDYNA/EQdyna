SUBROUTINE prop3d(e,pois,c)
  use globalvar
  implicit none
  include 'mpif.h'
  !
  !...program to read and store properties for 3D continuum
  !	elements. See book pp.83
  !	B.D. 8/19/05
  !
  integer(kind=4) :: m
  real(kind=8) :: temp,temp1,temp2
   !...material properties
  real (kind=8),dimension(numat) :: e,pois
  real (kind=8),dimension(nrowc,nrowc,numat) :: c
  !...initialze c to avoid many o assignment. This should be done 
  !	outside of the loop. B.D. 1/29/07
  c = 0.0
  do m=1,numat
    !...temp variable to faciliate calculation
    temp = e(m)*(1-pois(m))/((1+pois(m))*(1-2*pois(m)))
    temp1 = temp * pois(m) / (1-pois(m))
    temp2 = temp * (1-2*pois(m)) / (2*(1-pois(m)))
    !NOTE: the following initialization should be outside of 
    ! the loop. This cause a big problem for when I try to 
    ! use more than 2 types of material in Jan 2007 for
    ! scec code validation! B.D. 1/29/07
    !...initialze c to avoid many o assignment
    !c = 0.0
    !...calculate non-zero top diagonal elements of C
    c(1,1,m) = temp
    c(1,2,m) = temp1
    c(1,3,m) = temp1
    c(2,2,m) = temp
    c(2,3,m) = temp1
    c(3,3,m) = temp
    c(4,4,m) = temp2
    c(5,5,m) = temp2
    c(6,6,m) = temp2
    !...use symmetry for other non-zero elements
    c(2,1,m) = c(1,2,m)
    c(3,1,m) = c(1,3,m)
    c(3,2,m) = c(2,3,m)
  enddo
end subroutine prop3d
