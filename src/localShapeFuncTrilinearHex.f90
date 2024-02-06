! Copyright(C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine localShapeFuncTrilinearHex
  use globalvar
  implicit none
  ! reduced order one-point Gaussian rule.

  integer(kind = 4) :: i,j
  real (kind = dp) :: cst = 1.0d0/8.0d0    !constant used intensively
  real (kind = dp), dimension(3,8) :: acoor = reshape((/ -1.0d0, -1.0d0, -1.0d0, &
                                                     1.0d0, -1.0d0, -1.0d0, &
                                                     1.0d0,  1.0d0, -1.0d0,  &
                                                    -1.0d0,  1.0d0, -1.0d0, &
                                                    -1.0d0, -1.0d0,  1.0d0, &
                                                     1.0d0, -1.0d0,  1.0d0,  &
                                                     1.0d0,  1.0d0,  1.0d0,   &
                                                    -1.0d0,  1.0d0,  1.0d0/),(/3,8/))

  w = 8.0d0    !weight for 1-point Gaussian rule
  do i = 1, 8
    shl(4,i) = cst    !the 4th is shape function itself
    do j = 1, 3
      shl(j,i) = cst * acoor(j,i)    !derivatives
    enddo
  enddo
end subroutine localShapeFuncTrilinearHex
