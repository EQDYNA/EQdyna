!/* Copyright (C) 2006-2020, Earthquake Modeling Lab @ Texas A&M University. 
! * All Rights Reserved.
! * This code is part of software EQdyna, please see EQdyna License Agreement
! * attached before you copy, download, install or use EQdyna./
subroutine qdcshl
  use globalvar
  implicit none
  !
  !...the local shape function subroutine for element group of
  !	8-node trilinear hexahedral (brick) element: 
  !	only one-point Gaussian rule is considered here. 
  ! 	B.D. 8/19/05
  !
  integer(kind=4) :: i,j
  real(kind = dp) :: cst = 1.0d0/8.0d0	!constant used intensively
  real(kind = dp), dimension(3,8) :: acoor = reshape((/ -1.0d0, -1.0d0, -1.0d0, &
													 1.0d0, -1.0d0, -1.0d0, &
													 1.0d0,  1.0d0, -1.0d0,  &
													-1.0d0,  1.0d0, -1.0d0, &
													-1.0d0, -1.0d0,  1.0d0, &
													 1.0d0, -1.0d0,  1.0d0,  &
													 1.0d0,  1.0d0,  1.0d0,   &
													-1.0d0,  1.0d0,  1.0d0/),(/3,8/))
  !
  w = 8.0d0	!weight for 1-point Gaussian rule
  !
  do i=1,8
    shl(4,i) = cst	!the 4th is shape function itself
    do j=1,3
      shl(j,i) = cst * acoor(j,i)	!derivatives
    enddo
  enddo
  !
end subroutine qdcshl
