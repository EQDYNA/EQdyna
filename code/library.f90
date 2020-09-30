! Table of Contents of functions and subroutines. 
! #
! # vlm
! -----------------------------------------------
subroutine fb1(x,ww,w,res)
  implicit none
  real (kind = 8) :: x, ww, w, res
  ! ww == W;
  if (abs(x)<=ww) then 
	res = 1.0d0
  elseif (abs(x)>ww.and.abs(x)<ww+w) then 
	res = 0.5d0*(1.0d0 + dtanh(w/(abs(x)-ww-w) + w/(abs(x)-ww)))
  elseif (abs(x)>=ww+w) then 
    res = 0.0d0
  endif
end subroutine fb1
! -----------------------------------------------
subroutine fb2(y,ww,w,res)
  implicit none
  real (kind = 8) :: y, ww, w, res
  ! ww == W;  
  if (y<0.0d0) then
	write(*,*) 'z coordinates should be positive for B3'
	stop
  endif   
  if (y<=w) then 
	res = 0.5d0*(1.0d0 + dtanh(w/(w-y) - w/y))
  elseif (y>=w.and.y<=ww) then 
	res = 1.0d0 
  elseif (y>ww.and.y<ww+w) then 
	res = 0.5d0*(1.0d0 + dtanh(w/(y-ww-w) + w/(y-ww)))
  elseif (y>=ww+w) then 
    res = 0.0d0
  endif
end subroutine fb2
!------------------------------------------------------
subroutine fb3(y,ww,w,res)
  implicit none
  real (kind = 8) :: y, ww, w, res
  ! ww == W;  
  if (y<0.0d0) then
	write(*,*) 'z coordinates should be positive for B3'
	stop
  endif 
  if (y<=ww) then 
	res = 1.0d0 
  elseif (y>ww.and.y<ww+w) then 
	res = 0.5d0*(1.0d0 + dtanh(w/(y-ww-w) + w/(y-ww)))
  elseif (y>=ww+w) then 
    res = 0.0d0
  endif
end subroutine fb3
! -----------------------------------------------
subroutine vlm(xl,volume)
  use globalvar
  implicit none
  !
  !...program to calculate volume of a hexlahedron from
  ! nodal coordinates. See Belytschko et al.(1984) for 
  ! reference.
  ! B.D. 8/21/05
  !
  integer(kind=4) :: i
  real(kind=8) :: volume
  real(kind=8),dimension(nesd,nen) :: xl
  integer(kind=4),dimension(8,8) :: it = reshape((/ &
    1,2,3,4,5,6,7,8, 2,3,4,1,6,7,8,5, 3,4,1,2,7,8,5,6, &
    4,1,2,3,8,5,6,7, 5,8,7,6,1,4,3,2, 6,5,8,7,2,1,4,3, &
    7,6,5,8,3,2,1,4, 8,7,6,5,4,3,2,1/),(/8,8/))
  real(kind=8),dimension(8) :: bb
  !
  do i=1,8
    bb(i) = xl(2,it(2,i))*(xl(3,it(6,i))-xl(3,it(3,i))+xl(3,it(5,i)) &
      -xl(3,it(4,i))) + xl(2,it(3,i))*(xl(3,it(2,i))-xl(3,it(4,i))) &
      + xl(2,it(4,i))*(xl(3,it(3,i))-xl(3,it(8,i))+xl(3,it(2,i)) &
      -xl(3,it(5,i))) + xl(2,it(5,i))*(xl(3,it(8,i))-xl(3,it(6,i)) &
      +xl(3,it(4,i))-xl(3,it(2,i))) + xl(2,it(6,i))*(xl(3,it(5,i)) &
      -xl(3,it(2,i))) + xl(2,it(8,i))*(xl(3,it(4,i))-xl(3,it(5,i)))
  enddo
  !
  volume = 0.0
  do i=1,nen
    volume = volume + xl(1,i) * bb(i)
  enddo
  volume = volume/12.0
  !  
end subroutine vlm
  
  
    
  
