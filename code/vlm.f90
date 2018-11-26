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
  
  
    
  
