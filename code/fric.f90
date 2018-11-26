SUBROUTINE slip_weak(slip,fricsgl,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear slip-weakening
  ! friction law for fault dynamics. B.D. 8/19/06
  !...revised for ExGM 100runs. B.D. 8/10/10
  !...revised for SCEC TPV19. B.D. 1/8/12
  !  fricsgl(i,*),i=1 mus, 2 mud, 3 do, 4 cohesion, 
  !  5 time for fixed rutpure, 6 for pore pressure
  !
  real (kind=8) :: xmu,slip
  real (kind=8),dimension(6) :: fricsgl
  !
  if(slip == 0.0) then
    xmu = fricsgl(1)	!xmu is frictional coefficient, node by node on fault
  elseif(slip < fricsgl(3)) then
    xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*slip/fricsgl(3)
  endif
  if(slip >= fricsgl(3)) then
    xmu = fricsgl(2)
  endif
  !
end SUBROUTINE slip_weak

!================================================

SUBROUTINE time_weak(trupt,fricsgl,xmu)
  use globalvar
  implicit none
  !
  !### subroutine to implement linear time-weakening
  ! friction law for fault dynamics. B.D. 8/19/06
  !
  real (kind=8) :: xmu,trupt
  real (kind=8),dimension(2) :: fricsgl
  !
  if(trupt <= 0.0) then
    xmu = fricsgl(1)
  elseif(trupt < critt0) then
    xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*trupt/critt0
  else
    xmu = fricsgl(2)
  endif
  !
end SUBROUTINE time_weak

