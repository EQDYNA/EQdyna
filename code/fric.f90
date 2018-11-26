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
  real (kind=8),dimension(20) :: fricsgl
  !
  if(abs(slip).lt.1.0e-10) then
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
  real (kind=8),dimension(20) :: fricsgl
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

SUBROUTINE rate_state_ageing_law(V,theta,fricsgl,xmu,dxmudv)
  use globalvar
  implicit none
  !
  !### subroutine to implement rate- and state- 
  ! friction law for fault dynamics. Bin Luo 4/9/2014
  !
  real (kind=8) :: xmu, dxmudv
  real (kind=8) :: V,theta
  real (kind=8) :: A,B,L,f0,V0
  real (kind=8),dimension(20) :: fricsgl
  real (kind=8) :: tmp, tmpc
  !
  A  = fricsgl(9)
  B  = fricsgl(10)
  L  = fricsgl(11)
  f0 = fricsgl(13)
  V0 = fricsgl(12)

  tmpc = 1.0d0 / (2.0d0 * V0) * dexp((f0 + B * dlog(V0*theta/L)) / A)
  tmp = (V+1.d-30) * tmpc
  xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
  dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2) ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
  theta = L/V + (theta - L/V)*dexp(-V*dt/L)
  !
end SUBROUTINE rate_state_ageing_law

!================================================

SUBROUTINE rate_state_slip_law(V,psi,fricsgl,xmu,dxmudv)
  use globalvar
  implicit none
  !
  !### subroutine to implement rate- and state- 
  ! friction law for fault dynamics. Bin Luo 4/9/2014
  !
  real (kind=8) :: xmu, dxmudv
  real (kind=8) :: V,psi,psiss,fLV,fss
  real (kind=8) :: A,B,L,f0,V0,fw,Vw
  real (kind=8),dimension(20) :: fricsgl
  real (kind=8) :: tmp, tmpc
  !
  A  = fricsgl(9)
  B  = fricsgl(10)
  L  = fricsgl(11)
  f0 = fricsgl(13)
  V0 = fricsgl(12)
  fw = fricsgl(14)
  Vw = fricsgl(15)

  tmpc = 1.0d0 / (2.0d0 * V0) * dexp(psi/A)
  tmp = (V+1.d-30) * tmpc
  xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
  dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2)  ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
  fLV = f0 - (B - A) * dlog(V/V0)
  fss = fw + (fLV - fw) / ((1.0d0 + (V/Vw)**8)**0.125d0)
  psiss = A * dlog(2.0d0 * V0 / V * dsinh(fss/A))
  psi = psiss + (psi - psiss) * dexp(-V*dt/L)
  !
end SUBROUTINE rate_state_slip_law



