! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine slip_weak(slip,fricsgl,xmu)
  use globalvar
  implicit none

  real (kind = dp) :: xmu,slip
  real (kind = dp),dimension(20) :: fricsgl
  !
  if(abs(slip).lt.1.0e-10) then
    xmu = fricsgl(1)    !xmu is frictional coefficient, node by node on fault
  elseif(slip < fricsgl(3)) then
    xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*slip/fricsgl(3)
  endif
  if(slip >= fricsgl(3)) then
    xmu = fricsgl(2)
  endif
  !
end subroutine slip_weak

subroutine time_weak(trupt,fricsgl,xmu)

    use globalvar
    implicit none

    real (kind = dp) :: xmu,trupt
    real (kind = dp),dimension(20) :: fricsgl
    
    if(trupt <= 0.0d0) then
        xmu = fricsgl(1)
    elseif(trupt < fricsgl(5)) then
        xmu = fricsgl(1) - (fricsgl(1) - fricsgl(2))*trupt/fricsgl(5)
    else
        xmu = fricsgl(2)
    endif

end subroutine time_weak

subroutine rate_state_ageing_law(V2,theta,fricsgl,xmu,dxmudv)
  use globalvar
  implicit none

  real (kind = dp) :: xmu, dxmudv
  real (kind = dp) :: V2,theta
  real (kind = dp) :: A,B,L,f0,V0
  real (kind = dp),dimension(100) :: fricsgl
  real (kind = dp) :: tmp, tmpc
  !
  A  = fricsgl(9)
  B  = fricsgl(10)
  L  = fricsgl(11)
  f0 = fricsgl(13)
  V0 = fricsgl(12)

  tmpc = 1.0d0 / (2.0d0 * V0) * dexp((f0 + B * dlog(V0*theta/L)) / A)
  tmp = (V2+1.d-30) * tmpc
  xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
  dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2) ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
  theta = L/V2 + (theta - L/V2)*dexp(-V2*dt/L)
  !
end subroutine rate_state_ageing_law

subroutine rate_state_slip_law(V2,psi,fricsgl,xmu,dxmudv)
  use globalvar
  implicit none

  real (kind = dp) :: xmu, dxmudv
  real (kind = dp) :: V2,psi,psiss,fLV,fss,fssa
  real (kind = dp) :: A,B,L,f0,V0,fw,Vw
  real (kind = dp),dimension(100) :: fricsgl
  real (kind = dp) :: tmp, tmpc
  !
  A  = fricsgl(9)
  B  = fricsgl(10)
  L  = fricsgl(11)
  f0 = fricsgl(13)
  V0 = fricsgl(12)
  fw = fricsgl(14)
  Vw = fricsgl(15)

  tmpc = 1.0d0 / (2.0d0 * V0) * dexp(psi/A)
  tmp = (V2+1.d-30) * tmpc
  xmu = A * dlog(tmp + sqrt(tmp**2 + 1.0d0)) !arcsinh(z)= ln(z+sqrt(z^2+1))
  dxmudv = A * tmpc / sqrt(1.0d0 + tmp**2)  ! d(arcsinh(z))/dz = 1/sqrt(1+z^2)
  fLV = f0 - (B - A) * dlog(V2/V0)
  !fLV = max(1.0d-8, fLV)
  fss = fw + (fLV - fw) / ((1.0d0 + (V2/Vw)**8)**0.125d0)
  fssa = fss/A
  !fssa = max(1.0d-8, fssa)
  ! Using sinh(x) = (exp(x) - exp(-x))/2
  !psiss = A * dlog(2.0d0 * V0 / V2 * dsinh(fss/A))
  psiss = A * dlog(2.0d0 * V0 / V2 * (dexp(fssa) - dexp(-fssa))/2.0d0)
  psi = psiss + (psi - psiss) * dexp(-V2*dt/L)
  !
end subroutine rate_state_slip_law



