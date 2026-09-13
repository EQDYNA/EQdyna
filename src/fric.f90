! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine slip_weak(slip,fricsgl,xmu)
  use globalvar
  implicit none

  real (kind = dp) :: xmu,slip
  real (kind = dp),dimension(20) :: fricsgl
  !
  if(abs(slip).lt.1.0e-10) then
    xmu = fricsgl(FRIC_SLOT_SW_FS)    !xmu is frictional coefficient, node by node on fault
  elseif(slip < fricsgl(FRIC_SLOT_SW_D0)) then
    xmu = fricsgl(FRIC_SLOT_SW_FS) - (fricsgl(FRIC_SLOT_SW_FS) - fricsgl(FRIC_SLOT_SW_FD))*slip/fricsgl(FRIC_SLOT_SW_D0)
  endif
  if(slip >= fricsgl(FRIC_SLOT_SW_D0)) then
    xmu = fricsgl(FRIC_SLOT_SW_FD)
  endif
  !
end subroutine slip_weak

subroutine time_weak(trupt,fricsgl,xmu)

    use globalvar
    implicit none

    real (kind = dp) :: xmu,trupt
    real (kind = dp),dimension(20) :: fricsgl
    
    if(trupt <= 0.0d0) then
        xmu = fricsgl(FRIC_SLOT_SW_FS)
    elseif(trupt < fricsgl(FRIC_SLOT_TW_T0)) then
        xmu = fricsgl(FRIC_SLOT_SW_FS) - (fricsgl(FRIC_SLOT_SW_FS) - fricsgl(FRIC_SLOT_SW_FD))*trupt/fricsgl(FRIC_SLOT_TW_T0)
    else
        xmu = fricsgl(FRIC_SLOT_SW_FD)
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
  A  = fricsgl(FRIC_SLOT_RSF_A)
  B  = fricsgl(FRIC_SLOT_RSF_B)
  L  = fricsgl(FRIC_SLOT_RSF_DC)
  f0 = fricsgl(FRIC_SLOT_RSF_R0)
  V0 = fricsgl(FRIC_SLOT_RSF_V0)

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
  A  = fricsgl(FRIC_SLOT_RSF_A)
  B  = fricsgl(FRIC_SLOT_RSF_B)
  L  = fricsgl(FRIC_SLOT_RSF_DC)
  f0 = fricsgl(FRIC_SLOT_RSF_R0)
  V0 = fricsgl(FRIC_SLOT_RSF_V0)
  fw = fricsgl(FRIC_SLOT_RSF_FW)
  Vw = fricsgl(FRIC_SLOT_RSF_VW)

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



