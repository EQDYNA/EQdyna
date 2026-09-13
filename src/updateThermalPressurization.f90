subroutine updateThermalPressurization

use globalvar
implicit none
!
! This subroutine is used to calculate pore pressure change due 
! to thermopressurization based on equation 12 and 13 in the benchmark
! description of TPV105 3D. 
! DL and BL, 20200901
!
integer(kind = 4)  :: i, ift, j, k
real   (kind = dp) :: htp, rouctp, lamta, gama, omega, kapa, tmp, tmp2, ker

do ift = 1, ntotft
    do i = 1, nftnd(ift)
    
        gama  = fric(FRIC_SLOT_TP_LAMBDA,i,ift)/fric(FRIC_SLOT_TP_ROUC,i,ift) ! tp_lambda/tp_rouc
        omega = fric(FRIC_SLOT_TP_A_HY,i,ift) ! tp_a_hy
        kapa  = fric(FRIC_SLOT_TP_A_TH,i,ift) ! tp_a_th

        tmp = 0.0d0
        do j = 1, nt-1
            ker = -kapa/(omega - kapa)/(4.0d0*kapa*(nt-j)*dt + 2.0d0*fric(FRIC_SLOT_TP_H,i,ift)**2)**0.5
            ker = ker + omega/(omega - kapa)/(4.0d0*omega*(nt-j)*dt + 2.0d0*fric(FRIC_SLOT_TP_H,i,ift)**2)**0.5
            tmp = tmp + abs(onFaultTPHist(2,i,j,ift))*onFaultTPHist(1,i,j,ift)*ker*dt
        enddo
        patnode(i,ift) = tmp*gama/(pi)**0.5

        tmp = 0.0d0
        do j = 1, nt-1
            ker = 1.0d0/(4.0d0*kapa*(nt-j)*dt + 2.0d0*fric(FRIC_SLOT_TP_H,i,ift)**2)**0.5
            tmp = tmp + abs(onFaultTPHist(2,i,j,ift))*onFaultTPHist(1,i,j,ift)*ker*dt
        enddo
        Tatnode(i,ift) = tmp/fric(FRIC_SLOT_TP_ROUC,i,ift)/(pi)**0.5
        fric(FRIC_SLOT_TP_NORM_TP,i,ift) = patnode(i,ift)
        fric(FRIC_SLOT_TP_TEMP,i,ift) = Tatnode(i,ift) + fric(FRIC_SLOT_TP_TINI,i,ift) ! + Tini
    enddo   
enddo

end subroutine updateThermalPressurization
  
  
    
  
