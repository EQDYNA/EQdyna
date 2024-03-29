subroutine thermop

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
	
		gama  = fric(19,i,ift)/fric(18,i,ift) ! tp_lambda/tp_rouc
		omega = fric(16,i,ift) ! tp_a_hy
		kapa  = fric(17,i,ift) ! tp_a_th 

		tmp = 0.0d0 
		do j = 1, nt-1
			ker = -kapa/(omega - kapa)/(4.0d0*kapa*(nt-j)*dt + 2.0d0*fric_tp_h**2)**0.5
			ker = ker + omega/(omega - kapa)/(4.0d0*omega*(nt-j)*dt + 2.0d0*fric_tp_h**2)**0.5
			tmp = tmp + abs(onFaultTPHist(2,i,j,ift))*onFaultTPHist(1,i,j,ift)*ker*dt  
		enddo 
		patnode(i,ift) = tmp*gama/(pi)**0.5

		tmp = 0.0d0 
		do j = 1, nt-1
			ker = 1.0d0/(4.0d0*kapa*(nt-j)*dt + 2.0d0*fric_tp_h**2)**0.5
			tmp = tmp + abs(onFaultTPHist(2,i,j,ift))*onFaultTPHist(1,i,j,ift)*ker*dt
		enddo 
		Tatnode(i,ift) = tmp/fric(18,i,ift)/(pi)**0.5
		fric(51,i,ift) = patnode(i,ift) 
		fric(52,i,ift) = Tatnode(i,ift) + fric(41,i,ift) ! + Tini
	enddo 	
enddo

end subroutine thermop
  
  
    
  
