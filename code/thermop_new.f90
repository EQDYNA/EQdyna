subroutine thermop(nsmp, Tatnode, patnode, ift, nftnd, fric, nt)
  use globalvar
  implicit none
  !
  ! This subroutine is used to calculate pore pressure change due 
  ! to thermopressurization based on equation 12 and 13 in the benchmark
  ! description of TPV105 3D. 
  ! DL and BL, 20200901
  !
  integer(kind=4) :: i, ift, j, k, nt, nsmp(2,nftnd), nftnd
  real(kind=8) :: Tatnode(nftnd, numtp), Tatnode1(nftnd,numtp), patnode(nftnd,numtp), patnode1(nftnd,numtp), Tpartt(nftnd,numtp), ppartt(nftnd,numtp), &
	fric(100,nftnd)
  real(kind=8) :: htp, rouctp, lamta
    
  do i = 1, nftnd
    do j = 2,numtp-1
		Tpartt(i,j) = fric_tp_a_th* (Tatnode(i,j+1) - 2.0d0*Tatnode(i,j) + Tatnode(i,j-1))/dxtp/dxtp + &
					dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)/fric_tp_rouc/fric_tp_h/(2.0d0*pi)**0.5*abs(fric(50,i))*fric(49,i)
					
		Tatnode1(i,j) = Tatnode(i,j) + dt * Tpartt(i,j)
	enddo 
	Tatnode1(i,1) = Tatnode(i,2)
	Tatnode1(i,numtp) = fric_tp_Tini
	do j = 1, numtp
		Tatnode1(i,j) = (Tatnode1(i,j) + Tatnode(i,j))/2.0d0
	enddo 
    do j = 2,numtp-1
		Tpartt(i,j) = fric_tp_a_th* (Tatnode1(i,j+1) - 2.0d0*Tatnode1(i,j) + Tatnode1(i,j-1))/dxtp/dxtp + &
					dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)/fric_tp_rouc/fric_tp_h/(2.0d0*pi)**0.5*abs(fric(50,i))*fric(49,i)					
		Tatnode(i,j) = Tatnode(i,j) + dt * Tpartt(i,j)
	enddo 		
	Tpartt(i,1) = 1.0d0/fric_tp_rouc/fric_tp_h/(2.0d0*pi)**0.5*abs(fric(50,i))*fric(49,i)
	Tatnode(i,1) = Tatnode(i,1) + dt * Tpartt(i,1)
	Tatnode(i,numtp) = fric_tp_Tini
	
	do j = 2, numtp-1
		ppartt(i,j) = fric(20,i)* (patnode(i,j+1) - 2.0d0*patnode(i,j) + patnode(i,j-1))/dxtp/dxtp + &
					fric_tp_lambda * Tpartt(i,j)
					!fric(20,i): alpha_hy
		patnode1(i,j) = patnode(i,j) + dt * ppartt(i,j)
	enddo 
	patnode1(i,1) = patnode(i,2)
	patnode1(i,numtp) = fric_tp_pini
	do j = 1, numtp
		patnode1(i,j) = (patnode1(i,j) + patnode(i,j))/2.0d0
	enddo 	
	do j = 2, numtp-1
		ppartt(i,j) = fric(20,i)* (patnode1(i,j+1) - 2.0d0*patnode1(i,j) + patnode1(i,j-1))/dxtp/dxtp + &
					fric_tp_lambda * Tpartt(i,j)
					!fric(20,i): alpha_hy
		patnode(i,j) = patnode(i,j) + dt * ppartt(i,j)
	enddo 	
	
	ppartt(i,1) = fric_tp_lambda * Tpartt(i,1)
	patnode(i,1) = patnode(i,1) + dt * ppartt(i,1)
	patnode(i,numtp) = fric_tp_pini
	
	fric(51,i) = patnode(i,1) - fric_tp_pini
	fric(52,i) = Tatnode(i,1)
	!write(*,*) 'ftnode,i,p,T',fric(51,i)/1.0d6,fric(52,i)
  enddo


end subroutine thermop
  
  
    
  
