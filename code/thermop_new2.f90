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
  real(kind=8) :: Tatnode(nftnd, numtp), patnode(nftnd,numtp), Tpartt(nftnd,numtp), ppartt(nftnd,numtp), &
	fric(100,nftnd), Tatnode1(nftnd, numtp), patnode1(nftnd,numtp)
  real(kind=8) :: htp, rouctp, lamta, aa0, aa1, aa2, aa3 
  
  
  do i = 1, nftnd
	aa0 = fric_tp_a_th/dxtp/dxtp/2.0d0
	aa1 = 1.0d0/fric_tp_rouc/fric_tp_h/(2.0d0*pi)**0.5*abs(fric(50,i))*fric(49,i)
	aa2 = 1.0d0/dt + fric_tp_a_th/dxtp/dxtp 
	aa3 = 1.0d0/dt
	
	Tatnode1(i,numtp) = fric_tp_Tini
	Tatnode1(i,numtp-1) = aa2*Tatnode1(i,numtp) - aa3*Tatnode(i,numtp) - aa0*(fric_tp_Tini + Tatnode(i,numtp-1) - 2.0d0*fric_tp_Tini) - aa1*dexp(-((numtp-1)*dxtp)**2/fric_tp_h**2/2.0d0) - aa0*fric_tp_Tini 
	Tatnode1(i,numtp-1) = Tatnode1(i,numtp-1)/aa0 	
	do j = numtp-1, 2, -1 
		Tatnode1(i,j-1) = aa2*Tatnode1(i,j) - aa3*Tatnode(i,j) - aa0*(Tatnode(i,j+1) + Tatnode(i,j-1) - 2.0d0*Tatnode(i,j)) - aa1*dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0) - aa0*Tatnode1(i,j+1) 		
		Tatnode1(i,j-1) = Tatnode1(i,j-1)/aa0
	enddo 
	! The temperature for t+1 at nodes from j = 1, numtp are solved. 
	do j = 1, numtp
		if (j == 1) then
			Tpartt(i,j) = aa0*(Tatnode1(i,j+1) + Tatnode1(i,j+1) - 2.0d0*Tatnode1(i,j)) + aa0*(Tatnode(i,j+1) + Tatnode(i,j+1) - 2.0d0*Tatnode(i,j)) + aa1*dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)
		elseif (j>1 .and. j<numtp) then 
			Tpartt(i,j) = aa0*(Tatnode1(i,j+1) + Tatnode1(i,j-1) - 2.0d0*Tatnode1(i,j)) + aa0*(Tatnode(i,j+1) + Tatnode(i,j-1) - 2.0d0*Tatnode(i,j)) + aa1*dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)		
		elseif (j == numtp) then
			Tpartt(i,j) = 0.0d0 
		endif 
	enddo 
    if (i==1) then 
		write(*,*) 'T profile at t', Tatnode(i,1), Tatnode(i,2), Tatnode(i,3),Tatnode(i,4), Tatnode(i,5), Tatnode(i,6), Tatnode(i,numtp-5), Tatnode(i,numtp-4), Tatnode(i,numtp-3),Tatnode(i,numtp-2), Tatnode(i,numtp-1),Tatnode(i,numtp)
		write(*,*) 'T profile at t+1', Tatnode1(i,1), Tatnode1(i,2), Tatnode1(i,3), Tatnode1(i,4), Tatnode1(i,5), Tatnode1(i,6), Tatnode1(i,numtp-5), Tatnode1(i,numtp-4),Tatnode1(i,numtp-3), Tatnode1(i,numtp-2), Tatnode1(i,numtp-1),Tatnode1(i,numtp)
	endif 
	
	aa0 = fric(20,i)/dxtp/dxtp/2.0d0
	aa1 = fric_tp_lambda
	aa2 = 1.0d0/dt + fric(20,i)/dxtp/dxtp 
	aa3 = 1.0d0/dt  
 	do j = numtp, 2, -1
		if (j == numtp) then 
			patnode1(i,numtp) = fric_tp_pini
			patnode1(i,j-1) = aa2*patnode1(i,j) - aa3*patnode(i,j) - aa0*(fric_tp_pini + patnode(i,j-1) - 2.0d0*fric_tp_pini) - aa1*Tpartt(i,j) - aa0*fric_tp_pini 
			patnode1(i,j-1) = patnode1(i,j-1)/aa0
		elseif (j>2 .and. j<numtp) then 
			patnode1(i,j-1) = aa2*patnode1(i,j) - aa3*patnode(i,j) - aa0*(patnode(i,j+1) + patnode(i,j-1) - 2.0d0*patnode(i,j)) - aa1*Tpartt(i,j) - aa0*patnode1(i,j+1) 	
			patnode1(i,j-1) = patnode1(i,j-1)/aa0			
		elseif (j == 2) then 
			patnode1(i,j-1) = aa2*patnode1(i,j) - aa3*patnode(i,j) - aa0*(patnode(i,j+1) + patnode(i,j-1) - 2.0d0*patnode(i,j)) - aa1*Tpartt(i,j) - aa0*patnode1(i,j+1) 
			patnode1(i,j-1) = patnode1(i,j-1)/aa0			
		endif 
	enddo 
	! The pore pressure for t+1 at nodes from j = 1, numtp are solved. 
	! do j = 1, numtp
		! if (j == 1) then
			! ppartt(j) = aa0*(patnode1(j+1) + patnode1(j+1) - 2.0d0*patnode1(j)) + aa0*(patnode(j+1) + patnode(j+1) - 2.0d0*patnode(j)) + aa1*dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)
		! elseif (j>1 .and. j<numtp) then 
			! ppartt(j) = aa0*(patnode1(j+1) + patnode1(j-1) - 2.0d0*patnode1(j)) + aa0*(patnode(j+1) + patnode(j-1) - 2.0d0*patnode(j)) + aa1*dexp(-((j-1)*dxtp)**2/fric_tp_h**2/2.0d0)		
		! elseif (j == numtp) then
			! ppartt(j) = 0.0d0 
		! endif 
	! enddo  
    if (i==1) then 
		write(*,*) 'p profile at t', patnode(i,1), patnode(i,2), patnode(i,3), patnode(i,numtp-2), patnode(i,numtp-1),patnode(i,numtp)
		write(*,*) 'p profile at t+1', patnode1(i,1), patnode1(i,2), patnode1(i,3), patnode1(i,numtp-2), patnode1(i,numtp-1),patnode1(i,numtp)
	endif 
	do j = 1, numtp
		Tatnode(i,j) = Tatnode1(i,j) 
		patnode(i,j) = patnode1(i,j)
	enddo
	
	fric(51,i) = patnode(i,1) - fric_tp_pini
	fric(52,i) = Tatnode(i,1)
	!write(*,*) 'ftnode,i,p,T',fric(51,i)/1.0d6,fric(52,i)
	
  enddo


end subroutine thermop
  
  
    
  
