subroutine wedge(cenx, ceny, cenz, nelement, ntags, iy, iz, nftndtmp) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: nelement, ntags, iy, iz, nftndtmp, k, i, neworder(nen)
	real (kind = dp) :: cenx, ceny, cenz

	if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 	

		! Degenerate the brick element into two wedge elements.
		!     8
		!  5       7     
		! 	    6		 ! Brick
		!     4
		!  1       3
		! 	    2

		!     
		!  5       7(8)     
		! 	    6		 ! Wedge one
		!  1       3(4)
		! 	    2		
		et(nelement) = 11 ! et ==11 indicates a wedge element.
		ien(1,nelement) = plane1(iy-1,iz-1)
		ien(2,nelement) = plane2(iy-1,iz-1)
		ien(3,nelement) = plane2(iy,iz-1)
		ien(4,nelement) = ien(3,nelement)
		ien(5,nelement) = plane1(iy-1,iz)
		ien(6,nelement) = plane2(iy-1,iz)
		ien(7,nelement) = plane2(iy,iz)
		ien(8,nelement) = ien(7,nelement)	
		
		neworder = (/1,2,3,3,5,6,7,7/)
		call reorder(neworder, nelement, iy, iz)
		
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		!       6     
		!  7(8)      5      Wedge two
		!   	2       
		!  3(4)      1				
		nelement = nelement + 1
		et(nelement) = 11
		ien(1,nelement) = plane2(iy,iz-1)
		ien(2,nelement) = plane1(iy,iz-1)
		ien(3,nelement) = plane1(iy-1,iz-1)
		ien(4,nelement) = ien(3,nelement)
		ien(5,nelement) = plane2(iy,iz)
		ien(6,nelement) = plane1(iy,iz)	
		ien(7,nelement) = plane1(iy-1,iz)
		ien(8,nelement) = ien(7,nelement)
		
		neworder = (/3,4,1,1,7,8,5,5/)
		call reorder(neworder, nelement, iy, iz)
		
		stressCompIndexArr(nelement)=ntags
		ntags = ntags + 12

		mat(nelement,1)=material(1,1)
		mat(nelement,2)=material(1,2)
		mat(nelement,3)=material(1,3)				
		mat(nelement,5)=mat(nelement,2)**2*mat(nelement,3)!miu=vs**2*rho
		mat(nelement,4)=mat(nelement,1)**2*mat(nelement,3)-2*mat(nelement,5)!lam=vp**2*rho-2*miu	

		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement) == nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 		
	endif 
end subroutine

subroutine wedge4num(cenx, ceny, cenz, nelement) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: nelement	
	real (kind = dp) :: cenx, ceny, cenz
	
	if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 
		nelement = nelement + 1
	endif 
end subroutine

subroutine tetra(cenx, ceny, cenz, nelement, ntags, iy, iz, nftndtmp) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: nelement, ntags, iy, iz, nftndtmp, k, i, neworder(nen)
	real (kind = dp) :: cenx, ceny, cenz

		if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 	
		! Orders in hexahedrons. 
		
		! ien(1,nelement) = plane1(iy-1,iz-1)
		! ien(2,nelement) = plane2(iy-1,iz-1)
		! ien(3,nelement) = plane2(iy,iz-1)
		! ien(4,nelement) = plane1(iy,iz-1)
		! ien(5,nelement) = plane1(iy-1,iz)
		! ien(6,nelement) = plane2(iy-1,iz)
		! ien(7,nelement) = plane2(iy,iz)
		! ien(8,nelement) = plane1(iy,iz)	
		
		!2-1-3-6; 
		!5-6-8-1; 
		!7-8-6-3;
		!4-1-8-3;
		!6-3-8-1;
		 
		 et(nelement)=12
		 ien(1,nelement) = plane1(iy-1,iz-1)!#1
		 ien(2,nelement) = plane2(iy,iz-1)!#3
		 ien(3,nelement) = plane2(iy-1,iz)!#6
		 ien(5,nelement) = plane2(iy-1,iz-1)!#2
		 ien(4,nelement) = ien(3,nelement)
		 ien(6,nelement) = ien(5,nelement)
		 ien(7,nelement) = ien(5,nelement)
		 ien(8,nelement) = ien(7,nelement)
		 
		neworder = (/1,3,6,6,2,2,2,2/)
		call reorder(neworder, nelement, iy, iz)	
		
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
		 !STEP2: Create the second one.
		 nelement=nelement+1
		 ien(1,nelement) = plane2(iy-1,iz)!#6   5-6-8-1
		 ien(2,nelement) = plane1(iy,iz)!#8				  
		 ien(3,nelement) = plane1(iy-1,iz-1)!#1
		 ien(5,nelement) = plane1(iy-1,iz)!#5
		 ien(4,nelement) = ien(3,nelement)
		 ien(6,nelement) = ien(5,nelement)
		 ien(7,nelement) = ien(5,nelement)
		 ien(8,nelement) = ien(7,nelement)

		 
		neworder = (/6,8,1,1,5,5,5,5/)
		call reorder(neworder, nelement, iy, iz)	
				 
		 mat(nelement,1)=mat(nelement-1,1)
		 mat(nelement,2)=mat(nelement-1,2)
		 mat(nelement,3)=mat(nelement-1,3)
		 mat(nelement,4)=mat(nelement-1,4)
		 mat(nelement,5)=mat(nelement-1,5)
		 et(nelement)=12 
		 stressCompIndexArr(nelement)= ntags
		 ntags = ntags + 12					
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
		!STEP3: Create the third one.
		nelement=nelement+1	 
		ien(1,nelement) = plane1(iy,iz)!#8   7-8-6-3
		ien(2,nelement) = plane2(iy-1,iz)!6				  
		ien(3,nelement) = plane2(iy,iz-1)!#3
		ien(5,nelement) = plane2(iy,iz)!#7			
		 ien(4,nelement) = ien(3,nelement)
		 ien(6,nelement) = ien(5,nelement)
		 ien(7,nelement) = ien(5,nelement)
		 ien(8,nelement) = ien(7,nelement)
		 
		neworder = (/8,6,3,3,7,7,7,7/)
		call reorder(neworder, nelement, iy, iz)	
				 
		mat(nelement,1)=mat(nelement-1,1)
		mat(nelement,2)=mat(nelement-1,2)
		mat(nelement,3)=mat(nelement-1,3)
		mat(nelement,4)=mat(nelement-1,4)
		mat(nelement,5)=mat(nelement-1,5)	
		et(nelement)=12 
		stressCompIndexArr(nelement)= ntags
		ntags = ntags + 12					
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
			
		!STEP4: Create the 4th one.
		nelement = nelement + 1
		ien(1,nelement) = plane1(iy-1,iz-1)!#1 4-1-8-3;
		ien(2,nelement) = plane1(iy,iz)!#8					  
		ien(3,nelement) = plane2(iy,iz-1)!#3
		ien(5,nelement) = plane1(iy,iz-1)!#4			
		 ien(4,nelement) = ien(3,nelement)
		 ien(6,nelement) = ien(5,nelement)
		 ien(7,nelement) = ien(5,nelement)
		 ien(8,nelement) = ien(7,nelement)
		 
		neworder = (/1,8,3,3,4,4,4,4/)
		call reorder(neworder, nelement, iy, iz)	
				 
		mat(nelement,1)=mat(nelement-1,1)
		mat(nelement,2)=mat(nelement-1,2)
		mat(nelement,3)=mat(nelement-1,3)
		mat(nelement,4)=mat(nelement-1,4)
		mat(nelement,5)=mat(nelement-1,5)
		et(nelement)=12 
		stressCompIndexArr(nelement)= ntags
		ntags = ntags + 12
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 	

		!STEP5: Create the 5th one.					
		nelement = nelement + 1
		ien(1,nelement) = plane2(iy,iz-1)!#3  6-3-8-1;  
		ien(2,nelement) = plane1(iy,iz)!#8
		ien(3,nelement) = plane1(iy-1,iz-1)!#1	
		ien(5,nelement) = plane2(iy-1,iz)!#6		
		 ien(4,nelement) = ien(3,nelement)
		 ien(6,nelement) = ien(5,nelement)
		 ien(7,nelement) = ien(5,nelement)
		 ien(8,nelement) = ien(7,nelement)
		 
		neworder = (/3,8,1,1,6,6,6,6/)
		call reorder(neworder, nelement, iy, iz)	
				 
		mat(nelement,1)=mat(nelement-1,1)
		mat(nelement,2)=mat(nelement-1,2)
		mat(nelement,3)=mat(nelement-1,3)
		mat(nelement,4)=mat(nelement-1,4)
		mat(nelement,5)=mat(nelement-1,5)
		et(nelement)=12 
		stressCompIndexArr(nelement)=ntags
		ntags = ntags + 12
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,nelement)==nsmp(1,i,1)) then
						ien(k,nelement) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
	endif		
end subroutine

subroutine tetra4num(cenx, ceny, cenz, nelement) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: nelement	
	real (kind = dp) :: cenx, ceny, cenz
	
	if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 	
		nelement=nelement+4
	endif 
end subroutine

subroutine reorder(neworder, nelement, iy, iz)

	use globalvar
	implicit none
	
	integer (kind = 4) :: i, neworder(nen), ientmp(nen), nelement, iz, iy
	
		ientmp(1) = plane1(iy-1,iz-1)
		ientmp(2) = plane2(iy-1,iz-1)
		ientmp(3) = plane2(iy,iz-1)
		ientmp(4) = plane1(iy,iz-1)
		ientmp(5) = plane1(iy-1,iz)
		ientmp(6) = plane2(iy-1,iz)
		ientmp(7) = plane2(iy,iz)
		ientmp(8) = plane1(iy,iz)	
		
		do i = 1, nen
			ien(i,nelement) = ientmp(neworder(i))
		enddo 
		
end subroutine reorder