subroutine wedge(cenx, ceny, cenz, elemCount, ntags, iy, iz, nftndtmp) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: elemCount, ntags, iy, iz, nftndtmp, k, i, neworder(nen)
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
		et(elemCount) = 11 ! et ==11 indicates a wedge element.
		ien(1,elemCount) = plane1(iy-1,iz-1)
		ien(2,elemCount) = plane2(iy-1,iz-1)
		ien(3,elemCount) = plane2(iy,iz-1)
		ien(4,elemCount) = ien(3,elemCount)
		ien(5,elemCount) = plane1(iy-1,iz)
		ien(6,elemCount) = plane2(iy-1,iz)
		ien(7,elemCount) = plane2(iy,iz)
		ien(8,elemCount) = ien(7,elemCount)	
		
		neworder = (/1,2,3,3,5,6,7,7/)
		call reorder(neworder, elemCount, iy, iz)
		
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		!       6     
		!  7(8)      5      Wedge two
		!   	2       
		!  3(4)      1				
		elemCount = elemCount + 1
		et(elemCount) = 11
		ien(1,elemCount) = plane2(iy,iz-1)
		ien(2,elemCount) = plane1(iy,iz-1)
		ien(3,elemCount) = plane1(iy-1,iz-1)
		ien(4,elemCount) = ien(3,elemCount)
		ien(5,elemCount) = plane2(iy,iz)
		ien(6,elemCount) = plane1(iy,iz)	
		ien(7,elemCount) = plane1(iy-1,iz)
		ien(8,elemCount) = ien(7,elemCount)
		
		neworder = (/3,4,1,1,7,8,5,5/)
		call reorder(neworder, elemCount, iy, iz)
		
		stressCompIndexArr(elemCount)=ntags
		ntags = ntags + 12

		mat(elemCount,1)=material(1,1)
		mat(elemCount,2)=material(1,2)
		mat(elemCount,3)=material(1,3)				
		mat(elemCount,5)=mat(elemCount,2)**2*mat(elemCount,3)!miu=vs**2*rho
		mat(elemCount,4)=mat(elemCount,1)**2*mat(elemCount,3)-2*mat(elemCount,5)!lam=vp**2*rho-2*miu	

		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount) == nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 		
	endif 
end subroutine

subroutine wedge4num(cenx, ceny, cenz, elemCount) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: elemCount	
	real (kind = dp) :: cenx, ceny, cenz
	
	if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 
		elemCount = elemCount + 1
	endif 
end subroutine

subroutine tetra(cenx, ceny, cenz, elemCount, ntags, iy, iz, nftndtmp) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: elemCount, ntags, iy, iz, nftndtmp, k, i, neworder(nen)
	real (kind = dp) :: cenx, ceny, cenz

		if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 	
		! Orders in hexahedrons. 
		
		! ien(1,elemCount) = plane1(iy-1,iz-1)
		! ien(2,elemCount) = plane2(iy-1,iz-1)
		! ien(3,elemCount) = plane2(iy,iz-1)
		! ien(4,elemCount) = plane1(iy,iz-1)
		! ien(5,elemCount) = plane1(iy-1,iz)
		! ien(6,elemCount) = plane2(iy-1,iz)
		! ien(7,elemCount) = plane2(iy,iz)
		! ien(8,elemCount) = plane1(iy,iz)	
		
		!2-1-3-6; 
		!5-6-8-1; 
		!7-8-6-3;
		!4-1-8-3;
		!6-3-8-1;
		 
		 et(elemCount)=12
		 ien(1,elemCount) = plane1(iy-1,iz-1)!#1
		 ien(2,elemCount) = plane2(iy,iz-1)!#3
		 ien(3,elemCount) = plane2(iy-1,iz)!#6
		 ien(5,elemCount) = plane2(iy-1,iz-1)!#2
		 ien(4,elemCount) = ien(3,elemCount)
		 ien(6,elemCount) = ien(5,elemCount)
		 ien(7,elemCount) = ien(5,elemCount)
		 ien(8,elemCount) = ien(7,elemCount)
		 
		neworder = (/1,3,6,6,2,2,2,2/)
		call reorder(neworder, elemCount, iy, iz)	
		
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
		 !STEP2: Create the second one.
		 elemCount=elemCount+1
		 ien(1,elemCount) = plane2(iy-1,iz)!#6   5-6-8-1
		 ien(2,elemCount) = plane1(iy,iz)!#8				  
		 ien(3,elemCount) = plane1(iy-1,iz-1)!#1
		 ien(5,elemCount) = plane1(iy-1,iz)!#5
		 ien(4,elemCount) = ien(3,elemCount)
		 ien(6,elemCount) = ien(5,elemCount)
		 ien(7,elemCount) = ien(5,elemCount)
		 ien(8,elemCount) = ien(7,elemCount)

		 
		neworder = (/6,8,1,1,5,5,5,5/)
		call reorder(neworder, elemCount, iy, iz)	
				 
		 mat(elemCount,1)=mat(elemCount-1,1)
		 mat(elemCount,2)=mat(elemCount-1,2)
		 mat(elemCount,3)=mat(elemCount-1,3)
		 mat(elemCount,4)=mat(elemCount-1,4)
		 mat(elemCount,5)=mat(elemCount-1,5)
		 et(elemCount)=12 
		 stressCompIndexArr(elemCount)= ntags
		 ntags = ntags + 12					
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
		!STEP3: Create the third one.
		elemCount=elemCount+1	 
		ien(1,elemCount) = plane1(iy,iz)!#8   7-8-6-3
		ien(2,elemCount) = plane2(iy-1,iz)!6				  
		ien(3,elemCount) = plane2(iy,iz-1)!#3
		ien(5,elemCount) = plane2(iy,iz)!#7			
		 ien(4,elemCount) = ien(3,elemCount)
		 ien(6,elemCount) = ien(5,elemCount)
		 ien(7,elemCount) = ien(5,elemCount)
		 ien(8,elemCount) = ien(7,elemCount)
		 
		neworder = (/8,6,3,3,7,7,7,7/)
		call reorder(neworder, elemCount, iy, iz)	
				 
		mat(elemCount,1)=mat(elemCount-1,1)
		mat(elemCount,2)=mat(elemCount-1,2)
		mat(elemCount,3)=mat(elemCount-1,3)
		mat(elemCount,4)=mat(elemCount-1,4)
		mat(elemCount,5)=mat(elemCount-1,5)	
		et(elemCount)=12 
		stressCompIndexArr(elemCount)= ntags
		ntags = ntags + 12					
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
			
		!STEP4: Create the 4th one.
		elemCount = elemCount + 1
		ien(1,elemCount) = plane1(iy-1,iz-1)!#1 4-1-8-3;
		ien(2,elemCount) = plane1(iy,iz)!#8					  
		ien(3,elemCount) = plane2(iy,iz-1)!#3
		ien(5,elemCount) = plane1(iy,iz-1)!#4			
		 ien(4,elemCount) = ien(3,elemCount)
		 ien(6,elemCount) = ien(5,elemCount)
		 ien(7,elemCount) = ien(5,elemCount)
		 ien(8,elemCount) = ien(7,elemCount)
		 
		neworder = (/1,8,3,3,4,4,4,4/)
		call reorder(neworder, elemCount, iy, iz)	
				 
		mat(elemCount,1)=mat(elemCount-1,1)
		mat(elemCount,2)=mat(elemCount-1,2)
		mat(elemCount,3)=mat(elemCount-1,3)
		mat(elemCount,4)=mat(elemCount-1,4)
		mat(elemCount,5)=mat(elemCount-1,5)
		et(elemCount)=12 
		stressCompIndexArr(elemCount)= ntags
		ntags = ntags + 12
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 	

		!STEP5: Create the 5th one.					
		elemCount = elemCount + 1
		ien(1,elemCount) = plane2(iy,iz-1)!#3  6-3-8-1;  
		ien(2,elemCount) = plane1(iy,iz)!#8
		ien(3,elemCount) = plane1(iy-1,iz-1)!#1	
		ien(5,elemCount) = plane2(iy-1,iz)!#6		
		 ien(4,elemCount) = ien(3,elemCount)
		 ien(6,elemCount) = ien(5,elemCount)
		 ien(7,elemCount) = ien(5,elemCount)
		 ien(8,elemCount) = ien(7,elemCount)
		 
		neworder = (/3,8,1,1,6,6,6,6/)
		call reorder(neworder, elemCount, iy, iz)	
				 
		mat(elemCount,1)=mat(elemCount-1,1)
		mat(elemCount,2)=mat(elemCount-1,2)
		mat(elemCount,3)=mat(elemCount-1,3)
		mat(elemCount,4)=mat(elemCount-1,4)
		mat(elemCount,5)=mat(elemCount-1,5)
		et(elemCount)=12 
		stressCompIndexArr(elemCount)=ntags
		ntags = ntags + 12
		if (ceny>0.0d0) then 
			do i=1,nftndtmp
				do k=1,nen
					if(ien(k,elemCount)==nsmp(1,i,1)) then
						ien(k,elemCount) = nsmp(2,i,1)  !use master node for the node!
					endif
				enddo
			enddo		
		endif 
		
	endif		
end subroutine

subroutine tetra4num(cenx, ceny, cenz, elemCount) 
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: elemCount	
	real (kind = dp) :: cenx, ceny, cenz
	
	if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and.ceny<dx.and.ceny>-dx.and.cenz>fltxyz(1,3,1)) then 	
		elemCount=elemCount+4
	endif 
end subroutine

subroutine reorder(neworder, elemCount, iy, iz)

	use globalvar
	implicit none
	
	integer (kind = 4) :: i, neworder(nen), ientmp(nen), elemCount, iz, iy
	
		ientmp(1) = plane1(iy-1,iz-1)
		ientmp(2) = plane2(iy-1,iz-1)
		ientmp(3) = plane2(iy,iz-1)
		ientmp(4) = plane1(iy,iz-1)
		ientmp(5) = plane1(iy-1,iz)
		ientmp(6) = plane2(iy-1,iz)
		ientmp(7) = plane2(iy,iz)
		ientmp(8) = plane1(iy,iz)	
		
		do i = 1, nen
			ien(i,elemCount) = ientmp(neworder(i))
		enddo 
		
end subroutine reorder