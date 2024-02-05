! Subroutine list:
! 1. output_onfault_st
! 2. output_offfault_st
! 3. output_frt
! 4. output_timeanalysis
! 5. output_plastic_strain
!#1
subroutine output_onfault_st

	use globalvar
	implicit none
	
	integer (kind = 4) :: i, j 
	
	if(n4onf>0) then
		do i=1,n4onf
			j=anonfs(3,i)
			if(j==1)  then  !main fault stations
				sttmp = '      '
				dptmp = '      '
				write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)/100.d0) 
				write(dptmp,'(i4.3)') int(abs(xonfs(2,anonfs(2,i),j))/100.d0) 
				open(51,file='faultst'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

				sttmp = '      '
				dptmp = '      '
				write(sttmp,'(f5.1)') xonfs(1,anonfs(2,i),j)/1000.d0 
				write(dptmp,'(f5.1)') abs(xonfs(2,anonfs(2,i),j)/1000.d0) 
				loca = '# location = on fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km down-dip'		
			endif
			write(51,*) '# ',projectname
			write(51,*) '# Author=',author
			call date_and_time(values=time_array)
			write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
				'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
				':',time_array(7)
			write(51,*) '# code = EQdyna3D'
			write(51,*) '# element_size =',dx
			write(51,'( a14,f8.4,a3)') '# time_step =', dt, ' s'
			write(51,'( a19,i6)') '# num_time_steps =', locplt-1
			!write(51,*) loca !Disable the loca to avoid writting in two lines.
			write(51,*) '# Time series in 11 columns in format E15.7'
			write(51,*) '# Column #1 = Time (s)'
			write(51,*) '# Column #2 = horizontal slip (m)'
			write(51,*) '# Column #3 = horizontal slip rate (m/s)'
			write(51,*) '# Column #4 = horizontal shear stress (MPa)'
			write(51,*) '# Column #5 = down-dip slip (m)'
			write(51,*) '# Column #6 = down-dip slip rate (m/s)'
			write(51,*) '# Column #7 = down-dip shear stress (MPa)'
			write(51,*) '# Column #8 = normal stress (MPa)'
			write(51,*) '# Column #9 = state variable psi (dimensionless)'
			write(51,*) '# Column #10 = Temperature (degrees Kelvin)'
			write(51,*) '# Column #11 = Pore pressure (MPa)'		
			write(51,*) '# The line below lists the names of the data fields:'
			write(51,'(1X,103A)') 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress psi temperature pressure'
			write(51,*) '#'
					! fltsta(1,locplt-1,j) = time
					! fltsta(2,locplt-1,j) = sliprates
					! fltsta(3,locplt-1,j) = sliprated
					! fltsta(4,locplt-1,j) = state(i)
					! fltsta(5,locplt-1,j) = slips
					! fltsta(6,locplt-1,j) = slipd
					! fltsta(7,locplt-1,j) = slipn
					! fltsta(8,locplt-1,j) = tstk
					! fltsta(9,locplt-1,j) = tdip
					! fltsta(10,locplt-1,j) = tnrm
					! fltsta(11,locplt-1,j) = fric(51,i)
					! fltsta(12,locplt-1,j) = fric(52,i)		
			do j=1,locplt-1
				write(51,'( E21.13,10E15.7)') fltsta(1,j,i),fltsta(5,j,i),fltsta(2,j,i),fltsta(8,j,i)/1.0d6,&
						-fltsta(6,j,i),-fltsta(3,j,i),-fltsta(9,j,i)/1.0d6, -fltsta(10,j,i)/1.0d6, fltsta(4,j,i), fltsta(12,j,i), fltsta(11,j,i)/1.0d6  
			enddo
			close(51)
		enddo
	endif 
end subroutine output_onfault_st
!#2
subroutine output_offfault_st
	
	use globalvar
	implicit none
	
	integer (kind = 4) :: i, j 	
	
	if(n4out>0) then
		do i=1,n4out
			bodytmp = '      '
			sttmp = '      '
			dptmp = '      '
			write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i))/100.d0) 
			write(sttmp,'(i4.3)') int(x4nds(1,an4nds(1,i))/100.d0) 
			write(dptmp,'(i4.3)') int(abs(x4nds(3,an4nds(1,i)))/100.d0) 

			open(51,file='body'//trim(adjustl(bodytmp))//'st'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

			bodytmp = '      '
			sttmp = '      '
			dptmp = '      '
			write(bodytmp,'(f5.1)') x4nds(2,an4nds(1,i))/1000. 
			write(sttmp,'(f5.1)') x4nds(1,an4nds(1,i))/1000. 
			write(dptmp,'(f5.1)') abs(x4nds(3,an4nds(1,i)))/1000. 
			loca = '# location = '//trim(adjustl(bodytmp))//' km off fault, '//trim(adjustl(sttmp))//' km along strike'//trim(adjustl(dptmp))//' km depth'
			write(51,*) '# ',projectname
			write(51,*) '# Author=',author
			call date_and_time(values=time_array)
			write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
					'/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
					':',time_array(7)
			write(51,*) '# code = EQdyna3D'
			write(51,*) '# element_size =',dx
			write(51,'( a14,f8.4,a3)') '# time_step=', dt, ' s'
			write(51,'( a19,i6)') '# num_time_steps=',locplt
			!write(51,*) loca
			!write(51,*) '# Time series in 11 columns in format E15.7 for data'
			write(51,*) '# Column #1 = Time (s)'
			write(51,*) '# Column #2 = horizontal displacement (m)'
			write(51,*) '# Column #3 = horizontal displacement (m)'
			write(51,*) '# Column #3 = horizontal velocity (m/s)'
			write(51,*) '# Column #4 = vertical displacement (m)'
			write(51,*) '# Column #5 = vertical velocity (m/s)'
			write(51,*) '# Column #6 = normal displacement (m)'
			write(51,*) '# Column #7 = normal velocity (m/s)'
			write(51,*) '#'
			write(51,*) '# The line below lists the names of the data fields:'
			write(51,*) 't h-disp h-vel v-disp v-vel n-disp n-vel'
			do j=1,locplt
				write(51,'( E21.13,6E15.7)') dout(1,j),dout((i-1)*6+2,j), &
				dout((i-1)*6+3,j),-dout((i-1)*6+6,j),-dout((i-1)*6+7,j), &
				dout((i-1)*6+4,j),dout((i-1)*6+5,j)
			enddo
			close(51)
		enddo
	endif
end subroutine output_offfault_st

!#3
subroutine output_frt
	! The subroutine output_frt generates frt.txt* files for each MPI process.
	! frt.txt* contain on-fault variables for visualization 
	!   and restart files for the next deformation phase. 
	use globalvar
	implicit none
	
	integer (kind = 4) :: i, j 	
	
	if(nftnd(1) > 0) then
		open(unit=10004+me,file='frt.txt'//mm,status='unknown')
		
		write(10004+me,'(1x,22e18.7e4)')    &
				! 3 coordinates of the fault nodes.
			((x(j,nsmp(1,i,1)), j = 1,3), & 
				! rupture time 
			fnft(i,1),                    & 
				! 71-73: final slips, slipd, slipn
				! 74-76: final sliprates, sliprated, slipraten
			(fric(j,i,1), j = 71,76),     &
				! 47: final slip rate
			fric(47,i,1),                 &
				! 78: final effective normal stress, tnrm 
				! 79: final shear strike stress, tstk
				! 80: final shear dip stress, tdip 
			fric(78,i,1),                 &
			fric(79,i,1),                 &
			fric(80,i,1),                 &
				! 31-33, vxm, vym, vzm, 3 vel components of master nodes.
				! 34-36, vxs, vys, vzs, 3 vel components of slave nodes.
			fric(31,i,1),                 &
			fric(32,i,1),                 &
			fric(33,i,1),                 &
			fric(34,i,1),                 &
			fric(35,i,1),                 &
			fric(36,i,1),                 &
				! 20: state variable in RSF
			fric(20,i,1),                 &
				! 21: state variable for normal stress variation (Shi and Day)
			fric(23,i,1),                 &
				!
			i=1,nftnd(1)) ! Finish the write(10004,me, ...) line.
			
		close(10004+me)
	endif
end subroutine output_frt

!#4
subroutine output_timeanalysis

	use globalvar
	implicit none
	
	integer (kind = 4) :: i, j 	
	
	open(unit=14,file='timeinfo'//mm,status='unknown')	!rupture time
		write(14,'(1x,10e18.7e4,2i10)') (timeused(i),i=1,9),btime,totalNumOfElements,totalNumOfEquations
	close(14)
end subroutine output_timeanalysis

!#5
subroutine output_plastic_strain
	use globalvar
	implicit none
	integer (kind = 4) :: i, j 	
	real (kind = dp) :: sc(3)
	if (output_plastic == 1) then	
		
		do i=1,totalNumOfElements
			if ((pstrain(i)>1.0d-4).and.(abs(x(1,ien(1,i)))<5.0d3).and.(abs(x(2,ien(1,i)))<2.0d3).and.(abs(x(3,ien(1,i)))<8.0d3)) then 
				open(unit=10007+me,file='pstr.txt'//mm,status='unknown',position='append')
				sc=0.0d0
				do j=1,8
					sc(1)=sc(1)+x(1,ien(j,i))
					sc(2)=sc(2)+x(2,ien(j,i))
					sc(3)=sc(3)+x(3,ien(j,i))
				enddo
				sc(1)=sc(1)/8.0d0
				sc(2)=sc(2)/8.0d0
				sc(3)=sc(3)/8.0d0			
				write(10007+me,'(1x,16e18.7e4)') sc(1),sc(2),sc(3),pstrain(i),(s1(stressCompIndexArr(i)+j),j=1,12)	
			endif
		enddo
	endif
end subroutine output_plastic_strain

!#6
subroutine find_surface_node_id
	use globalvar
	implicit none
	integer (kind = 4) :: i, j 	
	real (kind = dp) :: sc(3)
	if (outputGroundMotion == 1) then	
		do i=1,numnp
			if ((x(1,i)<fltxyz(2,1,1)+20.0d3) .and. (x(1,i)>fltxyz(1,1,1)-20.0d3) .and. (abs(x(2,i))<20.0d3) .and. (abs(x(3,i))<dx/1000)) then 
				surface_nnode = surface_nnode + 1
				surface_node_id(surface_nnode) = i
				open(unit=10008+me,file='surface_coor.txt'//mm,status='unknown',position='append')		
					write(10008+me,'(1x,3e18.7e4)') x(1,i),x(2,i),x(3,i)	
			endif
		enddo
	endif
end subroutine find_surface_node_id

!#7
subroutine output_gm
	use globalvar
	implicit none
	integer (kind = 4) :: i, j, itag

	! if (output_ground_motion == 1 .and. surface_nnode > 0) then	
		! open(unit=1009,file='gm.txt'//mm,status='unknown',position='append')		
			! do i=1,surface_nnode
				! itag = surface_node_id(i)
				! write(1009,'(1x,3e15.7)') v(1,itag), v(2,itag), v(3,itag)
		! enddo
	! endif
	if (outputGroundMotion == 1 .and. surface_nnode > 0) then	
		open(unit=10009+me,file='gm'//mm,status='unknown',position='append', access='stream')		
			do i=1,surface_nnode
				itag = surface_node_id(i)
				write(10009+me) v(1,itag), v(2,itag), v(3,itag)
			enddo
	endif	
end subroutine output_gm
