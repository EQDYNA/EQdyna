! subroutine Read_Input_Files(C_elastic, C_nuclea, friclaw, ntotft, nucfault, &
							! npx, npy, npz, &
							! term, dt, &
							! xmin, xmax, ymin, ymax, zmin, zmax, &
							! dis4uniF, dis4uniB, rat, dx, 
							! fxmin, fxmax, fymin, fymax, fzmin, fzmax, &
							! nmat, n2mat material, &
							! me)
! ! This subroutine is to read controllable parameters and model information from input files. 
! ! Input files include 

	! use globalvar
	! implicit none
	! include 'mpif.h'
	! integer(kind=4)::C_elastic, C_nuclea, friclaw, ntotft, nucfault, &
								! npx, npy, npz, &
								! nmat, n2mat, &
								! me
	! real(kind=8):: term, dt, &
					! xmin, xmax, ymin, ymax, zmin, zmax, &
					! dis4uniF, dis4uniB, rat, dx
	! real(kind=8),allocatable:: fxmin, fxmax, fymin, fymax, fzmin, fzmax, &
								! material
								
	! allocate(fxmin(ntotft),fxmax(ntotft),fymin(ntotft),fymax(ntotft),fzmin(ntotft),fzmax(ntotft),material(nmat,n2mat))

	! call readglobal(C_elastic, C_nuclea, friclaw, ntotft, nucfault, &
								! npx, npy, npz, &
								! term, dt, me)
	! call readmodelgeometry(xmin, xmax, ymin, ymax, zmin, zmax, &
								! dis4uniF, dis4uniB, rat, dx, me)	
	! call readfaultgeometry(fxmin, fxmax, fymin, fymax, fzmin, fzmax, me)
	! call readmaterial(material, me)

! end subroutine Read_Input_Files
! #1 readgloabl -----------------------------------------------------------
subroutine readglobal(me)
! This subroutine is read information from FE_global.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: me
	
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
		write(*,*) 'Checking FE_global.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Global.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Global.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1001, file = 'FE_Global.txt', form = 'formatted', status = 'old')
		read(1001,*) C_elastic
		read(1001,*) C_nuclea
		read(1001,*) friclaw
		read(1001,*) ntotft
		read(1001,*) nucfault
		read(1001,*) 
		read(1001,*) npx, npy, npz
		read(1001,*)
		read(1001,*) term
		read(1001,*) dt 
		read(1001,*)
		read(1001,*) nmat, n2mat
	close(1001)
end subroutine readglobal 
! #2 readmodelgeometry -------------------------------------------------
subroutine readmodelgeometry(me)
! This subroutine is read information from FE_global.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: me
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Model_Geometry.txt", EXIST=file_exists)
		write(*,*) 'Checking FE_Model_Geometry.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Model_Geometry.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Model_Geometry.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Model_Geometry.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1002, file = 'FE_Model_Geometry.txt', form = 'formatted', status = 'old')
		read(1002,*) xmin, xmax
		read(1002,*) ymin, ymax
		read(1002,*) zmin, zmax
		read(1002,*) 
		read(1002,*) dis4uniF, dis4uniB
		read(1002,*) rat
		read(1002,*) dx 
	close(1002)
end subroutine readmodelgeometry

! #3 readfaultgeometry
subroutine readfaultgeometry(me)
! This subroutine is read information from FE_Fault_Geometry.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4)::me,i
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Fault_Geometry.txt", EXIST=file_exists)
		write(*,*) 'Checking FE_Fault_Geometry.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Fault_Geometry.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Fault_Geometry.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Fault_Geometry.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1003, file = 'FE_Fault_Geometry.txt', form = 'formatted', status = 'old')
		do i = 1, ntotft
			read(1003,*) 
			read(1003,*) fxmin(i), fxmax(i)
			read(1003,*) fymin(i), fymax(i)
			read(1003,*) fzmin(i), fzmax(i)
		enddo
	close(1003)
end subroutine readfaultgeometry

! #4 readmaterial --------------------------------------------------------
subroutine readmaterial(me)
! This subroutine is read information from FE_Material.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: me, i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Material.txt", EXIST=file_exists)
		write(*,*) 'Checking FE_Material.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Material.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Material.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Material.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1004, file = 'FE_Material.txt', form = 'formatted', status = 'old')
		do i = 1, nmat
			read(1004,*) (material(i,j), j = 1, n2mat)
		enddo 
	close(1002)
end subroutine readmaterial
