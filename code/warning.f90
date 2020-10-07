subroutine warning

	use globalvar
	implicit none
	
	if (C_elastic==0.and.C_Q==1) then
		write(*,*) 'Q model can only work with elastic code'
		stop 1001 
	endif
	if (C_Q==1.and.rat>1) then
		write(*,*) 'Q model can only work with uniform element size'
		write(*,*) 'rat should be 1.0'
		stop 1002
	endif	

end subroutine warning