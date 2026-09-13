! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT

! Table of Contents of functions and subroutines.
! #
! # memory_estimate
! -----------------------------------------------
! fb1, fb2, fb3, vlm moved to func_lib.f90.

subroutine memory_estimate
    use globalvar
    implicit none
    
    real(kind = dp) :: memory = 0.0d0 ! in bytes
    
    if (me == 0) then
        !write(*,*) memory, 'GB memory would be used on me=0'
        ! An estimate of 2GB/million elements memory is needed; 
        ! Test done on Ubuntu docker for test.drv.a6 with 4 cores;
        write(*,*) 1.54d0*totalNumOfElements/1.0e6*npx*npy*npz, ' GB memory is expected for EQdyna ... ...'
    endif
end subroutine memory_estimate
