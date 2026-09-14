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
    integer (kind = 4) :: nranks
    
    if (me == 0) then
        ! 1.54 GB per million elements, measured on Ubuntu docker for
        ! test.drv.a6 with 4 cores. totalNumOfElements is THIS rank's count,
        ! so the job total is that times the rank count. Both are reported
        ! because sizing a job needs the total while judging per-step cost
        ! needs the per-rank figure (performance scales with cells per rank).
        nranks = npx*npy*npz
        write(*,'(A,I12)')    ' Cells per rank (rank 0)     : ', totalNumOfElements
        write(*,'(A,I12)')    ' Cells total (all ranks)     : ', totalNumOfElements*nranks
        write(*,'(A,I12)')    ' MPI ranks (npx*npy*npz)     : ', nranks
        write(*,'(A,F12.3,A)')' Memory per rank             : ', 1.54d0*totalNumOfElements/1.0d6, ' GB'
        write(*,'(A,F12.3,A)')' Memory total (all ranks)    : ', 1.54d0*totalNumOfElements/1.0d6*nranks, ' GB'
    endif
end subroutine memory_estimate
