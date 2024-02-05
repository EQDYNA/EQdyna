! Table of Contents of functions and subroutines. 
! #
! # vlm
! -----------------------------------------------
subroutine fb1(xtmp,ww,wtmp,res)
    use globalvar
    implicit none
    real (kind = dp) :: xtmp, ww, wtmp, res
    ! ww == W;
    if (abs(xtmp)<=ww) then 
        res = 1.0d0
    elseif (abs(xtmp)>ww.and.abs(xtmp)<ww+wtmp) then 
        res = 0.5d0*(1.0d0 + dtanh(wtmp/(abs(xtmp)-ww-wtmp) + wtmp/(abs(xtmp)-ww)))
    elseif (abs(xtmp)>=ww+wtmp) then 
        res = 0.0d0
    endif
end subroutine fb1
! -----------------------------------------------
subroutine fb2(ytmp,ww,wtmp,res)
    use globalvar
    implicit none
    real (kind = dp) :: ytmp, ww, wtmp, res
    ! ww == W;  
    if (ytmp<0.0d0) then
        write(*,*) 'z coordinates should be positive for B3'
        stop
    endif   
    if (ytmp<=wtmp) then 
        res = 0.5d0*(1.0d0 + dtanh(wtmp/(wtmp-ytmp) - wtmp/ytmp))
    elseif (ytmp>=wtmp.and.ytmp<=ww) then 
        res = 1.0d0 
    elseif (ytmp>ww.and.ytmp<ww+wtmp) then 
        res = 0.5d0*(1.0d0 + dtanh(wtmp/(ytmp-ww-wtmp) + wtmp/(ytmp-ww)))
    elseif (ytmp>=ww+wtmp) then 
        res = 0.0d0
    endif
end subroutine fb2
!------------------------------------------------------
subroutine fb3(ytmp,ww,wtmp,res)
    use globalvar
    implicit none
    real (kind = dp) :: ytmp, ww, wtmp, res
    ! ww == W;  
    if (ytmp<0.0d0) then
        write(*,*) 'z coordinates should be positive for B3'
        stop
    endif 
    if (ytmp<=ww) then 
        res = 1.0d0 
    elseif (ytmp>ww.and.ytmp<ww+wtmp) then 
        res = 0.5d0*(1.0d0 + dtanh(wtmp/(ytmp-ww-wtmp) + wtmp/(ytmp-ww)))
    elseif (ytmp>=ww+wtmp) then 
        res = 0.0d0
    endif
end subroutine fb3
! -----------------------------------------------
subroutine vlm(xl,volume)
    use globalvar
    implicit none
    !
    !...program to calculate volume of a hexlahedron from
    ! nodal coordinates. See Belytschko et al.(1984) for 
    ! reference.
    ! B.D. 8/21/05
    !
    integer(kind=4) :: i
    real(kind = dp) :: volume
    real(kind = dp),dimension(nesd,nen) :: xl
    integer(kind=4),dimension(8,8) :: it = reshape((/ &
        1,2,3,4,5,6,7,8, 2,3,4,1,6,7,8,5, 3,4,1,2,7,8,5,6, &
        4,1,2,3,8,5,6,7, 5,8,7,6,1,4,3,2, 6,5,8,7,2,1,4,3, &
        7,6,5,8,3,2,1,4, 8,7,6,5,4,3,2,1/),(/8,8/))
    real(kind = dp),dimension(8) :: bb
    !
    do i=1,8
        bb(i) = xl(2,it(2,i))*(xl(3,it(6,i))-xl(3,it(3,i))+xl(3,it(5,i)) &
            -xl(3,it(4,i))) + xl(2,it(3,i))*(xl(3,it(2,i))-xl(3,it(4,i))) &
            + xl(2,it(4,i))*(xl(3,it(3,i))-xl(3,it(8,i))+xl(3,it(2,i)) &
            -xl(3,it(5,i))) + xl(2,it(5,i))*(xl(3,it(8,i))-xl(3,it(6,i)) &
            +xl(3,it(4,i))-xl(3,it(2,i))) + xl(2,it(6,i))*(xl(3,it(5,i)) &
            -xl(3,it(2,i))) + xl(2,it(8,i))*(xl(3,it(4,i))-xl(3,it(5,i)))
    enddo
    !
    volume = 0.0d0
    do i=1,nen
        volume = volume + xl(1,i) * bb(i)
    enddo
    volume = volume/12.0d0
    !  
end subroutine vlm

subroutine memory_estimate
    use globalvar
    implicit none
    
    real(kind = dp) :: memory = 0.0d0 ! in bytes
    
    if (me == 0) then
        !write(*,*) memory, 'GB memory would be used on me=0'
        ! An estimate of 2GB/million elements memory is needed; 
        ! Test done on Ubuntu docker for test.drv.a6 with 4 cores;
        write(*,*) 1.54d0*numel/1.0e6*npx*npy*npz, ' GB memory is expected for EQdyna ... ...'
    endif
end subroutine memory_estimate
