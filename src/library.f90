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
    volume = 0.0
    do i=1,nen
        volume = volume + xl(1,i) * bb(i)
    enddo
    volume = volume/12.0
    !  
end subroutine vlm
  
subroutine insert_rough_fault(xcoor, ycoor, zcoor, ycoort, pfx, pfz)
    ! This subroutine is to modify ycoor if a rough_fault interface is inserted.
    
    use globalvar
    implicit none
    real (kind = dp) :: xcoor, ycoor, zcoor, peak, ycoort, pfx, pfz
    real (kind = dp) :: fx1, fx2, fz1
    integer (kind = 4) :: ixx, izz
    
    fx1 = rough_fx_min
    fx2 = rough_fx_max
    fz1 = rough_fz_min
    if ((xcoor < fx2 + tol) .and. (xcoor > fx1 - tol) .and. (zcoor > fz1 - tol)) then 
        ixx = (xcoor - fx1)/dx + 1
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor < fx1 - tol) .and. (zcoor > fz1 - tol) ) then
        ixx = 1
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor > fx2 + tol) .and. (zcoor > fz1 - tol)) then 
        ixx = nnx
        izz = (zcoor - fz1)/dx + 1
    elseif ((xcoor < fx2 + tol) .and. (xcoor > fx1 - tol) .and. (zcoor < fz1 - tol)) then 
        ixx = (xcoor - fx1)/dx + 1
        izz = 1
    elseif ((xcoor < fx1 - tol) .and. (zcoor < fz1 - tol)) then 
        ixx = 1
        izz = 1 
    elseif ((xcoor > fx2 + tol) .and. (zcoor < fz1 - tol)) then 
        ixx = nnx
        izz = 1
    endif 
    
    peak = rough_geo(1,nnz*(ixx-1)+izz)
    pfx = rough_geo(2,nnz*(ixx-1)+izz)
    pfz = rough_geo(3,nnz*(ixx-1)+izz)    
    
    if (ycoor > -tol) then
        ycoort = ycoor*(ymax - peak)/ymax + peak
    elseif (ycoor < -tol) then 
        ycoort = ycoor*(peak - ymin)/(-ymin) + peak 
    endif 
    
end subroutine insert_rough_fault

subroutine memory_estimate
    use globalvar
    implicit none
    
    real(kind = dp) :: memory = 0.0d0 ! in bytes
    
    memory = memory + 4*(maxm+(5+2*ndof)*numnp)
    memory = memory + 4*(nen+2*nee+1+2*14+nen*4+2*(nrowsh-1)*nen)*numel
    memory = memory + 4*(2+120*2)*nftmx*ntotft
    memory = memory + 4*numel + 8*5*maxm
    memory = memory + 8*4*neq + 8*2*ndof*numnp
    memory = memory/1024/1024
    
    if (me == 0) then
        write(*,*) memory, 'GB memory would be used on me=0'
    endif
end subroutine memory_estimate
    
subroutine plastic_set_mat_stress(depthz, nelement)
    use globalvar
    implicit none
    
    real(kind = dp) :: depthz, vstmp, vptmp, routmp, tmp2, norm_thres ! in positive meters
    integer(kind = 4) :: nelement
    
    depthz = depthz/1.0d3 ! in km.
    norm_thres = -100.0d6
    
    if (depthz<0.03d0) then 
        vstmp = 2.206d0*depthz**0.272
    elseif (depthz<0.19d0) then 
        vstmp = 3.542d0*depthz**0.407
    elseif (depthz<4.0d0) then 
        vstmp = 2.505d0*depthz**0.199
    elseif (depthz<8.0d0) then 
        vstmp = 2.927d0*depthz**0.086
    else
        vstmp = 2.927d0*8.0d0**0.086
    endif 
    vptmp = max(1.4d0+1.14d0*vstmp, 1.68d0*vstmp)
    routmp = 2.4405d0 + 0.10271d0*vstmp
    !roumax = 2.4405d0 + 0.10271d0*2.927d0*8.0d0**0.086
    vstmp = vstmp*1.0d3
    vptmp = vptmp*1.0d3
    routmp = routmp*1.0d3
    !roumax = roumax*1.0d3
    
    mat(nelement,1)=vptmp
    mat(nelement,2)=vstmp
    mat(nelement,3)=routmp                
    mat(nelement,5)=mat(nelement,2)**2*mat(nelement,3)!miu=vs**2*rho
    mat(nelement,4)=mat(nelement,1)**2*mat(nelement,3)-2*mat(nelement,5)!lam=vp**2*rho-2*miu    
    
    depthz = depthz*1.0d3 
    !tmp2 = min(depthz, 5.0d3) * grav
    tmp2 = depthz * grav
    
    if(et(nelement)==1)then
        eleporep(nelement)= 0.0d0!rhow*tmp2*gama  !pore pressure>0
        s1(ids(nelement)+3)=-(roumax- rhow*(gamar+1.0d0))*tmp2  !vertical, comp<0    
        s1(ids(nelement)+1)=s1(ids(nelement)+3)
        s1(ids(nelement)+2)=s1(ids(nelement)+3) 
        s1(ids(nelement)+6)=-s1(ids(nelement)+3) *0.33d0
        ! if (s1(ids(nelement)+3)<norm_thres) then 
            ! s1(ids(nelement)+3) = norm_thres
            ! s1(ids(nelement)+1) = s1(ids(nelement)+3)
            ! s1(ids(nelement)+2) = s1(ids(nelement)+3) 
            ! s1(ids(nelement)+6) = (roumax -rhow)*0.33d0*tmp2*abs(norm_thres/(roumax*tmp2))
        ! endif         
    elseif(et(nelement)==2)then
        eleporep(nelement)= 0.0d0!rhow*tmp2*gama  !pore pressure>0
        s1(ids(nelement)+3+15)=-(roumax-rhow*(gamar+1.0d0))*tmp2  !vertical, comp<0
        s1(ids(nelement)+1+15)=s1(ids(nelement)+3+15)
        s1(ids(nelement)+2+15)=s1(ids(nelement)+3+15) 
        s1(ids(nelement)+6+15)=-s1(ids(nelement)+3+15)*0.33d0        
        ! if (s1(ids(nelement)+3+15)<norm_thres) then 
            ! s1(ids(nelement)+3+15) = norm_thres
            ! s1(ids(nelement)+1+15) = s1(ids(nelement)+3+15)
            ! s1(ids(nelement)+2+15) = s1(ids(nelement)+3+15) 
            ! s1(ids(nelement)+6+15) = (roumax -rhow)*0.33d0*tmp2*abs(norm_thres/(roumax*tmp2))
        ! endif         
    endif        
    
end subroutine plastic_set_mat_stress
