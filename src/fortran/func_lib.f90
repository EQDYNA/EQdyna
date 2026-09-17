! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT

! Shared helper subroutines used by more than one caller across src/.
! Table of Contents of functions and subroutines.
! #1 pmlRegionDistance
! #2 insertFaultInterface
! #3 fb1
! #4 fb2
! #5 fb3
! #6 vlm
! #7 checkPMLAlignment
! #8 devStrDepthTaper

! #1 pmlRegionDistance
subroutine pmlRegionDistance(x, y, z, xmax0, xmin0, ymax0, ymin0, zmin0, damp)
    ! Classifies a point (x,y,z) into the PML region cascade and returns the
    ! raw (undamped) distances damp(1:3) to the xmax0/xmin0, ymax0/ymin0, and
    ! zmin0 boundaries. Consolidates the four identical region-cascade copies
    ! previously duplicated in computePMLDampingVector.f90 (twice) and
    ! assembleGlobalKU.f90's calcPMLElemKU (twice).
    !
    ! Boundary tests are INCLUSIVE (>=/<=), standardized 2026-09-13
    ! (pathway_forward item 13): a point exactly on a PML bound classifies
    ! into the region with zero damping distance (harmless) instead of
    ! falling out of the cascade. Element centers additionally may never
    ! lie exactly on a PML bound -- enforced at mesh time by
    ! checkPMLAlignment (called from meshgen).
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: x, y, z, xmax0, xmin0, ymax0, ymin0, zmin0
    real (kind = dp), dimension(3) :: damp
    logical :: xHi, xLo, yHi, yLo

    xHi = (x>=xmax0)
    xLo = (x<=xmin0)
    yHi = (y>=ymax0)
    yLo = (y<=ymin0)

    if (z<=zmin0) then
        damp(3) = abs(z-zmin0)
    else
        damp(3) = 0.0d0
    endif

    if (xHi.and.yHi) then !region 11
        damp(1)=abs(x-xmax0)
        damp(2)=abs(y-ymax0)
    elseif (xHi.and.yLo) then !region 12
        damp(1)=abs(x-xmax0)
        damp(2)=abs(y-ymin0)
    elseif (xLo.and.yLo) then !region 13
        damp(1)=abs(x-xmin0)
        damp(2)=abs(y-ymin0)
    elseif (xLo.and.yHi) then !region 14
        damp(1)=abs(x-xmin0)
        damp(2)=abs(y-ymax0)
    elseif (xHi.and.y>ymin0.and.y<ymax0) then !region 1_12
        damp(1)=abs(x-xmax0)
        damp(2)=0.0d0
    elseif (yLo.and.x>xmin0.and.x<xmax0) then !region 1_23
        damp(1)=0.0d0
        damp(2)=abs(y-ymin0)
    elseif (xLo.and.y>ymin0.and.y<ymax0) then !region 1_34
        damp(1)=abs(x-xmin0)
        damp(2)=0.0d0
    elseif (yHi.and.x>xmin0.and.x<xmax0) then !region 1_41
        damp(1)=0.0d0
        damp(2)=abs(y-ymax0)
    else
    !Middle area 9 missing previously.
    !Feb.18.2016/D.Liu
        damp(1)=0.0d0
        damp(2)=0.0d0
    endif
end subroutine pmlRegionDistance

! #2 insertFaultInterface
subroutine insertFaultInterface(nodeCoor, ycoort, pfx, pfz)
    ! This subroutine is to modify ycoor to insert a rough/& dipping fault interface.
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: nodeCoor(10), peak, ycoort, pfx, pfz, fx1, fx2, fz1
    integer (kind = 4) :: ixx, izz
    
    fx1 = rough_fx_min
    fx2 = rough_fx_max
    fz1 = rough_fz_min
    ! Index (ixx, izz) are counted from the fault corner (rough_fx_min, rough_fz_min)
    if ((nodeCoor(1) < fx2 + tol) .and. (nodeCoor(1) > fx1 - tol) .and. (nodeCoor(3) > fz1 - tol)) then 
        ixx = nint((nodeCoor(1) - fx1)/dx) + 1
        izz = nint((nodeCoor(3) - fz1)/dz) + 1
    elseif ((nodeCoor(1) < fx1 - tol) .and. (nodeCoor(3) > fz1 - tol) ) then
        ixx = 1
        izz = nint((nodeCoor(3) - fz1)/dz) + 1
    elseif ((nodeCoor(1) > fx2 + tol) .and. (nodeCoor(3) > fz1 - tol)) then 
        ixx = nnx
        izz = nint((nodeCoor(3) - fz1)/dz) + 1
    elseif ((nodeCoor(1) < fx2 + tol) .and. (nodeCoor(1) > fx1 - tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = nint((nodeCoor(1) - fx1)/dx) + 1
        izz = 1
    elseif ((nodeCoor(1) < fx1 - tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = 1
        izz = 1 
    elseif ((nodeCoor(1) > fx2 + tol) .and. (nodeCoor(3) < fz1 - tol)) then 
        ixx = nnx
        izz = 1
    endif 
    
    peak = rough_geo(1,nnz*(ixx-1)+izz)
    pfx  = rough_geo(2,nnz*(ixx-1)+izz)
    pfz  = rough_geo(3,nnz*(ixx-1)+izz)    
    
    if (nodeCoor(2) > -tol) then
        ycoort = nodeCoor(2)*(ymax - peak)/ymax + peak
    elseif (nodeCoor(2) < -tol) then 
        ycoort = nodeCoor(2)*(peak - ymin)/(-ymin) + peak 
    endif 
    
end subroutine insertFaultInterface

! #3 fb1
subroutine fb1(xtmp,ww,wtmp,res)
    use globalvar
    use errorCodes
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

! #4 fb2
subroutine fb2(ytmp,ww,wtmp,res)
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: ytmp, ww, wtmp, res
    ! ww == W;  
    if (ytmp<0.0d0) then
        call abortRun(ERR_NUM_NEGATIVE_DEPTH, &
            'Negative depth passed to fb2; z coordinates must be positive here.')
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

! #5 fb3
subroutine fb3(ytmp,ww,wtmp,res)
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: ytmp, ww, wtmp, res
    ! ww == W;  
    if (ytmp<0.0d0) then
        call abortRun(ERR_NUM_NEGATIVE_DEPTH, &
            'Negative depth passed to fb3; z coordinates must be positive here.')
    endif 
    if (ytmp<=ww) then 
        res = 1.0d0 
    elseif (ytmp>ww.and.ytmp<ww+wtmp) then 
        res = 0.5d0*(1.0d0 + dtanh(wtmp/(ytmp-ww-wtmp) + wtmp/(ytmp-ww)))
    elseif (ytmp>=ww+wtmp) then 
        res = 0.0d0
    endif
end subroutine fb3

! #6 vlm
subroutine vlm(xl,volume)
    use globalvar
    use errorCodes
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

subroutine checkPMLAlignment(elemCenter)
    ! Mesh-time precheck (PROJECT_RULES.md rule 2; pathway_forward item 13):
    ! the PML region cascade assumes no element center lies exactly on a PML
    ! boundary plane. This held by construction on structured meshes but was
    ! never enforced; a stretched or degenerated mesh could violate it
    ! silently. Fail loudly at mesh time instead.
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: elemCenter(3)

    if (nPML <= 0) return
    if (elemCenter(1) == PMLb(1) .or. elemCenter(1) == PMLb(2) .or. &
        elemCenter(2) == PMLb(3) .or. elemCenter(2) == PMLb(4) .or. &
        elemCenter(3) == PMLb(5)) then
        write(*,*) 'checkPMLAlignment: element center exactly on a PML bound at', &
            elemCenter(1), elemCenter(2), elemCenter(3)
        call abortRun(ERR_NUM_PML_ALIGNMENT, &
            'An element centre lies exactly on a PML bound, making region classification ambiguous. Adjust the mesh or nPML.')
    endif
end subroutine checkPMLAlignment

! #8 devStrDepthTaper
function devStrDepthTaper(depth) result(omega)
    ! Depth taper applied to the OFF-FAULT deviatoric pre-stress built by
    ! meshgen.f90's setPlasticStress. Mirrors SCEC TPV29/TPV30's Omega(depth)
    ! (spec part 4, "Initial Stress Tensor"): the deviatoric component of the
    ! initial stress tapers linearly to zero between two depths, while the
    ! vertical (lithostatic) component keeps growing. TPV30's own values are
    ! 17000 m and 22000 m.
    !
    ! Without it, devStr stays a FIXED fraction of |strVert| at every depth,
    ! so the deviatoric stress grows without bound with depth -- fine for the
    ! cases gated before v5.9.0 (all of which leave the taper inactive and are
    ! bit-for-bit unchanged), wrong for TPV30.
    !
    ! The taper is INACTIVE when devStrTaperDepthEnd <= devStrTaperDepthStart
    ! (case.setup writes 0.0 0.0 for a case that does not configure it, and
    ! refuses a half-configured pair), and then returns exactly 1.0d0 -- an
    ! IEEE-exact multiplicative identity, which is what keeps the pre-taper
    ! result bit-for-bit.
    !
    ! depth is positive DOWN, in m -- the same argument setPlasticStress
    ! already takes.
    use globalvar
    implicit none
    real (kind = dp) :: depth, omega

    omega = 1.0d0
    if (devStrTaperDepthEnd > devStrTaperDepthStart) then
        if (depth >= devStrTaperDepthEnd) then
            omega = 0.0d0
        elseif (depth > devStrTaperDepthStart) then
            omega = (devStrTaperDepthEnd - depth) &
                    /(devStrTaperDepthEnd - devStrTaperDepthStart)
        endif
    endif
end function devStrDepthTaper
