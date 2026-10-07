! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine wedge(cenx, ceny, cenz, elemCount, stressDofCount, iy, iz, nftndtmp, iFault)
    ! Row 153 checkpoint 1: iFault added (was hardcoded "1" throughout this
    ! subroutine, via fltxyz(...,1) and bare C_degen); every existing caller
    ! passes iFault=1, so this is a pure refactor (see meshgen.f90 call site).
    use globalvar
    implicit none

    integer (kind = 8) :: elemCount, stressDofCount
    integer (kind = 4) :: iy, iz, nftndtmp, k, i, neworder(nen), iFault
    real (kind = dp) :: cenx, ceny, cenz, pointToFaultDist, tangentDip, x0, pNorm
    logical :: onDegenPlane
        onDegenPlane = .false.
        if (faultDegenStyle(iFault) == 1) then
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            pointToFaultDist = abs(ceny*tangentDip+cenz)/(1.d0+tangentDip**2)**0.5
            if (cenx>fltxyz(1,1,iFault).and.cenx<fltxyz(2,1,iFault).and. &
                    ceny>fltxyz(1,2,iFault).and.ceny<fltxyz(2,2,iFault).and. &
                    cenz>fltxyz(1,3,iFault).and. &
                    pointToFaultDist<tol) onDegenPlane = .true.
        elseif (faultDegenStyle(iFault) == 2) then
            ! Row 153 checkpoint 2a: style 2 is the TPV24/25 branch fault --
            ! a strike-tilted plane in x-y, extruded along z (vertical fault,
            ! dipping along-strike rather than in dip). x0 is the junction's
            ! x-coordinate, derived from this fault's own box lower edge
            ! (fltxyz(1,1,iFault) = x0 + dx, one cell gap so the junction
            ! column itself is excluded -- the same convention
            ! probe_branch_mesh.f90 (PR #143) used and the owner's "junction
            ! node owned by the main fault only" decision requires). Verified
            ! against that probe: neworder11/12 and the (cenx,ceny) distance
            ! test below are exactly its derived/confirmed values.
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            pNorm = (1.d0+tangentDip**2)**0.5
            x0 = fltxyz(1,1,iFault) - dx
            pointToFaultDist = abs(ceny+(cenx-x0)*tangentDip)/pNorm
            if (cenx>fltxyz(1,1,iFault).and.cenx<fltxyz(2,1,iFault).and. &
                    ceny>fltxyz(1,2,iFault).and.ceny<fltxyz(2,2,iFault).and. &
                    cenz>fltxyz(1,3,iFault).and. &
                    pointToFaultDist<tol) onDegenPlane = .true.
        endif
    if (onDegenPlane) then
        if (faultDegenStyle(iFault) == 1) then
        ! Degenerate the brick element into two wedge elements.
        !       8
        !  5        7
        !      6         ! Brick
        !       4
        !  1        3
        !      2
        ! Collapse 4->3, 8->3
        !       4
        !  1        8
        !      5
        !       3
        !  2        7
        !      6
        elemTypeArr(elemCount) = 11 ! wedge below fault
        neworder = (/5,1,4,4,6,2,3,3/)
        call reorder(neworder, elemCount, iy, iz)
        ! Collapse 4->3, 8->3
        !       2
        !  3        6
        !      7
        !       1
        !  4        5
        !      8
        elemCount = elemCount + 1
        elemTypeArr(elemCount) = 12
        neworder = (/4,8,5,5,3,7,6,6/)
        call reorder(neworder, elemCount, iy, iz)
        else
        ! Style 2 (x-y tilt, z-extruded): neworder11/neworder12, verified
        ! by probe_branch_mesh.f90 (positive Jacobian on every wedge
        ! element it builds).
        elemTypeArr(elemCount) = 11
        neworder = (/8,5,6,6,4,1,2,2/)
        call reorder(neworder, elemCount, iy, iz)
        elemCount = elemCount + 1
        elemTypeArr(elemCount) = 12
        neworder = (/6,7,8,8,2,3,4,4/)
        call reorder(neworder, elemCount, iy, iz)
        endif

        ! The above node uses existing array spaces for stress and mat.
        stressCompIndexArr(elemCount)=stressDofCount
        stressDofCount = stressDofCount + 12

        mat(elemCount,1)=material(1,1)
        mat(elemCount,2)=material(1,2)
        mat(elemCount,3)=material(1,3)                
        mat(elemCount,5)=mat(elemCount,2)**2*mat(elemCount,3)!miu=vs**2*rho
        mat(elemCount,4)=mat(elemCount,1)**2*mat(elemCount,3)-2*mat(elemCount,5)!lam=vp**2*rho-2*miu    
    endif                 
end subroutine

subroutine wedge4num(cenx, ceny, cenz, elemCount, iFault)
    ! Row 153 checkpoint 1: iFault added (see wedge() above); every existing
    ! caller passes iFault=1, so this is a pure refactor.
    use globalvar
    implicit none

    integer (kind = 8) :: elemCount
    integer (kind = 4) :: iFault
    real (kind = dp) :: cenx, ceny, cenz, tangentDip, pointToFaultDist, x0, pNorm
    logical :: onDegenPlane

        onDegenPlane = .false.
        if (faultDegenStyle(iFault) == 1) then
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            pointToFaultDist = abs(ceny*tangentDip+cenz)/(1.d0+tangentDip**2)**0.5
            if (cenx>fltxyz(1,1,iFault).and.cenx<fltxyz(2,1,iFault).and. &
                    ceny>fltxyz(1,2,iFault).and.ceny<fltxyz(2,2,iFault).and. &
                    cenz>fltxyz(1,3,iFault).and. &
                    pointToFaultDist<dx/100.d0) onDegenPlane = .true.
        elseif (faultDegenStyle(iFault) == 2) then
            ! Mirrors wedge()'s style-2 test above -- same x0/pNorm derivation.
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            pNorm = (1.d0+tangentDip**2)**0.5
            x0 = fltxyz(1,1,iFault) - dx
            pointToFaultDist = abs(ceny+(cenx-x0)*tangentDip)/pNorm
            if (cenx>fltxyz(1,1,iFault).and.cenx<fltxyz(2,1,iFault).and. &
                    ceny>fltxyz(1,2,iFault).and.ceny<fltxyz(2,2,iFault).and. &
                    cenz>fltxyz(1,3,iFault).and. &
                    pointToFaultDist<dx/100.d0) onDegenPlane = .true.
        endif

    if (onDegenPlane) then
        elemCount = elemCount + 1
    endif
end subroutine

subroutine reorder(neworder, elemCount, iy, iz)

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, neworder(nen), iz, iy
    integer (kind = 8) :: nodeIdPerElem(nen), elemCount
    
        nodeIdPerElem(1) = plane1(iy-1,iz-1)
        nodeIdPerElem(2) = plane2(iy-1,iz-1)
        nodeIdPerElem(3) = plane2(iy,iz-1)
        nodeIdPerElem(4) = plane1(iy,iz-1)
        nodeIdPerElem(5) = plane1(iy-1,iz)
        nodeIdPerElem(6) = plane2(iy-1,iz)
        nodeIdPerElem(7) = plane2(iy,iz)
        nodeIdPerElem(8) = plane1(iy,iz)    
        
        do i = 1, nen
            nodeElemIdRelation(i,elemCount) = nodeIdPerElem(neworder(i))
        enddo 
        
end subroutine reorder
