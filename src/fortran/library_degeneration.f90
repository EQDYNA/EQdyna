! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine wedge(cenx, ceny, cenz, elemCount, stressDofCount, iy, iz, nftndtmp) 
    
    use globalvar
    implicit none
    
    integer (kind = 4) :: elemCount, stressDofCount, iy, iz, nftndtmp, k, i, neworder(nen)
    real (kind = dp) :: cenx, ceny, cenz, pointToFaultDist, tangentDip
        tangentDip = dtan(C_degen/180.d0*pi)
        pointToFaultDist = abs(ceny*tangentDip+cenz)/(1.d0+tangentDip**2)**0.5
    if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and. &
            ceny>fltxyz(1,2,1).and.ceny<fltxyz(2,2,1).and. &
            cenz>fltxyz(1,3,1).and. &
            pointToFaultDist<tol) then
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

subroutine wedge4num(cenx, ceny, cenz, elemCount) 
    
    use globalvar
    implicit none
    
    integer (kind = 4) :: elemCount    
    real (kind = dp) :: cenx, ceny, cenz, tangentDip, pointToFaultDist
    
        tangentDip = dtan(C_degen/180.d0*pi)
        pointToFaultDist = abs(ceny*tangentDip+cenz)/(1.d0+tangentDip**2)**0.5

    if (cenx>fltxyz(1,1,1).and.cenx<fltxyz(2,1,1).and. &
            ceny>fltxyz(1,2,1).and.ceny<fltxyz(2,2,1).and. &
            cenz>fltxyz(1,3,1).and. &
            pointToFaultDist<dx/100.d0) then

        elemCount = elemCount + 1
    endif 
end subroutine

subroutine reorder(neworder, elemCount, iy, iz)

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, neworder(nen), nodeIdPerElem(nen), elemCount, iz, iy
    
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
