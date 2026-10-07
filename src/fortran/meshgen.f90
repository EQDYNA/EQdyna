! Copyright (C) 2006 Benchun Duan <bduan@tamu.edu>, Dunyu Liu <dliu@ig.utexas.edu>
! MIT
subroutine meshgen
    ! Create regular node grids and hexahedral elements for this MPI process.
    ! Create split nodes.
    ! Distort the regular mesh for complex fault geometries.
    ! Set material propertiesa, initial stress, and other element-wise properties. 
    
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    ! incremental variables 
    ! Item 143: counts, ids and offsets are 64-bit (globalvar.f90); the
    ! per-axis grid sizes and loop indices stay kind=4.
    integer (kind = 8) :: nodeCount=0, elemCount=0, equationNumCount=0, eqNumIndexArrLocTag=0, stressDofCount=0
    integer (kind = 8) :: msnode, n1,n2,n3,n4,m1,m2,m3,m4
    integer (kind = 4) :: nxt,nyt,nzt,nx,ny,nz,ix,iy,iz, &
                       edgex1,edgey1,edgez1,i,j,k,i1,j1,k1,edgezn, &
                       nxuni,nyuni,nzuni,ift, &
                       mex,mey,mez,itmp1,&
                       numOfDofPerNodeTmp, nodeXyzIndex(10),isOnFt
    integer (kind = 4), dimension(ntotft) :: nftnd0,ixfi,izfi,ifs,ifd
    integer (kind = 8), allocatable :: fltrc(:,:,:,:)
    ! Temporary real variables
    real (kind = dp) :: nodeCoor(10), elementCenterCoor(3), &
                       a,b,area,aa1,bb1,cc1,dd1,p1,q1, ycoort, pfx = 0.0d0, pfz = 0.0d0
    real (kind = dp) :: xline(10000), yline(10000), zline(10000), modelBoundCoor(3,2)
    ! pathway_forward.md item 94 (owner ruling 2026-09-24): clamp off-fault
    ! station depth to the nearest z-grid node instead of dropping a station
    ! whose requested depth is not exactly a grid z-plane. xGridFull/
    ! yGridFull/zGridFull are the FULL (unsliced) 1D coordinate arrays
    ! getLocalOneDimCoorArrAndSize builds internally before slicing to this
    ! rank's local piece -- identical bit-for-bit on every rank (a pure
    ! function of dz/fltxyz/rat/nPML/zmin/zmax, all read identically
    ! everywhere), so no MPI communication is needed to find the nearest
    ! node. Only z needs this: x and y already snap to the nearest interior
    ! node inside setSurfaceStation itself.
    real (kind = dp) :: xGridFull(10000), yGridFull(10000), zGridFull(10000)
    integer (kind = 4) :: xGridFullSize, yGridFullSize, zGridFullSize
    real (kind = dp) :: x4ndsSnapZ(max(1,totalNumOfOffSt))
    logical :: x4ndsZValid(max(1,totalNumOfOffSt))
    integer (kind = 4) :: iSt2, kNearest
    real (kind = dp) :: distBest, distCur

    call calcXyzMPIId(mex, mey, mez)
    call getLocalOneDimCoorArrAndSize(nxt, nxuni, edgex1, mex, nx, xline, modelBoundCoor, 1, xGridFull, xGridFullSize)
    call getLocalOneDimCoorArrAndSize(nyt, nyuni, edgey1, mey, ny, yline, modelBoundCoor, 2, yGridFull, yGridFullSize)
    call getLocalOneDimCoorArrAndSize(nzt, nzuni, edgezn, mez, nz, zline, modelBoundCoor, 3, zGridFull, zGridFullSize)
    ! Finding 1 (row 94 audit): "clamp to the nearest node INSIDE the mesh"
    ! does not mean "nearest node of the raw array", because the raw array
    ! includes the PML -- the absorbing-boundary layer setNumDof (below)
    ! itself marks non-physical by testing nodeCoor(3) < PMLb(5) (12-dof
    ! absorbing formulation instead of the normal ndof). PMLb(5) is set by
    ! the dimId==3 call just above to globalOneDimCoorArr(nPML+1), the
    ! shallowest PML node's neighbour -- the code's OWN boundary between
    ! "real material" and "absorbing layer". A station is a request for a
    ! PHYSICAL measurement point, so the valid clamp band is
    ! [PMLb(5), modelBoundCoor(3,2)] (modelBoundCoor(3,2) is the free
    ! surface, z=0) -- NOT [zmin, 0], which would let a station land inside
    ! the PML. A requested depth outside that band is left UNCLAMPED
    ! (x4ndsZValid=.false.): setSurfaceStation below never matches it on
    ! z, so it stays a genuine, loud DROP (report_dropped_offfault_st names
    ! the checked reason), not a snap into physically meaningless territory.
    !
    ! For an in-band request, nearest-node search is restricted to k with
    ! zGridFull(k) >= PMLb(5)-tol -- redundant given the band check above
    ! (the array is monotonic, so an in-band request's nearest node is
    ! always in-band too), but kept explicit rather than relied upon. A tie
    ! -- exactly half a grid spacing away from two nodes -- keeps the lower
    ! index (deeper/more negative z, whichever the ascending array reaches
    ! first with the strict `<`); ties are not expected at any station
    ! coordinate this project uses.
    x4ndsZValid = .false.
    do iSt2 = 1, totalNumOfOffSt
        if (x4nds(3,iSt2) < PMLb(5)-tol .or. x4nds(3,iSt2) > modelBoundCoor(3,2)+tol) cycle
        distBest = huge(distBest)
        kNearest = 1
        do k = 1, zGridFullSize
            if (zGridFull(k) < PMLb(5)-tol) cycle
            distCur = abs(zGridFull(k) - x4nds(3,iSt2))
            if (distCur < distBest) then
                distBest = distCur
                kNearest = k
            endif
        enddo
        x4ndsSnapZ(iSt2) = zGridFull(kNearest)
        x4ndsZValid(iSt2) = .true.
    enddo
    ! Finding 6 (row 94 audit): persist the band check so
    ! report_dropped_offfault_st (library_output.f90) can distinguish a
    ! CHECKED depth-out-of-band drop from an unconfirmed one. Identical on
    ! every rank (x4ndsZValid depends only on the global z grid and the
    ! request), so a plain overwrite each call is fine -- no reduction needed.
    x4ndsZValidPersist = x4ndsZValid
    xmin = modelBoundCoor(1,1)
    xmax = modelBoundCoor(1,2)
    ymin = modelBoundCoor(2,1)
    ymax = modelBoundCoor(2,2)
    zmin = modelBoundCoor(3,1)
    zmax = modelBoundCoor(3,2)
    ! Move adjacent plane1 and plane2 to create hexahedral meshes
    allocate(plane1(ny+ntotft,nz),plane2(ny+ntotft,nz),fltrc(2,nxuni,nzuni,ntotft))
    plane1 = 0
    plane2 = 0
    
    numcount    = 0
    numcount(1) = nx
    numcount(2) = ny
    numcount(3) = nz
    
    ! Initialize scalars

    msnode   = gridNodeCount(nx, ny, nz)
    numOfOnFaultStCount    = 0
    numOfOffFaultStCount    = 0
    OffFaultStNodeIdIndex   = 0
    ! Initialize arrays
    nftnd0   = 0
    
    ixfi = 0
    izfi = 0
    
    allocate(n4yn(totalNumOfOffSt))
    n4yn = 0
    
    ! Loop over x,y,z grids to create nodes and elements
    ! Create node
    do ix = 1, nx
        do iz = 1, nz
            do iy = 1, ny
                call initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
                call createNode(nodeCoor, xline(ix), yline(iy), zline(iz), nodeCount, nodeXyzIndex)
                if (insertFaultType > 0) then 
                    call insertFaultInterface(nodeCoor, ycoort, pfx, pfz)
                    meshCoor(2,nodeCount) = ycoort
                endif 
                
                call setNumDof(nodeCoor, numOfDofPerNodeTmp)

                eqNumStartIndexLoc(nodeCount) = eqNumIndexArrLocTag
                numOfDofPerNodeArr(nodeCount) = numOfDofPerNodeTmp
                
                call setEquationNumber(nodeXyzIndex, nodeCoor, eqNumIndexArrLocTag, equationNumCount, numOfDofPerNodeTmp)
                call setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nodeCount, x4ndsSnapZ, x4ndsZValid)
                call createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nodeCount, msnode, nftnd0, equationNumCount, eqNumIndexArrLocTag, &
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
                
                ! Create element
                if(ix>=2 .and. iy>=2 .and. iz>=2) then
                    call createElement(elemCount, stressDofCount, iy, iz, elementCenterCoor)
                    call checkPMLAlignment(elementCenterCoor)
                    call setElementMaterial(elemCount, elementCenterCoor)
                    
                    ! Row 153 checkpoint 1: the old C_degen>3.0d0 block only
                    ! ever tested/collapsed against fault 1 (hardcoded "1"
                    ! below). Generalized to loop over every fault and use
                    ! THAT fault's own faultDegenStyle/box, exiting on the
                    ! first fault that actually collapses the element (an
                    ! element can only belong to one fault's wedge pair). No
                    ! existing case combines ntotft>1 with degeneration, so
                    ! this loop executes exactly once, for fault 1, on every
                    ! committed case -- bit-identical to the old hardcoded call.
                    do ift = 1, ntotft
                        if (faultDegenStyle(ift) > 0) then
                            call wedge(elementCenterCoor(1), elementCenterCoor(2), elementCenterCoor(3), elemCount, stressDofCount, iy, iz, nftnd0(ift), ift)
                            isOnFt=0
                            call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(1,elemCount)), ift, isOnFt)
                            if (isOnFt==1 .and. elemTypeArr(elemCount)==1) elemTypeArr(elemCount) = 13

                            isOnFt=0
                            call checkIsOnFault(meshCoor(1:3,nodeElemIdRelation(2,elemCount)), ift, isOnFt)
                            if (isOnFt==1 .and. elemTypeArr(elemCount)==1) elemTypeArr(elemCount) = 13

                            ! Row 153 checkpoint 2a audit fix: exit only on
                            ! ACTUAL degeneration (type 11/12), not on a mere
                            ! type-13 "brick adjacent to this fault's plane"
                            ! match -- aligns this loop's stop condition with
                            ! countMeshEntities.f90's (which only ever exits
                            ! when wedge4num incremented elementCount, i.e.
                            ! actual degeneration). No registered case has
                            ! more than one fault with faultDegenStyle>0, so
                            ! this loop still executes exactly one matching
                            ! iteration everywhere it used to: bit-identical.
                            if (elemTypeArr(elemCount)==11 .or. elemTypeArr(elemCount)==12) exit
                        endif
                    enddo
                    
                    call replaceSlaveWithMasterNode(nodeCoor, elemCount, nftnd0) 
                    if (C_elastic == 0) call setPlasticStress(-0.5d0*(zline(iz)+zline(iz-1)) + 7.3215d0, elemCount)          
                 endif!if element
            enddo!iy
        enddo!iz
        plane1 = plane2
    enddo!ix
    
    sizeOfStressDofIndexArr=stressDofCount
    
    call meshGenError(nx, ny, nz, nodeCount, msnode, elemCount, equationNumCount, eqNumIndexArrLocTag, nftnd0)
    
    ! Row 17 (multi-fault): fltMPI reset ONCE here, before the per-fault
    ! MPI4arn loop below -- not inside MPI4arn itself any more (see that
    ! subroutine's header comment). No-op at ntotft==1.
    fltMPI = .false.

    ! compute on-fault area associated with each fault node pair and distance from source
    do ift=1,ntotft
        if(nftnd0(ift)>0) then
        !...element areas and distribute evenly to its four nodes
        do i=2,ifd(ift)
            do j=2,ifs(ift)
            !...4 nodes of quadrilateral
                n1 = fltrc(1,j,i,ift) !use nodal number in nodeCount
                n2 = fltrc(1,j-1,i,ift)
                n3 = fltrc(1,j-1,i-1,ift)
                n4 = fltrc(1,j,i-1,ift)
                m1 = fltrc(2,j,i,ift) !use nodal number in nftnd0
                m2 = fltrc(2,j-1,i,ift)
                m3 = fltrc(2,j-1,i-1,ift)
                m4 = fltrc(2,j,i-1,ift)
                !...calculate area of the quadrilateral
                !......if fault is not in coor axes plane
                aa1=sqrt((meshCoor(1,n2)-meshCoor(1,n1))**2 + (meshCoor(2,n2)-meshCoor(2,n1))**2 &
                + (meshCoor(3,n2)-meshCoor(3,n1))**2)
                bb1=sqrt((meshCoor(1,n3)-meshCoor(1,n2))**2 + (meshCoor(2,n3)-meshCoor(2,n2))**2 &
                + (meshCoor(3,n3)-meshCoor(3,n2))**2)
                cc1=sqrt((meshCoor(1,n4)-meshCoor(1,n3))**2 + (meshCoor(2,n4)-meshCoor(2,n3))**2 &
                + (meshCoor(3,n4)-meshCoor(3,n3))**2)
                dd1=sqrt((meshCoor(1,n1)-meshCoor(1,n4))**2 + (meshCoor(2,n1)-meshCoor(2,n4))**2 &
                + (meshCoor(3,n1)-meshCoor(3,n4))**2)
                p1=sqrt((meshCoor(1,n4)-meshCoor(1,n2))**2 + (meshCoor(2,n4)-meshCoor(2,n2))**2 &
                + (meshCoor(3,n4)-meshCoor(3,n2))**2)
                q1=sqrt((meshCoor(1,n3)-meshCoor(1,n1))**2 + (meshCoor(2,n3)-meshCoor(2,n1))**2 &
                + (meshCoor(3,n3)-meshCoor(3,n1))**2)
                area = 0.25d0 * sqrt(4*p1*p1*q1*q1 - &
               (bb1*bb1 + dd1*dd1 - aa1*aa1 -cc1*cc1)**2) 
                !...distribute above area to 4 nodes evenly
                area = 0.25d0 * area
                arn(m1,ift) = arn(m1,ift) + area
                arn(m2,ift) = arn(m2,ift) + area
                arn(m3,ift) = arn(m3,ift) + area
                arn(m4,ift) = arn(m4,ift) + area
                enddo
            enddo
        endif

        call MPI4arn(nx, ny, nz, mex, mey, mez, nftnd0(ift), ift)
    enddo  
end subroutine meshgen
!==================================================================================================
!**************************************************************************************************
!==================================================================================================

subroutine setElementMaterial(elemCount, elementCenterCoor)
! Subroutine velocityStructure will asign Vp, Vs and rho 
!   based on input from bMaterial.txt, which is created by 
!   case input file user_defined_param.py.

    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8) :: elemCount
    integer (kind = 4) :: i, sideOfFault, idx(3)
    real (kind = dp) :: elementCenterCoor(3), vptmp, vstmp, rhotmp

    if (nmat == 1 .and. n2mat == 3) then
    ! homogenous material
        mat(elemCount,1)  = material(1,1)
        mat(elemCount,2)  = material(1,2)
        mat(elemCount,3)  = material(1,3)
    elseif (nmat > 1 .and. n2mat == 4) then 
        ! 1D velocity structure
        ! material = material(i,j), i=1,nmat, and j=1,4
        ! for j
        !   1: bottom depth of a layer (should be positive in m).
        !   2: Vp, m/s
        !   3: Vs, m/s
        !   4: rho, kg/m3
        if (abs(elementCenterCoor(3)) < material(1,1)) then
            mat(elemCount,1)  = material(1,2)
            mat(elemCount,2)  = material(1,3)
            mat(elemCount,3)  = material(1,4)
        else
            do i = 2, nmat
                if (abs(elementCenterCoor(3)) < material(i,1) &
                    .and. abs(elementCenterCoor(3)) >= material(i-1,1)) then
                    
                    mat(elemCount,1)  = material(i,2)
                    mat(elemCount,2)  = material(i,3)
                    mat(elemCount,3)  = material(i,4)
                endif
            enddo
        endif
    elseif (nmat > 1 .and. n2mat == 5) then
        ! Two-sided 1D velocity structure (SCEC TPV35: a different 1D
        ! profile on each side of the fault). Columns:
        !   1: bottom depth of a layer (positive, m)   2: Vp   3: Vs   4: rho
        !   5: side, -1 for element centres at y < the fault y-plane,
        !      +1 for y > it.
        ! The fault y-plane is the ONE vertical plane every fault shares
        ! (readmaterial's checkTwoSidedMaterialTable refuses anything else),
        ! so minval over the fault axis is that plane, not a fault-1 read.
        ! Within a side the rule is the n2mat==4 one above (ascending
        ! bottoms, first match). The table is validated once in
        ! readmaterial (ERR_CFG_MATERIAL_TABLE_INVALID): coplanar vertical
        ! faults, side column -1/+1, per-side bottoms strictly ascending.
        if (elementCenterCoor(2) < minval(fltxyz(1,2,:))) then
            sideOfFault = -1
        else
            sideOfFault = 1
        endif
        do i = 1, nmat
            if (nint(material(i,5)) /= sideOfFault) cycle
            if (abs(elementCenterCoor(3)) < material(i,1)) then
                mat(elemCount,1)  = material(i,2)
                mat(elemCount,2)  = material(i,3)
                mat(elemCount,3)  = material(i,4)
                exit
            endif
        enddo
    elseif (nmat > 1 .and. n2mat == 6) then
        ! 3D structured material grid (SCEC TPV34: CVM-H sampled at the
        ! uniform element-centre spacing). bMaterial.txt rows are
        ! [x y z vp vs rho]; readmaterial's buildMaterialGrid3D validated
        ! them into matGrid3D with matGridOrigin/Spacing/Count. Each element
        ! reads the NEAREST grid cell to its centre, clamped to the grid
        ! (piecewise constant, never interpolated): a uniform-belt element
        ! centre lies ON the grid and reads its own sample exactly; a
        ! stretched or PML element off the grid reads the nearest sample.
        ! floor(off + 0.5) is the tie rule the Python port uses too (nint
        ! and numpy.rint differ at exact .5 offsets).
        do i = 1, 3
            idx(i) = floor((elementCenterCoor(i) - matGridOrigin(i))/matGridSpacing(i) + 0.5d0) + 1
            idx(i) = max(1, min(matGridCount(i), idx(i)))
        enddo
        mat(elemCount,1) = matGrid3D(1, idx(1), idx(2), idx(3))
        mat(elemCount,2) = matGrid3D(2, idx(1), idx(2), idx(3))
        mat(elemCount,3) = matGrid3D(3, idx(1), idx(2), idx(3))
    endif

    ! calculate lambda and mu from Vp, Vs and rho.
    ! mu = Vs**2*rho
    mat(elemCount,5)  = mat(elemCount,2)**2*mat(elemCount,3)
    ! lambda = Vp**2*rho-2*mu
    mat(elemCount,4)  = mat(elemCount,1)**2*mat(elemCount,3)-2.0d0*mat(elemCount,5)

end subroutine setElementMaterial

subroutine MPI4arn(nx, ny, nz, mex, mey, mez, totalNumFaultNode, iFault)
! Add up arn from neighbor MPI blocks.
    use globalvar
    use errorCodes
    implicit none 
    include 'mpif.h'
    
    integer (kind = 4) :: bndl,bndr,bndf,bndb,bndd,bndu, nx, ny, nz, mex, mey, mez
    integer (kind = 4) :: totalNumFaultNode, iFault, i

    ! Row 17 (multi-fault): fltMPI is reset ONCE by the caller (meshgen.f90,
    ! before its `do ift=1,ntotft` loop that calls this subroutine once per
    ! fault) -- NOT here any more. Resetting it on every call, as before,
    ! wiped out a direction flag an EARLIER fault had already set whenever a
    ! LATER fault did not touch that same boundary direction. fltMPI(k) is
    ! still only ever set to .true. (never back to .false.) below, in
    ! syncArnBoundary, so it correctly accumulates an OR across every fault's
    ! calls. No-op at ntotft==1 (one call either way).
    !
    ! fltgm/fltl/fltr/fltf/fltb/fltd/fltu/fltnum are all (.., ntotft) now,
    ! pre-allocated once in allocInit (eqdyna3d.f90) -- this subroutine fills
    ! column iFault in place instead of deallocating/reallocating a
    ! single-fault-sized array every call (which is what used to make every
    ! fault but the last invisible to assembleGlobalMass.f90's
    ! MPI4NodalQuant/addFaultBoundaryTerm downstream). Reduces to the old
    ! single-fault fill at ntotft==1 (column 1 only, same values).
    fltnum(:,iFault) = 0
    do i = 1, totalNumFaultNode
        if(mod(fltgm(i,iFault),10)==1 ) then
            fltnum(1,iFault) = fltnum(1,iFault) + 1
            fltl(fltnum(1,iFault),iFault) = i
        endif
        if(mod(fltgm(i,iFault),10)==2 ) then
            fltnum(2,iFault) = fltnum(2,iFault) + 1
            fltr(fltnum(2,iFault),iFault) = i
        endif
        if(mod(fltgm(i,iFault),100)-mod(fltgm(i,iFault),10)==10 ) then
            fltnum(3,iFault) = fltnum(3,iFault) + 1
            fltf(fltnum(3,iFault),iFault) = i
        endif
        if(mod(fltgm(i,iFault),100)-mod(fltgm(i,iFault),10)==20 ) then
            fltnum(4,iFault) = fltnum(4,iFault) + 1
            fltb(fltnum(4,iFault),iFault) = i
        endif
        if(fltgm(i,iFault)-mod(fltgm(i,iFault),100)==100 ) then
            fltnum(5,iFault) = fltnum(5,iFault) + 1
            fltd(fltnum(5,iFault),iFault) = i
        endif
        if(fltgm(i,iFault)-mod(fltgm(i,iFault),100)==200 ) then
            fltnum(6,iFault) = fltnum(6,iFault) + 1
            fltu(fltnum(6,iFault),iFault) = i
        endif
    enddo
    ! FIX (pathway_forward.md, "fault-plane-on-MPI-boundary halving"): arn is
    ! accumulated per fault node from the LOCAL ifs x ifd quad grid (lines
    ! ~115-152 above), independent of y -- so whenever this rank's boundary
    ! in a direction the fault has ZERO nominal extent in (per fltxyz; always
    ! y for EQdyna's x-z-plane faults) carries fault nodes, this rank has
    ! already built the FULL local fault-node arn there (both ranks rebuild
    ! the whole ifs x ifd grid identically), so summing across that boundary
    ! DUPLICATES rather than DIVIDES the area. Directions the fault has real
    ! extent in (x, z -- or any future non-degenerate direction) are genuine
    ! divisions: each rank only owns a slice of the fault there, and summing
    ! is correct (verified: scratch/mira_faultmpi_guard, x-split and z-split
    ! reproduce the serial arn exactly; y-split alone doubled it before this
    ! fix, matched after).
    !
    ! A prior attempt skipped the ENTIRE `call syncArnBoundary` for such
    ! directions and was reverted: fltMPI(k) is set at the top of
    ! syncArnBoundary and is ALSO the gate `addFaultBoundaryTerm`
    ! (assembleGlobalMass.f90) uses to decide whether fnms/nodalMassArr need
    ! a cross-rank correction on that same boundary -- skipping the whole
    ! call silently skipped that correction too and broke a rough-fault
    ! 2,2,1 run (zero mass / NaN velocity). This fix instead still calls
    ! syncArnBoundary (fltMPI(k) still gets set, the exchange still happens,
    ! fnms/nodalMassArr's correction is untouched) and only skips arn's OWN
    ! add-back when the direction is degenerate (see syncArnBoundary below).
    ! Direct evidence fnms is NOT affected by this bug (unlike arn): dumped
    ! fnms at the same physical fault node is bit-identical between serial
    ! and a fault-splitting decomposition (scratch/mira_faultmpi_guard) --
    ! fnms's volume-element accumulation is genuinely owned exclusively by
    ! one side's rank per split-node copy, so its cross-rank correction on
    ! this boundary is a no-op there in the cases checked; left unchanged.
    if (npx > 1) then
        bndl=1  !left boundary
        bndr=nx ! right boundary
        if (mex == masterProcsId) then
            bndl=0
        elseif (mex == npx-1) then
            bndr=0
        endif

        if (bndl/=0) then
            if(fltnum(1,iFault)>0 ) call syncArnBoundary(fltl(:,iFault), fltnum(1,iFault), 1, me-npy*npz, 1000, iFault, 1)
        endif !if bhdl/=0

        if (bndr/=0) then
            if(fltnum(2,iFault)>0 ) call syncArnBoundary(fltr(:,iFault), fltnum(2,iFault), 2, me+npy*npz, 1000, iFault, 1)
        endif !bndr/=0
    endif !npx>1
!*****************************************************************************************
    if (npy > 1) then   !MPI:no fault message passing along y due to the x-z vertical fault.
        bndf=1  !front boundary
        bndb=ny ! back boundary
        if (mey == masterProcsId) then
            bndf=0
        elseif (mey == npy-1) then
            bndb=0
        endif

        if (bndf/=0) then
            if(fltnum(3,iFault)>0) call syncArnBoundary(fltf(:,iFault), fltnum(3,iFault), 3, me-npz, 2000, iFault, 2)
        endif !bhdf/=0

        if (bndb/=0) then
            if(fltnum(4,iFault)>0) call syncArnBoundary(fltb(:,iFault), fltnum(4,iFault), 4, me+npz, 2000, iFault, 2)
        endif !bndb/=0
     endif !npy>1
!*****************************************************************************************
    if (npz > 1) then
        bndd=1  !lower(down) boundary
        bndu=nz ! upper boundary
        if (mez == masterProcsId) then
            bndd=0
        elseif (mez == npz-1) then
            bndu=0
        endif
        if (bndd/=0) then
            if(fltnum(5,iFault)>0) call syncArnBoundary(fltd(:,iFault), fltnum(5,iFault), 5, me-1, 3000, iFault, 3)
        endif !bhdd/=0

        if (bndu/=0) then
            if(fltnum(6,iFault)>0) call syncArnBoundary(fltu(:,iFault), fltnum(6,iFault), 6, me+1, 3000, iFault, 3)
        endif !bndu/=0
    endif !npz>1
contains
    subroutine syncArnBoundary(idxArr, n, k, neighbor, tagBase, ift, dimId)
    ! Send this rank's arn values for the given boundary's fault nodes to
    ! neighbor, receive neighbor's values for the same nodes, and accumulate
    ! -- UNLESS the fault has zero NOMINAL GRID extent in the dimension
    ! (dimId: 1=x, 2=y, 3=z) this boundary lies along, in which case this
    ! rank's local arn is already the fault's full local contribution and
    ! the neighbor's value is a duplicate, not a partial sum.
    !
    ! "NOMINAL GRID extent" is fltxyz, and it is NOT the physical dip. The
    ! two ways checkIsOnFault selects fault nodes give opposite answers here:
    !
    !   insertFaultType > 0, C_degen == 0  (e.g. test.tpv10, dip 60)
    !       nodes chosen by nodeCoor(2) == 0.0d0, exact equality on the
    !       UNBLENDED grid y. The fault is ONE y-index plane however steeply
    !       it dips; insertFaultInterface displaces only the physical y into
    !       ycoort/meshCoor. fltxyz y-extent is 0.0/0.0, so a y boundary
    !       COINCIDES with the fault and both ranks hold the COMPLETE
    !       tributary area  ->  DUPLICATE, skip the add-back.
    !
    !   C_degen > 3  (wedge degeneration, e.g. test.tpv36/37, dip 15)
    !       nodes chosen by |z + y*tan(C_degen)| < dx/100. The fault SPANS a
    !       range of grid y, fltxyz y-extent is faultWidth*cos(dip), so a y
    !       boundary CROSSES it and each rank holds a partial area
    !       ->  DIVIDE, add back.
    !
    ! A steeply dipping inserted fault therefore behaves exactly like a
    ! vertical one for this decision. Anyone tempted to branch on the dip
    ! angle here would break test.tpv10.
    !
    ! Evidence: DUPLICATE -- test_fault_mpi_boundary_arn.py (serial == xsplit
    ! == zsplit == ysplit, exactly). DIVIDE -- test_dipping_fault_y_split.py
    ! (tpv36 at (2,2,1) vs (2,1,2), tractions equal to 1.0e-08, ratio
    ! 1.000000 at 3416 nodes; audited 2026-09-16). fltMPI(k) is still set and the
    ! exchange still happens either way: addFaultBoundaryTerm
    ! (assembleGlobalMass.f90) depends on fltMPI(k)/the send-recv having run,
    ! independent of what this subroutine does with arn.
        integer (kind = 4), intent(in) :: idxArr(:), n, k, neighbor, tagBase, ift, dimId
        integer (kind = 4) :: j, jMPIstatus(MPI_STATUS_SIZE), jMPIerr
        real (kind = dp), allocatable, dimension(:) :: sendBuf, recvBuf

        fltMPI(k)=.true.
        allocate(sendBuf(n),recvBuf(n))
        sendBuf = 0.
        recvBuf = 0.
        do j = 1, n
            sendBuf(j)=arn(idxArr(j),ift)
        enddo
        call requireValidNeighbor(neighbor, 'syncArnBoundary (meshgen/MPI4arn)', &
            dimId, k, 'fault-node arn boundary exchange; entered only when fltnum(k)>0 on THIS rank, which is a local condition')
        call mpi_sendrecv(sendBuf, n, MPI_DOUBLE_PRECISION, neighbor, tagBase+me, &
            recvBuf, n, MPI_DOUBLE_PRECISION, neighbor, tagBase+neighbor, &
            MPI_COMM_WORLD, jMPIstatus, jMPIerr)
        if (fltxyz(2,dimId,ift) /= fltxyz(1,dimId,ift)) then
            do j = 1, n
                arn(idxArr(j),ift) = arn(idxArr(j),ift) + recvBuf(j)
            enddo
        endif
        deallocate(sendBuf,recvBuf)
    end subroutine syncArnBoundary
end subroutine MPI4arn

subroutine meshGenError(nx, ny, nz, nodeCount, msnode, elemCount, equationNumCount, eqNumIndexArrLocTag, nftnd0)
! Check consistency between mesh4 and meshgen
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: nx, ny, nz, nftnd0(ntotft)
    integer (kind = 8) :: nodeCount, msnode, elemCount, equationNumCount, eqNumIndexArrLocTag
    integer (kind = 4) :: i
    if (sizeOfStressDofIndexArr>=(5*sizeOfEqNumIndexArr)) then
        write(*,*) '5*sizeOfEqNumIndexArr',sizeOfEqNumIndexArr,'is not enough for sizeOfStressDofIndexArr',sizeOfStressDofIndexArr
        call abortRun(ERR_MESH_STRESS_ARR_SMALL, &
            'stressArr is too small for this mesh; 5*sizeOfEqNumIndexArr does not cover sizeOfStressDofIndexArr.')
    endif
    if(nodeCount/=gridNodeCount(nx, ny, nz).or.msnode/=totalNumOfNodes.or.elemCount/=totalNumOfElements.or.equationNumCount/=totalNumOfEquations) then
        write(*,*) 'Inconsistency in node/element/equation/between meshgen and countMeshEntities: stop!',me
        write(*,*) 'nodeCount&totalNumOfNodes=',nodeCount,totalNumOfNodes
        write(*,*) 'elemCount,totalNumOfElements=',elemCount,totalNumOfElements
        write(*,*) 'equationNumCount,totalNumOfEquations=',equationNumCount,totalNumOfEquations
        write(*,*) 'nodeCount,nx,ny,nz',nodeCount,nx,ny,nz
        write(*,*) 'msnode,totalNumOfNodes',msnode,totalNumOfNodes
        call abortRun(ERR_MESH_COUNT_MISMATCH, &
            'meshgen node/element/equation tallies disagree with countMeshEntities (counts printed above).')
    endif
    if(eqNumIndexArrLocTag/=sizeOfEqNumIndexArr) then
        write(*,*) 'Inconsistency in eqNumIndexArrLocTag and sizeOfEqNumIndexArr: stop!',me
        write(*,*) eqNumIndexArrLocTag,sizeOfEqNumIndexArr
        call abortRun(ERR_MESH_EQNUM_MISMATCH, &
            'eqNumIndexArrLocTag /= sizeOfEqNumIndexArr (values printed above).')
    endif
    do i=1,ntotft
        if(nftnd0(i)/=nftnd(i)) then
            write(*,*) 'Inconsistency in fault between meshgen and countMeshEntities:',me,i
            call abortRun(ERR_MESH_FAULT_MISMATCH, &
                'nftnd0 from meshgen disagrees with nftnd from countMeshEntities (rank and fault printed above).')
        endif
    enddo
end subroutine meshGenError

subroutine calcXyzMPIId(mex, mey, mez)
    use globalvar 
    use errorCodes
    implicit none
    integer (kind = 4) :: mex, mey, mez
     
    mex=int(me/(npy*npz))
    mey=int((me-mex*npy*npz)/npz)
    mez=int(me-mex*npy*npz-mey*npz)
end subroutine calcXyzMPIId

subroutine getLocalOneDimCoorArrAndSize(globalOneDimCoorArrSize, numOfNodesWithUniformGridsize, &
    frontEdgeNodeId, MPIXyzId, localOneDimCoorArrSize, localOneDimCoorArr, modelBoundCoor, dimId, &
    globalOneDimCoorArrOut, globalOneDimCoorArrOutSize)
    ! globalOneDimCoorArrOut/globalOneDimCoorArrOutSize (pathway item 94,
    ! owner ruling 2026-09-24): the FULL 1D grid this subroutine already
    ! builds before slicing it down to this rank's local piece, exposed so a
    ! caller can find the nearest node to an arbitrary requested coordinate
    ! without any MPI communication -- every rank builds this same array
    ! from the same inputs (dz/fltxyz/rat/nPML/zmin/zmax etc.), so it is
    ! identical bit-for-bit everywhere. Fixed-size (10000), matching
    ! localOneDimCoorArr's own existing convention, so the implicit
    ! (no-interface) call sites need only pass a same-shaped actual argument.
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: globalOneDimCoorArrSize, numOfNodesWithUniformGridsize, dimId
    integer (kind = 4) :: frontEdgeNodeId, localOneDimCoorArrSize, MPIXyzId
    integer (kind = 4) :: numOfNodesPerMPI, residualNumOfNodes
    integer (kind = 4) :: i, numOfMPIXyz, ift
    integer (kind = 4) :: globalOneDimCoorArrOutSize
    real (kind = dp) :: gridSize, frontEdgeCoor, backEdgeCoor, &
            minCoor, maxCoor, coorTmp, gridSizeTmp, localOneDimCoorArr(10000), &
            modelBoundCoor(3,2), globalOneDimCoorArrOut(10000)
    real (kind = dp) :: fltMin1, fltMax1, commensOffset
    real (kind = dp), allocatable :: globalOneDimCoorArr(:)
    character(len=300) :: reasonMsg

    ! Row 17 rebased (restore per-fault mesh extent, invent nothing): the
    ! uniform x/z belt (dimId==1/3) and the uniform-y belt (dimId==2) used to
    ! be anchored on FAULT 1's box alone (fltxyz(.,.,1)) -- every other
    ! fault's box was never consulted here at all, so a case with a fault
    ! outside fault 1's x/z extent, or outside a hand-widened fixed y-margin,
    ! either refused to run (checkInputConsistency's old guards) or meshed
    ! wrong. fltMin1/fltMax1 are now the UNION (minval/maxval) over every
    ! fault's own box on this axis, mirroring eqquasi's b6010e2 (x/z) and
    ! 2769c73 (y). ntotft==1 collapses fltMin1==fltMax1==fault 1's own single
    ! value on every axis, so every formula below is bit-identical to the old
    ! fault-1-only one in that case.
    if (dimId == 1) then
        fltMin1 = minval(fltxyz(1,1,1:ntotft))
        fltMax1 = maxval(fltxyz(2,1,1:ntotft))
        numOfNodesWithUniformGridsize = nint((fltMax1 - fltMin1)/dx) + 1
        gridSize      = dx
        frontEdgeCoor = fltMin1
        backEdgeCoor  = fltMax1
        minCoor       = xmin
        maxCoor       = xmax
        numOfMPIXyz   = npx
    elseif (dimId == 2) then
        if (C_degen == 0.0d0) then
            fltMin1 = minval(fltxyz(1,2,1:ntotft))
            fltMax1 = maxval(fltxyz(2,2,1:ntotft))
        else
            ! C_degen>3 (wedge-degeneration dipping fault, e.g. tpv36/tpv37):
            ! a DIFFERENT, pre-existing, single-fault-only mechanism
            ! (orthogonal to this work -- see checkInputConsistency.f90's own
            ! "C_degen>3 ... legitimate fymin /= fymax" comment) where
            ! fltxyz(1,2,1)/fltxyz(2,2,1) describe a dip-projection y-RANGE,
            ! not a fault y-PLANE the uniform belt should be anchored to --
            ! the belt was never built from fltxyz's y at all before this
            ! mission (frontEdgeCoor/backEdgeCoor were -/+ dis4uniF/B*dy
            ! regardless), and generalizing it to the union here broke
            ! test.tpv36 (measured: test_rank_local_mesh.py's 8-rank fault
            ! ownership count changed). Left exactly as before, bit-identical
            ! -- this mission's multi-fault scope is ntotft==1-only for
            ! C_degen>3 regardless, so no multi-fault y-belt case is lost.
            fltMin1 = 0.0d0
            fltMax1 = 0.0d0
        endif
        numOfNodesWithUniformGridsize = nint((fltMax1 - fltMin1)/dy) + dis4uniF + dis4uniB + 1
        gridSize      = dy
        frontEdgeCoor = fltMin1 - dis4uniF*dy
        backEdgeCoor  = fltMax1 + dis4uniB*dy
        minCoor       = ymin
        maxCoor       = ymax
        numOfMPIXyz   = npy
    elseif (dimId == 3) then
        fltMin1 = minval(fltxyz(1,3,1:ntotft))
        fltMax1 = maxval(fltxyz(2,3,1:ntotft))
        numOfNodesWithUniformGridsize = nint((fltMax1 - fltMin1)/dz) + 1
        gridSize      = dz
        frontEdgeCoor = fltMin1
        backEdgeCoor  = fltMax1
        minCoor       = zmin
        maxCoor       = zmax
        numOfMPIXyz   = npz
    endif

    ! Hard refuse (not a silent clamp or a widened margin) when a fault's
    ! bound on this axis is not an integer number of gridSize steps from
    ! fltMin1, the belt origin just computed -- eqquasi's eeac6f9/a761f33
    ! finding, ported here: a fault edge that falls between node lines meshes
    ! with fewer fault nodes than declared, or none, SILENTLY. Checked once
    ! per fault per bound (lower and upper -- a planar y-fault has both
    ! bounds equal, so this subsumes the old "fymin not a multiple of dy"
    ! check with the now-correct, union-based origin instead of a hardcoded
    ! 0). ntotft==1 makes fltMin1 this very fault's own bound, so the offset
    ! is always exactly 0 and this can never fire for a single-fault case
    ! (bit-identical no-op). Same abortRun-with-the-actual-numbers-named style
    ! as the overflow refusal below.
    !
    ! Tolerance: `tol` (globalvar.f90, 1.0d-5 m), NOT gridSize/100 (victor-reyes
    ! audit, PR #76 MAJOR 1). eqquasi's own eeac6f9/a761f33 use dx/100 for this
    ! check, but EQdyna's own downstream consumer, checkIsOnFault, only ever
    ! matches a node to a fault plane within `tol`=1e-5 m regardless of
    ! gridSize -- so a gridSize/100 pass band (2 m at dy=200) let a fault
    ! bound land up to 2 m off a node line, this check pass silently, and
    ! checkIsOnFault then fail to match ANY node (1e-5 m << 2 m), i.e. the
    ! exact silent zero-fault-node failure this guard exists to prevent. Using
    ! `tol` ties the two checks to the same constant by construction. The
    ! prior (pre-this-PR) single-fault check in checkInputConsistency.f90 used
    ! 1.0d-6 relative to dy (abs(fymin/dy - nint(fymin/dy)) > 1.0d-6, i.e.
    ! 1e-6*dy absolute -- 2e-4 m at dy=200), itself looser than `tol`=1e-5 m
    ! at any dy>10; `tol` is the tighter, defensible, already-shared constant.
    ! Measured on TPV22 (dy=200, fault y-offset 1600 m) and TPV23 (dy=250,
    ! fault y-offset 1000 m): both offsets are exact integer multiples of dy
    ! in double precision (1600.0d0/200.0d0 and 1000.0d0/250.0d0 both exact
    ! integers, no representable remainder), so commensOffset==0.0d0 exactly
    ! for both cases -- tightening this check to `tol` does not newly refuse
    ! either mesh.
    !
    ! dimId==2 skips entirely when C_degen/=0: the belt for that axis is the
    ! untouched fixed-margin one (see the dimId==2 branch above), not built
    ! from fltxyz's y at all, so testing fltxyz(.,2,ift)'s dip-projection
    ! y-range for commensurability with dy would be checking a quantity this
    ! guard has nothing to do with (and tpv36/tpv37's own y-range is not
    ! generally dy-commensurate from 0 -- it never needed to be).
    do ift = 1, ntotft
        if (dimId == 2 .and. C_degen /= 0.0d0) cycle
        do i = 1, 2
            commensOffset = fltxyz(i,dimId,ift) - fltMin1 - &
                dble(nint((fltxyz(i,dimId,ift) - fltMin1)/gridSize))*gridSize
            if (abs(commensOffset) > tol) then
                write(reasonMsg,'(a,i0,a,i0,a,i0,a,f0.3,a,f0.3,a,f0.3,a)') &
                    'getLocalOneDimCoorArrAndSize: dimId=', dimId, ', fault ', ift, &
                    ' bound ', i, ' = ', fltxyz(i,dimId,ift), &
                    ' is not an integer multiple of gridSize=', gridSize, &
                    ' from this axis'' belt origin=', fltMin1, &
                    '; it would fall between mesh node lines and mesh with too few ' // &
                    'fault nodes, or none, silently. Fix par.faultgeom or dx/dy/dz.'
                if (dimId == 2) then
                    call abortRun(ERR_GEOM_MULTIFAULT_Y_BAD, trim(reasonMsg))
                else
                    call abortRun(ERR_GEOM_MULTIFAULT_XZ_BAD, trim(reasonMsg))
                endif
            endif
        enddo
    enddo

    coorTmp = frontEdgeCoor
    gridSizeTmp = gridSize
    do i = 1, np
        gridSizeTmp = gridSizeTmp * rat
        coorTmp = coorTmp - gridSizeTmp
        if (coorTmp <= minCoor) exit
    enddo 
    frontEdgeNodeId = i + nPML
  
    coorTmp = backEdgeCoor
    gridSizeTmp = gridSize
    do i = 1, np
        gridSizeTmp = gridSizeTmp * rat
        coorTmp = coorTmp + gridSizeTmp
        if (coorTmp >= maxCoor) exit
    enddo
    if (dimId == 3) i = -nPML
    globalOneDimCoorArrSize = numOfNodesWithUniformGridsize + frontEdgeNodeId + i + nPML
    ! Finding 2 (row 94 audit): globalOneDimCoorArrOut/xline/yline/zline/
    ! localOneDimCoorArr are ALL fixed-size(10000) below and in the caller
    ! (meshgen); the copy at the end of this subroutine
    ! (globalOneDimCoorArrOut(1:globalOneDimCoorArrSize) = ...) is an
    ! unchecked write into that fixed buffer. A grid fine/large enough to
    ! exceed 10000 nodes used to overflow it silently -- hard-refuse instead
    ! (rule 2: no silent overflow), before the allocate below, so the message
    ! names the actual oversized count rather than segfaulting downstream.
    if (globalOneDimCoorArrSize > 10000) then
        write(reasonMsg,'(a,i0,a,i0,a)') &
            'getLocalOneDimCoorArrAndSize: dimId=', dimId, &
            ' built a full 1D grid of ', globalOneDimCoorArrSize, &
            ' nodes, exceeding the fixed-size 10000 buffer (xline/yline/zline/' // &
            'globalOneDimCoorArrOut/localOneDimCoorArr); reduce nx/ny/nz (dx/dy/dz) for ' // &
            'this dimension, or raise the fixed buffer size (10000) in meshgen.f90 and ' // &
            'this subroutine.'
        call abortRun(ERR_MESH_GRID_TOO_LARGE, trim(reasonMsg))
    endif
    allocate(globalOneDimCoorArr(globalOneDimCoorArrSize))

    numOfNodesPerMPI = int((globalOneDimCoorArrSize+numOfMPIXyz-1)/numOfMPIXyz)
    residualNumOfNodes = (globalOneDimCoorArrSize+numOfMPIXyz-1) - numOfNodesPerMPI * numOfMPIXyz

    if (MPIXyzId<(numOfMPIXyz-residualNumOfNodes)) then
        localOneDimCoorArrSize = numOfNodesPerMPI
    else
        localOneDimCoorArrSize = numOfNodesPerMPI + 1
    endif
    ! Row 132(1): a rank whose local axis holds fewer than 2 nodes owns no
    ! element on that axis (a 1-node slab sits on both of its own faces at
    ! once) and, separately, breaks setSurfaceStation's elseif chain (ix==1/
    ! ix==nx, iy==1/iy==ny are mutually exclusive branches -- when nx==1 (or
    ! ny==1) `ix==1` is tested first and always wins, so the `ix==nx` branch,
    ! which is ALWAYS an ownership candidate per row 127's rule, becomes
    ! unreachable and any station this rank was supposed to own as the lower
    ! side of a seam is silently dropped, on no rank). The python port
    ! refuses this same shape loudly already (`check_partition_1d`,
    ! src/python/eqdyna/meshgen.py) -- Fortran did not. Every rank computes
    ! its own local size here, so whichever rank(s) are too thin abort (and
    ! MPI_Abort tears down the whole communicator), which is equivalent to a
    ! global check.
    if (localOneDimCoorArrSize < 2) then
        write(reasonMsg,'(a,i0,a,i0,a,i0,a,i0,a)') &
            'getLocalOneDimCoorArrAndSize: dimId=', dimId, ', MPIXyzId=', MPIXyzId, &
            ' of ', numOfMPIXyz, ' holds a local 1D slice of only ', localOneDimCoorArrSize, &
            ' node(s); every rank needs >=2 nodes on each axis to own at least one element ' // &
            'and for setSurfaceStation''s ownership branches to be reachable -- reduce ' // &
            'npx/npy/npz for this axis, or use a finer grid (dx/dy/dz).'
        call abortRun(ERR_MPI_AXIS_TOO_THIN, trim(reasonMsg))
    endif
    ! Finding 2: the per-rank local slice is copied into localOneDimCoorArr,
    ! also fixed-size(10000), by the loops just below -- this used to only
    ! PRINT a warning and continue straight into the overflowing write.
    if (localOneDimCoorArrSize > 10000) then
        write(reasonMsg,'(a,i0,a,i0,a,i0,a,i0,a)') &
            'getLocalOneDimCoorArrAndSize: dimId=', dimId, ', MPIXyzId=', MPIXyzId, &
            ' has a local slice of ', localOneDimCoorArrSize, &
            ' nodes (', 10000, &
            '-node fixed buffer localOneDimCoorArr); reduce nx/ny/nz per rank or increase npx/npy/npz.'
        call abortRun(ERR_MESH_GRID_TOO_LARGE, trim(reasonMsg))
    endif

    globalOneDimCoorArr(frontEdgeNodeId+1) = frontEdgeCoor
    gridSizeTmp = gridSize
    do i = frontEdgeNodeId, 1, -1
        gridSizeTmp = gridSizeTmp * rat 
        globalOneDimCoorArr(i) = globalOneDimCoorArr(i+1) - gridSizeTmp
    enddo 
    do i = frontEdgeNodeId+2, frontEdgeNodeId + numOfNodesWithUniformGridsize
        globalOneDimCoorArr(i) = globalOneDimCoorArr(i-1) + gridSize
    enddo 
    if (dimId < 3) then 
        gridSizeTmp = gridSize
        do i = frontEdgeNodeId+numOfNodesWithUniformGridsize+1, globalOneDimCoorArrSize
            gridSizeTmp = gridSizeTmp * rat
            globalOneDimCoorArr(i) = globalOneDimCoorArr(i-1) + gridSizeTmp
        enddo 
    endif 
    if (dimId == 1) then 
        modelBoundCoor(1,1) = globalOneDimCoorArr(1)
        modelBoundCoor(1,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)
        PMLb(1) = globalOneDimCoorArr(globalOneDimCoorArrSize-nPML)
        PMLb(2) = globalOneDimCoorArr(nPML+1)
        PMLb(6) = globalOneDimCoorArr(globalOneDimCoorArrSize) - globalOneDimCoorArr(globalOneDimCoorArrSize-1)
    elseif (dimId == 2) then 
        modelBoundCoor(2,1) = globalOneDimCoorArr(1)
        modelBoundCoor(2,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)     
        PMLb(3) = globalOneDimCoorArr(globalOneDimCoorArrSize-nPML)
        PMLb(4) = globalOneDimCoorArr(nPML+1) 
        PMLb(7) = globalOneDimCoorArr(globalOneDimCoorArrSize) - globalOneDimCoorArr(globalOneDimCoorArrSize-1)
    elseif (dimId == 3) then
        modelBoundCoor(3,1) = globalOneDimCoorArr(1) 
        modelBoundCoor(3,2) = globalOneDimCoorArr(globalOneDimCoorArrSize)    
        PMLb(5) = globalOneDimCoorArr(nPML+1)
        PMLb(8) = globalOneDimCoorArr(2) - globalOneDimCoorArr(1)
    endif 

    if (MPIXyzId <= (numOfMPIXyz - residualNumOfNodes)) then
        do i = 1, localOneDimCoorArrSize
            localOneDimCoorArr(i) = globalOneDimCoorArr((numOfNodesPerMPI-1)*MPIXyzId+i)
        enddo
    else
        do i = 1, localOneDimCoorArrSize
            localOneDimCoorArr(i) = globalOneDimCoorArr((numOfNodesPerMPI-1)*MPIXyzId+i+(MPIXyzId-numOfMPIXyz+residualNumOfNodes))
        enddo
    endif

    globalOneDimCoorArrOutSize = globalOneDimCoorArrSize
    globalOneDimCoorArrOut(1:globalOneDimCoorArrSize) = globalOneDimCoorArr(1:globalOneDimCoorArrSize)
end subroutine getLocalOneDimCoorArrAndSize

subroutine setNumDof(nodeCoor, numOfDofPerNodeTmp)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: numOfDofPerNodeTmp
    real (kind = dp) :: nodeCoor(10)
    numOfDofPerNodeTmp = ndof !Default
    if (nodeCoor(1)>PMLb(1) .or. nodeCoor(1)<PMLb(2) .or. nodeCoor(2)>PMLb(3) &
        .or. nodeCoor(2)<PMLb(4) .or. nodeCoor(3)<PMLb(5)) then
        numOfDofPerNodeTmp = 12 ! Modify if inside PML
    endif   
end subroutine setNumDof

subroutine setSurfaceStation(nodeXyzIndex, nodeCoor, xline, yline, nodeCount, x4ndsSnapZ, x4ndsZValid)
    ! pathway_forward.md item 94 (owner ruling 2026-09-24): the depth test
    ! below compares against x4ndsSnapZ(i), the requested depth CLAMPED to
    ! the nearest node of the full z grid (computed once in meshgen, before
    ! the node loop), not the raw request x4nds(3,i) -- so a station whose
    ! requested depth is not itself a grid z-plane still matches, at the
    ! nearest node, instead of matching nothing. x and y are unchanged: they
    ! already snap to the nearest interior node below.
    !
    ! Finding 1 (row 94 audit): x4ndsZValid(i) is .false. for a station whose
    ! requested depth falls outside the physical (non-PML) clamp band (see
    ! meshgen's computation of x4ndsSnapZ/x4ndsZValid, above) -- for such a
    ! station x4ndsSnapZ(i) was never assigned a meaningful value, so the
    ! depth test below is gated on x4ndsZValid(i) first (short-circuit
    ! .and.): it can never match, and the station stays a genuine DROP.
    !
    ! Row 127 (this ruling): getLocalOneDimCoorArrAndSize overlaps every
    ! rank's local grid with its neighbour by exactly one node on EVERY axis
    ! (local index 1 of rank r == local index n of rank r-1, for x, y, and,
    ! when npz>1, z too -- confirmed algebraically from that subroutine's
    ! index formula). Two defects follow from treating that shared node
    ! asymmetrically:
    !   (a) DROP: the old Part1/2/3 below required "iy>1 .and. iy<ny" in
    !       ALL THREE x-branches -- there was no branch at all for iy==1 or
    !       iy==ny, so a station on a y-partition seam matched on NO rank
    !       (tpv8 station 11 at (2,2,1): 14 body files instead of 15).
    !   (b) DUPLICATE: ix==1/ix==nx already had branches, but neither
    !       excluded the case where the OTHER rank sharing that seam node
    !       also matches -- two ranks writing the same station file.
    ! Fix, applied identically on x, y and z: the rank with the LOWER MPI
    ! coordinate along an axis owns a station that lands on the seam it
    ! shares with its higher-coordinate neighbour. Concretely: this rank's
    ! local index 1 on an axis is a candidate only when this rank has no
    ! lower neighbour on that axis (its MPI coordinate is 0); this rank's
    ! local index n (nx/ny/nz) is ALWAYS a candidate -- it is either the
    ! true global high edge, or the seam with a higher-coordinate neighbour,
    ! which this rank owns as the lower of the pair. Interior nodes (1<idx<n)
    ! are never shared and need no gate. z has no separate index branch
    ! (depth is matched by value against x4ndsSnapZ, not by iz), so its
    ! ownership gate is a single early return covering the whole node.
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: nodeXyzIndex(10), ix, iy, iz, nx, ny, i
    integer (kind = 8) :: nodeCount
    integer (kind = 4) :: mex, mey, mez
    real (kind = dp) :: nodeCoor(10), xline(nodeXyzIndex(4)), yline(nodeXyzIndex(5))
    real (kind = dp) :: x4ndsSnapZ(max(1,totalNumOfOffSt))
    logical :: x4ndsZValid(max(1,totalNumOfOffSt))
    logical :: xMatch, yMatch

    ix = nodeXyzIndex(1)
    iy = nodeXyzIndex(2)
    iz = nodeXyzIndex(3)
    nx = nodeXyzIndex(4)
    ny = nodeXyzIndex(5)

    call calcXyzMPIId(mex, mey, mez)

    ! z-axis ownership gate (see header comment): skip this whole node if a
    ! lower-mez neighbour owns the z-seam it sits on.
    if (iz==1 .and. mez/=0) return

    do i=1,totalNumOfOffSt
        if(n4yn(i)/=0) cycle
        if (.not. (x4ndsZValid(i) .and. abs(nodeCoor(3)-x4ndsSnapZ(i))<tol)) cycle

        !x-axis match, mirroring the old Part1/2/3, gated on x-ownership.
        xMatch = .false.
        if (ix>1 .and. ix<nx) then
            if(abs(nodeCoor(1)-x4nds(1,i))<tol .or.&
            (x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
            (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1))) .or. &
            (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
            (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) xMatch = .true.
        elseif (ix==1) then
            if (mex==0) then
                if(abs(nodeCoor(1)-x4nds(1,i))<tol .or. &
                (x4nds(1,i)>nodeCoor(1).and.x4nds(1,i)<xline(ix+1).and. &
                (x4nds(1,i)-nodeCoor(1))<(xline(ix+1)-x4nds(1,i)))) xMatch = .true.
            endif
        elseif (ix==nx) then
            if(abs(nodeCoor(1)-x4nds(1,i))<tol .or. &
            (x4nds(1,i)>xline(ix-1).and.x4nds(1,i)<nodeCoor(1).and. &
            (nodeCoor(1)-x4nds(1,i))<(x4nds(1,i)-xline(ix-1)))) xMatch = .true.
        endif
        if (.not. xMatch) cycle

        !y-axis match: same structure as x, the row-127 fix (was missing
        !entirely for iy==1/iy==ny), gated on y-ownership.
        yMatch = .false.
        if (iy>1 .and. iy<ny) then
            if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
            (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
            (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1))) .or. &
            (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
            (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) yMatch = .true.
        elseif (iy==1) then
            if (mey==0) then
                if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
                (x4nds(2,i)>nodeCoor(2).and.x4nds(2,i)<yline(iy+1).and. &
                (x4nds(2,i)-nodeCoor(2))<(yline(iy+1)-x4nds(2,i)))) yMatch = .true.
            endif
        elseif (iy==ny) then
            if(abs(nodeCoor(2)-x4nds(2,i))<tol .or. &
            (x4nds(2,i)>yline(iy-1).and.x4nds(2,i)<nodeCoor(2).and. &
            (nodeCoor(2)-x4nds(2,i))<(x4nds(2,i)-yline(iy-1)))) yMatch = .true.
        endif
        if (.not. yMatch) cycle

        n4yn(i) = 1
        numOfOffFaultStCount = numOfOffFaultStCount + 1
        OffFaultStNodeIdIndex(1,numOfOffFaultStCount) = i
        OffFaultStNodeIdIndex(2,numOfOffFaultStCount) = nodeCount
        exit     !if node found, jump out the loop
    enddo
end subroutine setSurfaceStation

subroutine reduceOffFaultStationCoor(actualCoorHere, actualCoorGlobal)
    ! pathway_forward.md item 94 (owner ruling 2026-09-24: clamp station
    ! depth to the nearest node). The REAL(dp) MPI_Allreduce for the snap
    ! report (checkOffFaultStationCoverage/report_dropped_offfault_st,
    ! eqdyna3d.f90/library_output.f90) lives HERE, in a file with no other
    ! MPI_Allreduce call, deliberately: eqdyna3d.f90's own MPI_Allreduce
    ! calls are all LOGICAL (checkFaultMPIAlignment/checkOffFaultStation-
    ! Coverage/checkOnFaultStationCoverage), and gfortran's no-explicit-
    ! interface argument check refuses two calls to the SAME external name
    ! IN ONE FILE whose buffers differ in type or rank (see that
    ! subroutine's own header comment). It cannot live in library_output.f90
    ! either: testsys/regression/test_station_header_*.py compiles
    ! globalvar.f90+library_output.f90 alone with plain gfortran (no MPI
    ! wrapper, no mpif.h on the include path) to check the station writers
    ! in isolation, and an `include 'mpif.h'` there breaks that compile.
    ! meshgen.f90 carries no such isolated-compile test and already
    ! `include`s mpif.h in this same file (the `meshgen` subroutine, for its
    ! MPI_sendrecv calls), so this is a safe, uncontested home for it.
    ! MAX over a very-negative sentinel (set by the caller on every rank
    ! that did not match a given station) picks up whichever rank(s) did;
    ! safe because no real model coordinate is anywhere near -1e30, and a
    ! station matched by more than one rank (the documented MPI-boundary
    ! caveat, testsys/matrix.py) lands on the identical bit value on each of
    ! them, so MAX picks it either way.
    use globalvar
    implicit none
    include 'mpif.h'
    real (kind = dp), intent(in)  :: actualCoorHere(3,totalNumOfOffSt)
    real (kind = dp), intent(out) :: actualCoorGlobal(3,totalNumOfOffSt)
    integer (kind = 4) :: iMPIerr

    call MPI_Allreduce(actualCoorHere(1,1), actualCoorGlobal(1,1), 3*totalNumOfOffSt, &
        MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, iMPIerr)
end subroutine reduceOffFaultStationCoor

subroutine setEquationNumber(nodeXyzIndex, nodeCoor, eqNumIndexArrLocTag, equationNumCount, numOfDofPerNodeTmp)
    use globalvar 
    use errorCodes
    implicit none
    integer (kind = 4) :: iDof, numOfDofPerNodeTmp, nodeXyzIndex(10)
    integer (kind = 8) :: eqNumIndexArrLocTag, equationNumCount
    real (kind = dp) :: nodeCoor(10)
    
    do iDof = 1,numOfDofPerNodeTmp
        if(abs(nodeCoor(1)-xmin)<tol.or.abs(nodeCoor(1)-xmax)<tol.or.abs(nodeCoor(2)-ymin)<tol &
        .or.abs(nodeCoor(2)-ymax)<tol.or.abs(nodeCoor(3)-zmin)<tol) then
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag+1
            eqNumIndexArr(eqNumIndexArrLocTag) = -1 
            ! Dof = -1 for fixed boundary nodes, no eq #. 
        else
            equationNumCount      = equationNumCount + 1
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag+1
            eqNumIndexArr(eqNumIndexArrLocTag) = equationNumCount
            
            !Count # of DOF on MPI boundaries
            if (nodeXyzIndex(1)==1) then !Left
                numcount(4)=numcount(4)+1
            endif
            if (nodeXyzIndex(1)==nodeXyzIndex(4)) then !Right
                numcount(5)=numcount(5)+1
            endif
            if (nodeXyzIndex(2)==1) then !Front
                numcount(6)=numcount(6)+1
            endif
            if (nodeXyzIndex(2)==nodeXyzIndex(5)) then !Back
                numcount(7)=numcount(7)+1
            endif
            if (nodeXyzIndex(3)==1) then !Down
                numcount(8)=numcount(8)+1
            endif
            if (nodeXyzIndex(3)==nodeXyzIndex(6)) then !Up
                numcount(9)=numcount(9)+1
            endif
        endif
    enddo
    
end subroutine setEquationNumber


subroutine createElement(elemCount, stressDofCount, iy, iz, elementCenterCoor)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8) :: elemCount, stressDofCount
    integer (kind = 4) :: iy, iz, i, j 
    real (kind = dp) :: elementCenterCoor(3)
    
    elementCenterCoor = 0.0d0 
    
    elemCount        = elemCount + 1
    elemTypeArr(elemCount)    = 1 
    nodeElemIdRelation(1,elemCount) = plane1(iy-1,iz-1)
    nodeElemIdRelation(2,elemCount) = plane2(iy-1,iz-1)
    nodeElemIdRelation(3,elemCount) = plane2(iy,iz-1)
    nodeElemIdRelation(4,elemCount) = plane1(iy,iz-1)
    nodeElemIdRelation(5,elemCount) = plane1(iy-1,iz)
    nodeElemIdRelation(6,elemCount) = plane2(iy-1,iz)
    nodeElemIdRelation(7,elemCount) = plane2(iy,iz)
    nodeElemIdRelation(8,elemCount) = plane1(iy,iz)
    stressCompIndexArr(elemCount)   = stressDofCount
    
    do i=1,nen
        do j=1,3
            elementCenterCoor(j) = elementCenterCoor(j) + meshCoor(j,nodeElemIdRelation(i,elemCount))
        enddo
    enddo
    elementCenterCoor = elementCenterCoor/8.0d0
    
    if (elementCenterCoor(1)>PMLb(1).or.elementCenterCoor(1)<PMLb(2) &
    .or.elementCenterCoor(2)>PMLb(3).or.elementCenterCoor(2)<PMLb(4) &
    .or.elementCenterCoor(3)<PMLb(5)) then
        elemTypeArr(elemCount) = 2
        stressDofCount        = stressDofCount+15+6
    else
        stressDofCount        = stressDofCount+12
    endif
    
    ! assign nodeElemIdRelation(elemCount), ids(elemCount), et(elemCount)
    ! return elementCenterCoor, 
    ! update elemCount, stressDofCount
end subroutine createElement

 subroutine replaceSlaveWithMasterNode(nodeCoor, elemCount, nftnd0)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8) :: elemCount
    integer (kind = 4) :: iFault, iFaultNodePair, nftnd0(ntotft), k
    real (kind = dp) :: nodeCoor(10)
    logical :: isAboveSomeFaultPlane
    ! The default grids only contain slave nodes.
    ! This subroutine will replace slave nodes with corresponding master nodes.

    ! Row 17 (multi-fault): the y-test used to be `nodeCoor(2)>0.0d0 .and.
    ! abs(nodeCoor(2)-dy)<tol` -- one cell above fault 1's plane (y=0) alone.
    ! For ntotft>1 an element one cell above fault 2's plane (or any fault's)
    ! never matched, so its slave-node references were never swapped for
    ! master-node ones -- the mesh stayed disconnected across fault 2's
    ! split, which showed up downstream as a zero/singular mass entry and a
    ! NaN velocity at a fault-2 node within the first 2 steps (measured on
    ! test.multifault2, rank 1/3, nodes at fault 2's z-edges). Generalized to
    ! ANY fault's plane-plus-one-cell.
    !
    ! Row 17 rebased: x/z bounds USED to stay keyed to fault 1's box alone
    ! (checkInputConsistency.f90 used to require every fault share it) --
    ! now getLocalOneDimCoorArrAndSize meshes the UNION of every fault's x/z
    ! extent, so a fault whose x/z differs from fault 1's (TPV22/23's two
    ! faults, e.g.) needs its OWN x/z tested alongside its OWN y-plane, one
    ! fault at a time, not fault 1's x/z paired with any(...) fault's y.
    ! Reduces to the old any(...)-across-y test bit-for-bit at ntotft==1 (one
    ! fault, so "test fault i's own x/z and y" IS "test fault 1's x/z and y").
    isAboveSomeFaultPlane = .false.
    do iFault = 1, ntotft
        if (nodeCoor(1)>(fltxyz(1,1,iFault)-tol) .and. nodeCoor(1)<(fltxyz(2,1,iFault)+dx+tol) .and. &
            nodeCoor(3)>(fltxyz(1,3,iFault)-tol) .and. &
            abs(nodeCoor(2) - (fltxyz(1,2,iFault) + dy)) < tol) then
            isAboveSomeFaultPlane = .true.
            exit
        endif
    enddo
    if ((elemTypeArr(elemCount) == 1 .and. isAboveSomeFaultPlane) &
         .or. elemTypeArr(elemCount)==12 .or. elemTypeArr(elemCount)==13 ) then
        do iFault = 1, ntotft
            do iFaultNodePair = 1, nftnd0(iFault)
                do k = 1,nen
                    if(nodeElemIdRelation(k,elemCount)==nsmp(1,iFaultNodePair,iFault)) then
                        nodeElemIdRelation(k,elemCount) = nsmp(2,iFaultNodePair,iFault)  !use master node for the node!
                    endif
                enddo
            enddo
        enddo
    endif      
    
end subroutine replaceSlaveWithMasterNode

subroutine checkIsOnFault(nodeCoor, iFault, isOnFault)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: isOnFault, iFault
    real (kind = dp) :: nodeCoor(10), distToFault, tangentDip, x0, pNorm
    isOnFault = 0
 
    if(nodeCoor(1)>=(fltxyz(1,1,iFault)-tol).and.nodeCoor(1)<=(fltxyz(2,1,iFault)+tol).and. &
        nodeCoor(2)>=(fltxyz(1,2,iFault)-tol).and.nodeCoor(2)<=(fltxyz(2,2,iFault)+tol).and. &
        nodeCoor(3)>=(fltxyz(1,3,iFault)-tol) .and. nodeCoor(3)<=(fltxyz(2,3,iFault)+tol)) then
        ! Row 17 (multi-fault): this used to test nodeCoor(2)==0.0d0 by EXACT
        ! float equality -- correct only for a single fault pinned at the
        ! global y=0 plane. Generalized to the SAME test against this
        ! fault's own y-plane (fltxyz(1,2,iFault), which is fymin(iFault);
        ! a planar vertical fault has fymin==fymax so either bound names the
        ! same plane). abs(...)<tol is a superset of the old exact test (0 is
        ! trivially within tol of itself) and tol=1e-5 is far smaller than any
        ! dy in use, so this is bit-identical at ntotft=1, fault 1 at y=0.
        ! The rough-fault y-blend (ycoort/meshCoor) is untouched -- this test
        ! runs against the UNBLENDED nodeCoor(2), exactly as before.
        ! Row 153 checkpoint 1: per-fault faultDegenStyle(iFault)/
        ! faultDegenAngle(iFault) replace the bare C_degen test -- derived
        ! from C_degen identically for every fault (readfaultgeometry), so
        ! this branch is bit-identical to the old C_degen==0.0d0/C_degen>3.
        ! test on every existing case.
        if (faultDegenStyle(iFault)==0 .and. abs(nodeCoor(2) - fltxyz(1,2,iFault)) < tol) then
            isOnFault = 1
        elseif (faultDegenStyle(iFault)==1) then
            if (fltxyz(1,2,iFault)>=fltxyz(2,2,iFault)) write(*,*) 'ymax should be > ymin. Wrong geo, exit'
            ! Row 153 checkpoint 2a audit fix: offset by THIS fault's own
            ! y-origin (fltxyz(1,2,iFault) = fymin) before applying the dip
            ! tilt -- previously hardcoded to assume the dip plane's trace
            ! sits at y=0. Every committed dipping-fault case (TPV36/37) has
            ! fymin==0.0d0 (see tpv36_37_common.py), so subtracting it is a
            ! no-op there: bit-identical. A second dipping fault with
            ! fymin!=0 (none registered yet) is now handled correctly too.
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            distToFault = abs(nodeCoor(3)+(nodeCoor(2)-fltxyz(1,2,iFault))*tangentDip)
            distToFault = distToFault/(1.d0+tangentDip**2)**0.5
            if (distToFault < dx/100.d0) isOnFault = 1
        elseif (faultDegenStyle(iFault)==2) then
            ! Style 2: TPV24/25 branch fault, x-y strike tilt, z-extruded.
            ! x0 is derived from this fault's own box (fxmin - dx, the one
            ! cell gap that excludes the junction column -- see wedge()'s
            ! matching comment in library_degeneration.f90).
            tangentDip = dtan(faultDegenAngle(iFault)/180.d0*pi)
            pNorm = (1.d0+tangentDip**2)**0.5
            x0 = fltxyz(1,1,iFault) - dx
            distToFault = abs(nodeCoor(2)+(nodeCoor(1)-x0)*tangentDip)/pNorm
            if (distToFault < dx/100.d0) isOnFault = 1
    endif
    endif
end subroutine checkIsOnFault

subroutine createMasterNode(nodeXyzIndex, nxuni, nzuni, nodeCoor, ycoort, nodeCount, msnode, nftnd0, equationNumCount, eqNumIndexArrLocTag,&
                            pfx, pfz, ixfi, izfi, ifs, ifd, fltrc)
use globalvar 
use errorCodes
implicit none
integer (kind = 4) :: iFault, iFaultNodePair, isOnFault, nftnd0(ntotft), i, nxuni, nzuni
integer (kind = 8) :: nodeCount, msnode, equationNumCount, eqNumIndexArrLocTag
integer (kind = 8) :: fltrc(2,nxuni,nzuni,ntotft)
integer (kind = 4) :: ixfi(ntotft), izfi(ntotft), ifs(ntotft), ifd(ntotft), nodeXyzIndex(10)
integer (kind = 4) :: mex, mey, mez
logical :: isOnFaultStationOwner
real (kind = dp) :: nodeCoor(10), ycoort, pfx, pfz

! Row 131: setOnFaultStation (below) matches a fault node by VALUE only,
! with no ix/iz-keyed gate at all. getLocalOneDimCoorArrAndSize overlaps
! every rank's local grid with its neighbour by exactly one node on every
! axis (same fact row 127 used for setSurfaceStation), so a fault node
! sitting exactly on a shared x seam (npx>1) or z seam (npz>1) is
! physically present, and isOnFault==1, on BOTH ranks -- and both write
! the same faultst* file (last-writer-wins race; live on tpv8 at (2,2,1),
! ranks 0 and 2 both write rank 0's four faultst files). This must NOT
! change which rank owns the split-node PAIR itself (nsmp/nftnd0/msnode/
! fltgm/un/us/ud below are untouched -- every rank still builds its own
! local copy for the solve/frt output); it only gates who additionally
! counts/writes the STATION record for a shared-seam node. Same ownership
! rule as row 127: the lower-MPI-coordinate rank of a seam pair owns it
! (local index 1 on an axis is owned only when this rank has no lower
! neighbour there; local index n is always owned).
call calcXyzMPIId(mex, mey, mez)

do iFault = 1, ntotft

    isOnFault = 0 
    call checkIsOnFault(nodeCoor, iFault, isOnFault)

    if (isOnFault == 1) then
        nftnd0(iFault)                = nftnd0(iFault) + 1 ! # of split-node pairs + 1
        nsmp(1,nftnd0(iFault),iFault) = nodeCount              ! set Slave node nodeID to nsmp
        ! Row 17 (multi-fault): msnode used to be
        ! nx*ny*nz + nftnd0(iFault) -- a PER-FAULT sequence number reused as
        ! a GLOBAL node id. At ntotft>1 fault 2's msnode range (1..nftnd0(2))
        ! collided with fault 1's (1..nftnd0(1)) whenever nftnd0(2) <=
        ! nftnd0(1), aliasing two different physical master nodes onto one
        ! id (ERR 45, exit 45, previously a hard refusal below).
        !
        ! A per-fault BLOCK offset of a FIXED size (e.g. (iFault-1)*nftmx,
        ! tried and reverted here) does not work: nftmx is this RANK's max
        ! fault-node count over ALL faults, so whichever fault has the most
        ! nodes ON THIS RANK sets the block size for every fault's block --
        ! including faults that come before it. When a LATER fault has MORE
        ! local nodes than an EARLIER one, the earlier fault's block is
        ! oversized and the later fault's range overruns totalNumOfNodes
        ! (measured: SIGSEGV on test.multifault2 rank 1, msnode 66749 vs
        ! totalNumOfNodes 65788 -- fault 2 locally outnumbered fault 1).
        !
        ! Fixed with a running total across every fault's split-node pairs
        ! created SO FAR (sum(nftnd0), which already includes THIS
        ! increment): strictly increasing by exactly 1 per master node
        ! created, in MESH-SCAN order (which interleaves faults, since faults
        ! differ only in y and the scan visits all y for each x,z), and so
        ! never exceeds totalNumOfNodes - nx*ny*nz (= sum of every fault's
        ! FINAL local count). The one consumer that needs to invert this back
        ! to a fault+local-index pair -- addFaultBoundaryTerm,
        ! assembleGlobalMass.f90 -- does not reconstruct the formula; it
        ! looks the already-recorded id up from nsmp(2,i,ift) instead (the
        ! same mapping this line writes to nsmp two lines below), so the
        ! interleaving is invisible to it. Reduces to the old formula
        ! bit-for-bit at ntotft==1 (sum over a length-1 array).
        ! Item 143: gridNodeCount (globalvar.f90) does the nx*ny*nz product
        ! in 64-bit; written inline with kind=4 operands it wraps past 2^31.
        msnode                        = gridNodeCount(nodeXyzIndex(4), nodeXyzIndex(5), nodeXyzIndex(6)) + sum(nftnd0) ! create Master node at the end of regular grids

        eqNumStartIndexLoc(msnode) = eqNumIndexArrLocTag
        numOfDofPerNodeArr(msnode) = 3         
        nsmp(2,nftnd0(iFault),iFault) = msnode !set Master node nodeID to nsmp
        plane2(nodeXyzIndex(5)+iFault,nodeXyzIndex(3)) = msnode
        
        meshCoor(1,msnode) = nodeCoor(1)
        meshCoor(2,msnode) = nodeCoor(2)
        meshCoor(3,msnode) = nodeCoor(3)
        if (insertFaultType > 0) then 
            meshCoor(2,msnode) = ycoort
        endif
        
        !set Equation Numbers for the newly created Master Node.
        do i = 1, ndof
            equationNumCount      = equationNumCount + 1
            eqNumIndexArrLocTag      = eqNumIndexArrLocTag + 1
            eqNumIndexArr(eqNumIndexArrLocTag) = equationNumCount
        enddo
        
        ! Count split-node pair # for MPI. fltgm is a global var, now
        ! per-fault (Row 17: fltgm(nftnd0(iFault)) alone collided across
        ! faults whenever nftnd0(iFault) ran over the same range for two
        ! faults -- indexed by (localFaultNode, iFault) instead).
        !
        ! The fltnum(1..6) accumulation that used to run here in parallel is
        ! REMOVED, not just re-indexed: MPI4arn (below, called once per fault
        ! right after this loop finishes) unconditionally resets fltnum and
        ! recomputes it from fltgm over the same nodes, so this running
        ! total was never read before being overwritten -- true before this
        ! change too (see MPI4arn's own header comment) and unaffected by
        ! moving to a per-fault fltgm.
        if(nodeXyzIndex(1) == 1) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 1
        endif
        if(nodeXyzIndex(1) == nodeXyzIndex(4)) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 2
        endif
        if(nodeXyzIndex(2) == 1) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 10
        endif
        if(nodeXyzIndex(2) == nodeXyzIndex(5)) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 20
        endif
        if(nodeXyzIndex(3) == 1) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 100
        endif
        if(nodeXyzIndex(3) == nodeXyzIndex(6)) then
            fltgm(nftnd0(iFault),iFault) = fltgm(nftnd0(iFault),iFault) + 200
        endif

        ! setOnFaultStation -- gated (row 131) so a fault node on a shared
        ! x, y or z seam is claimed for STATION OUTPUT by exactly one rank
        ! (the lower-MPI-coordinate side), while still creating its own
        ! full split-node pair above, unconditionally, on both ranks.
        ! The y term matters when an npy boundary lands ON the fault plane
        ! (checkFaultMPIAlignment's DUPLICATE case): both ranks then hold
        ! every fault node, and without it both write every faultst* file.
        isOnFaultStationOwner = .not. ((nodeXyzIndex(1)==1 .and. mex/=0) .or. &
                                        (nodeXyzIndex(2)==1 .and. mey/=0) .or. &
                                        (nodeXyzIndex(3)==1 .and. mez/=0))
        if (isOnFaultStationOwner) then
        do i = 1, nonfs(iFault)
            if(abs(nodeCoor(1)-xonfs(1,i,iFault))<tol .and. &
                abs(nodeCoor(3)-xonfs(2,i,iFault))<tol) then
                numOfOnFaultStCount = numOfOnFaultStCount + 1
                anonfs(1,numOfOnFaultStCount) = nftnd0(iFault)
                anonfs(2,numOfOnFaultStCount) = i
                anonfs(3,numOfOnFaultStCount) = iFault
                exit
            endif
        enddo
        endif  

        ! set unit vectors to split-node pair    
        un(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(2,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))*dsin(fltxyz(2,4,iFault))
        un(3,nftnd0(iFault),iFault) = dcos(fltxyz(2,4,iFault))        
        us(1,nftnd0(iFault),iFault) = -dsin(fltxyz(1,4,iFault))
        us(2,nftnd0(iFault),iFault) = -dcos(fltxyz(1,4,iFault))
        us(3,nftnd0(iFault),iFault) = 0.0d0
        ud(1,nftnd0(iFault),iFault) = dcos(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(2,nftnd0(iFault),iFault) = dsin(fltxyz(1,4,iFault))*dcos(fltxyz(2,4,iFault))
        ud(3,nftnd0(iFault),iFault) = dsin(fltxyz(2,4,iFault))
        
        if (insertFaultType >0) then
            un(1,nftnd0(iFault),iFault) = -pfx/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(2,nftnd0(iFault),iFault) = 1.0d0/(pfx**2 + 1.0d0 + pfz**2)**0.5
            un(3,nftnd0(iFault),iFault) = -pfz/(pfx**2 + 1.0d0 + pfz**2)**0.5    
            us(1,nftnd0(iFault),iFault) = 1.0d0/(1.0d0 + pfx**2)**0.5
            us(2,nftnd0(iFault),iFault) = pfx/(1.0d0 + pfx**2)**0.5
            us(3,nftnd0(iFault),iFault) = 0.0d0
            ud(1,nftnd0(iFault),iFault) = us(2,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault) &
                - us(3,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault)
            ud(2,nftnd0(iFault),iFault) = us(3,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault) &
                - us(1,nftnd0(iFault),iFault)*un(3,nftnd0(iFault),iFault)
            ud(3,nftnd0(iFault),iFault) = us(1,nftnd0(iFault),iFault)*un(2,nftnd0(iFault),iFault) &
                - us(2,nftnd0(iFault),iFault)*un(1,nftnd0(iFault),iFault)
        endif                             
        
        !...prepare for area calculation
        if(ixfi(iFault)==0) ixfi(iFault)=nodeXyzIndex(1)
        if(izfi(iFault)==0) izfi(iFault)=nodeXyzIndex(3)
        ifs(iFault)=nodeXyzIndex(1)-ixfi(iFault)+1
        ifd(iFault)=nodeXyzIndex(3)-izfi(iFault)+1
        fltrc(1,ifs(iFault),ifd(iFault),iFault) = msnode    !master node
        fltrc(2,ifs(iFault),ifd(iFault),iFault) = nftnd0(iFault) !fault node num in sequence
    endif 
enddo 

end subroutine createMasterNode

subroutine createNode(nodeCoor, xcoor, ycoor, zcoor, nodeCount, nodeXyzIndex)
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 8) :: nodeCount
    integer (kind = 4) :: nodeXyzIndex(10), iy, iz
    real (kind = dp) :: nodeCoor(10), xcoor, ycoor, zcoor
    iy = nodeXyzIndex(2)
    iz = nodeXyzIndex(3)
    nodeCoor(1)   = xcoor
    nodeCoor(2)   = ycoor
    nodeCoor(3)   = zcoor
    nodeCount         = nodeCount + 1        
    plane2(iy,iz) = nodeCount
    meshCoor(1,nodeCount) = nodeCoor(1)
    meshCoor(2,nodeCount) = nodeCoor(2)
    meshCoor(3,nodeCount) = nodeCoor(3) 
end subroutine createNode

! insertFaultInterface moved to func_lib.f90.

subroutine initializeNodeXyzIndex(ix, iy, iz, nx, ny, nz, nodeXyzIndex)
    use globalvar 
    use errorCodes
    implicit none 
    integer (kind = 4) :: ix, iy, iz, nx, ny, nz, nodeXyzIndex(10)
    nodeXyzIndex(1) = ix
    nodeXyzIndex(2) = iy
    nodeXyzIndex(3) = iz
    nodeXyzIndex(4) = nx
    nodeXyzIndex(5) = ny
    nodeXyzIndex(6) = nz
end subroutine initializeNodeXyzIndex

subroutine setPlasticStress(depth, elemCount)
    use globalvar
    use errorCodes
    implicit none
    
    real(kind = dp) :: depth, vstmp, vptmp, routmp, strVert, devStr
    real(kind = dp) :: devStrDepthTaper
    integer(kind = 8) :: elemCount
    integer(kind = 4) :: etTag

    etTag = 0
    if (elemTypeArr(elemCount)==2) etTag = 1 ! adjustment for PML elements
    
    eleporep(elemCount) = 0.0d0  !rhow*tmp2*gama  !pore pressure>0
    strVert            = -(roumax- rhow*(gamar+1.0d0))*depth*grav ! should be negative   
    ! devStrDepthTaper (func_lib.f90) is SCEC TPV29/30's Omega(depth): the
    ! deviatoric component tapers to zero over a depth interval while the
    ! vertical component keeps growing. It returns exactly 1.0d0 when the
    ! taper is not configured, so this line is bit-for-bit the pre-v5.9.0
    ! `abs(strVert)*devStrToStrVertRatio` for every such case.
    devStr             = abs(strVert)*devStrToStrVertRatio*devStrDepthTaper(depth) ! positive
    
    stressArr(stressCompIndexArr(elemCount)+3+15*etTag) = strVert
    stressArr(stressCompIndexArr(elemCount)+1+15*etTag) = strVert - devStr*dcos(2.0d0*str1ToFaultAngle)
    stressArr(stressCompIndexArr(elemCount)+2+15*etTag) = strVert + devStr*dcos(2.0d0*str1ToFaultAngle)
    stressArr(stressCompIndexArr(elemCount)+6+15*etTag) = devStr*dsin(2.0d0*str1ToFaultAngle)
    if (stressArr(stressCompIndexArr(elemCount)+2+15*etTag) >= 0.0d0) write(*,*) 'WARNING: positive Sigma3 ... ...'
    
end subroutine setPlasticStress
