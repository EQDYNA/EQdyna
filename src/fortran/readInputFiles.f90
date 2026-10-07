
subroutine readglobal
! This subroutine is read information from bglobal.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, ios


    call requireInputFile("bGlobal.txt")
    
    open(unit = 1001, file = 'bGlobal.txt', form = 'formatted', status = 'old')
        read(1001,*) mode
        read(1001,*) C_elastic
        read(1001,*) C_nuclea
        read(1001,*) C_degen
        ! Row 153 checkpoint 1 audit fix: checkIsOnFault (meshgen.f90) takes
        ! neither its C_degen==0 nor its C_degen>3.0d0 branch for 0<C_degen<=3
        ! -- isOnFault stays 0 for every node, which this refactor's
        ! faultDegenStyle/faultDegenAngle derivation (readfaultgeometry,
        ! below) would otherwise silently map to style=0 (vertical planar
        ! fault) instead of refusing. Python's meshgen.py:_check_is_on_fault_vec
        ! raises NotImplementedError for the same range; match it here.
        if (C_degen > 0.0d0 .and. C_degen <= 3.0d0) call abortRun(ERR_GEOM_DEGEN_UNSUPPORTED, &
            'bGlobal.txt: C_degen is in (0, 3], a value checkIsOnFault takes neither its ' // &
            '==0 nor its >3 branch for -- no fault node would be found. Use C_degen=0 ' // &
            '(vertical planar fault) or C_degen>3 (wedge-degeneration dip angle in degrees).')
        read(1001,*) insertFaultType
        read(1001,*) friclaw
        read(1001,*) ntotft
        read(1001,*) nucfault
        read(1001,*) TPV
        read(1001,*) output_plastic
        read(1001,*) outputGroundMotion
        read(1001,*) outputFinalSurfDisp
        read(1001,*) 
        read(1001,*) npx, npy, npz
        read(1001,*)
        read(1001,*) totalSimuTime
        read(1001,*) dt 
        read(1001,*)
        read(1001,*) nmat, n2mat
        read(1001,*) roumax, rhow, gamar
        read(1001,*) rdampk, vmaxPML
        read(1001,*) 
        read(1001,*) xsource, ysource, zsource
        read(1001,*) nucR, nucRuptVel, nucdtau0, nucT
        read(1001,*) str1ToFaultAngle, devStrToStrVertRatio
        read(1001,*) bulk, coheplas
        read(1001,*) fstrike, fdip
        read(1001,*) slipRateThres

        ! Viscoplastic / plastic-output block. Every entry here used to be a
        ! constant compiled into the solver: tv was 2*dz/NUC_VS_FIXED
        ! (readInputFiles.f90), the deviatoric pre-stress had no depth taper at
        ! all (meshgen.f90's setPlasticStress), and the plastic-strain output
        ! window was a fixed 5x2x8 km box (library_output.f90). They are read,
        ! never defaulted here: case.setup writes all three unconditionally, so
        ! a file without them is a STALE file and says so (rule 2).
        read(1001,*,iostat=ios)
        if (ios /= 0) call stopStaleGlobal('the viscoplastic/plastic-output block separator')
        read(1001,*,iostat=ios) tv
        if (ios /= 0) call stopStaleGlobal('the viscoplastic relaxation time Tv, s (par.viscoplasticRelaxTime)')
        read(1001,*,iostat=ios) devStrTaperDepthStart, devStrTaperDepthEnd
        if (ios /= 0) call stopStaleGlobal('the deviatoric pre-stress taper depths, m positive down (par.devStrTaperDepthStart/End)')
        read(1001,*,iostat=ios) (plasticOutputHalfWidth(i), i = 1, 3)
        if (ios /= 0) call stopStaleGlobal('the plastic-strain output window half-widths, m (par.plasticOutputHalfWidth)')
        ! Station normal-stress sign convention (board row 22a). SCEC specs do
        ! not agree on it -- TPV29/30, TPV10/11 and TPV36/37 say "Positive means
        ! extension", TPV103/104 and TPV105-3D say "Positive means
        ! compression" -- so it is a per-case input, not a constant. Until
        ! row 22a library_output.f90 negated column 8 for EVERY case, i.e.
        ! wrote compression-positive everywhere.
        read(1001,*,iostat=ios) nStressOutSign
        if (ios /= 0) call stopStaleGlobal('the station normal-stress sign convention, +1 extension / -1 compression (par.faultStNormalStressSign)')
        if (nStressOutSign /= 1 .and. nStressOutSign /= -1) call abortRun(ERR_CFG_NSTRESS_SIGN_INVALID, &
            'bGlobal.txt station normal-stress sign must be +1 (extension) or -1 (compression); re-run case.setup.')

    close(1001)
    str1ToFaultAngle = str1ToFaultAngle*pi/180.0d0 !convert degrees to radian

end subroutine readglobal

subroutine stopStaleGlobal(what)
! Stop loudly on a bGlobal.txt written by an older case.setup than this binary
! (PROJECT_RULES.md rule 2: no substituted default for absent input). Every
! rank reads the same file and reaches the same verdict; only the master
! prints, then all ranks stop.
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    character (len=*) :: what
    integer (kind = 4) :: iMPIerr

    if (me == masterProcsId) then
        write(*,*) 'readglobal: bGlobal.txt ends before ', trim(what), '.'
        write(*,*) '  This file was written by an older case.setup than this binary.'
        write(*,*) '  Re-run case.setup in this directory to regenerate it.'
    endif
    call MPI_Barrier(MPI_COMM_WORLD, iMPIerr)
    call abortRun(ERR_INPUT_FILE_STALE, &
        'bGlobal.txt ends before '//trim(what)//'; re-run case.setup.')
end subroutine stopStaleGlobal
! #2 readmodelgeometry -------------------------------------------------
subroutine readmodelgeometry
! This subroutine is read information from bglobal.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    
    call requireInputFile("bModelGeometry.txt")
    
    open(unit = 1002, file = 'bModelGeometry.txt', form = 'formatted', status = 'old')
        read(1002,*) xmin, xmax
        read(1002,*) ymin, ymax
        read(1002,*) zmin, zmax
        read(1002,*) 
        read(1002,*) dis4uniF, dis4uniB
        read(1002,*) rat
        read(1002,*) dx, dy, dz 
    close(1002)
end subroutine readmodelgeometry

! #3 readfaultgeometry
subroutine readfaultgeometry
! This subroutine is read information from bFault_Geometry.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4)::i
    
    call requireInputFile("bFaultGeometry.txt")
    
    open(unit = 1003, file = 'bFaultGeometry.txt', form = 'formatted', status = 'old')
        do i = 1, ntotft
            read(1003,*) 
            read(1003,*) fxmin(i), fxmax(i)
            read(1003,*) fymin(i), fymax(i)
            read(1003,*) fzmin(i), fzmax(i)
        enddo
    close(1003)

    do i = 1, ntotft
        fltxyz(1,1,i)=fxmin(i)
        fltxyz(2,1,i)=fxmax(i)
        fltxyz(1,2,i)=fymin(i)
        fltxyz(2,2,i)=fymax(i)
        fltxyz(1,3,i)=fzmin(i)
        fltxyz(2,3,i)=fzmax(i)
        fltxyz(1,4,i)=fstrike*pi/180.0d0
        if (C_degen>3.0d0) then
            fltxyz(2,4,i) = C_degen*pi/180.0d0
        else
            fltxyz(2,4,i) = 90.d0*pi/180.d0
        endif
        ! Row 153 checkpoint 1: per-fault degeneration style/angle, derived
        ! from C_degen exactly as fltxyz(2,4,i) above -- every fault gets the
        ! SAME style/angle C_degen already gave it (uniform test), so this is
        ! a pure refactor, not a behavior change.
        if (C_degen>3.0d0) then
            faultDegenStyle(i) = 1
            faultDegenAngle(i) = C_degen
        else
            faultDegenStyle(i) = 0
            faultDegenAngle(i) = 0.d0
        endif
    enddo
    
end subroutine readfaultgeometry

! #4 readmaterial --------------------------------------------------------
subroutine readmaterial
! This subroutine is read information from bMaterial.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    call requireInputFile("bMaterial.txt")
    
    open(unit = 1004, file = 'bMaterial.txt', form = 'formatted', status = 'old')
        do i = 1, nmat
            read(1004,*) (material(i,j), j = 1, n2mat)
        enddo
    close(1004)

    if (n2mat == 5) call checkTwoSidedMaterialTable
    if (n2mat == 6) call buildMaterialGrid3D

    ccosphi = coheplas*dcos(atan(bulk))
    sinphi  = dsin(atan(bulk))
    nstep   = idnint(totalSimuTime/dt)
    rdampk  = rdampk*dt
    ! tv is no longer derived here -- it is read from bGlobal.txt (readglobal).
end subroutine readmaterial

subroutine checkTwoSidedMaterialTable
! Validate the n2mat==5 two-sided 1D material table (meshgen.f90's
! setElementMaterial third branch, introduced for SCEC TPV35): the side
! column decides a layer set by the sign of (element-centre y - the fault
! y-plane), which only means something when every fault lies on ONE common
! vertical plane (fymin == fymax, identical across faults); the
! per-side first-match lookup is only reproducible for strictly ascending
! bottoms; and a side with no rows would leave every element on it without
! a material (caught later as ERR_MESH_MATERIAL_UNSET, but with no hint why).
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: i, s, nrows
    real (kind = dp) :: lastBottom

    if (nmat < 2) call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
        'bMaterial.txt: a two-sided (n2mat=5) table needs nmat >= 2 rows (at least one per side).')
    if (maxval(fltxyz(2,2,:)) - minval(fltxyz(1,2,:)) > tol) call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
        'bMaterial.txt: the two-sided (n2mat=5) material table takes sides against ONE vertical y-plane; every fault must have fymin == fymax and share the same y.')
    if (C_degen /= 0) call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
        'bMaterial.txt: the two-sided (n2mat=5) material table needs a vertical planar fault (C_degen=0); a dipping/degenerate fault has no single y-plane to take sides against.')
    do s = -1, 1, 2
        nrows = 0
        lastBottom = -1.0d0
        do i = 1, nmat
            if (nint(material(i,5)) /= -1 .and. nint(material(i,5)) /= 1) &
                call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
                'bMaterial.txt: column 5 of a two-sided (n2mat=5) table must be -1 (y below the fault plane) or +1 (y above it).')
            if (nint(material(i,5)) /= s) cycle
            nrows = nrows + 1
            if (material(i,1) <= lastBottom) call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
                'bMaterial.txt: layer bottoms within one side of a two-sided (n2mat=5) table must be strictly ascending and positive.')
            lastBottom = material(i,1)
        enddo
        if (nrows == 0) call abortRun(ERR_CFG_MATERIAL_TABLE_INVALID, &
            'bMaterial.txt: a two-sided (n2mat=5) table must have at least one row for each side (-1 and +1).')
    enddo
end subroutine checkTwoSidedMaterialTable

subroutine buildMaterialGrid3D
! Build the n2mat==6 3D structured material grid (SCEC TPV34: CVM-H
! sampled at the mesh's uniform element-centre spacing) from bMaterial.txt
! rows [x y z vp vs rho] (m, m/s, kg/m3; EQdyna frame, z <= 0 underground).
! The grid is self-describing: per axis the origin is the smallest
! coordinate, the spacing the smallest positive offset from it (1 m for a
! single-plane axis), the count nint((max-min)/spacing)+1. nmat must equal
! nx*ny*nz, every cell must be filled exactly once by an on-grid row, and
! vp/vs/rho must be positive (ERR_CFG_MATERIAL_GRID_INVALID otherwise).
! meshgen.f90's setElementMaterial then reads the NEAREST cell to each
! element centre, clamped to the grid: piecewise constant, never
! interpolated (rule 17 step 2) -- an element centre on the grid (every
! uniform-belt element when the grid spacing equals dx) reads its own
! sample exactly.
    use globalvar
    use errorCodes
    implicit none
    integer (kind = 4) :: i, k, idx(3)
    integer (kind = 4), allocatable :: filled(:,:,:)
    real (kind = dp) :: cmin(3), cmax(3), off

    if (nmat < 2) call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
        'bMaterial.txt: a 3D material grid (n2mat=6) needs nmat >= 2 rows.')
    do k = 1, 3
        cmin(k) = minval(material(:,k))
        cmax(k) = maxval(material(:,k))
        matGridSpacing(k) = huge(1.0d0)
        do i = 1, nmat
            off = material(i,k) - cmin(k)
            if (off > tol .and. off < matGridSpacing(k)) matGridSpacing(k) = off
        enddo
        if (matGridSpacing(k) == huge(1.0d0)) matGridSpacing(k) = 1.0d0   ! single-plane axis
        matGridOrigin(k) = cmin(k)
        matGridCount(k) = nint((cmax(k) - cmin(k))/matGridSpacing(k)) + 1
    enddo
    if (product(matGridCount) /= nmat) call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
        'bMaterial.txt: the 3D material grid (n2mat=6) rows do not form a complete uniform nx*ny*nz block (nmat /= nx*ny*nz).')

    ! The grid must cover the DECLARED mesh box (xmin/xmax/ymin/ymax/zmin/
    ! zmax from bModelGeometry.txt, already read by readmodelgeometry before
    ! this subroutine runs) at cell-centre granularity: the grid's
    ! nearest-neighbour reach is [origin - spacing/2, origin +
    ! spacing*(count-1) + spacing/2] per axis, and that reach must contain
    ! [box_min, box_max]. A stretched/PML element legitimately lands OFF the
    ! grid and clamps to the nearest sample -- setElementMaterial's n2mat==6
    ! branch does that deliberately and is unaffected by this check. What
    ! this check catches is the grid not even covering the declared domain,
    ! which would silently clamp core/interior elements (the ones the
    ! physics depends on) to an edge cell with no indication.
    if (matGridOrigin(1) - 0.5d0*matGridSpacing(1) > xmin + tol .or. &
        matGridOrigin(1) + matGridSpacing(1)*(matGridCount(1)-1) + 0.5d0*matGridSpacing(1) < xmax - tol .or. &
        matGridOrigin(2) - 0.5d0*matGridSpacing(2) > ymin + tol .or. &
        matGridOrigin(2) + matGridSpacing(2)*(matGridCount(2)-1) + 0.5d0*matGridSpacing(2) < ymax - tol .or. &
        matGridOrigin(3) - 0.5d0*matGridSpacing(3) > zmin + tol .or. &
        matGridOrigin(3) + matGridSpacing(3)*(matGridCount(3)-1) + 0.5d0*matGridSpacing(3) < zmax - tol) &
        call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
        'bMaterial.txt: the 3D material grid (n2mat=6) does not cover the mesh box '// &
        '(xmin/xmax/ymin/ymax/zmin/zmax from bModelGeometry.txt) -- a core element '// &
        'outside the grid would silently clamp to the edge cell instead of reading '// &
        'the right material.')
    allocate(matGrid3D(3, matGridCount(1), matGridCount(2), matGridCount(3)))
    allocate(filled(matGridCount(1), matGridCount(2), matGridCount(3)))
    filled = 0
    do i = 1, nmat
        do k = 1, 3
            off = (material(i,k) - matGridOrigin(k))/matGridSpacing(k)
            idx(k) = nint(off) + 1
            if (abs(off - nint(off)) > 1.0d-6 .or. idx(k) < 1 .or. idx(k) > matGridCount(k)) &
                call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
                'bMaterial.txt: a 3D material grid (n2mat=6) row has a coordinate that is not on the uniform grid.')
        enddo
        if (filled(idx(1),idx(2),idx(3)) /= 0) call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
            'bMaterial.txt: a 3D material grid (n2mat=6) cell is given twice.')
        if (material(i,4) <= 0.0d0 .or. material(i,5) <= 0.0d0 .or. material(i,6) <= 0.0d0) &
            call abortRun(ERR_CFG_MATERIAL_GRID_INVALID, &
            'bMaterial.txt: a 3D material grid (n2mat=6) row has vp, vs or rho <= 0.')
        filled(idx(1),idx(2),idx(3)) = 1
        matGrid3D(1:3, idx(1), idx(2), idx(3)) = material(i,4:6)
    enddo
    deallocate(filled)
end subroutine buildMaterialGrid3D

! #6 readstations --------------------------------------------------------
subroutine readstations1
! This subroutine is read information from bStations.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    call requireInputFile("bStations.txt")
    
    open(unit = 1006, file = 'bStations.txt', form = 'formatted', status = 'old')
        read(1006,*) totalNumOfOffSt
        read(1006,*) (nonfs(i), i = 1, ntotft)
        !write(*,*) 'totalNumOfOffSt,nonfs',totalNumOfOffSt, (nonfs(i), i = 1, ntotft), me
    close(1006)
end subroutine readstations1
! #7 readstations2 --------------------------------------------------------
subroutine readstations2
! This subroutine is read information from bStations.txt
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    call requireInputFile("bStations.txt")
    
    open(unit = 1006, file = 'bStations.txt', form = 'formatted', status = 'old')
        read(1006,*) 
        read(1006,*) 
        read(1006,*)
        do i = 1, ntotft
            do j = 1, nonfs(i)
                read(1006,*) xonfs(1,j,i), xonfs(2,j,i)
            enddo 
        enddo
        read(1006,*)
        do i = 1, totalNumOfOffSt
            read(1006,*) x4nds(1,i), x4nds(2,i), x4nds(3,i)
        enddo 
    close(1006)
    
    xonfs=xonfs*1000.0d0  !convert from km to m
    x4nds=x4nds*1000.0d0        
        
end subroutine readstations2

! #8 read_rough_geometry ------------------------------------------------
subroutine read_fault_rough_geometry
! This subroutine is read information from bFault_Rough_Geometry.txt
!
! The header/row checks below are defense in depth behind the case.setup-time
! validator (scripts/lib.py:validateFaultRoughGeometry, called from
! scripts/case.setup for every insertFaultType > 0). They are not redundant:
! these files get hand-edited, and a case is routinely set up on one machine
! and copied to an HPC where case.setup never runs again. Until v5.6.0 this
! reader took the header entirely on faith -- it read nnx/nnz, allocated
! rough_geo(3,nnx*nnz) and read exactly that many rows -- so a file that was
! short, or whose header disagreed with the mesh being built, either read past
! the end of the file or silently morphed every node onto the wrong surface.
! PROJECT_RULES.md rule 2: fail loudly instead.
!
! What the mesh actually requires (src/func_lib.f90:insertFaultInterface):
!   ixx = nint((x - rough_fx_min)/dx) + 1,  1 <= ixx <= nnx
!   izz = nint((z - rough_fz_min)/dz) + 1,  1 <= izz <= nnz
! i.e. the file's grid must BE the fault grid of the mesh being built: same
! corner, same spacing, same node counts. Note the indexing uses the MESH dx/dz
! while rough_fx_max is derived from the HEADER's dx, so a spacing mismatch is
! invisible at read time and stretches the surface at mesh time.
    use globalvar
    use errorCodes
    implicit none
    real (kind = dp) :: nnxTmp, nnzTmp, spanx, spanz, extraTmp
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j, ios, nnxExpect, nnzExpect

    call requireInputFile("bFault_Rough_Geometry.txt")

    open(unit = 1008, file = 'bFault_Rough_Geometry.txt', form = 'formatted', status = 'old')
        read(1008,*,iostat=ios) nnxTmp, nnzTmp
        if (ios /= 0) call stopRoughGeometry('header row 1 must read "nnx nnz" but is missing or non-numeric.')
        read(1008,*,iostat=ios) dxtmp, rough_fx_min, rough_fz_min
        if (ios /= 0) call stopRoughGeometry('header row 2 must read "dx fxmin fzmin" but is missing or non-numeric.')
    close(1008)
    nnx = nint(nnxTmp)
    nnz = nint(nnzTmp)

    if (abs(nnxTmp - nnx) > tol .or. abs(nnzTmp - nnz) > tol .or. nnx < 2 .or. nnz < 2) then
        if (me == masterProcsId) write(*,*) 'header nnx, nnz read as ', nnxTmp, nnzTmp
        call stopRoughGeometry('header nnx and nnz must be integers >= 2.')
    endif

    ! spacing: the header's dx must be the mesh's dx, because insertFaultInterface
    ! indexes the rough grid with the mesh dx/dz.
    if (abs(dxtmp - dx) > tol) then
        if (me == masterProcsId) write(*,*) 'file dx = ', dxtmp, ' but the mesh dx = ', dx
        call stopRoughGeometry('the rough-geometry file was sampled at a different cell size than this mesh.')
    endif

    ! corner: every index is counted from (rough_fx_min, rough_fz_min).
    if (abs(rough_fx_min - fltxyz(1,1,1)) > tol .or. abs(rough_fz_min - fltxyz(1,3,1)) > tol) then
        if (me == masterProcsId) then
            write(*,*) 'file corner (fxmin, fzmin) = ', rough_fx_min, rough_fz_min
            write(*,*) 'mesh fault corner          = ', fltxyz(1,1,1), fltxyz(1,3,1)
        endif
        call stopRoughGeometry('the rough-geometry file starts at a different fault corner than this mesh.')
    endif

    ! node counts: the fault must land on a whole number of cells, and the file
    ! must carry exactly that many nodes.
    spanx = (fltxyz(2,1,1) - fltxyz(1,1,1))/dx
    spanz = (fltxyz(2,3,1) - fltxyz(1,3,1))/dz
    if (abs(spanx - nint(spanx)) > tol .or. abs(spanz - nint(spanz)) > tol) then
        if (me == masterProcsId) write(*,*) 'fault extent / cell size = ', spanx, spanz, ' (both must be whole numbers)'
        call stopRoughGeometry('the fault edges do not land on rough-geometry nodes for this dx/dz.')
    endif
    nnxExpect = nint(spanx) + 1
    nnzExpect = nint(spanz) + 1
    if (nnx /= nnxExpect .or. nnz /= nnzExpect) then
        if (me == masterProcsId) then
            write(*,*) 'file   nnx, nnz = ', nnx, nnz
            write(*,*) 'mesh   nnx, nnz = ', nnxExpect, nnzExpect, ' from the fault extent and dx, dz = ', dx, dz
        endif
        call stopRoughGeometry('the rough-geometry grid is not the fault grid of the mesh being built.')
    endif

    rough_fx_max = (nnx - 1)*dxtmp + rough_fx_min

    allocate(rough_geo(3,nnx*nnz))

    open(unit = 1008, file = 'bFault_Rough_Geometry.txt', form = 'formatted', status = 'old')
        read(1008,*)
        read(1008,*)
        do i = 1, nnx*nnz
            read(1008,*,iostat=ios) rough_geo(1,i), rough_geo(2,i), rough_geo(3,i)
            if (ios /= 0) then
                if (me == masterProcsId) then
                    write(*,*) 'stopped at data row ', i, ' of ', nnx*nnz, ' (file row ', i+2, ')'
                endif
                call stopRoughGeometry('the file is short of, or has a malformed, "y dy/dx dy/dz" data row.')
            endif
            ! NaN/Inf here propagates straight into the mesh coordinates, where
            ! it is far harder to trace back to this file.
            do j = 1, 3
                if (rough_geo(j,i) /= rough_geo(j,i) .or. abs(rough_geo(j,i)) > 1.0d30) then
                    if (me == masterProcsId) write(*,*) 'non-finite value in column ', j, ' of data row ', i
                    call stopRoughGeometry('the rough-geometry file holds a NaN or Inf.')
                endif
            enddo
        enddo
        ! Per-cell fault-normal offset: the surface climbs |dy/dx|*dx from one
        ! x column to the next and |dy/dz|*dz from one z row to the next, while
        ! the nearest off-fault node layer is dy away. At an offset of one full
        ! cell a column's fault node passes its neighbour's first off-fault
        ! layer and insertFaultInterface tangles the elements. Checked on the
        ! RATIO, not the raw slope: a planar fault dipping at 45 degrees has
        ! dy/dz = cot(45) = 1 exactly while its true offset is dx*cos(45) < dy,
        ! so a raw-slope test would refuse every fault dipping 45 degrees or
        ! less (test.tpv10 dips 60 and carries |dy/dz| = 0.577).
        if (max(maxval(abs(rough_geo(2,:)))*dx/dy, &
                maxval(abs(rough_geo(3,:)))*dz/dy) >= 1.0d0) then
            if (me == masterProcsId) then
                write(*,*) 'max |dy/dx|*dx/dy = ', maxval(abs(rough_geo(2,:)))*dx/dy
                write(*,*) 'max |dy/dz|*dz/dy = ', maxval(abs(rough_geo(3,:)))*dz/dy
                write(*,*) 'dx, dy, dz = ', dx, dy, dz
            endif
            call stopRoughGeometry('the fault surface climbs a full fault-normal cell or more between adjacent fault nodes; the inserted elements would tangle.')
        endif

        ! A file with MORE rows than the header advertises was written for a
        ! different grid; the tail would be silently ignored.
        read(1008,*,iostat=ios) extraTmp
        if (ios == 0) then
            if (me == masterProcsId) write(*,*) 'expected exactly ', nnx*nnz, ' data rows after the 2 header rows'
            call stopRoughGeometry('the rough-geometry file has more data rows than its header declares.')
        endif
    close(1008)

end subroutine read_fault_rough_geometry

subroutine stopRoughGeometry(reason)
! Stop loudly on a bad bFault_Rough_Geometry.txt, naming the file, the reason
! and the fix (PROJECT_RULES.md rule 2). Every rank reads the same file and so
! reaches the same verdict; only the master prints, then all ranks stop.
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    character (len=*) :: reason
    integer (kind = 4) :: iMPIerr

    if (me == masterProcsId) then
        write(*,'(1X,A)') 'read_fault_rough_geometry: bFault_Rough_Geometry.txt is not usable for this mesh.'
        write(*,*) '  ', reason
        write(*,*) '  Regenerate it for this case (case.setup, which validates it, or'
        write(*,*) '  scripts/convertFaultGeometry for a supplied surface), or fix par so'
        write(*,*) '  the fault grid matches the file.'
    endif
    call MPI_Barrier(MPI_COMM_WORLD, iMPIerr)
    call abortRun(ERR_GEOM_ROUGH_INVALID, &
        'bFault_Rough_Geometry.txt is not usable for this mesh: '//trim(reason))
end subroutine stopRoughGeometry

subroutine requireInputFile(fileName)
! Stop with a clear message if a required input file is missing.
    use globalvar
    use errorCodes
    implicit none
    include 'mpif.h'
    character (len=*) :: fileName
    logical :: file_exists

    if (me == 0) then
        INQUIRE(FILE=fileName, EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            call abortRun(ERR_INPUT_FILE_MISSING, &
                trim(fileName)//' is required but missing. Run case.setup in this directory.')
        endif
    endif
end subroutine requireInputFile
