
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

    ccosphi = coheplas*dcos(atan(bulk))
    sinphi  = dsin(atan(bulk))
    nstep   = idnint(totalSimuTime/dt)
    rdampk  = rdampk*dt
    ! tv is no longer derived here -- it is read from bGlobal.txt (readglobal).
end subroutine readmaterial

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
        write(*,*) 'read_fault_rough_geometry: bFault_Rough_Geometry.txt is not usable for this mesh.'
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
