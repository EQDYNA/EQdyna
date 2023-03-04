
subroutine readglobal
! This subroutine is read information from bglobal.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i 
    
    
    if (me == 0) then 
        INQUIRE(FILE="bGlobal.txt", EXIST=file_exists)
        !write(*,*) 'Checking bGlobal.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bGlobal.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bGlobal.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bGlobal.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
    open(unit = 1001, file = 'bGlobal.txt', form = 'formatted', status = 'old')
        read(1001,*) mode
        read(1001,*) C_elastic
        read(1001,*) C_nuclea
        read(1001,*) C_degen
        read(1001,*) rough_fault
        read(1001,*) friclaw
        read(1001,*) ntotft
        read(1001,*) nucfault
        read(1001,*) TPV
        read(1001,*) output_plastic
        read(1001,*) outputGroundMotion
        read(1001,*) 
        read(1001,*) npx, npy, npz
        read(1001,*)
        read(1001,*) term
        read(1001,*) dt 
        read(1001,*)
        read(1001,*) nmat, n2mat
        read(1001,*) 
        read(1001,*) xsource, ysource, zsource
    close(1001)
end subroutine readglobal 
! #2 readmodelgeometry -------------------------------------------------
subroutine readmodelgeometry
! This subroutine is read information from bglobal.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    
    if (me == 0) then 
        INQUIRE(FILE="bModelGeometry.txt", EXIST=file_exists)
        !write(*,*) 'Checking bModel_Geometry.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bModelGeometry.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bModelGeometry.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bModelGeometry.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
    open(unit = 1002, file = 'bModelGeometry.txt', form = 'formatted', status = 'old')
        read(1002,*) xmin, xmax
        read(1002,*) ymin, ymax
        read(1002,*) zmin, zmax
        read(1002,*) 
        read(1002,*) dis4uniF, dis4uniB
        read(1002,*) rat
        read(1002,*) dx 
    close(1002)
end subroutine readmodelgeometry

! #3 readfaultgeometry
subroutine readfaultgeometry
! This subroutine is read information from bFault_Geometry.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4)::i
    
    if (me == 0) then 
        INQUIRE(FILE="bFaultGeometry.txt", EXIST=file_exists)
        !write(*,*) 'Checking bFault_Geometry.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bFaultGeometry.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bFaultGeometry.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bFaultGeometry.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
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
        fltxyz(2,4,i)=fdip*pi/180.0d0
    enddo
    
end subroutine readfaultgeometry

! #4 readmaterial --------------------------------------------------------
subroutine readmaterial
! This subroutine is read information from bMaterial.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    if (me == 0) then 
        INQUIRE(FILE="bMaterial.txt", EXIST=file_exists)
        !write(*,*) 'Checking bMaterial.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bMaterial.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bMaterial.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bMaterial.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
    open(unit = 1004, file = 'bMaterial.txt', form = 'formatted', status = 'old')
        do i = 1, nmat
            read(1004,*) (material(i,j), j = 1, n2mat)
        enddo 
    close(1004)

    ccosphi=coheplas*dcos(atan(bulk))
    sinphi=dsin(atan(bulk))
    nstep=idnint(term/dt)
    rdampk=rdampk*dt    
    tv = 2.0d0*dx/3464.0d0
end subroutine readmaterial

! ! #5 readfric --------------------------------------------------------
! 20221128. Not used. Fric params are input from netCDF read on the whole map.
! subroutine readfric
! ! This subroutine is read information from bMaterial.txt
    ! use globalvar
    ! implicit none
    ! include 'mpif.h'

    ! logical::file_exists
    ! integer(kind=4):: i, j 
    
    ! if (me == 0) then 
        ! INQUIRE(FILE="bFric.txt", EXIST=file_exists)
        ! !write(*,*) 'Checking bFric.txt by the master procs', me
        ! if (file_exists .eqv. .FALSE.) then
            ! write(*,*) 'bFric.txt is required but missing ...'
            !     
        ! endif 
    ! endif 
    ! if (me == 0) then 
        ! INQUIRE(FILE="bFric.txt", EXIST=file_exists)
        ! if (file_exists .eqv. .FALSE.) then
            ! write(*,*) 'bFric.txt is still missing, so exiting EQdyna'
            ! stop
        ! endif 
    ! endif     
    
    ! open(unit = 1005, file = 'bFric.txt', form = 'formatted', status = 'old')
    ! if (friclaw == 1) then ! slip-weakening
        ! read(1005,*) fric_sw_fs
        ! read(1005,*) fric_sw_fd
        ! read(1005,*) fric_sw_D0
    ! endif 
    ! if (friclaw == 3) then ! aging law
        ! do i = 1, 4
            ! read(1005,*)
        ! enddo
        ! read(1005,*) fric_rsf_a
        ! read(1005,*) fric_rsf_deltaa0
        ! read(1005,*) fric_rsf_b
        ! read(1005,*) fric_rsf_Dc
        ! read(1005,*) fric_rsf_r0
        ! read(1005,*) fric_rsf_v0
        ! read(1005,*)
        ! read(1005,*) fric_rsf_vinix    
        ! read(1005,*) fric_rsf_viniz    
    ! endif     
    ! if (friclaw == 4) then ! strong-rate weakening
        ! do i = 1, 4
            ! read(1005,*)
        ! enddo
        ! read(1005,*) fric_rsf_a
        ! read(1005,*) fric_rsf_deltaa0
        ! read(1005,*) fric_rsf_b
        ! read(1005,*) fric_rsf_Dc
        ! read(1005,*) fric_rsf_r0
        ! read(1005,*) fric_rsf_v0
        ! read(1005,*)
        ! read(1005,*) fric_rsf_vinix    
        ! read(1005,*) fric_rsf_viniz    
        ! read(1005,*)
        ! read(1005,*) fric_rsf_fw
        ! read(1005,*) fric_rsf_vw
        ! read(1005,*) fric_rsf_deltavw0    
    ! endif         
    ! if (friclaw == 5) then ! strong-rate weakening + termop
        ! do i = 1, 4
            ! read(1005,*)
        ! enddo
        ! read(1005,*) fric_rsf_a
        ! read(1005,*) fric_rsf_deltaa0
        ! read(1005,*) fric_rsf_b
        ! read(1005,*) fric_rsf_Dc
        ! read(1005,*) fric_rsf_r0
        ! read(1005,*) fric_rsf_v0
        ! read(1005,*)
        ! read(1005,*) fric_rsf_vinix    
        ! read(1005,*) fric_rsf_viniz    
        ! read(1005,*)
        ! read(1005,*) fric_rsf_fw
        ! read(1005,*) fric_rsf_vw
        ! read(1005,*) fric_rsf_deltavw0        
        ! read(1005,*)
        ! read(1005,*) fric_tp_a_th
        ! read(1005,*) fric_tp_rouc
        ! read(1005,*) fric_tp_lambda
        ! read(1005,*) fric_tp_h
        ! read(1005,*) fric_tp_a_hy
        ! read(1005,*) fric_tp_deltaa_hy0
        ! read(1005,*) fric_ww
        ! read(1005,*) fric_w
        ! read(1005,*) 
        ! read(1005,*) fric_ini_sliprate
        ! read(1005,*) fric_tp_Tini
        ! read(1005,*) fric_tp_pini
        ! ! read(1005,*) dxtp
        ! ! read(1005,*) tpw
    ! endif     
    ! close(1005)
! end subroutine readfric

! #6 readstations --------------------------------------------------------
subroutine readstations1
! This subroutine is read information from bStations.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    if (me == 0) then 
        INQUIRE(FILE="bStations.txt", EXIST=file_exists)
        !write(*,*) 'Checking bStations.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bStations.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bStations.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bStations.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
    open(unit = 1006, file = 'bStations.txt', form = 'formatted', status = 'old')
        read(1006,*) n4nds
        read(1006,*) (nonfs(i), i = 1, ntotft)
        !write(*,*) 'n4nds,nonfs',n4nds, (nonfs(i), i = 1, ntotft), me
    close(1006)
end subroutine readstations1
! #7 readstations2 --------------------------------------------------------
subroutine readstations2
! This subroutine is read information from bStations.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j 
    
    if (me == 0) then 
        INQUIRE(FILE="bStations.txt", EXIST=file_exists)
        !write(*,*) 'Checking bStations.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bStations.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bStations.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bStations.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
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
        do i = 1, n4nds
            read(1006,*) x4nds(1,i), x4nds(2,i), x4nds(3,i)
        enddo 
    close(1006)
    
    xonfs=xonfs*1000.0d0  !convert from km to m
    x4nds=x4nds*1000.0d0        
        
end subroutine readstations2

! #8 read_rough_geometry ------------------------------------------------
subroutine read_fault_rough_geometry
! This subroutine is read information from bFault_Rough_Geometry.txt
    use globalvar
    implicit none
    include 'mpif.h'

    logical::file_exists
    integer(kind=4):: i, j
    
    if (me == 0) then 
        INQUIRE(FILE="bFault_Rough_Geometry.txt", EXIST=file_exists)
        !write(*,*) 'Checking bFault_Rough_Geometry.txt by the master procs', me
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bFault_Rough_Geometry.txt is required but missing ...'
                
        endif 
    endif 
    if (me == 0) then 
        INQUIRE(FILE="bFault_Rough_Geometry.txt", EXIST=file_exists)
        if (file_exists .eqv. .FALSE.) then
            write(*,*) 'bFault_Rough_Geometry.txt is still missing, so exiting EQdyna'
            stop
        endif 
    endif     
    
    open(unit = 1008, file = 'bFault_Rough_Geometry.txt', form = 'formatted', status = 'old')
        read(1008,*) nnx, nnz
        read(1008,*) dxtmp 
    close(1008)
    
    allocate(rough_geo(3,nnx*nnz))
    
    open(unit = 1008, file = 'bFault_Rough_Geometry.txt', form = 'formatted', status = 'old')
        read(1008,*)
        read(1008,*)
        do i = 1, nnx*nnz
            read(1008,*) rough_geo(1,i), rough_geo(2,i), rough_geo(3,i)
        enddo
    close(1008)    
        
end subroutine read_fault_rough_geometry