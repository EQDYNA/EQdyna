! netcdf_io reads and writes out arrays into netcdf format.
! A Fortran example code is shown here https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_wr.f90
! A companion exmaple code to read the data written by sfc_pres_temp_wr.f90 could be found via 
! https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_rd.f90.
! Note that if a 6 X 12 data on a lat-lon grid is to be created, the variable array should be in 
! this format - on_fault_vars(lon, lat). 
! In our cases, lat = z/dip and lon = x/strike. So, on_fault_vars should be on_fault_vars(nx,nz)

! NOTE that FORTRAN 90 NetCDF API requires variable dimensions in a reverse order from Python. 

! subroutines contained in this file includes: 
! - #1: netcdf_read_on_fault_eqdyna
! - #2: netcdf_read_on_fault_eqdyna_restart
! - #A1: check

! Subroutine #1.
! netcdf_read_on_fault reads in on-fault quantities from netcdf files created by case.setup.
subroutine netcdf_read_on_fault_eqdyna
    use netcdf
    use globalvar
    implicit none 
    character (len = 50 ) :: infile
    integer (kind = 4) :: ncid,  var_id(40), i, j, nvar, ii, jj, fnx, fnz, ift
    real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars
    real (kind = dp)   :: xcord, zcord
    
    fnx  = (fxmax(1) - fxmin(1))/dx+1
    fnz  = (fzmax(1) - fzmin(1))/dx+1
    
    infile = "on_fault_vars_input.nc"
    
    nvar = 24
    allocate(on_fault_vars(fnx,fnz,nvar))
    
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file. 
    call check( nf90_open(infile, NF90_NOWRITE, ncid))
    
    ! Get the varid of the data variables, based on their names.
    call check( nf90_inq_varid(ncid, "sw_fs",    var_id(1)))
    call check( nf90_inq_varid(ncid, "sw_fd",    var_id(2)))
    call check( nf90_inq_varid(ncid, "sw_D0",    var_id(3)))
    call check( nf90_inq_varid(ncid, "rsf_a",    var_id(4)))
    call check( nf90_inq_varid(ncid, "rsf_b",    var_id(5)))
    call check( nf90_inq_varid(ncid, "rsf_Dc",   var_id(6)))
    call check( nf90_inq_varid(ncid, "rsf_v0",   var_id(7)))
    call check( nf90_inq_varid(ncid, "rsf_r0",   var_id(8)))
    call check( nf90_inq_varid(ncid, "rsf_fw",   var_id(9)))
    call check( nf90_inq_varid(ncid, "rsf_vw",   var_id(10)))
    call check( nf90_inq_varid(ncid, "tp_a_hy",  var_id(11)))
    call check( nf90_inq_varid(ncid, "tp_a_th",  var_id(12)))
    call check( nf90_inq_varid(ncid, "tp_rouc",  var_id(13)))
    call check( nf90_inq_varid(ncid, "tp_lambda",var_id(14)))
    call check( nf90_inq_varid(ncid, "tp_h",     var_id(15)))
    call check( nf90_inq_varid(ncid, "tp_Tini",  var_id(16)))
    call check( nf90_inq_varid(ncid, "tp_pini",  var_id(17)))
    call check( nf90_inq_varid(ncid, "init_slip_rate",     var_id(18)))
    call check( nf90_inq_varid(ncid, "init_strike_shear",  var_id(19)))
    call check( nf90_inq_varid(ncid, "init_normal_stress", var_id(20)))
    call check( nf90_inq_varid(ncid, "init_state",         var_id(21)))
    call check( nf90_inq_varid(ncid, "tw_t0",          var_id(22)))
    call check( nf90_inq_varid(ncid, "cohesion",       var_id(23)))
    call check( nf90_inq_varid(ncid, "init_dip_shear", var_id(24)))
    ! Read the data
    do i = 1, nvar
        call check( nf90_get_var(ncid, var_id(i), on_fault_vars(:,:,i)))
    enddo        

    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            xcord          = x(1, nsmp(1,i,ift))
            zcord          = x(3, nsmp(1,i,ift))
            ii             = (xcord - fxmin(ift))/dx + 1
            jj             = (zcord - fzmin(ift))/dx + 1
            fric(1,  i, 1) = on_fault_vars(ii,jj,1) ! sw_fs
            fric(2,  i, 1) = on_fault_vars(ii,jj,2) ! sw_fd
            fric(3,  i, 1) = on_fault_vars(ii,jj,3) ! sw_D0
            fric(9,  i, 1) = on_fault_vars(ii,jj,4) ! rsf_a
            fric(10, i, 1) = on_fault_vars(ii,jj,5) ! rsf_b
            fric(11, i, 1) = on_fault_vars(ii,jj,6) ! rsf_Dc
            fric(12, i, 1) = on_fault_vars(ii,jj,7) ! rsf_v0
            fric(13, i, 1) = on_fault_vars(ii,jj,8) ! rsf_r0
            fric(14, i, 1) = on_fault_vars(ii,jj,9) ! rsf_fw
            fric(15, i, 1) = on_fault_vars(ii,jj,10)! rsf_vw
            fric(16, i, 1) = on_fault_vars(ii,jj,11)! tp_a_hy
            fric(17, i, 1) = on_fault_vars(ii,jj,12)! tp_a_th
            fric(18, i, 1) = on_fault_vars(ii,jj,13)! tp_rouc
            fric(19, i, 1) = on_fault_vars(ii,jj,14)! tp_lambda
            fric(40, i, 1) = on_fault_vars(ii,jj,15)! tp_h
            fric(41, i, 1) = on_fault_vars(ii,jj,16)! tp_Tini 
            fric(42, i, 1) = on_fault_vars(ii,jj,17)! tp_pini 
            fric(46, i, 1) = on_fault_vars(ii,jj,18)! creeping slip rate, lower bound
            fric(8,  i, 1) = on_fault_vars(ii,jj,19)! init_strike_shear
            fric(7,  i, 1) = on_fault_vars(ii,jj,20)! init_norm
            fric(20, i, 1) = on_fault_vars(ii,jj,21)! init_state variable
            fric(47, i, 1) = fric(46, i, 1)         ! peak slip rate
            fric(25, i, 1) = 0.0d0                  ! vini_norm 
            fric(26, i, 1) = fric(46, i, 1)         ! vinix
            fric(27, i, 1) = 0.0d0                  ! viniz
            fric(5,  i, 1) = on_fault_vars(ii,jj,22)! tw_t0
            fric(4,  i, 1) = on_fault_vars(ii,jj,23)! cohesion
            fric(49, i, 1) = on_fault_vars(ii,jj,24)! init_dip_shear
            
            fric(23, i, 1) = abs(fric(7, i, 1))     ! initialize theta_pc as abs(normal stress)
        enddo 
    enddo
    
    ! Close the file, freeing all resources.
    call check( nf90_close(ncid))

end subroutine netcdf_read_on_fault_eqdyna

! Subroutine #2.
! netcdf_read_on_fault_eqdyna_restart reads in additional on-fault quantities from restart *.r.nc netcdf files created by previous cycles.
subroutine netcdf_read_on_fault_eqdyna_restart
    use netcdf
    use globalvar
    implicit none 
    character (len = 50 ) :: infile
    integer (kind = 4) :: ncid,  var_id(20), i, j, nvar, fnx, fnz, ii, jj, ift
    real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars
    real (kind = dp)   :: xcord, zcord
    
    infile = "fault.r.nc"
    
    fnx  = (fxmax(1) - fxmin(1))/dx+1
    fnz  = (fzmax(1) - fzmin(1))/dx+1
    
    ! Read in initial conditions from restart files fault.r.nc spun off by EQquasi. 
    nvar = 12
    allocate(on_fault_vars(fnx,fnz,nvar))    
    
    ! NOTE. the array structure is different than loading python generated nc file.
    ! here we follow the structure of subroutine netcdf_write_on_fault.
    ! on_fault_vars is now nxt by nzt!!!
    ! allocate(on_fault_vars(nxt,nzt,nvar))
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file. 
    call check( nf90_open(infile, NF90_NOWRITE, ncid))

    ! Get the varid of the data variables, based on their names.
    ! 'shear_strike', 'shear_dip', 'effective_normal', 'slip_rate' , 'state_variable', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs'
    call check( nf90_inq_varid(ncid, "shear_strike",     var_id(1)))
    call check( nf90_inq_varid(ncid, "shear_dip",        var_id(2)))
    call check( nf90_inq_varid(ncid, "effective_normal", var_id(3)))
    call check( nf90_inq_varid(ncid, "slip_rate",        var_id(4)))
    call check( nf90_inq_varid(ncid, "state_variable",   var_id(5)))
    call check( nf90_inq_varid(ncid, "state_normal",     var_id(6)))
    call check( nf90_inq_varid(ncid, "vxm",              var_id(7)))
    call check( nf90_inq_varid(ncid, "vym",              var_id(8)))
    call check( nf90_inq_varid(ncid, "vzm",              var_id(9)))
    call check( nf90_inq_varid(ncid, "vxs",              var_id(10)))
    call check( nf90_inq_varid(ncid, "vys",              var_id(11)))
    call check( nf90_inq_varid(ncid, "vzs",              var_id(12)))
    ! Read the data
    do i = 1, nvar
        call check( nf90_get_var(ncid, var_id(i), on_fault_vars(:,:,i)))
    enddo         
    
    do ift = 1, ntotft
        do i = 1, nftnd(ift)
            xcord            = x(1, nsmp(1,i,ift))
            zcord            = x(3, nsmp(1,i,ift))
            ii               = (xcord - fxmin(ift))/dx + 1
            jj               = (zcord - fzmin(ift))/dx + 1
            fric(8,  i, ift) = on_fault_vars(ii,jj, 1) ! tstk0
            fric(49, i, ift) = on_fault_vars(ii,jj, 2) ! tdip0
            fric(7,  i, ift) = on_fault_vars(ii,jj, 3) ! tnorm0
            fric(47, i, ift) = on_fault_vars(ii,jj, 4) ! sliprate
            fric(20, i, ift) = on_fault_vars(ii,jj, 5) ! state
            fric(23, i, ift) = on_fault_vars(ii,jj, 6) ! state_normal
            fric(31, i, ift) = on_fault_vars(ii,jj, 7) ! vxm
            fric(32, i, ift) = on_fault_vars(ii,jj, 8) ! vym
            fric(33, i, ift) = on_fault_vars(ii,jj, 9) ! vzm
            fric(34, i, ift) = on_fault_vars(ii,jj, 10)! vxs
            fric(35, i, ift) = on_fault_vars(ii,jj, 11)! vys
            fric(36, i, ift) = on_fault_vars(ii,jj, 12)! vzs
            !fric(23, (i-1)*nzt+j, 1) = abs(fric(7, (i-1)*nzt+j, 1))! initialize theta_pc as abs(normal stress)
        enddo 
    enddo 
    ! Close the file, freeing all resources.
    call check( nf90_close(ncid))
    
    deallocate(on_fault_vars)

end subroutine netcdf_read_on_fault_eqdyna_restart


! #A1
subroutine check(status)
    use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, nf90_strerror(status)
      stop "Stopped"
    end if
end subroutine check  
