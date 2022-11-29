! netcdf_io reads and writes out arrays into netcdf format.
! A Fortran example code is shown here https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_wr.f90

! subroutines contained in this file includes: 
! - netcdf_write
! - netcdf_write_on_fault

! subroutine netcdf_write(outfile, outtype)
	! use netcdf
	! use globalvar
	! implicit none 
	! character (len = 50 ) :: outfile, outtype, lat_name, lon_name, lat_units, lon_units, var_1_name, var_1_unit, UNITS
	! integer (kind = 4) :: ncid, lat_dimid, lon_dimid, lat_varid, lon_varid, var_1_id, ilat, ilon, i, j, nlat, nlon
	! integer (kind = 4) :: dimids(2),  nsmp_tmp(2,nftmx)
	! integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index

	! UNITS = 'units'
	! if (outtype == 'disp') then 
		! lat_name = 'node_id'
		! lon_name = 'disp_dim_id'
		! lat_units = 'unit'
		! lon_units = 'unit'
		! var_1_name = 'disp'
		! var_1_unit = 'meters'
		! nlat = numnp
		! nlon = ndof
	! elseif (outtype == 'coor') then 
		! lat_name = 'node_id'
		! lon_name = 'coor_dim_id'
		! lat_units = 'unit'
		! lon_units = 'unit'
		! var_1_name = 'coordinates'
		! var_1_unit = 'meters'
		! nlat = numnp
		! nlon = ndof
	! elseif (outtype == 'ien') then 
		! lat_name = 'element_id'
		! lon_name = 'node_index_id'
		! lat_units = 'unit'
		! lon_units = 'unit'
		! var_1_name = 'node_id'
		! var_1_unit = 'unit'
		! nlat = numel
		! nlon = nen
	! elseif (outtype == 'nsmp') then
		! lat_name = 'fault_node_id'
		! lon_name = 'slave_master_id'
		! lat_units = 'unit'
		! lon_units = 'unit'
		! var_1_name = 'node_id'
		! var_1_unit = 'unit'
		! nlat = nftmx
		! nlon = 2
		! do i = 1, nftmx
			! do j = 1, 2
				! nsmp_tmp(j,i) = nsmp(j,i,1)
			! enddo 
		! enddo 
	! endif 
	! allocate(lat_index(nlat))
	! allocate(lon_index(nlon))
	! lat_index = (/ (i, i = 1, nlat) /)
	! lon_index = (/ (i, i = 1, nlon) /)
	! ! Create the netCDF file.
	! call check(nf90_create(outfile, NF90_CLOBBER, ncid))
	
	! ! Define the dimensions.
	! call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
	! call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
	
	! ! Define coordiante variables. They will hold the coordinate 
	! ! information, that is, the latitudes (y), and longitudes (x). A varid is 
	! ! returned for each.
	! call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
	! call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
	
	! ! Assign units attributes to coordinate var data. This attaches a
	! ! text attribute to each of the coordinate variables, containing the
	! ! units.
	! call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
	! call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )	

	! ! Define the netcdf variables. The dimids array is used to pass the 
	! ! dimids of the dimensions of the netCDF variables.
	! dimids = (/ lon_dimid, lat_dimid /)
	! if (outtype == 'disp' .or. outtype == 'coor') then 
		! call check(nf90_def_var(ncid, var_1_name, NF90_REAL, dimids, var_1_id)) 
	! elseif (outtype == 'nsmp' .or. outtype == 'ien') then 
		! call check(nf90_def_var(ncid, var_1_name, NF90_INT, dimids, var_1_id)) 
	! endif 
	! ! Assign units attributes to the pressure and temperature netCDF
	! ! variables.
	! call check(nf90_put_att(ncid, var_1_id, UNITS, var_1_unit))
	
	! ! End definitions.
	! call check(nf90_enddef(ncid))
	! ! Write data.
	! ! Write the coordinate variable data. This will put the x, and y 
	! ! of our data grid into the netCDF file.
	! call check(nf90_put_var(ncid, lat_varid, lat_index))
	! call check(nf90_put_var(ncid, lon_varid, lon_index))
	
	! ! Write the data. This will write our displacement fields "cons" which is defined 
	! ! globally. Its dimension is 3 by numnp (row by column).
	! if (outtype == 'disp') then 
		! call check(nf90_put_var(ncid, var_1_id, cons))
	! elseif (outtype == 'coor') then 
		! call check(nf90_put_var(ncid, var_1_id, x))
	! elseif (outtype == 'ien') then 
		! call check(nf90_put_var(ncid, var_1_id, ien))
	! elseif (outtype == 'nsmp') then 
		! call check(nf90_put_var(ncid, var_1_id, nsmp_tmp))
	! endif 
	! ! Close the file.
	! call check(nf90_close(ncid))
	! !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

! end subroutine netcdf_write

! ! netcdf_write_on_fault writes on-fault quantities to netcdf files.
! subroutine netcdf_write_on_fault(outfile)
	! use netcdf
	! use globalvar
	! implicit none 
	! character (len = 50 ) :: outfile, outtype, lat_name, lon_name, lat_units, lon_units, UNITS
	! character (len = 50), allocatable, dimension(:) :: var_name, var_unit
	! integer (kind = 4) :: ncid, lat_dimid, lon_dimid, lat_varid, lon_varid, var_id(20), ilat, ilon, i, j, nlat, nlon, nvar
	! integer (kind = 4) :: dimids(2)
	! integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index
	! real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars

	! nvar = 12 ! nvar variables.
	! allocate(var_name(nvar), var_unit(nvar), on_fault_vars(nxt, nzt, nvar))
	! on_fault_vars = 0.0d0
	
	! UNITS = 'units'
	! lat_name = 'nid_dip'
	! lon_name = 'nid_strike'
	! lat_units = 'unit'
	! lon_units = 'unit'
	! var_name = [ character(len=20) :: 'shear_strike', 'shear_dip', 'effective_normal', 'slip_rate' , 'state_variable', 'state_normal', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs']
	! var_unit = [ character(len=20) :: 'Pa'          , 'Pa'       , 'Pa'              , 'm/s'       , 'unit'          , 'Pa'          , 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s']
	! nlat = nzt ! Total nodes along dip. 
	! nlon = nxt ! Total nodes along strike.  
	
	! allocate(lat_index(nlat))
	! allocate(lon_index(nlon))
	
	! lat_index = (/ (i, i = 1, nlat) /)
	! lon_index = (/ (i, i = 1, nlon) /)
	
	! do i = 1, nxt
		! do j = 1, nzt
			! on_fault_vars(i,j, 1) = fric(28, (i-1)*nzt+j, 1) ! tstk0
			! on_fault_vars(i,j, 2) = fric(29, (i-1)*nzt+j, 1) ! tdip0
			! on_fault_vars(i,j, 3) = fric(30, (i-1)*nzt+j, 1) ! tnorm0
			! on_fault_vars(i,j, 4) = fric(26, (i-1)*nzt+j, 1) ! sliprate
			! on_fault_vars(i,j, 5) = fric(20, (i-1)*nzt+j, 1) ! state
			! on_fault_vars(i,j, 6) = fric(23, (i-1)*nzt+j, 1) ! state variable for normal stress, theta_pc
			! on_fault_vars(i,j, 7) = fric(31, (i-1)*nzt+j, 1) ! vxm
			! on_fault_vars(i,j, 8) = fric(32, (i-1)*nzt+j, 1) ! vym
			! on_fault_vars(i,j, 9) = fric(33, (i-1)*nzt+j, 1) ! vzm
			! on_fault_vars(i,j, 10) = fric(34, (i-1)*nzt+j, 1) ! vxs
			! on_fault_vars(i,j, 11) = fric(35, (i-1)*nzt+j, 1) ! vys
			! on_fault_vars(i,j, 12) = fric(36, (i-1)*nzt+j, 1) ! vzs
		! enddo 
	! enddo 
	! ! Create the netCDF file.
	! call check(nf90_create(outfile, NF90_CLOBBER, ncid))
	
	! ! Define the dimensions.
	! call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
	! call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
	
	! ! Define coordiante variables. They will hold the coordinate 
	! ! information, that is, the latitudes (y), and longitudes (x). A varid is 
	! ! returned for each.
	! call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
	! call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
	
	! ! Assign units attributes to coordinate var data. This attaches a
	! ! text attribute to each of the coordinate variables, containing the
	! ! units.
	! call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
	! call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )	

	! ! Define the netcdf variables. The dimids array is used to pass the 
	! ! dimids of the dimensions of the netCDF variables.
	! dimids = (/ lon_dimid, lat_dimid /)
	
	! do i = 1, nvar
		! call check(nf90_def_var(ncid, var_name(i), NF90_REAL, dimids, var_id(i))) 
		! ! Assign units attributes to the pressure and temperature netCDF
		! ! variables.
		! call check(nf90_put_att(ncid, var_id(i), UNITS, var_unit(i)))
	! enddo 

	! ! End definitions.
	! call check(nf90_enddef(ncid))
	
	! ! Write data.
	! ! Write the coordinate variable data. This will put the x, and y 
	! ! of our data grid into the netCDF file.
	! call check(nf90_put_var(ncid, lat_varid, lat_index))
	! call check(nf90_put_var(ncid, lon_varid, lon_index))
	
	! ! Write the data. This will write our displacement fields "cons" which is defined 
	! ! globally. Its dimension is 3 by numnp (row by column).
	! do i = 1,nvar
		! call check(nf90_put_var(ncid, var_id(i), on_fault_vars(:,:,i)))
	! enddo
	
	! ! Close the file.
	! call check(nf90_close(ncid))
	! !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

! end subroutine netcdf_write_on_fault

! subroutine netcdf_write_roughness(outfile)
	! use netcdf
	! use globalvar
	! implicit none 
	! character (len = 50 ) :: outfile, outtype, lat_name, lon_name, lat_units, lon_units, UNITS
	! character (len = 50), allocatable, dimension(:) :: var_name, var_unit
	! integer (kind = 4) :: ncid, lat_dimid, lon_dimid, lat_varid, lon_varid, var_id(20), ilat, ilon, i, j, nlat, nlon, nvar
	! integer (kind = 4) :: dimids(2)
	! integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index
	! real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars

	! nvar = 3 ! nvar variables.
	! allocate(var_name(nvar), var_unit(nvar), on_fault_vars(nnx, nnz, nvar))
	! on_fault_vars = 0.0d0
	
	! UNITS = 'units'
	! lat_name = 'nid_dip'
	! lon_name = 'nid_strike'
	! lat_units = 'unit'
	! lon_units = 'unit'
	! var_name = [ character(len=20) :: 'peak', 'pypx', 'pypz' ]
	! var_unit = [ character(len=20) :: 'm'          , 'unit'       , 'unit']
	! nlat = nnz ! Total nodes along dip. 
	! nlon = nnx ! Total nodes along strike.  
	
	! allocate(lat_index(nlat))
	! allocate(lon_index(nlon))
	
	! lat_index = (/ (i, i = 1, nlat) /)
	! lon_index = (/ (i, i = 1, nlon) /)
	
	! do i = 1, nnx
		! do j = 1, nnz
			! on_fault_vars(i,j, 1) = rough_geo(1, (i-1)*nnz+j) ! tstk0
			! on_fault_vars(i,j, 2) = rough_geo(2, (i-1)*nnz+j) ! tdip0
			! on_fault_vars(i,j, 3) = rough_geo(3, (i-1)*nnz+j) ! tnorm0
		! enddo 
	! enddo 
	! ! Create the netCDF file.
	! call check(nf90_create(outfile, NF90_CLOBBER, ncid))
	
	! ! Define the dimensions.
	! call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
	! call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
	
	! ! Define coordiante variables. They will hold the coordinate 
	! ! information, that is, the latitudes (y), and longitudes (x). A varid is 
	! ! returned for each.
	! call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
	! call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
	
	! ! Assign units attributes to coordinate var data. This attaches a
	! ! text attribute to each of the coordinate variables, containing the
	! ! units.
	! call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
	! call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )	

	! ! Define the netcdf variables. The dimids array is used to pass the 
	! ! dimids of the dimensions of the netCDF variables.
	! dimids = (/ lon_dimid, lat_dimid /)
	
	! do i = 1, nvar
		! call check(nf90_def_var(ncid, var_name(i), NF90_REAL, dimids, var_id(i))) 
		! ! Assign units attributes to the pressure and temperature netCDF
		! ! variables.
		! call check(nf90_put_att(ncid, var_id(i), UNITS, var_unit(i)))
	! enddo 

	! ! End definitions.
	! call check(nf90_enddef(ncid))
	
	! ! Write data.
	! ! Write the coordinate variable data. This will put the x, and y 
	! ! of our data grid into the netCDF file.
	! call check(nf90_put_var(ncid, lat_varid, lat_index))
	! call check(nf90_put_var(ncid, lon_varid, lon_index))
	
	! ! Write the data. This will write our displacement fields "cons" which is defined 
	! ! globally. Its dimension is 3 by numnp (row by column).
	! do i = 1,nvar
		! call check(nf90_put_var(ncid, var_id(i), on_fault_vars(:,:,i)))
	! enddo
	
	! ! Close the file.
	! call check(nf90_close(ncid))
	! !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

! end subroutine netcdf_write_roughness

! Subroutine #4.
! netcdf_read_on_fault reads in on-fault quantities from netcdf files created by case.setup.
subroutine netcdf_read_on_fault_eqdyna(infile)
	use netcdf
	use globalvar
	implicit none 
	character (len = 50 ) :: infile
	integer (kind = 4) :: ncid,  var_id(40), i, j, nvar, ii, jj, fnx, fnz, ift
	real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars
	real (kind = dp)   :: xcord, zcord
	
	fnx  = (fxmax(1) - fxmin(1))/dx+1
	fnz  = (fzmax(1) - fzmin(1))/dx+1
	
	nvar = 21
	allocate(on_fault_vars(fnz,fnx,nvar))
	
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
	call check( nf90_inq_varid(ncid, "init_shear_stress",  var_id(19)))
	call check( nf90_inq_varid(ncid, "init_normal_stress", var_id(20)))
	call check( nf90_inq_varid(ncid, "init_state",         var_id(21)))
	
	! Read the data
	do i = 1, nvar
		call check( nf90_get_var(ncid, var_id(i), on_fault_vars(:,:,i)))
	enddo		

	do ift = 1, ntotft
		do i = 1, nftnd(ift)
			xcord = x(1, nsmp(1,i,ift))
			zcord = x(3, nsmp(1,i,ift))
			ii    = (xcord - fxmin(ift))/dx + 1
			jj    = (zcord - fzmin(ift))/dx + 1
			fric(1,  i, 1) = on_fault_vars(jj,ii,1)! sw_fs
			fric(2,  i, 1) = on_fault_vars(jj,ii,2)! sw_fd
			fric(3,  i, 1) = on_fault_vars(jj,ii,3)! sw_D0
			fric(9,  i, 1) = on_fault_vars(jj,ii,4)! rsf_a
			fric(10, i, 1) = on_fault_vars(jj,ii,5)! rsf_b
			fric(11, i, 1) = on_fault_vars(jj,ii,6)! rsf_Dc
			fric(12, i, 1) = on_fault_vars(jj,ii,7)! rsf_v0
			fric(13, i, 1) = on_fault_vars(jj,ii,8)! rsf_r0
			fric(14, i, 1) = on_fault_vars(jj,ii,9) ! rsf_fw
			fric(15, i, 1) = on_fault_vars(jj,ii,10)! rsf_vw
			fric(16, i, 1) = on_fault_vars(jj,ii,11)! tp_a_hy
			fric(17, i, 1) = on_fault_vars(jj,ii,12)! tp_a_th
			fric(18, i, 1) = on_fault_vars(jj,ii,13)! tp_rouc
			fric(19, i, 1) = on_fault_vars(jj,ii,14)! tp_lambda
			fric(40, i, 1) = on_fault_vars(jj,ii,15)! tp_h
			fric(41, i, 1) = on_fault_vars(jj,ii,16)! tp_Tini 
			fric(42, i, 1) = on_fault_vars(jj,ii,17)! tp_pini 
			fric(46, i, 1) = on_fault_vars(jj,ii,18)! init_slip rate
			fric(8,  i, 1) = on_fault_vars(jj,ii,19)! init_shear
			fric(7,  i, 1) = on_fault_vars(jj,ii,20)! init_norm
			fric(20, i, 1) = on_fault_vars(jj,ii,21)! init_state variable
			fric(47, i, 1) = fric(46, i, 1)! peak slip rate
			fric(25, i, 1) = 0.0d0! vini_norm 
			fric(26, i, 1) = fric(46, i, 1) ! vinix
			fric(27, i, 1) = 0.0d0! viniz
			
			fric(23, i, 1) = abs(fric(7, i, 1))! initialize theta_pc as abs(normal stress)
		enddo 
	enddo
	
	! Close the file, freeing all resources.
	call check( nf90_close(ncid))

end subroutine netcdf_read_on_fault_eqdyna

! ! Subroutine #5.
! ! netcdf_read_on_fault_restart reads in additional on-fault quantities from netcdf files created by previous cycles.
! subroutine netcdf_read_on_fault_restart(infile1, infile2)
	! use netcdf
	! use globalvar
	! implicit none 
	! character (len = 50 ) :: infile1, infile2
	! integer (kind = 4) :: ncid,  var_id(20), i, j, nvar
	! real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars
	
	! ! Read in 5 variables a, b, Dc, v0, r0 from .
	! nvar = 5
	! allocate(on_fault_vars(nzt,nxt,nvar))
	
	! ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file. 
	! call check( nf90_open(infile1, NF90_NOWRITE, ncid))
	
	! ! Get the varid of the data variables, based on their names.
	! call check( nf90_inq_varid(ncid, "a", var_id(1)))
	! call check( nf90_inq_varid(ncid, "b", var_id(2)))
	! call check( nf90_inq_varid(ncid, "Dc", var_id(3)))
	! call check( nf90_inq_varid(ncid, "v0", var_id(4)))
	! call check( nf90_inq_varid(ncid, "r0", var_id(5)))
	! !call check( nf90_inq_varid(ncid, "init_slip_rate", var_id(6)))
	! !call check( nf90_inq_varid(ncid, "init_shear_stress", var_id(7)))
	! !call check( nf90_inq_varid(ncid, "init_normal_stress", var_id(8)))
	! !call check( nf90_inq_varid(ncid, "init_state", var_id(9)))
	
	! ! Read the data
	! do i = 1, nvar
		! call check( nf90_get_var(ncid, var_id(i), on_fault_vars(:,:,i)))
	! enddo		 
	! do i = 1, nxt
		! do j = 1, nzt
			! fric(9, (i-1)*nzt+j, 1) = on_fault_vars(j,i,1) ! a
			! fric(10, (i-1)*nzt+j, 1) = on_fault_vars(j,i,2)! b
			! fric(11, (i-1)*nzt+j, 1) = on_fault_vars(j,i,3)! Dc
			! fric(12, (i-1)*nzt+j, 1) = on_fault_vars(j,i,4)! v0
			! fric(13, (i-1)*nzt+j, 1) = on_fault_vars(j,i,5)! r0
			! !fric(46, (i-1)*nzt+j, 1) = on_fault_vars(j,i,6)! init_slip_rate
			! !fric(8, (i-1)*nzt+j, 1) = on_fault_vars(j,i,7)! shear
			! !fric(7, (i-1)*nzt+j, 1) = on_fault_vars(j,i,8)! norm
			! !fric(20, (i-1)*nzt+j, 1) = on_fault_vars(j,i,9)! norm
			! !fric(47, (i-1)*nzt+j, 1) = fric(46, (i-1)*nzt+j, 1)! peak slip rate
		! enddo 
	! enddo 
	! ! Close the file, freeing all resources.
	! call check( nf90_close(ncid))
	
	! deallocate(on_fault_vars)
	
	! ! Phase two, read in initial conditions from restart files fault.r.nc 
	! nvar = 12 
	! ! NOTE. the array structure is different than loading python generated nc file.
	! ! here we follow the structure of subroutine netcdf_write_on_fault.
	! ! on_fault_vars is now nxt by nzt!!!
	! allocate(on_fault_vars(nxt,nzt,nvar))
	! ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file. 
	! call check( nf90_open(infile2, NF90_NOWRITE, ncid))

	! ! Get the varid of the data variables, based on their names.
	! ! 'shear_strike', 'shear_dip', 'effective_normal', 'slip_rate' , 'state_variable', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs'
	! call check( nf90_inq_varid(ncid, "shear_strike", var_id(1)))
	! call check( nf90_inq_varid(ncid, "shear_dip", var_id(2)))
	! call check( nf90_inq_varid(ncid, "effective_normal", var_id(3)))
	! call check( nf90_inq_varid(ncid, "slip_rate", var_id(4)))
	! call check( nf90_inq_varid(ncid, "state_variable", var_id(5)))
	! call check( nf90_inq_varid(ncid, "state_normal", var_id(6)))
	! call check( nf90_inq_varid(ncid, "vxm", var_id(7)))
	! call check( nf90_inq_varid(ncid, "vym", var_id(8)))
	! call check( nf90_inq_varid(ncid, "vzm", var_id(9)))
	! call check( nf90_inq_varid(ncid, "vxs", var_id(10)))
	! call check( nf90_inq_varid(ncid, "vys", var_id(11)))
	! call check( nf90_inq_varid(ncid, "vzs", var_id(12)))
	! ! Read the data
	! do i = 1, nvar
		! call check( nf90_get_var(ncid, var_id(i), on_fault_vars(:,:,i)))
	! enddo		 
	
	! do i = 1, nxt
		! do j = 1, nzt
			! fric(8, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 1)! tstk0
			! fric(49, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 2)! tdip0
			! fric(7, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 3)! tnorm0
			! fric(46, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 4)! sliprate
			! fric(20, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 5)! state
			! fric(23, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 6)! state_normal
			! fric(31, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 7)! vxm
			! fric(32, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 8)! vym
			! fric(33, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 9)! vzm
			! fric(34, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 10)! vxs
			! fric(35, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 11)! vys
			! fric(36, (i-1)*nzt+j, 1) = on_fault_vars(i,j, 12)! vzs
			! fric(47, (i-1)*nzt+j, 1) = fric(46, (i-1)*nzt+j, 1)! peak slip rate
			! !fric(23, (i-1)*nzt+j, 1) = abs(fric(7, (i-1)*nzt+j, 1))! initialize theta_pc as abs(normal stress)
		! enddo 
	! enddo 
	! ! Close the file, freeing all resources.
	! call check( nf90_close(ncid))
	
	! deallocate(on_fault_vars)

! end subroutine netcdf_read_on_fault_restart

subroutine check(status)
	use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, nf90_strerror(status)
      stop "Stopped"
    end if
end subroutine check  
