"""Port of src/fortran/netcdf_io.f90's netcdf_read_on_fault_eqdyna_restart
(subroutine #2, the mode == 2 / EQsimu restart read of fault.r.nc).

The mode-1 reader (netcdf_read_on_fault_eqdyna, on_fault_vars_input.nc) stays
in readInputFiles.read_on_fault_vars, where it already lived.

Checkpoints, in the Fortran's order (netcdf_io.f90:129-202):
  C1  ntotft > 1 and dz /= dx (exact float equality)  -> ERR_CFG_RESTART_DZ_NE_DX (19)
  C2  open fault.r.nc, look up the 12 variables in the Fortran order
      (a missing file or variable -> ERR_NETCDF (71), 'NetCDF: ...')
  C3  fault count nftFile = len(dim 3) of shear_strike if it is rank 3, else 1;
      nftFile < ntotft or > ntotft -> ERR_INPUT_FILE_STALE (22), the two messages
  C4  per fault: (fnx, fnz) = nint(extent/spacing)+1 corner of its own slice
      (start (1,1,ift), count (fnx,fnz,1)); a rank-2 file is read as the
      corner of the whole variable
  C5  per fault node: ii/jj by round-half-away nint of the SLAVE node's x/z,
      then the 12 fric slots below

Slot map (globalvar.f90): shear_strike -> 8 INIT_STRIKE_SHEAR, shear_dip -> 49
INIT_DIP_SHEAR, effective_normal -> 7 INIT_NORM, slip_rate -> 47 PEAK_SLIPRATE,
state_variable -> 20 STATE, state_normal -> 23 THETA_PC, vxm/vym/vzm -> 31/32/33,
vxs/vys/vzs -> 34/35/36. After the read the Fortran flow runs init_vel
(eqdyna3d.f90:143), which seeds v1 at the fault nodes from slots 31-36;
assembleGlobalMass.init_vel is that, and driver.run now starts from it.

Deliberate differences from the Fortran, all HARD FAILURES where the Fortran
has undefined behaviour or a library error:
  - the (ii, jj) of a node outside its fault's (fnx, fnz) grid raises (Fortran
    indexes out of bounds, no check);
  - a fault's (fnx, fnz) larger than the file's corner raises 'NetCDF: Index
    exceeds dimension bound', as nf90_get_var does (numpy would truncate);
  - the dz/dx message prints repr() floats, Fortran G0 (same numbers, not the
    same digits);
  - a variable of rank other than 2 or 3 raises (Fortran would hand it to
    nf90_get_var and fail there).
"""
import numpy as np

from .checkInputConsistency import (ERR_CFG_RESTART_DZ_NE_DX, ERR_INPUT_FILE_STALE,
                                    ERR_NETCDF, InputConsistencyError)

# (variable name, fric slot) in netcdf_read_on_fault_eqdyna_restart's order.
RESTART_VARS = (
    ('shear_strike', 8), ('shear_dip', 49), ('effective_normal', 7),
    ('slip_rate', 47), ('state_variable', 20), ('state_normal', 23),
    ('vxm', 31), ('vym', 32), ('vzm', 33),
    ('vxs', 34), ('vys', 35), ('vzs', 36),
)


def _nint(x):
    """Fortran nint: round half away from zero."""
    return int(np.sign(x) * np.floor(np.abs(x) + 0.5)) if x != 0 else 0


def _netcdf_error(what):
    return InputConsistencyError(ERR_NETCDF, 'NetCDF: ' + what)


def read_on_fault_restart(nc_path, fric, ntotft, fxmin, fxmax, fzmin, fzmax,
                          dx, dz, meshCoor, nsmp):
    """Return a copy of `fric` ((nftnd+1, 101), 1-indexed rows and slots, as
    readInputFiles.read_on_fault_vars returns it) with the restart slots filled.

    fxmin/fxmax/fzmin/fzmax: per-fault sequences (index ift-1), length ntotft.
    meshCoor: (N+1, 3), 1-indexed rows. nsmp: (nftnd, 3) int64 [slave id,
    master id, 0-indexed fault], 1-indexed node ids, as build_node_coordinates
    returns it.
    """
    import netCDF4

    if ntotft > 1 and dz != dx:
        raise InputConsistencyError(
            ERR_CFG_RESTART_DZ_NE_DX,
            'multi-fault restart (mode = 2) refused: dz = %r but dx = %r; EQquasi '
            'spaces fault.r.nc dip rows by dx, EQdyna by dz, so the file would be '
            'read misaligned -- set dz equal to dx' % (dz, dx))
    for name, seq in (('fxmin', fxmin), ('fxmax', fxmax), ('fzmin', fzmin), ('fzmax', fzmax)):
        if len(seq) != ntotft:
            raise ValueError('read_on_fault_restart: %s has %d entries, ntotft = %d'
                             % (name, len(seq), ntotft))

    try:
        ds = netCDF4.Dataset(nc_path, 'r')
    except OSError as exc:
        raise _netcdf_error('cannot open %s: %s' % (nc_path, exc)) from None
    out = np.array(fric, copy=True)
    try:
        ds.set_auto_maskandscale(False)
        vars_ = []
        for vname, _slot in RESTART_VARS:
            if vname not in ds.variables:
                raise _netcdf_error('Variable not found (%s)' % vname)
            vars_.append(ds.variables[vname])
        ndims_var = vars_[0].ndim
        if ndims_var not in (2, 3):
            raise _netcdf_error('shear_strike has rank %d, expected 2 or 3' % ndims_var)
        nftFile = vars_[0].shape[0] if ndims_var == 3 else 1
        if nftFile < ntotft:
            raise InputConsistencyError(
                ERR_INPUT_FILE_STALE,
                'fault.r.nc holds %d fault(s) but ntotft = %d; write it with an '
                'EQquasi that carries the nid_fault dimension' % (nftFile, ntotft))
        if nftFile > ntotft:
            raise InputConsistencyError(
                ERR_INPUT_FILE_STALE,
                'fault.r.nc holds %d fault(s) but ntotft = %d; it was written for a '
                'case with more faults than this one' % (nftFile, ntotft))

        fault_of = np.asarray(nsmp[:, 2], dtype=np.int64)
        for ift in range(1, ntotft + 1):
            fnx = _nint((fxmax[ift - 1] - fxmin[ift - 1]) / dx) + 1
            fnz = _nint((fzmax[ift - 1] - fzmin[ift - 1]) / dz) + 1
            on_fault = np.empty((fnx, fnz, len(vars_)))   # Fortran (ii, jj, var)
            for k, v in enumerate(vars_):
                # netCDF4 shape is (fault, dip, strike) / (dip, strike): the
                # reverse of Fortran's (strike, dip, fault).
                if v.shape[-1] < fnx or v.shape[-2] < fnz:
                    raise _netcdf_error('Index exceeds dimension bound (%s: file %s, '
                                        'fault %d needs %d x %d)' % (RESTART_VARS[k][0],
                                                                      v.shape, ift, fnz, fnx))
                blk = v[ift - 1, :fnz, :fnx] if ndims_var == 3 else v[:fnz, :fnx]
                on_fault[:, :, k] = np.asarray(blk).T
            for i in np.nonzero(fault_of == ift - 1)[0]:
                row = int(i) + 1
                slave = int(nsmp[i, 0])
                ii = _nint((meshCoor[slave, 0] - fxmin[ift - 1]) / dx) + 1
                jj = _nint((meshCoor[slave, 2] - fzmin[ift - 1]) / dz) + 1
                if not (1 <= ii <= fnx and 1 <= jj <= fnz):
                    raise _netcdf_error('Index exceeds dimension bound (fault %d node row %d '
                                        'maps to (%d, %d) outside its %d x %d grid)'
                                        % (ift, row, ii, jj, fnx, fnz))
                for k, (_vname, slot) in enumerate(RESTART_VARS):
                    out[row, slot] = on_fault[ii - 1, jj - 1, k]
    finally:
        ds.close()
    return out
