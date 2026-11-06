#! /usr/bin/env python3
"""
Crop a real EQquasi fault.r.nc to the small restart fixture committed at
case_input/test.restart2/fault.r.nc (board row 16c).

EQdyna writes no restart dump of its own, so the only real source of a
fault.r.nc is an EQquasi run. The Fortran reader (netcdf_read_on_fault_eqdyna_
restart) and the python reader (eqdyna.netcdf_io) both take the (fnx, fnz)
corner of the file starting at index (1,1) with dip index 1 = fzmin, strike
index 1 = fxmin. This script cuts a rectangle [DIP0:DIP1, STRIKE0:STRIKE1] of
the source's (nid_dip, nid_strike) arrays and writes it as a new rank-2 file,
so the cut's first row/column becomes the reader's (1,1).

No value is altered: the 12 variables the reader needs are copied as stored
(float32), in the reader's names; the three slip variables the reader does not
read are not carried. Deterministic: no timestamp or history attribute is
written, so the same source and window give byte-identical output.

Usage: crop_restart_fixture.py SOURCE DEST   (source is opened read-only)
"""
import sys

import netCDF4
import numpy as np

VARS = ('shear_strike', 'shear_dip', 'effective_normal', 'slip_rate',
        'state_variable', 'state_normal', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs')
# Source: EQquasi bp1001.qd.250 run, Q3/fault.r.nc: 161 (dip) x 321 (strike),
# dx = 250 m, dip index 0 = 40 km deep, index 160 = free surface, strike index
# 160 = x = 0. The window is z in [-4 km, 0], x in [-10 km, +10 km]: 17 x 81
# fault nodes, centred on the source run's fastest-slipping surface patch.
# Dip rows 100-140 of the source (z -15 to -5 km) are locked at slip rates of
# 1e-28 to 1e-34 m/s, where EQdyna's aging law (fric.f90
# rate_state_ageing_law, theta = L/V + (theta - L/V)*exp(-V dt/L)) cancels to
# theta = 0 and goes NaN at step 2-5; the window stays above that band.
DIP0, DIP1 = 144, 161
STRIKE0, STRIKE1 = 120, 201


def crop(src_path, dst_path, dip0=DIP0, dip1=DIP1, strike0=STRIKE0, strike1=STRIKE1):
    src = netCDF4.Dataset(src_path, 'r')
    try:
        for v in VARS:
            if v not in src.variables:
                raise KeyError('source %s lacks variable %s' % (src_path, v))
        if not (0 <= dip0 < dip1 <= src.dimensions['nid_dip'].size and
                0 <= strike0 < strike1 <= src.dimensions['nid_strike'].size):
            raise ValueError('window [%d:%d, %d:%d] outside the source %dx%d' % (
                dip0, dip1, strike0, strike1, src.dimensions['nid_dip'].size,
                src.dimensions['nid_strike'].size))
        dst = netCDF4.Dataset(dst_path, 'w', format='NETCDF3_CLASSIC')
        try:
            dst.createDimension('nid_dip', dip1 - dip0)
            dst.createDimension('nid_strike', strike1 - strike0)
            for dim, lo, hi in (('nid_dip', dip0, dip1), ('nid_strike', strike0, strike1)):
                cv = dst.createVariable(dim, 'i4', (dim,))
                cv.units = src.variables[dim].units
                cv[:] = np.asarray(src.variables[dim][lo:hi])
            for v in VARS:
                sv = src.variables[v]
                sv.set_auto_maskandscale(False)
                dv = dst.createVariable(v, sv.dtype, ('nid_dip', 'nid_strike'))
                dv.units = sv.units
                dv[:, :] = sv[dip0:dip1, strike0:strike1]
        finally:
            dst.close()
    finally:
        src.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    crop(sys.argv[1], sys.argv[2])
