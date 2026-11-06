"""Tests for the mode = 2 restart fixture and its crop script.

The committed case_input/test.restart2/fault.r.nc must hold the 12 variables
the Fortran and python restart readers consume, as float32 rank-2 arrays whose
shape equals the case's fault node grid, and testsys/crop_restart_fixture.py
must reproduce a source window exactly and deterministically.
"""
import os

import netCDF4
import numpy as np
import pytest

from conftest import REPO_ROOT
from testsys import crop_restart_fixture as crf

FIXTURE = os.path.join(REPO_ROOT, 'case_input', 'test.restart2', 'fault.r.nc')


def test_committed_fixture_shape_vars_and_size():
    assert os.path.getsize(FIXTURE) <= 1_000_000
    with netCDF4.Dataset(FIXTURE) as ds:
        assert ds.dimensions['nid_dip'].size == crf.DIP1 - crf.DIP0 == 17
        assert ds.dimensions['nid_strike'].size == crf.STRIKE1 - crf.STRIKE0 == 81
        for v in crf.VARS:
            a = np.asarray(ds.variables[v][:])
            assert a.dtype == np.float32, v
            assert a.shape == (17, 81), v
            assert np.all(np.isfinite(a)), v
        assert float(np.min(ds.variables['slip_rate'][:])) > 0.0


def _write_source(path, nd=6, ns=7):
    rng = np.random.default_rng(7)
    with netCDF4.Dataset(path, 'w', format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('nid_dip', nd)
        ds.createDimension('nid_strike', ns)
        for dim, n in (('nid_dip', nd), ('nid_strike', ns)):
            cv = ds.createVariable(dim, 'i4', (dim,))
            cv.units = 'index'
            cv[:] = np.arange(n)
        for v in crf.VARS + ('extra_slip',):
            x = ds.createVariable(v, 'f4', ('nid_dip', 'nid_strike'))
            x.units = 'u'
            x[:, :] = rng.standard_normal((nd, ns)).astype(np.float32)


def test_crop_copies_window_exactly_and_deterministically(tmp_path):
    src, d1, d2 = (str(tmp_path / n) for n in ('s.nc', 'a.nc', 'b.nc'))
    _write_source(src)
    crf.crop(src, d1, 1, 5, 2, 6)
    crf.crop(src, d2, 1, 5, 2, 6)
    assert open(d1, 'rb').read() == open(d2, 'rb').read()
    with netCDF4.Dataset(src) as s, netCDF4.Dataset(d1) as d:
        assert 'extra_slip' not in d.variables
        for v in crf.VARS:
            assert np.array_equal(d.variables[v][:], s.variables[v][1:5, 2:6]), v


def test_crop_refuses_missing_variable_and_bad_window(tmp_path):
    src = str(tmp_path / 's.nc')
    _write_source(src)
    with pytest.raises(ValueError):
        crf.crop(src, str(tmp_path / 'x.nc'), 0, 99, 0, 3)
    bad = str(tmp_path / 'bad.nc')
    with netCDF4.Dataset(bad, 'w', format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('nid_dip', 2)
        ds.createDimension('nid_strike', 2)
    with pytest.raises(KeyError):
        crf.crop(bad, str(tmp_path / 'y.nc'), 0, 1, 0, 1)
