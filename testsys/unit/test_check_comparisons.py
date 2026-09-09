"""Unit tests for check.test.py's compare_txt_files / compare_nc_files.

These are pure functions (no import-time side effects; see check.test.py's
module docstring), so we import the real module from the repo root and
call them directly against small synthetic fixtures built in tmp_path --
never against the live test/ or test.reference.results/ trees.
"""
import importlib.machinery
import importlib.util
import os

import numpy as np
import xarray as xr

from conftest import REPO_ROOT


def _load_check_test():
    path = os.path.join(REPO_ROOT, 'check.test.py')
    loader = importlib.machinery.SourceFileLoader('check_test_under_test', path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


check_test = _load_check_test()


def _write_frt(path, values):
    with open(path, 'w') as f:
        f.write(' '.join(str(v) for v in values))


def test_compare_txt_files_identical_reports_success(tmp_path):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    _write_frt(a, [1.0, 2.0, 3.0])
    _write_frt(b, [1.0, 2.0, 3.0])

    result = check_test.compare_txt_files(str(a), str(b))

    assert result.startswith('SUCCESS')


def test_compare_txt_files_value_beyond_threshold_reports_fail(tmp_path):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    _write_frt(a, [1.0, 2.0, 3.0])
    _write_frt(b, [1.0, 2.0, 3.01])  # 0.01 > THRESHOLD (1e-3)

    result = check_test.compare_txt_files(str(a), str(b))

    assert result.startswith('FAIL')


def test_compare_txt_files_within_threshold_reports_success(tmp_path):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    _write_frt(a, [1.0, 2.0, 3.0])
    _write_frt(b, [1.0, 2.0, 3.0 + 5e-4])  # under THRESHOLD (1e-3)

    result = check_test.compare_txt_files(str(a), str(b))

    assert result.startswith('SUCCESS')


def test_compare_txt_files_length_mismatch_reports_fail(tmp_path):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    _write_frt(a, [1.0, 2.0, 3.0])
    _write_frt(b, [1.0, 2.0])

    result = check_test.compare_txt_files(str(a), str(b))

    assert result.startswith('FAIL')


def _make_dataset(value):
    return xr.Dataset({'slip': (('node',), np.array([value, value, value]))})


def test_compare_nc_files_identical_reports_success(tmp_path):
    fn1 = str(tmp_path / 'ref.nc')
    fn2 = str(tmp_path / 'test.nc')
    _make_dataset(1.23).to_netcdf(fn1)
    _make_dataset(1.23).to_netcdf(fn2)

    result = check_test.compare_nc_files(fn1, fn2)

    assert result.startswith('SUCCESS')


def test_compare_nc_files_value_beyond_threshold_reports_fail(tmp_path):
    fn1 = str(tmp_path / 'ref.nc')
    fn2 = str(tmp_path / 'test.nc')
    _make_dataset(1.0).to_netcdf(fn1)
    _make_dataset(2.0).to_netcdf(fn2)  # far beyond rtol/atol=1e-3

    result = check_test.compare_nc_files(fn1, fn2)

    assert result.startswith('FAIL')
