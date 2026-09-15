"""Unit tests for scripts/convertFaultGeometry.

The utility's contract is narrow and checkable: given an arbitrary (x, z, y)
surface it must land EXACTLY on the fault grid a case asks for, keep the
supplied values where no interpolation is needed, and never emit a file that
scripts/lib.py:validateFaultRoughGeometry would reject. Expected values here
come from the analytic surface the fixture is built from, never from the
utility's own output.
"""
import importlib.machinery
import importlib.util
import os

import numpy as np
import pytest

import lib
from conftest import REPO_ROOT


def _load():
    # convertFaultGeometry has no .py suffix, like create.newcase.
    path = os.path.join(REPO_ROOT, 'scripts', 'convertFaultGeometry')
    loader = importlib.machinery.SourceFileLoader('convertFaultGeometry_uut', path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


cfg = _load()

# A smooth analytic surface on a 100 m grid covering 8 km x 4 km.
SRC_DX = 100.0
FXMIN, FXMAX = -4000.0, 4000.0
FZMIN, FZMAX = -4000.0, 0.0


def _analytic(x, z, amp=60.0):
    kx, kz = 2.0*np.pi/8000.0, 2.0*np.pi/4000.0
    return amp*np.sin(kx*x)*np.cos(kz*z)


def _writePoints(tmp_path, dx=SRC_DX, name='surface.txt', subsample=None,
                 seed=5):
    xs = FXMIN + dx*np.arange(round((FXMAX - FXMIN)/dx) + 1)
    zs = FZMIN + dx*np.arange(round((FZMAX - FZMIN)/dx) + 1)
    ZZ, XX = np.meshgrid(zs, xs, indexing='ij')
    rows = np.column_stack([XX.ravel(), ZZ.ravel(), _analytic(XX, ZZ).ravel()])
    if subsample is not None:
        rng = np.random.default_rng(seed)
        keep = rng.choice(rows.shape[0], size=subsample, replace=False)
        rows = rows[keep]
    path = str(tmp_path / name)
    np.savetxt(path, rows, fmt='%.9e')
    return path


def _validateKw(dx, **over):
    # dy is required by the validator and never defaulted (rule 2): it sets the
    # units of the element-tangling check.
    kw = dict(dx=dx, dz=dx, dy=dx, fxmin=FXMIN, fxmax=FXMAX, fzmin=FZMIN,
              fzmax=FZMAX, verbose=False)
    kw.update(over)
    return kw


def test_exact_decimation_keeps_the_supplied_values_bit_for_bit(tmp_path):
    src = _writePoints(tmp_path)
    out = str(tmp_path / 'g.txt')
    cfg.convertFaultGeometry(src, dx=500.0, dy=500.0, out=out, verbose=False)

    _, y, _, _ = lib.readFaultRoughGeometry(out)
    xs = FXMIN + 500.0*np.arange(y.shape[1])
    zs = FZMIN + 500.0*np.arange(y.shape[0])
    ZZ, XX = np.meshgrid(zs, xs, indexing='ij')
    # 500 m nodes are a subset of the 100 m input, so no interpolation may
    # happen: the values must be the analytic ones to write precision.
    assert y == pytest.approx(_analytic(XX, ZZ), abs=1e-6)
    assert y.shape == (9, 17)


def test_the_written_file_passes_the_case_setup_validator(tmp_path):
    src = _writePoints(tmp_path)
    out = str(tmp_path / 'g.txt')
    d = cfg.convertFaultGeometry(src, dx=500.0, dy=500.0, out=out, verbose=False)
    assert (d['nnx'], d['nnz']) == (17, 9)
    # and independently, run the validator again on the file on disk
    lib.validateFaultRoughGeometry(out, **_validateKw(500.0))


def test_derivative_columns_track_the_analytic_derivative(tmp_path):
    src = _writePoints(tmp_path)
    out = str(tmp_path / 'g.txt')
    cfg.convertFaultGeometry(src, dx=200.0, dy=200.0, out=out, verbose=False)
    _, _, dydx, _ = lib.readFaultRoughGeometry(out)
    xs = FXMIN + 200.0*np.arange(dydx.shape[1])
    zs = FZMIN + 200.0*np.arange(dydx.shape[0])
    ZZ, XX = np.meshgrid(zs, xs, indexing='ij')
    kx, kz = 2.0*np.pi/8000.0, 2.0*np.pi/4000.0
    exact = 60.0*kx*np.cos(kx*XX)*np.cos(kz*ZZ)
    # interior only: the boundary lines carry a first-order one-sided stencil
    assert dydx[:, 1:-1] == pytest.approx(exact[:, 1:-1], abs=2e-4)


def test_a_dx_the_input_was_never_sampled_at_is_interpolated(tmp_path):
    src = _writePoints(tmp_path)
    out = str(tmp_path / 'g.txt')
    # 250 m nodes are NOT a subset of the 100 m input grid.
    d = cfg.convertFaultGeometry(src, dx=250.0, dy=250.0, out=out, verbose=False)
    assert (d['nnx'], d['nnz']) == (33, 17)
    _, y, _, _ = lib.readFaultRoughGeometry(out)
    xs = FXMIN + 250.0*np.arange(y.shape[1])
    zs = FZMIN + 250.0*np.arange(y.shape[0])
    ZZ, XX = np.meshgrid(zs, xs, indexing='ij')
    # bilinear on a 100 m grid: error ~ (h^2/8)*|y''| = 0.15 m for this surface
    assert y == pytest.approx(_analytic(XX, ZZ), abs=0.3)


def test_scattered_points_are_resampled_onto_the_grid(tmp_path):
    src = _writePoints(tmp_path, subsample=1200)
    out = str(tmp_path / 'g.txt')
    d = cfg.convertFaultGeometry(src, dx=500.0, dy=500.0,
                                 fxmin=-3500.0, fxmax=3500.0,
                                 fzmin=-3500.0, fzmax=0.0,
                                 out=out, verbose=False)
    assert (d['nnx'], d['nnz']) == (15, 8)
    _, y, _, _ = lib.readFaultRoughGeometry(out)
    xs = -3500.0 + 500.0*np.arange(y.shape[1])
    zs = -3500.0 + 500.0*np.arange(y.shape[0])
    ZZ, XX = np.meshgrid(zs, xs, indexing='ij')
    assert y == pytest.approx(_analytic(XX, ZZ), abs=2.0)


def test_a_fault_larger_than_the_supplied_surface_is_refused(tmp_path):
    # rule 2: extrapolating a rough surface is not a thing to do quietly.
    src = _writePoints(tmp_path)
    with pytest.raises(lib.FaultGeometryError) as e:
        cfg.convertFaultGeometry(src, dx=500.0, dy=500.0, fxmin=-6000.0, fxmax=6000.0,
                                 out=str(tmp_path / 'g.txt'), verbose=False)
    assert 'outside the supplied surface' in str(e.value)


def test_a_dx_that_does_not_divide_the_fault_is_refused(tmp_path):
    src = _writePoints(tmp_path)
    with pytest.raises(lib.FaultGeometryError) as e:
        cfg.convertFaultGeometry(src, dx=300.0, dy=300.0,
                                 out=str(tmp_path / 'g.txt'), verbose=False)
    assert 'whole number' in str(e.value)


def test_non_finite_input_points_are_refused(tmp_path):
    src = _writePoints(tmp_path)
    rows = np.loadtxt(src)
    rows[17, 2] = np.nan
    np.savetxt(src, rows, fmt='%.9e')
    with pytest.raises(lib.FaultGeometryError) as e:
        cfg.convertFaultGeometry(src, dx=500.0, dy=500.0,
                                 out=str(tmp_path / 'g.txt'), verbose=False)
    assert 'NaN' in str(e.value)


def test_column_order_can_be_remapped(tmp_path):
    # Asserting only y.shape here was tautological: the shape comes from the
    # x/z columns' min/max, so an implementation that ignored `cols` entirely
    # and hard-coded y = data[:, 2] produced the same shape and passed. The
    # real claim is that remapping the columns recovers the SAME SURFACE, so
    # compare the values against the unswapped conversion.
    src = _writePoints(tmp_path)
    ref = str(tmp_path / 'ref.txt')
    cfg.convertFaultGeometry(src, dx=500.0, dy=500.0, out=ref, verbose=False)
    _, yRef, dxRef, dzRef = lib.readFaultRoughGeometry(ref)

    rows = np.loadtxt(src)
    swapped = str(tmp_path / 'swapped.txt')
    np.savetxt(swapped, rows[:, [0, 2, 1]], fmt='%.9e')   # x, y, z
    out = str(tmp_path / 'g.txt')
    cfg.convertFaultGeometry(swapped, dx=500.0, dy=500.0, cols=(0, 2, 1), out=out,
                             verbose=False)
    _, y, dydx, dydz = lib.readFaultRoughGeometry(out)

    assert y.shape == yRef.shape == (9, 17)
    np.testing.assert_allclose(y, yRef, rtol=0, atol=1e-6)
    np.testing.assert_allclose(dydx, dxRef, rtol=0, atol=1e-9)
    np.testing.assert_allclose(dydz, dzRef, rtol=0, atol=1e-9)
    # And the surface must not be trivially flat, or the comparison above
    # would hold for any implementation at all.
    assert np.abs(y).max() > 1.0
