"""Unit tests for scripts/lib.py (fault-loading and boxcar-taper helpers).

Every expected value below is computed independently of the function under
test -- either an exact closed form (asinh(1) = ln(1+sqrt(2))) or a
documented boundary/symmetry property of the boxcar tapers -- never a copy
of the function's own output.
"""
import math
import types

import numpy as np
import pytest

import lib


# ---- B1: even (symmetric) boxcar in x, full width 2*ww, taper width w ----

def test_B1_is_flat_one_inside_plateau():
    assert lib.B1(0.0, 10.0, 2.0) == 1.0
    assert lib.B1(10.0, 10.0, 2.0) == 1.0  # boundary is included (<=ww)


def test_B1_is_zero_beyond_taper():
    assert lib.B1(12.0, 10.0, 2.0) == 0.0  # exactly at ww+w
    assert lib.B1(50.0, 10.0, 2.0) == 0.0


def test_B1_is_even_in_x():
    # B1 only ever reads abs(x); a sign-handling bug would break this.
    for x in (0.5, 5.0, 10.5, 11.9, 20.0):
        assert lib.B1(x, 10.0, 2.0) == lib.B1(-x, 10.0, 2.0)


def test_B1_taper_is_strictly_between_bounds_inside_transition():
    res = lib.B1(11.0, 10.0, 2.0)  # strictly inside (ww, ww+w)
    assert 0.0 < res < 1.0


# ---- B2: one-sided (y>=0) boxcar with a taper to 0 at the free surface ----

def test_B2_negative_y_exits_cleanly_not_with_nameerror():
    # Regression: lib.py used to do `from sys import *`, which does not bind
    # the name `sys` itself -- the guard clause's `sys.exit()` therefore
    # raised NameError instead of exiting, on this exact code path.
    with pytest.raises(SystemExit):
        lib.B2(-1.0, 15.0, 3.0)


def test_B3_negative_y_exits_cleanly_not_with_nameerror():
    with pytest.raises(SystemExit):
        lib.B3(-1.0, 15.0, 3.0)


def test_B2_tapers_to_zero_at_the_free_surface():
    # y=0 sits inside the near-surface taper branch (y<w): documented to
    # go to 0, not 1 -- this is what distinguishes B2 from B3.
    assert lib.B2(0.0, 15.0, 3.0) == 0.0


def test_B2_is_flat_one_on_the_plateau():
    assert lib.B2(3.0, 15.0, 3.0) == 1.0   # y == w
    assert lib.B2(10.0, 15.0, 3.0) == 1.0  # w < y < ww
    assert lib.B2(15.0, 15.0, 3.0) == 1.0  # y == ww


def test_B2_is_zero_beyond_far_taper():
    assert lib.B2(18.0, 15.0, 3.0) == 0.0  # y == ww+w
    assert lib.B2(100.0, 15.0, 3.0) == 0.0


# ---- B3: no near-surface taper; flat 1 from y=0 up to ww ----

def test_B3_is_flat_one_from_surface_to_plateau_edge():
    assert lib.B3(0.0, 15.0, 3.0) == 1.0
    assert lib.B3(15.0, 15.0, 3.0) == 1.0


def test_B3_is_zero_beyond_far_taper():
    assert lib.B3(18.0, 15.0, 3.0) == 0.0


# ---- linear1: piecewise-linear boxcar; exact midpoint is hand-computable ----

def test_linear1_plateau_and_zero_regions():
    assert lib.linear1(0.0, 10.0, 2.0) == 1.0
    assert lib.linear1(12.0, 10.0, 2.0) == 0.0


def test_linear1_taper_midpoint_is_exactly_half():
    # At x = ww + w/2, linear1 = 1 - (w/2)/w = 0.5 exactly, by construction
    # of a piecewise-*linear* ramp -- independent of lib.py's own formula.
    res = lib.linear1(11.0, 10.0, 2.0)
    assert res == 0.5


# ---- shear_steady_state: exact closed form via asinh(1) = ln(1+sqrt(2)) ----

def test_shear_steady_state_matches_closed_form_asinh_one():
    # Choose load_rate == v0 (kills the b*log(v0/load_rate) term) and
    # r0 == 0 (kills the remaining exponent), so the asinh argument
    # collapses to slip_rate/(2*v0); choosing slip_rate = 2*v0 makes that
    # argument exactly 1.
    a, norm = 0.01, -1.0
    v0 = 1e-6
    res = lib.shear_steady_state(a=a, b=0.02, v0=v0, r0=0.0,
                                  load_rate=v0, norm=norm, slip_rate=2 * v0)
    expected = -norm * a * math.log(1 + math.sqrt(2))  # asinh(1), closed form
    assert math.isclose(res, expected, rel_tol=1e-12)


# ---- state_steady_state: friclaw 3 and 4/5 exact closed forms ----

def test_state_steady_state_friclaw3_matches_hand_computed_value():
    # shear/norm/a = 1 and r0 chosen to cancel log(2*sinh(1)) exactly,
    # slip_rate == v0 kills the remaining log term, so state = d0/v0.
    r0 = math.log(2 * math.sinh(1.0))
    state = lib.state_steady_state(a=1.0, b=1.0, d0=1.0, v0=2.0, r0=r0,
                                    shear=1.0, norm=1.0, slip_rate=2.0,
                                    friclaw=3)
    assert math.isclose(state, 0.5, rel_tol=1e-12)


def test_state_steady_state_friclaw5_matches_hand_computed_value():
    # a=1, v0=slip_rate=1, shear/norm/a=1 -> state = log(2*sinh(1)).
    state = lib.state_steady_state(a=1.0, b=1.0, d0=1.0, v0=1.0, r0=0.0,
                                    shear=1.0, norm=1.0, slip_rate=1.0,
                                    friclaw=5)
    expected = math.log(2 * math.sinh(1.0))
    assert math.isclose(state, expected, rel_tol=1e-12)


# ---- sort_nicely / alphanum_key: numeric-aware filename ordering ----
# (moved from scripts/plot_on_fault_vars and plot_on_fault_vars2, which
#  duplicated these verbatim; plot_on_fault_vars now imports them from lib)

def test_tryint_converts_digits_and_passes_through_other_strings():
    assert lib.tryint("42") == 42
    assert lib.tryint("abc") == "abc"


def test_sort_nicely_orders_embedded_numbers_numerically_not_lexically():
    # Plain string sort would put "fault.10.nc" before "fault.2.nc";
    # sort_nicely must not.
    names = ["fault.10.nc", "fault.2.nc", "fault.1.nc"]
    assert lib.sort_nicely(names) == ["fault.1.nc", "fault.2.nc", "fault.10.nc"]


def test_sort_nicely_sorts_in_place_and_returns_the_same_list():
    names = ["b2", "b10", "b1"]
    result = lib.sort_nicely(names)
    assert result is names
    assert result == ["b1", "b2", "b10"]


# ---- loadFrtData: frt.txt* loader shared by plotRuptureDynamics/plotSlipAndRPT ----

def _make_par(tmp_path):
    """Minimal par-like object covering exactly what loadFrtData reads."""
    fxmin, fxmax, dx = 0.0, 1.0, 1.0
    fzmin, fzmax, dz = 0.0, 0.0, 1.0
    na = round((fxmax - fxmin) / dx + 1)  # 2
    ma = round((fzmax - fzmin) / dz + 1)  # 1
    return types.SimpleNamespace(
        nx=1, ny=1, nz=1,
        fxmin=fxmin, fxmax=fxmax, dx=dx,
        fzmin=fzmin, fzmax=fzmax, dz=dz,
        dip=90.0,  # sin(90 deg) == 1, keeps along-dip arithmetic trivial
        fx=np.linspace(fxmin, fxmax, na),
        fz=np.linspace(fzmin, fzmax, ma),
    )


def test_loadFrtData_grids_two_nodes_from_a_synthetic_frt_file(tmp_path, monkeypatch):
    par = _make_par(tmp_path)
    # Columns (0-indexed, matching lib.loadFrtData's docstring, 1-indexed there):
    # 0 x, 1 y, 2 z, 3 rupture time, 4 slip_s, 5 slip_d, 6-8 unused,
    # 9 peak slip rate, 10 final slip rate, 11 normal, 12 shear, 13 dip shear,
    # 14-19 vxm,vym,vzm,vxs,vys,vzs, 20 state, 21 state_normal.
    row0 = [0, 0, 0, 1.0, 3.0, 4.0, 0, 0, 0, 5.0, 6.0, 7.0, 8.0, 9.0,
            10, 11, 12, 13, 14, 15, 16.0, 17.0]
    row1 = [1, 0, 0, 2.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0, 0, 0, 0, 0.0, 0.0]
    monkeypatch.chdir(tmp_path)
    np.savetxt(tmp_path / "frt.txt0", np.array([row0, row1]))

    xx, zz, rupt, rupt2d, fVarArr, magnitude = lib.loadFrtData(par)

    assert rupt2d.shape == (1, 2, 100)
    assert fVarArr.shape == (1, 2, 100)
    assert rupt2d[0, 0, 1] == pytest.approx(5.0)   # slip magnitude = hypot(3,4)
    assert rupt2d[0, 0, 0] == pytest.approx(1.0)   # rupture time
    assert fVarArr[0, 0, 4] == pytest.approx(16.0)  # state_variable (col 20)
    assert fVarArr[0, 0, 5] == pytest.approx(17.0)  # state_normal (col 21)
    assert rupt[0, 2] == pytest.approx(1.0)         # rupture time for ii=0
    assert rupt[1, 2] == pytest.approx(2.0)         # rupture time for ii=1
    assert magnitude != 0.0  # moment accumulated from the non-zero node


# ---- bFault_Rough_Geometry.txt validation --------------------------------
#
# The "good" fixture is an ANALYTIC surface -- y and both derivative columns
# come from closed-form expressions, never from the finite differences the
# validator itself uses -- and every failure case below is that same file with
# one thing perturbed, so no test can pass by echoing the implementation.
# The shipped official SCEC TPV29 surface is exercised too: v5.6.0's validator
# must accept the gated cases' files unchanged.

import os

GEOM_FXMIN, GEOM_FXMAX = -5000.0, 5000.0
GEOM_FZMIN, GEOM_FZMAX = -5000.0, 0.0
GEOM_DX = 500.0
TPV29_SHIPPED = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    'case_input', 'test.tpv29', 'bFault_Rough_Geometry.tpv29.50m.txt')


def _analyticSurface(amp=100.0, dx=GEOM_DX):
    """y = amp*sin(kx x)*cos(kz z) with kx, kz one full period across the
    fault, and its exact derivatives. Peak |slope| = amp*kx."""
    nnx = round((GEOM_FXMAX - GEOM_FXMIN)/dx) + 1
    nnz = round((GEOM_FZMAX - GEOM_FZMIN)/dx) + 1
    x = GEOM_FXMIN + dx*np.arange(nnx)
    z = GEOM_FZMIN + dx*np.arange(nnz)
    ZZ, XX = np.meshgrid(z, x, indexing='ij')
    kx = 2.0*np.pi/(GEOM_FXMAX - GEOM_FXMIN)
    kz = 2.0*np.pi/(GEOM_FZMAX - GEOM_FZMIN)
    y = amp*np.sin(kx*XX)*np.cos(kz*ZZ)
    dydx = amp*kx*np.cos(kx*XX)*np.cos(kz*ZZ)
    dydz = -amp*kz*np.sin(kx*XX)*np.sin(kz*ZZ)
    return y, dydx, dydz


def _geomKwargs(dx=GEOM_DX, **over):
    # dy is supplied explicitly. It used to be omitted, and the validator
    # silently substituted dx for it -- so every offset assertion below was
    # exercising the fallback, not the check. The fallback is gone (the
    # validator now SKIPS the offset check and says so when dy is missing),
    # which is what surfaced this.
    kw = dict(dx=dx, dz=dx, dy=dx, fxmin=GEOM_FXMIN, fxmax=GEOM_FXMAX,
              fzmin=GEOM_FZMIN, fzmax=GEOM_FZMAX,
              ymin=-20000.0, ymax=20000.0, verbose=False)
    kw.update(over)
    return kw


def _writeGoodGeometry(tmp_path, amp=100.0, dx=GEOM_DX,
                       fxmin=GEOM_FXMIN, fzmin=GEOM_FZMIN):
    y, dydx, dydz = _analyticSurface(amp, dx)
    path = str(tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE)
    lib.writeFaultRoughGeometry(y, dydx, dydz, dx, fxmin, fzmin, fname=path)
    return path


def _perturbLines(path, mutate):
    """Rewrite the file with `mutate` applied to its list of lines."""
    with open(path) as f:
        lines = f.readlines()
    with open(path, 'w') as f:
        f.writelines(mutate(lines))
    return path


def test_validate_accepts_a_clean_analytic_surface(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    d = lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert (d['nnx'], d['nnz']) == (21, 11)
    assert d['warnings'] == []          # peak slope 0.126, under the 0.2 warn
    # maxSlope is the peak of the analytic derivative fields over the grid
    _, adx, adz = _analyticSurface()
    assert d['maxSlope'] == pytest.approx(
        max(np.abs(adx).max(), np.abs(adz).max()), rel=1e-6)
    assert 0.1 < d['maxSlope'] < 0.2


def test_validate_accepts_the_official_tpv29_surface_unchanged():
    # The gated test.tpv29 ships this file (the 50 m spec sampling);
    # validation must not reject it. It carries SCEC's ANALYTIC
    # derivatives, not finite differences, so this is also the real-world
    # check on the derivative tolerance: at 50 m the residual against a
    # central difference is 7e-4, against a bound of 4x the curvature the
    # derivative column itself implies.
    d = lib.validateFaultRoughGeometry(
        TPV29_SHIPPED, dx=50.0, dz=50.0, dy=50.0,
        fxmin=-20000.0, fxmax=20000.0, fzmin=-20000.0, fzmax=0.0,
        ymin=-30000.0, ymax=30000.0, verbose=False)
    assert (d['nnx'], d['nnz']) == (801, 401)
    assert d['maxSlope'] == pytest.approx(0.6615, abs=1e-3)
    assert len(d['warnings']) == 1      # 0.66 > 0.2: a warning, never a failure
    # With dy supplied the hard offset check actually runs on the real file.
    # cot(dip) style reasoning does not apply here -- this is a rough surface,
    # and its steepest per-cell climb must stay below one fault-normal cell.
    assert d['maxOffset'] < 1.0
    assert d['maxOffset'] == pytest.approx(0.6615, abs=1e-3)   # dx == dy here


def test_reader_unpacks_the_z_fastest_row_order(tmp_path):
    # rough_geo(:, nnz*(ix-1) + iz) in insertFaultInterface: z runs fastest.
    y = np.arange(6, dtype=float).reshape(2, 3)      # (nnz=2, nnx=3)
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, np.zeros_like(y), np.zeros_like(y),
                                100.0, 0.0, 0.0, fname=path)
    raw = np.loadtxt(path, skiprows=2)
    assert list(raw[:, 0]) == [0.0, 3.0, 1.0, 4.0, 2.0, 5.0]   # (ix, iz) order
    header, back, _, _ = lib.readFaultRoughGeometry(path)
    assert (header['nnx'], header['nnz']) == (3, 2)
    assert np.array_equal(back, y)


def test_missing_file_is_rejected_with_the_way_to_produce_one():
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.readFaultRoughGeometry('no_such_bFault_Rough_Geometry.txt')
    assert 'insertFaultType' in str(e.value)


def test_wrong_nx_in_the_header_is_rejected_naming_both_counts(tmp_path):
    path = _writeGoodGeometry(tmp_path)

    def mutate(lines):
        lines[0] = '22\t11\t0\n'          # was 21 x 11
        return lines
    _perturbLines(path, mutate)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    msg = str(e.value)
    assert '22' in msg and '21' in msg


def test_wrong_nz_in_the_header_is_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: ['21\t12\t0\n'] + l[1:])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert 'nnz' in str(e.value)


def test_wrong_dx_in_the_header_is_rejected(tmp_path):
    # A surface sampled at 250 m dropped into a 500 m case: same row count,
    # same corner, silently stretched by 2x if nobody looks.
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l[:1] + ['250.0\t-5000.0\t-5000.0\n'] + l[2:])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert '250.0' in str(e.value) and '500.0' in str(e.value)


def test_wrong_origin_is_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l[:1] + ['500.0\t-4000.0\t-6000.0\n'] + l[2:])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    msg = str(e.value)
    assert 'fxmin' in msg and 'fzmin' in msg


def test_short_file_is_rejected(tmp_path):
    # The exact shape of the original hazard: the Fortran reads nnx*nnz rows
    # from a file that does not have them.
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l[:-1])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    msg = str(e.value)
    assert '230' in msg and '231' in msg      # 21*11 = 231 rows expected


def test_long_file_is_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l + [l[-1]])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert '232' in str(e.value)


def test_ragged_rows_are_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l[:-1] + ['1.0\t2.0\n'])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert '3' in str(e.value)


def test_nan_in_the_surface_column_is_rejected_with_its_location(tmp_path):
    path = _writeGoodGeometry(tmp_path)

    def mutate(lines):
        lines[7] = 'nan\t0.0\t0.0\n'          # data row 6 -> ix=0, iz=5
        return lines
    _perturbLines(path, mutate)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    msg = str(e.value)
    assert 'non-finite' in msg and 'ix=0' in msg and 'iz=5' in msg


def test_inf_in_a_derivative_column_is_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: l[:5] + ['0.0\tinf\t0.0\n'] + l[6:])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert 'dy/dx' in str(e.value)


def test_derivative_column_scaled_by_two_is_rejected(tmp_path):
    # The classic stale/mis-scaled-derivative bug: the surface column is fine,
    # the derivative column is not the derivative of it.
    y, dydx, dydz = _analyticSurface()
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, 2.0*dydx, dydz, GEOM_DX,
                                GEOM_FXMIN, GEOM_FZMIN, fname=path)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert 'dy/dx' in str(e.value) and 'central difference' in str(e.value)


def test_derivative_column_from_a_different_surface_is_rejected(tmp_path):
    # Surface updated, derivatives left behind (both are individually
    # plausible fields -- only their relationship is wrong).
    y, _, dydz = _analyticSurface(amp=100.0)
    _, dydxOther, _ = _analyticSurface(amp=80.0)
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, dydxOther, dydz, GEOM_DX,
                                GEOM_FXMIN, GEOM_FZMIN, fname=path)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert 'dy/dx' in str(e.value)


def test_swapped_derivative_columns_are_rejected(tmp_path):
    y, dydx, dydz = _analyticSurface()
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, dydz, dydx, GEOM_DX,
                                GEOM_FXMIN, GEOM_FZMIN, fname=path)
    with pytest.raises(lib.FaultGeometryError):
        lib.validateFaultRoughGeometry(path, **_geomKwargs())


def test_a_correct_but_rough_surface_only_warns(tmp_path):
    # peak slope amp*kz = 0.63: past generateFaultInterface's 0.2 warning
    # but nowhere near element tangling. test.drv.a6 (0.38) and test.tpv29
    # (0.65) both live here, so this MUST stay a warning.
    path = _writeGoodGeometry(tmp_path, amp=500.0)   # amp*kz = 0.628
    d = lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert len(d['warnings']) == 1
    assert 0.2 < d['maxSlope'] < 1.0


def test_excessive_per_cell_offset_is_rejected(tmp_path):
    # amp*kz = 2.5 with dz == dy: adjacent fault nodes are offset by more than
    # a full fault-normal cell, so the inserted element layer tangles.
    path = _writeGoodGeometry(tmp_path, amp=2000.0)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs(dy=GEOM_DX))
    assert 'tangle' in str(e.value)


def test_a_dipping_planar_fault_is_not_mistaken_for_a_tangled_mesh(tmp_path):
    # THE false positive a raw-slope test would produce. A planar fault dipping
    # at 45 deg has dz = dx*sin(45) and dy/dz = cot(45) = 1 exactly -- while its
    # actual per-row offset is dx*cos(45) = 0.71 of a dy = dx cell, and nothing
    # tangles at any dip. test.tpv10 is the gated case in this family (60 deg,
    # raw |dy/dz| = 0.577). Build the surface the way generateFaultInterface
    # does for insertFaultType=1 and check every dip from 80 down to 20 deg.
    for dipDeg in (80.0, 60.0, 45.0, 30.0, 20.0):
        dip = math.radians(dipDeg)
        dx = 500.0
        dz = dx*math.sin(dip)
        nnx, nnz = 11, 9
        fzmin = -dz*(nnz - 1)
        iz = np.arange(nnz)[:, None]*np.ones(nnx)
        y = (nnz - 1 - iz)*dx*math.cos(dip)
        dydx = np.zeros_like(y)
        dydz = np.full_like(y, -dx*math.cos(dip)/dz)      # = -cot(dip)
        path = str(tmp_path / f'dip{dipDeg:g}.txt')
        lib.writeFaultRoughGeometry(y, dydx, dydz, dx, -2500.0, fzmin,
                                    fname=path)
        d = lib.validateFaultRoughGeometry(
            path, dx=dx, dz=dz, dy=dx, fxmin=-2500.0, fxmax=2500.0,
            fzmin=fzmin, fzmax=0.0, ymin=-20000.0, ymax=20000.0, verbose=False)
        # rel 1e-7: the file carries 7 significant digits ('%.7e')
        assert d['maxSlope'] == pytest.approx(1.0/math.tan(dip), rel=1e-7)
        assert d['maxOffset'] == pytest.approx(math.cos(dip), rel=1e-7)
        assert d['maxOffset'] < 1.0


def test_surface_outside_the_model_domain_is_rejected(tmp_path):
    # Roughness supplied in km instead of m is the usual way this happens.
    path = _writeGoodGeometry(tmp_path, amp=100.0)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path,
                                       **_geomKwargs(ymin=-50.0, ymax=50.0))
    assert 'model domain' in str(e.value)


def test_fault_extent_not_a_whole_number_of_cells_is_rejected(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs(dx=300.0, dz=300.0))
    assert 'whole number of cells' in str(e.value)


def test_case_nfx_inconsistent_with_its_own_extent_is_reported(tmp_path):
    path = _writeGoodGeometry(tmp_path)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs(nfx=20))
    assert 'par.nfx' in str(e.value)


def test_every_problem_is_reported_in_one_pass(tmp_path):
    # One case.setup should tell the whole story, not one item per re-run.
    path = _writeGoodGeometry(tmp_path)
    _perturbLines(path, lambda l: ['22\t11\t0\n', '250.0\t-4000.0\t-5000.0\n'] + l[2:])
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert str(e.value).count('\n  - ') >= 3


def test_forCase_wrapper_reads_the_grid_off_par(tmp_path, monkeypatch):
    path = _writeGoodGeometry(tmp_path)
    monkeypatch.chdir(tmp_path)
    par = types.SimpleNamespace(
        dx=GEOM_DX, dz=GEOM_DX, dy=GEOM_DX,
        fxmin=GEOM_FXMIN, fxmax=GEOM_FXMAX,
        fzmin=GEOM_FZMIN, fzmax=GEOM_FZMAX,
        ymin=-20000.0, ymax=20000.0, nfx=21, nfz=11)
    d = lib.validateFaultRoughGeometryForCase(par, verbose=False)
    assert d['nnx'] == 21
    par.dx = 250.0                       # case moved, file did not
    with pytest.raises(lib.FaultGeometryError):
        lib.validateFaultRoughGeometryForCase(par, verbose=False)
    assert os.path.basename(path) == lib.FAULT_ROUGH_GEOMETRY_FILE


# ---- geometry-source resolution guard (pathway_forward item 31) -----------
#
# test.tpv29 shipped its surface at 100 m while its FULL_SPECS entry asks for
# the 50 m spec resolution, so the full tier could not even be set up: the
# compset's own decimator raised "dx=50.0 is not an integer multiple of the
# shipped geometry spacing 100.0 m", which is true and useless. The general
# form is "the case wants a dx its geometry source cannot supply".

def test_a_dx_finer_than_the_supplied_source_is_refused():
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.requireFaultGeometryResolution(25.0, 50.0, 'the shipped surface')
    msg = str(e.value)
    assert '25.0' in msg and '50.0' in msg      # requested and native spacing
    assert 'never refined' in msg               # and why
    assert '100.0' in msg                       # and a dx that would work


def test_a_dx_that_is_not_a_multiple_of_the_source_spacing_is_refused():
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.requireFaultGeometryResolution(75.0, 50.0, 'the shipped surface')
    msg = str(e.value)
    # The message must name the requested dx, the source spacing, and a dx
    # that WOULD work. The old second assertion was `'100.0' in msg or
    # '50.0' in msg`, which the first assertion already guarantees -- it
    # could not fail.
    assert '75.0' in msg, 'message must name the requested dx'
    assert '50.0' in msg, 'message must name the source spacing'
    assert '100.0' in msg, 'message must suggest a dx that actually works'


def test_exact_multiples_of_the_source_spacing_are_allowed():
    for dx in (50.0, 100.0, 200.0, 250.0, 500.0):
        lib.requireFaultGeometryResolution(dx, 50.0, 'the shipped surface')


def test_no_declared_source_means_no_resolution_constraint():
    # Generated geometry (insertFaultType 1/2) has no fixed source sampling.
    lib.requireFaultGeometryResolution(37.0, None)


def test_the_case_validator_applies_the_source_resolution_guard(tmp_path,
                                                                monkeypatch):
    _writeGoodGeometry(tmp_path)
    monkeypatch.chdir(tmp_path)
    par = types.SimpleNamespace(
        dx=GEOM_DX, dz=GEOM_DX, dy=GEOM_DX,
        fxmin=GEOM_FXMIN, fxmax=GEOM_FXMAX,
        fzmin=GEOM_FZMIN, fzmax=GEOM_FZMAX,
        ymin=-20000.0, ymax=20000.0, nfx=21, nfz=11,
        faultGeometrySourceDx=500.0,
        faultGeometrySourceName='the shipped surface')
    lib.validateFaultRoughGeometryForCase(par, verbose=False)   # dx == source
    par.faultGeometrySourceDx = 750.0        # source coarser than the case dx
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometryForCase(par, verbose=False)
    assert 'never refined' in str(e.value)


def test_defaultParameters_declares_the_source_spacing_slot():
    import defaultParameters
    par = defaultParameters.parameters()
    assert par.faultGeometrySourceDx is None
    assert isinstance(par.faultGeometrySourceName, str)


# ---- robustness properties of the supplied-geometry path -------------------
#
# One test per property, each written so that REMOVING the property makes it
# fail: idempotent re-setup, no clobbering of a deliberately placed file,
# provenance cross-check for a file copied in from elsewhere, and the
# too-fine-dx refusal. Bad inputs are perturbed good ones throughout.

def _caseWithWriter(tmp_path, dx=GEOM_DX, amp=100.0, calls=None):
    """A `par` stand-in whose faultGeometryWriter records every call."""
    y, dydx, dydz = _analyticSurface(amp, dx)

    def writer(out):
        if calls is not None:
            calls.append(out)
        lib.writeFaultRoughGeometry(y, dydx, dydz, dx, GEOM_FXMIN, GEOM_FZMIN,
                                    fname=out,
                                    provenance=dict(source='unit-test writer',
                                                    sourceDx=f'{dx:g}'))
    return types.SimpleNamespace(
        dx=dx, dz=dx, dy=dx,
        fxmin=GEOM_FXMIN, fxmax=GEOM_FXMAX,
        fzmin=GEOM_FZMIN, fzmax=GEOM_FZMAX,
        ymin=-20000.0, ymax=20000.0,
        nfx=round((GEOM_FXMAX - GEOM_FXMIN)/dx) + 1,
        nfz=round((GEOM_FZMAX - GEOM_FZMIN)/dx) + 1,
        faultGeometrySourceDx=dx,
        faultGeometrySourceName='the unit-test surface',
        faultGeometryWriter=writer)


def test_ensure_generates_the_file_when_it_is_missing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = []
    par = _caseWithWriter(tmp_path, calls=calls)
    d = lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    assert calls == [lib.FAULT_ROUGH_GEOMETRY_FILE]
    assert d['nnx'] == 21


def test_running_setup_twice_is_idempotent(tmp_path, monkeypatch):
    # Byte-for-byte identical file, and the writer is not called the second
    # time -- the property that makes case.setup safe to re-run.
    monkeypatch.chdir(tmp_path)
    calls = []
    par = _caseWithWriter(tmp_path, calls=calls)
    lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    first = (tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE).read_bytes()
    lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    assert (tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE).read_bytes() == first
    assert len(calls) == 1          # second call left the correct file alone


def test_a_correct_file_placed_by_the_case_author_is_not_clobbered(
        tmp_path, monkeypatch):
    # The 50 m TPV29 incident: a deliberately placed, CORRECT surface was
    # silently replaced by whatever the compset's default path produced.
    monkeypatch.chdir(tmp_path)
    calls = []
    par = _caseWithWriter(tmp_path, calls=calls)
    # A different but equally valid surface for the same grid.
    y, dydx, dydz = _analyticSurface(amp=137.0)
    lib.writeFaultRoughGeometry(y, dydx, dydz, GEOM_DX, GEOM_FXMIN, GEOM_FZMIN,
                                fname=lib.FAULT_ROUGH_GEOMETRY_FILE)
    mine = (tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE).read_bytes()
    lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    assert (tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE).read_bytes() == mine
    assert calls == []


def test_a_file_that_does_not_match_the_case_is_regenerated(tmp_path,
                                                            monkeypatch):
    # Copied in from a coarser case: wrong grid, so it must be replaced, not
    # kept -- the other half of the no-clobber rule.
    monkeypatch.chdir(tmp_path)
    calls = []
    par = _caseWithWriter(tmp_path, calls=calls)
    y, dydx, dydz = _analyticSurface(dx=1000.0)
    lib.writeFaultRoughGeometry(y, dydx, dydz, 1000.0, GEOM_FXMIN, GEOM_FZMIN,
                                fname=lib.FAULT_ROUGH_GEOMETRY_FILE)
    d = lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    assert calls == [lib.FAULT_ROUGH_GEOMETRY_FILE]
    assert d['nnx'] == 21 and d['dx'] == GEOM_DX


def test_a_mismatched_file_with_no_writer_is_refused_not_silently_used(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    par = _caseWithWriter(tmp_path)
    par.faultGeometryWriter = None
    y, dydx, dydz = _analyticSurface(dx=1000.0)
    lib.writeFaultRoughGeometry(y, dydx, dydz, 1000.0, GEOM_FXMIN, GEOM_FZMIN,
                                fname=lib.FAULT_ROUGH_GEOMETRY_FILE)
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    assert '1000.0' in str(e.value) and '500.0' in str(e.value)


def test_provenance_sidecar_round_trips_and_is_deterministic(tmp_path):
    y, dydx, dydz = _analyticSurface()
    path = str(tmp_path / 'g.txt')
    fields = dict(source='the unit-test surface', sourceDx='500')
    lib.writeFaultRoughGeometry(y, dydx, dydz, GEOM_DX, GEOM_FXMIN, GEOM_FZMIN,
                                fname=path, provenance=fields)
    side = lib.faultGeometryProvenancePath(path)
    first = open(side).read()
    prov = lib.readFaultGeometryProvenance(path)
    assert prov['source'] == 'the unit-test surface'
    assert float(prov['dx']) == GEOM_DX and int(prov['nnx']) == 21
    lib.writeFaultRoughGeometry(y, dydx, dydz, GEOM_DX, GEOM_FXMIN, GEOM_FZMIN,
                                fname=path, provenance=fields)
    assert open(side).read() == first        # no timestamps: byte-identical


def test_a_provenance_sidecar_that_disagrees_with_the_file_is_rejected(
        tmp_path):
    # "Copied the wrong file in": the surface and its provenance record no
    # longer describe the same thing.
    y, dydx, dydz = _analyticSurface()
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, dydx, dydz, GEOM_DX, GEOM_FXMIN, GEOM_FZMIN,
                                fname=path, provenance=dict(source='src'))
    lib.validateFaultRoughGeometry(path, **_geomKwargs())     # consistent: OK
    lib.writeFaultGeometryProvenance(path, dict(source='src', dx='250',
                                                nnx='41', nnz='21'))
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.validateFaultRoughGeometry(path, **_geomKwargs())
    msg = str(e.value)
    assert 'provenance' in msg and '250' in msg and '500.0' in msg


def test_validation_reports_the_provenance_it_found(tmp_path):
    y, dydx, dydz = _analyticSurface()
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, dydx, dydz, GEOM_DX, GEOM_FXMIN, GEOM_FZMIN,
                                fname=path,
                                provenance=dict(source='official 25 m data',
                                                sourceDx='25'))
    d = lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert d['provenance']['source'] == 'official 25 m data'


def test_no_sidecar_is_not_an_error(tmp_path):
    # Hand-written files have no provenance; that is allowed, it is advisory.
    path = _writeGoodGeometry(tmp_path)
    d = lib.validateFaultRoughGeometry(path, **_geomKwargs())
    assert d['provenance'] is None


def test_ensure_refuses_a_dx_the_source_cannot_reach_before_touching_the_file(
        tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    calls = []
    par = _caseWithWriter(tmp_path, calls=calls)
    par.dx = par.dz = par.dy = 250.0                 # finer than the 500 m source
    par.faultGeometrySourceAvailableDx = [500.0, 1000.0]
    with pytest.raises(lib.FaultGeometryError) as e:
        lib.ensureFaultRoughGeometryForCase(par, verbose=False)
    msg = str(e.value)
    assert '250.0' in msg and '500.0' in msg
    assert 'Sources available here: 500 m, 1000 m.' in msg
    assert calls == []                               # nothing was written
    assert not (tmp_path / lib.FAULT_ROUGH_GEOMETRY_FILE).exists()


def test_provenance_keeps_full_precision_for_a_dipping_fault_corner(tmp_path):
    # Regression: the sidecar wrote its numbers with '%g' (6 significant
    # digits). A dipping fault's fzmin is fzmin*sin(dip) -- test.tpv10's is
    # -12990.381056766579 -- which '%g' rounds to -12990.4, 2 cm away from the
    # header. The cross-check then fired on the writer's OWN rounding and
    # case.setup refused a perfectly good generated file.
    fzmin = -12990.381056766579
    # 6x6 nodes, not 4x3: the interior derivative check needs at least 5 along
    # each axis for the y''' stencil and now REFUSES a grid smaller than that
    # rather than falling back to a weaker bound (rule 2). This test is about
    # the provenance sidecar's precision, so it just needs a grid the validator
    # can actually clear.
    y = np.zeros((6, 6))
    path = str(tmp_path / 'g.txt')
    lib.writeFaultRoughGeometry(y, y, y, 500.0, -15000.0, fzmin, fname=path,
                                provenance=dict(source='dipping-fault test'))
    prov = lib.readFaultGeometryProvenance(path)
    assert abs(float(prov['fzmin']) - fzmin) <= lib.FAULT_GEOM_COORD_TOL
    lib.validateFaultRoughGeometry(
        path, dx=500.0, dz=500.0, dy=500.0, fxmin=-15000.0, fxmax=-12500.0,
        fzmin=fzmin, fzmax=fzmin + 2500.0, verbose=False)


def test_every_full_spec_entry_has_the_keys_its_consumer_reads():
    """FULL_SPECS entries must carry exactly the keys run_e2e_full.py reads.

    Guards a real defect: the test.tpv29 entry was written with decomp=(4,1,4)
    and source=... while the consumer reads spec['nx'], spec['ny'], spec['nz']
    and spec['citation']. Being first in the dict, it made the entire full tier
    die with KeyError before any case launched -- and nothing caught it, because
    the fast e2e tier never imports full_specs.py at all.
    """
    import re
    import sys
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    specsDir = os.path.join(root, 'testsys', 'e2e')
    sys.path.insert(0, specsDir)
    try:
        from full_specs import FULL_SPECS
    finally:
        sys.path.remove(specsDir)

    # The required keys are read out of the consumer, not hard-coded here, so
    # this test follows run_e2e_full.py if it starts reading something new.
    consumer = open(os.path.join(specsDir, 'run_e2e_full.py'),
                    errors='replace').read()
    required = set(re.findall(r"spec\[['\"](\w+)['\"]\]", consumer))
    assert required, 'could not find any spec[...] reads in run_e2e_full.py'

    for case, spec in FULL_SPECS.items():
        missing = sorted(required - set(spec))
        assert not missing, (
            '%s is missing %s; run_e2e_full.py reads those keys directly '
            'and will raise KeyError' % (case, missing))
