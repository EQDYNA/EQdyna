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
