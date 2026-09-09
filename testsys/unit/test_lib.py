"""Unit tests for scripts/lib.py (fault-loading and boxcar-taper helpers).

Every expected value below is computed independently of the function under
test -- either an exact closed form (asinh(1) = ln(1+sqrt(2))) or a
documented boundary/symmetry property of the boxcar tapers -- never a copy
of the function's own output.
"""
import math

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
