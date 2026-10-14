"""Unit tests for the rate-and-state normal-stress inputs (issues #189, #200)
and the non-compressive stop (#201).

  - lib.resolveNormalStressInputs: the defaults and the refusals case.setup
    applies before it writes bGlobal.txt;
  - defaultParameters: the default values, which must reproduce the
    behaviour before the inputs existed (L = Dc, caps iff insertFaultType>0
    and C_elastic==1, the old -10/-40 MPa constants);
  - faulting.normal_stress_caps_on: the tri-state switch;
  - fric.rate_state_normal_stress: L_pc <= 0 is bitwise the Dc form,
    L_pc > 0 is the closed form with that L;
  - driver.latch_normal_stress / latched_normal_stress_error: the latch keeps
    the first event and the message names time, node, coordinates and value.

The end-to-end parse (bGlobal.txt -> both readers, stale and invalid lines)
is testsys/regression/test_nstress_sign_convention.py; the end-to-end stop
is testsys/regression/test_rsf_noncompressive_stop.py.
"""
import types

import numpy as np
import pytest

import defaultParameters as dp
import lib
from eqdyna import driver, faulting, fric
from eqdyna import globalvar as gv
from eqdyna.checkInputConsistency import InputConsistencyError


def _par(**kw):
    return types.SimpleNamespace(**kw)


# --- defaults --------------------------------------------------------------

def test_resolve_defaults_on_an_empty_namespace():
    assert lib.resolveNormalStressInputs(_par()) == (0.0, -1, -10.0e6, -40.0e6)


def test_default_parameters_equal_the_resolver_defaults():
    p = dp.parameters
    assert (p.fric_pc_L, p.C_normal_stress_caps, p.min_norm, p.max_norm) == \
        lib.resolveNormalStressInputs(_par())


def test_defaults_equal_the_old_python_constants():
    # The constants the caps used before #200 (globalvar.f90's old
    # parameters, mirrored in the Python globalvar).
    assert dp.parameters.min_norm == gv.MIN_NORM == -10.0e6
    assert dp.parameters.max_norm == gv.MAX_NORM == -40.0e6


def test_resolve_passes_valid_values_through():
    got = lib.resolveNormalStressInputs(
        _par(fric_pc_L=0.25, C_normal_stress_caps=1, min_norm=-5e6, max_norm=-5e6))
    assert got == (0.25, 1, -5e6, -5e6)
    assert isinstance(got[1], int)


# --- refusals --------------------------------------------------------------

@pytest.mark.parametrize('kw', [
    dict(fric_pc_L=float('nan')),
    dict(C_normal_stress_caps=2),
    dict(C_normal_stress_caps=-2),
    dict(C_normal_stress_caps=0.5),
    dict(min_norm=-40e6, max_norm=-10e6),      # min < max
    dict(min_norm=0.0),                        # not compressive
    dict(min_norm=1e6, max_norm=-40e6),
])
def test_resolve_refuses_invalid(kw):
    with pytest.raises(ValueError):
        lib.resolveNormalStressInputs(_par(**kw))


# --- caps switch -----------------------------------------------------------

@pytest.mark.parametrize('sw,ift,ce,on', [
    (-1, 0, 1, False), (-1, 1, 1, True), (-1, 2, 1, True), (-1, 1, 0, False),
    (-1, 0, 0, False),
    (0, 0, 1, False), (0, 1, 1, False), (0, 2, 0, False),
    (1, 0, 1, True), (1, 0, 0, True), (1, 2, 0, True),
])
def test_caps_switch_truth_table(sw, ift, ce, on):
    finv = dict(C_normal_stress_caps=sw, insertFaultType=ift, C_elastic=ce)
    assert faulting.normal_stress_caps_on(finv) is on


def test_caps_switch_default_when_absent():
    assert faulting.normal_stress_caps_on(dict(insertFaultType=1, C_elastic=1))
    assert not faulting.normal_stress_caps_on(dict(insertFaultType=0, C_elastic=1))


# --- theta_pc relaxation length --------------------------------------------

def _fric(n, dc):
    f = np.zeros((n, 100))
    f[:, gv.RSF_DC] = dc
    return f


def test_L_default_is_bitwise_the_dc_form():
    rng = np.random.default_rng(189)
    n = 64
    dc = rng.uniform(0.1, 0.6, n)
    V2, th, tn = rng.uniform(0, 3, n), rng.uniform(5e7, 2e8, n), -rng.uniform(5e7, 2e8, n)
    dt = 0.0125
    f = _fric(n, dc)
    ref_dot = -V2 / dc * (th - np.abs(tn))
    ref = th + ref_dot * dt
    for L_pc in (0.0, -1.0):
        got, got_dot = fric.rate_state_normal_stress(np, V2, th, tn, f, dt, L_pc)
        assert np.array_equal(got, ref) and np.array_equal(got_dot, ref_dot)
    got, _ = fric.rate_state_normal_stress(np, V2, th, tn, f, dt)
    assert np.array_equal(got, ref)


def test_L_positive_uses_that_length_for_every_node():
    n = 5
    f = _fric(n, np.linspace(0.1, 0.5, n))
    V2 = np.full(n, 2.0); th = np.full(n, 1.0e8); tn = np.full(n, -1.2e8); dt = 0.01
    got, dot = fric.rate_state_normal_stress(np, V2, th, tn, f, dt, 0.4)
    # -V/L*(theta - |tn|) = -2/0.4*(1e8 - 1.2e8) = +1e8 per s
    assert np.allclose(dot, 1.0e8, rtol=0, atol=1e-6)
    assert np.allclose(got, 1.0e8 + 1.0e8 * dt, rtol=0, atol=1e-6)


# --- the non-compressive latch ---------------------------------------------

def test_latch_keeps_the_first_event_and_lowest_row():
    z = np.zeros(driver.NSVIOL_LEN)
    ok = np.array([-1e6, -2e6, -3e6, -4e6])
    assert np.array_equal(driver.latch_normal_stress(np, z, ok, 0.1), z)
    # the boundary value 0 counts (Fortran stops at >= 0)
    s1 = driver.latch_normal_stress(np, z, np.array([-1e6, 0.0, 2e6, -4e6]), 0.2)
    assert s1.tolist() == [1.0, 0.2, 1.0, 0.0]
    # a later step with another bad row does not overwrite it
    s2 = driver.latch_normal_stress(np, s1, np.array([3e6, -2e6, -3e6, -4e6]), 0.3)
    assert np.array_equal(s2, s1)


def _fake_S():
    # two faults, rows 0-2 on fault 0 and rows 3-4 on fault 1
    meshCoor = np.arange(30, dtype=float).reshape(10, 3) * 100.0
    return dict(fault_of=np.array([0, 0, 0, 1, 1]),
                nsmp1=np.array([2, 4, 6, 7, 9]),
                meshCoor=meshCoor, friclaw=5)


def test_latched_error_none_when_nothing_latched():
    assert driver.latched_normal_stress_error(_fake_S(), np.zeros(4)) is None


def test_latched_error_names_time_node_and_value():
    err = driver.latched_normal_stress_error(
        _fake_S(), np.array([1.0, 0.75, 4.0, 3.5e6]), rank=2)
    assert isinstance(err, InputConsistencyError)
    assert err.code == 65
    msg = str(err)
    assert 'rank 2, friclaw 5' in msg
    assert 'time 7.50000000e-01 s' in msg
    assert 'fault 2, fault node 2' in msg          # row 4 is fault 1's 2nd row
    x, y, z = _fake_S()['meshCoor'][9]
    assert 'node x, y, z (m) %.8e %.8e %.8e' % (x, y, z) in msg
    assert 'normal traction 3.50000000e+00 MPa' in msg
