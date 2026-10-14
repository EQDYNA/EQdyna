"""Unit tests for the rate-and-state normal-stress inputs (issues #189, #200)
and the non-compressive stop (#201).

  - lib.resolveNormalStressInputs: the defaults and the refusals case.setup
    applies before it writes bGlobal.txt;
  - defaultParameters: the default values, which must reproduce the
    behaviour before the inputs existed (L = Dc, caps iff insertFaultType>0
    and C_elastic==1, the old -10/-40 MPa constants);
  - faulting.normal_stress_caps_on: the tri-state switch;
  - fric.rate_state_normal_stress, the three fric_pc_L modes (EQquasi's
    meaning): unset (None) is bitwise the Dc form; <= 0 is no relaxation,
    theta_pc = |tnrm| exactly with no division (L = 0 included); > 0 is the
    update with that L;
  - driver.latch_normal_stress / latched_normal_stress_error: the latch keeps
    the first event and the message names time, node, coordinates and value.

The end-to-end parse (bGlobal.txt -> both readers, stale and invalid lines)
is testsys/regression/test_nstress_sign_convention.py; the end-to-end stop
is testsys/regression/test_rsf_noncompressive_stop.py; the three modes end
to end in both codes are testsys/regression/test_fric_pc_L_modes.py.
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
    assert lib.resolveNormalStressInputs(_par()) == (None, -1, -10.0e6, -40.0e6)


def test_default_fric_pc_L_is_unset():
    assert dp.parameters.fric_pc_L is None


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


@pytest.mark.parametrize('L,want', [(None, None), (0, 0.0), (0.0, 0.0),
                                    (-1, -1.0), (-0.5, -0.5), (0.004, 0.004),
                                    (2, 2.0), (np.float64(0.1), 0.1)])
def test_resolve_fric_pc_L_modes(L, want):
    got = lib.resolveNormalStressInputs(_par(fric_pc_L=L))[0]
    assert got == want
    assert got is None if want is None else type(got) is float


# --- refusals --------------------------------------------------------------

@pytest.mark.parametrize('kw', [
    dict(fric_pc_L=float('nan')),
    dict(fric_pc_L=float('inf')),
    dict(fric_pc_L=-float('inf')),
    dict(fric_pc_L='nan'),
    dict(fric_pc_L='x'),
    dict(fric_pc_L=[0.1]),
    dict(fric_pc_L=True),
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


def _state(n, seed):
    rng = np.random.default_rng(seed)
    dc = rng.uniform(0.1, 0.6, n)
    V2, th, tn = rng.uniform(0, 3, n), rng.uniform(5e7, 2e8, n), -rng.uniform(5e7, 2e8, n)
    return dc, V2, th, tn


def test_L_unset_is_bitwise_the_dc_form():
    n = 64
    dc, V2, th, tn = _state(n, 189)
    dt = 0.0125
    f = _fric(n, dc)
    ref_dot = -V2 / dc * (th - np.abs(tn))
    ref = th + ref_dot * dt
    for args in ((None,), ()):          # explicit None and the default
        got, got_dot = fric.rate_state_normal_stress(np, V2, th, tn, f, dt, *args)
        assert np.array_equal(got, ref) and np.array_equal(got_dot, ref_dot)


def test_L_unset_equals_L_given_as_a_uniform_dc():
    # unset reads each node's Dc; giving that same Dc as the number is the
    # same arithmetic, so the same bits.
    n = 16
    _, V2, th, tn = _state(n, 7)
    f = _fric(n, np.full(n, 0.4))
    a = fric.rate_state_normal_stress(np, V2, th, tn, f, 0.01, None)
    b = fric.rate_state_normal_stress(np, V2, th, tn, f, 0.01, 0.4)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


@pytest.mark.parametrize('L_pc', [0.0, -0.0, -1.0, -1e30])
def test_L_nonpositive_is_no_relaxation(L_pc):
    # EQquasi's meaning: theta_pc = |tnrm| at once, whatever theta_pc and the
    # slip rate were; a branch, so L = 0 never divides (no warning raised).
    n = 32
    dc, V2, th, tn = _state(n, 200)
    V2[0] = 0.0                          # V = 0 and L = 0: would be 0/0
    dt = 0.0125
    f = _fric(n, dc)
    with np.errstate(all='raise'):
        got, got_dot = fric.rate_state_normal_stress(np, V2, th, tn, f, dt, L_pc)
    assert np.array_equal(got, np.abs(tn))
    assert np.all(np.isfinite(got_dot))
    assert np.array_equal(got_dot, (np.abs(tn) - th) / dt)


def test_L_nonpositive_does_not_depend_on_dc_or_slip_rate():
    n = 8
    _, V2, th, tn = _state(n, 5)
    a, _ = fric.rate_state_normal_stress(np, V2, th, tn, _fric(n, np.full(n, 0.1)), 0.01, 0.0)
    b, _ = fric.rate_state_normal_stress(np, 3 * V2, th, tn, _fric(n, np.full(n, 0.9)), 0.01, -2.0)
    assert np.array_equal(a, b)


def test_L_positive_uses_that_length_for_every_node():
    n = 5
    f = _fric(n, np.linspace(0.1, 0.5, n))
    V2 = np.full(n, 2.0); th = np.full(n, 1.0e8); tn = np.full(n, -1.2e8); dt = 0.01
    got, dot = fric.rate_state_normal_stress(np, V2, th, tn, f, dt, 0.4)
    # -V/L*(theta - |tn|) = -2/0.4*(1e8 - 1.2e8) = +1e8 per s
    assert np.allclose(dot, 1.0e8, rtol=0, atol=1e-6)
    assert np.allclose(got, 1.0e8 + 1.0e8 * dt, rtol=0, atol=1e-6)


def test_L_modes_under_jax_match_numpy():
    jnp = pytest.importorskip('jax.numpy')
    import jax
    jax.config.update('jax_enable_x64', True)
    n = 16
    dc, V2, th, tn = _state(n, 31)
    f = _fric(n, dc)
    for L_pc in (None, 0.0, -1.0, 0.25):
        a = fric.rate_state_normal_stress(np, V2, th, tn, f, 0.01, L_pc)
        b = fric.rate_state_normal_stress(jnp, jnp.asarray(V2), jnp.asarray(th),
                                          jnp.asarray(tn), jnp.asarray(f), 0.01, L_pc)
        assert np.allclose(a[0], np.asarray(b[0]), rtol=1e-15, atol=0)
        assert np.allclose(a[1], np.asarray(b[1]), rtol=1e-13, atol=1e-9)


# --- the non-compressive latch ---------------------------------------------

def test_latch_keeps_the_first_event_and_lowest_row():
    z = np.zeros(driver.NSVIOL_LEN)
    ok = np.array([-1e6, -2e6, -3e6, -4e6])
    assert np.array_equal(driver.latch_normal_stress(np, z, ok, 0.1, 1), z)
    # the boundary value 0 counts (Fortran stops at >= 0)
    s1 = driver.latch_normal_stress(np, z, np.array([-1e6, 0.0, 2e6, -4e6]), 0.2, 2)
    assert s1.tolist() == [1.0, 0.2, 2.0, 1.0, 0.0]
    # a later step with another bad row does not overwrite it
    s2 = driver.latch_normal_stress(np, s1, np.array([3e6, -2e6, -3e6, -4e6]), 0.3, 3)
    assert np.array_equal(s2, s1)


def _fake_S():
    # two faults, rows 0-2 on fault 0 and rows 3-4 on fault 1
    meshCoor = np.arange(30, dtype=float).reshape(10, 3) * 100.0
    return dict(fault_of=np.array([0, 0, 0, 1, 1]),
                nsmp1=np.array([2, 4, 6, 7, 9]),
                meshCoor=meshCoor, friclaw=5)


def test_latched_error_none_when_nothing_latched():
    assert driver.latched_normal_stress_error(_fake_S(), np.zeros(driver.NSVIOL_LEN)) is None


def test_latched_error_names_time_node_and_value():
    err = driver.latched_normal_stress_error(
        _fake_S(), np.array([1.0, 0.75, 6.0, 4.0, 3.5e6]), rank=2)
    assert isinstance(err, InputConsistencyError)
    assert err.code == 65
    msg = str(err)
    assert 'rank 2, friclaw 5' in msg
    assert 'time 7.50000000e-01 s, step 6,' in msg
    assert 'fault 2, fault node 2' in msg          # row 4 is fault 1's 2nd row
    x, y, z = _fake_S()['meshCoor'][9]
    assert 'node x, y, z (m) %.8e %.8e %.8e' % (x, y, z) in msg
    assert 'normal traction 3.50000000e+00 MPa' in msg


def test_latched_error_refuses_S_without_fault_of():
    S = _fake_S()
    del S['fault_of']
    with pytest.raises(ValueError, match='fault_of'):
        driver.latched_normal_stress_error(
            S, np.array([1.0, 0.75, 6.0, 4.0, 3.5e6]))
