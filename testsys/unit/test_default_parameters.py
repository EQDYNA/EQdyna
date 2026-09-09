"""Unit tests for scripts/defaultParameters.py invariants.

defaultParameters.parameters is a class whose body runs the fault-grid and
on-fault-array setup at import time -- so "importing it" already exercises
the behaviour we want to check. Every expected value here is either an
integer computed by hand from the documented default geometry, or a ratio
taken verbatim from the source comment ("initial shear stress = -0.41 *
initial normal stress") -- not a copy of whatever the array happens to
contain.
"""
import numpy as np

import defaultParameters as dp
from defaultParameters import mod4dip


def test_fault_grid_node_counts_match_hand_computed_defaults():
    # nfx = (fxmax-fxmin)/dx + 1 = (22e3 - (-22e3))/500 + 1 = 89
    # nfz = (fzmax-fzmin)/dz + 1 = (0 - (-22e3))/500 + 1   = 45
    assert dp.parameters.nfx == 89
    assert dp.parameters.nfz == 45


def test_fault_coordinate_arrays_span_the_documented_extent():
    p = dp.parameters
    assert len(p.fx) == p.nfx
    assert len(p.fz) == p.nfz
    assert p.fx[0] == p.fxmin
    assert p.fx[-1] == p.fxmax
    assert p.fz[0] == p.fzmin
    assert p.fz[-1] == p.fzmax


def test_on_fault_vars_shape_matches_grid_and_100_variable_slots():
    p = dp.parameters
    assert p.on_fault_vars.shape == (p.nfz, p.nfx, 100)


def test_initial_shear_stress_is_minus_0p41_times_normal_stress_everywhere():
    # Documented default: on_fault_vars[...,8] = -0.41 * on_fault_vars[...,7]
    # (fault_sw_fs-like static friction assumption baked into the default
    # case). Checked over the whole grid, not a single cell, so an indexing
    # bug in the nested ix/iz loop would also be caught.
    p = dp.parameters
    shear = p.on_fault_vars[:, :, 8]
    normal = p.on_fault_vars[:, :, 7]
    assert np.allclose(shear, -0.41 * normal)


def test_mod4dip_vertical_fault_is_the_asymptotic_no_dip_limit():
    # dip=90 degrees is a vertical strike-slip fault: the dip rotation
    # should be a no-op on the along-dip cell size and leave the source
    # on the y=0 plane (asymptotic-limit check, PROJECT_RULES.md-style).
    dx, fzmin, srcDowndip, vp = 500.0, -22.0e3, 7.5e3, 6.0e3
    dz, fzmin_out, ysource, zsource, dt = mod4dip(90.0, dx, fzmin, srcDowndip, vp)
    assert dz == dx
    assert abs(ysource) < 1e-9  # cos(90 deg) == 0, up to float error
    assert zsource == srcDowndip
    assert fzmin_out == fzmin
    assert dt == 0.5 * dz / vp


def test_mod4dip_45_degrees_matches_hand_computed_trig():
    import math
    dx, fzmin, srcDowndip, vp = 500.0, -22.0e3, 7.5e3, 6.0e3
    dz, fzmin_out, ysource, zsource, dt = mod4dip(45.0, dx, fzmin, srcDowndip, vp)
    s = math.sqrt(2) / 2  # sin(45) == cos(45) == sqrt(2)/2, exact
    assert math.isclose(dz, dx * s, rel_tol=1e-12)
    assert math.isclose(ysource, -srcDowndip * s, rel_tol=1e-12)
    assert math.isclose(zsource, srcDowndip * s, rel_tol=1e-12)
    assert math.isclose(fzmin_out, fzmin * s, rel_tol=1e-12)
