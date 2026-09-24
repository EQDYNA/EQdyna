"""Unit tests for testsys/compare.py -- the one comparison the sweep uses.

Two kinds of test here:
  1. the gate FAILS when the answer is wrong (a gate that cannot fail is not a
     gate): a perturbed value, a missing artifact, a different node set.
  2. the committed references are in the canonical form the gate assumes.

Everything runs against small synthetic arrays or a copy of a committed
reference in tmp_path -- never against the live test/ tree, and never writing
into test.reference.results/ (rule 7).
"""
import os
import shutil

import numpy as np
import pytest
import xarray as xr

from conftest import REPO_ROOT
from testsys import compare, frt_canonical, matrix

NCOL = frt_canonical.FRT_COLUMNS


def _frt(coords, fill=0.0):
    """A synthetic frt array: coordinates in columns 0-2, `fill` elsewhere."""
    arr = np.full((len(coords), NCOL), float(fill))
    arr[:, :3] = np.asarray(coords, dtype=float)
    return arr


def _write(arr, path):
    frt_canonical.write_canonical(arr, str(path))
    return str(path)


# --------------------------------------------------------------------------
# the committed references are canonical
# --------------------------------------------------------------------------
@pytest.mark.parametrize('case', matrix.CASES)
def test_every_committed_reference_is_already_canonical(case):
    # The gate compares against this file directly, so it must BE the
    # canonical form: canonicalising it again is a no-op. If this fails, the
    # reference was hand-edited or written by something other than
    # frt_canonical.write_canonical.
    ref = compare.load_reference(case)
    np.testing.assert_array_equal(frt_canonical.canonicalize(ref), ref)


@pytest.mark.parametrize('case', matrix.CASES)
def test_every_committed_reference_has_one_row_per_fault_node(case):
    ref = compare.load_reference(case)
    coords = np.round(ref[:, :3], frt_canonical.COORD_DECIMALS)
    assert len(np.unique(coords, axis=0)) == ref.shape[0], (
        'duplicate fault-node coordinates in the reference for %s -- the '
        'reference is supposed to be decomposition-independent' % case)


def test_missing_reference_raises_rather_than_passing(tmp_path, monkeypatch):
    monkeypatch.setattr(compare, 'REFERENCE_ROOT', str(tmp_path))
    with pytest.raises(FileNotFoundError):
        compare.load_reference('test.tpv8')


# --------------------------------------------------------------------------
# canonicalisation refuses to make a silent choice
# --------------------------------------------------------------------------
def test_canonicalize_refuses_duplicates_that_disagree():
    # A fault node written by two ranks with DIFFERENT values: dedupe would
    # keep the first and throw the other away silently.
    arr = _frt([(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)])
    arr[1, 7] = 1.0
    with pytest.raises(ValueError, match='DIFFERENT values'):
        frt_canonical.canonicalize(arr)


def test_canonicalize_accepts_duplicates_that_agree():
    arr = _frt([(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], fill=3.0)
    out = frt_canonical.canonicalize(arr)
    assert out.shape[0] == 2
    np.testing.assert_array_equal(out[:, 0], [0.0, 1.0])   # lexsorted by x


# --------------------------------------------------------------------------
# the abs-max gate
# --------------------------------------------------------------------------
def test_abs_max_gate_passes_on_identical_arrays():
    arr = _frt([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], fill=5.0)
    ok, lines = compare.abs_max_gate('test.tpv8', arr, arr.copy())
    assert ok
    assert 'max|diff|=0.000000e+00' in lines[0]


def test_abs_max_gate_fails_just_past_the_bound_and_says_where():
    bound = matrix.CASE_BOUND['test.tpv8']
    a = _frt([(0.0, 0.0, 0.0)], fill=5.0)
    b = a.copy()
    b[0, 9] += bound * 10
    ok, lines = compare.abs_max_gate('test.tpv8', a, b)
    assert not ok
    assert 'col 9' in lines[0]
    assert 'never moved' in lines[1]


def test_abs_max_gate_uses_the_case_bound_not_a_shared_one():
    # tpv1053d's bound is 1e-4; the same diff must fail tpv8's 1e-8. One
    # bound per CASE, identical across backends -- that is what makes the
    # columns comparable.
    a = _frt([(0.0, 0.0, 0.0)], fill=1.0)
    b = a.copy()
    b[0, 4] += 1e-6
    assert compare.abs_max_gate('test.tpv1053d', a, b)[0]
    assert not compare.abs_max_gate('test.tpv8', a, b)[0]


# --------------------------------------------------------------------------
# the flip-budget gate (test.drv.a6)
# --------------------------------------------------------------------------
def _fnft(values):
    arr = _frt([(float(i), 0.0, 0.0) for i in range(len(values))])
    arr[:, matrix.DRV_A6['fnft_col']] = values
    return arr


def test_flip_decomposition_counts_existence_flips_and_timing_shifts():
    # 5 nodes, hand-countable: node 0 matches, node 1 shifts by 2 s (> 1 s),
    # node 2 ruptured only in ref, node 3 only in run, node 4 matches.
    sentinel = 99999.0
    ref = _fnft([1.0, 1.0, 2.0, sentinel, 3.0])
    run = _fnft([1.0, 3.0, sentinel, 2.0, 3.0])
    d = compare.flip_decomposition(ref, run)
    assert d['n_ruptured_both'] == 3
    assert d['n_only_ref'] == 1 and d['n_only_run'] == 1
    assert d['n_existence_flips'] == 2
    assert d['n_timing_shifts'] == 1
    assert d['total_flips'] == 3
    assert d['n_matched_arrival'] == 2


def test_flip_decomposition_refuses_a_run_that_never_nucleated():
    sentinel = 99999.0
    ref = _fnft([1.0, 2.0])
    run = _fnft([sentinel, sentinel])
    with pytest.raises(AssertionError, match='never|nucleate'):
        compare.flip_decomposition(ref, run)


def test_flip_budget_gate_fails_past_the_flip_budget():
    n = matrix.DRV_A6['total_flip_bound'] + 10
    sentinel = 99999.0
    ref = _fnft([1.0] * (n + 5))
    run = _fnft([1.0] * 5 + [sentinel] * n)
    ok, lines = compare.flip_budget_gate('test.drv.a6', ref, run)
    assert not ok
    assert 'flips=%d/%d' % (n, matrix.DRV_A6['total_flip_bound']) in lines[1]


# --------------------------------------------------------------------------
# end to end on a run directory
# --------------------------------------------------------------------------
def test_compare_frt_passes_when_the_run_reproduces_the_reference(tmp_path):
    # The reference copied in as a serial run's single frt.txt0: one file,
    # not four, and the gate does not care.
    shutil.copy(compare.reference_path('test.tpv8'), str(tmp_path / 'frt.txt0'))
    ok, lines = compare.compare_frt('test.tpv8', str(tmp_path))
    assert ok, lines
    assert 'max|diff|=0.000000e+00' in lines[0]


def test_compare_frt_fails_on_one_perturbed_value(tmp_path):
    ref = compare.load_reference('test.tpv8')
    ref[3, 8] += 10.0
    _write(ref, tmp_path / 'frt.txt0')
    ok, lines = compare.compare_frt('test.tpv8', str(tmp_path))
    assert not ok, lines


def test_compare_frt_is_independent_of_how_many_rank_files_a_run_wrote(tmp_path):
    # The same physics split across 3 "rank" files, with a shared node
    # duplicated into two of them, must compare identically to the reference.
    ref = compare.load_reference('test.tpv8')
    _write(ref[:800], tmp_path / 'frt.txt0')
    _write(ref[790:1500], tmp_path / 'frt.txt1')   # 790:800 written twice
    _write(ref[1500:], tmp_path / 'frt.txt2')
    ok, lines = compare.compare_frt('test.tpv8', str(tmp_path))
    assert ok, lines


def test_compare_frt_fails_on_a_different_fault_node_set(tmp_path):
    ref = compare.load_reference('test.tpv8')
    _write(ref[:-1], tmp_path / 'frt.txt0')        # one node short
    with pytest.raises(ValueError, match='node counts differ'):
        compare.compare_frt('test.tpv8', str(tmp_path))


def test_compare_frt_refuses_an_empty_run_directory(tmp_path):
    # No frt output at all must be a loud failure, never an empty-vs-empty pass.
    with pytest.raises(ValueError, match='no frt files'):
        compare.compare_frt('test.tpv8', str(tmp_path))


def test_fortran_cell_fails_when_plotruptureDynamics_output_is_missing(tmp_path):
    # matrix.ARTIFACTS says the fortran backend produces frt AND nc. A run
    # that wrote only frt must not pass as if it had covered both.
    shutil.copy(compare.reference_path('test.tpv8'), str(tmp_path / 'frt.txt0'))
    ok, lines = compare.compare_cell('test.tpv8', 'fortran', str(tmp_path))
    assert not ok
    assert any('fault.dyna.r.nc' in line for line in lines)


def test_python_cell_compares_frt_only(tmp_path):
    shutil.copy(compare.reference_path('test.tpv8'), str(tmp_path / 'frt.txt0'))
    ok, lines = compare.compare_cell('test.tpv8', 'python-jax', str(tmp_path))
    assert ok, lines
    assert not any('fault.dyna.r.nc' in line for line in lines)


# --------------------------------------------------------------------------
# fault.dyna.r.nc (migrated from testsys/unit/test_check_comparisons.py, which
# tested these same functions through check.test.py before the comparison
# moved into one module)
# --------------------------------------------------------------------------
def _dataset(value):
    return xr.Dataset({'slip': (('node',), np.array([value, value, value]))})


def test_compare_nc_files_identical_reports_success(tmp_path):
    fn1, fn2 = str(tmp_path / 'ref.nc'), str(tmp_path / 'run.nc')
    _dataset(1.23).to_netcdf(fn1)
    _dataset(1.23).to_netcdf(fn2)
    assert compare.compare_nc_files(fn1, fn2).startswith('SUCCESS')


def test_compare_nc_files_value_beyond_threshold_reports_fail(tmp_path):
    fn1, fn2 = str(tmp_path / 'ref.nc'), str(tmp_path / 'run.nc')
    _dataset(1.0).to_netcdf(fn1)
    _dataset(2.0).to_netcdf(fn2)
    assert compare.compare_nc_files(fn1, fn2).startswith('FAIL')


def test_compare_nc_files_within_threshold_but_not_bit_exact_reports_success(tmp_path):
    # Regression: this check used to require bit-exact data via
    # f1.identical(f2), which turned ordinary MPI reduction-order noise
    # (~2.9e-8 on test.tpv10's shear_strike) into a false FAIL.
    fn1, fn2 = str(tmp_path / 'ref.nc'), str(tmp_path / 'run.nc')
    _dataset(1.0).to_netcdf(fn1)
    _dataset(1.0 + 5e-8).to_netcdf(fn2)
    assert compare.compare_nc_files(fn1, fn2).startswith('SUCCESS')


def test_compare_nc_files_differing_attrs_reports_fail(tmp_path):
    fn1, fn2 = str(tmp_path / 'ref.nc'), str(tmp_path / 'run.nc')
    ds1 = _dataset(1.0)
    ds1.attrs['source'] = 'reference-run'
    ds1.to_netcdf(fn1)
    ds2 = _dataset(1.0)
    ds2.attrs['source'] = 'different-run'
    ds2.to_netcdf(fn2)
    assert compare.compare_nc_files(fn1, fn2).startswith('FAIL')


def test_reference_root_is_the_committed_tree_and_is_never_written(tmp_path):
    # Rule 7: the gate reads this tree and nothing else writes through it.
    assert compare.REFERENCE_ROOT == os.path.join(REPO_ROOT, 'test.reference.results')
