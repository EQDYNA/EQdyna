#! /usr/bin/env python3
"""Unit tests for testsys/frt_canonical.py -- the decomposition-independent
canonical form for frt output.

What these pin down, and why each matters:

  * ORDER INDEPENDENCE. The canonical array must not depend on how many ranks
    wrote the input or in what order the files were concatenated. That is the
    whole property the canonical form exists to provide, so it is asserted by
    reversing the file order and demanding a bit-identical result -- not
    allclose.

  * THE DEDUPE IS REAL, NOT COSMETIC. Partition-boundary nodes are written by
    every rank that owns them. The committed references duplicate 1.0-3.6% of
    rows, so a test that only checked "output has some rows" would pass on a
    no-op dedupe. These assert the exact post-dedupe row counts.

  * FAILURES ARE LOUD. An empty file list, a wrong column count, and a node-set
    mismatch each raise. An empty array compares equal to another empty array,
    so "no files found" must never be able to look like a pass (rule 2).

  * ALIGNMENT ACROSS RANK COUNTS. A serial 1-file output and a 2-rank
    reference must reduce to the same node set with coordinates agreeing
    exactly. This is the property that lets every backend be compared against
    one artifact.

Cheap (rule 9): reads committed reference files, no build, no simulation.
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'testsys'))

import frt_canonical as F  # noqa: E402

REF = os.path.join(REPO, 'test.reference.results')

# Node counts of the committed canonical references, measured when they were
# generated from the per-rank output of the runs that produced them. They pin
# the fixtures as well as the dedupe: a canonical reference that changed node
# count is a different fault discretisation, not a tolerance question.
EXPECTED_ROWS = {
    'test.tpv8': 1891,
    'test.tpv10': 1891,
    'test.tpv104': 2701,
    'test.tpv1053d': 4005,
    'test.drv.a6': 5151,
    'test.tpv29': 3321,
    'test.meng2023a': 651,
    'test.meng2023cb': 651,
    'test.tpv36': 3477,
}

# How the synthetic multi-rank fixtures below are built. OVERLAP_ROWS is the
# number of rows deliberately written into two "rank" files, standing in for
# the partition-boundary nodes every rank that owns them writes.
RANKS = 3
OVERLAP_ROWS = 7


def _canonical(name):
    """The committed canonical reference for a case. A missing reference is a
    FAILURE, not a skip: "I could not check this" and "this is fine" are
    different answers (rule 2)."""
    p = os.path.join(REF, name, 'frt.canonical.txt')
    if not os.path.isfile(p):
        raise AssertionError('missing canonical reference %s -- the reference '
                             'tree is incomplete' % p)
    return np.loadtxt(p)


def _split_into_rank_files(arr, out_dir):
    """Write `arr` as RANKS per-rank files with OVERLAP_ROWS shared rows
    between neighbours, and the rows shuffled within each file.

    This is what a real MPI run produces -- boundary nodes duplicated across
    ranks, node order arbitrary -- built deterministically from the reference
    so the expected canonical result is known exactly.
    """
    edges = np.linspace(0, arr.shape[0], RANKS + 1).astype(int)
    rng = np.random.default_rng(20260915)
    paths = []
    for r in range(RANKS):
        lo = max(0, edges[r] - (OVERLAP_ROWS if r else 0))
        block = arr[lo:edges[r + 1]]
        block = block[rng.permutation(block.shape[0])]
        p = os.path.join(str(out_dir), 'frt.txt%d' % r)
        F.write_canonical(block, p)
        paths.append(p)
    return paths


@pytest.mark.parametrize('name,rows', sorted(EXPECTED_ROWS.items()))
def test_canonical_reference_has_the_expected_node_count(name, rows):
    assert _canonical(name).shape == (rows, F.FRT_COLUMNS)


@pytest.mark.parametrize('name,rows', sorted(EXPECTED_ROWS.items()))
def test_multi_rank_output_reduces_to_the_reference(name, rows, tmp_path):
    # The property the canonical form exists for: split across ranks, with
    # duplicated boundary rows and arbitrary order, and it reduces back to the
    # exact same table -- bit-identical, not allclose.
    ref = _canonical(name)
    _split_into_rank_files(ref, tmp_path)
    got = F.canonical_from_case(str(tmp_path))
    assert got.shape == (rows, F.FRT_COLUMNS)
    assert np.array_equal(got, ref)


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_form_is_independent_of_file_order(name, tmp_path):
    files = _split_into_rank_files(_canonical(name), tmp_path)
    forward = F.canonicalize(F.load_frt(files))
    reverse = F.canonicalize(F.load_frt(list(reversed(files))))
    assert np.array_equal(forward, reverse)


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_rows_are_sorted_by_position(name):
    a = _canonical(name)
    order = np.lexsort((a[:, 2], a[:, 1], a[:, 0]))
    assert np.array_equal(order, np.arange(a.shape[0]))


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_coordinates_are_unique(name):
    a = _canonical(name)
    key = np.round(a[:, :3], F.COORD_DECIMALS)
    assert np.unique(key, axis=0).shape[0] == a.shape[0]


def test_the_multi_rank_fixtures_actually_contain_duplicates(tmp_path):
    """Guards against the dedupe becoming a no-op.

    If the fixtures had no duplicated boundary rows, every dedupe assertion
    above would pass on an identity function and this module would be
    vacuous -- the same tautology trap that let a source scan report PASS over
    thirteen live offenders earlier in this project.
    """
    ref = _canonical('test.tpv8')
    files = _split_into_rank_files(ref, tmp_path)
    raw = F.load_frt(files).shape[0]
    assert raw == ref.shape[0] + (RANKS - 1) * OVERLAP_ROWS
    assert raw > ref.shape[0]


def test_empty_file_list_raises_rather_than_returning_empty():
    # An empty array compares equal to another empty array; "no files" must not
    # be able to look like a pass.
    with pytest.raises(ValueError, match='no frt files'):
        F.load_frt([])


def test_wrong_column_count_raises(tmp_path):
    bad = tmp_path / 'frt.txt0'
    np.savetxt(str(bad), np.zeros((3, 5)))
    with pytest.raises(ValueError, match='columns'):
        F.load_frt([str(bad)])


def test_align_rejects_a_different_node_count():
    a = _canonical('test.tpv8')
    with pytest.raises(ValueError, match='node counts differ'):
        F.align(a, a[:-5])


def test_align_rejects_shifted_coordinates():
    a = _canonical('test.tpv8')
    b = a.copy()
    b[:, 0] += 500.0          # one grid spacing: a real misalignment
    with pytest.raises(ValueError, match='same nodes'):
        F.align(a, b)


def test_rank_files_are_ordered_numerically_not_lexically(tmp_path):
    for r in (0, 1, 2, 10):
        np.savetxt(str(tmp_path / f'frt.txt{r}'), np.zeros((1, F.FRT_COLUMNS)))
    got = [os.path.basename(p) for p in F.frt_rank_files(str(tmp_path))]
    assert got == ['frt.txt0', 'frt.txt1', 'frt.txt2', 'frt.txt10']


def test_serial_output_aligns_with_multi_rank_output(tmp_path):
    """The property the whole design exists for: a 1-file serial result and a
    3-file MPI result of the same physics align node for node."""
    ref = _canonical('test.tpv8')
    ranks = tmp_path / 'ranks'
    ranks.mkdir()
    raw = F.load_frt(_split_into_rank_files(ref, ranks))
    one = tmp_path / 'frt.txt0'
    F.write_canonical(ref, str(one))

    a, b = F.align(F.load_frt([str(one)]), raw)
    assert a.shape == b.shape
    assert np.abs(a[:, :3] - b[:, :3]).max() == 0.0


def test_canonicalize_refuses_to_choose_between_disagreeing_duplicates():
    """The dedupe drops a row; that is only lossless if the duplicates agree.

    Measured on all eight committed references: 358 duplicate groups, max
    intra-duplicate spread 0.000e+00 across all 22 columns. So exact equality
    is the right test, and silently keeping the first row would be a choice
    between two answers (rule 2).
    """
    arr = np.zeros((2, F.FRT_COLUMNS))
    arr[1, 5] = 1.0                      # same node, different physics
    with pytest.raises(ValueError, match='DIFFERENT values'):
        F.canonicalize(arr)
