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

# Measured from the committed references. These are the post-dedupe counts, so
# they fail if the dedupe silently becomes a no-op.
EXPECTED_ROWS = {
    'test.tpv8': 1891,
    'test.tpv10': 1891,
    'test.tpv104': 2701,
    'test.tpv1053d': 4005,
    'test.drv.a6': 5151,
    'test.tpv29': 3321,
    'test.meng2023a': 651,
    'test.meng2023cb': 651,
}


def _case(name):
    d = os.path.join(REF, name)
    if not F.frt_rank_files(d):
        pytest.skip(f'{name} reference not present')
    return d


@pytest.mark.parametrize('name,rows', sorted(EXPECTED_ROWS.items()))
def test_canonical_row_count_matches_the_deduped_node_set(name, rows):
    assert F.canonical_from_case(_case(name)).shape == (rows, F.FRT_COLUMNS)


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_form_is_independent_of_file_order(name):
    d = _case(name)
    files = F.frt_rank_files(d)
    forward = F.canonicalize(F.load_frt(files))
    reverse = F.canonicalize(F.load_frt(list(reversed(files))))
    # bit-identical, not allclose: order independence is exact or it is not a
    # canonical form.
    assert np.array_equal(forward, reverse)


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_rows_are_sorted_by_position(name):
    a = F.canonical_from_case(_case(name))
    order = np.lexsort((a[:, 2], a[:, 1], a[:, 0]))
    assert np.array_equal(order, np.arange(a.shape[0]))


@pytest.mark.parametrize('name', sorted(EXPECTED_ROWS))
def test_canonical_coordinates_are_unique(name):
    a = F.canonical_from_case(_case(name))
    key = np.round(a[:, :3], F.COORD_DECIMALS)
    assert np.unique(key, axis=0).shape[0] == a.shape[0]


def test_multi_rank_references_actually_contain_duplicates():
    """Guards against the dedupe becoming a no-op.

    If no reference had duplicate boundary nodes, every dedupe assertion above
    would pass on an identity function and this whole module would be
    vacuous -- the same tautology trap that let a source scan report PASS over
    thirteen live offenders earlier in this project.
    """
    found = []
    for name in EXPECTED_ROWS:
        d = os.path.join(REF, name)
        if not F.frt_rank_files(d):
            continue
        raw = F.load_frt(F.frt_rank_files(d)).shape[0]
        if raw > EXPECTED_ROWS[name]:
            found.append((name, raw, EXPECTED_ROWS[name]))
    assert found, ('no reference has duplicate boundary rows, so the dedupe is '
                   'untested by the assertions in this file')


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
    a = F.canonical_from_case(_case('test.tpv8'))
    with pytest.raises(ValueError, match='node counts differ'):
        F.align(a, a[:-5])


def test_align_rejects_shifted_coordinates():
    a = F.canonical_from_case(_case('test.tpv8'))
    b = a.copy()
    b[:, 0] += 500.0          # one grid spacing: a real misalignment
    with pytest.raises(ValueError, match='same nodes'):
        F.align(a, b)


def test_rank_files_are_ordered_numerically_not_lexically(tmp_path):
    for r in (0, 1, 2, 10):
        np.savetxt(str(tmp_path / f'frt.txt{r}'), np.zeros((1, F.FRT_COLUMNS)))
    got = [os.path.basename(p) for p in F.frt_rank_files(str(tmp_path))]
    assert got == ['frt.txt0', 'frt.txt1', 'frt.txt2', 'frt.txt10']


def test_serial_output_aligns_with_a_multi_rank_reference(tmp_path):
    """The property the whole design exists for.

    Builds a synthetic 'serial' output by canonicalising a multi-rank
    reference, then aligns it back against the raw per-rank files. A 1-file
    result and a 2-file result must reduce to the same node set.
    """
    d = _case('test.tpv8')
    raw = F.load_frt(F.frt_rank_files(d))
    serial = F.canonicalize(raw)
    one = tmp_path / 'frt.txt0'
    F.write_canonical(serial, str(one))

    a, b = F.align(F.load_frt([str(one)]), raw)
    assert a.shape == b.shape
    assert np.abs(a[:, :3] - b[:, :3]).max() == 0.0
