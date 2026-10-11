"""Unit tests for scripts/readFaultTimeSeries.py (the src_ts* reader) and for
the field table it documents.

The synthetic files are written in the exact layout library_output.f90's
output_fault_ts writes (header lines copied from a real run's
src_ts_hdr0); every expected value is the known array the test wrote,
never the reader's own output.
"""
import os
import re

import numpy as np
import pytest

import readFaultTimeSeries as rfts

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NAMES = ['sliprate', 'sliprate_strike', 'sliprate_dip', 'slip_strike', 'slip_dip',
         'slip_normal', 'slip_path', 'sliprate_peak', 'traction_normal',
         'traction_strike', 'traction_dip', 'traction_shear']


def write_case(d, tag, rank, data, xyz, dt=0.04, nstep=None, friclaw=1, names=NAMES):
    """data: (nrec, nnode, nfield) array -> src_ts<tag><rank> + header."""
    nrec, nnode, nfield = data.shape
    with open(os.path.join(d, 'src_ts_hdr%s%d' % (tag, rank)), 'w') as f:
        f.write('# EQdyna whole-fault time series header\n')
        f.write('# binary = src_ts%s%d\n' % (tag, rank))
        f.write('# format = float64, native byte order, no record markers\n')
        f.write('# layout = [step][node][field], one record per time step\n')
        f.write('# fault = %d\n# rank = %d\n# friclaw = %d\n' % (2 if tag else 1, rank, friclaw))
        f.write('# nnode = %d\n# nfield = %d\n' % (nnode, nfield))
        f.write('# dt =  %.16E\n# nstep = %d\n' % (dt, nstep or nrec))
        f.write('# time of record k (0-based) = (k + 1) * dt\n')
        for k in range(nfield):
            f.write('# field %d %s m some meaning, with spaces\n' % (k, names[k]))
        f.write('# rows below: x y z (m) of each fault node, in binary node order\n')
        for row in xyz:
            f.write(' %24.16E %24.16E %24.16E\n' % tuple(row))
    data.astype('<f8').tofile(os.path.join(d, 'src_ts%s%d' % (tag, rank)))


def test_layout_is_step_node_field(tmp_path):
    nrec, nnode, nfield = 5, 3, len(NAMES)
    data = np.arange(nrec * nnode * nfield, dtype=float).reshape(nrec, nnode, nfield)
    xyz = np.array([[0.0, 0.0, -100.0], [100.0, 0.0, -100.0], [200.0, 0.0, -100.0]])
    write_case(str(tmp_path), '', 0, data, xyz, dt=0.25)
    ts = rfts.read_fault_ts(str(tmp_path / 'src_ts_hdr0'))
    for k, name in enumerate(NAMES):
        assert ts[name].shape == (nrec, nnode)
        np.testing.assert_array_equal(ts[name], data[:, :, k])
    np.testing.assert_array_equal(ts['xyz'], xyz)
    np.testing.assert_allclose(ts['t'], [0.25, 0.5, 0.75, 1.0, 1.25])
    assert ts['meta']['nfield'] == nfield and ts['fields'][3]['desc'] == 'some meaning, with spaces'
    # a field subset and the memory-mapped path give the same values
    sub = rfts.read_fault_ts(str(tmp_path / 'src_ts_hdr0'), fields=['slip_dip'], mmap=True)
    assert 'slip_strike' not in sub
    np.testing.assert_array_equal(sub['slip_dip'], data[:, :, 4])


def test_partial_run_reads_written_steps_only(tmp_path):
    data = np.random.default_rng(0).normal(size=(4, 2, len(NAMES)))
    write_case(str(tmp_path), '', 0, data, np.zeros((2, 3)), nstep=100)
    ts = rfts.read_fault_ts(str(tmp_path / 'src_ts_hdr0'))
    assert ts['t'].size == 4 and ts['meta']['nstep'] == 100


def test_truncated_binary_is_refused(tmp_path):
    data = np.ones((3, 2, len(NAMES)))
    write_case(str(tmp_path), '', 0, data, np.zeros((2, 3)))
    with open(tmp_path / 'src_ts0', 'ab') as f:
        f.write(b'\0' * 8)
    with pytest.raises(ValueError, match='whole number'):
        rfts.read_fault_ts(str(tmp_path / 'src_ts_hdr0'))


def test_read_case_concatenates_ranks_per_fault(tmp_path):
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=(6, 2, len(NAMES))), rng.normal(size=(6, 3, len(NAMES)))
    c = rng.normal(size=(6, 4, len(NAMES)))
    write_case(str(tmp_path), '', 0, a, np.full((2, 3), 0.0))
    write_case(str(tmp_path), '', 3, b, np.full((3, 3), 3.0))
    write_case(str(tmp_path), 'ft2_', 1, c, np.full((4, 3), 1.0))
    f1 = rfts.read_case(str(tmp_path), '')
    np.testing.assert_array_equal(f1['traction_dip'], np.concatenate([a, b], axis=1)[:, :, 10])
    np.testing.assert_array_equal(f1['ranks'], [0, 0, 3, 3, 3])
    f2 = rfts.read_case(str(tmp_path), 'ft2_')
    np.testing.assert_array_equal(f2['sliprate'], c[:, :, 0])


def _fortran_fields():
    """(name, unit) of every `call add(...)` in faultTsFields, file order."""
    src = open(os.path.join(REPO_ROOT, 'src', 'fortran', 'library_output.f90')).read()
    body = src[src.index('subroutine faultTsFields'):src.index('end subroutine faultTsFields')]
    return re.findall(r"call add\(\s*\w+\s*,\s*'(\w+)'\s*,\s*'([^']*)'", body)


def test_fortran_common_fields_match_reader_docs():
    names = [n for n, _ in _fortran_fields()]
    assert names[:12] == NAMES
    # every field the solver can write is documented for users
    doc = open(os.path.join(REPO_ROOT, 'docs', 'user', 'outputs.md')).read()
    missing = [n for n in dict.fromkeys(names)
               if '`%s`' % n not in doc and not re.match(r'vel_(master|slave)_[xyz]$', n)]
    assert not missing, missing
    assert 'vel_master_x/y/z' in doc and 'vel_slave_x/y/z' in doc
