"""Regression cover for `eqdyna3d.serial_buckets` (master CI run 36082084976).

run_case once built `io` from the 'write frt' phase alone, so the 'write
stations' phase added with station output (row 114) fell into unaccounted_s
and failed profile_schema's 5% sum check on a CI runner (0.80 s of 15.40 s).
These tests assert the behaviour: every recorded phase lands in a bucket, and
an unmapped phase is refused instead of silently vanishing.
"""
import pytest

from eqdyna import eqdyna3d as E
from eqdyna import profile_emit


def _prof(**phases):
    p = E.Profile('numpy')
    for name, secs in phases.items():
        p[name.replace('_', ' ')] = secs
    return p


def test_every_phase_is_bucketed_including_station_writes():
    p = E.Profile('numpy')
    p.update({'setup (mesh+input)': 2.0, 'resolve solver': 0.5, 'solve': 10.0,
              'write frt': 0.25, 'write stations': 0.8})
    b = E.serial_buckets(p, fault_s=3.0)
    assert set(b) == set(profile_emit.BUCKET_KEYS)
    assert sum(b.values()) == pytest.approx(sum(p.values()))
    assert b['io'] == pytest.approx(1.05)
    assert b['setup'] == pytest.approx(2.5)
    assert b['element'] == pytest.approx(7.0) and b['fault'] == 3.0
    assert b['exchange'] == 0.0 and b['wait'] == 0.0


def test_unmapped_phase_is_refused():
    p = E.Profile('numpy')
    p.update({'solve': 1.0, 'write something new': 0.3})
    with pytest.raises(ValueError, match='write something new'):
        E.serial_buckets(p, fault_s=0.0)
