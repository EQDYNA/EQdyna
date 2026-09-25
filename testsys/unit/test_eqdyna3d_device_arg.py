"""Regression cover for board row 56 (owner ruling 2026-09-24): no `--device
auto`. The default is cpu for serial AND --mpi runs; GPU only on an explicit
`--device gpu`; `--device auto` is refused naming the two valid choices."""
import os
import subprocess
import sys

import pytest

from eqdyna import eqdyna3d as E


class _Stop(Exception):
    pass


def _device_seen(monkeypatch, argv):
    seen = []

    def record(device):
        seen.append(device)
        raise _Stop()
    monkeypatch.setattr(E, '_select_device', record)
    monkeypatch.setattr(sys, 'argv', ['eqdyna'] + argv)
    with pytest.raises(_Stop):
        E.main()
    return seen


def test_default_device_is_cpu_serial_and_mpi(monkeypatch):
    assert _device_seen(monkeypatch, ['case']) == ['cpu']
    assert _device_seen(monkeypatch, ['case', '--mpi']) == ['cpu']


def test_explicit_gpu_is_passed_through_serial_and_mpi(monkeypatch):
    assert _device_seen(monkeypatch, ['case', '--device', 'gpu']) == ['gpu']
    assert _device_seen(monkeypatch, ['case', '--mpi', '--device', 'gpu']) == ['gpu']


def test_device_auto_is_refused_naming_both_choices():
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(E.__file__)))
    env = dict(os.environ, PYTHONPATH=pkg_root)
    r = subprocess.run([sys.executable, '-m', 'eqdyna', 'no_such_case', '--device', 'auto'],
                       capture_output=True, text=True, timeout=60, env=env)
    assert r.returncode == 2
    assert "'auto' is not a device" in r.stderr
    assert 'cpu' in r.stderr and 'gpu' in r.stderr
    assert 'row 56' in r.stderr


def test_numpy_refuses_explicit_gpu(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['eqdyna', 'case', '--backend', 'numpy', '--device', 'gpu'])
    with pytest.raises(SystemExit, match='meaningless with --backend numpy'):
        E.main()
