"""Regression cover for board row 56 (owner ruling 2026-09-24): no `--device
auto`. The default is cpu for serial AND --mpi runs; GPU only on an explicit
`--device cuda`; `--device auto` is refused naming the two valid choices."""
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


def test_explicit_cuda_is_passed_through_serial_and_mpi(monkeypatch):
    assert _device_seen(monkeypatch, ['case', '--device', 'cuda']) == ['cuda']
    assert _device_seen(monkeypatch, ['case', '--mpi', '--device', 'cuda']) == ['cuda']


def test_device_auto_is_refused_naming_both_choices():
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(E.__file__)))
    env = dict(os.environ, PYTHONPATH=pkg_root)
    r = subprocess.run([sys.executable, '-m', 'eqdyna', 'no_such_case', '--device', 'auto'],
                       capture_output=True, text=True, timeout=60, env=env)
    assert r.returncode == 2
    assert "'auto' is not a device" in r.stderr
    assert 'cpu' in r.stderr and 'cuda' in r.stderr
    assert 'row 56' in r.stderr


def test_numpy_refuses_explicit_gpu(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['eqdyna', 'case', '--backend', 'numpy', '--device', 'cuda'])
    with pytest.raises(SystemExit, match='meaningless with --backend numpy'):
        E.main()


def test_old_gpu_spelling_is_refused_with_hint():
    import argparse
    with pytest.raises(argparse.ArgumentTypeError, match="spelled 'cuda'"):
        E._device_arg('gpu')


def test_sweep_passes_device_explicitly(monkeypatch, tmp_path):
    """The critical row-56 audit finding: run_e2e's GPU sweep set only
    JAX_PLATFORMS=cuda, which the cpu default then overwrote -> a green sweep
    on the CPU. The launched command must carry --device itself."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'e2e'))
    import run_e2e
    seen = []

    def fake_call(cmd, cwd, env, prefix):
        seen.append(cmd)
        return 0, 'o', 'e'
    monkeypatch.setattr(run_e2e, '_call_kept', fake_call)
    monkeypatch.setattr(run_e2e, 'gpu_env', lambda env, idx: env)
    for dev, idx in (('cpu', None), ('cuda', 0)):
        seen.clear()
        try:
            run_e2e.run_standalone(str(tmp_path), 'python-jax', device=dev, gpu_index=idx)
        except Exception:
            pass   # post-run artifact checks fail on the fake case dir; the launch is what is tested
        assert seen, 'run_standalone never launched'
        cmd = seen[0]
        assert cmd[cmd.index('--device') + 1] == dev
