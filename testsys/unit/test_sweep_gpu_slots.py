#! /usr/bin/env python3
"""
The sweep's GPU device pool (pathway item 58).

`test.tpv36` x python-jax failed the 2026-09-22 GPU sweep with
RESOURCE_EXHAUSTED and passed at 7.06e-09 when run alone: XLA preallocates a
fraction of the card per process, the sweep runs cells concurrently, and two
cells cannot both hold 75% of one device. That is a harness defect, not a
parity one, and the harness fix is one concurrent python-jax cell per visible
card with the preallocation stated rather than inherited.

These tests cover the scheduling and the refusals, not the CUDA runtime: no
GPU is required to run them, which is the point -- CI has no card and must
still catch a change that lets two cells share one.
"""
import importlib.util
import os
import sys
import threading

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTSYS)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

_spec = importlib.util.spec_from_file_location(
    'run_e2e_under_test', os.path.join(TESTSYS, 'e2e', 'run_e2e.py'))
E = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(E)


def test_slots_never_exceed_the_device_count_and_stay_distinct():
    """Six cells, two cards: at most two in flight, and no two holders ever
    hold the same index. A semaphore without the index would allow the
    second, which is the collision itself."""
    slots = E.GpuSlots([0, 1])
    lock = threading.Lock()
    held = set()
    peak = [0]
    collisions = []

    def worker():
        idx = slots.acquire()
        try:
            with lock:
                if idx in held:
                    collisions.append(idx)
                held.add(idx)
                peak[0] = max(peak[0], len(held))
            for _ in range(200):           # hold it long enough to overlap
                threading.Event().wait(0.0005)
        finally:
            with lock:
                held.discard(idx)
            slots.release(idx)

    threads = [threading.Thread(target=worker) for _ in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=60)
    assert not any(t.is_alive() for t in threads), 'a cell never got a device'
    assert collisions == []
    assert peak[0] == 2, 'peak concurrency %d, want 2' % peak[0]
    assert sorted(slots._free) == [0, 1], 'devices were not all returned'


def test_gpu_env_pins_the_card_and_states_the_fraction():
    env = E.gpu_env({'PATH': '/usr/bin'}, 3)
    assert env['CUDA_VISIBLE_DEVICES'] == '3'
    assert env['XLA_PYTHON_CLIENT_MEM_FRACTION'] == E.GPU_MEM_FRACTION
    assert env['PATH'] == '/usr/bin', 'gpu_env must not drop the caller env'


def test_run_standalone_refuses_a_gpu_cell_with_no_slot():
    """No fallback. Without a reserved card this used to launch anyway and
    race whatever else the sweep had started."""
    with pytest.raises(RuntimeError, match='needs a GPU slot'):
        E.run_standalone('/nonexistent/test.tpv36.python-jax', 'python-jax',
                         device='cuda', env={}, gpu_index=None)


def test_run_cell_refuses_a_gpu_cell_with_no_pool(monkeypatch):
    monkeypatch.setattr(E, 'make_serial_case', lambda *a, **k: None)
    with pytest.raises(RuntimeError, match='no GpuSlots pool'):
        E.run_cell('test.tpv36', 'python-jax', '/nonexistent', None, {},
                   'cuda', gpu_slots=None)


def test_cpu_sweeps_are_untouched(monkeypatch):
    """A cpu sweep must not acquire, pin, or even ASK about a card -- a
    CPU-only box would otherwise fail on nvidia-smi."""
    def boom():
        raise AssertionError('visible_gpu_indices called on a cpu sweep')
    monkeypatch.setattr(E, 'visible_gpu_indices', boom)
    calls = {}

    def fake_standalone(case_dir, backend, device='cpu', env=None, gpu_index=None):
        calls['gpu_index'] = gpu_index
        calls['device'] = device
        return 'frt'

    monkeypatch.setattr(E, 'make_serial_case', lambda *a, **k: None)
    monkeypatch.setattr(E, 'run_standalone', fake_standalone)
    E.run_cell('test.tpv8', 'python-jax', '/tmp', None, {}, 'cpu')
    assert calls == {'gpu_index': None, 'device': 'cpu'}


def test_visible_gpu_indices_honours_the_override(monkeypatch):
    monkeypatch.setenv('EQDYNA_E2E_GPUS', '0,2,3')
    assert E.visible_gpu_indices() == [0, 2, 3]


def test_visible_gpu_indices_refuses_an_empty_answer(monkeypatch):
    """`nvidia-smi` present but reporting nothing must not degrade to one
    slot or to the CPU: a --device cuda sweep that ran on the CPU would
    report a GPU measurement it never made."""
    monkeypatch.delenv('EQDYNA_E2E_GPUS', raising=False)

    class R:
        returncode = 0
        stdout = ''
        stderr = ''

    monkeypatch.setattr(E.subprocess, 'run', lambda *a, **k: R())
    with pytest.raises(RuntimeError, match='no CUDA devices'):
        E.visible_gpu_indices()
