"""run_mpi_scaling.main keeps going past a refused point (2026-09-27).

A refusal (a non-positive per-step figure) used to propagate out of main()
and drop every point of the pass: one refused 8-rank point lost 2 of 3
passes. Now the refused repeat is recorded on its row, never as a number,
the other points are measured and saved, and the run exits non-zero.

Driven end to end through main() with the box and the solvers stubbed out:
no MPI, jax, mesh or Fortran binary runs, and the snapshot and ledger go to
tmp_path, not docs/."""
import importlib.util
import json
import os
import sys

import pytest

from conftest import REPO_ROOT

PERF = os.path.join(REPO_ROOT, 'testsys', 'perf')


def _load(monkeypatch):
    monkeypatch.setenv('EQDYNAROOT', REPO_ROOT)
    sys.path.insert(0, PERF)
    spec = importlib.util.spec_from_file_location(
        'run_mpi_scaling', os.path.join(PERF, 'run_mpi_scaling.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _jax_point(ms):
    return dict(ms_per_step=ms, ms_per_step_wall=None, compile_s=1.0,
                solve_lo_s=1.0, solve_hi_s=2.0, wall_lo_s=3.0, wall_hi_s=4.0,
                rank_ms=[ms], rank_ms_max=ms, rank_ms_mean=ms, straggler=1.0,
                mpi_ms=[0.1], exchange_ms_min=0.1, exchange_excess_ms_est=[0.0],
                wait_ms=[0.0], eff=[1.0], threads=[1], cpus_allowed=[1],
                Ei=[10], Ep=[5], halo=[3], platform='cpu', gpu_mem=None,
                device_peak_gb=[None])


def _stub_box(ms, monkeypatch, tmp_path, appended):
    rs = ms.rs
    import ledger
    monkeypatch.setattr(ms, 'OUT', str(tmp_path / 'snap.json'))
    monkeypatch.setattr(ms.numa, 'numa_topology', lambda: {0: [2, 3, 4, 5]})
    monkeypatch.setattr(ms, 'spread_cpus',
                        lambda nodes, k, exclude=(), max_busy=1.0:
                        (list(range(2, 2 + k)), {c: 0.0 for c in range(2, 2 + k)}))
    monkeypatch.setattr(rs, 'build_py_case', lambda case: str(tmp_path))
    monkeypatch.setattr(rs, 'make_case', lambda *a, **k: None)
    monkeypatch.setattr(rs, 'read_term_dt', lambda d: (1.0, 0.01))
    monkeypatch.setattr(ledger, 'box_tenancy',
                        lambda b: dict(busy=0, total=4))
    real_rows = ledger.rows_from_mpi_scaling_snapshot
    monkeypatch.setattr(ledger, 'rows_from_mpi_scaling_snapshot',
                        lambda meta, snap, ten: real_rows(
                            dict(meta, sha='abc1234'),
                            'docs/perf_snapshots/fake.json', ten))
    monkeypatch.setattr(ledger, 'append_rows',
                        lambda rows: appended.extend(rows) or len(rows))


def test_refused_point_does_not_drop_the_pass(monkeypatch, tmp_path, capsys):
    ms = _load(monkeypatch)
    rs = ms.rs
    out = str(tmp_path / 'snap.json')
    appended = []
    _stub_box(ms, monkeypatch, tmp_path, appended)

    def fake_jax(case_dir, cpus, n, n_lo, n_hi, sync, platform, warmup):
        if n == 2:
            raise ms.numa.NonPositiveMeasurement(
                'per-step by difference over rank solve time came out -1.000 '
                'ms at 2 ranks')
        return _jax_point(100.0 / n)
    monkeypatch.setattr(ms, 'per_step_jax_mpi', fake_jax)
    monkeypatch.setattr(rs, 'per_step_fortran',
                        lambda work, n, *a: (0.1 / n, 1.0, 2.0, 3.0, '', ''))
    monkeypatch.setattr(sys, 'argv', ['run_mpi_scaling', '--ranks', '1,2,4',
                                      '--placement', 'spread'])

    with pytest.raises(SystemExit) as ei:
        ms.main()
    assert 'refused' in str(ei.value) and '(2, ' in str(ei.value)

    rows = json.load(open(out))['rows']
    by = {r['ranks']: r for r in rows}
    assert sorted(by) == [1, 2, 4]
    # the points either side of the refusal were measured and kept
    assert by[1]['jax']['ms_per_step'] == 100.0
    assert by[4]['jax']['ms_per_step'] == 25.0
    # the refused point is recorded as a refusal, never as a jax number,
    # and its Fortran column was still measured
    assert 'jax' not in by[2] and 'jax_halo' not in by[2]
    assert by[2]['refused'][0]['backend'] == 'jax_halo'
    assert 'solve time' in by[2]['refused'][0]['reason']
    assert by[2]['fortran']['ms_per_step'] == pytest.approx(50.0)
    assert 'REFUSED' in capsys.readouterr().out
    # the ledger got the measured points only: 2 jax + 3 fortran
    assert sorted((r['backend'], r['ranks']) for r in appended) == [
        ('fortran', 1), ('fortran', 2), ('fortran', 4),
        ('python-jax-mpi', 1), ('python-jax-mpi', 4)]


def test_partial_refusal_failed_mpirun_and_fortran_refusal(monkeypatch, tmp_path):
    """--repeats 2: jax at 1 rank loses one repeat to a failed mpirun (None)
    and keeps the other; Fortran at 2 ranks is refused on both repeats. The
    counted-repeat numbers are recorded, both failures are refusals, and the
    run exits non-zero."""
    ms = _load(monkeypatch)
    appended = []
    _stub_box(ms, monkeypatch, tmp_path, appended)
    calls = {'jax': 0}

    def fake_jax(case_dir, cpus, n, n_lo, n_hi, sync, platform, warmup):
        calls['jax'] += 1
        return None if (n == 1 and calls['jax'] == 1) else _jax_point(100.0 / n)

    def fake_fortran(work, n, *a):
        if n == 2:
            raise ms.numa.NonPositiveMeasurement('per-step by difference came '
                                                 'out -0.1 s/step')
        return (0.1 / n, 1.0, 2.0, 3.0, '', '')
    monkeypatch.setattr(ms, 'per_step_jax_mpi', fake_jax)
    monkeypatch.setattr(ms.rs, 'per_step_fortran', fake_fortran)
    monkeypatch.setattr(sys, 'argv', ['run_mpi_scaling', '--ranks', '1,2',
                                      '--placement', 'spread', '--repeats', '2'])
    with pytest.raises(SystemExit) as ei:
        ms.main()
    assert 'refused' in str(ei.value)
    by = {r['ranks']: r for r in json.load(open(tmp_path / 'snap.json'))['rows']}
    assert by[1]['jax']['repeats_counted'] == 1
    assert by[1]['refused'] == [dict(backend='jax_halo',
                                     reason='mpirun returned non-zero')]
    assert by[1]['fortran']['repeats_counted'] == 2
    assert 'fortran' not in by[2]
    assert [x['backend'] for x in by[2]['refused']] == ['fortran', 'fortran']
    assert by[2]['jax']['repeats_counted'] == 2
    assert ('fortran', 2) not in [(r['backend'], r['ranks']) for r in appended]


def test_a_non_refusal_error_still_stops_the_run(monkeypatch, tmp_path):
    """Only numa.NonPositiveMeasurement is a per-point refusal. A guard such
    as the wrong-device check raises a plain RuntimeError and must stop the
    sweep at once, not be logged as REFUSED while every later point runs."""
    ms = _load(monkeypatch)
    _stub_box(ms, monkeypatch, tmp_path, [])
    seen = []

    def fake_jax(case_dir, cpus, n, n_lo, n_hi, sync, platform, warmup):
        seen.append(n)
        raise RuntimeError("asked for platform 'cpu'; the ranks reported gpu")
    monkeypatch.setattr(ms, 'per_step_jax_mpi', fake_jax)
    monkeypatch.setattr(sys, 'argv', ['run_mpi_scaling', '--ranks', '1,2',
                                      '--placement', 'spread'])
    with pytest.raises(RuntimeError, match='asked for platform'):
        ms.main()
    assert seen == [1]
