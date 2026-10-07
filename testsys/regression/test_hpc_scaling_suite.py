#! /usr/bin/env python3
"""
Regression guard: testsys/perf/hpc_scaling_suite.py (board item 142)'s
collect/analyze arithmetic and refusal logic.

WHAT THIS DOES NOT DO: launch a real HPC submission, or run the actual
`--submit --test size` until-OOM sweep -- this box is memory-constrained
(measured full swap in a prior session) and that real run is deliberately
deferred, per the task spec, to real HPC hardware. This file is a pure
unit/regression guard over the --collect/--analyze arithmetic, using real
committed per-rank profile fixtures
(testsys/regression/fixtures/profile_guard/fortran/, a real 4-rank
test.tpv8 sweep) reshaped into synthetic lo/hi pairs -- not a fresh
solver launch.

WHY THIS MATTERS: the per-step-by-difference arithmetic
(run_numa_scaling.per_step_and_fixed) breaks if you difference a FIXED
cost (setup, io -- one-time, does not scale with step count) the same way
as a PER-STEP cost (element, fault, exchange, wait). A first draft of this
script's per_rank_metrics() did exactly that and tripped
NonPositiveMeasurement on synthetic data where the fixed buckets came out
to ~0 when differenced over nsteps -- a bug in the suite's own code, not a
genuine measurement. This file pins that split so it cannot regress:
PER_STEP_BUCKETS are differenced, FIXED_BUCKETS are read undifferenced off
the hi run.

Also pins: `level=big` refuses `backend=python-jax-mpi` IN CODE (the
owner's spec), and `--analyze` never crashes the whole report when one
point's parity check raises (e.g. a run directory with no frt output at
all) -- it declares that point 'fail' with the exception text and keeps
going.

Cheap (rule 9): no solver launch, no subprocess, no real sbatch -- pure
in-memory/tempdir arithmetic against committed fixtures and synthetic
data. Exits non-zero on any failure.
"""
import copy
import json
import os
import sys
import tarfile
import tempfile

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
for _p in (TESTSYS, REPO_ROOT, os.path.join(REPO_ROOT, 'src', 'python')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from testsys.perf import hpc_scaling_suite as h  # noqa: E402
import run_numa_scaling as numa                   # noqa: E402

FIXTURE = os.path.join(TESTSYS, 'regression', 'fixtures', 'profile_guard',
                       'fortran', 'profile.rank0.json')


def _load_fixture_row():
    with open(FIXTURE) as f:
        return json.load(f)


def _scaled_row(base, nsteps, per_step_scale):
    """A synthetic row derived from the real fixture: PER_STEP_BUCKETS scaled
    by (nsteps/base_nsteps)*per_step_scale (so lo/hi differ only in step
    count, as a real run would), FIXED_BUCKETS held at the fixture's own
    value (a one-time cost does not grow with nsteps)."""
    row = copy.deepcopy(base)
    base_nsteps = base['nsteps']
    buckets = dict(row['buckets_s'])
    per_step_total = 0.0
    for key in h.PER_STEP_BUCKETS:
        per_rank_step = base['buckets_s'][key] / base_nsteps
        buckets[key] = per_rank_step * nsteps * per_step_scale
        per_step_total += buckets[key]
    fixed_total = sum(base['buckets_s'][key] for key in h.FIXED_BUCKETS)
    row['nsteps'] = nsteps
    row['buckets_s'] = buckets
    row['total_s'] = per_step_total + fixed_total
    row['loop_s'] = per_step_total
    row['unaccounted_s'] = 0.0
    return row


def check_bucket_partition():
    assert set(h.PER_STEP_BUCKETS) | set(h.FIXED_BUCKETS) == set(h.profile_schema.BUCKET_KEYS)
    assert set(h.PER_STEP_BUCKETS) & set(h.FIXED_BUCKETS) == set()
    print('  PASS  PER_STEP_BUCKETS/FIXED_BUCKETS partition profile_schema.BUCKET_KEYS exactly')


def check_big_refuses_python_jax_mpi():
    for test in h.TEST_TYPES:
        try:
            h._validate_combo(test, 'big', 'python-jax-mpi')
        except h.SuiteError:
            pass
        else:
            raise AssertionError(
                'level=big, backend=python-jax-mpi, test=%r was NOT refused' % test)
    # fortran at big must NOT be refused by this particular rule.
    for test in h.TEST_TYPES:
        h._validate_combo(test, 'big', 'fortran')
    print('  PASS  level=big refuses backend=python-jax-mpi in code, for every test type')


def check_per_rank_metrics_fixed_vs_per_step():
    """The core arithmetic fix: FIXED_BUCKETS must come back as the hi run's
    own raw value (undifferenced); PER_STEP_BUCKETS must come back as a
    positive per-step figure equal to the per-rank-step rate used to build
    the synthetic rows."""
    base = _load_fixture_row()
    lo = _scaled_row(base, 20, 1.0)
    hi = _scaled_row(base, 60, 1.0)
    per_rank = h.per_rank_metrics([lo], [hi])
    r = per_rank[0]
    for key in h.FIXED_BUCKETS:
        expected = base['buckets_s'][key]
        got = r['buckets_fixed_s'][key]
        assert abs(got - expected) < 1e-9, (
            'FIXED_BUCKETS[%r] was differenced instead of read raw: got %r, expected %r'
            % (key, got, expected))
    for key in h.PER_STEP_BUCKETS:
        expected_rate = base['buckets_s'][key] / base['nsteps']
        got = r['buckets_per_step_s'][key]
        assert got > 0, 'PER_STEP_BUCKETS[%r] per-step rate came out non-positive: %r' % (key, got)
        assert abs(got - expected_rate) < 1e-6, (
            'PER_STEP_BUCKETS[%r] per-step rate: got %r, expected %r'
            % (key, got, expected_rate))
    assert r['per_step_s'] > 0
    print('  PASS  per_rank_metrics: FIXED_BUCKETS read raw off hi, '
          'PER_STEP_BUCKETS differenced and positive')


def check_fixed_buckets_would_trip_guard_if_differenced():
    """Reproduce the exact bug this file exists to pin against: differencing
    a FIXED bucket (constant across lo/hi, as a one-time cost is) over
    nsteps must hit run_numa_scaling's own non-positive guard -- proving the
    guard the repo already has is a real backstop, and that
    per_rank_metrics must route FIXED_BUCKETS around it, not through it."""
    base = _load_fixture_row()
    lo = _scaled_row(base, 20, 1.0)
    hi = _scaled_row(base, 60, 1.0)
    key = h.FIXED_BUCKETS[0]
    try:
        numa.per_step_and_fixed(lo['buckets_s'][key], hi['buckets_s'][key],
                                lo['nsteps'], hi['nsteps'])
    except numa.NonPositiveMeasurement:
        pass
    else:
        raise AssertionError(
            'differencing a FIXED bucket over nsteps did not trip the '
            'non-positive guard -- the fixture no longer demonstrates the bug '
            'this test pins against')
    print('  PASS  differencing a FIXED bucket over nsteps still trips '
          'NonPositiveMeasurement (confirms why FIXED_BUCKETS must bypass it)')


def check_build_points_strong_weak_size():
    machine_cfg = dict(cores_per_node=128)
    for test in ('strong', 'weak', 'size'):
        points = h.build_points(test, 'medium', machine_cfg, max_nodes=4)
        assert len(points) >= 2, '%s: expected >=2 scaling points, got %d' % (test, len(points))
        ids = [p['point_id'] for p in points]
        assert len(ids) == len(set(ids)), '%s: duplicate point_id(s): %r' % (test, ids)
        for p in points:
            assert p['ranks'] >= 1
            assert len(p['decomp']) == 3
            assert p['decomp'][0] * p['decomp'][1] * p['decomp'][2] == p['ranks'], (
                '%s point %s: decomp %r does not multiply to ranks %d'
                % (test, p['point_id'], p['decomp'], p['ranks']))
    strong = h.build_points('strong', 'medium', machine_cfg, max_nodes=4)
    total_elems = [round(h.elements_for_dx(p['dx'])) for p in strong]
    spread = (max(total_elems) - min(total_elems)) / max(total_elems)
    assert spread < 0.05, (
        'strong scaling must hold total element count ~fixed while ranks grow: '
        'got %r (spread %.1f%%)' % (total_elems, spread * 100))
    weak = h.build_points('weak', 'medium', machine_cfg, max_nodes=4)
    per_core = [h.elements_for_dx(p['dx']) / p['ranks'] for p in weak]
    spread = (max(per_core) - min(per_core)) / max(per_core)
    assert spread < 0.05, (
        'weak scaling must hold elements-per-core ~fixed while ranks grow: '
        'got %r (spread %.1f%%)' % (per_core, spread * 100))
    size = h.build_points('size', 'medium', machine_cfg, max_nodes=4)
    ranks = {p['ranks'] for p in size}
    assert len(ranks) == 1, 'size test must hold core count fixed: got rank counts %r' % ranks
    dxs = [p['dx'] for p in size]
    assert len(set(dxs)) == len(dxs), 'size test must sweep dx, got duplicates: %r' % dxs
    print('  PASS  build_points: strong holds elements fixed, weak holds '
          'elements/core fixed, size holds ranks fixed and sweeps dx '
          '(checked %d/%d/%d points)' % (len(strong), len(weak), len(size)))


def _write_manifest_and_profiles(suite_dir, point, lo_nsteps, hi_nsteps,
                                 nranks, base_row, case='test.tpv8', backend='fortran'):
    os.makedirs(suite_dir, exist_ok=True)
    point = dict(point)
    point['runs'] = {}
    for tag, nsteps in (('lo', lo_nsteps), ('hi', hi_nsteps)):
        rd = os.path.join(suite_dir, point['point_id'], tag)
        os.makedirs(rd)
        for r in range(nranks):
            row = _scaled_row(base_row, nsteps, 1.0)
            row['rank'] = r
            row['nranks'] = nranks
            with open(os.path.join(rd, 'profile.rank%d.json' % r), 'w') as f:
                json.dump(row, f)
        point['runs'][tag] = dict(run_dir=rd, nsteps=nsteps)
    manifest = dict(schema=h.SUITE_SCHEMA_ID, machine='ubuntu', backend=backend,
                    test='strong', level='medium', case=case, account='x',
                    max_nodes=1, nsteps_lo=lo_nsteps, nsteps_hi=hi_nsteps,
                    created_utc='test', points=[point])
    with open(os.path.join(suite_dir, h.MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f)
    return manifest


def check_collect_and_analyze_round_trip():
    """cmd_collect -> analyze on two synthetic points derived from the real
    fixture: point 0 (case/backend the matrix table gates, no real frt
    output present) must come back 'fail' with a readable reason rather than
    raising; point 1 (not the suite's first point) must come back
    'declared', per the no-silent-skip contract, and must never touch
    compare_cell at all."""
    base_row = _load_fixture_row()
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_suite_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    point0 = dict(point_id='strong_n0001', ranks=2, nnodes=1, dx=170.3,
                 decomp=[2, 1, 1], backend='fortran')
    point1 = dict(point_id='strong_n0002', ranks=4, nnodes=1, dx=120.0,
                 decomp=[2, 2, 1], backend='fortran')
    manifest = _write_manifest_and_profiles(
        suite_dir, point0, 20, 60, 2, base_row, case='test.tpv8', backend='fortran')
    # add the second point onto the same manifest/suite_dir
    rd_lo = os.path.join(suite_dir, point1['point_id'], 'lo')
    rd_hi = os.path.join(suite_dir, point1['point_id'], 'hi')
    for rd, nsteps in ((rd_lo, 20), (rd_hi, 60)):
        os.makedirs(rd)
        for r in range(4):
            row = _scaled_row(base_row, nsteps, 1.0)
            row['rank'] = r
            row['nranks'] = 4
            with open(os.path.join(rd, 'profile.rank%d.json' % r), 'w') as f:
                json.dump(row, f)
    point1 = dict(point1, runs=dict(lo=dict(run_dir=rd_lo, nsteps=20),
                                    hi=dict(run_dir=rd_hi, nsteps=60)))
    manifest['points'].append(point1)
    with open(os.path.join(suite_dir, h.MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f)

    import argparse
    args = argparse.Namespace(suite_dir=suite_dir, out=None)
    h.cmd_collect(args)
    tgz = os.path.join(suite_dir, 'hpc_scaling_ubuntu_fortran_strong_medium.tgz')
    assert os.path.isfile(tgz), 'cmd_collect did not produce the expected tarball'
    with tarfile.open(tgz, 'r:gz') as tar:
        names = tar.getnames()
    assert any('profile.rank' in n for n in names), (
        'collected tarball has no profile.rank*.json members: %r' % names)

    report = h.analyze(tgz)
    assert report['schema'] == 'eqdyna-hpc-scaling-report/1'
    assert len(report['points']) == 2
    p0, p1 = report['points']
    assert p0['point_id'] == 'strong_n0001'
    assert p0['parity_status'] == 'fail', (
        'point 0 (gated cell, no real frt output) should be a declared FAIL, '
        'not a silent pass or a crash: got %r' % p0['parity_status'])
    assert any('compare_cell raised' in ln or 'frt' in ln.lower()
              for ln in p0['parity_lines']), (
        'point 0 fail reason is not legible: %r' % p0['parity_lines'])
    assert p1['parity_status'] == 'declared', (
        'point 1 (not the suite\'s first point) must be DECLARED incomparable, '
        'got %r' % p1['parity_status'])
    for key in h.PER_STEP_BUCKETS:
        assert key in p0['bucket_breakdown_per_step_s'], (
            'bucket_breakdown_per_step_s missing PER_STEP_BUCKETS key %r '
            '(regression of the KeyError bug this suite\'s per_rank_metrics '
            'fix addressed)' % key)
    for key in h.FIXED_BUCKETS:
        assert key not in p0['bucket_breakdown_per_step_s'], (
            'bucket_breakdown_per_step_s must not contain FIXED bucket %r '
            '-- those are one-time costs, not a per-step rate' % key)
    assert p0['mean_per_step_s'] > 0
    assert p0['setup_and_io_fixed_s'] > 0
    print('  PASS  cmd_collect + analyze round trip: tarball has profile '
          'files, point 0 declared FAIL with a legible reason (no crash), '
          'point 1 declared incomparable, bucket breakdown has exactly the '
          'PER_STEP_BUCKETS keys')


def main():
    print('Regression guard: testsys/perf/hpc_scaling_suite.py (item 142)')
    checks = [check_bucket_partition,
              check_big_refuses_python_jax_mpi,
              check_per_rank_metrics_fixed_vs_per_step,
              check_fixed_buckets_would_trip_guard_if_differenced,
              check_build_points_strong_weak_size,
              check_collect_and_analyze_round_trip]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append(str(e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_hpc_scaling_suite (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_hpc_scaling_suite')
    return 0


if __name__ == '__main__':
    sys.exit(main())
