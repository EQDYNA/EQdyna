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


def _scaled_row(base, nsteps, per_step_scale, mpi_init=0.5):
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
    s = base['buckets_s']['setup']   # the four setup parts partition the setup bucket exactly
    row['setup_detail_s'] = dict(mpi_init=mpi_init, input_read=0.4 * s, meshing=0.3 * s,
                                 netcdf_open=0.2 * s, mass_assembly=0.1 * s)
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

FAKE_ACCOUNT = 'FAKE-ACCT-7731'


def _write_suite(suite_dir, specs, repeats=2, case='test.tpv8', backend='fortran',
                 nlo=100, nhi=200):
    """A synthetic suite directory + manifest (as --submit writes it) over
    `specs` = [(point_id, ranks, term_s)]. Every run dir carries real-fixture-
    derived profile files and a bGlobal.txt whose lines 15/16 are term and dt."""
    base = _load_fixture_row()
    points = []
    for pid, ranks, term in specs:
        point = dict(point_id=pid, ranks=ranks, nnodes=1, dx=500.0, decomp=[ranks, 1, 1],
                     backend=backend, runs={'lo': [], 'hi': []})
        for k in range(repeats):
            for tag, nsteps in (('lo', nlo), ('hi', nhi)):
                name = h.run_name(tag, k)
                rd = os.path.join(suite_dir, pid, name)
                os.makedirs(rd)
                with open(os.path.join(rd, 'bGlobal.txt'), 'w') as f:
                    f.write('0\n' * 15 + '%r\n0.005\n' % term)
                for r in range(ranks):
                    row = _scaled_row(base, nsteps, 1.0 + 0.1 * k, mpi_init=0.25 * ranks)
                    row['rank'], row['nranks'], row['backend'] = r, ranks, backend
                    with open(os.path.join(rd, 'profile.rank%d.json' % r), 'w') as f:
                        json.dump(row, f)
                point['runs'][tag].append(dict(name=name, run_dir=rd, nsteps=nsteps))
        points.append(point)
    manifest = dict(schema=h.SUITE_SCHEMA_ID, machine='ubuntu', backend=backend,
                    test='strong', level='medium', case=case, max_nodes=1,
                    nsteps_lo=nlo, nsteps_hi=nhi, repeats=repeats, created_utc='test',
                    provenance=dict(git_sha='0' * 40, git_dirty=False, machine='ubuntu',
                                    loaded_modules=[], nsteps_lo=nlo, nsteps_hi=nhi,
                                    repeats=repeats, rank_layout=[], date_utc='test'),
                    points=points)
    with open(os.path.join(suite_dir, h.MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f)
    return manifest


def check_step_pair_cli_defaults_and_validation():
    args = h.build_parser().parse_args(['--submit'])
    assert (args.nsteps_lo, args.nsteps_hi, args.repeats) == (100, 200, 1), (
        'default step pair/repeats are %r' % ((args.nsteps_lo, args.nsteps_hi, args.repeats),))
    args = h.build_parser().parse_args(['--submit', '--nsteps-lo', '7', '--nsteps-hi', '19'])
    assert (args.nsteps_lo, args.nsteps_hi) == (7, 19)
    for lo, hi, rep in ((20, 20, 1), (30, 20, 1), (0, 20, 1), (5, 15, 0)):
        try:
            h._check_steps(lo, hi, rep)
        except h.SuiteError:
            continue
        raise AssertionError('step pair %r repeats %r was not refused' % ((lo, hi), rep))
    print('  PASS  --nsteps-lo/-hi default 100/200, settable, a bad pair or repeats refused')


def check_job_id_is_last_all_digit_line():
    banner = (b'-----------------------------------------------------------------\n'
              b'          Welcome to the Vista Supercomputer, project 12345 here\n'
              b'-> Submission accepted\n12345\n  \n4242\n')
    assert h.parse_job_id(banner) == '4242', 'banner/digit-line output parsed to %r' % h.parse_job_id(banner)
    assert h.parse_job_id(b'Submitted batch job 987\n') == '987'
    assert h.parse_job_id(b'555;vista\n') == '555'
    try:
        h.parse_job_id(b'sbatch: error: invalid account FAKE\n')
    except h.SuiteError as e:
        assert 'FAKE' not in str(e), 'error text echoed sbatch output'
    else:
        raise AssertionError('output with no job id was accepted')
    print('  PASS  parse_job_id takes the last all-digit line, then the last all-digit token, else raises')


class _FakeSbatch:
    """Stands in for subprocess.check_output: `git` passes through to the real
    one; `sbatch` is recorded and answers with a banner that ECHOES the account
    (the worst case for the account-leak guard) before the id."""
    def __init__(self, real):
        self.real, self.sbatch_argv = real, []

    def __call__(self, cmd, *a, **kw):
        if cmd[0] == 'sbatch':
            self.sbatch_argv.append(list(cmd))
            return ('Welcome to Vista. Charging account %s\n7001\n' % FAKE_ACCOUNT).encode()
        return self.real(cmd, *a, **kw)


def _all_text(root):
    out = []
    for dp, _dn, fns in os.walk(root):
        for fn in fns:
            with open(os.path.join(dp, fn), errors='replace') as f:
                out.append((os.path.join(dp, fn), f.read()))
    return out


def check_submit_never_stores_the_account():
    import contextlib
    import io
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_submit_')
    fake = _FakeSbatch(h.subprocess.check_output)
    real_co, real_build = h.subprocess.check_output, h.build_point_case
    h.subprocess.check_output = fake
    h.build_point_case = lambda dst, case, dx, decomp, nsteps, env=None: os.makedirs(dst)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            rc = h.main(['--submit', '--machine', 'vista', '--account', FAKE_ACCOUNT,
                         '--backend', 'fortran', '--test', 'strong', '--level', 'medium',
                         '--max-nodes', '1', '--repeats', '2', '--out-dir', tmpd,
                         '--module', 'netcdf', '--mpirun', 'mpirun'])
    finally:
        h.subprocess.check_output, h.build_point_case = real_co, real_build
    assert rc == 0, 'submit failed rc=%r: %s' % (rc, buf.getvalue())
    assert len(fake.sbatch_argv) == 1 and fake.sbatch_argv[0][:4] == ['sbatch', '--parsable', '-A', FAKE_ACCOUNT], (
        'the account must reach sbatch on its command line: %r' % fake.sbatch_argv)
    leaks = [p for p, t in _all_text(tmpd) if FAKE_ACCOUNT in t]
    assert not leaks, 'account string written into %r' % leaks
    assert FAKE_ACCOUNT not in buf.getvalue(), 'account string echoed on stdout'
    with open(os.path.join(tmpd, h.MANIFEST_NAME)) as f:
        manifest = json.load(f)
    assert 'account' not in manifest
    assert manifest['points'][0]['job_id'] == '7001', 'job_id is %r' % manifest['points'][0]['job_id']
    prov = manifest['provenance']
    assert len(prov['git_sha']) == 40 and isinstance(prov['git_dirty'], bool), prov
    assert (prov['nsteps_lo'], prov['nsteps_hi'], prov['repeats'], prov['machine']) == (100, 200, 2, 'vista')
    assert prov['rank_layout'] == [dict(point_id='strong_n0001', ranks=144, nnodes=1,
                                        decomp=list(h._decompose(144)))], prov['rank_layout']
    assert prov['date_utc'].endswith('Z')
    runs = manifest['points'][0]['runs']
    assert [r['name'] for r in runs['lo']] == ['lo_r0', 'lo_r1'] and [r['nsteps'] for r in runs['hi']] == [200, 200]
    with open(os.path.join(tmpd, 'strong_n0001', 'submit.sh')) as f:
        script = f.read()
    assert script.count('-np 144') == 4, 'two repeats x lo/hi = 4 run lines: %r' % script
    print('  PASS  submit: account only on the sbatch line (not in manifest, script, stdout), '
          'job id parsed from the banner, provenance stamped, repeats built')


def check_vista_dry_run_script():
    import contextlib
    import io
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_dry_')
    argv = ['--submit', '--dry-run', '--machine', 'vista', '--backend', 'fortran',
            '--test', 'strong', '--level', 'medium', '--max-nodes', '1', '--out-dir', tmpd,
            '--module', 'netcdf', '--mpirun', 'mpirun', '--walltime', '00:15:00',
            '--bin', '$SCRATCH/eqdyna-vista/bin/eqdyna', '--work-root', '$SCRATCH/eqdyna-scaling/v1']
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = h.main(argv)
    assert rc == 0, buf.getvalue()
    with open(os.path.join(tmpd, 'strong_n0001', 'submit.sh')) as f:
        script = f.read()
    for want in ('#SBATCH -N 1\n', '#SBATCH -n 144\n', '#SBATCH -p gg\n', '#SBATCH -t 00:15:00\n',
                 'module load netcdf\n',
                 'cd $SCRATCH/eqdyna-scaling/v1/strong_n0001/lo_r0 && mpirun -np 144 $SCRATCH/eqdyna-vista/bin/eqdyna\n',
                 'cd $SCRATCH/eqdyna-scaling/v1/strong_n0001/hi_r0 && mpirun -np 144 $SCRATCH/eqdyna-vista/bin/eqdyna\n'):
        assert want in script, 'job script lacks %r:\n%s' % (want, script)
    assert '#SBATCH -A' not in script
    assert 'sbatch -A <ALLOC> %s' % os.path.join(tmpd, 'strong_n0001', 'submit.sh') in buf.getvalue()
    assert [n for n in os.listdir(os.path.join(tmpd, 'strong_n0001'))] == ['submit.sh'], 'dry run built a case'
    assert not os.path.exists(os.path.join(tmpd, h.MANIFEST_NAME)), 'dry run wrote a manifest'
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = h.main(argv + ['--account', FAKE_ACCOUNT])
    assert rc == 1 and FAKE_ACCOUNT not in buf.getvalue(), 'dry run with --account was not refused quietly'
    print('  PASS  vista dry run: one-node script (gg, 144 ranks, 00:15:00, module load netcdf, '
          '$SCRATCH paths, no -A), no case built, --account refused')


def check_per_step_from_loop_s_and_nonpositive_raises():
    base = _load_fixture_row()
    lo, hi = _scaled_row(base, 100, 1.0), _scaled_row(base, 200, 1.0)
    lo['loop_s'], hi['loop_s'] = 10.0, 20.0
    lo['total_s'], hi['total_s'] = 50.0, 40.0     # a wrapper-style total that does NOT difference
    r = h.per_rank_metrics([lo], [hi])[0]
    assert abs(r['per_step_s'] - 0.1) < 1e-12, 'per-step must difference loop_s (0.1 s), got %r' % r['per_step_s']
    hi['loop_s'] = 10.0
    try:
        h.per_rank_metrics([lo], [hi])
    except numa.NonPositiveMeasurement:
        pass
    else:
        raise AssertionError('a non-positive loop_s difference was reported as a number')
    print('  PASS  per-step = difference of each rank\'s own loop_s; non-positive raises')


def check_overhead_and_setup_split_required():
    base = _load_fixture_row()
    lo, hi = _scaled_row(base, 100, 1.0), _scaled_row(base, 200, 1.0, mpi_init=1.5)
    r = h.per_rank_metrics([lo], [hi])[0]
    want = 1.5 + base['buckets_s']['setup'] + base['buckets_s']['io']
    assert abs(r['overhead']['total_s'] - want) < 1e-9, 'overhead %r != init+setup+io %r' % (r['overhead'], want)
    assert r['overhead']['init_s'] == 1.5 and set(r['setup_detail_s']) == set(h.DETAIL_KEYS['fortran'])
    missing = copy.deepcopy(hi)
    del missing['setup_detail_s']
    short = copy.deepcopy(hi)
    short['setup_detail_s']['mass_assembly'] += 1.0 + base['buckets_s']['setup']   # parts exceed setup
    for bad, why in ((missing, 'no setup_detail_s'), (short, 'parts do not sum to setup')):
        try:
            h.per_rank_metrics([lo], [bad])
        except h.SuiteError:
            continue
        raise AssertionError('%s was accepted' % why)
    print('  PASS  overhead = mpi_init + setup + io per rank; a missing or non-summing setup split raises')


def check_profile_schema_setup_detail_shape():
    row = _scaled_row(_load_fixture_row(), 100, 1.0)
    h.profile_schema.validate(row)
    for bad in ({}, {'mpi_init': -1.0}, {'mpi_init': 'x'}, [1.0]):
        row2 = dict(row, setup_detail_s=bad)
        try:
            h.profile_schema.validate(row2)
        except ValueError:
            continue
        raise AssertionError('setup_detail_s %r passed the schema' % (bad,))
    print('  PASS  profile_schema: setup_detail_s optional, shape checked')


def check_collect_analyze_overhead_repeats_parity():
    """cmd_collect -> analyze on two points (2 and 4 ranks) x 2 repeats. Point 0
    runs at term 0.5 s != the reference's 5 s -> DECLARED with the measured
    mismatch; point 1 claims the reference term and mesh, so compare_cell runs
    (stubbed to FAIL) -> the point is flagged and --analyze exits non-zero."""
    import argparse
    import contextlib
    import io
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_suite_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    _write_suite(suite_dir, [('strong_n0001', 2, 0.5), ('strong_n0002', 4, 5.0)])
    real_cc = h.compare.compare_cell
    h.compare.compare_cell = lambda case, backend, run_dir: (False, ['frt: STUB MISMATCH'])
    try:
        h.cmd_collect(argparse.Namespace(suite_dir=suite_dir, out=None))
    finally:
        h.compare.compare_cell = real_cc
    tgz = os.path.join(suite_dir, 'hpc_scaling_ubuntu_fortran_strong_medium.tgz')
    with tarfile.open(tgz, 'r:gz') as tar:
        names = tar.getnames()
    assert 'strong_n0002/hi_r1/profile.rank3.json' in names, 'repeat run dirs not packed: %r' % names[:6]
    report = h.analyze(tgz)
    p0, p1 = report['points']
    assert report['provenance']['git_sha'] == '0' * 40
    assert p0['parity_status'] == 'declared' and 'timed term 0.5' in ' '.join(p0['parity_lines']), p0['parity_lines']
    assert p1['parity_status'] == 'fail' and 'STUB MISMATCH' in ' '.join(p1['parity_lines']), p1['parity_lines']
    assert abs(p0['overhead_s']['init']['mean'] - 0.5) < 1e-9 and abs(p1['overhead_s']['init']['mean'] - 1.0) < 1e-9, (
        'overhead init per rank count: %r / %r' % (p0['overhead_s']['init'], p1['overhead_s']['init']))
    assert set(p0['setup_detail_mean_s']) == set(h.DETAIL_KEYS['fortran'])
    for p in (p0, p1):
        assert p['repeats'] == 2 and len(p['repeat_per_step_s']) == 2
        assert p['repeat_rel_range'] > 0.05, 'repeats scaled 1.0/1.1 must show spread: %r' % p['repeat_rel_range']
        assert p['overhead_s']['total']['max'] >= p['overhead_s']['total']['mean'] >= p['overhead_s']['total']['min'] > 0
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = h.cmd_analyze(argparse.Namespace(tgz=tgz, out=None))
    assert rc == 1 and 'strong_n0002' in buf.getvalue().splitlines()[-1], 'a failed parity point did not fail --analyze'
    try:
        h.point_parity(dict(point_id='x', runs={'hi': [dict(name='hi_r0')]}))
    except h.SuiteError:
        pass
    else:
        raise AssertionError('a timed run with no parity record was accepted')
    print('  PASS  collect+analyze: repeats packed, per-rank-count overhead and spread reported, '
          'declared vs failed parity, failed parity exits non-zero, unchecked run refused')


def main():
    print('Regression guard: testsys/perf/hpc_scaling_suite.py (item 142)')
    checks = [check_bucket_partition,
              check_big_refuses_python_jax_mpi,
              check_per_rank_metrics_fixed_vs_per_step,
              check_fixed_buckets_would_trip_guard_if_differenced,
              check_build_points_strong_weak_size,
              check_step_pair_cli_defaults_and_validation,
              check_job_id_is_last_all_digit_line,
              check_submit_never_stores_the_account,
              check_vista_dry_run_script,
              check_per_step_from_loop_s_and_nonpositive_raises,
              check_overhead_and_setup_split_required,
              check_profile_schema_setup_detail_shape,
              check_collect_analyze_overhead_repeats_parity]
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
