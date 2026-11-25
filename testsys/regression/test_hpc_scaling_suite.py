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
owner's spec); `--analyze` on a python-jax-mpi suite (fault/wait declared
unmeasured, reported "not measured", never differenced); a failing sbatch is a
SuiteError that never carries the allocation; the job id is the last all-digit
line or an error; and the PARITY block: one parity point per suite compared
with compare_cell, every timed point declared, `--analyze` non-zero unless the
parity point passed.

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


def _scaled_row(base, nsteps, per_step_scale, mpi_init=0.5, backend='fortran'):
    """A synthetic row derived from the real fixture: PER_STEP_BUCKETS scaled
    by (nsteps/base_nsteps)*per_step_scale (so lo/hi differ only in step
    count, as a real run would), FIXED_BUCKETS held at the fixture's own
    value (a one-time cost does not grow with nsteps). backend='python-jax-mpi'
    reshapes it the way that emitter writes: fault and wait are 0.0 (declared
    unmeasured), the setup split has the python keys and no mpi_init."""
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
    s = base['buckets_s']['setup']   # the setup parts partition the setup bucket exactly
    if backend == 'fortran':
        row['setup_detail_s'] = dict(mpi_init=mpi_init, input_read=0.4 * s, meshing=0.3 * s,
                                     netcdf_open=0.2 * s, mass_assembly=0.1 * s)
    else:
        row['backend'] = backend
        buckets['fault'] = 0.0
        buckets['wait'] = 0.0
        row['loop_s'] = buckets['element'] + buckets['exchange']
        row['total_s'] = row['loop_s'] + fixed_total
        row['setup_detail_s'] = dict(input_mesh=0.5 * s, setup_exchange=0.3 * s, driver_setup=0.2 * s)
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
                 nlo=100, nhi=200, parity='pass'):
    """A synthetic suite directory + manifest (as --submit writes it) over
    `specs` = [(point_id, ranks)]. Every run dir carries real-fixture-derived
    profile files. `parity` = 'absent' writes no parity_point; anything else
    writes one (its profile files under parity/gate, no result yet: --collect
    adds it)."""
    base = _load_fixture_row()
    points = []
    for pid, ranks in specs:
        point = dict(point_id=pid, ranks=ranks, nnodes=1, dx=500.0, decomp=[ranks, 1, 1],
                     backend=backend, parity_declaration=h.TIMED_PARITY_DECLARATION,
                     runs={'lo': [], 'hi': []})
        for k in range(repeats):
            for tag, nsteps in (('lo', nlo), ('hi', nhi)):
                name = h.run_name(tag, k)
                rd = os.path.join(suite_dir, pid, name)
                os.makedirs(rd)
                for r in range(ranks):
                    row = _scaled_row(base, nsteps, 1.0 + 0.1 * k, mpi_init=0.25 * ranks,
                                      backend=backend)
                    row['rank'], row['nranks'], row['backend'] = r, ranks, backend
                    with open(os.path.join(rd, 'profile.rank%d.json' % r), 'w') as f:
                        json.dump(row, f)
                point['runs'][tag].append(dict(name=name, run_dir=rd, nsteps=nsteps))
        points.append(point)
    manifest = dict(schema=h.SUITE_SCHEMA_ID, machine='ubuntu', backend=backend,
                    test='strong', level='medium', case=case, max_nodes=1,
                    nsteps_lo=nlo, nsteps_hi=nhi, repeats=repeats, created_utc='test',
                    provenance=dict(git_sha='0' * 40, git_dirty=False, solver_git_sha='0' * 40,
                                    machine='ubuntu', loaded_modules=[], nsteps_lo=nlo,
                                    nsteps_hi=nhi, repeats=repeats, rank_layout=[],
                                    date_utc='test'),
                    points=points)
    if parity != 'absent':
        rd = os.path.join(suite_dir, 'parity', 'gate')
        os.makedirs(rd)
        for r in range(4):
            row = _scaled_row(base, 50, 1.0, backend=backend)
            row['rank'], row['nranks'], row['backend'] = r, 4, backend
            with open(os.path.join(rd, 'profile.rank%d.json' % r), 'w') as f:
                json.dump(row, f)
        manifest['parity_point'] = dict(point_id='parity', ranks=4, nnodes=1, dx=500.0,
                                        decomp=[2, 2, 1], term_s=5.0, backend=backend,
                                        host_point=specs[0][0],
                                        run=dict(name='gate', run_dir=rd))
    with open(os.path.join(suite_dir, h.MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f)
    return manifest


def _collect(suite_dir, ok=True, raises=None):
    """cmd_collect with compare_cell stubbed; returns the list of
    (case, backend, run_dir) it was called with."""
    import argparse
    calls = []

    def stub(case, backend, run_dir):
        calls.append((case, backend, run_dir))
        if raises:
            raise raises
        return ok, ([] if ok else ['frt: STUB MISMATCH']) + ['frt: STUB COMPARED']
    real = h.compare.compare_cell
    h.compare.compare_cell = stub
    try:
        h.cmd_collect(argparse.Namespace(suite_dir=suite_dir, out=None))
    finally:
        h.compare.compare_cell = real
    return calls


def _tar_suite(suite_dir, tgz):
    """Pack manifest + profile files as they are (no --collect)."""
    with tarfile.open(tgz, 'w:gz') as tar:
        tar.add(os.path.join(suite_dir, h.MANIFEST_NAME), arcname=h.MANIFEST_NAME)
        for dp, _dn, fns in os.walk(suite_dir):
            for fn in sorted(fns):
                if fn.startswith('profile.rank'):
                    full = os.path.join(dp, fn)
                    tar.add(full, arcname=os.path.relpath(full, suite_dir))




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


def _all_text(root):
    out = []
    for dp, _dn, fns in os.walk(root):
        for fn in fns:
            with open(os.path.join(dp, fn), errors='replace') as f:
                out.append((os.path.join(dp, fn), f.read()))
    return out


def check_job_id_is_last_all_digit_line_and_nothing_else():
    banner = (b'-----------------------------------------------------------------\n'
              b'          Welcome to the Vista Supercomputer, project 12345 here\n'
              b'-> Submission accepted\n12345\n  \n4242\n')
    assert h.parse_job_id(banner) == '4242', 'banner/digit-line output parsed to %r' % h.parse_job_id(banner)
    assert h.parse_job_id(b'7\n') == '7'
    for bad in (b'Welcome to Vista, valid through 2026 for this allocation\n',
                b'Submitted batch job 987\n', b'555;vista\n', b'', b'12 34\n'):
        try:
            h.parse_job_id(bad)
        except h.SuiteError as e:
            assert bad.decode().strip() in str(e) or bad == b'', 'error does not quote the output: %s' % e
        else:
            raise AssertionError('output %r with no all-digit line was accepted' % bad)
    try:
        h.parse_job_id(('acct %s nope\n' % FAKE_ACCOUNT).encode(), redact=FAKE_ACCOUNT)
    except h.SuiteError as e:
        assert FAKE_ACCOUNT not in str(e) and h.ALLOC_PLACEHOLDER in str(e), str(e)
    else:
        raise AssertionError('a no-id output was accepted')
    print('  PASS  parse_job_id: last all-digit line only, no word fallback, error quotes the output '
          '(allocation redacted)')


def _fake_sbatch_dir(root, mode):
    """A directory with an executable `sbatch` to put first on PATH. Modes:
    ok (banner echoing the account, then the id), fail (exit 3, stderr echoes
    the account), banner (exit 0, a sentence containing 2026, no id line)."""
    d = os.path.join(root, 'fakebin_' + mode)
    os.makedirs(d)
    body = {
        'ok': 'echo "Welcome. Charging account %s"\necho 7001\n' % FAKE_ACCOUNT,
        'fail': 'echo "sbatch: error: invalid account %s" >&2\nexit 3\n' % FAKE_ACCOUNT,
        'banner': 'echo "Welcome to Vista, valid through 2026, account %s"\n' % FAKE_ACCOUNT,
    }[mode]
    path = os.path.join(d, 'sbatch')
    with open(path, 'w') as f:
        f.write('#! /bin/bash\n' + body)
    os.chmod(path, 0o755)
    return d


def _run_submit(tmpd, mode, extra=(), default_bin=False):
    """h.main --submit on vista with a fake sbatch (a real executable on PATH,
    so the real subprocess code runs) and the case builders stubbed to mkdir.
    Returns (rc, stdout, stderr)."""
    import contextlib
    import io
    fake = _fake_sbatch_dir(tmpd, mode)
    real_build, real_pbuild, real_path = h.build_point_case, h.build_parity_case, os.environ['PATH']
    h.build_point_case = lambda dst, case, dx, decomp, nsteps, env=None: os.makedirs(dst)
    h.build_parity_case = lambda dst, case, decomp, term, env=None: os.makedirs(dst)
    os.environ['PATH'] = fake + os.pathsep + real_path
    out, err = io.StringIO(), io.StringIO()
    argv = ['--submit', '--machine', 'vista', '--account', FAKE_ACCOUNT,
            '--backend', 'fortran', '--test', 'strong', '--level', 'medium',
            '--max-nodes', '1', '--repeats', '2', '--out-dir', os.path.join(tmpd, 'suite'),
            '--module', 'netcdf', '--mpirun', 'mpirun'] + list(extra)
    if not default_bin:
        argv += ['--bin', os.path.join(tmpd, 'nogit', 'bin', 'eqdyna')]
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = h.main(argv)
    finally:
        h.build_point_case, h.build_parity_case, os.environ['PATH'] = real_build, real_pbuild, real_path
    return rc, out.getvalue(), err.getvalue()


def check_submit_never_stores_the_account():
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_submit_')
    os.makedirs(os.path.join(tmpd, 'nogit', 'bin'))
    rc, out, err = _run_submit(tmpd, 'ok')
    assert rc == 0, 'submit failed rc=%r: %s %s' % (rc, out, err)
    suite = os.path.join(tmpd, 'suite')
    leaks = [p for p, t in _all_text(suite) if FAKE_ACCOUNT in t]
    assert not leaks, 'account string written into %r' % leaks
    assert FAKE_ACCOUNT not in out + err, 'account string echoed on stdout/stderr'
    with open(os.path.join(suite, h.MANIFEST_NAME)) as f:
        manifest = json.load(f)
    assert 'account' not in manifest
    assert manifest['points'][0]['job_id'] == '7001', 'job_id is %r' % manifest['points'][0]['job_id']
    prov = manifest['provenance']
    assert len(prov['git_sha']) == 40 and isinstance(prov['git_dirty'], bool), prov
    assert prov['solver_git_sha'] == h.UNKNOWN_BIN_SHA == 'unknown (bin not in a git tree)', prov['solver_git_sha']
    assert (prov['nsteps_lo'], prov['nsteps_hi'], prov['repeats'], prov['machine']) == (100, 200, 2, 'vista')
    assert prov['rank_layout'] == [dict(point_id='strong_n0001', ranks=144, nnodes=1,
                                        decomp=list(h._decompose(144)), dx=prov['rank_layout'][0]['dx'])], prov['rank_layout']
    assert prov['date_utc'].endswith('Z')
    runs = manifest['points'][0]['runs']
    assert [r['name'] for r in runs['lo']] == ['lo_r0', 'lo_r1'] and [r['nsteps'] for r in runs['hi']] == [200, 200]
    pp = manifest['parity_point']
    assert (pp['ranks'], pp['decomp'], pp['term_s'], pp['dx'], pp['host_point']) == (
        4, [2, 2, 1], 5.0, 500.0, 'strong_n0001'), pp
    assert all(p['parity_declaration'] == 'declared: dx/term differ from the reference'
               for p in manifest['points'])
    with open(os.path.join(suite, 'strong_n0001', 'submit.sh')) as f:
        script = f.read()
    assert script.count('-np 144') == 4, 'two repeats x lo/hi = 4 run lines: %r' % script
    assert 'cd %s/parity/gate && mpirun -np 4 ' % suite in script, 'parity run not in the first job: %r' % script
    print('  PASS  submit: account only on the sbatch line (not in manifest, script, stdout, stderr), '
          'job id parsed from the banner, provenance + parity point stamped, parity run in the same job')


def check_solver_sha_of_the_bin():
    sha = h._git('rev-parse', 'HEAD')
    got = h.solver_git_sha('fortran', os.path.join(TESTSYS, 'perf', 'eqdyna'))   # a path inside this git tree
    assert got == sha, 'bin inside the tree: %r != HEAD %r' % (got, sha)
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_sha_')
    assert h.solver_git_sha('fortran', os.path.join(tmpd, 'eqdyna')) == h.UNKNOWN_BIN_SHA
    assert h.solver_git_sha('fortran', '$SCRATCH_UNSET_XYZ/bin/eqdyna') == h.UNKNOWN_BIN_SHA
    print('  PASS  solver_git_sha: the bin\'s tree HEAD, else "unknown (bin not in a git tree)"')


def check_failed_sbatch_raises_suite_error_without_account():
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_sbfail_')
    os.makedirs(os.path.join(tmpd, 'nogit', 'bin'))
    rc, out, err = _run_submit(tmpd, 'fail')
    suite = os.path.join(tmpd, 'suite')
    assert rc == 1, 'a failing sbatch must exit 1 (SuiteError), got %r: %s %s' % (rc, out, err)
    assert 'exited 3' in out and h.ALLOC_PLACEHOLDER in out, 'exit code / redacted stderr missing: %r' % out
    leaks = [p for p, t in _all_text(tmpd) if FAKE_ACCOUNT in t and 'fakebin_' not in p]
    assert FAKE_ACCOUNT not in out + err and not leaks, 'account leaked: out/err %r, files %r' % (out + err, leaks)
    with open(os.path.join(suite, h.MANIFEST_NAME)) as f:
        manifest = json.load(f)
    assert manifest['points'][0]['job_id'] is None, 'manifest must exist BEFORE any submit (job_id still null)'
    print('  PASS  sbatch failure: SuiteError with exit code + redacted stderr, no argv, no account anywhere; '
          'manifest already written')


def check_banner_only_sbatch_reply_raises():
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_banner_')
    os.makedirs(os.path.join(tmpd, 'nogit', 'bin'))
    rc, out, err = _run_submit(tmpd, 'banner')
    assert rc == 1 and 'no all-digit job id line' in out and 'valid through 2026' in out, (
        'a banner-only reply was not refused with the output quoted: rc=%r %r' % (rc, out))
    assert FAKE_ACCOUNT not in out + err, 'account echoed from the banner'
    print('  PASS  banner-only sbatch reply (a sentence containing 2026) raises, output quoted, account redacted')


def check_work_root_cases_built_where_the_job_cds():
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_wr_')
    os.makedirs(os.path.join(tmpd, 'nogit', 'bin'))
    wr = os.path.join(tmpd, 'workroot')
    rc, out, err = _run_submit(tmpd, 'ok', extra=['--work-root', wr])
    assert rc == 0, out + err
    suite = os.path.join(tmpd, 'suite')
    with open(os.path.join(suite, h.MANIFEST_NAME)) as f:
        manifest = json.load(f)
    run_dir = manifest['points'][0]['runs']['lo'][0]['run_dir']
    assert run_dir == os.path.join(wr, 'strong_n0001', 'lo_r0') and os.path.isdir(run_dir), run_dir
    assert manifest['parity_point']['run']['run_dir'] == os.path.join(wr, 'parity', 'gate')
    with open(os.path.join(suite, 'strong_n0001', 'submit.sh')) as f:
        assert 'cd %s && mpirun -np 144' % run_dir in f.read(), 'the job does not cd where the case was built'
    tmp2 = tempfile.mkdtemp(prefix='test_hpc_scaling_wr2_')
    os.makedirs(os.path.join(tmp2, 'nogit', 'bin'))
    rc, out, err = _run_submit(tmp2, 'ok', extra=['--work-root', '$NO_SUCH_VAR_HPC_142/x'])
    assert rc == 1 and 'does not resolve' in out, 'an unresolvable --work-root was accepted: %r' % out
    print('  PASS  --work-root: cases built (and recorded) where the job cds; an unresolved variable refused')




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
                 'module load netcdf || exit 1\n',
                 'cd $SCRATCH/eqdyna-scaling/v1/strong_n0001/lo_r0 && mpirun -np 144 $SCRATCH/eqdyna-vista/bin/eqdyna\n',
                 'cd $SCRATCH/eqdyna-scaling/v1/strong_n0001/hi_r0 && mpirun -np 144 $SCRATCH/eqdyna-vista/bin/eqdyna\n',
                 'cd $SCRATCH/eqdyna-scaling/v1/parity/gate && mpirun -np 4 $SCRATCH/eqdyna-vista/bin/eqdyna\n'):
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


def check_collect_analyze_overhead_repeats_and_parity_point():
    """cmd_collect -> analyze on two points (2 and 4 ranks) x 2 repeats plus
    the parity point: repeats packed, per-rank-count overhead and spread
    reported, compare_cell called ONCE (the parity run, never a timed run),
    every timed point declared."""
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_suite_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    _write_suite(suite_dir, [('strong_n0001', 2), ('strong_n0002', 4)])
    calls = _collect(suite_dir, ok=True)
    assert calls == [('test.tpv8', 'fortran', os.path.join(suite_dir, 'parity', 'gate'))], (
        'compare_cell must run on exactly the parity run: %r' % calls)
    tgz = os.path.join(suite_dir, 'hpc_scaling_ubuntu_fortran_strong_medium.tgz')
    with tarfile.open(tgz, 'r:gz') as tar:
        names = tar.getnames()
    assert 'strong_n0002/hi_r1/profile.rank3.json' in names and 'parity/gate/profile.rank3.json' in names, names[:8]
    report = h.analyze(tgz)
    p0, p1 = report['points']
    assert report['provenance']['git_sha'] == '0' * 40
    assert abs(p0['overhead_s']['init']['mean'] - 0.5) < 1e-9 and abs(p1['overhead_s']['init']['mean'] - 1.0) < 1e-9, (
        'overhead init per rank count: %r / %r' % (p0['overhead_s']['init'], p1['overhead_s']['init']))
    assert set(p0['setup_detail_mean_s']) == set(h.DETAIL_KEYS['fortran'])
    for p in (p0, p1):
        assert p['repeats'] == 2 and len(p['repeat_per_step_s']) == 2
        assert p['repeat_rel_range'] > 0.05, 'repeats scaled 1.0/1.1 must show spread: %r' % p['repeat_rel_range']
        assert p['overhead_s']['total']['max'] >= p['overhead_s']['total']['mean'] >= p['overhead_s']['total']['min'] > 0
        assert p['parity_declaration'] == 'declared: dx/term differ from the reference'
        assert p['buckets_not_measured'] == [], 'fortran measures every bucket: %r' % p['buckets_not_measured']
    print('  PASS  collect+analyze: repeats packed, overhead/spread per rank count, compare_cell on the '
          'parity run only, every timed point declared')


def check_analyze_python_jax_mpi_suite():
    """M1: --analyze on a synthetic python-jax-mpi suite (fault=0.0, wait=0.0
    in every profile, setup split input_mesh/setup_exchange/driver_setup)."""
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_pyjax_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    _write_suite(suite_dir, [('strong_n0001', 4), ('strong_n0002', 8)], repeats=1,
                 backend='python-jax-mpi')
    calls = _collect(suite_dir, ok=True)
    assert calls[0][1] == 'python-jax', 'parity must compare against the python-jax cell: %r' % calls
    tgz = os.path.join(suite_dir, 'hpc_scaling_ubuntu_python-jax-mpi_strong_medium.tgz')
    report = h.analyze(tgz)
    for p in report['points']:
        assert p['mean_per_step_s'] > 0, 'per-step must come from loop_s: %r' % p['mean_per_step_s']
        assert p['buckets_not_measured'] == ['fault', 'wait'], p['buckets_not_measured']
        assert p['bucket_breakdown_per_step_s']['fault'] is None and p['bucket_breakdown_per_step_s']['wait'] is None
        assert p['bucket_breakdown_per_step_s']['element'] > 0 and p['bucket_breakdown_per_step_s']['exchange'] > 0
        assert p['overhead_s']['init'] is None, 'python-jax-mpi mpi_init is declared not measured'
        assert set(p['setup_detail_mean_s']) == {'input_mesh', 'setup_exchange', 'driver_setup'}
        assert abs(sum(p['setup_detail_mean_s'].values()) - p['overhead_s']['setup']['mean']) < 1e-9
    import argparse
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = h.cmd_analyze(argparse.Namespace(tgz=tgz, out=None))
    assert rc == 0, 'a clean python-jax-mpi suite must analyse green: %s' % buf.getvalue()[-400:]
    print('  PASS  --analyze on a python-jax-mpi suite: fault/wait reported "not measured", per-step from '
          'loop_s, setup split incl. driver_setup sums to setup')


def check_unmeasured_bucket_only_while_zero_and_setup_sum_two_sided():
    base = _load_fixture_row()
    lo = _scaled_row(base, 100, 1.0, backend='python-jax-mpi')
    hi = _scaled_row(base, 200, 1.0, backend='python-jax-mpi')
    lo['buckets_s']['wait'] = 0.5    # a nonzero wait is a measurement: it is differenced (0.5 -> 0.0 is non-positive) -> below noise, NOT "not measured" (None)
    got = h.per_rank_metrics([lo], [hi])[0]['buckets_per_step_s']['wait']
    assert got == h.BELOW_NOISE, 'a nonzero wait was waved through as %r' % (got,)
    s = base['buckets_s']['setup']
    for backend in ('fortran', 'python-jax-mpi'):
        lo = _scaled_row(base, 100, 1.0, backend=backend)
        hi = _scaled_row(base, 200, 1.0, backend=backend)
        hi['buckets_s']['setup'] += 1.0 + s     # parts now SHORT of setup: the sum check is two-sided
        try:
            h.per_rank_metrics([lo], [hi])
        except h.SuiteError:
            continue
        raise AssertionError('%s: setup parts short of setup were accepted' % backend)
    print('  PASS  unmeasured only while zero in both runs; setup-split sum is two-sided for both backends')


def check_parity_block_verdicts_and_exit_codes():
    import argparse
    import contextlib
    import io

    def analyse(parity, ok=True, raises=None, collect=True):
        tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_par_')
        suite_dir = os.path.join(tmpd, 'suite')
        os.makedirs(suite_dir)
        _write_suite(suite_dir, [('strong_n0001', 2), ('strong_n0002', 4)], repeats=1, parity=parity)
        if collect:
            _collect(suite_dir, ok=ok, raises=raises)
            tgz = os.path.join(suite_dir, 'hpc_scaling_ubuntu_fortran_strong_medium.tgz')
        else:
            tgz = os.path.join(tmpd, 'raw.tgz')
            _tar_suite(suite_dir, tgz)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = h.cmd_analyze(argparse.Namespace(tgz=tgz, out=None))
        return rc, buf.getvalue()

    rc, text = analyse('pass')
    assert rc == 0 and 'PARITY: PASS -- compared against the reference: 1 point(s); declared not comparable: 2' in text, text[-700:]
    for pid in ('strong_n0001', 'strong_n0002'):
        assert '  %s: declared: dx/term differ from the reference' % pid in text, 'declaration not printed for %s' % pid
    rc, text = analyse('pass', ok=False)
    assert rc == 1 and 'PARITY: FAIL' in text and 'STUB MISMATCH' in text, 'a failed parity point did not fail --analyze'
    rc, text = analyse('pass', raises=RuntimeError('no frt'))
    assert rc == 1 and 'PARITY: FAIL' in text and 'no frt' in text, 'a compare_cell that raises must be a failed parity point'
    rc, text = analyse('absent', collect=False)
    assert rc == 1 and 'PARITY: ABSENT -- compared against the reference: 0 point(s); declared not comparable: 2' in text, text[-500:]
    rc, text = analyse('missing', collect=False)
    assert rc == 1 and 'PARITY: MISSING' in text and 'compared against the reference: 0 point(s)' in text, text[-500:]

    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_par2_')
    os.makedirs(os.path.join(tmpd, 's'))
    _write_suite(os.path.join(tmpd, 's'), [('strong_n0001', 2)], repeats=1, parity='absent')
    try:
        _collect(os.path.join(tmpd, 's'))
    except h.SuiteError:
        pass
    else:
        raise AssertionError('--collect accepted a manifest with no parity point')
    try:
        h.parity_block(dict(points=[dict(point_id='x')], parity_point=None))
    except h.SuiteError:
        pass
    else:
        raise AssertionError('a timed point without its declaration label was accepted')
    print('  PASS  PARITY block: pass exits 0; fail, raising compare, absent and missing exit 1; '
          'counts and per-point declarations printed; collect refuses a suite with no parity point')


# --------------------------------------------------------------------------
# Row 142 Vista smoke findings: dx must mesh, job script, netCDF4, installer
# --------------------------------------------------------------------------
RAW_VISTA_DX = 170.30165587605515   # the dx the first Vista smoke wrote; solver exit 33 on every rank
LEVEL_MAX_NODES = {'medium': 16, 'large': 64, 'big': 256}
ELEMENT_ERROR_BOUND = 0.015         # measured max 0.98% (large); bound = 1.5x margin


def _point_targets(test, level, points, cores):
    tgt = h.LEVEL_TARGET_ELEMENTS[level]
    if test == 'strong':
        return [tgt for _ in points]
    if test == 'weak':
        return [tgt * p['ranks'] / float(cores) for p in points]
    return [tgt * f for f in (0.1, 0.3, 1.0, 3.0, 10.0)]


def check_every_point_dx_meshes_and_stays_near_target():
    cfg = dict(cores_per_node=144)
    for level, max_nodes in LEVEL_MAX_NODES.items():
        for test in ('strong', 'weak', 'size'):
            points = h.build_points(test, level, cfg, max_nodes)
            for p, target in zip(points, _point_targets(test, level, points, 144)):
                try:
                    h.grid_line_check(h.DEFAULT_CASE, p['dx'])
                except h.SuiteError as e:
                    raise AssertionError('%s/%s point %s dx=%r does not mesh: %s'
                                         % (level, test, p['point_id'], p['dx'], e))
                err = abs(h.elements_for_dx(p['dx']) / target - 1.0)
                assert err < ELEMENT_ERROR_BOUND, (
                    '%s/%s point %s: snapped dx %r is %.2f%% off the target element count '
                    '(bound %.1f%%)' % (level, test, p['point_id'], p['dx'], 100 * err,
                                        100 * ELEMENT_ERROR_BOUND))
    try:
        h.grid_line_check(h.DEFAULT_CASE, RAW_VISTA_DX)
    except h.SuiteError as e:
        assert 'not an integer multiple' in str(e), e
    else:
        raise AssertionError('the solver rule accepted the raw dx %r that Vista refused: the '
                             'check cannot fail' % RAW_VISTA_DX)
    assert h.dx_for_target_elements(h.LEVEL_TARGET_ELEMENTS['medium']) == RAW_VISTA_DX, (
        'RAW_VISTA_DX no longer is the unsnapped medium dx; update this test')
    print('  PASS  every point of every test x level (%s) meshes under the solver rule; snapped dx '
          'stays within %.1f%% of target; the raw Vista dx is refused' % (
              ', '.join(LEVEL_MAX_NODES), 100 * ELEMENT_ERROR_BOUND))


def check_medium_smallest_points_build_and_mesh_lines():
    """create.newcase + case.setup (build_point_case) for the smallest point of each
    medium test type, then the Python setup/consistency check on the BUILT case (the
    solver's per-axis belt rule on the files the solver will read). large/big smallest
    points are 1e9+ elements: their dx is covered by the pure check above, not built."""
    from testsys.common import mpi_child_env
    env = mpi_child_env(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([os.path.join(REPO_ROOT, 'bin'), os.path.join(REPO_ROOT, 'scripts'),
                                   env.get('PATH', '')])
    cfg = dict(cores_per_node=144)
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_build_')
    for test in ('strong', 'weak', 'size'):
        pts = h.build_points(test, 'medium', cfg, 1)
        p = pts[0]
        try:
            case = h.build_point_case(os.path.join(tmpd, p['point_id']), h.DEFAULT_CASE,
                                      p['dx'], p['decomp'], 100, env)
        except (h.SuiteError, RuntimeError) as e:
            raise AssertionError('medium/%s smallest point %s (dx=%r) did not build/mesh: %s'
                                 % (test, p['point_id'], p['dx'], e))
        with open(os.path.join(case, 'bModelGeometry.txt')) as f:
            assert repr(p['dx']) in f.read(), 'built case does not carry the point dx %r' % p['dx']
    print('  PASS  medium smallest points (strong/weak/size) build and pass the solver axis rule')


def check_job_script_parity_first_no_errexit_and_rc_recorded():
    import subprocess
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_job_')
    d1, d2, d3 = (os.path.join(tmpd, n) for n in ('parity', 'lo', 'hi'))
    for d in (d1, d2, d3):
        os.makedirs(d)
    path = os.path.join(tmpd, 'submit.sh')
    h.write_job_script(path, jobname='j', nnodes=1, ranks=4, partition='gg', walltime='00:05:00',
                       modules=[], commands=['cd %s && true' % d1, 'cd %s && false' % d2,
                                             'cd %s && true' % d3], host_label='ubuntu')
    with open(path) as f:
        script = f.read()
    assert 'set -e' not in script and 'set -euo' not in script, 'errexit hides every point after a failure'
    r = subprocess.run(['bash', path], capture_output=True, text=True)
    out = r.stdout.splitlines()
    want = ['%s %s 0' % (h.POINT_RC_TAG, d1), '%s %s 1' % (h.POINT_RC_TAG, d2),
            '%s %s 0' % (h.POINT_RC_TAG, d3)]
    assert out == want, 'per-point rc lines %r != %r (stderr %r)' % (out, want, r.stderr)
    assert r.returncode == 1, 'a failed point must fail the job, got rc %d' % r.returncode
    # order: the suite's first job runs the parity line before the timed lines
    tmpd2 = tempfile.mkdtemp(prefix='test_hpc_scaling_order_')
    argv = ['--submit', '--dry-run', '--machine', 'vista', '--backend', 'fortran', '--test', 'strong',
            '--level', 'medium', '--max-nodes', '1', '--out-dir', tmpd2, '--bin', '/b/eqdyna']
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        assert h.main(argv) == 0
    with open(os.path.join(tmpd2, 'strong_n0001', 'submit.sh')) as f:
        lines = f.read().splitlines()
    runs = [ln for ln in lines if ln.startswith('cd ')]
    assert '/parity/gate && ' in runs[0], 'parity point is not the first run line: %r' % runs
    assert len(runs) == 3, runs
    assert 'mpirun -np 144' in runs[1], runs
    print('  PASS  job script: parity line first, no errexit, every point rc recorded, a failed point '
          'fails the job (rc 1) after the rest ran')


def check_missing_netcdf4_is_a_clear_hard_failure():
    saved = sys.modules.get('netCDF4', None)
    had = 'netCDF4' in sys.modules
    sys.modules['netCDF4'] = None          # import netCDF4 -> ImportError
    try:
        try:
            h.require_netcdf4()
        except h.SuiteError as e:
            assert 'netCDF4' in str(e) and 'install-eqdyna.sh -e vista' in str(e), e
        else:
            raise AssertionError('a missing netCDF4 did not raise')
    finally:
        if had:
            sys.modules['netCDF4'] = saved
        else:
            del sys.modules['netCDF4']
    h.require_netcdf4()                    # present here: must not raise
    print('  PASS  missing netCDF4 raises SuiteError naming `install-eqdyna.sh -e vista`')


def check_installer_vista_branch_creates_pinned_venv():
    import subprocess
    path = os.path.join(REPO_ROOT, 'install-eqdyna.sh')
    r = subprocess.run(['bash', '-n', path], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(path) as f:
        text = f.read()
    block = text[text.index('        vista)'):text.index('        grace)')]
    for want in ('python3 -m venv "$VENV"', 'numpy==2.0.2', 'netCDF4==1.7.2', 'xarray==2024.7.0',
                 '${EQDYNA_VENV:-$SCRATCH/eqdyna-venv}', 'if [ -n "$ENV" ]'):
        assert want in block, 'install-eqdyna.sh vista branch lacks %r' % want
    print('  PASS  install-eqdyna.sh -e vista: venv under $SCRATCH with numpy/netCDF4/xarray pinned')


def check_noisy_bucket_reported_below_noise_not_refused():
    """Row 142 (Vista smoke, 2026-11-24): one rank's tiny sub-bucket came out
    lo > hi and --analyze refused the WHOLE report. A sub-bucket is now
    reported below noise with the count of ranks affected and differenced on
    the rank-MEAN; loop_s keeps its hard raise (check_per_step_from_loop_s_...)."""
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_noise_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    _write_suite(suite_dir, [('strong_n0001', 4)], repeats=1)
    lo_path = os.path.join(suite_dir, 'strong_n0001', 'lo_r0', 'profile.rank0.json')
    hi_path = os.path.join(suite_dir, 'strong_n0001', 'hi_r0', 'profile.rank0.json')
    with open(lo_path) as f:
        lo = json.load(f)
    with open(hi_path) as f:
        hi = json.load(f)
    delta = 2.0 * hi['buckets_s']['exchange']       # lo exchange now exceeds hi: negative difference
    lo['buckets_s']['exchange'] += delta
    lo['loop_s'] += delta
    lo['total_s'] += delta
    with open(lo_path, 'w') as f:
        json.dump(lo, f)
    assert lo['loop_s'] < hi['loop_s'], 'fixture must keep loop_s positive: %r %r' % (lo['loop_s'], hi['loop_s'])
    try:
        numa.per_step_and_fixed(lo['buckets_s']['exchange'], hi['buckets_s']['exchange'],
                                lo['nsteps'], hi['nsteps'])
    except numa.NonPositiveMeasurement:
        pass
    else:
        raise AssertionError('fixture rank 0 exchange must difference non-positive')
    tgz = os.path.join(tmpd, 'noisy.tgz')
    with tarfile.open(tgz, 'w:gz') as tar:
        tar.add(suite_dir, arcname='.')
    report = h.analyze(tgz)
    pt = report['points'][0]
    bn = pt['buckets_below_noise']
    assert list(bn) == ['exchange'] and bn['exchange']['ranks_below_noise'] == 1 \
        and bn['exchange']['ranks_total'] == 4, bn
    assert bn['exchange']['rank_mean_per_step_s'] > 0 and \
        pt['bucket_breakdown_per_step_s']['exchange'] == bn['exchange']['rank_mean_per_step_s'], (
        bn, pt['bucket_breakdown_per_step_s'])
    assert 'exchange' not in pt['buckets_not_measured'], pt['buckets_not_measured']
    print('  PASS  noisy sub-bucket: reported below noise (1 of 4 ranks), rank-mean differenced, '
          'report not refused')


def check_job_script_sets_host_label_to_machine():
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_host_')
    path = os.path.join(tmpd, 'submit.sh')
    h.write_job_script(path, jobname='j', nnodes=1, ranks=4, partition='gg', walltime='00:05:00',
                       modules=[], commands=['cd %s && true' % tmpd], host_label='vista')
    with open(path) as f:
        text = f.read()
    assert 'export EQDYNA_HOST_LABEL=vista\n' in text, text
    print('  PASS  job script exports EQDYNA_HOST_LABEL=<machine> (profile host is not "workstation")')


def check_provenance_dirty_ignores_untracked():
    import subprocess
    repo = tempfile.mkdtemp(prefix='test_hpc_scaling_prov_')
    g = ['git', '-C', repo, '-c', 'user.email=t@t', '-c', 'user.name=t']
    subprocess.check_call(g + ['init', '-q'])
    with open(os.path.join(repo, 'a.txt'), 'w') as f:
        f.write('x\n')
    subprocess.check_call(g + ['add', 'a.txt'])
    subprocess.check_call(g + ['commit', '-q', '-m', 'm'])
    pts = [dict(point_id='strong_n0001', ranks=4, nnodes=1, dx=500.0, decomp=[4, 1, 1])]
    old = h.REPO_ROOT
    h.REPO_ROOT = repo
    try:
        with open(os.path.join(repo, 'build.log'), 'w') as f:
            f.write('artifact\n')
        clean = h.provenance('ubuntu', 100, 200, 1, pts, 'fortran', '/nonexistent/bin/eqdyna')
        with open(os.path.join(repo, 'a.txt'), 'w') as f:
            f.write('changed\n')
        dirty = h.provenance('ubuntu', 100, 200, 1, pts, 'fortran', '/nonexistent/bin/eqdyna')
    finally:
        h.REPO_ROOT = old
    assert clean['git_dirty'] is False and clean['git_untracked_count'] == 1, clean
    assert dirty['git_dirty'] is True, dirty
    print('  PASS  provenance: untracked build artifacts do not make git_dirty (counted instead); '
          'a tracked edit does')


def check_mixed_unmeasured_and_noisy_bucket_raises():
    """Audit MAJOR: python-jax-mpi declares `wait` unmeasured (0.0/0.0). One rank
    that measured it nonzero and differenced non-positive, beside ranks that read
    0.0/0.0, must raise -- never be folded into "not measured"."""
    tmpd = tempfile.mkdtemp(prefix='test_hpc_scaling_mixed_')
    suite_dir = os.path.join(tmpd, 'suite')
    os.makedirs(suite_dir)
    _write_suite(suite_dir, [('strong_n0001', 4)], repeats=1, backend='python-jax-mpi')
    lo_path = os.path.join(suite_dir, 'strong_n0001', 'lo_r0', 'profile.rank0.json')
    with open(lo_path) as f:
        lo = json.load(f)
    lo['buckets_s']['wait'] = 0.5
    lo['total_s'] += 0.5          # keep unaccounted_s == total_s - sum(buckets_s)
    with open(lo_path, 'w') as f:
        json.dump(lo, f)
    tgz = os.path.join(tmpd, 'mixed.tgz')
    _tar_suite(suite_dir, tgz)
    try:
        h.analyze(tgz)
    except h.SuiteError as e:
        assert 'declared unmeasured on some ranks' in str(e), str(e)
    else:
        raise AssertionError('a nonzero non-positive wait beside 0/0 ranks was reported as not measured')
    print('  PASS  mixed unmeasured/noisy bucket: raises, not folded into "not measured"')


def main():
    print('Regression guard: testsys/perf/hpc_scaling_suite.py (item 142)')
    checks = [check_bucket_partition,
              check_big_refuses_python_jax_mpi,
              check_per_rank_metrics_fixed_vs_per_step,
              check_fixed_buckets_would_trip_guard_if_differenced,
              check_build_points_strong_weak_size,
              check_step_pair_cli_defaults_and_validation,
              check_job_id_is_last_all_digit_line_and_nothing_else,
              check_submit_never_stores_the_account,
              check_solver_sha_of_the_bin,
              check_failed_sbatch_raises_suite_error_without_account,
              check_banner_only_sbatch_reply_raises,
              check_work_root_cases_built_where_the_job_cds,
              check_vista_dry_run_script,
              check_per_step_from_loop_s_and_nonpositive_raises,
              check_overhead_and_setup_split_required,
              check_profile_schema_setup_detail_shape,
              check_collect_analyze_overhead_repeats_and_parity_point,
              check_analyze_python_jax_mpi_suite,
              check_unmeasured_bucket_only_while_zero_and_setup_sum_two_sided,
              check_parity_block_verdicts_and_exit_codes,
              check_every_point_dx_meshes_and_stays_near_target,
              check_medium_smallest_points_build_and_mesh_lines,
              check_job_script_parity_first_no_errexit_and_rc_recorded,
              check_missing_netcdf4_is_a_clear_hard_failure,
              check_installer_vista_branch_creates_pinned_venv,
              check_noisy_bucket_reported_below_noise_not_refused,
              check_job_script_sets_host_label_to_machine,
              check_provenance_dirty_ignores_untracked,
              check_mixed_unmeasured_and_noisy_bucket_raises]
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
