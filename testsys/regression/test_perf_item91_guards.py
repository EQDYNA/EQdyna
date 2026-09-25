#! /usr/bin/env python3
"""Behavioural regression guards for pathway item 91's six perf-tool
defects (testsys/perf/run_shard_scaling.py, run_numa_scaling.py,
run_perf.py, probe_scatter_bandwidth.py, run_setup_probe.py,
run_jaxmpi_ab.py -- sub-item e only). Full defect text: docs/BOARD_HISTORY.md
section "### Item 91".

Every check drives the real function with inputs that produce the bad case
and asserts it raises or labels -- never a grep over source text (rule 10a).
No MPI, no box-state dependence, no jax/numactl/git subprocess actually
runs: every subprocess boundary is monkeypatched so this guard is seconds
and machine-independent.

  (a) run_shard_scaling snapshot filename carries seconds, not just the
      date, and its rows carry the engine/policy keys the shared ledger
      reader (ledger.rows_from_scaling_snapshot) requires.
  (b) three per-step-by-difference sites (run_numa_scaling.per_step_and_fixed,
      probe_scatter_bandwidth.per_iter, run_perf.steady_state_per_step)
      raise on a non-positive result instead of recording it.
  (c) run_numa_scaling.time_one binds memory (--membind) as well as cpu
      (--physcpubind), via numactl, not bare taskset.
  (e) run_jaxmpi_ab.require_git_sha / run_setup_probe.require_git_sha raise
      on a git failure instead of returning ''.
  (f) run_perf.require_affinity_pinned raises when the pin did not take
      effect, instead of warning and proceeding.
  (g) run_shard_scaling.record_point / run_numa_scaling.record_baseline
      label the speedup baseline explicitly and never silently swap it
      when the intended first point was skipped.

Anti-vacuous-green discipline (papercuts): every verdict prints the content
property it asserted on.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, 'testsys', 'perf')
sys.path.insert(0, PERF)
sys.path.insert(0, ROOT)

import run_shard_scaling as rss    # noqa: E402
import run_numa_scaling as rns     # noqa: E402
import run_perf as rp              # noqa: E402
import probe_scatter_bandwidth as psb  # noqa: E402
import run_setup_probe as rsp      # noqa: E402
import run_jaxmpi_ab as rjab       # noqa: E402
import ledger                      # noqa: E402

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


def raises(fn, exc_types, what, must_contain=()):
    try:
        fn()
    except exc_types as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        check(not absent,
              '%s -> %s but message lacks %s: %s'
              % (what, type(e).__name__, absent, msg[:220])
              if absent else
              '%s -> %s; message: %s'
              % (what, type(e).__name__, msg[:220].replace('\n', ' ')))
        return
    check(False, '%s -> nothing raised' % what)


# ---------------------------------------------------------------- (a) ----

def check_91a_snapshot_naming_and_ledger():
    """Reloads the REAL run_shard_scaling module twice, under two faked wall
    times 5 minutes apart on the same calendar day, and compares its own
    module-level OUT global -- not a re-implementation of the naming
    expression. This is what actually catches the pre-fix `%Y-%m-%d`-only
    defect: with only the date, out1 == out2 here."""
    import importlib
    real_strftime = time.strftime

    def fake(t):
        def f(fmt, *a):
            return real_strftime(fmt, *a) if a else real_strftime(fmt, t)
        return f

    t1 = time.localtime(1_800_000_000)          # some day, some HH:MM:SS
    t2 = time.localtime(1_800_000_000 + 300)     # same day, 5 min later
    try:
        time.strftime = fake(t1)
        importlib.reload(rss)
        out1 = rss.OUT
        time.strftime = fake(t2)
        importlib.reload(rss)
        out2 = rss.OUT
    finally:
        time.strftime = real_strftime
        importlib.reload(rss)   # leave the module in its normal state
    check(out1 != out2,
          '91a: run_shard_scaling.OUT, reloaded at two wall times 5 min '
          'apart on the SAME calendar day, differ (%r vs %r)'
          % (out1, out2))

    # 91a's ledger half (2026-09-24): rows_from_scaling_snapshot(tool=...)
    # must file a run_shard_scaling point under its OWN tool name and the
    # EXISTING 'threads' parallelism value, never mislabelled as a
    # run_scaling row -- driven with a synthetic shard-shaped meta, no jax
    # or box-state dependence.
    meta = dict(sha='abc1234', host='h', date='2026-09-24 00:00',
               case='test.tpv104', busy_ceiling=0.2,
               rows=[dict(mode='element', n=2, ms_per_step=12.5,
                         cpus=[0, 1], n_lo=20, n_hi=60,
                         busy=(0.1, 0.1, 0.1, [], {0: 0.05, 1: 0.05}))])
    rows = ledger.rows_from_scaling_snapshot(
        meta, 'docs/perf_snapshots/fake_shard.json', dict(busy=0, total=4),
        tool='run_shard_scaling')
    check(len(rows) == 1 and rows[0]['tool'] == 'run_shard_scaling'
          and rows[0]['backend'] == 'python-jax'
          and rows[0]['parallelism'] == 'threads'
          and rows[0]['mode'] == 'element',
          '91a: a shard-shaped point files as tool=%r backend=%r '
          'parallelism=%r mode=%r (want run_shard_scaling/python-jax/'
          'threads/element)'
          % (rows[0].get('tool'), rows[0].get('backend'),
             rows[0].get('parallelism'), rows[0].get('mode')))
    check('run_shard_scaling' in open(os.path.join(
            PERF, 'run_shard_scaling.py')).read()
          and 'ledger.append_rows' in open(os.path.join(
              PERF, 'run_shard_scaling.py')).read(),
          '91a: run_shard_scaling.py itself calls ledger.append_rows (the '
          'wiring this check above proves the SHAPE of)')


def check_91_wall_clock_guard():
    """Row 91, second half (2026-09-24): run_mpi_scaling.per_step_jax_mpi's
    ms_per_step_wall difference (line ~409) must raise on a non-positive
    result, the SAME pattern PR #15 already applied to the rank-solve-time
    difference two lines below it (`ps_solve <= 0`) and to
    run_numa_scaling.per_step_and_fixed / probe_scatter_bandwidth.per_iter /
    run_perf.steady_state_per_step (check_91b above). Driven with
    jax_mpi_once faked to return a DECREASING wall clock as nsteps grows --
    every other dependency (profile_record.capture_run) fails past a
    nonexistent case_dir and is caught by per_step_jax_mpi's own
    warn-only try/except, so no real MPI/jax ever runs."""
    import run_mpi_scaling as rms

    def fake_decreasing(case_dir, n, cpus, ranks, sync, platform):
        return (15.0, []) if n == 20 else (5.0, [])

    orig = rms.jax_mpi_once
    try:
        rms.jax_mpi_once = fake_decreasing
        raises(lambda: rms.per_step_jax_mpi('/does/not/exist', [0], 1, 20,
                                            60, 'halo'),
               (RuntimeError,),
               '91 run_mpi_scaling.per_step_jax_mpi: wall clock DECREASED '
               'from n_lo=20 to n_hi=60',
               must_contain=('per-step', 'wall clock'))
    finally:
        rms.jax_mpi_once = orig

    def _rank(ms):
        return dict(ms_per_step=ms, mpi_ms_per_step=0.1, wait_ms_per_step=0.05,
                   effective_cores=1.0, Ei=10, Ep=5, halo_eqs=3,
                   carry_bytes_total=100, N_local=50, NEQ_local=30,
                   threads=1, cpus_allowed=1)

    def fake_increasing(case_dir, n, cpus, ranks, sync, platform):
        return ((5.0, [_rank(1.0)]) if n == 20 else (15.0, [_rank(2.0)]))

    orig = rms.jax_mpi_once
    try:
        rms.jax_mpi_once = fake_increasing
        r = rms.per_step_jax_mpi('/does/not/exist', [0], 1, 20, 60, 'halo')
        check(r is not None and r['ms_per_step_wall'] > 0,
              '91 run_mpi_scaling.per_step_jax_mpi: a genuine increasing-wall '
              'case still returns (ms_per_step_wall=%r)'
              % (r or {}).get('ms_per_step_wall'))
    finally:
        rms.jax_mpi_once = orig


# ---------------------------------------------------------------- (b) ----

def check_91b_nonpositive_raises():
    raises(lambda: rns.per_step_and_fixed(10.0, 9.0, 20, 60), (RuntimeError,),
           '91b run_numa_scaling.per_step_and_fixed: t_hi < t_lo',
           must_contain=('per-step',))
    ps, fixed = rns.per_step_and_fixed(1.0, 2.0, 10, 20)
    check(ps > 0, '91b run_numa_scaling.per_step_and_fixed: a genuine '
                  'positive case still returns (ps=%.4f)' % ps)

    orig_time_one = psb.time_one
    try:
        psb.time_one = lambda op, iters, cpus, node_map, dump_hlo_dir=None: (
            {30: 5.0, 100: 4.0}[iters])
        raises(lambda: psb.per_iter('control', [], {}, 30, 100),
               (RuntimeError,),
               '91b probe_scatter_bandwidth.per_iter: t_hi < t_lo',
               must_contain=('per-iteration',))
    finally:
        psb.time_one = orig_time_one

    orig_time_python = rp.time_python
    try:
        rp.time_python = lambda engine, n: (
            (0.0, 5.0) if n == 10 else (0.0, 4.0))
        raises(lambda: rp.steady_state_per_step('numpy', 10, 20),
               (RuntimeError,),
               '91b run_perf.steady_state_per_step: solve(n_hi) < solve(n_lo)',
               must_contain=('per-step',))
    finally:
        rp.time_python = orig_time_python


# ---------------------------------------------------------------- (c) ----

def check_91c_membind():
    got = {}

    def fake_run(cmd, env=None, capture_output=None, text=None):
        got['cmd'] = cmd
        class R:
            returncode = 0
            stdout = 'WALL 1.23\n'
            stderr = ''
        return R()

    orig_run = rns.subprocess.run
    try:
        rns.subprocess.run = fake_run
        node_map = {0: [0, 1, 2, 3, 4, 5, 6, 7],
                    1: [8, 9, 10, 11, 12, 13, 14, 15]}
        w = rns.time_one('/case', 20, [8, 9], node_map)
    finally:
        rns.subprocess.run = orig_run
    cmd = got['cmd']
    check(w == 1.23, '91c: time_one still returns the wall seconds (%.2f)' % w)
    check(cmd[0] == 'numactl',
          '91c: time_one launches via numactl, not bare taskset (cmd[0]=%r)'
          % cmd[0])
    check(any(c == '--physcpubind=8,9' for c in cmd),
          '91c: --physcpubind carries the exact cpu list (cmd=%s)' % cmd)
    check(any(c == '--membind=1' for c in cmd),
          '91c: --membind carries the NUMA node(s) those cpus belong to '
          '(node 1 for cpus 8,9) -- the fix run_scaling.py:18-26 documents '
          '(cmd=%s)' % cmd)
    check(rns._nodes_of(node_map, [8, 9]) == [1],
          '91c: _nodes_of([8,9]) == [1] (pure helper, no subprocess)')
    check(rns._nodes_of(node_map, [0, 8]) == [0, 1],
          '91c: _nodes_of([0,8]) == [0,1] spans both nodes when cpus do')


# ---------------------------------------------------------------- (e) ----

def check_91e_git_sha():
    for mod, name in ((rjab, 'run_jaxmpi_ab'), (rsp, 'run_setup_probe')):
        raises(lambda mod=mod: mod.require_git_sha(
                   1, '', 'fatal: bad revision', 'git rev-parse'),
               (SystemExit,),
               '91e %s.require_git_sha: git returncode=1' % name,
               must_contain=('git rev-parse',))
        raises(lambda mod=mod: mod.require_git_sha(0, '', '', 'git rev-parse'),
               (SystemExit,),
               '91e %s.require_git_sha: returncode=0 but empty stdout '
               '(the exact defect: recording sha="")' % name)
        sha = mod.require_git_sha(0, ' abcd123\n', '', 'git rev-parse')
        check(sha == 'abcd123',
              '91e %s.require_git_sha: success returns the stripped sha '
              '(%r)' % (name, sha))


# ---------------------------------------------------------------- (f) ----

def check_91f_affinity_gate():
    try:
        rp.require_affinity_pinned('[0]', '0')
        ok = True
    except SystemExit:
        ok = False
    check(ok, '91f: an affinity match ([0] for core 0) does not raise')
    raises(lambda: rp.require_affinity_pinned('[0,1]', '0'), (SystemExit,),
           '91f: an affinity mismatch ([0,1] for core 0) refuses instead of '
           'warning-and-proceeding')


# ---------------------------------------------------------------- (g) ----

def check_91g_run_scaling_loop_labels_baseline():
    """91g residual (item 91 PR #15 NOTES): run_scaling.py's own python loop,
    driven FOR REAL through main() with every measurement and box probe
    stubbed (no case, no jax, no numactl, nothing written to the repo): with
    th=1 SKIPPED, the th=2 row must carry baseline_n=2 and a NOTE must be
    printed -- before, base.setdefault silently made th=2 read 1.00x."""
    import io, contextlib, json as _json, tempfile as _tf, types
    import run_scaling as rscal
    tmp = _tf.mkdtemp(prefix='item91g_rs.')
    saved = {k: getattr(rscal, k) for k in
             ('free_node_map', 'select_cpus', 'build_py_case', 'per_step_py',
              'nodes_of', 'sh', 'OUT', 'ROOT')}
    saved_numa = (rscal.numa.numa_topology, rscal.numa.require_idle)
    saved_led = (ledger.append_rows, ledger.box_tenancy)
    saved_argv = sys.argv
    try:
        rscal.numa.numa_topology = lambda: {0: [0, 1, 2, 3]}
        # A real require_idle() 5-tuple (load1, load5, load15, other_users,
        # busy_dict), not the old placeholder `{}` -- since row 76 (2026-09-24)
        # rows_from_scaling_snapshot reads this row's 'busy' field for real
        # (ledger.contention_from_check), an empty/wrong-shaped stub here
        # would make the REAL converter raise "named no cpus", not the
        # vacuous no-op it used to be.
        rscal.numa.require_idle = lambda cpus, ceil, ov: (
            0.1, 0.1, 0.1, [], {c: 0.05 for c in cpus})
        rscal.free_node_map = lambda nodes, ceil, ov: {0: [0, 1, 2, 3]}
        rscal.select_cpus = lambda free, k, policy: None if k == 1 else list(range(k))
        rscal.build_py_case = lambda case: tmp
        rscal.per_step_py = lambda d, cpus, nodes, backend, lo, hi: (0.2 / len(cpus), 1.0, 2.0, 3.0)
        rscal.nodes_of = lambda nodes, cpus: [0]
        rscal.sh = lambda cmd, **kw: types.SimpleNamespace(stdout='abc1234\n', returncode=0)
        rscal.OUT = os.path.join(tmp, 'scaling_last.json')
        rscal.ROOT = tmp
        ledger.append_rows = lambda rows: len(rows)
        ledger.box_tenancy = lambda ceil: dict(busy=0, total=4)
        sys.argv = ['run_scaling.py', '--skip-fortran', '--backends', 'numpy',
                    '--py-threads', '1,2,4', '--policies', 'compact']
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rscal.main()
        rows = _json.load(open(rscal.OUT))['rows']
    finally:
        for k, v in saved.items():
            setattr(rscal, k, v)
        rscal.numa.numa_topology, rscal.numa.require_idle = saved_numa
        ledger.append_rows, ledger.box_tenancy = saved_led
        sys.argv = saved_argv
    got = [(r['n'], r.get('baseline_n'), round(r['speedup'], 2)) for r in rows]
    check(got == [(2, 2, 1.0), (4, 2, 2.0)] and 'NOTE: baseline for mode' in buf.getvalue(),
          '91g run_scaling main(): th=1 skipped -> rows (n, baseline_n, speedup) '
          '%r want [(2, 2, 1.0), (4, 2, 2.0)] and a printed NOTE' % (got,))


def check_91g_baseline_labelling():
    base, base_n = {}, {}
    _, bn1, note1 = rss.record_point(base, base_n, 'element', 2, 100.0, 1)
    check(bn1 == 2 and note1 is not None and '1' in note1 and '2' in note1,
          '91g run_shard_scaling.record_point: n=1 was skipped, n=2 sets '
          'the baseline and gets a non-None NOTE naming both n=%r vs the '
          'devices[0]=1 it silently used to be read against (note=%r)'
          % (bn1, note1))
    speedup2, bn2, note2 = rss.record_point(base, base_n, 'element', 2, 100.0, 1)
    check(speedup2 == 1.0 and note2 is None,
          '91g: the already-set baseline point itself gets no repeat NOTE '
          '(speedup=%.2f, note=%r)' % (speedup2, note2))

    base2, base_n2 = {}, {}
    _, bn3, note3 = rss.record_point(base2, base_n2, 'element', 1, 100.0, 1)
    check(bn3 == 1 and note3 is None,
          '91g: when n=1 (devices[0]) DOES set the baseline, no NOTE fires '
          '(bn=%r note=%r)' % (bn3, note3))

    state = {}
    speedup, blabel, note = rns.record_baseline(state, 'spread-8', 50.0,
                                                 'within-node-8')
    check(blabel == 'spread-8' and note is not None
          and 'within-node-8' in note and 'spread-8' in note,
          '91g run_numa_scaling.record_baseline: first config skipped, '
          'second config (%r) becomes baseline with a non-None NOTE (%r)'
          % (blabel, (note or '')[:80]))


def check_91a_backfill_detects_shard_schema():
    """Audit finding (PR #34): _rows_from_snapshot_file sent a shard-shaped
    snapshot to rows_from_scaling_snapshot with the default tool, filing its
    points as run_scaling. It must detect the shard schema."""
    import json, tempfile
    snapdir = os.path.join(ROOT, 'docs', 'perf_snapshots')
    meta = dict(case='test.tpv8', sha='0000000', host='h', date='2026-09-25 00:00',
                n_lo=5, n_hi=10, busy_ceiling=0.5, overridden=False,
                fortran_quoted={}, element_baseline_n=1, skipped=[],
                rows=[dict(n=1, mode='element', ms_per_step=10.0, n_lo=5, n_hi=10,
                           cpus=[0], busy={'0': 0.0})])
    fd, path = tempfile.mkstemp(suffix='.json', prefix='shard_test_', dir=snapdir)
    try:
        with os.fdopen(fd, 'w') as fh:
            json.dump(meta, fh)
        rows, _ = ledger._rows_from_snapshot_file(path)
    finally:
        os.remove(path)
    tools = {r['tool'] for r in rows}
    if tools != {'run_shard_scaling'}:
        FAILURES.append('91a-backfill: shard snapshot backfilled as %r' % tools)
        print('FAIL -- 91a-backfill: shard snapshot backfilled as %r, want run_shard_scaling' % tools)
        return
    bad = dict(meta['rows'][0]); bad.pop('mode')
    try:
        ledger.rows_from_scaling_snapshot(dict(meta, rows=[dict(bad, engine='python-jax')]),
                                          'x', None, tool='run_shard_scaling')
        FAILURES.append('91a-backfill: a shard row carrying engine was accepted')
        print('FAIL -- 91a-backfill: a shard row carrying engine was accepted')
        return
    except ValueError:
        pass
    print('PASS -- 91a-backfill: a shard-schema snapshot backfills as tool=run_shard_scaling; '
          'a shard row carrying engine is refused')


CHECKS = [
    ('91a', check_91a_snapshot_naming_and_ledger),
    ('91a-backfill', check_91a_backfill_detects_shard_schema),
    ('91-wall-clock', check_91_wall_clock_guard),
    ('91b', check_91b_nonpositive_raises),
    ('91c', check_91c_membind),
    ('91e', check_91e_git_sha),
    ('91f', check_91f_affinity_gate),
    ('91g', check_91g_baseline_labelling),
    ('91g-run_scaling', check_91g_run_scaling_loop_labels_baseline),
]


def main():
    for tag, fn in CHECKS:
        print('\n-- item %s --' % tag)
        fn()
    print()
    if FAILURES:
        print('FAIL test_perf_item91_guards: %d check(s) failed' % len(FAILURES))
        return 1
    print('SUCCESS test_perf_item91_guards: all item-91 perf-tool defects '
          'guarded behaviourally (%d checks)' % len(CHECKS))
    return 0


if __name__ == '__main__':
    sys.exit(main())


