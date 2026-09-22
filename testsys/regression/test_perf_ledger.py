#! /usr/bin/env python3
"""Guard for the append-only perf ledger (docs/perf_ledger.jsonl, written by
testsys/perf/ledger.py).

Incident guarded: PROJECT_RULES rule 19 (2026-09-21) -- a committed shared
mutable perf artifact (scaling_last.json) was nearly overwritten, destroying
the data a board row cited, with two agents holding it dirty at once. The
ledger is the positive counterpart; this guard makes sure it stays one.

Anti-vacuous-green discipline (papercuts: four gates went green while testing
NOTHING): every verdict below prints the CONTENT property it asserted on --
row counts, parsed field values, byte counts -- never just "file exists".

Checks:
  1. round-trip     -- an appended row parses back with the same values.
  2. append-only    -- N appends -> N lines; after each append the previous
                       bytes are a byte-identical prefix.
  3. rejection      -- 10 malformed rows each raise ValueError (rule 2).
  4. concurrency    -- 6 processes x 40 rows appended simultaneously ->
                       exactly 240 parseable, distinct lines (flock +
                       single O_APPEND write actually holds).
  5. committed ledger -- every line of docs/perf_ledger.jsonl validates.
  6. history immutability -- HEAD's ledger content is a byte-prefix of the
                       working copy (append-only against git history).
  7. platform identity -- (incident: 812.78 ms/step, a CPU measurement
                       recorded under a cuda label, 2026-09-21 gpu1 snapshot)
                       a row APPENDED without platform fields is refused; a
                       'gpu' claim without positive per-device memory
                       evidence is refused; legacy platform-less rows still
                       validate when READ; unknown is never cpu.
  8. evidence over label -- on the REAL committed GPU snapshots: gpu1 (label
                       cuda, 0 MiB delta on all four A100s) converts to
                       platform 'unknown', never cpu/gpu; gpu_final converts
                       to 'gpu' with devices counted from the deltas; a
                       gpu_final with its evidence stripped converts to
                       'unknown', not to the label.
  9. cell-wall-clock -- a cell row must carry wall_s and an explicitly-null
                       ms_per_step (the metric-mixing that once produced a
                       spurious 33% JAX regression must be structurally
                       impossible); the e2e converter emits no row for a
                       FAILED cell.
 10. sweep capture never gates -- capture_e2e_cells writes snapshot + rows on
                       a clean path; with tenancy unreadable the _or_warn
                       wrapper WARNS, appends nothing, and raises nothing.
"""
import copy
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, 'testsys', 'perf')
sys.path.insert(0, PERF)
import ledger  # noqa: E402

FAILURES = []


def check(ok, what):
    print('%s -- %s' % ('PASS' if ok else 'FAIL', what))
    if not ok:
        FAILURES.append(what)


def good_row(i=0, proc=0):
    return dict(ts_utc='2026-09-21T12:00:00Z', snapshot_date_local='2026-09-21 12:00',
                sha='0123abc', host='cotopaxi', tool='run_scaling',
                case='test.tpv104', backend='fortran', ranks=4,
                ms_per_step=100.0 + i, rank_ms_min=90.0, rank_ms_max=110.0,
                rank_ms_mean=100.0, effective_cores=3.8, threads_per_rank=1,
                busy_ceiling=0.5, tenancy_busy=20, tenancy_total=64,
                metric='per-step-by-difference', n_lo=20, n_hi=60,
                snapshot='docs/perf_snapshots/scaling_2026-09-21_x.json',
                platform='cpu', devices=None,
                platform_evidence='synthetic test row: pinned to one cpu',
                proc=proc, i=i)


def main():
    tmp = tempfile.mkdtemp(prefix='perf_ledger_test.')
    led = os.path.join(tmp, 'ledger.jsonl')

    # 1. round-trip: parsed values, not just presence.
    ledger.append(good_row(i=7), path=led)
    lines = open(led, 'rb').read().split(b'\n')[:-1]
    row = json.loads(lines[0])
    check(len(lines) == 1, 'round-trip: 1 append -> 1 line (got %d)' % len(lines))
    check(row['ms_per_step'] == 107.0 and row['ranks'] == 4
          and row['snapshot'].startswith('docs/perf_snapshots/'),
          'round-trip: parsed ms_per_step=%r ranks=%r snapshot=%r match input'
          % (row['ms_per_step'], row['ranks'], row['snapshot']))

    # 2. append-only: prefix bytes never change.
    n_appends = 5
    prefix_ok = True
    prev = open(led, 'rb').read()
    for i in range(1, n_appends):
        ledger.append(good_row(i=i), path=led)
        now = open(led, 'rb').read()
        if not now.startswith(prev) or len(now) <= len(prev):
            prefix_ok = False
        prev = now
    nlines = len(prev.split(b'\n')) - 1
    check(nlines == n_appends,
          'append-only: %d appends -> %d lines' % (n_appends, nlines))
    check(prefix_ok, 'append-only: after each of %d appends the earlier '
                     'content was a byte-identical prefix' % n_appends)

    # 3. rejection battery: every malformed row must raise (rule 2).
    def broken(**kw):
        r = copy.deepcopy(good_row())
        for k, v in kw.items():
            if v is Ellipsis:
                del r[k]
            else:
                r[k] = v
        return r
    cases = [
        ('missing sha', broken(sha=Ellipsis)),
        ('missing snapshot', broken(snapshot=Ellipsis)),
        ('absolute snapshot path', broken(snapshot='/etc/passwd')),
        ('snapshot outside docs/perf_snapshots/', broken(snapshot='docs/x.json')),
        ('non-positive ms_per_step', broken(ms_per_step=-41.67)),
        ('malformed sha', broken(sha='not-a-sha')),
        ('n_lo >= n_hi', broken(n_lo=60, n_hi=60)),
        ('unknown backend', broken(backend='python-torch')),
        ('null tenancy on a live (non-backfilled) row',
         broken(tenancy_busy=None, tenancy_total=None)),
        ('unknown metric', broken(metric='total-wall-clock')),
    ]
    for what, r in cases:
        try:
            ledger.append(r, path=led)
            check(False, 'rejection: %s raised ValueError' % what)
        except ValueError:
            check(True, 'rejection: %s raised ValueError' % what)
    lines_after = sum(1 for _ in open(led))
    check(lines_after == n_appends,
          'rejection: no rejected row reached the ledger (still %d lines, '
          'expected %d)' % (lines_after, n_appends))

    # 3b. a backfilled row MAY carry null tenancy (old snapshots lack it).
    r = broken(tenancy_busy=None, tenancy_total=None)
    r['backfilled_from'] = 'docs/perf_snapshots/scaling_2026-09-21_x.json'
    try:
        ledger.validate(r)
        check(True, 'backfilled row with null tenancy validates')
    except ValueError as e:
        check(False, 'backfilled row with null tenancy validates (%s)' % e)

    # 4. concurrency: 6 processes x 40 rows, all at once.
    nproc, per = 6, 40
    led2 = os.path.join(tmp, 'concurrent.jsonl')
    prog = ('import sys; sys.path.insert(0, %r); import ledger\n'
            'sys.path.insert(0, %r); from test_perf_ledger import good_row\n'
            'p = int(sys.argv[1])\n'
            'for i in range(%d): ledger.append(good_row(i=i, proc=p), path=%r)\n'
            % (PERF, HERE, per, led2))
    procs = [subprocess.Popen([sys.executable, '-c', prog, str(p)])
             for p in range(nproc)]
    rcs = [p.wait() for p in procs]
    check(all(rc == 0 for rc in rcs),
          'concurrency: all %d appender processes exited 0 (rcs %s)'
          % (nproc, rcs))
    raw = open(led2, 'rb').read()
    lines = raw.split(b'\n')[:-1]
    parsed, ids = [], set()
    torn = 0
    for ln in lines:
        try:
            d = json.loads(ln)
            parsed.append(d)
            ids.add((d['proc'], d['i']))
        except ValueError:
            torn += 1
    check(len(lines) == nproc * per,
          'concurrency: %d lines for %d expected appends' % (len(lines), nproc * per))
    check(torn == 0, 'concurrency: 0 torn/unparseable lines (got %d)' % torn)
    check(len(ids) == nproc * per,
          'concurrency: %d distinct (proc,i) payloads of %d expected -- no '
          'append silently clobbered another' % (len(ids), nproc * per))

    # 4b. snapshot->ledger converters, on synthetic recorded-shape input
    # (the committed mpi snapshot is all-SKIPPED, so without this the jax-mpi
    # converter path would be green-by-absence).
    mpi_meta = dict(case='test.tpv104', sha='abc1234', host='cotopaxi',
                    date='2026-09-21 16:37', n_lo=40, n_hi=160, max_busy=0.6,
                    rows=[dict(ranks=2, skipped=True, cpus=[2, 3]),
                          dict(ranks=4, cpus=[8, 9, 16, 17], n_lo=40, n_hi=160,
                               jax_halo=dict(ms_per_step=605.0,
                                             rank_ms=[600.0, 605.0, 590.0, 601.0],
                                             eff=[3.9, 3.8, 4.0, 3.7],
                                             threads=[1, 1, 1, 1]),
                               fortran=dict(ms_per_step=101.5))])
    conv = ledger.rows_from_mpi_scaling_snapshot(
        mpi_meta, 'docs/perf_snapshots/mpi_scaling_x.json',
        dict(busy=30, total=64))
    check(len(conv) == 2,
          'converter(mpi): 1 skipped + 1 measured row -> 2 ledger rows '
          '(got %d)' % len(conv))
    jrow = next(r for r in conv if r['backend'] == 'python-jax-mpi')
    frow = next(r for r in conv if r['backend'] == 'fortran')
    check(jrow['rank_ms_min'] == 590.0 and jrow['rank_ms_max'] == 605.0
          and abs(jrow['rank_ms_mean'] - 599.0) < 1e-9
          and jrow['sync'] == 'halo' and jrow['tenancy_busy'] == 30,
          'converter(mpi): jax row min/max/mean=%r/%r/%r sync=%r tenancy=%r/%r'
          % (jrow['rank_ms_min'], jrow['rank_ms_max'], jrow['rank_ms_mean'],
             jrow['sync'], jrow['tenancy_busy'], jrow['tenancy_total']))
    check(frow['ms_per_step'] == 101.5 and frow['rank_ms_min'] is None,
          'converter(mpi): fortran row ms_per_step=%r, per-rank fields null'
          % frow['ms_per_step'])
    sc_meta = dict(case='test.tpv104', sha='abc1234', host='cotopaxi',
                   date='2026-09-21 11:35', busy_ceiling=0.45,
                   rows=[dict(engine='python-numpy', n=8, policy='spread',
                              ms_per_step=250.25, n_lo=20, n_hi=60,
                              cpus=[0, 8, 16, 24, 32, 40, 48, 56])])
    conv2 = ledger.rows_from_scaling_snapshot(
        sc_meta, 'docs/perf_snapshots/scaling_x.json', dict(busy=10, total=64))
    check(len(conv2) == 1 and conv2[0]['backend'] == 'python-numpy'
          and conv2[0]['ms_per_step'] == 250.25 and conv2[0]['ranks'] == 8
          and conv2[0]['policy'] == 'spread',
          'converter(scaling): row backend=%r ms_per_step=%r ranks=%r policy=%r'
          % (conv2[0]['backend'], conv2[0]['ms_per_step'], conv2[0]['ranks'],
             conv2[0]['policy']))

    # 7. platform identity. Incident: 812.78 ms/step -- a CPU measurement
    # recorded under a cuda label (eqdyna3d --device auto overriding
    # JAX_PLATFORMS=cuda; 0 MiB delta on all four A100s).
    legacy = good_row()
    for k in ('platform', 'platform_evidence', 'devices'):
        del legacy[k]
    try:
        ledger.validate(legacy)
        check(True, 'platform: legacy row (no platform key) validates when READ')
    except ValueError as e:
        check(False, 'platform: legacy row (no platform key) validates when '
                     'READ (%s)' % e)
    try:
        ledger.append(legacy, path=led)
        check(False, 'platform: APPENDING a platform-less row raises ValueError')
    except ValueError:
        check(True, 'platform: APPENDING a platform-less row raises ValueError')
    naked_gpu = broken(platform='gpu', devices=1,
                       platform_evidence='claims a GPU with no evidence')
    try:
        ledger.append(naked_gpu, path=led)
        check(False, "platform: 'gpu' claim without gpu_delta_mib evidence "
                     'raises ValueError')
    except ValueError:
        check(True, "platform: 'gpu' claim without gpu_delta_mib evidence "
                    'raises ValueError')
    evidenced_gpu = broken(platform='gpu', devices=2,
                           platform_evidence='nvidia-smi delta on devices 0,1')
    evidenced_gpu['gpu_delta_mib'] = {'0': 30842, '1': 30842}
    try:
        ledger.validate(evidenced_gpu, appending=True)
        check(True, "platform: 'gpu' claim WITH per-device delta evidence "
                    'validates')
    except ValueError as e:
        check(False, "platform: 'gpu' claim WITH per-device delta evidence "
                     'validates (%s)' % e)
    try:
        ledger.append(broken(devices=3), path=led)
        check(False, 'platform: a cpu row carrying devices=3 raises ValueError')
    except ValueError:
        check(True, 'platform: a cpu row carrying devices=3 raises ValueError')
    unk = broken(platform='unknown',
                 platform_evidence='snapshot records no device evidence')
    try:
        ledger.validate(unk, appending=True)
        check(unk['platform'] not in ('cpu', 'gpu'),
              "platform: 'unknown' validates and is distinguishable from cpu "
              'and gpu (%r)' % unk['platform'])
    except ValueError as e:
        check(False, "platform: 'unknown' validates (%s)" % e)
    n_after_7 = sum(1 for _ in open(led))
    check(n_after_7 == n_appends,
          'platform: no refused row reached the ledger (still %d lines)'
          % n_after_7)

    # 8. evidence over label, on the REAL committed GPU snapshots. gpu1 is
    # the incident itself: label cuda, 0 MiB delta everywhere, 812.78 ms/step.
    snapdir = os.path.join(ROOT, 'docs', 'perf_snapshots')
    gpu1_rel = 'docs/perf_snapshots/mpi_scaling_2026-09-21_tpv104_gpu1.json'
    gpu1 = json.load(open(os.path.join(snapdir, os.path.basename(gpu1_rel))))
    rows1 = ledger.rows_from_mpi_scaling_snapshot(gpu1, gpu1_rel, None,
                                                  backfilled_from=gpu1_rel)
    check(len(rows1) == 1 and rows1[0]['platform'] == 'unknown'
          and abs(rows1[0]['ms_per_step'] - 812.784) < 0.01,
          'evidence-over-label: gpu1 (the 812.78 incident, label=cuda, all '
          'deltas 0) converts to platform=%r, NOT cpu and NOT gpu'
          % rows1[0]['platform'])
    final_rel = ('docs/perf_snapshots/'
                 'mpi_scaling_2026-09-21_tpv104_gpu_final.json')
    final = json.load(open(os.path.join(snapdir, os.path.basename(final_rel))))
    rowsf = ledger.rows_from_mpi_scaling_snapshot(final, final_rel, None,
                                                  backfilled_from=final_rel)
    devs = [r['devices'] for r in rowsf]
    check(len(rowsf) == 3 and all(r['platform'] == 'gpu' for r in rowsf)
          and devs == [1, 2, 4]
          and all(any(v > 0 for v in r['gpu_delta_mib'].values())
                  for r in rowsf)
          and all(r['device_peak_gb'] for r in rowsf),
          'evidence-over-label: gpu_final converts to 3 platform=gpu rows, '
          'devices=%s counted from the nvidia-smi deltas, allocator '
          'device_peak_gb carried' % devs)
    check(all(r['ranks'] == r_meta['ranks'] for r, r_meta
              in zip(rowsf, final['rows'])),
          'evidence-over-label: ranks still counts MPI processes (%s); '
          'devices is the separate GPU count'
          % [r['ranks'] for r in rowsf])
    # the NEGATIVE: strip/zero the device evidence of a genuine GPU snapshot;
    # the converter must NOT record the label.
    stripped = copy.deepcopy(final)
    for r in stripped['rows']:
        for k in list(r):
            if k.startswith('jax'):
                r[k]['gpu_mem']['delta'] = {d: 0 for d in r[k]['gpu_mem']['delta']}
                r[k].pop('device_peak_gb', None)
    rows_s = ledger.rows_from_mpi_scaling_snapshot(stripped, final_rel, None,
                                                   backfilled_from=final_rel)
    check(all(r['platform'] == 'unknown' for r in rows_s)
          and all(r['platform'] != 'cpu' for r in rows_s)
          and all(r['platform'] != 'gpu' for r in rows_s),
          'evidence-over-label (NEGATIVE): gpu_final with its device evidence '
          'zeroed converts to platforms %s -- the cuda label alone is not '
          'recorded as a measurement' % sorted({r['platform'] for r in rows_s}))
    # and a label claiming cpu over live GPU evidence is equally distrusted.
    lied = copy.deepcopy(final)
    for r in lied['rows']:
        for k in list(r):
            if k.startswith('jax'):
                r[k]['platform'] = 'cpu'
    rows_l = ledger.rows_from_mpi_scaling_snapshot(lied, final_rel, None,
                                                   backfilled_from=final_rel)
    check(all(r['platform'] == 'unknown' for r in rows_l),
          'evidence-over-label (NEGATIVE): a "cpu" label over positive GPU '
          'deltas converts to %s, not cpu'
          % sorted({r['platform'] for r in rows_l}))

    # 9. cell-wall-clock is a DIFFERENT metric and cannot be read as per-step.
    def cell_row(**kw):
        r = good_row()
        r.update(metric='cell-wall-clock', ms_per_step=None, n_lo=None,
                 n_hi=None, wall_s=123.4, tool='run_e2e',
                 rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
                 effective_cores=None, threads_per_rank=None)
        r.update(kw)
        return r
    try:
        ledger.validate(cell_row(), appending=True)
        check(True, 'cell-wall-clock: a well-formed cell row validates '
                    '(wall_s=123.4, ms_per_step null)')
    except ValueError as e:
        check(False, 'cell-wall-clock: a well-formed cell row validates (%s)'
                     % e)
    for what, r in [
            ('cell row carrying a numeric ms_per_step (the averaging trap)',
             cell_row(ms_per_step=123.4)),
            ('cell row without wall_s', cell_row(wall_s=None)),
            ('cell row with a fake n_lo/n_hi pair',
             cell_row(n_lo=20, n_hi=60))]:
        try:
            ledger.validate(r, appending=True)
            check(False, 'cell-wall-clock: %s raises ValueError' % what)
        except ValueError:
            check(True, 'cell-wall-clock: %s raises ValueError' % what)
    e2e_meta = dict(tool='run_e2e', sha='abc1234', host='cotopaxi',
                    date='2026-09-22 07:00', label='default', device='cpu',
                    jobs_budget=8, tenancy_ceiling=0.5,
                    cells=[dict(case='test.tpv8', backend='fortran', ok=True,
                                seconds=41.2, ranks=4, platform='cpu',
                                platform_evidence='fortran: no GPU path'),
                           dict(case='test.tpv8', backend='python-jax',
                                ok=False, seconds=3.1, ranks=1,
                                platform='cpu',
                                platform_evidence='JAX_PLATFORMS=cpu pinned'),
                           dict(case='test.tpv8', backend='python-numpy',
                                ok=True, seconds=97.6, ranks=1,
                                platform='cpu',
                                platform_evidence='numpy: no GPU path')])
    erows = ledger.rows_from_e2e_results(
        e2e_meta, 'docs/perf_snapshots/e2e_cells_x.json',
        dict(busy=12, total=64))
    check(len(erows) == 2
          and sorted(r['backend'] for r in erows) == ['fortran',
                                                      'python-numpy']
          and all(r['metric'] == 'cell-wall-clock' for r in erows)
          and all(r['ms_per_step'] is None for r in erows)
          and {r['backend']: r['wall_s'] for r in erows}
              == {'fortran': 41.2, 'python-numpy': 97.6},
          'e2e converter: 2 PASSING cells -> 2 rows (wall_s 41.2/97.6, '
          'ms_per_step null); the FAILED python-jax cell produced no row')

    # 10. sweep capture: real write path, then the degrade-to-warning path.
    import contextlib
    import io
    root2 = os.path.join(tmp, 'e2e_root')
    real_tenancy = ledger.box_tenancy
    try:
        ledger.box_tenancy = lambda ceiling: dict(busy=9, total=64)
        n = ledger.capture_e2e_cells_or_warn(e2e_meta, root=root2)
        led3 = os.path.join(root2, 'docs', 'perf_ledger.jsonl')
        got = [json.loads(ln) for ln in open(led3)]
        snaps = os.listdir(os.path.join(root2, 'docs', 'perf_snapshots'))
        check(n == 2 and len(got) == 2
              and all(g['tenancy_busy'] == 9 for g in got)
              and len(snaps) == 1
              and all(g['snapshot'] == 'docs/perf_snapshots/' + snaps[0]
                      for g in got),
              'capture: clean path appended %r row(s), wrote 1 snapshot (%s) '
              'that every row points at' % (n, snaps))

        def _no_tenancy(ceiling):
            raise SystemExit('FAIL: numactl --hardware gave no topology')
        ledger.box_tenancy = _no_tenancy
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            n2 = ledger.capture_e2e_cells_or_warn(e2e_meta, root=root2)
        still = sum(1 for _ in open(led3))
        check(n2 is None and still == 2 and 'WARNING' in buf.getvalue()
              and 'verdict' in buf.getvalue(),
              'capture: with tenancy unreadable the wrapper returned %r, '
              'appended nothing (still %d rows), raised nothing, and '
              'PRINTED a warning' % (n2, still))
    finally:
        ledger.box_tenancy = real_tenancy

    # 5. the COMMITTED ledger validates row by row.
    real = os.path.join(ROOT, 'docs', 'perf_ledger.jsonl')
    if not os.path.exists(real):
        check(False, 'committed ledger docs/perf_ledger.jsonl exists')
    else:
        bad, nrows, nbackfilled = [], 0, 0
        for lineno, ln in enumerate(open(real), 1):
            try:
                d = json.loads(ln)
                ledger.validate(d)
                nrows += 1
                nbackfilled += 1 if 'backfilled_from' in d else 0
            except ValueError as e:
                bad.append('line %d: %s' % (lineno, e))
        check(not bad and nrows >= 1,
              'committed ledger: %d row(s) all validate (%d backfilled); '
              'problems: %s' % (nrows, nbackfilled, bad or 'none'))

        # 6. append-only against git history: HEAD content is a byte-prefix.
        rel = 'docs/perf_ledger.jsonl'
        exists = subprocess.run(['git', 'cat-file', '-e', 'HEAD:%s' % rel],
                                cwd=ROOT, capture_output=True)
        staged = subprocess.run(['git', 'ls-files', '--cached', '--', rel],
                                cwd=ROOT, capture_output=True, text=True)
        if exists.returncode == 0:
            head = subprocess.run(['git', 'show', 'HEAD:%s' % rel], cwd=ROOT,
                                  capture_output=True).stdout
            work = open(real, 'rb').read()
            check(work.startswith(head),
                  'history: HEAD ledger (%d bytes) is a byte-prefix of the '
                  'working copy (%d bytes) -- append-only holds against git '
                  'history' % (len(head), len(work)))
        elif staged.stdout.strip():
            check(True, 'history: ledger staged but not yet in HEAD '
                        '(first-commit state); prefix check trivially true')
        else:
            check(False, 'history: docs/perf_ledger.jsonl is neither in HEAD '
                         'nor staged -- an untracked ledger is not evidence '
                         '(rule 19)')

    print()
    if FAILURES:
        print('FAIL test_perf_ledger: %d check(s) failed' % len(FAILURES))
        return 1
    print('SUCCESS test_perf_ledger: ledger append-only, validated, '
          'concurrency-safe')
    return 0


if __name__ == '__main__':
    sys.exit(main())
