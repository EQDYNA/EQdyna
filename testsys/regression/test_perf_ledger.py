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
