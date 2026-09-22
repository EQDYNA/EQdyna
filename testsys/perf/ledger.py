#! /usr/bin/env python3
"""Append-only perf ledger: docs/perf_ledger.jsonl (PROJECT_RULES rule 19's
positive counterpart).

WHY. A perf number produced today is quoted in a report and lost, or lands in
a shared mutable `*_last.*` file that the next run overwrites (rule 19's
2026-09-21 incident nearly destroyed the data item 33's row cites). Three
findings of this campaign -- the run_scaling Fortran metric biased in jax's
favour for days (91.25 vs a true 64.74 ms/step at 16 ranks), a 47% tenancy
drift on a single-device jax reading (611.00 -> 895.98 ms/step in two days),
and item 47's 2.68x per-rank spread on identical work -- are all invisible in
one-off numbers and obvious in a time series. The ledger is that time series.

FORMAT: JSONL, one JSON object per line. Chosen over TSV because rows carry
provenance the schema will grow (per-rank lists, cpu sets, sync mode) that TSV
can only hold via an ad-hoc second encoding, while a JSONL append is
line-atomic, self-describing, and never forces a column migration that
rewrites history.

CONTRACT.
- The ledger is APPEND-ONLY. Never rewritten, never sorted in place, exactly
  one line per measurement. `testsys/regression/test_perf_ledger.py` enforces
  this mechanically: the committed (HEAD) content must be a byte-prefix of the
  working copy.
- The ledger is the INDEX over `docs/perf_snapshots/`, not a replacement:
  every row points (relative path) at the dated, immutable snapshot file that
  holds the full run record.
- Rule 19 applies to the ledger itself: a row is evidence because it carries
  the SHA that produced it, and because the appended line can never be
  silently replaced.
- CONCURRENT APPENDS: `append` takes fcntl.flock(LOCK_EX) on the ledger file
  and issues the whole line as a SINGLE os.write on an O_APPEND descriptor.
  flock serialises writers that use this module; O_APPEND + single-write means
  even a writer that somehow bypassed the lock could interleave lines, not
  bytes. Proven by the concurrency check in test_perf_ledger.py.
- Rule 2: a row missing a required field, or carrying an uninterpretable
  value, RAISES. Nothing is defaulted in. Fields a tool genuinely does not
  measure are explicit nulls from an allowed list, never absent keys.

Row schema (one measurement = one line):
  ts_utc            ISO-8601 Z -- when the line was APPENDED (for live rows,
                    the end of the run; for backfilled rows, the backfill).
  snapshot_date_local  the producing run's own 'date' field, box-local time.
  sha, host, tool, case, backend, ranks, ms_per_step
  rank_ms_min/max/mean, effective_cores, threads_per_rank   (nullable: not
                    every tool measures per-rank quantities)
  busy_ceiling      the ceiling the tool selected/refused cpus against.
  tenancy_busy/tenancy_total  whole-box cpus over that ceiling at run time
                    (nullable ONLY on backfilled rows -- old snapshots did not
                    record whole-box tenancy).
  metric, n_lo, n_hi   'per-step-by-difference' over (n_lo, n_hi) steps.
  snapshot          RELATIVE path under docs/perf_snapshots/ to the dated file.
  backfilled_from   present iff the row was seeded from a committed snapshot
                    file rather than written by the run itself.
Extra keys (policy, sync, cpus, effective_cores_per_rank, ...) are allowed;
required keys are not negotiable.

Backfill (ONLY from committed dated snapshot files, never from prose):
    python3 testsys/perf/ledger.py backfill docs/perf_snapshots/<file>.json ...
"""
import fcntl
import json
import os
import re
import sys
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(TESTSYS))
LEDGER_RELPATH = os.path.join('docs', 'perf_ledger.jsonl')
LEDGER = os.path.join(ROOT, LEDGER_RELPATH)
SNAPSHOT_DIR_RELPATH = os.path.join('docs', 'perf_snapshots')

BACKENDS = ('fortran', 'python-numpy', 'python-jax', 'python-jax-mpi')
METRICS = ('per-step-by-difference',)

REQUIRED = ('ts_utc', 'snapshot_date_local', 'sha', 'host', 'tool', 'case',
            'backend', 'ranks', 'ms_per_step', 'rank_ms_min', 'rank_ms_max',
            'rank_ms_mean', 'effective_cores', 'threads_per_rank',
            'busy_ceiling', 'tenancy_busy', 'tenancy_total', 'metric',
            'n_lo', 'n_hi', 'snapshot')
# Nullable = "this tool does not measure that quantity", recorded as an
# explicit null so absence-of-measurement is distinguishable from a bug.
NULLABLE = ('rank_ms_min', 'rank_ms_max', 'rank_ms_mean', 'effective_cores',
            'threads_per_rank')

_TS_RE = re.compile(r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$')
_SHA_RE = re.compile(r'^[0-9a-f]{7,40}$')


def validate(row):
    """Raise ValueError on the first uninterpretable thing about `row`.
    A row that passes here is meant to be readable YEARS later with no access
    to the session that produced it (rule 6: provenance travels with the
    number)."""
    if not isinstance(row, dict):
        raise ValueError('row must be a dict, got %r' % type(row))
    missing = [k for k in REQUIRED if k not in row]
    if missing:
        raise ValueError('row missing required field(s) %s -- a ledger row '
                         'without them is not interpretable later' % missing)
    backfilled = 'backfilled_from' in row
    for k in REQUIRED:
        if row[k] is None:
            if k in NULLABLE:
                continue
            if backfilled and k in ('tenancy_busy', 'tenancy_total',
                                    'snapshot_date_local'):
                continue
            raise ValueError('required field %r is null and not in the '
                             'nullable set %s' % (k, NULLABLE))
    if not _TS_RE.match(str(row['ts_utc'])):
        raise ValueError('ts_utc %r is not ISO-8601 UTC (YYYY-MM-DDTHH:MM:SSZ)'
                         % row['ts_utc'])
    if not _SHA_RE.match(str(row['sha'])):
        raise ValueError('sha %r is not a git sha -- rule 19: a perf row is '
                         'evidence only when pinned to a commit' % row['sha'])
    if row['backend'] not in BACKENDS:
        raise ValueError('backend %r not in %s' % (row['backend'], BACKENDS))
    if row['metric'] not in METRICS:
        raise ValueError('metric %r not in %s -- rule 5: never invent a '
                         'metric silently; widen METRICS deliberately'
                         % (row['metric'], METRICS))
    if not (isinstance(row['ranks'], int) and row['ranks'] >= 1):
        raise ValueError('ranks %r must be an int >= 1' % row['ranks'])
    ms = row['ms_per_step']
    if not (isinstance(ms, (int, float)) and ms > 0):
        raise ValueError('ms_per_step %r must be a positive number (a '
                         'non-positive per-step figure is an invalid '
                         'measurement, not a slow one)' % ms)
    for k in ('n_lo', 'n_hi'):
        if not (isinstance(row[k], int) and row[k] >= 1):
            raise ValueError('%s %r must be an int >= 1' % (k, row[k]))
    if not row['n_lo'] < row['n_hi']:
        raise ValueError('n_lo %r must be < n_hi %r for a by-difference '
                         'metric' % (row['n_lo'], row['n_hi']))
    bc = row['busy_ceiling']
    if not (isinstance(bc, (int, float)) and 0 <= bc <= 1):
        raise ValueError('busy_ceiling %r must be a fraction in [0, 1]' % bc)
    tb, tt = row['tenancy_busy'], row['tenancy_total']
    if (tb is None) != (tt is None):
        raise ValueError('tenancy_busy/tenancy_total must be both set or '
                         'both null, got %r/%r' % (tb, tt))
    if tb is not None:
        if not (isinstance(tb, int) and isinstance(tt, int)
                and 0 <= tb <= tt and tt > 0):
            raise ValueError('tenancy %r/%r must be ints, 0 <= busy <= total'
                             % (tb, tt))
    snap = row['snapshot']
    if (not isinstance(snap, str) or os.path.isabs(snap) or '..' in snap
            or not snap.replace(os.sep, '/').startswith('docs/perf_snapshots/')):
        raise ValueError('snapshot %r must be a relative path under '
                         'docs/perf_snapshots/ -- the ledger is the index '
                         'over the dated snapshot files' % snap)
    mn, mx, mean = row['rank_ms_min'], row['rank_ms_max'], row['rank_ms_mean']
    if mn is not None and mx is not None and mean is not None:
        if not (mn <= mean <= mx):
            raise ValueError('rank_ms min/mean/max %r/%r/%r are not ordered'
                             % (mn, mean, mx))


def append(row, path=LEDGER):
    """Validate, then append ONE line: flock(LOCK_EX) around a single
    os.write on an O_APPEND descriptor. Raises on any short write."""
    validate(row)
    line = json.dumps(row, sort_keys=True, separators=(',', ':'))
    if '\n' in line or '\r' in line:
        raise ValueError('row serialised with an embedded newline -- refusing '
                         'to corrupt a line-oriented ledger')
    data = (line + '\n').encode('utf-8')
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        n = os.write(fd, data)
        if n != len(data):
            raise OSError('short write to %s: %d of %d bytes -- the ledger '
                          'may now hold a torn line; do not retry blindly, '
                          'inspect it' % (path, n, len(data)))
    finally:
        os.close(fd)   # releases the flock


def append_rows(rows, path=LEDGER):
    for r in rows:
        append(r, path)
    return len(rows)


def utc_now():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def box_tenancy(ceiling):
    """Whole-box tenancy right now: how many cpus are over `ceiling`, out of
    how many read. Reuses run_numa_scaling's measured busy fractions (one bug
    fixed once). Raises if per-cpu utilisation cannot be read at all -- a
    tenancy field defaulted in would be exactly the silent fallback rule 2
    forbids."""
    sys.path.insert(0, TESTSYS)
    import run_numa_scaling as numa
    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology -- cannot '
                         'record box tenancy for the perf ledger.')
    all_cpus = sorted(c for cs in nodes.values() for c in cs)
    busy = numa.cpu_busy_fractions(all_cpus)
    vals = [b for b in (busy or {}).values() if b is not None]
    if not vals:
        raise SystemExit('FAIL: could not read per-cpu utilisation from '
                         '/proc/stat -- cannot record box tenancy.')
    return dict(busy=sum(1 for b in vals if b > ceiling), total=len(vals))


def _base(meta, tool, snapshot, tenancy, backfilled_from):
    row = dict(ts_utc=utc_now(), snapshot_date_local=meta.get('date'),
               sha=meta['sha'], host=meta['host'], tool=tool,
               case=meta['case'], metric='per-step-by-difference',
               snapshot=snapshot,
               tenancy_busy=None if tenancy is None else tenancy['busy'],
               tenancy_total=None if tenancy is None else tenancy['total'])
    if backfilled_from is not None:
        row['backfilled_from'] = backfilled_from
    return row


def rows_from_scaling_snapshot(meta, snapshot, tenancy, backfilled_from=None):
    """Ledger rows from a run_scaling.py snapshot dict (its `meta`).
    One row per measured (engine, n, policy) point; skipped configs produce
    no row -- they are recorded in the snapshot itself."""
    out = []
    for r in meta['rows']:
        row = _base(meta, 'run_scaling', snapshot, tenancy, backfilled_from)
        row.update(backend=r['engine'], ranks=r['n'],
                   ms_per_step=r['ms_per_step'],
                   # run_scaling measures aggregate per-step only; per-rank
                   # quantities are explicit nulls, not zeros.
                   rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
                   effective_cores=None, threads_per_rank=None,
                   busy_ceiling=meta['busy_ceiling'],
                   n_lo=r['n_lo'], n_hi=r['n_hi'],
                   policy=r['policy'], cpus=r['cpus'])
        validate(row)
        out.append(row)
    return out


def rows_from_mpi_scaling_snapshot(meta, snapshot, tenancy, backfilled_from=None):
    """Ledger rows from a run_mpi_scaling.py snapshot dict. One row per
    (ranks, sync) jax-mpi point plus one per same-session fortran point."""
    out = []
    for r in meta['rows']:
        if r.get('skipped'):
            continue
        for k in sorted(r):
            if not k.startswith('jax_'):
                continue
            b = r[k]
            rank_ms = b['rank_ms']
            row = _base(meta, 'run_mpi_scaling', snapshot, tenancy,
                        backfilled_from)
            row.update(backend='python-jax-mpi', ranks=r['ranks'],
                       ms_per_step=b['ms_per_step'],
                       rank_ms_min=min(rank_ms), rank_ms_max=max(rank_ms),
                       rank_ms_mean=sum(rank_ms) / len(rank_ms),
                       effective_cores=sum(b['eff']) / len(b['eff']),
                       effective_cores_per_rank=b['eff'],
                       threads_per_rank=max(b['threads']),
                       busy_ceiling=meta['max_busy'],
                       n_lo=r['n_lo'], n_hi=r['n_hi'],
                       sync=k[len('jax_'):], cpus=r['cpus'])
            validate(row)
            out.append(row)
        if 'fortran' in r:
            f = r['fortran']
            row = _base(meta, 'run_mpi_scaling', snapshot, tenancy,
                        backfilled_from)
            row.update(backend='fortran', ranks=r['ranks'],
                       ms_per_step=f['ms_per_step'],
                       rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
                       effective_cores=None, threads_per_rank=None,
                       busy_ceiling=meta['max_busy'],
                       n_lo=r['n_lo'], n_hi=r['n_hi'], cpus=r['cpus'])
            validate(row)
            out.append(row)
    return out


def backfill(snapshot_path, path=LEDGER):
    """Seed ledger rows from ONE committed dated snapshot file. The snapshot
    file is the provenance; prose (session logs, board rows) is never a
    backfill source. Returns the number of rows appended (0 is a valid
    answer: a snapshot whose configs were all SKIPPED holds no measurement)."""
    rel = os.path.relpath(os.path.abspath(snapshot_path), ROOT)
    if not rel.replace(os.sep, '/').startswith('docs/perf_snapshots/'):
        raise SystemExit('FAIL: backfill source %s is not under '
                         'docs/perf_snapshots/ -- only committed dated '
                         'snapshots have the provenance a ledger row needs.'
                         % snapshot_path)
    meta = json.load(open(snapshot_path))
    if 'max_busy' in meta and 'busy_ceiling' not in meta:
        rows = rows_from_mpi_scaling_snapshot(meta, rel, None,
                                              backfilled_from=rel)
    elif 'busy_ceiling' in meta and 'max_busy' not in meta:
        rows = rows_from_scaling_snapshot(meta, rel, None,
                                          backfilled_from=rel)
    else:
        raise SystemExit('FAIL: cannot tell which tool wrote %s (looked for '
                         'exactly one of max_busy/busy_ceiling); refusing to '
                         'guess a schema.' % snapshot_path)
    return append_rows(rows, path)


def main(argv):
    if len(argv) >= 3 and argv[1] == 'backfill':
        total = 0
        for p in argv[2:]:
            n = backfill(p)
            print('%s: %d row(s) appended (backfilled) to %s' %
                  (p, n, LEDGER_RELPATH))
            total += n
        print('backfill done: %d row(s) total' % total)
        return 0
    print('usage: python3 testsys/perf/ledger.py backfill '
          'docs/perf_snapshots/<file>.json [...]')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
