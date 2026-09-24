#! /usr/bin/env python3
"""Small CLI over docs/run_profiles.jsonl (testsys/profile_record.py's
append-only record). Three questions, three subcommands -- no general query
language, because the three are all this landing was asked for:

  buckets   bucket breakdown for one case, grouped by sha (did a change move
            time between buckets, or just noise).
  compare   backend-vs-backend bucket breakdown at one rank count.
  diff      per-step-by-difference of buckets between two step counts, for
            one (case, backend, ranks) -- the project's standard metric
            (see ledger.py / run_mpi_scaling.py's own differencing, and
            CLAUDE.md's "Per-step-by-difference breaks at >=4 MPI ranks"
            note: this differences EACH RANK'S OWN total_s, never a wrapper
            wall clock, and refuses a non-positive result rather than
            reporting it).

Every subcommand reads the record file read-only; none of them writes to it.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))          # testsys/perf/
TESTSYS = os.path.dirname(HERE)
sys.path.insert(0, TESTSYS)
import profile_record  # noqa: E402

RECORD_PATH = profile_record.RECORD_PATH


def load_rows(path=RECORD_PATH, exclude_dirty=False):
    """Every profile_record.py row carries `tree_dirty` (required on that
    schema, unlike the legacy-tolerant perf ledger). `exclude_dirty=True`
    drops any row whose working tree was locally modified at capture time --
    a row from an uncommitted tree is not evidence about the sha it happens
    to be stamped with."""
    if not os.path.exists(path):
        return []
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    if exclude_dirty:
        rows = [r for r in rows if not r.get('tree_dirty', False)]
    return rows


def filter_rows(rows, **kw):
    out = rows
    for k, v in kw.items():
        if v is None:
            continue
        out = [r for r in out if r.get(k) == v]
    return out


def _bucket_means(rows):
    keys = sorted(set().union(*(r['buckets_s'] for r in rows))) if rows else []
    n = len(rows)
    return {k: sum(r['buckets_s'].get(k, 0.0) for r in rows) / n for k in keys}


def cmd_buckets(a):
    rows = filter_rows(load_rows(a.path, a.exclude_dirty), case=a.case,
                       backend=a.backend, ranks=a.ranks)
    if a.sha:
        rows = [r for r in rows if r['sha'] in a.sha]
    if not rows:
        raise SystemExit('no rows match case=%r backend=%r ranks=%r sha=%r'
                         % (a.case, a.backend, a.ranks, a.sha))
    by_sha = {}
    for r in rows:
        by_sha.setdefault(r['sha'], []).append(r)
    for sha in sorted(by_sha):
        rs = by_sha[sha]
        means = _bucket_means(rs)
        print('sha=%-10s n_rows=%-3d ' % (sha, len(rs))
              + '  '.join('%s=%.4f' % (k, means[k]) for k in sorted(means)))
    return 0


def cmd_compare(a):
    rows_a = filter_rows(load_rows(a.path, a.exclude_dirty), case=a.case,
                         backend=a.backend_a, ranks=a.ranks)
    rows_b = filter_rows(load_rows(a.path, a.exclude_dirty), case=a.case,
                         backend=a.backend_b, ranks=a.ranks)
    if not rows_a or not rows_b:
        raise SystemExit(
            'missing rows for case=%r ranks=%r: backend=%r has %d row(s), '
            'backend=%r has %d row(s)' % (a.case, a.ranks, a.backend_a,
                                          len(rows_a), a.backend_b, len(rows_b)))
    ma, mb = _bucket_means(rows_a), _bucket_means(rows_b)
    keys = sorted(set(ma) | set(mb))
    print('%-10s %-16s %-16s %s' % ('bucket', a.backend_a, a.backend_b,
                                    'ratio(%s/%s)' % (a.backend_b, a.backend_a)))
    for k in keys:
        va, vb = ma.get(k, 0.0), mb.get(k, 0.0)
        ratio = (vb / va) if va else float('nan')
        print('%-10s %-16.4f %-16.4f %.3f' % (k, va, vb, ratio))
    return 0


def cmd_diff(a):
    rows = filter_rows(load_rows(a.path, a.exclude_dirty), case=a.case,
                       backend=a.backend, ranks=a.ranks)
    lo = {r['rank']: r for r in rows if r['nsteps'] == a.n_lo}
    hi = {r['rank']: r for r in rows if r['nsteps'] == a.n_hi}
    common = sorted(set(lo) & set(hi))
    if not common:
        raise SystemExit(
            'need rows at both nsteps=%d (%d found) and nsteps=%d (%d found) '
            'sharing at least one rank id, for case=%r backend=%r ranks=%r'
            % (a.n_lo, len(lo), a.n_hi, len(hi), a.case, a.backend, a.ranks))
    dn = a.n_hi - a.n_lo
    for rank in common:
        lo_r, hi_r = lo[rank], hi[rank]
        total_ms = (hi_r['total_s'] - lo_r['total_s']) * 1000.0 / dn
        if total_ms <= 0:
            raise SystemExit(
                'rank %d: non-positive per-step total (%.6g ms/step) between '
                'nsteps=%d and nsteps=%d -- per CLAUDE.md convention this is a '
                'bug to raise on, not a number to report' % (rank, total_ms,
                                                              a.n_lo, a.n_hi))
        keys = sorted(set(lo_r['buckets_s']) | set(hi_r['buckets_s']))
        parts = ['%s=%.4fms' % (k, (hi_r['buckets_s'].get(k, 0.0)
                                    - lo_r['buckets_s'].get(k, 0.0))
                                 * 1000.0 / dn)
                 for k in keys]
        print('rank=%-3d total=%.4fms/step  %s' % (rank, total_ms,
                                                    ' '.join(parts)))
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--path', default=RECORD_PATH,
                   help='docs/run_profiles.jsonl to read (default: the '
                        'checked-out repo\'s own)')
    p.add_argument('--exclude-dirty', action='store_true',
                   help='drop rows captured from a locally modified '
                        'src/testsys tree -- such a row is evidence about '
                        'an uncommitted tree, not about the sha it is '
                        'stamped with (2026-09-23 incident)')
    sub = p.add_subparsers(dest='cmd', required=True)

    b = sub.add_parser('buckets', help='bucket breakdown for a case, by sha')
    b.add_argument('case')
    b.add_argument('--backend')
    b.add_argument('--ranks', type=int)
    b.add_argument('--sha', action='append',
                   help='restrict to this sha (repeatable); default: all shas present')
    b.set_defaults(func=cmd_buckets)

    c = sub.add_parser('compare', help='backend-vs-backend buckets at N ranks')
    c.add_argument('case')
    c.add_argument('ranks', type=int)
    c.add_argument('backend_a')
    c.add_argument('backend_b')
    c.set_defaults(func=cmd_compare)

    d = sub.add_parser('diff', help='per-step-by-difference of buckets '
                                     'between two step counts')
    d.add_argument('case')
    d.add_argument('backend')
    d.add_argument('ranks', type=int)
    d.add_argument('n_lo', type=int)
    d.add_argument('n_hi', type=int)
    d.set_defaults(func=cmd_diff)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == '__main__':
    sys.exit(main())
