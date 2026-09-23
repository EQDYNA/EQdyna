#! /usr/bin/env python3
"""Append-only, durable record of every run's per-rank profile:
docs/run_profiles.jsonl (sibling to docs/perf_ledger.jsonl).

WHY A SEPARATE FILE, NOT A ROW SHAPE ON THE EXISTING LEDGER.
`testsys/perf/ledger.py` is one measurement per line: a single per-step or
wall-clock NUMBER, with the provenance that number needs. A profiled run
emits `nranks` per-rank BUCKET BREAKDOWNS at once, a different cardinality
and a different question ("where did the time go on THIS rank" vs "how fast
was this cell"). Folding one into the other would either lose the per-rank
rows or force every ledger reader to filter out a shape it does not expect.
The two files share exactly what already applies to both: append-only,
line-atomic, validated-before-write, and requiring the sha/host/timestamp
provenance rule 6 asks of every number in this repo.

CONTRACT (identical to ledger.py's, restated because this module does not
import ledger.append/validate directly -- the row shapes differ):
  - APPEND-ONLY. `append()`/`append_rows()` take fcntl.flock(LOCK_EX) on the
    record file and issue the whole line as a SINGLE os.write on an
    O_APPEND descriptor -- the same concurrency-safety argument ledger.py
    makes, proved the same way in the regression test.
  - Rule 2: a row missing a required field, or carrying an uninterpretable
    value, RAISES. Nothing is defaulted in.
  - A run that finishes without its rows appended is a FAILURE: `capture_run`
    raises on any problem (missing/invalid profile files, a tenancy read
    that cannot be measured, a short write) rather than degrading to a
    warning. Unlike `ledger.capture_e2e_cells_or_warn`, there is no `_or_warn`
    wrapper here -- profiling data is the deliverable of this landing, not a
    free side-measurement of a gate that must stay green regardless, so the
    harness is meant to turn this raise into a failed cell.

REUSED FROM ledger.py, READ-ONLY (no edit to that module):
  BACKENDS, TENANCY_REFERENCE_CEILING, box_tenancy(), utc_now(),
  the sha/timestamp regexes. One definition of "what is a git sha" and one
  definition of "what counts as whole-box tenancy" between the two records,
  not two that can drift apart.

Row schema, one per (case, backend, run, rank):
  schema            'eqdyna-run-profile/1'
  ts_utc, sha       provenance, same regexes ledger.py enforces.
  case, backend     as in ledger.py.
  term              the term-axis value this run belongs to (opaque here --
                     the term axis itself lands separately; this module only
                     requires it be present and non-null, never guessed).
  ranks             total MPI ranks in this run (== nranks reported by the
                     profile files).
  rank              this row's own rank id, 0 <= rank < ranks.
  host, pid         from the rank's own profile.rank<r>.json.
  cpus_allowed, numa_nodes   this rank's affinity, from its own profile file.
  ranks_per_node    rank COUNT per host, derived across every rank file of
                     THIS run (not per-rank data -- one shared list on every
                     row of the run, so a reader can see the placement
                     without re-deriving it from the other rows).
  nsteps, sampling_every, buckets_s, loop_s, total_s, unaccounted_s
                     carried straight from the validated profile.rank<r>.json
                     (bucket TOTALS only -- no per-step trace, by design: the
                     per-rank bucket file is already the fine-grained record;
                     this ledger is the durable cross-run index over it).
  tenancy_busy, tenancy_total, busy_ceiling
                     whole-box tenancy at capture time, via ledger.box_tenancy
                     -- "host tenancy" reusing exactly what the perf ledger
                     already measures, not a second implementation of it.
"""
import fcntl
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))          # testsys/
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(HERE)
PERF = os.path.join(HERE, 'perf')
sys.path.insert(0, PERF)
import ledger  # noqa: E402  -- read-only reuse: BACKENDS, box_tenancy,
                              # utc_now, TENANCY_REFERENCE_CEILING, the
                              # sha/ts regexes. Not edited by this module.
sys.path.insert(0, HERE)
import profile_schema  # noqa: E402

RECORD_SCHEMA_ID = 'eqdyna-run-profile/1'
RECORD_RELPATH = os.path.join('docs', 'run_profiles.jsonl')
RECORD_PATH = os.path.join(ROOT, RECORD_RELPATH)

REQUIRED = ('schema', 'ts_utc', 'sha', 'case', 'backend', 'term', 'ranks',
            'rank', 'host', 'pid', 'cpus_allowed', 'numa_nodes',
            'ranks_per_node', 'nsteps', 'sampling_every', 'buckets_s',
            'loop_s', 'total_s', 'unaccounted_s', 'tenancy_busy',
            'tenancy_total', 'busy_ceiling')


def _row_identity(row):
    return ('case=%r backend=%r ranks=%r rank=%r term=%r'
            % (row.get('case'), row.get('backend'), row.get('ranks'),
               row.get('rank'), row.get('term')))


def validate_row(row):
    """Raise ValueError on the first uninterpretable thing about one
    docs/run_profiles.jsonl row. Re-checks the bucket/sum arithmetic
    independently (defense in depth: this is the DURABLE copy, and it must
    not simply trust that profile_schema.validate() was actually called on
    the source file before this row was built)."""
    if not isinstance(row, dict):
        raise ValueError('row must be a dict, got %r' % type(row))
    missing = [k for k in REQUIRED if k not in row]
    if missing:
        raise ValueError('row missing required field(s) %s -- a record row '
                         'without them is not interpretable later' % missing)

    if row['schema'] != RECORD_SCHEMA_ID:
        raise ValueError('schema %r != %r' % (row['schema'], RECORD_SCHEMA_ID))
    if not ledger._TS_RE.match(str(row['ts_utc'])):
        raise ValueError('ts_utc %r is not ISO-8601 UTC (YYYY-MM-DDTHH:MM:SSZ)'
                         % row['ts_utc'])
    if not ledger._SHA_RE.match(str(row['sha'])):
        raise ValueError('sha %r is not a git sha -- a profile row is '
                         'evidence only when pinned to a commit' % row['sha'])
    if not (isinstance(row['case'], str) and row['case'].strip()):
        raise ValueError('case %r must be a non-empty string' % row['case'])
    if row['backend'] not in ledger.BACKENDS:
        raise ValueError('backend %r not in %s' % (row['backend'], ledger.BACKENDS))
    if row['term'] is None:
        raise ValueError('term is null on row from %s -- every row must say '
                         'which term-axis value produced it; no value is '
                         'defaulted in' % _row_identity(row))

    profile_schema.check_nonneg_int('run_profiles row', 'ranks', row['ranks'], 1)
    profile_schema.check_nonneg_int('run_profiles row', 'rank', row['rank'], 0)
    if row['rank'] >= row['ranks']:
        raise ValueError('rank %r must be < ranks %r on row from %s'
                         % (row['rank'], row['ranks'], _row_identity(row)))
    if not (isinstance(row['host'], str) and row['host'].strip()):
        raise ValueError('host %r must be a non-empty string' % row['host'])
    profile_schema.check_nonneg_int('run_profiles row', 'pid', row['pid'], 1)

    profile_schema.check_int_list('run_profiles row', 'cpus_allowed', row['cpus_allowed'])
    profile_schema.check_int_list('run_profiles row', 'numa_nodes', row['numa_nodes'])
    # ranks_per_node is a list of COUNTS, one per host -- repeats are
    # expected (two hosts each running 2 ranks is [2, 2]), so unique=False.
    profile_schema.check_int_list('run_profiles row', 'ranks_per_node',
                                  row['ranks_per_node'], unique=False)

    profile_schema.check_nonneg_int('run_profiles row', 'nsteps', row['nsteps'], 1)
    profile_schema.check_nonneg_int('run_profiles row', 'sampling_every',
                                    row['sampling_every'], 1)
    if not (isinstance(row['loop_s'], (int, float)) and not isinstance(row['loop_s'], bool)
            and row['loop_s'] >= 0):
        raise ValueError('loop_s %r must be a number >= 0' % row['loop_s'])
    profile_schema.check_buckets('run_profiles row', row['total_s'],
                                 row['buckets_s'], row['unaccounted_s'])

    tb, tt = row['tenancy_busy'], row['tenancy_total']
    if not (isinstance(tb, int) and isinstance(tt, int) and 0 <= tb <= tt and tt > 0):
        raise ValueError('tenancy %r/%r must be ints, 0 <= busy <= total, total > 0'
                         % (tb, tt))
    bc = row['busy_ceiling']
    if not (isinstance(bc, (int, float)) and 0 <= bc <= 1):
        raise ValueError('busy_ceiling %r must be a fraction in [0, 1]' % bc)


def append(row, path=RECORD_PATH):
    """Validate, then append ONE line: flock(LOCK_EX) around a single
    os.write on an O_APPEND descriptor. Raises on any short write. Mirrors
    ledger.append's concurrency argument exactly (same lock discipline, same
    single-write discipline) -- see that module's docstring for why each
    piece is load-bearing."""
    validate_row(row)
    line = json.dumps(row, sort_keys=True, separators=(',', ':'))
    if '\n' in line or '\r' in line:
        raise ValueError('row serialised with an embedded newline -- refusing '
                         'to corrupt a line-oriented record')
    data = (line + '\n').encode('utf-8')
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        n = os.write(fd, data)
        if n != len(data):
            raise OSError('short write to %s: %d of %d bytes -- the record '
                          'may now hold a torn line; do not retry blindly, '
                          'inspect it' % (path, n, len(data)))
    finally:
        os.close(fd)   # releases the flock


def append_rows(rows, path=RECORD_PATH):
    for r in rows:
        append(r, path)
    return len(rows)


def _ranks_per_node(per_rank_rows):
    """[int, ...]: rank count per host, one entry per host, sorted by host
    name for a deterministic list regardless of rank-file iteration order.
    Repeats are expected and kept (two hosts each running 2 ranks -> [2, 2])."""
    counts = {}
    for row in per_rank_rows:
        counts[row['host']] = counts.get(row['host'], 0) + 1
    return [counts[h] for h in sorted(counts)]


def capture_run(run_dir, *, case, backend, ranks, term, sha):
    """Validate every profile.rank<r>.json under `run_dir`, then append one
    docs/run_profiles.jsonl row per rank. Returns the list of rows appended.

    RAISES on any problem -- missing/invalid profile files, a caller-claimed
    `ranks`/`backend` that disagrees with what the files themselves report, an
    unmeasurable tenancy read, or a short write -- so that a run finishing
    without its rows appended is a FAILURE the harness can turn into a failed
    cell, exactly as the calling contract requires. There is no `_or_warn`
    variant: unlike a free perf-ledger side-measurement, this data IS the
    deliverable."""
    per_rank = profile_schema.validate_run_dir(run_dir)   # raises on any problem
    nranks = per_rank[0]['nranks']
    if nranks != ranks:
        raise ValueError(
            'capture_run: run_dir %r profiles report nranks=%d but caller '
            'passed ranks=%r for %s x %s -- the harness and the profile files '
            'disagree on how many ranks this run had' % (run_dir, nranks,
                                                          ranks, case, backend))
    file_backend = per_rank[0]['backend']
    if file_backend != backend:
        raise ValueError(
            'capture_run: run_dir %r profiles report backend=%r but caller '
            'passed backend=%r for case %r' % (run_dir, file_backend, backend, case))
    if not (isinstance(sha, str) and ledger._SHA_RE.match(sha)):
        raise ValueError('capture_run: sha %r is not a git sha' % (sha,))
    if term is None:
        raise ValueError('capture_run: term must not be None -- every row '
                         'must say which term-axis value produced it')

    tenancy = ledger.box_tenancy(ledger.TENANCY_REFERENCE_CEILING)  # raises if
                                                                     # unmeasurable
    rpn = _ranks_per_node(per_rank)
    ts = ledger.utc_now()

    rows = []
    for p in per_rank:
        row = dict(
            schema=RECORD_SCHEMA_ID, ts_utc=ts, sha=sha, case=case,
            backend=backend, term=term, ranks=nranks, rank=p['rank'],
            host=p['host'], pid=p['pid'],
            cpus_allowed=list(p['cpus_allowed']),
            numa_nodes=list(p['numa_nodes']),
            ranks_per_node=rpn,
            nsteps=p['nsteps'], sampling_every=p['sampling_every'],
            buckets_s=dict(p['buckets_s']),
            loop_s=p['loop_s'], total_s=p['total_s'],
            unaccounted_s=p['unaccounted_s'],
            tenancy_busy=tenancy['busy'], tenancy_total=tenancy['total'],
            busy_ceiling=ledger.TENANCY_REFERENCE_CEILING,
        )
        validate_row(row)
        rows.append(row)

    n = append_rows(rows)
    if n != len(per_rank):
        raise RuntimeError(
            'capture_run: appended %d row(s) but expected %d for %s x %s -- '
            'a run finishing without ALL its rows appended is a failure'
            % (n, len(per_rank), case, backend))
    return rows
