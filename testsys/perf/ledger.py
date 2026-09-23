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
  metric            'per-step-by-difference' (n_lo/n_hi required ints) or
                    'cell-wall-clock' (one e2e cell's whole wall time: setup +
                    solve + compare + XLA compile). THE TWO ARE NOT COMPARABLE
                    -- mixing compile time into a per-step verdict once
                    produced a spurious 33% JAX regression -- so a
                    cell-wall-clock row carries its value in `wall_s` and MUST
                    carry ms_per_step/n_lo/n_hi as explicit nulls: there is no
                    ms_per_step number on it to average against a per-step row.
  n_lo, n_hi        the two step counts of a by-difference metric.
  snapshot          RELATIVE path under docs/perf_snapshots/ to the dated file.
  backfilled_from   present iff the row was seeded from a committed snapshot
                    file rather than written by the run itself.
  platform          'cpu' | 'gpu' | 'unknown' -- what ACTUALLY ran, derived
                    from the snapshot's own recorded device evidence (nvidia-smi
                    per-device memory deltas), NEVER from the requested label
                    and never defaulted. The 2026-09-21 812.78 ms/step row is
                    why: a CPU measurement recorded under a cuda label
                    (eqdyna3d --device auto overriding JAX_PLATFORMS=cuda),
                    0 MiB delta on all four A100s. 'unknown' is for a snapshot
                    that does not say, or whose label contradicts its evidence
                    -- distinguishable from 'cpu' by construction.
                    REQUIRED on every row this module APPENDS from now on.
                    The 11 rows appended before 2026-09-22 lack it (the ledger
                    is append-only and is never rewritten); a row WITHOUT the
                    key is a legacy row and must be read as platform UNKNOWN,
                    never as cpu.
  parallelism       'mpi' | 'threads' | 'serial' -- WHAT THE `ranks` NUMBER
                    ON THIS ROW COUNTS. Without it one column carried two
                    incommensurable quantities: on a run_scaling python row
                    `ranks` is the cpus in ONE process's affinity mask (XLA and
                    OpenBLAS auto-size their pools from the numactl pin;
                    run_scaling.time_one_py launches a single process), while
                    on a run_scaling fortran row or any run_mpi_scaling row it
                    is `mpirun -np n` SEPARATE PROCESSES. Reading the two as
                    one column published "Fortran 16.24x vs jax 2.59x at 16"
                    on 2026-09-23 -- largely 16 MPI ranks against 16 threads in
                    one process; the real-ranks measurement put jax at 1.58x of
                    Fortran, not 4.19x.
                      'mpi'     ranks = MPI processes, separate address spaces.
                      'threads' ranks = cpus made available to ONE process
                                (its affinity mask); library thread pools
                                auto-size from it. No MPI.
                      'serial'  one process, no parallel axis selected or
                                measured; ranks is 1 by construction and says
                                NOTHING about how many cores the run used.
                    REQUIRED on every row this module APPENDS from now on, on
                    the same terms as `platform`: rows appended before
                    2026-09-23 lack the key (the ledger is append-only and is
                    never rewritten), and a row WITHOUT the key is a legacy row
                    whose parallelism is UNKNOWN -- never read as 'mpi', which
                    is exactly the defect. So a consumer tests
                    `'parallelism' in row`, not `row.get('parallelism', ...)`:
                    presence of the key is the old/new discriminator, and the
                    two populations must not be plotted on one axis.
  platform_evidence what proved `platform` (a sentence naming the measurement).
  devices           int count of GPUs shown busy by the memory evidence, on
                    'gpu' rows; explicit null otherwise. This is the companion
                    that ends the `ranks` ambiguity: `ranks` always counts MPI
                    processes; `devices` counts GPUs those processes provably
                    touched.
  gpu_peak_mib / gpu_delta_mib   nvidia-smi per-device peak / delta-over-idle
                    (MiB) sampled during the run; device_peak_gb the per-rank
                    jax allocator peak. Carried onto gpu rows (capacity claims
                    are answered from here), null/absent elsewhere.
  wall_s            the value field of a cell-wall-clock row (seconds).
  supersedes        {line, ts_utc, reason} -- this row RESTATES, with honest
                    identity, the measurement at that 1-based line of this
                    ledger. Append-only supersession: the old line stays
                    byte-identical, the newest row for a (snapshot, ranks,
                    sync) wins. Written by `reissue`, which refuses unless the
                    old line's snapshot/case/backend/ranks/sync/ms_per_step
                    match the re-derived row exactly.
Extra keys (policy, sync, cpus, effective_cores_per_rank, verdict, ...) are
allowed; required keys are not negotiable.

Backfill (ONLY from committed dated snapshot files, never from prose):
    python3 testsys/perf/ledger.py backfill docs/perf_snapshots/<file>.json ...
Reissue (append corrected identities for rows that predate `platform`):
    python3 testsys/perf/ledger.py reissue docs/perf_snapshots/<file>.json 4,5
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
METRICS = ('per-step-by-difference', 'cell-wall-clock')
PLATFORMS = ('cpu', 'gpu', 'unknown')
# What the `ranks` field of a row COUNTS. There is no 'unknown' member on
# purpose: a tool always knows how it launched its own run, so a row it cannot
# classify is a bug in the tool, not a third state. Legacy rows (before
# 2026-09-23) carry no key at all -- that absence, and only that absence, is
# read as unknown.
PARALLELISM = ('mpi', 'threads', 'serial')
# What chose the cpus on a run_mpi_scaling row. 'packed' and 'spread' are
# `run_mpi_scaling.py --placement`; 'explicit' is `--cpus` (the caller named
# the cpus, so neither selection algorithm ran). Scoped to tool='run_mpi_scaling'
# only -- run_scaling has its own, unrelated 'policy' field, and e2e cells
# select no cpus at all. Required on every run_mpi_scaling row APPENDED from
# 2026-09-23; absent on everything written before, and on every row from a
# different tool. See `_check_placement`.
PLACEMENT = ('packed', 'spread', 'explicit')
# The ceiling the e2e cell capture measures whole-box tenancy against. e2e
# selects no cpus, so unlike the scaling tools it has no selection ceiling;
# this is a declared reference for the tenancy statistic only, recorded on
# the row as busy_ceiling so the number stays interpretable.
TENANCY_REFERENCE_CEILING = 0.5

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


def _check_platform(row):
    """The platform identity checks. A row that CLAIMS a GPU must show the
    per-device memory evidence; a cpu/unknown row must not carry a device
    count. Guarded incident: 812.78 ms/step, a CPU measurement recorded under
    a cuda label (2026-09-21, mpi_scaling_..._gpu1.json)."""
    for k in ('platform', 'platform_evidence', 'devices'):
        if k not in row:
            raise ValueError('row carries platform identity but is missing '
                             '%r -- platform, platform_evidence and devices '
                             'travel together or not at all' % k)
    p = row['platform']
    if p not in PLATFORMS:
        raise ValueError('platform %r not in %s -- "what was requested" is '
                         'not a platform; only what the device evidence '
                         'proves, or unknown' % (p, PLATFORMS))
    ev = row['platform_evidence']
    if not (isinstance(ev, str) and ev.strip()):
        raise ValueError('platform_evidence %r must be a non-empty string '
                         'naming the measurement that proved platform=%r'
                         % (ev, p))
    if p == 'gpu':
        delta = row.get('gpu_delta_mib') or {}
        pos = [v for v in delta.values()
               if isinstance(v, (int, float)) and v > 0]
        if not pos:
            raise ValueError(
                "platform 'gpu' claimed without positive per-device memory "
                "evidence (gpu_delta_mib) -- a number that cannot show the "
                "GPU it ran on is the 812.78 ms/step incident in reverse; "
                "record it as 'unknown' or not at all")
        if not (isinstance(row['devices'], int) and row['devices'] >= 1):
            raise ValueError('devices %r must be an int >= 1 on a gpu row '
                             '(the count of devices the memory evidence '
                             'shows busy)' % row['devices'])
    elif row['devices'] is not None:
        raise ValueError('devices %r must be null on a %r row -- it counts '
                         'GPUs proven by memory evidence, nothing else'
                         % (row['devices'], p))


def _row_identity(row):
    """A short 'who wrote this' string for a refusal message, built only from
    fields already required, so it is safe to call on a row that failed."""
    return ('tool=%r case=%r backend=%r ranks=%r snapshot=%r'
            % (row.get('tool'), row.get('case'), row.get('backend'),
               row.get('ranks'), row.get('snapshot')))


def _check_parallelism(row):
    """Fail closed on the `ranks` discriminator. Guarded incident
    (2026-09-23): `ranks` meant MPI processes on fortran/run_mpi_scaling rows
    and cpus-in-one-affinity-mask on run_scaling python rows, with nothing on
    the row saying which, and the two were divided against each other into a
    published 4.19x port deficit that was really 1.58x.

    There is no default. A defaulted 'mpi' would reproduce the defect exactly:
    the threaded rows are the ones that would inherit it."""
    if 'parallelism' not in row:
        raise ValueError(
            'row is missing required field %r -- every appended row must say '
            'what its `ranks` number counts, one of %s; the emitter that '
            'built this row (%s) has not been wired. No value is defaulted '
            'in: an assumed "mpi" on a threaded row is the 2026-09-23 defect '
            'this field exists to end.'
            % ('parallelism', PARALLELISM, _row_identity(row)))
    p = row['parallelism']
    if p not in PARALLELISM:
        raise ValueError(
            'parallelism %r not in %s on row from %s -- widen PARALLELISM '
            'deliberately, with the reader contract updated, or fix the '
            'emitter' % (p, PARALLELISM, _row_identity(row)))
    if p == 'serial' and row['ranks'] != 1:
        raise ValueError(
            "parallelism 'serial' with ranks=%r on row from %s -- serial "
            'means one process and no parallel axis, so ranks is 1 by '
            'construction' % (row['ranks'], _row_identity(row)))


def _check_placement(row):
    """Fail closed on the `placement` discriminator for run_mpi_scaling rows
    (2026-09-23). Guarded incident: least_loaded_cpus sorts by (busy, node,
    cpu); on a quiet box the node tiebreak dominates and it PACKS ranks onto
    the fewest NUMA nodes, and a whole day of 8/16-rank tables carried that
    bias with nothing on the row saying so. test.tpv104 16 ranks, same
    session: packed 96.0/56.0 ms/step = 1.71x jax/Fortran, spread 52.5/63.3 =
    0.83x (docs/perf_ledger.jsonl lines 361-364) -- opposite verdicts, same
    code, same rank count.

    Scoped to tool='run_mpi_scaling': that is the only tool this field
    describes (run_scaling has its own 'policy' field; e2e cells select no
    cpus). No value is defaulted in -- an assumed 'packed' would reproduce
    the exact defect this field exists to end."""
    if 'placement' not in row:
        raise ValueError(
            'row is missing required field %r -- every run_mpi_scaling row '
            'must say how its cpus were chosen, one of %s (row: %s). No '
            'value is defaulted in.' % ('placement', PLACEMENT,
                                        _row_identity(row)))
    p = row['placement']
    if p not in PLACEMENT:
        raise ValueError('placement %r not in %s on row from %s'
                         % (p, PLACEMENT, _row_identity(row)))
    rpn = row.get('ranks_per_node')
    if not (isinstance(rpn, list) and rpn
            and all(isinstance(x, int) and x >= 0 for x in rpn)):
        raise ValueError(
            'ranks_per_node %r must be a non-empty list of ints on a row '
            'carrying placement=%r (row: %s)'
            % (rpn, p, _row_identity(row)))


def validate(row, appending=False):
    """Raise ValueError on the first uninterpretable thing about `row`.
    A row that passes here is meant to be readable YEARS later with no access
    to the session that produced it (rule 6: provenance travels with the
    number).

    `appending=True` (what `append` uses) additionally REQUIRES the platform
    identity fields (from 2026-09-22) and `parallelism` (from 2026-09-23):
    every row written from then on says what actually ran and what its `ranks`
    number counts. With the default (reading history), a row without a `platform` key is
    accepted as a legacy row -- the 11 pre-2026-09-22 lines are append-only
    and will never be rewritten to gain the key -- and must be read as
    platform UNKNOWN, never as cpu. A legacy row that DOES carry the key is
    held to the full checks."""
    if not isinstance(row, dict):
        raise ValueError('row must be a dict, got %r' % type(row))
    missing = [k for k in REQUIRED if k not in row]
    if missing:
        raise ValueError('row missing required field(s) %s -- a ledger row '
                         'without them is not interpretable later' % missing)
    backfilled = 'backfilled_from' in row
    is_cell = row.get('metric') == 'cell-wall-clock'
    for k in REQUIRED:
        if row[k] is None:
            if k in NULLABLE:
                continue
            if backfilled and k in ('tenancy_busy', 'tenancy_total',
                                    'snapshot_date_local'):
                continue
            if is_cell and k in ('ms_per_step', 'n_lo', 'n_hi'):
                continue     # required-NULL there; enforced below
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
    if appending or 'platform' in row:
        _check_platform(row)
    # Same legacy contract as platform: REQUIRED when appending, tolerated as
    # absent (= UNKNOWN, never a default) when reading the 11+ rows that
    # predate the field, and fully checked on any row that does carry it.
    if appending or 'parallelism' in row:
        _check_parallelism(row)
    # Scoped to run_mpi_scaling: 'placement' describes how ITS cpus were
    # chosen and has no meaning for run_scaling (which has 'policy') or e2e
    # (which pins none). Required on append for that tool only; tolerated as
    # absent on read (legacy rows, and rows from every other tool).
    if row.get('tool') == 'run_mpi_scaling' and (appending
                                                 or 'placement' in row):
        _check_placement(row)
    if row['metric'] == 'per-step-by-difference':
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
    else:  # cell-wall-clock
        if row['ms_per_step'] is not None:
            raise ValueError(
                'a cell-wall-clock row must carry ms_per_step as an explicit '
                'null, got %r -- wall clock mixes XLA compile and setup into '
                'the number (the spurious 33%% JAX regression), so it must be '
                'STRUCTURALLY impossible to average against a per-step row; '
                'its value lives in wall_s' % row['ms_per_step'])
        w = row.get('wall_s')
        if not (isinstance(w, (int, float)) and w > 0):
            raise ValueError('cell-wall-clock row needs wall_s as a positive '
                             'number of seconds, got %r' % w)
        for k in ('n_lo', 'n_hi'):
            if row[k] is not None:
                raise ValueError('%s %r must be null on a cell-wall-clock '
                                 'row: there is no by-difference pair, and a '
                                 'fake one would invite the per-step reading'
                                 % (k, row[k]))
    sup = row.get('supersedes')
    if sup is not None:
        if not (isinstance(sup, dict) and isinstance(sup.get('line'), int)
                and sup['line'] >= 1 and isinstance(sup.get('reason'), str)
                and sup['reason'].strip()):
            raise ValueError('supersedes %r must be {line: int >= 1, ts_utc, '
                             'reason: non-empty} naming the ledger line it '
                             'restates' % (sup,))
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
    """Validate (strictly: a row being APPENDED must carry its platform
    identity), then append ONE line: flock(LOCK_EX) around a single os.write
    on an O_APPEND descriptor. Raises on any short write."""
    validate(row, appending=True)
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


FORTRAN_CPU_EVIDENCE = ('fortran solver: CPU-only, src/fortran has no GPU '
                        'path')


def mpi_block_platform(block):
    """(platform, evidence, devices, gpu_peak_mib, gpu_delta_mib,
    device_peak_gb) for one run_mpi_scaling jax_<sync> block, derived from the
    block's own recorded DEVICE EVIDENCE -- the nvidia-smi per-device memory
    delta sampled during the run -- never from the requested label alone.

    Why the label is not trusted: mpi_scaling_2026-09-21_tpv104_gpu1.json
    carries platform="cuda" and an 812.78 ms/step figure, and its own
    nvidia-smi record shows a 0 MiB delta on all four A100s (eqdyna3d's
    --device auto override forced CPU under JAX_PLATFORMS=cuda). That number
    is a CPU measurement; this function returns 'unknown' for it, which can
    never be mistaken for 'cpu' or averaged as a GPU point.
    A cpu label IS trusted: JAX_PLATFORMS=cpu pins the jax backend to host
    devices and cannot silently land on a GPU."""
    label = block.get('platform')
    gm = block.get('gpu_mem') or {}
    delta = dict(gm.get('delta') or {})
    peak = dict(gm['peak']) if gm.get('peak') else None
    dpg = block.get('device_peak_gb')
    if label is None:
        return ('unknown', 'snapshot block records neither a platform label '
                'nor device-memory evidence (predates the GPU mode)',
                None, None, None, None)
    busy_devs = sorted(k for k, v in delta.items() if v > 0)
    if label == 'cpu':
        if busy_devs:
            return ('unknown', 'label "cpu" contradicted by the evidence: '
                    'nvidia-smi read a positive memory delta on device(s) %s '
                    'during the run; recording neither cpu nor gpu'
                    % ','.join(busy_devs), None, peak, delta, dpg)
        return ('cpu', 'requested cpu: JAX_PLATFORMS=cpu pins the jax backend '
                'to host devices; it cannot land on a GPU', None, None, None,
                None)
    if busy_devs:
        return ('gpu', 'nvidia-smi memory delta > 0 MiB on device(s) %s '
                'sampled during the run' % ','.join(busy_devs),
                len(busy_devs), peak, delta, dpg)
    return ('unknown', 'label %r contradicted by the evidence: nvidia-smi '
            'read a 0 MiB memory delta on every device during the run (the '
            '2026-09-21 812.78 ms/step incident: a CPU measurement under a '
            'GPU label); recording neither cpu nor gpu' % label,
            None, peak, delta or None, dpg)


_SCALING_PARALLELISM_BY_ENGINE = {
    # run_scaling.py:325-356 builds `mpirun -np n` for the fortran engine:
    # n SEPARATE PROCESSES.
    'fortran': 'mpi',
    # run_scaling.py:453-478 (time_one_py) launches ONE `sys.executable -c`
    # under numactl with n cpus in its affinity mask; OpenBLAS (numpy) and XLA
    # (jax) auto-size their thread pools from that mask. The threaded mode is
    # a real execution mode and stays measured -- only the label was missing.
    'python-numpy': 'threads',
    'python-jax': 'threads',
}


def _SCALING_PARALLELISM(engine):
    """The discriminator for one run_scaling engine. Raises on an engine this
    function has not been taught -- a new engine must declare what its `ranks`
    counts before its first number is written down, not after someone divides
    it by a Fortran row."""
    try:
        return _SCALING_PARALLELISM_BY_ENGINE[engine]
    except KeyError:
        raise ValueError(
            'run_scaling engine %r has no declared parallelism (known: %s) '
            '-- refusing to guess what its `ranks` counts'
            % (engine, sorted(_SCALING_PARALLELISM_BY_ENGINE))) from None


def rows_from_scaling_snapshot(meta, snapshot, tenancy, backfilled_from=None):
    """Ledger rows from a run_scaling.py snapshot dict (its `meta`).
    One row per measured (engine, n, policy) point; skipped configs produce
    no row -- they are recorded in the snapshot itself.

    Platform is 'cpu' for every engine BY CONSTRUCTION of the tool, not by
    default: fortran has no GPU path, numpy has no GPU backend, and
    run_scaling pins JAX_PLATFORMS=cpu on its python runs (run_scaling.py's
    env setup), which jax cannot override onto a GPU."""
    _CPU_EV = {'fortran': FORTRAN_CPU_EVIDENCE,
               'python-numpy': 'numpy backend: host arrays only, no GPU path',
               'python-jax': 'run_scaling pins JAX_PLATFORMS=cpu; jax cannot '
                             'land on a GPU under that pin'}
    out = []
    for r in meta['rows']:
        row = _base(meta, 'run_scaling', snapshot, tenancy, backfilled_from)
        row.update(backend=r['engine'], ranks=r['n'],
                   ms_per_step=r['ms_per_step'],
                   parallelism=_SCALING_PARALLELISM(r['engine']),
                   platform='cpu', devices=None,
                   platform_evidence=_CPU_EV.get(
                       r['engine'], 'run_scaling is a CPU-only tool'),
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
        # Propagated only if the producing row carries it: a snapshot from
        # before 2026-09-23 has no 'placement' key at all, and that absence
        # must reach the ledger row as absence (validate()'s `appending=True`
        # then refuses it, same contract as platform/parallelism), never as
        # a guessed 'packed'.
        placement_fields = {k: r[k] for k in ('placement', 'ranks_per_node')
                            if k in r}
        for k in sorted(r):
            if not k.startswith('jax_'):
                continue
            b = r[k]
            rank_ms = b['rank_ms']
            plat, ev, ndev, peak, delta, dpg = mpi_block_platform(b)
            row = _base(meta, 'run_mpi_scaling', snapshot, tenancy,
                        backfilled_from)
            row.update(backend='python-jax-mpi', ranks=r['ranks'],
                       ms_per_step=b['ms_per_step'],
                       # every run_mpi_scaling point is `mpirun -np ranks`,
                       # one process per rank, on both backends.
                       parallelism='mpi',
                       platform=plat, platform_evidence=ev, devices=ndev,
                       gpu_peak_mib=peak, gpu_delta_mib=delta,
                       device_peak_gb=dpg,
                       rank_ms_min=min(rank_ms), rank_ms_max=max(rank_ms),
                       rank_ms_mean=sum(rank_ms) / len(rank_ms),
                       effective_cores=sum(b['eff']) / len(b['eff']),
                       effective_cores_per_rank=b['eff'],
                       threads_per_rank=max(b['threads']),
                       busy_ceiling=meta['max_busy'],
                       n_lo=r['n_lo'], n_hi=r['n_hi'],
                       sync=k[len('jax_'):], cpus=r['cpus'],
                       **placement_fields)
            validate(row)
            out.append(row)
        if 'fortran' in r:
            f = r['fortran']
            row = _base(meta, 'run_mpi_scaling', snapshot, tenancy,
                        backfilled_from)
            row.update(backend='fortran', ranks=r['ranks'],
                       ms_per_step=f['ms_per_step'],
                       parallelism='mpi',
                       platform='cpu', devices=None,
                       platform_evidence=FORTRAN_CPU_EVIDENCE,
                       rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
                       effective_cores=None, threads_per_rank=None,
                       busy_ceiling=meta['max_busy'],
                       n_lo=r['n_lo'], n_hi=r['n_hi'], cpus=r['cpus'],
                       **placement_fields)
            validate(row)
            out.append(row)
    return out


_E2E_PARALLELISM_BY_BACKEND = {
    # run_e2e.py:475-478 -- these two cells launch `mpirun -np ranks`
    # (matrix.FORTRAN_RANKS / matrix.PY_MPI_RANKS), one process per rank.
    'fortran': 'mpi',
    'python-jax-mpi': 'mpi',
    # run_e2e.py:483-490 -- one unpinned process, ranks fixed at 1. The cell
    # selects no cpus and sets no thread-count env, so `ranks`=1 counts THE
    # PROCESS and is not a statement about cores used; 'serial' records
    # exactly that and must not be read as "one core".
    'python-numpy': 'serial',
    'python-jax': 'serial',
}


def _E2E_PARALLELISM(backend):
    """The discriminator for one e2e cell. Raises on an unknown backend: a new
    column must declare how it launches before its wall clock is recorded."""
    try:
        return _E2E_PARALLELISM_BY_BACKEND[backend]
    except KeyError:
        raise ValueError(
            'e2e backend %r has no declared parallelism (known: %s) -- '
            'refusing to guess what its `ranks` counts'
            % (backend, sorted(_E2E_PARALLELISM_BY_BACKEND))) from None


def rows_from_e2e_results(meta, snapshot, tenancy):
    """Ledger rows from one e2e sweep (testsys/e2e/run_e2e.py): metric
    'cell-wall-clock', one row per PASSING cell.

    wall_s is the cell's whole wall time -- create.newcase, case.setup, the
    solve, canonicalisation, comparison, and (for jax cells) XLA compile --
    which is exactly why this is a DIFFERENT metric from
    'per-step-by-difference' and carries no ms_per_step at all (validate
    enforces the explicit null). Failed cells produce no row: their wall time
    measures the failure path, not the solver, and a red cell must not leave
    a green-looking timing behind."""
    out = []
    for c in meta['cells']:
        if not c['ok']:
            continue
        row = dict(ts_utc=utc_now(), snapshot_date_local=meta['date'],
                   sha=meta['sha'], host=meta['host'], tool='run_e2e',
                   case=c['case'], backend=c['backend'],
                   metric='cell-wall-clock', snapshot=snapshot,
                   tenancy_busy=tenancy['busy'],
                   tenancy_total=tenancy['total'],
                   ranks=c['ranks'], ms_per_step=None,
                   parallelism=_E2E_PARALLELISM(c['backend']),
                   wall_s=c['seconds'],
                   rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
                   effective_cores=None, threads_per_rank=None,
                   busy_ceiling=meta['tenancy_ceiling'], n_lo=None, n_hi=None,
                   platform=c['platform'],
                   platform_evidence=c['platform_evidence'],
                   devices=None, verdict='SUCCESS', selection=meta['label'])
        validate(row)
        out.append(row)
    return out


def capture_e2e_cells(meta, root=ROOT):
    """Write one e2e sweep's per-cell wall clocks as a dated snapshot file
    plus ledger rows. Returns (n_rows_appended, snapshot_relpath). RAISES on
    any problem -- the sweep itself must call capture_e2e_cells_or_warn."""
    name = 'e2e_cells_%s_%d.json' % (time.strftime('%Y-%m-%d_%H%M%S'),
                                     os.getpid())
    snap_dir = os.path.join(root, SNAPSHOT_DIR_RELPATH)
    os.makedirs(snap_dir, exist_ok=True)
    snap_path = os.path.join(snap_dir, name)
    if os.path.exists(snap_path):
        raise RuntimeError('snapshot %s already exists -- refusing to '
                           'overwrite a dated file (rule 19)' % snap_path)
    with open(snap_path, 'w') as f:
        json.dump(meta, f, indent=1)
    rel = os.path.join(SNAPSHOT_DIR_RELPATH, name).replace(os.sep, '/')
    tenancy = box_tenancy(meta['tenancy_ceiling'])
    rows = rows_from_e2e_results(meta, rel, tenancy)
    return append_rows(rows, os.path.join(root, LEDGER_RELPATH)), rel


def capture_e2e_cells_or_warn(meta, root=ROOT):
    """capture_e2e_cells, degraded to a LOUD WARNING on any error. The sweep
    is a parity gate; a ledger problem must never change a cell's verdict or
    the sweep's exit code, in either direction -- the timing data point is
    lost, not the gate. SystemExit is included deliberately: box_tenancy
    raises it where numactl is absent (CI runners record no tenancy and
    therefore no row). Returns the row count, or None on failure."""
    try:
        n, rel = capture_e2e_cells(meta, root)
        print('perf ledger: %d cell-wall-clock row(s) appended to %s '
              '(snapshot %s)' % (n, LEDGER_RELPATH, rel))
        return n
    except (Exception, SystemExit) as exc:          # noqa: BLE001
        print('WARNING: perf-ledger capture failed (%s: %s) -- the sweep '
              'verdict above is unaffected; this run\'s timing data point is '
              'lost, not the gate.' % (type(exc).__name__, exc))
        return None


def _rows_from_snapshot_file(snapshot_path):
    """(rows, rel) from ONE committed dated snapshot file, tool detected from
    the snapshot's own schema. The snapshot file is the provenance; prose
    (session logs, board rows) is never a source."""
    rel = os.path.relpath(os.path.abspath(snapshot_path), ROOT)
    if not rel.replace(os.sep, '/').startswith('docs/perf_snapshots/'):
        raise SystemExit('FAIL: source %s is not under docs/perf_snapshots/ '
                         '-- only committed dated snapshots have the '
                         'provenance a ledger row needs.' % snapshot_path)
    meta = json.load(open(snapshot_path))
    if 'max_busy' in meta and 'busy_ceiling' not in meta:
        return rows_from_mpi_scaling_snapshot(meta, rel, None,
                                              backfilled_from=rel), rel
    if 'busy_ceiling' in meta and 'max_busy' not in meta:
        return rows_from_scaling_snapshot(meta, rel, None,
                                          backfilled_from=rel), rel
    raise SystemExit('FAIL: cannot tell which tool wrote %s (looked for '
                     'exactly one of max_busy/busy_ceiling); refusing to '
                     'guess a schema.' % snapshot_path)


def backfill(snapshot_path, path=LEDGER):
    """Seed ledger rows from ONE committed dated snapshot file. Returns the
    number of rows appended (0 is a valid answer: a snapshot whose configs
    were all SKIPPED holds no measurement)."""
    rows, _rel = _rows_from_snapshot_file(snapshot_path)
    return append_rows(rows, path)


def reissue(snapshot_path, lines, path=LEDGER):
    """Append CORRECTED identities for existing ledger lines that predate the
    platform field, re-derived from the SAME snapshot. Append-only
    supersession: the old lines stay byte-identical; each new row carries
    supersedes={line, ts_utc, reason} and the newest row for a (snapshot,
    ranks, sync) wins. Refuses unless every superseded line matches the
    re-derived row on snapshot/case/backend/ranks/sync AND on ms_per_step --
    a reissue changes a row's identity, never its number."""
    rows, rel = _rows_from_snapshot_file(snapshot_path)
    if len(rows) != len(lines):
        raise SystemExit('FAIL: %s re-derives %d row(s) but %d superseded '
                         'line(s) were named -- the mapping must be exact, '
                         'in converter order.' % (snapshot_path, len(rows),
                                                  len(lines)))
    existing = open(path).read().splitlines()
    for row, ln in zip(rows, lines):
        if not 1 <= ln <= len(existing):
            raise SystemExit('FAIL: line %d is not in the ledger (%d lines).'
                             % (ln, len(existing)))
        old = json.loads(existing[ln - 1])
        for k in ('snapshot', 'case', 'backend', 'ranks', 'sync',
                  'ms_per_step'):
            if old.get(k) != row.get(k):
                raise SystemExit(
                    'FAIL: ledger line %d has %s=%r but the re-derived row '
                    'has %r -- refusing to mark supersession across two '
                    'different measurements.' % (ln, k, old.get(k),
                                                 row.get(k)))
        if 'platform' in old:
            raise SystemExit('FAIL: ledger line %d already carries a '
                             'platform (%r) -- reissue exists for legacy '
                             'rows only.' % (ln, old['platform']))
        row['supersedes'] = dict(
            line=ln, ts_utc=old['ts_utc'],
            reason='original row predates the platform field; identity '
                   're-derived from the snapshot\'s own device evidence, '
                   'values unchanged')
    n = append_rows(rows, path)
    for row in rows:
        print('  line %d (%s ranks=%s sync=%s) -> platform=%r devices=%r'
              % (row['supersedes']['line'], row['case'], row['ranks'],
                 row.get('sync'), row['platform'], row['devices']))
    return n


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
    if len(argv) == 4 and argv[1] == 'reissue':
        lines = [int(x) for x in argv[3].split(',') if x]
        n = reissue(argv[2], lines)
        print('reissue done: %d corrected row(s) appended to %s; the '
              'superseded lines stay byte-identical (append-only)'
              % (n, LEDGER_RELPATH))
        return 0
    print('usage: python3 testsys/perf/ledger.py backfill '
          'docs/perf_snapshots/<file>.json [...]\n'
          '       python3 testsys/perf/ledger.py reissue '
          'docs/perf_snapshots/<file>.json <line,line,...>')
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
