#! /usr/bin/env python3
"""Schema + validation for the per-rank runtime profile every backend now
emits (owner: wei-lin, item 2026-09-23; capture emitter is mira-volkov's
parallel mission -- this module only reads and judges what it wrote).

`<run_dir>/profile.rank<r>.json`, one file per MPI rank (`nranks`=1 for a
serial cell):
    {"schema": "eqdyna-profile/1",
     "backend": "fortran|python-numpy|python-jax|python-jax-mpi",
     "rank": int, "nranks": int, "host": str, "pid": int,
     "cpus_allowed": [int...], "numa_nodes": [int...],
     "nsteps": int, "sampling_every": int,
     "buckets_s": {"setup","element","fault","exchange","wait","io"},
     "loop_s": f, "total_s": f, "unaccounted_s": f}

Fields may be ADDED to this schema; an existing field is never renamed --
`validate()` checks for the REQUIRED set and tolerates (and, for buckets_s,
individually validates) extras.

THE SUM CHECK, AND WHY IT IS THE ONE THAT MATTERS.
`total_s` and `loop_s` are measured by their OWN timers; no bucket is a
remainder computed FROM total_s. That independence is exactly what lets this
validator ask a real question: do the independently-timed buckets, plus the
profiler's own `unaccounted_s`, actually add up to the independently-timed
`total_s`? A "yes" is not guaranteed by construction -- it is the property a
correct profiler has and a broken one does not.

Background this guards against (CLAUDE.md, "things that will bite you"):
`writeCompTime` initialised to 0 and read from nowhere let a 511x error in
`compTimeInSeconds(2)` survive completely undetected, because nothing ever
summed that bucket against anything. Two checks close that gap:

  1. `unaccounted_s` must equal `total_s - sum(buckets_s.values())` TO FLOAT
     PRECISION (FLOAT_TOL, below). The narrow check: the emitter's own
     bookkeeping of "what fell outside every bucket" must agree with what a
     remainder computed from its OTHER two independent numbers says it
     should be. A profiler whose `unaccounted_s` field is stale, mis-summed,
     or copy-pasted from the wrong run fails this immediately, even if the
     numbers involved are individually plausible.
  2. `|total_s - sum(buckets_s.values())| <= max(SUM_TOLERANCE * total_s,
     SUM_FLOOR_S)`. The
     wide check, and the one a x511 bucket error cannot survive: it bounds
     how much of `total_s` is allowed to sit in that gap at all.

SUM_TOLERANCE = 0.05 (5%). Justification: buckets_s times independently-timed
sections of a step loop that ALSO contains code between/around those
sections -- Python/MPI call overhead, GC, OS scheduling jitter, clock
resolution -- none of which is itself bucketed. That gap is real and is not a
bug; 5% sits comfortably above ordinary jitter (this repo's own run-to-run
noise floor is measured elsewhere at 1-2%, e.g. run_scaling.py's "reproducible
to under 1%") while being far too tight for a gross accounting bug to hide in:
a single bucket wrong by 511x (the compTimeInSeconds(2) incident), or simply
omitted, moves the gap by orders of magnitude past 5%, not by a fraction of
it. It is a named module constant, not a magic number, so a future
measurement of the real jitter floor can widen it deliberately, in one place,
with the justification updated here.

SUM_FLOOR_S = 2.0 s, added 2026-09-24 because the 5% alone was a FLAKY gate
on short runs. The gap is mostly a FIXED per-run overhead (import, output,
the code around the timed sections), not a proportional one. Measured over
409 committed `docs/run_profiles.jsonl` rows: the largest absolute gap is
2.73 s, on a 188 s fortran tpv36 run (1.4%). But test.tpv8 x python-jax, the
shortest cell (~15 s), reads up to 4.8% locally, and on the 2-core CI runner
it read 0.80 s of 15.40 s = 5.2% and failed master (run 36082084976) with physics
well inside its bound. A 2 s floor admits that fixed overhead on a short run.
It changes nothing for runs of 40 s and longer, where 5% already exceeds 2 s.
A x511 bucket still moves the gap by orders of magnitude past either limit.
"""
import json
import math
import os
import re

SCHEMA_ID = 'eqdyna-profile/1'
BACKENDS = ('fortran', 'python-numpy', 'python-jax', 'python-jax-mpi')
BUCKET_KEYS = ('setup', 'element', 'fault', 'exchange', 'wait', 'io')

REQUIRED = ('schema', 'backend', 'rank', 'nranks', 'host', 'pid',
            'cpus_allowed', 'numa_nodes', 'nsteps', 'sampling_every',
            'buckets_s', 'loop_s', 'total_s', 'unaccounted_s')

# See module docstring for the justification of both constants.
SUM_TOLERANCE = 0.05
SUM_FLOOR_S = 2.0
FLOAT_TOL = 1e-9

_RANK_FILE_RE = re.compile(r'^profile\.rank(\d+)\.json$')


def check_nonneg_int(context, name, v, min_=0):
    """Shared by profile_record.py (rule 1: one copy of a check every row
    this schema's data feeds gets held to)."""
    if not (isinstance(v, int) and not isinstance(v, bool) and v >= min_):
        raise ValueError('%s: %s %r must be an int >= %d' % (context, name, v, min_))


def check_int_list(context, name, v, unique=True, allow_empty=False):
    """A list of non-negative ints. `unique=False` for a field like
    `ranks_per_node` where the SAME count legitimately repeats across hosts;
    `unique=True` (default) for a set-like field such as `cpus_allowed` or
    `numa_nodes`, where a repeated entry means the same cpu/node was listed
    twice and is itself a bug worth catching."""
    if not isinstance(v, list) or (not v and not allow_empty):
        raise ValueError('%s: %s %r must be a non-empty list' % (context, name, v))
    if unique and len(set(v)) != len(v):
        raise ValueError('%s: %s %r contains duplicate entries' % (context, name, v))
    for x in v:
        if not (isinstance(x, int) and not isinstance(x, bool) and x >= 0):
            raise ValueError('%s: %s entry %r must be a non-negative int'
                             % (context, name, x))


def check_buckets(context, total_s, buckets_s, unaccounted_s,
                  bucket_keys=BUCKET_KEYS):
    """The sum check described in the module docstring. Shared by
    profile_record.py so a row re-validates the SAME arithmetic the emitted
    profile file was already held to, rather than trusting the file's own
    say-so a second time down the pipeline."""
    if not isinstance(buckets_s, dict):
        raise ValueError('%s: buckets_s %r must be a dict' % (context, buckets_s))
    missing = [k for k in bucket_keys if k not in buckets_s]
    if missing:
        raise ValueError('%s: buckets_s missing required bucket(s) %s'
                         % (context, missing))
    for k, v in buckets_s.items():
        if not (isinstance(v, (int, float)) and not isinstance(v, bool) and v >= 0):
            raise ValueError('%s: buckets_s[%r] = %r must be a non-negative number'
                             % (context, k, v))
    for name, v in (('total_s', total_s), ('unaccounted_s', unaccounted_s)):
        if not (isinstance(v, (int, float)) and not isinstance(v, bool)):
            raise ValueError('%s: %s %r must be a number' % (context, name, v))
    if total_s <= 0:
        raise ValueError('%s: total_s %r must be > 0' % (context, total_s))

    bucket_sum = sum(buckets_s.values())
    gap = total_s - bucket_sum
    if not math.isclose(unaccounted_s, gap, rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL):
        raise ValueError(
            '%s: unaccounted_s %.9g does not equal total_s - sum(buckets_s) = '
            '%.9g (to float precision, tol=%.1e) -- unaccounted_s must be the '
            'SAME remainder this check computes independently, or an emitter '
            'bug in the unaccounted timer goes undetected exactly like the '
            '511x compTimeInSeconds(2) error this check exists to catch'
            % (context, unaccounted_s, gap, FLOAT_TOL))
    allowed = max(SUM_TOLERANCE * total_s, SUM_FLOOR_S)
    if abs(gap) > allowed:
        raise ValueError(
            '%s: |total_s - sum(buckets_s)| = %.6g exceeds max(%.0f%% of total_s, '
            '%.1f s) = %.6g (total_s %.6g) -- buckets_s + unaccounted_s must '
            'account for essentially '
            'all of total_s; a bucket inflated or deflated by orders of '
            'magnitude (the 511x compTimeInSeconds(2) incident) is exactly '
            'what this catches' % (context, abs(gap), SUM_TOLERANCE * 100,
                                   SUM_FLOOR_S, allowed, total_s))


def validate(path_or_dict):
    """Raise ValueError on the first uninterpretable thing about one
    profile.rank<r>.json row. Returns the validated dict (loaded from disk if
    a path was given)."""
    if isinstance(path_or_dict, dict):
        row, context = path_or_dict, '<dict>'
    elif isinstance(path_or_dict, (str, os.PathLike)):
        path = str(path_or_dict)
        try:
            with open(path) as f:
                row = json.load(f)
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError('cannot read/parse profile file %s: %s' % (path, exc)) from exc
        context = path
    else:
        raise ValueError('validate() expects a path or dict, got %r' % type(path_or_dict))

    missing = [k for k in REQUIRED if k not in row]
    if missing:
        raise ValueError('%s: missing required field(s) %s' % (context, missing))

    if row['schema'] != SCHEMA_ID:
        raise ValueError('%s: schema %r != %r' % (context, row['schema'], SCHEMA_ID))
    if row['backend'] not in BACKENDS:
        raise ValueError('%s: backend %r not in %s' % (context, row['backend'], BACKENDS))

    check_nonneg_int(context, 'rank', row['rank'], 0)
    check_nonneg_int(context, 'nranks', row['nranks'], 1)
    if row['rank'] >= row['nranks']:
        raise ValueError('%s: rank %r must be < nranks %r'
                         % (context, row['rank'], row['nranks']))
    if not (isinstance(row['host'], str) and row['host'].strip()):
        raise ValueError('%s: host %r must be a non-empty string' % (context, row['host']))
    check_nonneg_int(context, 'pid', row['pid'], 1)

    check_int_list(context, 'cpus_allowed', row['cpus_allowed'])
    check_int_list(context, 'numa_nodes', row['numa_nodes'])

    check_nonneg_int(context, 'nsteps', row['nsteps'], 1)
    check_nonneg_int(context, 'sampling_every', row['sampling_every'], 1)

    if not (isinstance(row['loop_s'], (int, float)) and not isinstance(row['loop_s'], bool)
            and row['loop_s'] >= 0):
        raise ValueError('%s: loop_s %r must be a number >= 0' % (context, row['loop_s']))

    check_buckets(context, row['total_s'], row['buckets_s'], row['unaccounted_s'])

    return row


def validate_run_dir(run_dir):
    """Validate every profile.rank<r>.json under `run_dir` and their mutual
    consistency: same nranks, same backend, and rank files present for EXACTLY
    0..nranks-1 (missing -> ValueError naming which rank; unexpected extra ->
    ValueError naming which). Returns the validated rows, sorted by rank."""
    if not os.path.isdir(run_dir):
        raise ValueError('run_dir %r is not a directory' % (run_dir,))
    found = {}
    for name in sorted(os.listdir(run_dir)):
        m = _RANK_FILE_RE.match(name)
        if not m:
            continue
        found[int(m.group(1))] = os.path.join(run_dir, name)
    if not found:
        raise ValueError('no profile.rank<N>.json files found under %r' % (run_dir,))

    rows = {}
    for r, path in found.items():
        row = validate(path)
        if row['rank'] != r:
            raise ValueError('%s: filename declares rank %d but the row carries rank=%r'
                             % (path, r, row['rank']))
        rows[r] = row

    nranks_seen = {row['nranks'] for row in rows.values()}
    if len(nranks_seen) != 1:
        raise ValueError('run_dir %r: inconsistent nranks across rank files: %s'
                         % (run_dir, sorted(nranks_seen)))
    nranks = nranks_seen.pop()
    expected, present = set(range(nranks)), set(rows)
    missing = sorted(expected - present)
    extra = sorted(present - expected)
    if missing:
        raise ValueError(
            'run_dir %r: missing rank file(s) for rank(s) %s (nranks=%d, '
            'found ranks %s)' % (run_dir, missing, nranks, sorted(present)))
    if extra:
        raise ValueError(
            'run_dir %r: unexpected rank file(s) for rank(s) %s outside '
            '0..nranks-1=%d' % (run_dir, extra, nranks - 1))

    backends_seen = {row['backend'] for row in rows.values()}
    if len(backends_seen) != 1:
        raise ValueError('run_dir %r: inconsistent backend across rank files: %s'
                         % (run_dir, sorted(backends_seen)))

    return [rows[r] for r in range(nranks)]
