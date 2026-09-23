"""profile_emit.py -- no Fortran counterpart (like backend.py/__main__.py).

Writes the ALWAYS-ON per-rank `<run_dir>/profile.rank<r>.json`
(schema `eqdyna-profile/1`, testsys/profile_schema.py) for every Python
backend: python-numpy, python-jax (serial, nranks=1) and python-jax-mpi
(driver.run_mpi, real MPI). Mirrors library_output.f90's `output_profile`
(src/fortran/library_output.f90 #10) -- same schema, same six buckets, same
EQDYNA_PROFILE=0 off switch, same default-ON.

BUCKET DEPARTURE FROM FORTRAN, documented per rule 23 (Fortran is the
reference and DEFINES the buckets; the port may depart in structure only for
measured perf, documented): Fortran keeps `element` and `fault` disjoint
because assembleGlobalKU/calcHourglassResist and faulting.f90 are separate
subroutine calls with their own MPI_WTIME() brackets (driver.f90). Every
Python backend (numpy AND jax) fuses velDispUpdate + both element kernels +
faulting into ONE step function (driver.make_step / driver.run_mpi's
a_body/b_body) -- see driver.py's module docstring. Splitting `element` from
`fault` back apart needs a NEW timer INSIDE that fused step:
  - numpy executes that step eagerly and synchronously (no queue, no device
    to drain), so a `time.perf_counter()` pair around the `FLT.faulting`
    call inside `driver.make_step_parts`'s `part_b` costs no new sync and
    changes no arithmetic. Landed 2026-09-23: `driver.run` accumulates it
    into a `fault_timer` dict, built only `if not B.is_jax(xp)`, and returns
    it as `out['fault_s']`; `eqdyna3d.run_case` subtracts it back out of
    `Profile['solve']` to get `element`. **python-numpy now matches
    Fortran's `fault` boundary.**
  - jax: still impossible without a NEW `block_until_ready` between the
    element kernels and faulting, which the mission's "no new sync" rule
    forbids on the default path -- and a timer placed anyway would measure
    the ONE-TIME trace, not the per-step cost (wrong, not just imprecise).
    `driver.make_step_parts` refuses to build the timer at all when
    `B.is_jax(xp)`, so `fault_s` stays 0.0 there by construction.
So `fault_s` is folded into `element_s` and always reported as 0.0 in the
`fault` bucket on **python-jax and python-jax-mpi only**. `element` means
"element + fault, fused" for those two backends -- read it that way, not as
a like-for-like column against Fortran's `element`. On python-numpy,
`element` and `fault` are disjoint and match Fortran's boundary.

`exchange`/`wait` are genuinely 0.0 for the two serial backends (no MPI to
speak of) and are REAL measurements (not folded, not skipped) for
python-jax-mpi -- see driver.run_mpi's now-unconditional t_compute/t_mpi
accumulation, added by this landing at the EXISTING sync point
(`jax.block_until_ready(hv)`, driver.py, already unconditional before this
change) rather than a new one. `wait` stays 0.0 on python-jax-mpi's default
path for a different, real reason: the production exchange has no collective
barrier (nearest-neighbour Sendrecv only, MPI4NodalQuant.exchange) -- there is
nothing to measure without ADDING one, which the mission explicitly forbids
on the default path. Any per-step blocking-on-a-slow-neighbour is inside
Sendrecv itself and is therefore already inside `exchange`, exactly mirroring
how Fortran's own `wait` is a sub-timer NESTED inside (and subtracted back
out of) `exchange` -- the same nesting relationship, minus the sub-timer,
because there is no separate call here to subtract.
"""
import json
import os
import re

SCHEMA_ID = 'eqdyna-profile/1'
BUCKET_KEYS = ('setup', 'element', 'fault', 'exchange', 'wait', 'io')
OFF_ENV = 'EQDYNA_PROFILE'

_RANGE_RE = re.compile(r'(\d+)(?:-(\d+))?')


def _parse_range_list(text):
    """"48-51,60,62-63" -> [48,49,50,51,60,62,63]. Same format Fortran's
    output_profile parses from /proc and /sys (kernel cpu-range-list)."""
    out = []
    for tok in text.strip().split(','):
        tok = tok.strip()
        if not tok:
            continue
        m = _RANGE_RE.match(tok)
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) else lo
        out.extend(range(lo, hi + 1))
    return out


def cpus_allowed():
    """This process's cpu affinity mask, exactly as Fortran's
    output_profile/readCpusAllowed reads /proc/self/status -- but Python
    already has this as a syscall wrapper, so no /proc parse is needed."""
    return sorted(os.sched_getaffinity(0))


def numa_nodes(cpu_list):
    """NUMA node(s) backing `cpu_list`, from
    /sys/devices/system/node/node<k>/cpulist -- same source and overlap test
    as Fortran's computeNumaNodes (library_output.f90)."""
    cpu_set = set(cpu_list)
    nodes = []
    base = '/sys/devices/system/node'
    if not os.path.isdir(base):
        return nodes
    for name in sorted(os.listdir(base)):
        if not name.startswith('node'):
            continue
        try:
            k = int(name[4:])
        except ValueError:
            continue
        path = os.path.join(base, name, 'cpulist')
        try:
            with open(path) as f:
                node_cpus = set(_parse_range_list(f.read()))
        except OSError:
            continue
        if cpu_set & node_cpus:
            nodes.append(k)
    return sorted(nodes)


def enabled():
    """EQDYNA_PROFILE=0 is the OFF switch (A/B instrument only); default ON,
    same contract as the Fortran side's output_profile."""
    return os.environ.get(OFF_ENV, '') != '0'


def build_row(backend, rank, nranks, nsteps, buckets_s, loop_s, total_s,
              sampling_every=1):
    """Build one validated-shape profile row. Raises (does not default) on a
    missing bucket key -- callers must supply all six, 0.0 where genuinely
    zero or folded elsewhere (documented at the call site and in this
    module's docstring), never omitted."""
    missing = [k for k in BUCKET_KEYS if k not in buckets_s]
    if missing:
        raise ValueError('profile_emit.build_row: buckets_s missing %s -- '
                         'every one of %s must be supplied explicitly'
                         % (missing, BUCKET_KEYS))
    cpus = cpus_allowed()
    unaccounted = total_s - sum(buckets_s[k] for k in BUCKET_KEYS)
    return {
        'schema': SCHEMA_ID,
        'backend': backend,
        'rank': int(rank),
        'nranks': int(nranks),
        'host': os.uname().nodename,
        'pid': os.getpid(),
        'cpus_allowed': cpus,
        'numa_nodes': numa_nodes(cpus),
        'nsteps': int(nsteps),
        'sampling_every': int(sampling_every),
        'buckets_s': {k: float(buckets_s[k]) for k in BUCKET_KEYS},
        'loop_s': float(loop_s),
        'total_s': float(total_s),
        'unaccounted_s': float(unaccounted),
    }


def write_profile(run_dir, backend, rank, nranks, nsteps, buckets_s, loop_s,
                  total_s, sampling_every=1):
    """Write `<run_dir>/profile.rank<rank>.json`. No-op (returns None) when
    EQDYNA_PROFILE=0, identical off-switch contract to the Fortran side."""
    if not enabled():
        return None
    row = build_row(backend, rank, nranks, nsteps, buckets_s, loop_s, total_s,
                    sampling_every=sampling_every)
    path = os.path.join(run_dir, 'profile.rank%d.json' % rank)
    with open(path, 'w') as f:
        json.dump(row, f, indent=None, separators=(',', ':'))
    return path
