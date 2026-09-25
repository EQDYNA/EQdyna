#! /usr/bin/env python3
"""Strong-scaling measurement for the jax backend's EXPLICIT domain
decomposition (backend.run_time_loop_sharded), report-only, never a gate.
pathway_forward item 43's settling experiment.

WHAT IT MEASURES, AND WHY THREE MODES INSTEAD OF ONE. The question item 43
left open is not "is sharded jax faster" but "WHAT bounds it". A single
ms/step column cannot answer that, so every device count is measured under
three configurations whose DIFFERENCES attribute the cost:

  element            the real thing: element arrays cut across N devices,
                     one psum of the whole nodal force array per step at
                     driver.f90:27's position. Correct physics.
  element-nosync     identical, collective removed (EQDYNA_SHARD_SYNC=none).
                     WRONG ANSWER, valid timing. element minus this IS the
                     collective's cost, measured rather than inferred.
  element+nodal      element arrays AND the node-axis arrays of
                     velDispUpdate cut across devices (EQDYNA_SHARD_MODE).
                     WRONG ANSWER (no halo: each device's nodal arrays are
                     stale where another device owns the node), valid timing.
                     This is the CEILING a correct halo implementation could
                     reach, measured before paying for the halo -- it differs
                     from a correct one only by the boundary exchange, which
                     is a psum over |boundary| words instead of |NEQ|.

So: element+nodal tells you whether the halo is worth building, and
element minus element-nosync tells you how much of the naive version's
plateau is the all-reduce as opposed to the replicated nodal work.

PER-STEP BY DIFFERENCE, two fresh processes per point (run_scaling.py's
`per_step_py` technique and the same reason: jax caches compiled code
in-process, so a second call in one interpreter pays no compile and would
cancel the wrong term). XLA compile is ~15% of a 114-step jax run; a total
wall clock would make the answer a property of the metric. `88f5227` fixed
exactly this bias on the Fortran side -- do not reintroduce it here.

PLACEMENT AND BUSY-CHECKING are not re-implemented: `numa_topology`,
`cpu_busy_fractions`, `require_idle` come from run_numa_scaling.py and
`compact_cpus`/`spread_cpus`/`free_node_map`/`select_cpus`/`numactl_prefix`/
`nodes_of` from run_scaling.py, so a fix to either is a fix here too. This
box carries ~10-25 foreign single-core tenants smeared across all 8 NUMA
nodes, so a point whose cpus are not free is SKIPPED and recorded, never
silently included. `--wait-for-room S` retries the free-cpu probe for up to
S seconds before giving up on a point, because on this box occupancy moves on
a minutes timescale and a sweep that gives up instantly measures nothing.

FORTRAN IS NOT RE-MEASURED HERE. Its 1/2/4/8/16-core numbers on this case and
box are recorded in pathway_forward item 43 (931.65/466.64/237.62/117.07/
64.74 ms/step, by the same n_lo/n_hi difference); `--fortran` re-runs them
through run_scaling.py's own path if the comparison needs refreshing rather
than a second implementation of it.

Usage:
    python3 testsys/perf/run_shard_scaling.py [--case test.tpv104]
        [--devices 1,2,4,8,16] [--modes element,element-nosync,element+nodal]
        [--n-lo 20] [--n-hi 60] [--busy-ceiling 0.2] [--wait-for-room 600]
        [--repeats 1] [--i-know-the-box-is-busy]
"""
import argparse
import json
import os
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(TESTSYS))
PYTHON_PKG = os.path.join(ROOT, 'src', 'python')
OUT = os.path.join(ROOT, 'docs', 'perf_snapshots',
                   'shard_scaling_%s_%s.json' % (time.strftime('%Y-%m-%d_%H%M%S'),
                                          os.environ.get('EQDYNA_SNAPSHOT_TAG', 'tpv104')))
# PROJECT_RULES rule 19: dated and immutable, never a shared *_last.* file.
# A path whose name says "last" cites whatever ran most recently, so it is
# not evidence a board row can point at; two tools writing one such file
# silently discard each other. Item 91a: the timestamp must include time of
# day (not just the date), or two runs on the same day overwrite each other
# -- exactly the failure mode this comment already warned about.

sys.path.insert(0, TESTSYS)
import run_numa_scaling as numa      # noqa: E402
import run_scaling as rs             # noqa: E402
import ledger                        # noqa: E402

MODES = {
    'element':        dict(EQDYNA_SHARD_MODE='element', EQDYNA_SHARD_SYNC='psum'),
    'element-nosync': dict(EQDYNA_SHARD_MODE='element', EQDYNA_SHARD_SYNC='none'),
    'element+nodal':  dict(EQDYNA_SHARD_MODE='element+nodal', EQDYNA_SHARD_SYNC='psum'),
}
# Recorded Fortran MPI baseline, same case/box/metric (pathway item 43,
# 2026-09-19, sha 737358f). Quoted for the table, NOT re-derived here.
FORTRAN_MS = {1: 931.65, 2: 466.64, 4: 237.62, 8: 117.07, 16: 64.74}

CHILD = r"""
import os, sys, time
sys.path.insert(0, %(pkg)r)
from eqdyna import eqdyna3d
t0 = time.time()
eqdyna3d.run_case(%(case)r, nsteps=%(n)d, verbose=False, backend='jax')
print('WALL', time.time() - t0)
import jax
print('DEVICES', len(jax.devices()))
print('THREADS', len(os.listdir('/proc/%%d/task' %% os.getpid())))
"""


def time_one(case_dir, nsteps, cpus, node_map, ndev, mode):
    """One fresh, pinned process. Returns (wall, ndev_seen, nthreads) or None.

    Does NOT set XLA_FLAGS/OMP_NUM_THREADS/OPENBLAS_NUM_THREADS: XLA auto-sizes
    its intra-op pool from the process's cpu affinity, which numactl has
    already set (run_scaling.py's measured finding). The device count travels
    as EQDYNA_JAX_DEVICES and backend.py turns it into
    --xla_force_host_platform_device_count before jax initialises -- setting
    XLA_FLAGS here as well would be a second authority for it and backend.py
    raises on that deliberately."""
    script = CHILD % dict(pkg=PYTHON_PKG, case=case_dir, n=nsteps)
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    env['EQDYNA_JAX_DEVICES'] = str(ndev)
    env.update(MODES[mode])
    for k in ('XLA_FLAGS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        env.pop(k, None)
    cmd = rs.numactl_prefix(cpus, node_map) + [sys.executable, '-c', script]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1200:]); print(r.stderr[-1500:])
        return None
    got = {}
    for line in r.stdout.splitlines():
        p = line.split()
        if p and p[0] in ('WALL', 'DEVICES', 'THREADS'):
            got[p[0]] = float(p[1])
    if len(got) != 3:
        raise RuntimeError('child did not report WALL/DEVICES/THREADS: %r'
                           % r.stdout[-500:])
    if int(got['DEVICES']) != ndev:
        raise RuntimeError(
            'child ran on %d jax devices, %d were requested -- a mislabelled '
            'scaling point is worse than a missing one'
            % (int(got['DEVICES']), ndev))
    return got['WALL'], int(got['DEVICES']), int(got['THREADS'])


# Item 91g's baseline bookkeeping lives in run_scaling.py (one copy, shared
# with that tool's two loops); this name is kept for callers and the guard.
record_point = rs.record_point


def per_step(case_dir, cpus, node_map, ndev, mode, n_lo, n_hi):
    lo = time_one(case_dir, n_lo, cpus, node_map, ndev, mode)
    hi = time_one(case_dir, n_hi, cpus, node_map, ndev, mode)
    if lo is None or hi is None:
        return None
    ps, fixed = numa.per_step_and_fixed(lo[0], hi[0], n_lo, n_hi)
    return dict(ms_per_step=ps * 1e3, fixed_s=fixed,
                wall_lo_s=lo[0], wall_hi_s=hi[0], threads=hi[2])


def room(nodes, k, ceiling, override, wait_s):
    """(cpus, busy) for k free cpus, compact over currently-free nodes, or
    (None, None) after `wait_s` seconds of retrying. Occupancy on this box
    moves over minutes; a single probe is a coin flip, not a verdict."""
    deadline = time.time() + wait_s
    while True:
        free = rs.free_node_map(nodes, ceiling, override)
        cpus = rs.select_cpus(free, k, 'compact')
        if cpus is not None:
            busy = numa.require_idle(cpus, ceiling, override)
            if busy is not None:
                return cpus, busy
        if time.time() >= deadline:
            return None, None
        time.sleep(20)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--devices', default='1,2,4,8,16')
    ap.add_argument('--modes', default='element,element-nosync,element+nodal')
    ap.add_argument('--n-lo', type=int, default=20)
    ap.add_argument('--n-hi', type=int, default=60)
    ap.add_argument('--repeats', type=int, default=1)
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    ap.add_argument('--wait-for-room', type=float, default=600.0)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    a = ap.parse_args()
    devices = [int(x) for x in a.devices.split(',') if x]
    modes = [m for m in a.modes.split(',') if m]
    for m in modes:
        if m not in MODES:
            raise SystemExit('FAIL: unknown mode %r (have %r)' % (m, sorted(MODES)))
    for n in devices:
        if n > 16:
            raise SystemExit('FAIL: %d devices requested -- this experiment is '
                             'capped at 16 (owner constraint).' % n)

    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology; this '
                         'experiment is about placement and cannot run blind.')
    rs.CASE = a.case
    case_dir = rs.build_py_case(a.case)
    sha = rs.sh('git -C %s rev-parse --short HEAD' % ROOT).stdout.strip()
    rows, skipped = [], []
    print('case %s  sha %s  host %s  n_lo/n_hi %d/%d'
          % (a.case, sha, os.uname().nodename, a.n_lo, a.n_hi), flush=True)

    base = {}
    base_n = {}
    # Item 91g: the speedup column baselines on the first NON-SKIPPED n for
    # each mode. If devices[0] (the intended n=1 baseline) is skipped or
    # fails, base_n[mode] silently becomes whatever n measured first, and a
    # speedup of 1.00x at that n would be misread as "no scaling" rather than
    # "this IS the baseline". So the baseline point is recorded explicitly
    # (baseline_n on every row, in both the printed table and the snapshot)
    # and a mismatch against devices[0] is printed loudly, never silent.
    for mode in modes:
        print('\n-- %s --' % mode, flush=True)
        for n in devices:
            label = '%s dev=%-3d' % (mode, n)
            cpus, busy = room(nodes, n, a.busy_ceiling, a.i_know_the_box_is_busy,
                              a.wait_for_room)
            if cpus is None:
                print('  %-26s SKIPPED (no %d free cpus within %.0fs)'
                      % (label, n, a.wait_for_room), flush=True)
                skipped.append(dict(label=label, n=n, mode=mode))
                continue
            best = None
            for _ in range(a.repeats):
                r = per_step(case_dir, cpus, nodes, n, mode, a.n_lo, a.n_hi)
                if r is None:
                    continue
                if best is None or r['ms_per_step'] < best['ms_per_step']:
                    best = r
            if best is None:
                print('  %-26s FAILED' % label, flush=True)
                continue
            speedup, bn, note = record_point(base, base_n, mode, n,
                                             best['ms_per_step'], devices[0])
            if note:
                print('  NOTE: ' + note, flush=True)
            best.update(mode=mode, n=n, cpus=cpus, nodes=rs.nodes_of(nodes, cpus),
                        speedup=speedup, baseline_n=bn, busy=busy,
                        loadavg=os.getloadavg(), n_lo=a.n_lo, n_hi=a.n_hi)
            rows.append(best)
            print('  %-26s %9.2f ms/step  speedup %5.2fx (vs n=%d)  fixed %6.2fs  '
                  'threads %3d  wall %.1f/%.1fs  cpus=%s nodes=%s'
                  % (label, best['ms_per_step'], best['speedup'], bn,
                     best['fixed_s'], best['threads'], best['wall_lo_s'],
                     best['wall_hi_s'], cpus, best['nodes']), flush=True)

    el = {r['n']: r['ms_per_step'] for r in rows if r['mode'] == 'element'}
    el_base_n = min(el) if el else None
    print('\n%-14s %10s %10s %8s %8s' % ('devices/ranks', 'fortran', 'jax-shard',
                                         'f-up', 'j-up (vs n=%s)' % el_base_n))
    if el_base_n is not None and el_base_n != devices[0]:
        print('  NOTE: j-up baselines on n=%d (n=%d, the requested first '
              'point, has no "element"-mode measurement) -- not a silent '
              'baseline swap.' % (el_base_n, devices[0]))
    for n in devices:
        f = FORTRAN_MS.get(n)
        j = el.get(n)
        print('%-14d %10s %10s %8s %8s'
              % (n, '%.2f' % f if f else '-', '%.2f' % j if j else '-',
                 '%.2fx' % (FORTRAN_MS[1] / f) if f else '-',
                 '%.2fx' % (el[el_base_n] / j) if j else '-'))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    meta = dict(case=a.case, sha=sha, host=os.uname().nodename,
               date=time.strftime('%Y-%m-%d %H:%M'), n_lo=a.n_lo,
               n_hi=a.n_hi, busy_ceiling=a.busy_ceiling,
               overridden=bool(a.i_know_the_box_is_busy),
               fortran_quoted=FORTRAN_MS, rows=rows, skipped=skipped,
               element_baseline_n=el_base_n)
    json.dump(meta, open(OUT, 'w'), indent=1)
    print('\nsaved %s' % OUT)
    if skipped:
        print('%d point(s) SKIPPED as busy (not measured, not silently dropped): %s'
              % (len(skipped), [s['label'] for s in skipped]))

    # Item 91a's ledger half (2026-09-24): ledger.rows_from_scaling_snapshot
    # now takes `tool` as a parameter and declares a shard point's `ranks` as
    # the EXISTING 'threads' parallelism value (see that function's
    # docstring) -- shard rows file under tool='run_shard_scaling', never
    # mislabelled as run_scaling's own points.
    snap_rel = os.path.relpath(OUT, ROOT).replace(os.sep, '/')
    tenancy = ledger.box_tenancy(a.busy_ceiling)
    nledger = ledger.append_rows(
        ledger.rows_from_scaling_snapshot(meta, snap_rel, tenancy,
                                          tool='run_shard_scaling'))
    print('%d ledger row(s) appended to %s (box tenancy %d/%d cpus over %.2f)'
          % (nledger, ledger.LEDGER_RELPATH, tenancy['busy'], tenancy['total'],
             a.busy_ceiling))


if __name__ == '__main__':
    main()
