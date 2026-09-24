#! /usr/bin/env python3
"""
Straggler-spread probe: IDENTICAL memory-streaming work on k individually
pinned cpus, run concurrently; per-process ms/iter and EFFECTIVE_CORES
(cpu_used / wall). Pathway item 47a -- a REWRITE, 2026-09-24, of the
2026-09-21 session-scratchpad probe behind item 47's 2.68x record, which was
never committed and no longer exists on disk. The original's kernel is known
(docs/SESSION_LOG_2026-09-21_jax-scaling.md:212: `(a*c+0.5).sum()` over a
30.5 MB array, one pinned process per rank); its exact harness is not, so a
number from this tool is NOT the 2026-09-21 number, and item 47's figure was
itself taken under a tenant load (corrected 2026-09-22). No MPI: the
question is the box's per-core throughput spread, not the solver's.

Placement reuses testsys/perf/run_scaling.py's free_node_map / select_cpus /
numactl_prefix verbatim (cpu AND memory bound). Each worker reports
`ms_per_iter` from a timed loop AFTER a warm-up; a worker whose cpu time
fell under --min-effective-cores of its wall time is REPORTED, not dropped
(it is the straggler signal, rule 2). A non-positive timing raises.

    python3 testsys/perf/threadprobe.py --k 16 --iters 200
Writes docs/perf_snapshots/threadprobe_<YYYY-mm-dd_HHMMSS>.json (rule 19).
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

WORKER = r'''
import os, sys, time, json
import numpy as np
n = int(sys.argv[1]); iters = int(sys.argv[2]); warm = int(sys.argv[3])
a = np.ones(n); c = np.full(n, 1.0000001)
for _ in range(warm):
    (a * c + 0.5).sum()
w0 = time.perf_counter(); c0 = time.process_time()
for _ in range(iters):
    (a * c + 0.5).sum()
wall = time.perf_counter() - w0; cpu = time.process_time() - c0
print(json.dumps(dict(wall_s=wall, cpu_s=cpu, iters=iters,
                      cpus_allowed=len(os.sched_getaffinity(0)))))
'''


def summarize(results):
    """Per-worker rows -> spread statistics. Raises (rule 2) on any
    non-positive timing or an empty result set; returns a dict."""
    if not results:
        raise ValueError('threadprobe: no worker results to summarize')
    ms = []
    for r in results:
        if r['wall_s'] <= 0 or r['iters'] <= 0:
            raise ValueError('threadprobe: non-positive timing %r -- an invalid '
                             'measurement, not a fast one' % (r,))
        ms.append(1e3 * r['wall_s'] / r['iters'])
    mean = sum(ms) / len(ms)
    return dict(n=len(ms), ms_min=min(ms), ms_max=max(ms), ms_mean=mean,
                spread_max_over_min=max(ms) / min(ms),
                straggler_max_over_mean=max(ms) / mean,
                effective_cores=[r['cpu_s'] / r['wall_s'] for r in results])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--k', type=int, default=16, help='pinned worker processes')
    ap.add_argument('--mb', type=float, default=30.5, help='array size per worker (MB)')
    ap.add_argument('--iters', type=int, default=200)
    ap.add_argument('--warm', type=int, default=20)
    ap.add_argument('--policy', default='compact', choices=('compact', 'spread'))
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    ap.add_argument('--min-effective-cores', type=float, default=0.95)
    a = ap.parse_args()
    import run_scaling as rs
    import run_numa_scaling as numa
    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology')
    free = rs.free_node_map(nodes, a.busy_ceiling, a.i_know_the_box_is_busy)
    cpus = rs.select_cpus(free, a.k, a.policy)
    if cpus is None:
        raise SystemExit('FAIL: no %d-cpu %s placement among cpus under busy ceiling '
                         '%.2f right now (free nodes %s) -- not measured'
                         % (a.k, a.policy, a.busy_ceiling, sorted(free)))
    n = int(a.mb * 1e6 / 8)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    procs = [subprocess.Popen(rs.numactl_prefix([c], nodes) +
                              [sys.executable, '-c', WORKER, str(n), str(a.iters), str(a.warm)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
             for c in cpus]
    results = []
    for c, p in zip(cpus, procs):
        out, err = p.communicate()
        if p.returncode != 0:
            raise SystemExit('FAIL: worker on cpu %d exited %d: %s' % (c, p.returncode, err[-300:]))
        results.append(dict(json.loads(out.strip().splitlines()[-1]), cpu=c))
    s = summarize(results)
    sha_r = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'],
                           capture_output=True, text=True)
    if sha_r.returncode != 0 or not sha_r.stdout.strip():
        raise SystemExit('FAIL: could not read the git sha (rule 19): %s' % sha_r.stderr.strip())
    low = [r['cpu'] for r, e in zip(results, s['effective_cores']) if e < a.min_effective_cores]
    for r, e in zip(results, s['effective_cores']):
        print('cpu %3d  %8.3f ms/iter  EFFECTIVE_CORES %.2f%s'
              % (r['cpu'], 1e3 * r['wall_s'] / r['iters'], e,
                 '  <-- below %.2f' % a.min_effective_cores if e < a.min_effective_cores else ''))
    print('k=%d %s: spread max/min %.2fx, straggler max/mean %.2fx, loadavg %s'
          % (s['n'], a.policy, s['spread_max_over_min'], s['straggler_max_over_mean'],
             os.getloadavg()))
    out = os.path.join(ROOT, 'docs', 'perf_snapshots',
                       'threadprobe_%s.json' % time.strftime('%Y-%m-%d_%H%M%S'))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(dict(sha=sha_r.stdout.strip(), host=os.uname().nodename,
                   date=time.strftime('%Y-%m-%d %H:%M'), loadavg=os.getloadavg(),
                   args=vars(a), cpus=cpus, results=results, summary=s,
                   below_min_effective_cores=low), open(out, 'w'), indent=1)
    print('saved %s' % os.path.relpath(out, ROOT))


if __name__ == '__main__':
    main()
