#!/usr/bin/env python3
"""Spread-vs-packed placement arm, test.tpv104, 16 ranks.

The question: the ~50 ms/step "exchange" the 1D jax-MPI slab reports at 16
ranks -- is it halo cost (the surface-to-volume argument item 43 was priced
on), or is it memory-bandwidth starvation from least_loaded_cpus packing all
16 ranks onto 2-3 of this box's 8 NUMA nodes?

Discriminator: same tool, same case, same code path, ONLY the cpu set differs.
  spread = 2 ranks on each of the 8 NUMA nodes (4x the memory bandwidth)
  packed = the tool's own default selection
Halo cost predicts spread is the same or WORSE (more inter-node hops).
Bandwidth starvation predicts the exchange figure collapses.

Fortran runs in both arms as the control: it is 3D-decomposed and was not
bandwidth-limited at this size, so if spread moves Fortran too the effect is
box-level, not port-specific.
"""
import json
import os
import subprocess
import sys
import time

W = '/home/utig5/dliu/EQdyna/.claude/worktrees/wei-s3'
PERF = os.path.join(W, 'testsys', 'perf')
OUT = os.path.dirname(os.path.abspath(__file__))
NODES = {n: list(range(n * 8, n * 8 + 8)) for n in range(8)}
NEVER = {0, 1}          # measured: read busy 0.00, deliver 0.39 effective cores


def busy_now(dur=3.0):
    def snap():
        d = {}
        for ln in open('/proc/stat'):
            if ln.startswith('cpu') and ln[3].isdigit():
                f = ln.split()
                c = int(f[0][3:])
                v = [int(x) for x in f[1:]]
                d[c] = (sum(v), v[3] + v[4])
        return d
    a = snap()
    time.sleep(dur)
    b = snap()
    out = {}
    for c in a:
        dt = b[c][0] - a[c][0]
        di = b[c][1] - a[c][1]
        out[c] = (dt - di) / dt if dt else 0.0
    return out


def spread_set(per_node=2, ceil=0.10):
    """The `per_node` least-loaded cpus on EVERY node. Fails loudly if a node
    cannot supply them under `ceil` -- a spread arm with a hole in it is a
    different experiment, not a slightly worse one."""
    b = busy_now()
    chosen, short = [], []
    for n, cs in sorted(NODES.items()):
        free = sorted((b[c], c) for c in cs if c not in NEVER and b[c] <= ceil)
        if len(free) < per_node:
            short.append((n, len(free)))
        chosen += [c for _, c in free[:per_node]]
    if short:
        return None, short, b
    return sorted(chosen), None, b


def run(tag, extra, ranks, repeats=2):
    log = os.path.join(OUT, 'spread_arm_%s.log' % tag)
    cmd = [sys.executable, os.path.join(PERF, 'run_mpi_scaling.py'),
           '--case', 'test.tpv104', '--ranks', str(ranks),
           '--n-lo', '20', '--n-hi', '60', '--max-busy', '0.2',
           '--syncs', 'halo', '--warmup', '--repeats', str(repeats)] + extra
    print('\n===== %s =====\n%s' % (tag, ' '.join(cmd)), flush=True)
    env = dict(os.environ)
    env['EQDYNAROOT'] = W
    env['PATH'] = '%s/bin:%s/scripts:%s' % (W, W, env.get('PATH', ''))
    t0 = time.time()
    with open(log, 'w') as fh:
        fh.write(' '.join(cmd) + '\n')
        fh.flush()
        r = subprocess.run(cmd, cwd=PERF, env=env, stdout=fh,
                           stderr=subprocess.STDOUT)
    print('%s rc=%d  %.1f min  -> %s' % (tag, r.returncode,
                                         (time.time() - t0) / 60.0, log),
          flush=True)
    return r.returncode


def arm_skipped(tag):
    """Did the run SKIP on the busy ceiling rather than measure?"""
    log = os.path.join(OUT, 'spread_arm_%s.log' % tag)
    return 'SKIPPED' in open(log).read()


def main():
    print('start %s' % time.strftime('%H:%M:%S'), flush=True)

    # --- spread arm, retried on a ceiling skip with a freshly measured set ---
    done = False
    for attempt in range(1, 9):
        cpus, short, b = spread_set()
        if cpus is None:
            print('attempt %d: nodes short of 2 free cpus: %s  '
                  'tenancy>0.2 %d/64  median busy %.3f  max busy %.3f '
                  '-- waiting 120s'
                  % (attempt, short,
                     sum(1 for v in b.values() if v > 0.2),
                     sorted(b.values())[len(b) // 2], max(b.values())),
                  flush=True)
            time.sleep(120)
            continue
        print('attempt %d: spread cpus %s  busy %s  tenancy %d/64'
              % (attempt, cpus, {c: round(b[c], 3) for c in cpus},
                 sum(1 for v in b.values() if v > 0.2)), flush=True)
        tag = 'spread16_a%d' % attempt
        run(tag, ['--cpus', ','.join(str(c) for c in cpus)], 16)
        if not arm_skipped(tag):
            done = True
            break
        print('attempt %d SKIPPED on the ceiling -- re-measuring' % attempt,
              flush=True)
        time.sleep(120)
    if not done:
        print('SPREAD ARM DID NOT MEASURE -- not reporting a number', flush=True)

    # --- packed arm, same window, the tool's own default selection ---
    run('packed16', ['--exclude-cpus', '0,1'], 16)

    print('done %s' % time.strftime('%H:%M:%S'), flush=True)


if __name__ == '__main__':
    main()
