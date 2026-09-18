#! /usr/bin/env python3
"""
Diagnose the jax-CPU compact-pinned scaling plateau seen past 16 cores
(pathway_forward item 33 continuation). Report-only; never a gate; never
wired into any test tier.

QUESTION. Measured strong-scaling for jax-compact peaks 2.52x at 16 cores and
FALLS to 1.98x at 32. Is that:
  (a) the duplicate-index scatter-add op specifically failing to
      thread-parallelize on the CPU backend, or
  (b) memory bandwidth saturating for EVERY op at this problem size (so
      scatter is not special), or
  (c) neither -- pin it to a different op via HLO.

This script isolates ONE control elementwise op (no scatter semantics, same
total bytes moved) against ONE scatter-add op (jax's `.at[idx].add`, the CPU
lowering of the same primitive the port's `assembleGlobalKU`/
`calcHourglassResist` use for `np.add.at` under the jax backend), both
jitted, at compact-pinned 8 and 32 cores, and compares each op's own
8->32 speedup ratio. If the control also fails to scale, the plateau is
shared (verdict b); if only scatter fails, scatter is the specific
bottleneck (verdict a).

ARRAY SIZING. The real per-step HLO dump for test.tpv8 (cited in the item-33
mission) showed f64 arrays up to shape (3818584,) among the 556 scatter ops
in one step -- that is this script's N_SOURCE. Per-node-sharing multiplicity
in a hex mesh (`assembleGlobalKU`'s scatter-add of nodal force/mass
contributions) is ~8 elements touching each interior node, so N_TARGET =
N_SOURCE // 8 = 477323 unique destination indices, each duplicated exactly 8
times in a FIXED, once-shuffled index array (`np.random.default_rng(0)`) --
duplicated-but-not-contiguous, which is the harder case for a scatter
implementation (contiguous runs of the same index can sometimes be
partially handled by a segmented reduction; a shuffled pattern cannot).
This is a stated assumption, not a measurement of the real element-to-node
connectivity table -- disclosed here rather than silently presented as exact.

CONTROL SIZING. Same *total bytes* as the scatter's SOURCE array (the array
that dominates memory traffic each op: N_SOURCE * 8 bytes ~= 30.5 MB), so
the two benchmarks move a comparable volume of memory per call and a
bandwidth ceiling (if it exists) should show up in both equally.

METHOD. Per-iteration cost by DIFFERENCE over two iteration counts (n_lo,
n_hi), same discipline as `run_scaling.py`/`run_numa_scaling.py`: each op
loops internally via `lax.fori_loop` (a real HLO while-loop, not
unrolled/algebraically-simplified by XLA) for `iters` steps, one fresh
process per (op, cores, iters) so a compiled function is never reused
across a process boundary and JIT compile time cancels in the difference,
not by accident.  n_lo/n_hi are recorded in every output, exactly to avoid
item 40's failure mode (a from-memory "4.0x" figure with no iteration counts
attached and never reproduced).

PINNING. Reuses `run_scaling.py`'s free-node-aware compact pinning
(`free_node_map`, `select_cpus`, `numactl_prefix`) verbatim -- no second,
divergent pinning implementation. Only COMPACT is measured (the mission
does not ask about locality here, item 33's own compact/spread contrast
already exists in `run_scaling.py`/`run_numa_scaling.py`).

Usage:
    JAX_PLATFORMS=cpu /home/utig5/dliu/gns/gns/venv_cotopaxi/bin/python3 \\
        testsys/perf/probe_scatter_bandwidth.py [--cores 8,32]
        [--n-lo 30] [--n-hi 100] [--busy-ceiling 0.2]
        [--i-know-the-box-is-busy] [--dump-hlo]

Writes testsys/perf/scatter_bandwidth_last.json.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(TESTSYS))
OUT = os.path.join(TESTSYS, 'scatter_bandwidth_last.json')

sys.path.insert(0, TESTSYS)
import run_scaling as rs   # noqa: E402  (free_node_map, select_cpus, numactl_prefix, compact_cpus)
import run_numa_scaling as numa  # noqa: E402  (numa_topology, require_idle)

N_SOURCE = 3818584          # real HLO dump, test.tpv8, one step: f64 array up to this shape
MULTIPLICITY = 8            # stated assumption: ~8 hex elements share an interior node
N_TARGET = N_SOURCE // MULTIPLICITY
assert N_TARGET * MULTIPLICITY == N_SOURCE, 'N_SOURCE must be an exact multiple of MULTIPLICITY'

WORKER = r'''
import os, sys, time
os.environ['JAX_ENABLE_X64'] = '1'
import jax
import jax.numpy as jnp
import numpy as np
from functools import partial

N_SOURCE = %d
N_TARGET = %d
MULT = %d
op = %r
iters = %d

@partial(jax.jit, static_argnames=('iters',))
def run_control(x0, iters):
    def body(i, x):
        return x * 2.0 + 1.0
    return jax.lax.fori_loop(0, iters, body, x0)

@partial(jax.jit, static_argnames=('iters',))
def run_scatter(carry0, idx, src, iters):
    def body(i, carry):
        return carry.at[idx].add(src)
    return jax.lax.fori_loop(0, iters, body, carry0)

if op == 'control':
    x0 = jnp.arange(N_SOURCE, dtype=jnp.float64) * 1e-6
    r = run_control(x0, iters)
    r.block_until_ready()
else:
    rng = np.random.default_rng(0)
    idx_np = np.repeat(np.arange(N_TARGET, dtype=np.int32), MULT)
    rng.shuffle(idx_np)                       # fixed, once-shuffled -- see module docstring
    idx = jnp.asarray(idx_np)
    src = jnp.arange(N_SOURCE, dtype=jnp.float64) * 1e-6
    carry0 = jnp.zeros(N_TARGET, dtype=jnp.float64)
    r = run_scatter(carry0, idx, src, iters)
    r.block_until_ready()

t0 = time.time()
if op == 'control':
    r = run_control(x0, iters)
else:
    r = run_scatter(carry0, idx, src, iters)
r.block_until_ready()
print('WALL', time.time() - t0)
'''
# NOTE: the block above times a SECOND call after a warm-up call above it --
# see time_one() docstring for why the warmup is INSIDE this same process
# rather than a cancel-by-difference across two processes for the whole
# thing: fori_loop trip count is a static (jit) argument, so a second call
# with the SAME iters in the SAME process hits the compilation cache and
# times pure execution once compiled. Combined with the OUTER n_lo/n_hi
# difference (two separate fresh processes at different iters), this
# cancels both (1) this process's own interpreter/import/compile cost via
# the outer difference AND (2) leaves no ambiguity about whether the timed
# region includes a compile. Both n_lo and n_hi runs pay their own compile
# during the (untimed) warmup call.


def time_one(op, iters, cpus, node_map, dump_hlo_dir=None):
    script = WORKER % (N_SOURCE, N_TARGET, MULTIPLICITY, op, iters)
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    env.pop('XLA_FLAGS', None)
    env.pop('OMP_NUM_THREADS', None)
    env.pop('OPENBLAS_NUM_THREADS', None)
    if dump_hlo_dir:
        os.makedirs(dump_hlo_dir, exist_ok=True)
        env['XLA_FLAGS'] = '--xla_dump_to=%s' % dump_hlo_dir
    cmd = rs.numactl_prefix(cpus, node_map) + [sys.executable, '-c', script]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        return None
    for line in r.stdout.splitlines():
        if line.startswith('WALL'):
            return float(line.split()[1])
    return None


def per_iter(op, cpus, node_map, n_lo, n_hi, dump_hlo_dir=None):
    t_lo = time_one(op, n_lo, cpus, node_map)
    t_hi = time_one(op, n_hi, cpus, node_map, dump_hlo_dir=dump_hlo_dir)
    if t_lo is None or t_hi is None:
        return None, None
    per = (t_hi - t_lo) / float(n_hi - n_lo)
    fixed = t_lo - n_lo * per
    return per, fixed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cores', default='8,32')
    ap.add_argument('--n-lo', type=int, default=30)
    ap.add_argument('--n-hi', type=int, default=100)
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    ap.add_argument('--dump-hlo', action='store_true',
                     help='also dump XLA HLO for the scatter op at the largest '
                          'core count requested (nice-to-have corroboration)')
    a = ap.parse_args()
    cores = [int(x) for x in a.cores.split(',') if x]
    for c in cores:
        if c > 32:
            raise SystemExit('FAIL: %d cores requested -- owner constraint is '
                             'at most 32 (one socket), never 64.' % c)

    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology.')
    print('topology: %d NUMA node(s) x %d cpu(s)' % (len(nodes), len(nodes[min(nodes)])))
    print('N_SOURCE=%d  N_TARGET=%d  MULTIPLICITY=%d  n_lo=%d n_hi=%d'
          % (N_SOURCE, N_TARGET, MULTIPLICITY, a.n_lo, a.n_hi))

    rows = {}
    for op in ('control', 'scatter'):
        rows[op] = {}
        print('\n-- op=%s --' % op)
        for k in cores:
            free = rs.free_node_map(nodes, a.busy_ceiling, a.i_know_the_box_is_busy)
            cpus = rs.select_cpus(free, k, 'compact')
            label = 'op=%-7s cores=%-3d' % (op, k)
            if cpus is None:
                print('  %s SKIPPED (no %d-cpu compact placement among free nodes %s)'
                      % (label, k, sorted(free)))
                rows[op][k] = dict(skipped='no_free_node_placement', free_nodes=sorted(free))
                continue
            chk = numa.require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
            if chk is None:
                print('  %s SKIPPED (cpus busy, not overridden)' % label)
                rows[op][k] = dict(skipped='busy', cpus=cpus)
                continue
            dump_dir = None
            if a.dump_hlo and op == 'scatter' and k == max(cores):
                dump_dir = tempfile.mkdtemp(prefix='hlo_scatter_%dc_' % k,
                                            dir=os.environ.get('TMPDIR', '/tmp'))
            per, fixed = per_iter(op, cpus, nodes, a.n_lo, a.n_hi, dump_hlo_dir=dump_dir)
            if per is None:
                print('  %s FAILED' % label)
                rows[op][k] = dict(failed=True, cpus=cpus)
                continue
            used_nodes = rs.nodes_of(nodes, cpus)
            rows[op][k] = dict(us_per_iter=per * 1e6, fixed_s=fixed, cpus=cpus,
                                nodes=used_nodes, hlo_dump=dump_dir)
            print('  %s %10.2f us/iter  fixed %6.2fs  cpus=%s nodes=%s%s'
                  % (label, per * 1e6, fixed, cpus, used_nodes,
                     ('  hlo->%s' % dump_dir) if dump_dir else ''))

    print('\n-- speedup 8->32 (per-op, higher is better; ~4x would be linear) --')
    verdicts = {}
    for op in ('control', 'scatter'):
        r8 = rows[op].get(8)
        r32 = rows[op].get(32)
        if r8 and r32 and 'us_per_iter' in r8 and 'us_per_iter' in r32:
            ratio = r8['us_per_iter'] / r32['us_per_iter']
            verdicts[op] = ratio
            print('  %-8s 8c=%9.2f us/iter  32c=%9.2f us/iter  speedup=%.2fx'
                  % (op, r8['us_per_iter'], r32['us_per_iter'], ratio))
        else:
            print('  %-8s inconclusive (missing 8c or 32c measurement)' % op)

    if 'control' in verdicts and 'scatter' in verdicts:
        c, s = verdicts['control'], verdicts['scatter']
        print('\nVERDICT:')
        if c < 1.5 and s < 1.5:
            print('  (b) memory bandwidth: BOTH ops fail to scale past 8 cores '
                  '(control=%.2fx, scatter=%.2fx) -- scatter is not special.' % (c, s))
        elif c >= 1.5 and s < c * 0.6:
            print('  (a) scatter-specific: control scales (%.2fx) but scatter does '
                  'not keep pace (%.2fx) -- the scatter-add op is the CPU bottleneck.'
                  % (c, s))
        else:
            print('  (c) inconclusive by this test alone: control=%.2fx scatter=%.2fx '
                  '-- does not cleanly separate; see notes for a follow-up.' % (c, s))

    sha = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'],
                         capture_output=True, text=True).stdout.strip()
    meta = dict(sha=sha, host=os.uname().nodename, date=time.strftime('%Y-%m-%d %H:%M'),
               n_source=N_SOURCE, n_target=N_TARGET, multiplicity=MULTIPLICITY,
               n_lo=a.n_lo, n_hi=a.n_hi, busy_ceiling=a.busy_ceiling,
               overridden=bool(a.i_know_the_box_is_busy), rows=rows,
               speedup_8_to_32=verdicts)
    json.dump(meta, open(OUT, 'w'), indent=1)
    print('\nprovenance: sha=%s host=%s; saved %s' % (sha, os.uname().nodename, OUT))


if __name__ == '__main__':
    main()
