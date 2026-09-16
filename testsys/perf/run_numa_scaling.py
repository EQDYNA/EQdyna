#! /usr/bin/env python3
"""
JAX-CPU core scaling, separating CORE COUNT from NUMA LOCALITY.
Report-only. Never a gate. pathway_forward item 33.

WHY THIS EXISTS. Every core-scaling number this repo quotes is stale,
unverified, or taken under load:

  * `2.00x / 3.96x / 7.48x on 2/4/8 cores` and `~2.5x plateau when
    memory-bound` PREDATE the Python restructure and the HLO-literal fix. They
    have not been reproduced since and were quoted in 2026-09-16's session as
    if current.
  * The one figure measured on current code is a single ratio: test.tpv8,
    193 ms/step pinned to 1 core vs 13.8 ms/step unpinned = 14x. One point is
    not a curve.

AND THE OLD CURVE IS MECHANICALLY SUSPECT. This box is 8 NUMA nodes of 8 cores
(2 sockets x 32, AMD EPYC 7532, 1 thread/core, ~128 GB/node). The inherited
knee sits at EXACTLY 8 cores -- the NUMA node size. So it may be measuring
remote-memory cost rather than any property of the solver, and nobody has
separated the two. That is what the SPREAD configuration below is for: the same
8 cores, one per node, maximum NUMA distance. If 8-within-a-node is fast and
8-across-nodes is slow, the knee is locality. If both are the same, the knee is
the solver running out of parallel work.

WHY IT REFUSES TO RUN ON A BUSY BOX. A memory-bandwidth measurement taken while
ten other users are streaming memory measures them, not us. This is not a
caveat to note in the output -- it invalidates the number, so the script FAILS
(rule 2: a check that cannot run must not produce a result that looks like
one). Override only if you know why you are doing it.

METHOD. Per-step cost by DIFFERENCE over two step counts, so the fixed cost
(interpreter, case load, XLA compile) cancels exactly -- the same technique
run_perf.py uses, and for the same reason: at 114 steps compile is 15% of a
jax run, so a total-wall-clock number drifts when XLA changes and the solver
does not.

CASE CHOICE MATTERS. Default is test.tpv104 (2701 fault nodes), NOT test.tpv8
(1891). tpv8 is small enough that running out of parallel work is a plausible
alternative explanation for the knee, which would confound the whole point.

Usage:
    python3 testsys/perf/run_numa_scaling.py [--case test.tpv104]
                                             [--steps 114] [--factor 3]
                                             [--load-ceiling 1.0]
                                             [--i-know-the-box-is-busy]
"""
import argparse
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
PYTHON_PKG = os.path.join(REPO_ROOT, 'src', 'python')
OUT = os.path.join(TESTSYS, 'perf', 'numa_scaling_last.json')


def numa_topology():
    """{node: [cpu, ...]} from numactl, or {} if unavailable."""
    try:
        txt = subprocess.run(['numactl', '--hardware'], capture_output=True,
                             text=True).stdout
    except FileNotFoundError:
        return {}
    nodes = {}
    for m in re.finditer(r'^node (\d+) cpus: (.+)$', txt, re.M):
        nodes[int(m.group(1))] = [int(c) for c in m.group(2).split()]
    return nodes


def require_idle(ceiling, override):
    """Refuse to measure bandwidth on a contended box."""
    load1, load5, load15 = os.getloadavg()
    others = set()
    ps = subprocess.run(['ps', '-eo', 'user='], capture_output=True, text=True)
    me = os.environ.get('USER', '')
    for u in ps.stdout.split():
        if u != me and not u.startswith('_') and u not in (
                'root', 'daemon', 'systemd+', 'message+', 'syslog', 'nobody',
                'dbus', 'chrony', 'polkitd', 'uuidd', 'statd', 'rtkit',
                'colord', 'avahi', 'kernoops', 'whoopsie', 'gdm', 'lightdm',
                'libstor+', 'systemd-network', 'systemd-resolve'):
            others.add(u)
    if load1 <= ceiling:
        print('  load %.2f (1 min) <= %.2f -- box is quiet enough' % (load1, ceiling))
        return load1, load5, load15, sorted(others)
    msg = ('REFUSING TO MEASURE: load average is %.2f / %.2f / %.2f and the '
           'ceiling is %.2f.\n'
           '  Other users with processes: %s\n'
           '  A core-scaling measurement is a MEMORY BANDWIDTH measurement. '
           'Taken under someone\n'
           '  else\'s memory traffic it measures them, and the number is not '
           'wrong-but-usable,\n'
           '  it is meaningless. pathway item 33 exists because the previous '
           'figures were taken\n'
           '  exactly this way.\n'
           '  Wait for an idle box, or pass --i-know-the-box-is-busy to record '
           'a number that\n'
           '  must NOT be quoted as this repo\'s scaling curve.'
           % (load1, load5, load15, ceiling, ', '.join(sorted(others)) or 'none'))
    if not override:
        raise SystemExit('FAIL: ' + msg)
    print('  WARNING, OVERRIDDEN: ' + msg)
    return load1, load5, load15, sorted(others)


def build_case(case_name):
    sys.path.insert(0, os.path.join(REPO_ROOT, 'testsys', 'e2e'))
    import run_e2e                                        # noqa: E402
    d = os.path.join(TESTSYS, 'perf', 'numa_case', case_name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(os.path.dirname(d), exist_ok=True)
    run_e2e.make_serial_case(case_name, d, run_e2e.base_env())
    return d


def time_one(case_dir, nsteps, cpus):
    """Wall seconds for `nsteps` on exactly `cpus`, in a fresh process.

    Fresh process on purpose: jax caches compiled functions in-process, so a
    second call in the same interpreter pays no compile and the difference
    below would cancel the wrong term.
    """
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from eqdyna import eqdyna3d\n"
        "import time; t0 = time.time()\n"
        "eqdyna3d.run_case(%r, nsteps=%d, verbose=False, backend='jax')\n"
        "print('WALL', time.time() - t0)\n" % (PYTHON_PKG, case_dir, nsteps))
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    # Let XLA use the cores taskset gives it; do NOT force single-threaded.
    env.pop('XLA_FLAGS', None)
    env.pop('OMP_NUM_THREADS', None)
    cpulist = ','.join(str(c) for c in cpus)
    r = subprocess.run(['taskset', '-c', cpulist, sys.executable, '-c', script],
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(r.stderr[-1500:])
        return None
    for line in r.stdout.splitlines():
        if line.startswith('WALL'):
            return float(line.split()[1])
    return None


def per_step(case_dir, cpus, n_lo, n_hi):
    t_lo = time_one(case_dir, n_lo, cpus)
    t_hi = time_one(case_dir, n_hi, cpus)
    if t_lo is None or t_hi is None:
        return None, None
    ps = (t_hi - t_lo) / float(n_hi - n_lo)
    fixed = t_lo - n_lo * ps
    return ps, fixed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--steps', type=int, default=114)
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--load-ceiling', type=float, default=1.0)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    a = ap.parse_args()

    print('JAX-CPU core scaling: core count vs NUMA locality (item 33)')
    nodes = numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology; this '
                         'experiment is about NUMA and cannot run blind.')
    per_node = len(nodes[min(nodes)])
    print('  topology: %d NUMA node(s) x %d cpu(s)' % (len(nodes), per_node))
    load = require_idle(a.load_ceiling, a.i_know_the_box_is_busy)

    node0 = nodes[min(nodes)]
    configs = []
    k = 1
    while k <= per_node:
        configs.append(('within-node-%d' % k, node0[:k]))
        k *= 2
    # same core count as the largest within-node point, maximum NUMA distance
    spread = [nodes[n][0] for n in sorted(nodes)][:per_node]
    if len(spread) == per_node:
        configs.append(('spread-%d-one-per-node' % per_node, spread))
    # one full socket, if the topology has at least 4 nodes
    if len(nodes) >= 4:
        sock = [c for n in sorted(nodes)[:len(nodes) // 2] for c in nodes[n]]
        configs.append(('socket-%d' % len(sock), sock))
        allc = [c for n in sorted(nodes) for c in nodes[n]]
        configs.append(('all-%d' % len(allc), allc))

    print('  case    : %s (built fresh, serial)' % a.case)
    case_dir = build_case(a.case)
    n_lo, n_hi = a.steps, a.steps * a.factor
    print('  method  : per-step by difference, %d vs %d steps, fresh process each'
          % (n_lo, n_hi))
    print()

    results = {}
    base = None
    print('  %-26s %6s  %12s  %10s  %s'
          % ('configuration', 'cpus', 'ms/step', 'speedup', 'fixed cost'))
    for label, cpus in configs:
        ps, fixed = per_step(case_dir, cpus, n_lo, n_hi)
        if ps is None:
            print('  %-26s %6d  FAILED' % (label, len(cpus)))
            continue
        if base is None:
            base = ps
        results[label] = dict(cpus=cpus, ms_per_step=ps * 1e3,
                              speedup=base / ps, fixed_s=fixed)
        print('  %-26s %6d  %12.2f  %9.2fx  %8.2f s'
              % (label, len(cpus), ps * 1e3, base / ps, fixed))

    w = results.get('within-node-%d' % per_node)
    sp = results.get('spread-%d-one-per-node' % per_node)
    print('\n  VERDICT on the knee at %d cores:' % per_node)
    if w and sp:
        ratio = sp['ms_per_step'] / w['ms_per_step']
        print('    %d cores within one node : %8.2f ms/step' % (per_node, w['ms_per_step']))
        print('    %d cores one per node     : %8.2f ms/step  (%.2fx slower)'
              % (per_node, sp['ms_per_step'], ratio))
        if ratio > 1.25:
            print('    -> LOCALITY. The same core count is materially slower when '
                  'spread across\n       NUMA nodes, so the knee is remote memory, '
                  'not the solver running out of work.')
        else:
            print('    -> NOT locality. Same core count performs the same spread '
                  'or packed, so the\n       knee is the solver running out of '
                  'parallel work at this problem size.')
    else:
        print('    inconclusive: one of the two configurations did not run.')

    payload = dict(
        item='pathway_forward 33', case=a.case, steps=[n_lo, n_hi],
        topology={str(k): v for k, v in nodes.items()},
        load_avg='%.2f %.2f %.2f' % load[:3], other_users=load[3],
        host=socket.gethostname(), platform=platform.platform(),
        timestamp=time.strftime('%Y-%m-%d %H:%M:%S'),
        overridden=bool(a.i_know_the_box_is_busy), results=results)
    with open(OUT, 'w') as f:
        json.dump(payload, f, indent=2)
    print('\n  wrote %s' % OUT)
    if a.i_know_the_box_is_busy:
        print('  NOTE: taken on a busy box by override -- do NOT quote as the '
              'scaling curve.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
