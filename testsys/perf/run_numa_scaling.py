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
another process is streaming memory ON THE SAME CPUS measures that process,
not us. This is not a caveat to note in the output -- it invalidates the
number, so the script FAILS (rule 2: a check that cannot run must not produce
a result that looks like one). Override only if you know why you are doing it.

THE CEILING IS PER-CPU, NOT A WHOLE-BOX LOAD AVERAGE (F3, 2026-09-17). The
first version of this ceiling was `os.getloadavg()[0] <= 1.0`, which on this
64-core box refused at load 1.28 -- one core busy (the owner's own `train.py`,
pinned to a single core for 1d21h) and 62 free. That is not conservative, it
is the wrong QUESTION: load average is an absolute queue length, not a
fraction of capacity, so a fixed ceiling of 1.0 is unmeasurable by
construction on any machine with even one long-running background process,
regardless of how quiet everything else is. What this experiment actually
needs to know is "are the SPECIFIC cpus I am about to pin to free" -- a
foreign process on a different NUMA node is irrelevant to a measurement that
never touches it. So the gate now samples `/proc/stat` for exactly the cpu
list a given configuration is about to use (see `cpu_busy_fractions`) and
refuses per-configuration, not once globally: a busy core on socket 1 skips
only the configurations that would have used it, not the whole run.

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
                                             [--busy-ceiling 0.2]
                                             [--i-know-the-box-is-busy]
"""
import argparse
import json
import os
import platform
import re
import socket
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
PYTHON_PKG = os.path.join(REPO_ROOT, 'src', 'python')
OUT = os.path.join(TESTSYS, 'perf', 'numa_scaling_last.json')

sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(TESTSYS, 'perf'))
import perflib  # noqa: E402  (acquire_case_lock, rebuild_serial_case)

# What a second concurrent invocation costs, printed by the refusal (item 77).
LOCK_CONSEQUENCE = [
    'A second invocation in this checkout rmtrees and rebuilds that SAME case'
    ' directory',
    'while the first is TIMING jax out of it. The first does not crash: it'
    ' reports',
    'seconds, and the collision arrives as a NUMA locality effect that is'
    ' really a',
    'rebuilt case (rule 21a, pathway item 77).',
    '',
    'NOT waiting. NOT rebuilding anyway. NOT falling back to a second case'
    ' directory --',
    'each of those is the silent fallback rule 2 forbids. Run your perf tool in'
    ' its own',
    'git worktree (rule 21a), or wait for the holder above to finish.']

# One lock per process; see run_scaling._case_lock for why the memo exists.
_case_lock = None


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


def _other_users():
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
    return sorted(others)


def cpu_busy_fractions(cpus, sample_s=0.3):
    """Busy fraction (0..1) for EXACTLY the requested cpu ids, from a short
    /proc/stat sample -- not the whole-box load average (see module
    docstring, F3). Returns {} if per-cpu stats could not be read for these
    ids at all; the caller must treat that as a hard failure, not as
    'idle' (rule 2: a check that cannot evaluate must fail, not default to a
    pass)."""
    def read():
        rows = {}
        try:
            with open('/proc/stat') as f:
                for line in f:
                    if not line.startswith('cpu') or line[3] == ' ':
                        continue
                    parts = line.split()
                    cpu_id = int(parts[0][3:])
                    nums = [int(x) for x in parts[1:]]
                    idle = nums[3] + nums[4]        # idle + iowait
                    rows[cpu_id] = (idle, sum(nums))
        except (OSError, ValueError, IndexError):
            return None
        return rows
    t0 = read()
    time.sleep(sample_s)
    t1 = read()
    if t0 is None or t1 is None:
        return {}
    out = {}
    for c in cpus:
        if c not in t0 or c not in t1:
            continue
        d_idle = t1[c][0] - t0[c][0]
        d_total = t1[c][1] - t0[c][1]
        out[c] = 0.0 if d_total <= 0 else max(0.0, min(1.0, 1.0 - d_idle / d_total))
    return out


def require_idle(cpus, busy_ceiling, override):
    """Refuse to measure THESE cpus if any one of them is busier than
    `busy_ceiling` (a 0..1 fraction), using per-cpu utilisation rather than
    a whole-box load average -- a foreign process on cpus this configuration
    never touches does not disqualify it.

    Returns (load1, load5, load15, other_users, busy_dict) if free enough (or
    overridden). Returns None if refused and NOT overridden -- the caller
    skips just this configuration; item 33 stays runnable on every OTHER cpu
    set even while one is busy. Raises SystemExit only if utilisation could
    not be measured AT ALL for these cpus -- that is a hard failure, distinct
    from 'measured, and busy'.
    """
    load1, load5, load15 = os.getloadavg()
    others = _other_users()
    busy = cpu_busy_fractions(cpus)
    if not busy:
        raise SystemExit(
            'FAIL: could not read per-cpu utilisation for cpus %s from '
            '/proc/stat -- a check that cannot evaluate must fail, not '
            'default to "idle" (rule 2).' % cpus)
    worst_cpu = max(busy, key=busy.get)
    worst = busy[worst_cpu]
    print('    cpus %s busy: %s  (whole-box load %.2f/%.2f/%.2f, other users: %s)'
          % (cpus, {c: '%.0f%%' % (f * 100) for c, f in sorted(busy.items())},
             load1, load5, load15, ', '.join(others) or 'none'))
    if worst <= busy_ceiling:
        print('    cpu %d worst at %.0f%% <= %.0f%% ceiling -- free enough'
              % (worst_cpu, worst * 100, busy_ceiling * 100))
        return load1, load5, load15, others, busy
    msg = ('REFUSING TO MEASURE cpus %s: cpu %d is %.0f%% busy, ceiling is '
           '%.0f%%.\n'
           '    Other users with processes: %s\n'
           '    A core-scaling measurement is a MEMORY BANDWIDTH measurement '
           'on THESE cpus.\n'
           '    Taken while one of them is busy the number measures that '
           'process, not us --\n'
           '    not wrong-but-usable, meaningless. pathway item 33 exists '
           'because earlier\n'
           '    figures were taken exactly this way.\n'
           '    Wait for these cpus to free up, or pass '
           '--i-know-the-box-is-busy to record\n'
           '    a number that must NOT be quoted as this repo\'s scaling '
           'curve.'
           % (cpus, worst_cpu, worst * 100, busy_ceiling * 100,
              ', '.join(others) or 'none'))
    if override:
        print('    WARNING, OVERRIDDEN: ' + msg)
        return load1, load5, load15, others, busy
    print('    ' + msg)
    return None


def build_case(case_name):
    """create.newcase + forced serial + case.setup, under testsys/perf/.

    GATE 0 (item 77): the exclusive lock on the directory holding the case,
    taken before anything is imported, created or deleted. It lives here, in
    the function that does the rmtree, rather than in main(), for the same
    reason `run_scaling.build_py_case` does: this module is IMPORTED by
    `run_scaling`, `run_mpi_scaling`, `run_shard_scaling`,
    `probe_scatter_bandwidth` and `run_jaxmpi_ab`, so a future caller reaching
    `numa.build_case` directly must be guarded by the same acquire and not by
    whatever this file's main() happens to do.

    `perflib.acquire_case_lock` derives the resource from the directory
    actually about to be destroyed. The lock is KEPT after this returns (the
    sweep then times jax out of this directory) and memoised per process, for
    the reason `run_scaling.build_py_case` records: flock is per open file
    description, so a second acquire here would refuse against our own pid.
    """
    global _case_lock
    d = os.path.join(TESTSYS, 'perf', 'numa_case', case_name)
    if _case_lock is None:
        _case_lock = perflib.acquire_case_lock(d, LOCK_CONSEQUENCE)
    return perflib.rebuild_serial_case(case_name, d)


def _nodes_of(node_map, cpus):
    """Sorted NUMA node ids that any of `cpus` belongs to. Small local copy
    of `run_scaling.nodes_of` -- not imported from there because
    run_scaling.py itself imports THIS module (`import run_numa_scaling as
    numa`), so the reverse import would be circular."""
    want = set(cpus)
    return sorted(n for n, cs in node_map.items() if want & set(cs))


def time_one(case_dir, nsteps, cpus, node_map):
    """Wall seconds for `nsteps` on exactly `cpus`, in a fresh process.

    Fresh process on purpose: jax caches compiled functions in-process, so a
    second call in the same interpreter pays no compile and the difference
    below would cancel the wrong term.

    Item 91c: pinned with `numactl --physcpubind=... --membind=...`, not bare
    `taskset -c`. Measured by `run_scaling.py` (module docstring there,
    :18-26): taskset binds the cpu mask only, not memory -- first-touch
    allocation can still land on a remote node, which defeats the entire
    point of a tool whose job is separating core count from NUMA locality.
    `node_map` is this process's own `numa_topology()` result; the nodes
    bound are exactly the ones `cpus` belongs to, never guessed.
    """
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from eqdyna import eqdyna3d\n"
        "import time; t0 = time.time()\n"
        "eqdyna3d.run_case(%r, nsteps=%d, verbose=False, backend='jax')\n"
        "print('WALL', time.time() - t0)\n" % (PYTHON_PKG, case_dir, nsteps))
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    # Let XLA use the cores numactl gives it; do NOT force single-threaded.
    env.pop('XLA_FLAGS', None)
    env.pop('OMP_NUM_THREADS', None)
    nodes = _nodes_of(node_map, cpus)
    cmd = ['numactl', '--physcpubind=%s' % ','.join(str(c) for c in cpus),
           '--membind=%s' % ','.join(str(n) for n in nodes),
           sys.executable, '-c', script]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(r.stderr[-1500:])
        return None
    for line in r.stdout.splitlines():
        if line.startswith('WALL'):
            return float(line.split()[1])
    return None


def per_step_and_fixed(t_lo, t_hi, n_lo, n_hi):
    """(seconds per step, fixed seconds) from two wall times over two step
    counts -- the family's per-step-by-difference arithmetic, in ONE place.

    The two lines below were written out independently at four call sites
    (`per_step` here, `run_scaling.per_step_py`, `run_scaling.per_step_fortran`
    and `run_shard_scaling.per_step`), so a correction to either of them landed
    in one tool and stayed missing in the other three. The semantics are
    untouched: `ps` is the slope over the step-count difference and `fixed` the
    intercept, i.e. everything that does NOT scale with the step count
    (interpreter start, case load, XLA compile, MPI init, netCDF open) and that
    the difference is taken in order to cancel.

    ITEM 91b (2026-09-24): this function USED to not judge the result at all,
    on the theory that only `run_mpi_scaling.per_step_jax_mpi`'s per-rank
    solve-time differencing needed the check. That theory did not survive
    contact with the two-wall-clock case: a wall-clock difference can go
    non-positive for exactly the same reason (fixed cost dominating, box
    noise) and it is measured here too, by `run_shard_scaling.per_step` and
    this module's own `per_step`, above. `run_mpi_scaling.per_step_jax_mpi`
    keeps its own copy of the check because it differences the RANKS' OWN
    SOLVE TIME rather than two wall clocks -- a different quantity, so a
    second raise here is not a duplicate of that one.
    """
    ps = (t_hi - t_lo) / float(n_hi - n_lo)
    if ps <= 0:
        raise RuntimeError(
            'per-step by difference came out %.6f s/step (t_lo=%.3fs at '
            'n_lo=%d, t_hi=%.3fs at n_hi=%d). That is not a slow '
            'measurement, it is an invalid one -- treat a non-positive '
            'result as a bug to raise on, not a number to report.'
            % (ps, t_lo, n_lo, t_hi, n_hi))
    return ps, t_lo - n_lo * ps


def record_baseline(state, label, ps, first_label):
    """Pure bookkeeping for the single speedup baseline across all
    configurations (item 91g), factored out of main()'s loop so it is
    directly unit-testable with no jax/numactl involved.

    `state` is a dict this function fills in place the first time it is
    called (keys 'value', 'label'). Returns (speedup, baseline_label,
    note); `note` is a non-None warning string exactly when the baseline
    being set is NOT `first_label` -- i.e. the intended first configuration
    was skipped or failed and every speedup on this run is grounded on a
    later configuration instead. The caller must not silently drop `note`."""
    note = None
    if 'value' not in state:
        state['value'] = ps
        state['label'] = label
        if label != first_label:
            note = ('baseline is %r (the first configuration, %r, was '
                    'skipped or failed) -- speedup values are relative to '
                    '%r, not %r' % (label, first_label, label, first_label))
    return state['value'] / ps, state['label'], note


def per_step(case_dir, cpus, node_map, n_lo, n_hi):
    t_lo = time_one(case_dir, n_lo, cpus, node_map)
    t_hi = time_one(case_dir, n_hi, cpus, node_map)
    if t_lo is None or t_hi is None:
        return None, None
    return per_step_and_fixed(t_lo, t_hi, n_lo, n_hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--steps', type=int, default=114)
    ap.add_argument('--factor', type=int, default=3)
    ap.add_argument('--busy-ceiling', type=float, default=0.2,
                     help='max per-cpu busy fraction (0..1) tolerated on the '
                          'exact cpus a configuration is about to use')
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    a = ap.parse_args()

    print('JAX-CPU core scaling: core count vs NUMA locality (item 33)')
    nodes = numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology; this '
                         'experiment is about NUMA and cannot run blind.')
    per_node = len(nodes[min(nodes)])
    print('  topology: %d NUMA node(s) x %d cpu(s)' % (len(nodes), per_node))

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
    # Socket-fill region. Item 33 asks for BOTH 16 and 32 within one socket --
    # jumping straight from spread-8 to a full socket left a 4x gap with no
    # point in between (F2, found 2026-09-17). The owner's constraint is "at
    # most 32, not 64": a single process has no reason to ever span both
    # sockets, so there is deliberately no all-cores config here (F1, same
    # date -- the previous version generated one, contradicting that
    # constraint).
    sorted_nodes = sorted(nodes)
    if len(sorted_nodes) >= 2:
        two_node = [c for n in sorted_nodes[:2] for c in nodes[n]]
        configs.append(('multi-node-%d' % len(two_node), two_node))
    if len(nodes) >= 4:
        sock = [c for n in sorted_nodes[:len(nodes) // 2] for c in nodes[n]]
        configs.append(('socket-%d' % len(sock), sock))

    print('  case    : %s (built fresh, serial)' % a.case)
    case_dir = build_case(a.case)
    n_lo, n_hi = a.steps, a.steps * a.factor
    print('  method  : per-step by difference, %d vs %d steps, fresh process each'
          % (n_lo, n_hi))
    print()

    results = {}
    skipped = {}
    base_state = {}
    # Item 91g: the baseline used to become whichever configuration
    # measured first, silently, if `configs[0]` (the intended baseline) was
    # skipped or failed. `record_baseline` (above) makes that explicit:
    # baseline_label on every result and in the payload, plus a loud NOTE if
    # it is not configs[0] -- never a silent baseline swap.
    print('  %-26s %6s  %12s  %10s  %s'
          % ('configuration', 'cpus', 'ms/step', 'speedup', 'fixed cost'))
    for label, cpus in configs:
        print('  %-26s %6d  checking...' % (label, len(cpus)))
        chk = require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
        if chk is None:
            print('  %-26s %6d  SKIPPED (cpus busy, not overridden)'
                  % (label, len(cpus)))
            skipped[label] = dict(cpus=cpus)
            continue
        load1, load5, load15, others, busy = chk
        ps, fixed = per_step(case_dir, cpus, nodes, n_lo, n_hi)
        if ps is None:
            print('  %-26s %6d  FAILED' % (label, len(cpus)))
            continue
        speedup, base_label, note = record_baseline(base_state, label, ps,
                                                     configs[0][0])
        if note:
            print('  NOTE: ' + note)
        results[label] = dict(cpus=cpus, ms_per_step=ps * 1e3,
                              speedup=speedup, baseline_label=base_label,
                              fixed_s=fixed,
                              load_avg='%.2f %.2f %.2f' % (load1, load5, load15),
                              other_users=others, cpu_busy=busy)
        print('  %-26s %6d  %12.2f  %9.2fx (vs %s)  %8.2f s'
              % (label, len(cpus), ps * 1e3, speedup, base_label, fixed))

    if not results:
        raise SystemExit(
            'FAIL: every configuration was skipped or failed -- nothing '
            'measured. %d skipped as busy: %s. Wait for those cpus to free '
            'up, or pass --i-know-the-box-is-busy.'
            % (len(skipped), sorted(skipped)))
    if skipped:
        print('\n  %d configuration(s) SKIPPED as busy (not measured, not a '
              'silent gap): %s' % (len(skipped), sorted(skipped)))

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
        busy_ceiling=a.busy_ceiling, skipped=skipped,
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
