#! /usr/bin/env python3
"""Box-state busy probe: per-cpu busy fractions and a free-cpu-by-NUMA-node
map, extracted from testsys/perf/run_scaling.py's `free_node_map` (row 92's
probe half, 2026-09-24) so every perf tool that needs to know which cpus are
free reads ONE implementation, not a second copy that can drift from it
(PROJECT_RULES rule 1). `run_scaling.py` now imports `free_node_map` from
here rather than defining it a second time.

Reuses `run_numa_scaling.numa_topology`/`cpu_busy_fractions` -- the SAME
two-read /proc/stat idle-time delta every perf tool already measures
busy-ness with -- rather than a third sampling method. `ledger.py`'s
`contention` field (row 76) is built from each tool's own pre-flight busy
check (this sampler, already called), not from a fresh call here -- see
`ledger.contention_from_check`.

CLI: `python3 testsys/perf/busy_probe.py [--busy-ceiling 0.2]` prints the
live per-cpu busy map and a per-NUMA-node free/total summary, for a human
checking box state before launching a sweep.
"""
import argparse
import os
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TESTSYS)
import run_numa_scaling as numa   # noqa: E402


def free_node_map(all_nodes, busy_ceiling, override):
    """{node: cpus} restricted to the individual cpus currently under
    `busy_ceiling`, so `compact_cpus`/`spread_cpus` (run_scaling.py) build
    placements out of room that actually exists right now instead of always
    starting at node 0. Probes fresh on every call -- occupancy moves over
    the course of a multi-minute sweep, so a configuration built early from a
    stale free-set could target a cpu that has since gone busy (the per-cpu
    `require_idle` check downstream still catches that case, but asking for
    the wrong cpu in the first place is the defect this function fixes).

    PER-CPU, NOT WHOLE-NODE (2026-09-18 fix, item 33 thread-scaling
    investigation). The original version required EVERY cpu on a node to be
    idle before offering ANY of that node's cpus -- correct for a box where
    interference clusters by node, wrong for a box with foreign single-core
    jobs smeared across every node, where no node is ever seen fully idle in
    a snapshot even though most cores are free at any instant. Cpu-level
    filtering keeps `compact_cpus` filling node-by-node in cpu-id order
    exactly as before (so it stays as compact as the actually-free cpus
    allow) but no longer demands more idle room than a configuration
    actually needs.

    `override` mirrors `--i-know-the-box-is-busy`: if the busy ceiling itself
    is being bypassed, there is nothing to filter FOR, so this reverts to the
    full topology in node-number order -- identical to pre-fix behaviour.

    Raises SystemExit if per-cpu utilisation could not be read AT ALL (same
    hard-failure discipline as `cpu_busy_fractions`/`require_idle`: a check
    that cannot evaluate must fail, not silently call every cpu free). An
    individual cpu whose utilisation could not be read is dropped from the
    free set (not assumed free, not used to invalidate cpus that WERE
    read)."""
    if override:
        return dict(all_nodes)
    all_cpus = sorted(c for cs in all_nodes.values() for c in cs)
    busy = numa.cpu_busy_fractions(all_cpus)
    if not busy:
        raise SystemExit(
            'FAIL: could not read per-cpu utilisation for cpus %s from '
            '/proc/stat -- cannot tell which cpus are free (rule 2: a '
            'check that cannot evaluate must fail).' % all_cpus)
    free = {}
    for n, cpus in all_nodes.items():
        idle_cpus = [c for c in cpus
                    if busy.get(c) is not None and busy[c] <= busy_ceiling]
        if idle_cpus:
            free[n] = idle_cpus
    return free


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    a = ap.parse_args()
    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology; cannot '
                         'probe box state.')
    all_cpus = sorted(c for cs in nodes.values() for c in cs)
    busy = numa.cpu_busy_fractions(all_cpus)
    if not busy:
        raise SystemExit('FAIL: could not read per-cpu utilisation from '
                         '/proc/stat.')
    print('per-cpu busy fraction (ceiling %.2f):' % a.busy_ceiling)
    for c in all_cpus:
        b = busy.get(c)
        flag = 'free' if (b is not None and b <= a.busy_ceiling) else 'BUSY'
        print('  cpu %3d  busy=%-6s  %s'
              % (c, ('%.3f' % b) if b is not None else 'unreadable', flag))
    print('\nper-NUMA-node summary:')
    for n in sorted(nodes):
        cpus = nodes[n]
        free = [c for c in cpus
               if busy.get(c) is not None and busy[c] <= a.busy_ceiling]
        print('  node %2d: %d/%d cpus free' % (n, len(free), len(cpus)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
