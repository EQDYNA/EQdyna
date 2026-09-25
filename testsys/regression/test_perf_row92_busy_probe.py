#! /usr/bin/env python3
"""Behavioural regression guard for board row 92 (2026-09-24): the box-busy
probe (free_node_map, per-cpu busy fraction) must live in ONE place --
testsys/perf/busy_probe.py -- and run_scaling.py must import it rather than
keep a second, independently-drifting copy.

No /proc/stat read, no numactl: run_numa_scaling.cpu_busy_fractions is
monkeypatched, so this file is seconds and machine-independent.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, "testsys", "perf")
sys.path.insert(0, PERF)

import busy_probe  # noqa: E402
import run_scaling as rs  # noqa: E402
import run_numa_scaling as numa  # noqa: E402

FAILURES = []


def check(ok, what):
    print("%s -- %s" % ("PASS" if ok else "FAIL", what))
    if not ok:
        FAILURES.append(what)


def check_one_main_function():
    n = sum(1 for line in open(os.path.join(PERF, "busy_probe.py"))
           if line.startswith("def main"))
    check(n == 1, "row92: busy_probe.py defines exactly one main() (got %d)"
                  % n)


def check_run_scaling_imports_not_redefines():
    check(rs.free_node_map is busy_probe.free_node_map,
          "row92: run_scaling.free_node_map IS busy_probe.free_node_map "
          "(same function object, not a second copy)")


def check_free_node_map_filters_by_busy_ceiling():
    """The real free_node_map, driven behaviourally: node 0's cpu 1 is over
    ceiling and must be dropped; node 1 is entirely idle and kept whole."""
    nodes = {0: [0, 1], 1: [2, 3]}
    orig = numa.cpu_busy_fractions
    try:
        numa.cpu_busy_fractions = lambda cpus: {0: 0.05, 1: 0.90,
                                                2: 0.02, 3: 0.03}
        free = busy_probe.free_node_map(nodes, 0.2, False)
    finally:
        numa.cpu_busy_fractions = orig
    check(free == {0: [0], 1: [2, 3]},
          "row92: free_node_map(ceiling=0.2) -> %r (want cpu 1 dropped from "
          "node 0, node 1 kept whole)" % (free,))


def check_override_bypasses_filter():
    nodes = {0: [0, 1]}
    orig = numa.cpu_busy_fractions
    try:
        numa.cpu_busy_fractions = lambda cpus: (_ for _ in ()).throw(
            AssertionError("cpu_busy_fractions must not be called when "
                           "override=True"))
        free = busy_probe.free_node_map(nodes, 0.2, True)
    finally:
        numa.cpu_busy_fractions = orig
    check(free == nodes,
          "row92: free_node_map(override=True) returns the full topology "
          "unfiltered, and never reads /proc/stat at all -> %r" % (free,))


def check_unreadable_utilisation_raises():
    orig = numa.cpu_busy_fractions
    try:
        numa.cpu_busy_fractions = lambda cpus: {}
        try:
            busy_probe.free_node_map({0: [0, 1]}, 0.2, False)
        except SystemExit as e:
            check("cannot tell which cpus are free" in str(e),
                  "row92: free_node_map raises SystemExit naming the "
                  "failure when /proc/stat is unreadable (got: %s)" % e)
        else:
            check(False, "row92: free_node_map did not raise when "
                        "cpu_busy_fractions returned {}")
    finally:
        numa.cpu_busy_fractions = orig


CHECKS = [
    ("one-main", check_one_main_function),
    ("no-duplicate", check_run_scaling_imports_not_redefines),
    ("filters-busy-ceiling", check_free_node_map_filters_by_busy_ceiling),
    ("override-bypasses", check_override_bypasses_filter),
    ("unreadable-raises", check_unreadable_utilisation_raises),
]


def main():
    for tag, fn in CHECKS:
        print("\n-- %s --" % tag)
        fn()
    print()
    if FAILURES:
        print("FAIL test_perf_row92_busy_probe: %d check(s) failed" % len(FAILURES))
        return 1
    print("SUCCESS test_perf_row92_busy_probe: all %d checks passed" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
