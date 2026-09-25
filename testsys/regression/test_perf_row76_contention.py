#! /usr/bin/env python3
"""Behavioural regression guards for board row 76 (2026-09-24): every NEW
docs/perf_ledger.jsonl row must carry a contention field -- the busy
fraction of the SPECIFIC cpus a measurement used, complementing the
existing whole-box tenancy_busy/tenancy_total (see ledger.py's module
docstring for the field's precise definition and sampling window).

Every check drives the real ledger functions with synthetic, in-memory
rows -- no /proc/stat read, no numactl, no git subprocess -- so this file is
seconds and machine-independent (same discipline as
test_perf_item91_guards.py).

Also covers item 6ii's ledger.tree_dirty_once: memoized per process, so a
tool that captures many times per invocation reads the tree state ONCE.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
PERF = os.path.join(ROOT, "testsys", "perf")
sys.path.insert(0, PERF)

import ledger  # noqa: E402
import run_numa_scaling as numa  # noqa: E402

FAILURES = []


def check(ok, what):
    print("%s -- %s" % ("PASS" if ok else "FAIL", what))
    if not ok:
        FAILURES.append(what)


def raises(fn, exc_types, what, must_contain=()):
    try:
        fn()
    except exc_types as e:
        msg = str(e)
        absent = [f for f in must_contain if f not in msg]
        if absent:
            check(False, "%s -> %s but message lacks %s: %s"
                  % (what, type(e).__name__, absent, msg[:220]))
        else:
            check(True, "%s -> %s; message: %s"
                  % (what, type(e).__name__, msg[:220].replace("\n", " ")))
        return
    check(False, "%s -> nothing raised" % what)


def _row(**overrides):
    row = dict(ts_utc="2026-09-24T00:00:00Z", snapshot_date_local="2026-09-24",
              sha="abc1234", host="h", tool="run_scaling", case="test.tpv104",
              backend="fortran", ranks=1, ms_per_step=10.0,
              rank_ms_min=None, rank_ms_max=None, rank_ms_mean=None,
              effective_cores=None, threads_per_rank=None,
              busy_ceiling=0.2, tenancy_busy=0, tenancy_total=4,
              metric="per-step-by-difference", n_lo=20, n_hi=60,
              snapshot="docs/perf_snapshots/fake.json",
              platform="cpu", platform_evidence="fortran: CPU-only",
              devices=None, parallelism="mpi",
              contention=dict(cpus=[0, 1], busy=0, total=2))
    row.update(overrides)
    return row


def check_contention_required_on_append():
    row = _row()
    del row["contention"]
    raises(lambda: ledger.validate(row, appending=True), (ValueError,),
           "row76: validate(appending=True) on a row with NO contention key",
           must_contain=("contention",))
    ledger.validate(row, appending=False)
    check(True, "row76: validate(appending=False) tolerates a legacy row "
                "with no contention key")


def check_contention_shape_enforced():
    bad_cases = [
        (dict(cpus=[], busy=0, total=0), "empty cpus list"),
        (dict(cpus=[0, 1], busy=3, total=2), "busy > total"),
        (dict(cpus=[0, 1], busy=0, total=3), "total != len(cpus)"),
        ("not-a-dict", "contention is not a dict"),
    ]
    for bad, what in bad_cases:
        raises(lambda bad=bad: ledger.validate(_row(contention=bad),
                                               appending=True),
               (ValueError,), "row76: contention=%r (%s) refused" % (bad, what))
    ledger.validate(_row(), appending=True)
    check(True, "row76: a well-shaped contention field passes "
                "validate(appending=True)")


def check_contention_from_check_two_shapes():
    five_tuple = (0.1, 0.2, 0.3, [], {0: 0.05, 1: 0.90})
    c1 = ledger.contention_from_check(five_tuple, 0.2)
    check(c1 == dict(cpus=[0, 1], busy=1, total=2),
          "row76: 5-tuple shape -> %r (want cpus=[0,1] busy=1 total=2)" % (c1,))

    bare_dict = {"2": 0.05, "3": 0.10}
    c2 = ledger.contention_from_check(bare_dict, 0.2)
    check(c2 == dict(cpus=[2, 3], busy=0, total=2),
          "row76: bare str-keyed dict shape -> %r (want cpus=[2,3] busy=0 "
          "total=2)" % (c2,))

    raises(lambda: ledger.contention_from_check({}, 0.2), (ValueError,),
           "row76: contention_from_check({}, ...) refuses an empty busy check")
    raises(lambda: ledger.contention_from_check("nonsense", 0.2),
           (ValueError,),
           "row76: contention_from_check on an uninterpretable shape refuses")


def check_box_tenancy_carries_cpus_and_contention_from_tenancy():
    orig_topo = numa.numa_topology
    orig_busy = numa.cpu_busy_fractions
    try:
        numa.numa_topology = lambda: {0: [0, 1], 1: [2, 3]}
        numa.cpu_busy_fractions = lambda cpus: {
            c: (0.9 if c == 0 else 0.05) for c in cpus}
        t = ledger.box_tenancy(0.2)
    finally:
        numa.numa_topology = orig_topo
        numa.cpu_busy_fractions = orig_busy
    check(t["cpus"] == [0, 1, 2, 3] and t["busy"] == 1 and t["total"] == 4,
          "row76: box_tenancy(0.2) on a synthetic 4-cpu box -> %r "
          "(want cpus=[0,1,2,3] busy=1 total=4)" % (t,))
    c = ledger.contention_from_tenancy(t)
    check(c == dict(cpus=[0, 1, 2, 3], busy=1, total=4),
          "row76: contention_from_tenancy(box_tenancy(...)) relabels the "
          "SAME reading -> %r" % (c,))
    raises(lambda: ledger.contention_from_tenancy(dict(busy=0, total=1)),
           (ValueError,),
           "row76: contention_from_tenancy on an old tenancy dict with no "
           "cpus key refuses rather than guessing")


def check_run_jaxmpi_ab_wires_contention():
    src = open(os.path.join(PERF, "run_jaxmpi_ab.py")).read()
    check("contention=ledger.contention_from_check" in src,
          "row76: run_jaxmpi_ab.py builds its ledger rows with a "
          "contention= field (it hand-builds rows outside "
          "rows_from_mpi_scaling_snapshot, so the universal check is the "
          "only thing that would have caught a miss here)")


def check_tree_dirty_once_memoizes():
    calls = []

    def fake(paths=("src", "testsys"), root=ledger.ROOT):
        calls.append(1)
        return len(calls) == 1

    orig = ledger.tree_dirty
    orig_cache = dict(ledger._tree_dirty_cache)
    try:
        ledger.tree_dirty = fake
        ledger._tree_dirty_cache.clear()
        r1 = ledger.tree_dirty_once()
        r2 = ledger.tree_dirty_once()
        r3 = ledger.tree_dirty_once()
    finally:
        ledger.tree_dirty = orig
        ledger._tree_dirty_cache.clear()
        ledger._tree_dirty_cache.update(orig_cache)
    check(len(calls) == 1,
          "row112.ii: tree_dirty_once() called 3x invoked the real "
          "tree_dirty() %d time(s), want exactly 1" % len(calls))
    check(r1 is True and r2 is True and r3 is True,
          "row112.ii: every call returns the FIRST call's value (%r,%r,%r)"
          % (r1, r2, r3))


def check_capture_sites_use_tree_dirty_once():
    for fname in ("run_perf.py", "run_scaling.py", "run_mpi_scaling.py"):
        src = open(os.path.join(PERF, fname)).read()
        check("profile_record.ledger.tree_dirty()" not in src,
              "row112.ii: %s no longer calls the raw tree_dirty() directly"
              % fname)
        check("profile_record.ledger.tree_dirty_once()" in src,
              "row112.ii: %s calls the memoized tree_dirty_once()" % fname)


CHECKS = [
    ("contention-required", check_contention_required_on_append),
    ("contention-shape", check_contention_shape_enforced),
    ("contention-from-check", check_contention_from_check_two_shapes),
    ("box-tenancy-cpus", check_box_tenancy_carries_cpus_and_contention_from_tenancy),
    ("jaxmpi-ab-wired", check_run_jaxmpi_ab_wires_contention),
    ("tree-dirty-once", check_tree_dirty_once_memoizes),
    ("tree-dirty-call-sites", check_capture_sites_use_tree_dirty_once),
]


def main():
    for tag, fn in CHECKS:
        print("\n-- %s --" % tag)
        fn()
    print()
    if FAILURES:
        print("FAIL test_perf_row76_contention: %d check(s) failed" % len(FAILURES))
        return 1
    print("SUCCESS test_perf_row76_contention: all %d checks passed" % len(CHECKS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
