#! /usr/bin/env python3
"""
Regression guard: every regression script is assigned to EXACTLY ONE CI
shard (testsys/ci_shard.py), never zero and never more than one.

Background (owner request, 2026-09-23): unit-regression was CI's critical
path (3.4 of a 4.1 min run). It is now split across 3 parallel CI jobs
(.github/workflows/test.yml), each running `python3 testsys/ci_shard.py run
<N>`. A script silently dropped from every shard -- e.g. a future
test_*.py added to testsys/regression/ without a matching edit to
ci_shard.SHARDS -- would still show green in CI (nothing FAILED; it simply
never RAN), which is exactly the "green result that tested nothing" shape
this repo's own papercuts.md warns about. This file is that guard, wired
into the regression tier itself (test_ci_shard.SHARDS's own module) so it
runs on every `run.py regression` invocation, in local dev AND inside
whichever shard it lands in -- it does not depend on a human remembering to
run `ci_shard.py verify` separately.

WHAT THIS PINS:
  1. every test_*.py that exists on disk under testsys/regression/ is
     assigned to at least one shard (no orphans);
  2. no script is assigned to more than one shard (no double-run, which
     would silently widen the critical path back out and mask a shard
     that is secretly missing something else);
  3. no shard names a script that is not actually on disk (a stale entry
     left behind by a rename/delete, quietly making that shard shorter
     than it looks).

This is intentionally NOT a check on shard balance (rule: scope split, not
oversight -- see ci_shard.py's own module docstring). An imbalanced-but-
complete partition is a performance issue for a human to rebalance, not a
correctness issue this guard exists to catch.

Cheap (rule 9): pure Python, one glob, one dict walk. Milliseconds.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from testsys import ci_shard


def main():
    print("Regression guard: every regression script is in exactly one CI shard")
    ok, msg = ci_shard.verify()
    print(("  PASS  " if ok else "  FAIL  ") + msg)
    if not ok:
        print("\nFAIL test_ci_shard_coverage")
        return 1
    print("\nSUCCESS test_ci_shard_coverage")
    return 0


if __name__ == "__main__":
    sys.exit(main())
