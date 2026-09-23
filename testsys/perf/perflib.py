#! /usr/bin/env python3
"""Shared helpers for the tools that rebuild a FIXED in-repo case directory.

WHY THIS FILE EXISTS. Four builders --

    run_perf.build_perf_case                (testsys/perf/perf_case, and
                                             perf_case_tpv29 when
                                             run_tpv29_pinned_compare.py
                                             repoints PERF_CASE)
    run_scaling.build_py_case               (testsys/perf/scaling_case)
    run_numa_scaling.build_case             (testsys/perf/numa_case)
    probe_plastic_traction.build_case       (testsys/parity/probe_case)

-- carried the SAME body in four copies: derive the lock resource from the
directory about to be destroyed, acquire it non-blocking, turn a refusal into
`SystemExit('FAIL: ...')`, import the e2e sweep's own case builder, rmtree,
makedirs, make_serial_case. That duplication is not incidental to this repo's
history: the unlocked-rmtree defect had to be found and repaired FOUR separate
times (pathway items 70, 74, 77, 81), once per copy, because a fix to one copy
could not reach the others. Two functions here, one call each, is the state in
which that defect has one place to be wrong.

WHAT IS DELIBERATELY NOT SHARED. The LOCK LIFETIME differs per tool and is a
real difference, not an accident to be flattened into a flag:

  * run_perf, run_scaling, run_numa_scaling KEEP the lock after the builder
    returns -- the caller then TIMES out of the directory it just built, and a
    lock released at the end of the build would let a second invocation rebuild
    the tree mid-measurement, which is the item-74/77 collision itself.
  * run_scaling and run_numa_scaling additionally memoise their lock in a
    module global, because their builder has several callers and flock is
    per open file description: a second acquire from the same process would
    refuse against its OWN pid, a false collision.
  * probe_plastic_traction RELEASES as soon as the write is done -- its
    documented contract, since its measurement phase never touches the tree
    destructively.

So each tool keeps those three lines at its own call site, where the lifetime
is visible, and only the two steps that were genuinely identical live here.
`testsys/regression/test_perf_tool_locks.py` observes all of it by RUNNING each
builder with its destroyers replaced by recorders, including which builders
still hold the lock when they return.

THE ROOT IS THIS CHECKOUT'S, NEVER $EQDYNAROOT (item 77). `run_scaling.ROOT`
honours $EQDYNAROOT while the case directory it rebuilds is derived from
__file__; a lock rooted in $EQDYNAROOT would guard one tree and destroy
another. REPO_ROOT below is derived from this file's own location for that
reason, and all four callers already derived theirs the same way.
"""
import os
import shutil
import sys

PERF_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(PERF_DIR))
E2E_DIR = os.path.join(REPO_ROOT, 'testsys', 'e2e')

sys.path.insert(0, REPO_ROOT)
from testsys import runlock  # noqa: E402


def acquire_case_lock(case_dir, consequence):
    """GATE 0: take the exclusive lock on the directory HOLDING `case_dir`,
    before anything is imported, created or deleted.

    The resource is derived from the directory actually about to be destroyed,
    so the lock follows the path rather than restating it -- which is what
    keeps `run_perf.build_perf_case` guarded when `run_tpv29_pinned_compare.py`
    repoints PERF_CASE at perf_case_tpv29/.

    `consequence` is the refusal's closing paragraph: what THIS tool's
    collision costs. It stays with the tool, since the cost does.

    Returns the held `runlock.RunTreeLock` -- the caller decides whether to
    keep it (it is still held if the reference is dropped: the fd stays open
    and release is registered atexit) or to release it. Raises SystemExit on
    refusal, non-blocking, never waiting and never falling back (rule 2).
    """
    resource = os.path.relpath(os.path.dirname(case_dir), REPO_ROOT)
    try:
        return runlock.acquire(REPO_ROOT, resource, consequence=consequence)
    except runlock.RunTreeLocked as exc:
        raise SystemExit('FAIL: %s' % exc)


def rebuild_serial_case(case_name, case_dir):
    """rmtree `case_dir` and rebuild it: create.newcase + forced serial +
    case.setup, through the e2e sweep's OWN builder (rule 1) rather than a
    second copy of case setup that can drift from scripts/.

    THE CALLER MUST ALREADY HOLD THE LOCK on the parent directory
    (`acquire_case_lock`). This function destroys; it does not decide when it
    is safe to.
    """
    sys.path.insert(0, E2E_DIR)
    import run_e2e  # noqa: E402
    if os.path.isdir(case_dir):
        shutil.rmtree(case_dir)
    os.makedirs(os.path.dirname(case_dir), exist_ok=True)
    run_e2e.make_serial_case(case_name, case_dir, run_e2e.base_env())
    return case_dir
