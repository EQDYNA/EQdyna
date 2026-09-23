#! /usr/bin/env python3
"""
CI sharding for testsys' unit+regression tier (owner request, 2026-09-23:
"make CI faster without weakening any guard").

This file belongs to testsys/regression's OWNER, not to testsys/run.py,
testsys/matrix.py or testsys/e2e/run_e2e.py -- it never imports any of them
and none of them import it. run.py's own 'unit regression' invocation
(local dev, `python3 testsys/run.py unit regression`) is UNCHANGED and stays
the full, unsharded run; sharding is a CI-workflow concern layered on top,
not a new test tier and not a rewrite of how a developer runs the suite
locally.

WHAT THIS IS.
  SHARDS            -- a static partition of testsys/regression/test_*.py's
                       basenames into N groups, hand-balanced against a real
                       timing run on this box (see below for the numbers).
  UNIT_PYTEST_SHARD -- which shard ALSO runs the full `pytest testsys/unit`
                       tier (the unit tier is fast enough as a whole -- ~12 s
                       -- that splitting it further than "one shard runs all
                       of it" was not worth a second pytest invocation).
  `run <shard>`     -- runs that shard's regression scripts (each a
                       standalone script with its own SUCCESS/FAIL banner,
                       same contract as run.py's run_regression) plus the
                       unit tier if this is UNIT_PYTEST_SHARD, and returns
                       the OR of every exit code.
  `verify`          -- the partition guard: every test_*.py that actually
                       exists on disk is in EXACTLY ONE shard. Prints the
                       counts it compared (files on disk, files assigned,
                       files duplicated/missing/stale) -- never a bare
                       pass/fail (this repo's own papercuts.md: "a green
                       result that tested nothing").

WHY STATIC RATHER THAN AUTO-BALANCED AT RUNTIME. A shard assignment that
silently re-balances itself when a script's runtime drifts would also
silently swallow a NEW script nobody added to any shard (it would just land
wherever the balancer put it, indistinguishable from a script placed there
deliberately) -- exactly the failure mode `verify` exists to catch. A static
list makes "this script is not in any shard" a loud, mechanical fact instead
of an emergent property of a scheduler. Rebalancing is a human edit to
SHARDS, same as adding a new regression script is.

test_ci_shard_coverage.py (testsys/regression/) is `verify` wired into the
regression tier itself, so the guard runs on every `run.py regression`
invocation and inside whichever shard it is assigned to -- not just when a
human remembers to run this file's own CLI.
"""
import glob
import os
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
REGRESSION_DIR = os.path.join(TESTSYS, "regression")

# Hand-balanced against a real timed run of every regression script on this
# box (numactl-pinned, mpirun cells at 4 ranks), 2026-09-23. Totals: shard 1
# 41.5 s, shard 2 41.5 s, shard 3 41.4 s (regression only; shard 3 also carries
# the ~12.4 s unit pytest tier -- see UNIT_PYTEST_SHARD). Rebalance by hand if
# a future change shifts one script's cost substantially; `verify` will not
# catch an IMBALANCED-but-complete partition, only an incomplete one -- that
# is a deliberate scope split, not an oversight.
SHARDS = {
    "1": [
        "test_ci_dependencies.py",
        "test_dipping_fault_y_split.py",
        "test_multifault_item10_fix.py",
        "test_multifault_item9_fix.py",
        "test_multifault_refused.py",
        "test_perf_parallelism_discriminator.py",
        "test_rough_fault_normal_consistency.py",
        "test_rsfNucleation_tpv2802_td.py",
        "test_station_header_column_count.py",
        "test_stop_exit_status.py",
        "test_symlink_integrity.py",
        "test_term_axis.py",
    ],
    "2": [
        "test_ci_board_separation_step.py",
        "test_dev_str_depth_taper.py",
        "test_e2e_run_tree_lock.py",
        "test_equilibrium_dump.py",
        "test_fault_geometry_guard.py",
        "test_fault_mpi_boundary_arn.py",
        "test_history_table.py",
        "test_make_default_goal.py",
        "test_perf_mpi_placement.py",
        "test_pretag_ci_negative.py",
        "test_readme_commands.py",
        "test_sweep_core_budget.py",
        "test_sweep_tenancy_budget_2026_09_23.py",
    ],
    "3": [
        "test_sweep_speed_2026_09_23.py",  # added at merge (wei-lin): landed d488dae after this partition was timed; 0.2 s
        "test_ci_shard_coverage.py",
        "test_ci_workflow_coverage.py",
        "test_create_newcase.py",
        "test_docker_guide_no_pinned_version.py",
        "test_drucker_prager_kernel.py",
        "test_fractal_fault_geometry_derivatives.py",
        "test_multifault_item7_fix.py",
        "test_perf_ledger.py",
        "test_perf_tool_locks.py",
        "test_pml_region_axes.py",
        "test_precommit_board_separation_guard.py",
        "test_precommit_main_checkout_guard.py",
        "test_pretag_sweep_negative.py",
        "test_publish_image_fetch_depth.py",
        "test_release_complete.py",
        "test_station_header_location_stamp.py",
        "test_stress_i0_carry_aliasing.py",
        "test_version_banner.py",
    ],
}
UNIT_PYTEST_SHARD = "3"


def discover_regression_scripts():
    """Same discovery rule as run.py's run_regression: every test_*.py
    directly under testsys/regression/, sorted. Deliberately reimplemented
    (not imported from run.py) -- this file must stay import-independent of
    run.py/matrix.py/run_e2e.py."""
    return sorted(
        os.path.basename(p)
        for p in glob.glob(os.path.join(REGRESSION_DIR, "test_*.py"))
    )


def verify():
    """Partition guard: every script on disk is in EXACTLY ONE shard.
    Returns (ok: bool, message: str). message always states the counts
    compared, per this repo's own papercuts.md rule."""
    on_disk = set(discover_regression_scripts())
    assigned_lists = [name for names in SHARDS.values() for name in names]
    assigned = set(assigned_lists)

    problems = []
    if len(assigned_lists) != len(assigned):
        dupes = sorted({n for n in assigned_lists if assigned_lists.count(n) > 1})
        problems.append("duplicated across shards (%d): %s" % (len(dupes), dupes))

    missing = sorted(on_disk - assigned)
    if missing:
        problems.append("on disk but in NO shard (%d): %s" % (len(missing), missing))

    stale = sorted(assigned - on_disk)
    if stale:
        problems.append("assigned to a shard but no longer on disk (%d): %s" % (len(stale), stale))

    msg = ("on-disk=%d assigned=%d(across %d shards, %d unique) missing=%d stale=%d"
           % (len(on_disk), len(assigned_lists), len(SHARDS), len(assigned),
              len(missing), len(stale)))
    if problems:
        return False, msg + " -- " + "; ".join(problems)
    return True, msg + " -- partition OK (every on-disk script in exactly one shard)"


def run_shard(shard_id):
    if shard_id not in SHARDS:
        print("ci_shard: FAIL - unknown shard %r (known: %s)" % (shard_id, sorted(SHARDS)))
        return 1
    print("\n==== testsys: ci_shard %s ====" % shard_id)
    overall = 0
    for name in SHARDS[shard_id]:
        path = os.path.join(REGRESSION_DIR, name)
        print("-- %s --" % name)
        rc = subprocess.call([sys.executable, path], cwd=REPO_ROOT)
        print("%s %s (exit %d)" % ("SUCCESS" if rc == 0 else "FAIL", name, rc))
        overall = overall or rc
    if shard_id == UNIT_PYTEST_SHARD:
        print("\n==== testsys: unit (carried by shard %s) ====" % shard_id)
        rc = subprocess.call(
            [sys.executable, "-m", "pytest", "-v", os.path.join(TESTSYS, "unit")],
            cwd=REPO_ROOT,
        )
        overall = overall or rc
    return overall


def main(argv):
    if len(argv) >= 1 and argv[0] == "verify":
        ok, msg = verify()
        print(("PASS" if ok else "FAIL") + " ci_shard verify: " + msg)
        return 0 if ok else 1
    if len(argv) >= 2 and argv[0] == "run":
        return run_shard(argv[1])
    print("usage: ci_shard.py verify | ci_shard.py run <shard-id>")
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
