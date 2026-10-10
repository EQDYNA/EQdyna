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
                       exists on disk, MINUS the explicit
                       regression_sweep_exclusions.EXCLUDED_FROM_SWEEP set
                       (see below), is in EXACTLY ONE shard, and no excluded
                       script is assigned to any shard at all. Prints the
                       counts it compared (files on disk, files excluded,
                       files assigned, files duplicated/missing/stale) --
                       never a bare pass/fail (this repo's own
                       papercuts.md: "a green result that tested nothing").

EXCLUDED_FROM_SWEEP (2026-09-30, owner requirement (c)):
  `testsys/regression_sweep_exclusions.py`'s `EXCLUDED_FROM_SWEEP` names the
  regression script(s) that are release-PROCESS checks, not per-commit
  regression guards (currently just `test_release_complete.py` -- see that
  module's own docstring for the 2026-09-30 incident this closes). Sharding
  this file into a CI shard would put it back in a job that runs on every
  commit and every PR, which is exactly the failure mode being fixed:
  `test_release_complete.py`'s assertions are meaningful only at the moment
  of auditing a release ALREADY tagged, and permanently fail between
  releases. `verify()` therefore checks it is in NO shard, and
  `discover_regression_scripts()` still reports it as present on disk (it
  is a real file) so a human reading `verify`'s printed counts sees it
  accounted for, not silently vanished.

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

# Single-file import, no transitive dependency on run.py/matrix.py/run_e2e.py
# -- see that module's own docstring and this file's module docstring above.
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
from testsys import regression_sweep_exclusions  # noqa: E402

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
        "test_ci_pr_policy_step.py",  # added 2026-09-23 (Iris, pr-enforce): parses test.yml, one subprocess, <1 s
        "test_ci_fast_lane_unit_regression.py",  # added 2026-10-07 (fast lane runs unit-regression): parses test.yml twice, <1 s
        "test_ci_push_trigger_filter.py",  # added 2026-09-24 (wei-lin, row 87): parses test.yml, <1 s
        "test_pr_policy_guard.py",  # added 2026-09-23 (Iris, pr-enforce): scratch git repo, no network, <1 s
        "test_prepush_pr_policy_guard.py",  # added 2026-09-23 (Iris, pr-enforce): bare remote + real pushes, ~1-2 s
        "test_dipping_fault_y_split.py",
        "test_multifault_item10_fix.py",
        "test_multifault_item9_fix.py",
        "test_multifault_refused.py",
        "test_multifault_no_ntotft1_specialcase.py",  # added row 17: pure text scan, sub-second
        "test_tpv2223_multifault_routing.py",  # replaces test_multifault_two_fault_smoke.py (item 17 section A, test.multifault2 retired): serial (1 rank, both faults) Fortran routing check on test.tpv22's own real compset, ~2.3 min
        "test_perf_parallelism_discriminator.py",
        "test_perf_item91_guards.py",  # added 2026-09-24 (Iris, item 91 perf-
                                 # tool defects): pure monkeypatched
                                 # behavioural checks, no MPI/subprocess/box
                                 # dependence, well under 1 s.
        "test_perf_row76_contention.py",  # added 2026-09-24 (Iris, row 76
                                 # contention field + row 112.ii tree_dirty_
                                 # once): pure in-memory ledger checks, no
                                 # /proc/stat/numactl/git subprocess, <1 s.
        "test_perf_row92_busy_probe.py",  # added 2026-09-24 (Iris, row 92
                                 # busy_probe extraction): monkeypatched
                                 # cpu_busy_fractions, no real /proc/stat
                                 # read, <1 s.
        "test_perf_row112_e2e_ranks.py",  # added 2026-09-24 (Iris, row 112.i
                                 # cell_cost/profile_ranks non-positive-rank
                                 # guard): pure monkeypatched matrix table,
                                 # no case build/solver run, <1 s.
        "test_release_evidence_tree_clean.py",  # added 2026-09-23 (Iris, rule-24 tree_clean fix): 2 sandbox git-init scenarios, well under 1 s
        "test_content_key_sweep_evidence.py",  # added 2026-09-30 (PR #63): content-key evidence, 6 sandbox scenarios + mutation, ~1 s
        "test_evidence_tpv12_both_ways.py",  # added 2026-10-07 (victor-reyes
                                 # audit, PR #147): 3 subprocess invocations
                                 # of testsys/parity/evidence_tpv12_scec_
                                 # comparison.py against already-committed
                                 # station text (real data, a scaled-column
                                 # copy, a copy missing frt.canonical.txt),
                                 # no network/build, well under a few seconds.
        "test_evidence_tpv13_both_ways.py",  # added 2026-10-07 (row 149
        "test_tpv12_tpv13_border_rupture_eligible.py",  # added 2026-10-09 (defect B): imports case params, no solver, <1 s
                                 # follow-on): same pattern as test_evidence_
                                 # tpv12_both_ways.py, 4 subprocess invocations
                                 # of testsys/parity/evidence_tpv13_scec_
                                 # comparison.py against already-committed
                                 # station/frt text (real data, a scaled-
                                 # column copy, a copy missing frt.canonical.
                                 # txt, a collapsed-rupture fraction-bound
                                 # fixture), no network/build, well under a
                                 # few seconds.
        "test_profile_guard.py",  # added 2026-09-23 (item 3, profile-guard
                                 # testsys half); ~0.03 s, pure fixture/schema
                                 # checks, no subprocess -- negligible to
                                 # shard 1's timed 41.5 s total.
        "test_hpc_scaling_suite.py",  # added 2026-10-06 (item 142, HPC scaling
                                 # suite): pure in-memory collect/analyze
                                 # arithmetic against a committed profile_guard
                                 # fixture, no subprocess/sbatch/solver launch,
                                 # <1 s -- negligible to shard 1's timed total.
        "test_profile_env_strict.py",  # added 2026-09-23 (profile-fix audit,
                                 # item 3 follow-up): EQDYNA_PROFILE strict
                                 # parse, both languages; ~1 s (one mpirun -np
                                 # 1 launch that aborts at env-parse, before
                                 # any input file) -- negligible to shard 1.
        "test_profile_ranks_helper.py",  # added 2026-09-23 (Iris, combo-fix
                                 # item 3 follow-up): profile_ranks() vs
                                 # cell_cost() conflation guard; pure-Python
                                 # monkeypatched calls, no subprocess,
                                 # negligible to shard 1's timed 41.5 s total.
        "test_rough_fault_normal_consistency.py",
        "test_rsfNucleation_tpv2802_td.py",
        "test_station_header_column_count.py",
        "test_int64_index_width.py",  # added 2026-10-03 (item 143, int64 indices): 3 tiny gfortran builds of globalvar.f90 + a driver, <3 s
        "test_offfault_station_dropped_report.py",  # added 2026-09-24 (wei-lin, item 94): 2 tiny gfortran builds, <2 s
        "test_onfault_station_dropped_report.py",  # added 2026-09-24 (wei-lin, item 116): 2 tiny gfortran builds, <2 s
        "test_offfault_station_header_actual_node.py",  # added 2026-09-25 (mira-volkov, row 94 audit finding 4): 1 tiny gfortran build, <1 s
        "test_offfault_station_header_actual_node_python.py",  # added 2026-09-25 (mira-volkov, row 94 audit finding 4): pure Python, <1 s
        "test_stop_exit_status.py",
        "test_symlink_integrity.py",
        "test_shared_compset_file_loud.py",  # added 2026-09-24 (wei-lin, item 44): 3 create.newcase/case.setup subprocesses, ~3 s
        "test_term_axis.py",
        "test_station_gate.py",  # added 2026-09-24 (wei-lin, owner gate design): tempdir copies of committed refs, 10 scenarios, ~1 s
        "test_nstress_sign_convention.py",  # added 2026-09-24 (wei-lin, row 22a): 11 case-param subprocesses + tempdir gate checks, ~1 s
        "test_threadprobe.py",  # added 2026-09-24 (wei-lin, item 47a): pure function, <1 s
        "test_no_hardcoded_paths.py",  # added 2026-10-07 (path scrub, owner decision 2026-10-06): one git grep subprocess, <1 s
        "test_defect_c_depth_offset_removed.py",  # added 2026-10-09 (mira-volkov,
                                 # board PR #164, Defect C): pure regex/text
                                 # checks over already-committed meshgen.f90/
                                 # meshgen.py, plus 2 in-memory reverted-text
                                 # fixtures (rule 14a both-ways); no
                                 # subprocess/build, <1 s.
    ],
    "2": [
        "test_check_input_consistency.py",  # added 2026-09-24 (mira-volkov,
                                 # checkInputConsistency port): one build +
                                 # 2 tiny 1-rank runs each on Fortran and
                                 # python-numpy, well under shard 2's other
                                 # subprocess-heavy scripts.
        "test_material_grid3d_coverage.py",  # added 2026-10-06 (mira-volkov,
                                 # TPV34 PR #93 pre-merge audit Medium): one
                                 # build + 2 tiny 1-rank runs each on Fortran
                                 # and python-numpy, same cost profile as
                                 # test_check_input_consistency.py above.
        "test_offfault_station_depth_selection.py",  # added 2026-09-25
                                 # (mira-volkov, row 94 audit finding 3): one
                                 # serial case.setup + one 3-step 1-rank
                                 # mpirun eqdyna run + one direct Python
                                 # build_station_matching call on the same
                                 # case dir -- comparable cost to
                                 # test_check_input_consistency.py above.
        "test_ci_board_separation_step.py",
        "test_dev_str_depth_taper.py",
        "test_e2e_run_tree_lock.py",
        "test_e2e_python_cell_logs_kept.py",  # added 2026-09-24 (wei-lin, item 109): one tee child + one failed jax import, ~3 s
        "test_e2e_cell_timeout.py",  # added 2026-10-01 (Iris, per-cell timeout
                                 # feature, the 11h45m test.tpv1053d hang):
                                 # two plain `sleep` children (2s, 5s, the
                                 # second killed at a 1s deadline) plus pure
                                 # function checks against the real
                                 # docs/perf_ledger.jsonl -- no build, no
                                 # solver, measured ~3.5 s total.
        "test_equilibrium_dump.py",
        "test_fault_geometry_guard.py",
        "test_fault_mpi_boundary_arn.py",
        "test_history_table.py",
        "test_make_default_goal.py",
        "test_perf_mpi_placement.py",
        "test_readme_executes.py",  # added 2026-09-25 (Iris, README-executes
                                 # gate): default FAST mode only (parser +
                                 # synthetic mutation self-test, no clone, no
                                 # network); well under 1 s. Its FULL mode
                                 # (real clone + real run, ~2-4 min) only
                                 # fires under EQDYNA_README_GATE=full, set
                                 # by `run.py readme`, never by this shard.
        "test_rank_local_mesh.py",  # added 2026-09-24 (mira-volkov, item 64):
                                 # 3 serial case builds + 6 decompositions'
                                 # rank-local builds in one process, no MPI;
                                 # the heaviest script in this shard.
        "test_row114_station_output.py",  # added 2026-09-24 (mira-volkov,
                                 # row 114 station output): one serial case
                                 # build + a 3-step run on numpy AND jax,
                                 # well under this shard's other
                                 # subprocess-heavy scripts.
        "test_row120_mpi_station_output.py",  # added 2026-09-25 (mira-volkov,
                                 # row 120 python-jax-mpi station output):
                                 # incremental src/fortran build (no-op if
                                 # already built) + one 5-step 4-rank mpirun
                                 # + 4x in-process build_solver_state (no
                                 # mpi/jax on the python side) -- measured
                                 # ~4.9 s total, comparable to
                                 # test_row114_station_output.py above.
        "test_row127_station_ownership.py",  # added 2026-09-25 (mira-volkov,
                                 # row 127 station ownership at a shared
                                 # MPI-partition boundary): shares
                                 # test_row120's incremental src/fortran
                                 # build, runs two 4-rank mpirun cases (5
                                 # steps each, one real tpv8 (2,2,1), one
                                 # synthetic (2,1,2) exercising npz>1) with
                                 # per-rank output directories, plus 8
                                 # in-process build_solver_state calls --
                                 # measured ~9.5 s total.
        "test_row132_axis_too_thin.py",  # added 2026-09-25 (mira-volkov, row
                                 # 132(1) axis-too-thin guard): shares
                                 # test_row127's incremental src/fortran
                                 # build; one ~20-line standalone program
                                 # compile+link (no mpirun) + one Python
                                 # call -- sub-second beyond the shared build.
        "test_row153_degen_range_refused.py",  # added 2026-10-07 (row 153
                                 # checkpoint 1 audit fix): reuses the already-
                                 # built bin/eqdyna, 7 one-rank mpirun launches
                                 # against a 4-line bGlobal.txt (aborts at
                                 # readglobal, before any mesh work) + one
                                 # subprocess launch of case.setup -- well
                                 # under this shard's other mpirun-heavy
                                 # scripts.
        "test_pretag_ci_negative.py",
        "test_readme_commands.py",
        "test_user_docs_commands.py",  # added 2026-09-25 (board row 126, docs
                                 # site): reuses test_readme_commands.py's own
                                 # checks over docs/user/**/*.md, filesystem
                                 # and import checks only, well under 1 s.
        "test_params_reference_freshness.py",  # added 2026-09-25 (board row
                                 # 126, docs site): one AST parse of
                                 # scripts/defaultParameters.py plus a string
                                 # compare, well under 1 s.
        "test_user_docs_style.py",  # added 2026-09-24 (wei-lin, user-facing docs rule): pure text checks, <1 s
        "test_user_docs_coverage.py",  # added 2026-10-03 (owner: "keep docs
                                 # synced", pathway_forward.md section A
                                 # header note): two file reads + a list
                                 # membership check, well under 1 s.
        "test_sweep_core_budget.py",
        "test_sweep_tenancy_budget_2026_09_23.py",
        "test_slot47_peak_sliprate.py",  # added 2026-10-02 (item 17b, PR #76):
                                 # one serial test.tpv8 build+run at a
                                 # shortened 1 s term (vs its 5 s gate term),
                                 # comparable cost to
                                 # test_offfault_station_depth_selection.py
                                 # above.
    ],
    "3": [
        "test_tenancy_without_numa.py",  # added 2026-09-23 (wei-lin, PR #6 CI smoke fix): <2 s
        "test_backend_axis.py",  # added 2026-09-23 (wei-lin, numpy out of the gates): <1 s
        "test_perf_meta_imports.py",  # added 2026-09-23 (wei-lin): <1 s
        "test_sweep_speed_2026_09_23.py",  # added at merge (wei-lin): landed 90514cc after this partition was timed; 0.2 s
        "test_src_stamp.py",  # added 2026-09-24 (wei-lin, source stamp): synthetic binaries + one refused run_e2e, <5 s
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
        "test_pretag_release_docs_negative.py",  # added 2026-09-30 (owner
                                 # requirement (b)/(d)): 5 sandbox git-commit
                                 # scenarios, well under 1 s, same shape as
                                 # test_pretag_sweep_negative.py above.
        "test_publish_image_fetch_depth.py",
        "test_station_header_location_stamp.py",
        "test_stress_i0_carry_aliasing.py",
        "test_version_banner.py",
        "test_root_allowlist.py",  # added 2026-09-24 (wei-lin, root-notes move): one git ls-files, <1 s
        "test_text_line_format.py",  # 2026-09-29: static scan of src/fortran, <1 s
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
    """Partition guard: every script on disk, MINUS the explicit
    regression_sweep_exclusions.EXCLUDED_FROM_SWEEP set, is in EXACTLY ONE
    shard -- and no excluded script is assigned to ANY shard (owner
    requirement (c): a release-process check like test_release_complete.py
    must not be sharded back into a per-commit CI job). Returns (ok: bool,
    message: str). message always states the counts compared, per this
    repo's own papercuts.md rule -- excluded=N printed explicitly, not
    folded silently into missing/stale."""
    on_disk = set(discover_regression_scripts())
    excluded = set(regression_sweep_exclusions.EXCLUDED_FROM_SWEEP)
    swept_on_disk = on_disk - excluded
    assigned_lists = [name for names in SHARDS.values() for name in names]
    assigned = set(assigned_lists)

    problems = []
    if len(assigned_lists) != len(assigned):
        dupes = sorted({n for n in assigned_lists if assigned_lists.count(n) > 1})
        problems.append("duplicated across shards (%d): %s" % (len(dupes), dupes))

    missing = sorted(swept_on_disk - assigned)
    if missing:
        problems.append("on disk (swept) but in NO shard (%d): %s" % (len(missing), missing))

    stale = sorted(assigned - swept_on_disk)
    if stale:
        problems.append("assigned to a shard but no longer on disk, or excluded from the "
                        "sweep (%d): %s" % (len(stale), stale))

    wrongly_assigned = sorted(excluded & assigned)
    if wrongly_assigned:
        problems.append("EXCLUDED_FROM_SWEEP script(s) assigned to a shard anyway (%d): %s "
                        "-- this would re-run a release-process check on every commit's CI"
                        % (len(wrongly_assigned), wrongly_assigned))

    not_on_disk_excluded = sorted(excluded - on_disk)
    if not_on_disk_excluded:
        problems.append("EXCLUDED_FROM_SWEEP names a file no longer on disk (%d): %s -- "
                        "stale exclusion entry" % (len(not_on_disk_excluded), not_on_disk_excluded))

    msg = ("on-disk=%d excluded=%d swept=%d assigned=%d(across %d shards, %d unique) "
           "missing=%d stale=%d"
           % (len(on_disk), len(excluded), len(swept_on_disk), len(assigned_lists),
              len(SHARDS), len(assigned), len(missing), len(stale)))
    if problems:
        return False, msg + " -- " + "; ".join(problems)
    return True, msg + (" -- partition OK (every swept on-disk script in exactly "
                        "one shard, every excluded script in none)")


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
