# Item 77 -- the three remaining perf-tool holes (iris, 2026-09-23)

Branch `iris/item77-perf-locks`, worktree
`/home/utig5/dliu/EQdyna/.claude/worktrees/iris-item77`, from `4ee171b`.
Finishes item 74. Nothing here runs a perf measurement or launches `mpirun`.

## 1. `run_scaling.build_py_case` -> lock on `testsys/perf/scaling_case`

Acquire sits INSIDE `build_py_case`, first statement after `d` is formed, and
before the `run_e2e` import, the `rmtree` and `make_serial_case`
(`testsys/perf/run_scaling.py:428-447`). Not in `main()`: main() is the path
this builder is LEAST often reached by. Five tools call it --
`run_scaling.main`, `run_mpi_scaling.py:408`, `run_jaxmpi_ab.py:169`,
`run_setup_probe.py:142`, `run_shard_scaling.py:195` -- and the last four call
`rs.build_py_case` directly and never reach this module's `main()`. A lock in
`main()` would guard one of five callers and LOOK like it guarded all five.

The resource is `os.path.relpath(os.path.dirname(d), REPO_ROOT)`, derived from
the directory actually about to be destroyed (item 74's precedent in
`run_perf.build_perf_case`), so it follows the path instead of restating it.

`REPO_ROOT` is NEW and file-derived. This module's existing `ROOT` honours
`$EQDYNAROOT`; `d` does not (it is built from `TESTSYS`). A lock rooted at
`ROOT` would be taken in one checkout while the rebuild happened in another.
Pinned by `check_the_builder_lock_ignores_a_foreign_eqdynaroot`.

`_case_lock` memo: flock is per open file description, so a second acquire in
the SAME process refuses against its own pid. Nothing stops a caller calling
the builder twice. The memo is not a fallback -- the lock is genuinely held.

## 2. `run_numa_scaling.build_case` -> lock on `testsys/perf/numa_case`

Same shape, same depth (`testsys/perf/run_numa_scaling.py:223-245`). This
module is IMPORTED by five other perf tools, so a future direct call to
`numa.build_case` must be guarded by the builder's own acquire.

## 3. `$EQDYNAROOT` mismatch in `run_jaxmpi_ab.py` -> REFUSE

`require_root_is_this_checkout()`, first statement of `main()`, before argparse
and before Gate 0 (`--notes` defaults to a path under ROOT, and the lock would
otherwise be CREATED in the foreign tree before the refusal).

REFUSE, not WARN. Evidence, checked rather than assumed -- every setter of this
variable in the repo sets it to its OWN tree:

    install-eqdyna.sh:119,152          export EQDYNAROOT=$(pwd)
    .github/workflows/test.yml (x6)    export EQDYNAROOT=$(pwd)
    testsys/e2e/run_e2e.py:174         env['EQDYNAROOT'] = REPO_ROOT  (file-derived)
    testsys/e2e/run_e2e_full.py:153    env['EQDYNAROOT'] = REPO_ROOT  (file-derived)
    4 x testsys/regression/*.py        env['EQDYNAROOT'] = ROOT       (own sandbox)

No tool in `testsys/perf/` sets it at all -- they only READ it. So no supported
workflow produces a mismatch, and refusing costs nothing. The only way to get
one is a stale export (install-eqdyna.sh exports `$(pwd)` and the export
survives a `cd` into a worktree), which is exactly the damage case.

A warning would not do: `stage()` OVERWRITES two tracked files in a tree this
session does not own, and the `finally` clause writes them again on the way
out. The operator is by construction not watching that tree, so a warning in a
4-arm measurement log is a silent failure with extra text (rule 2).

## Guard

`testsys/regression/test_perf_tool_locks.py`: 10 -> 17 checks. The three new
entry-point checks hold the real lock, invoke the real code path, and assert
non-zero exit + a refusal naming the holder + the rule-2 shape + the protected
directory UNCHANGED (a sentinel file survives). The builders are invoked as the
four other tools invoke them (`rs.build_py_case(case)` in a fresh interpreter),
NOT through their `main()`, which probes NUMA topology and cpu busy-ness first
and would make the guard pass or fail on this box's tenancy.

`check_the_mismatch_refusal_does_not_fire_when_the_roots_agree` is the negative
control: a check that refused unconditionally would make every other
`run_jaxmpi_ab` check pass for the wrong reason.

Mutation-verified, mutating ONLY `testsys/runlock.py` and reverting:
  M1  holder pid line -> 'MUTANT'          -> all 3 new checks FAIL (and 4 old)
  M2  LOCK_EX -> LOCK_SH (2nd acquire wins) -> run_scaling builder exits 0,
      rebuilds the case, sentinel gone; check FAILS with
      "exited 0 while another invocation held the lock".

## Still unguarded, NOT in this item's scope

`testsys/parity/probe_plastic_traction.py:80` rmtrees and rebuilds the fixed
in-repo path `testsys/parity/probe_case/test.drv.a6` through the same
`run_e2e.make_serial_case`. Same class, one caller, no lock. Reported, not
fixed here.

`run_scaling.ROOT` (and `run_mpi_scaling`, `run_shard_scaling`,
`run_setup_probe`, `ledger`, `probe_scatter_bandwidth`) still honour
`$EQDYNAROOT` for `PYTHON_PKG` -- i.e. WHICH tree's solver source gets timed.
That is a wrong-label hazard, not a write, so item 77's refusal was applied
only to the tool that WRITES.
