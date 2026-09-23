# Item 74 -- lock the two testsys/perf/ tools that rebuild fixed in-repo paths

Worktree `.claude/worktrees/iris-item74`, branch `iris/item74-perf-locks`,
branched from `51f20e0`. Not merged, not pushed. 2026-09-23.

## API shape chosen

`runlock.acquire(repo_root, resource)` keeps ONE `repo_root`; `resource`
becomes a relative PATH instead of a single directory name. The lockfile stays
BESIDE the guarded directory -- in that directory's parent, `.<basename>.lock`
-- which for a one-component resource is byte-for-byte the old rule
(`test` -> `REPO_ROOT/.test.lock`), and for `src/python/eqdyna` is
`src/python/.eqdyna.lock`.

Rejected alternative: a per-call deeper `repo_root`
(`acquire(ROOT + "/src/python", "eqdyna")`). It needs no change to runlock,
but it moves path arithmetic into every caller under an argument named
`repo_root` that is then not the repo root; a caller that joins wrong puts a
lockfile somewhere plausible and guards nothing, silently. One validated place
beats N unvalidated ones. `lock_path` now raises `ValueError` on absolute,
empty, `.`, `..` or trailing-separator forms -- guarded by the new test.

Second API change: `acquire(consequence=[...])`. What a collision COSTS
differs by tool, and a refusal that describes the wrong incident teaches the
reader to skip it. The tail of the refusal is now the caller's; the default is
the e2e text, unchanged, so item 70's guard passes untouched.

## Where the lock is taken

- `run_perf.build_perf_case()`, first statement, resource derived from
  `os.path.dirname(PERF_CASE)`. In the FUNCTION that does the rmtree, not in
  `main()`: `run_tpv29_pinned_compare.py` reassigns `run_perf.PERF_CASE` to
  `perf_case_tpv29/` and calls this function directly, never reaching
  `main()`. A lock in `main()` would guard one of two callers and look like it
  guarded both. Deriving the resource also makes the lock follow whichever
  directory is actually about to be destroyed.
- `run_jaxmpi_ab.main()`, immediately after `parse_args`, resource
  `src/python/eqdyna` -- the PACKAGE `stage()` copies arm files into, not the
  `__pycache__` beneath it that the tool happens to delete. Before the vault,
  before any staging, before the notes file is opened.

Board row 74 says the resource is `testsys/perf_case`; the real path is
`testsys/perf/perf_case` (`run_perf.TESTSYS` is `testsys/perf`, a misnomer).

## Verification (fresh, this worktree)

- `testsys/regression/test_perf_tool_locks.py` -- 10 checks, exit 0, 0.17 s.
- `testsys/regression/test_e2e_run_tree_lock.py` -- UNCHANGED
  (md5 572a735d7c25d7b59f23ee361b77a2a5, empty `git diff 51f20e0`), 9 checks,
  exit 0.
- `python3 testsys/run.py unit regression` -- exit 0, 73 s.
- Mutations, on the LOCK HELPER only (never an entry point's refusal path):
  (1) lockfile moved INSIDE the guarded directory -> 2 checks fail, including
  the one pinning item 70's unchanged paths; (2) run_perf acquires AFTER the
  rmtree -> "THE DEFECT ITSELF: rmtree-d ... despite refusing" plus the
  source-order check. Both reverted; no measurement and no mpirun was ever run.

## Still unguarded after this change (for the board)

- `run_scaling.build_py_case` -> `testsys/perf/scaling_case/<case>` (rmtree +
  rebuild). Called by `run_scaling.py`, `run_mpi_scaling.py` and
  `run_jaxmpi_ab.py`, so a jaxmpi A/B and an MPI scaling run still collide.
- `run_numa_scaling.build_case` -> `testsys/perf/numa_case/<case>`, same shape.
- `run_jaxmpi_ab` resolves ROOT from `EQDYNAROOT` when set: run it from a
  worktree with EQDYNAROOT pointing at the main checkout and it stages arm
  files into the MAIN checkout's package (rule 21b). The lock follows ROOT, so
  the two are consistent -- but it then locks and writes the wrong tree.
