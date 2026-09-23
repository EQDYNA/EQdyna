# item 70 -- a lock on the rotated run tree (working notes, NOT committed)

Branch `fix/item70-e2e-test-lock`, worktree `/home/utig5/dliu/EQdyna-wt-item70`.

## H1 -- flock beats an O_EXCL lockfile here. HELD.
Both refuse a second holder; they differ on holder death. flock is released by
the kernel, so a SIGKILLed sweep leaves nothing to clean up. An O_EXCL file
outlives its holder and forces every acquisition to guess liveness from a pid,
which is wrong under pid reuse in both directions -- and a lock that needs a
human to unstick it after a crash gets worked around. `testsys/perf/ledger.py`
already uses flock, so this is the existing convention, not a new one.
Consequence accepted and documented: flock cannot name the holder, so the
holder writes a pid/started/cwd/argv record INTO the locked file and a refuser
reads it. Record = provenance; flock = the lock. Missing record -> refuse
anyway, saying the record is missing.

## H2 -- take the lock before Gate 1, not just before the rotation. HELD.
Brief said "before the rotation". Taking it at Gate 0 (ahead of the
create.newcase guard and the Fortran build) is strictly stronger, costs
nothing, and means a doomed invocation is refused before it spends a build.
It also incidentally serialises `bin/`.

## H3 -- the guard must assert the tree was NOT renamed, not just the exit code.
HELD, and mutation B below is the proof: a tool that rotates and THEN refuses
exits non-zero and prints the refusal, so an exit-code-only guard passes on the
exact defect. The sentinel-file assertion is the one with teeth.

## Mutation verification (contained: no mutant reached a build or a solver)
- MUTATION A, `fcntl.flock(...)` deleted from `acquire()`:
  KILLED `check_a_second_acquire_refuses`, KILLED `check_the_refusal_does_not_block`.
  `check_a_dead_holder_is_taken_over...` SURVIVED -- correctly: it tests the
  stale-record path, which does not depend on flock's exclusion. Noted, not
  papered over.
- MUTATION B, the pre-item-70 unguarded rotation reinstated ahead of the lock:
  guard exit 1, KILLED `check_run_e2e_refuses_and_does_not_rotate` (sentinel
  gone from the live tree) and KILLED
  `check_the_lock_is_taken_before_the_rotation_in_both_tools` (acquire at line
  579, rotate at line 567).
- Clean tree: 9/9 PASS, exit 0.

## Incident during this session, recorded because it cost real time
My FIRST mutation script ran the whole guard under MUTATION A. With flock
neutralised, `run_e2e.py` did not refuse -- it proceeded through Gate 1 into
`install-eqdyna.sh` (a real Fortran build) and `run_e2e_full.py` then found the
freshly installed `bin/eqdyna` and launched `mpirun -np 16 eqdyna` on a box
reserved for a 32-rank perf measurement. Killed within ~2 min; `bin/`,
`test*/`, the perf-ledger append and its snapshot file were all reverted.
LESSON, general: a mutant that disables a REFUSAL turns a cheap guard into the
expensive thing the guard exists to avoid running. Mutate the refusal only
against in-process checks; mutate the ORDERING when you need the entry point.

## Two-invocation demonstration (both invocations real)
#1 `run_e2e.py --cases test.tpv8 --backends python-numpy` (no fortran cell ->
no build, one process). #2 the same command while #1 was live: exit 1, refusal
naming pid/started/cwd/cmd, no rotation. #1 then SIGKILLed; the next acquire
printed the stale-takeover NOTICE naming the dead pid.

## Out of scope, confirmed untouched
`testsys/perf/run_perf.py:172`, `testsys/perf/run_jaxmpi_ab.py:88`. Same defect
class; `runlock.py` is generic over the resource name and fits both unchanged.
