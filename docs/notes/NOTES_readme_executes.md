# Session notes: README-executes gate (mission: stranger-clone, mechanical)

Worktree `agent-a54ea4e3f9ec4f45d`, branch `wei/readme-accurate`, base
`366a88d` (conductor's WIP README rewrite).

## Checkpoint 1

Landed:

* `testsys/regression/test_readme_executes.py` -- new. One parser + one
  execution engine, two depths:
  - FAST (default; what `run.py regression`/ci_shard run): parses
    README.md's fenced blocks under Requirements/Install/Quick
    start/Python-solver, asserts all four exist, asserts `e2e-full` never
    appears inside a fenced block, asserts the ONLY reason any line is
    skipped is a literal `# needs root` suffix, and runs a built-in
    mutation self-test against a tiny SYNTHETIC README (not the real one)
    -- intact copy passes, two independently-mutated copies (an export
    dropped mid-session; a command typo'd) each fail AND are attributed to
    the exact section/line. Sub-second, no clone, no network.
  - FULL (`EQDYNA_README_GATE=full`, set by the new `run.py readme`
    OPTIONAL_TIERS entry): clones THIS commit into a temp dir and executes
    the real README.md's selected fenced blocks verbatim, in ONE bash
    session, in an `env -i`-equivalent environment (fresh HOME,
    `PATH=/usr/bin:/bin`, nothing else -- no PYTHONPATH/EQDYNAROOT/venv
    inherited). Then checks the Quick start bullet list's own documented
    outputs (parsed from the README text itself) exist and are non-empty,
    and that the Python-solver step wrote a non-empty `frt.txt0`.
* Design decisions and why, all in the file's own docstring (not repeated
  here): pip lines run for real (idempotent, network reachable, and
  marking them would hide the exact bug class this gate exists to catch);
  `set -u` is load-bearing (an unbound-variable abort from a removed
  export is a REAL failure mode, found by the mutation self-test itself --
  it kills bash BEFORE the `__rc`-check wrapper runs, so `attribute_failure`
  falls back to the last `RUN README.md:N` tag echoed, not just the
  explicit `GATE-FAIL` tag).
* `testsys/run.py` -- added `run_readme()` (sets `EQDYNA_README_GATE=full`,
  calls the same file above), registered in `RUNNERS`/`OPTIONAL_TIERS`.
  `release` now expands to `('release', 'readme')` via a new
  `RELEASE_EXPANDS` tuple (mission: "include it in run.py release") -- the
  release gate now also re-proves the README a new user would actually
  follow, at release time (not on every `unit regression` invocation).
* `testsys/ci_shard.py` -- registered `test_readme_executes.py` in shard 2
  (FAST mode only; sub-second; the FULL mode never fires here since no
  shard sets `EQDYNA_README_GATE`).
* `testsys/regression/test_perf_item91_guards.py` --
  `check_91a_backfill_detects_shard_schema` fixed: it wrote its temp
  snapshot under this TEST FILE's own repo root while
  `ledger._rows_from_snapshot_file` resolves against `ledger.ROOT`
  (follows `$EQDYNAROOT`); any process with EQDYNAROOT pointed elsewhere
  (which is exactly what building the README gate does, on purpose) made
  this test die via a bare `SystemExit` instead of reporting a clean FAIL.
  Fixed to write under `ledger.ROOT`'s own `docs/perf_snapshots`, and to
  catch `SystemExit` from the function under test as a FAIL, not a crash.
  Mutation-verified: reproduced the old crash directly against `ledger`
  with a foreign `EQDYNAROOT`, confirmed the fixed test passes clean under
  the same foreign `EQDYNAROOT`.

Verified so far: `python3 testsys/run.py unit regression` -- SUCCESS
(both tiers), including the new/edited files, against a real
`./install-eqdyna.sh -m ubuntu` build on this branch.

## Checkpoint 2 -- ran the FULL gate for real against README.md

`python3 testsys/run.py readme` (== `EQDYNA_README_GATE=full python3
testsys/regression/test_readme_executes.py`), real clone of HEAD (`28fbff2`
at the time), real `env -i`-equivalent environment: **FAILS**, correctly.
Two independent README defects found, in order:

1. **Blocking, first hit.** README.md:55 `create.newcase ~/runs/tpv8
   test.tpv8` -- `scripts/create.newcase:49` calls `os.mkdir(case_path)`,
   which requires the PARENT directory (`~/runs`) to already exist. On a
   genuinely fresh clone/`$HOME`, `~/runs` does not exist, so this is
   `FileNotFoundError: [Errno 2] No such file or directory:
   '<HOME>/runs/tpv8'`. This is either a missing `mkdir -p ~/runs` step in
   the README or a `create.newcase` fix (`os.mkdir` -> `os.makedirs`) --
   both are outside this session's remit (README prose / production code).
   Gate output (`GATE-FAIL section=Quick_start line=55 rc=1`) names it
   exactly; verified with `run.py readme` invoked for real a second time
   (a fresh clone reproduces it identically at the same line).

2. **Found by manual continuation past #1** (same kept work dir,
   `create.newcase` run by hand with `~/runs` pre-created, purely to check
   whether a SECOND defect exists downstream -- not part of the committed
   gate's own recorded run, and no README/production file was touched to
   do this). Once past #1, the Fortran Quick start (README.md:55-58) runs
   clean end to end: `case.setup`, `bash run.sh` (4 MPI ranks) all exit 0,
   and every one of the 5 documented outputs
   (`cRuptureDynamics.png`, `frt.txt0`/`frt.txt2`, 8 `faultst*.txt`,
   11 `body*.txt`, `fault.dyna.r.nc`) exists and is non-empty. The Python
   solver section (README.md:80-84) then fails:
   `ModuleNotFoundError: No module named 'jax'` at
   `src/python/eqdyna/eqdyna3d.py:895` (`_select_device` -> `import jax`).
   Cause: `pip install jax` (README, "Python solver" paragraph, prose
   sentence, single-backtick inline code) is NOT inside a fenced code
   block -- a user (or gate) that runs only the fenced blocks under
   "Python solver" never installs jax. This is the SAME bug class this
   whole gate was built to catch (see this file's own module docstring's
   motivating incident: `python -m eqdyna.standalone`/ModuleNotFoundError).

Per mission instruction ("If the README fails your gate, REPORT the
failing step and the exact error; do not rewrite the README yourself"):
neither defect was fixed here. Kept evidence: `/tmp/eqdyna_readme_gate_
vp6p06ht/work` (the `run.py readme` invocation whose output is quoted
above); earlier kept dirs from the same investigation:
`/tmp/eqdyna_readme_gate_okbfceyj/work` (manually continued past #1 to
find #2, in the SAME dir -- so its `home/runs/` contains my manual
continuation, not only the gate's own output).

The FAST tier (parser + mutation self-test, what `run.py regression` and
CI actually run today) is unaffected by either defect and stays green --
by design, since FULL mode is opt-in and never runs at seconds-scale.

## Gate status at this checkpoint

* `python3 testsys/run.py unit regression` -- SUCCESS (real build via
  `./install-eqdyna.sh -m ubuntu`, this branch, this worktree).
* `python3 testsys/run.py readme` -- **FAILS**, correctly, on a real
  README defect (#1 above) that predates this session (conductor's WIP
  commit `366a88d`). This is the gate doing its job, not a gate bug --
  see the mutation-test evidence above and in the file's own checks for
  why a false failure here is implausible (the same harness, on a
  synthetic README, distinguishes an intact copy from two independently
  broken ones cleanly). Reported to the conductor/owner rather than
  patched.
