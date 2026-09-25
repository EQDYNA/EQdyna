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

Not yet done at this checkpoint: the actual FULL-mode run (`EQDYNA_README_
GATE=full` / `run.py readme`) against the real README.md, on a real fresh
clone. That is next -- expect it to surface whatever gap prompted this
mission (the docstring already names the one candidate found by inspection:
`pip install jax` is prose/inline in the "Python solver" section, not
inside a fenced block, so a gate that runs fenced blocks only will not
install jax before `python3 -m eqdyna . --backend jax` runs).
