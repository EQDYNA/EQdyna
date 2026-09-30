#! /usr/bin/env python3
"""ONE shared change classifier (owner, 2026-09-30: "fewer full sweeps and
fewer releases"; victor-reyes design review, same day, built to THAT, not
the earlier sketch this docstring's history still shows in git blame).

Every tracked path is exactly one of PHYSICS, USER_FACING, or INTERNAL,
decided in THIS PRECEDENCE ORDER: USER_FACING's short, explicit list is
checked first, then INTERNAL's short, explicit list, and everything else
-- DEFAULT-PHYSICS, deliberately, not merely as an unclassifiable fallback:
a path this module has never heard of, or a testsys/ file nobody has
carved out as internal yet, gets the full-sweep-and-PR-required treatment
rather than silently riding along (rule 2, fail toward the most expensive
answer).

TWO SEPARATE QUESTIONS, TWO SEPARATE FUNCTIONS (review item 6). Editing
testsys/matrix.py or testsys/e2e/ changes what the GATE checks -- it needs
a fresh sweep to prove the gate itself still passes (classify_path says
PHYSICS, so needs_sweep.py and pr_policy.py both require it) -- but it does
NOT by itself mean the SOLVER'S OUTPUT changed, so it must not make rule
27's release-due-AT-ONCE trigger fire the way an actual solver-code or
reference change does. `classify_path` answers "does this need a sweep /a
PR"; `is_output_change` answers rule 27's narrower "did output actually
change", reusing the ORIGINAL, unchanged check_release_due.py definition
(CODE_DIRS with its line-level comment/blank filter, OUTPUT_FILES,
REFERENCE_DIR) -- moved here so there is exactly one copy, never two.

REUSED BY (no second copy of any of this):
  - testsys/regression/check_release_due.py -- rule 27's DUE-at-once test
    is is_output_change(); rule 27 (b)'s PR-count now counts only
    classify_path(path) == USER_FACING paths per PR (review item 6: an
    INTERNAL-only PR never counts, and PHYSICS alone -- a matrix.py-only
    PR, say -- needs a sweep but does not itself count toward (b) either;
    see check_release_due.py's own docstring for the exact rule).
  - testsys/needs_sweep.py -- the advisory "does this range need
    testsys/run.py e2e" tool (classify_path == PHYSICS on any changed path).
  - testsys/pr_policy.py -- rule 25's gated-path set is WIDENED (never
    narrowed) by classify_path == PHYSICS, alongside the historical
    GATED_PREFIXES (src/, testsys/, .github/) it already enforced; see
    pr_policy.py's own module docstring for why the union, not a
    replacement.

WHY EACH NON-OBVIOUS CALL, review items 2 and 8:
  - Default PHYSICS for testsys/*.py not on the internal list: run.py
    dispatches every tier; common.py/runlock.py/profile_record.py/
    profile_schema.py are imported by testsys/e2e/run_e2e.py and
    testsys/perf/ledger.py (profile_record's own append path); ledger.py
    and run_numa_scaling.py are imported BY profile_record/run_e2e in
    turn. A silent bug in any of these can produce a wrong-but-green sweep
    without a single src/ line changing. testsys/perf/ as a WHOLE is
    likewise NOT on the internal list (only testsys/unit/, testsys/
    regression/, testsys/hooks/ are) -- ledger.py and run_numa_scaling.py
    are the two review-named reasons, and the rest of testsys/perf/ (the
    opt-in, report-only perf tools) defaults PHYSICS too rather than
    inventing a THIRD, narrower carve-out with its own failure mode.
  - scripts/ defaults PHYSICS (case.setup, lib.py, defaultParameters.py,
    machines.py, src_hash.py -- explicitly named by the review -- and
    every other script, INCLUDING ones loaded by file path rather than
    package import, e.g. case.setup's `from lib import ...` and
    testsys/run.py's `_load_by_path` for scripts/machines.py and
    scripts/src_hash.py: a path-based load is not exempt just because it
    is not a package import) EXCEPT the legacy MATLAB ground-motion
    post-processing scripts (*.m) and scripts/figures/, which are
    classified INTERNAL: dead/legacy analysis tooling, invoked by nothing
    this repo's own gates or CLIs load, and not a build/mesh/case/plot
    tool that shapes what a RUN of the solver does.
  - Dockerfile/.dockerignore/.github/workflows/publish.yml/VERSION/
    README.md/ubuntu.env.sh are USER_FACING, not PHYSICS (review item 8,
    correcting an earlier draft that had install-eqdyna.sh and
    ubuntu.env.sh classified the same way): install-eqdyna.sh COMPILES the
    solver (a botched build changes what binary every cell runs), so it
    stays PHYSICS; ubuntu.env.sh only bootstraps OS package dependencies
    (apt-get) and never touches the compiler invocation, so it is
    USER_FACING like the rest of the install/distribution surface.
  - docs/user/ is USER_FACING (the owner's own "?" -- decided yes: this is
    what a user reads to install and run a case); the REST of docs/ is
    INTERNAL (notes, evidence, session logs, the perf ledger/snapshots,
    board history).
  - Any `check_*` file (by basename, wherever it lives) is INTERNAL: it is
    a checker ABOUT the release/board process, not part of what it checks.
  - This file itself is PHYSICS (review's explicit instruction): it drives
    pr_policy's gate and needs_sweep's verdict, so a bug in IT is exactly
    as dangerous as a bug in the gate it feeds.

.github/ is classified INTERNAL, decided deliberately (the owner's own
request to say why): a workflow file changes HOW CI runs -- job order,
shard count, timeouts, which branch triggers what -- never WHAT the solver
computes. Rule 15b already draws exactly this line for CI generally ("CI
never re-verifies physics"). The one way a `.github/` edit could matter
physics-wise is if it changed which CELLS CI runs, and that membership
lives in testsys/matrix.py's CI_CELLS -- a file already classified PHYSICS
-- so a workflow-only edit cannot silently narrow physics coverage without
ALSO touching a path this classifier separately catches. This does not
weaken rule 25's protection for .github/ itself: pr_policy.py keeps its
own historical GATED_PREFIXES entry for '.github/' regardless of this
classifier's answer (see pr_policy.py's docstring) -- INTERNAL here means
"no sweep needed", not "no PR needed".
"""
import os

PHYSICS = 'PHYSICS'
USER_FACING = 'USER_FACING'
INTERNAL = 'INTERNAL'

# --- rule 27's NARROW "output actually changed" set (is_output_change),
# unchanged from the original check_release_due.py; deliberately NOT the
# same thing as "needs a sweep" (classify_path == PHYSICS) -- see the
# module docstring's "TWO SEPARATE QUESTIONS" section.
CODE_DIRS = ('src/fortran/', 'src/python/eqdyna/')
OUTPUT_FILES = ('src/fortran/library_output.f90',
                'src/python/eqdyna/library_output.py',
                'scripts/plotRuptureDynamics')
REFERENCE_DIR = 'test.reference.results/'

# --- USER_FACING: a short, explicit list, checked FIRST (docs/user/ must
# win over the broader docs/ -> INTERNAL rule below).
USER_FACING_EXACT = frozenset({
    'VERSION', 'README.md', 'Dockerfile', '.dockerignore', 'ubuntu.env.sh',
    '.github/workflows/publish.yml',
})
USER_FACING_PREFIXES = ('docs/user/',)

# --- INTERNAL: a short, explicit list, checked SECOND. Everything else is
# PHYSICS by default (see module docstring).
INTERNAL_PREFIXES = ('testsys/unit/', 'testsys/regression/', 'testsys/hooks/',
                     'docs/', '.github/')
INTERNAL_EXACT = frozenset({
    'testsys/ci_shard.py', 'testsys/ci_status.py', 'testsys/pr_policy.py',
    'testsys/check_board_separation.py',
    'testsys/check_git_config_no_credential.py',
    'pathway_forward.md', 'PROJECT_RULES.md', 'CLAUDE.md',
})


def _is_check_star(path):
    return os.path.basename(path).startswith('check_')


def _is_legacy_scripts_asset(path):
    """*.m (legacy MATLAB ground-motion post-processing) and scripts/figures/
    -- dead/reference material, loaded by nothing this repo's gates or CLIs
    import or exec."""
    if not path.startswith('scripts/'):
        return False
    return path.endswith('.m') or path == 'scripts/figures' or path.startswith('scripts/figures/')


def classify_path(path):
    """PHYSICS, USER_FACING, or INTERNAL for `path` (repo-relative, forward
    slashes, no leading './'). See the module docstring for precedence and
    reasoning; unclassifiable -> PHYSICS by construction (the final
    `return PHYSICS` below, not a special case)."""
    path = path.replace(os.sep, '/')
    if path in USER_FACING_EXACT or path.startswith(USER_FACING_PREFIXES):
        return USER_FACING
    if (path in INTERNAL_EXACT or path.startswith(INTERNAL_PREFIXES)
            or _is_check_star(path) or _is_legacy_scripts_asset(path)):
        return INTERNAL
    return PHYSICS


def is_comment_or_blank(line, path):
    """True only for a changed line we can PROVE is non-code (rule 27
    review item 7: an OpenMP sentinel, `!$` / `!$OMP`, is semantically CODE
    -- it can add or remove parallelism -- even though it starts with a
    Fortran comment character, so it must NOT be swept up as a comment)."""
    s = line.strip()
    if not s:
        return True
    if path.endswith(('.f90', '.F90', '.f', '.F')):
        if s[:2].upper() == '!$':
            return False
        return s.startswith('!')
    if path.endswith('.py'):
        return s.startswith('#')
    return False                       # unknown language: count it IN


def changed_lines(diff_text):
    """The +/- body lines of a unified diff, headers excluded."""
    out = []
    for l in diff_text.splitlines():
        if l.startswith(('+++', '---', '@@', 'diff ', 'index ', 'new file',
                         'deleted file', 'similarity', 'rename ', 'Binary')):
            if l.startswith('Binary'):
                out.append(l)          # binary change: cannot classify -> IN
            continue
        if l.startswith(('+', '-')):
            out.append(l[1:])
    return out


def is_output_change(path, diff_text=None):
    """True for rule 27's ORIGINAL, narrow "physics or output change"
    criteria 1-3 ONLY -- never widened by classify_path's broader PHYSICS
    bucket (review item 6): editing testsys/matrix.py or scripts/lib.py
    needs a sweep (classify_path says PHYSICS) but is not, by itself, an
    output change.

    `diff_text` is the unified diff (`-U0`) of `path` in one commit; if
    omitted, a CODE_DIRS path counts IN unconditionally (conservative:
    no diff to prove otherwise, rule 2)."""
    if path.startswith(REFERENCE_DIR):
        return True
    if path in OUTPUT_FILES:
        return True
    if path.startswith(CODE_DIRS):
        if diff_text is None:
            return True
        lines = changed_lines(diff_text)
        if not lines:
            return True
        return not all(is_comment_or_blank(l, path) for l in lines)
    return False
