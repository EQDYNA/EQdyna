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
PHYSICS, so pr_policy.py requires a PR for it) -- but it does
NOT by itself mean the SOLVER'S OUTPUT changed, so it must not make rule
27's release-due-AT-ONCE trigger fire the way an actual solver-code or
reference change does. `classify_path` answers "does this need a sweep /a
PR"; `is_output_change` answers rule 27's narrower "did output actually
change", reusing the ORIGINAL, unchanged check_release_due.py definition
(CODE_DIRS with its line-level comment/blank filter, OUTPUT_FILES,
REFERENCE_DIR) -- moved here so there is exactly one copy, never two.

REUSED BY (no second copy of any of this):
  - testsys/regression/check_release_due.py -- rule 27's DUE-at-once test
    is is_output_change(); rule 27 (b)'s PR-count excludes only a PR whose
    EVERY changed path classifies INTERNAL (review item 6, corrected from
    an earlier draft that narrowed the count to USER_FACING-only): a
    PHYSICS-only PR -- a matrix.py-only PR, say -- needs a sweep AND still
    counts toward (b)'s threshold; classify_path == PHYSICS on its own just
    does not trip the AT-ONCE trigger (is_output_change is narrower than
    that); see check_release_due.py's own docstring for the exact rule.
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
  - Any `check_*` file, by basename, under `testsys/` (root or any
    subdirectory) is INTERNAL: it is a checker ABOUT the release/board
    process, not part of what it checks. Scoped to `testsys/` (victor-reyes
    audit of PR #62, finding 4, corrected from "wherever it lives"): a
    `check_*` name outside `testsys/` is not this repo's release-process
    tooling and must not be exempted from review just by that name.
  - This file itself is PHYSICS (review's explicit instruction): it drives
    pr_policy's gate and the release sweep carry-forward, so a bug in IT is exactly
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
import re

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
# win over the broader docs/ -> INTERNAL rule below). `Docker.guide.md`
# added 2026-09-30 (victor-reyes audit of PR #62, finding 2): a root-level
# user install guide, sibling to Dockerfile/README.md/VERSION, not physics.
USER_FACING_EXACT = frozenset({
    'VERSION', 'README.md', 'Dockerfile', '.dockerignore', 'ubuntu.env.sh',
    '.github/workflows/publish.yml', 'Docker.guide.md',
})
USER_FACING_PREFIXES = ('docs/user/',)

# --- INTERNAL: a short, explicit list, checked SECOND. Everything else is
# PHYSICS by default (see module docstring). `LICENSE`, `pastReleaseNotes.md`,
# `.gitignore` added 2026-09-30 (victor-reyes audit of PR #62, finding 2):
# plain doc/meta files outside src/testsys/.github that the PHYSICS-default
# was newly gating for a PR with no physics content to review at all.
INTERNAL_PREFIXES = ('testsys/unit/', 'testsys/regression/', 'testsys/hooks/',
                     'docs/', '.github/')
INTERNAL_EXACT = frozenset({
    'testsys/ci_shard.py', 'testsys/ci_status.py', 'testsys/pr_policy.py',
    'testsys/check_board_separation.py',
    'testsys/check_git_config_no_credential.py',
    'pathway_forward.md', 'PROJECT_RULES.md', 'CLAUDE.md',
    'LICENSE', 'pastReleaseNotes.md', '.gitignore',
})

# --- default-PHYSICS, but not UNREVIEWED (victor-reyes audit of PR #62,
# finding on the totality test): every path that falls all the way through
# to the bare `return PHYSICS` at the end of classify_path is expected to
# live under one of these prefixes or be one of these exact scripts --
# reviewed, deliberately, same day as the classifier itself. A tracked path
# that classifies PHYSICS WITHOUT matching one of these (or the USER_FACING/
# INTERNAL lists above) is a review gap: some brand-new top-level file or
# directory nobody has looked at yet, not something this module already
# reasoned about. `explicit_bucket` below is what a totality test asserts
# against -- classify_path's own return value never depends on it, so this
# list growing narrower cannot silently change any gate's answer, only the
# totality test's opinion of whether the answer was reviewed.
DEFAULT_PHYSICS_REVIEWED_PREFIXES = ('src/', 'testsys/', 'case_input/',
                                     'test.reference.results/', 'scripts/')
DEFAULT_PHYSICS_REVIEWED_EXACT = frozenset({'testNameList.py', 'install-eqdyna.sh'})


# --- ITEM 3 (owner, 2026-09-30, "fewer full sweeps and fewer releases"): a
# THIRD question, again answered by its own function (see "TWO SEPARATE
# QUESTIONS" above -- now three): not "does a PR need review" (classify_path
# / PHYSICS, which fails toward the expensive answer for anything testsys/
# or scripts/ has not been explicitly carved out of) and not "did rule 27's
# release-due-at-once trigger fire" (is_output_change), but "does the LOCAL
# RELEASE SWEEP need to run again before this exact tree can be tagged, or
# can the last swept release's evidence carry forward" -- the ONE answer to
# "does this change need a full sweep" (owner 2026-09-30: "no sweep when no
# physics change"; testsys/needs_sweep.py, a second, broader answer, deleted) -- consumed by
# testsys/content_key.py's compute_release_physics,
# testsys/regression/check_pretag_ci.py's release_physics_key fallback, and
# testsys/run.py's run_release(). This is the owner's OWN closed, exact
# list, deliberately NARROWER than classify_path's PHYSICS bucket: it is not
# "everything under testsys/ or scripts/ classify_path would call PHYSICS",
# it is exactly the set the owner named -- solver source, case inputs,
# reference oracles, the case name list, the three pass/fail-judging
# modules, and the case-pipeline scripts that shape what a run computes.
# Do not add to or remove from this list without the owner's sign-off: too
# WIDE defeats item 3's whole point (fewer full sweeps); too NARROW ships a
# stale sweep for a tree that actually changed physics.
RELEASE_PHYSICS_PREFIXES = ('src/', 'case_input/', 'test.reference.results/')
RELEASE_PHYSICS_EXACT = frozenset({
    'testNameList.py',
    'testsys/matrix.py', 'testsys/compare.py', 'testsys/frt_canonical.py',
    # The case-pipeline inputs (owner's own phrase): case.setup and
    # create.newcase THEMSELVES, plus every scripts/ file either one
    # actually `import`s or reads FOR ITS OWN LOGIC -- not the *.m /
    # plotting/analysis utilities create.newcase merely byte-copies
    # alongside them into every new case directory (those are already
    # INTERNAL above via _is_legacy_scripts_asset, same "shapes nothing a
    # run computes" reasoning). scripts/lib.py and scripts/machines.py are
    # case.setup's own direct imports (case.setup:23-24, `from lib import
    # ensureFaultRoughGeometryForCase, resolveNormalStressSign,
    # resolveViscoplasticParams` / `from machines import
    # write_slurm_header`); scripts/defaultParameters.py is not imported by
    # case.setup directly but by every case's own user_defined_params.py,
    # and is named explicitly by the owner.
    'scripts/case.setup', 'scripts/create.newcase',
    'scripts/lib.py', 'scripts/machines.py', 'scripts/defaultParameters.py',
})

# The one line item 3 carves back OUT of `src/`: a version-only bump to
# this exact runtime banner must not, by itself, demand a fresh release
# sweep -- the banner never changes what the solver computes (rule 11 keeps
# it in sync with VERSION precisely because it is otherwise inert). An
# exact-shape regex, not a loose one: a banner line that changes SHAPE (not
# just the version token) is exactly the case that must NOT be silently
# swallowed by too loose a pattern.
VERSION_BANNER_FILE = 'src/fortran/eqdyna3d.f90'
VERSION_BANNER_RE = re.compile(
    r"^\s*write\(\*,\*\)\s*'=+\s*Welcome to EQdyna\s+\S+\s+=+'\s*$")


def is_release_physics_path(path):
    """True for a path in item 3's closed, owner-named release-physics set
    (RELEASE_PHYSICS_PREFIXES/EXACT above) -- narrower than classify_path's
    PHYSICS bucket by design; see the block comment above for why the two
    must not be conflated or one derived from the other."""
    path = path.replace(os.sep, '/')
    return path in RELEASE_PHYSICS_EXACT or path.startswith(RELEASE_PHYSICS_PREFIXES)


def is_release_physics_line(path, line):
    """False only for a line of VERSION_BANNER_FILE that IS the version
    banner (item 3's one named exclusion); True for every other line of
    every other release-physics path, so a caller filtering line-by-line
    counts everything else IN (rule 2)."""
    if path.replace(os.sep, '/') != VERSION_BANNER_FILE:
        return True
    return not VERSION_BANNER_RE.match(line)


def _is_check_star(path):
    """A `check_*` basename is INTERNAL only under `testsys/` (root or any
    subdirectory, e.g. `testsys/regression/`) -- victor-reyes audit of PR
    #62, finding 4: scoped down from "wherever it lives", which would have
    silently swallowed a hypothetical `scripts/check_something.py` or
    `case_input/.../check_whatever.py` that is not a checker ABOUT this
    repo's release/board process at all."""
    return path.startswith('testsys/') and os.path.basename(path).startswith('check_')


def _is_legacy_scripts_asset(path):
    """*.m (legacy MATLAB ground-motion post-processing) and scripts/figures/
    -- dead/reference material, loaded by nothing this repo's gates or CLIs
    import or exec."""
    if not path.startswith('scripts/'):
        return False
    return path.endswith('.m') or path == 'scripts/figures' or path.startswith('scripts/figures/')


def explicit_bucket(path):
    """Which bucket `path` matches EXPLICITLY -- USER_FACING, INTERNAL, or a
    REVIEWED default-PHYSICS prefix/exact name -- or None if it falls
    through to the bare, UNREVIEWED PHYSICS default. Exists only so a
    totality test can assert every tracked path was actually thought about
    once; classify_path below never consults this and its answer for a
    given path never changes because of it (rule 2: fail toward the
    expensive gating answer regardless of review status)."""
    path = path.replace(os.sep, '/')
    if path in USER_FACING_EXACT or path.startswith(USER_FACING_PREFIXES):
        return USER_FACING
    if (path in INTERNAL_EXACT or path.startswith(INTERNAL_PREFIXES)
            or _is_check_star(path) or _is_legacy_scripts_asset(path)):
        return INTERNAL
    if (path in DEFAULT_PHYSICS_REVIEWED_EXACT
            or path.startswith(DEFAULT_PHYSICS_REVIEWED_PREFIXES)):
        return PHYSICS
    return None


def classify_path(path):
    """PHYSICS, USER_FACING, or INTERNAL for `path` (repo-relative, forward
    slashes, no leading './'). See the module docstring for precedence and
    reasoning; unclassifiable -> PHYSICS by construction (the final
    `return PHYSICS` below, not a special case -- see explicit_bucket above
    for the separate question of whether that default was ever reviewed)."""
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
