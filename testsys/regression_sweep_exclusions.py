#! /usr/bin/env python3
"""The regression script(s) that are RELEASE-PROCESS checks, not per-commit
regression guards -- excluded from `testsys/run.py`'s `run_regression()`
sweep and from `testsys/ci_shard.py`'s CI shards, while staying exactly
where they are: real `test_*.py` files under `testsys/regression/`, still
directly runnable, still the release process's own explicit check (paired
with `testsys/regression/check_pretag_ci.py`, PROJECT_RULES rule 15/25).

THE INCIDENT THIS FIXES (2026-09-30). `test_release_complete.py` evaluates
its assertions against whatever `VERSION` currently says, for a tag that
ALREADY exists -- meaningful only at the moment of auditing a just-cut
release, never at an arbitrary later commit. `v5.20.2`'s own board
Tasks-done row and its GitHub Release were both created in commits/actions
AFTER the tag (a sequencing mistake), so two of that file's checks
(`check_pathway_tasks_done_row`, and what was then `check_network_side` --
split 2026-09-30 into `check_tag_pushed_to_remote` and
`check_github_release_published`, see that file's own module docstring) can
never pass again for that exact tagged sha -- a fact about the FROZEN,
IMMUTABLE tagged tree (rule 8), not about any commit that comes after it.

`testsys/run.py`'s `run_regression()` swept every `test_*.py` under
`testsys/regression/` unconditionally (`glob.glob('test_*.py')`), so it
swept this one file into the generic tier that runs on EVERY commit and
EVERY PR's `unit-regression` CI shard, and -- via that same command --
into `.github/workflows/publish.yml`'s `build-and-push` gate, which runs
`testsys/run.py unit regression` INSIDE the built image before the image is
pushed. Every ordinary commit's CI, and every release's Docker publish,
therefore failed for a reason that has nothing to do with that commit's own
content: `v5.20.2`'s image was never published (confirmed 404 on the GHCR
manifest) because the very first gate step after the tag push permanently
fails this way.

Excluded here, explicitly, rather than silently renamed, moved, or hidden:
the file stays exactly where it is, still named `test_*.py`, still a real
regression script -- it is simply not part of the automatic sweep every
OTHER `test_*.py` file here joins by just existing on disk. Both
`testsys/run.py`'s `run_regression()` and `testsys/ci_shard.py`'s shard
partition (`SHARDS` + `verify()`) import this ONE list so the two can never
independently drift on what is excluded -- exactly the "one shared copy"
discipline `testsys/content_key.py`'s `ALLOWED_EXACT_PATHS` /
`testsys/change_class.py`'s classifier already use elsewhere in this repo.

This module has NO other imports (stdlib or otherwise) so both
`testsys/run.py` and `testsys/ci_shard.py` -- the latter deliberately
import-independent of `run.py`, `matrix.py` and `run_e2e.py` per its own
module docstring -- can import it without acquiring any new dependency.
"""

# Basenames only, matched against os.path.basename(<test_*.py path>).
EXCLUDED_FROM_SWEEP = frozenset({
    'test_release_complete.py',
})
