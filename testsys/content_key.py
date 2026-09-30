#! /usr/bin/env python3
"""Content-key sweep evidence (item 2b, 2026-09-30 classifier work): ONE
copy of "what counts as the tree a release sweep swept", shared by the
writer (`testsys/e2e/run_e2e.py`'s `write_release_evidence`) and the two
readers (`testsys/regression/check_pretag_ci.py`'s pre-tag guard and, via
its own re-use of that module, `test_release_complete.py`'s post-hoc check).

WHY. The pre-existing acceptance test (`_sweep_is_ancestor` +
`_sweep_disallowed_paths` in check_pretag_ci.py, unchanged by this module)
required the swept sha to be a real git ancestor of the tag commit, walked
by `git diff --name-only swept_sha tag_sha`. Two real situations break that:

  1. A squash-merge (rule 25's ONLY merge strategy) mints a brand-new commit
     for the PR; the swept sha (a commit on the PR branch, pre-squash) is
     often NOT an ancestor of that new commit at all, even though the
     merged commit's TREE is byte-identical to the branch tip's tree.
  2. A history rewrite (2026-09-30, scec_archive/ dropped, see pathway item
     36) drops old commit OBJECTS outright. A `summary.json` committed
     before the rewrite names a sha `git cat-file -e` can no longer prove
     exists at all -- `_sweep_is_ancestor` has no object to compare.

Content-key evidence answers a narrower, more robust question: not "is the
swept sha reachable from the tag", but "does every tracked path OUTSIDE the
narrow evidence/ledger/board allow-list have the SAME BLOB at the swept
commit as at the tag commit". That question needs no ancestry at all, and
needs the swept commit's git OBJECT to exist only if you insist on
recomputing its key from git -- which nothing here does: the swept side's
key is read straight out of the JSON string `summary.json` already
recorded, at write time, when the object certainly still existed. Two trees
with the same PHYSICS-and-code content (identical outside the allow-list)
are the same evidence, regardless of how they are related in the commit
graph or whether the older one survived a rewrite.

The allow-list is intentionally the SAME set check_pretag_ci.py's older
ancestor-diff check already used (docs/evidence/, docs/perf_snapshots/,
docs/perf_ledger.jsonl, docs/run_profiles.jsonl, pathway_forward.md) -- not
widened, not narrowed: still "evidence, ledgers, and the board", the paths
rule 15d's own two-commit release order already expects to change after a
sweep with no re-sweep required. Anything else -- including `.github/`
(INTERNAL to `testsys/change_class.py`, but that classifier answers a
different question: whether a PR is needed, not whether a sweep is still
valid) and `scripts/machines.py` -- changing so much as one byte moves the
key and invalidates the evidence.
"""
import hashlib
import subprocess

ALLOWED_EXACT_PATHS = ('docs/perf_ledger.jsonl', 'docs/run_profiles.jsonl',
                        'pathway_forward.md')
ALLOWED_PATH_PREFIXES = ('docs/evidence/', 'docs/perf_snapshots/')


def path_allowed(path):
    """True for a path the content key deliberately ignores (evidence,
    ledgers, the board) -- see the module docstring for why this list is
    narrow and unchanged from check_pretag_ci.py's pre-existing one."""
    if path in ALLOWED_EXACT_PATHS:
        return True
    return any(path.startswith(prefix) for prefix in ALLOWED_PATH_PREFIXES)


def compute(repo_root, sha):
    """sha256 hex digest over every tracked (path, blob-sha) pair at `sha`
    NOT covered by `path_allowed`, sorted so the result depends only on
    CONTENT, never on `git ls-tree`'s own (already-sorted, but not
    guaranteed stable across git versions for every locale) output order.

    Raises `subprocess.CalledProcessError` if `sha` is not a readable git
    object in `repo_root` -- this function is only ever called with a sha
    this checkout is known to have (the tag commit itself); a caller that
    wants to compare against a swept sha that might be GONE (a rewritten
    history) reads that side's key from the committed JSON string instead
    of calling this a second time (see check_pretag_ci.py)."""
    out = subprocess.run(
        ['git', '-C', repo_root, 'ls-tree', '-r', '--full-tree', sha],
        capture_output=True, text=True, check=True).stdout
    entries = []
    for line in out.splitlines():
        if not line.strip():
            continue
        meta, path = line.split('\t', 1)
        blob = meta.split()[2]
        if not path_allowed(path):
            entries.append('%s\0%s' % (path, blob))
    entries.sort()
    h = hashlib.sha256()
    for entry in entries:
        h.update(entry.encode('utf-8', errors='surrogateescape'))
        h.update(b'\n')
    return h.hexdigest()
