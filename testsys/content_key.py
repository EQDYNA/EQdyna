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

ITEM 3 (owner, 2026-09-30, "fewer full sweeps and fewer releases") adds a
SECOND, narrower key alongside `compute` above: `compute_release_physics`,
over only `change_class.is_release_physics_path`'s closed, owner-named set
(solver source, case inputs, reference oracles, the case name list, the
pass/fail-judging modules, the case-pipeline scripts) rather than "every
tracked path outside the small evidence/ledger allow-list". A docs/ edit,
a README change, or a board update moves `compute`'s key (forcing, before
this item, a full re-sweep for content that never touched physics) but
does NOT move `compute_release_physics`'s -- see that function's own
docstring and testsys/regression/check_pretag_ci.py's release_physics_key
fallback for how the two keys are used together, `compute` first (strict,
existing behaviour unchanged) and `compute_release_physics` only as a
fallback when `compute` no longer matches.
"""
import hashlib
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import change_class  # noqa: E402

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
    """sha256 hex digest over every tracked (path, mode, object-sha) triple
    at `sha`
    NOT covered by `path_allowed` (mode included so a mode-only change,
    e.g. a lost exec bit on a script, moves the key exactly as it shows up
    in the ancestor rule's `git diff --name-only`), sorted so the result depends only on
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
        mode, _otype, blob = meta.split()
        if not path_allowed(path):
            entries.append('%s\0%s\0%s' % (path, mode, blob))
    entries.sort()
    h = hashlib.sha256()
    for entry in entries:
        h.update(entry.encode('utf-8', errors='surrogateescape'))
        h.update(b'\n')
    return h.hexdigest()


def compute_release_physics(repo_root, sha):
    """sha256 hex digest over every tracked (path, mode, blob-or-normalized-
    content) triple at `sha` that testsys.change_class.is_release_physics_path
    selects (item 3, owner 2026-09-30) -- the owner's narrower, closed
    release-physics set, NOT `compute`'s broad "everything outside the
    evidence/ledger allow-list". This is what lets check_pretag_ci.py accept
    a PRIOR release's swept evidence for a tag whose tree differs from that
    swept tree in only non-physics-bearing ways (a docs/ edit, a README
    change, a board update) -- `compute`'s stricter key would already have
    moved for any such edit and forced a needless full re-sweep.

    change_class.VERSION_BANNER_FILE is special-cased: its content is read
    with `git cat-file` and every line that is the version banner
    (change_class.is_release_physics_line) is dropped before hashing, so a
    version-only bump to that one line does not move this key either (item
    3's one explicit exclusion) -- every other line of that file, and every
    byte of every other release-physics path, still counts via its blob sha,
    exactly as `compute` does for its own broader set.

    Raises `subprocess.CalledProcessError` under the same conditions as
    `compute` above (sha not a readable object in repo_root)."""
    out = subprocess.run(
        ['git', '-C', repo_root, 'ls-tree', '-r', '--full-tree', sha],
        capture_output=True, text=True, check=True).stdout
    entries = []
    for line in out.splitlines():
        if not line.strip():
            continue
        meta, path = line.split('\t', 1)
        mode, _otype, blob = meta.split()
        if not change_class.is_release_physics_path(path):
            continue
        if path.replace(os.sep, '/') == change_class.VERSION_BANNER_FILE:
            content = subprocess.run(
                ['git', '-C', repo_root, 'cat-file', '-p', blob],
                capture_output=True, text=True, check=True).stdout
            kept = [l for l in content.splitlines()
                   if change_class.is_release_physics_line(path, l)]
            normalized_digest = hashlib.sha256(
                '\n'.join(kept).encode('utf-8', errors='surrogateescape')).hexdigest()
            entries.append('%s\0%s\0N:%s' % (path, mode, normalized_digest))
        else:
            entries.append('%s\0%s\0%s' % (path, mode, blob))
    entries.sort()
    h = hashlib.sha256()
    for entry in entries:
        h.update(entry.encode('utf-8', errors='surrogateescape'))
        h.update(b'\n')
    return h.hexdigest()
