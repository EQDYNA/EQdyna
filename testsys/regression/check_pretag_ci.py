#! /usr/bin/env python3
"""
Pre-`git tag` gate (rule 15 step 6/7): "wait for CI to go green on the
pushed commit, then tag" is a sentence, not a mechanism. Nothing enforced it
before. Run this BEFORE `git tag`; it exits non-zero unless it is safe to.

THE INCIDENT (2026-09-21): `v5.13.1` was tagged at `dfee14d` while the only
CI run for that exact SHA was the run the tag push itself triggered -- at
`git tag` time there was no completed run for `dfee14d` at all. It happened
to conclude green (run 35667535066, 7/7); nothing would have noticed had it
gone red, and a pushed tag cannot be re-pointed (rule 8).

THE STRUCTURAL TRAP: `dfee14d` touches ONLY `pathway_forward.md`, one of
`.github/workflows/test.yml`'s own `paths-ignore` entries, so it could never
have had a pre-tag run of its OWN -- waiting for one would wait forever. The
workflow's comment claims "a release commit always also touches non-md files
... so it still triggers CI regardless"; true of the bundled release commit,
false of the follow-up Tasks-done-row commit rule 15 step 4 itself asks for
when it lands as its own commit rather than folded into step 5's "commit
everything above together". This is why `--ack-paths-ignored-parent`
exists below, rather than the check silently falling back to the parent on
its own: a releaser must choose to accept parent evidence, not have it
handed to them unasked.

Deliberately not auto-run by testsys/run.py's regression tier (its filename
does not match `test_*.py`): its exit codes are not "0 good, nonzero bad" --
PENDING is not a bug in the tree, it means "CI has not finished, ask again
in a minute", and gating every ordinary commit's fast tier on that would
fail CI runs for reasons that have nothing to do with the code under test.
This is a manual, human-invoked gate at exactly one moment: right before
`git tag`.

Exit codes (rule 2: situations that mean different things get different
codes; a script that prints a paragraph and exits 1 either way just means
the reader has to open the log to find out which):

  0  PASS            a completed, successful CI run exists for the SHA
                      being asked about (or, with --ack-paths-ignored-parent,
                      for the nearest ancestor that could have one).
  1  FAIL             a completed CI run for the relevant SHA finished
                      WITHOUT success -- CI ran and failed. Do not tag.
  2  PENDING          no completed run yet -- either none has started, or
                      one is in progress. Wait; do not tag.
  3  PATHS_IGNORED    this SHA touches only paths-ignore'd files and so can
                      never trigger its own CI run. Re-run with
                      --ack-paths-ignored-parent to accept the nearest
                      ancestor's evidence instead, or tag the ancestor.
  4  UNVERIFIED       gh is missing/unauthenticated, or a network call
                      failed. Not evidence of failure (rule 2) -- re-run
                      where gh can actually reach the API.
  5  SWEEP_INSUFFICIENT   CI is green for the SHA (the check above passed),
                      but no committed local release sweep justifies
                      tagging it. Added 2026-09-23 (release-path guard,
                      owner-approved design: CI stops running the e2e sweep,
                      so the LOCAL release sweep -- every supported cell,
                      everyday + RELEASE_ONLY, all at the ONE GATE_TERM_S --
                      becomes the science gate and a tag must be refused
                      without one). See `evaluate_sweep_evidence` below.

Only checked when the CI check above is PASS -- CI-not-green and
sweep-not-found are different situations (rule 2), and the CI check is more
fundamental, so it is reported alone first.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(TESTSYS)
sys.path.insert(0, ROOT)
from testsys import ci_status, matrix  # noqa: E402

EXIT_CODES = {
    ci_status.PASS: 0,
    ci_status.FAIL: 1,
    ci_status.PENDING: 2,
    ci_status.PATHS_IGNORED: 3,
    ci_status.UNVERIFIED: 4,
}
SWEEP_INSUFFICIENT_EXIT = 5

# --- release-path guard: a committed local RELEASE sweep -- every supported
# cell (everyday cells + matrix.RELEASE_ONLY), all at the ONE GATE_TERM_S --
# not just green CI (2026-09-23) ---------------------------------------------
EVIDENCE_GLOB = os.path.join('docs', 'evidence', 'sweep-*', 'summary.json')
# Exactly the allowed set for rule 15d's two-commit release order: the
# board-only Tasks-done-row commit touches only pathway_forward.md, and the
# evidence itself lands as its own commit touching only these evidence/ledger
# paths. VERSION/README/pastReleaseNotes.md are deliberately NOT here -- if
# those changed after the sweep, the sweep no longer describes the tagged
# tree's physics-relevant content and must be re-run.
# docs/run_profiles.jsonl added 2026-09-24: the release sweep appends its own
# per-rank profile rows there (profile_record.py, PR #6), exactly as it appends
# docs/perf_ledger.jsonl; without it the v5.17.0 sweep's own rows could not
# land before the tag. It is a run record, not physics input.
SWEEP_ALLOWED_EXACT_PATHS = ('docs/perf_ledger.jsonl', 'docs/run_profiles.jsonl',
                             'pathway_forward.md')
SWEEP_ALLOWED_PATH_PREFIXES = ('docs/evidence/', 'docs/perf_snapshots/')
SWEEP_REQUIRED_FIELDS = ('sha', 'tree_clean', 'term', 'n_runnable', 'n_success',
                         'cells', 'started_utc', 'finished_utc')


def _sweep_path_allowed(path):
    if path in SWEEP_ALLOWED_EXACT_PATHS:
        return True
    return any(path.startswith(prefix) for prefix in SWEEP_ALLOWED_PATH_PREFIXES)


def _sweep_git(repo_root, *args):
    r = subprocess.run(('git', '-C', repo_root) + args, capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def _sweep_is_ancestor(repo_root, ancestor, sha):
    if ancestor == sha:
        return True
    rc, _out, _err = _sweep_git(repo_root, 'merge-base', '--is-ancestor', ancestor, sha)
    return rc == 0


def _sweep_disallowed_paths(repo_root, sha_from, sha_to):
    """Paths changed between the two commits that are outside the
    evidence/ledger allow-list. Empty result required for the swept SHA to
    justify tagging a later commit."""
    if sha_from == sha_to:
        return []
    rc, out, err = _sweep_git(repo_root, 'diff', '--name-only', sha_from, sha_to)
    if rc != 0:
        raise ValueError('git diff --name-only %s %s failed: %s'
                         % (sha_from, sha_to, err))
    changed = [p for p in out.splitlines() if p.strip()]
    return sorted(p for p in changed if not _sweep_path_allowed(p))


def default_full_runnable_count():
    """The full runnable cell count DERIVED from this checkout's
    testsys/matrix.py (rule 15 step 1: never a constant in a guard), imported
    from ROOT so this always reads the checkout the guard itself is running
    against, not whatever repo a caller's `repo_root` argument points at."""
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    from testsys import matrix as _matrix
    runnable, _declared_unsupported = _matrix.cells()
    return len(runnable)


def find_sweep_summaries(repo_root):
    """[(path, data_or_None, parse_error_or_None)] for every
    docs/evidence/sweep-*/summary.json under repo_root."""
    paths = sorted(glob.glob(os.path.join(repo_root, EVIDENCE_GLOB)), reverse=True)
    out = []
    for p in paths:
        try:
            with open(p) as fh:
                data = json.load(fh)
        except (OSError, ValueError) as exc:
            out.append((p, None, 'could not parse %s: %s' % (p, exc)))
            continue
        missing = [f for f in SWEEP_REQUIRED_FIELDS if f not in data]
        if missing:
            out.append((p, None, '%s is missing required field(s) %r' % (p, missing)))
            continue
        out.append((p, data, None))
    return out


def evaluate_sweep_candidate(path, data, tag_sha, repo_root, full_runnable_count):
    """Reasons this ONE summary.json fails to justify tagging `tag_sha`
    (empty list means it qualifies)."""
    reasons = []
    sha = data.get('sha')
    if not (isinstance(sha, str) and re.fullmatch(r'[0-9a-fA-F]{40}', sha)):
        return ['%s: "sha" is not a 40-char commit sha (%r)' % (path, sha)]
    if data.get('term') != matrix.GATE_TERM_S:
        reasons.append('%s: term is %r, want %r (matrix.GATE_TERM_S -- there '
                       'is only one term, 2026-09-23)'
                       % (path, data.get('term'), matrix.GATE_TERM_S))
    if data.get('tree_clean') is not True:
        reasons.append('%s: tree_clean is %r, want true' % (path, data.get('tree_clean')))
    cells = data.get('cells')
    if not isinstance(cells, list):
        reasons.append('%s: "cells" is %r, want a list' % (path, cells))
        cells = []
    expected = full_runnable_count()
    n_runnable = data.get('n_runnable')
    n_success = data.get('n_success')
    if n_runnable != expected:
        reasons.append(
            '%s: n_runnable=%r but the tagged tree\'s testsys/matrix.py '
            'currently declares %d runnable cells' % (path, n_runnable, expected))
    if len(cells) != n_runnable:
        reasons.append('%s: %d cell record(s) but n_runnable=%r'
                       % (path, len(cells), n_runnable))
    if n_success != n_runnable:
        reasons.append('%s: n_success=%r != n_runnable=%r' % (path, n_success, n_runnable))
    bad = [c for c in cells if not isinstance(c, dict) or c.get('verdict') != 'SUCCESS']
    if bad:
        names = ['%s/%s=%s' % (c.get('case'), c.get('backend'), c.get('verdict'))
                 if isinstance(c, dict) else repr(c) for c in bad]
        reasons.append('%s: %d cell(s) not SUCCESS: %s' % (path, len(bad), ', '.join(names)))
    try:
        if not _sweep_is_ancestor(repo_root, sha, tag_sha):
            reasons.append(
                '%s: swept sha %s is neither %s itself nor an ancestor of it'
                % (path, sha, tag_sha))
            return reasons
    except Exception as exc:  # pragma: no cover - defensive, git itself failing
        return reasons + ['%s: could not check ancestry of %s: %s' % (path, sha, exc)]
    try:
        disallowed = _sweep_disallowed_paths(repo_root, sha, tag_sha)
    except ValueError as exc:
        return reasons + [str(exc)]
    if disallowed:
        reasons.append(
            '%s: commit(s) between swept sha %s and tagged sha %s touch '
            'file(s) outside the evidence/ledger allow-list: %s'
            % (path, sha, tag_sha, ', '.join(disallowed)))
    return reasons


def evaluate_sweep_evidence(tag_sha, repo_root=None, full_runnable_count=None):
    """(ok, message). PASS iff some committed docs/evidence/sweep-*/summary.json
    justifies tagging `tag_sha`: term==matrix.GATE_TERM_S (the one term every
    cell runs at), tree_clean, every declared cell SUCCESS,
    n_runnable==n_success==len(cells)== the CURRENT matrix.py's runnable count
    (matrix.cells(), i.e. everyday cells + RELEASE_ONLY -- the release
    selection), and the swept sha is `tag_sha` itself or an ancestor whose
    only intervening changes are evidence/ledger paths (rule 15d)."""
    if repo_root is None:
        repo_root = ROOT
    if full_runnable_count is None:
        full_runnable_count = default_full_runnable_count
    found = find_sweep_summaries(repo_root)
    if not found:
        return False, ('no %s found -- a release tag needs a committed local '
                       'release sweep (matrix.GATE_TERM_S) for this exact SHA '
                       '(or a permitted ancestor per rule 15d)' % EVIDENCE_GLOB)
    all_reasons = []
    for path, data, err in found:
        if err:
            all_reasons.append(err)
            continue
        reasons = evaluate_sweep_candidate(path, data, tag_sha, repo_root,
                                           full_runnable_count)
        if not reasons:
            return True, ('%s: swept sha %s justifies tagging %s (n_runnable=%d, '
                         'all SUCCESS)' % (path, data['sha'], tag_sha, data['n_runnable']))
        all_reasons.extend(reasons)
    return False, ('no committed sweep evidence justifies tagging %s:\n  %s'
                   % (tag_sha, '\n  '.join(all_reasons)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pre-tag', metavar='SHA', required=True,
                   help='the commit you are about to `git tag`')
    p.add_argument('--ack-paths-ignored-parent', action='store_true',
                   help='if SHA cannot trigger its own CI run, accept the '
                        'nearest ancestor commit\'s green run as evidence '
                        'instead (prints exactly which SHA that is)')
    args = p.parse_args()

    try:
        requested_full_sha = ci_status.resolve_sha(args.pre_tag)
    except ValueError as exc:
        print('FAIL  %s' % exc)
        return EXIT_CODES[ci_status.FAIL]

    result = ci_status.evaluate_pretag(
        requested_full_sha, ack_paths_ignored_parent=args.ack_paths_ignored_parent)

    print('%s  %s' % (result.status, result.message))
    if result.evidence_sha and result.evidence_sha != requested_full_sha:
        print('  (evidence sha: %s -- NOT %s itself)'
              % (result.evidence_sha, requested_full_sha))
    if result.status != ci_status.PASS:
        return EXIT_CODES[result.status]

    # CI is green for this exact SHA. Additionally required (release-path
    # guard, 2026-09-23): a committed local release sweep (matrix.GATE_TERM_S) for this SHA or a
    # permitted ancestor (rule 15d). Different situation, different code
    # (rule 2) -- do not fold this into the CI outcome above.
    ok, sweep_msg = evaluate_sweep_evidence(requested_full_sha)
    print('%s  %s' % ('PASS' if ok else 'SWEEP_INSUFFICIENT', sweep_msg))
    if not ok:
        return SWEEP_INSUFFICIENT_EXIT
    return 0


if __name__ == '__main__':
    sys.exit(main())
