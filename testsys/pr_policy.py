#! /usr/bin/env python3
"""
Decision logic for the owner-approved hybrid PR workflow (2026-09-23):
changes under src/ or testsys/ reach master ONLY through a merged pull
request (squash-merge on GitHub). Everything else (docs, board, evidence,
session logs, rule text, reference artifacts) may still be pushed directly.

ONE COPY of the decision logic, called from two places, on purpose:

  - .github/workflows/test.yml's pr-policy-gate job runs this module in
    ci-check mode on every push to master. It is the AUTHORITATIVE gate --
    it can call GitHub's own commits-pulls API, so it can tell a direct push
    apart from a squash-merge landing.
  - testsys/hooks/pre-push (installed the same way as testsys/hooks/
    pre-commit, via core.hooksPath) runs this module in push-guard mode
    before a local push reaches the remote at all. It needs NO API call: a
    squash/rebase/merge-commit PR landing on GitHub arrives in a local clone
    via fetch, never via push from a developer's own machine, so a LOCAL
    push that carries new src/ or testsys/ content is by construction a
    direct push. push-guard therefore refuses unconditionally the moment a
    new commit in the push range touches a gated path -- no evidence to
    weigh, nothing to call out to.

Both modes share: which paths are "gated" (touches_gated_paths), how a push
range is resolved from a before/after sha pair (resolve_push_range), how a
commit's changed paths are read (commit_files), and the per-commit / whole-
range evaluation (evaluate_commit_gate / evaluate_range). ci-check adds one
thing on top: PR evidence (commit_pr_evidence) that can turn a gated commit
from a violation into an accounted-for PR merge.

EVIDENCE SOURCE FOR ci-check, AND WHY. GitHub's "list pull requests
associated with a commit" endpoint (GET /repos/OWNER/REPO/commits/SHA/pulls)
is used as the AUTHORITATIVE signal: it is GitHub's own record of which
PR(s) a commit belongs to, computed the same way regardless of which of the
three merge strategies (merge commit, squash, rebase) closed the PR, so it
does not need to guess at a commit-subject convention. A squash-merge
subject's "(#NNN)" suffix is used ONLY as a FALLBACK, and only when the API
call itself fails (network/auth/rate-limit) -- never to override what the
API said. It is fallback-only, not primary, for two reasons: (1) it only
exists for GitHub's default squash-merge subject format, so it is blind to
the other two merge strategies and to a squash performed with a custom
subject; (2) a commit subject is developer-writable text with no
verification behind it, so trusting it as primary evidence would let a
hand-typed subject spoof the gate. When the API fails AND no fallback
suffix matches, this module RAISES (PolicyCheckUnavailable) rather than
defaulting to green OR red by guess -- an unverifiable check is a
CI-environment problem to fix, not a policy verdict to report.

push-guard needs none of this: see above, it never calls the API at all.
"""
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request
from collections import namedtuple

# '.github/' added 2026-09-23 by the conductor, loudly: without it one direct
# push could delete the pr-policy-gate job itself, so the owner's "mechanical"
# enforcement would be deletable. The owner named src/ and testsys/; this is
# the one constant to edit to reverse it.
GATED_PREFIXES = ('src/', 'testsys/', '.github/')
MASTER_REF = 'master'
ZERO_SHA = '0' * 40
SQUASH_SUFFIX_RE = re.compile(r'\(#(\d+)\)\s*$')

Evidence = namedtuple('Evidence', 'via_pr source detail')
Decision = namedtuple('Decision', 'sha ok gated_files reason')


class GitCommandError(RuntimeError):
    """A git invocation this module needed failed. Raised, never swallowed
    -- a range this module cannot read is a range it must not silently pass."""


class PolicyCheckUnavailable(RuntimeError):
    """The gate could not be evaluated at all (unresolvable push range, or
    PR evidence unobtainable with no fallback match). Distinct from a policy
    VIOLATION: this says the check itself is broken, not that a rule was
    broken -- but it is still refused, never defaulted to a pass."""


def run_git(args, cwd=None):
    r = subprocess.run(['git'] + args, cwd=cwd, capture_output=True,
                       text=True, timeout=120)
    if r.returncode != 0:
        raise GitCommandError('git %s failed in %r (exit %d): %s'
                              % (' '.join(args), cwd, r.returncode, r.stderr.strip()))
    return r.stdout


def commit_files(sha, cwd=None):
    """Paths a commit changes against its (first) parent; root commit: all."""
    # --no-renames: a rename out of a gated prefix must report the OLD path
    # (git's rename detection prints only the new one, so src/a -> docs/a
    # used to read as ungated). -z: no quoting of non-ASCII paths.
    out = run_git(['show', '--pretty=format:', '--name-only', '--no-renames',
                   '-z', '-m', '--first-parent', sha], cwd=cwd)
    return [p for p in (x.strip() for x in out.split('\0')) if p]


def touches_gated_paths(files, gated_prefixes=GATED_PREFIXES):
    """Sorted subset of files that starts with a gated prefix. Empty means
    the commit is not gated at all (docs, board, evidence, ...)."""
    return sorted({f for f in files
                  for prefix in gated_prefixes if f.startswith(prefix)})


def resolve_push_range(before, head, have_commit):
    """Turn a (before, head) sha pair into a rev-list range spec.

    have_commit(sha) -> bool is injectable (tests pass a stub; real callers
    pass something that shells out to a commit-existence check).

    - before is ZERO_SHA (or empty): the ref is being created for the first
      time. Every ancestor of head is new, so the range is head itself
      (a plain rev-list of head walks the whole history) -- a safe
      INCLUDE-EVERYTHING answer, never narrower than reality.
    - before is unreachable in this clone (force-push rewrote history, or an
      overly shallow fetch): there is NO safe superset fallback available
      here, unlike check_board_separation.py's CI step. That step can fall
      back to origin/master..HEAD because it also runs on branches other
      than master, where origin/master is a genuinely different reference
      point. This check runs ONLY on a push to master, so by the time it
      runs, origin/master already equals head (the push already landed)
      and that fallback would silently resolve to an empty range. Raise
      instead of guessing.
    - otherwise: the ordinary case, before..head.
    """
    if not before or before == ZERO_SHA:
        return head, 'ROOT (ref created; range is every ancestor of head)'
    if not have_commit(before):
        raise PolicyCheckUnavailable(
            "cannot resolve push range: before-sha %r is not reachable in "
            "this clone (force-push, or a too-shallow fetch). There is no "
            "safe superset fallback here (origin/master already equals head "
            "by the time this check runs), so refusing rather than silently "
            "checking nothing." % before)
    return '%s..%s' % (before, head), 'NORMAL (before..head)'


def commit_pr_evidence(sha, subject, api_fetcher):
    """Decide whether sha arrived via a merged PR.

    api_fetcher(sha) -> list[dict] mirrors the commits-pulls API's JSON
    array (each dict a PR; merged_at non-null means merged). It must raise
    on failure (any exception) -- caught here and, if the subject carries a
    squash-merge "(#NNN)" suffix, treated as the fallback signal; otherwise
    re-raised as PolicyCheckUnavailable.
    """
    try:
        prs = api_fetcher(sha)
    except Exception as e:
        m = SQUASH_SUFFIX_RE.search(subject)
        if m:
            return Evidence(True, 'squash-subject-fallback',
                            'commits-api unavailable (%s); subject %r ends '
                            'with (#%s)' % (e, subject, m.group(1)))
        raise PolicyCheckUnavailable(
            "cannot determine PR association for %s: commits-api failed (%s) "
            "and subject %r carries no (#NNN) squash-merge suffix to "
            "fall back on" % (sha, e, subject))
    # A PR counts only if it was merged INTO master AND GitHub's own merge
    # commit for it is THIS commit AND this commit is not the PR's head. The
    # last clause catches the accident this module exists for: an agent
    # fast-forwarding master to an open PR's head, which GitHub then marks
    # "merged" with that head as merge_commit_sha.
    merged, rejected = [], []
    for pr in prs:
        why = []
        if not pr.get('merged_at'):
            why.append('not merged')
        if (pr.get('base') or {}).get('ref') != MASTER_REF:
            why.append('base %r != %r' % ((pr.get('base') or {}).get('ref'), MASTER_REF))
        if pr.get('merge_commit_sha') != sha:
            why.append('merge_commit_sha %r != this commit' % pr.get('merge_commit_sha'))
        if (pr.get('head') or {}).get('sha') == sha:
            why.append('this commit IS the PR head (a direct push of the branch, not a merge)')
        (rejected if why else merged).append((pr.get('number'), why))
    if merged:
        return Evidence(True, 'commits-api',
                        'merged into %s as this commit: PR(s) %s'
                        % (MASTER_REF, sorted(n for n, _ in merged)))
    return Evidence(False, 'commits-api',
                    'no PR merged into %s as this commit (checked %d associated '
                    'PR(s)%s)' % (MASTER_REF, len(prs),
                                 ''.join('; #%s: %s' % (n, ', '.join(w)) for n, w in rejected)))


def evaluate_commit_gate(sha, subject, files, verify_pr,
                         gated_prefixes=GATED_PREFIXES):
    """One commit's verdict. verify_pr(sha, subject) -> Evidence, or None
    to skip PR verification entirely (push-guard mode: a gated commit is
    always a violation, no evidence to weigh -- see module docstring)."""
    gated = touches_gated_paths(files, gated_prefixes)
    if not gated:
        return Decision(sha, True, gated, 'touches no gated path')
    if verify_pr is None:
        return Decision(sha, False, gated,
                        'touches gated path(s) %s; a LOCAL PUSH can never be '
                        'a PR merge (a squash/rebase/merge-commit PR lands '
                        'via fetch from GitHub, not push), so this is a '
                        'direct push and is refused' % gated)
    evidence = verify_pr(sha, subject)
    if evidence.via_pr:
        return Decision(sha, True, gated,
                        'touches gated path(s) %s; PR-merged (%s: %s)'
                        % (gated, evidence.source, evidence.detail))
    return Decision(sha, False, gated,
                    'touches gated path(s) %s WITHOUT an associated merged '
                    'PR (%s: %s)' % (gated, evidence.source, evidence.detail))


def evaluate_range(range_spec, cwd, verify_pr):
    """Every commit on range_spec's first-parent chain, merges included (see
    the comment below). Returns (ok: bool, decisions: list[Decision])."""
    # --first-parent, NOT --no-merges: every commit on master's own chain is
    # evaluated, merges included, each diffed against its first parent. A
    # merge commit carrying its own gated edits (an "evil merge") used to be
    # dropped entirely. Commits reachable only through a merge's second
    # parent are covered by that merge's first-parent diff.
    shas = run_git(['rev-list', '--first-parent', range_spec], cwd=cwd).split()
    decisions = []
    for sha in shas:
        subject = run_git(['log', '-1', '--pretty=format:%s', sha], cwd=cwd).strip()
        files = commit_files(sha, cwd=cwd)
        decisions.append(evaluate_commit_gate(sha, subject, files, verify_pr))
    return all(d.ok for d in decisions), decisions


def github_commits_pulls_fetcher(repo_slug, token):
    """Real api_fetcher for ci-check: GitHub's commits-pulls endpoint,
    stdlib-only (urllib + json -- no requests dependency to add)."""
    if not repo_slug or not token:
        raise PolicyCheckUnavailable(
            'ci-check needs GITHUB_REPOSITORY and a token (GITHUB_TOKEN or '
            'GH_TOKEN) to query the commits-pulls API; repo=%r token_present=%s'
            % (repo_slug, bool(token)))

    def fetch(sha):
        url = 'https://api.github.com/repos/%s/commits/%s/pulls' % (repo_slug, sha)
        req = urllib.request.Request(url, headers={
            'Accept': 'application/vnd.github+json',
            'Authorization': 'Bearer %s' % token,
            'X-GitHub-Api-Version': '2022-11-28',
            'User-Agent': 'EQdyna-pr-policy',
        })
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                if resp.status != 200:
                    raise RuntimeError('HTTP %d from commits-pulls endpoint' % resp.status)
                return json.loads(resp.read().decode('utf-8'))
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError) as e:
            raise RuntimeError('commits-pulls API call failed: %s' % e)
    return fetch


def have_commit(sha, cwd=None):
    r = subprocess.run(['git', 'cat-file', '-e', '%s^{commit}' % sha],
                       cwd=cwd, capture_output=True, timeout=30)
    return r.returncode == 0


def main(argv, cwd=None, fetcher_factory=github_commits_pulls_fetcher):
    if len(argv) != 3 or argv[0] not in ('ci-check', 'push-guard'):
        print(__doc__.strip())
        print('\npr_policy: REFUSED -- usage: pr_policy.py '
             '{ci-check|push-guard} <before-sha> <after-sha>')
        return 2
    mode, before, after = argv
    try:
        range_spec, note = resolve_push_range(
            before, after, have_commit=lambda s: have_commit(s, cwd=cwd))
        print('pr_policy %s: range=%s (%s)' % (mode, range_spec, note))
        verify_pr = None
        if mode == 'ci-check':
            repo_slug = os.environ.get('GITHUB_REPOSITORY')
            token = os.environ.get('GITHUB_TOKEN') or os.environ.get('GH_TOKEN')
            fetch = fetcher_factory(repo_slug, token)
            verify_pr = lambda sha, subject: commit_pr_evidence(sha, subject, fetch)
        ok, decisions = evaluate_range(range_spec, cwd=cwd, verify_pr=verify_pr)
    except (PolicyCheckUnavailable, GitCommandError) as e:
        print('pr_policy %s: REFUSED -- %s' % (mode, e))
        return 1

    gated = [d for d in decisions if d.gated_files]
    print('pr_policy %s: %d commit(s) in range, %d touching a gated path %s'
         % (mode, len(decisions), len(gated), GATED_PREFIXES))
    for d in gated:
        print('  %s %s -- %s' % ('OK ' if d.ok else 'RED', d.sha[:12], d.reason))

    if not ok:
        bad = [d for d in decisions if not d.ok]
        print('\nFAIL pr_policy %s: %d commit(s) touch %s without an '
             'associated merged PR:' % (mode, len(bad), GATED_PREFIXES))
        for d in bad:
            print('  %s  %s' % (d.sha[:12], d.reason))
        if mode == 'push-guard':
            print('\nPush changes under src/ or testsys/ through a pull '
                 'request (squash-merge) instead of directly to master.')
        return 1
    print('PASS pr_policy %s: every gated commit in range (%d total, %d '
         'gated) is a merged pull request, or the range carries no gated '
         'commit at all' % (mode, len(decisions), len(gated)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
