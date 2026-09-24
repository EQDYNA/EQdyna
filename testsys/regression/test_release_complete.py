#! /usr/bin/env python3
"""
Regression guard: a released version must be released EVERYWHERE (rules 11, 15).

Rule 15's checklist has been followed from memory and skipped steps (the
Tasks-done row, `gh release create`) more than once with a green gate
throughout -- the argument for checking it instead of reading it.

WHAT THIS PINS, for the version currently in VERSION:

  1. the Fortran runtime banner matches VERSION (also covered by
     test_version_banner.py -- kept here because a release is the moment the
     two can drift);
  2. README.md's leading `# News` block is THIS version, not a previous one;
  3. pathway_forward.md has a Tasks-done row naming it;
  4. an annotated git tag `vX.Y.Z` exists locally;
  5. a completed, SUCCESSFUL CI run exists for the exact SHA the tag points
     at (added 2026-09-21, see check_ci_green_for_tagged_sha below -- this is
     the post-hoc half of the pre-tag gate in
     testsys/regression/check_pretag_ci.py; that one runs BEFORE `git tag`
     and blocks it, this one runs AFTER and catches a tag that landed on an
     un-green or never-run SHA anyway);
  6. a committed full-term local sweep (docs/evidence/sweep-*/summary.json)
     justifies the exact tagged SHA (added 2026-09-23, see
     check_sweep_evidence_for_tagged_sha below -- the release-path guard:
     CI stopped running the e2e sweep, so this is now the only mechanical
     physics check at release time, required IN ADDITION to #5, not instead
     of it; post-hoc half of check_pretag_ci.py's --pre-tag guard, exit code
     5, SWEEP_INSUFFICIENT).

Steps that need the network -- the tag being pushed, CI being green on it, and
the GitHub Release existing -- are checked ONLY when `gh` is available and
authenticated, and are reported as UNVERIFIED (not passed) otherwise. "I could
not check this" and "this is fine" must not share an exit code (rule 2), but
nor should a regression tier fail because a laptop is offline: an unverifiable
network check prints UNVERIFIED and does not affect the exit code, while every
local check is mandatory.

A version still in development trips none of the checks above: they SKIP unless
a tag for VERSION already exists, so they fire on released versions and stay
quiet while VERSION is ahead of the tags.

THAT SKIP IS WHAT LET v5.13.0 THROUGH, so one check now runs before it and
always. Keying everything off VERSION tested exactly one direction -- VERSION
ahead of the tags -- and the direction that failed for real was the mirror
image: v5.13.0 was tagged, annotated and pushed 21 commits past v5.12.0 with
NO commit touching VERSION, so the tagged tree printed "Welcome to EQdyna
5.12.0" and this guard was green the whole time, cheerfully re-verifying
v5.12.0, a release that was already complete. A guard that reports on the
version named in a file cannot see a release that the file was never told
about. `check_version_not_behind_newest_tag` compares VERSION against the
highest semver tag REACHABLE FROM HEAD and runs unconditionally, so a dropped
bump fails the cheap tier at the first commit after the tag instead of
surviving to the next release.

Cheap (rule 9): file reads plus at most two short `gh` calls.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from testsys import ci_status  # noqa: E402

import importlib.util as _importlib_util


def _load_check_pretag_ci():
    """The pre-tag guard module, imported from its path (it is not a
    package) -- reused here so the post-hoc sweep-evidence check below and
    the pre-hoc one in check_pretag_ci.py --pre-tag share ONE implementation
    of `evaluate_sweep_evidence` rather than two that could drift apart."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'check_pretag_ci.py')
    spec = _importlib_util.spec_from_file_location(
        'check_pretag_ci_for_release_complete', path)
    mod = _importlib_util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


check_pretag_ci = _load_check_pretag_ci()


def version():
    return open(os.path.join(ROOT, 'VERSION'), errors='replace').read().strip()


def _git(*args):
    r = subprocess.run(('git', '-C', ROOT) + args, capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def tag_exists(v):
    rc, out = _git('tag', '-l', 'v' + v)
    return rc == 0 and out.strip() == 'v' + v


def newest_reachable_tag():
    """The highest semver `vX.Y.Z` tag reachable from HEAD, as a tuple, or
    None if this history carries no release tag at all.

    `--merged HEAD` deliberately, not `git tag -l`: a tag on some branch this
    tree is not descended from says nothing about whether THIS tree declares
    its own version correctly, and failing on one would be a false positive on
    any maintenance branch.
    """
    rc, out = _git('tag', '--merged', 'HEAD', '-l', 'v*')
    if rc != 0:
        return None
    best = None
    for line in out.splitlines():
        m = re.match(r'^v(\d+)\.(\d+)\.(\d+)$', line.strip())
        if not m:
            continue                       # -rc/-dev and other non-releases
        t = tuple(int(g) for g in m.groups())
        if best is None or t > best:
            best = t
    return best


def check_version_not_behind_newest_tag(v):
    """VERSION must not be BEHIND a tag that is already released.

    The direction every other check in this file is blind to. v5.13.0 was
    annotated, pushed and 21 commits past v5.12.0 with no commit touching
    VERSION, so the tagged tree announced itself as 5.12.0 at runtime while
    this guard read VERSION, found v5.12.0 released and complete, and passed.

    Runs before -- and independently of -- the `tag_exists(VERSION)` skip,
    because that skip is precisely the hole: "no tag for VERSION" was read as
    "in development" without ever asking whether a LATER tag exists.
    """
    newest = newest_reachable_tag()
    if newest is None:
        print('  SKIP  no vX.Y.Z tag reachable from HEAD -- nothing released '
              'for VERSION to be behind')
        return
    m = re.match(r'^(\d+)\.(\d+)\.(\d+)$', v)
    if not m:
        raise AssertionError('VERSION is %r, not a bare X.Y.Z semver -- every '
                             'check here and the release tag name derive from '
                             'it' % v)
    mine = tuple(int(g) for g in m.groups())
    newest_s = '%d.%d.%d' % newest
    if mine < newest:
        raise AssertionError(
            'VERSION says %s but v%s is already released and reachable from '
            'HEAD. That release\'s tree declares itself as an OLDER version: '
            'the runtime banner, README and VERSION all name %s, so a v%s '
            'binary misattributes every run to %s. This is the v5.13.0 '
            'failure -- the bump was dropped, not deliberately held, and a '
            'pushed tag cannot be re-pointed (rule 8). Bump VERSION, the '
            'banner and README together (rule 11) and cut the next patch.'
            % (v, newest_s, v, newest_s, v))
    print('  PASS  VERSION %s is not behind the newest reachable tag (v%s)'
          % (v, newest_s))


def check_banner_matches(v):
    p = os.path.join(ROOT, 'src', 'fortran', 'eqdyna3d.f90')
    text = open(p, errors='replace').read()
    m = re.search(r'Welcome to EQdyna ([0-9.]+)', text)
    if not m:
        raise AssertionError('no version banner found in %s' % p)
    if m.group(1) != v:
        raise AssertionError('banner says %s, VERSION says %s -- the banner is '
                             'the first line of every run and the number users '
                             'quote in bug reports' % (m.group(1), v))
    print('  PASS  runtime banner == VERSION (%s)' % v)


def check_readme_news_leads_with_this_version(v):
    p = os.path.join(ROOT, 'README.md')
    text = open(p, errors='replace').read()
    m = re.search(r'^\* \d{8} v([0-9.]+) release notes', text, re.M)
    if not m:
        raise AssertionError("README.md has no '* YYYYMMDD vX.Y.Z release "
                             "notes' block (rule 15 step 3)")
    if m.group(1) != v:
        raise AssertionError(
            "README.md's LEADING News block is v%s but VERSION is %s. Rule 15 "
            "step 3: the current release leads README, and the previous one "
            "moves to pastReleaseNotes.md." % (m.group(1), v))
    print('  PASS  README News block leads with v%s' % v)


def check_pathway_tasks_done_row(v):
    p = os.path.join(ROOT, 'pathway_forward.md')
    text = open(p, errors='replace').read()
    if not re.search(r'v%s\b' % re.escape(v), text):
        raise AssertionError(
            'pathway_forward.md never mentions v%s. Rule 15 step 4 wants a '
            'Tasks-done row for the release; this was skipped on BOTH v5.8.0 '
            'and v5.8.1.' % v)
    # a Tasks-done row, not just a passing mention inside some item body
    if not re.search(r'^\| 20\d\d-\d\d-\d\d \|.*v%s' % re.escape(v), text, re.M):
        raise AssertionError(
            'pathway_forward.md mentions v%s but has no dated Tasks-done row '
            'for it (rule 15 step 4).' % v)
    print('  PASS  pathway_forward.md has a Tasks-done row for v%s' % v)


def check_tag_is_annotated(v):
    """Local cat-file first (works on a normal developer clone). Falls back
    to a remote peel check when that says 'commit' -- confirmed 2026-09-16
    (v5.8.2's tag-triggered Actions run, id 35122388271) that this is NOT
    always a real lightweight tag: actions/checkout fetches a tag-triggered
    event via `git fetch --no-tags ... +<sha>:refs/tags/vX.Y.Z` (its own log
    line), which creates a LOCAL ref pointing straight at the commit
    regardless of what the tag object on the remote actually is. v5.8.2 was
    verified annotated on the remote (`git ls-remote --tags origin` shows
    both `refs/tags/v5.8.2` -> tag object and `refs/tags/v5.8.2^{}` -> its
    peeled commit) at the exact moment this local check said 'commit'. So
    the local answer is checkout-mechanics-dependent, not authoritative;
    the remote's peeled-ref count is what an annotated tag actually is."""
    rc, out = _git('cat-file', '-t', 'v' + v)
    if rc == 0 and out == 'tag':
        print('  PASS  v%s is an annotated tag' % v)
        return
    if rc == 0 and out == 'commit':
        rc2, out2 = _git('ls-remote', '--tags', 'origin',
                         'refs/tags/v' + v, 'refs/tags/v' + v + '^{}')
        if rc2 != 0:
            # NOT the same as "confirmed lightweight" -- rc2!=0 means the
            # remote round-trip itself failed (found 2026-09-16, v5.8.6's
            # verify-published-image run: a container built via `docker
            # build` then run on a LATER, separate job/runner has whatever
            # git credential state got baked into its .git config at build
            # time, which is not guaranteed valid for that later job's
            # network context -- `git ls-remote origin` errored outright,
            # rc=128, not "0 lines" from a real answer). Treating a failed
            # network call as proof of "not annotated" was the bug: this
            # is exactly the class of check check_network_side already
            # handles below with UNVERIFIED, not FAIL, for the identical
            # reason (rule 2: "I could not check this" and "this is fine"
            # must not share an exit code, but an unreachable network is
            # not evidence of failure either).
            print('  UNVERIFIED  v%s tag-annotation status (git ls-remote failed, '
                  'rc=%d -- local cat-file said "commit", ambiguous without the '
                  'remote peel; not a pass, not a fail) -- re-run where the '
                  'remote is reachable with valid credentials' % (v, rc2))
            return
        lines = [l for l in out2.splitlines() if l.strip()]
        if len(lines) == 2:
            print('  PASS  v%s is an annotated tag (local checkout fetched it as '
                  'a bare ref -- confirmed via remote peel instead)' % v)
            return
        raise AssertionError(
            'v%s is not an ANNOTATED tag -- confirmed via the remote: '
            '`git ls-remote --tags origin` returned %d line(s) for it, not the 2 '
            'an annotated tag always shows (the ref plus its peeled `^{}` '
            'commit). Release notes live in the tag body; a lightweight '
            'tag carries none.' % (v, len(lines)))
    raise AssertionError(
        'v%s is not an ANNOTATED tag (git cat-file -t says %r). Release '
        'notes live in the tag body; a lightweight tag carries none.'
        % (v, out))


SWEEP_EVIDENCE_FIRST_VERSION = (5, 17, 0)


def check_sweep_evidence_for_tagged_sha(v):
    """The release-path guard's post-hoc half (2026-09-23, companion to
    check_ci_green_for_tagged_sha immediately below): CI no longer runs the
    e2e sweep, so a committed full-term LOCAL sweep is the only mechanical
    check of the physics for a release, and the tagged sha must be backed by
    one -- shares its logic with check_pretag_ci.py's --pre-tag guard
    (exit code 5, SWEEP_INSUFFICIENT) via `_load_check_pretag_ci` above, so
    the pre-hoc and post-hoc halves cannot silently drift apart."""
    # PROSPECTIVE, by the owner's 2026-09-23 decision: the requirement was
    # adopted after v5.16.2 was tagged and "v5.17.0 becomes the first release
    # gated under the new rule". A tag older than that cannot have been cut
    # against a rule that did not exist; every tag from 5.17.0 on is bound.
    if tuple(int(x) for x in v.split('.')[:3]) < SWEEP_EVIDENCE_FIRST_VERSION:
        print('  N/A   sweep evidence: v%s predates the requirement (binds '
              'from v%s on)' % (v, '.'.join(map(str, SWEEP_EVIDENCE_FIRST_VERSION))))
        return
    try:
        sha = ci_status.resolve_sha('v' + v)
    except ValueError as exc:
        raise AssertionError('cannot resolve tag v%s to a commit: %s' % (v, exc))
    ok, msg = check_pretag_ci.evaluate_sweep_evidence(sha)
    if not ok:
        raise AssertionError(
            'v%s is tagged at %s but no committed full-term local sweep '
            'justifies it: %s' % (v, sha, msg))
    print('  PASS  %s' % msg)


def check_network_side(v):
    """Pushed tag, CI, and the GitHub Release. UNVERIFIED when gh is absent."""
    unverified = []
    rc, out = _git('ls-remote', '--tags', 'origin', 'refs/tags/v' + v)
    if rc != 0:
        unverified.append('tag on remote (git ls-remote failed; offline?)')
    elif ('refs/tags/v' + v) not in out:
        raise AssertionError('v%s is not on the remote -- rule 15 step 6/7' % v)
    else:
        print('  PASS  v%s is pushed to origin' % v)

    if subprocess.run(['which', 'gh'], capture_output=True).returncode != 0:
        unverified.append('GitHub Release (gh not installed)')
    else:
        r = subprocess.run(['gh', 'release', 'view', 'v' + v, '--json', 'tagName'],
                           cwd=ROOT, capture_output=True, text=True)
        if r.returncode != 0:
            if 'release not found' in (r.stderr or '').lower():
                raise AssertionError(
                    'no GitHub Release for v%s. The tag exists but the Releases '
                    'page -- the thing users look at -- does not show it. '
                    'Rule 15 step 9: gh release create v%s --notes-file <notes> '
                    '--latest. This was skipped on both v5.8.0 and v5.8.1.'
                    % (v, v))
            unverified.append('GitHub Release (gh call failed: %s)'
                              % (r.stderr or '').strip()[:80])
        else:
            print('  PASS  GitHub Release exists for v%s' % v)
    for u in unverified:
        print('  UNVERIFIED  %s -- not a pass; re-run where it can be checked' % u)


def check_every_workflow_for_tagged_sha(v, sha, workflow_name, self_run_id):
    """Item 78, post-tag half: EVERY workflow with a run at the tagged sha --
    tag-triggered runs INCLUDED, since after the tag they exist and are the
    only runs publish.yml ever makes -- must be green. v5.16.0 is the case:
    its test.yml runs were green, and 'Publish EQdyna Docker image' run
    35819232914, triggered by the v5.16.0 tag push itself, failed, so no
    image was published and nothing in the release checks said so. The
    pre-tag gate cannot see that run (it does not exist before the tag);
    this is the check that can. A run still in progress is UNVERIFIED, as
    above -- inside the tag's own CI job, publish.yml is usually still
    running."""
    try:
        runs = ci_status.find_all_workflow_runs(sha)
        failed, pending, names = ci_status.classify_every_workflow(
            runs, workflow_name, exclude_run_id=self_run_id)
    except ci_status.GhUnavailable as exc:
        print('  UNVERIFIED  every-workflow CI status for tagged sha %s (%s) -- '
              'not a pass, not a fail' % (sha, exc))
        return
    if failed:
        raise AssertionError(
            'v%s is tagged at %s and workflow(s) %s ran for that sha and '
            'finished WITHOUT success -- %s being green is not CI being green '
            '(item 78; v5.16.0 shipped with no Docker image this way). A pushed '
            'tag cannot be re-pointed (rule 8); this needs a human decision.'
            % (v, sha, ', '.join(repr(n) for n in failed), workflow_name))
    if pending:
        print('  UNVERIFIED  workflow(s) %s for tagged sha %s (v%s) still in '
              'progress -- re-run after they finish'
              % (', '.join(repr(n) for n in pending), sha, v))
        return
    print('  PASS  every workflow with a run at tagged sha %s is green (%d: %s)'
          % (sha, len(names), ', '.join(names)))


def check_ci_green_for_tagged_sha(v):
    """A completed, successful CI run must exist for the exact SHA `vX.Y.Z`
    points at (rule 15 step 6/7's post-hoc half; testsys/regression/
    check_pretag_ci.py is the pre-hoc half that runs BEFORE `git tag`).

    Self-reference: when THIS check runs inside the very CI run that the tag
    push just triggered, `gh run list` will show that run for the tagged SHA
    with status=in_progress -- it cannot be "completed" because it is, right
    now, the thing executing this line. Waiting for it would deadlock the job
    waiting on itself; silently treating "no OTHER completed run yet" as PASS
    would rubber-stamp exactly the ordering the 2026-09-21 incident showed is
    unsafe. So: identify that run by GITHUB_RUN_ID and exclude only it before
    judging the rest -- neither fail nor pass by accident on it.
    """
    try:
        sha = ci_status.resolve_sha('v' + v)
    except ValueError as exc:
        raise AssertionError('cannot resolve tag v%s to a commit: %s' % (v, exc))

    self_run_id = None
    github_run_id = os.environ.get('GITHUB_RUN_ID')
    github_sha = os.environ.get('GITHUB_SHA')
    if github_run_id and github_sha == sha:
        self_run_id = int(github_run_id)

    try:
        workflow_name = ci_status.parse_workflow_name()
        runs = ci_status.find_ci_runs(sha, workflow_name)
    except ci_status.GhUnavailable as exc:
        print('  UNVERIFIED  CI status for tagged sha %s (%s) -- not a pass, '
              'not a fail; re-run where gh can reach the API' % (sha, exc))
        return

    status = ci_status.classify_runs(runs, exclude_run_id=self_run_id)
    if status == 'PASS':
        print('  PASS  a completed, successful %s run exists for tagged sha '
              '%s (v%s)' % (workflow_name, sha, v))
        check_every_workflow_for_tagged_sha(v, sha, workflow_name, self_run_id)
        return
    if status == 'FAIL':
        raise AssertionError(
            'v%s is tagged at %s but the completed %s run(s) for that exact '
            'sha did not succeed -- this is the ordering rule 15 step 6/7 '
            'exists to prevent: a released tag pointing at a red commit. A '
            'pushed tag cannot be re-pointed (rule 8); this needs a human '
            'decision, not a silent pass.' % (v, sha, workflow_name))
    if self_run_id is not None and status in ('IN_PROGRESS', 'NONE'):
        print('  UNVERIFIED  CI status for tagged sha %s (v%s) -- this check '
              'is itself running inside run %d for that sha, which by '
              'definition has not completed yet; excluding it leaves no '
              'other run to judge. Re-run after this job finishes.'
              % (sha, v, self_run_id))
        return
    if status == 'IN_PROGRESS':
        print('  UNVERIFIED  a %s run for tagged sha %s (v%s) is still in '
              'progress -- not provably green yet, not a fail either'
              % (workflow_name, sha, v))
        return
    raise AssertionError(
        'v%s is tagged at %s but no %s run exists for that exact sha at all '
        '-- rule 15 step 6/7 requires CI to have run and gone green on it.'
        % (v, sha, workflow_name))


def main():
    v = version()
    print('Regression guard: release completeness for VERSION %s' % v)
    failures = []

    # UNCONDITIONAL, and before the skip. Every other check reports on the
    # version VERSION names; this one is the only one that can see a release
    # VERSION was never told about, which is what v5.13.0 was.
    try:
        check_version_not_behind_newest_tag(v)
    except AssertionError as e:
        failures.append('check_version_not_behind_newest_tag: %s' % e)
        print('  FAIL  check_version_not_behind_newest_tag: %s' % e)

    if not tag_exists(v):
        print('  SKIP  no tag v%s yet -- VERSION is ahead of every reachable '
              'tag, so this is a version in development, not a half-finished '
              'release. The check above already ruled out the other '
              'direction.' % v)
    else:
        for c in (check_banner_matches,
                  check_readme_news_leads_with_this_version,
                  check_pathway_tasks_done_row, check_tag_is_annotated,
                  check_ci_green_for_tagged_sha,
                  check_sweep_evidence_for_tagged_sha, check_network_side):
            try:
                c(v)
            except AssertionError as e:
                failures.append('%s: %s' % (c.__name__, e))
                print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_release_complete (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_release_complete')
    return 0


if __name__ == '__main__':
    sys.exit(main())
