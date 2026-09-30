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
  4. a completed, SUCCESSFUL CI run exists for the exact SHA the tag points
     at (added 2026-09-21, see check_ci_green_for_tagged_sha below -- this is
     the post-hoc half of the pre-tag gate in
     testsys/regression/check_pretag_ci.py; that one runs BEFORE `git tag`
     and blocks it, this one runs AFTER and catches a tag that landed on an
     un-green or never-run SHA anyway);
  5. a committed full-term local sweep (docs/evidence/sweep-*/summary.json)
     justifies the exact tagged SHA (added 2026-09-23, see
     check_sweep_evidence_for_tagged_sha below -- the release-path guard:
     CI stopped running the e2e sweep, so this is now the only mechanical
     physics check at release time, required IN ADDITION to #4, not instead
     of it; post-hoc half of check_pretag_ci.py's --pre-tag guard, exit code
     5, SWEEP_INSUFFICIENT).

Dropped 2026-09-30 (owner): "an annotated git tag exists" was a pinned check
(`check_tag_is_annotated`) but `gh release create --target <sha>` -- the new
one-step tag+Release command (rule 15 step 7) -- mints a LIGHTWEIGHT tag, not
an annotated one, so the old check would fail by construction under the new
process. Removed entirely, not loosened: there is no replacement check for
tag annotation.

Steps that need the network -- the tag being pushed, CI being green on it, and
the GitHub Release existing -- are checked ONLY when `gh` is available and
authenticated, and are reported as UNVERIFIED (not passed) otherwise. "I could
not check this" and "this is fine" must not share an exit code (rule 2), but
nor should a regression tier fail because a laptop is offline: an unverifiable
network check prints UNVERIFIED and does not affect the exit code, while every
local check is mandatory.

**check_tag_pushed_to_remote vs. check_github_release_published (split
2026-09-30, owner course correction on top of the same day's
regression_sweep_exclusions.py fix).** These used to be one function,
`check_network_side`, called unconditionally by `main()`. That was already
wrong by construction for the GitHub-Release half: any automatic run of
this file (the tag's own CI, or `publish.yml`'s build-and-push image gate)
fires from the SAME `gh release create --target <sha>` command (rule 15
step 7/8, rewritten 2026-09-30) that is what creates the Release in the
first place -- there is no guaranteed ordering between "the Release object
is queryable via the API" and "the ref-creation webhook fired CI" within
that one call, so treating it as certain either way would be a guess, not
a check. Excluding this whole file from the
automatic sweep (`testsys/regression_sweep_exclusions.py`) already stops it
from firing automatically at all, but the function split is kept anyway, as
the mechanical guarantee: `check_tag_pushed_to_remote` runs in `main()`'s
normal check list (a tag being on the remote is true by definition any time
after the tag push, including inside the tag's own CI), while
`check_github_release_published` runs ONLY when `main()` is invoked with
`--post-publish` -- a flag meant to be passed by hand (or by a scripted
step), once, immediately after `gh release create` actually completes.
Default invocation (no flag) SKIPS it outright rather than failing or
passing it by accident.

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
import argparse
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from testsys import ci_status, release_docs  # noqa: E402

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
    """Reuses testsys.release_docs.readme_news_leads_with (ONE copy, shared
    with check_pretag_ci.py's pre-hoc evaluate_release_docs_ready, so the
    two can never drift apart -- see release_docs.py's module docstring)."""
    p = os.path.join(ROOT, 'README.md')
    text = open(p, errors='replace').read()
    ok, msg = release_docs.readme_news_leads_with(text, v)
    if not ok:
        raise AssertionError("README.md %s" % msg)
    print('  PASS  README News block %s' % msg)


def check_pathway_tasks_done_row(v):
    """Reuses testsys.release_docs.tasks_done_row_present (ONE copy, shared
    with check_pretag_ci.py's pre-hoc evaluate_release_docs_ready -- see
    release_docs.py's module docstring)."""
    p = os.path.join(ROOT, 'pathway_forward.md')
    text = open(p, errors='replace').read()
    ok, msg = release_docs.tasks_done_row_present(text, v)
    if not ok:
        raise AssertionError('pathway_forward.md %s' % msg)
    print('  PASS  pathway_forward.md %s' % msg)


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


def check_tag_pushed_to_remote(v):
    """Half of the old check_network_side (split 2026-09-30): is v%s pushed
    to origin. Fine to run ANY time after the tag push, automatically
    included -- unlike check_github_release_published below, this one is
    true by definition the instant the tag push completes, which by
    definition has already happened by the time any CI for that tag fires."""
    rc, out = _git('ls-remote', '--tags', 'origin', 'refs/tags/v' + v)
    if rc != 0:
        print('  UNVERIFIED  tag on remote (git ls-remote failed; offline?) '
              '-- not a pass; re-run where it can be checked')
        return
    if ('refs/tags/v' + v) not in out:
        raise AssertionError('v%s is not on the remote -- rule 15 step 6/7' % v)
    print('  PASS  v%s is pushed to origin' % v)


def check_github_release_published(v):
    """The OTHER half of the old check_network_side (split 2026-09-30, owner
    course correction): does a GitHub Release exist for v%s.

    THIS CHECK IS NOT SAFE to run automatically inside the tag's own CI run,
    or inside `.github/workflows/publish.yml`'s build-and-push image gate:
    both fire off the SAME `gh release create --target <sha>` command
    (rule 15 step 7/8, rewritten 2026-09-30) that mints the Release, and
    there is no guaranteed ordering between "the ref-creation webhook fired
    this CI run" and "the Release object is queryable via the API" within
    one API call -- folding this into the same check list main() always
    runs (the old check_network_side did exactly that, back when a
    separate, later `gh release create` truly could not have run yet) meant
    the Release half failed every time it ran automatically at tag-creation
    time, for a reason that had nothing to do with the tagged commit's
    content (this is what happened to v5.20.2, see
    testsys/regression_sweep_exclusions.py's module docstring).

    Deliberately NOT called from main()'s default check list. Invoke it via
    `--post-publish` on the command line, by hand or from a scripted step,
    ONCE, immediately after `gh release create` has actually completed --
    never automatically on tag push."""
    if subprocess.run(['which', 'gh'], capture_output=True).returncode != 0:
        print('  UNVERIFIED  GitHub Release (gh not installed) -- not a '
              'pass; re-run where it can be checked')
        return
    r = subprocess.run(['gh', 'release', 'view', 'v' + v, '--json', 'tagName'],
                       cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        if 'release not found' in (r.stderr or '').lower():
            raise AssertionError(
                'no GitHub Release for v%s. Rule 15 step 7/8: '
                'gh release create v%s --target <sha> --notes-file <notes> '
                '--latest -- this mints the tag AND the Release together; '
                'if the tag exists but this fails, the one command never '
                'ran (or ran against a different sha). This class of gap '
                '(gh release create skipped) previously hit v5.8.0/v5.8.1.'
                % (v, v))
        print('  UNVERIFIED  GitHub Release (gh call failed: %s) -- not a '
              'pass; re-run where it can be checked'
              % (r.stderr or '').strip()[:80])
        return
    print('  PASS  GitHub Release exists for v%s' % v)


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
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        '--post-publish', action='store_true',
        help='Also run check_github_release_published. Pass this ONLY by '
             'hand (or from a scripted step), ONCE, immediately after '
             '`gh release create --target <sha>` has actually completed '
             '(rule 15 step 7/8). Never pass it automatically from the '
             'tag\'s own CI run or publish.yml\'s build-and-push image gate: '
             'both fire off that SAME command, with no guaranteed ordering '
             'against the Release object becoming queryable, so the '
             'GitHub-Release check is not safe to treat as a pass or a '
             'fail there. Default (flag absent): that check is SKIPPED, not '
             'failed and not passed.')
    args = p.parse_args()

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
                  check_pathway_tasks_done_row,
                  check_ci_green_for_tagged_sha,
                  check_sweep_evidence_for_tagged_sha,
                  check_tag_pushed_to_remote):
            try:
                c(v)
            except AssertionError as e:
                failures.append('%s: %s' % (c.__name__, e))
                print('  FAIL  %s: %s' % (c.__name__, e))

        # check_github_release_published is POST-PUBLICATION ONLY (owner
        # course correction, 2026-09-30): it can never pass automatically at
        # tag-push time (see its own docstring and --post-publish's help
        # text above), so it is never in the loop above. Run it here, by
        # itself, only when explicitly asked.
        if args.post_publish:
            try:
                check_github_release_published(v)
            except AssertionError as e:
                failures.append('check_github_release_published: %s' % e)
                print('  FAIL  check_github_release_published: %s' % e)
        else:
            print('  SKIP  check_github_release_published -- post-publication'
                  '-only check; re-run with --post-publish right after '
                  '`gh release create` completes')
    if failures:
        print('\nFAIL test_release_complete (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_release_complete')
    return 0


if __name__ == '__main__':
    sys.exit(main())
