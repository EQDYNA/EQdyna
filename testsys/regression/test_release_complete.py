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
  4. an annotated git tag `vX.Y.Z` exists locally.

Steps that need the network -- the tag being pushed, CI being green on it, and
the GitHub Release existing -- are checked ONLY when `gh` is available and
authenticated, and are reported as UNVERIFIED (not passed) otherwise. "I could
not check this" and "this is fine" must not share an exit code (rule 2), but
nor should a regression tier fail because a laptop is offline: an unverifiable
network check prints UNVERIFIED and does not affect the exit code, while every
local check is mandatory.

A version still in development trips none of this: the guard SKIPS entirely
unless a tag for VERSION already exists, so it fires on released versions and
stays quiet while VERSION is ahead of the tags.

Cheap (rule 9): file reads plus at most two short `gh` calls.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def version():
    return open(os.path.join(ROOT, 'VERSION'), errors='replace').read().strip()


def _git(*args):
    r = subprocess.run(('git', '-C', ROOT) + args, capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def tag_exists(v):
    rc, out = _git('tag', '-l', 'v' + v)
    return rc == 0 and out.strip() == 'v' + v


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


def main():
    v = version()
    print('Regression guard: release completeness for VERSION %s' % v)
    if not tag_exists(v):
        print('  SKIP  no tag v%s yet -- VERSION is ahead of the tags, so this '
              'is a version in development, not a half-finished release.' % v)
        print('\nSUCCESS test_release_complete (nothing released to check)')
        return 0
    failures = []
    for c in (check_banner_matches, check_readme_news_leads_with_this_version,
              check_pathway_tasks_done_row, check_tag_is_annotated,
              check_network_side):
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
