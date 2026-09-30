#! /usr/bin/env python3
"""
Rule 21c / pathway item 75 -- the MERGE-BOUNDARY half of the enforcement.

`testsys/hooks/pre-commit` refuses a commit whose STAGED set mixes
PROJECT_RULES.md or pathway_forward.md with any other path. That check is
blind in exactly one way, measured 2026-09-23 and not assumed:

    `git commit --amend` compares the index against the commit it REPLACES,
    so amending code into an existing board-only commit (or a board file into
    an existing code commit) shows the hook a clean, unmixed staged set. The
    commit that lands is mixed anyway.

Nothing a pre-commit hook can see distinguishes an amend from a first commit
-- git passes the hook no arguments and sets no marker -- so the gap is closed
HERE instead, over finished history, which is also where the original incident
was actually caught (the conductor reading a merge diff by hand).

Usage, range REQUIRED -- there is no default range, because a wrong default
would silently check the empty set and report success:

    python3 testsys/check_board_separation.py origin/master..HEAD
    python3 testsys/check_board_separation.py origin/master..HEAD --squash-check

Exit 0 if every non-merge commit in the range is a permitted mix (see below)
or board-free; exit 1 naming each offender and the other files it carries.

MERGE COMMITS ARE SKIPPED, and this is a deliberate, stated limit: an
automatic merge introduces no content of its own, and its first-parent diff is
the whole of the other branch, which would flag every merge of a board branch.
An "evil merge" -- a merge commit hand-edited to carry content -- is therefore
not caught by this check. Committed here so nobody reads it as total.

CORRECTED 2026-09-30 (the v5.20.2 release-gate incident, PR #69 -> squash
`ff06534`). PR #69 bundled a code fix (`ci_status.py`, `check_pretag_ci.py`,
`test.yml`) with its OWN matching `PROJECT_RULES.md` rule text (rules 15a/
15g/24/24a) and a `pathway_forward.md` board row (item 139), split across
FOUR individually-clean commits on the branch (code, rules, code, rules).
Each commit, and the whole PR range, passed this check (`pull_request` event,
per-commit). Squash-merging the PR flattened all four into ONE commit on
master, `ff06534`, whose single diff mixes `PROJECT_RULES.md` AND
`pathway_forward.md` with the code files -- and THAT is what the next
`push`-triggered run correctly caught and failed on. The repo's squash-merge-
only setting (rule 25) means a clean multi-commit PR branch is not evidence
the eventual master commit will be clean; only checking the PR's own commits
individually, as before, cannot see this. Two changes, owner-authorized
2026-09-30, close it:

  (1) `PROJECT_RULES.md` may now share a commit with non-board files -- this
      is consilium's own "a rule ships with its refusing check in the same
      change" convention, and rule 21c's separation exists to protect
      REVERTABILITY and single-writer discipline, not to forbid a rule and
      its own enforcing code landing together. `pathway_forward.md` gets NO
      such exception: it is "the board" rule 21c's rationale is really about
      (two sessions editing status rows at once), and it stays solo -- mixing
      it with ANY other file, `PROJECT_RULES.md` included, is still refused.
  (2) A NEW `--squash-check` mode (see `check_squash_aggregate` below):
      treats the WHOLE range as if it were squash-merged into one commit --
      exactly what will happen on merge -- and applies the same separation
      rule to that aggregate path set. Wired into `.github/workflows/test.yml`
      for `pull_request` events only (a `push` is never squashed by this
      repo's settings, so the existing per-commit check already covers it).
      This is what makes the NEXT PR like #69 fail its OWN CI before merge,
      not only master's push-triggered run after.

`ff06534` itself is a one-time, NAMED exemption (`EXEMPT_SHAS` below) rather
than perpetually red: it predates this correction, and even under the
corrected rule it would still be flagged (it carries `pathway_forward.md`
mixed with code, which stays forbidden under (1) above). Recorded, not
silently passed -- the exemption prints as EXEMPT in the output, not SUCCESS.

GENERALIZED THE SAME DAY (owner course correction, after (1)/(2) above had
already landed): the real defect was narrower than "this one check is
blind" -- it is that ANY check judging a PR's shape by its individual
commits, rather than by the ONE diff a squash-merge will actually produce,
can pass a PR whose landed result it never looked at. So `--squash-check`'s
role changed from an ADDITION alongside the per-commit check to the SOLE
GATING result for `pull_request` events (see `main` below): the aggregate
`base...head` diff -- exactly what lands on master -- is what decides
pass/fail; the per-commit result still prints, but as a diagnostic only,
because a squash-merged PR's individual commits never reach master at all.
This is what makes "green PR" and "green master" the SAME claim rather than
two claims that can disagree, which is what actually happened with `ff06534`.
`push` events (never squashed, rule 25) keep the plain per-commit check as
sole and gating, because for THEM the individual commits genuinely are what
lands.

The other two commit/PR-shape checks in this repo were read against this
same principle and need no change: `testsys/pr_policy.py` (the
`pr-policy-gate` job) only runs on `push` to master (never squashed), and
`testsys/change_class.py` classifies a PATH/diff-line pair, not a commit
range, so it has no per-commit-vs-aggregate distinction to get wrong.
"""
import subprocess
import sys

RULE_BOOK_FILE = 'PROJECT_RULES.md'
BOARD_FILE = 'pathway_forward.md'
BOARD_FILES = (RULE_BOOK_FILE, BOARD_FILE)

# One-time, named exemptions for historical commits that predate the
# 2026-09-30 correction and cannot be un-mixed (history is immutable once
# pushed, rule 8-adjacent). Each entry names the incident; this is a ledger,
# not a loophole -- adding to it is exactly as reviewable as any other rule
# change, since it lives in this file under version control.
EXEMPT_SHAS = {
    'ff06534dd741cc53665ad025ece89853ba9f0c60':
        "PR #69 squash-merge incident, 2026-09-30: mixes PROJECT_RULES.md, "
        "pathway_forward.md (item 139), and code in one commit. Predates "
        "this correction; the pathway_forward.md mix would still be refused "
        "under the corrected rule (board stays solo), so this is recorded "
        "as a named exemption, not covered by the new PROJECT_RULES.md "
        "allowance.",
}


def git(args):
    r = subprocess.run(['git'] + args, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise SystemExit('check_board_separation: `git %s` failed (exit %d): %s'
                         % (' '.join(args), r.returncode, r.stderr.strip()))
    return r.stdout


def commit_files(sha):
    """Paths a commit changes against its single parent (root commit: all)."""
    out = git(['show', '--pretty=format:', '--name-only', '-m',
               '--first-parent', sha])
    return [p for p in (line.strip() for line in out.splitlines()) if p]


def range_files(lo, hi):
    """The UNION of paths changed across lo..hi as ONE diff -- what a
    squash-merge of this range into a single commit will actually produce.
    Distinct from the per-commit check: a range whose individual commits are
    each clean can still union to a mixed set (PR #69 / `ff06534`).

    THREE-DOT (`lo...hi`, merge-base diff), never two-dot (victor-reyes
    audit, PR #70, 2026-09-30): for a `pull_request` event, CI passes `lo` =
    `github.event.pull_request.base.sha`, which is MASTER'S CURRENT TIP at
    CI-run time, not the commit the PR branched from -- it moves every time
    something else lands on master while this PR is open. A two-dot diff
    (`git diff lo hi`) compares the two TREES directly, so it would include
    every board-only push that reaches master after the branch point,
    re-appearing as a false "board file in this PR's union". Three-dot
    diffs from `git merge-base lo hi` instead, which is exactly this PR's
    OWN changes regardless of what has since landed on master. Confirmed:
    a code-only branch behind a master-only board push showed
    `pathway_forward.md` in the two-dot union and did NOT in the three-dot
    one."""
    out = git(['diff', '--name-only', '%s...%s' % (lo, hi)])
    return [p for p in (line.strip() for line in out.splitlines()) if p]


def improper_mix(paths):
    """(is_violation, board_paths, other_paths) for one path set.

    `pathway_forward.md` mixed with anything else is ALWAYS a violation --
    the board stays solo, no exception. `PROJECT_RULES.md` mixed with other
    files is NOT a violation on its own (2026-09-30 correction: a rule may
    share a commit with its own enforcing code) -- only flagged when
    `pathway_forward.md` is also present alongside other files.
    """
    board = [p for p in paths if p in BOARD_FILES]
    other = [p for p in paths if p not in BOARD_FILES]
    if not board or not other:
        return False, board, other
    if BOARD_FILE in paths:
        return True, board, other
    # Only RULE_BOOK_FILE among the board files, alongside other files:
    # permitted since 2026-09-30.
    return False, board, other


def check_per_commit(rng):
    """(exit_code) -- the original per-commit check, corrected per
    `improper_mix` above and exempting `EXEMPT_SHAS`."""
    shas = git(['rev-list', '--no-merges', rng]).split()
    offenders = []
    exempted = []
    for sha in shas:
        paths = commit_files(sha)
        violation, board, other = improper_mix(paths)
        if not violation:
            continue
        if sha in EXEMPT_SHAS:
            exempted.append((sha, EXEMPT_SHAS[sha]))
            continue
        subject = git(['log', '-1', '--pretty=format:%s', sha]).strip()
        offenders.append((sha, subject, board, other))

    for sha, reason in exempted:
        print('EXEMPT check_board_separation: %s is a named historical '
              'exemption, not a pass: %s' % (sha[:12], reason))

    if offenders:
        print('FAIL check_board_separation: %d commit(s) in %s MIX '
              'pathway_forward.md with other files (PROJECT_RULES.md rule '
              '21c, pathway item 75) -- PROJECT_RULES.md alone with other '
              'files is permitted since 2026-09-30, pathway_forward.md is '
              'not' % (len(offenders), rng))
        for sha, subject, board, other in offenders:
            print('  %s  %s' % (sha[:12], subject))
            print('      board: %s' % ', '.join(board))
            print('      other: %s' % ', '.join(other))
        print('  Split into two commits -- pathway_forward.md alone, then '
              'the rest -- so a conductor can drop the board-row change '
              'with one `git revert`.')
        return 1
    print('SUCCESS check_board_separation: %d non-merge commit(s) in %s, '
          'none improperly mixing pathway_forward.md with other files (%d '
          'named exemption(s) applied)' % (len(shas), rng, len(exempted)))
    return 0


def check_squash_aggregate(lo, hi):
    """(exit_code) -- treats lo..hi as ONE diff, exactly what this repo's
    squash-merge-only setting (rule 25) will produce on merge. Catches the
    PR #69 shape: every individual commit clean, the union mixed."""
    paths = range_files(lo, hi)
    violation, board, other = improper_mix(paths)
    if violation:
        print('FAIL check_board_separation --squash-check: squashing %s..%s '
              'into one commit (what merging this PR will do) would MIX '
              'pathway_forward.md with other files, even though no '
              'individual commit in the range does -- this is the PR #69 / '
              '`ff06534` incident (2026-09-30). Split the board-row change '
              'out of this PR; land it as its own direct push per rule 25.'
              % (lo, hi))
        print('      board: %s' % ', '.join(board))
        print('      other: %s' % ', '.join(other))
        return 1
    print('SUCCESS check_board_separation --squash-check: %s..%s would not '
          'mix pathway_forward.md with other files if squashed into one '
          'commit' % (lo, hi))
    return 0


def main(argv):
    if len(argv) not in (2, 3) or (len(argv) == 3 and argv[2] != '--squash-check'):
        print(__doc__.strip())
        print('\ncheck_board_separation: REFUSED -- exactly one commit range '
              '(and optionally --squash-check) is required; there is no '
              'default.')
        return 2
    rng = argv[1]
    if len(argv) == 2:
        # push event (never squashed by this repo's settings, rule 25):
        # each commit lands on master exactly as committed, so per-commit
        # purity IS what master will show. Gating.
        return check_per_commit(rng)

    # --squash-check: pull_request event. GENERALIZED 2026-09-30 (owner):
    # "every check about commit or PR shape must judge the PR as it will
    # land: one squashed diff from base...head, evaluated in PR CI. Then a
    # green PR means a green master by construction." The aggregate check
    # is the ONLY thing that gates here -- it is exactly what the post-merge
    # push-triggered check will see, because a squash-merged PR becomes ONE
    # commit on master whose file set is this same union. Per-commit purity
    # among the PR's OWN, soon-to-be-discarded commits is informational
    # only: a PR with messy intermediate commits but a clean squash result
    # must not be refused for a shape nobody will ever see on master, and a
    # PR whose intermediate commits are each clean but whose union is mixed
    # (PR #69 / `ff06534`) must NOT pass just because each commit looked
    # fine on its own.
    if '..' not in rng:
        print('check_board_separation: REFUSED -- --squash-check needs a '
              'lo..hi range, got %r' % rng)
        return 2
    per_commit_code = check_per_commit(rng)
    if per_commit_code != 0:
        print('NOTE check_board_separation: the per-commit result above is '
              'INFORMATIONAL ONLY under --squash-check -- this repo '
              'squash-merges every PR (rule 25), so individual branch '
              'commits never reach master; only the aggregate result below '
              'gates this PR.')
    lo, hi = rng.split('..', 1)
    return check_squash_aggregate(lo, hi)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
