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

Exit 0 if every non-merge commit in the range is board-only or board-free;
exit 1 naming each offender and the other files it carries.

MERGE COMMITS ARE SKIPPED, and this is a deliberate, stated limit: an
automatic merge introduces no content of its own, and its first-parent diff is
the whole of the other branch, which would flag every merge of a board branch.
An "evil merge" -- a merge commit hand-edited to carry content -- is therefore
not caught by this check. Committed here so nobody reads it as total.
"""
import subprocess
import sys

BOARD_FILES = ('PROJECT_RULES.md', 'pathway_forward.md')


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


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip())
        print('\ncheck_board_separation: REFUSED -- exactly one commit range '
              'is required; there is no default.')
        return 2
    rng = argv[1]
    shas = git(['rev-list', '--no-merges', rng]).split()
    offenders = []
    for sha in shas:
        paths = commit_files(sha)
        board = [p for p in paths if p in BOARD_FILES]
        other = [p for p in paths if p not in BOARD_FILES]
        if board and other:
            subject = git(['log', '-1', '--pretty=format:%s', sha]).strip()
            offenders.append((sha, subject, board, other))

    if offenders:
        print('FAIL check_board_separation: %d commit(s) in %s MIX the rule '
              'book / board with other files (PROJECT_RULES.md rule 21c, '
              'pathway item 75)' % (len(offenders), rng))
        for sha, subject, board, other in offenders:
            print('  %s  %s' % (sha[:12], subject))
            print('      board: %s' % ', '.join(board))
            print('      other: %s' % ', '.join(other))
        print('  Split each into two commits -- board alone, then the rest --')
        print('  so a conductor can drop the board change with one `git revert`.')
        return 1
    print('SUCCESS check_board_separation: %d non-merge commit(s) in %s, none '
          'mixing PROJECT_RULES.md/pathway_forward.md with other files'
          % (len(shas), rng))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
