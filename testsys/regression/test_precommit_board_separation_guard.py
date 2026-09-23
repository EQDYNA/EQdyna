#! /usr/bin/env python3
"""
Regression guard for pathway_forward.md item 75 (P3) -- PROJECT_RULES.md
rule 21c.

INCIDENT (2026-09-23): two subagents dispatched in one night's autopilot each
edited PROJECT_RULES.md and pathway_forward.md unasked, from inside their own
CODE missions. The conductor stripped both files back before landing, so
nothing reached master -- but catching it required reading the merge diff by
hand, and the next conductor may not.

Rule 21c says those two files have exactly one writer. Git cannot enforce
that: a hunk carries no session identity, and nothing in a commit proves who
typed it. So `testsys/hooks/pre-commit` enforces the strongest thing git CAN
see -- SEPARATION. A commit whose staged set contains either board file
TOGETHER WITH any path outside those two is refused, which makes every board
change arrive as its own reviewable commit that one `git revert` removes.

This guard asserts the four behaviours that matter, each driven by a REAL
`git commit` in a throwaway `tempfile.mkdtemp()` sandbox repo -- never against
this repository, and never touching this repository's shared `core.hooksPath`
(it is shared by every live worktree, and other sessions are working here):

  1. a BOARD-ONLY commit is ALLOWED;
  2. a CODE-ONLY commit is ALLOWED;
  3. a MIXED commit is REFUSED, and the message NAMES the offending file and
     says to split it into two commits;
  4. the main-checkout refusal (item 68, rules 21/21b) STILL FIRES -- the new
     check is second, and a reordering that let a main-checkout commit through
     would be the older and worse regression.

It also MEASURES, rather than assuming, the three git operations that would
jam this project's merge path if the hook fired on them. Each is asserted, so
a future git version that changes the answer fails this guard instead of
silently blocking merges:

  5. `git merge --ff-only` creates no commit -> hook never runs;
  6. an AUTOMATIC merge commit from `git merge` of two branches -- measured
     below. If git ran `pre-commit` for it, merging a board branch into a code
     branch would be REFUSED, because the merge commit's staged set is the
     union of both sides. The measurement is the point of the test;
  7. `git commit --amend` on a board-only commit stays allowed.

MEASURED 2026-09-23 with git in this environment: (5) exit 0, no commit
created, hook never invoked; (6) exit 0 -- git does NOT run `pre-commit` for
an automatic merge commit, so merging a board branch with a code branch is
UNAFFECTED, and that is asserted rather than assumed because if it ever
changes it jams the project's merge path; (7) exit 0.

The one BLIND SPOT, also measured and asserted: an amend that ADDS code to an
existing board-only commit is allowed, because `git diff --cached` compares
against the commit being replaced. Git gives pre-commit no way to detect an
amend, so that gap is closed at the merge boundary by
`testsys/check_board_separation.py`, and this guard asserts that checker
catches exactly the commit the amend produced.

Cheap (rule 9): a handful of `git init`s and one-line commits, no build.
Under ~2 s. Exits non-zero on any failure (rule 2); no skips.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HOOK_RELPATH = 'testsys/hooks/pre-commit'
HOOK_PATH = os.path.join(ROOT, *HOOK_RELPATH.split('/'))
HOOKS_DIRNAME = 'testsys/hooks'
BOARD_FILES = ('PROJECT_RULES.md', 'pathway_forward.md')

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='item75 guard', GIT_AUTHOR_EMAIL='item75@example.invalid',
    GIT_COMMITTER_NAME='item75 guard', GIT_COMMITTER_EMAIL='item75@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def git(args, cwd, check=True):
    r = subprocess.run(['git'] + args, cwd=cwd, capture_output=True,
                       text=True, timeout=60, env=SANDBOX_ENV)
    if check and r.returncode != 0:
        raise RuntimeError('git %s failed in %s (exit %d)\nstdout: %s\nstderr: %s'
                           % (' '.join(args), cwd, r.returncode, r.stdout, r.stderr))
    return r


def write(cwd, name, text):
    path = os.path.join(cwd, name)
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, 'w') as fh:
        fh.write(text)


def build_sandbox(tmp):
    """Throwaway repo with the REAL hook active, plus a linked worktree.

    The sandbox carries both board files and one source file, so every case
    below is a commit of the same SHAPE as the real thing.
    """
    main = os.path.join(tmp, 'main')
    os.makedirs(main)
    git(['init', '-q', '.'], cwd=main)
    os.makedirs(os.path.join(main, *HOOKS_DIRNAME.split('/')))
    shutil.copyfile(HOOK_PATH, os.path.join(main, *HOOK_RELPATH.split('/')))
    os.chmod(os.path.join(main, *HOOK_RELPATH.split('/')), 0o755)
    write(main, 'PROJECT_RULES.md', 'seed rules\n')
    write(main, 'pathway_forward.md', 'seed board\n')
    write(main, 'src/solver.txt', 'seed code\n')
    git(['add', '-A'], cwd=main)
    # Seed with the hook deliberately inactive, so the sandbox can be built.
    git(['-c', 'core.hooksPath=%s' % os.path.join(tmp, 'no-such-hooks'),
         'commit', '-q', '-m', 'seed'], cwd=main)
    git(['config', 'core.hooksPath', HOOKS_DIRNAME], cwd=main)
    linked = os.path.join(tmp, 'linked-wt')
    git(['worktree', 'add', '-q', '-b', 'wt', linked], cwd=main)
    return main, linked


def stage_and_commit(cwd, edits, message):
    for name, text in edits.items():
        write(cwd, name, text)
    git(['add'] + list(edits), cwd=cwd)
    return git(['commit', '-m', message], cwd=cwd, check=False)


def out(r):
    return r.stdout + r.stderr


def check_board_only_allowed(wt, fails, log):
    r = stage_and_commit(wt, {'pathway_forward.md': 'board edit 1\n'},
                         'board: item 75 row')
    log.append(('board-only commit', r))
    if r.returncode != 0:
        fails.append(
            'a BOARD-ONLY commit was REFUSED (exit %d): %r. Rule 21c permits '
            'the owning session to update the board; a hook that blocks it '
            'blocks the only legitimate way to write these files.'
            % (r.returncode, out(r)[:600]))


def check_both_board_files_allowed(wt, fails, log):
    r = stage_and_commit(
        wt, {'PROJECT_RULES.md': 'rules edit\n',
             'pathway_forward.md': 'board edit 2\n'},
        'board: rule 21c plus its row')
    log.append(('board-only commit, BOTH files', r))
    if r.returncode != 0:
        fails.append(
            'a commit touching BOTH board files and nothing else was REFUSED '
            '(exit %d): %r. A rule and the row that tracks it are one change; '
            'refusing that pair would force a nonsense split.'
            % (r.returncode, out(r)[:600]))


def check_code_only_allowed(wt, fails, log):
    r = stage_and_commit(wt, {'src/solver.txt': 'code edit\n'},
                         'code: touch the solver')
    log.append(('code-only commit', r))
    if r.returncode != 0:
        fails.append('a CODE-ONLY commit was REFUSED (exit %d): %r'
                     % (r.returncode, out(r)[:600]))


def check_mixed_refused(wt, fails, log):
    head_before = git(['rev-parse', 'HEAD'], cwd=wt).stdout.strip()
    r = stage_and_commit(
        wt, {'pathway_forward.md': 'board edit 3\n',
             'src/solver.txt': 'code edit 2\n'},
        'mixed: must be refused')
    log.append(('MIXED commit', r))
    text = out(r)
    if r.returncode == 0:
        fails.append(
            'the hook ALLOWED a commit mixing pathway_forward.md with '
            'src/solver.txt. That is exactly the 2026-09-23 incident shape '
            'and the whole point of item 75.')
    if 'REFUSED' not in text:
        fails.append('the mixing refusal does not contain "REFUSED": %r'
                     % text[:600])
    if 'src/solver.txt' not in text:
        fails.append(
            'the mixing refusal does not NAME the offending file '
            '(src/solver.txt), so the committer cannot tell what to remove: %r'
            % text[:600])
    if 'pathway_forward.md' not in text:
        fails.append('the mixing refusal does not name the board file: %r'
                     % text[:600])
    if 'two commits' not in text:
        fails.append(
            'the mixing refusal does not say what to do (no "two commits"): %r'
            % text[:600])
    head_after = git(['rev-parse', 'HEAD'], cwd=wt).stdout.strip()
    if head_after != head_before:
        fails.append('HEAD MOVED despite the mixing refusal (%s -> %s)'
                     % (head_before, head_after))
    # leave the sandbox clean for the cases that follow
    git(['restore', '--staged', '--worktree',
         'pathway_forward.md', 'src/solver.txt'], cwd=wt)


def check_main_checkout_still_refused(main, fails, log):
    """Item 68's refusal must still fire -- and must fire FIRST.

    A mixed commit attempted from the main checkout is refused for the OLDER
    reason; if the new check had been put ahead of it, a main-checkout
    board-only commit would sail through.
    """
    head_before = git(['rev-parse', 'HEAD'], cwd=main).stdout.strip()
    r = stage_and_commit(main, {'pathway_forward.md': 'from main checkout\n'},
                         'board-only, but from the main checkout')
    log.append(('main-checkout commit (board-only)', r))
    text = out(r)
    if r.returncode == 0:
        fails.append(
            'the hook ALLOWED a commit from the MAIN CHECKOUT. Item 68\'s '
            'refusal has regressed -- probably because the rule-21c check was '
            'placed ahead of it and returned 0.')
    if 'MAIN CHECKOUT' not in text:
        fails.append(
            'a main-checkout commit was refused, but NOT by the item 68 '
            'check (message lacks "MAIN CHECKOUT"): %r' % text[:600])
    if git(['rev-parse', 'HEAD'], cwd=main).stdout.strip() != head_before:
        fails.append('HEAD moved in the sandbox main checkout despite refusal')
    git(['restore', '--staged', '--worktree', 'pathway_forward.md'], cwd=main)


def measure_merge_and_amend(tmp, fails, log, measured):
    """Three git operations that must NOT be jammed by this hook.

    Each is MEASURED in its own sandbox and then asserted, so a git version
    that starts firing pre-commit for merges fails here loudly instead of
    blocking the project's merge path silently.
    """
    main, linked = build_sandbox(os.path.join(tmp, 'merge'))

    # (5) fast-forward: creates no commit, so pre-commit cannot run.
    git(['checkout', '-q', '-b', 'ffsrc'], cwd=linked)
    r = stage_and_commit(linked, {'src/solver.txt': 'ff code\n'}, 'ff: code')
    if r.returncode != 0:
        raise RuntimeError('sandbox setup: ff branch commit failed: %s' % out(r))
    r = git(['merge', '--ff-only', 'ffsrc'], cwd=main, check=False)
    log.append(('git merge --ff-only in the MAIN CHECKOUT', r))
    measured['ff_only_exit'] = r.returncode
    if r.returncode != 0:
        fails.append(
            'git merge --ff-only was REFUSED (exit %d): %r. A true '
            'fast-forward creates NO commit, so pre-commit must never see it; '
            'this is the only permitted HEAD move in the main checkout.'
            % (r.returncode, out(r)[:600]))

    # (6) an AUTOMATIC merge commit whose union of sides MIXES board and code.
    # Built in a fresh sandbox so it runs in a linked worktree (item 68's
    # refusal would otherwise mask the answer).
    main2, wt2 = build_sandbox(os.path.join(tmp, 'merge2'))
    git(['checkout', '-q', '-b', 'boardbr'], cwd=wt2)
    r = stage_and_commit(wt2, {'pathway_forward.md': 'board branch\n'},
                         'board: branch edit')
    if r.returncode != 0:
        raise RuntimeError('sandbox setup: board branch commit failed: %s' % out(r))
    git(['checkout', '-q', '-b', 'codebr', 'wt'], cwd=wt2)
    r = stage_and_commit(wt2, {'src/solver.txt': 'code branch\n'},
                         'code: branch edit')
    if r.returncode != 0:
        raise RuntimeError('sandbox setup: code branch commit failed: %s' % out(r))
    r = git(['merge', '--no-ff', '--no-edit', 'boardbr'], cwd=wt2, check=False)
    log.append(('git merge (AUTOMATIC merge commit, board branch into code '
                'branch)', r))
    measured['auto_merge_exit'] = r.returncode
    merged_ok = (r.returncode == 0)
    if merged_ok:
        parents = git(['rev-list', '--parents', '-n', '1', 'HEAD'],
                      cwd=wt2).stdout.split()
        measured['auto_merge_parents'] = len(parents) - 1
        if len(parents) - 1 != 2:
            fails.append('the merge did not produce a 2-parent commit '
                         '(parents=%d); the case measured nothing'
                         % (len(parents) - 1))
    else:
        measured['auto_merge_parents'] = 0
        fails.append(
            'BLOCKER: git ran pre-commit for an AUTOMATIC merge commit and '
            'the hook REFUSED it (exit %d): %r. A merge whose two sides are a '
            'board branch and a code branch has a union staged set, so this '
            'hook would block the project\'s entire merge path. The hook must '
            'detect MERGE_HEAD and exit 0 before the rule-21c check.'
            % (r.returncode, out(r)[:600]))

    # (7) `git commit --amend` on a board-only commit. Amend DOES run
    # pre-commit; the staged set it sees must still be board-only.
    git(['checkout', '-q', 'boardbr'], cwd=wt2)
    r = git(['commit', '--amend', '--no-edit'], cwd=wt2, check=False)
    log.append(('git commit --amend on a BOARD-ONLY commit', r))
    measured['amend_board_exit'] = r.returncode
    if r.returncode != 0:
        fails.append(
            'git commit --amend on a BOARD-ONLY commit was REFUSED (exit %d): '
            '%r. Amending a board commit is a board-only operation and must '
            'stay allowed.' % (r.returncode, out(r)[:600]))

    # and the mirror, which is the hook's ONE measured blind spot: amending
    # code into an existing board-only commit. `git diff --cached` compares
    # the index against the commit being REPLACED, so the hook sees a clean
    # code-only staged set and allows it, while the commit that lands is
    # mixed. Git passes pre-commit no argument and sets no marker that
    # distinguishes an amend, so this cannot be closed in the hook -- it is
    # closed at the merge boundary by testsys/check_board_separation.py, and
    # BOTH halves of that statement are asserted here.
    git(['checkout', '-q', 'boardbr'], cwd=wt2)
    base = git(['rev-parse', 'HEAD~1'], cwd=wt2).stdout.strip()
    write(wt2, 'src/solver.txt', 'amended in code\n')
    git(['add', 'src/solver.txt'], cwd=wt2)
    r = git(['commit', '--amend', '--no-edit'], cwd=wt2, check=False)
    log.append(('git commit --amend turning a board commit MIXED', r))
    measured['amend_mixed_exit'] = r.returncode
    if r.returncode != 0:
        fails.append(
            'the hook REFUSED an amend that turns a board-only commit mixed '
            '(exit %d). That is a BETTER outcome than measured on 2026-09-23, '
            'but this guard and check_board_separation.py document the '
            'opposite; re-measure and update both before trusting either.'
            % r.returncode)
    return base, wt2


def check_range_checker_catches_amend(base, wt2, fails, log, measured):
    """The merge-boundary half: it must FLAG the mixed commit the amend made,
    PASS a range that holds only clean commits, and REFUSE to run with no
    range rather than silently checking the empty set."""
    checker = os.path.join(ROOT, 'testsys', 'check_board_separation.py')
    if not os.path.exists(checker):
        fails.append('testsys/check_board_separation.py is missing, so the '
                     'amend gap measured above is not closed anywhere')
        return

    def run(args):
        return subprocess.run([sys.executable, checker] + args, cwd=wt2,
                              capture_output=True, text=True, timeout=120,
                              env=SANDBOX_ENV)

    r = run(['%s..HEAD' % base])
    log.append(('check_board_separation over the amended MIXED commit', r))
    measured['range_mixed_exit'] = r.returncode
    if r.returncode == 0:
        fails.append(
            'check_board_separation PASSED a range containing a commit that '
            'mixes pathway_forward.md with src/solver.txt. The amend gap is '
            'then closed nowhere and rule 21c is one `--amend` from bypassed.')
    if 'src/solver.txt' not in out(r):
        fails.append('check_board_separation did not NAME the offending other '
                     'file: %r' % out(r)[:600])

    r = run(['%s..%s' % (base, base)])
    log.append(('check_board_separation over an empty/clean range', r))
    measured['range_clean_exit'] = r.returncode
    if r.returncode != 0:
        fails.append('check_board_separation failed on a clean range (exit '
                     '%d): %r' % (r.returncode, out(r)[:600]))

    r = run([])
    log.append(('check_board_separation with NO range', r))
    measured['range_norange_exit'] = r.returncode
    if r.returncode == 0:
        fails.append(
            'check_board_separation with no range exited 0. A missing range '
            'must REFUSE, not fall back to checking nothing and reporting '
            'success (PROJECT_RULES rule 2).')


def check_hook_mentions_board_files(fails):
    """Item 75's own verification command: the hook must know these names."""
    if not os.path.exists(HOOK_PATH):
        fails.append('%s does not exist' % HOOK_RELPATH)
        return
    with open(HOOK_PATH) as fh:
        text = fh.read()
    for name in BOARD_FILES:
        if name not in text:
            fails.append(
                '%s does not mention %s, so rule 21c is not mechanical '
                '(this is item 75\'s own verification command).'
                % (HOOK_RELPATH, name))


def main():
    fails = []
    log = []
    measured = {}
    check_hook_mentions_board_files(fails)
    tmp = tempfile.mkdtemp(prefix='item75_board_sep_')
    try:
        os.makedirs(os.path.join(tmp, 'merge'))
        os.makedirs(os.path.join(tmp, 'merge2'))
        main_ck, wt = build_sandbox(tmp)
        check_board_only_allowed(wt, fails, log)
        check_both_board_files_allowed(wt, fails, log)
        check_code_only_allowed(wt, fails, log)
        check_mixed_refused(wt, fails, log)
        check_main_checkout_still_refused(main_ck, fails, log)
        base, wt2 = measure_merge_and_amend(tmp, fails, log, measured)
        check_range_checker_catches_amend(base, wt2, fails, log, measured)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if os.environ.get('ITEM75_VERBOSE'):
        for label, r in log:
            print('--- %s -> exit %d' % (label, r.returncode))
            print(out(r).rstrip())

    if fails:
        print('FAIL test_precommit_board_separation_guard')
        for f in fails:
            print(' -', f)
        return 1
    print('SUCCESS test_precommit_board_separation_guard: in a throwaway '
          'sandbox the hook ALLOWED board-only (either file, or both), '
          'ALLOWED code-only, REFUSED a mixed commit naming src/solver.txt '
          'and HEAD unmoved, still REFUSED a main-checkout commit via item '
          '68\'s check; measured git merge --ff-only exit %s, AUTOMATIC merge '
          'commit exit %s (%s parents), commit --amend board-only exit %s, '
          'amend-to-mixed ALLOWED by the hook (exit %s, git\'s one blind '
          'spot) and CAUGHT by check_board_separation.py (exit %s, clean '
          'range %s, no-range %s)'
          % (measured.get('ff_only_exit'), measured.get('auto_merge_exit'),
             measured.get('auto_merge_parents'),
             measured.get('amend_board_exit'),
             measured.get('amend_mixed_exit'),
             measured.get('range_mixed_exit'),
             measured.get('range_clean_exit'),
             measured.get('range_norange_exit')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
