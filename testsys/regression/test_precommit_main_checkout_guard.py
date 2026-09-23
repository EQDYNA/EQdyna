#! /usr/bin/env python3
"""
Regression guard for pathway_forward.md item 68 (P2) -- PROJECT_RULES.md
rules 21 and 21b.

INCIDENT (relayed 2026-09-22): twice in one day an agent committed directly
into the MAIN checkout /home/utig5/dliu/EQdyna instead of its own worktree,
and a third session had HEAD move underneath it mid-run as a result. The
failure is silent by construction -- the victim's next read of a source file
returns different bytes, no error raised, nothing recording that the tree
moved. Rules 21/21b were hortatory: enforced by being read.

The mechanism that makes them mechanical has THREE parts, and it is dead if
any one of them goes missing, so all three are asserted here:

  1. testsys/hooks/pre-commit is TRACKED by git. An untracked hook exists on
     one box and on no clone.
  2. It is EXECUTABLE -- git index mode 100755 AND the working-tree bit.
     Git silently ignores a non-executable hook: no error, no refusal, and
     the guard reads as passing because nothing happens.
  3. It ACTUALLY REFUSES when --git-dir equals --git-common-dir, and
     ACTUALLY ALLOWS a commit from a linked worktree. The second half is
     load-bearing in the opposite direction: a hook that blocked worktree
     commits would stop every agent in this project, which is a worse
     outage than the one being prevented.
  4. install-eqdyna.sh points core.hooksPath at the tracked directory.
     Git installs no hooks on clone, so an installer that forgets this
     ships a hook that never runs.

Parts 3's behaviour is exercised in a THROWAWAY `git init` sandbox under
tempfile.mkdtemp(), never against this repository: the refusal path has to
be driven by a real `git commit`, and driving it here would mean creating
commits in a live checkout. The sandbox gets its own core.hooksPath; this
repository's shared config is never touched (it is shared by every live
worktree).

Cheap (rule 9): two `git init`s, four commits of a one-line file, no build.
Well under 1 s. Exits non-zero on any failure (rule 2); no skips.
"""
import os
import shutil
import stat
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HOOK_RELPATH = 'testsys/hooks/pre-commit'
HOOK_PATH = os.path.join(ROOT, *HOOK_RELPATH.split('/'))
INSTALLER = os.path.join(ROOT, 'install-eqdyna.sh')
EXPECTED_HOOKSPATH = 'testsys/hooks'

# git env for the sandbox: no user identity is configured in a fresh
# `git init`, and an unset user.email makes `git commit` fail for a reason
# that has nothing to do with the hook -- which would make the refusal test
# pass for the wrong reason.
SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='item68 guard', GIT_AUTHOR_EMAIL='item68@example.invalid',
    GIT_COMMITTER_NAME='item68 guard', GIT_COMMITTER_EMAIL='item68@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def git(args, cwd, check=True):
    r = subprocess.run(['git'] + args, cwd=cwd, capture_output=True,
                       text=True, timeout=60, env=SANDBOX_ENV)
    if check and r.returncode != 0:
        raise RuntimeError('git %s failed in %s (exit %d)\nstdout: %s\nstderr: %s'
                           % (' '.join(args), cwd, r.returncode, r.stdout, r.stderr))
    return r


def check_tracked(fails):
    r = subprocess.run(['git', 'ls-files', '-s', '--', HOOK_RELPATH],
                       cwd=ROOT, capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        fails.append('git ls-files -s failed (exit %d): %s' % (r.returncode, r.stderr))
        return None
    if not r.stdout.strip():
        fails.append(
            '%s is NOT tracked by git. Rules 21/21b are mechanical only if the '
            'hook ships with the repository; an untracked hook exists on one '
            'box and on no clone.' % HOOK_RELPATH)
        return None
    # format: <mode> <sha> <stage>\t<path>
    return r.stdout.split()[0]


def check_executable(index_mode, fails):
    if index_mode is not None and index_mode != '100755':
        fails.append(
            '%s has git index mode %s, not 100755. Git SILENTLY ignores a '
            'non-executable hook -- no error, no refusal -- so a clone would '
            'get an inert guard that reads as installed.'
            % (HOOK_RELPATH, index_mode))
    if not os.path.exists(HOOK_PATH):
        fails.append('%s does not exist in the working tree' % HOOK_RELPATH)
        return
    mode = os.stat(HOOK_PATH).st_mode
    if not (mode & stat.S_IXUSR):
        fails.append('%s is not executable in the working tree (mode %o)'
                     % (HOOK_RELPATH, stat.S_IMODE(mode)))


def check_installer_sets_hookspath(fails):
    if not os.path.exists(INSTALLER):
        fails.append('install-eqdyna.sh not found at %s' % INSTALLER)
        return
    with open(INSTALLER) as fh:
        text = fh.read()
    wanted = 'git config core.hooksPath %s' % EXPECTED_HOOKSPATH
    if wanted not in text:
        fails.append(
            'install-eqdyna.sh does not contain %r. Git installs NO hooks on '
            'clone, so without this the tracked hook never runs and rules '
            '21/21b are back to being enforced by reading them.' % wanted)


def build_sandbox(tmp):
    """A throwaway repo with the REAL hook installed via core.hooksPath.

    Returns (main_checkout, linked_worktree).
    """
    main = os.path.join(tmp, 'main')
    os.makedirs(os.path.join(main, 'testsys', 'hooks'))
    git(['init', '-q', '.'], cwd=main)
    shutil.copyfile(HOOK_PATH, os.path.join(main, *HOOK_RELPATH.split('/')))
    os.chmod(os.path.join(main, *HOOK_RELPATH.split('/')), 0o755)
    git(['add', '-A'], cwd=main)
    # Seed commit with the hook deliberately NOT active yet (hooksPath at a
    # path that does not exist), so the sandbox can be built at all.
    git(['-c', 'core.hooksPath=%s' % os.path.join(tmp, 'no-such-hooks'),
         'commit', '-q', '-m', 'seed'], cwd=main)
    git(['config', 'core.hooksPath', EXPECTED_HOOKSPATH], cwd=main)
    linked = os.path.join(tmp, 'linked-wt')
    git(['worktree', 'add', '-q', '-b', 'wt', linked], cwd=main)
    return main, linked


def commit_a_file(cwd, name, message):
    with open(os.path.join(cwd, name), 'w') as fh:
        fh.write('item68\n')
    git(['add', name], cwd=cwd)
    return git(['commit', '-m', message], cwd=cwd, check=False)


def check_behaviour(fails):
    tmp = tempfile.mkdtemp(prefix='item68_hook_')
    try:
        main, linked = build_sandbox(tmp)

        # --- the refusal: main checkout, --git-dir == --git-common-dir ---
        dirs = git(['rev-parse', '--git-dir', '--git-common-dir'], cwd=main)
        a, b = [os.path.realpath(os.path.join(main, p))
                for p in dirs.stdout.split()]
        if a != b:
            fails.append('sandbox main checkout has DIFFERING git dirs (%s vs '
                         '%s) -- the sandbox is malformed and the refusal test '
                         'below would prove nothing' % (a, b))
            return
        head_before = git(['rev-parse', 'HEAD'], cwd=main).stdout.strip()
        r = commit_a_file(main, 'in_main.txt', 'must be refused')
        if r.returncode == 0:
            fails.append(
                'the hook ALLOWED a commit from the main checkout (git-dir == '
                'git-common-dir). That is exactly the incident rules 21/21b '
                'exist to stop.')
        combined = r.stdout + r.stderr
        if 'REFUSED' not in combined:
            fails.append('the refusal message does not contain "REFUSED": %r'
                         % combined[:400])
        if 'worktree' not in combined:
            fails.append(
                'the refusal message does not tell the committer what to do '
                'instead (no mention of a worktree): %r' % combined[:400])
        head_after = git(['rev-parse', 'HEAD'], cwd=main).stdout.strip()
        if head_after != head_before:
            fails.append('HEAD MOVED in the sandbox main checkout despite the '
                         'refusal (%s -> %s)' % (head_before, head_after))

        # --- the permission: a linked worktree must still be able to commit.
        # A hook that blocks this stops every agent in the project.
        r = commit_a_file(linked, 'in_worktree.txt', 'must succeed')
        if r.returncode != 0:
            fails.append(
                'the hook REFUSED a commit from a LINKED WORKTREE (exit %d): '
                '%r. This is the load-bearing case -- blocking it would stop '
                'every agent session in this project.'
                % (r.returncode, (r.stdout + r.stderr)[:400]))

        # --- and a true fast-forward in the main checkout creates no commit,
        # so it never reaches pre-commit and must not be blocked.
        r = git(['merge', '--ff-only', 'wt'], cwd=main, check=False)
        if r.returncode != 0:
            fails.append(
                'git merge --ff-only was refused in the sandbox main checkout '
                '(exit %d): %r. A true fast-forward creates no commit and must '
                'stay possible -- it is the only permitted HEAD move there.'
                % (r.returncode, (r.stdout + r.stderr)[:400]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    fails = []
    index_mode = check_tracked(fails)
    check_executable(index_mode, fails)
    check_installer_sets_hookspath(fails)
    if index_mode is not None and os.path.exists(HOOK_PATH):
        check_behaviour(fails)
    else:
        fails.append('hook missing or untracked -- its behaviour was not '
                     'exercised, so this run proves nothing about the refusal')

    if fails:
        print('FAIL test_precommit_main_checkout_guard')
        for f in fails:
            print(' -', f)
        return 1
    print('SUCCESS test_precommit_main_checkout_guard: %s is tracked (mode '
          '%s) and executable, install-eqdyna.sh sets core.hooksPath=%s, and '
          'in a throwaway sandbox the hook REFUSED a main-checkout commit '
          '(HEAD unmoved), ALLOWED a linked-worktree commit, and left '
          '`git merge --ff-only` working'
          % (HOOK_RELPATH, index_mode, EXPECTED_HOOKSPATH))
    return 0


if __name__ == '__main__':
    sys.exit(main())
