#! /usr/bin/env python3
"""
Regression guard for testsys/hooks/pre-push -- the LOCAL half of the
owner-approved hybrid PR workflow (2026-09-23): src/ and testsys/ reach
master only through a merged pull request.

Driven end-to-end through a REAL bare remote (`git init --bare`) and a real
local clone in a throwaway tempfile.mkdtemp() sandbox, never against this
repository or its shared core.hooksPath (other sessions are working here).
`git push` is invoked for real; the hook runs for real; only the outcome
(exit code, and the bare remote's actual ref sha) is asserted.

WHAT THIS PINS:
  1. a DOCS-only push to master succeeds, and the bare remote's master
     really advances to the pushed commit;
  2. a push adding a commit that touches ONLY testsys/ directly is REFUSED
     (git push exits non-zero), and the bare remote's master is UNCHANGED --
     not "the hook printed something", the actual side effect never happened;
  3. a push adding a MIXED docs+src commit is likewise refused;
  4. the SAME gated commit, pushed to a branch OTHER than master, succeeds --
     the policy is scoped to master, not every branch;
  5. deleting the master ref (`git push origin :master`) succeeds even
     though the history being removed contains gated commits -- a delete
     introduces no new content, so there is nothing to refuse.

Mutation check for this guard (done while writing it, stated so the next
person does not re-derive it): temporarily changed the hook's
`[ "$remote_ref" = "refs/heads/master" ]` guard to compare against a
nonexistent ref name, ran this file, and cases 2/3 both flipped from
REFUSED to ALLOWED (their remote-sha-unchanged assertions failed) while
case 1/4/5 stayed green -- confirming the guard is exercising the master-
scoping branch and not passing for an unrelated reason. Reverted before
committing.

Cheap (rule 9): git init --bare plus a handful of commits and pushes in a
temp dir, no build, no network. Under ~2 s.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HOOK_RELPATH = 'testsys/hooks/pre-push'
HOOK_PATH = os.path.join(ROOT, *HOOK_RELPATH.split('/'))
POLICY_RELPATH = 'testsys/pr_policy.py'
POLICY_PATH = os.path.join(ROOT, *POLICY_RELPATH.split('/'))
HOOKS_DIRNAME = 'testsys/hooks'

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='prepush guard', GIT_AUTHOR_EMAIL='prepush@example.invalid',
    GIT_COMMITTER_NAME='prepush guard', GIT_COMMITTER_EMAIL='prepush@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def sh(args, cwd, check=True):
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


def commit(cwd, edits, message):
    for name, text in edits.items():
        write(cwd, name, text)
    sh(['add'] + list(edits), cwd=cwd)
    sh(['commit', '-q', '-m', message], cwd=cwd)
    return sh(['rev-parse', 'HEAD'], cwd=cwd).stdout.strip()


def remote_master_sha(bare_dir):
    r = sh(['--git-dir', bare_dir, 'rev-parse', '--verify', '-q', 'refs/heads/master'],
          cwd=bare_dir, check=False)
    return r.stdout.strip() if r.returncode == 0 else None


def build_sandbox(tmp):
    """A bare remote plus a local clone with the REAL hook (and its
    testsys/pr_policy.py dependency) installed via core.hooksPath, seeded
    with one docs-only commit already pushed to master."""
    bare = os.path.join(tmp, 'remote.git')
    os.makedirs(bare)
    sh(['init', '-q', '--bare', '.'], cwd=bare)
    # Allow deleting the bare repo's current branch -- git's OWN default
    # (receive.denyDeleteCurrent=refuse) would otherwise refuse case 5 for a
    # reason that has nothing to do with our hook, and that refusal must not
    # be mistaken for the hook working.
    sh(['config', 'receive.denyDeleteCurrent', 'ignore'], cwd=bare)

    local = os.path.join(tmp, 'local')
    os.makedirs(local)
    sh(['init', '-q', '.'], cwd=local)
    os.makedirs(os.path.join(local, *HOOKS_DIRNAME.split('/')))
    shutil.copyfile(HOOK_PATH, os.path.join(local, *HOOK_RELPATH.split('/')))
    os.chmod(os.path.join(local, *HOOK_RELPATH.split('/')), 0o755)
    os.makedirs(os.path.join(local, 'testsys'), exist_ok=True)
    shutil.copyfile(POLICY_PATH, os.path.join(local, *POLICY_RELPATH.split('/')))
    sh(['config', 'core.hooksPath', HOOKS_DIRNAME], cwd=local)
    sh(['remote', 'add', 'origin', bare], cwd=local)

    seed = commit(local, {'README.md': 'seed\n'}, 'seed: docs only')
    r = sh(['push', 'origin', 'master'], cwd=local, check=False)
    if r.returncode != 0:
        raise RuntimeError('sandbox setup: seed push to master failed:\n%s\n%s'
                           % (r.stdout, r.stderr))
    return bare, local, seed


def out(r):
    return r.stdout + r.stderr


def check_docs_only_push_allowed(bare, local, fails, log):
    before = remote_master_sha(bare)
    sha = commit(local, {'README.md': 'docs edit 2\n'}, 'docs: second edit')
    r = sh(['push', 'origin', 'master'], cwd=local, check=False)
    log.append(('docs-only push', r))
    if r.returncode != 0:
        fails.append('a DOCS-only push to master was REFUSED (exit %d): %r'
                     % (r.returncode, out(r)[:600]))
    after = remote_master_sha(bare)
    if after != sha:
        fails.append('a DOCS-only push claimed success (exit %d) but the bare '
                     'remote master is %r, not the pushed commit %r'
                     % (r.returncode, after, sha))
    if after == before:
        fails.append('the bare remote master did not move at all after an '
                     'allowed push (still %r)' % before)


def check_testsys_push_refused(bare, local, fails, log):
    before = remote_master_sha(bare)
    sha = commit(local, {'testsys/scratch_gate.py': 'x = 1\n'},
                'add a testsys script directly')
    r = sh(['push', 'origin', 'master'], cwd=local, check=False)
    log.append(('testsys/-only push', r))
    text = out(r)
    if r.returncode == 0:
        fails.append('a push adding a testsys/-only commit was ALLOWED '
                     '(exit 0): %r' % text[:600])
    if 'REFUSED' not in text:
        fails.append('the testsys/-only refusal does not contain "REFUSED": %r'
                     % text[:600])
    after = remote_master_sha(bare)
    if after != before:
        fails.append('the bare remote master MOVED (%r -> %r) despite the '
                     'push being refused' % (before, after))
    # leave local in sync with the (unchanged) remote for the next case
    sh(['reset', '-q', '--hard', before], cwd=local)


def check_mixed_push_refused(bare, local, fails, log):
    before = remote_master_sha(bare)
    sha = commit(local, {'README.md': 'mixed docs\n',
                         'src/fortran/scratch.f90': '! stub\n'},
                'mixed: docs plus src, no PR')
    r = sh(['push', 'origin', 'master'], cwd=local, check=False)
    log.append(('mixed docs+src push', r))
    if r.returncode == 0:
        fails.append('a push adding a MIXED docs+src commit was ALLOWED '
                     '(exit 0): %r' % out(r)[:600])
    after = remote_master_sha(bare)
    if after != before:
        fails.append('the bare remote master MOVED (%r -> %r) despite the '
                     'mixed push being refused' % (before, after))
    sh(['reset', '-q', '--hard', before], cwd=local)


def check_non_master_branch_allowed(bare, local, fails, log):
    sh(['checkout', '-q', '-b', 'feature'], cwd=local)
    commit(local, {'testsys/scratch_gate2.py': 'x = 2\n'},
          'add a testsys script on a feature branch')
    r = sh(['push', 'origin', 'feature'], cwd=local, check=False)
    log.append(('testsys/-touching push to a NON-master branch', r))
    if r.returncode != 0:
        fails.append('a testsys/-touching push to refs/heads/feature (not '
                     'master) was REFUSED (exit %d): %r -- the policy is '
                     'scoped to master, not every branch'
                     % (r.returncode, out(r)[:600]))
    sh(['checkout', '-q', 'master'], cwd=local)


def check_master_delete_allowed(bare, local, fails, log):
    r = sh(['push', 'origin', '--delete', 'master'], cwd=local, check=False)
    log.append(('delete refs/heads/master', r))
    if r.returncode != 0:
        fails.append('deleting refs/heads/master was REFUSED (exit %d): %r '
                     '-- a delete introduces no new content and must never '
                     'be refused' % (r.returncode, out(r)[:600]))
    after = remote_master_sha(bare)
    if after is not None:
        fails.append('refs/heads/master still resolves to %r on the bare '
                     'remote after a claimed-successful delete' % after)


def check_hook_mentions_dependency(fails):
    if not os.path.exists(HOOK_PATH):
        fails.append('%s does not exist' % HOOK_RELPATH)
        return
    with open(HOOK_PATH) as fh:
        text = fh.read()
    if POLICY_RELPATH not in text:
        fails.append('%s does not invoke %s -- the decision logic would be '
                     're-implemented in shell instead of shared with CI'
                     % (HOOK_RELPATH, POLICY_RELPATH))
    if not os.access(HOOK_PATH, os.X_OK):
        fails.append('%s is not executable (rule 13)' % HOOK_RELPATH)


def main():
    fails = []
    log = []
    check_hook_mentions_dependency(fails)
    tmp = tempfile.mkdtemp(prefix='prepush_pr_policy_guard_')
    try:
        bare, local, seed = build_sandbox(tmp)
        check_docs_only_push_allowed(bare, local, fails, log)
        check_testsys_push_refused(bare, local, fails, log)
        check_mixed_push_refused(bare, local, fails, log)
        check_non_master_branch_allowed(bare, local, fails, log)
        check_master_delete_allowed(bare, local, fails, log)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if os.environ.get('PREPUSH_VERBOSE'):
        for label, r in log:
            print('--- %s -> exit %d' % (label, r.returncode))
            print(out(r).rstrip())

    if fails:
        print('FAIL test_prepush_pr_policy_guard')
        for f in fails:
            print(' -', f)
        return 1
    print('SUCCESS test_prepush_pr_policy_guard: against a real bare remote, '
         'a docs-only push to master succeeded and the remote sha advanced, '
         'a testsys/-only push and a mixed docs+src push were both REFUSED '
         'with the remote sha UNCHANGED, the same gated content pushed to a '
         'non-master branch succeeded, and deleting refs/heads/master '
         'succeeded despite gated history behind it')
    return 0


if __name__ == '__main__':
    sys.exit(main())
