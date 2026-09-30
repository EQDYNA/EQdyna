#! /usr/bin/env python3
"""
Negative test for the release-path guard's RELEASE-DOCS half
(`check_pretag_ci.evaluate_release_docs_ready`, added 2026-09-30 -- owner
requirement (b)/(d): "everything a post-tag check reads (board Tasks-done
row, notes, Release draft) lands before the tag" and "the tag's own CI must
be a re-run of checks already green, never a first run").

RULE 14a, MECHANICAL: a board/gate command must be able to come out BOTH
ways. This file shows the guard RED once (a candidate commit missing its
own Tasks-done row, or its own README News lead, is REFUSED) and then GREEN
(the identical commit, with both documents already in their post-release
shape, is ACCEPTED) -- never asserting only the direction that happens to
be convenient.

WHY A SEPARATE FILE, not folded into `test_pretag_sweep_negative.py` or
`test_pretag_ci_negative.py`. Those two already partition the guard's other
two gates (CI-outcome, sweep-evidence) into dedicated files with their own
fixtures; `test_pretag_ci_negative.py`'s `run_guard` STUBS both
`evaluate_sweep_evidence` and now `evaluate_release_docs_ready` to a fixed
`(True, ...)` for its CI-outcome-only cases (its real fixture SHAs are
mid-development commits with no Tasks-done row for their own VERSION, which
this gate would otherwise correctly refuse -- exactly why it must be
stubbed there). This file is the dedicated coverage for the function that
stub replaces.

Every case drives a REAL `git commit` in a throwaway `tempfile.mkdtemp()`
sandbox repo (same technique as `test_pretag_sweep_negative.py` and
`test_precommit_board_separation_guard.py`) -- never this repository.

Five scenarios:

  1. pathway_forward.md never mentions the commit's own version  -> refuse
  2. pathway_forward.md mentions it, but with no DATED Tasks-done
     row (rule 15 step 4's narrower requirement)                 -> refuse
  3. the Tasks-done row is present, but README's leading News
     block still names the PREVIOUS version                      -> refuse
  4. both documents are in their final, post-release shape        -> ACCEPT
  5. mutation self-check (rule 6): with the row's own regex check
     neutered to always agree, scenario 1's commit is accepted --
     proving scenario 1 refuses BECAUSE of the row check, not by
     accident.

Cheap (rule 9): a handful of `git init`s and commits, well under 1 s. No
network, no build.
"""
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTSYS)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import release_docs  # noqa: E402

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='release-docs guard', GIT_AUTHOR_EMAIL='release-docs@example.invalid',
    GIT_COMMITTER_NAME='release-docs guard', GIT_COMMITTER_EMAIL='release-docs@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)

PATHWAY_WITH_ROW = (
    '# Pathway forward\n\n'
    '| Date | Item | Notes |\n'
    '|---|---|---|\n'
    '| 2026-09-30 | release | v1.2.3 tagged and pushed |\n'
)
PATHWAY_MENTION_NO_ROW = (
    '# Pathway forward\n\n'
    'Some prose that happens to say v1.2.3 without being a dated row.\n'
)
PATHWAY_NO_MENTION = (
    '# Pathway forward\n\nNothing about any release here.\n'
)
README_LEADS_1_2_3 = '# News\n* 20260930 v1.2.3 release notes\n'
README_LEADS_1_2_2 = '# News\n* 20260901 v1.2.2 release notes\n'


def load_guard():
    """The real module, imported fresh from its path (it is not a package)."""
    path = os.path.join(HERE, 'check_pretag_ci.py')
    spec = importlib.util.spec_from_file_location(
        'check_pretag_ci_under_test_release_docs', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def git(args, cwd, check=True):
    r = subprocess.run(['git'] + args, cwd=cwd, capture_output=True, text=True,
                       timeout=60, env=SANDBOX_ENV)
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


def commit_all(cwd, message):
    git(['add', '-A'], cwd=cwd)
    git(['commit', '-q', '-m', message], cwd=cwd)
    return git(['rev-parse', 'HEAD'], cwd=cwd).stdout.strip()


def new_repo(tmp, name, version, pathway_text, readme_text):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    git(['init', '-q', '.'], cwd=d)
    write(d, 'VERSION', version + '\n')
    write(d, 'pathway_forward.md', pathway_text)
    write(d, 'README.md', readme_text)
    sha = commit_all(d, 'candidate commit for %s' % name)
    return d, sha


def check_no_mention_refused(guard, tmp, fails, log):
    d, sha = new_repo(tmp, 'case1_no_mention', '1.2.3',
                      PATHWAY_NO_MENTION, README_LEADS_1_2_2)
    ok, msg = guard.evaluate_release_docs_ready(sha, repo_root=d)
    log.append(('1 no mention of v1.2.3 at all', ok, msg))
    if ok:
        fails.append('1: pathway_forward.md never mentioning v1.2.3 was ACCEPTED')
    if 'pathway_forward.md' not in msg:
        fails.append('1: refusal does not name pathway_forward.md: %r' % msg)


def check_mention_no_row_refused(guard, tmp, fails, log):
    d, sha = new_repo(tmp, 'case2_mention_no_row', '1.2.3',
                      PATHWAY_MENTION_NO_ROW, README_LEADS_1_2_2)
    ok, msg = guard.evaluate_release_docs_ready(sha, repo_root=d)
    log.append(('2 mentions v1.2.3 but no dated row', ok, msg))
    if ok:
        fails.append('2: a bare mention with no dated Tasks-done row was ACCEPTED')
    if 'dated Tasks-done row' not in msg:
        fails.append('2: refusal does not name the missing dated row: %r' % msg)


def check_readme_stale_refused(guard, tmp, fails, log):
    d, sha = new_repo(tmp, 'case3_readme_stale', '1.2.3',
                      PATHWAY_WITH_ROW, README_LEADS_1_2_2)
    ok, msg = guard.evaluate_release_docs_ready(sha, repo_root=d)
    log.append(('3 Tasks-done row present, README still leads v1.2.2', ok, msg))
    if ok:
        fails.append('3: README leading with the PREVIOUS version was ACCEPTED')
    if 'README.md' not in msg or 'v1.2.2' not in msg:
        fails.append('3: refusal does not name the stale README lead: %r' % msg)


def check_both_ready_accepted(guard, tmp, fails, log):
    d, sha = new_repo(tmp, 'case4_both_ready', '1.2.3',
                      PATHWAY_WITH_ROW, README_LEADS_1_2_3)
    ok, msg = guard.evaluate_release_docs_ready(sha, repo_root=d)
    log.append(('4 both documents already post-release', ok, msg))
    if not ok:
        fails.append('4: a commit with both documents already in post-release '
                     'shape was REFUSED: %r' % msg)
    if sha not in msg:
        fails.append('4: acceptance does not name the commit it checked: %r' % msg)


def mutation_self_check(guard, tmp, fails, log):
    """Rule 6: prove case 1 actually fails BECAUSE of
    release_docs.tasks_done_row_present, not because of something else --
    neuter it to always agree and confirm the identical commit then passes,
    then restore it and confirm the refusal returns. This is the file's own
    red-then-green pair for rule 14a, driven twice against ONE commit."""
    d, sha = new_repo(tmp, 'case5_mutation', '1.2.3',
                      PATHWAY_NO_MENTION, README_LEADS_1_2_3)

    saved = release_docs.tasks_done_row_present
    guard.release_docs.tasks_done_row_present = (
        lambda pathway_text, v: (True, 'neutered: always agrees'))
    try:
        ok_neutered, msg_neutered = guard.evaluate_release_docs_ready(sha, repo_root=d)
    finally:
        guard.release_docs.tasks_done_row_present = saved
    log.append(('5a mutation: row check neutered (want PASS)', ok_neutered, msg_neutered))
    if not ok_neutered:
        fails.append('5a: neutering tasks_done_row_present to always agree still '
                     'refused (%r) -- case 1 is not isolating that check' % msg_neutered)

    ok_restored, msg_restored = guard.evaluate_release_docs_ready(sha, repo_root=d)
    log.append(('5b mutation: row check restored (want REFUSE)', ok_restored, msg_restored))
    if ok_restored:
        fails.append('5b: restoring the real row check no longer refused the same '
                     'commit (%r) -- the mutation leaked between calls' % msg_restored)


def main():
    print('Negative test: check_pretag_ci.evaluate_release_docs_ready, 5 scenarios')
    guard = load_guard()
    fails = []
    log = []
    tmp = tempfile.mkdtemp(prefix='release_docs_neg_')
    try:
        check_no_mention_refused(guard, tmp, fails, log)
        check_mention_no_row_refused(guard, tmp, fails, log)
        check_readme_stale_refused(guard, tmp, fails, log)
        check_both_ready_accepted(guard, tmp, fails, log)
        mutation_self_check(guard, tmp, fails, log)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for label, ok, msg in log:
        print('  %-55s ok=%s' % (label, ok))
    if os.environ.get('RELEASE_DOCS_NEG_VERBOSE'):
        for label, ok, msg in log:
            print('--- %s ---\n%s' % (label, msg))

    if fails:
        print('\nFAIL test_pretag_release_docs_negative (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_pretag_release_docs_negative: 3 refusals (no mention, '
          'mention-without-a-dated-row, stale README lead), 1 acceptance (both '
          'documents post-release), and a mutation pair proving the refusal is '
          'the row check\'s (rule 14a: shown red, then green, on the same commit)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
