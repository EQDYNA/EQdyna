#! /usr/bin/env python3
"""
Unit coverage for testsys/check_board_separation.py's 2026-09-30 correction
(PR #69 / commit `ff06534` incident).

Three things this file proves, none covered by the existing regression
guards (`test_precommit_board_separation_guard.py` tests the pre-commit
HOOK; `test_ci_board_separation_step.py` tests that test.yml WIRES the step
up -- neither exercises this script's own per-commit/aggregate logic):

  1. `PROJECT_RULES.md` mixed with other files is no longer an offender;
     `pathway_forward.md` mixed with other files still is, unconditionally
     (rule 21c's post-correction text: "the board stays separate").
  2. `EXEMPT_SHAS` removes a named historical commit from the offender list
     without silently passing it -- the output says EXEMPT, not SUCCESS.
  3. `--squash-check` reproduces the PR #69 shape: a range whose individual
     commits are each clean (PROJECT_RULES.md-only, code-only,
     pathway_forward.md-only, code-only) still fails once unioned into one
     diff, because the union mixes pathway_forward.md with code -- exactly
     what squash-merging the range would produce. The per-commit check on
     the SAME range must still pass, proving the two modes catch different
     things.

Real git sandboxes (tempfile.mkdtemp), never this repository's own tree or
its shared core.hooksPath.
"""
import os
import subprocess
import sys
import tempfile

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'testsys'))
import check_board_separation as cbs  # noqa: E402

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='cbs unit', GIT_AUTHOR_EMAIL='cbs@example.invalid',
    GIT_COMMITTER_NAME='cbs unit', GIT_COMMITTER_EMAIL='cbs@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def git(args, cwd, check=True):
    r = subprocess.run(['git'] + args, cwd=cwd, capture_output=True,
                       text=True, timeout=60, env=SANDBOX_ENV)
    if check and r.returncode != 0:
        raise RuntimeError('git %s failed in %s (exit %d)\n%s\n%s'
                           % (' '.join(args), cwd, r.returncode, r.stdout, r.stderr))
    return r


def write(cwd, name, text):
    with open(os.path.join(cwd, name), 'w') as fh:
        fh.write(text)


def commit(cwd, edits, message):
    for name, text in edits.items():
        write(cwd, name, text)
    git(['add'] + list(edits), cwd=cwd)
    git(['commit', '-q', '-m', message], cwd=cwd)
    return git(['rev-parse', 'HEAD'], cwd=cwd).stdout.strip()


@pytest.fixture
def repo():
    tmp = tempfile.mkdtemp(prefix='cbs_unit_')
    git(['init', '-q', '.'], cwd=tmp)
    write(tmp, 'PROJECT_RULES.md', 'seed rules\n')
    write(tmp, 'pathway_forward.md', 'seed board\n')
    os.makedirs(os.path.join(tmp, 'src'), exist_ok=True)
    write(tmp, 'src/solver.txt', 'seed code\n')
    git(['add', '-A'], cwd=tmp)
    git(['commit', '-q', '-m', 'seed'], cwd=tmp)
    yield tmp


def run_checker(args, cwd):
    return subprocess.run([sys.executable, os.path.join(ROOT, 'testsys',
                          'check_board_separation.py')] + args,
                         cwd=cwd, capture_output=True, text=True,
                         timeout=60, env=SANDBOX_ENV)


# --------------------------------------------------- (1) improper_mix() ---

def test_rulebook_alone_with_other_is_not_a_violation():
    violation, board, other = cbs.improper_mix(['PROJECT_RULES.md', 'src/solver.txt'])
    assert violation is False
    assert board == ['PROJECT_RULES.md']
    assert other == ['src/solver.txt']


def test_board_row_alone_with_other_is_a_violation():
    violation, board, other = cbs.improper_mix(['pathway_forward.md', 'src/solver.txt'])
    assert violation is True


def test_both_board_files_with_other_is_a_violation_via_board_row():
    violation, _board, _other = cbs.improper_mix(
        ['PROJECT_RULES.md', 'pathway_forward.md', 'src/solver.txt'])
    assert violation is True, (
        'pathway_forward.md present alongside other files must still refuse, '
        'even with PROJECT_RULES.md also in the set')


def test_board_files_alone_no_other_is_never_a_violation():
    violation, _b, _o = cbs.improper_mix(['PROJECT_RULES.md', 'pathway_forward.md'])
    assert violation is False


def test_no_board_files_is_never_a_violation():
    violation, _b, _o = cbs.improper_mix(['src/a.txt', 'src/b.txt'])
    assert violation is False


# ------------------------------------------- (1b) per-commit, real git ---

def test_per_commit_rulebook_plus_code_passes(repo):
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    commit(repo, {'PROJECT_RULES.md': 'new rule 99\n', 'src/solver.txt': 'v2\n'},
          'rules: new rule plus its enforcing code')
    r = run_checker(['%s..HEAD' % base], cwd=repo)
    assert r.returncode == 0, r.stdout + r.stderr
    assert 'SUCCESS' in r.stdout


def test_per_commit_board_row_plus_code_fails(repo):
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    commit(repo, {'pathway_forward.md': 'row 200\n', 'src/solver.txt': 'v2\n'},
          'code plus a smuggled board row')
    r = run_checker(['%s..HEAD' % base], cwd=repo)
    assert r.returncode == 1, r.stdout + r.stderr
    assert 'pathway_forward.md' in r.stdout
    assert 'src/solver.txt' in r.stdout


# --------------------------------------------------------- (2) exempt ---

def test_exempt_sha_is_reported_not_silently_passed(repo, monkeypatch):
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    bad_sha = commit(repo, {'pathway_forward.md': 'row 201\n', 'src/solver.txt': 'v3\n'},
                     'a historical mixed commit, to be exempted')
    monkeypatch.chdir(repo)
    monkeypatch.setitem(cbs.EXEMPT_SHAS, bad_sha, 'unit-test exemption')
    code = cbs.check_per_commit('%s..HEAD' % base)
    assert code == 0, 'an exempted sha must not fail the range'
    # Re-run WITHOUT the exemption to prove this sha is only passing because
    # of the exemption, not by accident of the (fixed) improper_mix logic.
    monkeypatch.delitem(cbs.EXEMPT_SHAS, bad_sha)
    code_unexempted = cbs.check_per_commit('%s..HEAD' % base)
    assert code_unexempted == 1, (
        'the same commit must fail once the exemption is removed -- '
        'otherwise the exemption proved nothing')


def test_exempt_sha_prints_exempt_not_success(repo, capsys, monkeypatch):
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    bad_sha = commit(repo, {'pathway_forward.md': 'row 202\n', 'src/solver.txt': 'v4\n'},
                     'another historical mixed commit')
    monkeypatch.chdir(repo)
    monkeypatch.setitem(cbs.EXEMPT_SHAS, bad_sha, 'unit-test exemption 2')
    cbs.check_per_commit('%s..HEAD' % base)
    captured = capsys.readouterr()
    assert 'EXEMPT' in captured.out
    assert captured.out.count('SUCCESS') <= captured.out.count('EXEMPT')


# --------------------------------------------------- (3) --squash-check ---

def test_squash_check_catches_union_mix_per_commit_check_does_not(repo):
    """The PR #69 shape: four commits, each individually clean (rulebook,
    code, board-row, code), whose UNION mixes pathway_forward.md with code.
    --squash-check must fail this range; the plain per-commit check on the
    identical range must pass."""
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    commit(repo, {'PROJECT_RULES.md': 'rule X\n'}, 'rules: X (rule-only)')
    commit(repo, {'src/solver.txt': 'code A\n'}, 'code: A (code-only)')
    commit(repo, {'pathway_forward.md': 'row Y\n'}, 'board: row Y (board-only)')
    commit(repo, {'src/solver.txt': 'code B\n'}, 'code: B (code-only)')

    plain = run_checker(['%s..HEAD' % base], cwd=repo)
    assert plain.returncode == 0, (
        'every individual commit is clean; the per-commit check must pass: %s'
        % (plain.stdout + plain.stderr))

    squashed = run_checker(['%s..HEAD' % base, '--squash-check'], cwd=repo)
    assert squashed.returncode == 1, (
        'the union of this range mixes pathway_forward.md with src/solver.txt '
        '-- squashing it (what this repo\'s merge settings do) would produce '
        'exactly the ff06534 incident, and --squash-check must catch it: %s'
        % (squashed.stdout + squashed.stderr))
    assert 'pathway_forward.md' in squashed.stdout


def test_squash_check_passes_a_genuinely_clean_range(repo):
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    commit(repo, {'PROJECT_RULES.md': 'rule X\n'}, 'rules: X')
    commit(repo, {'src/solver.txt': 'code A\n'}, 'code: A')
    r = run_checker(['%s..HEAD' % base, '--squash-check'], cwd=repo)
    assert r.returncode == 0, r.stdout + r.stderr


def test_squash_check_requires_a_dotdot_range(repo):
    r = run_checker(['HEAD', '--squash-check'], cwd=repo)
    assert r.returncode == 2


def test_no_args_refuses_with_exit_2(repo):
    r = run_checker([], cwd=repo)
    assert r.returncode == 2


def test_squash_check_passes_despite_a_messy_intermediate_commit(repo):
    """GENERALIZED 2026-09-30: an intermediate branch commit that mixes
    pathway_forward.md with code is informational-only under --squash-check
    if the FINAL union (what squash-merge actually lands) is clean -- e.g. a
    later commit in the same PR removes the board-row edit again before
    merge. Individual commits never reach master; only the aggregate does.
    """
    base = git(['rev-parse', 'HEAD'], cwd=repo).stdout.strip()
    commit(repo, {'pathway_forward.md': 'wip row\n', 'src/solver.txt': 'wip\n'},
          'wip: messy intermediate commit (mixed)')
    commit(repo, {'pathway_forward.md': 'seed board\n'},  # revert the board edit
          'fixup: drop the board-row edit before merge')

    per_commit_only = run_checker(['%s..HEAD' % base], cwd=repo)
    assert per_commit_only.returncode == 1, (
        'the intermediate commit really is mixed; the plain per-commit '
        'check must still see it')

    squashed = run_checker(['%s..HEAD' % base, '--squash-check'], cwd=repo)
    assert squashed.returncode == 0, (
        'the AGGREGATE union of this range restores pathway_forward.md to '
        'its original content and touches only src/solver.txt net -- '
        '--squash-check must pass, matching what squash-merging this PR '
        'would actually put on master: %s' % (squashed.stdout + squashed.stderr))
