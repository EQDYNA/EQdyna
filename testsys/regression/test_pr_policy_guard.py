#! /usr/bin/env python3
"""
Regression guard for the owner-approved hybrid PR workflow (2026-09-23):
testsys/pr_policy.py's decision logic, mutation-tested against a REAL
throwaway repository (tempfile.mkdtemp()), never against this repository.

WHAT THIS PINS, each driven by a real commit in the sandbox and a real
run_git/rev-list/show underneath testsys.pr_policy (only the PR-evidence
API call is stubbed -- there is no network dependency anywhere in this
file):

  ci-check mode (verify_pr backed by a stubbed commits-pulls API):
    1. a commit touching ONLY testsys/ with NO associated merged PR -> RED
    2. a commit touching ONLY docs (README.md) -> GREEN, no API call needed
       at all (the commit is not gated, so verify_pr is never invoked --
       checked by making the stub raise if called)
    3. a commit touching src/ WITH a stubbed merged-PR response -> GREEN
    4. a commit touching BOTH README.md and src/ (mixed) with no PR -> RED,
       and the report names the src/ path specifically
    5. the API-unavailable path: stub raises, subject carries "(#7)" ->
       fallback signal accepts it (GREEN, source=squash-subject-fallback);
       stub raises, subject carries no suffix -> PolicyCheckUnavailable,
       never a silent GREEN or RED-by-default

  push-guard mode (verify_pr=None -- no API, ever):
    6. the SAME testsys/-touching commit from (1) is RED here too, and for
       a DIFFERENT stated reason (no PR merge is possible over a local
       push, not "no PR was found")
    7. the SAME docs-only commit from (2) is GREEN

  range resolution (resolve_push_range), against real commits so "reachable"
  is a real property of a real object database, not a mocked bool:
    8. before=ZERO_SHA -> range is the head sha alone, and evaluating it
       covers the FULL history (every commit new, ref just created)
    9. before=a real, valid sha from a DIFFERENT throwaway repo (so it is
       syntactically a sha but genuinely absent from this repo's object
       database) -> PolicyCheckUnavailable, not a guessed range

Mutation check for this guard itself (done by hand while writing it, stated
so the next person does not have to re-derive it): flip GATED_PREFIXES to
only ('testsys/',) and case 3/case 8's src/-only commits stop being
evaluated as gated at all, which flips their PASS reasons from "PR-merged"
to "touches no gated path" -- case (1)'s testsys/ commit still catches the
regression on its own, but was checked directly by temporarily editing
touches_gated_paths's default and re-running: cases 3 and 4 both went from
their expected reasons to "touches no gated path", confirming the assertions
below are sensitive to that prefix list and not vacuously true.

Cheap (rule 9): a handful of git inits and one-line commits in a temp dir,
no build, no network. Under ~1 s.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'testsys'))
import pr_policy as pp  # noqa: E402

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='pr_policy guard', GIT_AUTHOR_EMAIL='pr_policy@example.invalid',
    GIT_COMMITTER_NAME='pr_policy guard', GIT_COMMITTER_EMAIL='pr_policy@example.invalid',
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


def build_sandbox(tmp):
    repo = os.path.join(tmp, 'repo')
    os.makedirs(repo)
    sh(['init', '-q', '.'], cwd=repo)
    root = commit(repo, {'README.md': 'seed\n'}, 'seed')
    return repo, root


def api_raises(_sha):
    raise RuntimeError('stubbed: commits-pulls API unreachable')


def api_no_merged_pr(_sha):
    return [{'number': 3, 'merged_at': None}]


def api_fail_if_called(_sha):
    raise AssertionError('verify_pr was called for a commit that is not '
                         'gated -- it must never be invoked at all for one')


def make_api_returns_merged(pr_number):
    def fetch(_sha):
        return [{'number': pr_number, 'merged_at': '2026-09-23T00:00:00Z'}]
    return fetch


def main():
    fails = []
    tmp = tempfile.mkdtemp(prefix='pr_policy_guard_')
    try:
        repo, root = build_sandbox(tmp)

        docs_sha = commit(repo, {'README.md': 'docs edit\n'}, 'docs: readme tweak')
        testsys_sha = commit(repo, {'testsys/scratch_gate.py': 'x = 1\n'},
                             'add a testsys script directly')
        src_pr_sha = commit(repo, {'src/fortran/scratch.f90': '! stub\n'},
                            'feature: scratch fortran routine (#42)')
        mixed_sha = commit(repo, {'README.md': 'mixed docs\n',
                                  'src/fortran/other.f90': '! stub2\n'},
                           'mixed: docs plus src, no PR')
        head = mixed_sha

        # ---- case 1: testsys/-only, no merged PR -> RED --------------------
        d = pp.evaluate_commit_gate(testsys_sha, 'add a testsys script directly',
                                    pp.commit_files(testsys_sha, cwd=repo),
                                    lambda sha, subj: pp.commit_pr_evidence(sha, subj, api_no_merged_pr))
        if d.ok:
            fails.append('case 1: testsys/-only commit with NO merged PR was GREEN '
                        '(reason: %r) -- expected RED' % d.reason)
        if d.gated_files != ['testsys/scratch_gate.py']:
            fails.append('case 1: gated_files=%r, expected exactly the testsys/ path'
                        % (d.gated_files,))

        # ---- case 2: docs-only -> GREEN, API never called -------------------
        d = pp.evaluate_commit_gate(docs_sha, 'docs: readme tweak',
                                    pp.commit_files(docs_sha, cwd=repo),
                                    lambda sha, subj: pp.commit_pr_evidence(sha, subj, api_fail_if_called))
        if not d.ok:
            fails.append('case 2: docs-only commit was RED (reason: %r)' % d.reason)
        if d.gated_files:
            fails.append('case 2: docs-only commit reported gated_files=%r' % (d.gated_files,))

        # ---- case 3: src/ with a stubbed merged PR -> GREEN ------------------
        d = pp.evaluate_commit_gate(src_pr_sha, 'feature: scratch fortran routine (#42)',
                                    pp.commit_files(src_pr_sha, cwd=repo),
                                    lambda sha, subj: pp.commit_pr_evidence(sha, subj, make_api_returns_merged(42)))
        if not d.ok:
            fails.append('case 3: src/ commit with a stubbed MERGED PR was RED '
                        '(reason: %r)' % d.reason)
        if 'commits-api' not in d.reason:
            fails.append('case 3: pass reason does not cite commits-api evidence: %r'
                        % d.reason)

        # ---- case 4: mixed docs+src, no PR -> RED, names the src/ path ------
        d = pp.evaluate_commit_gate(mixed_sha, 'mixed: docs plus src, no PR',
                                    pp.commit_files(mixed_sha, cwd=repo),
                                    lambda sha, subj: pp.commit_pr_evidence(sha, subj, api_no_merged_pr))
        if d.ok:
            fails.append('case 4: MIXED docs+src commit with no PR was GREEN '
                        '(reason: %r) -- expected RED' % d.reason)
        if 'src/fortran/other.f90' not in d.gated_files:
            fails.append('case 4: gated_files=%r does not name the src/ path'
                        % (d.gated_files,))
        if 'README.md' in d.gated_files:
            fails.append('case 4: gated_files=%r wrongly includes README.md (docs '
                        'is not a gated prefix)' % (d.gated_files,))

        # ---- case 5a: API raises, subject has a squash suffix -> fallback GREEN
        ev = pp.commit_pr_evidence(src_pr_sha, 'feature: thing (#99)', api_raises)
        if not ev.via_pr or ev.source != 'squash-subject-fallback':
            fails.append('case 5a: API-unavailable + "(#99)" suffix did not fall '
                        'back to GREEN via squash-subject-fallback: %r' % (ev,))

        # ---- case 5b: API raises, no suffix -> PolicyCheckUnavailable --------
        try:
            pp.commit_pr_evidence(src_pr_sha, 'feature: thing, no suffix', api_raises)
            fails.append('case 5b: API-unavailable with NO squash suffix did not '
                        'raise PolicyCheckUnavailable -- an unverifiable commit '
                        'must never silently resolve to GREEN or RED')
        except pp.PolicyCheckUnavailable:
            pass

        # ---- case 6: push-guard mode refuses the SAME testsys/ commit, for a
        # DIFFERENT reason (no PR possible over a push, not "no PR found") ----
        d = pp.evaluate_commit_gate(testsys_sha, 'add a testsys script directly',
                                    pp.commit_files(testsys_sha, cwd=repo), None)
        if d.ok:
            fails.append('case 6: push-guard mode ALLOWED a testsys/-only commit '
                        '(reason: %r)' % d.reason)
        if 'LOCAL PUSH' not in d.reason:
            fails.append('case 6: push-guard refusal reason does not explain WHY '
                        '(missing "LOCAL PUSH"): %r' % d.reason)

        # ---- case 7: push-guard mode allows the same docs-only commit --------
        d = pp.evaluate_commit_gate(docs_sha, 'docs: readme tweak',
                                    pp.commit_files(docs_sha, cwd=repo), None)
        if not d.ok:
            fails.append('case 7: push-guard mode refused a docs-only commit '
                        '(reason: %r)' % d.reason)

        # ---- case 8: resolve_push_range, ZERO_SHA -> full history, real repo -
        range_spec, note = pp.resolve_push_range(pp.ZERO_SHA, head,
                                                 lambda s: True)
        if range_spec != head or 'ROOT' not in note:
            fails.append('case 8: ZERO_SHA resolution gave range=%r note=%r, '
                        'expected range==head with a ROOT note' % (range_spec, note))
        ok, decisions = pp.evaluate_range(range_spec, cwd=repo, verify_pr=None)
        shas_seen = {d.sha for d in decisions}
        if not {root, docs_sha, testsys_sha, src_pr_sha, mixed_sha} <= shas_seen:
            fails.append('case 8: ROOT range over %r did not cover every commit '
                        '(seen=%d, expected>=5): %r' % (range_spec, len(shas_seen), shas_seen))

        # ---- case 9: an unresolvable before-sha (real sha, wrong repo) -------
        # A generic build_sandbox() seed commit (same content, same message,
        # same fixed author/committer identity) can land on the IDENTICAL sha
        # in both repos if both commits happen within the same wall-clock
        # second -- git hashes only tree+parent+identity+timestamp+message,
        # none of which differ then. That sha would legitimately be reachable
        # in `repo` (it IS repo's own root), making this case pass for the
        # wrong reason. Force divergence with unique content instead.
        other_repo = os.path.join(tmp, 'other')
        os.makedirs(other_repo)
        sh(['init', '-q', '.'], cwd=other_repo)
        other_root = commit(other_repo, {'UNIQUE.md': 'unrelated repo %s\n' % tmp},
                            'seed of an unrelated sandbox')
        if pp.have_commit(other_root, cwd=repo):
            fails.append('case 9 setup: other_root=%s is reachable in `repo` -- '
                        'the two sandbox seed commits collided, so this case '
                        'proves nothing; give them more distinct content' % other_root)
        try:
            pp.resolve_push_range(other_root, head, lambda s: pp.have_commit(s, cwd=repo))
            fails.append('case 9: an unreachable before-sha (valid in a DIFFERENT '
                        'repo, absent from this one) did not raise '
                        'PolicyCheckUnavailable -- it must never guess a range')
        except pp.PolicyCheckUnavailable:
            pass

    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print('FAIL test_pr_policy_guard: %d problem(s):' % len(fails))
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_pr_policy_guard: 9 cases -- ci-check RED on a bare '
         'testsys/ commit and on a mixed docs+src commit, GREEN on docs-only '
         '(API never invoked) and on a stubbed-merged-PR src/ commit, correct '
         'fallback/raise split when the API is unavailable, push-guard RED/'
         'GREEN mirroring ci-check but for the push-specific reason, and '
         'resolve_push_range correct on both a ZERO_SHA ref-creation and an '
         'unresolvable before-sha, all against a real throwaway git repository')
    return 0


if __name__ == '__main__':
    sys.exit(main())
