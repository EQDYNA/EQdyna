#! /usr/bin/env python3
"""
Negative test for the release-path guard's CONTENT-KEY sweep-evidence path
(`check_pretag_ci.evaluate_sweep_candidate`'s `content_key` branch, item 2b,
2026-09-30). See `testsys/content_key.py`'s module docstring for why this
path exists: a squash-merge (rule 25's only merge strategy) or a history
rewrite can leave the swept sha unreachable from, or entirely absent from,
the tagged commit's object graph even though the tagged TREE is physics-
identical to what was swept. `content_key` evidence answers "does every
tracked path outside content_key.ALLOWED_* have the same blob at the tag
commit as it did when the sweep ran" -- recomputed from the TAG commit
alone, never from the swept sha.

WHY A SEPARATE FILE. `test_pretag_sweep_negative.py` is the dedicated
coverage for evidence with NO `content_key` field -- the original
ancestor-based rule, unchanged. This file is the dedicated coverage for the
new field; the two partition the guard's two evidence-acceptance paths
instead of one file re-deriving the other's fixtures.

Five scenarios (mutation-tested; see `mutation_self_check`):
  1. same content_key; swept sha is NOT an ancestor of the tag sha (two
     sibling commits with independently-identical tracked content)  -> ACCEPT
  2. swept sha IS an ancestor; only docs/evidence/ added afterward   -> ACCEPT
  3. one byte changed in scripts/machines.py after the swept sha     -> REFUSE
  4. .github/ edited after the swept sha -- content-key evidence does
     NOT exempt .github/ the way testsys/change_class.py's classifier
     does (that classifier answers "does a PR need review", not "is
     this sweep still valid")                                       -> REFUSE
  5. the swept sha recorded in the evidence is not a real object in
     this sandbox's store at all -- content-key evidence is decided
     from the recorded string alone, never by trying to read the
     swept commit, so this is still a clean ACCEPT, never a crash    -> ACCEPT

Every case drives a REAL `git commit` in a throwaway `tempfile.mkdtemp()`
sandbox repo (same technique as test_pretag_sweep_negative.py) -- never this
repository. `full_runnable_count`/`runnable_cells` are injected fixed
callables for the same reason that file injects them: `testsys.matrix`
describes THIS checkout's cases, not the sandbox's.

Cheap (rule 9): a handful of `git init`s and one JSON write per case, well
under 1 s. No network, no build.
"""
import importlib.util
import json
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
from testsys import matrix  # noqa: E402

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='content key guard', GIT_AUTHOR_EMAIL='ck@example.invalid',
    GIT_COMMITTER_NAME='content key guard', GIT_COMMITTER_EMAIL='ck@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def load_guard():
    """The real module, imported fresh from its path (it is not a package),
    exactly as test_pretag_sweep_negative.py does -- `guard.content_key` is
    then the same `testsys.content_key` module it imported."""
    path = os.path.join(HERE, 'check_pretag_ci.py')
    spec = importlib.util.spec_from_file_location(
        'check_pretag_ci_under_test_content_key', path)
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


def new_repo(tmp, name):
    d = os.path.join(tmp, name)
    os.makedirs(d)
    git(['init', '-q', '.'], cwd=d)
    write(d, 'README.md', 'seed\n')
    commit_all(d, 'seed')
    return d


FULL_CELLS = [
    {'case': 'test.tpv8', 'backend': 'fortran', 'verdict': 'SUCCESS',
     'max_diff': 0.0, 'wall_s': 1.0},
    {'case': 'test.tpv8', 'backend': 'python-numpy', 'verdict': 'SUCCESS',
     'max_diff': 0.0, 'wall_s': 1.0},
]


def base_summary(sha, content_key_val, **overrides):
    d = dict(sha=sha, tree_clean=True, term=matrix.GATE_TERM_S, n_runnable=2,
             n_success=2, cells=[dict(c) for c in FULL_CELLS],
             started_utc='2026-09-23T00:00:00Z',
             finished_utc='2026-09-23T00:01:00Z', content_key=content_key_val)
    d.update(overrides)
    return d


def write_summary(cwd, sha, data):
    write(cwd, 'docs/evidence/sweep-%s/summary.json' % sha[:7], json.dumps(data))


def fixed_count():
    """Matches base_summary's n_runnable=2."""
    return 2


def fixed_cells():
    """Matches FULL_CELLS' (case, backend) set."""
    return {(c['case'], c['backend']) for c in FULL_CELLS}


def check_same_key_non_ancestor_accepted(guard, tmp, fails, log):
    d = new_repo(tmp, 'case1_non_ancestor')
    seed = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    write(d, 'scripts/machines.py', 'shared content\n')
    swept_sha = commit_all(d, 'swept: add scripts/machines.py')
    swept_key = guard.content_key.compute(d, swept_sha)
    # A sibling commit off the SAME parent, with independently-written but
    # byte-identical tracked content -> same content_key, but neither
    # commit is an ancestor of the other.
    git(['checkout', '-q', seed], cwd=d)
    write(d, 'scripts/machines.py', 'shared content\n')
    write_summary(d, swept_sha, base_summary(swept_sha, swept_key))
    tag_sha = commit_all(d, 'tag tree: independently identical scripts/machines.py')
    is_anc = guard._sweep_is_ancestor(d, swept_sha, tag_sha)
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count,
                                            runnable_cells=fixed_cells)
    log.append(('1 same content_key, non-ancestor swept sha', ok, msg))
    if is_anc:
        fails.append('1: test setup bug -- swept_sha turned out to be an '
                     'ancestor of tag_sha, so this does not exercise the '
                     'non-ancestor path at all')
    if not ok:
        fails.append('1: same content_key from a non-ancestor swept sha was '
                     'REFUSED: %r' % msg)


def check_evidence_only_addition_accepted(guard, tmp, fails, log):
    d = new_repo(tmp, 'case2_evidence_only')
    write(d, 'scripts/machines.py', 'v1\n')
    swept_sha = commit_all(d, 'swept: add scripts/machines.py')
    swept_key = guard.content_key.compute(d, swept_sha)
    write_summary(d, swept_sha, base_summary(swept_sha, swept_key))
    tag_sha = commit_all(d, 'evidence: record the sweep (docs/evidence/ only)')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count,
                                            runnable_cells=fixed_cells)
    log.append(('2 only docs/evidence/ added after the swept sha', ok, msg))
    if not ok:
        fails.append('2: a docs/evidence/-only commit after the swept sha '
                     'was REFUSED: %r' % msg)


def check_machines_py_byte_change_refused(guard, tmp, fails, log):
    d = new_repo(tmp, 'case3_machines_byte')
    write(d, 'scripts/machines.py', 'v1\n')
    swept_sha = commit_all(d, 'swept: scripts/machines.py = v1')
    swept_key = guard.content_key.compute(d, swept_sha)
    write_summary(d, swept_sha, base_summary(swept_sha, swept_key))
    commit_all(d, 'evidence: record the sweep')
    write(d, 'scripts/machines.py', 'v2\n')  # one byte different
    tag_sha = commit_all(d, 'code: change scripts/machines.py after the sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count,
                                            runnable_cells=fixed_cells)
    log.append(('3 scripts/machines.py changed after the swept sha', ok, msg))
    if ok:
        fails.append('3: scripts/machines.py changed after the swept sha but '
                     'was ACCEPTED (content_key did not move)')
    if 'content_key' not in msg:
        fails.append('3: refusal does not name content_key: %r' % msg)


def check_github_edit_refused(guard, tmp, fails, log):
    d = new_repo(tmp, 'case4_github')
    write(d, '.github/workflows/test.yml', 'jobs: a\n')
    swept_sha = commit_all(d, 'swept: .github/workflows/test.yml = a')
    swept_key = guard.content_key.compute(d, swept_sha)
    write_summary(d, swept_sha, base_summary(swept_sha, swept_key))
    commit_all(d, 'evidence: record the sweep')
    write(d, '.github/workflows/test.yml', 'jobs: b\n')
    tag_sha = commit_all(d, 'ci: edit .github/workflows/test.yml after the sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count,
                                            runnable_cells=fixed_cells)
    log.append(('4 .github/ edited after the swept sha', ok, msg))
    if ok:
        fails.append('4: .github/workflows/test.yml changed after the swept '
                     'sha but was ACCEPTED -- content_key must not treat '
                     '.github/ as excluded the way change_class.py does')
    if 'content_key' not in msg:
        fails.append('4: refusal does not name content_key: %r' % msg)


def check_swept_sha_absent_from_object_store(guard, tmp, fails, log):
    d = new_repo(tmp, 'case5_absent_swept_sha')
    write(d, 'scripts/machines.py', 'v1\n')
    tag_sha = commit_all(d, 'tag tree: scripts/machines.py = v1')
    real_key = guard.content_key.compute(d, tag_sha)
    fake_sha = 'f' * 40  # a well-formed sha that is never an object here
    data = base_summary(fake_sha, real_key)
    reasons = guard.evaluate_sweep_candidate(
        'sandbox-summary.json', data, tag_sha, d,
        full_runnable_count=fixed_count, runnable_cells=fixed_cells)
    log.append(('5 swept sha absent from the object store (want ACCEPT, no crash)',
               not reasons, reasons))
    if reasons:
        fails.append('5: content_key matched the tag commit but evidence '
                     'recording a swept sha this sandbox never committed was '
                     'still refused: %r' % reasons)


def mutation_self_check(guard, tmp, fails, log):
    """Rule 6: prove the content_key equality check is load-bearing by
    neutering it (a `content_key.compute` stand-in that always agrees with
    whatever the evidence claims) and showing scenario 3's refusal
    disappears, then restoring the real function and showing it returns."""
    d = new_repo(tmp, 'case_mutation')
    write(d, 'scripts/machines.py', 'v1\n')
    swept_sha = commit_all(d, 'swept: scripts/machines.py = v1')
    swept_key = guard.content_key.compute(d, swept_sha)
    write(d, 'scripts/machines.py', 'v2\n')
    tag_sha = commit_all(d, 'code: change scripts/machines.py after the sweep')
    data = base_summary(swept_sha, swept_key)

    real_compute = guard.content_key.compute
    guard.content_key.compute = lambda repo_root, sha: swept_key  # neutered
    try:
        reasons_broken = guard.evaluate_sweep_candidate(
            'sandbox-summary.json', data, tag_sha, d,
            full_runnable_count=fixed_count, runnable_cells=fixed_cells)
    finally:
        guard.content_key.compute = real_compute
    log.append(('mutation: content_key check neutered (want PASS)',
               not reasons_broken, reasons_broken))
    if reasons_broken:
        fails.append('mutation check: neutering content_key.compute to always '
                     'agree with the recorded value still produced reasons '
                     '(%r) -- scenario 3 is not isolating the content_key '
                     'comparison' % reasons_broken)

    reasons_fixed = guard.evaluate_sweep_candidate(
        'sandbox-summary.json', data, tag_sha, d,
        full_runnable_count=fixed_count, runnable_cells=fixed_cells)
    log.append(('mutation: content_key check restored (want REFUSE)',
               bool(reasons_fixed), reasons_fixed))
    if not reasons_fixed:
        fails.append('mutation check: restoring the real content_key.compute '
                     'against a scripts/machines.py change produced no '
                     'reasons -- scenario 3\'s assertion is vacuous')


def main():
    print('Negative test: check_pretag_ci content_key sweep evidence, 5 scenarios')
    guard = load_guard()
    fails = []
    log = []
    tmp = tempfile.mkdtemp(prefix='content_key_neg_')
    try:
        check_same_key_non_ancestor_accepted(guard, tmp, fails, log)
        check_evidence_only_addition_accepted(guard, tmp, fails, log)
        check_machines_py_byte_change_refused(guard, tmp, fails, log)
        check_github_edit_refused(guard, tmp, fails, log)
        check_swept_sha_absent_from_object_store(guard, tmp, fails, log)
        mutation_self_check(guard, tmp, fails, log)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for label, ok, msg in log:
        print('  %-60s ok=%s' % (label, ok))
    if os.environ.get('CONTENT_KEY_NEG_VERBOSE'):
        for label, ok, msg in log:
            print('--- %s ---\n%s' % (label, msg))

    if fails:
        print('\nFAIL test_content_key_sweep_evidence (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_content_key_sweep_evidence: 2 acceptances (same '
          'content_key from a non-ancestor swept sha; a docs/evidence/-only '
          'commit after the swept sha), 2 refusals (scripts/machines.py '
          'changed; .github/ edited), and one clean accept for a swept sha '
          'absent from the object store; mutation check confirmed the '
          'content_key comparison is load-bearing')
    return 0


if __name__ == '__main__':
    sys.exit(main())
