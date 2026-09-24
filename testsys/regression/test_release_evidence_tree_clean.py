#! /usr/bin/env python3
"""
Guard for run_e2e.py's release-evidence `tree_clean` field (2026-09-23 fix).

THE DEFECT THIS GUARDS. `write_release_evidence` used to compute
`tree_clean` from `git status --porcelain` taken at the END of the sweep,
after the sweep itself had already appended a row to
docs/perf_ledger.jsonl, written docs/perf_snapshots/e2e_cells_*.json and
created docs/evidence/sweep-<sha>/ -- none of which is gitignored. So
`tree_clean` read False for EVERY real sweep: proof is
docs/evidence/sweep-679d8ad/summary.json (committed in f9c055e), which reads
tree_clean: False from a run whose `git status --porcelain` was 0 lines when
it STARTED. `testsys/regression/check_pretag_ci.py` requires tree_clean is
True, so a real release sweep could never satisfy the pre-tag guard.

THE FIX. `capture_start_tree_state()` is called once, at the top of
run_e2e.main(), before Gate 0, before the build and before the test/
rotation -- i.e. before this process writes anything. Its (tree_clean,
dirty_lines) is threaded through to `write_release_evidence`, which no
longer recomputes `git status` itself.

WHAT THIS FILE CHECKS, on REAL git state in throwaway sandbox repos (never
this repository, never `core.hooksPath` -- same technique as
test_pretag_sweep_negative.py):
  1. clean start + the sweep's OWN artifacts written afterward -> tree_clean
     True in the written evidence payload (the exact case that was always
     False before this fix).
  2. a real tracked-file edit (testsys/ or src/) present BEFORE the sweep
     starts -> tree_clean False, and the file is NAMED in dirty_at_start.
     Reverting that edit restores tree_clean True and disturbs nothing else.
  3. mutation check (rule 6): calling capture_start_tree_state() again AFTER
     the sweep's own writes (i.e. what the reverted, buggy end-of-run
     capture point would see) flips the SAME sandbox to tree_clean False --
     proof the guard is actually sensitive to WHERE it is called, not just
     to whether a bug is present.

Cheap (rule 9): two `git init` sandboxes, a handful of file writes, no
network, no build -- well under 1 s.
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

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='tree-clean guard', GIT_AUTHOR_EMAIL='guard@example.invalid',
    GIT_COMMITTER_NAME='tree-clean guard', GIT_COMMITTER_EMAIL='guard@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def load_run_e2e():
    """The real module, imported fresh from its path -- testsys/e2e is not a
    package, and importing it fresh (rather than via sys.path + import)
    means test files that also monkeypatch REPO_ROOT cannot interfere with
    each other."""
    path = os.path.join(TESTSYS, 'e2e', 'run_e2e.py')
    spec = importlib.util.spec_from_file_location(
        'run_e2e_under_test_tree_clean', path)
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


def new_repo(tmp, name):
    """A minimal sandbox that mimics the shape this guard cares about: a
    tracked src/ file and a tracked testsys/ file, so 'a modified src/ or
    testsys/ file at start' (the task's own wording) is a real tracked-file
    edit, not a synthetic dict. docs/perf_ledger.jsonl and a docs/perf_
    snapshots/ placeholder are ALSO seeded as tracked (matching the real
    repo, where both are already committed and .gitignore does not exclude
    either) so that a later append/new-file inside docs/ shows up per-path
    in `git status --porcelain` (`M docs/perf_ledger.jsonl`, `?? docs/
    perf_snapshots/e2e_cells_fake.json`) instead of collapsing to a single
    `?? docs/` line for a directory that did not exist at all before --
    an artifact of an under-seeded sandbox, not of the code under test."""
    d = os.path.join(tmp, name)
    os.makedirs(d)
    git(['init', '-q', '.'], cwd=d)
    write(d, 'README.md', 'seed\n')
    write(d, 'src/fortran/placeholder.f90', '! sandbox stand-in\n')
    write(d, 'testsys/matrix.py', '# sandbox stand-in\n')
    write(d, 'docs/perf_ledger.jsonl', '')
    write(d, 'docs/perf_snapshots/.keep', '')
    commit_all(d, 'seed')
    return d


FAKE_RESULTS = [('test.tpv8', 'python-jax', True, 1.2, ['max|diff|=0.0'])]


def check_clean_start_survives_sweep_writes(run_e2e, tmp, fails, log):
    """Scenario 1 + scenario 3 (mutation) together: this is the exact shape
    of a real release sweep."""
    d = new_repo(tmp, 'clean_start')
    run_e2e.REPO_ROOT = d

    tc, dirty = run_e2e.capture_start_tree_state()
    log.append(('1a capture before any sweep write', tc, dirty))
    if not (tc is True and dirty == []):
        fails.append('1a: sandbox was freshly committed but capture_start_'
                     'tree_state() returned tree_clean=%r dirty=%r (want '
                     'True, [])' % (tc, dirty))
        return

    # Exactly what a real sweep writes, none of it gitignored: an appended
    # ledger row, a new perf snapshot, a new evidence directory.
    os.makedirs(os.path.join(d, 'docs', 'perf_snapshots'), exist_ok=True)
    with open(os.path.join(d, 'docs', 'perf_ledger.jsonl'), 'a') as f:
        f.write(json.dumps({'fake': 'row'}) + '\n')
    with open(os.path.join(d, 'docs', 'perf_snapshots', 'e2e_cells_fake.json'), 'w') as f:
        json.dump({'fake': 'snapshot'}, f)

    run_e2e.write_release_evidence(
        FAKE_RESULTS, True, False, '2026-09-23T00:00:00Z',
        '2026-09-23T00:05:00Z', tc, dirty)
    sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    summary_path = os.path.join(d, 'docs', 'evidence', 'sweep-%s' % sha[:7],
                                'summary.json')
    if not os.path.exists(summary_path):
        fails.append('1b: write_release_evidence did not write %s' % summary_path)
        return
    payload = json.load(open(summary_path))
    log.append(('1b written payload tree_clean', payload.get('tree_clean'),
               payload.get('dirty_at_start')))
    if payload.get('tree_clean') is not True:
        fails.append(
            '1b: THE DEFECT IS BACK -- a release sweep that started with a '
            'clean tree wrote tree_clean=%r after writing its own ledger/'
            'snapshot/evidence artifacts (dirty_at_start=%r); this is '
            'exactly the f9c055e/sweep-679d8ad shape'
            % (payload.get('tree_clean'), payload.get('dirty_at_start')))

    # Scenario 3 -- mutation check (rule 6): prove the guard is sensitive to
    # WHERE it is called, by calling it again now, after the sweep's own
    # writes -- i.e. what the reverted (buggy) end-of-run capture point
    # would see, on this SAME sandbox and these SAME artifacts.
    tc_end, dirty_end = run_e2e.capture_start_tree_state()
    log.append(('1c capture AFTER sweep writes (simulated end-of-run bug)',
               tc_end, dirty_end))
    if tc_end is not False:
        fails.append(
            '1c: mutation check failed -- capturing status AFTER the sweep '
            'wrote its own artifacts should read tree_clean=False (proving '
            'the fix actually depends on WHEN it is called), but got %r'
            % tc_end)
    expect_substrings = ('docs/perf_ledger.jsonl', 'docs/perf_snapshots',
                        'docs/evidence')
    for sub in expect_substrings:
        if not any(sub in ln for ln in dirty_end):
            fails.append('1c: end-of-run capture did not name %r among %r'
                         % (sub, dirty_end))


def check_dirty_tracked_file_at_start(run_e2e, tmp, fails, log):
    """Scenario 2: a real edit to a tracked testsys/ file, present BEFORE the
    sweep starts, must be named and must fail the guard; reverting it must
    restore a clean reading and disturb nothing else."""
    d = new_repo(tmp, 'dirty_start')
    run_e2e.REPO_ROOT = d

    tc0, dirty0 = run_e2e.capture_start_tree_state()
    log.append(('2a capture before any edit', tc0, dirty0))
    if not (tc0 is True and dirty0 == []):
        fails.append('2a: freshly committed sandbox read tree_clean=%r (want '
                     'True)' % tc0)
        return

    target = os.path.join(d, 'testsys', 'matrix.py')
    with open(target, 'a') as f:
        f.write('\n# guard-test mutation\n')

    tc1, dirty1 = run_e2e.capture_start_tree_state()
    log.append(('2b capture with testsys/matrix.py modified', tc1, dirty1))
    if tc1 is not False:
        fails.append('2b: a modified tracked testsys/ file did not flip '
                     'tree_clean to False (got %r)' % tc1)
    if not any('testsys/matrix.py' in ln for ln in dirty1):
        fails.append('2b: dirty_at_start does not NAME testsys/matrix.py: %r'
                     % dirty1)

    git(['checkout', '--', 'testsys/matrix.py'], cwd=d)
    tc2, dirty2 = run_e2e.capture_start_tree_state()
    log.append(('2c capture after revert', tc2, dirty2))
    if not (tc2 is True and dirty2 == []):
        fails.append('2c: reverting the edit did not restore tree_clean=True, '
                     'dirty=[] (got %r, %r) -- something else moved' % (tc2, dirty2))


def main():
    print('Guard: run_e2e.capture_start_tree_state / write_release_evidence '
          'tree_clean, 2 scenarios (real sandbox git state)')
    run_e2e = load_run_e2e()
    real_repo_root = run_e2e.REPO_ROOT
    fails = []
    log = []
    tmp = tempfile.mkdtemp(prefix='release_evidence_tree_clean_')
    try:
        check_clean_start_survives_sweep_writes(run_e2e, tmp, fails, log)
        check_dirty_tracked_file_at_start(run_e2e, tmp, fails, log)
    finally:
        run_e2e.REPO_ROOT = real_repo_root  # never leave this pointed at a deleted sandbox
        shutil.rmtree(tmp, ignore_errors=True)

    for label, a, b in log:
        print('  %-55s -> %r %r' % (label, a, b))

    if fails:
        print('\nFAIL test_release_evidence_tree_clean (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_release_evidence_tree_clean')
    return 0


if __name__ == '__main__':
    sys.exit(main())
