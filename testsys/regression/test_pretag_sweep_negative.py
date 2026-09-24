#! /usr/bin/env python3
"""
Negative test for the release-path guard's SWEEP-EVIDENCE half
(`check_pretag_ci.evaluate_sweep_evidence`, added 2026-09-23 -- owner-approved
design: CI stops running the e2e sweep, so a committed local RELEASE sweep
(every supported cell -- everyday cells + matrix.RELEASE_ONLY -- all at the
ONE matrix.GATE_TERM_S) becomes the only mechanical check of the physics for
a release, required in addition to green CI).

WHY A SEPARATE FILE. `test_pretag_ci_negative.py` drives the guard's CI
dimension across its six documented outcomes and, since this new gate
landed, deliberately STUBS `evaluate_sweep_evidence` to `(True, ...)` in every
one of those cases so its exit-0 assertions stay about CI alone (see that
file's "One more thing IS stubbed" note). This file is the dedicated coverage
for the function it stubs -- the two files partition the guard's two gates
instead of one file re-deriving the other's fixtures.

Every case drives a REAL `git commit` in a throwaway `tempfile.mkdtemp()`
sandbox repo (same technique as test_precommit_board_separation_guard.py) --
never this repository, never `core.hooksPath`. `full_runnable_count` is
injected to a fixed callable rather than importing the real `testsys.matrix`:
that module describes THIS checkout's cases, not the sandbox's, and asking it
about a fake repo would test nothing.

Six scenarios, five refusals and one acceptance -- rule 15d's own two-commit
release order is what makes the sixth (docs/evidence/-only commit after the
swept SHA) a real, expected shape rather than an edge case:

  1. no summary.json committed at all                    -> refuse
  2. term != matrix.GATE_TERM_S (a stale/wrong value)     -> refuse
  3. one cell's verdict != SUCCESS                        -> refuse
  4. n_runnable short of the (injected) matrix count      -> refuse
  5. swept SHA is an ancestor; a LATER commit touches src/ -> refuse, names it
  6. swept SHA is an ancestor; the only later commit touches
     docs/evidence/ only (an evidence-recording commit)   -> ACCEPT

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
# GATE_TERM_S only -- a scalar constant of THIS checkout, not a description of
# the sandbox's (fake) cases, so importing it here does not repeat the
# full_runnable_count injection rationale above (docstring, "WHY A SEPARATE
# FILE" / `fixed_count`).
from testsys import matrix  # noqa: E402

SANDBOX_ENV = dict(
    os.environ,
    GIT_AUTHOR_NAME='sweep guard', GIT_AUTHOR_EMAIL='sweep@example.invalid',
    GIT_COMMITTER_NAME='sweep guard', GIT_COMMITTER_EMAIL='sweep@example.invalid',
    GIT_CONFIG_NOSYSTEM='1', HOME=os.devnull + '-nonexistent',
)


def load_guard():
    """The real module, imported fresh from its path -- it is not a package."""
    path = os.path.join(HERE, 'check_pretag_ci.py')
    spec = importlib.util.spec_from_file_location(
        'check_pretag_ci_under_test_sweep', path)
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


def base_summary(sha, **overrides):
    d = dict(sha=sha, tree_clean=True, term=matrix.GATE_TERM_S, n_runnable=2,
             n_success=2, cells=[dict(c) for c in FULL_CELLS],
             started_utc='2026-09-23T00:00:00Z',
             finished_utc='2026-09-23T00:01:00Z')
    d.update(overrides)
    return d


def write_summary(cwd, sha, data):
    write(cwd, 'docs/evidence/sweep-%s/summary.json' % sha[:7], json.dumps(data))


def fixed_count():
    """Matches base_summary's n_runnable=2 -- see module docstring."""
    return 2


def check_missing_summary(guard, tmp, fails, log):
    d = new_repo(tmp, 'case1')
    tag_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('1 missing summary', ok, msg))
    if ok:
        fails.append('1: no summary.json exists anywhere, but the guard ACCEPTED')
    if 'summary.json' not in msg:
        fails.append('1: refusal does not name summary.json: %r' % msg)


def check_wrong_term(guard, tmp, fails, log):
    """A summary carrying a term other than matrix.GATE_TERM_S -- e.g. a
    stale 20 s/6 s per-case value from the retired two-term design, or any
    other drift -- must be refused: there is only ONE term now, and a sweep
    that ran at a different one says nothing about the gate this repo
    actually runs."""
    d = new_repo(tmp, 'case2')
    tag_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    wrong_term = matrix.GATE_TERM_S + 15.0
    write_summary(d, tag_sha, base_summary(tag_sha, term=wrong_term))
    tag_sha = commit_all(d, 'evidence: wrong-term sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('2 term!=GATE_TERM_S', ok, msg))
    if ok:
        fails.append('2: term=%r (!= matrix.GATE_TERM_S=%r) was ACCEPTED'
                     % (wrong_term, matrix.GATE_TERM_S))
    if 'term' not in msg or str(wrong_term) not in msg:
        fails.append('2: refusal does not name the bad term: %r' % msg)


def check_one_fail_cell(guard, tmp, fails, log):
    d = new_repo(tmp, 'case3')
    tag_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    cells = [dict(FULL_CELLS[0]), dict(FULL_CELLS[1], verdict='FAIL')]
    write_summary(d, tag_sha, base_summary(tag_sha, n_success=1, cells=cells))
    tag_sha = commit_all(d, 'evidence: one FAIL cell')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('3 one FAIL cell', ok, msg))
    if ok:
        fails.append('3: a cell with verdict=="FAIL" was ACCEPTED')
    if 'FAIL' not in msg:
        fails.append('3: refusal does not name the FAIL verdict: %r' % msg)


def check_n_runnable_short(guard, tmp, fails, log):
    d = new_repo(tmp, 'case4')
    tag_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    cells = [dict(FULL_CELLS[0])]
    write_summary(d, tag_sha, base_summary(tag_sha, n_runnable=1, n_success=1,
                                           cells=cells))
    tag_sha = commit_all(d, 'evidence: short sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('4 n_runnable short', ok, msg))
    if ok:
        fails.append('4: n_runnable=1 against a declared table of 2 was ACCEPTED')
    if 'n_runnable' not in msg:
        fails.append('4: refusal does not name n_runnable: %r' % msg)


def check_ancestor_src_change_refused(guard, tmp, fails, log):
    d = new_repo(tmp, 'case5')
    swept_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    write_summary(d, swept_sha, base_summary(swept_sha))
    swept_sha = commit_all(d, 'evidence: full sweep at this sha')
    write(d, 'src/solver.txt', 'a real code change after the sweep\n')
    tag_sha = commit_all(d, 'code: change after the sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('5 src/ change after swept ancestor', ok, msg))
    if ok:
        fails.append(
            '5: a src/ change landing AFTER the swept sha and BEFORE the tag '
            'was ACCEPTED -- the sweep no longer describes the tagged tree')
    if 'src/solver.txt' not in msg:
        fails.append('5: refusal does not NAME the offending file: %r' % msg)


def check_evidence_only_ancestor_accepted(guard, tmp, fails, log):
    d = new_repo(tmp, 'case6')
    # The tree actually swept -- recorded INSIDE the JSON below -- is the
    # seed commit, taken BEFORE the evidence commit that reports on it.
    swept_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    write_summary(d, swept_sha, base_summary(swept_sha))
    commit_all(d, 'evidence: full sweep at this sha')
    # rule 15d shape: a second evidence/ledger-only commit lands between the
    # swept sha and the tag (e.g. a perf ledger append), touching no physics.
    write(d, 'docs/perf_ledger.jsonl', '{"note": "post-sweep ledger append"}\n')
    # ...and the sweep's own per-rank profile record (2026-09-24, v5.17.0).
    write(d, 'docs/run_profiles.jsonl', '{"note": "post-sweep profile rows"}\n')
    tag_sha = commit_all(d, 'perf: ledger + run-profile append after the sweep')
    ok, msg = guard.evaluate_sweep_evidence(tag_sha, repo_root=d,
                                            full_runnable_count=fixed_count)
    log.append(('6 evidence-only change after swept ancestor', ok, msg))
    if not ok:
        fails.append(
            '6: two evidence/ledger-only commits after the swept sha (exactly '
            'rule 15d\'s shape) were REFUSED: %r' % msg)
    if swept_sha not in msg:
        fails.append('6: acceptance does not name the swept sha it relied on: %r'
                     % msg)


def mutation_self_check(guard, tmp, fails, log):
    """Rule 6's own demand applied to this file: prove case 4 actually fails
    when the n_runnable check it exercises is broken, not just when nothing
    is wired up. Calls evaluate_sweep_candidate directly with a
    full_runnable_count that always agrees with whatever n_runnable claims --
    the same defect shape as neutering the real check, without editing
    check_pretag_ci.py from inside a test run."""
    d = new_repo(tmp, 'case4_mutation')
    tag_sha = git(['rev-parse', 'HEAD'], cwd=d).stdout.strip()
    data = base_summary(tag_sha, n_runnable=1, n_success=1,
                        cells=[dict(FULL_CELLS[0])])
    reasons_broken = guard.evaluate_sweep_candidate(
        'sandbox-summary.json', data, tag_sha, d,
        full_runnable_count=lambda: data['n_runnable'])  # mutated: always agrees
    log.append(('mutation: n_runnable check neutered (want PASS)',
                not reasons_broken, reasons_broken))
    if reasons_broken:
        fails.append(
            'mutation check: neutering full_runnable_count to always agree '
            'with n_runnable still produced reasons (%r) -- case 4 is not '
            'isolating the n_runnable comparison the way it claims to'
            % reasons_broken)
    reasons_fixed = guard.evaluate_sweep_candidate(
        'sandbox-summary.json', data, tag_sha, d, full_runnable_count=fixed_count)
    log.append(('mutation: n_runnable check restored (want REFUSE)',
                bool(reasons_fixed), reasons_fixed))
    if not reasons_fixed:
        fails.append('mutation check: restoring fixed_count (2) against '
                     'n_runnable=1 produced no reasons -- the real case 4 '
                     'assertion is vacuous')


def main():
    print('Negative test: check_pretag_ci.evaluate_sweep_evidence, 6 scenarios')
    guard = load_guard()
    fails = []
    log = []
    tmp = tempfile.mkdtemp(prefix='sweep_evidence_neg_')
    try:
        check_missing_summary(guard, tmp, fails, log)
        check_wrong_term(guard, tmp, fails, log)
        check_one_fail_cell(guard, tmp, fails, log)
        check_n_runnable_short(guard, tmp, fails, log)
        check_ancestor_src_change_refused(guard, tmp, fails, log)
        check_evidence_only_ancestor_accepted(guard, tmp, fails, log)
        mutation_self_check(guard, tmp, fails, log)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    for label, ok, msg in log:
        print('  %-45s ok=%s' % (label, ok))
    if os.environ.get('SWEEP_NEG_VERBOSE'):
        for label, ok, msg in log:
            print('--- %s ---\n%s' % (label, msg))

    if fails:
        print('\nFAIL test_pretag_sweep_negative (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_pretag_sweep_negative: 5 refusals (missing summary, '
          'wrong term, one FAIL cell, n_runnable short, a src/ change after '
          'the swept ancestor) and 1 acceptance (an evidence/ledger-only '
          'commit after the swept ancestor, rule 15d\'s own shape); mutation '
          'check confirmed the n_runnable comparison is load-bearing')
    return 0


if __name__ == '__main__':
    sys.exit(main())
