#! /usr/bin/env python3
"""
Regression guard, board row 160: `run.py release`'s carry-forward must not
reuse a FAILED sweep (rules 2, 14a).

THE DEFECT: testsys/run.py _release_sweep_carry_forward returned skip=True
whenever a summary.json's release_physics_key matched this tree's, without
reading the verdict, so a sweep that FAILED cells was "carried forward" and
the release sweep silently skipped. check_pretag_ci's tag-time gate would
still refuse the tag, but only after the sweep was wrongly skipped.

Both outcomes are pinned, on a throwaway git repo + evidence dir (the repo
under test is never touched):
  1. key matches, summary FAILED (n_success < n_runnable, a non-SUCCESS cell)
     -> skip False, reason names "FAILED sweep: n_success/n_runnable".
  2. key matches, summary all-green -> skip True.
  3. key matches, verdict field malformed (cells not a list; n_success
     missing or a string) -> skip False, with a reason; no default.

Cheap (rule 9): a tmp git repo with one commit, a few small JSON files, no
solver, no MPI. Well under 1 s. Exits non-zero on any failure.

Usage: test_row160_carry_forward_verdict.py [path/to/run.py]  (the optional
path lets the guard be pointed at an unfixed copy to prove it goes red).
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUN_PY = os.path.join(ROOT, 'testsys', 'run.py')
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import content_key  # noqa: E402


def _load_run(path):
    spec = importlib.util.spec_from_file_location('run_under_test', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _git(repo, *args):
    subprocess.run(['git', '-C', repo] + list(args), check=True,
                   capture_output=True, text=True)


def _carry(run_py, summary):
    """(skip, message) from run_py's carry-forward over a tmp repo whose only
    evidence is `summary` with the tree's own release_physics_key stamped in."""
    with tempfile.TemporaryDirectory() as repo:
        _git(repo, 'init', '-q')
        _git(repo, 'config', 'user.email', 't@example.invalid')
        _git(repo, 'config', 'user.name', 'tmp')
        with open(os.path.join(repo, 'f.txt'), 'w') as fh:
            fh.write('x\n')
        _git(repo, 'add', 'f.txt')
        _git(repo, 'commit', '-q', '-m', 'init')
        head = subprocess.run(['git', '-C', repo, 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=True).stdout.strip()
        key = content_key.compute_release_physics(repo, head)
        d = os.path.join(repo, 'docs', 'evidence', 'sweep-abc')
        os.makedirs(d)
        summary = dict(summary, release_physics_key=key, sha=head)
        with open(os.path.join(d, 'summary.json'), 'w') as fh:
            json.dump(summary, fh)
        mod = _load_run(run_py)
        mod.REPO_ROOT = repo
        return mod._release_sweep_carry_forward()


def _cells(verdicts):
    return [{'case': 'test.c%d' % i, 'backend': 'fortran', 'verdict': v}
            for i, v in enumerate(verdicts)]


def check_failed_sweep_with_matching_key_is_not_carried(run_py=RUN_PY):
    skip, msg = _carry(run_py, {'n_runnable': 3, 'n_success': 2,
                                'cells': _cells(['SUCCESS', 'SUCCESS', 'FAIL'])})
    if skip is not False:
        raise AssertionError('a FAILED sweep with a matching key was carried '
                             'forward (skip=%r): %s' % (skip, msg))
    if 'FAILED sweep: 2/3' not in msg:
        raise AssertionError('refusal does not say "FAILED sweep: 2/3": %s' % msg)


def check_green_sweep_with_matching_key_is_carried(run_py=RUN_PY):
    skip, msg = _carry(run_py, {'n_runnable': 2, 'n_success': 2,
                                'cells': _cells(['SUCCESS', 'SUCCESS'])})
    if skip is not True:
        raise AssertionError('an all-green sweep with a matching key was not '
                             'carried forward (skip=%r): %s' % (skip, msg))


def check_malformed_verdict_is_not_carried(run_py=RUN_PY):
    for label, summ in (
            ('cells not a list', {'n_runnable': 1, 'n_success': 1, 'cells': 'ok'}),
            ('n_success missing', {'n_runnable': 1, 'cells': _cells(['SUCCESS'])}),
            ('n_success a string', {'n_runnable': 1, 'n_success': '1',
                                    'cells': _cells(['SUCCESS'])})):
        skip, msg = _carry(run_py, summ)
        if skip is not False:
            raise AssertionError('malformed verdict (%s) was carried forward '
                                 '(skip=%r): %s' % (label, skip, msg))
        if 'FAILED sweep' not in msg:
            raise AssertionError('malformed verdict (%s) gave no reason: %s'
                                 % (label, msg))


def main(run_py=RUN_PY):
    checks = [check_failed_sweep_with_matching_key_is_not_carried,
              check_green_sweep_with_matching_key_is_carried,
              check_malformed_verdict_is_not_carried]
    failures = []
    for c in checks:
        try:
            c(run_py)
        except AssertionError as e:
            failures.append(str(e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_row160_carry_forward_verdict (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_row160_carry_forward_verdict')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else RUN_PY))
