#! /usr/bin/env python3
"""
Regression guard: a fast-lane PR must still run `build` and
`unit-regression` (owner decision 2026-10-07).

Incident: the fast lane used to skip both, running only fast-lane-checks.
Repo-wide guards that live in unit-regression (test_no_hardcoded_paths.py
and the board/rule readers) scan the whole tree, so a docs-only PR that
added an offending line went green, merged, and turned master's push CI red
after the fact (41ffebf, 1515597; fixed by PR #142).

WHAT THIS PINS:
  1. the `if:` of `build` and of `unit-regression` does not read
     `detect-lane.outputs.lane` -- both run on either lane;
  2. `merge-gate` needs both jobs, and its `fast)` branch requires
     BUILD_RESULT and UNIT_REGRESSION_RESULT to be `success` (a skipped job
     never counts as a pass);
  3. rule 14a: the same check, run on a copy of test.yml with the old
     lane-gated `if:` restored, must FAIL -- proof it can come out both ways.

Cheap (rule 9): two parses of test.yml with the shared stdlib parser.
"""
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')

sys.path.insert(0, HERE)
import test_ci_board_separation_step as board_step  # noqa: E402 -- reuse its YAML parser

BOTH_LANE_JOBS = ('build', 'unit-regression')
RESULT_VARS = ('BUILD_RESULT', 'UNIT_REGRESSION_RESULT')


def check(path):
    """Return a list of problems with `path`; empty means the fast lane
    still runs build and unit-regression."""
    problems = []
    jobs = board_step.parse_workflow(path)['jobs']
    for name in BOTH_LANE_JOBS:
        job = jobs.get(name)
        if job is None:
            problems.append('job %r is missing' % name)
            continue
        job_if = str(job.get('if', ''))
        if 'detect-lane.outputs.lane' in job_if:
            problems.append('job %r is lane-gated again (if: %r) -- a fast-lane '
                            'PR would skip it' % (name, job_if))
    gate = jobs.get('merge-gate')
    if gate is None:
        problems.append('job merge-gate is missing')
        return problems
    needs = str(gate.get('needs', ''))
    for name in BOTH_LANE_JOBS:
        if not re.search(r'(^|[\s\[,])%s([\s\],]|$)' % re.escape(name), needs):
            problems.append('merge-gate does not need %r (needs: %r)' % (name, needs))
    body = '\n'.join(str(s.get('run', '')) for s in gate['steps'] if isinstance(s, dict))
    m = re.search(r'^\s*fast\)(.*?);;', body, re.S | re.M)
    if not m:
        problems.append('merge-gate has no `fast)` case branch')
        return problems
    for var in RESULT_VARS:
        if not re.search(r'"\$%s"\s*=\s*"success"' % var, m.group(1)):
            problems.append('merge-gate fast) branch does not require %s = success' % var)
    return problems


def main():
    print('Regression guard: fast-lane PRs still run build + unit-regression '
          '(owner decision 2026-10-07)')
    problems = check(WORKFLOW)
    jobs = board_step.parse_workflow(WORKFLOW)['jobs']
    print('  parsed %d job(s) from %s' % (len(jobs), WORKFLOW))
    for name in BOTH_LANE_JOBS:
        print('  %s if: %r' % (name, jobs.get(name, {}).get('if')))

    # rule 14a: the old lane-gated condition must make the check go red.
    text = open(WORKFLOW).read()
    old_if = ("if: ${{ !cancelled() && !failure() && (github.event_name != "
              "'pull_request' || needs.detect-lane.outputs.lane == 'full') }}")
    mutated, n = re.subn(r'(\n  unit-regression:\n    needs: [^\n]*\n)    if: [^\n]*\n',
                         lambda mm: mm.group(1) + '    ' + old_if + '\n', text)
    if n != 1:
        problems.append('could not build the rule-14a mutant (unit-regression '
                        'job header not found in the expected shape)')
    else:
        with tempfile.NamedTemporaryFile('w', suffix='.yml', delete=False) as f:
            f.write(mutated)
        try:
            red = check(f.name)
        finally:
            os.unlink(f.name)
        print('  mutant (unit-regression lane-gated again): %d problem(s)' % len(red))
        if not red:
            problems.append('the check passed a test.yml whose unit-regression '
                            'is lane-gated -- it cannot go red')

    if problems:
        print('\nFAIL: %d problem(s):' % len(problems))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nPASS: build and unit-regression run on both lanes, merge-gate '
          'requires both on the fast lane, and the lane-gated mutant goes red')
    return 0


if __name__ == '__main__':
    sys.exit(main())
