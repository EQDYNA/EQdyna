#! /usr/bin/env python3
"""
Regression guard: `.github/workflows/publish.yml` must check out at
`fetch-depth: 0`, because the image it builds ships this checkout's `.git`
and two in-image guards read real history from it.

THE INCIDENT. The root `Dockerfile` does `COPY . /opt/eqdyna`, so whatever
`.git` the publish job checked out is what lands in the image. The publish
job then GATES on `python3 testsys/run.py unit regression` run INSIDE that
image, and two of those guards need real history:

  * `test_history_table.py`      -- the 8 bulk-imported commits of 2018-11-26
  * `test_pretag_ci_negative.py` -- 5 referenced SHAs

At `actions/checkout`'s DEFAULT depth of 1 neither can see its evidence and
both FAIL -- correctly; a guard that cannot see its evidence must fail, not
skip. That took the publish workflow red on 894cdc1, and because the gate step
runs BEFORE the push step, `Push image` was skipped and **the v5.16.0 image
was never published** (run 35819232914's step list). Commit a99902f added
`fetch-depth: 0`; nothing under `testsys/` referenced publish.yml at all, so
that line could be reverted to 1, or the whole `with:` block deleted, and
every tier would still exit 0. This file closes that.

WHAT THIS PINS, in order:
  1. exactly one job in publish.yml builds the image (`docker build`), and it
     has at least one `actions/checkout` step -- an unchecked-out tree has no
     history to argue about;
  2. EVERY checkout step in that job carries an explicit `fetch-depth: 0`.
     Absent is not neutral: the action's default is 1, which is the broken
     state. So both mutations -- `fetch-depth: 1` and deleting the `with:`
     block -- fail here, and they fail with different values printed;
  3. the PREMISE still holds -- `Dockerfile` still copies the repo root
     (`COPY . /opt/eqdyna`), no `.dockerignore` exists to strip `.git` from
     that copy, and both history-reading guards still exist. If any of those
     changes, the depth requirement may no longer be the right one, and this
     guard says so instead of silently outliving its reason.

Structural parsing, not grep: publish.yml's own comment block explains
`fetch-depth 0` in prose, so a grep for `fetch-depth: 0` would still be
satisfied by the comment after the real setting was reverted. The parser is
the one already written for `test_ci_board_separation_step.py` (rule 10: one
guard per incident, but one parser -- a second copy is a second thing to keep
honest). That guard pins test.yml's `unit-regression` job; this one pins
publish.yml's build job. Different workflow, different incident, and every
FAIL line below names `publish.yml` by path so nobody debugs test.yml.

Cheap (rule 9): two file parses, no subprocess. Exits non-zero on any failure
(rule 2).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from test_ci_board_separation_step import parse_workflow  # noqa: E402

WORKFLOW_REL = os.path.join('.github', 'workflows', 'publish.yml')
WORKFLOW = os.path.join(ROOT, WORKFLOW_REL)
DOCKERFILE = os.path.join(ROOT, 'Dockerfile')
DOCKERIGNORE = os.path.join(ROOT, '.dockerignore')
COPY_LINE = 'COPY . /opt/eqdyna'
HISTORY_GUARDS = ('test_history_table.py', 'test_pretag_ci_negative.py')


def build_jobs(doc):
    """Jobs with a step whose `run:` body invokes `docker build` -- the job
    whose checkout becomes the image's `.git`."""
    hits = []
    for job_name, job in doc['jobs'].items():
        for step in job['steps']:
            if not isinstance(step, dict):
                continue
            body = step.get('run')
            if not isinstance(body, str):
                continue
            for raw in body.splitlines():
                if 'docker build' in raw.split('#', 1)[0]:
                    hits.append(job_name)
                    break
            else:
                continue
            break
    return hits


def checkout_steps(job):
    return [s for s in job['steps']
            if isinstance(s, dict)
            and str(s.get('uses', '')).startswith('actions/checkout')]


def main():
    print('Regression guard: %s must check out at fetch-depth 0 (the image '
          'ships this .git; two in-image guards read it)' % WORKFLOW_REL)
    problems = []

    if not os.path.isfile(WORKFLOW):
        print('\nFAIL: %s does not exist. The image publish path is gone or '
              'renamed; re-derive this guard in the same commit.' % WORKFLOW_REL)
        return 1

    doc = parse_workflow(WORKFLOW)
    print('  parsed %d job(s) from %s' % (len(doc['jobs']), WORKFLOW_REL))

    jobs = build_jobs(doc)
    if not jobs:
        print('\nFAIL: no job in %s runs `docker build`. The image is built '
              'somewhere this guard cannot see, so nothing here pins the '
              'depth of the checkout that becomes its .git.' % WORKFLOW_REL)
        return 1
    if len(jobs) > 1:
        problems.append('%d jobs in %s run `docker build` (%s) -- one '
                        'authoritative build job, or a reviewer cannot tell '
                        'which checkout lands in the published image'
                        % (len(jobs), WORKFLOW_REL, ', '.join(jobs)))

    job_name = jobs[0]
    job = doc['jobs'][job_name]
    cos = checkout_steps(job)
    print('  build job %r has %d actions/checkout step(s)' % (job_name, len(cos)))
    if not cos:
        problems.append('job %r of %s has NO actions/checkout step, so the '
                        'tree (and the .git) that `COPY . /opt/eqdyna` ships '
                        'comes from nowhere this guard can pin'
                        % (job_name, WORKFLOW_REL))
    for co in cos:
        with_block = co.get('with')
        depth = (with_block or {}).get('fetch-depth')
        if str(depth) != '0':
            problems.append(
                "%s: job %r checks out `%s` with fetch-depth=%s%s. The image "
                "ships THIS .git (Dockerfile: `%s`), and the in-image gate "
                "runs %s, which read real history -- at depth 1 both FAIL, the "
                "gate step stops the job, and `Push image` is skipped, so the "
                "tag's image is never published (this is incident 894cdc1 / "
                "v5.16.0). Set `with: fetch-depth: 0`."
                % (WORKFLOW_REL, job_name, co.get('uses'),
                   depth if depth is not None else 'absent',
                   '' if with_block else " (no `with:` block at all; "
                                         "actions/checkout's default is 1)",
                   COPY_LINE, ' and '.join(HISTORY_GUARDS)))
        else:
            print('  %r: fetch-depth 0 (explicit)' % co.get('uses'))

    # --- the premise: does the depth still matter for the reason claimed? ---
    if not os.path.isfile(DOCKERFILE):
        problems.append('Dockerfile is gone; %s cannot ship a .git and this '
                        'guard needs re-deriving' % WORKFLOW_REL)
    else:
        body = open(DOCKERFILE, errors='replace').read().splitlines()
        if not any(l.split('#', 1)[0].strip() == COPY_LINE for l in body):
            problems.append(
                'Dockerfile no longer contains `%s`, so the published image '
                'may not ship this checkout\'s .git at all. The fetch-depth 0 '
                'requirement asserted above rests on that line -- re-derive '
                'this guard in the same commit that changed it.' % COPY_LINE)
    if os.path.isfile(DOCKERIGNORE):
        ignored = [l.strip() for l in open(DOCKERIGNORE, errors='replace')
                   if l.strip() and not l.strip().startswith('#')]
        if any(p.rstrip('/') in ('.git', '**/.git') for p in ignored):
            problems.append(
                '.dockerignore excludes .git, so `%s` ships no history and the '
                'two in-image guards cannot pass at ANY fetch-depth. Depth 0 '
                'is then not the fix -- re-derive this guard.' % COPY_LINE)
    for guard in HISTORY_GUARDS:
        if not os.path.isfile(os.path.join(HERE, guard)):
            problems.append(
                'testsys/regression/%s is gone -- it is one of the two guards '
                'whose history reads make depth 0 necessary. If both are gone, '
                'this requirement should be re-derived, not left standing.'
                % guard)

    if problems:
        print('\nFAIL: %d problem(s) in %s:' % (len(problems), WORKFLOW_REL))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nPASS: %s job %r checks out at fetch-depth 0, and the premise '
          'holds (Dockerfile `%s`, no .dockerignore strip of .git, both '
          'history-reading guards present)' % (WORKFLOW_REL, job_name, COPY_LINE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
