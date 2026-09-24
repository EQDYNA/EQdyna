#! /usr/bin/env python3
"""
Regression guard: CI must still RUN the owner-approved hybrid PR workflow
gate (2026-09-23) -- the `pr-policy-gate` job of test.yml, which calls
testsys/pr_policy.py in `ci-check` mode -- and it must be wired correctly
enough to actually gate something. Mirrors
test_ci_board_separation_step.py's shape and its stated reason for existing
as a SIBLING rather than an extension of test_ci_workflow_coverage.py: one
guard per incident, and this one needs structural YAML parsing plus a
subprocess probe of a different artifact.

THREE EDITS would disarm it and leave every push exiting 0:
  1. delete the job (or the pr_policy.py invocation inside it);
  2. loosen or drop the job's `if:` so it also fires on a pull_request event,
     or drop it so it fires on push to EVERY branch -- either way it stops
     being the master-only gate the owner asked for;
  3. put `fetch-depth: 1` (or omit `fetch-depth`, whose default is 1) on the
     job's checkout -- github.event.before is then unreachable and
     testsys/pr_policy.py's resolve_push_range refuses rather than guesses
     (by design), which would make the JOB fail for the wrong reason on
     every real push, indistinguishable from the gate actually catching a
     violation unless someone reads the log.

WHAT THIS PINS, in order:
  1. exactly one step in the workflow invokes `testsys/pr_policy.py
     ci-check`, with two arguments (real YAML content, parsed structurally
     -- test.yml's own comments name this job in prose several times, so a
     grep would be satisfied by the comments alone after the step was
     deleted);
  2. the step carries no `if:` of its own (a per-step condition can disarm
     it while leaving it visible and the JOB green);
  3. the JOB's own `if:` mentions BOTH `github.event_name` and `github.ref`
     -- a job-level `if` that checked only one of those two would fire on
     every push to every branch (event_name-only) or on a pull_request
     event that happens to target master (ref-only, which pull_request
     does not set to refs/heads/master but is worth being explicit against
     rather than assuming);
  4. the job checks out with `fetch-depth: 0`;
  5. the step's env carries a `GITHUB_TOKEN`, without which
     github_commits_pulls_fetcher raises PolicyCheckUnavailable on every
     gated commit regardless of whether it was actually PR-merged;
  6. testsys/pr_policy.py exists, is executable, and still refuses a bad
     invocation with exit 2 (its documented CLI contract) -- proof the
     step's target is runnable, not just present.

Cheap (rule 9): one file parse (reusing test_ci_board_separation_step.py's
own parser -- imported, not copied, so the two guards can never disagree
about how to read test.yml) plus one sub-second subprocess.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')
POLICY = os.path.join(ROOT, 'testsys', 'pr_policy.py')
POLICY_REL = os.path.join('testsys', 'pr_policy.py')

sys.path.insert(0, HERE)
import test_ci_board_separation_step as board_step  # noqa: E402 -- reuse its YAML parser


def invocation_lines(run_body):
    hits = []
    for raw in run_body.splitlines():
        line = raw.split('#', 1)[0].strip()
        if POLICY_REL in line or 'pr_policy.py' in line:
            hits.append(line)
    return hits


def main():
    print('Regression guard: CI must still run pr_policy.py ci-check on '
         'push to master, at fetch-depth 0, with a token (owner-approved '
         'hybrid PR workflow, 2026-09-23)')
    problems = []
    doc = board_step.parse_workflow(WORKFLOW)

    found = []  # (job_name, job, step, invocation lines)
    for job_name, job in doc['jobs'].items():
        for step in job['steps']:
            if not isinstance(step, dict):
                continue
            body = step.get('run')
            if not isinstance(body, str):
                continue
            hits = invocation_lines(body)
            if hits:
                found.append((job_name, job, step, hits))

    print('  parsed %d job(s) from %s' % (len(doc['jobs']), WORKFLOW))
    if not found:
        print('\nFAIL: no step in %s runs `python3 %s ci-check`.' % (WORKFLOW, POLICY_REL))
        print('  The owner-approved hybrid PR workflow is then enforced by')
        print('  nothing on the CI side: a direct push touching src/ or')
        print('  testsys/ would land on master with a green workflow.')
        return 1
    if len(found) > 1:
        problems.append('%d steps invoke pr_policy.py (%s) -- one '
                        'authoritative invocation, or a reviewer cannot tell '
                        'which push is actually gated'
                        % (len(found), ', '.join(j for j, _, _, _ in found)))

    job_name, job, step, hits = found[0]
    print('  step %r in job %r invokes it: %r'
         % (step.get('name', '<unnamed>'), job_name, hits))

    ci_check_hits = [h for h in hits if re.search(r'pr_policy\.py\s+["\']?ci-check', h)]
    if not ci_check_hits:
        problems.append('the invocation %r does not call ci-check mode' % hits)

    with_two_args = [h for h in hits
                     if re.search(r'pr_policy\.py\s+["\']?ci-check["\']?\s+\S.*\S', h)]
    if not with_two_args:
        problems.append('the invocation %r does not appear to pass two '
                        'arguments (before-sha, after-sha); pr_policy.py '
                        'refuses with exit 2 if it does not get exactly '
                        'three CLI tokens' % hits)

    if 'if' in step:
        problems.append("the step carries its own `if: %s` -- a per-step "
                        "condition can disarm it while leaving it visible "
                        "and the job green" % step['if'])

    job_if = job.get('if')
    if not job_if:
        problems.append('job %r has no `if:` at all -- it would run on '
                        'every push to every branch and every pull_request, '
                        'not just a push landing on master' % job_name)
    else:
        if 'github.event_name' not in job_if:
            problems.append("job %r's `if: %s` does not mention "
                            "github.event_name" % (job_name, job_if))
        if 'github.ref' not in job_if:
            problems.append("job %r's `if: %s` does not mention github.ref"
                            % (job_name, job_if))
        if 'master' not in job_if:
            problems.append("job %r's `if: %s` does not name master"
                            % (job_name, job_if))

    checkouts = [s for s in job['steps']
                if isinstance(s, dict) and str(s.get('uses', '')).startswith('actions/checkout')]
    if not checkouts:
        problems.append('job %r has no actions/checkout step' % job_name)
    for co in checkouts:
        depth = (co.get('with') or {}).get('fetch-depth')
        if str(depth) != '0':
            problems.append(
                "job %r checks out with fetch-depth=%s (absent means the "
                "action's default, 1). pr_policy.py resolves a COMMIT RANGE "
                "from github.event.before; at depth 1 that sha is not in the "
                "clone and resolve_push_range refuses by design -- every "
                "real push would then fail this job for a CI-environment "
                "reason, not a policy reason."
                % (job_name, depth if depth is not None else 'absent'))

    env = step.get('env') or {}
    if not isinstance(env, dict) or 'GITHUB_TOKEN' not in env:
        problems.append(
            "the invoking step's env does not set GITHUB_TOKEN -- "
            "github_commits_pulls_fetcher raises PolicyCheckUnavailable for "
            "every gated commit without one, regardless of whether it was "
            "really PR-merged")

    if not os.path.isfile(POLICY):
        problems.append('%s does not exist -- the step invokes nothing' % POLICY)
    else:
        if not os.access(POLICY, os.X_OK):
            problems.append('%s is not executable (rule 13)' % POLICY)
        r = subprocess.run([sys.executable, POLICY], capture_output=True,
                           text=True, timeout=60)
        if r.returncode != 2:
            problems.append(
                'running pr_policy.py with no arguments exited %d, not the '
                '2 its own CLI contract documents -- it is present but no '
                'longer behaves as the step assumes. stderr: %r'
                % (r.returncode, r.stderr.strip()[:200]))
        else:
            print('  %s: present, executable, refuses a bad invocation with '
                 'exit 2' % POLICY_REL)

    if problems:
        print('\nFAIL: %d problem(s):' % len(problems))
        for p in problems:
            print('  - %s' % p)
        return 1
    print('\nPASS: one step (%r, job %r) runs pr_policy.py ci-check with two '
         'arguments and a token, unconditionally at the step level, gated at '
         'the job level to push-to-master only (if: %r), at fetch-depth 0; '
         'pr_policy.py is executable and refuses a bad invocation with exit 2'
         % (step.get('name', '<unnamed>'), job_name, job_if))
    return 0


if __name__ == '__main__':
    sys.exit(main())
