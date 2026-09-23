#! /usr/bin/env python3
"""
Regression guard: `.github/workflows/publish.yml` must check out at
`fetch-depth: 0`, because the image it builds ships this checkout's `.git`
and two in-image guards read real history from it. It must ALSO check out
with `persist-credentials: false` (added for the same incident class, below).

THE CREDENTIAL-IN-IMAGE INCIDENT (2026-09-23). `actions/checkout@v4`'s
DEFAULT is `persist-credentials: true`: with that default, checkout writes
the live `GITHUB_TOKEN` into the workspace's `.git/config` as an
`http.<origin>/.extraheader` entry, and the cleanup for that is a POST step
(`action.yml`: `post: dist/index.js`) that runs AFTER every main step in the
job -- including `docker build`. Because `Dockerfile` does
`COPY . /opt/eqdyna` and `.dockerignore` deliberately keeps `.git` (the exact
premise `git_exclusion_hits` below already checks, for the OTHER incident),
the default would bake that token into a layer of every image this workflow
publishes. The token is short-lived -- it expires when the job ends -- so
this is not a live-token leak, but the image outlives the job, so it is a
credential-in-published-artifact leak. Fixed by adding
`persist-credentials: false` to the checkout step; this guard pins that
value the same way it already pins `fetch-depth: 0`, for the same reason:
nothing under `testsys/` referenced this setting at all before this guard,
so it could be reverted (or the whole `with:` block deleted) and every tier
would still exit 0.

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
  2b. EVERY checkout step in that job carries an explicit
     `persist-credentials: false`. Absent is not neutral here either: the
     action's default is `true`, which is the broken state (see the
     credential-in-image incident above) -- deleting the whole `with:` block
     therefore fails BOTH 2 and 2b, not just one;
  3. the PREMISE still holds -- `Dockerfile` still copies the repo root
     (`COPY . /opt/eqdyna`) and both history-reading guards still exist;
  4. (pathway item 79) `.dockerignore` EXISTS -- a bare `COPY . /opt/eqdyna`
     with no `.dockerignore` copies every session/build artifact on a
     developer's box (`.claude/`, `test/`, `bin/`, ...) into the context --
     AND no pattern in it excludes `.git`. This is the trap item 79 called
     out explicitly: the obvious fix for the artifact problem is a
     `.dockerignore` that excludes `.git` (it looks like just more session
     cruft), which reproduces THIS SAME INCIDENT in a worse form -- no
     history at any fetch-depth, not just a shallow one. Checked by replaying
     the file's patterns, in order, against several representative paths
     under `.git` (the directory itself, and a few paths inside it), with
     dockerignore/.gitignore matching semantics: unanchored bare names and
     globs (`.git`, `.git/`, `.git*`) match at any depth and exclude
     everything below them, `**/.git` matches through any number of parent
     directories, and a `!`-prefixed line un-excludes a path a PRIOR pattern
     matched -- so a broad `*` with no matching `!.git` negation is caught
     too, not just the literal string `.git`. If any of this changes, the
     depth requirement may no longer be the right one, and this guard says so
     instead of silently outliving its reason.

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
import fnmatch
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


# ---------------------------------------------- .dockerignore vs `.git` ---
# A future editor's obvious fix for build-context bloat is a .dockerignore
# that excludes `.git` -- it looks like session cruft, not history. This
# replays a .dockerignore's patterns, IN ORDER, against representative paths
# under `.git`, using dockerignore/.gitignore matching rules: a bare name or
# glob with no leading '/' matches at ANY depth (not just the root); a
# pattern that matches a directory excludes everything below it too (an
# ancestor match); '**' matches zero or more WHOLE path segments; and a
# leading '!' un-excludes whatever a PRIOR pattern excluded. A naive
# substring/equality test on the raw pattern text would miss `.git*` and a
# broad `*` left un-negated for `.git` -- both real ways to exclude `.git`
# without the literal string `.git` ever appearing on a line by itself.
GIT_PROBE_PATHS = (
    '.git',
    '.git/HEAD',
    '.git/refs/heads/master',
    '.git/objects/pack/pack-0000000000000000000000000000000000000000.pack',
)


def _prefix_match(pat_segs, path_segs):
    """True if `pat_segs` matches `path_segs` itself, or matches a PREFIX of
    it (i.e. the pattern denotes an ancestor directory, which excludes
    everything below it) -- `**` in `pat_segs` may consume zero or more
    whole segments."""
    if not pat_segs:
        return True
    head = pat_segs[0]
    if head == '**':
        return (_prefix_match(pat_segs[1:], path_segs)
                or (bool(path_segs) and _prefix_match(pat_segs, path_segs[1:])))
    if not path_segs:
        return False
    return (fnmatch.fnmatchcase(path_segs[0], head)
            and _prefix_match(pat_segs[1:], path_segs[1:]))


def _pattern_matches_path(pattern, path):
    pattern = pattern.rstrip('/')
    anchored = pattern.startswith('/')
    pat_segs = pattern.lstrip('/').split('/') if pattern else []
    path_segs = path.split('/')
    if anchored:
        return _prefix_match(pat_segs, path_segs)
    # Unanchored: the pattern may start matching at ANY depth of the path
    # (gitignore semantics for a bare name), so try every suffix.
    return any(_prefix_match(pat_segs, path_segs[start:])
               for start in range(len(path_segs)))


def git_exclusion_hits(dockerignore_lines):
    """Replay `dockerignore_lines` in order against GIT_PROBE_PATHS. Returns
    a list of (pattern_line, probe_path) for every probe path still excluded
    after ALL lines (including `!` negations) are applied -- empty if none
    is."""
    excluded = dict.fromkeys(GIT_PROBE_PATHS, False)
    deciding_line = dict.fromkeys(GIT_PROBE_PATHS)
    for raw in dockerignore_lines:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        negated = line.startswith('!')
        body = line[1:].strip() if negated else line
        if not body:
            continue
        for probe in GIT_PROBE_PATHS:
            if _pattern_matches_path(body, probe):
                excluded[probe] = not negated
                deciding_line[probe] = line
    return [(deciding_line[p], p) for p in GIT_PROBE_PATHS if excluded[p]]


def main():
    print('Regression guard: %s must check out at fetch-depth 0 with '
          'persist-credentials false (the image ships this .git; two '
          'in-image guards read it, and the default credential persistence '
          'would ship the GITHUB_TOKEN in it)' % WORKFLOW_REL)
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

        persist = (with_block or {}).get('persist-credentials')
        if str(persist).lower() != 'false':
            problems.append(
                "%s: job %r checks out `%s` with persist-credentials=%s%s. "
                "The action's default (true) writes the live GITHUB_TOKEN "
                "into .git/config as an http.<origin>/.extraheader entry, "
                "the cleanup for which is a POST step that runs AFTER "
                "`docker build`; Dockerfile ships this .git into the image "
                "(`%s`, and `.dockerignore` deliberately keeps `.git`), so "
                "the default bakes the token into every published image "
                "layer. Set `with: persist-credentials: false`."
                % (WORKFLOW_REL, job_name, co.get('uses'),
                   persist if persist is not None else 'absent',
                   '' if with_block else " (no `with:` block at all; "
                                         "actions/checkout's default is true)",
                   COPY_LINE))
        else:
            print('  %r: persist-credentials false (explicit)' % co.get('uses'))

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
    if not os.path.isfile(DOCKERIGNORE):
        problems.append(
            '.dockerignore does not exist (pathway item 79). `%s` with no '
            '.dockerignore copies EVERY session/build artifact on a '
            "developer's box (.claude/, test/, bin/, ...) into the image "
            'build context.' % COPY_LINE)
    else:
        raw_lines = open(DOCKERIGNORE, errors='replace').read().splitlines()
        hits = git_exclusion_hits(raw_lines)
        if hits:
            detail = '; '.join(
                'pattern %r excludes probe path %r' % (ln, p) for ln, p in hits)
            problems.append(
                '.dockerignore excludes .git (%s), so `%s` ships no history and '
                'the two in-image guards (%s) cannot pass at ANY fetch-depth. '
                'Depth 0 is then not the fix -- re-derive this guard.'
                % (detail, COPY_LINE, ' and '.join(HISTORY_GUARDS)))
        else:
            print('  .dockerignore: present, no pattern excludes .git '
                  '(%d probe path(s) checked)' % len(GIT_PROBE_PATHS))
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
    print('\nPASS: %s job %r checks out at fetch-depth 0 with '
          'persist-credentials false, and the premise holds (Dockerfile '
          '`%s`, no .dockerignore strip of .git, both history-reading '
          'guards present)' % (WORKFLOW_REL, job_name, COPY_LINE))
    return 0


if __name__ == '__main__':
    sys.exit(main())
