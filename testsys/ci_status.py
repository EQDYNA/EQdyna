#! /usr/bin/env python3
"""
CI-run evidence for a commit SHA -- shared by the pre-tag gate
(testsys/regression/check_pretag_ci.py) and the post-hoc release guard
(testsys/regression/test_release_complete.py's check_ci_green_for_tagged_sha).

WHY THIS EXISTS (2026-09-21 incident): tag `v5.13.1` was pushed at `dfee14d`
while the only CI run for that exact SHA was the run the tag push itself
triggered -- at `git tag` time there was no completed run for `dfee14d` at
all. It happened to conclude green. Rule 15 step 6/7 says wait for CI to be
green on the pushed COMMIT, then tag; nothing mechanical enforced that
ordering, and nothing would have noticed had the tag-triggered run gone red,
because a pushed tag cannot be re-pointed (rule 8).

The commit that exposed this was `dfee14d`, which touches ONLY
`pathway_forward.md` -- one of `.github/workflows/test.yml`'s own
`paths-ignore` entries -- so it could never have had a pre-tag run of its
own; waiting for one would wait forever. That is a distinct situation from
"CI ran and failed" and from "CI is still running": rule 2 says those three
must not share an exit code, so this module classifies them separately
rather than folding the paths-ignore case into a generic "no evidence yet".

Everything network-touching here can fail to mean either PASS or FAIL --
`gh` missing, unauthenticated, or the network down. Per the precedent in
test_release_complete.py's check_network_side/check_tag_is_annotated,
"I could not check this" gets its own outcome (UNVERIFIED), never PASS and
never FAIL by default.

`gh run list --commit <sha> --json ...` was tested against this repo's real
history on gh 2.91.0 and returned `[]` for SHAs (`dfee14d`, `813b952`) that
demonstrably have runs (confirmed via `gh run view <id>`). So this module
never uses `--commit`: it lists runs per workflow name and filters on
`headSha` itself.
"""
import fnmatch
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKFLOW_PATH = os.path.join(ROOT, '.github', 'workflows', 'test.yml')

# The five outcomes a caller can get. PASS/UNVERIFIED are shared with the
# regression tier's PASS/UNVERIFIED convention; the other three are the
# ones rule 2 requires to stay apart: CI ran and failed, CI has not
# finished (or not started) yet, and "this SHA cannot trigger CI by design".
PASS = 'PASS'
FAIL = 'FAIL'
PENDING = 'PENDING'
PATHS_IGNORED = 'PATHS_IGNORED'
UNVERIFIED = 'UNVERIFIED'


class GhUnavailable(Exception):
    """`gh` is missing, unauthenticated, or the network call failed outright.
    Never treated as PASS or FAIL -- see module docstring."""


class Result:
    def __init__(self, status, message, evidence_sha=None, hops=None):
        self.status = status
        self.message = message
        self.evidence_sha = evidence_sha
        self.hops = hops or []

    def __repr__(self):
        return 'Result(%s, evidence_sha=%s)' % (self.status, self.evidence_sha)


def _git(*args):
    r = subprocess.run(('git', '-C', ROOT) + args, capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr.strip()


def resolve_sha(ref):
    """Full 40-char SHA for `ref`, or raise if it does not resolve."""
    rc, out, err = _git('rev-parse', '--verify', ref + '^{commit}')
    if rc != 0:
        raise ValueError('%r does not resolve to a commit: %s' % (ref, err))
    return out.strip()


def commit_parent(sha):
    """The single parent of `sha`, or None if it is a root commit.

    Raises on a merge commit (more than one parent): "diff against the
    parent" is ambiguous there, and this module would rather refuse than
    silently pick one and call a paths-ignore-only verdict on the wrong
    comparison.
    """
    rc, out, err = _git('rev-list', '--parents', '-n', '1', sha)
    if rc != 0:
        raise ValueError('cannot read parents of %s: %s' % (sha, err))
    parts = out.split()
    parents = parts[1:]
    if len(parents) > 1:
        raise ValueError(
            '%s is a merge commit (%d parents) -- paths-ignore-only diffing '
            'against a single parent is ambiguous for a merge; resolve by '
            'hand' % (sha, len(parents)))
    return parents[0] if parents else None


def commit_changed_paths(sha):
    """Paths this commit changes relative to its single parent."""
    parent = commit_parent(sha)
    if parent is None:
        raise ValueError('%s is a root commit -- no parent to diff against' % sha)
    rc, out, err = _git('diff-tree', '--no-commit-id', '--name-only', '-r', sha)
    if rc != 0:
        raise ValueError('git diff-tree failed for %s: %s' % (sha, err))
    return [p for p in out.splitlines() if p.strip()]


def parse_paths_ignore(workflow_text=None):
    """The `on.push.paths-ignore` glob list, read out of test.yml itself
    (rule 1: this is the only copy; a hardcoded second list would drift).
    Regex, not a YAML parser -- consistent with test_ci_dependencies.py and
    test_ci_workflow_coverage.py's existing treatment of this same file, and
    avoids adding a yaml dependency nothing else in testsys/ declares to CI.
    """
    if workflow_text is None:
        workflow_text = open(WORKFLOW_PATH, errors='replace').read()
    m = re.search(r'paths-ignore:\s*\n((?:\s*-\s*.+\n)+)', workflow_text)
    if not m:
        raise ValueError('no on.push.paths-ignore block found in %s -- did '
                         'the workflow change shape?' % WORKFLOW_PATH)
    patterns = re.findall(r"-\s*'([^']+)'", m.group(1))
    if not patterns:
        raise ValueError('paths-ignore block in %s parsed to zero patterns'
                         % WORKFLOW_PATH)
    return patterns


def parse_workflow_name(workflow_text=None):
    """The workflow's own `name:` -- used to ask `gh` for only this
    workflow's runs, so a push-triggered run of some OTHER workflow (this
    repo also has 'Publish EQdyna Docker image' on the same push event)
    never gets counted as CI evidence."""
    if workflow_text is None:
        workflow_text = open(WORKFLOW_PATH, errors='replace').read()
    m = re.search(r'^name:\s*(.+?)\s*$', workflow_text, re.M)
    if not m:
        raise ValueError('no top-level `name:` found in %s' % WORKFLOW_PATH)
    return m.group(1)


def path_is_ignored(path, patterns):
    for pat in patterns:
        if pat.endswith('/**'):
            prefix = pat[:-3]
            if path == prefix or path.startswith(prefix + '/'):
                return True
        elif fnmatch.fnmatch(path, pat):
            return True
    return False


def commit_is_paths_ignore_only(sha, patterns=None):
    """(bool, changed_paths) -- True iff every path this commit touches is
    covered by paths-ignore, i.e. this exact push could never have started
    its own CI run."""
    if patterns is None:
        patterns = parse_paths_ignore()
    changed = commit_changed_paths(sha)
    if not changed:
        raise ValueError('%s changes no paths relative to its parent' % sha)
    ignored = all(path_is_ignored(p, patterns) for p in changed)
    return ignored, changed


def _gh_run_list(workflow_name, limit=300):
    if subprocess.run(['which', 'gh'], capture_output=True).returncode != 0:
        raise GhUnavailable('gh is not installed')
    r = subprocess.run(
        ['gh', 'run', 'list', '--workflow', workflow_name, '--limit', str(limit),
         '--json', 'databaseId,headSha,conclusion,status,createdAt,headBranch'],
        cwd=ROOT, capture_output=True, text=True)
    if r.returncode != 0:
        raise GhUnavailable('gh run list failed (rc=%d): %s'
                            % (r.returncode, (r.stderr or '').strip()[:200]))
    import json
    try:
        return json.loads(r.stdout)
    except ValueError as exc:
        raise GhUnavailable('gh run list returned unparseable JSON: %s' % exc)


def find_ci_runs(sha, workflow_name=None, limit=300):
    """All runs of `workflow_name` (default: this repo's test workflow, read
    from test.yml) whose headSha is exactly `sha`. Deliberately NOT
    `gh run list --commit <sha>`: verified empty on gh 2.91.0 for SHAs that
    demonstrably have runs (dfee14d, 813b952, both confirmed via `gh run
    view <id>` on 2026-09-21) -- that flag does not work on this `gh`
    version against this repo, so this function lists by workflow and
    filters on headSha itself instead.
    """
    if workflow_name is None:
        workflow_name = parse_workflow_name()
    runs = _gh_run_list(workflow_name, limit=limit)
    return [r for r in runs if r.get('headSha') == sha]


_TAG_REF_CACHE = {}


def is_tag_ref(name):
    """True iff `name` is a TAG on origin, not a branch.

    Needed because `gh run list --json headBranch` reports the same short
    ref name for a tag push and a branch push -- for the tag `v5.13.1`
    pushed at `dfee14d`, the run's `headBranch` is literally `"v5.13.1"`,
    indistinguishable by name alone from a branch called that. Checked
    against the remote because that is where "is this ref a tag" is
    actually decided; `git ls-remote --exit-code origin refs/tags/<name>`
    exits 0 if it is, 2 if it is not, and anything else means the network
    call itself failed (raised as GhUnavailable, not treated as "not a
    tag" -- the same reasoning check_network_side already applies to a
    failed remote round-trip).
    """
    if name in _TAG_REF_CACHE:
        return _TAG_REF_CACHE[name]
    r = subprocess.run(
        ['git', '-C', ROOT, 'ls-remote', '--exit-code', 'origin', 'refs/tags/' + name],
        capture_output=True, text=True)
    if r.returncode == 0:
        result = True
    elif r.returncode == 2:
        result = False
    else:
        raise GhUnavailable('git ls-remote origin refs/tags/%s failed (rc=%d): %s'
                            % (name, r.returncode, (r.stderr or '').strip()[:200]))
    _TAG_REF_CACHE[name] = result
    return result


def drop_tag_triggered_runs(runs):
    """Runs whose headBranch is actually a tag ref, not a branch, removed.

    A tag push fires its own `push`-triggered run (test_release_complete.py's
    check_tag_is_annotated docstring documents the same mechanism from a
    different angle: v5.8.2's tag push produced run 35122388271, twelve
    minutes after and independent of the branch-push run). Pre-tag evidence
    that a SHA is safe to tag can never legitimately come from the run the
    tagging itself produced -- that run cannot exist yet at the moment this
    check needs to answer. Counting it anyway is precisely how `dfee14d` read
    as green: its ONLY run has headBranch `v5.13.1`, the tag being created,
    not a branch.
    """
    return [r for r in runs if not is_tag_ref(r.get('headBranch', ''))]


def classify_runs(runs, exclude_run_id=None):
    """'PASS' / 'FAIL' / 'IN_PROGRESS' / 'NONE' for a list of runs already
    filtered to one SHA and one workflow.

    A completed success anywhere in the list is enough for PASS even if an
    earlier attempt for the same SHA failed (a green re-run IS evidence).
    exclude_run_id drops a specific run (the self-referential in-flight run
    examining its own tag push) before classifying the rest.
    """
    if exclude_run_id is not None:
        runs = [r for r in runs if r.get('databaseId') != exclude_run_id]
    completed = [r for r in runs if r.get('status') == 'completed']
    if any(r.get('conclusion') == 'success' for r in completed):
        return 'PASS'
    if completed:
        return 'FAIL'
    if runs:
        return 'IN_PROGRESS'
    return 'NONE'


def evaluate_pretag(sha, ack_paths_ignored_parent=False, max_hops=10):
    """The pre-`git tag` gate. Returns a Result whose .status is one of
    PASS / FAIL / PENDING / PATHS_IGNORED / UNVERIFIED (never anything
    else -- rule 2: three genuinely different situations, three outcomes).

    Default (ack_paths_ignored_parent=False): a SHA that cannot trigger its
    own CI run is PATHS_IGNORED, full stop -- it does not silently fall back
    to the parent's evidence. Passing ack_paths_ignored_parent=True walks up
    parents (each one required to ALSO be paths-ignore-only, capped at
    max_hops) until it finds a commit with runs of its own, and reports
    PASS/FAIL/PENDING against THAT ancestor -- with evidence_sha naming
    exactly which commit the evidence came from, never the SHA that was
    asked about.

    Runs triggered by pushing a TAG at this SHA are excluded from evidence
    (see drop_tag_triggered_runs) -- that run is the one the tag creation
    itself produces and cannot exist yet at the moment this check needs to
    answer "is it safe to create the tag". Counting it is the exact defect
    this module exists to close.
    """
    try:
        full_sha = resolve_sha(sha)
    except ValueError as exc:
        return Result(FAIL, str(exc))

    try:
        workflow_name = parse_workflow_name()
        patterns = parse_paths_ignore()
    except (GhUnavailable, ValueError) as exc:
        return Result(UNVERIFIED, 'could not read %s: %s' % (WORKFLOW_PATH, exc))

    current = full_sha
    hops = []
    for _ in range(max_hops + 1):
        try:
            runs = find_ci_runs(current, workflow_name)
            runs = drop_tag_triggered_runs(runs)
        except GhUnavailable as exc:
            return Result(UNVERIFIED,
                          'could not query CI runs for %s: %s' % (current, exc))
        status = classify_runs(runs)
        if status == 'PASS':
            return Result(
                PASS,
                'a completed, successful %s run exists for %s'
                % (workflow_name, current) +
                ('' if not hops else ' (evidence for %s, reached by walking '
                                     'up %d paths-ignore-only commit(s): %s)'
                                     % (full_sha, len(hops),
                                        ' -> '.join(h[0] for h in hops))),
                evidence_sha=current, hops=hops)
        if status == 'FAIL':
            return Result(
                FAIL,
                'a completed %s run for %s finished WITHOUT success '
                '(conclusion != success) -- CI ran and failed'
                % (workflow_name, current),
                evidence_sha=current, hops=hops)
        if status == 'IN_PROGRESS':
            return Result(
                PENDING,
                'a %s run for %s exists but has not completed yet -- CI has '
                'not finished' % (workflow_name, current),
                evidence_sha=current, hops=hops)
        # status == 'NONE': no run at all for `current`.
        try:
            ignored, changed = commit_is_paths_ignore_only(current, patterns)
        except ValueError as exc:
            return Result(
                PENDING,
                'no %s run found yet for %s, and its paths-ignore status '
                'could not be determined (%s) -- push it and wait'
                % (workflow_name, current, exc))
        if not ignored:
            return Result(
                PENDING,
                'no %s run found yet for %s, which DOES touch non-ignored '
                'paths and so can trigger its own run -- push it and wait, '
                'do not tag' % (workflow_name, current))
        # This commit cannot trigger CI by design.
        if current == full_sha and not ack_paths_ignored_parent:
            return Result(
                PATHS_IGNORED,
                '%s touches ONLY paths-ignore\'d files (%s) and so could '
                'NEVER have its own CI run -- the workflow\'s paths-ignore '
                'list is %s. Tagging it would rest on its PARENT\'s run, '
                'which this check does not do by default. Re-run with '
                '--ack-paths-ignored-parent to accept the nearest ancestor '
                'commit\'s green run as the evidence instead.'
                % (current, ', '.join(changed), patterns))
        hops.append((current, changed))
        parent = commit_parent(current)
        if parent is None:
            return Result(
                FAIL,
                '%s is a root commit with no parent to fall back to, and it '
                'is itself paths-ignore-only -- there is no CI evidence '
                'this chain can ever produce' % current)
        current = parent
    return Result(
        FAIL,
        'walked %d paths-ignore-only ancestors from %s without finding a '
        'commit that could have its own CI run -- max_hops exceeded, this '
        'is almost certainly a workflow paths-ignore list that is too broad'
        % (max_hops, full_sha))
