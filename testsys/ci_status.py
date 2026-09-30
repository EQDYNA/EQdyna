#! /usr/bin/env python3
"""
CI-run evidence for a commit SHA -- shared by the pre-tag gate
(testsys/regression/check_pretag_ci.py) and the post-hoc release guard
(testsys/regression/test_release_complete.py's check_ci_green_for_tagged_sha).

WHY THIS EXISTS (2026-09-21 incident): tag `v5.13.1` was pushed at `b3697f8`
while the only CI run for that exact SHA was the run the tag push itself
triggered -- at `git tag` time there was no completed run for `b3697f8` at
all. It happened to conclude green. Rule 15 step 6/7 says wait for CI to be
green on the pushed COMMIT, then tag; nothing mechanical enforced that
ordering, and nothing would have noticed had the tag-triggered run gone red,
because a pushed tag cannot be re-pointed (rule 8).

The commit that exposed this was `b3697f8`, which touches ONLY
`pathway_forward.md` -- one of `.github/workflows/test.yml`'s own
`paths-ignore` entries -- so it could never have had a pre-tag run of its
own; waiting for one would wait forever. That is a distinct situation from
"CI ran and failed" and from "CI is still running": rule 2 says those three
must not share an exit code, so this module classifies them separately
rather than folding the paths-ignore case into a generic "no evidence yet".

Everything network-touching here can fail to mean either PASS or FAIL --
`gh` missing, unauthenticated, or the network down. Per the precedent in
test_release_complete.py's check_tag_pushed_to_remote (split, along with its
sibling `check_github_release_published`, out of what was one function,
`check_network_side`, 2026-09-30),
"I could not check this" gets its own outcome (UNVERIFIED), never PASS and
never FAIL by default.

`gh run list --commit <sha> --json ...` was tested against this repo's real
history on gh 2.91.0 and returned `[]` for SHAs (`b3697f8`, `8f6ff07`) that
demonstrably have runs (confirmed via `gh run view <id>`). So this module
never uses `--commit`: it lists runs per workflow name and filters on
`headSha` itself.

THE v5.20.2 INCIDENT (2026-09-30), closed by this module's second hardening
pass. `v5.20.2` was tagged at `30ae29e`, a commit touching ONLY
`docs/evidence/**`, `docs/perf_ledger.jsonl`, `docs/perf_snapshots/**` and
`docs/run_profiles.jsonl` -- every one of them then covered by test.yml's
blanket `docs/**` paths-ignore entry, so the commit could never trigger its
own push-run. `check_pretag_ci.py --ack-paths-ignored-parent` walked up to
`30ae29e`'s parent (`e0c407f`), found ITS green run, and let the tag land on
that ancestor's evidence. `.github/workflows/test.yml`'s own comment already
said "paths-ignore is not evaluated for tag pushes" -- true GitHub Actions
semantics -- so pushing the `v5.20.2` tag DID trigger a real, full CI run for
`30ae29e` itself, which then failed at `test_term_axis.py` (a fix for it
landed later, `850ce6a`). The docs-evidence commit was genuinely
test-relevant content (test_term_axis.py reads the evidence JSON its
release-sweep-carry-forward logic writes), misclassified as inert
"docs-only", and tagged on the strength of an ancestor's run that never
exercised its own tree.

Two changes close this permanently: (1) `.github/workflows/test.yml`'s
`paths-ignore` no longer blanket-ignores `docs/**` -- the exact evidence/
ledger paths above are narrowed OUT of it, so a commit touching only them
now triggers its own push-run, same as any other commit; and (2)
`evaluate_pretag`'s `--ack-paths-ignored-parent` no longer authorizes a
PASS under any circumstance -- an ancestor's green run is never evidence for
the SHA being tagged, full stop. The flag still walks ancestors and reports
what it finds, purely as a diagnostic for the human reading the refusal
message; see `_diagnose_paths_ignored_ancestor` and `evaluate_pretag` below.
`evaluate_image_workflow` is the OTHER half of the hardening: even with (1)
and (2), a completed green `test.yml` run for the exact SHA does not by
itself prove the Docker image workflow (`publish.yml`) has ever run for
it -- that workflow triggers only on a tag push or `workflow_dispatch`, so
nothing before this ran it pre-tag at all. Item 78's `classify_every_workflow`
already refuses a SHA where a workflow that DID run finished unsuccessfully;
`evaluate_image_workflow` refuses the stronger, previously-silent case of a
SHA where it never ran in the first place.
"""
import fnmatch
import os
import re
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKFLOW_PATH = os.path.join(ROOT, '.github', 'workflows', 'test.yml')
# The Docker image build-and-push workflow -- item (b)'s gate reads ITS
# `name:` the same way parse_workflow_name reads test.yml's, never a second,
# hardcoded copy of the literal string.
PUBLISH_WORKFLOW_PATH = os.path.join(ROOT, '.github', 'workflows', 'publish.yml')

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


def parse_workflow_name(workflow_text=None, path=None):
    """The workflow's own `name:` -- used to ask `gh` for only this
    workflow's runs, so a push-triggered run of some OTHER workflow (this
    repo also has 'Publish EQdyna Docker image' on the same push event)
    never gets counted as CI evidence.

    `path` (default WORKFLOW_PATH, i.e. test.yml): which workflow file to
    read `name:` out of when `workflow_text` is not supplied directly --
    `evaluate_image_workflow` below passes PUBLISH_WORKFLOW_PATH so the two
    callers share this one regex rather than keeping a second copy."""
    if path is None:
        path = WORKFLOW_PATH
    if workflow_text is None:
        workflow_text = open(path, errors='replace').read()
    m = re.search(r'^name:\s*(.+?)\s*$', workflow_text, re.M)
    if not m:
        raise ValueError('no top-level `name:` found in %s' % path)
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
    """Runs of `workflow_name`, or of EVERY workflow when it is None
    (item 78; `workflowName` is in the JSON so callers can group)."""
    if subprocess.run(['which', 'gh'], capture_output=True).returncode != 0:
        raise GhUnavailable('gh is not installed')
    cmd = ['gh', 'run', 'list', '--limit', str(limit), '--json',
           'databaseId,headSha,conclusion,status,createdAt,headBranch,workflowName']
    if workflow_name is not None:
        cmd[3:3] = ['--workflow', workflow_name]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
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
    demonstrably have runs (b3697f8, 8f6ff07, both confirmed via `gh run
    view <id>` on 2026-09-21) -- that flag does not work on this `gh`
    version against this repo, so this function lists by workflow and
    filters on headSha itself instead.
    """
    if workflow_name is None:
        workflow_name = parse_workflow_name()
    runs = _gh_run_list(workflow_name, limit=limit)
    return [r for r in runs if r.get('headSha') == sha]


def find_all_workflow_runs(sha, limit=300):
    """Runs of EVERY workflow (not only test.yml) whose headSha is `sha`."""
    return [r for r in _gh_run_list(None, limit=limit) if r.get('headSha') == sha]


def classify_every_workflow(runs, required_workflow, exclude_run_id=None):
    """Item 78: "test.yml is green" is not "CI is green for this SHA" --
    v5.16.0 passed rule 15a honestly while 'Publish EQdyna Docker image'
    failed on the same SHA (run 35819232914) and no image was published.

    `runs` is every workflow's runs for ONE sha. Groups them by workflowName
    and classifies each group with classify_runs. Returns
    (failed, pending, names): workflow names whose runs completed without
    any success, names with no completed run yet, and every name seen.
    Raises GhUnavailable when a run carries no workflowName, or when
    `required_workflow` is absent from the listing although its own
    per-workflow query found a run: the all-workflow listing's window is
    shorter than the per-workflow one, and an incomplete listing must not
    read as "no other workflow ran".
    """
    if exclude_run_id is not None:
        runs = [r for r in runs if r.get('databaseId') != exclude_run_id]
    groups = {}
    for r in runs:
        name = r.get('workflowName')
        if not name:
            raise GhUnavailable('run %s carries no workflowName -- cannot tell '
                                'which workflow it belongs to' % r.get('databaseId'))
        groups.setdefault(name, []).append(r)
    if required_workflow is not None and required_workflow not in groups:
        raise GhUnavailable('the all-workflow run listing holds no %r run for '
                            'this sha although the per-workflow query did -- '
                            'the listing window is too short to judge the '
                            'other workflows' % required_workflow)
    failed, pending = [], []
    for name in sorted(groups):
        status = classify_runs(groups[name])
        if status == 'FAIL':
            failed.append(name)
        elif status == 'IN_PROGRESS':
            pending.append(name)
    return failed, pending, sorted(groups)


_TAG_REF_CACHE = {}


def is_tag_ref(name):
    """True iff `name` is a TAG on origin, not a branch.

    Needed because `gh run list --json headBranch` reports the same short
    ref name for a tag push and a branch push -- for the tag `v5.13.1`
    pushed at `b3697f8`, the run's `headBranch` is literally `"v5.13.1"`,
    indistinguishable by name alone from a branch called that. Checked
    against the remote because that is where "is this ref a tag" is
    actually decided; `git ls-remote --exit-code origin refs/tags/<name>`
    exits 0 if it is, 2 if it is not, and anything else means the network
    call itself failed (raised as GhUnavailable, not treated as "not a
    tag" -- the same reasoning test_release_complete.py's
    check_tag_pushed_to_remote already applies to a failed remote
    round-trip).
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

    A tag push fires its own `push`-triggered run (v5.8.2's tag push produced
    run 35122388271, twelve minutes after and independent of the branch-push
    run -- see `test_release_complete.py`'s module docstring for that
    incident). Pre-tag evidence
    that a SHA is safe to tag can never legitimately come from the run the
    tagging itself produced -- that run cannot exist yet at the moment this
    check needs to answer. Counting it anyway is precisely how `b3697f8` read
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


def _diagnose_paths_ignored_ancestor(full_sha, workflow_name, patterns, max_hops):
    """Walk up parents from `full_sha` PURELY TO REPORT what the nearest
    ancestor's CI status is -- this NEVER authorizes anything; see the
    v5.20.2 incident in this module's docstring. Returns (ancestor_sha or
    None, human-readable note, hops), where `hops` is the ordered list of
    (sha, changed_paths) pairs walked over, oldest-asked-about first.

    Kept as a real, callable mechanism (not deleted) because
    --ack-paths-ignored-parent is still useful for a human deciding what to
    do next -- it is only the AUTHORIZING effect that is gone.
    """
    current = full_sha
    hops = []
    for _ in range(max_hops):
        parent = commit_parent(current)
        if parent is None:
            return None, ('%s is a root commit -- no ancestor to check' % current), hops
        try:
            ignored, changed = commit_is_paths_ignore_only(current, patterns)
        except ValueError as exc:
            return None, ('could not determine %s\'s paths-ignore status: %s'
                          % (current, exc)), hops
        hops.append((current, changed))
        current = parent
        try:
            runs = drop_tag_triggered_runs(find_ci_runs(current, workflow_name))
        except GhUnavailable as exc:
            return None, ('could not query CI runs for ancestor %s: %s'
                          % (current, exc)), hops
        status = classify_runs(runs)
        if status in ('PASS', 'FAIL', 'IN_PROGRESS'):
            return current, ('nearest ancestor with CI evidence is %s (%s)'
                             % (current, status)), hops
        # status == 'NONE': keep walking only if this ancestor is ALSO
        # paths-ignore-only; otherwise it could have its own run (just
        # hasn't yet) and is not itself part of the paths-ignore chain.
        try:
            ignored, _changed = commit_is_paths_ignore_only(current, patterns)
        except ValueError:
            return None, ('ancestor %s has no run and its paths-ignore '
                          'status could not be determined' % current), hops
        if not ignored:
            return None, ('ancestor %s has no run yet, but it CAN trigger '
                          'its own (it is not paths-ignore-only)' % current), hops
    return None, ('walked %d paths-ignore-only ancestors from %s without '
                 'finding one with CI evidence' % (max_hops, full_sha)), hops


def evaluate_pretag(sha, ack_paths_ignored_parent=False, max_hops=10):
    """The pre-`git tag` gate. Returns a Result whose .status is one of
    PASS / FAIL / PENDING / PATHS_IGNORED / UNVERIFIED (never anything
    else -- rule 2: five genuinely different situations, five outcomes).

    HARDENED 2026-09-30 (the v5.20.2 incident -- see module docstring):
    ancestor evidence NEVER authorizes a PASS for `sha`, under any
    circumstance. A completed, successful workflow run must exist FOR THE
    EXACT SHA being asked about, full stop. `ack_paths_ignored_parent` no
    longer changes the outcome of this function at all when `sha` is itself
    paths-ignore-only: the result is PATHS_IGNORED either way. What the flag
    still does is walk the ancestor chain and fold a diagnostic note about
    the nearest ancestor's own CI status into the message -- for a human
    deciding what to do next, never for the check to decide FOR them.

    Runs triggered by pushing a TAG at this SHA are excluded from evidence
    (see drop_tag_triggered_runs) -- that run is the one the tag creation
    itself produces and cannot exist yet at the moment this check needs to
    answer "is it safe to create the tag". Counting it is the exact defect
    this module exists to close (the ORIGINAL incident, `b3697f8`/v5.13.1).
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

    try:
        runs = find_ci_runs(full_sha, workflow_name)
        runs = drop_tag_triggered_runs(runs)
    except GhUnavailable as exc:
        return Result(UNVERIFIED,
                      'could not query CI runs for %s: %s' % (full_sha, exc))

    status = classify_runs(runs)
    if status == 'PASS':
        # Item 78: every OTHER workflow with a non-tag run at this exact sha
        # must be green too. Tag-triggered runs stay excluded here for the
        # same reason as above (they cannot exist before the tag);
        # test_release_complete.py reads those after the tag.
        try:
            every = drop_tag_triggered_runs(find_all_workflow_runs(full_sha))
            failed, pending, names = classify_every_workflow(every, workflow_name)
        except GhUnavailable as exc:
            return Result(UNVERIFIED,
                          'could not check every workflow for %s: %s'
                          % (full_sha, exc))
        if failed:
            return Result(
                FAIL,
                '%s is green for %s, but workflow(s) %s also ran for that '
                'sha and finished WITHOUT success -- every workflow '
                'triggered for the sha must be green (item 78)'
                % (workflow_name, full_sha, ', '.join(repr(n) for n in failed)),
                evidence_sha=full_sha)
        if pending:
            return Result(
                PENDING,
                '%s is green for %s, but workflow(s) %s for that sha have '
                'not completed yet -- wait' % (
                    workflow_name, full_sha, ', '.join(repr(n) for n in pending)),
                evidence_sha=full_sha)
        return Result(
            PASS,
            'a completed, successful %s run exists for %s, and every '
            'workflow with a non-tag run at that sha is green (%d: %s)'
            % (workflow_name, full_sha, len(names), ', '.join(names)),
            evidence_sha=full_sha)
    if status == 'FAIL':
        return Result(
            FAIL,
            'a completed %s run for %s finished WITHOUT success '
            '(conclusion != success) -- CI ran and failed'
            % (workflow_name, full_sha),
            evidence_sha=full_sha)
    if status == 'IN_PROGRESS':
        return Result(
            PENDING,
            'a %s run for %s exists but has not completed yet -- CI has '
            'not finished' % (workflow_name, full_sha),
            evidence_sha=full_sha)

    # status == 'NONE': no run at all for `full_sha` itself.
    try:
        ignored, changed = commit_is_paths_ignore_only(full_sha, patterns)
    except ValueError as exc:
        return Result(
            PENDING,
            'no %s run found yet for %s, and its paths-ignore status '
            'could not be determined (%s) -- push it and wait'
            % (workflow_name, full_sha, exc))
    if not ignored:
        return Result(
            PENDING,
            'no %s run found yet for %s, which DOES touch non-ignored '
            'paths and so can trigger its own run -- push it and wait, '
            'do not tag' % (workflow_name, full_sha))

    # full_sha cannot trigger its own CI run by design. This is PATHS_IGNORED
    # no matter what ack_paths_ignored_parent says -- an ancestor's evidence
    # is NEVER a substitute for this sha's own (the v5.20.2 fix).
    base_msg = (
        '%s touches ONLY paths-ignore\'d files (%s) and so could NEVER have '
        'its own CI run -- the workflow\'s paths-ignore list is %s. A real '
        'CI run at this EXACT sha is mandatory; there is no escape hatch '
        'that lets an ancestor\'s run stand in for it (owner decision, '
        '2026-09-30, closing the v5.20.2 incident: %s was tagged on its '
        'parent\'s green run and then failed CI for real on the tag push, '
        'because a tag push is NOT subject to paths-ignore). Push a commit '
        'that also touches a non-ignored path, or narrow the workflow\'s '
        'paths-ignore list so this one can trigger its own run.'
        % (full_sha, ', '.join(changed), patterns, full_sha))
    if not ack_paths_ignored_parent:
        return Result(PATHS_IGNORED, base_msg)

    ancestor_sha, note, hops = _diagnose_paths_ignored_ancestor(
        full_sha, workflow_name, patterns, max_hops)
    return Result(
        PATHS_IGNORED,
        base_msg + ' --ack-paths-ignored-parent no longer grants a PASS; it '
        'only reports a diagnostic: %s.' % note,
        evidence_sha=ancestor_sha, hops=hops)


def evaluate_image_workflow(sha, workflow_name=None):
    """Owner requirement (b), 2026-09-30 (v5.20.2 incident hardening): does
    the Docker image-publish workflow (.github/workflows/publish.yml) have a
    completed, GREEN run for the EXACT `sha` -- not "if it ran, is it green"
    (item 78's `classify_every_workflow`, which only judges workflows that
    DID run for a sha), but the stronger, previously-unchecked question: HAS
    it run at all. A sha where it never ran is refused exactly as if it had
    failed; there is no silent pass for "never scheduled".

    `publish.yml` triggers only on a tag push or `workflow_dispatch` (see
    that file's own `on:` block) -- it never fires on an ordinary push to
    master. So before any tag exists, the ONLY way a completed green run of
    this workflow can exist for a given sha is a `workflow_dispatch` run
    against it, done by hand ahead of time. This function does not care HOW
    the run got triggered, only that one exists and is green; tag-triggered
    runs are not excluded here (unlike evaluate_pretag's test.yml check)
    because pre-tag there cannot be one yet for this sha -- if one somehow
    already exists (re-verifying a sha that was tagged before), it is real
    evidence, not the same self-referential trap drop_tag_triggered_runs
    exists to close for test.yml.

    Returns a Result whose .status is PASS / FAIL / PENDING / UNVERIFIED
    (never PATHS_IGNORED -- publish.yml carries no paths-ignore filter of
    its own to trip)."""
    if workflow_name is None:
        workflow_name = parse_workflow_name(path=PUBLISH_WORKFLOW_PATH)
    try:
        runs = find_ci_runs(sha, workflow_name)
    except GhUnavailable as exc:
        return Result(UNVERIFIED,
                      'could not query %s runs for %s: %s' % (workflow_name, sha, exc))
    status = classify_runs(runs)
    if status == 'PASS':
        return Result(PASS,
                      'a completed, successful %s run exists for %s'
                      % (workflow_name, sha), evidence_sha=sha)
    if status == 'FAIL':
        return Result(
            FAIL,
            'a completed %s run for %s finished WITHOUT success -- the image '
            'build, its gate (the repo\'s own test tiers run INSIDE the '
            'built image), or the push to the registry failed'
            % (workflow_name, sha), evidence_sha=sha)
    if status == 'IN_PROGRESS':
        return Result(
            PENDING,
            'a %s run for %s exists but has not completed yet -- wait'
            % (workflow_name, sha), evidence_sha=sha)
    # status == 'NONE': the workflow has never run for this sha at all.
    return Result(
        FAIL,
        'no %s run exists for %s at all -- this workflow triggers only on a '
        'tag push or workflow_dispatch, so a sha that has never been '
        'dispatched has no evidence it can even build the image, let alone '
        'that the image passes its own gate. `gh workflow run` takes a '
        'BRANCH OR TAG name for --ref, never a commit sha (rule 24a audit, '
        '2026-09-30: an earlier draft of this message said `--ref %s`, '
        'which `gh` rejects outright) -- dispatch it while %s is still the '
        'tip of the branch you intend to tag, e.g. `gh workflow run %r '
        '--ref master` run immediately after merging, before anything else '
        'lands on master, and wait for it to go green before tagging -- a '
        'sha where the image workflow never ran is refused, not silently '
        'passed (the v5.20.2 gap this closes: that tag shipped with a 404 '
        'GHCR manifest and nothing here said so beforehand).'
        % (workflow_name, sha, sha, sha, workflow_name))
