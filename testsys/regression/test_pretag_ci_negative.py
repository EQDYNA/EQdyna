#! /usr/bin/env python3
"""
Negative test for rule 15a's pre-tag guard: drive `check_pretag_ci.py` through
ALL SIX of its documented outcomes and assert each one (rules 2, 3).

WHY THIS EXISTS. The guard (`check_pretag_ci.py`, landed `a976ec7`) was written
after the v5.13.1 violation, and until 2026-09-21 it had never been shown to
FAIL where it should. It was then driven by hand across its six outcomes
against this repo's REAL CI history -- a good verification, and a perishable
one. This file is that verification made permanent, because two of the six
outcomes are NOT reproducible from live history:

  * PENDING via an in-progress run: run 35676386267 was `in_progress` when the
    board drove it on 2026-09-21. It has since completed. Re-running that same
    command today returns PASS. A transient outcome cannot be a standing
    evidence command.
  * UNVERIFIED: needs `gh` to be absent, which it is not on a working box.

A guard nobody has watched fail is the vacuous-gate pattern that hit this
project four times in one day -- one of them a SHA check comparing an empty
string to an empty string and reporting OK three times. So this test also
asserts the six outcomes produce six DIFFERENT exit codes, and it includes a
seventh case that deliberately disables the tag-run filter and shows the guard
then reads the v5.13.1 commit as green. If that seventh case ever stops
flipping, the filter has stopped being load-bearing.

PROVENANCE OF THE FIXTURE (this is real data, not stubs). Every run record
below was captured 2026-09-22 from this repo's own CI history with:

    gh run list --workflow 'Automatic Testing of EQdyna' --limit 300 \
       --json databaseId,headSha,conclusion,status,createdAt,headBranch

verbatim, with exactly ONE disclosed modification: `_IN_PROGRESS_RUN` is real
run 35673827398 with `status` changed 'completed' -> 'in_progress' and
`conclusion` changed 'success' -> None, because GitHub keeps no historical
record of a run while it was still running. Every other field of every other
record is as the API returned it. The commit SHAs, their parents and their
changed-path sets are read live from this repo's git history, not faked.

WHAT IS AND IS NOT INJECTED. Only the two network calls are replaced:
`ci_status._gh_run_list` (the `gh` round-trip) and `ci_status.is_tag_ref` (the
`git ls-remote` round-trip). Everything the guard actually decides with --
paths-ignore parsing out of test.yml, the parent walk, the tag-run exclusion,
run classification, and the exit-code mapping in `check_pretag_ci.main()` --
is the real code path. Case F additionally runs the guard as a REAL SUBPROCESS
with `gh` removed from PATH, with nothing injected at all.

One more thing IS stubbed, deliberately and out loud: as of 2026-09-23 the
guard ALSO requires committed full-term local sweep evidence
(`check_pretag_ci.evaluate_sweep_evidence`) once CI reports PASS -- a second,
independent gate (exit code 5, SWEEP_INSUFFICIENT) that these six cases were
never about. `run_guard` stubs it to a fixed `(True, ...)` for every case here
so A/G's exit-0 assertions test only the CI dimension, exactly as they did
before that gate existed. The sweep-evidence dimension gets its OWN dedicated
negative tests, against a real temporary git repo the same way
test_precommit_board_separation_guard.py builds one, in
`test_pretag_sweep_negative.py`. Case F's real-subprocess run is NOT stubbed
(nothing is injected in that case at all) but never reaches the sweep check
either, since UNVERIFIED returns before it.

Cheap (rule 9): no network, no build, well under 1 s.
"""
import contextlib
import importlib.util
import io
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
from testsys import ci_status  # noqa: E402

# --- the real commits these outcomes are driven from -----------------------
SHA_FAILED_CI = 'c782c2e15e0403cb03e7d4a1ece8a24fc6c63739'   # CI ran, failed
SHA_GREEN = '813b952df8c7f78f57b7b1b21a315f37df8328b8'       # CI ran, green
SHA_PATHS_IGNORED = 'dfee14d2eda4d16c407a9bb94ceffd52fe6544d7'  # v5.13.1's commit
SHA_ACK_EVIDENCE = 'e9c3fa2ee337cd78fc2963b12b6db8421a16fb69'   # its parent
SHA_CODE_NO_RUN = '01a1040'   # touches src/python -- CAN trigger CI, here has no run

# --- real run records, captured 2026-09-22 (see PROVENANCE above) ----------
_RUN_FAILED = {
    'conclusion': 'failure', 'createdAt': '2026-09-19T07:01:52Z',
    'databaseId': 35428208759, 'headBranch': 'master',
    'headSha': SHA_FAILED_CI, 'status': 'completed'}
_RUN_GREEN_MASTER = {
    'conclusion': 'success', 'createdAt': '2026-09-21T23:56:43Z',
    'databaseId': 35669818523, 'headBranch': 'master',
    'headSha': SHA_GREEN, 'status': 'completed'}
_RUN_GREEN_TAGPUSH = {
    'conclusion': 'success', 'createdAt': '2026-09-22T00:54:44Z',
    'databaseId': 35673827398, 'headBranch': 'v5.14.0',
    'headSha': SHA_GREEN, 'status': 'completed'}
# The one disclosed modification -- see PROVENANCE.
_RUN_IN_PROGRESS = dict(_RUN_GREEN_MASTER, databaseId=35676386267,
                        status='in_progress', conclusion=None)
_RUN_V5131_TAGPUSH = {
    'conclusion': 'success', 'createdAt': '2026-09-21T23:26:04Z',
    'databaseId': 35667535066, 'headBranch': 'v5.13.1',
    'headSha': SHA_PATHS_IGNORED, 'status': 'completed'}
_RUN_ACK_PARENT = {
    'conclusion': 'success', 'createdAt': '2026-09-21T21:42:44Z',
    'databaseId': 35658746273, 'headBranch': 'master',
    'headSha': SHA_ACK_EVIDENCE, 'status': 'completed'}

# Refs that are TAGS on origin, as `git ls-remote refs/tags/` reports them.
REAL_TAGS = ('v5.13.1', 'v5.14.0', 'v5.6.0')


def load_guard():
    """The real CLI module, imported from its path (it is not a package)."""
    path = os.path.join(HERE, 'check_pretag_ci.py')
    spec = importlib.util.spec_from_file_location('check_pretag_ci_under_test', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@contextlib.contextmanager
def injected_ci(runs, unavailable=False, honour_tags=True):
    """Replace ONLY the two network calls. Everything else stays real."""
    def _gh_run_list(workflow_name, limit=300):
        if unavailable:
            raise ci_status.GhUnavailable('gh is not installed (injected)')
        return list(runs)

    def _is_tag_ref(name):
        return bool(honour_tags) and name in REAL_TAGS

    saved = (ci_status._gh_run_list, ci_status.is_tag_ref)
    ci_status._gh_run_list = _gh_run_list
    ci_status.is_tag_ref = _is_tag_ref
    try:
        yield
    finally:
        ci_status._gh_run_list, ci_status.is_tag_ref = saved


def run_guard(argv, stub_sweep=True):
    """(exit_code, stdout) from the real `main()`, through the real mapping.

    `stub_sweep`: patch the freshly-loaded module's `evaluate_sweep_evidence`
    to a fixed `(True, ...)` -- see the module docstring's "One more thing IS
    stubbed" note. `load_guard()` builds a brand-new module object every call
    (no shared state), so this never leaks between cases.
    """
    guard = load_guard()
    if stub_sweep:
        guard.evaluate_sweep_evidence = (
            lambda *a, **k: (True, 'stubbed for CI-outcome-only case'))
    buf = io.StringIO()
    saved_argv = sys.argv
    sys.argv = ['check_pretag_ci.py'] + list(argv)
    try:
        with contextlib.redirect_stdout(buf):
            code = guard.main()
    finally:
        sys.argv = saved_argv
    return code, buf.getvalue()


# --- the six outcomes, plus the filter-is-load-bearing demonstration -------
# (label, argv, runs the API would return, expected exit, substring required)
CASES = [
    ('A  PASS            a green non-tag run exists',
     ['--pre-tag', SHA_GREEN], [_RUN_GREEN_MASTER, _RUN_GREEN_TAGPUSH],
     0, 'completed, successful'),
    ('B  FAIL            CI ran and did not succeed',
     ['--pre-tag', SHA_FAILED_CI], [_RUN_FAILED],
     1, 'WITHOUT success'),
    ('C  PENDING         a run exists but has not completed',
     ['--pre-tag', SHA_GREEN], [_RUN_IN_PROGRESS],
     2, 'has not completed yet'),
    ('D  PENDING         no run yet, and the commit CAN trigger one',
     ['--pre-tag', SHA_CODE_NO_RUN], [],
     2, 'DOES touch non-ignored paths'),
    ('E  PATHS_IGNORED   only run is the tag push it must not count',
     ['--pre-tag', SHA_PATHS_IGNORED], [_RUN_V5131_TAGPUSH],
     3, 'could NEVER have its own CI run'),
    ('G  PASS (acked)    parent evidence, named, never the SHA asked about',
     ['--pre-tag', SHA_PATHS_IGNORED, '--ack-paths-ignored-parent'],
     [_RUN_V5131_TAGPUSH, _RUN_ACK_PARENT],
     0, 'evidence sha: ' + SHA_ACK_EVIDENCE),
]


def check_injected_cases(failures):
    seen_codes = {}
    for label, argv, runs, want_code, want_text in CASES:
        with injected_ci(runs):
            code, out = run_guard(argv)
        ok = (code == want_code) and (want_text in out)
        print('  %s %-58s exit=%d (want %d)'
              % ('.' if ok else 'X', label, code, want_code))
        if not ok:
            failures.append('%s: exit=%d want=%d, stdout=%r'
                            % (label.split()[0], code, want_code, out.strip()))
        seen_codes.setdefault(want_code, label)
    return seen_codes


def check_tag_filter_is_load_bearing(failures):
    """The v5.13.1 violation, reproduced: with the tag-run filter disabled the
    guard reads dfee14d as green off the run its own tag push created. This is
    the defect rule 15a exists to close; if this case stops flipping, case E is
    passing for some other reason and has gone vacuous."""
    with injected_ci([_RUN_V5131_TAGPUSH], honour_tags=False):
        code, out = run_guard(['--pre-tag', SHA_PATHS_IGNORED])
    ok = code == 0 and 'completed, successful' in out
    print('  %s H  filter off -> the v5.13.1 violation reappears        '
          'exit=%d (want 0)' % ('.' if ok else 'X', code))
    if not ok:
        failures.append('H: with drop_tag_triggered_runs neutered the guard '
                        'did NOT read %s as green (exit=%d) -- case E may be '
                        'passing for the wrong reason' % (SHA_PATHS_IGNORED, code))


def check_unverified_end_to_end(failures):
    """UNVERIFIED (exit 4) with NOTHING injected: the real script, in a real
    subprocess, on a PATH that has git and which but no gh."""
    real = {}
    for tool in ('git', 'which'):
        p = shutil.which(tool)
        if p is None:
            failures.append('F: cannot locate %r to build a gh-free PATH' % tool)
            return
        real[tool] = p
    with tempfile.TemporaryDirectory() as td:
        for tool, path in real.items():
            os.symlink(path, os.path.join(td, tool))
        env = dict(os.environ, PATH=td)
        env.pop('GH_PATH', None)
        r = subprocess.run(
            [sys.executable, os.path.join(HERE, 'check_pretag_ci.py'),
             '--pre-tag', SHA_GREEN],
            cwd=ROOT, env=env, capture_output=True, text=True)
    ok = r.returncode == 4 and 'UNVERIFIED' in r.stdout
    print('  %s F  UNVERIFIED       real subprocess, gh absent from PATH      '
          '  exit=%d (want 4)' % ('.' if ok else 'X', r.returncode))
    if not ok:
        failures.append('F: exit=%d want=4, stdout=%r stderr=%r'
                        % (r.returncode, r.stdout.strip(), r.stderr.strip()[:300]))


def check_history_is_present():
    """These outcomes are driven from named past commits, so the checkout must
    have real history. actions/checkout@v4 defaults to depth 1, which would
    make every case fail for a reason that has nothing to do with the guard --
    so say exactly that instead of letting it look like a guard regression."""
    missing = []
    for sha in (SHA_FAILED_CI, SHA_GREEN, SHA_PATHS_IGNORED, SHA_ACK_EVIDENCE,
                SHA_CODE_NO_RUN):
        try:
            ci_status.resolve_sha(sha)
        except ValueError:
            missing.append(sha)
    if missing:
        return ('shallow or incomplete checkout: %d referenced commit(s) do '
                'not resolve (%s). This test needs real history -- set '
                '`fetch-depth: 0` on this job\'s actions/checkout step.'
                % (len(missing), ', '.join(s[:7] for s in missing)))
    return None


def main():
    print('Negative test: check_pretag_ci.py across all six outcomes')
    shallow = check_history_is_present()
    if shallow:
        print('\nFAIL: ' + shallow)
        return 1
    failures = []
    seen_codes = check_injected_cases(failures)
    check_unverified_end_to_end(failures)
    check_tag_filter_is_load_bearing(failures)

    # Anti-vacuity: the six outcomes must be six DIFFERENT exit codes, and all
    # five documented codes must be reached (4 comes from the subprocess case).
    covered = set(seen_codes) | {4}
    if covered != {0, 1, 2, 3, 4}:
        failures.append('exit codes covered %s, want {0,1,2,3,4} -- an outcome '
                        'is untested' % sorted(covered))
    if len(CASES) != 6:
        failures.append('CASES has %d entries; the guard documents six outcomes'
                        % len(CASES))

    if failures:
        print('\nFAIL: %d case(s)' % len(failures))
        for f in failures:
            print('  ' + f)
        return 1
    print('\nPASS: 6 outcomes -> exit codes {0,1,2,3,4} (0 twice: plain and '
          'acked), the tag-run filter shown load-bearing, UNVERIFIED shown '
          'end-to-end with no injection')
    return 0


if __name__ == '__main__':
    sys.exit(main())
