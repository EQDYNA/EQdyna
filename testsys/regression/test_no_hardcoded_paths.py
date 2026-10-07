#! /usr/bin/env python3
"""
Regression guard: no hardcoded machine-local paths in the tracked tree
(owner decision 2026-10-06, PR #130 follow-on -- the repo is PUBLIC).

INCIDENT: a repo-wide audit on 2026-10-06 found absolute paths tied to one
developer's machine (home directories under two different usernames) and to
one agent session's ephemeral `/tmp` scratch directory leaked into
43 tracked files -- rule text, the board, session logs, notes, evidence, and
one piece of actually-executed code (a hardcoded scratch path in
`testsys/parity/evidence_drv_a6_marginal_population.py`). None of that is a
credential, but a PUBLIC repo should not carry a contributor's home
directory layout, and a hardcoded absolute path in functional code is a
portability bug independent of publicness. Those 43 were scrubbed by hand;
this guard is what keeps a 44th from landing unnoticed.

WHAT THIS CHECKS: `git grep -c` (not a hand-maintained file list, so it
catches anything NEW) for the same two patterns across every file git
currently tracks, this guard's own source excluded (see MECHANISM below).
Any hit is a hard FAIL -- there is no allowlist and no skip, per rule 2's
"no third state" discipline used by the rest of testsys/.

NO EXEMPTIONS: `docs/BOARD_HISTORY.md` carried four pre-existing leaked paths
inside archived cells, flagged for an owner call rather than scrubbed, since
the file declares itself a verbatim archive. The owner decided (PR #136,
"Scrub them (Recommended)") to scrub them with inline placeholders instead of
keeping a named exemption; that exemption is removed here (rule 14a: this
guard must be able to fail on a real leak in that file again, not just on
everything else).

MECHANISM (so this file does not trip its own check): the two regexes are
built from concatenated string fragments rather than spelled out whole, so
`git grep` over this file's own source does not match them.

Cheap (rule 9): two `git grep -c` subprocess calls, no build, no network.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SELF_RELPATH = 'testsys/regression/test_no_hardcoded_paths.py'

# No exemptions (see module docstring, rule 14a): every tracked file must be
# clean, including docs/BOARD_HISTORY.md since PR #136/#137.
EXEMPT = set()

# Built from fragments so this file's own source never matches its own
# check (a plain scan for the literal pattern strings would self-trip).
_TMP_SCRATCH = '/' + 'tmp/claude-'
_HOME_DIRS = '/' + 'home/(utig5|staff)/'
PATTERN = _TMP_SCRATCH + '|' + _HOME_DIRS


def git_grep_matches(pattern):
    """Return {relpath: count} for every tracked file matching pattern,
    via a real `git grep -c`, not a hand-maintained file list."""
    r = subprocess.run(
        ['git', 'grep', '-c', '-E', pattern, '--', '.'],
        cwd=ROOT, capture_output=True, text=True, timeout=60)
    # git grep exits 1 when there are zero matches (not an error here);
    # anything else (2+) is a real failure of the tool itself.
    if r.returncode not in (0, 1):
        raise RuntimeError('git grep failed (exit %d): %s'
                           % (r.returncode, r.stderr))
    hits = {}
    for line in r.stdout.splitlines():
        path, _, count = line.rpartition(':')
        if path:
            hits[path] = int(count)
    return hits


def main():
    hits = git_grep_matches(PATTERN)

    # This guard's own source legitimately contains the pattern fragments
    # assembled into a regex string -- never the literal pattern itself, so
    # it should not appear in git grep's results; if it ever does, that is
    # itself a bug in the fragment-splitting above, worth failing loudly on
    # rather than silently excluding.
    if SELF_RELPATH in hits:
        print('FAIL test_no_hardcoded_paths: this guard\'s own source '
              '(%s) matched the pattern it checks for -- the fragment '
              'splitting meant to prevent that has broken' % SELF_RELPATH)
        return 1

    unexpected = {p: c for p, c in hits.items() if p not in EXEMPT}

    if unexpected:
        print('FAIL test_no_hardcoded_paths: %d tracked file(s) carry a '
              'hardcoded machine-local or scratch path (pattern %r):'
              % (len(unexpected), PATTERN))
        for path in sorted(unexpected):
            print('  %s (%d match(es))' % (path, unexpected[path]))
        print('  Fix: replace with a placeholder (prose) or a relative '
              'path / env var (functional code) -- see pathway_forward.md '
              'for the 2026-10-06 scrub this guard enforces. No exemptions '
              '(rule 14a) -- every tracked file must be clean.')
        return 1

    exempted_present = sorted(p for p in EXEMPT if p in hits)
    print('  PASS  0 unexpected tracked files match %r via `git grep` '
          '(%d file(s) matched overall, all accounted for by %d documented '
          'exemption(s): %s)'
          % (PATTERN, len(hits), len(exempted_present),
             exempted_present or 'none'))
    print('SUCCESS test_no_hardcoded_paths')
    return 0


if __name__ == '__main__':
    sys.exit(main())
