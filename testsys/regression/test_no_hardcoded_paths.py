#! /usr/bin/env python3
"""
Regression guard: no hardcoded machine-local paths, this box's hostname, or
a maintainer's personal email in the tracked tree (owner decision 2026-10-06,
PR #130 follow-on; widened 2026-10-14, item 3a -- the repo is PUBLIC).

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

WIDENED 2026-10-14 (item 3a, board commit 9f723f2): a second audit found this
box's hostname baked into ~207 committed perf/profile artifacts (recorded for
provenance, never keyed on by code -- see the forward-commit PR for that
read) and a maintainer's personal email address in ~40 src/fortran headers,
`scripts/*.m`, and example config values. Both scrubbed by a forward commit;
this guard now also fails on either pattern reappearing.

ONE NAMED HOSTNAME EXEMPTION: `docs/perf_ledger.jsonl` is append-only by design
(`testsys/regression/test_perf_ledger.py` check 6, "history immutability") --
rewriting a past row's `host` field is editing committed history, which that
guard exists to catch, not a thing this guard should force. Scrubbed here by
NOT editing it (item 3a, 2026-10-14): every row recorded before this date
still carries the real hostname, a measurement artifact, not a credential;
every row appended from now on should prefer a neutral host label (see that
file's own writer for how). HOST_EXEMPT names this one file so the exemption
is visible and so a SECOND file carrying the hostname still fails loudly.

WHAT THIS CHECKS: `git grep -c` (not a hand-maintained file list, so it
catches anything NEW) for four patterns across every file git currently
tracks, this guard's own source excluded (see MECHANISM below). Any hit is a
hard FAIL -- there is no allowlist and no skip, per rule 2's "no third state"
discipline used by the rest of testsys/, EXCEPT the email pattern, which
carries a narrow, named exemption (EMAIL_EXEMPT below) for the three files
that are the project's single canonical author-contact declaration.

NO EXEMPTIONS (paths, hostname): `docs/BOARD_HISTORY.md` carried four
pre-existing leaked paths inside archived cells, flagged for an owner call
rather than scrubbed, since the file declares itself a verbatim archive. The
owner decided (PR #136, "Scrub them (Recommended)") to scrub them with inline
placeholders instead of keeping a named exemption; that exemption is removed
here (rule 14a: this guard must be able to fail on a real leak in that file
again, not just on everything else).

MECHANISM (so this file does not trip its own check): the regexes are built
from concatenated string fragments rather than spelled out whole, so
`git grep` over this file's own source does not match them. The hostname
fragment is this box's `hostname -s` output, kept to a single split literal
(not reproduced in this docstring) rather than a hash, because a hash would
not let a human read this file and know what it guards; the split already
keeps it out of a plain-text search of the tracked tree.

Cheap (rule 9): four `git grep -c` subprocess calls, no build, no network.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SELF_RELPATH = 'testsys/regression/test_no_hardcoded_paths.py'

# No exemptions (see module docstring, rule 14a): every tracked file must be
# clean, including docs/BOARD_HISTORY.md since PR #136/#137.
EXEMPT = set()

# The project's single canonical author-contact declaration: a maintainer
# email is expected here and must not be scrubbed (it is the one place a
# reader is told who to contact -- scrubbing it everywhere else, and keeping
# it only here, is the point of the 2026-10-14 scrub).
EMAIL_EXEMPT = {'LICENSE', 'README.md', 'docs/user/citing.md'}

# The one file the append-only perf-ledger guard (test_perf_ledger.py check 6)
# forbids rewriting past rows of -- see module docstring, "ONE NAMED HOSTNAME
# EXEMPTION" above.
HOST_EXEMPT = {'docs/perf_ledger.jsonl'}

# Built from fragments so this file's own source never matches its own
# check (a plain scan for the literal pattern strings would self-trip).
_TMP_SCRATCH = '/' + 'tmp/claude-'
_HOME_DIRS = '/' + 'home/(utig5|staff)/'
PATTERN = _TMP_SCRATCH + '|' + _HOME_DIRS

# This box's hostname (`hostname -s`), split so this file's own source does
# not match the pattern it checks for.
_HOST = 'coto' + 'paxi'
HOST_PATTERN = _HOST

# The two maintainer personal addresses scrubbed 2026-10-14, split the same
# way as the hostname above.
_EMAIL_1 = 'dliu' + '@ig.utexas.edu'
_EMAIL_2 = 'bduan' + '@tamu.edu'
EMAIL_PATTERN = _EMAIL_1 + '|' + _EMAIL_2


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
    host_hits = git_grep_matches(HOST_PATTERN)
    email_hits = git_grep_matches(EMAIL_PATTERN)

    # This guard's own source legitimately contains the pattern fragments
    # assembled into a regex string -- never the literal pattern itself, so
    # it should not appear in git grep's results; if it ever does, that is
    # itself a bug in the fragment-splitting above, worth failing loudly on
    # rather than silently excluding.
    for label, hs in (('path', hits), ('host', host_hits), ('email', email_hits)):
        if SELF_RELPATH in hs:
            print('FAIL test_no_hardcoded_paths: this guard\'s own source '
                  '(%s) matched the %s pattern it checks for -- the fragment '
                  'splitting meant to prevent that has broken'
                  % (SELF_RELPATH, label))
            return 1

    unexpected = {p: c for p, c in hits.items() if p not in EXEMPT}
    unexpected_host = {p: c for p, c in host_hits.items() if p not in HOST_EXEMPT}
    unexpected_email = {p: c for p, c in email_hits.items()
                         if p not in EMAIL_EXEMPT}

    problems = []
    if unexpected:
        problems.append(('path', PATTERN, unexpected))
    if unexpected_host:
        problems.append(('hostname', HOST_PATTERN, unexpected_host))
    if unexpected_email:
        problems.append(('email', EMAIL_PATTERN, unexpected_email))

    if problems:
        for label, pattern, bad in problems:
            print('FAIL test_no_hardcoded_paths: %d tracked file(s) carry a '
                  'leaked %s (pattern %r):' % (len(bad), label, pattern))
            for path in sorted(bad):
                print('  %s:%d match(es)' % (path, bad[path]))
        print('  Fix: replace with a placeholder (prose) or a relative '
              'path / env var (functional code) -- see pathway_forward.md '
              'for the 2026-10-06 scrub, and item 3a (2026-10-14) for the '
              'hostname/email scrub, this guard enforces. No exemptions for '
              'path/hostname (rule 14a); email is exempt only in %s.'
              % sorted(EMAIL_EXEMPT))
        return 1

    exempted_present = sorted(p for p in EXEMPT if p in hits)
    email_exempted_present = sorted(p for p in EMAIL_EXEMPT if p in email_hits)
    host_exempted_present = sorted(p for p in HOST_EXEMPT if p in host_hits)
    print('  PASS  0 unexpected tracked files match %r via `git grep` '
          '(%d file(s) matched overall, all accounted for by %d documented '
          'exemption(s): %s)'
          % (PATTERN, len(hits), len(exempted_present),
             exempted_present or 'none'))
    print('  PASS  0 unexpected tracked files match the hostname pattern '
          '(%d file(s) matched overall, all accounted for by %d documented '
          'exemption(s): %s)'
          % (len(host_hits), len(host_exempted_present),
             host_exempted_present or 'none'))
    print('  PASS  0 unexpected tracked files match the email pattern '
          '(%d file(s) matched overall, all accounted for by %d documented '
          'exemption(s): %s)'
          % (len(email_hits), len(email_exempted_present),
             email_exempted_present or 'none'))
    print('SUCCESS test_no_hardcoded_paths')
    return 0


if __name__ == '__main__':
    sys.exit(main())
