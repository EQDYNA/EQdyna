#! /usr/bin/env python3
"""
Regression guard: the repo ROOT holds only an allowlisted set of TRACKED
entries (PROJECT_RULES rule 1, owner 2026-09-24: "follow the rules").

Incident: ten mission notes (NOTES_*.md) were committed at the top level
between 2026-09-17 and 2026-09-24 and nothing caught it; they now live in
docs/notes/. A new tracked root entry fails here by name, with where it
belongs. Adding a genuinely new root entry is a deliberate edit to ALLOWED.

Both ways (rule 14a): before checking the real tree, the checker is run on a
synthetic listing with one stray NOTES file and one stray evidence file and
must flag exactly those two; a checker that cannot go red fails first.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ALLOWED = frozenset({
    'case_input', 'CLAUDE.md', 'Dockerfile', 'Docker.guide.md', '.dockerignore',
    'docs', '.github', '.gitignore', 'install-eqdyna.sh', 'LICENSE',
    'pastReleaseNotes.md', 'pathway_forward.md', 'PROJECT_RULES.md',
    'README.md', 'scripts', 'src', 'testNameList.py', 'test.reference.results',
    'testsys', 'ubuntu.env.sh', 'VERSION',
})


def where_it_belongs(name):
    if name.startswith('NOTES_'):
        return 'mission notes go to docs/notes/'
    if name.startswith('SESSION_LOG'):
        return 'session logs go to docs/'
    return ('evidence goes to docs/evidence/, snapshots to docs/perf_snapshots/, '
            'scratch stays untracked (scratch/ is gitignored)')


def stray(top_level_entries):
    return sorted(set(top_level_entries) - ALLOWED)


def tracked_top_level():
    out = subprocess.run(['git', 'ls-files', '-z'], cwd=ROOT, capture_output=True,
                         text=True, check=True).stdout
    return {p.split('/', 1)[0] for p in out.split('\0') if p}


def main():
    synthetic = set(ALLOWED) | {'NOTES_row999.md', 'evidence_dump.json'}
    got = stray(synthetic)
    if got != ['NOTES_row999.md', 'evidence_dump.json']:
        print('FAIL test_root_allowlist: negative control flagged %r, want the 2 '
              'synthetic strays -- the checker cannot go red' % got)
        return 1
    entries = tracked_top_level()
    if not entries:
        print('FAIL test_root_allowlist: git ls-files listed nothing -- not a checkout?')
        return 1
    bad = stray(entries)
    if bad:
        print('FAIL test_root_allowlist: %d tracked root entr%s outside the allowlist:'
              % (len(bad), 'y' if len(bad) == 1 else 'ies'))
        for b in bad:
            print(' - %s: %s' % (b, where_it_belongs(b)))
        return 1
    print('SUCCESS test_root_allowlist: %d tracked root entries, all allowlisted; '
          'negative control flagged 2 of 2 synthetic strays' % len(entries))
    return 0


if __name__ == '__main__':
    sys.exit(main())
