#! /usr/bin/env python3
"""
Regression guard: user-facing docs are written for users (owner 2026-09-24:
"README is user facing. This is too long to read"; "some parts format too
bad"; "need a rule for this").

Applies to README.md, Docker.guide.md, case_input/*/README.md and every page
under docs/user/. For each file:
  1. no internal markers in PROSE (outside fenced code): PR numbers, rule
     numbers, board references ("pathway", "board row", "item <n>"), agent
     names, source `file.py:<line>` / `file.f90:<line>` references;
     a backticked 7-40 hex SHA is refused anywhere, code included;
  2. exactly one H1, it is the first heading, and nothing deeper than H3;
  3. every prose bullet line is at most MAX_BULLET_CHARS long;
  4. README.md is at most README_MAX_LINES lines.
Commands and paths a user runs live in fenced code or inline code and are not
marker-checked (a `testsys/run.py` path is legitimate user content).

Both ways (rule 14a): each check first runs on synthetic documents, one per
violation, that must be flagged, and on a clean one that must pass.
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

README_MAX_LINES = 150
MAX_BULLET_CHARS = 200
AGENTS = ('mira', 'iris', 'lars', 'kai', 'haruto', 'nadia', 'sophia', 'zofia',
          'victor', 'wei-lin', 'wei lin', 'dunyu-liu', 'anya', 'marta', 'priya')
# Checked in prose, with inline code stripped (commands and paths a user runs
# are legitimate content). Agent names are matched case-SENSITIVELY in their
# lowercase handle form, so the seismic data service IRIS or the Argonne
# machine Mira in ordinary prose is not flagged.
PROSE_MARKERS = [
    ('PR number', re.compile(r'\bPR\s*#\d+|\(#\d+\)')),
    ('rule number', re.compile(r'\brules?\s+\d+[a-z]?\b', re.I)),
    ('board reference', re.compile(r'pathway_forward|board row|\bitem\s+\d+', re.I)),
    ('agent name', re.compile(r'\b(' + '|'.join(re.escape(a) for a in AGENTS) + r')\b')),
    ('bare commit SHA', re.compile(r'(?<![\w/.-])(?=[0-9a-f]*[a-f])(?=[0-9a-f]*[0-9])(?:[0-9a-f]{7,12}|[0-9a-f]{40})(?![\w/.-])')),
]
# Checked on the WHOLE line, inline code included: a file:line citation, a
# developer-notes reference or a gate-internal name is never user content,
# backticked or not.
ANYWHERE_MARKERS = [
    ('source line reference', re.compile(r'\.(py|f90|sh|md|yml|txt):\d+')),
    ('developer notes reference', re.compile(r'\bNOTES_\w+')),
    ('gate internals', re.compile(r'\b(CASE_BOUND|GATE_TERM_S|CI_CELLS|GATE_STATIONS|STATION_BOUND|RELEASE_ONLY)\b')),
]
SHA = re.compile(r'`[0-9a-f]{7,40}`')


def user_docs():
    files = ['README.md', 'Docker.guide.md']
    files += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, 'case_input', '*', 'README.md')))
    files += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, 'docs', 'user', '**', '*.md'), recursive=True))
    return [f for f in files if os.path.isfile(os.path.join(ROOT, f))]


def problems(name, text):
    bad = []
    lines = text.split('\n')
    if SHA.search(text):
        bad.append('%s: backticked commit SHA %s' % (name, SHA.search(text).group(0)))
    fence = False
    headings = []
    for n, line in enumerate(lines, 1):
        if line.lstrip().startswith('```'):
            fence = not fence
            continue
        if fence:
            continue
        prose = re.sub(r'`[^`]*`', '', line)       # inline code is user content
        for label, rx in PROSE_MARKERS:
            m = rx.search(prose)
            if m:
                bad.append('%s:%d: %s %r' % (name, n, label, m.group(0)))
        for label, rx in ANYWHERE_MARKERS:
            m = rx.search(line)
            if m:
                bad.append('%s:%d: %s %r' % (name, n, label, m.group(0)))
        h = re.match(r'(#{1,6})\s', line)
        if h:
            headings.append((n, len(h.group(1))))
        if re.match(r'\s*([*-]|\d+\.)\s', line) and len(line) > MAX_BULLET_CHARS:
            bad.append('%s:%d: bullet is %d chars (max %d)' % (name, n, len(line), MAX_BULLET_CHARS))
    h1 = [n for n, lvl in headings if lvl == 1]
    if len(h1) != 1:
        bad.append('%s: %d H1 headings (want exactly 1)' % (name, len(h1)))
    elif headings[0][1] != 1:
        bad.append('%s:%d: first heading is not the H1' % (name, headings[0][0]))
    for n, lvl in headings:
        if lvl > 3:
            bad.append('%s:%d: heading deeper than H3' % (name, n))
    if name == 'README.md' and len(lines) > README_MAX_LINES:
        bad.append('README.md: %d lines (max %d)' % (len(lines), README_MAX_LINES))
    return bad


def negative_controls():
    clean = '# Title\n\n## Install\n\n```\npython3 testsys/run.py unit\n```\n\n* Short bullet with `scripts/lib.py`.\n'
    cases = {
        'PR number': clean + '\nFixed in PR #12.\n',
        'rule number': clean + '\nSee rule 25.\n',
        'board reference': clean + '\nTracked in pathway_forward.md.\n',
        'source line': clean + '\nSee meshgen.f90:612 for details.\n',
        'agent name': clean + '\nPorted by mira.\n',
        'sha': clean + '\nReference `51b7649`.\n',
        'two H1': clean + '\n# Second\n',
        'H4': clean + '\n#### Deep\n',
        'long bullet': clean + '\n* ' + 'x' * MAX_BULLET_CHARS + '\n',
        'file:line in inline code': clean + '\nSee `meshgen.f90:612`.\n',
        'md:line': clean + '\nSee README.md:40.\n',
        'notes ref': clean + '\nSee `docs/notes/NOTES_tpv30_gate.md`.\n',
        'bare sha': clean + '\nFixed in 51b7649 last week.\n',
        'gate internals': clean + '\nCompared at `CASE_BOUND`.\n',
    }
    legit = clean + ('\nChecksum md5 f6df5ececc9d95d0c6db43c9ba0c3c6e.\nWaveforms can be fetched from IRIS; the run used a Mira '
                     'allocation and the rupture pathway was bilateral. Run '
                     '`python3 testsys/run.py e2e` in `docs/user/`.\n')
    if problems('synthetic-legit.md', legit):
        fails.append('legitimate user content was flagged: %r' % problems('synthetic-legit.md', legit))
    fails = []
    if problems('synthetic-clean.md', clean):
        fails.append('clean synthetic doc was flagged: %r' % problems('synthetic-clean.md', clean))
    for label, text in cases.items():
        if not problems('synthetic.md', text):
            fails.append('negative control %r was NOT flagged -- the check cannot go red' % label)
    long_readme = '# T\n' + 'line\n' * README_MAX_LINES
    if not problems('README.md', long_readme):
        fails.append('negative control README length cap was NOT flagged')
    return fails


def main():
    fails = negative_controls()
    if fails:
        print('FAIL test_user_docs_style (checker)')
        for f in fails:
            print(' -', f)
        return 1
    files = user_docs()
    bad = []
    for f in files:
        with open(os.path.join(ROOT, f), errors='replace') as fh:
            bad += problems(f, fh.read())
    if bad:
        print('FAIL test_user_docs_style: %d problem(s) in %d user-facing file(s)' % (len(bad), len(files)))
        for b in bad:
            print(' -', b)
        return 1
    print('SUCCESS test_user_docs_style: %d user-facing file(s) clean; 15 negative '
          'controls flagged, clean and legitimate-content controls passed' % len(files))
    return 0


if __name__ == '__main__':
    sys.exit(main())
