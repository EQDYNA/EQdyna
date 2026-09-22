#! /usr/bin/env python3
"""
Regression guard + generator for docs/HISTORY.md: the legacy timeline
recovered from the commit messages that carry it (item 49, rules 1, 2, 3).

THE PROBLEM. This repo's first eight commits are version drops whose MESSAGES
carry real authorship dates spanning 2014-09-12 to 2016-10-05, but all eight
are STAMPED 2018-11-26 -- the bulk-import date. Every date-ordered tool (git
log, blame, GitHub's history view, any release-cadence analysis) therefore
sees a flat wall on one day and the true timeline is unreachable.

WHY A TABLE AND NOT A FIX. Restoring the real author dates IS a history
rewrite: it moves every downstream SHA -- 723 commits, all 20 tags, both
forks. Item 36 is decided, no history rewrite, and that stands. `git notes`
was the other cheap option and was not taken: notes live in a ref most clones
never fetch and most viewers never show, so the recovered timeline would be
as unreachable as the problem it fixes. A committed Markdown table is
readable by anyone with the repo, survives clone and fork, and -- unlike
notes -- can be checked by a test.

NOTHING IS INVENTED HERE. Every field is parsed out of the commit message it
belongs to. Two of the eight state only a month ('Date 09/2015'), so the
table carries a PRECISION column and says `month` rather than fabricating a
day. One states no version at all -- `221b3ad` mentions "Version 3.2.1" only
as the version it was MERGED INTO ("3D MPI (Bin Luo) incorporated with
Version 3.2.1") -- so the version parser accepts a `Version` token only where
it is declared in the header position, and that row reads `(unstated)`.

Run with --regenerate to rewrite docs/HISTORY.md from the messages; run with
no arguments (what testsys does) to verify the committed table still matches
what the messages say, and that the eight are still exactly the commits the
import flattened.

Cheap (rule 9): git log only, well under 1 s.
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DOC = os.path.join(ROOT, 'docs', 'HISTORY.md')

# The bulk import: every commit git stamps with this AUTHOR date.
IMPORT_AUTHOR_DATE = '2018-11-26'

BEGIN = '<!-- BEGIN GENERATED TABLE -- testsys/regression/test_history_table.py -->'
END = '<!-- END GENERATED TABLE -->'


def git(*args):
    r = subprocess.run(('git', '-C', ROOT) + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('git %s failed: %s' % (' '.join(args), r.stderr.strip()))
    return r.stdout


def import_commits():
    """The commits the 2018-11-26 bulk import flattened, oldest first."""
    out = git('log', '--reverse', '--date=short', '--format=%H\t%ad\t%cd')
    rows = []
    for line in out.splitlines():
        sha, ad, cd = line.split('\t')
        if ad == IMPORT_AUTHOR_DATE:
            rows.append((sha, ad, cd))
    return rows


def message(sha):
    return git('log', '-1', '--format=%B', sha)


def parse_authored_date(text):
    """('YYYY-MM-DD'|'YYYY-MM', 'day'|'month') from the message's own Date."""
    m = re.search(r'\bDate:?\s+(\d{2})/(\d{2})/(\d{4})\b', text)
    if m:
        return '%s-%s-%s' % (m.group(3), m.group(1), m.group(2)), 'day'
    m = re.search(r'\bDate:?\s+(\d{2})/(\d{4})\b', text)
    if m:
        return '%s-%s' % (m.group(2), m.group(1)), 'month'
    return None, None


def parse_author(text):
    m = re.search(r'\bAuthor:?\s+(.+?)\s+Date\b', text)
    return m.group(1).strip() if m else None


def parse_version(text):
    """The version this commit DECLARES, or None.

    Header position only: the token must follow the Email field (or open the
    message). `221b3ad` names "Version 3.2.1" mid-sentence as the version its
    MPI work was incorporated into; reading that as its own version would be
    inventing a number the commit never claims.
    """
    m = re.search(r'\bEmail:?\s+\S+@\S+', text)
    tail = text[m.end():] if m else text
    m = re.match(r'\s*Version\s+(\d+(?:\.\d+)*)', tail)
    return m.group(1) if m else None


def parse_headline(text):
    """An EXCERPT of the message's own verification line -- never paraphrased
    and never sourced from outside the commit. It stops at the first line
    break or full stop, so a citation split across lines (7d603a1's "Ma\\nand
    Liu (2006)") shows as a truncation; the commit message remains the full
    text and this column is labelled an excerpt for that reason."""
    for pat in (r'verified against ([^.\n]+)', r'tested against ([^.\n]+)',
                r'validated against ([^.\n]+)'):
        m = re.search(pat, text, re.I)
        if m:
            return 'verified against %s' % m.group(1).strip()
    return ''


def rows():
    out = []
    for sha, ad, cd in import_commits():
        text = message(sha)
        date, precision = parse_authored_date(text)
        out.append({
            'sha': sha, 'short': sha[:7], 'stamp_author': ad, 'stamp_commit': cd,
            'date': date, 'precision': precision,
            'author': parse_author(text), 'version': parse_version(text),
            'note': parse_headline(text),
        })
    return out


def render_table(data):
    lines = [BEGIN, '',
             '| # | commit | TRUE authored date | precision | author(s) | version declared | git stamp (author / committer) | the commit\'s own verification line (excerpt) |',
             '|---|--------|--------------------|-----------|-----------|------------------|-------------------------------|----------------------------------------------|']
    for i, r in enumerate(data, 1):
        lines.append('| %d | `%s` | **%s** | %s | %s | %s | %s / %s | %s |' % (
            i, r['short'], r['date'] or '(unstated)', r['precision'] or '--',
            r['author'] or '(unstated)',
            r['version'] or '(unstated)',
            r['stamp_author'], r['stamp_commit'],
            r['note'] or '--'))
    lines += ['', END]
    return '\n'.join(lines)


PREAMBLE = """# EQdyna legacy history, 2014-2016

The first %(n)d commits of this repository are version drops of code written
between **%(first)s and %(last)s**, bulk-imported on **%(stamp)s**. Git
stamps all %(n)d with the import date, so `git log`, `git blame`, GitHub's
history view and any date-ordered analysis see a flat wall on one day. The
real dates were never lost -- they are written in the commit MESSAGES -- but
no tool reads them.

This file is that timeline, parsed out of those messages.

**It is a record, not a repair.** Restoring the true author dates would
rewrite every downstream SHA (%(total)d commits, %(tags)d tags, both forks);
item 36 decided against a history rewrite and that stands. `git notes` was
the other cheap option and was rejected: notes live in a ref most clones
never fetch and most viewers never show.

**Nothing below is inferred.** Every cell is parsed from the commit it
describes by `testsys/regression/test_history_table.py`, which also VERIFIES
this table on every regression run -- if a message and this table ever
disagree, the test fails. Two commits state only a month, so the precision
column says `month` rather than inventing a day. One declares no version at
all (`221b3ad` mentions 3.2.1 only as the version its MPI work was merged
INTO), so its version reads `(unstated)`.

One detail the earlier prose account got slightly wrong, corrected here from
the log: the eight share an AUTHOR date of %(stamp)s, but `d1dc019`'s
COMMITTER date is 2018-11-27 -- the import ran past midnight. This table keys
on the author date, and shows both stamps.

Regenerate with `python3 testsys/regression/test_history_table.py --regenerate`.

## The pre-2018 source itself

A legacy archive of the pre-2018 source exists on disk at
`scratch/legacy_import/Examples` (gitignored like every `scratch/` asset),
hash-manifested at `docs/evidence/legacy_import_archive.sha256`. It is NOT
needed to read this table: the dates come from the commit messages, which are
in the repo. The archive's own bundled `.git` covers only 4.1.1/4.1.2/4.1.3
and stamps all three 2018-11-26 -- a second snapshot import, not a timeline,
and importing it would add a second flat wall rather than fix the first.

## The timeline

"""


def build_doc(data):
    tags = len(git('tag').split())
    total = int(git('rev-list', '--count', 'HEAD').strip())
    dated = [r['date'] for r in data if r['date']]
    head = PREAMBLE % {'n': len(data), 'first': min(dated), 'last': max(dated),
                       'stamp': IMPORT_AUTHOR_DATE, 'total': total, 'tags': tags}
    return head + render_table(data) + '\n'


def main(argv):
    regen = '--regenerate' in argv
    data = rows()
    problems = []

    if len(data) != 8:
        problems.append('expected 8 bulk-imported commits stamped %s, found %d '
                        '-- history changed shape' % (IMPORT_AUTHOR_DATE, len(data)))
    for r in data:
        if not r['date']:
            problems.append('%s: no Date in its message -- cannot be tabulated '
                            'without inventing one' % r['short'])
        if not r['author']:
            problems.append('%s: no Author in its message' % r['short'])

    # The table is only about commits whose stamp is WRONG. If a stated date
    # ever agrees with the stamp, this file has drifted from its own purpose.
    for r in data:
        if r['date'] and r['date'].startswith(IMPORT_AUTHOR_DATE[:4]):
            problems.append('%s: states %s, same year as the import stamp -- '
                            'it may not belong in this table' % (r['short'], r['date']))

    if regen:
        os.makedirs(os.path.dirname(DOC), exist_ok=True)
        open(DOC, 'w').write(build_doc(data))
        print('regenerated %s (%d rows)' % (DOC, len(data)))
        return 1 if problems else 0

    print('Regression guard: docs/HISTORY.md matches the commit messages')
    if not os.path.exists(DOC):
        print('\nFAIL: %s does not exist -- run with --regenerate' % DOC)
        return 1
    have = open(DOC, errors='replace').read()
    want = build_doc(data)
    if have != want:
        problems.append('docs/HISTORY.md does not match what the commit '
                        'messages say -- regenerate it (--regenerate) and '
                        'review the diff; do NOT hand-edit the generated table')
        h = have.split(BEGIN)[-1].split(END)[0].strip().splitlines()
        w = want.split(BEGIN)[-1].split(END)[0].strip().splitlines()
        for i in range(max(len(h), len(w))):
            a = h[i] if i < len(h) else '(missing)'
            b = w[i] if i < len(w) else '(missing)'
            if a != b:
                problems.append('  line %d:\n    committed: %s\n    messages:  %s'
                                % (i + 1, a, b))

    for r in data:
        print('  %s  %-10s %-5s  %-22s v%-9s %s'
              % (r['short'], r['date'], r['precision'], r['author'],
                 r['version'] or '(unstated)', r['note'] or ''))

    if problems:
        print('\nFAIL: %d problem(s)' % len(problems))
        for p in problems:
            print('  ' + p)
        return 1
    print('\nPASS: %d legacy commits, dates %s..%s, table matches the messages'
          % (len(data), min(r['date'] for r in data), max(r['date'] for r in data)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
