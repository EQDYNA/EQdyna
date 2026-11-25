#! /usr/bin/env python3
"""
Regression guard (rule 28): no scheduler accounting dump and no literal
compute-allocation identifier in committed evidence, docs, board or logs.

INCIDENT CLASS: sacct/squeue/sinfo/scontrol output carries the account
(allocation) name of whoever ran it. Pasted into an evidence README, a log, a
board row or a PR it publishes that identifier in a PUBLIC repo. Record only
the wall-time/elapsed/state fields, extracted by an explicit field list.

WHAT THIS CHECKS, over every tracked file outside CODE_PREFIXES, with no
allowlist: (1) sacct/scontrol column headers or key=value fields that name the
account (`Account` next to another column header, an `Account` cell in a
pipe/markdown table, `sacct -o ...account...`, `AllocTRES`, `AllocCPUS`,
`Account=<value>`, a bare `account: <value>` line); (2) a literal value in
`SLURM_JOB_ACCOUNT=`/`SBATCH_ACCOUNT=`, after `--account`, or after
`sbatch/srun/salloc/sacct ... -A` / `#SBATCH -A`. Matching is case-insensitive
for the column/field forms. A value that starts with `<`, `$`, `%`, `{` or a
quote is a placeholder and is not a value.

SIZE/BINARY: every tracked file is scanned, including files over 5 MB and
NUL-bearing ones (decoded with errors='replace'); a file that cannot be opened
fails the test, never a silent skip. The hit line text is never printed (it may
be the value); only file:line and the pattern label.

SCOPE LIMIT, stated not hidden: code, tests and batch templates
(src/, scripts/, testsys/, .github/) legitimately WRITE `-A <variable>` and
their fixtures use fake names; a literal there is a code finding reported by
`--all`, routed to its owner, not gated here. Commit messages and PR text are
not files: that part of rule 28 is a norm.

INDEX SYNC: the stated "N numbered rules" and "<word> lettered sub-rules", the
"M `## ` headings total" count, and an index line for EVERY numbered and
lettered heading are checked against PROJECT_RULES.md's real headings.

Fragments are split so this file never matches its own patterns.
Every outcome is exercised below by synthetic text (rule 14a): each pattern
fires individually, placeholder text stays green, a doctored index goes red.
"""
import re
import subprocess
import sys

ROOT = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True,
                      text=True, check=True).stdout.strip()
SELF = 'testsys/regression/test_no_scheduler_dump.py'
CODE_PREFIXES = ('src/', 'scripts/', 'testsys/', '.github/')
_VAL = r'[A-Za-z0-9]'          # a real value, not <x>, $x, %s, {x}, "x"
_I = re.I
PATTERNS = [
    ('account column header', re.compile(r'\bAcc' + r'ount\s+(Cluster|User|State|Partition|Elapsed|JobName|AllocCPUS)\b', _I)),
    ('account column cell', re.compile(r'(?:^|\|)\s*Acc' + r'ount\s*(?:\||$)', _I)),
    ('sacct format account', re.compile(r'\bsacct\b[^\n]*(?:-o|--format)[ =]\S*\bAcc' + r'ount\b', _I)),
    ('allocation field', re.compile(r'\bAlloc' + r'(TRES|CPUS)\b', _I)),
    ('account key=value', re.compile(r'\bAcc' + r'ount=' + _VAL, _I)),
    ('account yaml line', re.compile(r'^\s*acc' + r'ount:\s*' + _VAL + r'[\w.-]*\s*$', _I)),
    ('env account', re.compile(r'\b(?:SLURM_JOB_ACC' + r'OUNT|SBATCH_ACC' + r'OUNT)=\s*[\'"]?' + _VAL)),
    ('--account value', re.compile(r'--acc' + r'ount[= ]\s*' + _VAL)),
    ('-A value', re.compile(r'(?:sbatch|srun|salloc|sacct|#SBATCH)\b[^\n]*\s-A[ =]?\s*' + _VAL)),
]
# One red sample per pattern label; each must fire that label itself.
RED = {
    'account column header': 'JobID Acc' + 'ount Cluster State',
    'account column cell': 'JobID|Acc' + 'ount|State',
    'sacct format account': 'sacct -P -o JobID,acc' + 'ount,Elapsed',
    'allocation field': 'Alloc' + 'TRES=cpu=4',
    'account key=value': 'Acc' + 'ount=PROJ1',
    'account yaml line': 'acc' + 'ount: PROJ1',
    'env account': 'SLURM_JOB_ACC' + 'OUNT=PROJ1',
    '--account value': 'srun --acc' + 'ount=PROJ1',
    '-A value': 'sbatch -A PROJ1 job.sh',
}
# Further spellings that must also be red (any pattern).
RED_ALSO = ['SBATCH_ACC' + 'OUNT=PROJ1', '| JobID | Acc' + 'ount | State |', '#SBATCH -A PROJ1',
            'sacct --format=JobID,Acc' + 'ount']
GREEN = ('sbatch -A <alloc> job.sh\n--acc' + 'ount $ACCT\n#SBATCH -A "%s"\n'
         'sacct -o JobID,Elapsed,State\ngrep -A 3 foo\nSLURM_JOB_ACC' + 'OUNT=$X\n'
         'SBATCH_ACC' + 'OUNT=<alloc>\nacc' + 'ount: <alloc>\nAcc' + 'ount=%s\n'
         'an acc' + 'ount of the problem\n')

_WORDS = {w: i for i, w in enumerate(
    'zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen '
    'fifteen sixteen seventeen eighteen nineteen'.split())}
_WORDS.update({'twenty': 20, 'thirty': 30, 'forty': 40, 'fifty': 50, 'sixty': 60,
               'seventy': 70, 'eighty': 80, 'ninety': 90})


def scan_text(text):
    """Return [(lineno, label, line)] for every pattern hit in text."""
    hits = []
    for n, line in enumerate(text.splitlines(), 1):
        for label, rx in PATTERNS:
            if rx.search(line):
                hits.append((n, label, line.strip()[:100]))
    return hits


def tracked_hits(include_code):
    """Return (hits, errors). Every file is scanned binary-safe; an unreadable file is an
    error, never a skip."""
    out = subprocess.run(['git', 'ls-files', '-z'], cwd=ROOT, capture_output=True,
                         check=True).stdout.decode().split('\0')
    hits, errors = [], []
    for rel in filter(None, out):
        if rel == SELF or (not include_code and rel.startswith(CODE_PREFIXES)):
            continue
        try:
            with open(ROOT + '/' + rel, 'rb') as f:
                raw = f.read()
        except OSError as e:
            errors.append((rel, type(e).__name__))
            continue
        for n, label, _line in scan_text(raw.decode('utf-8', 'replace')):
            hits.append((rel, n, label))
    return hits, errors


def check_index(rules):
    """Return a list of problems: stated counts vs real heading counts, and index lines."""
    problems = []
    heads = re.findall(r'^## (\d+)([a-z]?)\. ', rules, re.M)
    nnum = sum(1 for _, let in heads if not let)
    nlet = sum(1 for _, let in heads if let)
    ntot = len(re.findall(r'^## ', rules, re.M))
    if nnum + nlet != ntot:
        problems.append('%d `## ` headings are not of the form "N." or "Na."' % (ntot - nnum - nlet))
    m = re.search(r'\b(\d+) numbered rules', rules)
    if not m or int(m.group(1)) != nnum:
        problems.append('index states %s numbered rules; real numbered headings %d'
                        % (m.group(1) if m else 'none', nnum))
    m = re.search(r'\b([a-z]+(?:-[a-z]+)?) lettered sub-rules', rules)
    try:
        stated = sum(_WORDS[p] for p in m.group(1).split('-')) if m else None
    except KeyError:
        stated = None
    if stated != nlet:
        problems.append('index states %s lettered sub-rules; real lettered headings %d'
                        % (stated if stated is not None else 'unparseable', nlet))
    if not re.search(r'\b%d `## ` headings total' % ntot, rules):
        problems.append('index does not state "%d `## ` headings total"' % ntot)
    if not re.search(r'reads %d\b' % ntot, rules):
        problems.append('index verify line does not read %d' % ntot)
    first = re.search(r'^## 1\. ', rules, re.M)
    index = rules[:first.start()] if first else rules
    for num, let in heads:
        if not re.search(r'^%s%s\. ' % (num, let), index, re.M):
            problems.append('rule %s%s has no index line' % (num, let))
    return problems


def main():
    # Rule 14a: every pattern fires on its own sample; placeholders stay green.
    for label, _ in PATTERNS:
        assert label in RED, 'self-test: no red sample for pattern %r' % label
    for label, sample in RED.items():
        got = {lab for _, lab, _ in scan_text(sample)}
        assert label in got, 'self-test: pattern %r did not fire on its own sample' % label
    for sample in RED_ALSO:
        assert scan_text(sample), 'self-test: extra red sample did not fire: %r' % sample[:24]
    ngreen = len(scan_text(GREEN))
    assert ngreen == 0, 'self-test: placeholder text produced %d hit(s) (want 0)' % ngreen
    with open(ROOT + '/PROJECT_RULES.md', encoding='utf-8') as f:
        book = f.read()
    assert check_index(book.replace('28 numbered rules', '27 numbered rules', 1)), \
        'self-test: index check stayed green on a wrong numbered count'
    nolet = re.sub(r'^28\. .*\n', '', book, count=1, flags=re.M)
    assert any('no index line' in x for x in check_index(nolet)), \
        'self-test: a dropped index line was not caught'
    print('  PASS  self-test: %d patterns each fire on their own sample (red), placeholder text -> 0 (green), '
          'doctored index -> red' % len(PATTERNS))
    if '--all' in sys.argv:
        hits, errors = tracked_hits(True)
        for h in hits:
            print('  %s:%d [%s]' % h)
        for e in errors:
            print('  UNREADABLE %s (%s)' % e)
        return 0
    hits, errors = tracked_hits(False)
    if hits or errors:
        for h in hits:
            print('  FAIL  %s:%d [%s]' % h)
        for e in errors:
            print('  FAIL  unreadable %s (%s)' % e)
        print('FAIL test_no_scheduler_dump: %d hit(s), %d unreadable file(s); rule 28 forbids scheduler '
              'dumps and literal allocation identifiers in committed text. Report, do not paste the '
              'value (line text is deliberately not printed).' % (len(hits), len(errors)))
        return 1
    problems = check_index(book)
    if problems:
        for pr in problems:
            print('  FAIL  PROJECT_RULES.md: ' + pr)
        print('FAIL test_no_scheduler_dump: rule-book index out of sync (%d problem(s))' % len(problems))
        return 1
    print('  PASS  0 scheduler-dump / allocation-value hits in non-code tracked files (all files scanned, '
          'binary-safe); index counts and lines ok')
    print('SUCCESS test_no_scheduler_dump')
    return 0


if __name__ == '__main__':
    sys.exit(main())
