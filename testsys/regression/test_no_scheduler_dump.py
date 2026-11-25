#! /usr/bin/env python3
"""
Regression guard (rule 28): no scheduler accounting dump and no literal
compute-allocation identifier in committed evidence, docs, board or logs.

INCIDENT CLASS: sacct/squeue/sinfo/scontrol output carries the account
(allocation) name of whoever ran it. Pasted into an evidence README, a log, a
board row or a PR it publishes that identifier in a PUBLIC repo. Record only
the wall-time/elapsed/state fields, extracted by an explicit field list.

WHAT THIS CHECKS, over every tracked text file outside CODE_PREFIXES, with no
allowlist: (1) sacct/scontrol column headers or key=value fields that name the
account (`Account` next to another column header, `AllocTRES`, `AllocCPUS`,
`Account=<value>`); (2) a literal value after `--account` or `sbatch/srun/
salloc/sacct ... -A` or `#SBATCH -A`. A value that starts with `<`, `$`, `%`,
`{` or a quote is a placeholder and is not a value.

SCOPE LIMIT, stated not hidden: code, tests and batch templates
(src/, scripts/, testsys/, .github/) legitimately WRITE `-A <variable>` and
their fixtures use fake names; a literal there is a code finding reported by
`--all`, routed to its owner, not gated here. Commit messages and PR text are
not files: that part of rule 28 is a norm.

Fragments are split so this file never matches its own patterns.
Both outcomes are exercised below by synthetic text (rule 14a).
"""
import re
import subprocess
import sys

ROOT = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True,
                      text=True, check=True).stdout.strip()
SELF = 'testsys/regression/test_no_scheduler_dump.py'
CODE_PREFIXES = ('src/', 'scripts/', 'testsys/', '.github/')
MAX_BYTES = 5 * 1024 * 1024
_VAL = r'[A-Za-z0-9]'          # a real value, not <x>, $x, %s, {x}, "x"
PATTERNS = [
    ('account column header', re.compile(r'\bAcc' + r'ount\s+(Cluster|User|State|Partition|Elapsed|JobName|AllocCPUS)\b')),
    ('allocation field', re.compile(r'\bAlloc' + r'(TRES|CPUS)\b')),
    ('account key=value', re.compile(r'\bAcc' + r'ount=' + _VAL)),
    ('--account value', re.compile(r'--acc' + r'ount[= ]\s*' + _VAL)),
    ('-A value', re.compile(r'(?:sbatch|srun|salloc|sacct|#SBATCH)\b[^\n]*\s-A[ =]?\s*' + _VAL)),
]


def scan_text(text):
    """Return [(lineno, label, line)] for every pattern hit in text."""
    hits = []
    for n, line in enumerate(text.splitlines(), 1):
        for label, rx in PATTERNS:
            if rx.search(line):
                hits.append((n, label, line.strip()[:100]))
    return hits


def tracked_hits(include_code):
    out = subprocess.run(['git', 'ls-files', '-z'], cwd=ROOT, capture_output=True,
                         check=True).stdout.decode().split('\0')
    hits = []
    for rel in filter(None, out):
        if rel == SELF or (not include_code and rel.startswith(CODE_PREFIXES)):
            continue
        try:
            with open(ROOT + '/' + rel, 'rb') as f:
                raw = f.read(MAX_BYTES + 1)
        except OSError:
            continue
        if len(raw) > MAX_BYTES or b'\0' in raw:
            continue
        for n, label, line in scan_text(raw.decode('utf-8', 'replace')):
            hits.append((rel, n, label, line))
    return hits


def main():
    # Both outcomes, synthetic (rule 14a): a dump must go red, placeholders green.
    red = ('JobID Acc' + 'ount Cluster State\n12 abc123 x COMPLETED\n'
           'srun --acc' + 'ount=PROJ1\nsbatch -A PROJ1 job.sh\n#SBATCH -A PROJ1\n')
    green = ('sbatch -A <alloc> job.sh\n--acc' + 'ount $ACCT\n#SBATCH -A "%s"\n'
             'sacct -o JobID,Elapsed,State\ngrep -A 3 foo\n')
    nred, ngreen = len(scan_text(red)), len(scan_text(green))
    assert nred == 4 and ngreen == 0, 'self-test: red=%d (want 4) green=%d (want 0)' % (nred, ngreen)
    print('  PASS  self-test: dump text -> %d hits (red), placeholder text -> %d (green)' % (nred, ngreen))
    if '--all' in sys.argv:
        for h in tracked_hits(True):
            print('  %s:%d [%s] %s' % h)
        return 0
    hits = tracked_hits(False)
    if hits:
        for h in hits:
            print('  FAIL  %s:%d [%s] %s' % h)
        print('FAIL test_no_scheduler_dump: %d hit(s); rule 28 forbids scheduler dumps and literal '
              'allocation identifiers in committed text. Report, do not paste the value.' % len(hits))
        return 1
    # Rule-index sync: the stated heading count must equal the real one.
    with open(ROOT + '/PROJECT_RULES.md') as f:
        rules = f.read()
    total = len(re.findall(r'^## ', rules, re.M))
    assert re.search(r'\b%d `## ` headings total' % total, rules), \
        'PROJECT_RULES.md index states a heading count other than the real %d' % total
    assert re.search(r'^28\. ', rules, re.M), 'rule 28 missing from the index'
    print('  PASS  0 scheduler-dump / allocation-value hits in non-code tracked files; index count %d ok' % total)
    print('SUCCESS test_no_scheduler_dump')
    return 0


if __name__ == '__main__':
    sys.exit(main())
