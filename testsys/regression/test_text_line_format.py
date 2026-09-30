#! /usr/bin/env python3
"""
Regression guard: no `src/fortran` list-directed write of a string literal is
longer than 80 columns. Such text is written with an explicit '(1X,A)'.

Incident (2026-09-29, first ls6 sweep): `write(51,*) '# Time series in 8
columns; ...'` is 83 columns. gfortran writes it whole; ifort ends a
list-directed record at 80, so the station header split in two and the orphan
tail (` 6.7`, `E16.7`) carried no leading '#'. The station reader took it for a
data row, and all 11 Fortran e2e cells failed parsing their own output -- a
defect no gfortran-only CI could see. The same wrap cut the FATAL reason
mid-word (`...or C_elas` / `tic=0.`).

This is a static check: it needs no Fortran compiler, so it runs on every
machine, including the ones whose compiler does not wrap. What it cannot see
is a long CONCATENATED or variable-carrying list-directed write; those are the
cases '(1X,A)' fixes, and the fix for a hit is the same.

Both ways (rule 14a): the checker is first run on a synthetic 83-column write,
which it must flag, and on the '(1X,A)' form of the same text, which it
must not. A checker that cannot go red fails first.
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')

LIMIT = 80  # ifort's list-directed record length

# write(<unit>,*) '<string literal>'   -- and nothing else on the line
LIST_DIRECTED_LITERAL = re.compile(
    r"^\s*write\s*\(\s*[^,()]+\s*,\s*\*\s*\)\s*('(?:[^']|'')*')\s*(?:!.*)?$", re.I)


def offenders(lines):
    """[(lineno, columns)] for list-directed literal writes that ifort would
    wrap. The record is a leading blank plus the literal, un-doubled."""
    out = []
    for n, line in enumerate(lines, 1):
        m = LIST_DIRECTED_LITERAL.match(line)
        if m:
            cols = 1 + len(m.group(1)[1:-1].replace("''", "'"))
            if cols > LIMIT:
                out.append((n, cols))
    return out


def main():
    long_text = '# ' + 'x' * 81
    bad = offenders(["  write(51,*) '%s'" % long_text])
    good = offenders(["  write(51,'(1X,A)') '%s'" % long_text,
                      "  write(51,*) '# short'"])
    if bad != [(1, 84)] or good:
        print('FAIL test_text_line_format: negative control wrong -- flagged %r for '
              'the 83-col list-directed write (want [(1, 84)]) and %r for the '
              'the (1X,A) / short forms (want none); the checker cannot go red'
              % (bad, good))
        return 1

    files = sorted(glob.glob(os.path.join(FSRC, '*.f90')))
    if not files:
        print('FAIL test_text_line_format: no .f90 under %s' % FSRC)
        return 1
    fails, nwrites = [], 0
    for path in files:
        with open(path, errors='replace') as f:
            lines = f.read().split('\n')
        nwrites += sum(1 for l in lines if LIST_DIRECTED_LITERAL.match(l))
        for n, cols in offenders(lines):
            fails.append('%s:%d is %d columns as a list-directed write (ifort wraps '
                         'at %d) -- write it with (1X,A)'
                         % (os.path.relpath(path, ROOT), n, cols, LIMIT))
    if fails:
        print('FAIL test_text_line_format (%d)' % len(fails))
        for f in fails:
            print('  - ' + f)
        return 1
    print('PASS test_text_line_format: %d .f90 file(s), %d list-directed literal '
          'write(s) inspected, none over %d columns; negative control flags an '
          '83-col write' % (len(files), nwrites, LIMIT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
