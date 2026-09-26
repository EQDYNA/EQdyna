#!/usr/bin/env python3
"""Rule 27's release-cadence check (pathway_forward.md item 133).

A release is DUE when a physics-or-output change has landed since the last
tag AND either 7 days have passed since that tag or 10 PRs have merged since
it. This prints exactly one line and ALWAYS exits 0. It is an advisory board
Command (rule 14), never a regression-tier check. The file is named check_*,
not test_*, so `testsys/run.py regression` never picks it up.

A "physics or output change" is any of (rule 27):
  1. a file under src/fortran/ or src/python/eqdyna/ whose diff since the tag
     is not entirely comment or blank lines. A diff this cannot classify
     counts IN (rule 2: fail toward "yes"). In Python only `#` lines and
     blank lines count as comments; docstring text is not separated out, so
     a docstring-only change counts IN, the conservative side.
  2. any change to an output writer: library_output.f90 / library_output.py,
     or scripts/plotRuptureDynamics (writes fault.dyna.r.nc);
  3. any change under test.reference.results/.

Usage:  python3 testsys/regression/check_release_due.py
"""
import datetime
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DAYS_DUE = 7
PRS_DUE = 10
CODE_DIRS = ('src/fortran/', 'src/python/eqdyna/')
OUTPUT_FILES = ('src/fortran/library_output.f90',
                'src/python/eqdyna/library_output.py',
                'scripts/plotRuptureDynamics')
REFERENCE_DIR = 'test.reference.results/'


def git(*args):
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True,
                          text=True, check=True).stdout


def is_comment_or_blank(line, path):
    """True only for a changed line we can PROVE is non-code."""
    s = line.strip()
    if not s:
        return True
    if path.endswith(('.f90', '.F90', '.f', '.F')):
        return s.startswith('!')
    if path.endswith('.py'):
        return s.startswith('#')
    return False                       # unknown language: count it IN


def changed_lines(diff_text):
    """The +/- body lines of a unified diff, headers excluded."""
    out = []
    for l in diff_text.splitlines():
        if l.startswith(('+++', '---', '@@', 'diff ', 'index ', 'new file',
                         'deleted file', 'similarity', 'rename ', 'Binary')):
            if l.startswith('Binary'):
                out.append(l)          # binary change: cannot classify -> IN
            continue
        if l.startswith(('+', '-')):
            out.append(l[1:])
    return out


def classify(path, diff_text):
    """Why `path` is a physics/output change, or None if it is not one."""
    if path.startswith(REFERENCE_DIR):
        return 'reference changed (%s)' % path
    if path in OUTPUT_FILES:
        return 'output writer changed (%s)' % path
    if path.startswith(CODE_DIRS):
        lines = changed_lines(diff_text)
        if not lines:
            return 'code changed, diff unclassifiable (%s)' % path
        if all(is_comment_or_blank(l, path) for l in lines):
            return None
        return 'code changed (%s)' % path
    return None


def main():
    try:
        tag = git('describe', '--tags', '--abbrev=0').strip()
    except subprocess.CalledProcessError:
        print('RELEASE DUE: no tag found (cannot measure; failing toward due)')
        return 0
    tag_date = datetime.datetime.fromisoformat(
        git('log', '-1', '--format=%cI', tag).strip())
    days = (datetime.datetime.now(datetime.timezone.utc) - tag_date).days
    subjects = git('log', '--format=%s', '%s..HEAD' % tag)
    # squash subjects end '(#NN)'; a real merge commit says 'Merge pull request #NN'
    prs = len(set(re.findall(r'\(#(\d+)\)', subjects))
              | set(re.findall(r'Merge pull request #(\d+)', subjects)))

    reasons = []
    for path in [p for p in git('diff', '--name-only', '%s..HEAD' % tag).splitlines() if p]:
        why = classify(path, git('diff', '-U0', '%s..HEAD' % tag, '--', path))
        if why:
            reasons.append(why)

    counts = '%d PRs / %d days since %s' % (prs, days, tag)
    if reasons and (days >= DAYS_DUE or prs >= PRS_DUE):
        trip = 'days' if days >= DAYS_DUE else 'PRs'
        print('RELEASE DUE: physics/output change since tag (%s; %d file(s)) '
              'and the %s threshold tripped, %s'
              % (reasons[0], len(reasons), trip, counts))
    else:
        print('release not due: %s, physics/output change since tag: %s'
              % (counts, 'yes (%d file(s))' % len(reasons) if reasons else 'no'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
