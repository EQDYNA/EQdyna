#!/usr/bin/env python3
"""Rule 27's release-cadence check (pathway_forward.md item 133).

Owner, 2026-09-25: "release as soon as a physics or output change lands, or
at the latest after a week or ~10 PRs, whichever comes first." So a release
is DUE when EITHER
  (a) a physics-or-output change has landed since the last tag -- at once,
      no threshold; OR
  (b) at least one PR has merged since the tag AND 7 days have passed since
      it or 10 PRs have merged since it.
A docs/board-only stretch (no PR) never forces a release. This prints
exactly one line and ALWAYS exits 0. It is an advisory board Command
(rule 14), never a regression-tier check. The file is named check_*, not
test_*, so `testsys/run.py regression` never picks it up.

A "physics or output change" is, per commit since the tag, any of (rule 27):
  1. a file under src/fortran/ or src/python/eqdyna/ whose diff in that
     commit is not entirely comment or blank lines. A diff this cannot
     classify counts IN (rule 2: fail toward "yes"). In Python only `#`
     lines and blank lines count as comments; docstring text is not
     separated out, so a docstring-only change counts IN, the conservative
     side.
  2. any change to an output writer: library_output.f90 / library_output.py,
     or scripts/plotRuptureDynamics (writes fault.dyna.r.nc);
  3. any change under test.reference.results/.
EXCEPT a commit listed in docs/release_exempt.txt with its evidence that
the change is output-neutral (bit-identical gates -- a refactor, a
retirement, a perf change). An exemption is a written claim, one line
`<sha> <evidence>`; a line without evidence, or a sha that does not name a
commit since the tag, exempts nothing (fail toward "yes"). An exempt commit
still counts toward (b)'s PR count.

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
EXEMPT_PATH = os.path.join(ROOT, 'docs', 'release_exempt.txt')


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


def load_exemptions(path):
    """{sha_prefix: evidence} from `<sha> <evidence>` lines; '#' comments and
    blanks ignored. A line with no evidence, or a non-hex/short sha, is not
    an exemption -- dropped here, so the commit counts IN."""
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        sha, _, evidence = line.partition(' ')
        if re.fullmatch(r'[0-9a-f]{7,40}', sha) and evidence.strip():
            out[sha] = evidence.strip()
    return out


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

    exempt = load_exemptions(EXEMPT_PATH)
    reasons, n_exempt = [], 0
    for sha in git('log', '--format=%H', '%s..HEAD' % tag).split():
        if any(sha.startswith(e) for e in exempt):
            n_exempt += 1
            continue
        files = git('diff-tree', '--no-commit-id', '--name-only', '-r', '-m',
                    '--first-parent', sha).split('\n')
        for path in sorted(set(p for p in files if p)):
            why = classify(path, git('show', '-U0', '--format=',
                                     '--first-parent', '-m', sha, '--', path))
            if why:
                reasons.append('%s in %s' % (why, sha[:7]))

    counts = '%d PRs / %d days since %s' % (prs, days, tag)
    ex = ' (%d exempt commit(s), docs/release_exempt.txt)' % n_exempt if n_exempt else ''
    if reasons:
        print('RELEASE DUE: physics/output change since tag (%s; %d file(s)), '
              '%s%s' % (reasons[0], len(reasons), counts, ex))
    elif prs >= 1 and (days >= DAYS_DUE or prs >= PRS_DUE):
        trip = 'days' if days >= DAYS_DUE else 'PRs'
        print('RELEASE DUE: the %s threshold tripped with no physics/output '
              'change, %s%s' % (trip, counts, ex))
    else:
        print('release not due: %s, physics/output change since tag: no%s'
              % (counts, ex))
    return 0


if __name__ == '__main__':
    sys.exit(main())
