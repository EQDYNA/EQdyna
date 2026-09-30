#!/usr/bin/env python3
"""Rule 27's release-cadence check (pathway_forward.md item 133).

Owner, 2026-09-25: "release as soon as a physics or output change lands, or
at the latest after a week or ~10 PRs, whichever comes first"; the PR
threshold lowered to 5 by the owner on 2026-09-28 ("I think 5 PR is good
enough to warrant a release" -> "Yes"). So a release
is DUE when EITHER
  (a) a physics-or-output change has landed since the last tag -- at once,
      no threshold; OR
  (b) at least one NON-INTERNAL PR has merged since the tag AND 7 days have
      passed since it or 5 NON-INTERNAL PRs have merged since it.
A docs/board-only stretch (no PR) never forces a release. This prints
exactly one line and ALWAYS exits 0. It is an advisory board Command
(rule 14), never a regression-tier check. The file is named check_*, not
test_*, so `testsys/run.py regression` never picks it up.

A "physics or output change" is testsys/change_class.is_output_change (item
1/6, 2026-09-30 classifier work): per commit since the tag, any of (rule 27,
UNCHANGED from before the classifier existed) --
  1. a file under src/fortran/ or src/python/eqdyna/ whose diff in that
     commit is not entirely comment or blank lines. A diff this cannot
     classify counts IN (rule 2: fail toward "yes"). In Python only `#`
     lines and blank lines count as comments; docstring text is not
     separated out, so a docstring-only change counts IN, the conservative
     side. An OpenMP sentinel (`!$`/`!$OMP`) counts as CODE even though it
     starts with a Fortran comment character -- it can add or remove
     parallelism (review item 7).
  2. any change to an output writer: library_output.f90 / library_output.py,
     or scripts/plotRuptureDynamics (writes fault.dyna.r.nc);
  3. any change under test.reference.results/.
This is DELIBERATELY NARROWER than "needs a sweep": editing
testsys/matrix.py, testsys/compare.py, or anything under testsys/e2e/
classifies PHYSICS (change_class.classify_path) -- it needs a fresh sweep
to prove the gate still passes -- but it is not itself criterion 1-3 above,
so it does not make a release due at once (review item 6: the gate's own
comparison machinery changing is not the same claim as "the solver's output
changed").
EXCEPT a commit listed in docs/release_exempt.txt with its evidence that
the change is output-neutral (bit-identical gates -- a refactor, a
retirement, a perf change). An exemption is a written claim, one line
`<sha> <evidence>`; a line without evidence, or a sha that does not name a
commit since the tag, exempts nothing (fail toward "yes"). An exemption
waives DUE-AT-ONCE only (review item 7) -- it never waives the release
sweep carry-forward (change_class.is_release_physics_path), which does not
consult this file at all, NOR (b)'s PR count: an
exempted src/ commit still counts (see below).

(b)'s PR COUNT excludes ONLY INTERNAL-only PRs (owner's own framing:
"INTERNAL-only PRs... never trigger; they ride along" -- NOT narrowed to
USER_FACING specifically, corrected from an earlier draft of this file that
did narrow it that way): a PR whose EVERY changed path classifies INTERNAL
(docs, board, testsys/unit, testsys/regression, ...) never counts, at all,
toward either the day or the PR threshold. A PR touching even one PHYSICS
or USER_FACING path counts -- INCLUDING a PHYSICS-only PR that is exempted
from the at-once trigger above (a src/ refactor with proven bit-identical
output, say): the exemption says "this did not change output", it does not
say "this PR did not happen". Squash-merge (rule 25) is this repo's ONLY
merge strategy, so one commit == one PR; a PR number is read from that
commit's OWN subject.

Usage:  python3 testsys/regression/check_release_due.py
"""
import datetime
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import change_class  # noqa: E402

DAYS_DUE = 7
PRS_DUE = 5
# CODE_DIRS/OUTPUT_FILES/REFERENCE_DIR now LIVE in change_class.py (item 1,
# 2026-09-30: "fewer full sweeps and fewer releases" -- one shared
# classifier, no second copy). Re-exported here under their original names
# so nothing else that imports them from this module breaks.
CODE_DIRS = change_class.CODE_DIRS
OUTPUT_FILES = change_class.OUTPUT_FILES
REFERENCE_DIR = change_class.REFERENCE_DIR
DEFAULT_EXEMPT_PATH = os.path.join(ROOT, 'docs', 'release_exempt.txt')
EXEMPT_PATH = DEFAULT_EXEMPT_PATH


def git(*args):
    # errors='replace': a blob with a Latin-1 comment and no NUL byte is not
    # "Binary" to git but is not UTF-8 either; decoding it strictly crashed
    # the one-line/exit-0 contract (PR #49 audit, measured).
    return subprocess.run(['git', *args], cwd=ROOT, capture_output=True,
                          text=True, errors='replace', check=True).stdout


def watched(path):
    """True if classify() could return a reason for `path` -- so only these
    paths are diffed per commit (a per-file, per-commit `git show -U0` is
    not cheap; rule 9)."""
    return (path.startswith(REFERENCE_DIR) or path in OUTPUT_FILES
            or path.startswith(CODE_DIRS))


def classify(path, diff_text):
    """Why `path` is a physics/output change (rule 27's DUE-at-once
    criteria 1-3), or None if it is not one. Thin wrapper around
    change_class.is_output_change -- this function's own job is only the
    human-readable REASON text; the yes/no logic lives in ONE place."""
    if not change_class.is_output_change(path, diff_text):
        return None
    if path.startswith(REFERENCE_DIR):
        return 'reference changed (%s)' % path
    if path in OUTPUT_FILES:
        return 'output writer changed (%s)' % path
    if path.startswith(CODE_DIRS):
        return 'code changed (%s)' % path
    return 'output change (%s)' % path  # defensive; watched() names no other case


def load_exemptions(text):
    """({sha_prefix: evidence}, n_ignored) from `<sha> <evidence>` lines;
    '#' comments and blanks skipped. A line with no evidence, or a
    non-hex/short sha, is not an exemption -- counted as ignored, so the
    commit counts IN and the output says a line was dropped."""
    out, ignored = {}, 0
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        sha, _, evidence = line.partition(' ')
        if re.fullmatch(r'[0-9a-f]{7,40}', sha) and evidence.strip():
            out[sha] = evidence.strip()
        else:
            ignored += 1
    return out, ignored


def exemptions_text():
    """The exemption list AS COMMITTED at HEAD, not the working tree: an
    uncommitted local edit must not already exempt a commit. Tests patch
    EXEMPT_PATH to a file instead."""
    if EXEMPT_PATH != DEFAULT_EXEMPT_PATH:
        return open(EXEMPT_PATH).read() if os.path.exists(EXEMPT_PATH) else ''
    try:
        return git('show', 'HEAD:' + os.path.relpath(EXEMPT_PATH, ROOT))
    except subprocess.CalledProcessError:
        return ''


def _pr_number(subject):
    """The PR number a commit SUBJECT names, or None. Squash-merge subjects
    end '(#NN)'; a real merge commit says 'Merge pull request #NN'."""
    m = re.search(r'\(#(\d+)\)', subject) or re.search(r'Merge pull request #(\d+)', subject)
    return m.group(1) if m else None


def _counts_toward_pr_threshold(files):
    """True unless EVERY path this commit touches classifies INTERNAL (the
    owner's own framing: "INTERNAL-only PRs... never trigger; they ride
    along"). A PHYSICS commit counts too, even one exempted from the
    at-once trigger above (review's own required fixture: "exempt src
    counts toward 5") -- exemption waives DUE-AT-ONCE only, never this
    count; only a commit that touches NOTHING but INTERNAL paths (docs,
    board, testsys/unit, testsys/regression, ...) is excluded."""
    return any(change_class.classify_path(p) != change_class.INTERNAL
              for p in files if p)


def main():
    """One line, exit 0, whatever happens: an unexpected error prints a DUE
    line naming it (rule 2: fail toward "yes") instead of a traceback."""
    try:
        return _main()
    except Exception as exc:                      # noqa: BLE001
        print('RELEASE DUE: checker error, failing toward due (%s: %s)'
              % (type(exc).__name__, str(exc).splitlines()[0] if str(exc) else ''))
        return 0


def _main():
    try:
        tag = git('describe', '--tags', '--abbrev=0').strip()
    except subprocess.CalledProcessError:
        print('RELEASE DUE: no tag found (cannot measure; failing toward due)')
        return 0
    tag_date = datetime.datetime.fromisoformat(
        git('log', '-1', '--format=%cI', tag).strip())
    days = (datetime.datetime.now(datetime.timezone.utc) - tag_date).days

    exempt, n_ignored = load_exemptions(exemptions_text())
    reasons, n_exempt = [], 0
    counted_prs, total_prs = set(), set()
    # \x1f (ASCII unit separator) as the sha/subject delimiter: a commit
    # subject can contain almost anything BUT that.
    log = git('log', '--format=%H\x1f%s', '%s..HEAD' % tag)
    for line in log.splitlines():
        if not line.strip():
            continue
        sha, _, subject = line.partition('\x1f')
        pr_num = _pr_number(subject)
        files = [p for p in git('diff-tree', '--no-commit-id', '--name-only', '-r', '-m',
                                '--first-parent', sha).split('\n') if p]
        if pr_num:
            total_prs.add(pr_num)
            if _counts_toward_pr_threshold(files):
                counted_prs.add(pr_num)
        if any(sha.startswith(e) for e in exempt):
            n_exempt += 1
            continue
        for path in sorted(set(p for p in files if watched(p))):
            why = classify(path, git('show', '-U0', '--format=',
                                     '--first-parent', '-m', sha, '--', path))
            if why:
                reasons.append('%s in %s' % (why, sha[:7]))

    prs = len(counted_prs)
    counts = ('%d PR(s) counted (%d total PR(s), INTERNAL-only excluded) / '
             '%d days since %s' % (prs, len(total_prs), days, tag))
    ex = ' (%d exempt commit(s), docs/release_exempt.txt)' % n_exempt if n_exempt else ''
    if n_ignored:
        ex += ' (%d malformed exemption line(s) ignored)' % n_ignored
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
