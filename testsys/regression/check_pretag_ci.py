#! /usr/bin/env python3
"""
Pre-`git tag` gate (rule 15 step 6/7): "wait for CI to go green on the
pushed commit, then tag" is a sentence, not a mechanism. Nothing enforced it
before. Run this BEFORE `git tag`; it exits non-zero unless it is safe to.

THE INCIDENT (2026-09-21): `v5.13.1` was tagged at `dfee14d` while the only
CI run for that exact SHA was the run the tag push itself triggered -- at
`git tag` time there was no completed run for `dfee14d` at all. It happened
to conclude green (run 35667535066, 7/7); nothing would have noticed had it
gone red, and a pushed tag cannot be re-pointed (rule 8).

THE STRUCTURAL TRAP: `dfee14d` touches ONLY `pathway_forward.md`, one of
`.github/workflows/test.yml`'s own `paths-ignore` entries, so it could never
have had a pre-tag run of its OWN -- waiting for one would wait forever. The
workflow's comment claims "a release commit always also touches non-md files
... so it still triggers CI regardless"; true of the bundled release commit,
false of the follow-up Tasks-done-row commit rule 15 step 4 itself asks for
when it lands as its own commit rather than folded into step 5's "commit
everything above together". This is why `--ack-paths-ignored-parent`
exists below, rather than the check silently falling back to the parent on
its own: a releaser must choose to accept parent evidence, not have it
handed to them unasked.

Deliberately not auto-run by testsys/run.py's regression tier (its filename
does not match `test_*.py`): its exit codes are not "0 good, nonzero bad" --
PENDING is not a bug in the tree, it means "CI has not finished, ask again
in a minute", and gating every ordinary commit's fast tier on that would
fail CI runs for reasons that have nothing to do with the code under test.
This is a manual, human-invoked gate at exactly one moment: right before
`git tag`.

Exit codes (rule 2: situations that mean different things get different
codes; a script that prints a paragraph and exits 1 either way just means
the reader has to open the log to find out which):

  0  PASS            a completed, successful CI run exists for the SHA
                      being asked about (or, with --ack-paths-ignored-parent,
                      for the nearest ancestor that could have one).
  1  FAIL             a completed CI run for the relevant SHA finished
                      WITHOUT success -- CI ran and failed. Do not tag.
  2  PENDING          no completed run yet -- either none has started, or
                      one is in progress. Wait; do not tag.
  3  PATHS_IGNORED    this SHA touches only paths-ignore'd files and so can
                      never trigger its own CI run. Re-run with
                      --ack-paths-ignored-parent to accept the nearest
                      ancestor's evidence instead, or tag the ancestor.
  4  UNVERIFIED       gh is missing/unauthenticated, or a network call
                      failed. Not evidence of failure (rule 2) -- re-run
                      where gh can actually reach the API.
"""
import argparse
import os
import sys

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(TESTSYS))
from testsys import ci_status  # noqa: E402

EXIT_CODES = {
    ci_status.PASS: 0,
    ci_status.FAIL: 1,
    ci_status.PENDING: 2,
    ci_status.PATHS_IGNORED: 3,
    ci_status.UNVERIFIED: 4,
}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pre-tag', metavar='SHA', required=True,
                   help='the commit you are about to `git tag`')
    p.add_argument('--ack-paths-ignored-parent', action='store_true',
                   help='if SHA cannot trigger its own CI run, accept the '
                        'nearest ancestor commit\'s green run as evidence '
                        'instead (prints exactly which SHA that is)')
    args = p.parse_args()

    try:
        requested_full_sha = ci_status.resolve_sha(args.pre_tag)
    except ValueError as exc:
        print('FAIL  %s' % exc)
        return EXIT_CODES[ci_status.FAIL]

    result = ci_status.evaluate_pretag(
        requested_full_sha, ack_paths_ignored_parent=args.ack_paths_ignored_parent)

    print('%s  %s' % (result.status, result.message))
    if result.evidence_sha and result.evidence_sha != requested_full_sha:
        print('  (evidence sha: %s -- NOT %s itself)'
              % (result.evidence_sha, requested_full_sha))
    return EXIT_CODES[result.status]


if __name__ == '__main__':
    sys.exit(main())
