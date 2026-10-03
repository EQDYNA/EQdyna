#! /usr/bin/env python3
"""
Regression guard: user docs must not drift from shipped features (owner,
2026-10-03, "Go. But you should keep docs synced right? Why not?" --
pathway_forward.md section A header note). Checks two COVERAGE properties,
not style or command-resolution (test_user_docs_style.py and
test_user_docs_commands.py already own those):

  1. Every gated case in testNameList.nameList is named somewhere in
     docs/user/benchmarks.md -- a case with no mention at all is a case the
     user docs never told anyone exists.
  2. Every output-filename family this codebase's solver and
     post-processing actually emit for a multi-fault case
     (scripts/lib.py:faultTag, src/fortran/library_output.f90,
     scripts/plotRuptureDynamics) is named somewhere in docs/user/outputs.md.

Measured gap this guard closes (2026-10-03): test.tpv22/test.tpv23 (item 17,
landed weeks earlier) were entirely absent from benchmarks.md, and the
per-fault output-naming convention (faultTag/'ft<N>_') was undocumented in
outputs.md. Neither test_user_docs_commands.py nor test_user_docs_style.py
would have caught either gap -- both check STYLE and COMMAND RESOLUTION, not
whether a shipped feature is mentioned anywhere at all. Confirmed RED on the
pre-sync docs (session log / PR description carries the run), GREEN after.

Cheap (rule 9): two file reads and a list membership check, well under a
second.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import testNameList  # noqa: E402

BENCHMARKS_MD = os.path.join(ROOT, 'docs', 'user', 'benchmarks.md')
OUTPUTS_MD = os.path.join(ROOT, 'docs', 'user', 'outputs.md')

# The base (pre-faultTag-suffix) output-filename families this codebase's
# solver and post-processing actually emit -- hand-maintained, like
# test_user_docs_commands.py's own checks list. A brand new output family
# needs a line added HERE and in outputs.md; this guard only catches drift
# between the two already-named sides, not a family neither side has named.
OUTPUT_PATTERNS = (
    'faultst',
    'src_evol',
    'fault.dyna.r',
    'SCECRuptureTime',
    'cRuptureDynamics',
    'frt.txt',
)


def check_benchmarks_coverage():
    text = open(BENCHMARKS_MD).read()
    missing = [c for c in testNameList.nameList if c not in text]
    assert not missing, (
        'docs/user/benchmarks.md never mentions: %r -- a gated case with no '
        'user-doc coverage at all' % missing)


def check_outputs_coverage():
    text = open(OUTPUTS_MD).read()
    missing = [p for p in OUTPUT_PATTERNS if p not in text]
    assert not missing, (
        'docs/user/outputs.md never mentions these output-filename families: '
        '%r' % missing)


def main():
    print('Regression guard: user docs must cover every gated case and output family')
    checks = (
        ('benchmarks.md case coverage', check_benchmarks_coverage),
        ('outputs.md output-family coverage', check_outputs_coverage),
    )
    failures = []
    for label, fn in checks:
        try:
            fn()
            print('  PASS  %s' % label)
        except AssertionError as e:
            failures.append(str(e))
            print('  FAIL  %s: %s' % (label, e))
    if failures:
        print('\nFAIL test_user_docs_coverage (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_user_docs_coverage: %d case(s), %d output pattern(s) covered'
          % (len(testNameList.nameList), len(OUTPUT_PATTERNS)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
