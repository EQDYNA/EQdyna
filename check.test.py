#!/usr/bin/env python3
"""
Compare an ALREADY-COMPLETED Fortran run tree against test.reference.results/
(PROJECT_RULES.md rules 2, 3, 5, 7).

This is the fortran column of the e2e sweep's comparison, for a tree that has
already been run -- `python3 check.test.py` after a by-hand mpirun, without
re-running anything. It is NOT a second comparison implementation: it calls
testsys/compare.py's compare_cell, the same function testsys/e2e/run_e2e.py
calls for every cell of the sweep. One implementation, two entry points.

WHAT CHANGED, AND WHY IT MATTERED
This file used to carry

    fileNameList = ['fault.dyna.r.nc','frt.txt0','frt.txt1','frt.txt2','frt.txt3']

-- a gate whose file list was a statement about npx*npy*npz. It also skipped
(`continue`) any file absent from the reference, so a case whose reference had
2 rank files was silently compared on 2 of the 4 names and a case with 4 on all
4, with no line saying which. frt is now compared in its canonical,
decomposition-independent form (testsys/frt_canonical.py): ONE reference file
per case, one row per fault node, ordered by position. A serial run, a 2-rank
run and a 4-rank run are all compared the same way, to the same artifact.

fault.dyna.r.nc is still compared, deliberately. It is lossy relative to frt
(12 resampled dip x strike variables against 18 physics columns at native node
resolution), so it adds no physics coverage -- but it is the only thing that
exercises scripts/plotRuptureDynamics, which would otherwise be gated by
nothing. See matrix.ARTIFACTS.

"Pass" means every selected case compared clean and the process exits 0. A
FAIL, a missing artifact, or zero comparisons run is a non-zero exit (rule 2).
"""
import argparse
import os
import sys

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from testsys import compare, matrix  # noqa: E402

THRESHOLD = matrix.THRESHOLD  # rule 5's one outer sanity bound, defined once.


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--test-root', default='test',
                    help='tree holding the completed run(s) (default: test)')
    ap.add_argument('--cases', help='comma-separated subset of the case axis '
                                    '(default: every gated case)')
    args = ap.parse_args(argv)

    cases = args.cases.split(',') if args.cases else list(matrix.CASES)
    unknown = [c for c in cases if c not in matrix.CASES]
    if unknown:
        print('check.test.py: FAIL - unknown case(s) %s; gated cases are %s'
              % (unknown, ', '.join(matrix.CASES)))
        return 1

    print('check.test.py: comparing %d case(s) of the fortran column against '
          '%s: %s' % (len(cases), compare.REFERENCE_ROOT, ', '.join(cases)))
    print('check.test.py: artifacts=%s, one bound per case (rule 5), outer '
          'sanity threshold=%.0e' % ('+'.join(matrix.ARTIFACTS['fortran']),
                                     THRESHOLD))

    failures = []
    for case in cases:
        run_dir = os.path.join(args.test_root, case)
        print(' ')
        print('-- %s (%s, gate: %s) --'
              % (case, run_dir, matrix.gate_description(case)))
        if not os.path.isdir(run_dir):
            print('FAIL missing run directory %s' % run_dir)
            failures.append(case)
            continue
        try:
            ok, lines = compare.compare_cell(case, 'fortran', run_dir)
        except Exception as exc:            # noqa: BLE001 - reported, not swallowed
            ok, lines = False, ['%s: %s' % (type(exc).__name__, exc)]
        for line in lines:
            print('   ' + line)
        print('%s %s' % ('SUCCESS' if ok else 'FAIL', case))
        if not ok:
            failures.append(case)

    print(' ')
    if not cases:
        print('check.test.py: FAIL - no comparisons were run; a check that '
              'compared nothing must never exit green')
        return 1
    print('check.test.py: %d/%d case(s) SUCCESS'
          % (len(cases) - len(failures), len(cases)))
    if failures:
        print('check.test.py: FAIL -', len(failures), 'case(s) failed:',
              ', '.join(failures))
        return 1
    print('check.test.py: SUCCESS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
