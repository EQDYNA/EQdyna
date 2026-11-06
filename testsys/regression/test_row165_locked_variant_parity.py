#! /usr/bin/env python3
"""Row 165 reproducer (rule 10): fortran (4 ranks) vs python-jax (serial) on a synthetic
locked-node variant of the test.restart2 fixture (row165_locked_variant.py).

MEASURED 2026-11-05 (this box): at nt=1 the two agree to <= 1e-15 on every column (35 nodes differ, tdip
roundoff ~1e-18); tdip, pure roundoff noise, grows ~1e3 per step (nt=2 1e-9, nt=8 3e-3 on strike columns,
tdip 6.5e3 Pa); at 5 s max|diff| is 4.2e5 Pa on tstk (col 12) of 1.4e7. Controls: fortran -O3 == fortran
-O0 bit-for-bit (0.0), but the SAME fortran source built -O3 -ffast-math differs from fortran by 4.1e5 Pa
(jax: 4.2e5; numpy backend vs fortran 4.0e5), so a roundoff-only change of the reference itself moves the
answer as much as the port does: an ill-conditioned (roundoff-amplifying) variant, not a port divergence.
The test DECLARES that: it fails if the first-step agreement is lost (a real port divergence) or if the
amplification disappears (then the declaration must be removed and the bound tightened).
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import row165_locked_variant as R  # noqa: E402

DT = 5.0 / 240
FIRST_STEP_BOUND = 1e-15      # nt=1: roundoff agreement, any column
AMPLIFIED_AT_5S = 1e4         # declared: measured 4.2e5 Pa on tstk


def run_pair(root, fx, nt):
    term = (nt - 0.5) * DT
    fc = R.stage(root, "f%d" % nt, fx, term)
    jc = R.stage(root, "j%d" % nt, fx, term)
    return R.run_fortran(fc), R.run_python(jc)


def main():
    if shutil.which("mpirun") is None or shutil.which("eqdyna") is None:
        raise SystemExit("FAIL: mpirun/eqdyna not on PATH (install-eqdyna.sh, PATH=$EQDYNAROOT/bin)")
    fails = 0
    with tempfile.TemporaryDirectory() as root:
        fx = os.path.join(root, "variant.nc")
        R.make_variant(fx)
        for nt, name in ((1, "first step"), (240, "5 s")):
            f, j = run_pair(root, fx, nt)
            rep = sorted(R.diff_report(f, j, name), key=lambda r: -r[1])
            c, dmax, row, xyz, fv, jv, nd = rep[0]
            print("nt=%d worst col %d max|diff| %.3e at x=%g z=%g (f %.9e j %.9e), %d nodes differ in that col" % (nt, c, dmax, xyz[0], xyz[2], fv, jv, nd))
            if nt == 1 and dmax > FIRST_STEP_BOUND:
                print("FAIL first-step disagreement %.3e > %g: a real port divergence" % (dmax, FIRST_STEP_BOUND)); fails += 1
            if nt == 240:
                d12 = max(r[1] for r in rep if r[0] == 12)
                if d12 < AMPLIFIED_AT_5S:
                    print("FAIL tstk max|diff| %.3e < %g: amplification gone, tighten the bound" % (d12, AMPLIFIED_AT_5S)); fails += 1
    print("SUCCESS: 0 failure(s)" if not fails else "FAILURE: %d failure(s)" % fails)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
