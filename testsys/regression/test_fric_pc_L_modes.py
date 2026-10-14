#! /usr/bin/env python3
"""
Regression guard (issue #189): par.fric_pc_L means what EQquasi's fric_pc_L
means, for every value, in BOTH codes, end to end (case.setup -> bGlobal.txt
-> bin/eqdyna and `python3 -m eqdyna --backend numpy`):

  None (default, unset) -> L = each node's Dc, the behaviour before the input;
  <= 0                  -> no relaxation, theta_pc = |sigma_n| at once;
  > 0                   -> L in m.

The tiny friclaw-4 case of test_rsf_noncompressive_stop.py (10 steps), with
the caps on and max_norm = -60 MPa, so the -120 MPa initial normal traction
is capped to -60 MPa while theta_pc starts at 120 MPa: the modes then differ
in how fast the strength follows the capped traction. Checks:

  1. case.setup writes the "set L" line: unset -> "0 0.0", a number -> "1 L";
  2. every mode exits 0 in both codes;
  3. unset equals an explicit L = Dc (0.4 m): identical station values in each code;
  4. L = 0 equals L = -1 (both: no relaxation): identical, in each code;
  5. L = 0 differs from unset (the instantaneous strength drop moves the
     station shear by > 0.1 MPa); L = 1e-5 m differs from both;
  6. Fortran and python agree in every mode (station file, |diff| <= 1e-9).

Fails (never skips) when mpirun or bin/eqdyna is missing.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_rsf_noncompressive_stop as base  # noqa: E402  (case + runners)

NSTEPS = 10
STATION = 'faultst000dp005.txt'
CAPS = 'par.max_norm = -60.0e6\n'
MODES = [  # tag, user_defined_params line, the bGlobal.txt line case.setup writes
    ('unset', '', '0 0.0'),
    ('dc', 'par.fric_pc_L = par.fric_rsf_Dc\n', '1 0.4'),
    ('zero', 'par.fric_pc_L = 0.0\n', '1 0.0'),
    ('neg', 'par.fric_pc_L = -1.0\n', '1 -1.0'),
    ('L1e-5', 'par.fric_pc_L = 1e-5\n', '1 1e-05'),
]
TOL_CODES = 1e-9


def _load(path):
    # the data rows only: skip the '#' header and the column-name line
    rows = [ln for ln in open(path)
            if ln.strip() and ln.split()[0][0] in '+-.0123456789']
    return np.loadtxt(rows)


def run_modes(tmp):
    out = {}
    for tag, line, want_line in MODES:
        cd = base._make_case(tmp, tag, tension=False, caps=1, nsteps=NSTEPS,
                             extra=CAPS + line)
        got_line = open(os.path.join(cd, 'bGlobal.txt')).read().splitlines()[-2]
        if got_line != want_line:
            raise AssertionError('%s: bGlobal.txt L line %r, expected %r'
                                 % (tag, got_line, want_line))
        rc_f, o_f = base._run_fortran(cd)
        if rc_f != 0:
            raise AssertionError('%s: Fortran exited %d:\n%s' % (tag, rc_f, o_f[-2000:]))
        f = _load(os.path.join(cd, STATION))
        os.remove(os.path.join(cd, STATION))
        rc_p, o_p = base._run_python(cd, nsteps=NSTEPS)
        if rc_p != 0:
            raise AssertionError('%s: python exited %d:\n%s' % (tag, rc_p, o_p[-2000:]))
        p = _load(os.path.join(cd, STATION))
        out[tag] = (f, p)
    return out


def check(out):
    fails = []

    def same(a, b, code):
        i = 0 if code == 'fortran' else 1
        if not np.array_equal(out[a][i], out[b][i]):
            fails.append('%s: %s != %s (max %.3e)'
                         % (code, a, b, np.abs(out[a][i] - out[b][i]).max()))

    def differ(a, b, atleast):
        d = np.abs(out[a][0] - out[b][0]).max()
        if not d > atleast:
            fails.append('fortran: %s vs %s differ by %.3e, expected > %g' % (a, b, d, atleast))

    for code in ('fortran', 'python'):
        same('unset', 'dc', code)
        same('zero', 'neg', code)
    differ('unset', 'zero', 0.1)
    differ('unset', 'L1e-5', 0.0)
    differ('zero', 'L1e-5', 0.1)
    for tag, _, _ in MODES:
        f, p = out[tag]
        d = np.abs(f - p).max()
        if not d <= TOL_CODES:
            fails.append('%s: Fortran vs python max |diff| %.3e > %g' % (tag, d, TOL_CODES))
    return fails


def main():
    if shutil.which(base.MPIRUN) is None:
        print('FAIL: %s not found' % base.MPIRUN)
        return 1
    try:
        with tempfile.TemporaryDirectory() as tmp:
            fails = check(run_modes(tmp))
    except Exception as e:  # noqa: BLE001 -- report, never skip
        fails = ['%s' % e]
    for f in fails:
        print('  FAIL ', f)
    if not fails:
        print('  PASS  fric_pc_L unset == Dc, 0 == -1, 0 and 1e-5 differ from '
              'unset, Fortran == python in all %d modes' % len(MODES))
    print('%s: %d failure(s)' % ('SUCCESS' if not fails else 'FAILURE', len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
