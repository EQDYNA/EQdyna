#! /usr/bin/env python3
"""
Regression guard: the on-fault station n-stress SIGN is each case's SCEC spec
convention, not a global constant (board row 22a).

THE DEFECT: library_output.f90 wrote column 8 of every faultst*.txt as
`-onFaultQuantHistSCECForm(10,j,i)/1.0d6` -- compression-positive for EVERY
case. test.tpv29 read +182.07 MPa at dp120 where its spec says "Positive means
extension" and EQdyna's own 2015 SCEC submission reads -181.51. Only TPV104 and
TPV105-3D ask for compression-positive. The frt gate never saw it: frt carries
no station column.

What is pinned here, all cheap (no solver, no MPI):
  1. every gated case DECLARES (par.faultStNormalStressSign) exactly the
     convention testsys/matrix.py NSTRESS_CONVENTION says its spec states --
     the table and the declaration are independent, so either being edited
     alone fails here;
  2. scripts/lib.py resolveNormalStressSign maps extension -> +1,
     compression -> -1 and refuses anything else;
  3. compare.nstress_sign_gate passes a correctly-signed run, and FAILS it
     under each mutation: the written sign flipped, the table's convention
     flipped, a zero first-step value, and a run with only a surface station.
The written files themselves are checked by the same gate on every fortran e2e
cell (matrix.ARTIFACTS['fortran'] carries 'nsign').
Exits non-zero on any failure.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

from testsys import compare, matrix  # noqa: E402

HEADER = (' # location = on fault, 0.0 km along strike, 12.0 km down-dip\n'
          ' # Column #8 = normal stress (MPa)\n'
          ' t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress n-stress\n')


def write_station(d, name, nstress):
    with open(os.path.join(d, name), 'w') as f:
        f.write(HEADER)
        for k in range(3):
            f.write('  %.13E' % (0.01 * (k + 1)) + '  %.7E' * 7 % (
                0, 0, 70.0, 0, 0, 0, nstress) + '\n')


def declared(case):
    out = subprocess.run(
        [sys.executable, '-c',
         'from user_defined_params import par; print(par.faultStNormalStressSign)'],
        cwd=os.path.join(ROOT, 'case_input', case), capture_output=True, text=True,
        env=dict(os.environ, PYTHONPATH=os.path.join(ROOT, 'scripts')))
    if out.returncode != 0:
        return 'ERROR: %s' % out.stderr.strip().splitlines()[-1:]
    return out.stdout.strip().splitlines()[-1]


def main():
    fails = []

    # 1. declaration == spec table, for every gated case
    for case in matrix.CASES:
        want = matrix.NSTRESS_CONVENTION[case][0]
        got = declared(case)
        if got != want:
            fails.append('%s declares %r, spec table says %r (%s)'
                         % (case, got, want, matrix.NSTRESS_CONVENTION[case][1]))
    n_comp = sum(v[0] == 'compression' for v in matrix.NSTRESS_CONVENTION.values())

    # 2. the case.setup mapping
    import lib

    class P:
        pass
    for value, sign in (('extension', 1), ('compression', -1)):
        p = P(); p.faultStNormalStressSign = value
        if lib.resolveNormalStressSign(p) != sign:
            fails.append('resolveNormalStressSign(%r) != %+d' % (value, sign))
    for bad in ('positive', None, 1):
        p = P(); p.faultStNormalStressSign = bad
        try:
            lib.resolveNormalStressSign(p)
            fails.append('resolveNormalStressSign(%r) did not raise' % (bad,))
        except ValueError:
            pass

    # 3. the gate, unmutated and under four mutations
    def gate(case, files):
        with tempfile.TemporaryDirectory() as d:
            for name, v in files:
                write_station(d, name, v)
            return compare.nstress_sign_gate(case, d)[0]

    ext, comp = 'test.tpv29', 'test.tpv104'
    good_ext = [('faultst000dp000.txt', 0.0), ('faultst000dp120.txt', -182.07),
                ('faultst-050dp100.txt', -151.7)]
    good_comp = [('faultst000dp075.txt', 120.0)]
    checks = [
        ('unmutated extension case passes', gate(ext, good_ext), True),
        ('unmutated compression case passes', gate(comp, good_comp), True),
        ('MUTATION written sign flipped (the pre-22a output) fails',
         gate(ext, [(n, -v) for n, v in good_ext]), False),
        ('MUTATION one station flipped fails',
         gate(ext, good_ext[:2] + [('faultst-050dp100.txt', 151.7)]), False),
        ('MUTATION compression case written extension-positive fails',
         gate(comp, [('faultst000dp075.txt', -120.0)]), False),
        ('MUTATION zero first-step n-stress fails',
         gate(ext, [('faultst000dp120.txt', 0.0)]), False),
        ('MUTATION only a surface station (nothing checked) fails',
         gate(ext, [('faultst000dp000.txt', 0.0)]), False),
    ]
    saved = matrix.NSTRESS_CONVENTION[ext]
    matrix.NSTRESS_CONVENTION[ext] = ('compression', saved[1])
    try:
        checks.append(('MUTATION table convention flipped fails a correct run',
                       gate(ext, good_ext), False))
    finally:
        matrix.NSTRESS_CONVENTION[ext] = saved
    for label, got, want in checks:
        if got != want:
            fails.append('%s: gate returned %s' % (label, got))

    for f in fails:
        print('FAIL', f)
    print('%s: %d cases declared vs spec table (%d compression), %d gate checks '
          '(%d mutations), %d failures'
          % ('FAIL' if fails else 'SUCCESS', len(matrix.CASES), n_comp,
             len(checks), sum(c[0].startswith('MUTATION') for c in checks),
             len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
