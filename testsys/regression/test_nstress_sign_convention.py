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


BGLOBAL = ('1\n1\n1\n0\n0\n1\n1\n1\n8\n0\n0\n0\n \n2 2 1\n \n5.0\n'
           '0.043736878936319105\n \n1 3\n2670.0 1000.0 0.66\n0.1 6000.0\n \n'
           '0.0 0.0 -12000.0\n1500.0 -9999.0 50000000.0 1.0\n45.0 0.33\n'
           '0.75 5000000.0\n270.0 90.0\n0.001\n \n0.28868360277136257\n'
           '0.0 0.0\n5000.0 2000.0 8000.0\n1\n0\n0 0.0\n-1 -10000000.0 -40000000.0\n')
# The file's last two lines: the theta_pc relaxation length as "set L"
# (par.fric_pc_L; set 0 = unset, each node's Dc; set 1 = L, <= 0 no
# relaxation) and the normal-stress caps switch and values
# (par.C_normal_stress_caps, min_norm, max_norm), as case.setup writes the
# defaults.
TAIL = '0 0.0\n-1 -10000000.0 -40000000.0\n'
# Everything before the sign line; the file then has the sign line, the
# fault time-series switch (par.outputFaultTimeSeries) and TAIL, in that order.
HEAD = BGLOBAL[:-len('1\n0\n' + TAIL)]
# Everything before TAIL.
PRE_TAIL = HEAD + '1\n0\n'


def write_station(d, name, nstress):
    # nstress == 'zeros' writes the phantom-station shape: every value 0.
    shear = 0.0 if nstress == 'zeros' else 70.0
    nstress = 0.0 if nstress == 'zeros' else nstress
    with open(os.path.join(d, name), 'w') as f:
        f.write(HEADER)
        for k in range(3):
            f.write('  %.13E' % (0.01 * (k + 1)) + '  %.7E' * 7 % (
                0, 0, shear, 0, 0, 0, nstress) + '\n')


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
    checks.append(('MUTATION all-zero phantom station file fails',
                   gate(ext, good_ext + [('faultst000dp000.txt', 'zeros')]), False))
    checks.append(('MUTATION unparseable station name fails',
                   gate(ext, good_ext + [('faultst000dpXYZ.txt', -1.0)]), False))
    for tok, want in (('0.1234567-100', 1.234567e-101), ('-0.5000000+100', -0.5e100),
                      ('0.2544000E+02', 25.44)):
        got = compare._fortran_float(tok)
        checks.append(('Fortran value %s parses' % tok, abs(got - want) <= 1e-12 * abs(want), True))
    for label, got, want in checks:
        if got != want:
            fails.append('%s: gate returned %s' % (label, got))

    # 4. read_bglobal: the trailing lines are read, and a missing or invalid
    # one is REFUSED, not defaulted (rule 2). BGLOBAL is test.tpv8's real
    # case.setup output (5 s gate term), verbatim.
    sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
    from eqdyna import checkInputConsistency, readInputFiles
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'bGlobal.txt')
        # A missing line is a stale file (ValueError, as before); an invalid
        # value is refused with Fortran's numbered code 15 (rule 23).
        # The time-series switch after it is pinned the same way: missing is
        # stale, anything but 0/1 is refused with code 18.
        variants = [('nStressOutSign', 'as written (+1)', BGLOBAL, 1),
                    ('nStressOutSign', 'compression (-1)', HEAD + '-1\n0\n' + TAIL, -1),
                    ('nStressOutSign', 'MUTATION sign line onward missing', HEAD, 'stale'),
                    ('nStressOutSign', 'MUTATION value 0', HEAD + '0\n0\n' + TAIL, 'exit 15'),
                    ('nStressOutSign', 'MUTATION value 2', HEAD + '2\n0\n' + TAIL, 'exit 15'),
                    ('outputFaultTimeSeries', 'as written (0)', BGLOBAL, 0),
                    ('outputFaultTimeSeries', 'on (1)', HEAD + '1\n1\n' + TAIL, 1),
                    ('outputFaultTimeSeries', 'MUTATION switch line onward missing', HEAD + '1\n', 'stale'),
                    ('outputFaultTimeSeries', 'MUTATION value 2', HEAD + '1\n2\n' + TAIL, 'exit 18'),
                    ('outputFaultTimeSeries', 'MUTATION value -1', HEAD + '1\n-1\n' + TAIL, 'exit 18'),
                    # The rate-and-state normal-stress lines (issues 189/200):
                    # missing is stale, an invalid switch or caps order is
                    # refused with code 10, as readglobal does.
                    ('fric_pc_L', 'as written (unset, each node Dc)', BGLOBAL, None),
                    ('fric_pc_L', 'set (0.004 m)', PRE_TAIL + '1 0.004\n-1 -10000000.0 -40000000.0\n', 0.004),
                    ('fric_pc_L', 'set 0 (no relaxation)', PRE_TAIL + '1 0.0\n-1 -10000000.0 -40000000.0\n', 0.0),
                    ('fric_pc_L', 'set -1 (no relaxation)', PRE_TAIL + '1 -1.0\n-1 -10000000.0 -40000000.0\n', -1.0),
                    ('fric_pc_L', 'unset ignores the value', PRE_TAIL + '0 0.5\n-1 -10000000.0 -40000000.0\n', None),
                    ('fric_pc_L', 'MUTATION both lines missing', PRE_TAIL, 'stale'),
                    ('fric_pc_L', 'MUTATION not a number', PRE_TAIL + '1 x\n-1 -1e7 -4e7\n', 'stale'),
                    ('fric_pc_L', 'MUTATION old one-number line', PRE_TAIL + '0.004\n-1 -1e7 -4e7\n', 'stale'),
                    ('fric_pc_L', 'MUTATION set flag 2', PRE_TAIL + '2 0.1\n-1 -1e7 -4e7\n', 'stale'),
                    ('C_normal_stress_caps', 'as written (-1)', BGLOBAL, -1),
                    ('C_normal_stress_caps', 'on (1)', PRE_TAIL + '0 0.0\n1 -5e6 -60e6\n', 1),
                    ('min_norm', 'set (-5 MPa)', PRE_TAIL + '0 0.0\n1 -5e6 -60e6\n', -5e6),
                    ('max_norm', 'set (-60 MPa)', PRE_TAIL + '0 0.0\n1 -5e6 -60e6\n', -60e6),
                    ('max_norm', 'equal to min_norm', PRE_TAIL + '0 0.0\n1 -5e6 -5e6\n', -5e6),
                    ('C_normal_stress_caps', 'MUTATION caps line missing', PRE_TAIL + '0 0.0\n', 'stale'),
                    ('C_normal_stress_caps', 'MUTATION caps line short', PRE_TAIL + '0 0.0\n-1 -1e7\n', 'stale'),
                    ('C_normal_stress_caps', 'MUTATION switch 2', PRE_TAIL + '0 0.0\n2 -1e7 -4e7\n', 'exit 10'),
                    ('C_normal_stress_caps', 'MUTATION switch -2', PRE_TAIL + '0 0.0\n-2 -1e7 -4e7\n', 'exit 10'),
                    ('C_normal_stress_caps', 'MUTATION min_norm < max_norm', PRE_TAIL + '0 0.0\n-1 -4e7 -1e7\n', 'exit 10'),
                    ('C_normal_stress_caps', 'MUTATION min_norm 0', PRE_TAIL + '0 0.0\n-1 0.0 -4e7\n', 'exit 10'),
                    ('C_normal_stress_caps', 'MUTATION min_norm tensile', PRE_TAIL + '0 0.0\n0 1e6 -4e7\n', 'exit 10')]
        for key, label, text, want in variants:
            with open(path, 'w') as f:
                f.write(text)
            try:
                got = readInputFiles.read_bglobal(path)[key]
            except checkInputConsistency.InputConsistencyError as e:
                got = 'exit %d' % e.code
            except ValueError:
                got = 'stale'
            checks.append(('read_bglobal %s %s' % (key, label), got, want))
            if got != want:
                fails.append('read_bglobal %s %s: got %r, want %r' % (key, label, got, want))

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
