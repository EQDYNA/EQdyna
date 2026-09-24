#! /usr/bin/env python3
"""
Regression guard: the station time-series gate and the python-jax nc gate
(owner gate design, 2026-09-24).

THE DEFECT SHAPE: the e2e gate compared ONLY frt.canonical.txt, so a wrong
sign in a station column (board row 22a: n-stress compression-positive in 9
of 11 cases) passed every gate for years. Every cell is now also gated on
selected on/off-fault station series (matrix.GATE_STATIONS, normalized,
matrix.STATION_BOUND) and on fault.dyna.r.nc, python-jax included.

What is pinned, on BEHAVIOUR (rule 10a), no solver run: a run directory is
built from the committed references themselves (so unmutated e == 0 exactly),
then mutated, and compare.station_gate / compare.compare_nc_files must say:
  - unmutated                                            -> pass
  - one column's sign flipped (the 22a shape)            -> FAIL
  - one sample perturbed just ABOVE the case bound       -> FAIL
  - the same perturbation just BELOW the bound           -> pass
  - a selected station file missing                      -> FAIL
  - the time axis shifted                                -> FAIL
  - one row missing                                      -> FAIL
  - nonzero values in an identically-zero column         -> FAIL
  - one fault.dyna.r.nc variable perturbed               -> FAIL (nc)
and every case must resolve the station gate one way (a bound, or a declared
reason) with a committed reference for every selected file.
Exits non-zero on any failure.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from testsys import compare, matrix  # noqa: E402

CASE = 'test.tpv8'


def write_station(path, header_lines, names, data):
    with open(path, 'w') as f:
        f.writelines(header_lines)
        f.write(' ' + ' '.join(names) + '\n')
        for row in data:
            f.write(' '.join(repr(float(v)) for v in row) + '\n')


def header_of(path):
    return [l for l in open(path) if l.split() and l.split()[0].startswith('#')]


def build_run(tmp, case=CASE):
    run = os.path.join(tmp, 'run')
    os.makedirs(run)
    for kind in ('on', 'off'):
        for fn in matrix.GATE_STATIONS[case][kind]:
            shutil.copy(compare.station_reference_path(case, fn), os.path.join(run, fn))
    return run


def mutate(run, fn, fun):
    p = os.path.join(run, fn)
    names, data = compare.read_station_file(p)
    write_station(p, header_of(p), names, fun(names, data.copy()))


def scale(kind, col_name):
    """S_q exactly as compare.normalize_station_error computes it."""
    s = 0.0
    for fn in matrix.GATE_STATIONS[CASE][kind]:
        names, d = compare.read_station_file(compare.station_reference_path(CASE, fn))
        s = max(s, float(np.max(np.abs(d[:, names.index(col_name)]))))
    return max(s, matrix.STATION_ZERO_FLOOR)


def main():
    fails, checks = [], []

    # table completeness: every case one way, every selected file committed
    for c in matrix.CASES:
        if c not in matrix.STATION_BOUND:
            fails.append('%s has no STATION_BOUND' % c)
        for kind in ('on', 'off'):
            for fn in matrix.GATE_STATIONS[c][kind]:
                try:
                    compare.station_reference_path(c, fn)
                except FileNotFoundError as e:
                    fails.append(str(e))
    for b in ('fortran', 'python-jax'):
        for a in ('frt', 'nc', 'nsign', 'station'):
            if a not in matrix.ARTIFACTS[b]:
                fails.append('%s lost artifact %r' % (b, a))

    bound = matrix.STATION_BOUND[CASE]
    on0 = matrix.GATE_STATIONS[CASE]['on'][0]
    S = scale('on', 'h-shear-stress')

    def bump(factor):
        def f(names, d):
            c = names.index('h-shear-stress')
            d[len(d) // 2, c] += factor * bound * S
            return d
        return f

    def flip(names, d):
        c = names.index('n-stress')
        d[:, c] = -d[:, c]
        return d

    def nan_sample(names, d):
        d[len(d) // 2, names.index('h-slip-rate')] = float('nan')
        return d

    def tshift(names, d):
        d[:, 0] += 1e-3
        return d

    # The identically-zero column: test.tpv104's on-fault 'temperature' is 0.0
    # at every selected station in the reference (no gated tpv8 column is), so
    # S_q = 0 and the scale is STATION_ZERO_FLOOR: noise of 1e-3 must fail
    # and 1e-12 (e = 1e-6, under tpv104's 1e-5) must pass.
    ZCASE, ZCOL = 'test.tpv104', 'temperature'

    def zero_col(level):
        def f(names, d):
            d[:, names.index(ZCOL)] = level
            return d
        return f

    scenarios = [
        ('unmutated', None, True),
        ('MUTATION n-stress sign flipped', flip, False),
        ('MUTATION one sample +2x bound', bump(2.0), False),
        ('below bound: one sample +0.5x bound', bump(0.5), True),
        ('MUTATION time axis shifted 1 ms', tshift, False),
        ('MUTATION one row dropped', lambda n, d: d[:-1], False),
        ('MUTATION selected file missing', 'delete', False),
        ('MUTATION one NaN sample', nan_sample, False),
        ('MUTATION field-name line changed', 'rename', False),
        ('MUTATION 1e-3 in an identically-zero column', zero_col(1e-3), False, ZCASE),
        ('below floor: 1e-12 in an identically-zero column', zero_col(1e-12), True, ZCASE),
    ]
    # the zero-column scenario must really target a zero column
    vs = [compare.read_station_file(compare.station_reference_path(ZCASE, fn))
          for fn in matrix.GATE_STATIONS[ZCASE]['on']]
    if any(np.any(d[:, n.index(ZCOL)]) for n, d in vs):
        fails.append('%s is not identically zero in the %s references; the '
                     'zero-column scenarios would test nothing' % (ZCOL, ZCASE))
    for sc in scenarios:
        label, fun, want = sc[:3]
        case = sc[3] if len(sc) > 3 else CASE
        with tempfile.TemporaryDirectory() as tmp:
            run = build_run(tmp, case)
            target = matrix.GATE_STATIONS[case]['on'][0]
            if fun == 'delete':
                os.remove(os.path.join(run, target))
            elif fun == 'rename':
                p = os.path.join(run, target)
                names, d = compare.read_station_file(p)
                write_station(p, header_of(p), names[:-1] + ['n-stress-renamed'], d)
            elif fun is not None:
                mutate(run, target, fun)
            got, lines = compare.station_gate(case, run)
        checks.append(label)
        if got != want:
            fails.append('%s: gate returned %s, wanted %s:\n    %s'
                         % (label, got, want, '\n    '.join(lines[:3])))

    # nc: the jax column's comparison is compare_nc_files against the ONE
    # committed fault.dyna.r.nc; a perturbed variable must fail it.
    from netCDF4 import Dataset
    with tempfile.TemporaryDirectory() as tmp:
        ref = compare.reference_path(CASE, compare.NC_NAME)
        cp = os.path.join(tmp, compare.NC_NAME)
        shutil.copy(ref, cp)
        if not compare.compare_nc_files(ref, cp).startswith('SUCCESS'):
            fails.append('nc: an unmutated copy of the reference did not pass')
        with Dataset(cp, 'a') as ds:
            v = 'shear_strike'  # a physics variable, not a coordinate
            a = np.asarray(ds.variables[v][:], dtype=float)
            a.flat[a.size // 2] += 0.01 * float(np.max(np.abs(a))) + 1.0
            ds.variables[v][:] = a
        if compare.compare_nc_files(ref, cp).startswith('SUCCESS'):
            fails.append('nc: MUTATION variable %r perturbed by 1% of its max still passed' % v)
        checks += ['nc unmutated', 'MUTATION nc variable perturbed']

    for f in fails:
        print('FAIL', f)
    print('%s: station gate %d scenarios (%d mutations) on %s, bound %.0e, '
          'S_q(h-shear-stress)=%.4g; %d cases resolved; %d failures'
          % ('FAIL' if fails else 'SUCCESS', len(checks),
             sum(c.startswith('MUTATION') for c in checks), CASE, bound, S,
             len(matrix.CASES), len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
