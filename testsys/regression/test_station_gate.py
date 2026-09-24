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


def build_run(tmp):
    run = os.path.join(tmp, 'run')
    os.makedirs(run)
    for kind in ('on', 'off'):
        for fn in matrix.GATE_STATIONS[CASE][kind]:
            shutil.copy(compare.station_reference_path(CASE, fn), os.path.join(run, fn))
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
        if (c in matrix.STATION_BOUND) == (c in matrix.STATION_UNSUPPORTED_CASES):
            fails.append('%s resolves the station gate %s ways' % (
                c, 'two' if c in matrix.STATION_BOUND else 'zero'))
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

    def tshift(names, d):
        d[:, 0] += 1e-3
        return d

    def zero_col_noise(names, d):
        # the identically-zero column: v-slip on this vertical strike-slip
        # fault is zero at every selected station in the reference
        c = names.index('v-slip')
        d[:, c] = 1e-3
        return d

    scenarios = [
        ('unmutated', None, True),
        ('MUTATION n-stress sign flipped', flip, False),
        ('MUTATION one sample +2x bound', bump(2.0), False),
        ('below bound: one sample +0.5x bound', bump(0.5), True),
        ('MUTATION time axis shifted 1 ms', tshift, False),
        ('MUTATION one row dropped', lambda n, d: d[:-1], False),
        ('MUTATION nonzero in an identically-zero column', zero_col_noise, False),
        ('MUTATION selected file missing', 'delete', False),
    ]
    # the zero-column scenario must really target a zero column
    vs = [compare.read_station_file(compare.station_reference_path(CASE, fn))
          for fn in matrix.GATE_STATIONS[CASE]['on']]
    if any(np.any(d[:, n.index('v-slip')]) for n, d in vs):
        fails.append('v-slip is not identically zero in the %s references; the '
                     'zero-column scenario would test nothing' % CASE)
    for label, fun, want in scenarios:
        with tempfile.TemporaryDirectory() as tmp:
            run = build_run(tmp)
            if fun == 'delete':
                os.remove(os.path.join(run, on0))
            elif fun is not None:
                mutate(run, on0, fun)
            got, lines = compare.station_gate(CASE, run)
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
            v = next(iter(ds.variables))
            a = np.asarray(ds.variables[v][:], dtype=float)
            a.flat[a.size // 2] += 1.0
            ds.variables[v][:] = a
        if compare.compare_nc_files(ref, cp).startswith('SUCCESS'):
            fails.append('nc: MUTATION variable %r perturbed by 1.0 still passed' % v)
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
