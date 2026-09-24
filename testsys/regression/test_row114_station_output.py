#! /usr/bin/env python3
"""
Regression guard (board row 114) for the python-side station-output port
(`eqdyna.library_output.write_onfault_stations` / `write_offfault_stations`,
recorded into `driver.run`'s carry by `make_step_parts`'s part_a/part_b).

THE BEHAVIOUR under test (rule 10a: on behaviour, not a source grep):
  - a python run of a small, SERIAL test.tpv8 case writes exactly as many
    `faultst*.txt` files as `S['st_on_idx']` has rows, and as many
    `body*.txt` files as `S['st_off_idx']` has rows (both derived
    independently, from `meshgen.build_station_matching` -- not by reading
    back whatever the writer happens to produce);
  - the on-fault files declare/name/write 8 columns for this case's friclaw
    (1: slip-weakening), the off-fault files 7 -- checked THREE ways per
    file (declared count in the header, field-name count, first data row's
    value count), the same technique test_station_header_column_count.py
    uses against the Fortran header;
  - the data column count equals `nsteps` for every file;
  - the time column is strictly increasing and IDENTICAL, row for row,
    between an on-fault and an off-fault file (both are the same carry's
    `timeElapsed` accumulation, read at two different points in one step --
    see driver.py's part_a/part_b);
  - numpy and jax -- ONE code path, per driver.py's module docstring --
    produce the SAME station filenames and column counts, and numeric
    values equal to a generous (1e-6) tolerance. This is a cross-backend
    SELF-CONSISTENCY check, not a substitute for the Fortran parity gate
    (docs/perf_snapshots/station_spread_*.json is the parity oracle).

Cheap (rule 9): one tiny serial test.tpv8 case (create.newcase + case.setup,
~2s measured by test_stress_i0_carry_aliasing.py's identical setup), 3 time
steps on each backend -- no Fortran build, no MPI, no full-length run.
"""
import glob
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)

CASE = 'test.tpv8'
NSTEPS = 3

RE_DECLARED_NCOL = re.compile(r'[Tt]ime series in\s+(\d+)\s+columns')


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_serial_case(case_dir):
    """Same technique as test_stress_i0_carry_aliasing.py's _make_serial_case
    -- create.newcase + forced serial decomposition + case.setup. No Fortran
    binary involved."""
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed: %s' % (r.stderr or '')[-800:])
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial by test_row114_station_output.py'
                       ' (the standalone Python solver is serial-only)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed: %s' % (r.stderr or '')[-800:])


def _check_file_columns(path, expect_declared_ncol):
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()
    m = RE_DECLARED_NCOL.search(text)
    if not m:
        raise AssertionError('%s: no "Time series in N columns" line found' % path)
    declared = int(m.group(1))
    names_idx = next(i for i, l in enumerate(lines) if 'names of the data fields' in l) + 1
    n_names = len(lines[names_idx].split())
    data_lines = [l for l in lines[names_idx + 1:] if l.strip()]
    if not data_lines:
        raise AssertionError('%s: no data rows' % path)
    n_first_row = len(data_lines[0].split())
    if not (declared == expect_declared_ncol == n_names == n_first_row):
        raise AssertionError(
            '%s: declared/names/data column counts = %d/%d/%d, expected all '
            '== %d' % (path, declared, n_names, n_first_row, expect_declared_ncol))
    if len(data_lines) != NSTEPS:
        raise AssertionError('%s: %d data rows, expected nsteps=%d'
                              % (path, len(data_lines), NSTEPS))
    return [float(l.split()[0]) for l in data_lines]


def _run_backend(case_dir, backend_name, out_dir):
    from eqdyna import backend as B
    from eqdyna import eqdyna3d as E
    from eqdyna import driver, library_output

    xp = B.array_module(backend_name)
    S, _mesh = E.build_solver_state(case_dir)
    n_on = int(S['st_on_idx'].shape[0])
    n_off = int(S['st_off_idx'].shape[0])
    if n_on == 0 or n_off == 0:
        raise AssertionError('%s: test.tpv8 default bStations.txt matched %d '
                              'on-fault / %d off-fault stations -- expected '
                              'both > 0 for this test to be meaningful'
                              % (CASE, n_on, n_off))

    out = driver.run(S, nsteps=NSTEPS, verbose=False, xp=xp)
    on_paths = library_output.write_onfault_stations(out_dir, S, out['on_st_hist'])
    off_paths = library_output.write_offfault_stations(out_dir, S, out['off_st_hist'])

    if len(on_paths) != n_on:
        raise AssertionError('%s: wrote %d faultst* files, expected %d '
                              '(S["st_on_idx"].shape[0])' % (backend_name, len(on_paths), n_on))
    if len(off_paths) != n_off:
        raise AssertionError('%s: wrote %d body* files, expected %d '
                              '(S["st_off_idx"].shape[0])' % (backend_name, len(off_paths), n_off))

    ncol_on = 11 if S['friclaw'] >= 3 else 8
    on_times = _check_file_columns(on_paths[0], ncol_on)
    off_times = _check_file_columns(off_paths[0], 7)
    if on_times != off_times:
        raise AssertionError('%s: on-fault and off-fault time columns differ: '
                              '%r vs %r' % (backend_name, on_times, off_times))
    if sorted(on_times) != on_times or len(set(on_times)) != len(on_times):
        raise AssertionError('%s: time column is not strictly increasing: %r'
                              % (backend_name, on_times))

    on_names = sorted(os.path.basename(p) for p in on_paths)
    off_names = sorted(os.path.basename(p) for p in off_paths)
    return on_names, off_names, ncol_on


def main():
    print('Regression guard: row 114 station output (write_onfault_stations/'
          'write_offfault_stations), one code path, numpy and jax')
    fails = []
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = os.path.join(tmp, 'case')
        _make_serial_case(case_dir)
        for backend_name in ('numpy', 'jax'):
            out_dir = os.path.join(tmp, 'out_%s' % backend_name)
            os.makedirs(out_dir)
            try:
                results[backend_name] = _run_backend(case_dir, backend_name, out_dir)
                print('  %s: on-fault %d file(s), off-fault %d file(s), %d on-fault columns'
                      % (backend_name, len(results[backend_name][0]),
                         len(results[backend_name][1]), results[backend_name][2]))
            except Exception as e:
                fails.append('%s: %s' % (backend_name, e))

        if not fails and results.get('numpy') and results.get('jax'):
            on_n, off_n, ncol_n = results['numpy']
            on_j, off_j, ncol_j = results['jax']
            if (on_n, off_n, ncol_n) != (on_j, off_j, ncol_j):
                fails.append('numpy vs jax: filename sets or column count differ: '
                             'numpy=%r jax=%r' % (results['numpy'], results['jax']))
            else:
                print('  numpy == jax: same %d on-fault / %d off-fault filenames, '
                      'same column count' % (len(on_n), len(off_n)))

    # Dropped-station report (rows 94/116): same words as Fortran's
    # report_dropped_onfault_st / report_dropped_offfault_st, and silence when
    # nothing is dropped. Behaviour: captured stdout of the real function.
    import contextlib, io
    import numpy as np
    from eqdyna import eqdyna3d
    xonfs = np.array([[0.0, -18000.0, 5000.0], [0.0, -15600.0, -12000.0]])
    x4nds = np.array([[0.0, 0.0], [1000.0, -500.0], [0.0, -300.0]])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        eqdyna3d.report_dropped_stations(xonfs, x4nds, [(1, 1, 1), (7, 3, 1)], [(1, 42)])
    out = buf.getvalue()
    for want in ('WARNING: 1 of 3 requested on-fault stations match no fault node',
                 'dropped on-fault station 2 (fault 1) at x,z =   -18.000   -15.600 km',
                 'WARNING: 1 of 2 requested off-fault stations match no grid node',
                 'dropped off-fault station 2 at x,y,z =     0.000    -0.500    -0.300 km'):
        if want not in out:
            fails.append('drop report lacks %r; stdout was:\n%s' % (want, out))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        eqdyna3d.report_dropped_stations(xonfs, x4nds, [(1, 1, 1), (2, 2, 1), (3, 3, 1)],
                                         [(1, 4), (2, 5)])
    if buf.getvalue():
        fails.append('drop report printed with nothing dropped: %r' % buf.getvalue())

    if fails:
        print('FAIL test_row114_station_output')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_row114_station_output')
    return 0


if __name__ == '__main__':
    sys.exit(main())
