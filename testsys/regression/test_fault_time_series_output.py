#! /usr/bin/env python3
"""
Board item 158: `output_fault_ts` (`src/fortran/library_output.f90:901-988`,
PR #177, par.outputFaultTimeSeries) had no gated or regression test -- this
is that guard.

METHOD: a real, tiny (forced serial, few-step) test.tpv8 case.

  (a) flag OFF (default 0): no `src_ts*` file of any kind is written.
  (b) flag ON: `src_ts_hdr0` + `src_ts0` are written; the header's declared
      shape (nnode, nfield, nstep) matches what the binary actually holds
      (size in bytes == 8 * nnode * nfield * nstep_written, the same
      independent re-derivation `readFaultTimeSeries.read_fault_ts` uses to
      drop a partial trailing record), and the LAST record's `slip_strike`
      / `traction_strike` at the fault node nearest one of this case's own
      on-fault stations equal that station's own `faultst*.txt` file (SCEC
      columns 2 / 4, the writer's documented units m / MPa) to within half
      a unit of the printed digit -- the one quantity (MPa shear stress)
      measurably rules out a field-order swap, a unit-scale slip, or the
      writer picking the wrong fault node, the three defects this flag's
      correspondence-doc entry says to watch for. Not ported: the
      correspondence doc
      (docs/fortran_python_correspondence.md) records `output_fault_ts` as
      a Fortran-only writer, refused by the Python backend when the flag is
      set -- there is no python counterpart to test here.

Cheap (rule 9): forced single-rank (so every fault node is in one file, no
cross-rank concatenation needed), term = 6*dt.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

CASE = 'test.tpv8'
NSTEP = 6


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'), os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(case_dir, env, flag):
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=env, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed: %s' % (r.stderr or '')[-800:])
    params_path = os.path.join(case_dir, 'user_defined_params.py')
    with open(params_path) as f:
        text = f.read().rstrip('\n')
    with open(params_path, 'w') as f:
        f.write(text + '\n\n# item 158 probe: forced serial, tiny term\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n'
                       'par.term = par.dt * %d\n'
                       'par.outputFaultTimeSeries = %d\n' % (NSTEP, flag))
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise AssertionError('case.setup failed: %s' % (r.stderr or '')[-800:])


def _run(case_dir, env):
    r = subprocess.run(['mpirun', '-np', '1', 'eqdyna'], cwd=case_dir, env=env,
                       capture_output=True, text=True, timeout=180)
    if r.returncode != 0:
        raise AssertionError('eqdyna exited %d:\n%s' % (r.returncode, (r.stdout + r.stderr)[-2000:]))


def main():
    if shutil.which('mpirun') is None:
        print('FAIL test_fault_time_series_output: mpirun not on PATH '
              '(this guard requires a real MPI run, it does not skip silently)')
        return 1
    binary = os.path.join(ROOT, 'bin', 'eqdyna')
    if not os.path.exists(binary):
        print('FAIL test_fault_time_series_output: %s not built '
              '(run install-eqdyna.sh first -- this guard does not skip silently)' % binary)
        return 1
    env = _env()
    fails = []

    with tempfile.TemporaryDirectory() as tmp:
        # (a) flag off: no src_ts* file at all
        off_dir = os.path.join(tmp, 'off')
        _make_case(off_dir, env, 0)
        _run(off_dir, env)
        ts_files = sorted(f for f in os.listdir(off_dir) if f.startswith('src_ts'))
        if ts_files:
            fails.append('flag off (default) still wrote src_ts* files: %r' % ts_files)

        # (b) flag on
        on_dir = os.path.join(tmp, 'on')
        _make_case(on_dir, env, 1)
        _run(on_dir, env)

        hdr_path = os.path.join(on_dir, 'src_ts_hdr0')
        bin_path = os.path.join(on_dir, 'src_ts0')
        if not (os.path.exists(hdr_path) and os.path.exists(bin_path)):
            print('FAIL test_fault_time_series_output')
            print('  - flag on did not write src_ts_hdr0/src_ts0; dir has: %r'
                 % sorted(os.listdir(on_dir)))
            return 1

        sys.path.insert(0, os.path.join(ROOT, 'scripts'))
        from readFaultTimeSeries import read_fault_ts, read_header

        meta, fields, xyz = read_header(hdr_path)
        nnode, nfield = meta['nnode'], meta['nfield']
        nbytes = os.path.getsize(bin_path)
        if nnode <= 0 or nfield <= 0:
            fails.append('header declares nnode=%r nfield=%r, expected both > 0' % (nnode, nfield))
        if nbytes % (8 * nnode * nfield) != 0:
            fails.append('src_ts0 is %d bytes, not a whole multiple of 8*nnode*nfield=%d'
                         % (nbytes, 8 * nnode * nfield))
        else:
            nrec = nbytes // (8 * nnode * nfield)
            if nrec != NSTEP:
                fails.append('src_ts0 holds %d complete step records, expected nstep=%d '
                             '(term = %d*dt)' % (nrec, NSTEP, NSTEP))

        if fails:
            print('FAIL test_fault_time_series_output')
            for f in fails:
                print('  -', f)
            return 1

        ts = read_fault_ts(hdr_path)
        if ts['xyz'].shape != (nnode, 3):
            fails.append('xyz coordinate block shape %r, expected (%d, 3)' % (ts['xyz'].shape, nnode))

        # cross-check against this case's own on-fault station file: default
        # test.tpv8 stations include (0.0, -7.5) km (strike, down-dip) --
        # scripts/defaultParameters.py st_coor_on_fault -- which case.setup
        # names faultst000dp075.txt.
        station_file = os.path.join(on_dir, 'faultst000dp075.txt')
        if not os.path.exists(station_file):
            fails.append('expected on-fault station file faultst000dp075.txt not found; '
                         'dir has: %r' % sorted(f for f in os.listdir(on_dir) if f.startswith('faultst')))
        else:
            xyz = ts['xyz']
            import numpy as np
            k = int(np.argmin(xyz[:, 0] ** 2 + (xyz[:, 2] + 7500.0) ** 2))
            node_x, node_z = xyz[k, 0], xyz[k, 2]
            if abs(node_x) > 1.0 or abs(node_z + 7500.0) > 1.0:
                fails.append('nearest src_ts node to station (0, -7.5 km) is at (%.3f, %.3f) m, '
                             'not within 1 m of the requested station' % (node_x, node_z))
            rows = [l.split() for l in open(station_file)
                   if l.strip() and not l.lstrip().startswith('#') and not l.strip()[0].isalpha()]
            st = np.array([[float(x) for x in row] for row in rows], dtype=float)
            n = min(len(st), ts['t'].size)
            if n == 0:
                fails.append('faultst000dp075.txt has no data rows to compare against')
            else:
                slip_strike = ts['slip_strike'][:n, k]
                traction_strike = ts['traction_strike'][:n, k] * 1e-6  # Pa -> MPa, station column 4
                d_slip = np.abs(slip_strike - st[:n, 1]).max()
                d_trac = np.abs(traction_strike - st[:n, 3]).max()
                # station text is E16.7 (7 significant digits); half a unit
                # of the last printed digit at this run's magnitude is a
                # generous absolute 1e-6 (slip, m) / 1e-3 (traction, MPa) --
                # loose enough to tolerate the text round-trip, tight enough
                # that a field swap (slip vs sliprate, or a missing 1e-6
                # Pa->MPa scale) fails by orders of magnitude.
                if not (d_slip <= 1e-6):
                    fails.append('src_ts slip_strike vs faultst000dp075.txt column 2 (h-slip): '
                                'max|diff| = %.3e m (tolerance 1e-6 m)' % d_slip)
                if not (d_trac <= 1e-3):
                    fails.append('src_ts traction_strike vs faultst000dp075.txt column 4 '
                                '(h-shear-stress): max|diff| = %.3e MPa (tolerance 1e-3 MPa)' % d_trac)

    if fails:
        print('FAIL test_fault_time_series_output')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_fault_time_series_output '
          '(flag off writes nothing; flag on shape + station cross-check agree, '
          '%d nodes x %d fields x %d steps)' % (nnode, nfield, NSTEP))
    return 0


if __name__ == '__main__':
    sys.exit(main())
