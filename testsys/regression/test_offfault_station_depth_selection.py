#! /usr/bin/env python3
"""
Regression guard (row 94 audit finding 3): the off-fault DEPTH-CLAMP
algorithm itself -- nearest-node snap, the exact-tie tie-break, an
out-of-band drop, and Fortran/Python agreement on the same station set.
`test_offfault_station_dropped_report.py`/`test_row114_station_output.py`
pin the REPORT wording with synthetic (matched, actual) inputs; neither
exercises the real nearest-node SEARCH in `setSurfaceStation`/
`build_station_matching`, so a bug in the search itself (wrong tie-break,
wrong PML-band cutoff) could pass both while this algorithm is wrong.

METHOD: a real, tiny (forced serial, 3-step) test.tpv8 case with 4
off-fault stations chosen from that CASE's OWN z-grid (not hand-derived --
see the constants' provenance below), each isolated to a DIFFERENT
physical node so the "one station per node, first-match-wins" contention
in setSurfaceStation (a real, pre-existing property, unrelated to this
mission) cannot mask the result:
  station 1 (x=0, y=1.0, z=-12.0 km): sits exactly on a grid z-plane ->
      no snap, matches at its own requested depth.
  station 2 (z=-13.9 km): 100 m off the nearest plane (z=-14.0 km,
      400 m from its other neighbour z=-13.5 km) -> unambiguous SNAP to
      z=-14.0 km.
  station 3 (z=-10.25 km): the EXACT midpoint between z=-10.5 km and
      z=-10.0 km (250 m each) -> a TRUE TIE. The owner-documented rule
      (meshgen.f90's setSurfaceStation comment) is "keeps the lower index
      (deeper/more negative z)" -- this asserts the actual result (-10.5,
      not -10.0), not just "some node".
  station 4 (z=-25.0 km): below PMLb(5) (this case's z-grid puts it at
      about -22.07 km) but still inside the padded domain (zline[0] is
      about -26.47 km) -- finding 1's physical-band clamp, not a
      request that is simply off the mesh entirely -- so it must DROP
      with finding 6's "checked cause" wording, not snap into the PML.

These z-values were derived by running this SAME case's grid-builder
once (`meshgen.build_grid_lines`, Python) and reading its own `zline`/
`pmlb['zmin0']` -- not invented -- so this test's literals ARE that
case's real z-plane values. x=0, y=1.0 km sits exactly on this case's x/y
grid too (an x/y snap is pre-existing, unrelated, and orthogonal to this
mission -- keeping x/y exact isolates the z behaviour under test).

Runs the REAL `bin/eqdyna` (skips loudly, never silently, if absent) AND
calls the REAL Python `meshgen.build_station_matching` on the identical
case directory, and requires them to report the SAME matched depth for
every station that matches on either side (finding 3's 4th requirement:
Fortran and Python choose the same node for the same station set).

Cheap (rule 9): forced single-rank, term = 3*dt (a handful of ms of
simulated time) -- mesh generation and the coverage/report subroutines run
regardless of term (eqdyna3d.f90's call order: meshgen ->
checkOffFaultStationCoverage happens before any time stepping), so this is
dominated by mesh setup, not dynamics.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

CASE = 'test.tpv8'
STATIONS_KM = [[0.0, 1.0, -12.0], [0.0, 1.0, -13.9], [0.0, 1.0, -10.25], [0.0, 1.0, -25.0]]


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'), os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(case_dir, env):
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=env, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed: %s' % (r.stderr or '')[-800:])
    params_path = os.path.join(case_dir, 'user_defined_params.py')
    with open(params_path) as f:
        text = f.read().rstrip('\n')
    with open(params_path, 'w') as f:
        f.write(text + '\n\n# row94 finding-3 probe: forced serial, tiny term, fixed stations\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n'
                       'par.st_coor_off_fault = %r\n'
                       'par.n_off_fault = len(par.st_coor_off_fault)\n'
                       'par.term = par.dt * 3\n' % STATIONS_KM)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise AssertionError('case.setup failed: %s' % (r.stderr or '')[-800:])


def main():
    if shutil.which('mpirun') is None:
        print('FAIL test_offfault_station_depth_selection: mpirun not on PATH '
              '(this guard requires a real MPI run, it does not skip silently)')
        return 1
    env = _env()
    binary = os.path.join(ROOT, 'bin', 'eqdyna')
    if not os.path.exists(binary):
        print('FAIL test_offfault_station_depth_selection: %s not built '
              '(run install-eqdyna.sh first -- this guard does not skip silently)' % binary)
        return 1

    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = os.path.join(tmp, 'case')
        try:
            _make_case(case_dir, env)
        except AssertionError as e:
            print('FAIL test_offfault_station_depth_selection')
            print(' -', e)
            return 1
        r = subprocess.run(['mpirun', '-np', '1', 'eqdyna'], cwd=case_dir, env=env,
                           capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            print('FAIL test_offfault_station_depth_selection')
            print(' - eqdyna exited %d:\n%s' % (r.returncode, (r.stdout + r.stderr)[-2000:]))
            return 1
        out = r.stdout

        if 'station 1' in out:
            fails.append('station 1 (exact grid-plane depth) mentioned in the drop/snap '
                         'report; stdout tail:\n%s' % out[-1500:])
        body_files = sorted(f for f in os.listdir(case_dir) if f.startswith('body'))
        if 'body010st000dp120.txt' not in body_files:
            fails.append('station 1 (z=-12.0 km exact) did not write body010st000dp120.txt; '
                         'wrote %r' % body_files)

        m = re.search(r'snapped off-fault station 2 .*?actual x,y,z =\s*\S+\s+\S+\s+(\S+)\s+km, '
                      r'distance =\s*(\S+)\s+km', out)
        if not m or [float(v) for v in m.groups()] != [-14.0, 0.1]:
            fails.append('station 2 snap not reported as actual z=-14.0 km, distance=0.100 km; '
                         'stdout tail:\n%s' % out[-1500:])
        if 'body010st000dp139.txt' not in body_files:
            fails.append('station 2 did not write body010st000dp139.txt (snap to -14.0 km); '
                         'wrote %r' % body_files)

        m = re.search(r'snapped off-fault station 3 .*?actual x,y,z =\s*\S+\s+\S+\s+(\S+)\s+km, '
                      r'distance =\s*(\S+)\s+km', out)
        if not m or [float(v) for v in m.groups()] != [-10.5, 0.25]:
            fails.append('station 3 (exact tie) not resolved to the deeper plane z=-10.5 km, '
                         'distance=0.250 km; stdout tail:\n%s' % out[-1500:])
        if 'body010st000dp103.txt' not in body_files:
            fails.append('station 3 (tie) did not write body010st000dp103.txt (deeper plane, '
                         '-10.5 km); wrote %r' % body_files)
        if 'body010st000dp100.txt' in body_files:
            fails.append('station 3 (tie) wrote body010st000dp100.txt -- resolved to the '
                         'SHALLOWER plane (-10.0 km), not the documented deeper one')

        if not re.search(r'dropped off-fault station 4 at x,y,z =\s*\S+\s+\S+\s+-25\.000 km '
                         r'\(checked cause: requested depth is outside the physical, non-PML '
                         r'mesh band\)', out):
            fails.append('station 4 (out-of-band depth) not reported as a checked-cause drop; '
                         'stdout tail:\n%s' % out[-1500:])
        if any('dp250' in f for f in body_files) or any('dp220' in f for f in body_files):
            fails.append('station 4 (out-of-band) wrote a body* file -- it must DROP, not '
                         'snap into the PML; wrote %r' % body_files)

        from eqdyna import readInputFiles, meshgen
        params, _g = readInputFiles.build_params(case_dir)
        xonfs, x4nds = readInputFiles.read_bstations(os.path.join(case_dir, 'bStations.txt'))
        xline, yline, zline, pmlb, bounds = meshgen.build_grid_lines(params)
        _anonfs, off_matches, z_valid = meshgen.build_station_matching(
            xline, yline, zline, params, xonfs, x4nds, pmlb['zmin0'], bounds[2][1])
        meshCoor, _nftnd, _nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)
        py_matched_z = {sc: float(meshCoor[nc][2]) for sc, nc in off_matches}
        want_z = {1: -12000.0, 2: -14000.0, 3: -10500.0}
        if py_matched_z != want_z:
            fails.append('python build_station_matching matched z %r, expected %r '
                         '(must agree with the Fortran run above)' % (py_matched_z, want_z))
        if 4 in py_matched_z or z_valid[3]:
            fails.append('python matched or validated station 4 (out-of-band) -- '
                         'z_valid=%r, matched=%r' % (z_valid.tolist(), py_matched_z))

    if fails:
        print('FAIL test_offfault_station_depth_selection')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_offfault_station_depth_selection '
          '(exact/snap/tie/drop all correct, Fortran and Python agree)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
