#! /usr/bin/env python3
"""
Regression guard (board row 167, deliverable 1): the post-meshing screen
(src/python/eqdyna/meshscreen.py + src/fortran/screenMesh.f90, exit codes
81-85) is GREEN on every registered case except the wedge-mesh cases
test.tpv36/test.tpv37, whose footwall cells are wired to hanging-wall nodes,
and each of its five checks (meshscreen.SCREEN_CODES) can actually fire.

Four tests:

  1. ALL CASES. `python3 -m eqdyna <case> --mesh-only` for every name in
     testNameList.nameList (the list is iterated itself, never a shorter hand
     copy). Every case exits 0 and prints `mesh-screen: PASS`, except
     WEDGE_EXPECTED, which must exit 81 with exactly the measured failing
     checks. That expectation is an explicit expected failure: when the
     wiring is fixed the wedge cases turn green and THIS test goes red with a
     message telling the fixer to update WEDGE_EXPECTED and the row 167 step.
  2. PLANTED DEFECTS. One defect per check, applied to the in-memory mesh of
     test.tpv8 (the exact arrays the run hands to screen_mesh); each must trip
     its own check, and the set it trips is compared with PLANT_EXPECTED. The
     unplanted mesh passes (control). A plant that trips more than its own
     check is recorded in PLANT_EXPECTED, not softened: the table is exact in
     both directions.
  3. FORTRAN/PYTHON PARITY. `mpirun -np 1 <binary> --mesh-only` on tpv8
     (green) and tpv36 (red): same verdict, same failing checks, same exit
     code as the Python backend. A missing binary is a failure, never a skip.
     The binary is src/fortran/eqdyna, built by `make` there if absent
     (bin/eqdyna is never touched).
  4. THRESHOLDS. SCREEN_* limits in meshscreen.py equal those in
     screenMesh.f90 (both parsed), and so do the check names and order.
"""
import concurrent.futures
import os
import re
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

import testNameList  # noqa: E402
from eqdyna import eqdyna3d, meshscreen  # noqa: E402

MESHSCREEN_PY = os.path.join(ROOT, 'src', 'python', 'eqdyna', 'meshscreen.py')
SCREENMESH_F90 = os.path.join(ROOT, 'src', 'fortran', 'screenMesh.f90')
FORTRAN_DIR = os.path.join(ROOT, 'src', 'fortran')
FORTRAN_BIN = os.path.join(FORTRAN_DIR, 'eqdyna')

# Measured 2026-11-23 (python mesh-only, screen commit 937ceb9):
# {case: (exit code, failing checks in SCREEN_CODES order)}.
WEDGE_EXPECTED = {
@@WEDGE@@
}
ROW167 = ('board row 167 (wedge-mesh mis-wiring, pathway_forward.md): the footwall brick is '
          'wired to hanging-wall nodes in the wedge mesh')

# plant name -> set of checks it trips on the planted tpv8 mesh (measured).
PLANT_EXPECTED = {
@@PLANT@@
}

COARSE_EXTRA = ''


def env():
    e = dict(os.environ)
    e['EQDYNAROOT'] = ROOT
    e['PATH'] = os.pathsep.join([os.path.join(ROOT, 'scripts'), e.get('PATH', '')])
    e['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    e['HWLOC_COMPONENTS'] = '-gl'
    e['OMP_NUM_THREADS'] = e['OPENBLAS_NUM_THREADS'] = '4'
    return e


def make_case(case, case_dir, extra=''):
    """create.newcase + serial decomposition (par.nx=ny=nz=1, as
    test_rank_local_mesh.py) + case.setup. Raises on any failure."""
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, case], env=env(), capture_output=True, text=True)
    assert r.returncode == 0, 'create.newcase %s failed: %s' % (case, (r.stderr or r.stdout)[-800:])
    p = os.path.join(case_dir, 'user_defined_params.py')
    with open(p) as f:
        text = f.read().rstrip('\n')
    with open(p, 'w') as f:
        f.write(text + '\n\n# test_mesh_screening.py: serial mesh-only input\n'
                'par.nx = 1\npar.ny = 1\npar.nz = 1\n' + extra)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env(),
                       capture_output=True, text=True)
    assert r.returncode == 0, 'case.setup %s failed: %s' % (case, (r.stderr or r.stdout)[-800:])


def parse_verdict(out, rc, who):
    """(verdict, checks tuple) from the `mesh-screen:` line; raises if absent."""
    lines = [l for l in out.splitlines() if l.startswith('mesh-screen:')]
    assert len(lines) == 1, '%s: expected exactly one mesh-screen line (rc %d), got %d; tail:\n%s' % (
        who, rc, len(lines), out[-1200:])
    m = re.match(r'mesh-screen: (PASS|FAIL) checks=(\S+)', lines[0])
    assert m, '%s: unparsable mesh-screen line: %r' % (who, lines[0])
    checks = () if m.group(2) == 'none' else tuple(m.group(2).split(','))
    return m.group(1), checks


def run_python_mesh_only(case_dir):
    r = subprocess.run([sys.executable, '-m', 'eqdyna', case_dir, '--mesh-only'], env=env(),
                       capture_output=True, text=True, timeout=1500)
    return r.returncode, r.stdout + r.stderr


def run_fortran_mesh_only(case_dir):
    r = subprocess.run(['mpirun', '-np', '1', FORTRAN_BIN, '--mesh-only'], cwd=case_dir, env=env(),
                       capture_output=True, text=True, timeout=1500)
    return r.returncode, r.stdout + r.stderr


def ordered(checks):
    return tuple(k for k in meshscreen.SCREEN_CODES if k in checks)


# ---------------------------------------------------------------------------
def test_all_cases(tmp):
    names = list(testNameList.nameList)
    assert names, 'testNameList.nameList is empty'
    assert set(WEDGE_EXPECTED) <= set(names), 'WEDGE_EXPECTED names a case not in nameList'
    dirs = {}
    for c in names:
        dirs[c] = os.path.join(tmp, 'all_' + c)
        make_case(c, dirs[c])
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs = {c: ex.submit(run_python_mesh_only, dirs[c]) for c in names}
        res = {c: f.result() for c, f in futs.items()}
    print('all-cases: %d cases meshed in %.0f s' % (len(names), time.time() - t0))
    errs = []
    for c in names:
        rc, out = res[c]
        verdict, checks = parse_verdict(out, rc, c)
        line = [l for l in out.splitlines() if l.startswith('mesh-screen:')][0][:300]
        if c in WEDGE_EXPECTED:
            code, want = WEDGE_EXPECTED[c]
            if verdict == 'PASS' and rc == 0:
                errs.append('%s: WEDGE EXPECTATION FLIPPED -- the screen is now GREEN on %s. The wedge '
                            'wiring has been fixed: remove %s from WEDGE_EXPECTED in this file and '
                            'close the matching step of %s.' % (c, c, c, ROW167))
            elif (rc, ordered(checks)) != (code, want):
                errs.append('%s: expected FAIL rc=%d checks=%s (%s), got %s rc=%d checks=%s; %s' % (
                    c, code, ','.join(want), ROW167, verdict, rc, ','.join(checks), line))
            else:
                print('  %-18s RED as expected  rc=%d checks=%s' % (c, rc, ','.join(checks)))
        else:
            if not (verdict == 'PASS' and rc == 0):
                errs.append('%s: must be GREEN, got %s rc=%d checks=%s; %s' % (
                    c, verdict, rc, ','.join(checks), line))
            else:
                print('  %-18s GREEN' % c)
    assert not errs, 'all-cases failed:\n' + '\n'.join(errs)


# ---------------------------------------------------------------------------
class Mesh:
    """The arrays build_solver_state hands to screen_mesh, captured in memory."""

    def __init__(self, case_dir):
        real = meshscreen.screen_mesh
        got = {}

        def grab(*a):
            got['a'] = a
            return real(*a)
        meshscreen.screen_mesh = grab
        try:
            eqdyna3d.build_solver_state(case_dir, mesh_only=True)
        finally:
            meshscreen.screen_mesh = real
        self.a = got['a']

    def copy(self):
        a = list(self.a)
        for i in range(8):
            a[i] = np.array(a[i], copy=True)
        return a


def screen_set(args):
    """(set of failing check names, messages) for an argument list; never raises on a bad mesh."""
    m, fail = meshscreen._screen_measure(*args)
    return set(fail), fail


def interior_pair(a):
    """0-based pair row nearest the centre of the fault patch."""
    X, nsmp = a[0], a[3]
    xm, zm = X[nsmp[:, 1], 0], X[nsmp[:, 1], 2]
    xc, zc = 0.5 * (xm.min() + xm.max()), 0.5 * (zm.min() + zm.max())
    return int(np.argmin(np.abs(xm - xc) + np.abs(zm - zc)))


def far_node(a):
    """A node at the +y box face (far from the fault), 1-indexed."""
    return int(np.argmax(a[0][1:, 1]) + 1)


def plant(a, name):
    """Apply ONE defect to the argument list a (in place)."""
    X, conn, et, nsmp, un, us, ud, arn = a[:8]
    p = interior_pair(a)
    s, m = int(nsmp[p, 0]), int(nsmp[p, 1])
    if name == 'sides':
        # one hanging-wall cell holds the pair's SLAVE node in place of its master
        e = int(np.nonzero((conn == m).any(axis=1))[0][0])
        conn[e][conn[e] == m] = s
    elif name == 'bridging':
        # one cell holds BOTH nodes of the pair: a corner of it off the fault plane becomes the slave
        e = int(np.nonzero((conn == m).any(axis=1))[0][0])
        k = int(np.nonzero(np.abs(X[conn[e], 1] - X[m, 1]) > 0)[0][0])
        conn[e][k] = s
    elif name == 'fault_area':
        arn[p + 1] *= 1.001
    elif name == 'normals':
        # rotate the (un, us) pair of one node by 0.01 rad about ud: still orthonormal, no longer
        # perpendicular to the meshed facets
        c, sn = np.cos(0.01), np.sin(0.01)
        n0, s0 = un[p + 1].copy(), us[p + 1].copy()
        un[p + 1] = c * n0 + sn * s0
        us[p + 1] = -sn * n0 + c * s0
    elif name == 'normals_flip':
        un[p + 1] = -un[p + 1]
    elif name == 'volume_mass_coord':
        # a +y-face node dragged through its cell (inverted Jacobian)
        far = far_node(a)
        e = int(np.nonzero((conn == far).any(axis=1))[0][0])
        X[far] = X[far] + 1.5 * (X[conn[e]].mean(axis=0) - X[far])
    elif name == 'volume_mass_box':
        # a +y-face node pushed outward: cell volumes no longer sum to the rank box
        X[far_node(a), 1] += 0.1 * (X[1:, 1].max() - X[1:, 1].min()) / 100.0
    else:
        raise AssertionError('unknown plant ' + name)


PLANTS = ('sides', 'bridging', 'fault_area', 'normals', 'normals_flip',
          'volume_mass_coord', 'volume_mass_box')
OWN = {'sides': 'sides', 'bridging': 'bridging', 'fault_area': 'fault_area', 'normals': 'normals',
       'normals_flip': 'normals', 'volume_mass_coord': 'volume_mass', 'volume_mass_box': 'volume_mass'}


def test_planted_defects(tmp):
    d = os.path.join(tmp, 'plant_tpv8')
    make_case('test.tpv8', d, COARSE_EXTRA)
    mesh = Mesh(d)
    got, fail = screen_set(mesh.copy())
    assert got == set(), 'control (unplanted tpv8) must pass, tripped %s: %s' % (sorted(got), fail)
    errs = []
    for name in PLANTS:
        a = mesh.copy()
        plant(a, name)
        got, fail = screen_set(a)
        want = set(PLANT_EXPECTED[name])
        print('  plant %-18s trips %-34s expected %s' % (name, sorted(got), sorted(want)))
        if got != want:
            errs.append('plant %s tripped %s, expected %s; messages %s' % (name, sorted(got), sorted(want), fail))
        if OWN[name] not in got:
            errs.append('plant %s did not trip its own check %s' % (name, OWN[name]))
    assert not errs, 'planted defects:\n' + '\n'.join(errs)


# ---------------------------------------------------------------------------
def fortran_binary():
    if not os.path.exists(FORTRAN_BIN):
        r = subprocess.run(['make'], cwd=FORTRAN_DIR, env=dict(env(), MACHINE='ubuntu'),
                           capture_output=True, text=True, timeout=600)
        assert r.returncode == 0 and os.path.exists(FORTRAN_BIN), (
            'Fortran binary %s missing and `make` in src/fortran failed (rc %d):\n%s' % (
                FORTRAN_BIN, r.returncode, (r.stdout + r.stderr)[-1500:]))
    return FORTRAN_BIN


def test_fortran_python_parity(tmp):
    fortran_binary()
    errs = []
    for case in ('test.tpv8', 'test.tpv36'):
        d = os.path.join(tmp, 'parity_' + case)
        make_case(case, d)
        prc, pout = run_python_mesh_only(d)
        frc, fout = run_fortran_mesh_only(d)
        pv = parse_verdict(pout, prc, 'python ' + case)
        fv = parse_verdict(fout, frc, 'fortran ' + case)
        print('  %-12s python rc=%d %s | fortran rc=%d %s' % (case, prc, pv, frc, fv))
        want_red = case in WEDGE_EXPECTED
        if (pv[0] == 'FAIL') != want_red:
            errs.append('%s: python verdict %s, expected %s' % (case, pv[0], 'FAIL' if want_red else 'PASS'))
        if ordered(fv[1]) != ordered(pv[1]) or fv[0] != pv[0] or frc != prc:
            errs.append('%s: PARITY BROKEN python rc=%d %s checks=%s vs fortran rc=%d %s checks=%s' % (
                case, prc, pv[0], ','.join(pv[1]), frc, fv[0], ','.join(fv[1])))
    assert not errs, 'parity:\n' + '\n'.join(errs)


# ---------------------------------------------------------------------------
def test_thresholds_equal():
    with open(MESHSCREEN_PY) as f:
        py = dict((k, float(v)) for k, v in re.findall(
            r'^(SCREEN_[A-Z_]+)\s*=\s*([0-9.eE+-]+)', f.read(), re.M))
    with open(SCREENMESH_F90) as f:
        txt = f.read()
    fo = dict((k, float(v.lower().replace('d', 'e'))) for k, v in re.findall(
        r'parameter\s*::\s*(SCREEN_[A-Z_]+)\s*=\s*([0-9.dD+-]+)', txt))
    assert len(py) == 7 and len(fo) == 7, 'expected 7 limits each, python %s fortran %s' % (sorted(py), sorted(fo))
    assert py == fo, 'screen limits differ: python %s fortran %s' % (py, fo)
    names = [s.strip() for s in re.search(r"CNAME\(5\)\s*=\s*\(/(.*?)/\)", txt).group(1).replace("'", '').split(',')]
    assert list(meshscreen.SCREEN_CODES) == names, (
        'check names/order differ: python %s fortran %s' % (list(meshscreen.SCREEN_CODES), names))


def main():
    tmp = tempfile.mkdtemp(prefix='mesh_screening_')
    tests = [('thresholds_equal', lambda: test_thresholds_equal()),
             ('planted_defects', lambda: test_planted_defects(tmp)),
             ('fortran_python_parity', lambda: test_fortran_python_parity(tmp)),
             ('all_cases', lambda: test_all_cases(tmp))]
    sel = sys.argv[1:]
    bad = 0
    for name, fn in tests:
        if sel and name not in sel:
            continue
        t0 = time.time()
        try:
            fn()
            print('SUCCESS %s (%.1f s)' % (name, time.time() - t0))
        except AssertionError as exc:
            bad += 1
            print('FAIL %s (%.1f s): %s' % (name, time.time() - t0, exc))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
