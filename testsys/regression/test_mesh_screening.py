#! /usr/bin/env python3
"""
Regression guard (board row 167, deliverable 1): the post-meshing screen
(src/python/eqdyna/meshscreen.py + src/fortran/screenMesh.f90, exit codes
81-85) is GREEN on every registered case except the wedge-mesh cases
test.tpv36/test.tpv37, whose footwall cells are wired to hanging-wall nodes,
and each of its five checks (meshscreen.SCREEN_CODES) can actually fire.

Tests (every one a hard failure, none skips):

  1. ALL CASES, both backends. `python3 -m eqdyna <case> --mesh-only` and
     `mpirun -np 1 eqdyna --mesh-only` for every name in
     testNameList.nameList (the list is iterated itself, never a shorter hand
     copy). Every case exits 0 and prints `mesh-screen: PASS`, except the
     cases in WEDGE_EXPECTED_RED, which must exit with exactly the listed
     code and failing checks. Both backends must agree on verdict, exit code,
     failing checks and the element/node/fault-pair counts; every verdict
     line is printed with those counts and the number of (cell,node)
     attachments tested, and a verdict over zero of them is refused (a screen
     that compared nothing proves nothing).
     WEDGE_EXPECTED_RED is THE list step 2 of row 167 edits: when the wiring
     is fixed the wedge cases turn green and this test goes red with a message
     naming the list; deleting the entry (and the matching UNSUPPORTED cells in
     testsys/matrix.py) is the whole edit.
  2. JUDGE BOTH WAYS (rule 14a). The function that turns a run into a verdict
     (judge) is fed a green wedge case, a wedge case failing the wrong check,
     a red non-wedge case and a green non-wedge case; it must flag the first
     three and pass the fourth, and the same for the parsed measurement.
  3. PLANTED DEFECTS, Python. One defect per check, applied to the in-memory
     mesh of test.tpv8 (the exact arrays the run hands to screen_mesh); each
     trips its own check, the tripped set is compared with PLANT_EXPECTED
     exactly (a plant that trips two checks is recorded as two, not softened),
     and the unplanted mesh is the control.
  4. PLANTED DEFECTS, Fortran. The same planted meshes are dumped to a raw
     stream file and fed to testsys/parity/screen_mesh_dump_driver.f90, which
     links the PRODUCTION screenMesh.o (no hook in the production binary); the
     process exit status must equal the code of the first expected check.
  5. 2-RANK FORTRAN. `mpirun -np 2` on tpv8 (green) and tpv36 (red, 81): the
     screen runs per rank on its own box.
  6. THRESHOLDS. SCREEN_* limits in meshscreen.py equal those in
     screenMesh.f90 (both parsed), and so do the check names and order.

mpirun on this class of box needs `env -u DISPLAY HWLOC_COMPONENTS=-gl,...` (a
dead X server blocks hwloc GL discovery forever) and is always wrapped in
`timeout -s KILL`.
"""
import concurrent.futures
import os
import re
import shutil
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
DRIVER_F90 = os.path.join(ROOT, 'testsys', 'parity', 'screen_mesh_dump_driver.f90')
FORTRAN_DIR = os.path.join(ROOT, 'src', 'fortran')
FORTRAN_BIN = os.path.join(FORTRAN_DIR, 'eqdyna')
MPIRUN_TIMEOUT_S = 1500

# THE expected-red list (row 167 step 1). {case: (exit code, failing checks in
# SCREEN_CODES order)}. Measured 2026-11-24 on both backends: the wedge mesh
# wires 120 footwall (master-attached) attachments to the wrong side, and the
# master/slave volume share ratio of 59 fault-node pairs exceeds 2 (peak 2.33).
# Step 2 of row 167 deletes these entries when the wiring is fixed.
WEDGE_EXPECTED_RED = {
    'test.tpv36': (81, ('sides', 'volume_mass')),
    'test.tpv37': (81, ('sides', 'volume_mass')),
}
ROW167 = ('board row 167 (wedge-mesh mis-wiring, pathway_forward.md): the footwall brick is '
          'wired to hanging-wall nodes in the wedge mesh')

# plant name -> set of checks it trips on the planted tpv8 mesh (measured).
PLANT_EXPECTED = {
    'sides': {'sides'},
    'sides_swap_node': {'sides', 'fault_area'},
    'bridging': {'bridging'},
    'bridging_cell': {'bridging', 'sides', 'volume_mass'},
    'fault_area': {'fault_area'},
    'normals': {'normals'},
    'volume_mass_coord': {'volume_mass'},
    'volume_mass_conn': {'volume_mass'},
}

COARSE_EXTRA = ''


def env():
    e = dict(os.environ)
    e['EQDYNAROOT'] = ROOT
    e['PATH'] = os.pathsep.join([os.path.join(ROOT, 'scripts'), e.get('PATH', '')])
    e['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    e.pop('DISPLAY', None)
    e['HWLOC_COMPONENTS'] = '-gl,-opencl,-cuda,-nvml'
    e['OMP_NUM_THREADS'] = e['OPENBLAS_NUM_THREADS'] = e['MKL_NUM_THREADS'] = '1'
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
    """(verdict, checks tuple, measured dict, first raw line) from the
    `mesh-screen:` line(s); raises if there is none. One line per rank: the
    verdict is FAIL if any rank failed, the checks are the union, and the four
    COUNT_KEYS are summed over ranks (the other fields are the first rank's)."""
    lines = [l for l in out.splitlines() if l.startswith('mesh-screen:')]
    assert lines, '%s: no mesh-screen line (rc %d); tail:\n%s' % (who, rc, out[-1200:])
    verdict, checks, meas = 'PASS', set(), None
    for l in lines:
        m = re.match(r'mesh-screen: (PASS|FAIL) checks=(\S+) (.*)$', l)
        assert m, '%s: unparsable mesh-screen line: %r' % (who, l)
        if m.group(1) == 'FAIL':
            verdict = 'FAIL'
        if m.group(2) != 'none':
            checks |= set(m.group(2).split(','))
        kv = dict(x.split('=') for x in m.group(3).split())
        if meas is None:
            meas = kv
        else:
            for k in COUNT_KEYS:
                meas[k] = str(int(meas[k]) + int(kv[k]))
    return verdict, ordered(checks), meas, lines[0]


def run_python_mesh_only(case_dir):
    r = subprocess.run([sys.executable, '-m', 'eqdyna', case_dir, '--mesh-only'], env=env(),
                       capture_output=True, text=True, timeout=3000)
    return r.returncode, r.stdout + r.stderr


def run_fortran_mesh_only(case_dir, ranks=1):
    r = subprocess.run(['timeout', '-s', 'KILL', str(MPIRUN_TIMEOUT_S), 'mpirun', '-np', str(ranks),
                        FORTRAN_BIN, '--mesh-only'], cwd=case_dir, env=env(),
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def ordered(checks):
    return tuple(k for k in meshscreen.SCREEN_CODES if k in checks)


COUNT_KEYS = ('elements', 'nodes', 'fault_pairs', 'sides_attachments')


def judge(case, rc, verdict, checks, meas):
    """Error string for one backend's result on one case, or None. The ONLY
    place a run becomes a verdict; test_judge_both_ways feeds it known-bad
    inputs."""
    counts = {k: int(meas[k]) for k in COUNT_KEYS}
    empty = [k for k, v in counts.items() if v <= 0]
    if empty:
        return '%s: screen compared nothing (%s == 0): a vacuous verdict' % (case, ','.join(empty))
    if case in WEDGE_EXPECTED_RED:
        code, want = WEDGE_EXPECTED_RED[case]
        if verdict == 'PASS' and rc == 0:
            return ('%s: WEDGE EXPECTATION FLIPPED -- the screen is now GREEN on %s. The wiring has been '
                    'fixed: delete %s from WEDGE_EXPECTED_RED in testsys/regression/test_mesh_screening.py '
                    '(and its UNSUPPORTED cells in testsys/matrix.py) and close the matching step of %s.'
                    % (case, case, case, ROW167))
        if (rc, ordered(checks)) != (code, want):
            return '%s: expected FAIL rc=%d checks=%s (%s), got %s rc=%d checks=%s' % (
                case, code, ','.join(want), ROW167, verdict, rc, ','.join(checks))
        return None
    if not (verdict == 'PASS' and rc == 0):
        return '%s: must be GREEN, got %s rc=%d checks=%s' % (case, verdict, rc, ','.join(checks))
    return None


def test_judge_both_ways():
    meas = {'elements': '10', 'nodes': '20', 'fault_pairs': '3', 'sides_attachments': '12'}
    wedge = sorted(WEDGE_EXPECTED_RED)[0]
    code, want = WEDGE_EXPECTED_RED[wedge]
    green = [c for c in testNameList.nameList if c not in WEDGE_EXPECTED_RED][0]
    assert judge(wedge, code, 'FAIL', want, meas) is None, 'expected-red result must be accepted'
    assert judge(green, 0, 'PASS', (), meas) is None, 'green result on a green case must be accepted'
    bad = {
        'green wedge case': judge(wedge, 0, 'PASS', (), meas),
        'wedge case, wrong check': judge(wedge, 82, 'FAIL', ('bridging',), meas),
        'wedge case, wrong code': judge(wedge, 85, 'FAIL', want, meas),
        'red non-wedge case': judge(green, 81, 'FAIL', ('sides',), meas),
        'vacuous (no attachments)': judge(green, 0, 'PASS', (), dict(meas, sides_attachments='0')),
        'vacuous (no fault pairs)': judge(green, 0, 'PASS', (), dict(meas, fault_pairs='0')),
    }
    missed = [k for k, v in bad.items() if v is None]
    assert not missed, 'judge failed to flag: %s' % missed
    assert 'FLIPPED' in bad['green wedge case'] and 'WEDGE_EXPECTED_RED' in bad['green wedge case']


# ---------------------------------------------------------------------------
def fortran_binary():
    if not os.path.exists(FORTRAN_BIN):
        r = subprocess.run(['make'], cwd=FORTRAN_DIR, env=dict(env(), MACHINE='ubuntu'),
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0 and os.path.exists(FORTRAN_BIN), (
            'Fortran binary %s missing and `make` in src/fortran failed (rc %d):\n%s' % (
                FORTRAN_BIN, r.returncode, (r.stdout + r.stderr)[-1500:]))
    return FORTRAN_BIN


def test_all_cases(tmp):
    names = list(testNameList.nameList)
    assert names, 'testNameList.nameList is empty'
    assert set(WEDGE_EXPECTED_RED) <= set(names), 'WEDGE_EXPECTED_RED names a case not in nameList'
    fortran_binary()
    dirs = {}
    for c in names:
        dirs[c] = os.path.join(tmp, 'all_' + c)
        make_case(c, dirs[c])
    t0 = time.time()

    def both(c):
        return run_python_mesh_only(dirs[c]), run_fortran_mesh_only(dirs[c])
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        futs = {c: ex.submit(both, c) for c in names}
        res = {c: f.result() for c, f in futs.items()}
    print('all-cases: %d cases, python + fortran, in %.0f s' % (len(names), time.time() - t0))
    errs = []
    for c in names:
        (prc, pout), (frc, fout) = res[c]
        pv = parse_verdict(pout, prc, 'python ' + c)
        fv = parse_verdict(fout, frc, 'fortran ' + c)
        mp = pv[2]
        print('  %-18s %s rc=%-2d checks=%-18s elements=%s nodes=%s fault_pairs=%s attachments=%s | fortran %s rc=%d' % (
            c, pv[0], prc, ','.join(pv[1]) or 'none', mp['elements'], mp['nodes'], mp['fault_pairs'],
            mp['sides_attachments'], fv[0], frc))
        for who, rc, v in (('python', prc, pv), ('fortran', frc, fv)):
            e = judge(c, rc, v[0], v[1], v[2])
            if e:
                errs.append('[%s] %s; %s' % (who, e, v[3][:300]))
        if (pv[0], pv[1], prc) != (fv[0], fv[1], frc):
            errs.append('%s: PARITY BROKEN python %s rc=%d %s vs fortran %s rc=%d %s' % (
                c, pv[0], prc, ','.join(pv[1]), fv[0], frc, ','.join(fv[1])))
        for k in COUNT_KEYS[:3] + ('det_bad_cells', 'sides_wrong_master', 'sides_wrong_slave', 'bridging_cells',
                                   'area_nofacet_pairs', 'area_badcount_pairs', 'share_empty_pairs',
                                   'mass_ratio_out_pairs', 'sides_attachments'):
            if pv[2][k] != fv[2][k]:
                errs.append('%s: count %s differs, python %s fortran %s' % (c, k, pv[2][k], fv[2][k]))
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
    """Apply ONE defect to the argument list a (in place). The plain-named
    plant of each check is the single-cause one the PLANT_EXPECTED table says
    trips only that check; the *_swap_node / *_cell / *_conn plants are the
    realistic wiring defects, which unavoidably disturb a second quantity and
    are recorded as exactly the set they trip."""
    X, conn, et, nsmp, un, us, ud, arn = a[:8]
    p = interior_pair(a)
    s, m = int(nsmp[p, 0]), int(nsmp[p, 1])
    if name == 'sides':
        # the pair's frame points the other way: every attached cell is on the wrong side
        un[p + 1] = -un[p + 1]
    elif name == 'sides_swap_node':
        # one hanging-wall cell holds the pair's SLAVE node in place of its master (coincident, so
        # no geometry changes; the cell loses a master facet)
        e = int(np.nonzero((conn == m).any(axis=1))[0][0])
        conn[e][conn[e] == m] = s
    elif name == 'bridging':
        # the slave node is 1 m from its master in strike (a sheared vertex: volumes and the
        # master-side facets are untouched)
        X[s, 0] += 1.0
    elif name == 'bridging_cell':
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
    elif name == 'volume_mass_coord':
        # a +y-face node dragged through its cell (inverted Jacobian)
        far = far_node(a)
        e = int(np.nonzero((conn == far).any(axis=1))[0][0])
        X[far] = X[far] + 1.5 * (X[conn[e]].mean(axis=0) - X[far])
    elif name == 'volume_mass_conn':
        # a far cell's corner wired to its neighbour corner's node: the cell collapses, volumes no
        # longer sum to the rank box
        far = far_node(a)
        e = int(np.nonzero(conn == far)[0][0])
        k = int(np.nonzero(conn[e] == far)[0][0])
        base = 4 * int(k > 3)
        conn[e][base + (k + 2) % 4] = conn[e][base + (k + 1) % 4]
    else:
        raise AssertionError('unknown plant ' + name)


PLANTS = ('sides', 'sides_swap_node', 'bridging', 'bridging_cell', 'fault_area', 'normals',
          'volume_mass_coord', 'volume_mass_conn')
OWN = {'sides': 'sides', 'sides_swap_node': 'sides', 'bridging': 'bridging', 'bridging_cell': 'bridging',
       'fault_area': 'fault_area', 'normals': 'normals',
       'volume_mass_coord': 'volume_mass', 'volume_mass_conn': 'volume_mass'}


def dump_mesh(a, path):
    """Raw stream dump in the layout screen_mesh_dump_driver.f90 reads."""
    X, conn, et, nsmp, un, us, ud, arn, lines, tol, planar = a
    xline, yline, zline = (np.asarray(v, dtype='<f8') for v in lines)
    nft = nsmp.shape[0]
    with open(path, 'wb') as f:
        f.write(np.array([X.shape[0] - 1, conn.shape[0]], dtype='<i8').tobytes())
        f.write(np.array([len(xline), len(yline), len(zline)], dtype='<i4').tobytes())
        f.write(np.array([nft], dtype='<i8').tobytes())
        f.write(np.array([0 if planar[0] else 1], dtype='<i4').tobytes())
        f.write(np.array([tol], dtype='<f8').tobytes())
        for v in (xline, yline, zline):
            f.write(v.tobytes())
        f.write(np.ascontiguousarray(X[1:], dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(conn, dtype='<i8').tobytes())
        f.write(np.ascontiguousarray(et, dtype='<i4').tobytes())
        f.write(np.ascontiguousarray(nsmp[:, :2], dtype='<i8').tobytes())
        for v in (un, us, ud):
            f.write(np.ascontiguousarray(v[1:nft + 1], dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(arn[1:nft + 1], dtype='<f8').tobytes())


def build_driver(tmp):
    fortran_binary()          # builds globalvar.o / errorCodes.o / screenMesh.o / .mod files in src/fortran
    exe = os.path.join(tmp, 'screen_mesh_dump_driver')
    r = subprocess.run(['mpif90', '-fopenmp', '-ffree-line-length-none', '-O1', '-I' + FORTRAN_DIR, '-J' + tmp,
                        DRIVER_F90, os.path.join(FORTRAN_DIR, 'globalvar.o'),
                        os.path.join(FORTRAN_DIR, 'errorCodes.o'), os.path.join(FORTRAN_DIR, 'screenMesh.o'),
                        '-o', exe], capture_output=True, text=True, timeout=600, env=env())
    assert r.returncode == 0, 'driver build failed:\n' + (r.stdout + r.stderr)[-1500:]
    return exe


def mesh_args_ok(a):
    nsmp = a[3]
    assert a[1].shape[0] > 0 and a[0].shape[0] > 1 and nsmp.shape[0] > 0, 'empty mesh handed to the screen'


def test_planted_defects(tmp):
    d = os.path.join(tmp, 'plant_tpv8')
    make_case('test.tpv8', d, COARSE_EXTRA)
    mesh = Mesh(d)
    mesh_args_ok(mesh.a)
    exe = build_driver(tmp)
    # control: the unplanted mesh passes in BOTH implementations, driver exit 0
    got, fail = screen_set(mesh.copy())
    assert got == set(), 'control (unplanted tpv8) must pass, tripped %s: %s' % (sorted(got), fail)
    cpath = os.path.join(tmp, 'control.dump')
    dump_mesh(mesh.copy(), cpath)
    r = subprocess.run([exe, cpath], capture_output=True, text=True, env=env(), timeout=600)
    assert r.returncode == 0 and 'returned (no abort)' in r.stdout, (
        'Fortran control (unplanted tpv8 dump) must pass, rc %d:\n%s' % (r.returncode, r.stdout[-800:]))
    errs = []
    for name in PLANTS:
        a = mesh.copy()
        plant(a, name)
        got, fail = screen_set(a)
        want = set(PLANT_EXPECTED[name])
        path = os.path.join(tmp, name + '.dump')
        dump_mesh(a, path)
        r = subprocess.run([exe, path], capture_output=True, text=True, env=env(), timeout=600)
        want_code = meshscreen.SCREEN_CODES[ordered(want)[0]]
        print('  plant %-18s python trips %-14s expected %-14s | fortran rc=%d expected %d' % (
            name, ','.join(sorted(got)), ','.join(sorted(want)), r.returncode, want_code))
        if got != want:
            errs.append('plant %s tripped %s, expected %s; messages %s' % (name, sorted(got), sorted(want), fail))
        if OWN[name] not in got:
            errs.append('plant %s did not trip its own check %s' % (name, OWN[name]))
        if r.returncode != want_code:
            errs.append('plant %s: Fortran driver rc %d, expected %d; tail %s' % (
                name, r.returncode, want_code, r.stdout[-500:]))
        else:
            fchecks = re.search(r'mesh-screen: FAIL checks=(\S+)', r.stdout)
            if not fchecks or ordered(fchecks.group(1).split(',')) != ordered(want):
                errs.append('plant %s: Fortran failing checks %s, expected %s' % (
                    name, fchecks.group(1) if fchecks else None, ordered(want)))
    assert not errs, 'planted defects:\n' + '\n'.join(errs)


# ---------------------------------------------------------------------------
def test_two_rank_fortran(tmp):
    fortran_binary()
    errs = []
    for case in ('test.tpv8', 'test.tpv36'):
        d = os.path.join(tmp, 'two_' + case)
        make_case(case, d, 'par.nx = 2\n')
        rc, out = run_fortran_mesh_only(d, ranks=2)
        v = parse_verdict(out, rc, 'fortran %s 2 ranks' % case)
        print('  %-12s 2 ranks: %s rc=%d checks=%s elements=%s fault_pairs=%s attachments=%s' % (
            case, v[0], rc, ','.join(v[1]) or 'none', v[2]['elements'], v[2]['fault_pairs'], v[2]['sides_attachments']))
        e = judge(case, rc, v[0], v[1], v[2])
        if e:
            errs.append('2 ranks: ' + e)
    assert not errs, '2-rank:\n' + '\n'.join(errs)


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
    try:
        run_all(tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_all(tmp):
    tests = [('thresholds_equal', lambda: test_thresholds_equal()),
             ('judge_both_ways', lambda: test_judge_both_ways()),
             ('planted_defects', lambda: test_planted_defects(tmp)),
             ('two_rank_fortran', lambda: test_two_rank_fortran(tmp)),
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
