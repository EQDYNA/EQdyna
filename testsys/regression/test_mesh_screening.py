#! /usr/bin/env python3
"""
Regression guard (board row 167): the post-meshing screen
(src/python/eqdyna/meshscreen.py + src/fortran/screenMesh.f90; exit 81 sides,
82 bridging) is GREEN on every registered case, and each of its two checks can
actually fire.

Tests (every one a hard failure, none skips):

  1. ALL CASES, both backends: `python3 -m eqdyna <case> --mesh-only` and
     `mpirun -np 1 eqdyna --mesh-only` for every name in testNameList.nameList.
     Every case exits 0 and prints `mesh-screen: PASS`. Both backends must
     agree on verdict, exit code, failing checks and every count; a verdict
     over zero elements, nodes, fault pairs or attachments is refused (a screen
     that compared nothing proves nothing).
  2. JUDGE BOTH WAYS (rule 14a): judge() is fed good and bad results and must
     accept the good and flag every bad one.
  3. PLANTED DEFECTS on the in-memory tpv8 mesh, Python and Fortran (the dump
     goes to screen_mesh_dump_driver.f90, linked against the production
     screenMesh.o): one defect per check trips exactly that check; the
     unplanted mesh is the control.
  4. 2-RANK FORTRAN: tpv8 and tpv36 (wedge mesh) at mpirun -np 2.
  6. PASS LINE ON A NORMAL RUN (row 169 (e)): a solving run (python numpy, fortran
     1 and 2 ranks, a 1-2 step term) prints exactly one `mesh-screen: PASS`
     line (rank 0), the two backends with the same counts.
  7. jax-MPI ABORT (row 169 (f)): a real 2-rank mpi4py run in which the screen
     fails on rank 1 only must MPI_Abort with exit status 81, not hang.
  8. STRIKE -90 TIE CELLS (row 169 (h)): -90 moves only round-off tie cells;
     the error text must not say it flips the side test.
  5. INDEPENDENT SIDE ORACLE (python backend, the Fortran connectivity being
     held identical by test 1's counts): screen_mesh's `sides` check and the
     meshgen side test share one inequality and one `un`, so they cannot catch
     an error in either. oracle_sides() takes each fault's plane normal from the pair
     nodes' OWN coordinates (cross product of the in-plane spans, oriented to
     +y: the declared dip direction) and tests every (cell,node) attachment against THAT plane.
     Run on every registered case's gate mesh, and on the small C_degen
     fixture, where the same oracle must flag the wrong attachments of the
     SAME mesh built with fstrike=0 (un then points along +x: a wrong rule
     input the screen passes).

mpirun on this class of box needs `env -u DISPLAY HWLOC_COMPONENTS=-gl,...` (a
dead X server blocks hwloc GL discovery forever) and is wrapped in
`timeout -s KILL`.
"""
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

OBJS = ('globalvar.o', 'errorCodes.o', 'screenMesh.o')
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

import testNameList  # noqa: E402
from testsys.common import MACHINE  # noqa: E402
from eqdyna import eqdyna3d, meshscreen  # noqa: E402

FORTRAN_DIR = os.path.join(ROOT, 'src', 'fortran')
FORTRAN_BIN = os.path.join(ROOT, 'bin', 'eqdyna')
DRIVER_F90 = os.path.join(ROOT, 'testsys', 'parity', 'screen_mesh_dump_driver.f90')
MPIRUN_TIMEOUT_S = 1500
COUNT_KEYS = ('elements', 'nodes', 'fault_pairs', 'sides_attachments')
PARITY_KEYS = COUNT_KEYS + ('sides_wrong_master', 'sides_wrong_slave', 'bridging_cells')

# plant -> checks it must trip on the planted tpv8 mesh. bridging_cell is the
# realistic wiring defect: one of the two nodes is always on the wrong side.
PLANT_EXPECTED = {
    'sides_swap_node': ('sides',),
    'bridging_shear': ('bridging',),
    'bridging_cell': ('sides', 'bridging'),
    # a MASTER node attached to a cell beyond the fault patch on the footwall side (row 169 (b)):
    # the patch-edge exemption used to hide this; only slave attachments are exempt there now
    'sides_outside_master': ('sides',),
}


def env():
    e = dict(os.environ)
    e['EQDYNAROOT'] = ROOT
    e['PATH'] = os.pathsep.join([os.path.join(ROOT, 'scripts'), e.get('PATH', '')])
    e['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    e.pop('DISPLAY', None)
    e['HWLOC_COMPONENTS'] = '-gl,-opencl,-cuda,-nvml'
    e['OMP_NUM_THREADS'] = e['OPENBLAS_NUM_THREADS'] = e['MKL_NUM_THREADS'] = '1'
    return e


def ordered(checks):
    return tuple(k for k in meshscreen.SCREEN_CODES if k in checks)


def fortran_binary():
    assert os.path.exists(FORTRAN_BIN), (
        'Fortran binary %s missing: run ./install-eqdyna.sh -m ubuntu first' % FORTRAN_BIN)


def make_case(case, case_dir, extra=''):
    """create.newcase + serial decomposition + case.setup. Raises on failure."""
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'), case_dir, case],
                       env=env(), capture_output=True, text=True)
    assert r.returncode == 0, 'create.newcase %s failed: %s' % (case, (r.stderr or r.stdout)[-800:])
    p = os.path.join(case_dir, 'user_defined_params.py')
    with open(p) as f:
        text = f.read().rstrip('\n')
    with open(p, 'w') as f:
        f.write(text + '\n\npar.nx = 1\npar.ny = 1\npar.nz = 1\n' + extra)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env(), capture_output=True, text=True)
    assert r.returncode == 0, 'case.setup %s failed: %s' % (case, (r.stderr or r.stdout)[-800:])


def parse_verdict(out, rc, who):
    """(verdict, checks, counts) from the `mesh-screen:` line(s), one per rank:
    FAIL if any rank failed, checks the union, counts summed over ranks."""
    lines = [l for l in out.splitlines() if l.startswith('mesh-screen:')]
    assert lines, '%s: no mesh-screen line (rc %d); tail:\n%s' % (who, rc, out[-1200:])
    verdict, checks, counts = 'PASS', set(), dict.fromkeys(PARITY_KEYS, 0)
    for l in lines:
        m = re.match(r'mesh-screen: (PASS|FAIL) checks=(\S+) (.*)$', l)
        assert m, '%s: unparsable mesh-screen line: %r' % (who, l)
        if m.group(1) == 'FAIL':
            verdict = 'FAIL'
        if m.group(2) != 'none':
            checks |= set(m.group(2).split(','))
        kv = dict(x.split('=') for x in m.group(3).split())
        for k in PARITY_KEYS:
            counts[k] += int(kv[k])
    return verdict, ordered(checks), counts


def run_python(case_dir):
    r = subprocess.run([sys.executable, '-m', 'eqdyna', case_dir, '--mesh-only'], env=env(),
                       capture_output=True, text=True, timeout=3000)
    return r.returncode, r.stdout + r.stderr


def run_fortran(case_dir, ranks=1, args=('--mesh-only',)):
    r = subprocess.run(['timeout', '-s', 'KILL', str(MPIRUN_TIMEOUT_S), 'mpirun', '-np', str(ranks),
                        FORTRAN_BIN, *args], cwd=case_dir, env=env(), capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def judge(case, rc, verdict, checks, counts):
    """Error string for one backend's result on one case, or None. The ONLY
    place a run becomes a verdict; test_judge_both_ways feeds it known-bad input."""
    empty = [k for k in COUNT_KEYS if counts[k] <= 0]
    if empty:
        return '%s: screen compared nothing (%s == 0): a vacuous verdict' % (case, ','.join(empty))
    if not (verdict == 'PASS' and rc == 0):
        return '%s: must be GREEN, got %s rc=%d checks=%s' % (case, verdict, rc, ','.join(checks))
    return None


def test_judge_both_ways():
    counts = dict.fromkeys(PARITY_KEYS, 5)
    case = testNameList.nameList[0]
    assert judge(case, 0, 'PASS', (), counts) is None
    bad = {
        'red case (sides)': judge(case, 81, 'FAIL', ('sides',), counts),
        'red case (bridging)': judge(case, 82, 'FAIL', ('bridging',), counts),
        'vacuous (no attachments)': judge(case, 0, 'PASS', (), dict(counts, sides_attachments=0)),
        'vacuous (no fault pairs)': judge(case, 0, 'PASS', (), dict(counts, fault_pairs=0)),
    }
    missed = [k for k, v in bad.items() if v is None]
    assert not missed, 'judge failed to flag: %s' % missed


def test_all_cases(tmp):
    names = list(testNameList.nameList)
    assert names, 'nameList empty'
    fortran_binary()
    dirs = {c: os.path.join(tmp, 'all_' + c) for c in names}

    # Independent runs through one pool of half the cores (cap 8: a 3.5M-element case needs GBs); every
    # worker is single-threaded (env()), so a 4-core runner runs 2 at a time. Row 169 (g): the unit of work
    # is a (case, backend) run, not a case: the one-task-per-case pool ran a case's python and fortran
    # back to back, so the chain of the largest case set the floor of the shard. Python runs are queued
    # first (the numpy mesh build is the longer of the two). No case, backend or count comparison is dropped.
    workers = min(max(1, len(os.sched_getaffinity(0))), 8)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(lambda c: make_case(c, dirs[c]), names))
        futs = {(c, 'python'): ex.submit(run_python, dirs[c]) for c in names}
        futs.update({(c, 'fortran'): ex.submit(run_fortran, dirs[c]) for c in names})
        res = {c: (futs[(c, 'python')].result(), futs[(c, 'fortran')].result()) for c in names}
    errs = []
    for c in names:
        (prc, pout), (frc, fout) = res[c]
        pv = parse_verdict(pout, prc, 'python ' + c)
        fv = parse_verdict(fout, frc, 'fortran ' + c)
        print('  %-18s python %s rc=%-2d fortran %s rc=%-2d checks=%-8s %s' % (
            c, pv[0], prc, fv[0], frc, ','.join(pv[1]) or 'none',
            ' '.join('%s=%d' % (k, pv[2][k]) for k in PARITY_KEYS)))
        for who, rc, v in (('python', prc, pv), ('fortran', frc, fv)):
            e = judge(c, rc, *v)
            if e:
                errs.append('[%s] %s' % (who, e))
        if (pv[0], pv[1], prc, pv[2]) != (fv[0], fv[1], frc, fv[2]):
            errs.append('%s: PARITY BROKEN python %s rc=%d %s %s vs fortran %s rc=%d %s %s' % (
                c, pv[0], prc, pv[1], pv[2], fv[0], frc, fv[1], fv[2]))
    assert not errs, 'all-cases failed:\n' + '\n'.join(errs)


class Mesh:
    """The arrays build_solver_state hands to screen_mesh, captured in memory:
    (meshCoor, conn, elem_type, nsmp, un, tol)."""

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
        assert self.a[1].shape[0] > 0 and self.a[3].shape[0] > 0, 'empty mesh handed to the screen'

    def copy(self):
        return [np.array(v, copy=True) for v in self.a[:5]] + [self.a[5]]


def plant(a, name):
    """Apply ONE defect, in place, at the pair nearest the centre of the fault patch."""
    X, conn, et, nsmp, un, tol = a
    xm, zm = X[nsmp[:, 1], 0], X[nsmp[:, 1], 2]
    p = int(np.argmin(np.abs(xm - 0.5 * (xm.min() + xm.max())) + np.abs(zm - 0.5 * (zm.min() + zm.max()))))
    s, m = int(nsmp[p, 0]), int(nsmp[p, 1])
    e = int(np.nonzero((conn == m).any(axis=1))[0][0])
    if name == 'sides_swap_node':
        # a hanging-wall cell holds the SLAVE node in place of the master (coincident nodes)
        conn[e][conn[e] == m] = s
    elif name == 'bridging_shear':
        X[s, 0] += 1.0       # the slave node 1 m from its master
    elif name == 'bridging_cell':
        # one cell holds BOTH nodes: a corner off the fault plane becomes the slave
        k = int(np.nonzero(np.abs(X[conn[e], 1] - X[m, 1]) > 0)[0][0])
        conn[e][k] = s
    elif name == 'sides_outside_master':
        # a cell beyond the patch's left edge, below the plane, holds the MASTER node instead of the slave
        e, s, m = outside_cell(X, conn, nsmp, un)
        conn[e][conn[e] == s] = m
    else:
        raise AssertionError('unknown plant ' + name)


def outside_cell(X, conn, nsmp, un):
    """(cell, slave node, master node): a cell beyond the left edge of the patch in x,
    below the fault plane (its centre on the -un side), holding the slave node of a
    left-edge pair at mid depth."""
    xm, zm = X[nsmp[:, 1], 0], X[nsmp[:, 1], 2]
    edge = np.nonzero(xm == xm.min())[0]
    p = int(edge[np.argmin(np.abs(zm[edge] - 0.5 * (zm.min() + zm.max())))])
    s, m = int(nsmp[p, 0]), int(nsmp[p, 1])
    cells = np.nonzero((conn == s).any(axis=1))[0]
    cen = X[conn[cells]].sum(axis=1) / 8.0
    below = ((cen - X[s]) @ un[p + 1] < 0) & (cen[:, 0] < xm.min())
    assert below.any(), 'no cell beyond the patch edge below the plane holds slave node %d' % s
    return int(cells[np.nonzero(below)[0][0]]), s, m


def dump_mesh(a, path):
    """Raw stream dump in the layout screen_mesh_dump_driver.f90 reads."""
    X, conn, et, nsmp, un, tol = a
    nft = nsmp.shape[0]
    with open(path, 'wb') as f:
        f.write(np.array([X.shape[0] - 1, conn.shape[0], nft], dtype='<i8').tobytes())
        f.write(np.array([tol], dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(X[1:], dtype='<f8').tobytes())
        f.write(np.ascontiguousarray(conn, dtype='<i8').tobytes())
        f.write(np.ascontiguousarray(et, dtype='<i4').tobytes())
        f.write(np.ascontiguousarray(nsmp[:, :2], dtype='<i8').tobytes())
        f.write(np.ascontiguousarray(un[1:nft + 1], dtype='<f8').tobytes())


def build_driver(tmp):
    fortran_binary()
    mk = dict(env(), MACHINE=MACHINE)
    r = subprocess.run(['make', *OBJS], cwd=FORTRAN_DIR, env=mk, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, 'make %s failed:\n%s' % (' '.join(OBJS), (r.stdout + r.stderr)[-1500:])
    exe = os.path.join(tmp, 'screen_mesh_dump_driver')
    r = subprocess.run(['mpif90', '-fopenmp', '-ffree-line-length-none', '-O1', '-I' + FORTRAN_DIR, '-J' + tmp,
                        DRIVER_F90, *[os.path.join(FORTRAN_DIR, o) for o in OBJS],
                        '-o', exe], capture_output=True, text=True, timeout=600, env=env())
    assert r.returncode == 0, 'driver build failed:\n' + (r.stdout + r.stderr)[-1500:]
    return exe


def run_driver(exe, a, path):
    dump_mesh(a, path)
    r = subprocess.run([exe, path], capture_output=True, text=True, env=env(), timeout=600)
    got = re.search(r'mesh-screen: FAIL checks=(\S+)', r.stdout)
    return r.returncode, ordered(got.group(1).split(',')) if got else ()


def test_planted_defects(tmp):
    d = os.path.join(tmp, 'plant_tpv8')
    make_case('test.tpv8', d)
    mesh = Mesh(d)
    exe = build_driver(tmp)
    errs = []
    # control: the unplanted mesh passes in BOTH implementations
    m, fail = meshscreen._screen_measure(*mesh.copy())
    rc, fchecks = run_driver(exe, mesh.copy(), os.path.join(tmp, 'control.dump'))
    if fail or rc != 0:
        errs.append('control (unplanted tpv8) must pass: python %s, fortran rc %d' % (sorted(fail), rc))
    for name, want in PLANT_EXPECTED.items():
        a = mesh.copy()
        plant(a, name)
        m, fail = meshscreen._screen_measure(*a)
        rc, fchecks = run_driver(exe, a, os.path.join(tmp, name + '.dump'))
        want_code = meshscreen.SCREEN_CODES[want[0]]
        print('  plant %-16s python %-16s fortran rc=%d %s (expected %s, rc %d)' % (
            name, ','.join(ordered(fail)), rc, ','.join(fchecks), ','.join(want), want_code))
        if ordered(fail) != want:
            errs.append('plant %s: python tripped %s, expected %s: %s' % (name, ordered(fail), want, fail))
        if (rc, fchecks) != (want_code, want):
            errs.append('plant %s: fortran rc %d checks %s, expected %d %s' % (name, rc, fchecks, want_code, want))
    assert not errs, 'planted defects:\n' + '\n'.join(errs)


def oracle_sides(X, conn, nsmp, tol):
    """(attachments, wrong_master, wrong_slave, min_margin) of every distinct
    (cell,node) attachment to a pair node, judged by the plane spanned by the
    fault's own master nodes (never by un): every MASTER attachment, wherever
    the cell is (the one-cell margin past the patch, where the replace loop
    still swaps, included), and every SLAVE attachment whose cell centre is
    inside the patch window (slave attachments beyond it are the edge-weld
    convention, exempt like the screen's). X is 1-indexed (row 0 unused); conn (E,8) 1-indexed."""
    nft = nsmp.shape[0]
    role = np.zeros(X.shape[0], dtype=np.int8)
    pid = np.full(X.shape[0], -1, dtype=np.int64)
    role[nsmp[:, 0]], role[nsmp[:, 1]] = 1, 2
    pid[nsmp[:, 0]] = pid[nsmp[:, 1]] = np.arange(nft)
    normal, lo, hi = np.zeros((nft, 3)), np.zeros((nft, 3)), np.zeros((nft, 3))
    for f in np.unique(nsmp[:, 2]):
        sel = np.nonzero(nsmp[:, 2] == f)[0]
        P = X[nsmp[sel, 1]]
        assert np.abs(X[nsmp[sel, 0]] - P).max() <= tol, 'fault %d pair does not coincide' % f
        # local normal from the fault nodes' own in-plane spans (central differences
        # along the fault's x and z grid lines): exact for a plane, local for a rough fault
        ix, iz = (np.unique(P[:, k], return_inverse=True)[1].ravel() for k in (0, 2))
        G = np.full((ix.max() + 1, iz.max() + 1), -1)
        G[ix, iz] = np.arange(sel.size)

        def nb(dx_, dz_):
            j = G[np.clip(ix + dx_, 0, G.shape[0] - 1), np.clip(iz + dz_, 0, G.shape[1] - 1)]
            return np.where(j >= 0, j, np.arange(sel.size))
        n = np.cross(P[nb(1, 0)] - P[nb(-1, 0)], P[nb(0, 1)] - P[nb(0, -1)])
        n /= np.linalg.norm(n, axis=1)[:, None]
        n = np.where(n[:, 1:2] > 0, n, -n)         # declared dip direction: master side is +y
        normal[sel], lo[sel], hi[sel] = n, P.min(axis=0), P.max(axis=0)
    ce, sl = np.nonzero(role[conn] != 0)
    ce, nodes = np.unique(np.stack([ce, conn[ce, sl]], axis=1), axis=0).T
    pp, rl = pid[nodes], role[nodes]
    cen = X[conn[ce]].sum(axis=1) / 8.0
    inside = ((cen[:, [0, 2]] >= lo[pp][:, [0, 2]]) & (cen[:, [0, 2]] <= hi[pp][:, [0, 2]])).all(axis=1)
    ok = inside | (rl == 2)
    ce, nodes, pp, rl, cen, inside = ce[ok], nodes[ok], pp[ok], rl[ok], cen[ok], inside[ok]
    margin = np.where(rl == 2, 1.0, -1.0) * np.einsum('ij,ij->i', cen - X[nodes], normal[pp])
    # A master attachment beyond the patch whose cell centre is within the mesh tolerance of the plane is a
    # round-off tie (wedge cells at the patch edge: |margin| 1e-14 m against >= 100 m for every decided
    # attachment, measured on the gate meshes), decided by rounding, not a side: not judged. Inside the patch
    # the test stays strict.
    wrong = ~(margin > 0.0) & ~(~inside & (margin > -tol))
    return (int(margin.size), int((wrong & (rl == 2)).sum()), int((wrong & (rl == 1)).sum()),
            float(margin.min()) if margin.size else float('nan'))


def fixture_mesh(fstrike):
    """The small C_degen fixture (testsys/unit/test_meshgen_c_degen_wedge.py) meshed
    with the given fstrike: (X, conn, nsmp, un, tol)."""
    from eqdyna import meshgen
    sys.path.insert(0, os.path.join(ROOT, 'testsys', 'unit'))
    from test_meshgen_c_degen_wedge import _small_c_degen36_params
    p = dict(_small_c_degen36_params(), fstrike=fstrike)
    xl, yl, zl, pmlb, _ = meshgen.build_grid_lines(p)
    X, _, nsmp = meshgen.build_node_coordinates(xl, yl, zl, p)
    un = meshgen.build_fault_geometry(xl, yl, zl, p, nsmp)[0]
    conn = meshgen.build_elements(xl, yl, zl, p, pmlb, nsmp, np.array([[6000.0, 3464.0, 2670.0]]), X, un)[0]
    return X, conn, nsmp, un, p['tol']


def test_side_oracle(tmp):
    errs = []
    # fixture, both ways: fstrike 270 (valid) is clean, fstrike 0 must be flagged
    X, conn, nsmp, un, tol = fixture_mesh(270.0)
    good = oracle_sides(X, conn, nsmp, tol)
    X, conn, nsmp, un, tol = fixture_mesh(0.0)
    bad = oracle_sides(X, conn, nsmp, tol)
    scr = meshscreen._screen_measure(X, conn, np.zeros(conn.shape[0], dtype=int), nsmp, un, tol)[1]
    print('  fixture oracle: fstrike=270 %s ; fstrike=0 %s ; screen on the fstrike=0 mesh: %s' % (
        good, bad, 'FAIL ' + ','.join(scr) if scr else 'PASS (blind)'))
    if good[0] == 0 or good[1] + good[2] != 0:
        errs.append('fixture fstrike=270: oracle must see attachments and zero wrong, got %s' % (good,))
    if bad[1] + bad[2] == 0:
        errs.append('oracle did not flag the fixture meshed with fstrike=0 (it must: un points along x)')
    for c in testNameList.nameList:
        d = os.path.join(tmp, 'oracle_' + c)
        make_case(c, d)
        a = Mesh(d).a
        r = oracle_sides(a[0], a[1], a[3], a[5])
        print('  %-18s oracle attachments=%d wrong_master=%d wrong_slave=%d min_margin=%.3E' % ((c,) + r))
        if r[0] == 0 or r[1] + r[2] != 0:
            errs.append('%s: oracle attachments=%d wrong_master=%d wrong_slave=%d' % ((c,) + r[:3]))
        if c == 'test.tpv8':
            # the margin must be checked: a wrong-side master attachment beyond the patch edge is flagged
            pl = [np.array(v, copy=True) for v in a[:5]] + [a[5]]
            plant(pl, 'sides_outside_master')
            pr = oracle_sides(pl[0], pl[1], pl[3], pl[5])
            print('  %-18s oracle on the planted outside master: wrong_master=%d (clean mesh %d)' % (c, pr[1], r[1]))
            if pr[1] != 1:
                errs.append('oracle missed (or over-counted) the planted outside-patch master attachment: wrong_master=%d' % pr[1])
        shutil.rmtree(d, ignore_errors=True)
    assert not errs, 'side oracle:\n' + '\n'.join(errs)


def test_strike_refused(tmp):
    """par.strike drives un, hence the side test: only exactly 270 is accepted, both
    backends exit 35 (ERR_GEOM_STRIKE_UNSUPPORTED) naming par.strike on every other
    value, including -90 (equal mod 360, but it moves tie cells at the patch edge)."""
    fortran_binary()
    errs = []
    for strike, want in ((270.0, 0), (0.0, 35), (90.0, 35), (-90.0, 35)):
        d = os.path.join(tmp, 'strike_%g' % strike)
        make_case('test.tpv8', d, 'par.strike = %r\n' % strike)
        for who, (rc, out) in (('python', run_python(d)), ('fortran', run_fortran(d))):
            named = 'par.strike' in out
            print('  strike %6g %-8s rc=%d (expected %d)' % (strike, who, rc, want))
            if rc != want or (want and not named):
                errs.append('strike %g %s: rc %d, expected %d, par.strike named: %s' % (strike, who, rc, want, named))
    assert not errs, 'strike refusal:\n' + '\n'.join(errs)


def test_two_rank_fortran(tmp):
    fortran_binary()
    errs = []
    for case in ('test.tpv8', 'test.tpv36'):
        d = os.path.join(tmp, 'two_' + case)
        make_case(case, d, 'par.nx = 2\n')
        rc, out = run_fortran(d, ranks=2)
        v = parse_verdict(out, rc, 'fortran %s 2 ranks' % case)
        print('  %-12s 2 ranks: %s rc=%d checks=%s %s' % (
            case, v[0], rc, ','.join(v[1]) or 'none', ' '.join('%s=%d' % (k, v[2][k]) for k in PARITY_KEYS)))
        e = judge(case, rc, *v)
        if e:
            errs.append('2 ranks: ' + e)
    assert not errs, '2-rank:\n' + '\n'.join(errs)


def test_pass_line_normal_run(tmp):
    """Row 169 (e): a SOLVING run logs the screen too: exactly one PASS line, rank 0 only."""
    fortran_binary()
    d = os.path.join(tmp, 'pass_line')
    make_case('test.tpv8', d, 'par.term = 0.05\n')       # dt = 0.0417 s: 2 steps
    r = subprocess.run([sys.executable, '-m', 'eqdyna', d, '2', '--backend', 'numpy'], env=env(),
                       capture_output=True, text=True, timeout=3000)
    assert r.returncode == 0, 'python numpy run failed rc=%d:\n%s' % (r.returncode, (r.stdout + r.stderr)[-1500:])
    got = {'python': r.stdout}
    d2 = os.path.join(tmp, 'pass_line2')
    make_case('test.tpv8', d2, 'par.term = 0.05\npar.nx = 2\n')
    for who, dd, ranks in (('fortran 1 rank', d, 1), ('fortran 2 ranks', d2, 2)):
        rc, out = run_fortran(dd, ranks, args=())
        assert rc == 0, '%s normal run failed rc=%d:\n%s' % (who, rc, out[-1500:])
        got[who] = out
    counts = {}
    for who, out in got.items():
        n = sum(l.startswith('mesh-screen: PASS') for l in out.splitlines())
        print('  %-16s mesh-screen PASS lines: %d' % (who, n))
        assert n == 1, '%s: expected exactly one mesh-screen: PASS line on a normal run, got %d' % (who, n)
        counts[who] = parse_verdict(out, 0, who)
    assert counts['python'][2] == counts['fortran 1 rank'][2], 'python vs fortran counts differ: %s' % counts


MPI_ABORT_DRIVER = """
import sys
sys.path.insert(0, %(src)r)
from mpi4py import MPI
from eqdyna import eqdyna3d, meshscreen
real = meshscreen.screen_mesh
def failing_on_rank_1(*a):
    if MPI.COMM_WORLD.Get_rank() == 1:      # a rank-local verdict: the mesh each rank screens is its own box
        raise meshscreen.MeshScreenError({'sides': 'planted failure on rank 1'}, {})
    return real(*a)
meshscreen.screen_mesh = failing_on_rank_1
sys.argv = ['eqdyna', %(case)r, '--mpi', '--backend', 'jax']
eqdyna3d.main()
"""


def test_mpi_abort(tmp):
    """Row 169 (f): eqdyna3d.main's jax-MPI `except MeshScreenError` -> _abort(comm=comm) ->
    comm.Abort(81), driven by a real 2-rank mpi4py job. The stand-in is the failure itself
    (screen_mesh raises on rank 1 only; a real rank-local defect cannot be planted without
    editing meshgen); everything after it is the production code. Without the Abort rank 0
    waits in the setup exchange (a hang, killed by the timeout: rc 137) or the launcher
    reports a plain process exit instead of an MPI abort."""
    d = os.path.join(tmp, 'mpi_abort')
    make_case('test.tpv8', d, 'par.nx = 2\n')
    drv = os.path.join(tmp, 'mpi_abort_driver.py')
    with open(drv, 'w') as f:
        f.write(MPI_ABORT_DRIVER % dict(src=os.path.join(ROOT, 'src', 'python'), case=d))
    r = subprocess.run(['timeout', '-s', 'KILL', '600', 'mpirun', '-np', '2', sys.executable, drv],
                       cwd=d, env=env(), capture_output=True, text=True)
    out = r.stdout + r.stderr
    print('  mpirun rc=%d, MPI abort reported: %s, FATAL block: %s' % (
        r.returncode, bool(re.search(r'MPI_ABORT|MPI_Abort', out)), 'EQdyna: FATAL' in out))
    assert r.returncode == meshscreen.ERR_MESH_SCREEN_SIDES, (
        'expected exit %d (MPI_Abort with the check code), got %d (137 = hung and killed):\n%s' % (
            meshscreen.ERR_MESH_SCREEN_SIDES, r.returncode, out[-1500:]))
    assert 'EQdyna: FATAL' in out and 'planted failure on rank 1' in out, 'FATAL block missing:\n' + out[-1500:]
    assert re.search(r'MPI_ABORT|MPI_Abort', out), 'the job did not end in an MPI abort:\n' + out[-1500:]


def test_strike_minus90_ties():
    """Row 169 (h): -90 is equal to 270 mod 360 and moves only round-off tie cells of the
    side test (measured on the gate meshes: 0 differing elements on tpv8 and tpv23, 3 on
    tpv36), it does NOT flip the cell-side test. Every element whose connectivity differs
    between the two spellings must have a centre within round-off of the plane; the strike
    error text (errorCodes.f90, docs/user/troubleshooting.md) must not claim a flip."""
    X, c270, nsmp, un, tol = fixture_mesh(270.0)
    _, cm90, _, un90, _ = fixture_mesh(-90.0)
    diff = np.nonzero((c270 != cm90).any(axis=1))[0]
    pair = np.zeros(X.shape[0], dtype=np.int64)
    pair[nsmp[:, 0]] = np.arange(1, nsmp.shape[0] + 1)
    worst = 0.0
    for e in diff:
        for k in np.nonzero(c270[e] != cm90[e])[0]:
            n = c270[e][k]
            cen = X[c270[e]].sum(axis=0) / 8.0
            worst = max(worst, abs(float((cen - X[n]) @ un[pair[n]])) / max(float(np.linalg.norm(cen - X[n])), 1e-300))
    print('  fixture: un max |270 - (-90)| = %.3E ; %d differing elements, worst tie |cos| %.3E' % (
        np.abs(un - un90).max(), diff.size, worst))
    assert np.abs(un - un90).max() < 1e-12, 'un differs by more than round-off between 270 and -90'
    assert worst < 1e-9, 'a cell that differs between -90 and 270 is not a round-off tie (|cos| %.3E)' % worst
    for path in ('docs/user/troubleshooting.md', 'src/fortran/errorCodes.f90'):
        with open(os.path.join(ROOT, path)) as f:
            text = f.read()
        assert 'flips the cell-side' not in text, path + ' says -90 flips the cell-side test; measured: only round-off ties move'


def main():
    tmp = tempfile.mkdtemp(prefix='mesh_screening_')
    tests = [('judge_both_ways', lambda: test_judge_both_ways()),
             ('planted_defects', lambda: test_planted_defects(tmp)),
             ('strike_refused', lambda: test_strike_refused(tmp)),
             ('two_rank_fortran', lambda: test_two_rank_fortran(tmp)),
             ('side_oracle', lambda: test_side_oracle(tmp)),
             ('strike_minus90_ties', lambda: test_strike_minus90_ties()),
             ('pass_line_normal_run', lambda: test_pass_line_normal_run(tmp)),
             ('mpi_abort', lambda: test_mpi_abort(tmp)),
             ('all_cases', lambda: test_all_cases(tmp))]
    sel = sys.argv[1:]
    bad = 0
    try:
        for name, fn in tests:
            if sel and name not in sel:
                continue
            try:
                fn()
                print('SUCCESS %s' % name)
            except AssertionError as exc:
                bad += 1
                print('FAIL %s: %s' % (name, exc))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
