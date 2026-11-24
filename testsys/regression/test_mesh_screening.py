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


def run_fortran(case_dir, ranks=1):
    r = subprocess.run(['timeout', '-s', 'KILL', str(MPIRUN_TIMEOUT_S), 'mpirun', '-np', str(ranks),
                        FORTRAN_BIN, '--mesh-only'], cwd=case_dir, env=env(), capture_output=True, text=True)
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

    def one(c):
        make_case(c, dirs[c])
        return run_python(dirs[c]), run_fortran(dirs[c])
    # independent cases through a pool of half the cores (cap 8: a 3.5M-element
    # case needs GBs); every worker is single-threaded (env()), so a 4-core
    # runner runs 2 at a time. No case, backend or count comparison is dropped.
    workers = min(max(1, len(os.sched_getaffinity(0)) // 2), 8)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {c: ex.submit(one, c) for c in names}
        res = {c: f.result() for c, f in futs.items()}
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
    else:
        raise AssertionError('unknown plant ' + name)


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


def main():
    tmp = tempfile.mkdtemp(prefix='mesh_screening_')
    tests = [('judge_both_ways', lambda: test_judge_both_ways()),
             ('planted_defects', lambda: test_planted_defects(tmp)),
             ('two_rank_fortran', lambda: test_two_rank_fortran(tmp)),
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
