#! /usr/bin/env python3
"""
Wedge-element MASS + INTERNAL-FORCE parity probe: Fortran vs the Python port,
on synthetic element geometries -- the missing piece the C_degen dynamics
port needed before any full-case run could be trusted (see
testsys/parity/evidence_c_degen_port.py, whose item (b) only proved "5 numpy
steps ran without raising", not parity).

SAME SHAPE AS testsys/regression/test_drucker_prager_kernel.py, generalized
from one kernel (calcElemKU's D-P block) to the FULL wedge-relevant chain:

    calcLocalShapeFunc -> calcGlobalShapeFunc (the elemTypeArr==11/12 merge)
    -> contm (per-node lumped mass) -> calcSSPhi4Hrgls (ss/phi, calls vlm)
    -> calcElemMass (gravity/body-force, proven no-op here) -> calcElemKU
    (interior force) -> calcHourglassResist (KF78 hourglass force)

compiled from the REAL, UNMODIFIED src/fortran/*.f90 (testsys/parity/
probe_wedge_kernel.f90 is a NEW driver, calling those subroutines directly;
no existing src/fortran file is touched), and compared against the
corresponding calls into the Python port:

    eqdyna.assembleGlobalMass.compute_element_shape / contm /
    _contm_wedge_node_mass / compute_hourglass
    eqdyna.assembleGlobalKU.build / alloc_scratch / assembleGlobalKU /
    calcHourglassResist

on a synthetic ONE-ELEMENT solver state (S/inv), not a re-implementation of
the arithmetic -- these are the exact functions eqdyna3d.py's
build_solver_state / driver.py's time loop call in production.

THREE CASES: hex (elemType=1, REGRESSION GUARD -- untouched by the wedge
fix, must stay bit-identical to the already-shipped non-wedge port),
wedge_degenerate (elemType=11, generic skew, local nodes 3==4 and 7==8
coincide), wedge_shallow_dip15 (elemType=12, same collapse but built by
rotating a cube 15 degrees about the strike axis -- tpv36's par.dip -- so
the off-diagonal Jacobian terms are the ones a real shallow-dipping wedge
produces). elemType 11 vs 12 take IDENTICAL code in every kernel here (per
direct reading), so covering both is a complete check of the
"==11 .or. ==12" / ">10" predicates, not a claim they differ numerically.

IDENTITY EQUATION-NUMBER MAP: both sides map (node i, local dof j) to a
UNIQUE global equation slot with no sharing, so the scattered force lands in
exactly one slot per (node, dof) and the comparison is element-local.

C_elastic is fixed to 1 (elastic) for every case: the Drucker-Prager branch
is already isolated and verified by test_drucker_prager_kernel.py, and this
port's own refusal (C_degen>3 + C_elastic==0 -> NotImplementedError) means a
wedge element can never reach that branch in production.

TOLERANCE: rtol=1e-10, atol=1e-9. Same-order IEEE-754 double arithmetic on
both sides, operand-for-operand matched by direct reading of both the
Fortran and the Python. Not bit-identical because compute_element_shape uses
np.einsum for the Jacobian/cofactor contractions while the Fortran uses an
explicit nested-loop sum -- both correctly-rounded IEEE-754 but not
guaranteed to associate identically, the same allowance already implicit in
the already-shipped non-wedge path.

Cheap: compiles a handful of small .f90 files (plus a MACHINE=ubuntu clean
rebuild of the full src/fortran tree, needed because contm/calcSSPhi4Hrgls/
calcHourglassResist live in files that include mpif.h) and runs 3 element
cases -- seconds. Fails loudly if gfortran/mpif90 is unavailable or the
build fails -- never silently skips.
"""
import os
import shutil
import subprocess
import sys
import time

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
PROBE_SRC = os.path.join(TESTSYS, 'probe_wedge_kernel.f90')
FSRC = os.path.join(REPO_ROOT, 'src', 'fortran')

sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))
from eqdyna import assembleGlobalKU, assembleGlobalMass  # noqa: E402

RTOL, ATOL = 1e-10, 1e-9

# ---------------------------------------------------------------------------
# synthetic geometry / material / kinematics, IDENTICAL bytes to both sides
# ---------------------------------------------------------------------------
_ACOOR = np.array([
    [-1.0, -1.0, -1.0], [1.0, -1.0, -1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 1.0], [-1.0, 1.0, 1.0],
])  # (8,3), Fortran node order, LOCAL [-1,1]^3 coords (calcLocalShapeFunc.f90)

DX = 500.0
DIP_DEG = 15.0  # tpv36's par.dip
RNG = np.random.default_rng(20260916)


def _hex_from_T(T, origin, degenerate):
    """8 physical corners = origin + T @ local, local = _ACOOR rows.
    degenerate: collapse local node 4 -> node 3's coords, 8 -> 7's (the
    calcGlobalShapeFunc.f90:22-28 merge's own node pairing) -- physically a
    hex degenerated to a wedge along that edge."""
    phys = (T @ _ACOOR.T).T + origin[None, :]   # (8,3)
    if degenerate:
        phys[3] = phys[2]
        phys[7] = phys[6]
    return phys


def _case_hex():
    """Non-degenerate hex, mild generic shear -- NOT axis-aligned, so the
    off-diagonal cofactor terms are actually exercised (an axis-aligned box
    makes them exactly 0 on both sides and hides a transposed-cofactor bug --
    this port's own history names exactly that failure mode)."""
    h = DX / 2.0
    T = np.array([[h, 0.03 * h, 0.02 * h],
                  [0.04 * h, h, -0.03 * h],
                  [0.01 * h, -0.02 * h, h]])
    origin = np.array([1000.0, 2000.0, -5000.0])
    xl = _hex_from_T(T, origin, degenerate=False)
    return dict(name='hex', elemType=1, xl=xl)


def _case_wedge_degenerate():
    """Degenerate wedge (elemType=11), generic asymmetric skew -- distinct
    from _case_hex's T so this is an independent geometric check."""
    h = DX / 2.0
    T = np.array([[h, -0.05 * h, 0.03 * h],
                  [0.02 * h, h, 0.04 * h],
                  [-0.03 * h, 0.01 * h, h]])
    origin = np.array([-3000.0, 500.0, -8000.0])
    xl = _hex_from_T(T, origin, degenerate=True)
    return dict(name='wedge_degenerate', elemType=11, xl=xl)


def _case_wedge_shallow_dip():
    """Degenerate wedge (elemType=12) built by rotating a cube DIP_DEG about
    the strike (x) axis -- the off-diagonal-Jacobian shape a shallow
    dipping-fault wedge produces (tpv36's par.dip=15)."""
    dip = np.radians(DIP_DEG)
    c, s = np.cos(dip), np.sin(dip)
    Rx = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    T = Rx @ (DX / 2.0 * np.eye(3))
    origin = np.array([0.0, 9000.0, -4000.0])
    xl = _hex_from_T(T, origin, degenerate=True)
    return dict(name='wedge_shallow_dip15', elemType=12, xl=xl)


CASES_GEOM = (_case_hex(), _case_wedge_degenerate(), _case_wedge_shallow_dip())

VP, VS, RHO = 6000.0, 3464.0, 2670.0   # tpv36's material (par.vp/vs/rou)
MIU = RHO * VS ** 2
LAM = RHO * VP ** 2 - 2.0 * MIU
MATE = np.array([VP, VS, RHO, LAM, MIU])

DZ = DX * np.sin(np.radians(DIP_DEG))
DT = 0.5 * DZ / VP           # tpv36's own par.dt formula
RDAMPK = 0.1                 # defaultParameters.py's rdampk default

for _c in CASES_GEOM:
    _c['mate'] = MATE
    _c['vl'] = RNG.normal(scale=1.0e-2, size=(8, 3))     # m/s, node-major
    _c['dl'] = RNG.normal(scale=1.0e-3, size=(8, 3))     # m
    _c['stress0'] = RNG.normal(scale=5.0e6, size=6)      # Pa, prior stress
    _c['dt'] = DT
    _c['rdampk'] = RDAMPK


def _fmt(arr):
    return ' '.join(repr(float(v)) for v in np.asarray(arr).ravel())


def build_stdin():
    lines = [str(len(CASES_GEOM))]
    for c in CASES_GEOM:
        lines.append(str(c['elemType']))
        lines.append(_fmt(c['xl']))
        lines.append(_fmt(c['mate']))
        lines.append(_fmt(c['vl']))
        lines.append(_fmt(c['dl']))
        lines.append(_fmt(c['stress0']))
        lines.append('%r %r' % (float(c['dt']), float(c['rdampk'])))
    return '\n'.join(lines) + '\n'


# ---------------------------------------------------------------------------
# Fortran side: fresh build + probe compile/link/run
# ---------------------------------------------------------------------------
def build_fresh_fortran_and_probe():
    env = dict(os.environ, MACHINE='ubuntu')
    subprocess.run(['bash', '-c', 'rm -f *.o eqdyna'], cwd=FSRC)
    t0 = time.time()
    r = subprocess.run(['make', 'MACHINE=ubuntu'], cwd=FSRC, env=env,
                        capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: make MACHINE=ubuntu exited %d' % r.returncode)
    print('built src/eqdyna fresh in %.1fs' % (time.time() - t0))

    probe_o = os.path.join(FSRC, '_probe_wedge_kernel.o')
    probe_bin = os.path.join(FSRC, '_probe_wedge_kernel')
    r = subprocess.run(['mpif90', '-c', '-fopenmp', '-ffree-line-length-none', '-O3',
                         '-I/usr/include', PROBE_SRC, '-o', probe_o],
                        cwd=FSRC, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: compiling probe_wedge_kernel.f90 exited %d' % r.returncode)

    objs = sorted(f for f in os.listdir(FSRC)
                  if f.endswith('.o') and f not in ('eqdyna3d.o', os.path.basename(probe_o)))
    r = subprocess.run(['mpif90', '-fopenmp', '-ffree-line-length-none', '-O3', probe_o]
                        + objs + ['-o', probe_bin, '-L/usr/lib/x86_64-linux-gnu',
                                  '-lnetcdf', '-lnetcdff'],
                        cwd=FSRC, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-4000:])
        print(r.stderr[-4000:])
        raise SystemExit('FAIL: linking probe_wedge_kernel exited %d' % r.returncode)
    return probe_bin


def run_fortran(probe_bin):
    stdin_text = build_stdin()
    r = subprocess.run(['mpirun', '-np', '1', probe_bin], input=stdin_text,
                        capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        raise SystemExit('FAIL: probe_wedge_kernel exited %d\nSTDOUT:\n%s\nSTDERR:\n%s'
                         % (r.returncode, r.stdout[-4000:], r.stderr[-4000:]))
    blocks = []
    lines = r.stdout.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == 'BEGIN_CASE':
            j = i + 1
            nums = []
            while lines[j].strip() != 'END_CASE':
                nums.extend(float(t) for t in lines[j].replace('D', 'E').split())
                j += 1
            blocks.append(np.array(nums))
            i = j + 1
        else:
            i += 1
    if len(blocks) != len(CASES_GEOM):
        raise SystemExit('FAIL: expected %d cases from probe, got %d'
                         % (len(CASES_GEOM), len(blocks)))
    return blocks


def unpack_fortran_block(flat):
    """Slice one case's flat token stream into named arrays, in the EXACT
    order probe_wedge_kernel.f90's single list-directed WRITE emits them,
    with gfortran's column-major array linearisation for each piece."""
    p = 0

    def take(n, shape=None):
        nonlocal p
        v = flat[p:p + n]
        p += n
        return v.reshape(shape, order='F') if shape else v

    det = take(1)[0]
    gsf = take(32, (4, 8))          # globalShapeFunc(4,8): gsf[row,node]
    xs = take(9, (3, 3))            # xs(3,3)
    elmass = take(24)               # elmass(24): node-major, ned-minor
    ss = take(6)
    phi = take(32, (8, 4))          # phi(node,mode)
    elresf = take(24)
    hg_force = take(24)
    assert p == flat.shape[0], (p, flat.shape[0])
    node_mass = elmass[0::3]        # elmass(3*(j-1)+1) == node j's mass (x,y,z equal)
    return dict(det=det, gsf=gsf, xs=xs, node_mass=node_mass, ss=ss, phi=phi,
                elresf=elresf, hg_force=hg_force)


# ---------------------------------------------------------------------------
# Python side: the REAL production functions, one synthetic interior element
# ---------------------------------------------------------------------------
def run_python(c):
    etype = np.array([c['elemType']])
    xl = c['xl'][None, :, :]        # (1,8,3)
    mat = c['mate'][None, :]        # (1,5)

    det, eleshp3, xs = assembleGlobalMass.compute_element_shape(xl, etype)
    ss, phi48 = assembleGlobalMass.compute_hourglass(xl, xs, mat, eleshp3)

    is_wedge = c['elemType'] in (11, 12)
    node_mass = (assembleGlobalMass._contm_wedge_node_mass(det, mat[:, 2])
                 if is_wedge else
                 np.repeat(assembleGlobalMass.contm(det, mat[:, 2])[:, None], 8, axis=1))

    N, E = 8, 1
    conn = np.arange(8, dtype=np.int64)[None, :]           # (1,8) 0-indexed
    eq_ids = np.zeros((N, 12), dtype=np.int64)
    for i in range(N):
        eq_ids[i, 0:3] = [3 * i + 1, 3 * i + 2, 3 * i + 3]  # identity map, matches Fortran
    NEQ = 24
    ndof = np.full(N, 3, dtype=np.int64)
    meshCoor = c['xl']

    S = dict(
        N=N, dt=c['dt'], rdampk=c['rdampk'], w=assembleGlobalMass._W,
        conn=conn, elemType=etype, mat=mat,
        eledet=det, eleshp=np.transpose(eleshp3, (0, 2, 1)),
        ss=ss, phi=np.transpose(phi48, (0, 2, 1)),
        ndof=ndof, eq_ids=eq_ids, NEQ=NEQ, meshCoor=meshCoor,
        PMLb=np.array([1.0e4, -1.0e4, 1.0e4, -1.0e4, -1.0e4, 500.0, 500.0, 500.0]),
        nPML=6, vmaxPML=1500.0, R=0.01,
        C_elastic=1, grav=9.8, roumax=2670.0, gamar=0.0, rhow=1000.0,
    )
    inv = assembleGlobalKU.build(S)
    scratch = assembleGlobalKU.alloc_scratch(np, inv)

    velArr = c['vl']
    dispArr = c['dl']
    force = np.zeros(NEQ + 1)
    stress_i = c['stress0'][None, :].copy()
    s_p = np.zeros((0, 15))

    force, stress_i, s_p = assembleGlobalKU.assembleGlobalKU(
        np, inv, velArr, force, stress_i, s_p, c['dt'], c['rdampk'], scratch)
    elresf = force[1:25].copy()   # identity map -> local (node,dof) order directly

    force2 = np.zeros(NEQ + 1)
    force2 = assembleGlobalKU.calcHourglassResist(np, inv, dispArr, velArr, force2, c['rdampk'])
    hg_force = force2[1:25].copy()

    gsf_full = np.zeros((4, 8))
    gsf_full[0:3, :] = eleshp3[0]
    row4 = np.full(8, 1.0 / 8.0)
    if is_wedge:
        row4[2] = 0.25
        row4[3] = 0.0
        row4[6] = 0.25
        row4[7] = 0.0
    gsf_full[3, :] = row4

    return dict(det=float(det[0]), gsf=gsf_full, xs=xs[0], node_mass=node_mass[0],
                ss=ss[0], phi=phi48[0], elresf=elresf, hg_force=hg_force)


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def compare_field(name, f, p, fails, tag):
    f = np.asarray(f, dtype=np.float64)
    p = np.asarray(p, dtype=np.float64)
    if f.shape != p.shape:
        fails.append('%s.%s: SHAPE mismatch fortran=%r python=%r' % (tag, name, f.shape, p.shape))
        return
    diff = np.abs(f - p)
    tol = ATOL + RTOL * np.abs(f)
    bad = diff > tol
    worst = float(diff.max()) if diff.size else 0.0
    print('    %-16s max|diff|=%.3e  (rtol=%.0e atol=%.0e)  %s'
         % (name, worst, RTOL, ATOL, 'OK' if not bad.any() else 'MISMATCH'))
    if bad.any():
        idx = int(np.argmax(diff))
        fails.append('%s.%s: max|diff|=%.6e at flat index %d (fortran=%.9e python=%.9e)'
                     % (tag, name, worst, idx, f.flat[idx], p.flat[idx]))


def main():
    if shutil.which('gfortran') is None or shutil.which('mpif90') is None:
        print('FAIL evidence_wedge_kernel: gfortran/mpif90 not on PATH -- this '
              'probe requires a real Fortran compile, it does not skip silently')
        return 1

    probe_bin = build_fresh_fortran_and_probe()
    fort_blocks = run_fortran(probe_bin)

    all_fails = []
    for c, flat in zip(CASES_GEOM, fort_blocks):
        tag = '%s(elemType=%d)' % (c['name'], c['elemType'])
        print('')
        print('==== case: %s ====' % tag)
        fb = unpack_fortran_block(flat)
        pb = run_python(c)
        print('    det: fortran=%.9e python=%.9e' % (fb['det'], pb['det']))
        compare_field('det', [fb['det']], [pb['det']], all_fails, tag)
        compare_field('gsf(shape-func)', fb['gsf'], pb['gsf'], all_fails, tag)
        compare_field('xs', fb['xs'], pb['xs'], all_fails, tag)
        compare_field('node_mass', fb['node_mass'], pb['node_mass'], all_fails, tag)
        compare_field('ss(hourglass)', fb['ss'], pb['ss'], all_fails, tag)
        compare_field('phi(hourglass)', fb['phi'], pb['phi'], all_fails, tag)
        compare_field('elresf(interior)', fb['elresf'], pb['elresf'], all_fails, tag)
        compare_field('hg_force', fb['hg_force'], pb['hg_force'], all_fails, tag)

    print('')
    if all_fails:
        print('FAIL evidence_wedge_kernel')
        for f in all_fails:
            print('  -', f)
        return 1
    print('SUCCESS evidence_wedge_kernel (%d cases: %s), fortran oracle == python '
          'port to rtol=%g/atol=%g on det/shape-func/xs/node_mass/ss/phi/'
          'elresf/hg_force' % (len(CASES_GEOM), ', '.join(c['name'] for c in CASES_GEOM),
                               RTOL, ATOL))
    return 0


if __name__ == '__main__':
    sys.exit(main())
