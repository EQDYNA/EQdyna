#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the TPV29/TPV30 rough-fault
normal: a fault-normal direction must be the derivative of the surface AT THE
MESH'S OWN SPACING, never a subsample of a derivative taken at a finer one.

THE INCIDENT (2026-09-22). `case_input/test.tpv29/tpv29GeometryTools.py`
decimated the official SCEC surface correctly -- `y[::step]` keeps official
values, no interpolation, rule 17 step 2 -- and then decimated the SUPPLIED
DERIVATIVE COLUMNS the same way, in both call sites (`decimateToDx` and
`convertOfficial25m`). That is a category error. `y` is a sampling of a
surface; `dF/dx` is not a value of that surface but a property of a GRID AT A
SPACING. Verified against SCEC's own file: their dF/dx, dF/dy columns ARE
`np.gradient(F, 25)` -- median relative difference 0.0215% / 0.0304%, corr
1.000000 over all 1601x801 nodes -- which is precisely why they cannot be
carried to a coarser grid. At the gate's dx = 500 m the shipped normals came
from a grid sampled 20x finer than the mesh they were attached to.

WHY IT MATTERS, AND WHY NOTHING CAUGHT IT. EQdyna hands `dy/dx`, `dy/dz` to
`func_lib.f90:insertFaultInterface` (lines 112-114), which builds the per-node
un/us/ud frame at `meshgen.f90:892-905`. Under `C_elastic = 0` the on-fault
traction is RECONSTRUCTED by projecting assembled nodal forces onto that
frame (`meshgen.f90:926-937`), so a normal that is not perpendicular to the
mesh facets it bounds leaks the lithostatic MEAN stress into the SHEAR
channel, with gain 1/(R sin 2a) = 6.30 for TPV29/30's own R = 0.160739092,
a = 40.375111 deg, against 0.163 into the normal -- a 38.7:1 asymmetry.
`readInputFiles.f90:341-345` checks only |dy/dx|*dx/dy < 1, a mesh-VALIDITY
bound (does the surface climb less than one cell between nodes), which a
wrong-spacing derivative passes comfortably. Nothing asserted CONSISTENCY
between the normals and the facets they bound. That gap is what this guard
closes; broader coverage of the same question is pathway item 63.

WHAT THIS GUARD ASSERTS, in order of what it would have caught:

  1. FAITHFUL, NOT INVENTIVE. On a committed 61x61 block of the official 25 m
     file -- cut around the node carrying the largest |dF/dx| anywhere in it,
     so the identity is tested where the surface is roughest -- recomputing
     with `np.gradient(F, 25, edge_order=2)` reproduces SCEC's OWN supplied
     derivative columns. If it did not, the replacement would be an operation
     of our own invention and everything downstream of it would be wrong.
  2. THE TOOL'S OWN SELF-CHECK BITES. `checkDerivativesAgainstOfficial`
     accepts the correct frame mapping and REFUSES a sign-flipped one.
  3. EVERY SHIPPED SURFACE IS MESH-CONSISTENT. The derivative columns of
     `bFault_Rough_Geometry.tpv29.{100m,50m}.txt` (and TPV30's byte-identical
     copy of the 100 m file) equal `np.gradient` of their own `y` column at
     their own header dx.
  4. DECIMATION RECOMPUTES. `faultGridForCase(500)` is mesh-consistent, and
     the old subsampling behaviour is measurably different -- asserted, so
     this guard cannot pass vacuously against the bug it was written for.
  5. THE CONSEQUENCE IS GONE. The metric (rule 4b: an exclusion names the
     metric that produced it) is tau/strength at every fault node of the
     dx = 500 m gate grid, with the traction taken on the MESH FACET normal
     and projected onto the ASSIGNED un/us/ud frame -- the C_elastic = 0
     reconstruction path in miniature. Mesh-consistent normals must put ZERO
     nodes at or above strength; subsampled normals must put some there.

Exit 0 = pass, non-zero = fail, one SUCCESS/FAIL line per check.
"""
import contextlib
import importlib.util
import os
import sys

import numpy as np


@contextlib.contextmanager
def inDir(d):
    """faultGridForCase/shippedSourceForDx name the shipped surface by BARE
    filename -- the tool's contract is that it runs in its own compset
    directory, the way user_defined_params.py imports it."""
    prev = os.getcwd()
    os.chdir(d)
    try:
        yield
    finally:
        os.chdir(prev)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
TPV29_DIR = os.path.join(REPO_ROOT, 'case_input', 'test.tpv29')
TPV30_DIR = os.path.join(REPO_ROOT, 'case_input', 'test.tpv30')
FIXTURE = os.path.join(HERE, 'tpv29_official_25m_block.txt')

OFFICIAL_DX = 25.0
GATE_DX = 500.0

# Bound 1: the derivative identity against SCEC's own columns. Measured on the
# full official grid: 0.0215% (dF/dx) / 0.0304% (dF/dy) median relative;
# 0.0123% on this block's interior. 0.1% is ~5-8x the worst measurement and is
# the same bound the tool itself enforces
# (tpv29GeometryTools.DERIVATIVE_FAITHFULNESS_MEDIAN_REL).
MEDIAN_REL_BOUND = 1.0e-3

# Bound 3: the shipped files store y and its derivatives at %.7e, so a
# derivative recomputed FROM the rounded y differs from the stored one by
# roughly (1e-7 * max|y|) / dx ~ 1725e-7/50 = 3.5e-6 at worst. 1e-4 absolute
# is ~30x that and still ~4000x smaller than the 3.9e-3 the bug moved these
# columns by, so it separates "text rounding" from "wrong spacing" cleanly.
STORED_DERIV_ATOL = 1.0e-4

# TPV29/TPV30 initial-stress constants, spec TPV29_30_Description_v06 p.12-13,
# copied from case_input/test.tpv30/user_defined_params.py (the compset is the
# authority; these are here so the check is standalone).
ROU, GRAV, MU_S = 2670., 9.8, 0.18
B11, B33, B13 = 1.025837, 0.974162, -0.158649

failures = []


def report(name, ok, detail):
    print('%s: %s - %s' % (name, 'SUCCESS' if ok else 'FAIL', detail))
    if not ok:
        failures.append(name)


def loadGeoTools(caseDir):
    """Import the compset's own copy of the tool, with scripts/ on the path
    (it imports lib). Named per directory so tpv29's and tpv30's byte-identical
    copies are two distinct modules and neither shadows the other."""
    for p in (os.path.join(REPO_ROOT, 'scripts'), caseDir):
        if p not in sys.path:
            sys.path.insert(0, p)
    path = os.path.join(caseDir, 'tpv29GeometryTools.py')
    spec = importlib.util.spec_from_file_location(
        'geoTools_' + os.path.basename(caseDir).replace('.', '_'), path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def loadFixture():
    """The committed official block, returned in the EQDYNA frame at 25 m:
    (y, dy/dx, dy/dz) shaped (nz, nx) with iz increasing UPWARD.

    The TPV frame has depth increasing with row index, so flipping axis 0
    gives z_eq increasing with index; that flip also negates a gradient taken
    along it, which IS the dy/dz_eq = -dF/dy_tpv mapping -- so the supplied
    column is negated once and the array flipped once, and a gradient of the
    flipped surface then matches it directly with no further sign."""
    d = np.loadtxt(FIXTURE)
    n = int(round(np.sqrt(d.shape[0])))
    if n*n != d.shape[0]:
        raise AssertionError('%s: %d rows is not a square block' % (FIXTURE, d.shape[0]))
    F = d[:, 0].reshape(n, n)
    dFdx = d[:, 1].reshape(n, n)
    dFdy = d[:, 2].reshape(n, n)
    return F[::-1], dFdx[::-1], -dFdy[::-1]


def medianRel(a, b):
    m = np.abs(b) > 1e-12
    return float(np.median(np.abs(a[m] - b[m]) / np.abs(b[m])))


# --- 1. the recomputation reproduces SCEC's own supplied columns -------------
def check_scec_derivative_identity(g):
    y, supX, supZ = loadFixture()
    recX, recZ = g.meshConsistentDerivatives(y, OFFICIAL_DX)
    # Interior only: np.gradient's edges are one-sided over the BLOCK, while
    # SCEC's are centered over the full grid. Comparing them would measure the
    # block's boundary, not the scheme.
    inner = (slice(1, -1), slice(1, -1))
    rx = medianRel(recX[inner], supX[inner])
    rz = medianRel(recZ[inner], supZ[inner])
    report('scec_derivative_identity', max(rx, rz) <= MEDIAN_REL_BOUND,
           'np.gradient(F, 25, edge_order=2) reproduces the official columns to '
           'median relative %.4f%% (dy/dx) / %.4f%% (dy/dz) on %d interior '
           'nodes, bound %.3f%%. SCEC\'s derivatives ARE gradients of their own '
           'surface at 25 m, which is why subsampling them is wrong.'
           % (100*rx, 100*rz, recX[inner].size, 100*MEDIAN_REL_BOUND))


# --- 2. the tool's own self-check accepts the right frame and refuses a wrong one
def check_tool_self_check_bites(g):
    y, supX, supZ = loadFixture()
    from lib import FaultGeometryError
    try:
        worst = g.checkDerivativesAgainstOfficial(y, supX, supZ, OFFICIAL_DX)
    except FaultGeometryError as exc:
        report('tool_self_check_accepts_correct_frame', False,
               'the correct frame mapping was REFUSED: %s' % exc)
        return
    report('tool_self_check_accepts_correct_frame', worst <= MEDIAN_REL_BOUND,
           'checkDerivativesAgainstOfficial accepts the official block, worst '
           'median relative %.4f%%' % (100*worst))

    raised = False
    try:
        # dy/dz WITHOUT the z_eq = -y_tpv sign flip: the single most likely way
        # to get this wrong, and the one that would silently mirror every
        # fault normal about the horizontal.
        g.checkDerivativesAgainstOfficial(y, supX, -supZ, OFFICIAL_DX)
    except FaultGeometryError:
        raised = True
    report('tool_self_check_refuses_sign_flip', raised,
           'a dy/dz sign flip is %s' % ('refused loudly (rule 2)' if raised else
                                        'ACCEPTED -- the self-check is vacuous'))


# --- 3. every shipped surface is consistent with its own header dx -----------
def check_shipped_surfaces(g):
    shipped = [(TPV29_DIR, 'bFault_Rough_Geometry.tpv29.100m.txt'),
               (TPV29_DIR, 'bFault_Rough_Geometry.tpv29.50m.txt'),
               (TPV30_DIR, 'bFault_Rough_Geometry.tpv29.100m.txt')]
    worst, worstName = 0.0, None
    for d, name in shipped:
        path = os.path.join(d, name)
        dx, y, dydx, dydz = g.loadEQdynaGeometry(path)
        recX, recZ = g.meshConsistentDerivatives(y, dx)
        e = max(float(np.abs(recX - dydx).max()), float(np.abs(recZ - dydz).max()))
        if e > worst:
            worst, worstName = e, '%s/%s' % (os.path.basename(d), name)
    report('shipped_surfaces_mesh_consistent', worst <= STORED_DERIV_ATOL,
           'over %d shipped surfaces the stored dy/dx, dy/dz columns match '
           'np.gradient of their own y at their own header dx to max|diff| '
           '%.3e (worst: %s), bound %.1e'
           % (len(shipped), worst, worstName, STORED_DERIV_ATOL))


# --- 4. decimation to the gate dx recomputes rather than subsamples ----------
def check_decimation_recomputes(g):
    with inDir(TPV29_DIR):
        y, dydx, dydz = g.faultGridForCase(GATE_DX)
    recX, recZ = g.meshConsistentDerivatives(y, GATE_DX)
    e = max(float(np.abs(recX - dydx).max()), float(np.abs(recZ - dydz).max()))
    ok = e == 0.0
    report('decimation_recomputes', ok,
           'faultGridForCase(%g) returns derivatives identical to np.gradient '
           'of its own surface at %g m (max|diff| %.3e)' % (GATE_DX, GATE_DX, e))

    # The guard must be able to FAIL: show the old behaviour is far outside it.
    srcDx, fname = g.shippedSourceForDx(GATE_DX)
    _, sY, sX, sZ = g.loadEQdynaGeometry(os.path.join(TPV29_DIR, fname))
    step = int(round(GATE_DX/srcDx))
    subX, subZ = sX[::step, ::step], sZ[::step, ::step]
    gap = max(float(np.abs(subX - recX).max()), float(np.abs(subZ - recZ).max()))
    report('decimation_discriminates', gap > 100*STORED_DERIV_ATOL,
           'subsampling the %g m derivative columns by %d (the bug) differs '
           'from recomputing at %g m by max %.4f -- %.0fx the tolerance this '
           'guard uses, so a regression could not hide inside it'
           % (srcDx, step, GATE_DX, gap, gap/STORED_DERIV_ATOL))


# --- 5. the consequence: no fault node at or above strength at the gate dx ---
def _tauOverStrength(y, assignedX, assignedZ, facetX, facetZ, dx):
    """tau/strength with the traction taken on the MESH FACET normal and
    resolved in the ASSIGNED un/us/ud frame -- the C_elastic = 0 FE
    reconstruction in miniature. When the two frames agree this reduces to the
    ordinary analytic resolution; when they disagree the mean stress leaks
    into the shear channel, which is the defect."""
    nz, nx = y.shape
    depth = np.abs(-20000.0 + dx*np.arange(nz))[:, None]*np.ones((1, nx))
    dEff = np.maximum(depth, dx/2.)
    omega = np.where(dEff <= 17000., 1.0, np.maximum(0.0, (22000. - dEff)/5000.))
    s22 = -(ROU - 1000.)*GRAV*dEff
    S = np.zeros(y.shape + (3, 3))
    S[..., 0, 0] = (omega*B11 + (1. - omega))*s22
    S[..., 1, 1] = (omega*B33 + (1. - omega))*s22
    S[..., 2, 2] = s22
    S[..., 0, 1] = S[..., 1, 0] = omega*B13*s22

    def frame(px, pz):
        den = np.sqrt(px**2 + 1.0 + pz**2)
        nn = np.stack([-px/den, 1.0/den, -pz/den], -1)
        dens = np.sqrt(1.0 + px**2)
        ss = np.stack([np.ones_like(px)/dens, px/dens, np.zeros_like(px)], -1)
        return nn, ss, np.cross(ss, nn)

    na, sa, da = frame(assignedX, assignedZ)
    nm, _, _ = frame(facetX, facetZ)
    t = np.einsum('...ij,...j->...i', S, nm)
    Tn = np.einsum('...i,...i->...', na, t)
    Ts = np.einsum('...i,...i->...', sa, t)
    Td = np.einsum('...i,...i->...', da, t)
    C0 = np.where(depth < 4000., 0.4e6 + 200.0*(4000. - depth), 0.4e6)
    return np.hypot(Ts, Td)/(MU_S*np.maximum(-Tn, 0.0) + C0)


def check_no_node_at_strength(g):
    with inDir(TPV29_DIR):
        y, dydx, dydz = g.faultGridForCase(GATE_DX)
    facetX, facetZ = g.meshConsistentDerivatives(y, GATE_DX)
    srcDx, fname = g.shippedSourceForDx(GATE_DX)
    _, _, sX, sZ = g.loadEQdynaGeometry(os.path.join(TPV29_DIR, fname))
    step = int(round(GATE_DX/srcDx))

    fixed = _tauOverStrength(y, dydx, dydz, facetX, facetZ, GATE_DX)[1:-1, 1:-1]
    bug = _tauOverStrength(y, sX[::step, ::step], sZ[::step, ::step],
                           facetX, facetZ, GATE_DX)[1:-1, 1:-1]
    nFixed, nBug = int((fixed >= 1.0).sum()), int((bug >= 1.0).sum())
    report('no_fault_node_at_strength', nFixed == 0 and fixed.max() < 1.0,
           'metric = tau/strength on the facet normal resolved in the assigned '
           'frame, %d interior nodes of the %g m gate grid: max %.5f, %d at or '
           'above strength' % (fixed.size, GATE_DX, fixed.max(), nFixed))
    report('at_strength_check_discriminates', nBug > 0,
           'the same metric with SUBSAMPLED normals (the bug): max %.5f, %d at '
           'or above strength -- the check is not vacuous' % (bug.max(), nBug))


def main():
    if not os.path.isfile(FIXTURE):
        print('rough_fault_normal_consistency: FAIL - fixture %s is missing; a '
              'check that cannot run must not read as a pass (rule 2)' % FIXTURE)
        return 1
    g29 = loadGeoTools(TPV29_DIR)
    g30 = loadGeoTools(TPV30_DIR)
    check_scec_derivative_identity(g29)
    check_tool_self_check_bites(g29)
    check_shipped_surfaces(g29)
    check_decimation_recomputes(g29)
    check_no_node_at_strength(g29)
    # TPV30 ships a byte-identical copy of the tool and of the 100 m surface;
    # assert that rather than re-running every check against it.
    same = (open(os.path.join(TPV29_DIR, 'tpv29GeometryTools.py'), 'rb').read()
            == open(os.path.join(TPV30_DIR, 'tpv29GeometryTools.py'), 'rb').read())
    report('tpv30_copy_is_identical', same and g30.SHIPPED_DX == g29.SHIPPED_DX,
           'case_input/test.tpv30/tpv29GeometryTools.py is byte-identical to '
           'test.tpv29\'s, so the checks above speak for both compsets')

    if failures:
        print('rough_fault_normal_consistency: FAIL - %d check(s) failed: %s'
              % (len(failures), ', '.join(failures)))
        return 1
    print('rough_fault_normal_consistency: SUCCESS - all checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
