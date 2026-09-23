#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the GENERATED fault surface:
`scripts/generateFaultInterface` with `par.insertFaultType = 2` (fractal) and
= 1 (planar dipping) must store derivative columns that are the derivative of
the surface column IT WROTE, taken at THAT SURFACE'S OWN dx and dz.

WHY A SEPARATE FILE FROM test_rough_fault_normal_consistency.py. That guard is
the incident record for the SHIPPED SCEC surface (702b56b): its fixture, its
constants (OFFICIAL_DX, GATE_DX, TPV29/30's R, a, B11/B33/B13) and all six of
its checks are TPV29/TPV30's, and it reaches them through `case_input/` on
disk. This guard answers the same QUESTION -- is a derivative column the
derivative of its own grid? -- on the other path into that question, the one
where EQdyna itself produces the surface. It shares no fixture, no constant and
no import with the TPV29 guard (this one needs matplotlib and a module-level
`par`), and when one fails the file name alone says which path broke. Pathway
item 71.

WHAT WAS UNCOVERED. The generated path is exercised today only as a SIDE
EFFECT: `scripts/case.setup:379-380` calls `ensureFaultRoughGeometryForCase`,
which reaches `lib._checkDerivativeColumn` (lib.py:816-817) on both columns.
That is a real gate, but it fires only when a case is set up, and nothing in
`testsys/` asserted it. `test.drv.a6` and `test.drv.a6.v2` -- the two gated
cases carrying `insertFaultType = 2` -- therefore depended on an assertion that
lived entirely outside the test system.

WHAT THIS ASSERTS, and how it is stronger than the gate it backs up.
`lib._checkDerivativeColumn` compares against a TRUNCATION-ERROR bound implied
by the data (4x a neighbourhood maximum), because a supplied file may carry the
ANALYTIC derivative. A GENERATED file has no such freedom: the generator writes
`np.gradient` itself (generateFaultInterface:94-95), so the stored columns must
equal `np.gradient` of the stored surface to the file's write precision and
nothing looser. This guard asserts that equality directly, at an absolute
1e-4 -- the same bound the TPV29 guard uses for its shipped files, and ~100x
the %f write precision of these files.

  1. GENERATED DERIVATIVES ARE GRADIENTS AT THEIR OWN SPACING. Both columns of
     a freshly generated insertFaultType=2 surface equal np.gradient of its own
     y column at its own dx (axis 1) and dz (axis 0). dx != dz deliberately, so
     a swap of the two spacings cannot pass.
  2. THE FRACTAL BRANCH ACTUALLY PRODUCES ROUGHNESS. A plane satisfies check 1
     trivially, so check 1 alone would pass on a dead insertFaultType=2 branch.
     Asserted: type 2's surface is rough, type 1's (at dip 90) is flat, and the
     two differ.
  3. THE CHECK DISCRIMINATES. The same columns compared against a gradient
     taken at the WRONG spacing (2dx / 2dz -- the defect class 702b56b paid
     for) miss by orders of magnitude more than the bound, so this guard cannot
     pass vacuously.
  4. THE PRODUCTION GATE ACCEPTS THE GENERATED FILE, without case.setup:
     `lib.validateFaultRoughGeometryForCase` -- the exact call at
     case.setup:379-380 -- run directly on the generated file and its
     provenance sidecar.
  5. CHECK 1 BITES. One stored dy/dx value moved by 0.05 in a copy of that
     same file must put check 1's residual outside its bound. An assertion not
     shown to fail on the defect it names is not an assertion.
  6. THE PRODUCTION GATE BITES TOO, at its own looser bound. Measured by
     bisection on this surface (2026-09-23): the smallest single-node dy/dx
     error `validateFaultRoughGeometryForCase` refuses at this node is 0.0568,
     which is 12.6% of the surface's own max|slope| of 0.4505 -- the
     data-implied truncation bound behaving as lib.py:556 documents, on a
     surface deliberately near the generator's own 0.2 roughness warning. It
     therefore does NOT catch the 0.05 of check 5. That is the whole argument
     for this guard: the generated path can be held to an absolute 1e-4, ~570x
     tighter, because the generator writes np.gradient itself.

STAND-IN, disclosed: `generateFaultInterface` does `from user_defined_params
import par` at import time, so this guard installs a minimal synthetic case in
sys.modules -- a 32x16 fault grid, dx=100, dz=50, dip=90, seedId fixed at 1 by
the script itself. It is a real par object in every attribute the generator and
the validator read; it is synthetic only in being small. Everything is written
under tempfile.mkdtemp() and removed; no in-repo path is written (item 70).

Exit 0 = pass, non-zero = fail, one SUCCESS/FAIL line per check.
"""
import os
import shutil
import sys
import tempfile
import types

os.environ.setdefault('MPLBACKEND', 'Agg')  # headless; must precede pyplot

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')
GENERATOR = os.path.join(SCRIPTS, 'generateFaultInterface')
GEOM_FILE = 'bFault_Rough_Geometry.txt'

# The generator writes with fmt='%f' -- 6 decimals -- so a derivative recomputed
# from the rounded y differs from the stored, also-rounded, column by at most
# ~5e-7 (the column's own rounding) + 5e-7/min(dx,dz) (the surface's, propagated
# through the difference) ~ 5.1e-7 here. 1e-4 is ~200x that, and is the same
# bound test_rough_fault_normal_consistency.py uses on the shipped files.
STORED_DERIV_ATOL = 1.0e-4

# A small, cheap fault grid. nz <= nx: the generator cuts the (nx, nx) fractal
# field down to [:nz, :nx].
NX, NZ = 32, 16
DX, DZ = 100.0, 50.0        # deliberately unequal: a dx/dz swap must not pass
DIP = 90.0                  # vertical, so the dip ramp is ~0 and type 1 is flat
FXMIN, FZMIN = -1600.0, -800.0

failures = []


def report(name, ok, detail):
    print('%s: %s - %s' % (name, 'SUCCESS' if ok else 'FAIL', detail))
    if not ok:
        failures.append(name)


def makePar(insertFaultType):
    """The minimal case the generator and lib's validator read. Plain
    attributes, no defaults machinery: every name below is one the code under
    test actually reads, and a missing one must raise AttributeError there."""
    par = types.SimpleNamespace(
        insertFaultType=insertFaultType, dip=DIP,
        nfx=NX, nfz=NZ, dx=DX, dz=DZ, dy=DX,
        fxmin=FXMIN, fxmax=FXMIN + (NX - 1)*DX,
        fzmin=FZMIN, fzmax=FZMIN + (NZ - 1)*DZ,
        ymin=-5000.0, ymax=5000.0)
    return par


def loadGenerator(par):
    """Import scripts/generateFaultInterface (no .py suffix) with `par` in
    place, fresh each time so the module-level `from user_defined_params import
    par` binds to THIS case."""
    import importlib.machinery
    import importlib.util
    stub = types.ModuleType('user_defined_params')
    stub.par = par
    sys.modules['user_defined_params'] = stub
    sys.modules.pop('generateFaultInterface', None)
    if SCRIPTS not in sys.path:
        sys.path.insert(0, SCRIPTS)
    loader = importlib.machinery.SourceFileLoader('generateFaultInterface',
                                                  GENERATOR)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


def generateInto(workDir, insertFaultType):
    """Run the generator's own entry point (so the insertFaultType dispatch at
    generateFaultInterface:174 is what selects the branch) in workDir, and read
    the file back with lib's reader. Returns (par, header, y, dydx, dydz)."""
    import lib
    par = makePar(insertFaultType)
    prev = os.getcwd()
    os.chdir(workDir)
    try:
        loadGenerator(par)._main_func('regression guard')
        header, y, dydx, dydz = lib.readFaultRoughGeometry(GEOM_FILE)
    finally:
        os.chdir(prev)
    return par, header, y, dydx, dydz


def mutate(work, name, srcDir, row, delta):
    """Copy a generated case and move ONE stored dy/dx value by delta.
    Returns the copy's directory. The original is never touched."""
    dst = os.path.join(work, name)
    shutil.copytree(srcDir, dst)
    path = os.path.join(dst, GEOM_FILE)
    with open(path) as f:
        lines = f.readlines()
    cols = lines[row].split()
    cols[1] = '%f' % (float(cols[1]) + delta)
    lines[row] = '\t'.join(cols) + '\n'
    with open(path, 'w') as f:
        f.writelines(lines)
    return dst


def main():
    sys.path.insert(0, SCRIPTS)
    import lib
    from lib import FaultGeometryError

    if not os.path.isfile(GENERATOR):
        print('fractal_fault_geometry_derivatives: FAIL - %s is missing; a '
              'check that cannot run must not read as a pass (rule 2)'
              % GENERATOR)
        return 1

    work = tempfile.mkdtemp(prefix='eqdyna-item71-')
    try:
        fractalDir = os.path.join(work, 'fractal')
        planeDir = os.path.join(work, 'plane')
        os.makedirs(fractalDir)
        os.makedirs(planeDir)

        par, header, y, dydx, dydz = generateInto(fractalDir, 2)

        # --- 1. the columns are gradients of the surface at its own dx, dz ---
        recX = np.gradient(y, DX, axis=1)
        recZ = np.gradient(y, DZ, axis=0)
        eX = float(np.abs(dydx - recX).max())
        eZ = float(np.abs(dydz - recZ).max())
        report('generated_derivatives_at_own_spacing',
               max(eX, eZ) <= STORED_DERIV_ATOL,
               'insertFaultType=2 surface %dx%d at dx=%g, dz=%g: stored dy/dx '
               'matches np.gradient(y, dx, axis=1) to max|diff| %.3e and '
               'stored dy/dz matches np.gradient(y, dz, axis=0) to %.3e, '
               'bound %.1e (dx != dz, so a spacing swap cannot pass)'
               % (header['nnx'], header['nnz'], header['dx'], DZ, eX, eZ,
                  STORED_DERIV_ATOL))

        # --- 2. the fractal branch is alive; a plane would pass check 1 -----
        _, _, yPlane, pX, pZ = generateInto(planeDir, 1)
        rough = float(max(np.abs(dydx).max(), np.abs(dydz).max()))
        flat = float(max(np.abs(pX).max(), np.abs(pZ).max()))
        report('fractal_branch_produces_roughness',
               rough > 100*STORED_DERIV_ATOL and flat <= STORED_DERIV_ATOL
               and float(np.abs(y - yPlane).max()) > 1.0,
               'insertFaultType=2 gives max|slope| %.4g and peak |y| %.4g m, '
               'insertFaultType=1 at dip %g gives max|slope| %.3e and a flat '
               'surface -- so check 1 is not being satisfied by a dead branch '
               'returning a plane' % (rough, np.abs(y).max(), DIP, flat))

        # --- 3. the identity discriminates: wrong spacing is far outside ----
        gap = max(float(np.abs(np.gradient(y, 2*DX, axis=1) - dydx).max()),
                  float(np.abs(np.gradient(y, 2*DZ, axis=0) - dydz).max()))
        report('wrong_spacing_is_discriminated', gap > 100*STORED_DERIV_ATOL,
               'differentiating the same surface at TWICE its spacing (the '
               'wrong-spacing defect class) misses the stored columns by max '
               '%.4g = %.0fx this bound, so a regression could not hide inside '
               'it' % (gap, gap/STORED_DERIV_ATOL))

        # --- 4. the production gate accepts it, without case.setup ----------
        prev = os.getcwd()
        os.chdir(fractalDir)
        try:
            diag = lib.validateFaultRoughGeometryForCase(par, verbose=False)
            accepted, why = True, ('max|slope| %.4g, max per-cell offset %.4g '
                                   'cell' % (diag['maxSlope'], diag['maxOffset']))
        except FaultGeometryError as exc:
            accepted, why = False, str(exc)
        finally:
            os.chdir(prev)
        report('production_gate_accepts_generated_surface', accepted,
               'lib.validateFaultRoughGeometryForCase (the call at '
               'case.setup:379-380) accepts the generated file and its '
               'provenance sidecar with no case.setup run: %s' % why)

        # --- 5. the identity itself bites on one perturbed stored value -----
        # Mutation baked in: an assertion that has not been shown to fail on
        # the defect it names is not an assertion.
        row = 2 + (NX//2)*NZ + NZ//2          # an interior node, z-fastest
        mut05 = mutate(work, 'mut05', fractalDir, row, 0.05)
        _, my, mdydx, _ = lib.readFaultRoughGeometry(
            os.path.join(mut05, GEOM_FILE))
        resid = float(np.abs(mdydx - np.gradient(my, DX, axis=1)).max())
        report('identity_refuses_perturbed_derivative',
               resid > STORED_DERIV_ATOL,
               'one stored dy/dx value moved +0.05 at data row %d puts the '
               'check-1 residual at %.4g = %.0fx the %.1e bound, so check 1 '
               'fails on it' % (row - 2, resid, resid/STORED_DERIV_ATOL,
                                STORED_DERIV_ATOL))

        # --- 6. the production gate bites too, at its own (looser) bound ----
        # MEASURED on this surface, by bisection, 2026-09-23: the smallest
        # single-node dy/dx error validateFaultRoughGeometryForCase refuses at
        # this node is 0.0568 -- 12.6% of the surface's max|slope| of 0.4505.
        # That is the data-implied truncation bound doing its job (lib.py:556
        # explains why a SUPPLIED file needs it), and it is exactly why check 1
        # exists: a GENERATED file writes np.gradient itself, so it can be held
        # to an absolute 1e-4 instead, ~570x tighter here. 0.2 is 3.5x the
        # measured threshold.
        mut20 = mutate(work, 'mut20', fractalDir, row, 0.2)
        os.chdir(mut20)
        try:
            lib.validateFaultRoughGeometryForCase(par, verbose=False)
            raisedMsg = None
        except FaultGeometryError as exc:
            raisedMsg = str(exc).splitlines()[0]
        finally:
            os.chdir(prev)
        report('production_gate_refuses_perturbed_derivative',
               raisedMsg is not None,
               'one stored dy/dx value moved +0.2 at data row %d is %s'
               % (row - 2,
                  ('refused loudly (rule 2): %s' % raisedMsg) if raisedMsg
                  else 'ACCEPTED -- the gate is vacuous on generated files'))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if failures:
        print('fractal_fault_geometry_derivatives: FAIL - %d check(s) failed: %s'
              % (len(failures), ', '.join(failures)))
        return 1
    print('fractal_fault_geometry_derivatives: SUCCESS - all checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
