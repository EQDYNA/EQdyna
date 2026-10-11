#! /usr/bin/env python3
"""
Defect A (resolution README, docs/evidence/barall-overlays-2026-10-09/
resolution/, owner decision 2026-10-11 board PR #178): an element must be
typed PML (type 2) by the centre of its GRID-LINE cell, the same unblended
coordinates setNumDof uses to give a node 12 dofs -- not by the mean of its
y-blended meshCoor.

insertFaultInterface's y-blend (func_lib.f90) stretches the whole footwall
from the y-min PML face to a dipping fault whose y-offset (`peak`) reaches
7.5 km on test.tpv10/12/13. Typing elements by the blended centre made
regular (type 1) elements in the y-min PML interface own 12-dof nodes: their
y/z force went into split x-dofs and the Method-2 initial stress there was
unbalanced. Measured at dx 500 m, nt==1 residual (EQDYNA_DUMP_EQUIL=1): up
to ~90 m/s^2 on one slanted node column before, 2.9e-2 m/s^2 (the generic
free-surface value) after. Outcome: TPV13 vs FaultMod at 8 s, slip ratio
0.61 -> 0.95, body peak-velocity ratio 4.24 -> 1.01, p90 |dt| 2.84 -> 0.53 s.

RULE 14a, both ways:
  PASS  the port's build_elements types a coarse test.tpv13 mesh with ZERO
        regular elements holding a 12-dof node, and meshgen.f90's
        createElement tests gridCenterCoor against PMLb.
  FAIL  the same invariant checker, fed the historical blended-centre typing
        of the SAME mesh, must find > 0 such elements; and the Fortran
        checker, fed a synthetic copy of createElement testing
        elementCenterCoor, must reject it.

Cheap (rule 9): no solver run. One case.setup at dx 1500 m plus the port's
mesh builders, a few seconds.
"""
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
FORTRAN_SRC = os.path.join(ROOT, 'src', 'fortran', 'meshgen.f90')
CASE = 'test.tpv13'
COARSE_DX = 1500.0

FORTRAN_TEST_RE = re.compile(r'if\s*\(\s*(\w+)\(1\)\s*>\s*PMLb\(1\)')


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'), env.get('PATH', '')])
    return env


def _coarse_case(tmp):
    case_dir = os.path.join(tmp, 'case')
    subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                    case_dir, CASE], env=_env(), check=True, capture_output=True)
    p = os.path.join(case_dir, 'user_defined_params.py')
    text = open(p).read()
    text, n = re.subn(r'^par\.dx = [0-9.]+', 'par.dx = %r' % COARSE_DX, text, flags=re.M)
    if n != 1:
        raise RuntimeError('could not override par.dx in %s' % p)
    open(p, 'w').write(text)
    subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=_env(),
                   check=True, capture_output=True)
    return case_dir


def mixed_elements(conn, is_pml_elem, num_dof):
    """Number of regular elements that hold at least one 12-dof node."""
    has12 = (num_dof[conn] == 12).any(axis=1)
    return int((~is_pml_elem & has12).sum())


def check_python():
    from eqdyna import readInputFiles, meshgen
    with tempfile.TemporaryDirectory() as tmp:
        cd = _coarse_case(tmp)
        params, g = readInputFiles.build_params(cd)
        material = readInputFiles.read_bmaterial(os.path.join(cd, 'bMaterial.txt'),
                                                 g['nmat'], g['n2mat'])
    xl, yl, zl, pmlb, _ = meshgen.build_grid_lines(params)
    mc, _, nsmp = meshgen.build_node_coordinates(xl, yl, zl, params)
    conn, et, _, _ = meshgen.build_elements(xl, yl, zl, params, pmlb, nsmp, material, mc)
    nd = np.asarray(meshgen.build_equation_numbers(xl, yl, zl, params, pmlb)[0])
    if nd.shape[0] == mc.shape[0] - 1:
        nd = np.concatenate([[0], nd])
    conn = np.asarray(conn)
    n12 = int((nd == 12).sum())
    new_bad = mixed_elements(conn, np.asarray(et) == 2, nd)
    c = mc[conn].mean(axis=1)  # historical: blended-centre typing
    old_pml = ((c[:, 0] > pmlb['xmax0']) | (c[:, 0] < pmlb['xmin0']) |
               (c[:, 1] > pmlb['ymax0']) | (c[:, 1] < pmlb['ymin0']) |
               (c[:, 2] < pmlb['zmin0']))
    old_bad = mixed_elements(conn, old_pml, nd)
    print('  python: %d elements, %d 12-dof nodes; regular elements holding a '
          '12-dof node: fixed=%d, historical blended-centre=%d'
          % (len(et), n12, new_bad, old_bad))
    fails = []
    if n12 == 0:
        fails.append('no 12-dof nodes at all -- the check tested nothing')
    if new_bad != 0:
        fails.append('build_elements leaves %d regular elements holding a 12-dof node' % new_bad)
    if old_bad == 0:
        fails.append('the checker found 0 mixed elements under the historical typing too '
                     '-- it cannot come out both ways on this mesh')
    return fails


def check_fortran_text(text):
    names = FORTRAN_TEST_RE.findall(text)  # setNumDof's nodeCoor test also matches
    if 'elementCenterCoor' in names:
        return False, 'an element is typed PML by its blended elementCenterCoor'
    if 'gridCenterCoor' not in names:
        return False, 'createElement grid-centre PML test not found (PMLb tests on: %s)' % names
    return True, 'ok (PMLb tests on: %s)' % ', '.join(names)


def check_fortran():
    text = open(FORTRAN_SRC).read()
    fails = []
    ok, msg = check_fortran_text(text)
    print('  fortran: live meshgen.f90 -> %s' % msg)
    if not ok:
        fails.append(msg)
    bad = text.replace('if (gridCenterCoor(1)>PMLb(1)', 'if (elementCenterCoor(1)>PMLb(1)')
    ok_bad, msg_bad = check_fortran_text(bad)
    print('  fortran: synthetic blended-centre copy -> %s' % msg_bad)
    if ok_bad:
        fails.append('the Fortran checker accepted the synthetic blended-centre copy')
    return fails


def main():
    fails = check_python() + check_fortran()
    if fails:
        print('\nFAIL test_pml_element_type_grid_centre (%d)' % len(fails))
        for f in fails:
            print(' -', f)
        return 1
    print('\nSUCCESS test_pml_element_type_grid_centre')
    return 0


if __name__ == '__main__':
    sys.exit(main())
