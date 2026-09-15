#! /usr/bin/env python3
"""Bit-for-bit oracle test for the VECTORIZED standalone setup path.

The setup builders in python/eqdyna/standalone/{meshgen,assembleGlobalMass}.py
were originally verbatim scalar ports of the Fortran and are now vectorized
for speed. Each vectorized builder keeps its original scalar loop alongside
it (`_build_elements_scalar`, `_build_equation_numbers_scalar`,
`_assemble_mass_scalar`); this test asserts the two agree BYTE FOR BYTE
(`np.array_equal` on ints, `tobytes()` equality on floats -- NOT allclose)
on a real mesh read from the real parity fixture, not a toy grid.

This is deliberately NOT a substitute for
testsys/parity/test_standalone_meshgen.py, which anchors the same builders
on the Fortran pydump oracle. This one catches vectorization drift; that one
catches divergence from Fortran. Both are required.

Run: python3 testsys/unit/test_standalone_setup_vectorized.py [fixture_name]
Raises loudly (non-zero exit) on any mismatch or missing fixture.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python'))

from eqdyna import meshgen, assembleGlobalMass
from eqdyna.eqdyna3d.readInputFiles import build_params, read_bmaterial

FIXTURE_NAME = sys.argv[1] if len(sys.argv) > 1 else 'test_tpv8_serial'
FIXTURE = os.path.join(HERE, 'fixtures', FIXTURE_NAME)
if not os.path.isdir(FIXTURE):
    raise FileNotFoundError(
        'missing fixture dir %s -- run testsys/parity/make_fixtures.py first' % FIXTURE)


def same_bytes(name, a, b):
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    if a.shape != b.shape:
        raise AssertionError('%s: shape %r (vectorized) vs %r (scalar oracle)'
                              % (name, a.shape, b.shape))
    if a.dtype != b.dtype:
        raise AssertionError('%s: dtype %r (vectorized) vs %r (scalar oracle)'
                              % (name, a.dtype, b.dtype))
    if a.tobytes() != b.tobytes():
        bad = np.nonzero(np.asarray(a != b))
        first = tuple(int(i[0]) for i in bad)
        raise AssertionError('%s: NOT bit-identical to the scalar oracle -- %d of %d '
                              'entries differ, first at %r: vectorized=%r scalar=%r'
                              % (name, int(np.count_nonzero(np.asarray(a != b))), a.size,
                                 first, a[first], b[first]))
    print('  %-22s bit-identical  %s %s' % (name, a.shape, a.dtype))


def main():
    params, g = build_params(FIXTURE)
    material = read_bmaterial(os.path.join(FIXTURE, 'bMaterial.txt'),
                               g['nmat'], g['n2mat'])
    xline, yline, zline, pmlb, _ = meshgen.build_grid_lines(params)
    nx, ny, nz = len(xline), len(yline), len(zline)
    print('fixture %s: grid %dx%dx%d (%d nodes, %d elements)'
          % (FIXTURE_NAME, nx, ny, nz, nx * ny * nz,
             (nx - 1) * (ny - 1) * (nz - 1)))

    # ---- on_fault_grid_mask vs the scalar is_on_fault, on EVERY node ----
    mask = meshgen.on_fault_grid_mask(xline, yline, zline, params)
    scalar_mask = np.array(
        [meshgen.is_on_fault(xline[ix], yline[iy], zline[iz],
                              params['fxmin'], params['fxmax'], params['fymin'],
                              params['fymax'], params['fzmin'], params['fzmax'],
                              params['tol'])
         for ix in range(nx) for iz in range(nz) for iy in range(ny)])
    same_bytes('on_fault mask', mask, scalar_mask)

    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)

    # ---- build_elements ----
    conn, etype, mat, depth = meshgen.build_elements(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    s_conn, s_etype, s_mat, s_depth = meshgen._build_elements_scalar(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    same_bytes('conn', conn, s_conn)
    same_bytes('elemType', etype, s_etype)
    same_bytes('mat', mat, s_mat)
    same_bytes('elem depth', depth, s_depth)

    # Element CENTERS are internal to build_elements (they only decide
    # elemType/material), so an ulp shift there would not show up in the
    # returned arrays on this mesh -- check them directly, since that is
    # exactly the reduction whose order numpy is free to choose.
    centers_vec = meshCoor[conn].mean(axis=1)
    centers_scalar = np.empty_like(centers_vec)
    for e in range(s_conn.shape[0]):
        centers_scalar[e] = meshCoor[s_conn[e]].mean(axis=0)
    same_bytes('element centers', centers_vec, centers_scalar)

    # ---- build_equation_numbers ----
    num_dof, eq_start, eq_nums, total_eqs = meshgen.build_equation_numbers(
        xline, yline, zline, params, pmlb)
    s_num_dof, s_eq_start, s_eq_nums, s_total = meshgen._build_equation_numbers_scalar(
        xline, yline, zline, params, pmlb)
    if total_eqs != s_total:
        raise AssertionError('totalNumOfEquations: vectorized %d vs scalar %d'
                              % (total_eqs, s_total))
    same_bytes('num_dof', num_dof, s_num_dof)
    same_bytes('eq_start', eq_start, s_eq_start)
    same_bytes('eq_nums (flat)',
               np.concatenate(eq_nums[1:]), np.concatenate(s_eq_nums[1:]))
    print('  %-22s vectorized %d == scalar %d' % ('total_eqs', total_eqs, s_total))

    # ---- assemble_mass ----
    xl = meshCoor[conn]
    det = assembleGlobalMass.compute_element_det(xl)
    n_nodes = meshCoor.shape[0] - 1
    nm, fn = assembleGlobalMass.assemble_mass(conn, mat, det, num_dof, eq_start, eq_nums,
                                          total_eqs, n_nodes)
    s_nm, s_fn = assembleGlobalMass._assemble_mass_scalar(
        s_conn, s_mat, det, s_num_dof, s_eq_start, s_eq_nums, s_total, n_nodes)
    same_bytes('fnms', fn, s_fn)
    same_bytes('nodalMassArr', nm, s_nm)

    # ---- pack_eq_ids vs the per-node loop it replaced ----
    ndof0, eq_ids = meshgen.pack_eq_ids(num_dof, eq_nums, n_nodes, ncols=12)
    s_eq_ids = np.zeros((n_nodes, 12), dtype=np.int64)
    s_ndof0 = np.zeros(n_nodes, dtype=np.int64)
    for node in range(1, n_nodes + 1):
        nd = int(s_num_dof[node])
        s_ndof0[node - 1] = nd
        eqs = np.array(s_eq_nums[node], dtype=np.int64)
        s_eq_ids[node - 1, :nd] = np.where(eqs > 0, eqs, 0)
    same_bytes('ndof (0-indexed)', ndof0, s_ndof0)
    same_bytes('eq_ids', eq_ids, s_eq_ids)

    print('SUCCESS test_standalone_setup_vectorized (%s): every vectorized setup '
          'builder is bit-identical to its scalar oracle' % FIXTURE_NAME)


if __name__ == '__main__':
    main()
