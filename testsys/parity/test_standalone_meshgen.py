#! /usr/bin/env python3
"""
Milestone 1 + 2 + 3 + 4 + 5 + 6 + 6.5 + 7 + 7.5 parity check for
python/eqdyna/standalone/{meshgen,native_input,mass_assembly,frt_writer}.py
against the fresh tpv8 serial fixtures in
testsys/parity/fixtures/test_tpv8_serial/ (pydump_meshCoor.txt,
pydump_conn.txt, pydump_nodeinfo.txt, pydump_fault.txt,
pydump_stations.txt, pydump_nodalmass.txt, pydump_fnms.txt,
pydump_elemgeo.txt, pydump_v1.txt, frt.txt0, bGlobal.txt,
bModelGeometry.txt, bFaultGeometry.txt, bMaterial.txt, bStations.txt,
on_fault_vars_input.nc). M7.5 adds eledet/eleshp/ss/phi
(compute_element_shape/compute_hourglass) and init_vel.

Milestone 6 changes WHERE the case parameters/material/station tables/
initial fric array come from: instead of a hand-transcribed PARAMS dict
(as M1-M5 used), everything is read natively from the case-input files via
python/eqdyna/standalone/native_input.py -- the SAME M1-M5 builders and
fixtures are then re-run/re-checked against these natively-read inputs, so
a bug in native_input.py's parsing would show up as an M1-M5 regression,
not pass silently because "the builders were already proven".

Run: python3 testsys/parity/make_fixtures.py   (fresh fixtures first, rule 4)
     python3 testsys/parity/test_standalone_meshgen.py
     python3 testsys/parity/test_standalone_meshgen.py test_drv_a6_serial   (or any other
     testsys/parity/fixtures/<name> directory make_fixtures.py --case <case> produced --
     optional positional arg, defaults to test_tpv8_serial for the wired-into-`parity`-tier
     invocation above)

Raises loudly (non-zero exit) on any mismatch or missing fixture -- no
silent pass, per PROJECT_RULES rule 2/3.
"""
import os
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python'))

from eqdyna.standalone.meshgen import (
    build_grid_lines, build_node_coordinates, build_elements, build_equation_numbers,
    build_fault_geometry, build_station_matching)
from eqdyna.standalone.native_input import (
    build_params, read_bmaterial, read_bstations, read_on_fault_vars)
from eqdyna.standalone.mass_assembly import (
    compute_element_det, assemble_mass, compute_element_shape, compute_hourglass, init_vel)
from eqdyna.standalone.frt_writer import read_frt, format_frt_row

FIXTURE = os.path.join(TESTSYS, 'fixtures',
                        sys.argv[1] if len(sys.argv) > 1 else 'test_tpv8_serial')

# ---- M6: everything below is read NATIVELY from the case-input files ----
PARAMS, GLOBALS = build_params(FIXTURE)
MATERIAL = read_bmaterial(os.path.join(FIXTURE, 'bMaterial.txt'),
                           GLOBALS['nmat'], GLOBALS['n2mat'])


def load_pydump_meshCoor(path, n):
    arr = np.loadtxt(path)
    assert arr.shape == (n, 3), arr.shape
    return arr


def load_pydump_stations(path):
    with open(path) as f:
        n_onf = int(f.readline())
        n_off = int(f.readline())
        anonfs = [tuple(int(v) for v in f.readline().split()) for _ in range(n_onf)]
        off = [tuple(int(v) for v in f.readline().split()) for _ in range(n_off)]
    return anonfs, off


def load_pydump_fault(path, n):
    nsmp = np.zeros((n, 2), dtype=np.int64)
    un = np.zeros((n, 3))
    us = np.zeros((n, 3))
    ud = np.zeros((n, 3))
    arn = np.zeros(n)
    fric = np.zeros((n, 101))  # column 0 unused, 1-indexed FRIC_SLOT_* columns
    with open(path) as f:
        for i, line in enumerate(f):
            vals = line.split()
            nsmp[i] = [int(vals[0]), int(vals[1])]
            un[i] = [float(v) for v in vals[2:5]]
            us[i] = [float(v) for v in vals[5:8]]
            ud[i] = [float(v) for v in vals[8:11]]
            arn[i] = float(vals[11])
            fric[i, 1:] = [float(v) for v in vals[12:112]]
    return nsmp, un, us, ud, arn, fric


def load_pydump_elemgeo(path, n):
    """pydump_elemgeo.txt: eledet(i), ((eleshp(j,k,i),j=1,3),k=1,8) [24],
    (ss(j,i),j=1,6) [6], ((phi(j,k,i),j=1,8),k=1,4) [32] per element.
    Returns (det, eleshp, ss, phi) with eleshp (E,8,3) [node,deriv] and
    phi (E,4,8) [mode,node], matching python/eqdyna/loading.py's load()
    convention (so this fixture format's own provenance is cross-checked
    against the already-verified Fortran-port loader's reshape logic)."""
    raw = np.loadtxt(path)
    assert raw.shape == (n, 63), raw.shape
    det = raw[:, 0]
    eleshp = raw[:, 1:25].reshape(n, 8, 3)
    ss = raw[:, 25:31]
    phi = raw[:, 31:63].reshape(n, 4, 8)
    return det, eleshp, ss, phi


def load_pydump_conn(path, n):
    conn = np.zeros((n, 8), dtype=np.int64)
    etype = np.zeros(n, dtype=np.int64)
    mat = np.zeros((n, 5))
    with open(path) as f:
        for i, line in enumerate(f):
            vals = line.split()
            conn[i] = [int(v) for v in vals[0:8]]
            etype[i] = int(vals[8])
            mat[i] = [float(v) for v in vals[9:14]]
    return conn, etype, mat


def main():
    for name in ('pydump_meshCoor.txt', 'pydump_conn.txt', 'pydump_header.txt'):
        p = os.path.join(FIXTURE, name)
        if not os.path.isfile(p):
            raise FileNotFoundError(
                'missing fixture %s -- run testsys/parity/make_fixtures.py first' % p)

    with open(os.path.join(FIXTURE, 'pydump_header.txt')) as f:
        header_lines = f.readlines()
    n_nodes_fortran = int(header_lines[0])
    n_elem_fortran = int(header_lines[1])

    xline, yline, zline, pmlb, _ = build_grid_lines(PARAMS)
    meshCoor, nftnd, nsmp = build_node_coordinates(xline, yline, zline, PARAMS)

    # ---- M1: node coordinates ----
    fortran_coor = load_pydump_meshCoor(
        os.path.join(FIXTURE, 'pydump_meshCoor.txt'), n_nodes_fortran)
    py_coor = meshCoor[1:]  # drop the unused 0 row -> 1-indexed node i is row i-1
    if py_coor.shape[0] != n_nodes_fortran:
        raise AssertionError('M1 node count mismatch: python %d vs Fortran %d'
                              % (py_coor.shape[0], n_nodes_fortran))
    m1_diff = np.max(np.abs(py_coor - fortran_coor))
    print('M1 meshCoor max abs diff: %e (n=%d)' % (m1_diff, n_nodes_fortran))
    if m1_diff > 1e-9:
        raise AssertionError('M1 meshCoor parity FAILED: max abs diff %e' % m1_diff)

    # ---- M2: element connectivity + material ----
    conn, elem_type, mat, _depth = build_elements(xline, yline, zline, PARAMS, pmlb, nsmp, MATERIAL, meshCoor)
    if conn.shape[0] != n_elem_fortran:
        raise AssertionError('M2 element count mismatch: python %d vs Fortran %d'
                              % (conn.shape[0], n_elem_fortran))
    fconn, fetype, fmat = load_pydump_conn(
        os.path.join(FIXTURE, 'pydump_conn.txt'), n_elem_fortran)

    conn_mismatches = np.sum(np.any(conn != fconn, axis=1))
    etype_mismatches = np.sum(elem_type != fetype)
    # mat's 5 columns span wildly different scales depending on the case's
    # velocity model (tpv8: vp/vs/rho ~1e3; drv.a6's per-layer moduli columns
    # ~1e10) -- an absolute tolerance calibrated against tpv8's scale is
    # meaningless for drv.a6's (same audit pattern already applied to
    # M6.5's nodalMassArr/fnms and M7.5's ss below: relative diff, not
    # absolute, for a scale-dependent quantity).
    mat_absdiff = np.abs(mat - fmat)
    mat_reldiff = mat_absdiff / np.where(np.abs(fmat) > 0, np.abs(fmat), 1.0)
    mat_diff = np.max(mat_absdiff)
    mat_reldiff_max = np.max(mat_reldiff)
    print('M2 conn mismatched elements: %d / %d' % (conn_mismatches, n_elem_fortran))
    print('M2 elemType mismatched elements: %d / %d' % (etype_mismatches, n_elem_fortran))
    print('M2 mat max abs/rel diff: %e / %e' % (mat_diff, mat_reldiff_max))

    if conn_mismatches:
        bad = np.argmax(np.any(conn != fconn, axis=1))
        raise AssertionError('M2 connectivity parity FAILED at element %d: '
                              'python=%s fortran=%s' % (bad, conn[bad], fconn[bad]))
    if etype_mismatches:
        bad = np.argmax(elem_type != fetype)
        raise AssertionError('M2 elemType parity FAILED at element %d: '
                              'python=%d fortran=%d' % (bad, elem_type[bad], fetype[bad]))
    if mat_reldiff_max > 1e-9:
        bad = np.unravel_index(np.argmax(mat_reldiff), fmat.shape)
        raise AssertionError('M2 material parity FAILED at element %d slot %d: '
                              'python=%r fortran=%r' % (bad[0], bad[1], mat[bad], fmat[bad]))

    # ---- M3: equation numbering ----
    num_dof, eq_start, eq_nums, total_eqs = build_equation_numbers(xline, yline, zline, PARAMS, pmlb)
    with open(os.path.join(FIXTURE, 'pydump_header.txt')) as f:
        header_lines = f.readlines()
    n_equations_fortran = int(header_lines[2])
    if total_eqs != n_equations_fortran:
        raise AssertionError('M3 totalNumOfEquations mismatch: python %d vs Fortran %d'
                              % (total_eqs, n_equations_fortran))

    dof_mismatches = 0
    start_mismatches = 0
    eqnum_mismatches = 0
    first_bad = None
    with open(os.path.join(FIXTURE, 'pydump_nodeinfo.txt')) as f:
        for node_id, line in enumerate(f, start=1):
            vals = [int(v) for v in line.split()]
            f_ndof, f_start = vals[0], vals[1]
            f_eqs = vals[2:2 + f_ndof]
            if num_dof[node_id] != f_ndof:
                dof_mismatches += 1
                if first_bad is None:
                    first_bad = ('dof', node_id, num_dof[node_id], f_ndof)
                continue
            if eq_start[node_id] != f_start:
                start_mismatches += 1
                if first_bad is None:
                    first_bad = ('start', node_id, eq_start[node_id], f_start)
            if list(eq_nums[node_id]) != f_eqs:
                eqnum_mismatches += 1
                if first_bad is None:
                    first_bad = ('eqnums', node_id, list(eq_nums[node_id]), f_eqs)
    print('M3 numOfDofPerNodeArr mismatches: %d / %d' % (dof_mismatches, n_nodes_fortran))
    print('M3 eqNumStartIndexLoc mismatches: %d / %d' % (start_mismatches, n_nodes_fortran))
    print('M3 eqNumIndexArr mismatches: %d / %d' % (eqnum_mismatches, n_nodes_fortran))
    if dof_mismatches or start_mismatches or eqnum_mismatches:
        raise AssertionError('M3 equation-numbering parity FAILED, first bad: %r' % (first_bad,))

    # ---- M4: split-node fault geometry (un/us/ud, arn) ----
    fault_fixture = os.path.join(FIXTURE, 'pydump_fault.txt')
    if not os.path.isfile(fault_fixture):
        raise FileNotFoundError(
            'missing fixture %s -- run testsys/parity/make_fixtures.py first' % fault_fixture)
    if nftnd != nsmp.shape[0]:
        raise AssertionError('M4 nftnd mismatch: build_node_coordinates %d vs nsmp rows %d'
                              % (nftnd, nsmp.shape[0]))
    f_nsmp, f_un, f_us, f_ud, f_arn, f_fric = load_pydump_fault(fault_fixture, nftnd)

    nsmp_mismatches = np.sum(np.any(nsmp != f_nsmp, axis=1))
    if nsmp_mismatches:
        bad = np.argmax(np.any(nsmp != f_nsmp, axis=1))
        raise AssertionError('M4 nsmp parity FAILED at fault node %d: python=%s fortran=%s'
                              % (bad, nsmp[bad], f_nsmp[bad]))

    un, us, ud, arn = build_fault_geometry(xline, yline, zline, PARAMS, nsmp)
    un_diff = np.max(np.abs(un[1:] - f_un))
    us_diff = np.max(np.abs(us[1:] - f_us))
    ud_diff = np.max(np.abs(ud[1:] - f_ud))
    arn_diff = np.max(np.abs(arn[1:] - f_arn))
    print('M4 un max abs diff: %e' % un_diff)
    print('M4 us max abs diff: %e' % us_diff)
    print('M4 ud max abs diff: %e' % ud_diff)
    print('M4 arn max abs diff: %e (n=%d)' % (arn_diff, nftnd))
    if un_diff > 1e-9 or us_diff > 1e-9 or ud_diff > 1e-9:
        raise AssertionError('M4 un/us/ud parity FAILED: max abs diff %e / %e / %e'
                              % (un_diff, us_diff, ud_diff))
    if arn_diff > 1e-6:
        raise AssertionError('M4 arn parity FAILED: max abs diff %e' % arn_diff)

    # ---- M5: on-/off-fault station matching ----
    stations_fixture = os.path.join(FIXTURE, 'pydump_stations.txt')
    bstations_fixture = os.path.join(FIXTURE, 'bStations.txt')
    for p in (stations_fixture, bstations_fixture):
        if not os.path.isfile(p):
            raise FileNotFoundError(
                'missing fixture %s -- run testsys/parity/make_fixtures.py first' % p)

    xonfs, x4nds = read_bstations(bstations_fixture)
    f_anonfs, f_off = load_pydump_stations(stations_fixture)

    py_anonfs, py_off = build_station_matching(xline, yline, zline, PARAMS, xonfs, x4nds)

    print('M5 numOfOnFaultStCount: python %d vs fortran %d' % (len(py_anonfs), len(f_anonfs)))
    print('M5 numOfOffFaultStCount: python %d vs fortran %d' % (len(py_off), len(f_off)))
    if py_anonfs != f_anonfs:
        raise AssertionError('M5 on-fault station parity FAILED: python=%r fortran=%r'
                              % (py_anonfs, f_anonfs))
    if py_off != f_off:
        raise AssertionError('M5 off-fault station parity FAILED: python=%r fortran=%r'
                              % (py_off, f_off))

    # ---- M6: native netCDF4 read of on_fault_vars_input.nc -> initial fric ----
    nc_fixture = os.path.join(FIXTURE, 'on_fault_vars_input.nc')
    if not os.path.isfile(nc_fixture):
        raise FileNotFoundError(
            'missing fixture %s -- run testsys/parity/make_fixtures.py first' % nc_fixture)
    fric = read_on_fault_vars(nc_fixture, PARAMS['fxmin'], PARAMS['fzmin'],
                               PARAMS['dx'], PARAMS['dz'], meshCoor, nsmp)
    fric_diff = np.max(np.abs(fric[1:] - f_fric))
    print('M6 fric (initial on-fault state) max abs diff: %e (n=%d, 100 slots)'
          % (fric_diff, nftnd))
    if fric_diff > 1e-6:
        bad = np.unravel_index(np.argmax(np.abs(fric[1:] - f_fric)), f_fric.shape)
        raise AssertionError('M6 fric parity FAILED at fault node %d slot %d: '
                              'python=%r fortran=%r' % (bad[0], bad[1],
                                                         fric[1:][bad], f_fric[bad]))

    print('SUCCESS test_standalone_meshgen (M1 + M2 + M3 + M4 + M5 + M6 parity)')

    # ---- M6.5: lumped-mass assembly (nodalMassArr, fnms) ----
    nodalmass_fixture = os.path.join(FIXTURE, 'pydump_nodalmass.txt')
    fnms_fixture = os.path.join(FIXTURE, 'pydump_fnms.txt')
    for p in (nodalmass_fixture, fnms_fixture):
        if not os.path.isfile(p):
            raise FileNotFoundError(
                'missing fixture %s -- run testsys/parity/make_fixtures.py first' % p)

    xl = meshCoor[conn]  # (E,8,3): conn is 1-indexed, meshCoor row-0 unused
    det = compute_element_det(xl)
    nodalMassArr, fnms = assemble_mass(conn, mat, det, num_dof, eq_start, eq_nums,
                                        total_eqs, n_nodes_fortran)

    f_nodalmass = np.loadtxt(nodalmass_fixture)
    f_fnms = np.loadtxt(fnms_fixture)
    if f_nodalmass.shape[0] != total_eqs:
        raise AssertionError('M6.5 nodalMassArr size mismatch: python %d vs fortran %d'
                              % (total_eqs, f_nodalmass.shape[0]))
    if f_fnms.shape[0] != n_nodes_fortran:
        raise AssertionError('M6.5 fnms size mismatch: python %d vs fortran %d'
                              % (n_nodes_fortran, f_fnms.shape[0]))

    # Nodal mass values here are O(1e12) (kg-scale lumped mass accumulated
    # over ~1e5 m^3 elements at rho~2700 kg/m^3) -- an absolute tolerance
    # calibrated for Pa/m-scale quantities elsewhere in this port would be
    # meaningless at this magnitude; use max relative diff instead (same
    # "roundoff-only" bar, just expressed correctly for this scale).
    nodalmass_absdiff = np.abs(nodalMassArr[1:] - f_nodalmass)
    fnms_absdiff = np.abs(fnms[1:] - f_fnms)
    nodalmass_reldiff = np.max(nodalmass_absdiff / np.abs(f_nodalmass))
    fnms_reldiff = np.max(fnms_absdiff / np.abs(f_fnms))
    print('M6.5 nodalMassArr max abs/rel diff: %e / %e (n=%d)'
          % (np.max(nodalmass_absdiff), nodalmass_reldiff, total_eqs))
    print('M6.5 fnms max abs/rel diff: %e / %e (n=%d)'
          % (np.max(fnms_absdiff), fnms_reldiff, n_nodes_fortran))
    if nodalmass_reldiff > 1e-9:
        bad = np.argmax(nodalmass_absdiff / np.abs(f_nodalmass))
        raise AssertionError('M6.5 nodalMassArr parity FAILED at eq %d: python=%r fortran=%r'
                              % (bad + 1, nodalMassArr[1 + bad], f_nodalmass[bad]))
    if fnms_reldiff > 1e-9:
        bad = np.argmax(fnms_absdiff / np.abs(f_fnms))
        raise AssertionError('M6.5 fnms parity FAILED at node %d: python=%r fortran=%r'
                              % (bad + 1, fnms[1 + bad], f_fnms[bad]))

    print('SUCCESS test_standalone_meshgen (M1 + M2 + M3 + M4 + M5 + M6 + M6.5 parity)')

    # ---- M7: frt.txt writer, round-trip byte-exact vs a real Fortran file ----
    frt_fixture = os.path.join(FIXTURE, 'frt.txt0')
    if not os.path.isfile(frt_fixture):
        raise FileNotFoundError(
            'missing fixture %s -- run testsys/parity/make_fixtures.py first' % frt_fixture)

    with open(frt_fixture) as f:
        original_lines = [l.rstrip('\n') for l in f]
    parsed = read_frt(frt_fixture)
    reproduced_lines = [format_frt_row(row) for row in parsed]

    line_mismatches = sum(1 for a, b in zip(original_lines, reproduced_lines) if a != b)
    print('M7 frt.txt round-trip line mismatches: %d / %d' % (line_mismatches, len(original_lines)))
    if len(original_lines) != len(reproduced_lines):
        raise AssertionError('M7 frt.txt row count mismatch: original %d vs reproduced %d'
                              % (len(original_lines), len(reproduced_lines)))
    if line_mismatches:
        bad = next(i for i, (a, b) in enumerate(zip(original_lines, reproduced_lines)) if a != b)
        raise AssertionError('M7 frt.txt writer parity FAILED at row %d:\noriginal  =%r\n'
                              'reproduced=%r' % (bad, original_lines[bad], reproduced_lines[bad]))

    print('SUCCESS test_standalone_meshgen (M1 + M2 + M3 + M4 + M5 + M6 + M6.5 + M7 parity)')

    # ---- M7.5: eledet/eleshp (calcGlobalShapeFunc), ss/phi (calcSSPhi4Hrgls), init_vel ----
    elemgeo_fixture = os.path.join(FIXTURE, 'pydump_elemgeo.txt')
    v1_fixture = os.path.join(FIXTURE, 'pydump_v1.txt')
    for p in (elemgeo_fixture, v1_fixture):
        if not os.path.isfile(p):
            raise FileNotFoundError(
                'missing fixture %s -- run testsys/parity/make_fixtures.py first' % p)

    f_det, f_eleshp, f_ss, f_phi = load_pydump_elemgeo(elemgeo_fixture, n_elem_fortran)

    xl = meshCoor[conn]  # (E,8,3), same convention as M6.5's det computation
    det75, eleshp75, xs75 = compute_element_shape(xl)
    ss75, phi75 = compute_hourglass(xl, xs75, mat, eleshp75)
    eleshp75 = np.transpose(eleshp75, (0, 2, 1))  # (E,3,8) -> (E,8,3): [node,deriv]
    phi75 = np.transpose(phi75, (0, 2, 1))        # (E,8,4) -> (E,4,8): [mode,node]

    det_diff = np.max(np.abs(det75 - f_det))
    eleshp_diff = np.max(np.abs(eleshp75 - f_eleshp))
    phi_diff = np.max(np.abs(phi75 - f_phi))
    print('M7.5 eledet max abs diff: %e' % det_diff)
    print('M7.5 eleshp max abs diff: %e' % eleshp_diff)
    print('M7.5 phi max abs diff: %e' % phi_diff)
    if det_diff > 1e-6:
        raise AssertionError('M7.5 eledet parity FAILED: max abs diff %e' % det_diff)
    if eleshp_diff > 1e-9:
        bad = np.unravel_index(np.argmax(np.abs(eleshp75 - f_eleshp)), f_eleshp.shape)
        raise AssertionError('M7.5 eleshp parity FAILED at element %d node %d deriv %d: '
                              'python=%r fortran=%r' % (bad[0], bad[1], bad[2],
                                                         eleshp75[bad], f_eleshp[bad]))
    # ss(1..6) span ~12 orders of magnitude per element (diagonal terms
    # ss(1)/ss(4)/ss(6) ~1e12 for a PML box; cross terms ss(2)/ss(3)/ss(5)
    # are analytically ~0 for this port's axis-aligned rectangular
    # elements and are themselves catastrophic-cancellation noise on both
    # sides -- an absolute tolerance calibrated to the diagonal scale would
    # be meaningless for the cross terms, and vice versa (same pattern as
    # M6.5's nodalMassArr/fnms relative-diff note). Scale each element's
    # row by its own max magnitude (the "co" prefactor scale shared by all
    # 6 slots of that element) instead.
    ss_scale = np.max(np.abs(f_ss), axis=1, keepdims=True)
    ss_reldiff = np.abs(ss75 - f_ss) / ss_scale
    ss_diff = np.max(ss_reldiff)
    print('M7.5 ss max per-element-scaled relative diff: %e' % ss_diff)
    if ss_diff > 1e-9:
        bad = np.unravel_index(np.argmax(ss_reldiff), f_ss.shape)
        raise AssertionError('M7.5 ss parity FAILED at element %d slot %d: python=%r fortran=%r'
                              % (bad[0], bad[1], ss75[bad], f_ss[bad]))
    if phi_diff > 1e-9:
        bad = np.unravel_index(np.argmax(np.abs(phi75 - f_phi)), f_phi.shape)
        raise AssertionError('M7.5 phi parity FAILED at element %d mode %d node %d: '
                              'python=%r fortran=%r' % (bad[0], bad[1], bad[2],
                                                         phi75[bad], f_phi[bad]))

    f_v1 = np.loadtxt(v1_fixture)
    if f_v1.shape[0] != total_eqs:
        raise AssertionError('M7.5 v1 size mismatch: python %d vs fortran %d'
                              % (total_eqs, f_v1.shape[0]))
    v1_py = init_vel(nsmp, eq_nums, fric, total_eqs)
    v1_diff = np.max(np.abs(v1_py[1:] - f_v1))
    print('M7.5 init_vel (v1) max abs diff: %e (n=%d)' % (v1_diff, total_eqs))
    if v1_diff > 1e-9:
        bad = np.argmax(np.abs(v1_py[1:] - f_v1))
        raise AssertionError('M7.5 init_vel parity FAILED at eq %d: python=%r fortran=%r'
                              % (bad + 1, v1_py[1 + bad], f_v1[bad]))

    print('SUCCESS test_standalone_meshgen '
          '(M1 + M2 + M3 + M4 + M5 + M6 + M6.5 + M7 + M7.5 parity)')


if __name__ == '__main__':
    main()
