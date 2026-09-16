"""Unit cover for the C_degen>3 wedge-degeneration mesh port added to
meshgen.py (library_degeneration.f90's wedge()/reorder() +
meshgen.f90:91-101's type-13 retag + meshgen.f90:763-767's checkIsOnFault
C_degen>3 branch).

The full-scale, real-data parity check against a freshly-built Fortran
binary lives in testsys/parity/evidence_c_degen_port.py (test.tpv36,
report-only, not wired into any gate). This file is the FAST unit tier:
it never shells out, builds no case directory, and runs in well under a
second -- but it is NOT vacuous. `_small_c_degen36_params` reproduces
test.tpv36's own conforming relationship (dy=dx*cos(dip), dz=dx*sin(dip))
at a much smaller domain, and the tests below assert wedge/retag elements
actually occur, not just that the code runs without raising.

Rule 2 note: `test_build_elements_vectorized_matches_scalar_oracle`
asserts byte-equality (np.array_equal), not a tolerance -- this is the
scalar-vs-vectorized invariant, not a physics comparison.
"""
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))

from eqdyna import meshgen  # noqa: E402


def _small_c_degen36_params(dip=15.0, dx=500.0, n_dip_steps=6):
    """A small, genuinely conforming C_degen>3 mesh: same dx/dy/dz
    relationship test.tpv36's user_defined_params.py uses
    (dy=dx*cos(dip), dz=dx*sin(dip)), just a much smaller fault width so
    the whole mesh is a few thousand elements instead of ~950k. This
    relationship is what makes element centers land EXACTLY (to roundoff)
    on the dipping fault plane at a handful of brick positions -- confirmed
    by the assertions in the tests below, not assumed."""
    dy = dx * np.cos(dip * np.pi / 180.0)
    dz = dx * np.sin(dip * np.pi / 180.0)
    fault_width = n_dip_steps * dx
    fymax = fault_width * np.cos(dip * np.pi / 180.0)
    fzmin = -fault_width * np.sin(dip * np.pi / 180.0)
    dis4uniB = int(round(fault_width / dx)) + 3
    dis4uniF = 2
    params = dict(
        dx=dx, dy=dy, dz=dz,
        fxmin=-2 * dx, fxmax=2 * dx, fzmin=fzmin, fzmax=0.0,
        fymin=0.0, fymax=fymax,
        dis4uniF=dis4uniF, dis4uniB=dis4uniB,
        xmin=-3 * dx, xmax=3 * dx, ymin=-2 * dy, ymax=(dis4uniB + 2) * dy,
        zmin=fzmin - 3 * dz, zmax=0.0,
        rat=1.025, nPML=6, tol=1.0e-5, R=0.01,
        fstrike=0.0, C_degen=dip, insertFaultType=0, rough=None,
    )
    return params


@pytest.fixture(scope='module')
def small_mesh():
    params = _small_c_degen36_params()
    material = np.array([[6000.0, 3464.0, 2670.0]])
    xline, yline, zline, pmlb, _ = meshgen.build_grid_lines(params)
    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)
    return dict(params=params, material=material, xline=xline, yline=yline,
                zline=zline, pmlb=pmlb, meshCoor=meshCoor, nftnd=nftnd, nsmp=nsmp)


class TestIsOnFaultCDegenBranch:
    """Direct arithmetic checks of the new is_on_fault/`_check_is_on_fault_vec`
    C_degen>3 branch (meshgen.f90:763-767) -- no mesh construction needed."""

    def test_point_exactly_on_dipping_plane_is_on_fault(self):
        dip = 15.0
        dx = 500.0
        y, z = 300.0, -300.0 * np.tan(dip * np.pi / 180.0)
        assert meshgen.is_on_fault(
            0.0, y, z, fxmin=-1000, fxmax=1000, fymin=-1000, fymax=1000,
            fzmin=-1000, fzmax=1000, tol=1e-5, c_degen=dip, dx=dx) is True

    def test_point_far_from_dipping_plane_is_not_on_fault(self):
        dip = 15.0
        assert meshgen.is_on_fault(
            0.0, 300.0, 0.0, fxmin=-1000, fxmax=1000, fymin=-1000, fymax=1000,
            fzmin=-1000, fzmax=1000, tol=1e-5, c_degen=dip, dx=500.0) is False

    def test_outside_box_is_not_on_fault_even_on_plane(self):
        dip = 15.0
        y, z = 300.0, -300.0 * np.tan(dip * np.pi / 180.0)
        assert meshgen.is_on_fault(
            5000.0, y, z, fxmin=-1000, fxmax=1000, fymin=-1000, fymax=1000,
            fzmin=-1000, fzmax=1000, tol=1e-5, c_degen=dip, dx=500.0) is False

    def test_0_lt_c_degen_le_3_raises(self):
        with pytest.raises(NotImplementedError):
            meshgen.is_on_fault(0.0, 0.0, 0.0, -1, 1, -1, 1, -1, 1, 1e-5, c_degen=2.0)

    def test_c_degen_gt_3_without_dx_raises(self):
        with pytest.raises(ValueError):
            meshgen.is_on_fault(0.0, 0.0, 0.0, -1, 1, -1, 1, -1, 1, 1e-5, c_degen=15.0)


class TestOnFaultGridMaskAgreesWithScalar:
    """`on_fault_grid_mask`'s C_degen>3 branch must agree with `is_on_fault`
    on EVERY node of a real (small) grid, not a sample (Pattern 1 in this
    discipline's charter: self-consistency checks must cover the whole
    function space)."""

    def test_c_degen_gt_3_matches_scalar_everywhere(self, small_mesh):
        p = small_mesh['params']
        xline, yline, zline = small_mesh['xline'], small_mesh['yline'], small_mesh['zline']
        mask = meshgen.on_fault_grid_mask(xline, yline, zline, p)
        expect = np.array([
            meshgen.is_on_fault(xline[ix], yline[iy], zline[iz], p['fxmin'], p['fxmax'],
                                 p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], p['tol'],
                                 p['C_degen'], p['dx'])
            for ix in range(len(xline)) for iz in range(len(zline)) for iy in range(len(yline))
        ])
        assert np.array_equal(mask, expect)
        # non-vacuous: this mesh must actually contain fault nodes, or the
        # check above is trivially true for the wrong reason.
        assert mask.sum() > 0

    def test_c_degen_0_unaffected_by_the_refactor(self):
        """C_degen==0 (tpv8/tpv104/tpv10/drv.a6's branch) is bit-for-bit
        unchanged by routing on_fault_grid_mask through the shared
        `_check_is_on_fault_vec` helper."""
        p = dict(_small_c_degen36_params())
        p['C_degen'] = 0.0
        xline = np.linspace(-1000, 1000, 5)
        yline = np.array([-10.0, 0.0, 10.0])
        zline = np.linspace(-1000, 0, 5)
        mask = meshgen.on_fault_grid_mask(xline, yline, zline, p)
        expect = np.array([
            meshgen.is_on_fault(xline[ix], yline[iy], zline[iz], p['fxmin'], p['fxmax'],
                                 p['fymin'], p['fymax'], p['fzmin'], p['fzmax'], p['tol'])
            for ix in range(len(xline)) for iz in range(len(zline)) for iy in range(len(yline))
        ])
        assert np.array_equal(mask, expect)
        assert mask.sum() > 0


class TestBuildElementsWedgeSplit:
    """build_elements (vectorized) vs _build_elements_scalar (the bit-for-
    bit oracle) on the small conforming C_degen>3 mesh. This mirrors the
    file's own established pattern (every other milestone in meshgen.py
    is checked this way) for the ONE new branch that didn't have a unit
    test at all before this port (the removed parity tier that used to
    cover it, test_meshgen_vectorized.py, no longer exists -- see
    testsys/parity/evidence_c_degen_port.py for the real-data,
    freshly-built-Fortran-binary parity gate that supersedes it)."""

    def test_vectorized_matches_scalar_oracle_bit_for_bit(self, small_mesh):
        m = small_mesh
        conn_v, et_v, mat_v, depth_v = meshgen.build_elements(
            m['xline'], m['yline'], m['zline'], m['params'], m['pmlb'],
            m['nsmp'], m['material'], m['meshCoor'])
        conn_s, et_s, mat_s, depth_s = meshgen._build_elements_scalar(
            m['xline'], m['yline'], m['zline'], m['params'], m['pmlb'],
            m['nsmp'], m['material'], m['meshCoor'])

        assert np.array_equal(conn_v, conn_s)
        assert np.array_equal(et_v, et_s)
        assert np.array_equal(mat_v, mat_s)
        assert np.array_equal(depth_v, depth_s)

        # non-vacuous: the wedge split and the type-13 retag must BOTH
        # actually fire on this mesh, or the byte-equality check above
        # passes for the trivial reason that nothing new ran.
        assert int(np.sum(et_v == 11)) > 0
        assert int(np.sum(et_v == 12)) > 0
        assert int(np.sum(et_v == 11)) == int(np.sum(et_v == 12))
        assert int(np.sum(et_v == 13)) > 0

    def test_element_count_grows_by_one_per_wedge_trigger(self, small_mesh):
        """E = n_brick + n_wedge_triggers (one extra element per trigger,
        the array-level mirror of wedge() incrementing Fortran's elemCount
        by 2 instead of createElement's usual 1)."""
        m = small_mesh
        nx, ny, nz = len(m['xline']), len(m['yline']), len(m['zline'])
        n_brick = (nx - 1) * (ny - 1) * (nz - 1)
        conn_v, et_v, _, _ = meshgen.build_elements(
            m['xline'], m['yline'], m['zline'], m['params'], m['pmlb'],
            m['nsmp'], m['material'], m['meshCoor'])
        n_triggers = int(np.sum(et_v == 11))
        assert conn_v.shape[0] == n_brick + n_triggers

    def test_type11_wedge_has_matching_corner_3_4(self, small_mesh):
        """meshgen.f90's own sanity check (assembleGlobalMass.f90:45-50):
        a degenerate wedge's nodes 3 and 4 must coincide by construction."""
        m = small_mesh
        conn_v, et_v, _, _ = meshgen.build_elements(
            m['xline'], m['yline'], m['zline'], m['params'], m['pmlb'],
            m['nsmp'], m['material'], m['meshCoor'])
        wedge11 = conn_v[et_v == 11]
        assert wedge11.shape[0] > 0
        assert np.array_equal(wedge11[:, 2], wedge11[:, 3])

    def test_type11_never_gets_slave_master_substitution(self, small_mesh):
        """meshgen.f90:733-736's replaceSlaveWithMasterNode guard: type-11
        NEVER gets the substitution, type-12/13 ALWAYS do (unconditionally)
        -- this is the OR-clause that was dead code for C_degen==0 and is
        exercised for the first time by this port. A type-11 element's
        corner ids must therefore be disjoint from nsmp's MASTER column
        (nsmp[:,1]) unless a master id coincidentally also appears as a
        regular grid node id, which cannot happen since master ids are
        appended after all regular node ids."""
        m = small_mesh
        conn_v, et_v, _, _ = meshgen.build_elements(
            m['xline'], m['yline'], m['zline'], m['params'], m['pmlb'],
            m['nsmp'], m['material'], m['meshCoor'])
        master_ids = set(m['nsmp'][:, 1].tolist())
        wedge11_nodes = set(conn_v[et_v == 11].ravel().tolist())
        assert wedge11_nodes.isdisjoint(master_ids)


class TestUnsupportedCDegenRaisesLoudly:
    """0<C_degen<=3: the Fortran's checkIsOnFault itself takes neither
    if/elseif branch there (isOnFault stays 0 unconditionally) -- this
    port refuses rather than silently reproducing that degenerate
    behavior. Boundary test per this discipline's new-feature checklist."""

    def test_build_elements_raises_for_0_lt_c_degen_le_3(self, small_mesh):
        m = small_mesh
        bad_params = dict(m['params'])
        bad_params['C_degen'] = 2.0
        with pytest.raises(NotImplementedError):
            meshgen.build_elements(m['xline'], m['yline'], m['zline'], bad_params,
                                    m['pmlb'], m['nsmp'], m['material'], m['meshCoor'])

    def test_build_elements_scalar_raises_for_0_lt_c_degen_le_3(self, small_mesh):
        m = small_mesh
        bad_params = dict(m['params'])
        bad_params['C_degen'] = 2.0
        with pytest.raises(NotImplementedError):
            meshgen._build_elements_scalar(m['xline'], m['yline'], m['zline'], bad_params,
                                            m['pmlb'], m['nsmp'], m['material'], m['meshCoor'])


class TestEqdyna3dGuards:
    """eqdyna3d.py's build_solver_state guards for the C_degen>3 combos
    this port does NOT support end to end (Flag-interaction tests per this
    discipline's new-feature checklist: C_degen>3 interacting with
    C_elastic==0, and with a mesh that actually contains wedge elements)."""

    def test_c_degen_gt_3_plus_c_elastic_0_refused(self, monkeypatch):
        sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))
        from eqdyna import eqdyna3d, readInputFiles

        def fake_build_params(case_dir):
            return dict(_small_c_degen36_params(), insertFaultType=0), dict(
                ntotft=1, friclaw=1, npx=1, npy=1, npz=1, C_elastic=0)
        monkeypatch.setattr(readInputFiles, 'build_params', fake_build_params)
        with pytest.raises(NotImplementedError, match='C_elastic'):
            eqdyna3d.build_solver_state('/nonexistent')

    def test_wedge_elements_no_longer_refused_by_compute_element_shape(self, monkeypatch):
        """UPDATED: this port's dynamics guard (`if np.any(elem_type in
        (11,12)): raise NotImplementedError(...)`, formerly sitting right
        before assembleGlobalMass.compute_element_shape in
        build_solver_state) is GONE -- the wedge-degeneration dynamics are
        now ported (see assembleGlobalMass.py/assembleGlobalKU.py's own
        docstrings and testsys/parity/evidence_wedge_kernel.py, the
        Fortran-vs-Python kernel probe this is a cheap unit-level echo of).
        `compute_element_shape` must now ACCEPT elem_type 11/12 and return a
        finite, positive determinant, not raise."""
        import numpy as _np
        from eqdyna import assembleGlobalMass

        # A minimal degenerate wedge: local node 4 collapsed onto node 3,
        # node 8 onto node 7 (the same node pairing calcGlobalShapeFunc.f90's
        # merge acts on), otherwise a generic hex.
        local = _np.array([
            [-1.0, -1.0, -1.0], [1.0, -1.0, -1.0], [1.0, 1.0, -1.0], [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 1.0], [-1.0, 1.0, 1.0],
        ]) * 250.0
        local[3] = local[2]
        local[7] = local[6]
        xl = local[None, :, :]
        for etype in (11, 12):
            elem_type = _np.array([etype])
            det, eleshp, xs = assembleGlobalMass.compute_element_shape(xl, elem_type)
            assert _np.isfinite(det).all() and (det > 0).all()
            assert eleshp.shape == (1, 3, 8) and xs.shape == (1, 3, 3)
            # calcGlobalShapeFunc.f90:22-28's merge zeroes local columns 4/8
            # (0-indexed 3/7) of the derivative rows for a wedge element.
            assert _np.array_equal(eleshp[0, :, 3], _np.zeros(3))
            assert _np.array_equal(eleshp[0, :, 7], _np.zeros(3))
