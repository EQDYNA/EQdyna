"""Tests for the PY-ONLY parallel paths: the shard_map split
(EQDYNA_JAX_DEVICES > 1) and the rank bookkeeping of the python-jax-mpi
decomposition (driver.run_mpi under mpirun).

Neither path is reached by the serial sweep columns, so they are tested here.

WHAT IS AND IS NOT COVERED HERE. Pure index algebra and arithmetic: the
Fortran partition formulas (meshgen.partition_1d / calc_xyz_mpi_id against
hand-derived tables), bitwise line slicing, face ordering, ownership, and the
loud refusals. That a rank's rank-local MESH is the serial mesh's box is a
real-case check and lives in testsys/regression/test_rank_local_mesh.py; that
an N-rank run lands inside the case bound is the e2e python-jax-mpi cell. A
unit test that "passed" on internal consistency while the physics moved is
exactly the failure mode this repo's rule 1 exists for.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, os.pardir,
                                'src', 'python'))
from eqdyna import MPI4NodalQuant as MQ   # noqa: E402
from eqdyna import backend as B           # noqa: E402
from eqdyna import meshgen                 # noqa: E402


def synthetic(nx=6, ny=3, nz=3, npml=1):
    """A structured brick mesh with a PML shell and a two-sided fault, as the
    index arrays the partition consumes. NOT a physics case: every float is a
    placeholder, because the shard_map split only ever SELECTS rows and never
    reads a value. The shapes and the connectivity are what is under test."""
    nnx, nny, nnz = nx + 1, ny + 1, nz + 1
    N = nnx * nny * nnz
    nid = np.arange(N).reshape(nnx, nny, nnz)
    conn = []
    etype = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                conn.append([nid[i, j, k], nid[i + 1, j, k], nid[i + 1, j + 1, k],
                             nid[i, j + 1, k], nid[i, j, k + 1], nid[i + 1, j, k + 1],
                             nid[i + 1, j + 1, k + 1], nid[i, j + 1, k + 1]])
                edge = (i < npml or i >= nx - npml or j < npml or j >= ny - npml)
                etype.append(2 if edge else 1)
    conn = np.array(conn)
    elemType = np.array(etype)
    E = conn.shape[0]
    # 3 dof everywhere except PML nodes (12) -- the shape velDispUpdate branches on
    ndof = np.full(N, 3)
    pml_nodes = np.unique(conn[elemType == 2])
    ndof[pml_nodes] = 12
    eq_ids = np.zeros((N, 12), dtype=np.int64)
    nxt = 1
    for n in range(N):
        for d in range(ndof[n]):
            eq_ids[n, d] = nxt
            nxt += 1
    NEQ = nxt - 1
    Ei = int((elemType != 2).sum()); Ep = int((elemType == 2).sum())
    S = dict(N=N, NEQ=NEQ, conn=conn, elemType=elemType, eq_ids=eq_ids, ndof=ndof)

    def col(n, *shape):
        return np.arange(n * int(np.prod(shape or (1,)))).reshape((n,) + shape) * 1.0

    inv = {}
    for k, g in B._ELEM_GROUP.items():
        n = {'Ei': Ei, 'Ep': Ep, 'E': E}[g]
        if k in B._RAVELLED:
            # BUILT FROM eq_ids AND conn, exactly as assembleGlobalKU.build
            # does it -- NOT `arange % (NEQ+1)`, which is what these were.
            # An arbitrary equation number is not merely unrealistic: it
            # addresses an equation the rank does not hold, so it cannot be
            # renumbered into the local set at all, and decompose refuses it
            # (correctly). The fixture has to satisfy the same postcondition
            # the real mesh does: every index a rank's elements produce
            # belongs to a node that rank's elements touch.
            ci = conn[elemType != 2]
            cp = conn[elemType == 2]
            if k == 'idxP12':
                is12 = ndof[cp] == 12
                inv[k] = [np.where(is12, eq_ids[cp, jj], 0).ravel()
                          for jj in range(12)]
            elif k == 'idxP3':
                is3 = ndof[cp] == 3
                inv[k] = [np.where(is3, eq_ids[cp, gg], 0).ravel()
                          for gg in range(3)]
            elif k in ('idxIx', 'idxIy', 'idxIz'):
                inv[k] = eq_ids[ci, 'xyz'.index(k[-1])].ravel()
            else:                                   # idxH0 / idxH1 / idxH2
                d = int(k[-1])
                slot = np.where(ndof[conn] == 3, d, 9 + d)
                inv[k] = np.take_along_axis(eq_ids[conn], slot[:, :, None],
                                            axis=2)[:, :, 0].ravel()
        elif k in ('conn', 'conn_i', 'conn_p'):
            inv[k] = conn if g == 'E' else conn[elemType != 2] if g == 'Ei' \
                else conn[elemType == 2]
        elif k == 'phi':
            inv[k] = col(n, 4, 8)
        elif k == 'ss':
            inv[k] = col(n, 6)
        elif k in ('stress_i0', 'pml_init6'):
            inv[k] = col(n, 6)
        elif k in ('dNx_i', 'dNy_i', 'dNz_i', 'dNx_p', 'dNy_p', 'dNz_p',
                   'wx_p', 'wy_p', 'wz_p', 'neg_det_w_p'):
            inv[k] = col(n, 8) if not k.startswith('neg') else col(n, 1)
        else:
            inv[k] = col(n)
    inv.update(Ei=Ei, Ep=Ep, E=E,
               int_nodes_idx=np.nonzero(ndof == 3)[0],
               pml_nodes_idx=np.nonzero(ndof == 12)[0])
    inv['idx3_v'] = eq_ids[inv['int_nodes_idx'], 0:3]
    inv['idx12_v'] = eq_ids[inv['pml_nodes_idx'], :]
    inv['a9'] = col(inv['pml_nodes_idx'].shape[0], 9)
    inv['b9'] = col(inv['pml_nodes_idx'].shape[0], 9) + 1.0

    # a fault: every node on the i == nx//2 plane, paired with its neighbour
    m = nid[nx // 2].ravel()
    s = nid[nx // 2 + 1].ravel()
    finv = dict(nftnd=m.shape[0], nsmp1=s, nsmp2=m,
                un=col(m.shape[0], 3), arn=col(m.shape[0]),
                nuc_radius=col(m.shape[0]), friclaw=1, tr=col(m.shape[0]),
                idxF_s=[eq_ids[s, d] for d in range(3)],
                idxF_m=[eq_ids[m, d] for d in range(3)])
    return S, inv, finv


def test_unclassified_invariant_raises():
    """A new array in assembleGlobalKU.build's dict must be classified
    element-axis or replicated deliberately. Defaulting either way is a
    silent wrong answer, so the shard_map split refuses instead."""
    S, inv, finv = synthetic()
    bad = dict(inv, something_new=np.zeros(3))
    with pytest.raises(KeyError, match='classified neither'):
        B._split_sharded(bad, 2, B.pad_counts(bad, 2), B._ELEM_GROUP)


@pytest.mark.parametrize('bad', ['', 'psum', 'ring'])
def test_bad_sync_mode_raises(bad, monkeypatch):
    monkeypatch.setenv(MQ.SYNC_ENV, bad)
    with pytest.raises(ValueError, match='must be one of'):
        MQ.sync_mode()


@pytest.mark.parametrize('bad', ['', 'nodal', 'element+halo'])
def test_bad_shard_mode_raises(bad, monkeypatch):
    monkeypatch.setenv(B.MODE_ENV, bad)
    with pytest.raises(ValueError, match='must be one of'):
        B.shard_mode()


@pytest.mark.parametrize('bad', ['0', '-1', 'two'])
def test_bad_device_count_raises(bad, monkeypatch):
    monkeypatch.setenv(B._DEVICES_ENV, bad)
    with pytest.raises((ValueError, TypeError)):
        B.jax_device_count()


def test_per_rank_cache_dir_is_distinct_per_rank(tmp_path, monkeypatch):
    """The wedge fix: each rank must resolve a DIFFERENT compilation-cache
    directory. Sharing one is what let JAX's per-key lock hang the
    test.tpv8 x python-jax-mpi cell (3 runs lost), and it buys no reuse
    because every rank compiles its own element shapes."""
    monkeypatch.setenv(B._CACHE_ENV, str(tmp_path))
    paths = [B._resolved_cache_path('rank%d' % r) for r in range(4)]
    assert len(set(paths)) == 4
    assert all(p.startswith(str(tmp_path)) for p in paths)
    # and the serial path is unchanged -- no subdir, same directory as before
    assert B._resolved_cache_path() == str(tmp_path)
    assert B._resolved_cache_path(None) == str(tmp_path)


def test_cache_dir_off_is_honoured_for_both_forms(monkeypatch):
    """`off` must stay off even with a subdir -- a per-rank subdirectory of
    "off" would be a real directory named off/rank0 and would silently
    re-enable the cache the user turned off."""
    monkeypatch.setenv(B._CACHE_ENV, 'off')
    assert B._resolved_cache_path() == 'off'
    assert B._resolved_cache_path('rank3') == 'off'


def test_second_different_cache_request_raises(tmp_path, monkeypatch):
    """jax.config is process-global and the first compile wins, so a second
    request for a different directory cannot be honoured. Raise instead of
    returning silently, which would leave the MPI path sharing the directory
    it asked not to share while looking fixed."""
    monkeypatch.setenv(B._CACHE_ENV, str(tmp_path))
    monkeypatch.setattr(B, '_cache_enabled', True)
    monkeypatch.setattr(B, '_cache_path', str(tmp_path / 'rank0'))
    B.enable_compilation_cache('rank0')          # same request: no-op
    with pytest.raises(RuntimeError, match='already pointed at'):
        B.enable_compilation_cache('rank1')


def test_shard_padding_is_exact_multiple_and_pads_divisors_with_one():
    """Padded element rows must contribute EXACTLY 0.0, which needs every
    pad value to be 0 -- except the three PML divisors, where 0 would make
    calcPMLElemKU divide 0/0 and scatter a NaN into the sink."""
    S, inv, finv = synthetic()
    counts = B.pad_counts(inv, 7)
    for g in ('Ei', 'Ep', 'E'):
        assert counts[g] % 7 == 0 and counts[g] >= int(inv[g])
        assert counts[g] - int(inv[g]) < 7
    p = B._prep('b1', np.zeros(inv['Ep']), 7, counts, B._ELEM_GROUP)
    assert p.shape[0] == counts['Ep']
    assert np.all(p[inv['Ep']:] == 1.0)
    q = B._prep('lam_p', np.ones(inv['Ep']), 7, counts, B._ELEM_GROUP)
    assert np.all(q[inv['Ep']:] == 0.0)


def test_allreduce_sync_is_refused_by_name(monkeypatch):
    """The allreduce sync reduced the FULL nodal array and so required it
    replicated at global extent on every rank -- the replication the local
    renumbering removes. It must fail with a sentence, not silently take the
    halo path under an allreduce label. Checked BEFORE any mesh work, so an
    old command line fails in a second rather than after a mesh build."""
    import types
    from eqdyna import driver

    monkeypatch.setenv(MQ.SYNC_ENV, 'allreduce')
    comm = types.SimpleNamespace(Get_rank=lambda: 0, Get_size=lambda: 4)
    fake_jnp = types.SimpleNamespace(__name__='jax.numpy')
    with pytest.raises(ValueError, match='allreduce'):
        driver.run_mpi({'nstep': 1}, comm, MQ.Partition(0, 4, (2, 2, 1)), {},
                       xp=fake_jnp)


def test_mpi_path_still_refuses_numpy():
    """Unchanged contract, re-asserted here because the allreduce check now
    sits next to it: the ordering must keep the numpy refusal FIRST, so a
    numpy run is told it is numpy rather than told about a sync mode."""
    import types
    from eqdyna import driver
    comm = types.SimpleNamespace(Get_rank=lambda: 0, Get_size=lambda: 4)
    with pytest.raises(RuntimeError, match='jax'):
        driver.run_mpi({'nstep': 1}, comm, MQ.Partition(0, 4, (2, 2, 1)), {},
                       xp=np)


# ---------------------------------------------------------------------------
# THE FORTRAN DECOMPOSITION (meshgen.f90 getLocalOneDimCoorArrAndSize /
# calcXyzMPIId), U1/U2 of design 56d2401. Tables below are HAND-DERIVED from
# the Fortran formulas, not produced by the code under test.
# ---------------------------------------------------------------------------
# (global nodes, ranks) -> [(local size, 0-based offset) per rank]
PARTITION_TABLE = {
    (9, 1): [(9, 0)],
    (9, 2): [(5, 0), (5, 4)],            # per 5, resid 0: shared plane 4
    (10, 2): [(5, 0), (6, 4)],           # per 5, resid 1: rank 1 = P-resid gets +1
    (10, 4): [(3, 0), (3, 2), (3, 4), (4, 6)],   # per 3, resid 1: only rank 3 gets +1
    (11, 4): [(3, 0), (3, 2), (4, 4), (4, 7)],   # resid 2: rank 2 is the `<` vs `<=` case
    (89, 2): [(45, 0), (45, 44)],        # test.tpv8 x line
    (89, 4): [(23, 0), (23, 22), (23, 44), (23, 66)],
}


def _hand_partition(n, P, m):
    """The same formula written the long way, as the Fortran reads."""
    per = int((n + P - 1) / P)
    resid = (n + P - 1) - per * P
    size = per if m < (P - resid) else per + 1
    if m <= (P - resid):
        start1 = (per - 1) * m + 1                 # Fortran 1-based first index
    else:
        start1 = (per - 1) * m + 1 + (m - P + resid)
    return size, start1 - 1


def test_partition_table_is_the_fortran_formula():
    for (n, P), want in PARTITION_TABLE.items():
        assert [_hand_partition(n, P, m) for m in range(P)] == want, (n, P)
        assert [meshgen.partition_1d(n, P, m) for m in range(P)] == want, (n, P)


@pytest.mark.parametrize('n', range(4, 60))
@pytest.mark.parametrize('P', [1, 2, 3, 4])
def test_partition_shares_exactly_one_plane_and_covers(n, P):
    parts = [meshgen.partition_1d(n, P, m) for m in range(P)]
    if min(sz for sz, _ in parts) < 2:
        with pytest.raises(ValueError, match='holds'):
            meshgen.check_partition_1d(n, P)
        return
    assert meshgen.check_partition_1d(n, P) == parts
    assert parts[0][1] == 0 and parts[-1][0] + parts[-1][1] == n
    for m in range(P - 1):
        assert parts[m + 1][1] == parts[m][1] + parts[m][0] - 1
    # every ELEMENT (gap between consecutive planes) belongs to exactly one rank
    owners = np.zeros(n - 1, dtype=int)
    for sz, off in parts:
        owners[off:off + sz - 1] += 1
    assert np.all(owners == 1)


def test_calc_xyz_mpi_id_is_z_fastest():
    ids = [meshgen.calc_xyz_mpi_id(r, 4, 2, 2) for r in range(16)]
    assert ids[:5] == [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0)]
    assert ids[15] == (3, 1, 1)
    assert len(set(ids)) == 16


def test_local_line_is_a_bitwise_slice_never_a_reaccumulation():
    """U2: the local line must be the global line's own doubles. A
    re-accumulated geometric stretch differs in the last bits, which is 1e-8 m
    at 1e5 m -- above frt_canonical.align's 1e-9."""
    g = np.cumsum(np.concatenate(([-3.0e4], 500.0 * 1.025 ** np.arange(88))))
    for P in (1, 2, 4):
        for m in range(P):
            loc, off = meshgen.local_line(g, P, m)
            n, o = meshgen.partition_1d(g.size, P, m)
            assert off == o and loc.size == n
            assert loc.tobytes() == g[o:o + n].tobytes()


def test_decomp_table_is_fortrans_and_refuses_other_counts():
    assert MQ.DECOMP == {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2),
                         16: (4, 2, 2), 32: (4, 4, 2)}
    for n, dims in MQ.DECOMP.items():
        assert dims[0] * dims[1] * dims[2] == n
    with pytest.raises(ValueError, match='DECOMP'):
        MQ.Partition.for_size(0, 3)
    with pytest.raises(ValueError, match='multiply'):
        MQ.Partition(0, 4, (2, 1, 1))
    with pytest.raises(ValueError, match='out of range'):
        MQ.Partition(4, 4, (2, 2, 1))


def test_neighbours_and_model_edges():
    """me -/+ npy*npz, npz, 1 (assembleGlobalMass.f90's dest/source)."""
    p = MQ.Partition(5, 16, (4, 2, 2))            # (mex,mey,mez) = (1,0,1)
    assert p.mexyz == (1, 0, 1)
    assert [p.neighbour(0, 0), p.neighbour(0, 1)] == [1, 9]
    assert [p.neighbour(1, 0), p.neighbour(1, 1)] == [3, 7]
    assert [p.neighbour(2, 0), p.neighbour(2, 1)] == [4, 6]
    assert [p.at_model_edge(d, ib) for d in range(3) for ib in (0, 1)] == \
        [False, False, True, False, False, True]


def test_face_node_ids_follow_the_fortran_loop_order():
    nx, ny, nz = 3, 4, 2
    nid = lambda ix, iy, iz: (ix - 1) * ny * nz + (iz - 1) * ny + iy
    assert list(MQ.face_node_ids(0, 1, nx, ny, nz)) == \
        [nid(3, iy, iz) for iz in (1, 2) for iy in (1, 2, 3, 4)]
    assert list(MQ.face_node_ids(1, 0, nx, ny, nz)) == \
        [nid(ix, 1, iz) for ix in (1, 2, 3) for iz in (1, 2)]
    assert list(MQ.face_node_ids(2, 1, nx, ny, nz)) == \
        [nid(ix, iy, 2) for ix in (1, 2, 3) for iy in (1, 2, 3, 4)]


def test_every_shared_node_has_exactly_one_owner():
    """The lowest holder writes: over a (2,2,2) split of a 5x4x3 grid every
    GLOBAL node is owned by exactly one rank."""
    dims, n = (2, 2, 2), (5, 4, 3)
    count = np.zeros(n, dtype=int)
    for r in range(8):
        p = MQ.Partition(r, 8, dims)
        sz = [meshgen.partition_1d(n[d], dims[d], p.mexyz[d]) for d in range(3)]
        nx, ny, nz = (s[0] for s in sz)
        ids = np.arange(1, nx * ny * nz + 1)
        own = MQ.owned_mask(p, ids, (nx, ny, nz))
        s = ids[own] - 1
        gx = sz[0][1] + s // (nz * ny)
        gz = sz[2][1] + (s % (nz * ny)) // ny
        gy = sz[1][1] + s % ny
        np.add.at(count, (gx, gy, gz), 1)
    assert np.all(count == 1)


def test_fault_boundary_lists_decode_fltgm():
    """createMasterNode's fltgm codes and MPI4arn's six lists."""
    nx, ny, nz = 3, 3, 3
    nid = lambda ix, iy, iz: (ix - 1) * ny * nz + (iz - 1) * ny + iy
    slaves = [nid(1, 2, 1), nid(2, 2, 2), nid(3, 2, 3), nid(1, 2, 3)]
    nsmp = np.array([[s, 100 + i] for i, s in enumerate(slaves)])
    lists = meshgen.fault_boundary_lists(nsmp, nx, ny, nz)
    assert [list(x) for x in lists] == [[1, 4], [3], [], [], [1], [3, 4]]

