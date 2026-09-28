"""Tests for the PY-ONLY parallel path: the rank bookkeeping of the
python-jax-mpi decomposition (driver.run_mpi under mpirun). (The shard_map
explicit-decomposition path this file used to also cover was retired
2026-09-28 -- owner decision: jax-MPI is the parallel path now.)

This path is not reached by the serial sweep columns, so it is tested here.

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


@pytest.mark.parametrize('bad', ['', 'psum', 'ring'])
def test_bad_sync_mode_raises(bad, monkeypatch):
    monkeypatch.setenv(MQ.SYNC_ENV, bad)
    with pytest.raises(ValueError, match='must be one of'):
        MQ.sync_mode()


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

