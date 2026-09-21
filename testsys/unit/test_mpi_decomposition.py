"""Tests for the two PY-ONLY parallel paths added for pathway item 43.

Neither path exists in the Fortran, so the 30-cell sweep cannot cover them:
the sweep runs the serial python column, and every one of these code paths is
reached only when EQDYNA_JAX_DEVICES > 1 (shard_map) or under mpirun
(driver.run_mpi). They are therefore tested here, in the same change that
adds them.

WHAT IS AND IS NOT COVERED HERE. These are the PARTITION's invariants --
completeness, disjointness, symmetry, and the loud failures -- on a synthetic
mesh, because they are pure index algebra and do not need a solver. The
numerical question (does an N-rank run land inside the case bound) is not a
unit test: it is a full-length gated run, and it is
testsys/perf/run_mpi_scaling.py's parity companion, documented in
NOTES_item43_mpi.md. A unit test that "passed" on internal consistency while
the physics moved is exactly the failure mode this repo's rule 1 exists for.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, os.pardir,
                                'src', 'python'))
from eqdyna import MPI4NodalQuant as MQ   # noqa: E402
from eqdyna import backend as B           # noqa: E402


def synthetic(nx=6, ny=3, nz=3, npml=1):
    """A structured brick mesh with a PML shell and a two-sided fault, as the
    index arrays the partition consumes. NOT a physics case: every float is a
    placeholder, because decompose() only ever SELECTS rows and never reads a
    value. The shapes and the connectivity are what is under test."""
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
            inv[k] = ([np.arange(n * 8) % (NEQ + 1) for _ in range(12)]
                      if k == 'idxP12' else
                      [np.arange(n * 8) % (NEQ + 1) for _ in range(3)]
                      if k == 'idxP3' else np.arange(n * 8) % (NEQ + 1))
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


@pytest.mark.parametrize('nranks', [1, 2, 3, 4])
def test_partition_is_complete_and_disjoint(nranks):
    """Every element assigned exactly once, every fault node owned exactly
    once. A partition that drops elements still runs and still produces
    plausible output -- it just solves a smaller problem."""
    S, inv, finv = synthetic()
    got_i = got_p = got_e = 0
    owned = []
    for r in range(nranks):
        d = MQ.decompose(S, inv, finv, r, nranks)
        got_i += d['inv']['Ei']; got_p += d['inv']['Ep']; got_e += d['inv']['E']
        owned.append(d['fault_rows'])
        # the restricted arrays must agree with the declared counts
        assert d['inv']['lam_i'].shape[0] == d['inv']['Ei']
        assert d['inv']['idxIx'].shape[0] == d['inv']['Ei'] * 8
        assert d['inv']['conn'].shape[0] == d['inv']['E']
        assert d['inv']['idxP12'][0].shape[0] == d['inv']['Ep'] * 8
    assert (got_i, got_p, got_e) == (inv['Ei'], inv['Ep'], inv['E'])
    allowned = np.concatenate(owned)
    assert np.array_equal(np.sort(allowned), np.arange(finv['nftnd']))


@pytest.mark.parametrize('nranks', [2, 3, 4])
def test_halo_is_symmetric_and_same_order(nranks):
    """Rank r's shared-equation list with s must be s's with r, element for
    element. Sendrecv pairs the buffers positionally: a different order on the
    two sides adds neighbours' partials to the WRONG equations, which is a
    wrong answer with no error anywhere."""
    S, inv, finv = synthetic()
    dec = [MQ.decompose(S, inv, finv, r, nranks) for r in range(nranks)]
    for r in range(nranks):
        halo_r = dec[r]['halo_idx']
        for s, pos_rs in dec[r]['neighbours']:
            pos_sr = dict(dec[s]['neighbours'])[r]
            assert np.array_equal(halo_r[pos_rs], dec[s]['halo_idx'][pos_sr])


@pytest.mark.parametrize('nranks', [2, 3, 4])
def test_every_touched_equation_is_local_or_exchanged(nranks):
    """The postcondition the whole design rests on: after the exchange, this
    rank's force must be complete at every equation it will READ. So an
    equation another rank also touches must be in the halo."""
    S, inv, finv = synthetic()
    conn = S['conn']
    dec = [MQ.decompose(S, inv, finv, r, nranks) for r in range(nranks)]
    touch = []
    for r in range(nranks):
        lo, hi = dec[r]['report']['elem_lo'], dec[r]['report']['elem_hi']
        t = np.zeros(S['N'], dtype=bool); t[conn[lo:hi].ravel()] = True
        e = S['eq_ids'][t].ravel()
        touch.append(np.unique(e[e > 0]))
    for r in range(nranks):
        others = np.unique(np.concatenate([touch[s] for s in range(nranks) if s != r]))
        shared = np.intersect1d(touch[r], others)
        assert np.array_equal(np.sort(dec[r]['halo_idx'].astype(np.int64)), shared)


def test_one_rank_is_the_identity_restriction():
    """1 rank must reduce to the serial arrays exactly -- same rows, same
    order, empty halo. This is what makes `mpirun -np 1` a usable
    equivalence check against the serial column."""
    S, inv, finv = synthetic()
    d = MQ.decompose(S, inv, finv, 0, 1)
    assert d['halo_idx'].size == 0 and d['neighbours'] == []
    for k in B._ELEM_GROUP:
        a, b = inv[k], d['inv'][k]
        if isinstance(a, list):
            for x, y in zip(a, b):
                assert np.array_equal(np.asarray(x), np.asarray(y))
        else:
            assert np.array_equal(np.asarray(a), np.asarray(b))
    assert d['finv']['nftnd'] == finv['nftnd']


def test_more_ranks_than_elements_raises():
    """Loudly, rather than handing some rank an empty subdomain that
    contributes nothing while still counting in the speedup."""
    S, inv, finv = synthetic(nx=2, ny=1, nz=1, npml=0)
    with pytest.raises(ValueError, match='would get 0 of'):
        MQ.decompose(S, inv, finv, 0, inv['E'] + 1)


def test_rank_out_of_range_raises():
    S, inv, finv = synthetic()
    with pytest.raises(ValueError, match='out of range'):
        MQ.decompose(S, inv, finv, 4, 4)


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
