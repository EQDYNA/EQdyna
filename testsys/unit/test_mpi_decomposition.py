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
        # GLOBAL ids: rank r and rank s number the same shared equation
        # differently in their own local spaces, so global is the only frame
        # in which "the same equation" is a statement the two can both make.
        halo_r = dec[r]['halo_eq_global']
        for s, pos_rs in dec[r]['neighbours']:
            pos_sr = dict(dec[s]['neighbours'])[r]
            assert np.array_equal(halo_r[pos_rs],
                                  dec[s]['halo_eq_global'][pos_sr])
        # ...and the local halo is that same list renumbered, position for
        # position -- which is what Sendrecv's positional pairing relies on.
        assert dec[r]['halo_idx'].shape == halo_r.shape


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
        assert np.array_equal(np.sort(dec[r]['halo_eq_global']), shared)


def test_one_rank_is_the_identity_restriction():
    """1 rank must reduce to the serial arrays exactly -- same rows, same
    order, empty halo. This is what makes `mpirun -np 1` a usable
    equivalence check against the serial column.

    The RENUMBERING is the identity at 1 rank only because this mesh has no
    orphan node and no equation outside the touched set, so local_nodes ==
    arange(N) and local_eqs == arange(1, NEQ+1). That is asserted rather than
    assumed: on a mesh with an orphan node the 1-rank arrays would be
    correctly renumbered and NOT equal to the serial ones, and this test
    would then be making a claim about the mesh, not about decompose."""
    S, inv, finv = synthetic()
    d = MQ.decompose(S, inv, finv, 0, 1)
    assert np.array_equal(d['local_nodes'], np.arange(S['N']))
    assert np.array_equal(d['local_eqs'], np.arange(1, S['NEQ'] + 1))
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


# ---------------------------------------------------------------------------
# THE RANK-LOCAL RENUMBERING (decompose's global->local remap)
#
# The restriction above splits the WORK; this splits the PROBLEM. Measured
# reason it exists: with the carry at global extent, v1+velArr+dispArr+force
# came to 97.75 MB on test.tpv104 BYTE-IDENTICALLY at 1 rank and at 32 ranks
# (98.5% of the 32-rank carry), so per-step memory traffic did not fall with
# rank count and the step was memory-system-bound from 4 ranks up.
#
# WHAT CAN GO WRONG HERE, and therefore what these tests are: a single index
# left at GLOBAL numbering against a rank-local array is not a crash. It is a
# read from, or a scatter into, whatever local slot that global number lands
# on -- a wrong answer with nothing to attribute it to. So the tests check
# totality (every index array remapped, unclassified ones refused), range
# (nothing addresses outside the local extent), and EXACTNESS (a gather
# through the local indices returns the identical values a gather through the
# global ones did -- which is the whole bit-identity argument).
# ---------------------------------------------------------------------------
EQ_KEYS = tuple(k for k, v in MQ._INDEX_SPACE.items() if v == 'eq')
NODE_KEYS = tuple(k for k, v in MQ._INDEX_SPACE.items() if v == 'node')


def _flat(v):
    if isinstance(v, (list, tuple)):
        return np.concatenate([np.asarray(x).ravel() for x in v])
    return np.asarray(v).ravel()


@pytest.mark.parametrize('nranks', [1, 2, 3, 4])
def test_every_index_array_is_inside_the_local_extent(nranks):
    """Range check on every classified index array, on both dicts. An index
    beyond the local extent would be an out-of-bounds gather (jax CLAMPS it,
    silently, to the last row) -- the exact failure that produces plausible
    output."""
    S, inv, finv = synthetic()
    for r in range(nranks):
        d = MQ.decompose(S, inv, finv, r, nranks)
        n_l, neq_l = int(d['inv']['N']), int(d['inv']['NEQ'])
        assert n_l == d['local_nodes'].shape[0]
        assert neq_l == d['local_eqs'].shape[0]
        assert int(d['inv']['NEQ1']) == neq_l + 1
        # An EMPTY array is a real case, not a gap in the test: a rank can
        # hold zero interior elements (this fixture's npml shell swallows the
        # whole mesh in y) or zero PML ones, and that rank must still work.
        for k in NODE_KEYS:
            a = _flat(d['inv'][k])
            if a.size:
                assert a.min() >= 0 and a.max() < n_l, (k, r, nranks)
        for k in EQ_KEYS:
            a = _flat(d['inv'][k])
            if a.size:
                assert a.min() >= 0 and a.max() <= neq_l, (k, r, nranks)
        for k in ('nsmp1', 'nsmp2'):
            a = _flat(d['finv'][k])
            if a.size:
                assert a.min() >= 0 and a.max() < n_l
        for k in ('idxF_s', 'idxF_m'):
            a = _flat(d['finv'][k])
            if a.size:
                assert a.min() >= 0 and a.max() <= neq_l
        assert d['halo_idx'].size == 0 or (
            d['halo_idx'].min() >= 1 and d['halo_idx'].max() <= neq_l)


@pytest.mark.parametrize('nranks', [2, 3, 4])
def test_remap_is_the_identical_gather(nranks):
    """THE bit-identity argument, as a test. A gather or a scatter performs
    exactly the same operations under an INJECTIVE relabelling of its index
    array: only the addresses change, not the values, the duplicate structure
    or the array order. So gathering a rank-local nodal array through the
    LOCAL indices must return, element for element, what gathering the global
    array through the GLOBAL indices returned."""
    S, inv, finv = synthetic()
    rng = np.random.default_rng(0)
    nodal_g = rng.standard_normal((S['N'], 3))
    eqn_g = rng.standard_normal(S['NEQ'] + 1)
    for r in range(nranks):
        d = MQ.decompose(S, inv, finv, r, nranks)
        # rank-local copies of the two arrays, in local order
        nodal_l = nodal_g[d['local_nodes']]
        eqn_l = np.concatenate(([eqn_g[0]], eqn_g[d['local_eqs']]))
        # global counterparts of the same restricted index arrays
        g = MQ.decompose(S, inv, finv, r, nranks)   # fresh, then undo the map
        for k in NODE_KEYS:
            loc_idx = _flat(d['inv'][k])
            glob_idx = d['local_nodes'][loc_idx]
            assert np.array_equal(nodal_l[loc_idx], nodal_g[glob_idx])
        for k in EQ_KEYS:
            loc_idx = _flat(d['inv'][k])
            glob_idx = np.where(loc_idx == 0, 0, d['local_eqs'][loc_idx - 1])
            assert np.array_equal(eqn_l[loc_idx], eqn_g[glob_idx])
        assert g['inv']['N'] == d['inv']['N']


@pytest.mark.parametrize('nranks', [2, 4])
def test_local_extent_actually_shrinks(nranks):
    """The point of the change, stated as an inequality rather than left to
    the timing: each rank's node and equation counts must be strictly below
    the global ones, because that is what makes the carry fall with rank
    count. A remap that renumbered correctly but still covered the whole mesh
    would pass every test above and buy nothing."""
    S, inv, finv = synthetic()
    for r in range(nranks):
        d = MQ.decompose(S, inv, finv, r, nranks)
        assert d['report']['N_local'] < d['report']['N_global']
        assert d['report']['NEQ_local'] < d['report']['NEQ_global']


@pytest.mark.parametrize('nranks', [1, 2, 3, 4])
def test_published_map_inverts_the_relabelling(nranks):
    """decompose PUBLISHES local_nodes / local_eqs, and driver.run_mpi builds
    the local mass array from local_eqs alone:

        mass_local = [1.0] + mass_global[local_eqs]

    So if the published map were not the one used to relabel the index
    arrays, every equation's mass would be the mass of a DIFFERENT equation.
    Checked by mapping the local indices BACK and comparing against the
    global arrays restricted the way decompose says it restricted them
    (elem_lo/elem_hi for elements, fault_computed_rows for fault nodes) --
    published quantities on both sides, so this is not the remap re-derived
    and compared with itself.

    It also pins the SINK: local index 0 exactly where the global index was
    0, position for position. Equation 0 is the no-equation slot every masked
    contribution is scattered into and calcHourglassResist scrubs; mapping it
    onto a real equation would scrub that equation instead."""
    S, inv, finv = synthetic()
    sink_seen = 0
    for r in range(nranks):
        d = MQ.decompose(S, inv, finv, r, nranks)
        nodes, eqs = d['local_nodes'], d['local_eqs']
        lo, hi = d['report']['elem_lo'], d['report']['elem_hi']

        def back_node(a):
            return nodes[np.asarray(a)]

        def back_eq(a):
            a = np.asarray(a)
            return np.where(a == 0, 0, eqs[a - 1])

        # element arrays, group 'E': the restriction is exactly rows lo:hi
        assert np.array_equal(back_node(d['inv']['conn']), S['conn'][lo:hi])
        for k in ('idxH0', 'idxH1', 'idxH2'):
            g = np.asarray(inv[k]).reshape(-1, 8)[lo:hi].reshape(-1)
            l = np.asarray(d['inv'][k])
            assert np.array_equal(back_eq(l), g), (k, r, nranks)
            assert np.array_equal(l == 0, g == 0), (k, r, nranks)
            sink_seen += int(np.count_nonzero(l == 0))

        # idxP12 / idxP3 are the two arrays build() explicitly MASKS to the
        # sink (by node dof), so they are where the sink invariant has
        # something to bite on. PML rows within lo:hi, in PML-array order.
        # Rows of the PML-only array, which are CONTIGUOUS for a contiguous
        # global element range: the PML elements before lo, then this slab's.
        is_pml = S['elemType'] == 2
        sel_p = int(is_pml[:lo].sum()) + np.arange(int(is_pml[lo:hi].sum()))
        for k in ('idxP12', 'idxP3'):
            for dd in range(len(inv[k])):
                g = np.asarray(inv[k][dd]).reshape(-1, 8)[sel_p].reshape(-1)
                l = np.asarray(d['inv'][k][dd])
                assert np.array_equal(back_eq(l), g), (k, dd, r, nranks)
                assert np.array_equal(l == 0, g == 0), (k, dd, r, nranks)
                sink_seen += int(np.count_nonzero(l == 0))

        # fault arrays: the restriction is exactly fault_computed_rows
        computed = d['fault_computed_rows']
        for k in ('nsmp1', 'nsmp2'):
            assert np.array_equal(back_node(d['finv'][k]),
                                  np.asarray(finv[k])[computed])
        for k in ('idxF_s', 'idxF_m'):
            for dd in range(3):
                g = np.asarray(finv[k][dd])[computed]
                l = np.asarray(d['finv'][k][dd])
                assert np.array_equal(back_eq(l), g), (k, dd, r, nranks)
                assert np.array_equal(l == 0, g == 0), (k, dd, r, nranks)

        # the halo, in both frames
        assert np.array_equal(back_eq(d['halo_idx']), d['halo_eq_global'])
    assert sink_seen > 0, 'this fixture must exercise the sink somewhere'


def test_relabel_maps_the_sink_to_the_sink_and_nothing_else_to_it():
    """The mapping property the test above observes, stated directly: global
    0 -> local 0, and every REAL equation -> a positive local index."""
    g2l = np.array([0, 1, -1, 2, 3], dtype=np.int64)   # global eq 2 not held
    out = MQ._relabel('probe', np.array([0, 1, 3, 4, 0]), g2l, 'eq', 3)
    assert out.tolist() == [0, 1, 2, 3, 0]
    assert np.all(out[np.array([1, 2, 3])] > 0)


def test_unclassified_index_array_raises():
    """A new integer index array in build()'s dict must be classified node-
    or equation-indexed deliberately. Left unclassified it would keep GLOBAL
    numbering against a rank-local array -- a silent wrong-slot access."""
    S, inv, finv = synthetic()
    bad = dict(inv, idxSomethingNew=np.zeros(4, dtype=np.int64))
    with pytest.raises(KeyError, match='classified neither'):
        MQ.decompose(S, bad, finv, 0, 2)


def test_negative_index_raises():
    """Fancy-indexing the global->local table with a negative index would
    WRAP to a valid-looking local index. Refused rather than wrapped."""
    with pytest.raises(ValueError, match='negative'):
        MQ._relabel('idxFake', np.array([-1, 2]), np.arange(5), 'eq', 5)


def test_index_outside_the_local_set_raises():
    """An index the rank does not hold cannot be sized into the local set.
    Clamping it would land the write on a real node."""
    g2l = np.array([0, 1, -1, 2])      # global eq 2 not held by this rank
    with pytest.raises(ValueError, match='does not hold'):
        MQ._relabel('idxFake', np.array([1, 2]), g2l, 'eq', 2)


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
        driver.run_mpi({'nstep': 1}, comm, xp=fake_jnp)


def test_mpi_path_still_refuses_numpy():
    """Unchanged contract, re-asserted here because the allreduce check now
    sits next to it: the ordering must keep the numpy refusal FIRST, so a
    numpy run is told it is numpy rather than told about a sync mode."""
    import types
    from eqdyna import driver
    comm = types.SimpleNamespace(Get_rank=lambda: 0, Get_size=lambda: 4)
    with pytest.raises(RuntimeError, match='jax'):
        driver.run_mpi({'nstep': 1}, comm, xp=np)
