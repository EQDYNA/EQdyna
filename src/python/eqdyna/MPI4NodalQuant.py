"""MPI4NodalQuant.py <- assembleGlobalMass.f90:58-245 (the MPI4NodalQuant
subroutine) + the per-rank domain decomposition meshgen.f90 performs from
npx/npy/npz.

WHY THIS FILE HAS A SUBROUTINE'S NAME AND NOT A FILE'S. Every other module
here is named after the Fortran FILE it ports. MPI4NodalQuant is a subroutine
inside assembleGlobalMass.f90, and folding 300 lines of domain decomposition
into assembleGlobalMass.py (670 lines of lumped-mass assembly) would bury the
one thing this code is: the nodal exchange, and the partition it exchanges
across. It is named after the Fortran entity it ports, which is the rule the
convention is actually serving.

WHAT THIS IS FOR. The jax backend's per-step kernel is FASTER than Fortran's
at one core (611 vs 931 ms/step on test.tpv104) and does not scale: 2.76x at
16 cores against Fortran MPI's 14.39x. Two routes were measured before this
one was written, and both are recorded so this file is not mistaken for a
first guess:

  XLA AUTOMATIC partitioning -- 0 of 38 HLO scatters ever get a partition
    annotation, and fusions partition only to round(sqrt(N)) (item 43).
  EXPLICIT jax sharding (shard_map over host CPU devices) -- implemented and
    measured; see backend.run_time_loop_sharded and
    testsys/perf/run_shard_scaling.py. It works and it is bounded by two
    things a device mesh cannot fix: the nodal stages stay REPLICATED on
    every device, and the collective is an all-reduce of the WHOLE nodal
    array (O(NEQ)) where MPI moves O(boundary).

This file is the third route and it is not a jax idea at all: it is the
Fortran decomposition, one OS PROCESS per rank, with jax owning only the
local element kernel -- which is exactly where it already wins. Nothing in
XLA has to cooperate: there is no partition annotation to hope for, no
host-platform device fiction, the duplicate-index scatter-add becomes a
purely LOCAL operation, and the halo is O(boundary) by construction.

THE DECOMPOSITION, and how it differs from Fortran's. Fortran gives each rank
its own MESH: meshgen.f90 generates only that rank's nodes and elements, with
rank-local numbering. This port keeps the PROVEN serial mesh build on every
rank and then RESTRICTS it -- same global node and equation numbering on
every rank, each rank holding the full-length nodal arrays but touching only
its own elements' entries. That trade is deliberate:

  + every index array (eq_ids, conn, idx*, nsmp) comes from the serial code
    path that is already gated on 10 cases x 3 backends, so the partition
    cannot introduce a numbering bug -- it only SELECTS rows.
  + a rank's frt output is a subset of the serial rows, so the existing
    per-rank `frt.txt<rank>` canonicalisation (testsys/frt_canonical.py,
    which already globs frt.txt* for exactly this reason) compares an
    N-rank python run against the SAME committed reference as a 4-rank
    Fortran run and a serial run. No new comparison path.
  - every rank pays the serial mesh build (fixed cost, differenced out of
    any per-step number) and holds the full-length nodal arrays
    (NEQ+1 doubles = 30.5 MB on test.tpv104, times nranks).

THE EXCHANGE. After each rank assembles ITS elements, an equation on a
subdomain boundary holds a PARTIAL sum on each rank that touches it.
MPI4NodalQuant below sums those partials, exactly as
assembleGlobalMass.f90:58-245 does for Fortran, and it moves only the shared
equations -- not the array. Everything after the exchange (velDispUpdate,
faulting, the mass divide) is then plain local work on a force array that is
COMPLETE at every equation this rank will read, which is the same
postcondition Fortran relies on and the reason no velocity exchange is
needed: a node shared by two ranks is updated redundantly by both, from
identical complete forces, so the two agree without talking.

BIT-EXACTNESS. An N-rank run is NOT bit-identical to a serial one and cannot
be: a boundary equation's contributions are summed in a different order
(own-partial + received-partial, instead of element by element), and float
addition is not associative. It is gated the way every other backend column
is gated -- against the committed canonical reference at the case's bound in
testsys/matrix.py. Deterministic and reproducible run to run: the neighbour
loop walks ranks in ascending order, so the additions happen in a fixed
order for a given rank count.
"""
import numpy as np

from . import backend as B

# Relative cost of a PML element against an interior one, used ONLY to choose
# the cut points so the ranks get equal WORK rather than equal element counts.
# calcPMLElemKU computes 15 split stress components and 12 force blocks
# against calcElemKU's 6 and 3, so a PML element is several times an interior
# one; the exact figure is a load-balance heuristic, not physics, which is why
# decompose() REPORTS each rank's (Ei, Ep) so an imbalance is visible in the
# measurement instead of hiding inside it.
PML_WEIGHT = 3.0


def _cuts(weights, nranks):
    """Contiguous element ranges whose weight sums are as equal as the
    integer cuts allow. Contiguity is the point: serial element order is the
    mesh generator's nested loop order, so a contiguous range is a spatial
    slab and its shared-equation set is a surface. decompose() measures that
    surface and reports it, rather than trusting the claim."""
    c = np.concatenate(([0.0], np.cumsum(weights)))
    total = c[-1]
    edges = [0]
    for r in range(1, nranks):
        edges.append(int(np.searchsorted(c, total * r / nranks)))
    edges.append(len(weights))
    # A rank with zero elements would make its local kernel shapes zero and
    # its exchange meaningless -- refuse loudly rather than produce a rank
    # that contributes nothing while still being counted in the speedup.
    for r in range(nranks):
        if edges[r + 1] <= edges[r]:
            raise ValueError(
                'MPI4NodalQuant._cuts: rank %d would get 0 of %d elements at '
                'nranks=%d. Use fewer ranks than elements.'
                % (r, len(weights), nranks))
    return edges


def _take(a, sel):
    return np.asarray(a)[sel]


def _take_ravelled(a, sel):
    """The (n*8,) index arrays assembleGlobalKU.build already ravelled:
    reshaped, row-selected, re-ravelled, so the kernel receives exactly the
    flat layout it receives serially and the within-rank scatter order is the
    serial one."""
    return np.asarray(a).reshape(-1, 8)[sel].reshape(-1)


# WHICH INDEX SPACE each integer index array addresses -- the table that makes
# the global->local remap below total rather than a list somebody remembered to
# extend. Every integer ARRAY in assembleGlobalKU.build's dict must appear in
# _INDEX_SPACE and every one in faulting.build's in _FAULT_INDEX_SPACE, or
# _relabel_all raises: an unremapped index array is not a crash, it is a read
# from (or a scatter into) whatever local slot that global number happens to
# land on -- a wrong answer with nothing to attribute it to. Same discipline,
# and the same reason, as backend._ELEM_GROUP.
_INDEX_SPACE = {
    # node-indexed: they address velArr / dispArr, which are (N_local, 3)
    'conn': 'node', 'conn_i': 'node', 'conn_p': 'node',
    'int_nodes_idx': 'node', 'pml_nodes_idx': 'node',
    # equation-indexed: they address v1 / force / mass, which are (NEQ_local+1,)
    'idx3_v': 'eq', 'idx12_v': 'eq',
    'idxIx': 'eq', 'idxIy': 'eq', 'idxIz': 'eq',
    'idxP12': 'eq', 'idxP3': 'eq',
    'idxH0': 'eq', 'idxH1': 'eq', 'idxH2': 'eq',
}
_FAULT_INDEX_SPACE = {
    'nsmp1': 'node', 'nsmp2': 'node',
    'idxF_s': 'eq', 'idxF_m': 'eq',
}


def _is_index_array(v):
    """An integer ARRAY (or list of them), i.e. something that addresses a
    nodal or equation array. Excludes the integer SCALARS in the same dicts
    (N, NEQ, Ei, friclaw, TPV, ...), which are counts and trace-time
    branches, not indices."""
    if isinstance(v, (list, tuple)):
        return bool(len(v)) and all(_is_index_array(x) for x in v)
    return (isinstance(v, np.ndarray) and v.ndim >= 1
            and v.dtype.kind in 'iu')


def _relabel(name, v, g2l, space, n_local):
    """`v` with every global index replaced by its local one.

    A gather or a scatter reproduces exactly the same additions, in the same
    order, under an INJECTIVE relabelling of its index array: duplicate
    structure and array order are untouched, only the addresses change. That
    is why this remap is expected to be BIT-IDENTICAL to the global-extent
    run at the same rank count, and it is gated as such rather than at the
    case bound (see NOTES_item43_mpi.md).
    """
    if isinstance(v, (list, tuple)):
        return [_relabel(name, x, g2l, space, n_local) for x in v]
    a = np.asarray(v)
    if a.size and int(a.min()) < 0:
        raise ValueError(
            'MPI4NodalQuant._relabel: %r holds a negative %s index (min %d). '
            'Fancy-indexing the global->local table with it would wrap to a '
            'VALID-looking local index and scatter into the wrong slot.'
            % (name, space, int(a.min())))
    out = g2l[a]
    nbad = int(np.count_nonzero(out < 0))
    if nbad:
        raise ValueError(
            'MPI4NodalQuant._relabel: %d of %d entries of %r address a %s '
            'this rank does not hold, so it cannot be sized into the local '
            '%s set (%d entries). The element restriction and the touched-set '
            'are out of step; do not clamp -- the clamped write would land on '
            'a real node.' % (nbad, out.size, name, space, space, n_local))
    return out


def _relabel_all(d, table, g2l_node, g2l_eq, n_node, n_eq):
    """`d` with every classified index array remapped, in place."""
    for k in list(d):
        v = d[k]
        if not _is_index_array(v):
            continue
        space = table.get(k)
        if space is None:
            raise KeyError(
                'MPI4NodalQuant._relabel_all: %r is an integer index array '
                'that is classified neither node-indexed nor equation-indexed. '
                'Add it to _INDEX_SPACE / _FAULT_INDEX_SPACE. Leaving it at '
                'GLOBAL numbering against a rank-local array is a silent '
                'wrong-slot access, not an error.' % (k,))
        d[k] = _relabel(k, v, g2l_node if space == 'node' else g2l_eq,
                        space, n_node if space == 'node' else n_eq)
    return d


def decompose(S, inv, finv, rank, nranks):
    """Rank-local `inv`, `finv`, exchange plan and output row selection.

    Returns a dict with
        inv, finv      -- the same dicts, element rows and fault-node rows
                          restricted to this rank, and EVERY index array
                          renumbered into this rank's local node / equation
                          space (inv['N'], inv['NEQ'], inv['NEQ1'] are the
                          local extents, which is what driver.run_mpi sizes
                          the carry on)
        local_nodes    -- global node id of each local node, ascending
        local_eqs      -- global equation id of each local equation,
                          ascending; local equation j+1 is global
                          local_eqs[j], and local 0 is the sink, as globally
        halo_idx       -- LOCAL equation indices this rank exchanges (int32)
        halo_eq_global -- the same equations as GLOBAL ids; the only form in
                          which two ranks can agree on a shared equation, so
                          it is what the symmetry invariant is checked on
        neighbours     -- [(rank, positions-into-halo_idx), ...], ascending
        fault_rows     -- fault-node rows this rank OWNS and therefore writes
        report         -- counts for the measurement to print (rule: a
                          multi-rank number that does not state its per-rank
                          element count cannot be checked)
    """
    if not 0 <= rank < nranks:
        raise ValueError('decompose: rank %d out of range for nranks %d'
                         % (rank, nranks))
    elemType = S['elemType']
    eq_ids = S['eq_ids']
    conn = S['conn']
    E = conn.shape[0]

    # --- element partition: contiguous, work-weighted -----------------------
    w = np.where(elemType == 2, PML_WEIGHT, 1.0)
    edges = _cuts(w, nranks)
    lo, hi = edges[rank], edges[rank + 1]
    mine = np.zeros(E, dtype=bool)
    mine[lo:hi] = True

    # Positions within the interior / PML / all-element arrays. build() keeps
    # those arrays in ascending global-element order, so a contiguous global
    # range maps to a contiguous range in each -- computed, not assumed.
    is_int = (elemType == 1) | (elemType > 10)
    is_pml = elemType == 2
    sel_i = np.nonzero(mine[np.nonzero(is_int)[0]])[0]
    sel_p = np.nonzero(mine[np.nonzero(is_pml)[0]])[0]
    sel_e = np.arange(lo, hi)

    # The element-axis classification lives in ONE table (backend._ELEM_GROUP,
    # written for the shard_map route) so that a new array in build()'s dict
    # cannot be handled by one decomposition and forgotten by the other.
    inv_l = dict(inv)
    for k, group in B._ELEM_GROUP.items():
        sel = {'Ei': sel_i, 'Ep': sel_p, 'E': sel_e}[group]
        v = inv[k]
        if isinstance(v, (list, tuple)):
            inv_l[k] = [_take_ravelled(x, sel) if k in B._RAVELLED else _take(x, sel)
                        for x in v]
        else:
            inv_l[k] = (_take_ravelled(v, sel) if k in B._RAVELLED
                        else _take(v, sel))
    inv_l['Ei'] = int(sel_i.shape[0])
    inv_l['Ep'] = int(sel_p.shape[0])
    inv_l['E'] = int(sel_e.shape[0])

    # --- nodal restriction: the nodes this rank's elements touch ------------
    # velDispUpdate then integrates exactly those nodes. A node shared with
    # another rank is integrated by BOTH, from the same post-exchange force,
    # so the two agree with no further communication -- this is the reason
    # MPI4NodalQuant exchanges force and nothing exchanges velocity.
    touched = np.zeros(S['N'], dtype=bool)
    touched[conn[lo:hi].ravel()] = True
    mask_int = touched[inv['int_nodes_idx']]
    mask_pml = touched[inv['pml_nodes_idx']]
    for k, m in (('int_nodes_idx', mask_int), ('idx3_v', mask_int),
                 ('pml_nodes_idx', mask_pml), ('idx12_v', mask_pml),
                 ('a9', mask_pml), ('b9', mask_pml)):
        inv_l[k] = np.asarray(inv[k])[m]

    # --- exchange plan ------------------------------------------------------
    # Every rank can compute every other rank's touched set (it holds the full
    # mesh), so the plan needs no communication to build and both sides of a
    # pair derive the SAME shared-equation list in the same (sorted) order.
    my_eq = _eqs_of(eq_ids, touched)
    neighbours, halo_list = [], []
    for s in range(nranks):
        if s == rank:
            continue
        t = np.zeros(S['N'], dtype=bool)
        t[conn[edges[s]:edges[s + 1]].ravel()] = True
        shared = np.intersect1d(my_eq, _eqs_of(eq_ids, t), assume_unique=True)
        if shared.size:
            neighbours.append((s, shared))
            halo_list.append(shared)
    halo_idx = (np.unique(np.concatenate(halo_list)) if halo_list
                else np.zeros(0, dtype=np.int64))
    pos = {int(e): i for i, e in enumerate(halo_idx)}
    neighbours = [(s, np.array([pos[int(e)] for e in sh], dtype=np.int64))
                  for s, sh in neighbours]

    # --- fault rows this rank owns (lowest rank that touches both nodes) ----
    fault_touched = touched[finv['nsmp1']] & touched[finv['nsmp2']]
    owner = np.full(finv['nftnd'], -1, dtype=np.int64)
    for s in range(nranks):
        t = np.zeros(S['N'], dtype=bool)
        t[conn[edges[s]:edges[s + 1]].ravel()] = True
        cand = t[finv['nsmp1']] & t[finv['nsmp2']]
        owner = np.where((owner < 0) & cand, s, owner)
    if np.any(owner < 0):
        raise ValueError(
            'MPI4NodalQuant.decompose: %d of %d fault nodes are touched by no '
            'rank -- they would be missing from every frt.txt and the gate '
            'would compare a short file without saying so.'
            % (int(np.count_nonzero(owner < 0)), finv['nftnd']))
    finv_l = _restrict_fault(finv, fault_touched, int(finv['nftnd']))
    fault_rows = np.nonzero(owner == rank)[0]
    # Rows this rank OWNS, expressed as positions within the rows it COMPUTES
    # (finv_l), because that is what the solver will hand back.
    computed = np.nonzero(fault_touched)[0]
    own_in_computed = np.nonzero(np.isin(computed, fault_rows))[0]

    # --- global -> local renumbering ---------------------------------------
    # WHY THIS EXISTS. Restricting the index ARRAYS (above) shrinks the work;
    # it does not shrink the arrays those indices address. The carry stayed at
    # global extent -- v1+velArr+dispArr+force = 97.75 MB on test.tpv104,
    # BYTE-IDENTICAL at 1 rank and at 32, 98.5% of the 32-rank carry -- so
    # compute fell with rank count and memory traffic did not, and from 4
    # ranks up the step was memory-system-bound (measured: 32 concurrent
    # INDEPENDENT zero-communication processes reproduce the in-situ per-step
    # cost to 3.7%, which excludes MPI, barrier and exchange outright).
    # Fortran goes SUPERLINEAR to 37.33x on 32 ranks precisely because its
    # per-rank working set shrinks into cache. Renumbering into the local node
    # and equation sets is what makes that available here.
    #
    # The local sets are exactly what the restriction above already computed:
    # `touched` (the nodes this rank's elements reach, own plus halo) and
    # `my_eq` (their equations). Nothing new is derived.
    #
    # THE SINK STAYS AT INDEX 0. Equation 0 is the no-equation slot every
    # masked-out contribution is scattered into and calcHourglassResist
    # scrubs, mass[0] is its dummy 1.0, and driver's mass divide skips it by
    # slicing [1:]. Mapping it anywhere else would either scrub a real
    # equation or divide the sink by a real mass.
    nodes_l = np.nonzero(touched)[0]
    g2l_node = np.full(S['N'], -1, dtype=np.int64)
    g2l_node[nodes_l] = np.arange(nodes_l.shape[0])
    eqs_l = my_eq
    g2l_eq = np.full(int(S['NEQ']) + 1, -1, dtype=np.int64)
    g2l_eq[0] = 0
    g2l_eq[eqs_l] = np.arange(1, eqs_l.shape[0] + 1)

    n_node_l = int(nodes_l.shape[0]); n_eq_l = int(eqs_l.shape[0])
    _relabel_all(inv_l, _INDEX_SPACE, g2l_node, g2l_eq, n_node_l, n_eq_l)
    _relabel_all(finv_l, _FAULT_INDEX_SPACE, g2l_node, g2l_eq,
                 n_node_l, n_eq_l)
    inv_l['N'] = n_node_l
    inv_l['NEQ'] = n_eq_l
    inv_l['NEQ1'] = n_eq_l + 1
    halo_global = halo_idx
    halo_local = _relabel('halo_idx', halo_idx, g2l_eq, 'eq', n_eq_l)

    report = dict(rank=rank, nranks=nranks, elem_lo=lo, elem_hi=hi,
                  Ei=inv_l['Ei'], Ep=inv_l['Ep'], E=inv_l['E'],
                  work=float(w[lo:hi].sum()), work_total=float(w.sum()),
                  nodes=int(touched.sum()), eqs=int(my_eq.size),
                  halo_eqs=int(halo_idx.size),
                  halo_frac=float(halo_idx.size) / max(int(my_eq.size), 1),
                  neighbours=[int(s) for s, _ in neighbours],
                  fault_computed=int(computed.size), fault_owned=int(fault_rows.size),
                  N_local=n_node_l, NEQ_local=n_eq_l,
                  N_global=int(S['N']), NEQ_global=int(S['NEQ']))
    return dict(inv=inv_l, finv=finv_l, halo_idx=halo_local.astype(np.int32),
                halo_eq_global=halo_global, local_nodes=nodes_l,
                local_eqs=eqs_l,
                neighbours=neighbours, fault_rows=fault_rows,
                fault_computed_rows=computed,
                own_in_computed=own_in_computed, report=report)


def restrict_rows(d, rows, n):
    """`d` with every per-fault-node array reduced to `rows`. Used for the
    thermal-pressurization constants and history (friclaw 5), which are
    built per fault node by updateThermalPressurization.build and must follow
    the same restriction as finv or the two disagree on which node is which."""
    mask = np.zeros(n, dtype=bool)
    mask[rows] = True
    return _restrict_fault(d, mask, n)


def _eqs_of(eq_ids, node_mask):
    """Sorted unique equation indices of the flagged nodes, sink (0) dropped.
    Index 0 is the no-equation sink every masked-out contribution is scattered
    into and calcHourglassResist scrubs; exchanging it would sum garbage."""
    e = eq_ids[node_mask].ravel()
    return np.unique(e[e > 0])


def _restrict_fault(d, mask, n):
    """`d` with its per-fault-node rows (leading axis == n) reduced to `mask`.
    Scalars pass through. `tr` is either None (no nucleation for this case) or
    one value per fault node -- both handled, because a wrong guess here
    silently disables nucleation, the exact defect rule 17 step 3 records."""
    out = {}
    for k, v in d.items():
        if k == 'nftnd':
            out[k] = int(mask.sum())
        elif isinstance(v, (list, tuple)) and len(v) and \
                all(hasattr(x, 'shape') and x.shape[:1] == (n,) for x in v):
            out[k] = [np.asarray(x)[mask] for x in v]
        elif hasattr(v, 'shape') and getattr(v, 'shape', (0,))[:1] == (n,):
            out[k] = np.asarray(v)[mask]
        else:
            out[k] = v
    return out


SYNC_ENV = 'EQDYNA_MPI_SYNC'
SYNCS = ('halo', 'allreduce')


def sync_mode():
    """'halo' (default: O(boundary), MPI4NodalQuant's own pattern) or
    'allreduce' (O(NEQ), no ownership bookkeeping at all).

    'allreduce' IS NO LONGER RUNNABLE and driver.run_mpi refuses it by name.
    It reduced the FULL nodal array, which required every rank to hold that
    array at global extent -- exactly the replication the rank-local
    renumbering removed. Under local extent each rank's force has a different
    length, and an Allreduce over unequal buffers is not a smaller version of
    the same thing, it is undefined. The value is still PARSED, and still
    rejected when misspelled, so an old command line fails with a sentence
    instead of silently taking the halo path under an allreduce label.

    Both are CORRECT -- they compute the same sum, and the choice is purely
    how much is moved. They are not bit-identical to each other: 'halo' sums
    own + neighbours in ascending rank order, 'allreduce' lets the MPI
    implementation choose, so the last bits differ. NO FALLBACK on a
    misspelled value: a run that silently took the O(NEQ) path while labelled
    halo would make the whole comparison meaningless."""
    import os
    s = os.environ.get(SYNC_ENV, 'halo')
    if s not in SYNCS:
        raise ValueError('%s=%r: must be one of %r' % (SYNC_ENV, s, SYNCS))
    return s


STEP_PROFILE_ENV = 'EQDYNA_MPI_STEP_PROFILE'


def step_profile():
    """True when driver.run_mpi must attribute each step's wall time to jitted
    compute, the barrier, the MPI call, the halo device-to-host copy, and the
    host-side remainder. OFF by default, and the default must stay off: this is
    a measurement knob, not a feature.

    Why it exists: the per-step cost of this path sits ~103-126 ms ABOVE T1/N
    at every rank count >= 4, and transport (1.33-11.08 ms of a 132-172 ms
    step), placement (EFFECTIVE_CORES 1.00 on 31 of 32 ranks) and element
    balance (recut from a 3.350x work spread to 1.004x moved the 32-rank point
    by 1.4%) have each been measured and each failed to account for it. A
    residual that has survived three refutations gets attributed, not guessed
    at a fourth time.

    NO FALLBACK on a misspelled value: a run that silently profiled nothing
    while labelled as profiling would print an empty attribution that reads as
    "no overhead found"."""
    import os
    v = os.environ.get(STEP_PROFILE_ENV, '0')
    if v not in ('0', '1'):
        raise ValueError('%s=%r: must be "0" or "1"' % (STEP_PROFILE_ENV, v))
    return v == '1'


def exchange(comm, neighbours, halo_vals):
    """MPI4NodalQuant's nodal sum, on the HALO VALUES only.

    `halo_vals` is this rank's partial sums at its shared equations, already
    on the host. Returns the DELTA to add (neighbours' partials), so the
    caller adds once on the device and the local values are never sent back
    through a second conversion.

    Neighbours are walked in ASCENDING RANK ORDER and summed in that order,
    so the result is reproducible run to run for a given rank count. The
    ACCUMULATION order below is still exactly that -- the transfers are what
    became concurrent, not the additions -- so this is bit-identical to the
    Sendrecv version it replaces, at every rank count.

    WHY NOT Sendrecv, WHICH THIS USED TO BE. A blocking Sendrecv walked in
    ascending rank order SERIALISES a chain decomposition. Rank r's first
    exchange is with r-1, but r-1's first is with r-2, so r cannot start until
    r-1 has finished with r-2, and the chain completes as a wave of nranks-1
    sequential pair transfers instead of two concurrent ones. That is not a
    mechanism argued from the code: t_mpi at 32 ranks on test.tpv104 rose
    MONOTONICALLY with rank index -- rank 0 1.11 ms, rank 1 1.74, rank 15
    4.50, rank 16 4.51, rank 31 8.38, mean 4.66 of a 43.39 ms differenced
    step -- which is the wave, measured. The old docstring's 0.116 ms for a
    160 kB ring at 16 ranks was a per-PAIR cost and did not survive to 32.

    Irecvs are posted BEFORE any Isend, so no transfer can block on an
    unposted receive and the pattern cannot deadlock for any neighbour graph,
    not merely for the two-neighbour slab case."""
    delta = np.zeros_like(halo_vals)
    if not neighbours:
        return delta
    recvs = [np.empty(pos.shape[0], dtype=halo_vals.dtype)
             for _, pos in neighbours]
    reqs = [comm.Irecv(r, source=s) for r, (s, _) in zip(recvs, neighbours)]
    sends = [np.ascontiguousarray(halo_vals[pos]) for _, pos in neighbours]
    reqs += [comm.Isend(b, dest=s) for b, (s, _) in zip(sends, neighbours)]
    from mpi4py import MPI
    MPI.Request.Waitall(reqs)
    # ASCENDING RANK ORDER, as before: `neighbours` is built ascending by
    # decompose and this loop does not reorder it.
    for (s, pos), recv in zip(neighbours, recvs):
        delta[pos] += recv
    return delta
