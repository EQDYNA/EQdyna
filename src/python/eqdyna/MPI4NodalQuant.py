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
#
# THIS WAS 3.0 AND 3.0 WAS WRONG BY A FACTOR OF ~5. The reasoning behind it --
# calcPMLElemKU computes 15 split stress components and 12 force blocks against
# calcElemKU's 6 and 3, so a PML element must be several times an interior one
# -- counts arithmetic and ignores that on the cases this path is measured on
# the INTERIOR kernel additionally runs plasticity (test.tpv104 is viscoplastic;
# calcElemKU's plastic branch touches interior elements only). decompose() did
# report (Ei, Ep) exactly so this could be checked, and when it was checked it
# did not hold:
#
#   MEASURED, test.tpv104, 32 ranks, halo sync, per-rank by difference over
#   40/160 steps with the barrier wait and the exchange subtracted
#   (docs/perf_snapshots/mpi_scaling_2026-09-22_012512_tpv104_1to32_cacheoff.json):
#     the 4 all-PML ranks (Ei=0, Ep=12217)  did  36.5 ms/step of own work
#     the 26 mixed ranks (Ei~19227, Ep~5808) did 114.8 ms/step
#   -> b = 2.99e-3 ms/PML element, a = 5.07e-3 ms/interior element, b/a = 0.59.
#   A no-intercept least squares over all 32 (Ei, Ep, own) triples of the same
#   run gives 0.563. The same fit at 16 ranks gives 1.07-1.16 and at 8 ranks
#   1.08 -- the spread is real and is why this is stated as a RANGE (0.56-1.16,
#   i.e. "about one", never three) rather than as a 3-digit constant: a
#   no-intercept element model absorbs each rank's nodal work into the element
#   coefficients, and the nodal share per element grows with the rank count.
#
# The value below was then chosen by DIRECT MEASUREMENT of the candidates at
# 16 and 32 ranks, not by taking the fit at its word -- see NOTES_mpi_balance.md
# M3 for the table. It is still a load-balance heuristic and not physics, which
# is why decompose() reports it (report['pml_weight']) beside the (Ei, Ep) it
# produced, and why EQDYNA_PML_WEIGHT can override it for a recalibration on a
# case whose interior kernel is cheaper (a non-plastic case: friclaw 1 without
# viscoplasticity does less interior work, so its true ratio is higher).
PML_WEIGHT = 0.6
PML_WEIGHT_ENV = 'EQDYNA_PML_WEIGHT'


def pml_weight():
    """PML_WEIGHT, or EQDYNA_PML_WEIGHT when set. Strictly positive and
    finite. NO FALLBACK on an unparseable value: a silently ignored weight
    would produce a partition nobody asked for, and the only symptom would be
    a speedup number that does not reproduce."""
    import os
    raw = os.environ.get(PML_WEIGHT_ENV)
    if raw is None:
        return PML_WEIGHT
    try:
        w = float(raw)
    except ValueError:
        raise ValueError('%s=%r is not a number. It is the relative cost of a '
                         'PML element against an interior one and it chooses '
                         'the element cut points.' % (PML_WEIGHT_ENV, raw))
    if not (w > 0.0) or w == float('inf'):
        raise ValueError('%s=%r must be > 0 and finite; a non-positive or '
                         'infinite weight makes the cumulative weight '
                         'non-increasing and the cuts meaningless.'
                         % (PML_WEIGHT_ENV, raw))
    return w


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


def _require_interior_on_every_rank(elemType, edges, nranks, pw):
    """Refuse a partition in which some rank owns 0 INTERIOR elements.

    This is a different condition from `_cuts`'s "0 elements", and it was the
    whole defect. The mesh generator's nested loop emits the two x-face PML
    slabs as contiguous blocks at the ENDS of the element array (test.tpv104:
    elements 0..31955 and 703494..734999 are all PML, 63462 of 735000), so an
    over-weighted PML makes the first and last cuts land inside those blocks
    and hands those ranks a slab that is pure PML. At PML_WEIGHT=3.0 that was
    4 of 32 ranks (and 2 of 16), and those 4 ranks did 36.5 ms/step of own
    work against the other 26 ranks' 114.8 -- so they sat at the barrier for
    95-105 ms of a 172.6 ms step and set the whole run's speedup to 5.28x.

    WHY THIS IS AN INVARIANT AND NOT A WARNING. A zero-interior rank is not
    wrong, it is UNDER-LOADED, and under-loading is invisible in every output
    the run produces -- the physics is identical and only the wall clock moves.
    A warning in a 32-rank log is a warning nobody reads. The arithmetic is
    also exact and worth stating in the message, because it says which way to
    turn the knob: with a contiguous partition, rank 0 holds interior elements
    only if its weight share exceeds the leading all-PML block's weight, i.e.
    only if  n_lead * pw < (Ei + Ep*pw) / nranks.

    NOT SILENTLY REPAIRED. Nudging the cut to swallow one interior element
    would satisfy the letter of this check and leave the rank under-loaded by
    the same 3x, which is exactly the "green result that tested nothing" shape
    this repo keeps writing down. The weight is the thing that is wrong, so
    the weight is what the message points at."""
    is_int = (elemType == 1) | (elemType > 10)
    ci = np.concatenate(([0], np.cumsum(is_int)))
    empty = [r for r in range(nranks)
             if ci[edges[r + 1]] - ci[edges[r]] == 0]
    if not empty:
        return
    ei = int(is_int.sum())
    ep = int(elemType.shape[0]) - ei
    lead = int(np.argmax(is_int)) if is_int.any() else 0
    raise ValueError(
        'MPI4NodalQuant: rank(s) %s of %d own 0 INTERIOR elements at '
        'PML_WEIGHT=%g (Ei=%d, Ep=%d, %d leading all-PML elements). Such a '
        'rank still runs and still produces correct physics -- it is simply '
        'under-loaded, and every other rank waits for it at the barrier, '
        'which is invisible in the output. Rank 0 gets interior elements only '
        'while %d*%g = %g < (%d + %d*%g)/%d = %g, so lower %s (measured '
        '0.56-1.16 on test.tpv104, NOT 3) or use fewer ranks.'
        % (empty, nranks, pw, ei, ep, lead,
           lead, pw, lead * pw, ei, ep, pw, nranks, (ei + ep * pw) / nranks,
           PML_WEIGHT_ENV))


def _take(a, sel):
    return np.asarray(a)[sel]


def _take_ravelled(a, sel):
    """The (n*8,) index arrays assembleGlobalKU.build already ravelled:
    reshaped, row-selected, re-ravelled, so the kernel receives exactly the
    flat layout it receives serially and the within-rank scatter order is the
    serial one."""
    return np.asarray(a).reshape(-1, 8)[sel].reshape(-1)


def decompose(S, inv, finv, rank, nranks):
    """Rank-local `inv`, `finv`, exchange plan and output row selection.

    Returns a dict with
        inv, finv      -- the same dicts, element rows and fault-node rows
                          restricted to this rank
        halo_idx       -- equation indices this rank exchanges (int32)
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
    pw = pml_weight()
    w = np.where(elemType == 2, pw, 1.0)
    edges = _cuts(w, nranks)
    lo, hi = edges[rank], edges[rank + 1]
    _require_interior_on_every_rank(elemType, edges, nranks, pw)
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

    report = dict(rank=rank, nranks=nranks, elem_lo=lo, elem_hi=hi,
                  Ei=inv_l['Ei'], Ep=inv_l['Ep'], E=inv_l['E'],
                  pml_weight=pw,
                  work=float(w[lo:hi].sum()), work_total=float(w.sum()),
                  nodes=int(touched.sum()), eqs=int(my_eq.size),
                  halo_eqs=int(halo_idx.size),
                  halo_frac=float(halo_idx.size) / max(int(my_eq.size), 1),
                  neighbours=[int(s) for s, _ in neighbours],
                  fault_computed=int(computed.size), fault_owned=int(fault_rows.size))
    return dict(inv=inv_l, finv=finv_l, halo_idx=halo_idx.astype(np.int32),
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


def allreduce(comm, partial):
    """The simple sync: one MPI_Allreduce(SUM) over the FULL nodal array.

    Correctness needs no ownership information -- assembly is a sum, every
    rank's array is zero where it assembled nothing, and the sum over ranks
    is the complete array on every rank. That is the whole reason to try this
    before the halo: it has no boundary-node lists, no neighbour topology and
    no npx/npy/npz.

    REPRODUCIBILITY IS NOT GUARANTEED HERE and is checked rather than
    assumed: MPI does not promise a fixed reduction order, and float addition
    is not associative, so the last bits may move run to run. See
    testsys/unit/test_mpi_decomposition.py for the measured answer on this
    installation (mpi4py 4.1.2 / OpenMPI) and driver.run_mpi's docstring for
    which mode the gate uses."""
    out = np.empty_like(partial)
    comm.Allreduce(partial, out, op=_MPI().SUM)
    return out


def _MPI():
    from mpi4py import MPI     # ImportError deliberately uncaught
    return MPI


def exchange(comm, neighbours, halo_vals):
    """MPI4NodalQuant's nodal sum, on the HALO VALUES only.

    `halo_vals` is this rank's partial sums at its shared equations, already
    on the host. Returns the DELTA to add (neighbours' partials), so the
    caller adds once on the device and the local values are never sent back
    through a second conversion.

    Neighbours are walked in ASCENDING RANK ORDER and summed in that order,
    so the result is reproducible run to run for a given rank count. Sendrecv
    (not Isend/Irecv) because it cannot deadlock and, at the two neighbours a
    slab decomposition produces, has nothing to gain from overlap: measured
    0.116 ms for a 160 kB ring exchange at 16 ranks, 0.18% of a 64.74 ms
    step."""
    delta = np.zeros_like(halo_vals)
    for s, pos in neighbours:
        send = np.ascontiguousarray(halo_vals[pos])
        recv = np.empty_like(send)
        comm.Sendrecv(send, dest=s, recvbuf=recv, source=s)
        delta[pos] += recv
    return delta
