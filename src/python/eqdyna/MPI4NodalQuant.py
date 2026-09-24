"""MPI4NodalQuant.py <- assembleGlobalMass.f90:58-298 (the MPI4NodalQuant
and processNodalQuantArr subroutines) + the rank bookkeeping meshgen.f90 does
from npx/npy/npz (calcXyzMPIId, the per-face neighbour ranks).

WHY THIS FILE HAS A SUBROUTINE'S NAME AND NOT A FILE'S. Every other module
here is named after the Fortran FILE it ports. MPI4NodalQuant is a subroutine
inside assembleGlobalMass.f90, and the nodal exchange -- plus the partition it
exchanges across -- is the one thing this code is. It is named after the
Fortran entity it ports.

WHAT THIS IS FOR. python-jax-mpi is the Fortran decomposition, one OS PROCESS
per rank, with jax owning only the local element kernel. XLA automatic
partitioning and explicit shard_map were both measured first and both lose to
it (0 of 38 HLO scatters partitioned; the shard_map route replicates the nodal
stages and all-reduces O(NEQ)); see backend.run_time_loop_sharded.

THE DECOMPOSITION IS FORTRAN'S (pathway item 64, owner decisions 2026-09-24:
"follow Fortran"). Each rank builds ONLY its own (x,y,z) box, with rank-local
numbering, exactly as meshgen.f90 does:

  * the split is (npx,npy,npz) = DECOMP[nranks] below -- the same table the
    perf tools hand the Fortran binary (testsys/perf/run_scaling.py imports
    it from here, so the two cannot drift);
  * rank -> (mex,mey,mez) is calcXyzMPIId, z fastest (meshgen.calc_xyz_mpi_id);
  * each dimension's global grid line is split by getLocalOneDimCoorArrAndSize
    (meshgen.partition_1d): neighbouring boxes SHARE one node plane and never
    an element, and the local line is a bitwise SLICE of the global one;
  * every node in the box is local -- including both halves of a split-node
    pair on a shared face, which Fortran's createMasterNode creates on every
    rank that holds the plane. A fault node pair is therefore never split
    across ranks, and a y boundary lying ON the fault plane is legal;
  * elements are the box's own ix,iy,iz >= 2, so the element sets partition.

Nothing here builds or holds the global mesh. The only global-extent objects
a rank keeps are the three 1D grid lines (O(n^(1/3))) and two O(1) integer
censuses (meshgen.fault_census / equation_census) used to check that the
ranks' owned entities add up.

THE EXCHANGE, verbatim in ORDER. A node on a shared plane holds a PARTIAL sum
on each rank that holds it. MPI4NodalQuant sums them the way
assembleGlobalMass.f90:95-217 does: three phases, x then y then z; in each,
the minus face then the plus face; per face, fetch this rank's values at the
face's regular nodes (loop order below) and -- when fltMPI(k), i.e. this rank
has fault nodes on that face -- at the master nodes on it
(addFaultBoundaryTerm), Sendrecv with the face neighbour, add what came back.
An edge or corner contribution reaches its other sharers by RELAY (the
x-phase result is what the y-phase sends), not by a 26-neighbour exchange.
That relay is load-bearing, not incidental: with it every sharer ends with
the same bits (IEEE addition is commutative, so a+b == b+a and
(a+b)+(c+d) == (c+d)+(a+b)), which is what lets a shared node -- and a
split-node pair on a shared face -- be integrated REDUNDANTLY by every rank
that holds it, from identical complete forces, with no velocity exchange.
A one-round own+sum(neighbours) exchange does not have that property at 3+
sharers.

Used for nodalMassArr (numDof 3) and fnms (numDof 1) once at setup, as
assembleGlobalMass.f90:41-42, and for the nodal force every step, as
driver.f90:27. Fortran's mpi_barrier after each phase is a timing device and
not reproduced; the Sendrecv chain is deadlock-free without it (each phase
completes from rank 0 upward along every line).

BIT-EXACTNESS. An N-rank run is not bit-identical to a serial one and cannot
be: a shared equation's partials are summed in a different order. It is gated
against the committed canonical reference at the case bound, like every
other column, and is reproducible run to run for a given decomposition.
"""
import numpy as np

from . import meshgen

# (npx, npy, npz) per rank count -- the Fortran perf/e2e decomposition. ONE
# copy: testsys/perf/run_scaling.py imports this name.
DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2),
          16: (4, 2, 2), 32: (4, 4, 2)}


class Partition(object):
    """This rank's place in an (npx,npy,npz) decomposition. Pure bookkeeping:
    no communicator, so a test can build any rank's mesh in one process."""

    def __init__(self, rank, nranks, dims):
        dims = tuple(int(d) for d in dims)
        if len(dims) != 3 or min(dims) < 1 or dims[0] * dims[1] * dims[2] != nranks:
            raise ValueError('Partition: decomposition %r does not multiply to '
                             '%d ranks' % (dims, nranks))
        if not 0 <= rank < nranks:
            raise ValueError('Partition: rank %d out of range for %d ranks'
                             % (rank, nranks))
        self.rank, self.nranks, self.dims = int(rank), int(nranks), dims
        self.mexyz = meshgen.calc_xyz_mpi_id(self.rank, *dims)

    @classmethod
    def for_size(cls, rank, nranks):
        """The decomposition for `nranks` is DECOMP's -- never inferred, never
        a nearest match. A rank count the table does not name is refused."""
        if nranks not in DECOMP:
            raise ValueError(
                'python-jax-mpi runs at the rank counts MPI4NodalQuant.DECOMP '
                'names (%s), the same (npx,npy,npz) table the Fortran is run '
                'with; got %d ranks.' % (sorted(DECOMP), nranks))
        return cls(rank, nranks, DECOMP[nranks])

    def neighbour(self, d, ib):
        """Face neighbour in dimension d (0=x,1=y,2=z), ib 0 = minus side,
        1 = plus side: me -/+ npy*npz, npz, 1 (assembleGlobalMass.f90:77-90)."""
        stride = (self.dims[1] * self.dims[2], self.dims[2], 1)[d]
        return self.rank - stride if ib == 0 else self.rank + stride

    def at_model_edge(self, d, ib):
        """bnd(ib)==0 in the Fortran: the minus face of the first rank and the
        plus face of the last rank along d are model boundary, not halo."""
        return self.mexyz[d] == 0 if ib == 0 else self.mexyz[d] == self.dims[d] - 1

    def slice_lines(self, xline, yline, zline):
        """Local grid lines (bitwise slices of the global ones) and their
        0-based global offsets."""
        out, offs = [], []
        for d, line in enumerate((xline, yline, zline)):
            loc, off = meshgen.local_line(line, self.dims[d], self.mexyz[d])
            out.append(loc)
            offs.append(off)
        return out, tuple(offs)


def face_node_ids(d, ib, nx, ny, nz):
    """1-based LOCAL regular node ids on face (d, ib), in MPI4NodalQuant's own
    loop order (assembleGlobalMass.f90:152-172 / :190-210): x faces loop iz
    then iy, y faces ix then iz, z faces ix then iy. Both sides of a face
    walk it in this order over the same (shared) local extents, which is
    what makes the Sendrecv buffers line up -- checked at setup by
    `handshake`, not assumed."""
    bnd = 1 if ib == 0 else (nx, ny, nz)[d]
    if d == 0:
        iz, iy = np.meshgrid(np.arange(1, nz + 1), np.arange(1, ny + 1), indexing='ij')
        ids = (bnd - 1) * ny * nz + (iz - 1) * ny + iy
    elif d == 1:
        ix, iz = np.meshgrid(np.arange(1, nx + 1), np.arange(1, nz + 1), indexing='ij')
        ids = (ix - 1) * ny * nz + (iz - 1) * ny + bnd
    else:
        ix, iy = np.meshgrid(np.arange(1, nx + 1), np.arange(1, ny + 1), indexing='ij')
        ids = (ix - 1) * ny * nz + (bnd - 1) * ny + iy
    return ids.ravel().astype(np.int64)


def build_faces(part, n_local, ndof, eq_ids, flt_lists, flt_mpi):
    """The exchange plan: one entry per face this rank exchanges, in the
    Fortran's order (x-, x+, y-, y+, z-, z+; model-boundary faces and
    undivided dimensions skipped).

    Each entry carries `nodes` (1-based local node ids: the face's regular
    nodes, then -- when fltMPI(k) -- nx*ny*nz + its fault-node indices, i.e.
    addFaultBoundaryTerm's master nodes) and `eqs` (processNodalQuantArr's
    numDof==3 slots: every eq>0 of those nodes in dof order). `ndof`/`eq_ids`
    are the 0-indexed-by-node S['ndof']/S['eq_ids'] tables (sink 0 = no
    equation)."""
    nx, ny, nz = n_local
    n_reg = nx * ny * nz
    faces = []
    for d in range(3):
        if part.dims[d] <= 1:
            continue
        for ib in (0, 1):
            if part.at_model_edge(d, ib):
                continue
            k = 2 * d + ib
            nodes = face_node_ids(d, ib, nx, ny, nz)
            if flt_mpi[k]:
                nodes = np.concatenate((nodes, n_reg + flt_lists[k]))
            if np.unique(nodes).size != nodes.size:
                raise ValueError('build_faces: face %d holds a node twice' % k)
            e = eq_ids[nodes - 1]
            live = (np.arange(e.shape[1])[None, :] < ndof[nodes - 1][:, None]) & (e > 0)
            faces.append(dict(d=d, ib=ib, k=k, nb=part.neighbour(d, ib),
                              nodes=nodes, eqs=e[live].astype(np.int64),
                              eqs_per_node=live.sum(axis=1).astype(np.int64)))
    return faces


def _tags(part, face, num_dof):
    """MPI4NodalQuant's sendtag/recvtag (assembleGlobalMass.f90:127-141)."""
    base = (0, 10000, 20000)[face['d']] * num_dof
    return base + part.rank, base + face['nb']


def relay(comm, part, faces, key, arr, num_dof):
    """MPI4NodalQuant(arr, numDof) on a 1-indexed builder array, IN PLACE:
    for each face in order, fetch arr[face[key]], Sendrecv with the face
    neighbour, add the neighbour's values back (processNodalQuantArr's
    fetch/add). `key` is 'eqs' for numDof 3 (nodalMassArr, indexed by
    equation) and 'nodes' for numDof 1 (fnms, indexed by node)."""
    for f in faces:
        idx = f[key]
        send = np.ascontiguousarray(arr[idx])
        recv = np.empty_like(send)
        st, rt = _tags(part, f, num_dof)
        comm.Sendrecv(send, dest=f['nb'], sendtag=st, recvbuf=recv,
                      source=f['nb'], recvtag=rt)
        arr[idx] += recv
    return arr


def handshake(comm, part, faces, meshCoor):
    """Once per run, per face: exchange the face's length, per-node equation
    counts and node COORDINATES with the neighbour and require bitwise
    equality. This is the replacement for the global-id symmetry check the
    build-then-restrict design had: without it a transposed loop on one side
    would sum the wrong partials and nothing downstream could attribute it.
    `meshCoor` is the 1-indexed (N+1, 3) builder array. Tags 70000+ sit
    above every MPI4NodalQuant tag (at most 20000*3 + rank)."""
    for f in faces:
        mine = np.concatenate((np.array([f['nodes'].size, f['eqs'].size], dtype=np.float64),
                               f['eqs_per_node'].astype(np.float64),
                               meshCoor[f['nodes']].ravel()))
        n_theirs = np.empty(1)
        comm.Sendrecv(np.array([float(mine.size)]), dest=f['nb'], sendtag=70000 + part.rank,
                      recvbuf=n_theirs, source=f['nb'], recvtag=70000 + f['nb'])
        if int(n_theirs[0]) != mine.size:
            raise RuntimeError(
                'MPI4NodalQuant.handshake: rank %d face %d carries %d values, '
                'neighbour %d carries %d -- the two sides do not describe the '
                'same shared plane.' % (part.rank, f['k'] + 1, mine.size, f['nb'],
                                        int(n_theirs[0])))
        theirs = np.empty_like(mine)
        comm.Sendrecv(mine, dest=f['nb'], sendtag=71000 + part.rank, recvbuf=theirs,
                      source=f['nb'], recvtag=71000 + f['nb'])
        if not np.array_equal(mine, theirs):
            nbad = int(np.count_nonzero(mine != theirs))
            raise RuntimeError(
                'MPI4NodalQuant.handshake: rank %d face %d and neighbour %d '
                'disagree in %d of %d entries (node count, eq count, per-node '
                'eq counts, coordinates). A Sendrecv over these buffers would '
                'add partials of DIFFERENT nodes.'
                % (part.rank, f['k'] + 1, f['nb'], nbad, mine.size))


def owned_mask(part, ids, n_local):
    """True for each 1-based local REGULAR node id in `ids` that this rank
    OWNS: the lowest rank holding it. A node on a shared plane is held by the
    ranks on both sides; since rank = mex*npy*npz + mey*npz + mez is monotone
    in each coordinate, the lowest holder is the one that is lowest in every
    dimension, i.e. this rank owns it unless it sits on this rank's minus
    face of a dimension in which this rank is not first. Pure arithmetic, the
    same answer on every rank with no communication. Used ONLY to decide who
    WRITES a fault row and who COUNTS an equation -- never who computes one
    (every holder computes, as Fortran does)."""
    nx, ny, nz = n_local
    s = np.asarray(ids, dtype=np.int64) - 1
    loc = (s // (nz * ny), s % ny, (s % (nz * ny)) // ny)     # (ix, iy, iz)
    own = np.ones(s.shape, dtype=bool)
    for d in range(3):
        if part.mexyz[d] > 0:
            own &= loc[d] > 0
    return own


def exchange(comm, part, faces, vals):
    """The per-step nodal-force MPI4NodalQuant (driver.f90:27), on the HALO
    VALUES only. `vals` is this rank's force at every face equation (the
    sorted union `halo_idx`, already on the host); each face's `pos` indexes
    into it. Returns the COMPLETE values, which the caller writes back over
    the device array -- a set, not an add of a delta, because own+delta is
    not own+recv in floating point once a relay has summed three terms."""
    for f in faces:
        pos = f['pos']
        send = np.ascontiguousarray(vals[pos])
        recv = np.empty_like(send)
        st, rt = _tags(part, f, 3)
        comm.Sendrecv(send, dest=f['nb'], sendtag=st, recvbuf=recv,
                      source=f['nb'], recvtag=rt)
        vals[pos] += recv
    return vals


def setup_exchange(comm, part, S, mesh):
    """Everything the Fortran exchanges at SETUP, in its order, then the
    per-step plan:

      handshake                        -- the shared planes really are shared
      MPI4arn (meshgen.f90:156)        -- arn, with the DIVIDE/DUPLICATE rule
      MPI4NodalQuant(nodalMassArr, 3)  -- assembleGlobalMass.f90:41
      MPI4NodalQuant(fnms, 1)          -- assembleGlobalMass.f90:42

    and re-binds S['arn'] / S['nodalMassArr'] / S['fnms'] to the completed
    values. Returns the plan dict driver.run_mpi steps with."""
    n_local = mesh['n_local']
    arn1 = mesh['arn1']
    flt_lists = mesh['flt_lists']
    flt_mpi = meshgen.flt_mpi_flags(part, flt_lists)
    faces = build_faces(part, n_local, S['ndof'], S['eq_ids'], flt_lists, flt_mpi)
    handshake(comm, part, faces, mesh['meshCoor'])
    if meshgen.mpi4arn(comm, part, arn1, flt_lists, mesh['fault_box']) != flt_mpi:
        raise RuntimeError('MPI4NodalQuant.setup_exchange: MPI4arn set fltMPI '
                           'differently from the plan it was built against')
    mass1, fnms1 = mesh['mass1'], mesh['fnms1']
    relay(comm, part, faces, 'eqs', mass1, 3)
    relay(comm, part, faces, 'nodes', fnms1, 1)
    # assemble_mass's invariant, which the two relays must preserve because
    # they add the same neighbour values in the same order: nodalMassArr[eq]
    # == fnms[node] bit for bit for every live equation. A mis-paired eq list
    # breaks it; checked, not assumed.
    e = S['eq_ids']
    live = (np.arange(e.shape[1])[None, :] < S['ndof'][:, None]) & (e > 0)
    node_of = np.nonzero(live)[0] + 1
    if not np.array_equal(mass1[e[live]], fnms1[node_of]):
        raise RuntimeError('MPI4NodalQuant.setup_exchange: after the relay, '
                           'nodalMassArr and fnms disagree at %d equation(s).'
                           % int(np.count_nonzero(mass1[e[live]] != fnms1[node_of])))
    S['arn'] = arn1[1:]
    S['nodalMassArr'] = mass1[1:]
    S['fnms'] = fnms1[1:]

    halo_idx = (np.unique(np.concatenate([f['eqs'] for f in faces]))
                if faces else np.zeros(0, dtype=np.int64))
    for f in faces:
        f['pos'] = np.searchsorted(halo_idx, f['eqs'])
    # Rows this rank WRITES: every local fault row is COMPUTED here (both
    # halves of each pair are local), but each is written by exactly one
    # rank -- the lowest holder -- so the frt files partition the fault and
    # frt_canonical's duplicate branch never fires for this backend.
    nsmp = mesh['nsmp']
    owned_rows = (np.nonzero(owned_mask(part, nsmp[:, 0], n_local))[0]
                  if nsmp.shape[0] else np.zeros(0, dtype=np.int64))
    _check_censuses(comm, part, S, mesh, owned_rows, live)
    return dict(faces=faces, halo_idx=halo_idx.astype(np.int32),
                fault_rows=owned_rows, flt_mpi=flt_mpi)


def _check_censuses(comm, part, S, mesh, owned_rows, live):
    """The partition's conservation checks, O(1) communication each, run on
    every MPI run (design 56d2401 section 1.2, R1). The ranks' OWNED fault
    nodes must be the global set -- count AND sum of global grid keys, so a
    node owned twice plus one owned never cannot pass -- and their OWNED
    equations must add up to the global totalNumOfEquations. Both globals
    come from meshgen's streaming censuses over the global grid lines, not
    from any rank's mesh, so a numbering or boundary-classification defect in
    the rank-local build has nothing to agree with by construction."""
    nx, ny, nz = mesh['n_local']
    nxg, nyg, nzg = mesh['n_global']
    ox, oy, oz = mesh['offsets']
    s0 = mesh['nsmp'][owned_rows, 0].astype(np.int64) - 1
    keys = ((ox + s0 // (nz * ny)) * nzg * nyg + (oz + (s0 % (nz * ny)) // ny) * nyg
            + (oy + s0 % ny))
    got = (comm.allreduce(int(owned_rows.size)), comm.allreduce(int(keys.sum())))
    if got != tuple(mesh['fault_census']):
        raise RuntimeError(
            'MPI4NodalQuant: the ranks OWN fault nodes (count, key sum) = %r, '
            'the global grid has %r. frt.txt* would be short or doubled and '
            'the canonical comparison would not say so.'
            % (got, tuple(mesh['fault_census'])))
    n_reg = nx * ny * nz
    per_node = live.sum(axis=1)
    own_reg = owned_mask(part, np.arange(1, n_reg + 1), mesh['n_local'])
    owned_eqs = int(per_node[:n_reg][own_reg].sum()) + int(per_node[n_reg + owned_rows].sum())
    total = comm.allreduce(owned_eqs)
    if total != mesh['equation_census']:
        raise RuntimeError(
            'MPI4NodalQuant: the ranks OWN %d equations, the global mesh has '
            '%d (countMeshEntities). A subdomain face was classified as model '
            'boundary, or a shared node was dropped or doubled.'
            % (total, mesh['equation_census']))


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
