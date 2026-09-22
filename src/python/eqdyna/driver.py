"""driver.py <- src/driver.f90. THE time loop, and velDispUpdate.

One loop, for every friction law and every backend, in driver.f90's order:

    do nt = 1, nstep
        timeElapsed = timeElapsed + dt
        call velDispUpdate
        nodalForceArr = 0
        call assembleGlobalKU
        call calcHourglassResist
        call faulting                     <- friclaw dispatches INSIDE here
        nodalForceArr = nodalForceArr / nodalMassArr
    enddo

What this replaces: six loops (port.py, port_rsf.py, port_tp.py and a JAX
twin of each), all of them this loop with a different `faulting` inlined and
the friclaw branch hoisted to module level. The friclaw branch belongs where
faulting.f90:21-22 puts it -- inside faulting -- and because friclaw is a
static Python int it is resolved at TRACE time under jit and costs nothing.

MPI is absent by construction: this port is serial (npx=npy=npz=1, enforced
by eqdyna3d.build_solver_state, which raises NotImplementedError otherwise),
so driver.f90:27's MPI4NodalQuant has no work to do.
That is a scope limit, not an omission that could silently mislead -- a
multi-rank case is refused before it reaches here.
"""
import os
import sys
import time

import numpy as np

from . import assembleGlobalKU as KU
from . import backend as B
from . import faulting as FLT
from . import globalvar as gv
from . import updateThermalPressurization as TP


def velDispUpdate(xp, inv, v1, velArr, dispArr, force, dt):
    """driver.f90:92-164 -- integrate velocity and displacement.

    Two node kinds: interior (3 dof) and PML (12 dof, 9 split components
    plus 3 velocity components, with the split components damped by the
    a9/b9 coefficients precomputed in assembleGlobalKU.build).
    """
    idx3_v = inv['idx3_v']
    int_nodes = inv['int_nodes_idx']

    accel3 = force[idx3_v]
    v1 = B.setat(xp, v1, idx3_v, v1[idx3_v] + accel3 * dt)
    # Read back ONCE, AFTER the store, and reuse. The read-back must stay
    # after the store and must NOT be replaced by the stored expression:
    # idx3_v repeats index 0 (the -1 fixed-boundary sentinel is mapped onto
    # the sink slot), so for those entries the store is last-write-wins and
    # what comes back is not what went in. Preserving that is the point.
    v3 = v1[idx3_v]
    velArr = B.setat(xp, velArr, int_nodes, v3)
    dispArr = B.addat(xp, dispArr, int_nodes, v3 * dt)

    idx12_v = inv['idx12_v']
    if inv['pml_nodes_idx'].shape[0]:
        pnodes = inv['pml_nodes_idx']
        f9 = force[idx12_v[:, 0:9]]
        vold9 = v1[idx12_v[:, 0:9]]
        v1 = B.setat(xp, v1, idx12_v[:, 0:9],
                     (f9 + vold9 * inv['a9']) / inv['b9'])
        f3v9 = force[idx12_v[:, 9:12]]
        v1 = B.setat(xp, v1, idx12_v[:, 9:12], v1[idx12_v[:, 9:12]] + f3v9 * dt)
        v1 = B.setat(xp, v1, 0, 0.0)
        # ONE (n_pml,12) gather instead of twelve (n_pml,) gathers of the
        # same v1. No store happens between them, so V[:,k] is element for
        # element the array v1[idx12_v[:,k]] would return.
        V = v1[idx12_v]
        vA = V[:, 0] + V[:, 1] + V[:, 2] + V[:, 9]
        vB = V[:, 3] + V[:, 4] + V[:, 5] + V[:, 10]
        vC = V[:, 6] + V[:, 7] + V[:, 8] + V[:, 11]
        has_eq = idx12_v[:, 0] > 0
        for k, vk in enumerate((vA, vB, vC)):
            velArr = B.setat(xp, velArr, (pnodes, k), xp.where(has_eq, vk, 0.0))
            dispArr = B.setat(xp, dispArr, (pnodes, k),
                              xp.where(has_eq, dispArr[pnodes, k] + vk * dt, 0.0))
    return v1, velArr, dispArr


def _device_peak_gb(jax):
    """This rank's peak DEVICE memory, in GB, from the allocator itself.

    Why not nvidia-smi: XLA preallocates ~75% of the card by default, so
    nvidia-smi reports ~30 GB for a run whose working set is a fraction of
    that -- a figure that answers "is the card busy", not "does this case
    fit". memory_stats()['peak_bytes_in_use'] is the allocator's own
    high-water mark and is the number a capacity claim has to be made on.

    -1.0 means the backend exposes no memory_stats (the CPU backend does
    not). A distinguishable sentinel rather than 0.0, because a printed 0.0
    would read as "this run used no memory" and get quoted as one.
    """
    st = jax.devices()[0].memory_stats()
    if not st or 'peak_bytes_in_use' not in st:
        return -1.0
    return round(st['peak_bytes_in_use'] / 1e9, 3)


# Position of nodalForceArr in the carry tuple. Named because two callers
# reach into the carry to exchange exactly that entry (make_step's in-process
# nodal_sync and run_mpi's out-of-process MPI4NodalQuant), and the same
# integer literal in two places is one place for them to disagree.
FORCE = 3


def make_step_parts(xp, inv, finv, tp, mass, scratch):
    """The step, split at driver.f90:27 -- MPI4NodalQuant's position.

    part_a: timeElapsed, velDispUpdate, zero the force, both element kernels.
    part_b: thermal pressurization, faulting, the mass divide.

    ONE body, split rather than copied, because the two callers need the seam
    in a different place in the STACK, not in the code: make_step closes the
    seam with backend.nodal_sync (identity when serial, a device-mesh
    collective under shard_map) and keeps ONE jitted time loop, while run_mpi
    must leave the jit at the seam to make an MPI call and therefore jits the
    two halves separately. Both perform the same operations, in the same
    order, on the same operands."""
    dt = inv['dt']; rdampk = inv['rdampk']
    tr = finv['tr']
    friclaw = finv['friclaw']

    def part_a(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist) = carry
        timeElapsed = timeElapsed + dt                       # driver.f90:12

        v1, velArr, dispArr = velDispUpdate(xp, inv, v1, velArr, dispArr,
                                            force, dt)

        force = B.setat(xp, force, slice(None), 0.0)         # driver.f90:23
        force, stress_i, s_p = KU.assembleGlobalKU(
            xp, inv, velArr, force, stress_i, s_p, dt, rdampk, scratch)
        force = KU.calcHourglassResist(xp, inv, dispArr, velArr, force, rdampk)

        # driver.f90:27 -- MPI4NodalQuant(nodalForceArr, 3). Identity when the
        # run is serial (one device, one subdomain, nothing to exchange); an
        # all-reduce over the device mesh when backend.run_time_loop_sharded
        # has cut the element arrays across devices. Its POSITION is the
        # Fortran's: after both element kernels, before faulting, which is
        # what lets faulting and the mass divide be plain replicated nodal
        # work on a force array that is already complete.
        return (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft,
                timeElapsed, sliprate_hist, shear_hist)

    def part_b(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist) = carry

        if friclaw == 5:                                   # driver.f90:28
            fric = TP.updateThermalPressurization(
                xp, tp, fric, sliprate_hist, shear_hist, nt, dt)

        fric, fnft, force = FLT.faulting(xp, finv, fric, fnft, velArr, dispArr,
                                         force, dt, timeElapsed, tr, nt)

        if friclaw == 5:
            # onFaultTPHist(1|2, i, nt, ift) -- written AFTER faulting, so the
            # next step's integral sees steps 1..nt. Column nt-1, 0-based.
            sliprate_hist = B.setat(xp, sliprate_hist, (slice(None), nt - 1),
                                    fric[:, gv.PEAK_SLIPRATE])
            shear_hist = B.setat(xp, shear_hist, (slice(None), nt - 1),
                                 fric[:, gv.SHEAR_MAG])

        # driver.f90:30 -- a DIVISION by the lumped mass, not a multiply by a
        # precomputed reciprocal. The jax/rsf/tp ports this replaces used
        # `force * inv_mass_full`, which is a different rounding and a step
        # AWAY from the reference; unifying on the Fortran's division moves
        # those columns' bits and is reported as such rather than hidden.
        # Slot 0 is the no-equation sink and is scrubbed, never divided.
        force = B.setat(xp, force, slice(1, None), force[1:] / mass[1:])
        force = B.setat(xp, force, 0, 0.0)

        return (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft,
                timeElapsed, sliprate_hist, shear_hist)

    return part_a, part_b


def make_step(xp, inv, finv, tp, mass, scratch):
    """Build the per-step closure. Called by backend.run_time_loop INSIDE
    the jit on the jax path, so that `inv`'s arrays resolve to jit arguments
    rather than closed-over HLO literals."""
    part_a, part_b = make_step_parts(xp, inv, finv, tp, mass, scratch)

    def step(carry, nt):
        carry = part_a(carry, nt)
        # driver.f90:27 -- MPI4NodalQuant(nodalForceArr, 3). Identity when the
        # run is serial (one subdomain, nothing to exchange); an all-reduce
        # over the device mesh when backend.run_time_loop_sharded has cut the
        # element arrays across devices. Its POSITION is the Fortran's: after
        # both element kernels, before faulting, which is what lets faulting
        # and the mass divide be plain local work on a force array that is
        # already complete.
        force = B.nodal_sync(xp, inv, carry[FORCE])
        carry = carry[:FORCE] + (force,) + carry[FORCE + 1:]
        return part_b(carry, nt)

    return step


def run(S, nsteps=None, verbose=True, xp=np):
    """The whole solve. `xp` selects the backend: numpy or jax.numpy.

    Returns the same dict the six modules this replaces returned, so
    eqdyna3d.run_case and library_output's frt writer need no change.
    """
    nsteps = nsteps or S['nstep']

    inv = KU.build(S)
    finv = FLT.build(S)
    # Thermal-pressurization constants and history. The history is allocated at
    # the FULL step count for friclaw 5 (the Fortran allocates onFaultTPHist
    # the same way) and at width ZERO otherwise, so the carry has ONE shape for
    # every friction law and the jax loop does not need a second signature.
    tp = TP.build(S, nsteps) if S['friclaw'] == 5 else None
    hist_w = nsteps if S['friclaw'] == 5 else 0

    # Forced-rupture time is pure geometry -- computed once, not per step.
    # None means "swtwNucleation does nothing for this case", which is a
    # different statement from "tr is 1e9 everywhere" and is kept distinct so
    # a case that should nucleate and does not cannot look like a no-op.
    finv['tr'] = (FLT.forced_rupture_time(np, finv)
                  if FLT.nucleation_enabled(finv) else None)

    mass = np.concatenate(([1.0], S['nodalMassArr']))   # index 0 = unused sink
    # driver.f90:30 divides unconditionally. A zero lumped mass would give
    # inf/nan here; the ports this replaces silently masked it with
    # `where(mass>0, 1/mass, 0)`, which turns a broken mesh into a plausible
    # answer. Checked loudly instead.
    bad = int(np.count_nonzero(mass[1:] <= 0.0))
    if bad:
        raise ValueError(
            'driver.run: %d of %d lumped nodal masses are <= 0. driver.f90:30 '
            'divides by them unconditionally, so this mesh cannot be stepped. '
            'Fix the mass assembly rather than masking the divide.'
            % (bad, mass.shape[0] - 1))

    B.check_index_width(inv)
    inv = B.to_device(xp, inv)
    finv = B.to_device(xp, finv)
    if tp is not None:
        tp = B.to_device(xp, tp)
    mass = xp.asarray(mass)

    N = S['N']; NEQ = S['NEQ']; nftnd = S['nftnd']
    z = (lambda *a: xp.zeros(*a))
    carry0 = (z(NEQ + 1), z((N, 3)), z((N, 3)), z(NEQ + 1),
              xp.asarray(inv['stress_i0']).copy(), z((inv['Ep'], 15)),
              xp.asarray(S['fric_init'].copy()),
              xp.full(nftnd, gv.FNFT_SENTINEL),
              xp.asarray(0.0),
              z((nftnd, hist_w)), z((nftnd, hist_w)))

    # Which carry entries live on the ELEMENT axis, and therefore get cut
    # across devices under explicit decomposition. Stated here, beside
    # carry0, because this is where the shapes are; backend must not infer it
    # from a shape (Ei == Ep is possible on a small mesh).
    carry_shard = (None, None, None, None, 'Ei', 'Ep',
                   None, None, None, None, None)

    ndev = B.jax_device_count()
    if ndev > 1 and not B.is_jax(xp):
        raise RuntimeError(
            'EQDYNA_JAX_DEVICES=%d requests explicit domain decomposition, which '
            'exists only on the jax backend, but this run is numpy. Refusing '
            'rather than running serial under a %d-device label.' % (ndev, ndev))

    if ndev > 1 and B.timing_only():
        # Loud, on stderr, EVERY run -- one of the measurement knobs is set
        # and the answer this run produces is not the physics. eqdyna3d.run_case
        # additionally refuses to write it to the normal frt path.
        print('driver.run: *** TIMING-ONLY RUN, RESULT IS NOT VALID PHYSICS *** '
              '(%s=%s, %s=%s)' % (B.MODE_ENV, B.shard_mode(),
                                  B.SYNC_ENV, B.shard_sync()), file=sys.stderr)
    if verbose:
        print('driver.run: %d steps, backend=%s, friclaw=%d, devices=%d, mode=%s/%s'
              % (nsteps, xp.__name__, S['friclaw'], ndev,
                 B.shard_mode(), B.shard_sync()))
    t0 = time.perf_counter()
    scratch = KU.alloc_scratch(xp, inv)
    mk = lambda i: make_step(xp, i, finv, tp, mass, scratch)   # noqa: E731
    if ndev > 1:
        carry = B.run_time_loop_sharded(xp, mk, inv, carry0, nsteps, ndev,
                                        carry_shard)
    else:
        carry = B.run_time_loop(xp, mk, inv, carry0, nsteps)
    elapsed = time.perf_counter() - t0
    if verbose:
        print('driver.run: %.3f s, %.3f ms/step' % (elapsed, elapsed / nsteps * 1e3))

    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist) = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr),
                fnft=np.asarray(fnft), fric=np.asarray(fric),
                force=np.asarray(force))


def run_mpi(S, comm, nsteps=None, verbose=True, xp=np):
    """The whole solve, ONE PROCESS PER RANK -- Fortran's decomposition, with
    jax owning only the local element kernel.

    Structure, and why it is not `run` with a flag: `run` hands the entire
    time loop to XLA as one jitted fori_loop, which is exactly what makes the
    serial jax kernel beat Fortran's (611 vs 931 ms/step). An MPI call cannot
    happen inside that. So the step is jitted in TWO halves either side of
    driver.f90:27, and MPI4NodalQuant runs between them on the host:

        part_a (jit)   velDispUpdate, zero force, both element kernels
                       + gather this rank's halo equations, on device
        MPI            Sendrecv the halo values with each neighbour
        part_b (jit)   add the received partials, faulting, mass divide

    The two jits are built ONCE, before the loop -- a jax.jit constructed
    inside the loop recompiles every call. The halo gather is an extra OUTPUT
    of part_a and the received delta an extra ARGUMENT of part_b, so the only
    host traffic per step is the halo itself (O(boundary)), and the 30.5 MB
    force array is never copied to the host or re-scattered eagerly.

    Every rank builds the full serial mesh and then restricts it
    (MPI4NodalQuant.decompose), so global node and equation numbering is
    shared by all ranks and every index array comes from the gated serial
    path. See that module for the trade and for what is NOT bit-identical.

    Returns the same dict `run` returns, plus the fault-row selectors the
    caller needs to write this rank's `frt.txt<rank>`, plus `report` -- the
    per-rank counts any multi-rank measurement must print to be checkable."""
    import jax
    from . import MPI4NodalQuant as MQ

    rank = comm.Get_rank(); nranks = comm.Get_size()
    if not B.is_jax(xp):
        raise RuntimeError('driver.run_mpi: the MPI path exists for the jax '
                           'backend only (numpy is explicitly out of scope); '
                           'got xp=%s' % xp.__name__)
    nsteps = nsteps or S['nstep']

    inv = KU.build(S)
    finv = FLT.build(S)
    tp = TP.build(S, nsteps) if S['friclaw'] == 5 else None
    hist_w = nsteps if S['friclaw'] == 5 else 0
    finv['tr'] = (FLT.forced_rupture_time(np, finv)
                  if FLT.nucleation_enabled(finv) else None)

    loc = MQ.decompose(S, inv, finv, rank, nranks)
    inv_l, finv_l = loc['inv'], loc['finv']
    computed = loc['fault_computed_rows']
    if tp is not None:
        tp = MQ.restrict_rows(tp, computed, int(finv['nftnd']))

    # Every fault node must be written by exactly one rank, or the frt files
    # the gate reads are short and nothing says so.
    owned_total = comm.allreduce(int(loc['fault_rows'].shape[0]))
    if owned_total != int(finv['nftnd']):
        raise RuntimeError(
            'driver.run_mpi: the ranks together own %d of %d fault nodes. '
            'frt.txt* would be missing %d rows and the canonical comparison '
            'would silently compare a shorter file.'
            % (owned_total, int(finv['nftnd']), int(finv['nftnd']) - owned_total))

    mass = np.concatenate(([1.0], S['nodalMassArr']))
    bad = int(np.count_nonzero(mass[1:] <= 0.0))
    if bad:
        raise ValueError(
            'driver.run_mpi: %d of %d lumped nodal masses are <= 0. '
            'driver.f90:30 divides by them unconditionally.'
            % (bad, mass.shape[0] - 1))

    B.check_index_width(inv_l)
    inv_l = B.to_device(xp, inv_l)
    finv_l = B.to_device(xp, finv_l)
    if tp is not None:
        tp = B.to_device(xp, tp)
    mass = xp.asarray(mass)

    N = S['N']; NEQ = S['NEQ']; nftnd_l = int(finv_l['nftnd'])
    z = (lambda *a: xp.zeros(*a))
    carry = (z(NEQ + 1), z((N, 3)), z((N, 3)), z(NEQ + 1),
             xp.asarray(inv_l['stress_i0']).copy(), z((inv_l['Ep'], 15)),
             xp.asarray(S['fric_init'][computed].copy()),
             xp.full(nftnd_l, gv.FNFT_SENTINEL),
             xp.asarray(0.0),
             z((nftnd_l, hist_w)), z((nftnd_l, hist_w)))

    B.enable_compilation_cache()
    scratch = KU.alloc_scratch(xp, inv_l)
    dyn, sta = B.promote(xp, inv_l)
    halo = xp.asarray(loc['halo_idx'])
    nbrs = loc['neighbours']

    def a_body(dyn_arrays, c, nt, h):
        part_a, _ = make_step_parts(xp, {**sta, **dyn_arrays}, finv_l, tp,
                                    mass, scratch)
        c = part_a(c, nt)
        return c, c[FORCE][h]

    def b_body(dyn_arrays, c, nt, h, delta):
        _, part_b = make_step_parts(xp, {**sta, **dyn_arrays}, finv_l, tp,
                                    mass, scratch)
        force = B.addat(xp, c[FORCE], h, delta)     # MPI4NodalQuant's sum
        c = c[:FORCE] + (force,) + c[FORCE + 1:]
        return part_b(c, nt)

    # DONATE THE CARRY. Without donation each of the two calls allocates a
    # fresh buffer for every one of the 11 carry entries -- ~120 MB of copy
    # per step on test.tpv104 -- because XLA may not write into an input it
    # does not own. The fused fori_loop of the serial path updates its carry
    # in place and pays none of that, which is why breaking the loop open for
    # MPI costs 3x per step until the carry is donated. Safe here: `carry` is
    # rebound from the return value immediately and the previous value is
    # never read again.
    a_jit = jax.jit(a_body, donate_argnums=(1,))
    b_jit = jax.jit(b_body, donate_argnums=(1,))
    sync = MQ.sync_mode()
    # Under 'allreduce' the force is already complete when part_b runs, so
    # part_b's halo add is a no-op -- kept (rather than branching part_b) so
    # BOTH sync modes execute the identical jitted code and a difference
    # between their timings is the exchange and nothing else.
    zero_delta = xp.zeros(halo.shape[0])

    rep = loc['report']
    if verbose:
        print('driver.run_mpi rank %d/%d: %d steps, friclaw=%d, elements '
              'Ei=%d Ep=%d E=%d (work %.1f%% of total), nodes=%d eqs=%d '
              'halo=%d (%.2f%% of eqs), neighbours=%s, fault computed=%d '
              'owned=%d'
              % (rank, nranks, nsteps, S['friclaw'], rep['Ei'], rep['Ep'],
                 rep['E'], 100.0 * rep['work'] / rep['work_total'],
                 rep['nodes'], rep['eqs'], rep['halo_eqs'],
                 100.0 * rep['halo_frac'], rep['neighbours'],
                 rep['fault_computed'], rep['fault_owned']), flush=True)

    comm.Barrier()
    c0 = os.times()
    t0 = time.perf_counter()
    t_mpi = t_wait = 0.0
    for nt in range(1, nsteps + 1):
        carry, hv = a_jit(dyn, carry, nt, halo)
        # BLOCK BEFORE STARTING THE MPI CLOCK. jax dispatch is asynchronous,
        # so a_jit returns before part_a has run and the first thing that
        # touches hv absorbs the whole element kernel. Timing the exchange
        # without this reported 298 of 471 ms/step "in MPI4NodalQuant" on a
        # ONE-rank run with zero neighbours -- i.e. it was measuring the
        # solver, not the exchange.
        jax.block_until_ready(hv)
        # Barrier FIRST, timed separately. Without it the fastest rank's
        # "exchange" time is mostly waiting for the slowest rank, and on this
        # box the per-rank spread is real (cpu 61 measured EFFECTIVE_CORES
        # 0.50 against cpu 60's 0.97 on identical work): rank 0 reported 296
        # of 622 ms/step "in MPI4NodalQuant" while rank 1 reported 0.5 ms for
        # the same exchange. Charging that to the collective would be a
        # measurement error in the exact shape this campaign is trying to
        # avoid, so load imbalance is t_wait and the exchange is t_mpi.
        t_w = time.perf_counter()
        comm.Barrier()
        t_wait += time.perf_counter() - t_w
        t1 = time.perf_counter()
        if sync == 'allreduce':
            # The simple version: no ownership bookkeeping at all. Assembly
            # is a SUM, so every rank can hold the full-length partial array
            # and ONE Allreduce makes every rank's copy the complete one.
            # Moves O(NEQ) where the halo moves O(boundary); both are measured
            # side by side in testsys/perf/run_mpi_scaling.py.
            total = MQ.allreduce(comm, np.asarray(jax.device_get(carry[FORCE])))
            t_mpi += time.perf_counter() - t1
            carry = carry[:FORCE] + (xp.asarray(total),) + carry[FORCE + 1:]
            carry = b_jit(dyn, carry, nt, halo, zero_delta)
        else:
            delta = MQ.exchange(comm, nbrs, np.asarray(jax.device_get(hv)))
            t_mpi += time.perf_counter() - t1
            carry = b_jit(dyn, carry, nt, halo, xp.asarray(delta))
    jax.block_until_ready(carry)
    elapsed = time.perf_counter() - t0
    c1 = os.times()
    # EFFECTIVE_CORES: cpu seconds this process consumed per wall second. A
    # rank pinned to one cpu should read ~1.0; below that it was starved (it
    # shared its cpu with a foreign tenant, or it was waiting on the halo),
    # and that is invisible in a wall-clock number alone.
    eff = ((c1[0] - c0[0]) + (c1[1] - c0[1])) / elapsed if elapsed > 0 else 0.0

    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist) = carry
    rep = dict(rep, ms_per_step=elapsed / nsteps * 1e3, solve_s=elapsed,
               mpi_ms_per_step=t_mpi / nsteps * 1e3,
               wait_ms_per_step=t_wait / nsteps * 1e3,
               sync=sync, nsteps=nsteps, effective_cores=eff,
               device_peak_gb=_device_peak_gb(jax),
               threads=len(os.listdir('/proc/%d/task' % os.getpid())),
               cpus=sorted(os.sched_getaffinity(0)),
               cpus_allowed=len(os.sched_getaffinity(0)))
    if verbose:
        print('driver.run_mpi rank %d: %.3f s, %.3f ms/step (sync=%s: %.3f '
              'ms/step exchanging, %.3f ms/step waiting at the barrier), '
              'EFFECTIVE_CORES %.2f, threads=%d, cpus_allowed=%d %s'
              % (rank, elapsed, rep['ms_per_step'], sync,
                 rep['mpi_ms_per_step'], rep['wait_ms_per_step'],
                 eff, rep['threads'], rep['cpus_allowed'], rep['cpus']),
              flush=True)
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr),
                fnft=np.asarray(fnft), fric=np.asarray(fric),
                force=np.asarray(force),
                fault_rows=loc['fault_rows'],
                own_in_computed=loc['own_in_computed'], report=rep)
