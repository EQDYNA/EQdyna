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
from . import profile_emit as _profile_emit
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


def make_step_parts(xp, inv, finv, tp, mass, scratch, fault_timer=None):
    """The step, split at driver.f90:27 -- MPI4NodalQuant's position.

    part_a: timeElapsed, velDispUpdate, zero the force, both element kernels.
    part_b: thermal pressurization, faulting, the mass divide.

    ONE body, split rather than copied, because the two callers need the seam
    in a different place in the STACK, not in the code: make_step closes the
    seam with backend.nodal_sync (identity when serial, a device-mesh
    collective under shard_map) and keeps ONE jitted time loop, while run_mpi
    must leave the jit at the seam to make an MPI call and therefore jits the
    two halves separately. Both perform the same operations, in the same
    order, on the same operands.

    `fault_timer`, if given, is a mutable `{'s': 0.0}` this accumulates
    `FLT.faulting`'s wall time into, split out along Fortran's own boundary
    (`fault` = compTimeInSeconds(6), faulting.f90:28 -- rule 23). Callers must
    pass `None` (the default) unless `not B.is_jax(xp)`: numpy executes this
    call eagerly and synchronously, so timing it costs no new sync and moves
    no bit; under jax tracing a `perf_counter()` pair here would measure the
    ONE-TIME trace, not the per-step cost, which is wrong, not just imprecise
    -- see profile_emit.py's docstring for why `fault` stays folded into
    `element` on that backend."""
    dt = inv['dt']; rdampk = inv['rdampk']
    tr = finv['tr']
    friclaw = finv['friclaw']

    st_off_idx = finv['st_off_idx']
    n_off_st = st_off_idx.shape[0]

    def part_a(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist, on_st_hist, off_st_hist) = carry
        timeElapsed = timeElapsed + dt                       # driver.f90:12

        v1, velArr, dispArr = velDispUpdate(xp, inv, v1, velArr, dispArr,
                                            force, dt)

        # driver.f90:21 storeOffFaultStData -- recorded right after
        # velDispUpdate, on THIS step's just-updated dispArr/velArr, exactly
        # where the Fortran call sits (before the force array is zeroed,
        # which off-fault stations never read anyway). idhist's dispOrVel
        # is always 1 or 2 (never 3, see eqdyna3d.f90:287-298's
        # allocInitAfterMeshGen loop) -- nodalForceArr/eqNum is dead code in
        # storeOffFaultStData and is not reproduced here.
        # Columns, from output_offfault_st's write order (library_output.f90
        # :235-242): t, x-disp, x-vel, -z-disp, -z-vel, y-disp, y-vel.
        if n_off_st:
            col = nt - 1
            xh = st_off_idx
            vals = xp.stack([
                xp.full((n_off_st,), timeElapsed),
                dispArr[xh, 0], velArr[xh, 0],
                -dispArr[xh, 2], -velArr[xh, 2],
                dispArr[xh, 1], velArr[xh, 1],
            ], axis=1)
            off_st_hist = B.setat(xp, off_st_hist, (slice(None), slice(None), col), vals)

        force = B.setat(xp, force, slice(None), 0.0)         # driver.f90:23
        force, stress_i, s_p = KU.assembleGlobalKU(
            xp, inv, velArr, force, stress_i, s_p, dt, rdampk, scratch)
        force = KU.calcHourglassResist(xp, inv, dispArr, velArr, force, rdampk,
                                       scratch)

        # part_a ENDS HERE, at driver.f90:27 -- the seam each caller closes
        # with its own MPI4NodalQuant (make_step's backend.nodal_sync,
        # run_mpi's real MPI exchange). Nothing is exchanged inside part_a;
        # see make_step for what happens at the seam and why it sits here.
        return (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft,
                timeElapsed, sliprate_hist, shear_hist, on_st_hist, off_st_hist)

    st_on_idx = finv['st_on_idx']
    n_on_st = st_on_idx.shape[0]
    st_ncols_on = finv['st_ncols_on']
    st_sign = finv['st_sign']

    def part_b(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist, on_st_hist, off_st_hist) = carry

        if friclaw == 5:                                   # driver.f90:28
            fric = TP.updateThermalPressurization(
                xp, tp, fric, sliprate_hist, shear_hist, nt, dt)

        if fault_timer is not None:
            _t0 = time.perf_counter()
            fric, fnft, force = FLT.faulting(xp, finv, fric, fnft, velArr,
                                             dispArr, force, dt, timeElapsed,
                                             tr, nt)
            fault_timer['s'] += time.perf_counter() - _t0
        else:
            fric, fnft, force = FLT.faulting(xp, finv, fric, fnft, velArr,
                                             dispArr, force, dt, timeElapsed,
                                             tr, nt)

        # faulting.f90:24 storeOnFaultStationQuantSCEC -- recorded right
        # after faulting (so fric holds THIS step's post-friction-law
        # traction/state), before the mass divide below (which never touches
        # fric). See library_output.write_onfault_stations for the exact
        # column derivation this mirrors (getNsdSlipSliprateTraction writes
        # fric[SLIP_STRIKE]/[SLIP_DIP]/[SLIPRATE_STRIKE]/[SLIPRATE_DIP]
        # BEFORE solveRSF's own background-slip-rate add at faulting.f90:219
        # -- so for friclaw>=3 the creep term is added back here, at record
        # time, exactly as solveRSF's local nsdSlipVector/nsdSliprateVector
        # would carry it into storeOnFaultStationQuantSCEC; friclaw<=2 never
        # adds it, matching solveSWTW, which never touches those slots).
        if n_on_st:
            col = nt - 1
            xh = st_on_idx
            slipS = fric[xh, gv.SLIP_STRIKE]
            slipD = fric[xh, gv.SLIP_DIP]
            srS = fric[xh, gv.SLIPRATE_STRIKE]
            srD = fric[xh, gv.SLIPRATE_DIP]
            if friclaw >= 3:
                vx = fric[xh, gv.VINI_X]; vz = fric[xh, gv.VINI_Z]
                slipS = slipS + vx * timeElapsed
                slipD = slipD + vz * timeElapsed
                srS = srS + vx
                srD = srD + vz
            hShear = fric[xh, gv.TRACT_STRIKE] / 1.0e6
            vShear = -fric[xh, gv.TRACT_DIP] / 1.0e6
            nStress = st_sign * fric[xh, gv.TRACT_NORM] / 1.0e6
            cols = [xp.full((n_on_st,), timeElapsed), slipS, srS, hShear,
                    -slipD, -srD, vShear, nStress]
            if st_ncols_on == 11:
                cols += [fric[xh, gv.STATE], fric[xh, gv.TP_TEMP],
                         (fric[xh, gv.TP_NORM_TP] + fric[xh, gv.TP_PINI]) / 1.0e6]
            vals = xp.stack(cols, axis=1)
            on_st_hist = B.setat(xp, on_st_hist, (slice(None), slice(None), col), vals)

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
                timeElapsed, sliprate_hist, shear_hist, on_st_hist, off_st_hist)

    return part_a, part_b


def build_invariants(S, nsteps):
    """The loop-invariant state both entry points build, identically.

    Returns (inv, finv, tp, hist_w). ONE copy, because this is where a fix
    like forced_rupture_time's would otherwise have to be made twice -- the
    exact shape of the port.py/port_jax.py divergence faulting.py's docstring
    records. `tp` is None and the history width 0 unless friclaw==5, so the
    carry has one shape for every friction law (see run's carry0).
    """
    inv = KU.build(S)
    finv = FLT.build(S)
    tp = TP.build(S, nsteps) if S['friclaw'] == 5 else None
    hist_w = nsteps if S['friclaw'] == 5 else 0
    # Forced-rupture time is pure geometry -- computed once, not per step.
    # None means "swtwNucleation does nothing for this case", which is a
    # different statement from "tr is 1e9 everywhere" and is kept distinct so
    # a case that should nucleate and does not cannot look like a no-op.
    finv['tr'] = (FLT.forced_rupture_time(np, finv)
                  if FLT.nucleation_enabled(finv) else None)

    # Row 114 -- station output (library_output.f90's output_onfault_st /
    # output_offfault_st). Index arrays are host numpy here (as every other
    # finv array is, at this point) and travel through B.to_device with the
    # rest of finv exactly like nsmp1/idxF_s/etc -- no new plumbing path.
    # st_on_idx/st_off_idx are 0-indexed rows into the fault-node axis
    # (matches S['nsmp1']/fric) and the node axis (matches
    # velArr/dispArr) respectively; empty (shape (0,)) on the python-jax-mpi
    # path (eqdyna3d.build_solver_state does not run the matching there --
    # see run_case_mpi's loud refusal). st_ncols_on/st_sign are plain Python
    # scalars (static under jit, like `friclaw`), not device arrays.
    finv['st_on_idx'] = np.asarray(S['st_on_idx'],
                                   dtype=np.int64)
    finv['st_off_idx'] = np.asarray(S['st_off_idx'],
                                    dtype=np.int64)
    finv['st_ncols_on'] = 11 if S['friclaw'] >= 3 else 8
    # Column 8's sign is the case's spec convention (board row 22a): no
    # default -- a state without it is a caller bug, not a sign to guess.
    finv['st_sign'] = float(S['nStressOutSign'])
    return inv, finv, tp, hist_w


def make_step(xp, inv, finv, tp, mass, scratch, fault_timer=None):
    """Build the per-step closure. Called by backend.run_time_loop INSIDE
    the jit on the jax path, so that `inv`'s arrays resolve to jit arguments
    rather than closed-over HLO literals.

    `fault_timer`: see make_step_parts. Threaded through unchanged; this
    function adds no timing of its own."""
    part_a, part_b = make_step_parts(xp, inv, finv, tp, mass, scratch,
                                     fault_timer=fault_timer)

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

    Returns the same dict the six modules this replaces returned, plus one
    new key, `fault_s` -- the wall time spent inside `FLT.faulting` across
    all steps, measured only `not B.is_jax(xp)` (see make_step_parts) and
    0.0 on jax. eqdyna3d.run_case reads it to split `fault` out of `element`
    for python-numpy; library_output's frt writer needs no change.
    """
    nsteps = nsteps or S['nstep']

    # Thermal-pressurization history is allocated at the FULL step count for
    # friclaw 5 (the Fortran allocates onFaultTPHist the same way) and at
    # width ZERO otherwise -- see build_invariants.
    inv, finv, tp, hist_w = build_invariants(S, nsteps)

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

    # Row 114 station-history widths, read off finv BEFORE B.to_device (still
    # plain host numpy/int here, same as every other finv array at this point).
    n_on_st = int(finv['st_on_idx'].shape[0])
    n_off_st = int(finv['st_off_idx'].shape[0])
    st_ncols_on = finv['st_ncols_on']

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
              z((nftnd, hist_w)), z((nftnd, hist_w)),
              z((n_on_st, st_ncols_on, nsteps)), z((n_off_st, 7, nsteps)))

    # Which carry entries live on the ELEMENT axis, and therefore get cut
    # across devices under explicit decomposition. Stated here, beside
    # carry0, because this is where the shapes are; backend must not infer it
    # from a shape (Ei == Ep is possible on a small mesh).
    carry_shard = (None, None, None, None, 'Ei', 'Ep',
                   None, None, None, None, None, None, None)

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
    # EQDYNA_PROFILE read ONCE here, before the loop is built -- never per
    # step, never via getenv inside a traced function. `fault_timer` lives
    # HERE, outside `mk`, so it survives regardless of how many times `mk`
    # is invoked underneath `run_time_loop`/`run_time_loop_sharded` -- every
    # `make_step` call closes over this SAME dict. None on jax (unaffected
    # by the flag: see make_step_parts's docstring for why a timer must not
    # exist on that path, not just go unread) AND None when
    # EQDYNA_PROFILE=0 -- the per-step perf_counter() pair inside
    # make_step_parts's part_b is a profiler addition and must not run when
    # the switch is off.
    profile_on = _profile_emit.enabled()
    fault_timer = {'s': 0.0} if (profile_on and not B.is_jax(xp)) else None
    mk = lambda i: make_step(xp, i, finv, tp, mass, scratch,   # noqa: E731
                             fault_timer=fault_timer)
    if ndev > 1:
        carry = B.run_time_loop_sharded(xp, mk, inv, carry0, nsteps, ndev,
                                        carry_shard)
    else:
        carry = B.run_time_loop(xp, mk, inv, carry0, nsteps)
    elapsed = time.perf_counter() - t0
    if verbose:
        print('driver.run: %.3f s, %.3f ms/step' % (elapsed, elapsed / nsteps * 1e3))
    fault_s = fault_timer['s'] if fault_timer is not None else 0.0

    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist, on_st_hist, off_st_hist) = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr),
                fnft=np.asarray(fnft), fric=np.asarray(fric),
                force=np.asarray(force), fault_s=fault_s,
                on_st_hist=np.asarray(on_st_hist), off_st_hist=np.asarray(off_st_hist))


def run_mpi(S, comm, part, plan, nsteps=None, verbose=True, xp=np):
    """The whole solve, ONE PROCESS PER RANK -- Fortran's decomposition, with
    jax owning only the local element kernel.

    Structure, and why it is not `run` with a flag: `run` hands the entire
    time loop to XLA as one jitted fori_loop, which is exactly what makes the
    serial jax kernel beat Fortran's (611 vs 931 ms/step). An MPI call cannot
    happen inside that. So the step is jitted in TWO halves either side of
    driver.f90:27, and MPI4NodalQuant runs between them on the host:

        part_a (jit)   velDispUpdate, zero force, both element kernels
                       + gather this rank's halo equations, on device
        MPI            MPI4NodalQuant's x/y/z relay over the halo values
        part_b (jit)   write back the completed halo values, faulting,
                       mass divide

    The two jits are built ONCE, before the loop -- a jax.jit constructed
    inside the loop recompiles every call. The halo gather is an extra OUTPUT
    of part_a and the completed values an extra ARGUMENT of part_b, so the
    only host traffic per step is the halo itself (O(boundary)).

    `S` is THIS RANK'S BOX, built rank-locally (eqdyna3d.build_solver_state
    with a Partition) and completed by MPI4NodalQuant.setup_exchange, which
    also produced `plan` (the exchange faces, the halo equations, the fault
    rows this rank writes). Every array here is already at local extent;
    there is no global mesh, no restriction and no relabelling left to do.

    Returns the same dict `run` returns, plus the fault-row selectors the
    caller needs to write this rank's `frt.txt<rank>`, plus `report` -- the
    per-rank counts any multi-rank measurement must print to be checkable."""
    import jax
    from . import MPI4NodalQuant as MQ

    # PRE-LOOP setup timer (invariants, device transfer, jit function
    # CONSTRUCTION -- not execution, so no sync). Added by the profile-emitter
    # landing (2026-09-23): without it, eqdyna3d.run_case_mpi's outer
    # 'solve' Profile phase (which wraps this ENTIRE call) attributed only
    # the step-loop portion to any bucket, and this pre-loop/jit-build cost
    # -- MEASURED on test.tpv8 x 4 ranks: 2.4-3.8 s of a 12.0 s total_s, i.e.
    # roughly 30%, not a rounding error -- fell into unaccounted_s, which
    # blew past profile_schema's 5% SUM_TOLERANCE. This is real pre-loop
    # work (Fortran's analogue is meshgen+assembleGlobalMass, its own
    # `setup` bucket), so it belongs in `setup`, not in a gap. No new sync:
    # build_invariants/to_device/jax.jit(...) construction are synchronous host-side
    # calls already executing on this path; this only wraps them with two
    # perf_counter() calls.
    #
    # EQDYNA_PROFILE read ONCE here, before setup or the loop -- never per
    # step, never via getenv again below. When off, t_setup0 stays 0.0 and
    # every timer this landing added downstream (t_setup, the per-step
    # compute-dispatch timing, the tail-drain timer) is skipped so OFF is
    # exactly the pre-profile-emitter per-step code path.
    profile_on = _profile_emit.enabled()
    t_setup0 = time.perf_counter() if profile_on else 0.0

    rank = comm.Get_rank(); nranks = comm.Get_size()
    if (rank, nranks) != (part.rank, part.nranks):
        raise RuntimeError('driver.run_mpi: communicator is rank %d/%d but the '
                           'mesh was built for rank %d/%d'
                           % (rank, nranks, part.rank, part.nranks))
    if not B.is_jax(xp):
        raise RuntimeError('driver.run_mpi: the MPI path exists for the jax '
                           'backend only (numpy is explicitly out of scope); '
                           'got xp=%s' % xp.__name__)
    nsteps = nsteps or S['nstep']
    sync = MQ.sync_mode()
    if sync == 'allreduce':
        # NO FALLBACK to halo. See MQ.sync_mode: the allreduce sync reduced
        # the full nodal array and therefore needed it replicated at global
        # extent on every rank, which no rank holds any more.
        raise ValueError(
            'driver.run_mpi: %s=allreduce reduced the FULL nodal array and '
            'required every rank to hold it at global extent. Each rank now '
            'holds only its own subdomain, so the ranks\' force arrays have '
            'different lengths and there is nothing an Allreduce could '
            'reduce. Use %s=halo.' % (MQ.SYNC_ENV, MQ.SYNC_ENV))

    inv, finv, tp, hist_w = build_invariants(S, nsteps)
    nftnd_l = int(finv['nftnd'])

    # driver.f90:30 divides by every lumped mass unconditionally. Each rank
    # checks its own (complete, post-relay) masses and the COUNT is summed, so
    # a bad equation anywhere stops every rank -- every equation is held by at
    # least one rank, and no rank holds a global array to check instead.
    mass_h = np.concatenate(([1.0], S['nodalMassArr']))
    bad = comm.allreduce(int(np.count_nonzero(mass_h[1:] <= 0.0)))
    if bad:
        raise ValueError(
            'driver.run_mpi: %d lumped nodal mass(es) <= 0 across the ranks. '
            'driver.f90:30 divides by them unconditionally.' % bad)

    # Row 114: python-jax-mpi does not port station output (build_solver_state
    # hands this path empty st_on_idx/st_off_idx -- see eqdyna3d.run_case_mpi's
    # loud warning when a case actually names stations). n_on_st_l/n_off_st_l
    # are therefore always 0 here, but the carry tuple still carries the two
    # slots so make_step_parts's part_a/part_b (ONE implementation, shared with
    # `run`) can unpack the same-length tuple on both entry points.
    n_on_st_l = int(finv['st_on_idx'].shape[0])
    n_off_st_l = int(finv['st_off_idx'].shape[0])
    st_ncols_on_l = finv['st_ncols_on']

    B.check_index_width(inv)
    fric_init = S['fric_init'].copy()
    inv = B.to_device(xp, inv)
    finv = B.to_device(xp, finv)
    if tp is not None:
        tp = B.to_device(xp, tp)
    mass = xp.asarray(mass_h)

    N_l = int(inv['N']); NEQ_l = int(inv['NEQ'])
    z = (lambda *a: xp.zeros(*a))
    carry = (z(NEQ_l + 1), z((N_l, 3)), z((N_l, 3)), z(NEQ_l + 1),
             xp.asarray(inv['stress_i0']).copy(), z((inv['Ep'], 15)),
             xp.asarray(fric_init),
             xp.full(nftnd_l, gv.FNFT_SENTINEL),
             xp.asarray(0.0),
             z((nftnd_l, hist_w)), z((nftnd_l, hist_w)),
             z((n_on_st_l, st_ncols_on_l, nsteps)), z((n_off_st_l, 7, nsteps)))
    # Per-rank carry bytes, from the ALLOCATED arrays rather than recomputed
    # from shapes: the primary check that the subdomain really shrank.
    carry_bytes = {n: int(a.nbytes) for n, a in zip(
        ('v1', 'velArr', 'dispArr', 'force', 'stress_i', 's_p', 'fric',
         'fnft', 'timeElapsed', 'sliprate_hist', 'shear_hist',
         'on_st_hist', 'off_st_hist'), carry)}
    carry_bytes['mass'] = int(mass.nbytes)

    # PER-RANK CACHE DIRECTORY. Every rank traces part_a/part_b against its
    # own Ei/Ep/halo shapes, so no two ranks can share a cache entry -- a
    # shared directory buys no reuse and pays JAX's per-key lock, which is
    # what wedged this cell three times (see backend.enable_compilation_cache).
    B.enable_compilation_cache(subdir='rank%d' % rank)
    scratch = KU.alloc_scratch(xp, inv)
    dyn, sta = B.promote(xp, inv)
    halo = xp.asarray(plan['halo_idx'])
    faces = plan['faces']

    def a_body(dyn_arrays, c, nt, h):
        part_a, _ = make_step_parts(xp, {**sta, **dyn_arrays}, finv, tp,
                                    mass, scratch)
        c = part_a(c, nt)
        return c, c[FORCE][h]

    def b_body(dyn_arrays, c, nt, h, vals):
        _, part_b = make_step_parts(xp, {**sta, **dyn_arrays}, finv, tp,
                                    mass, scratch)
        # MPI4NodalQuant's completed values, WRITTEN over the partials (see
        # MQ.exchange for why a set and not an add of a delta).
        force = B.setat(xp, c[FORCE], h, vals)
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

    rep = dict(rank=rank, nranks=nranks, decomp=part.dims, mexyz=part.mexyz,
               Ei=int(inv['Ei']), Ep=int(inv['Ep']), E=int(inv['E']),
               nodes=N_l, eqs=NEQ_l, halo_eqs=int(plan['halo_idx'].size),
               halo_frac=float(plan['halo_idx'].size) / max(NEQ_l, 1),
               neighbours=[int(f['nb']) for f in faces],
               fault_computed=nftnd_l,
               fault_owned=int(plan['fault_rows'].shape[0]),
               N_local=N_l, NEQ_local=NEQ_l)
    if verbose:
        print('driver.run_mpi rank %d/%d: %d steps, friclaw=%d, decomposition '
              '%r box %r, elements Ei=%d Ep=%d E=%d, LOCAL N=%d NEQ=%d, '
              'halo=%d (%.2f%% of eqs), faces with %s, fault computed=%d '
              'owned=%d, carry %.2f MB'
              % (rank, nranks, nsteps, S['friclaw'], part.dims, part.mexyz,
                 rep['Ei'], rep['Ep'], rep['E'], N_l, NEQ_l, rep['halo_eqs'],
                 100.0 * rep['halo_frac'], rep['neighbours'],
                 rep['fault_computed'], rep['fault_owned'],
                 sum(carry_bytes.values()) / 1e6), flush=True)

    # STEP ATTRIBUTION, off by default (MPI4NodalQuant.step_profile). The two
    # accumulators below cost two perf_counter calls per step when on and
    # nothing when off, and they do NOT add a block_until_ready -- t_compute
    # brackets exactly the dispatch-plus-block the loop already performs, so
    # the profiled run measures the same step the production run does.
    #
    # WHAT t_compute MEANS, exactly, because it is easy to over-read: jax
    # dispatch is asynchronous and nothing blocks on part_b, so the
    # block_until_ready(hv) below drains part_b of the PREVIOUS step as well as
    # part_a of this one. t_compute is therefore ALL jitted compute per step,
    # a and b together, not part_a alone. That is the quantity the plateau
    # question needs (compute vs not-compute); splitting a from b costs an
    # extra block per step and is what probe_mpi_step_split.py is for.
    prof = MQ.step_profile()
    comm.Barrier()
    # t_setup stops HERE, right after this pre-existing barrier (present on
    # the default path already, unconditional -- not added by this change):
    # build_invariants/to_device/jit-construction plus the rendezvous that lines
    # every rank up before the loop's own clock starts. Reported as `setup`
    # (see the comment above this function's `t_setup0`), not folded into
    # `wait_ms_per_step`/t_wait -- that field's existing, documented meaning
    # ("0.0 BY CONSTRUCTION on the production path", the comment below) is
    # about the PER-STEP barrier and stays exactly as it was.
    #
    # Skipped (stays 0.0) when EQDYNA_PROFILE=0 -- `profile_on` was read
    # once, before this loop, at the top of this function.
    t_setup = (time.perf_counter() - t_setup0) if profile_on else 0.0
    c0 = os.times()
    t0 = time.perf_counter()
    t_mpi = t_wait = t_compute = t_d2h = 0.0
    for nt in range(1, nsteps + 1):
        # t_compute is accumulated whenever EITHER the pre-existing
        # EQDYNA_MPI_STEP_PROFILE knob (`prof`) OR the profile-emitter's
        # EQDYNA_PROFILE switch (`profile_on`, default ON) is set. The
        # block below (jax.block_until_ready(hv)) already runs on every
        # step of the PRODUCTION path regardless -- it always did, for the
        # reason in the comment just below -- so timing around an
        # already-mandatory sync adds no new synchronisation point. But the
        # two perf_counter() calls themselves are new profiler-added cost,
        # and EQDYNA_PROFILE=0 (with `prof` also unset) must skip them, not
        # just skip the eventual file write. See profile_emit.py's module
        # docstring for what this number means (part_a AND the previous
        # step's part_b, i.e. element kernels + faulting fused, per
        # driver.run_mpi's own async-pipeline comment below) and why
        # `wait`/barrier timing is NOT extended the same way (it requires
        # comm.Barrier(), a real new collective the default path must not
        # pay).
        _time_compute = prof or profile_on
        ta0 = time.perf_counter() if _time_compute else 0.0
        carry, hv = a_jit(dyn, carry, nt, halo)
        # BLOCK BEFORE STARTING THE MPI CLOCK. jax dispatch is asynchronous,
        # so a_jit returns before part_a has run and the first thing that
        # touches hv absorbs the whole element kernel. Timing the exchange
        # without this reported 298 of 471 ms/step "in MPI4NodalQuant" on a
        # ONE-rank run with zero neighbours -- i.e. it was measuring the
        # solver, not the exchange.
        jax.block_until_ready(hv)
        if _time_compute:
            t_compute += time.perf_counter() - ta0
        # THE BARRIER IS A MEASUREMENT DEVICE AND RUNS ONLY UNDER THE PROFILE.
        #
        # What it is for, unchanged: without it the fastest rank's "exchange"
        # time is mostly waiting for the slowest rank, and on this box the
        # per-rank spread is real (cpu 61 measured EFFECTIVE_CORES 0.50
        # against cpu 60's 0.97 on identical work): rank 0 reported 296 of 622
        # ms/step "in MPI4NodalQuant" while rank 1 reported 0.5 ms for the
        # same exchange. Charging that to the collective would be a
        # measurement error, so under the profile load imbalance is t_wait and
        # the exchange is t_mpi.
        #
        # Why it must NOT run otherwise: it was unconditional, so every
        # production step paid a 32-rank global rendezvous that the algorithm
        # does not need. The exchange is nearest-neighbour (at most six face
        # neighbours) and deadlock-free on its own, so a rank that is
        # momentarily slow should delay its neighbours, not all 31 others. A global barrier per step makes the run pay the MAXIMUM
        # over 32 ranks of each step's jitter instead of letting jitter
        # average out along the chain. Measured at 32 ranks on test.tpv104,
        # differenced over 200 vs 40 steps: 13.15 ms/step mean in this
        # barrier, 30% of a 43.39 ms step.
        #
        # THE PROFILED STEP IS THEREFORE NOT THE PRODUCTION STEP, and the
        # report says so (`barrier_in_step`) so a profiled TOTAL cannot be
        # quoted as a production per-step cost.
        if prof:
            t_w = time.perf_counter()
            comm.Barrier()
            t_wait += time.perf_counter() - t_w
        t1 = time.perf_counter()
        if prof:
            # Split the halo device-to-host copy out of the exchange, so
            # t_mpi is PURE MPI under the profile. The default path leaves
            # them fused exactly as they were, because moving the clock
            # would redefine an already-published number.
            tg0 = time.perf_counter()
            hv_host = np.array(jax.device_get(hv))   # writable copy: exchange adds in place
            t_d2h += time.perf_counter() - tg0
            t1 = time.perf_counter()
            hvals = MQ.exchange(comm, part, faces, hv_host)
        else:
            hvals = MQ.exchange(comm, part, faces,
                                np.array(jax.device_get(hv)))
        t_mpi += time.perf_counter() - t1
        # Dispatch-side timing only (b_jit/xp.asarray do not block; jax
        # queues them and the actual device execution is drained by NEXT
        # iteration's `jax.block_until_ready(hv)` above, or by the final
        # block below on the last step) -- added by the profile-emitter
        # landing so the always-on `element` bucket also counts the host
        # dispatch/H2D-queue cost of part_b's input prep, which previously
        # sat entirely in unaccounted_s. No new sync: neither call here
        # blocks; this only wraps calls that were already being made. But
        # the perf_counter() pair is profiler-added cost in its own right,
        # so it is skipped -- not just left unwritten -- when `profile_on`
        # is False (EQDYNA_PROFILE=0).
        if profile_on:
            tb0 = time.perf_counter()
            carry = b_jit(dyn, carry, nt, halo, xp.asarray(hvals))
            t_compute += time.perf_counter() - tb0
        else:
            carry = b_jit(dyn, carry, nt, halo, xp.asarray(hvals))
    # Drains the LAST step's part_b (every earlier step's part_b was already
    # drained, one step later, by the loop's own block_until_ready(hv)).
    # This block is pre-existing and unconditional; only the timing around
    # it is new, so folding its cost into `element` (element+fault fused,
    # see profile_emit.py) adds no sync, just attributes an already-paid one.
    # Skipped the same way when EQDYNA_PROFILE=0: block_until_ready(carry)
    # still runs (it always did, unconditionally), only its two
    # perf_counter() calls are profiler-added and go away.
    if profile_on:
        tf0 = time.perf_counter()
        jax.block_until_ready(carry)
        t_compute += time.perf_counter() - tf0
    else:
        jax.block_until_ready(carry)
    elapsed = time.perf_counter() - t0
    c1 = os.times()
    # EFFECTIVE_CORES: cpu seconds this process consumed per wall second. A
    # rank pinned to one cpu should read ~1.0; below that it was starved (it
    # shared its cpu with a foreign tenant, or it was waiting on the halo),
    # and that is invisible in a wall-clock number alone.
    eff = ((c1[0] - c0[0]) + (c1[1] - c0[1])) / elapsed if elapsed > 0 else 0.0

    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist, on_st_hist, off_st_hist) = carry
    rep = dict(rep, ms_per_step=elapsed / nsteps * 1e3, solve_s=elapsed,
               mpi_ms_per_step=t_mpi / nsteps * 1e3,
               wait_ms_per_step=t_wait / nsteps * 1e3,
               # Whole-loop TOTALS (not ms/step), unconditional since this
               # landing: docs/run_profile.md's always-on profile.rank<r>.json
               # buckets (element=compute folded with fault, exchange=mpi,
               # wait=wait -- see profile_emit.py) read these directly rather
               # than re-deriving seconds from a ms/step average.
               compute_s=t_compute, mpi_s=t_mpi, wait_s=t_wait,
               setup_s=t_setup,
               sync=sync, nsteps=nsteps, effective_cores=eff,
               # False on the production path: the per-step global barrier is
               # a profiling device now (see the loop). wait_ms_per_step is
               # then 0.0 BY CONSTRUCTION -- there is nothing to wait at --
               # and imbalance shows up inside mpi_ms_per_step instead. A
               # 0.00 wait must not be read as "perfectly balanced".
               barrier_in_step=bool(prof),
               carry_bytes=carry_bytes,
               carry_bytes_total=sum(carry_bytes.values()),
               device_peak_gb=_device_peak_gb(jax),
               threads=len(os.listdir('/proc/%d/task' % os.getpid())),
               cpus=sorted(os.sched_getaffinity(0)),
               cpus_allowed=len(os.sched_getaffinity(0)))
    if prof:
        # The residual is elapsed MINUS everything attributed, so the five
        # buckets sum to the step by construction and an unattributed cost
        # cannot hide in a gap between two clocks. It covers the Python loop
        # itself, pytree flatten/unflatten on both dispatches, and
        # xp.asarray(hvals).
        acc = t_compute + t_wait + t_mpi + t_d2h
        rep = dict(rep, step_profile=True, barrier_in_step=True,
                   compute_ms_per_step=t_compute / nsteps * 1e3,
                   d2h_ms_per_step=t_d2h / nsteps * 1e3,
                   host_ms_per_step=(elapsed - acc) / nsteps * 1e3)
        print('driver.run_mpi rank %d STEP PROFILE (ms/step): compute %.2f '
              '(a+b, see comment) | barrier %.2f | mpi %.2f | halo d2h %.2f '
              '| host residual %.2f | TOTAL %.2f'
              % (rank, rep['compute_ms_per_step'], rep['wait_ms_per_step'],
                 rep['mpi_ms_per_step'], rep['d2h_ms_per_step'],
                 rep['host_ms_per_step'], rep['ms_per_step']), flush=True)
    if verbose:
        print('driver.run_mpi rank %d: %.3f s, %.3f ms/step (sync=%s: %.3f '
              'ms/step exchanging, %.3f ms/step waiting at the barrier), '
              'EFFECTIVE_CORES %.2f, threads=%d, cpus_allowed=%d %s'
              % (rank, elapsed, rep['ms_per_step'], sync,
                 rep['mpi_ms_per_step'], rep['wait_ms_per_step'],
                 eff, rep['threads'], rep['cpus_allowed'], rep['cpus']),
              flush=True)
    # velArr / dispArr / force are RANK-LOCAL (this rank's box numbering), so
    # they are returned under different keys than driver.run's global ones: a
    # consumer indexing them by a serial node or equation id would read a real
    # but WRONG row, and must fail with a KeyError instead. fric/fnft are this
    # rank's fault rows; `fault_rows` are the ones it WRITES.
    return dict(velArr_local=np.asarray(velArr),
                dispArr_local=np.asarray(dispArr),
                force_local=np.asarray(force),
                fnft=np.asarray(fnft), fric=np.asarray(fric),
                fault_rows=plan['fault_rows'],
                own_in_computed=plan['fault_rows'], report=rep)
