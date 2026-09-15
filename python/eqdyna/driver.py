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
faulting.f90:17-18 puts it -- inside faulting -- and because friclaw is a
static Python int it is resolved at TRACE time under jit and costs nothing.

MPI is absent by construction: this port is serial (npx=npy=npz=1, enforced
by standalone/main.py), so driver.f90:27's MPI4NodalQuant has no work to do.
That is a scope limit, not an omission that could silently mislead -- a
multi-rank case is refused before it reaches here.
"""
import time

import numpy as np

from . import assembleGlobalKU as KU
from . import backend as B
from . import faulting as FLT
from . import globalvar as gv
from . import updateThermalPressurization as TP


def velDispUpdate(xp, inv, v1, velArr, dispArr, force, dt):
    """driver.f90:101-169 -- integrate velocity and displacement.

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


def make_step(xp, inv, finv, tp, mass, scratch):
    """Build the per-step closure. Called by backend.run_time_loop INSIDE
    the jit on the jax path, so that `inv`'s arrays resolve to jit arguments
    rather than closed-over HLO literals."""
    dt = inv['dt']; rdampk = inv['rdampk']
    tr = finv['tr']
    friclaw = finv['friclaw']

    def step(carry, nt):
        (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
         sliprate_hist, shear_hist) = carry
        timeElapsed = timeElapsed + dt                       # driver.f90:12

        v1, velArr, dispArr = velDispUpdate(xp, inv, v1, velArr, dispArr,
                                            force, dt)

        force = B.setat(xp, force, slice(None), 0.0)         # driver.f90:23
        force, stress_i, s_p = KU.assembleGlobalKU(
            xp, inv, velArr, force, stress_i, s_p, dt, rdampk, scratch)
        force = KU.calcHourglassResist(xp, inv, dispArr, velArr, force, rdampk)

        if friclaw == 5:                                     # driver.f90:28
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

    return step


def run(S, nsteps=None, verbose=True, xp=np):
    """The whole solve. `xp` selects the backend: numpy or jax.numpy.

    Returns the same dict the six modules this replaces returned, so
    standalone/main.py and the frt writer need no change.
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
              xp.asarray(inv['stress_i0']), z((inv['Ep'], 15)),
              xp.asarray(S['fric_init'].copy()),
              xp.full(nftnd, gv.FNFT_SENTINEL),
              xp.asarray(0.0),
              z((nftnd, hist_w)), z((nftnd, hist_w)))

    if verbose:
        print('driver.run: %d steps, backend=%s, friclaw=%d'
              % (nsteps, xp.__name__, S['friclaw']))
    t0 = time.perf_counter()
    scratch = KU.alloc_scratch(xp, inv)
    carry = B.run_time_loop(xp, lambda i: make_step(xp, i, finv, tp, mass, scratch),
                            inv, carry0, nsteps)
    elapsed = time.perf_counter() - t0
    if verbose:
        print('driver.run: %.3f s, %.3f ms/step' % (elapsed, elapsed / nsteps * 1e3))

    (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed,
     sliprate_hist, shear_hist) = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr),
                fnft=np.asarray(fnft), fric=np.asarray(fric),
                force=np.asarray(force))
