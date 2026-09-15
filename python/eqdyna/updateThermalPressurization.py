"""updateThermalPressurization.py <- src/updateThermalPressurization.f90.

Pore-pressure and temperature rise from thermal pressurization, TPV105 3D
benchmark equations 12 and 13. Called from driver.f90:28, BEFORE faulting,
so the history it integrates covers steps 1..nt-1 only: the current step's
slip rate and shear traction do not exist yet.

ONE implementation, both backends -- there is nothing backend-specific here
once the history sum is written over a fixed width (see below).

fric(FRIC_SLOT_TP_H) is a PER-NODE array, not the module-level scalar
globalvar.f90 used to carry. That scalar was never assigned anywhere in
src/*.f90 (confirmed by compiling a minimal program that links only
globalvar.f90 and prints it), so it was identically 0.0 and the whole
kernel degenerated. This port reproduces the FIXED semantics.

COST, stated because it bounds where this port can be used: the Fortran
allocates onFaultTPHist at (2, nftmx, nstep, ntotft) up front and the inner
sum runs j=1..nt-1 every step, so memory is O(nftnd*nstep) and compute is
O(nftnd*nstep^2). For test.tpv1053d (nftnd=4005, nstep=120) that is ~7.7 MB
and trivial. For a production run (nftnd~50,000, nstep~200,000) it is
~160 GB and quadratic time. Neither a bounded-lookback truncation nor a
recursive reformulation exists in the Fortran, so neither is invented here;
this is an open limitation, not a silent one.
"""
from . import backend as B
from . import globalvar as gv


def build(S, nsteps):
    """Per-node TP constants. Read-only for the whole run."""
    fric = S['fric_init']
    return dict(
        gama=fric[:, gv.TP_LAMBDA] / fric[:, gv.TP_ROUC],
        omega=fric[:, gv.TP_A_HY],
        kapa=fric[:, gv.TP_A_TH],
        rouc=fric[:, gv.TP_ROUC],
        Tini=fric[:, gv.TP_TINI],
        h=fric[:, gv.TP_H].copy(),
        nsteps=nsteps,
    )


def updateThermalPressurization(xp, tp, fric, sliprate_hist, shear_hist, nt, dt):
    """updateThermalPressurization.f90:14-38, vectorized over fault nodes AND
    over the history index j.

    THE FIXED-WIDTH SUM. The Fortran sums j = 1..nt-1. Here the sum runs over
    the FULL preallocated width every step, with columns j >= nt-1 forced to
    contribute exactly 0.0 by an explicit xp.where. That is what lets one
    expression serve both backends: nt is a traced value under jit, so a
    variable-length slice is not available there, and a variable-length slice
    on numpy plus a masked full-width sum on jax would be two different
    summation trees and therefore two different answers.

    The mask is NOT redundant with "the history is still zero there". `age`
    goes negative for those columns, which makes the kernel's sqrt argument
    negative (NaN) or zero (inf), and NaN*0.0 is NaN, not 0.0 -- an unmasked
    sum would be poisoned rather than merely wasteful. The where is the
    guard, and it is exact: those columns contribute a literal 0.0.
    """
    gama = tp['gama']; omega = tp['omega']; kapa = tp['kapa']
    rouc = tp['rouc']; h = tp['h']
    j1 = xp.arange(tp['nsteps'])              # 0-based; Fortran's j = j1 + 1
    valid = j1 < (nt - 1)
    age = (nt - (j1 + 1)) * dt                # (nt-j)*dt
    age = xp.where(valid, age, 0.0)

    h2 = 2.0 * (h ** 2)
    denom_k = 4.0 * kapa[:, None] * age[None, :] + h2[:, None]
    denom_o = 4.0 * omega[:, None] * age[None, :] + h2[:, None]
    dom = (omega - kapa)[:, None]
    ker1 = -kapa[:, None] / dom / xp.sqrt(denom_k) + omega[:, None] / dom / xp.sqrt(denom_o)
    ker2 = 1.0 / xp.sqrt(denom_k)

    hist = xp.abs(shear_hist) * sliprate_hist
    hist = xp.where(valid[None, :], hist, 0.0)

    patnode = (hist * ker1 * dt).sum(axis=1) * gama / xp.sqrt(xp.pi)
    Tatnode = (hist * ker2 * dt).sum(axis=1) / rouc / xp.sqrt(xp.pi)

    fric = B.setat(xp, fric, (slice(None), gv.TP_NORM_TP), patnode)
    fric = B.setat(xp, fric, (slice(None), gv.TP_TEMP), Tatnode + tp['Tini'])
    return fric
