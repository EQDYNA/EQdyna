"""Forced-rupture nucleation, written once for every backend.

Port of faulting.f90:swtwNucleation (called from solveSWTW, i.e. friclaw 1 and
2, whenever C_nuclea==1 and this is nucfault).

WHY THIS IS ITS OWN MODULE

The friclaw-1 solver exists twice -- port.py (NumPy) and port_jax.py (JAX) --
and so do the friclaw-4 and friclaw-5 solvers. Fortran has ONE loop with
`if (friclaw<=2) call solveSWTW` inside it; the port hoisted that dispatch to
module level and duplicated the loop three times, then again per backend.

That cost is not theoretical. The nucleation branch below was missing from the
port entirely: test.tpv29 nucleated 0 of 3321 fault nodes against a reference
of 2974, which read as a 7.1e7 physics failure. When the fix was written into
port.py alone, tpv29 x python-numpy passed and tpv29 x python-jax kept failing
at the identical value -- one bug needing two fixes.

So it lives here, parameterised by the array module, and both backends call it.
`xp` is numpy or jax.numpy: every operation used below (sqrt, where, full,
errstate-free arithmetic) has the same name and semantics in both.

THE PHYSICS

The smoothed forced-rupture time is the TPV29 spec formula
(TPV29_30_Description_v06, Part 6):

    T(r) = (r + 0.081 * r_crit * (1/(1 - (r/r_crit)^2) - 1)) / (0.7 * Vs)

with Vs fixed at 3464 m/s. TPV36/37/201 reuse it; TPV202 uses a plain
r/nucRuptVel instead. Inside r <= r_crit the forced rupture is imposed, and the
friction coefficient is driven from static to dynamic over t0 seconds
(FRIC_SLOT_TW_T0).

Outside that set of TPVs, tr stays at its 1e9 default and the whole thing
degenerates to `min(fs, fricCoeff)` -- which is what the port used to do
unconditionally, and why nothing nucleated for a case that needed the formula.
"""

# globalvar.f90:79-81
NUC_VS_FIXED = 3464.0      # fixed shear-wave speed used by the formula, m/s
NUC_TAPER_COEF = 0.081     # rupture-time taper coefficient
NUC_VR_TO_VS = 0.7         # nucleation rupture-speed-to-Vs ratio

# TPVs whose spec uses the smoothed forced-rupture time. 29 is the SOURCE of
# the formula and was missing from Fortran's own list, which forced
# case_input/test.tpv29 to declare par.tpv = 36 to reach its own physics.
SMOOTHED_TPVS = (29, 36, 37, 201)
CONSTANT_VELOCITY_TPVS = (202,)


def source_radius(xp, coords, xsource, ysource, zsource):
    """Distance from each fault node to the hypocentre.

    `coords` is the SLAVE node coordinate, matching faulting.f90:388's
    meshCoor(:, nsmp(1,...)).  Loop-invariant: compute once, outside the step.
    """
    return xp.sqrt((coords[:, 0] - xsource) ** 2
                   + (coords[:, 1] - ysource) ** 2
                   + (coords[:, 2] - zsource) ** 2)


def forced_rupture_time(xp, radius, TPV, nucR, nucRuptVel):
    """tr per fault node. 1e9 (i.e. never) outside the nucleation patch.

    Loop-invariant -- depends only on geometry -- so callers should compute it
    once rather than per step.

    At r == nucR exactly the taper divides by zero and Fortran yields +Inf, so
    that node never forces. Reproduced rather than guarded, so the edge matches
    bit for bit.
    """
    never = xp.full(radius.shape, 1.0e9)
    inside = radius <= nucR
    if TPV in SMOOTHED_TPVS:
        ratio = xp.where(inside, radius / nucR, 0.0)
        taper = 1.0 / (1.0 - ratio ** 2) - 1.0
        tr = (radius + NUC_TAPER_COEF * nucR * taper) / (
            NUC_VR_TO_VS * NUC_VS_FIXED)
        return xp.where(inside, tr, never)
    if TPV in CONSTANT_VELOCITY_TPVS:
        return xp.where(inside, radius / nucRuptVel, never)
    return never


def apply(xp, fricCoeff, fs, fd, tw_t0, tr, timeElapsed):
    """faulting.f90:399-405 -- drive mu from fs toward fd over t0 after tr.

    Returns the updated friction coefficient. With tr at its 1e9 default this
    reduces exactly to min(fs, fricCoeff), the no-forcing case.
    """
    tc = xp.where(timeElapsed < tr, 0.0,
                  xp.where(timeElapsed < tr + tw_t0,
                           (timeElapsed - tr) / tw_t0, 1.0))
    return xp.minimum(fs + (fd - fs) * tc, fricCoeff)


def enabled(S):
    """Whether swtwNucleation does anything for this case."""
    return (S.get('C_nuclea', 0) == 1 and S.get('nucfault', 1) == 1
            and S.get('nucR', 0.0) > 0.0)
