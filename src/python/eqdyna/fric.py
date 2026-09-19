"""fric.py <- src/fric.f90. The four friction laws, one function each.

Every function takes the array module `xp` (numpy or jax.numpy) as its
first argument and is written ONCE. There is no numpy copy and no jax copy:
the expressions below are valid in both, so there is nothing to keep in
sync.

Vectorized over fault nodes rather than called per node as the Fortran is.
The scalar Fortran's if/elseif chains become xp.where, which evaluates both
arms -- see slip_weak's note for why that is exact here and where it is
merely harmless.
"""
from . import globalvar as gv


def slip_weak(xp, slip, fric, xmu_unused=None):
    """fric.f90:3-19 -- slip-weakening friction coefficient.

    The Fortran is an if/elseif/endif followed by a SEPARATE `if`:

        if      (abs(slip) < 1e-10)  xmu = fs
        elseif  (slip < D0)          xmu = fs - (fs-fd)*slip/D0
        endif
        if      (slip >= D0)         xmu = fd

    Note what that structure does NOT do: when the first `if` chain falls
    through (abs(slip) >= 1e-10 and slip >= D0), `xmu` is left at whatever
    the second `if` assigns -- and the second `if` covers exactly that case,
    so xmu is always assigned. Reproduced with two xp.where in the same
    order, which is what makes the second one override the first.

    Evaluating both arms of the middle branch divides by D0 for every node,
    including nodes where the branch is not taken. D0 is a per-node input
    and is nonzero for every case in the suite; if it were zero the Fortran
    would not evaluate that arm and this would produce inf. That is not a
    silent difference -- the xp.where discards it -- but it is the reason
    this function does not guard the divide: guarding it would change the
    taken arm's value, and the taken arm is the one that has to match.
    """
    fs = fric[:, gv.SW_FS]
    fd = fric[:, gv.SW_FD]
    D0 = fric[:, gv.SW_D0]
    xmu = xp.where(xp.abs(slip) < 1.0e-10, fs, fs - (fs - fd) * slip / D0)
    return xp.where(slip >= D0, fd, xmu)


def time_weak(xp, trupt, fric):
    """fric.f90:21-37 -- time-weakening friction coefficient.

        if      (trupt <= 0)   xmu = fs
        elseif  (trupt < t0)   xmu = fs - (fs-fd)*trupt/t0
        else                   xmu = fd

    This is slip_weak with `slip/D0` replaced by `trupt/TW_T0` and with the
    first test on the SIGN of trupt rather than on |slip| vs 1e-10 -- an
    if/elseif/ELSE, so unlike slip_weak the three arms are exhaustive and
    ordered, and the last is a plain else rather than a second overriding
    `if`.

    `trupt = timeElapsed - fnft` (faulting.f90:151). An unruptured node has
    fnft at its 99999.0 sentinel, so trupt is hugely negative, so the first
    arm fires and xmu = fs. That is why unruptured nodes need no special
    case here -- the sentinel does the work. It also means TW_T0 is divided
    into a negative number for those nodes in the untaken middle arm, which
    xp.where discards.
    """
    fs = fric[:, gv.SW_FS]
    fd = fric[:, gv.SW_FD]
    t0 = fric[:, gv.TW_T0]
    return xp.where(trupt <= 0.0, fs,
                    xp.where(trupt < t0, fs - (fs - fd) * trupt / t0, fd))


def rate_state_ageing_law(xp, V2, theta, fric, dt):
    """fric.f90:39-61 -- ageing law (friclaw==3).

    Returns (xmu, dxmudv, theta_new). The Fortran mutates `theta` in place;
    here it is returned, because the caller (NewtonRaphson) re-seeds it from
    a frozen baseline on every iteration and the in-place form hid that.
    """
    A = fric[:, gv.RSF_A]; B = fric[:, gv.RSF_B]; L = fric[:, gv.RSF_DC]
    f0 = fric[:, gv.RSF_R0]; V0 = fric[:, gv.RSF_V0]

    tmpc = 1.0 / (2.0 * V0) * xp.exp((f0 + B * xp.log(V0 * theta / L)) / A)
    tmp = (V2 + 1.0e-30) * tmpc                      # the 1e-30 is the Fortran's, verbatim
    xmu = A * xp.log(tmp + xp.sqrt(tmp ** 2 + 1.0))  # arcsinh(z) = ln(z + sqrt(z^2+1))
    dxmudv = A * tmpc / xp.sqrt(1.0 + tmp ** 2)
    theta_new = L / V2 + (theta - L / V2) * xp.exp(-V2 * dt / L)
    return xmu, dxmudv, theta_new


def rate_state_slip_law(xp, V2, psi, fric, dt):
    """fric.f90:63-95 -- slip law with strong rate weakening (friclaw 4/5).

    Returns (xmu, dxmudv, psi_new).

    psiss uses the Fortran's EXPANDED sinh -- `(exp(x) - exp(-x))/2` -- not
    xp.sinh. The commented-out `dsinh` form sits directly above it in
    fric.f90:91 and the expanded form is what is compiled; the two differ in
    the last bits, and substituting the library sinh here is exactly the
    "library function that shares a name" trap.
    """
    A = fric[:, gv.RSF_A]; B = fric[:, gv.RSF_B]; L = fric[:, gv.RSF_DC]
    f0 = fric[:, gv.RSF_R0]; V0 = fric[:, gv.RSF_V0]
    fw = fric[:, gv.RSF_FW]; Vw = fric[:, gv.RSF_VW]

    tmpc = 1.0 / (2.0 * V0) * xp.exp(psi / A)
    tmp = (V2 + 1.0e-30) * tmpc
    xmu = A * xp.log(tmp + xp.sqrt(tmp ** 2 + 1.0))
    dxmudv = A * tmpc / xp.sqrt(1.0 + tmp ** 2)
    fLV = f0 - (B - A) * xp.log(V2 / V0)
    fss = fw + (fLV - fw) / ((1.0 + (V2 / Vw) ** 8) ** 0.125)
    fssa = fss / A
    psiss = A * xp.log(2.0 * V0 / V2 * (xp.exp(fssa) - xp.exp(-fssa)) / 2.0)
    psi_new = psiss + (psi - psiss) * xp.exp(-V2 * dt / L)
    return xmu, dxmudv, psi_new


def rate_state_normal_stress(xp, V2, theta_pc, tnrm, fric, dt):
    """faulting.f90:39-56 -- Shi & Day (2013) eq. B8 normal-stress state.

    Lives in faulting.f90 rather than fric.f90 in the Fortran, but it is a
    constitutive update like the three above and is only ever called from
    the RSF solver, so it is kept with them. Returns (theta_pc_new,
    theta_pc_dot); L is Dc, reused as L_pc (the Fortran's own comment).
    """
    L = fric[:, gv.RSF_DC]
    theta_pc_dot = -V2 / L * (theta_pc - xp.abs(tnrm))
    return theta_pc + theta_pc_dot * dt, theta_pc_dot
