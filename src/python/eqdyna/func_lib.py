"""func_lib.py <- src/func_lib.f90.

Two routines:
  `region_damp`          <- `pmlRegionDistance`, vectorized: classifies a
                            point into one of the PML regions and returns the
                            three damping coefficients.
  `dev_str_depth_taper`  <- `devStrDepthTaper`: SCEC TPV29/30's Omega(depth)
                            taper on the off-fault deviatoric pre-stress.

As of pathway_forward.md item 13, `pmlRegionDistance` classifies BOTH its
callers -- node-based (computePMLDampingVector.f90/velDispUpdate) and
element-center-based (assembleGlobalKU.f90's calcPMLElemKU) -- with
INCLUSIVE (>=/<=) boundary tests, and no longer takes a boundInclusive
argument at all. This matches that: always inclusive, no parameter.

Before that change this port had FOLLOWED an earlier Fortran where the two
call sites genuinely differed (element centers strict, nodes inclusive).
Carrying that distinction now would be reproducing a bug that was fixed
upstream. (Element centers are additionally guaranteed never to land exactly
on a PML bound by checkPMLAlignment, so it was already a no-op for that
caller; parity confirms it.)
"""
import numpy as np


def region_damp(x, y, z, PMLb, nPML, vmaxPML, R):
    """Vectorized port of src/func_lib.f90's `pmlRegionDistance` (current
    master: region-14 tests y against ymax0, not xmax0 -- pathway_forward.md
    item 8's bug, fixed upstream in a844e92; boundary tests are INCLUSIVE
    (>=/<=) for both callers -- pathway_forward.md item 13, standardized
    2026-09-13. There is no longer a bound_inclusive distinction to
    reproduce."""
    xmax0, xmin0, ymax0, ymin0, zmin0, maxdx, maxdy, maxdz = PMLb
    n = x.shape[0]
    xHi = x >= xmax0; xLo = x <= xmin0; yHi = y >= ymax0; yLo = y <= ymin0

    d3 = np.where(z <= zmin0, np.abs(z - zmin0), 0.0)

    d1 = np.zeros(n); d2 = np.zeros(n)
    taken = np.zeros(n, dtype=bool)
    for cond, a, b in [
        (xHi & yHi, np.abs(x - xmax0), np.abs(y - ymax0)),          # region 11
        (xHi & yLo, np.abs(x - xmax0), np.abs(y - ymin0)),          # region 12
        (xLo & yLo, np.abs(x - xmin0), np.abs(y - ymin0)),          # region 13
        (xLo & yHi, np.abs(x - xmin0), np.abs(y - ymax0)),          # region 14 (fixed)
    ]:
        sel = cond & ~taken
        d1 = np.where(sel, a, d1); d2 = np.where(sel, b, d2); taken |= sel
    sel = ~taken & xHi & (y > ymin0) & (y < ymax0)                  # region 1_12
    d1 = np.where(sel, np.abs(x - xmax0), d1); d2 = np.where(sel, 0.0, d2); taken |= sel
    sel = ~taken & yLo & (x > xmin0) & (x < xmax0)                  # region 1_23
    d1 = np.where(sel, 0.0, d1); d2 = np.where(sel, np.abs(y - ymin0), d2); taken |= sel
    sel = ~taken & xLo & (y > ymin0) & (y < ymax0)                  # region 1_34
    d1 = np.where(sel, np.abs(x - xmin0), d1); d2 = np.where(sel, 0.0, d2); taken |= sel
    sel = ~taken & yHi & (x > xmin0) & (x < xmax0)                  # region 1_41
    d1 = np.where(sel, 0.0, d1); d2 = np.where(sel, np.abs(y - ymax0), d2); taken |= sel
    # else (middle area 9): d1=d2=0, already the default.

    out = []
    for d, delta in ((d1, nPML * maxdx), (d2, nPML * maxdy), (d3, nPML * maxdz)):
        out.append(3.0 * vmaxPML / 2.0 / delta * np.log(1.0 / R) * (d / delta) ** 2.0)
    return out  # damp1, damp2, damp3


def dev_str_depth_taper(depth, taper_start, taper_end):
    """Port of src/fortran/func_lib.f90's `devStrDepthTaper` -- SCEC
    TPV29/TPV30's Omega(depth) (spec part 4, "Initial Stress Tensor"): the
    off-fault DEVIATORIC pre-stress tapers linearly to zero between
    `taper_start` and `taper_end` (depths in m, positive down) while the
    vertical lithostatic component keeps growing. TPV30's values are 17000 m
    and 22000 m.

    The taper is INACTIVE when taper_end <= taper_start (what case.setup
    writes as `0.0 0.0` for a case that does not configure it), and then this
    returns exactly 1.0 -- an IEEE-exact multiplicative identity, which is
    what keeps a case that does not ask for the taper bit-for-bit identical to
    the pre-v5.9.0 solver.

    `depth` may be a scalar or an array; the return has its shape.
    """
    if not taper_end > taper_start:
        return np.ones_like(np.asarray(depth, dtype=np.float64))
    omega = (taper_end - np.asarray(depth, dtype=np.float64)) / (taper_end - taper_start)
    return np.clip(omega, 0.0, 1.0)
