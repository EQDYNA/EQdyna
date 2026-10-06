#! /usr/bin/env python3
"""
TPV34 (SCEC, Imperial Fault Model 1: a planar fault in the 3D CVM-H
velocity structure) case tools -- TPV34_Description_v10 Parts 2 and 3.

The velocity structure is DATA: `tpv34_cvmh_grid_<dx>m.txt.gz`, extracted
once from CVM-H by extract_cvmh_grid.py (spec Part 2 recipe, provenance in
the file header and README.md) on the uniform grid of EQdyna ELEMENT
CENTRES of this case's box. It becomes par.mat, the n2mat == 6 3D material
grid (rows `x y z vp vs rho`), which the solver reads nearest-cell -- never
interpolated. Fault stresses follow from the same samples (spec Part 3):

    tau0   = 30.00 MPa * (mu/mu0)    (right-lateral, positive strike shear)
    sigma0 = 60.00 MPa * (mu/mu0)    (compressive; on_fault_vars slot 7 is
                                      signed negative)
    mu0    = 32.03812032 GPa,  mu = rho*Vs^2 at the fault node

where mu at a fault node is, as the spec suggests, the average of the
adjacent elements' mu on BOTH sides of the fault: the (up to) 8 grid cells
whose centres are at x +- dx/2, y +- dx/2, z +- dx/2 of the node (4 at the
free surface, where no element lies above). Nucleation is an ADDED shear
stress at t = 0 (not a forced-rupture-time formula, so swtwNucleation has
no TPV34 branch and par.C_nuclea = 0): with r = sqrt(x^2 + (depth-7500)^2),

    4.95 MPa * (mu/mu0)                                      r <= 1400 m
    2.475 MPa * (1 + cos(pi (r-1400)/600)) * (mu/mu0)        1400 < r <= 2000
    0                                                        r > 2000 m

Friction (linear slip-weakening): mu_s 0.580, mu_d 0.450, d0 0.18 m,
cohesion C0 = 0.000425 MPa/m * (2400 m - depth) for depth <= 2400 m, else 0
(1.02 MPa at the surface).

Grid use is DECIMATION-FREE by construction: the case dx must EQUAL the
shipped grid spacing (the grid is sampled at the solver's own spacing; a
finer dx would need a new extraction, a coarser one would need a new
extraction too -- subsampling a 500 m element-centre grid gives centres of
a 1000 m mesh only at every other point and would misplace half of them).
lib.requireFaultGeometryResolution refuses anything else.

Frame mapping (same proper rotation as tpv29/tpv35):
    x_eq = x_tpv (strike)   y_eq = z_tpv (fault-normal, + = far side)
    z_eq = -y_tpv (up)
"""
import gzip
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

SHIPPED_DX = 500.0
SHIPPED_NAME = 'the shipped CVM-H element-centre grid tpv34_cvmh_grid_500m.txt.gz'
GRID_FILE = os.path.join(HERE, 'tpv34_cvmh_grid_500m.txt.gz')

# Spec Part 3 constants (TPV34_Description_v10).
MU0 = 32.03812032e9          # Pa, reference shear modulus
TAU0_REF = 30.00e6           # Pa, initial shear stress at mu = mu0
SIGMA0_REF = 60.00e6         # Pa, initial normal stress (compressive) at mu = mu0
NUC_TAU = 4.95e6             # Pa, nucleation shear-stress increment at mu = mu0
NUC_R_IN, NUC_R_OUT = 1400.0, 2000.0   # m
HYPO_X, HYPO_DEPTH = 0.0, 7500.0       # m
MU_S, MU_D, D0 = 0.580, 0.450, 0.18
COHESION_GRAD = 0.000425e6   # Pa per m of (2400 - depth)
COHESION_DEPTH = 2400.0      # m


def loadGrid(path=GRID_FILE):
    """The shipped grid as an (nmat, 6) array [x y z vp vs rho] plus its
    header lines (provenance). Rows are placed by the solver's own grid
    rule, so file order is irrelevant."""
    header = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                header.append(line.rstrip('\n'))
            else:
                break
    data = np.loadtxt(path, comments='#')
    if data.ndim != 2 or data.shape[1] != 6:
        raise ValueError('%s: expected 6 columns x y z vp vs rho, got shape %r'
                         % (path, data.shape))
    return data, header


def materialGrid(dx, path=GRID_FILE):
    """par.mat for this case: the shipped grid, accepted only at the dx it
    was sampled at (see module docstring)."""
    from lib import requireFaultGeometryResolution
    requireFaultGeometryResolution(dx, SHIPPED_DX, SHIPPED_NAME, availableDx=[SHIPPED_DX])
    mat, _ = loadGrid(path)
    return mat


def gridLookup(mat):
    """dict of (x, y, z) -> mu (Pa) for every grid cell, plus the axis
    spacing, for the fault-node averaging below."""
    mu = mat[:, 5] * mat[:, 4] ** 2
    return {(round(r[0], 3), round(r[1], 3), round(r[2], 3)): m for r, m in zip(mat, mu)}


def faultNodeMu(mat, dx, fx, fz, fy=0.0):
    """mu (Pa) at every fault node [iz, ix] (iz = 0 at fzmin, EQdyna's
    par.fz order): the mean of the adjacent element-centre mu on both sides
    of the fault -- the spec's own suggested method. The grid holds exactly
    those centres (x +- dx/2, fy +- dx/2, z +- dx/2); cells at z > 0 do not
    exist (free surface) and are left out, so a surface node averages 4
    cells and an interior node 8."""
    lut = gridLookup(mat)
    h = 0.5 * dx
    out = np.empty((len(fz), len(fx)))
    for iz, z in enumerate(fz):
        for ix, x in enumerate(fx):
            vals = []
            for sx in (-h, h):
                for sy in (-h, h):
                    for sz in (-h, h):
                        key = (round(x + sx, 3), round(fy + sy, 3), round(z + sz, 3))
                        if key[2] > 0.0:
                            continue
                        if key not in lut:
                            raise ValueError('TPV34: fault node (%g, %g) needs grid cell %r, '
                                             'which the shipped grid does not hold' % (x, z, key))
                        vals.append(lut[key])
            out[iz, ix] = np.mean(vals)
    return out


def nucleationStress(x, depth):
    """Spec Part 3 added shear stress at mu = mu0 (Pa), r from the hypocentre."""
    r = np.sqrt((x - HYPO_X) ** 2 + (depth - HYPO_DEPTH) ** 2)
    out = np.zeros_like(r)
    out[r <= NUC_R_IN] = NUC_TAU
    ring = (r > NUC_R_IN) & (r <= NUC_R_OUT)
    out[ring] = 0.5 * NUC_TAU * (1.0 + np.cos(np.pi * (r[ring] - NUC_R_IN) / (NUC_R_OUT - NUC_R_IN)))
    return out


def cohesion(depth):
    """Spec Part 3 cohesion C0 (Pa): 0.000425 MPa/m * (2400 - depth) above 2400 m."""
    return np.where(depth <= COHESION_DEPTH, COHESION_GRAD * (COHESION_DEPTH - depth), 0.0)


def faultState(mat, dx, fx, fz):
    """(mu_ratio, tau0, sigma0, C0) on the fault grid [iz, ix]: mu/mu0,
    initial strike shear incl. nucleation (Pa, right-lateral positive),
    initial normal stress (Pa, compressive, positive here -- the caller
    negates for slot 7), cohesion (Pa)."""
    mu = faultNodeMu(mat, dx, fx, fz)
    ratio = mu / MU0
    X, Z = np.meshgrid(np.asarray(fx, float), np.asarray(fz, float))
    depth = -Z
    tau0 = (TAU0_REF + nucleationStress(X, depth)) * ratio
    sigma0 = SIGMA0_REF * ratio
    return ratio, tau0, sigma0, cohesion(depth)


def onFaultStationsKm(dz):
    """Spec Part 3's 35 on-fault stations: x in {-12,-6,0,6,12} km x
    down-dip {0, 1.0, 2.4, 5.0, 7.5, 10.0, 12.0} km. setOnFaultStation needs
    an exact fault node (it does not snap), so each depth is rounded to the
    nearest fault plane for THIS dz (2.4 km -> 2.5 km at the 500 m gate,
    file name dp025 instead of the spec's dp024)."""
    depths = (0.0, 1.0e3, 2.4e3, 5.0e3, 7.5e3, 10.0e3, 12.0e3)
    out = []
    for x in (-12.0, -6.0, 0.0, 6.0, 12.0):
        for d in depths:
            out.append([x, -round(d / dz) * dz / 1.0e3])
    return out


def offFaultStationsKm():
    """Spec Part 3's 56 off-fault stations body<Z>st<X>dp<D>: at depths 0
    and 2.4 km, fault-normal offsets z = +-3, +-9 km at x in {-20,-10,0,10,
    20} km, +-15 km at x in {-15, 0, 15} km, and 0 km (in the fault plane,
    beyond its ends) at x = +-20 km. [x_eq, y_eq, z_eq] in km, y_eq = z_tpv."""
    out = []
    for depth in (0.0, 2.4):
        for z in (-9.0, -3.0, 3.0, 9.0):
            for x in (-20.0, -10.0, 0.0, 10.0, 20.0):
                out.append([x, z, -depth])
        for z in (-15.0, 15.0):
            for x in (-15.0, 0.0, 15.0):
                out.append([x, z, -depth])
        for x in (-20.0, 20.0):
            out.append([x, 0.0, -depth])
    return out
