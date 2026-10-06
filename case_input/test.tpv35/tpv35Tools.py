#! /usr/bin/env python3
"""
TPV35 (SCEC, Parkfield 2004 M6 validation benchmark) official-data tools.

SCEC supplies TPV35's fault state and material as DATA, not formulas
(tpv35_data_files.zip from https://strike.scec.org/cvws/tpv35docs.html,
sha256 c8d4c28d...52c2f3, fetched 2026-10-05). The four files ship verbatim
in this directory:

    tpv35_input_data.txt                  401 x 156 grid at 100 m:
                                          `nx ny x y mu_s tau0(MPa)`
    tpv35_station_locations.txt           43 surface stations, `name code z x`
    tpv35_velocity_structure_near_side.txt 1D layers (thick Vp Vs rho), z<0 side
    tpv35_velocity_structure_far_side.txt  1D layers (thick Vp Vs rho), z>0 side

Frame mapping (proper rotation, det=+1, identical to tpv29GeometryTools):
    x_eq = x_tpv (strike)   y_eq = z_tpv (fault-normal)   z_eq = -y_tpv (up)
so the spec's "near side" (z_tpv < 0) is y_eq < 0, "far side" is y_eq > 0.

The fault-state grid is DECIMATED, never interpolated (PROJECT_RULES rule 17
step 2): `faultGridForCase` takes every n-th official value, so a case dx
must be an integer multiple of the shipped 100 m AND tile the 40 km x 15.5 km
fault exactly -- that leaves 100 m (spec) and 500 m (gate). Anything else is
refused by lib.requireFaultGeometryResolution / this module with the reason.
The spec says "interpolate for non-grid points"; this case never has any.

Nucleation is built into the data (spec Part 3): inside a ~1.4 km patch
around (x=0, depth 8.1 km) mu_s*sigma_n < tau0, so those nodes fail at t=0
with no artificial nucleation (par.C_nuclea = 0).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

SHIPPED_DX = 100.0
SHIPPED_NAME = 'the official TPV35 input grid tpv35_input_data.txt'
INPUT_DATA = os.path.join(HERE, 'tpv35_input_data.txt')
STATIONS = os.path.join(HERE, 'tpv35_station_locations.txt')
VELOCITY = {-1: os.path.join(HERE, 'tpv35_velocity_structure_near_side.txt'),
            +1: os.path.join(HERE, 'tpv35_velocity_structure_far_side.txt')}

# Spec constants (TPV35_Description_v05, Part 2/3), not in the data files.
NORMAL_STRESS = 60.0e6   # Pa, constant over the fault
MU_D = 0.30
D0 = 0.15                # m
COHESION = 0.0           # Pa


def loadInputGrid(path=INPUT_DATA):
    """The official grid as dict(x (nx,), depth (nd,), mu_s (nd,nx),
    tau0 (nd,nx) in Pa). Rows are placed by their own (nx, ny) indices, so
    file order is irrelevant; every cell must be filled exactly once."""
    with open(path) as f:
        head = f.readline().split()
    nxm1, ndm1 = int(head[0]), int(head[1])
    xmin, xmax, dmin, dmax = (float(v) for v in head[2:6])
    data = np.loadtxt(path, skiprows=1)
    if data.shape != ((nxm1 + 1) * (ndm1 + 1), 6):
        raise ValueError('%s: header promises %d x %d nodes with 6 columns, '
                         'found shape %r' % (path, nxm1 + 1, ndm1 + 1, data.shape))
    ix = data[:, 0].astype(int)
    idp = data[:, 1].astype(int)
    x = np.linspace(xmin, xmax, nxm1 + 1)
    depth = np.linspace(dmin, dmax, ndm1 + 1)
    if (np.abs(data[:, 2] - x[ix]).max() > 1e-6 or
            np.abs(data[:, 3] - depth[idp]).max() > 1e-6):
        raise ValueError('%s: (nx, ny) indices disagree with the x, y columns' % path)
    mu_s = np.full((ndm1 + 1, nxm1 + 1), np.nan)
    tau0 = np.full((ndm1 + 1, nxm1 + 1), np.nan)
    mu_s[idp, ix] = data[:, 4]
    tau0[idp, ix] = data[:, 5] * 1.0e6
    if np.isnan(mu_s).any() or len(np.unique(ix * 10 ** 6 + idp)) != data.shape[0]:
        raise ValueError('%s: grid not covered exactly once' % path)
    return dict(x=x, depth=depth, mu_s=mu_s, tau0=tau0)


def faultGridForCase(dx, fx, fz, path=INPUT_DATA):
    """mu_s and tau0 (Pa) on the case fault grid [iz, ix], iz=0 at fzmin
    (EQdyna's par.fz order, ascending z = deepest first), by exact
    integer-stride decimation of the official grid. Refuses any dx the
    source cannot give exactly and any fault box that is not the official
    one."""
    from lib import requireFaultGeometryResolution, FaultGeometryError
    requireFaultGeometryResolution(dx, SHIPPED_DX, SHIPPED_NAME,
                                   availableDx=[SHIPPED_DX])
    g = loadInputGrid(path)
    stride = int(round(dx / SHIPPED_DX))
    xs, ds = g['x'][::stride], g['depth'][::stride]
    fx = np.asarray(fx, dtype=float)
    fz = np.asarray(fz, dtype=float)
    if (xs.shape != fx.shape or np.abs(xs - fx).max() > 1e-6):
        raise FaultGeometryError(
            'TPV35: the case fault x grid (%d nodes, %g..%g) is not the %d-stride '
            'decimation of the official grid (%d nodes, %g..%g). The fault must '
            'span exactly x in [%g, %g] m.' % (fx.size, fx[0], fx[-1], stride, xs.size,
                                              xs[0], xs[-1], g['x'][0], g['x'][-1]))
    if (ds.shape != fz.shape or np.abs(ds[::-1] + fz).max() > 1e-6):
        raise FaultGeometryError(
            'TPV35: the case fault z grid (%d nodes, %g..%g) is not the %d-stride '
            'decimation of the official depth grid (%d nodes, 0..%g). The fault '
            'must span exactly z in [%g, 0] m.' % (fz.size, fz[0], fz[-1], stride,
                                                  ds.size, ds[-1], -g['depth'][-1]))
    mu_s = g['mu_s'][::stride, ::stride][::-1, :]
    tau0 = g['tau0'][::stride, ::stride][::-1, :]
    return mu_s, tau0


def readVelocityTable(path):
    """[(thickness or None for the halfspace, vp, vs, rho), ...] top-down."""
    rows = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0].lower().startswith('thick'):
                continue
            thick = None if s[0] == '-' else float(s[0])
            rows.append((thick, float(s[1]), float(s[2]), float(s[3])))
    if rows[-1][0] is not None or any(r[0] is None for r in rows[:-1]):
        raise ValueError('%s: expected exactly one trailing halfspace row' % path)
    return rows


def materialTable(halfspaceBottom):
    """The n2mat==5 two-sided material table par.mat, rows
    [bottom depth (m), vp, vs, rho, side] with side -1 (y_eq<0, spec near
    side) and +1 (y_eq>0, far side); per-side bottoms ascending, the
    halfspace closed at `halfspaceBottom` (any depth >= the model bottom)."""
    rows = []
    for side in (-1, +1):
        bottom = 0.0
        for thick, vp, vs, rho in readVelocityTable(VELOCITY[side]):
            bottom = halfspaceBottom if thick is None else bottom + thick
            rows.append([bottom, vp, vs, rho, float(side)])
    return np.array(rows)


def offFaultStationsKm(path=STATIONS):
    """[[x_eq, y_eq, z_eq] in km, ...] for the 43 official surface stations
    in file order; y_eq = z_tpv, all at the free surface."""
    out = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s or s[0] == 'station_name':
                continue
            z_tpv, x_tpv = float(s[2]), float(s[3])
            out.append([x_tpv / 1.0e3, z_tpv / 1.0e3, 0.0])
    return out
