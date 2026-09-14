#! /usr/bin/env python3
"""
TPV29 official rough-fault geometry tools.

The SCEC TPV29/30 benchmarks (strike.scec.org/cvws, TPV29_30_Description_v06)
supply the rough fault surface z_tpv = f(x, y_tpv) on a 25 m grid
(tpv29_tpv30_geometry_25m_data.txt), where in the TPV frame x is
along-strike, y_tpv is depth (positive down) and z_tpv is fault-normal
(positive toward the far side).

EQdyna frame mapping (proper rotation, det = +1):
    x_eq = x_tpv (strike),  y_eq = z_tpv (fault-normal),  z_eq = -y_tpv (up)
so the EQdyna fault surface is y_eq = f with
    dy/dx_eq = df/dx_tpv     and     dy/dz_eq = -df/dy_tpv.

This module
  1. loads an EQdyna-format geometry file (the compset ships the official
     surface at 100 m as bFault_Rough_Geometry.tpv29.100m.txt, produced by
     exact 4:1 decimation of the official 25 m file -- no interpolation, all
     values are the official ones);
  2. decimates it to the case's mesh spacing par.dx (which must be an integer
     multiple of the shipped spacing and divide the 40 km x 20 km fault); and
  3. writes bFault_Rough_Geometry.txt in the format read by
     src/readInputFiles.f90:read_fault_rough_geometry /
     src/func_lib.f90:insertFaultInterface:
         row 0: nnx nnz 0
         row 1: dx fxmin fzmin
         rows:  [y, dy/dx, dy/dz], z-fastest (index nnz*(ix) + iz, iz counted
                upward from fzmin), x from fxmin.

Run as a script (in a case directory, after editing par.dx in
user_defined_params.py) it regenerates bFault_Rough_Geometry.txt:
    python3 tpv29GeometryTools.py

convertOfficial25m() converts the raw SCEC 25 m download into the shipped
EQdyna-format file; it is retained for provenance and is not needed to run a
case.
"""
import numpy as np

FXMIN, FXMAX = -20000.0, 20000.0
FZMIN, FZMAX = -20000.0, 0.0
SHIPPED = 'bFault_Rough_Geometry.tpv29.100m.txt'


def loadEQdynaGeometry(fname):
    """Read an EQdyna bFault_Rough_Geometry-format file.
    Returns (dx, y, dydx, dydz) with arrays shaped (nz, nx), index [iz, ix],
    iz = 0 at fzmin (bottom), ix = 0 at fxmin."""
    with open(fname) as f:
        nnx, nnz = [int(float(v)) for v in f.readline().split()[:2]]
        dx, fxmin, fzmin = [float(v) for v in f.readline().split()[:3]]
    if abs(fxmin - FXMIN) > 1e-6 or abs(fzmin - FZMIN) > 1e-6:
        raise ValueError(f'{fname}: fault corner ({fxmin},{fzmin}) is not the '
                         f'TPV29 corner ({FXMIN},{FZMIN})')
    data = np.loadtxt(fname, skiprows=2)
    if data.shape[0] != nnx * nnz:
        raise ValueError(f'{fname}: {data.shape[0]} rows, expected {nnx*nnz}')
    # file order is z-fastest within each x column
    y = data[:, 0].reshape(nnx, nnz).T.copy()
    dydx = data[:, 1].reshape(nnx, nnz).T.copy()
    dydz = data[:, 2].reshape(nnx, nnz).T.copy()
    return dx, y, dydx, dydz


def decimateToDx(srcDx, y, dydx, dydz, dx):
    """Exact decimation from spacing srcDx to spacing dx (no interpolation).
    dx must be an integer multiple of srcDx and divide the fault extents."""
    step = dx / srcDx
    if abs(step - round(step)) > 1e-9:
        raise ValueError(f'dx={dx} is not an integer multiple of the shipped '
                         f'geometry spacing {srcDx} m')
    step = int(round(step))
    for extent in (FXMAX - FXMIN, FZMAX - FZMIN):
        if abs(extent / dx - round(extent / dx)) > 1e-9:
            raise ValueError(f'dx={dx} does not divide the fault extent {extent} m')
    return y[::step, ::step], dydx[::step, ::step], dydz[::step, ::step]


def writeBFault(y, dydx, dydz, dx, out='bFault_Rough_Geometry.txt'):
    """Write EQdyna's bFault_Rough_Geometry.txt (z-fastest ordering)."""
    nz, nx = y.shape
    with open(out, 'w') as f:
        f.write(f'{nx}\t{nz}\t0\n')
        f.write(f'{dx:.6f}\t{FXMIN:.6f}\t{FZMIN:.6f}\n')
        for ix in range(nx):
            for iz in range(nz):
                f.write(f'{y[iz, ix]:.7e}\t{dydx[iz, ix]:.7e}\t{dydz[iz, ix]:.7e}\n')


def faultGridForCase(dx):
    """Return (y, dydx, dydz) arrays shaped (nfz, nfx), [iz, ix] with iz=0 at
    fzmin, decimated from the shipped official geometry to spacing dx."""
    srcDx, y, dydx, dydz = loadEQdynaGeometry(SHIPPED)
    return decimateToDx(srcDx, y, dydx, dydz, dx)


def convertOfficial25m(official, dx, out):
    """Provenance converter: raw SCEC 25 m file -> EQdyna format at dx.
    official: path to tpv29_tpv30_geometry_25m_data.txt."""
    data = np.loadtxt(official, skiprows=2)
    if data.shape != (1601 * 801, 7):
        raise ValueError(f'{official}: unexpected shape {data.shape}')
    # official row order: depth index (ny, 0 at surface) outer, nx inner
    F = data[:, 4].reshape(801, 1601)
    dFdx = data[:, 5].reshape(801, 1601)
    dFdy = data[:, 6].reshape(801, 1601)
    step = int(round(dx / 25.0))
    if abs(dx / 25.0 - step) > 1e-9:
        raise ValueError(f'dx={dx} is not a multiple of 25 m')
    # EQdyna iz counts up from the bottom (depth 20000) -> flip depth axis;
    # dy/dz_eq = -dF/dy_tpv (z_eq = -y_tpv).
    y = F[::-1][::step][:, ::step]
    dydx = dFdx[::-1][::step][:, ::step]
    dydz = -dFdy[::-1][::step][:, ::step]
    writeBFault(y, dydx, dydz, float(dx), out=out)
    return y, dydx, dydz


if __name__ == '__main__':
    from user_defined_params import par
    y, dydx, dydz = faultGridForCase(par.dx)
    writeBFault(y, dydx, dydz, par.dx)
    print(f'bFault_Rough_Geometry.txt written at dx={par.dx} m: '
          f'nx={y.shape[1]}, nz={y.shape[0]}, '
          f'y range [{y.min():.1f}, {y.max():.1f}] m')
