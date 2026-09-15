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
  1. loads an EQdyna-format geometry file. The compset SHIPS the official
     surface at two resolutions (SHIPPED_SOURCES), both exact decimations of
     the official 25 m file -- no interpolation, every value is an official
     one -- so a clean checkout runs both tiers without the 66 MB download:
         50 m  (the spec resolution, full tier)     14 MB
         100 m (spec-acceptable, and every coarser
                gate/trial dx decimates from it)   3.5 MB
  2. decimates the coarsest shipped source that can reach the case's mesh
     spacing par.dx (which must be a multiple of 50 m and divide the
     40 km x 20 km fault); and
  3. writes bFault_Rough_Geometry.txt in the format read by
     src/readInputFiles.f90:read_fault_rough_geometry /
     src/func_lib.f90:insertFaultInterface:
         row 0: nnx nnz 0
         row 1: dx fxmin fzmin
         rows:  [y, dy/dx, dy/dz], z-fastest (index nnz*(ix) + iz, iz counted
                upward from fzmin), x from fxmin.

Command line:

    # in a case directory: regenerate bFault_Rough_Geometry.txt for par.dx
    python3 tpv29GeometryTools.py

    # same, for an explicit dx (must be reachable from a shipped source)
    python3 tpv29GeometryTools.py --dx 100

    # from the OFFICIAL 25 m download: produce any dx it can reach
    # (25, 50, 100, 200, 500, ... m) -- this is how the shipped files were made
    python3 tpv29GeometryTools.py --official tpv29_tpv30_geometry_25m_data.txt \\
            --dx 50 --out bFault_Rough_Geometry.tpv29.50m.txt

Every path above validates what it wrote with
scripts/lib.py:validateFaultRoughGeometry before reporting success.

CONVENTION for any compset that supplies its own fault geometry (TPV30 uses
this identical surface): ship the resolutions the compset's own tiers need,
ship or reference the converter that derives them from the authoritative
source, declare the finest shipped spacing as par.faultGeometrySourceDx so
lib.requireFaultGeometryResolution can refuse an unreachable par.dx, and state
the provenance chain in the compset README.
"""
import argparse
import sys

import numpy as np

from lib import (FaultGeometryError,
                 requireFaultGeometryResolution,
                 validateFaultRoughGeometry,
                 writeFaultGeometryProvenance)

FXMIN, FXMAX = -20000.0, 20000.0
FZMIN, FZMAX = -20000.0, 0.0
YMIN, YMAX = -30000.0, 30000.0      # the case's model domain, for validation

OFFICIAL_DX = 25.0
OFFICIAL_SHAPE = (801, 1601)        # (depth nodes, strike nodes) at 25 m
OFFICIAL_URL = ('https://strike.scec.org/cvws/tpv29_30docs.html -> '
                'download/tpv29_tpv30_geometry_25m_data.zip')

# Shipped sources, FINEST FIRST. Both are exact decimations of the official
# 25 m data, so decimating either one to a common dx gives identical values;
# the 100 m file is kept because every gate/trial dx (100, 200, 500 m) reaches
# its grid from 3.5 MB instead of 14 MB.
SHIPPED_SOURCES = (
    (50.0,  'bFault_Rough_Geometry.tpv29.50m.txt'),
    (100.0, 'bFault_Rough_Geometry.tpv29.100m.txt'),
)
SHIPPED_DX = SHIPPED_SOURCES[0][0]   # finest shipped spacing = 50 m
SHIPPED_NAME = ('the official SCEC TPV29/30 surface shipped with this compset '
                '(50 m and 100 m)')


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
    # Shared verdict and shared wording for "this dx is not reachable from the
    # supplied surface" (scripts/lib.py), so every compset fails the same way.
    requireFaultGeometryResolution(dx, srcDx, SHIPPED_NAME)
    step = int(round(dx / srcDx))
    for extent in (FXMAX - FXMIN, FZMAX - FZMIN):
        if abs(extent / dx - round(extent / dx)) > 1e-9:
            raise ValueError(f'dx={dx} does not divide the fault extent {extent} m')
    return y[::step, ::step], dydx[::step, ::step], dydz[::step, ::step]


def _sourceLabel(dx):
    """Name and spacing of the shipped source a dx would come from, for the
    provenance sidecar; falls back to the official data for a dx no shipped
    source reaches (i.e. a convertOfficial25m run)."""
    try:
        srcDx, fname = shippedSourceForDx(dx)
        return f'{fname} (shipped with this compset)', srcDx
    except (ValueError, FaultGeometryError):
        return (f'the official SCEC 25 m data file '
                f'(tpv29_tpv30_geometry_25m_data.txt)'), OFFICIAL_DX


def shippedSourceForDx(dx):
    """The COARSEST shipped source that can reach dx by exact decimation.

    Coarsest, not finest: reaching 500 m from the 100 m file reads 3.5 MB
    instead of 14 MB, and the values are identical either way because both
    files are decimations of the same official 25 m data.
    """
    for srcDx, fname in reversed(SHIPPED_SOURCES):
        ratio = dx/srcDx
        if dx >= srcDx - 1e-9 and abs(ratio - round(ratio)) <= 1e-9:
            return srcDx, fname
    # Nothing shipped can do it -- report against the FINEST source, which is
    # the one that decides what is reachable at all.
    requireFaultGeometryResolution(dx, SHIPPED_DX, SHIPPED_NAME)
    raise ValueError(f'no shipped TPV29 source reaches dx={dx} m')


def writeBFault(y, dydx, dydz, dx, out='bFault_Rough_Geometry.txt',
                validate=True, verbose=False, source=None, sourceDx=None):
    """Write EQdyna's bFault_Rough_Geometry.txt (z-fastest ordering) plus its
    provenance sidecar, then check it with the shared validator so this tool
    cannot emit a file case.setup or EQdyna would reject."""
    nz, nx = y.shape
    with open(out, 'w') as f:
        f.write(f'{nx}\t{nz}\t0\n')
        f.write(f'{dx:.6f}\t{FXMIN:.6f}\t{FZMIN:.6f}\n')
        for ix in range(nx):
            for iz in range(nz):
                f.write(f'{y[iz, ix]:.7e}\t{dydx[iz, ix]:.7e}\t{dydz[iz, ix]:.7e}\n')
    if source is None:
        source, sourceDx = _sourceLabel(dx)
    writeFaultGeometryProvenance(out, dict(
        benchmark='SCEC TPV29/TPV30 (TPV29_30_Description_v06)',
        tool='case_input/test.tpv29/tpv29GeometryTools.py',
        source=source, sourceDx=f'{sourceDx:.12g}',
        method='exact decimation (no interpolation); every value is official',
        officialData=OFFICIAL_URL,
        dx=f'{dx:.12g}', nnx=str(nx), nnz=str(nz),
        fxmin=f'{FXMIN:.12g}', fzmin=f'{FZMIN:.12g}'))
    if validate:
        validateFaultRoughGeometry(out, dx=dx, dz=dx, dy=dx,
                                   fxmin=FXMIN, fxmax=FXMAX,
                                   fzmin=FZMIN, fzmax=FZMAX,
                                   ymin=YMIN, ymax=YMAX, verbose=verbose)
    return out


def faultGridForCase(dx):
    """Return (y, dydx, dydz) arrays shaped (nfz, nfx), [iz, ix] with iz=0 at
    fzmin, decimated from a shipped official surface to spacing dx."""
    # Checked before any source is read, so an unreachable dx costs nothing
    # and reports the reason rather than a decimation failure.
    requireFaultGeometryResolution(dx, SHIPPED_DX, SHIPPED_NAME)
    _, fname = shippedSourceForDx(dx)
    srcDx, y, dydx, dydz = loadEQdynaGeometry(fname)
    return decimateToDx(srcDx, y, dydx, dydz, dx)


def convertOfficial25m(official, dx, out, verbose=False):
    """Convert the OFFICIAL SCEC 25 m download to EQdyna format at spacing dx.

    This is how both shipped surfaces were produced, and the entry point a case
    author uses for a resolution neither shipped file can reach.

    official : path to tpv29_tpv30_geometry_25m_data.txt (from OFFICIAL_URL;
               1601 x 801 nodes x 7 columns, 66 MB unzipped)
    dx       : target spacing, m -- any integer multiple of 25 that divides the
               40 km x 20 km fault (25, 50, 100, 200, 250, 500, ...)
    out      : output path, EQdyna bFault_Rough_Geometry format

    Returns (y, dydx, dydz) shaped (nz, nx), [iz, ix], iz = 0 at fzmin.
    The written file is put through lib.validateFaultRoughGeometry before this
    returns; a file that fails is left on disk for inspection and the error is
    raised (PROJECT_RULES.md rule 2).
    """
    requireFaultGeometryResolution(dx, OFFICIAL_DX,
                                   f'the official SCEC 25 m data file ({official})')
    data = np.loadtxt(official, skiprows=2)
    ny, nx = OFFICIAL_SHAPE
    if data.shape != (nx*ny, 7):
        raise ValueError(
            f'{official}: shape {data.shape}, expected ({nx*ny}, 7). Is this '
            f'the official 25 m data file from {OFFICIAL_URL}?')
    # official row order: depth index (ny, 0 at surface) outer, nx inner
    F = data[:, 4].reshape(ny, nx)
    dFdx = data[:, 5].reshape(ny, nx)
    dFdy = data[:, 6].reshape(ny, nx)
    step = int(round(dx / OFFICIAL_DX))
    # EQdyna iz counts up from the bottom (depth 20000) -> flip depth axis;
    # dy/dz_eq = -dF/dy_tpv (z_eq = -y_tpv).
    y = F[::-1][::step][:, ::step]
    dydx = dFdx[::-1][::step][:, ::step]
    dydz = -dFdy[::-1][::step][:, ::step]
    writeBFault(y, dydx, dydz, float(dx), out=out, verbose=verbose,
                source=f'the official SCEC 25 m data file ({official})',
                sourceDx=OFFICIAL_DX)
    return y, dydx, dydz


def _main_func():
    p = argparse.ArgumentParser(
        description='TPV29/30 official rough-fault geometry: regenerate '
                    'bFault_Rough_Geometry.txt for a case, or derive a shipped '
                    'surface from the official 25 m SCEC download.',
        epilog=f'Official 25 m data: {OFFICIAL_URL}')
    p.add_argument('--dx', type=float, default=None,
                   help='target spacing, m (default: par.dx from '
                        'user_defined_params.py in the current directory)')
    p.add_argument('--official', default=None,
                   help='path to tpv29_tpv30_geometry_25m_data.txt; derive '
                        'from the official source instead of a shipped one')
    p.add_argument('-o', '--out', default='bFault_Rough_Geometry.txt')
    args = p.parse_args()

    dx = args.dx
    if dx is None:
        from user_defined_params import par
        dx = par.dx

    if args.official:
        y, _, _ = convertOfficial25m(args.official, dx, args.out, verbose=True)
        src = f'the official 25 m data ({args.official})'
    else:
        srcDx, fname = shippedSourceForDx(dx)
        y, dydx, dydz = faultGridForCase(dx)
        writeBFault(y, dydx, dydz, dx, out=args.out, verbose=True)
        src = f'the shipped {srcDx:g} m surface ({fname})'
    print(f'{args.out} written at dx={dx} m from {src}: '
          f'nx={y.shape[1]}, nz={y.shape[0]}, '
          f'y range [{y.min():.1f}, {y.max():.1f}] m')
    return 0


if __name__ == '__main__':
    sys.exit(_main_func())
