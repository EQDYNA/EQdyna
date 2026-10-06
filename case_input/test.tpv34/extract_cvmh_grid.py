#! /usr/bin/env python3
"""
Extract the SCEC TPV34 (Imperial Fault, Model 1) velocity structure from
CVM-H with vx_lite, exactly as TPV34_Description_v10 Part 2 prescribes, on
a uniform grid of EQdyna ELEMENT CENTRES, and write it as the shipped
`tpv34_cvmh_grid_<dx>m.txt.gz` that tpv34Tools.py turns into par.mat (the
n2mat == 6 3D material grid, rows `x y z vp vs rho`).

This is the ONE place CVM-H is queried. The solver never interpolates: each
element reads the grid cell nearest its centre (meshgen.f90
setElementMaterial / meshgen.py build_elements, n2mat == 6), so a uniform
belt element whose centre is a grid point reads its own CVM-H sample
exactly, and the fault's tau0/sigma0 (prop. to mu/mu0, spec Part 3) are
computed in tpv34Tools from the SAME samples.

Spec recipe (Part 2, "Algorithm for Obtaining Velocities and Densities"):
  Step 1  (x, y, z)_tpv: x along strike, y depth (positive down), z normal
          to the fault, positive on the far side.  EQdyna frame mapping
          (same proper rotation as tpv29/tpv35): x_eq = x, y_eq = z,
          z_eq = -y.
  Step 2  UTM zone 11 NAD27:  X = 648446 - 0.5802386 x - 0.8144465 z
                              Y = 3625237 + 0.8144465 x - 0.5802386 z
                              Z = max(y, 100 m)
  Step 4  vx_lite -s -z dep -m <model> < infile > outfile
  Step 5  columns 17, 18, 19 = Vp, Vs, rho (m/s, m/s, kg/m3)
  Step 6  Vp, Vs or rho <= 0 is an error -- refused here
  Step 7  if Vp <= 2984 or Vs <= 1400: Vp = 2984, Vs = 1400, rho = 2220.34

Provenance of what this script ran against (recorded in the output header):
  model   CVM-H data files as distributed by UCVM (SCECcode/cvmh main,
          model/config -> hypocenter data host), kept once read-only in
          ~/shared_dataset/scec_cvmh.15.1.1/ (MANIFEST.json has per-file
          md5).  The spec names CVM-H 15.1.0; UCVM ships 15.1.1.  The
          difference is labelled, not hidden: 15.1.1 is a packaging/bug-fix
          release of the same model; the hypocentre sample reproduces the
          2016 SCEC-submitted EQdyna tau0 to 3 digits (see README.md).
  code    vx_lite built from SCECcode/cvmh main (reports Version 11.9.0).

Usage (defaults are the test.tpv34 case box at the 500 m gate spacing):
  python3 extract_cvmh_grid.py --vx-lite <path>/vx_lite --model <path>/cvmh \\
      [--dx 500] [--box -30e3 30e3 -16e3 16e3 -26e3 0] [--out <file>]
Requires LD_LIBRARY_PATH to reach libvxapi.so if vx_lite is not installed.
"""
import argparse
import datetime
import gzip
import os
import subprocess
import sys

import numpy as np

# Spec Part 2, Step 2 and Step 7 constants (TPV34_Description_v10).
UTM_X0, UTM_Y0 = 648446.0, 3625237.0
COS_A, SIN_A = 0.5802386, 0.8144465
MIN_QUERY_DEPTH = 100.0
VP_MIN, VS_MIN, RHO_AT_MIN = 2984.0, 1400.0, 2220.34

# test.tpv34's model box (user_defined_params.py) and gate spacing.
DEFAULT_BOX = (-30.0e3, 30.0e3, -16.0e3, 16.0e3, -26.0e3, 0.0)
DEFAULT_DX = 500.0


def cell_centres(box, dx):
    """Uniform cell centres (EQdyna frame) tiling `box` at spacing dx;
    refuses a box that dx does not tile exactly."""
    axes = []
    for lo, hi in zip(box[0::2], box[1::2]):
        n = (hi - lo) / dx
        if abs(n - round(n)) > 1e-9:
            raise SystemExit('box extent %g..%g is not a multiple of dx=%g' % (lo, hi, dx))
        axes.append(lo + dx * (np.arange(int(round(n))) + 0.5))
    X, Y, Z = np.meshgrid(*axes, indexing='ij')
    return np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])


def to_utm(pts_eq):
    """EQdyna (x, y, z) -> vx_lite (Easting, Northing, Depth) per Step 1-2."""
    x_tpv = pts_eq[:, 0]
    z_tpv = pts_eq[:, 1]
    depth = -pts_eq[:, 2]
    X = UTM_X0 - COS_A * x_tpv - SIN_A * z_tpv
    Y = UTM_Y0 + SIN_A * x_tpv - COS_A * z_tpv
    Z = np.maximum(depth, MIN_QUERY_DEPTH)
    return np.column_stack([X, Y, Z])


def run_vx_lite(vx_lite, model, utm):
    infile = '\n'.join('%.3f %.3f %.3f' % tuple(r) for r in utm) + '\n'
    proc = subprocess.run([vx_lite, '-s', '-z', 'dep', '-m', model], input=infile,
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit('vx_lite failed (%d):\n%s' % (proc.returncode, proc.stderr[-2000:]))
    rows = [ln.split() for ln in proc.stdout.splitlines() if ln.strip()]
    if len(rows) != utm.shape[0] or any(len(r) != 19 for r in rows):
        raise SystemExit('vx_lite returned %d lines (expected %d, 19 columns each)'
                         % (len(rows), utm.shape[0]))
    out = np.array([[float(r[16]), float(r[17]), float(r[18])] for r in rows])
    return out, proc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--vx-lite', required=True)
    ap.add_argument('--model', required=True, help='CVM-H model directory (-m)')
    ap.add_argument('--dx', type=float, default=DEFAULT_DX)
    ap.add_argument('--box', type=float, nargs=6, default=DEFAULT_BOX,
                    metavar=('XMIN', 'XMAX', 'YMIN', 'YMAX', 'ZMIN', 'ZMAX'))
    ap.add_argument('--out', default=None)
    ap.add_argument('--model-label', default='CVM-H 15.1.1 (UCVM distribution)')
    a = ap.parse_args(argv)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'tpv34_cvmh_grid_%dm.txt.gz' % int(round(a.dx)))

    pts = cell_centres(a.box, a.dx)
    utm = to_utm(pts)
    raw, proc = run_vx_lite(a.vx_lite, a.model, utm)
    bad = np.nonzero((raw <= 0.0).any(axis=1))[0]
    if bad.size:                                   # Step 6
        raise SystemExit('vx_lite returned a non-positive Vp/Vs/rho at %d points, first '
                         'at EQdyna %r -> UTM %r: %r' % (bad.size, pts[bad[0]].tolist(),
                                                          utm[bad[0]].tolist(), raw[bad[0]].tolist()))
    clamp = (raw[:, 0] <= VP_MIN) | (raw[:, 1] <= VS_MIN)   # Step 7
    props = raw.copy()
    props[clamp] = [VP_MIN, VS_MIN, RHO_AT_MIN]

    counts = [int(round((hi - lo) / a.dx)) for lo, hi in zip(a.box[0::2], a.box[1::2])]
    header = [
        '# SCEC TPV34 velocity structure, extracted from CVM-H by extract_cvmh_grid.py',
        '# spec: TPV34_Description_v10 Part 2 (UTM zone 11 NAD27 formula, Z = max(depth, 100 m),',
        '#       vx_lite -s -z dep, Vp/Vs/rho = columns 17-19, min-velocity clamp Step 7)',
        '# model: %s; spec names CVM-H 15.1.0 -- difference labelled, see README.md' % a.model_label,
        '# vx_lite: %s' % ' '.join(proc.args),
        '# extracted: %s UTC' % datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M'),
        '# grid: uniform ELEMENT-CENTRE grid, EQdyna frame (x along strike, y fault-normal '
        '(+ = spec far side), z up, m); dx = %g m; box x[%g,%g] y[%g,%g] z[%g,%g]' % ((a.dx,) + tuple(a.box)),
        '# counts: nx ny nz = %d %d %d (%d rows); clamped to Vp=2984/Vs=1400/rho=2220.34 at %d rows'
        % (counts[0], counts[1], counts[2], pts.shape[0], int(clamp.sum())),
        '# columns: x y z vp vs rho   (m m m m/s m/s kg/m3); row order x-major, then y, then z',
    ]
    with gzip.open(out, 'wt') as f:
        f.write('\n'.join(header) + '\n')
        for p, q in zip(pts, props):
            f.write('%.1f %.1f %.1f %.2f %.2f %.2f\n' % (p[0], p[1], p[2], q[0], q[1], q[2]))
    print('wrote %s: %d rows, %d clamped (%.1f%%), Vs range %.0f..%.0f m/s'
          % (out, pts.shape[0], int(clamp.sum()), 100.0 * clamp.mean(),
             props[:, 1].min(), props[:, 1].max()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
