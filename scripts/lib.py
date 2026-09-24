#! /usr/bin/env python3
from math import *
import sys
import glob
import re
import os
from os.path import exists

import numpy as np

# functions are defined in lib.py under scripts/
# function lists:
# - shear_steady_state
# - state_steady_state
# - B1, defined in TPV104 and TPV105
# - B2 and B3, defined in TPV105
# - resolveViscoplasticParams, case.setup's resolver+validator for the
#   viscoplastic/plastic-output block at the end of bGlobal.txt
# - shearModulusFromPar, mu = rho*Vs^2 from the case's own material config
# - loadFrtData, shared frt.txt* loader for plotRuptureDynamics/plotSlipAndRPT
# - tryint, alphanum_key, sort_nicely, generate_gif, seek_numbers_filename,
#   shared filename-sorting/gif helpers for plot_on_fault_vars
# - bFault_Rough_Geometry.txt support (see the block comment further down):
#   readFaultRoughGeometryHeader / readFaultRoughGeometry (reader),
#   validateFaultRoughGeometry / validateFaultRoughGeometryForCase (checks),
#   ensureFaultRoughGeometryForCase (case.setup's idempotent keep-or-generate),
#   requireFaultGeometryResolution (a dx the geometry source cannot supply),
#   writeFaultRoughGeometry (writer),
#   writeFaultGeometryProvenance / readFaultGeometryProvenance (sidecar)

def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) #+ rou*vs/2.0*slip_rate
  return res

def state_steady_state(a,b,d0,v0,r0,shear,norm,slip_rate, friclaw):
  # calculate the state variable at steady state
    if friclaw == 3:
        tmp   = a*log(2.*sinh(abs(shear/norm/a))) - r0 - log(slip_rate/v0)
        state = d0/v0*exp(tmp/b)
    elif friclaw == 4 or friclaw == 5:
        state = a*log(2.*v0/slip_rate*sinh(abs(shear/norm/a)))
    return state

# Mathematically smooth version of the boxcar function B1, B2, and B3
def B1(x,ww,w):
  if abs(x)<=ww:
    res = 1.0
  elif abs(x)>ww and abs(x)<ww+w: 
    res = 0.5*(1. + tanh(w/(abs(x)-ww-w) + w/(abs(x)-ww)))
  elif abs(x)>=ww+w:
    res = 0.0
  return res

def B2(y,ww,w):
  if y<0:
    print('z coordinates should be positive for B3')
    sys.exit()
  if y<w:
    res = 0.5*(1.+tanh(w/(w-y) - w/(y+1e-30)))
  elif y>=w and y<=ww:
    res = 1.0
  elif y>ww and y<ww+w: 
    res = 0.5*(1. + tanh(w/(y-ww-w) + w/(y-ww)))
  elif y>=ww+w:
    res = 0.0
  return res
def B3(y,ww,w):
  if y<0:
    print('z coordinates should be positive for B3')
    sys.exit()
  if y<=ww:
    res = 1.
  elif y>ww and y<ww+w:
    res = 0.5*(1. + tanh(w/(y-ww-w) + w/(y-ww)))
  elif y>=ww+w:
    res = 0.
  return res

def linear1(x,ww,w):
  if abs(x)<=ww:
    res = 1.0
  elif abs(x)>ww and abs(x)<ww+w:
    res = 1. - abs((abs(x)-ww))/w
  elif abs(x)>=ww+w:
    res = 0.0
  return res

# globalvar.f90's own constant, m/s. Used here for ONE thing: reproducing the
# pre-v5.9.0 derivation of the viscoplastic relaxation time, 2*dz/NUC_VS_FIXED,
# which readInputFiles.f90 applied to every case when Tv had no input slot.
NUC_VS_FIXED = 3464.0


def resolveNormalStressSign(par):
    """The integer case.setup writes as bGlobal.txt's last line: +1 when the
    case's station files report n-stress positive in extension, -1 when
    positive in compression (par.faultStNormalStressSign; board row 22a).
    Anything else RAISES (rule 2) -- the sign is the case's spec convention and
    there is no neutral value to fall back to."""
    signs = {'extension': 1, 'compression': -1}
    value = par.faultStNormalStressSign
    if value not in signs:
        raise ValueError(
            'case.setup: par.faultStNormalStressSign must be %s (got %r). It '
            'is the n-stress sign convention the case SCEC spec states for '
            'faultst*.txt column 8.' % (' or '.join(repr(k) for k in signs), value))
    print('case.setup: station n-stress sign = %+d (positive means %s)'
          % (signs[value], value))
    return signs[value]


def resolveViscoplasticParams(par):
    """Resolve the viscoplastic/plastic-output block case.setup writes to the
    end of bGlobal.txt, and print the values it resolved.

    Returns (tv, taperStart, taperEnd, halfWidths):
      tv         -- viscoplastic (Duvaut-Lions) relaxation time, s.
                    par.viscoplasticRelaxTime, or 2*par.dz/3464 when that is
                    None -- the formula readInputFiles.f90 hardcoded until
                    v5.9.0, so an unset case is unchanged bit-for-bit (Python
                    and Fortran evaluate the same two IEEE-754 operations on
                    the same dz, and str() round-trips the double exactly).
      taperStart,
      taperEnd   -- depths (m, positive down) of the deviatoric pre-stress
                    taper, SCEC TPV29/30's Omega(depth). (0.0, 0.0) means no
                    taper, which func_lib.f90's devStrDepthTaper turns into an
                    exact 1.0 multiplier.
      halfWidths -- (|x|,|y|,|z|) half-widths (m) of the plastic-strain output
                    window.

    Every invalid combination RAISES (PROJECT_RULES.md rule 2) -- a
    half-configured taper, a non-positive Tv, or a non-positive window is a
    case-configuration error, not something to fill in with a guess.
    """
    tv = getattr(par, 'viscoplasticRelaxTime', None)
    if tv is None:
        tv = 2.0*par.dz/NUC_VS_FIXED
        tvSource = 'derived as 2*dz/%g (the pre-v5.9.0 hardcode)' % NUC_VS_FIXED
    else:
        tv = float(tv)
        tvSource = 'par.viscoplasticRelaxTime'
    if not tv > 0.0:
        raise ValueError(
            'case.setup: par.viscoplasticRelaxTime must be a positive time in '
            'seconds (got %r). It is Tv in exp(-dt/Tv) (calcElemKU.f90).' % tv)

    taperStart = getattr(par, 'devStrTaperDepthStart', None)
    taperEnd   = getattr(par, 'devStrTaperDepthEnd', None)
    if (taperStart is None) != (taperEnd is None):
        raise ValueError(
            'case.setup: par.devStrTaperDepthStart and '
            'par.devStrTaperDepthEnd must be set together (got %r and %r).\n'
            '  They are the two depths of SCEC TPV29/30\'s Omega(depth) taper '
            'on the off-fault deviatoric pre-stress; one of them alone does '
            'not define a taper, so this refuses rather than inventing the '
            'other.' % (taperStart, taperEnd))
    if taperStart is None:
        taperStart, taperEnd = 0.0, 0.0
        taperSource = 'no taper (devStr is a fixed fraction of |strVert| at every depth)'
    else:
        taperStart, taperEnd = float(taperStart), float(taperEnd)
        if taperStart < 0.0 or taperEnd <= taperStart:
            raise ValueError(
                'case.setup: the deviatoric-stress taper needs '
                '0 <= par.devStrTaperDepthStart < par.devStrTaperDepthEnd '
                '(depths in m, positive down); got %r and %r.'
                % (taperStart, taperEnd))
        taperSource = 'par.devStrTaperDepthStart/End'

    halfWidths = getattr(par, 'plasticOutputHalfWidth', None)
    if halfWidths is None or len(halfWidths) != 3:
        raise ValueError(
            'case.setup: par.plasticOutputHalfWidth must be three half-widths '
            '(|x|,|y|,|z|) in m for the plastic-strain output window (got %r).'
            % (halfWidths,))
    halfWidths = tuple(float(v) for v in halfWidths)
    if min(halfWidths) <= 0.0:
        raise ValueError(
            'case.setup: every par.plasticOutputHalfWidth entry must be '
            'positive (got %r); a non-positive half-width writes no '
            'plastic-strain output at all.' % (halfWidths,))

    print('VISCOPLASTIC: Tv = %r s, %s' % (tv, tvSource))
    print('VISCOPLASTIC: deviatoric-stress depth taper %r -> %r m, %s'
          % (taperStart, taperEnd, taperSource))
    print('VISCOPLASTIC: plastic-strain output window half-widths %r m' % (halfWidths,))
    return tv, taperStart, taperEnd, halfWidths


def shearModulusFromPar(par, depth=0.0):
    """Shear modulus mu = rho*Vs^2 (Pa) of THIS case's configured material, at
    `depth` (m, positive down).

    par.nmat == 1: the single (vp, vs, rou) block case.setup writes to
    bMaterial.txt, so depth is irrelevant.
    par.nmat  > 1: the layered table par.mat, whose rows are
    [layer bottom depth, vp, vs, rou] -- selected with the SAME comparisons
    meshgen.f90:187-198 uses to assign a material to an element
    (`abs(z) < material(1,1)` for the top layer, then
    `material(i-1,1) <= abs(z) < material(i,1)`), so the modulus used here is
    the one the solver actually used there.

    Raises if the depth lies below the deepest layer -- the same condition the
    Fortran refuses with ERR_MESH_MATERIAL_UNSET. There is no fallback
    constant (PROJECT_RULES.md rule 2): a wrong modulus silently rescales
    every reported seismic moment.
    """
    depth = abs(depth)
    if par.nmat == 1:
        return par.rou*par.vs**2
    mat = par.mat
    if depth < mat[0, 0]:
        return mat[0, 3]*mat[0, 2]**2
    for i in range(1, par.nmat):
        if mat[i-1, 0] <= depth < mat[i, 0]:
            return mat[i, 3]*mat[i, 2]**2
    raise ValueError(
        'shearModulusFromPar: depth %g m lies below the deepest layer in '
        'par.mat (bottom at %g m), so this case defines no material there -- '
        'the same state meshgen.f90 refuses with ERR_MESH_MATERIAL_UNSET.'
        % (depth, mat[par.nmat-1, 0]))


def loadFrtData(par):
    """Load and grid the on-fault frt.txt* output written by EQdyna.

    Shared by plotRuptureDynamics and plotSlipAndRPT: loops over
    frt.txt{me} for me in range(nx*ny*nz), grids each row onto the
    (dip, strike) fault mesh, and accumulates seismic moment.

    frt.txt* file structure (1-indexed as in the on-fault output; see
    the FRIC_SLOT_* slot map in src/globalvar.f90 for the underlying
    per-node fault-variable layout):
      1-3,   coorx,y,z
      4,     rupture time
      5-9,   final slips,d,n, final sliprates,d.
      10,    peak slip rate.
      11,    final sliprate magnitude
      12-14, final tnrm, tstk, tdip
      15-20, vxm,vym,vzm, vxs,vys,vzs.
      21,    state variable
      22,    state var for normal stress variation (Shi and Day)

    Returns (xx, zz, rupt, rupt2d, fVarArr, magnitude):
      xx, zz    -- meshgrid of along-strike/along-dip coordinates, km
      rupt      -- (na*ma, 3) flat [xcoor, along-dip distance, rupture time]
      rupt2d    -- (ma, na, 100) gridded rupture-time/slip/stress panels
      fVarArr   -- (ma, na, 100) gridded fault-restart variables, passed
                   to generateNcRestart(faultVarArr)
      magnitude -- moment magnitude computed from summed slip*area*shearMod
    """
    nprocs = par.nx*par.ny*par.nz
    na     = round((par.fxmax-par.fxmin)/par.dx+1)
    ma     = round((par.fzmax-par.fzmin)/par.dz+1)
    rupt   = np.zeros((na*ma,3))
    rupt2d = np.zeros((ma,na,100))
    fVarArr= np.zeros((ma,na,100))

    [xx,zz] = np.meshgrid(par.fx,par.fz/sin(par.dip/180.*pi))
    xx = xx/1.e3
    zz = zz/1.e3#/sin(par.dip/180.*pi) # along dip distance
    moment = 0.

    for me in range(nprocs):
      fname = 'frt.txt' + str(me)
      if exists(fname):
        print('Post-processing ' + fname + ' ... ...')
        a = np.loadtxt(fname)
        n, m = a.shape
        for i in range(n):
            #!! use round() instead of int()!!
            ii = round((a[i,0] - par.fxmin)/par.dx)
            jj = round((a[i,2] - par.fzmin)/par.dz)

            rupt[jj*na+ii,0] = a[i,0]  # xcoor
            rupt[jj*na+ii,1] = -a[i,2]/sin(par.dip/180.*pi) # zcoor to along dip distance, reverse sign to positive numbers.
            rupt[jj*na+ii,2] = a[i,3]  # rupture time

            rupt2d[jj,ii,0]  = a[i,3]  # rupture time
            rupt2d[jj,ii,1]  = (a[i,4]**2 + a[i,5]**2)**0.5  # slip magnitude
            rupt2d[jj,ii,2]  = a[i,9]                        # peak slip rate
            rupt2d[jj,ii,3]  = a[i,10]                       # final slip rate
            # mu = rho*Vs^2 of THIS case's material at this node's depth, not
            # the 3464^2*2800 constant this line carried until v5.9.0 -- that
            # was a density this repo's own default case does not use (2670),
            # and it rescales every reported moment/Mw.
            shearMod = shearModulusFromPar(par, abs(a[i,2]))
            moment = moment + rupt2d[jj,ii,1]*par.dx*par.dx*shearMod
            rupt2d[jj,ii,4]  = a[i,12]/1.e6 # final shear stress
            rupt2d[jj,ii,5]  = a[i,11]/1.e6 # final normal stress
            rupt2d[jj,ii,6]  = a[i,13]/1.e6 # final dip shear
            rupt2d[jj,ii,7]  = a[i,4] # final slip s
            rupt2d[jj,ii,8]  = a[i,5] # final slip d

            #
            # fVarArr will be passed to the function generateNcRestart(faultVarArr):
            fVarArr[jj,ii,0]  = a[i,12] # shear_strike, Pa
            fVarArr[jj,ii,1]  = a[i,13] # shear_dip, Pa
            fVarArr[jj,ii,2]  = a[i,11] # effective_normal, Pa
            fVarArr[jj,ii,3]  = a[i,10] # slip_rate, m/s
            fVarArr[jj,ii,4]  = a[i,20] # state_variable
            fVarArr[jj,ii,5]  = a[i,21] # state_normal
            fVarArr[jj,ii,6]  = a[i,14] # vxm, m/s
            fVarArr[jj,ii,7]  = a[i,15] # vym
            fVarArr[jj,ii,8]  = a[i,16] # vzm
            fVarArr[jj,ii,9]  = a[i,17] # vxs
            fVarArr[jj,ii,10] = a[i,18] # vys
            fVarArr[jj,ii,11] = a[i,19] # vzs

    magnitude = 2/3*log10(moment*1.e7)-10.7

    return xx, zz, rupt, rupt2d, fVarArr, magnitude

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    return [tryint(c) for c in re.split('([0-9]+)', s)]

def sort_nicely(l):
    l.sort(key=alphanum_key, reverse=False)
    return l

def generate_gif():
  # Lazy AND guarded: imageio is optional and not installed in CI. Every
  # other plot is already written by the time this runs, so the GIF is the
  # only thing lost -- say that instead of raising a bare ImportError.
  try:
      import imageio
  except ImportError:
      raise RuntimeError(
          'writing the animated GIF needs imageio, which is not installed. '
          'Install imageio, or just use the PNG frames already written in '
          'this directory -- no other output is affected.')
  filenames = glob.glob('.//*.png')
  filenames = sort_nicely(filenames)
  with imageio.get_writer('./on_fault_vars.gif', mode='I') as writer:
      for filename in filenames:
          image = imageio.v2.imread(filename)
          writer.append_data(image)

def seek_numbers_filename(x):
    return (x[6:10])

# ---------------------------------------------------------------------------
# bFault_Rough_Geometry.txt: reader, validator, writer
# ---------------------------------------------------------------------------
# File format, as consumed by src/readInputFiles.f90:read_fault_rough_geometry
# and src/func_lib.f90:insertFaultInterface:
#
#   row 0:  nnx  nnz  (third field ignored)
#   row 1:  dx   fxmin  fzmin
#   rows 2 .. nnx*nnz+1:  y  dy/dx  dy/dz     <- z-fastest within each x column
#
# The Fortran reads the header, allocates rough_geo(3, nnx*nnz), then reads
# exactly nnx*nnz rows.  insertFaultInterface then indexes that array with
#     ixx = nint((x - fxmin)/dx) + 1        (mesh dx, NOT the header's dx)
#     izz = nint((z - fzmin)/dz) + 1        (mesh dz)
#     rough_geo(:, nnz*(ixx-1) + izz)
# so the file's grid must be *exactly* the case's fault grid: same node counts,
# same spacing, same corner.  A file that disagrees either walks off the end of
# the array or silently morphs the mesh onto the wrong surface.  Nothing in the
# format is self-describing enough for the Fortran to notice on its own, which
# is what these checks exist for (PROJECT_RULES.md rule 2).

FAULT_ROUGH_GEOMETRY_FILE = 'bFault_Rough_Geometry.txt'

# Coordinates (m) must agree to 1 mm: far below any cell size EQdyna runs, and
# far above the 1e-6 m quantization of the '%f' writer in generateFaultInterface.
FAULT_GEOM_COORD_TOL = 1.0e-3

# Derivative-consistency safety factor; see _checkDerivativeColumn for the
# measurement behind it.
FAULT_GEOM_DERIV_SAFETY = 4.0
# Absolute floor, so a planar (all-zero) surface and the 1e-6 quantization of
# the '%f' writer do not trip the check.
FAULT_GEOM_DERIV_FLOOR = 1.0e-6

# WARN at 0.2 on the raw |dy/dx|, |dy/dz|: the mesh-distortion threshold
# scripts/generateFaultInterface has always warned at.  All three gated cases
# that insert a fault interface sit above it -- test.drv.a6's fractal 0.38,
# test.tpv29's official SCEC surface 0.65, and test.tpv10's 60-degree dipping
# plane 0.577 (= cot 60) -- so it is a note, never a failure.
FAULT_GEOM_SLOPE_WARN = 0.2
# A supplied surface whose steepest slope is below this is suspiciously flat.
# No check derived from the file alone can see a GLOBAL rescale -- scale y and
# both derivative columns by the same k and every self-consistency test is
# invariant (verified: k = 1e-3 passed everything). That is exactly the shape
# of a units error, a surface supplied in km where the code wants m. The domain
# check catches the k > 1 direction (the surface leaves the model), so this
# catches the other one. It is a WARNING, not a failure: a genuinely
# near-planar supplied surface is legal, it is just unusual enough to say so.
# For reference the official TPV29 surface has max |slope| = 0.65.
FAULT_GEOM_SLOPE_FLAT_WARN = 1.0e-3

# FAIL is a hard geometric limit on the per-cell y OFFSET, not on the slope.
# insertFaultInterface puts the fault node of column ix at y = peak(ix), while
# the nearest off-fault node layer sits a cell away at y ~ peak(ix) + dy.  The
# surface climbs by |dy/dx|*dx from one x column to the next and by
# |dy/dz|*dz from one z row to the next, so the elements tangle once
#     |dy/dx|*dx/dy >= 1   or   |dy/dz|*dz/dy >= 1.
# The dx/dy and dz/dy factors are load-bearing, not decoration: a dipping
# planar fault has dz = dx*sin(dip) and |dy/dz| = cot(dip), which passes 1 at a
# 45-degree dip -- while its actual per-row offset is dx*cos(dip) < dy = dx and
# nothing tangles at any dip.  Testing the raw slope would refuse every fault
# dipping at 45 degrees or less.
FAULT_GEOM_OFFSET_FAIL = 1.0


class FaultGeometryError(ValueError):
    """bFault_Rough_Geometry.txt does not describe the case's fault grid."""


def requireFaultGeometryResolution(dx, sourceDx, sourceName='the supplied '
                                   'fault-geometry source', availableDx=None):
    """Fail loudly when a case asks for a dx its geometry source cannot give.

    A supplied fault surface is a FIXED sampling of a specific random
    realisation.  It can be coarsened exactly (take every n-th node -- still
    the official values), but it cannot be refined: interpolating a fractal
    surface invents roughness that is not in the benchmark, so the only honest
    answer to "I need a finer grid" is a finer source.  Both halves of that --
    "dx must be a multiple of sourceDx" and "dx must not be finer than
    sourceDx" -- are the same check, and it lives here rather than in any one
    compset's tooling so every case supplying geometry gets the same verdict
    and the same message.

    Prompted by a real defect: test.tpv29 shipped only a 100 m surface while
    its FULL_SPECS entry asks for the 50 m spec resolution, so the full tier
    died inside the compset's own decimator with "dx=50.0 is not an integer
    multiple of the shipped geometry spacing 100.0 m" -- a true statement that
    says nothing about what to do (pathway_forward item 31).
    """
    if sourceDx is None:
        return
    sourceDx, dx = float(sourceDx), float(dx)
    have = ('' if not availableDx else
            ' Sources available here: '
            + ', '.join(f'{float(v):g} m' for v in availableDx) + '.')
    if sourceDx <= 0.0:
        raise FaultGeometryError(
            f'{sourceName} declares a non-positive native spacing {sourceDx}.')
    ratio = dx/sourceDx
    if dx < sourceDx - 1.0e-9:
        raise FaultGeometryError(
            f'this case requests dx = {dx} m, finer than {sourceName}, which '
            f'is sampled at {sourceDx} m. A supplied fault surface can be '
            f'coarsened exactly but never refined -- interpolating it would '
            f'invent roughness the benchmark does not have. Supply a source '
            f'at {dx} m or finer, or set par.dx to a multiple of {sourceDx} m '
            f'(e.g. {sourceDx}, {2*sourceDx}, {5*sourceDx} m).{have}')
    if abs(ratio - round(ratio)) > 1.0e-9:
        raise FaultGeometryError(
            f'this case requests dx = {dx} m, which is not an integer multiple '
            f'of the {sourceDx} m spacing of {sourceName} ({dx}/{sourceDx} = '
            f'{ratio}). Exact decimation needs a whole-number stride, and '
            f'interpolating a supplied rough surface would change the '
            f'benchmark. Use a multiple of {sourceDx} m '
            f'(e.g. {int(round(ratio))*sourceDx} or '
            f'{(int(round(ratio))+1)*sourceDx} m).')


def _faultGeomClose(a, b):
    return abs(a - b) <= FAULT_GEOM_COORD_TOL + 1.0e-9*abs(b)


def _faultGeomNodeCount(lo, hi, h, axis, problems):
    """Node count of the fault grid along one axis, the way defaultParameters.py
    computes par.nfx/par.nfz (round((hi-lo)/h + 1)) and the way
    insertFaultInterface indexes rough_geo (nint((coord-lo)/h) + 1).  The two
    agree only when (hi-lo)/h is an integer, so say so when it is not."""
    if h <= 0.0:
        problems.append(f'{axis} spacing must be positive, got {h}')
        return None
    span = (hi - lo)/h
    if abs(span - round(span)) > 1.0e-6:
        problems.append(
            f'the fault does not land on a whole number of cells along {axis}: '
            f'({hi} - {lo})/{h} = {span}, not an integer. insertFaultInterface '
            f'indexes the rough grid with nint((coord - min)/spacing), so the '
            f'fault edge would fall between rough-geometry nodes.')
        return None
    return int(round(span)) + 1


def readFaultRoughGeometryHeader(fname=FAULT_ROUGH_GEOMETRY_FILE):
    """Read just the two header rows: dict(nnx, nnz, dx, fxmin, fzmin).

    Separate from readFaultRoughGeometry so validateFaultRoughGeometry can
    report a header that disagrees with the case BEFORE the row count it
    implies -- otherwise a wrong nnx surfaces only as "wrong number of rows",
    which points at the wrong end of the problem.
    """
    if not exists(fname):
        raise FaultGeometryError(
            f'{fname} is required for insertFaultType > 0 but does not exist. '
            f'Cases with insertFaultType 1 or 2 get it from '
            f'scripts/generateFaultInterface at case.setup time; with '
            f'insertFaultType 3 the case itself must supply it.')

    with open(fname) as f:
        head = [f.readline() for _ in range(2)]
    if not head[1].strip():
        raise FaultGeometryError(
            f'{fname}: expected two header rows ("nnx nnz" then '
            f'"dx fxmin fzmin") before the data; the file has fewer than two '
            f'non-empty lines.')
    try:
        fields0 = [float(v) for v in head[0].split()]
        fields1 = [float(v) for v in head[1].split()]
    except ValueError as exc:
        raise FaultGeometryError(f'{fname}: header rows are not numeric: {exc}')
    if len(fields0) < 2 or len(fields1) < 3:
        raise FaultGeometryError(
            f'{fname}: header row 0 needs "nnx nnz" (got {len(fields0)} '
            f'field(s)) and row 1 needs "dx fxmin fzmin" (got {len(fields1)}).')

    nnxTmp, nnzTmp = fields0[0], fields0[1]
    for name, val in (('nnx', nnxTmp), ('nnz', nnzTmp)):
        if abs(val - round(val)) > 1.0e-6 or round(val) < 2:
            raise FaultGeometryError(
                f'{fname}: header {name} = {val} is not an integer >= 2. '
                f'The Fortran does nint() on it and allocates nnx*nnz rows.')
    nnx, nnz = int(round(nnxTmp)), int(round(nnzTmp))
    return dict(nnx=nnx, nnz=nnz,
                dx=fields1[0], fxmin=fields1[1], fzmin=fields1[2])


def readFaultRoughGeometry(fname=FAULT_ROUGH_GEOMETRY_FILE, header=None):
    """Read a bFault_Rough_Geometry.txt-format file.

    Returns (header, y, dydx, dydz) where header is a dict with keys
    nnx, nnz, dx, fxmin, fzmin and the three arrays are shaped (nnz, nnx),
    indexed [iz, ix] with iz = 0 at fzmin and ix = 0 at fxmin -- the same
    [iz, ix] convention as par.on_fault_vars.

    Raises FaultGeometryError for anything that stops the file being read as
    that format at all (missing, truncated header, ragged rows, wrong row
    count).  Value-level checks live in validateFaultRoughGeometry.
    """
    if header is None:
        header = readFaultRoughGeometryHeader(fname)
    nnx, nnz = header['nnx'], header['nnz']

    try:
        data = np.loadtxt(fname, skiprows=2, ndmin=2)
    except ValueError as exc:
        raise FaultGeometryError(
            f'{fname}: data rows are not a clean N-by-3 table of numbers '
            f'("y dy/dx dy/dz" per row): {exc}')
    if data.ndim != 2 or data.shape[1] != 3:
        raise FaultGeometryError(
            f'{fname}: data rows have {data.shape[1] if data.ndim == 2 else "?"} '
            f'column(s), expected exactly 3 ("y dy/dx dy/dz").')
    if data.shape[0] != nnx*nnz:
        raise FaultGeometryError(
            f'{fname}: {data.shape[0] + 2} rows total ({data.shape[0]} data '
            f'rows), expected {nnx*nnz + 2} (2 header + nnx*nnz = {nnx}*{nnz} '
            f'= {nnx*nnz}). read_fault_rough_geometry reads exactly nnx*nnz '
            f'rows after the header, so a short file reads past the end of the '
            f'data and a long one silently ignores the tail.')

    # file order is z-fastest within each x column (rough_geo index nnz*(ix)+iz)
    y    = data[:, 0].reshape(nnx, nnz).T.copy()
    dydx = data[:, 1].reshape(nnx, nnz).T.copy()
    dydz = data[:, 2].reshape(nnx, nnz).T.copy()
    return header, y, dydx, dydz


def _checkDerivativeColumn(y, deriv, h, axis, label, problems,
                           derivSafety=FAULT_GEOM_DERIV_SAFETY):
    """Check that a derivative column is the derivative of the surface column.

    The reference is np.gradient(y, h, axis) -- second-order central in the
    interior, first-order one-sided on the two boundary lines -- which is
    exactly what scripts/generateFaultInterface writes, so generated files
    agree to their own write precision.

    A supplied file may instead carry the ANALYTIC derivative (SCEC's TPV29
    file does), and then the residual is just the stencil's truncation error,
    which is bounded by a quantity the file itself supplies:

        interior:  central(y) - y'  =  (h^2/6) y'''  ~  (1/6) * D2(y')
        boundary:  onesided(y) - y' = -(h/2) y''     ~  (1/2) * D1(y')

    where D2/D1 are the second/first differences of the DERIVATIVE column.
    Measured on the official SCEC TPV29 surface, max residual / max predicted
    bound is 1.00 at 25 m (its native sampling), 1.00 at 100 m, 1.01 at 200 m,
    1.09 at 500 m and 1.27-1.47 at 1000 m for the interior test, and 0.98-1.10
    for the boundary test over the same range.  FAULT_GEOM_DERIV_SAFETY = 4.0
    therefore clears every measured case by ~3x while still
    catching, say, a derivative column scaled by 2 (residual ~ max|y'|, which
    is 0.66 for TPV29, against a bound of 4*0.055 = 0.22 at 500 m).

    Comparing maxima rather than point by point is deliberate: the predicted
    bound is a field-level quantity, and a point-by-point test would divide by
    a locally vanishing curvature.
    """
    n = y.shape[axis]
    if n < 2:
        return
    g = np.gradient(y, h, axis=axis)
    resid = np.abs(deriv - g)
    take = (lambda a, sl: a[:, sl]) if axis == 1 else (lambda a, sl: a[sl])

    # The residual is compared to its bound ELEMENTWISE, not max-to-max.
    # Aggregating both sides with .max() first was a demonstrated false-accept:
    # the bound is then set by the roughest cell in the whole field, so an
    # error 3000x the legitimate local tolerance at the flattest cell passed
    # unnoticed (verified on an analytic surface, 2026-09-14). The floor keeps
    # a genuinely flat region from demanding an impossible zero tolerance.
    def _neighbourhoodMax(arr, axes):
        """Max of arr over +/-1 along `axes`, edge-replicated at the borders.

        The pointwise truncation bound (h^2/6)|y'''| collapses to ~0 wherever
        y''' has a zero crossing, while the actual residual there is set by the
        next term in the expansion. On the real TPV29 100 m surface that made
        the worst pointwise residual/bound ratio 162, against 1.05 once the
        bound is taken over a 3x3 neighbourhood -- so the smoothing is what
        makes an elementwise test possible at all. It stays local: the bound
        grows only to its immediate neighbours, not to the field maximum.
        """
        out = arr
        for ax in axes:
            if arr.shape[ax] < 3:
                continue
            pad = [(0, 0)]*arr.ndim
            pad[ax] = (1, 1)
            P = np.pad(out, pad, mode='edge')
            sl = [slice(None)]*arr.ndim
            acc = None
            for k in range(3):
                sl[ax] = slice(k, k + arr.shape[ax])
                v = P[tuple(sl)]
                acc = v if acc is None else np.maximum(acc, v)
            out = acc
        return out

    def _flag(residArr, boundArr, stencil, where, axes):
        boundArr = _neighbourhoodMax(boundArr, axes)
        boundArr = np.maximum(FAULT_GEOM_DERIV_FLOOR, derivSafety*boundArr)
        bad = residArr > boundArr
        if not bad.any():
            return
        ratio = residArr/boundArr
        k = np.unravel_index(np.argmax(ratio), ratio.shape)
        problems.append(
            f'{label} column disagrees with a {stencil} of the surface column '
            f'{where}: {int(bad.sum())} of {bad.size} nodes exceed the '
            f'truncation-error bound the derivative column itself implies '
            f'({derivSafety}x). Worst at interior index {tuple(int(v) for v in k)}: '
            f'|d - {stencil.split()[0]}(y)| = {residArr[k]:.4e} vs bound '
            f'{boundArr[k]:.4e} ({ratio[k]:.1f}x over). The derivative columns '
            f'look stale, mis-scaled, or locally corrupted relative to the '
            f'surface column.')

    # boundary lines: first-order one-sided stencil, residual ~ (h/2)|y"|
    edgeResid = np.stack([take(resid, 0), take(resid, -1)])
    edgeBound = np.stack([np.abs(take(deriv, 1) - take(deriv, 0)),
                          np.abs(take(deriv, -1) - take(deriv, -2))])/2.0
    _flag(edgeResid, edgeBound, 'onesided difference',
          f'on the {label[-1]}-boundary lines', axes=(1,))

    # interior: second-order central stencil, residual ~ (h^2/6)|y'''|
    #
    # The bound is estimated from the SURFACE column, not from the derivative
    # column. Deriving it from `deriv` lets a corrupted derivative inflate its
    # own tolerance: a spike in `deriv` raises the second difference of `deriv`
    # at exactly the node being tested, so the error hides behind the bound it
    # just created (verified -- a 3000x local error passed that way). y is
    # independent of the column under test, so that feedback is broken.
    #   central-difference error ~ (h^2/6)|y'''|, and
    #   y''' ~ (y[i+2] - 2y[i+1] + 2y[i-1] - y[i-2]) / (2h^3)
    # so the bound is |that 4-point difference| / (12h).
    if n >= 5:
        intResid = take(resid, slice(2, -2))
        intBound = np.abs(take(y, slice(4, None))
                          - 2.0*take(y, slice(3, -1))
                          + 2.0*take(y, slice(1, -3))
                          - take(y, slice(None, -4)))/(12.0*h)
        _flag(intResid, intBound, 'central difference', 'in the interior',
              axes=(0, 1))
    else:
        # No weaker-stencil fallback. The derivative-derived bound lets a
        # corrupted derivative column set its own tolerance, so running it
        # would report a pass that means nothing. A grid too small for the
        # 4-point y''' stencil cannot be checked, and an uncheckable file is
        # refused rather than waved through.
        problems.append(
            f'{label} column cannot be verified: the grid has only {n} nodes '
            f'along this axis and the interior derivative check needs at '
            f'least 5 for the y\'\'\' stencil. A fault grid this small is '
            f'almost certainly a mistake in fxmin/fxmax/fzmin/fzmax or dx.')


def validateFaultRoughGeometry(fname, dx, dz, fxmin, fxmax, fzmin, fzmax, dy,
                               ymin=None, ymax=None,
                               nfx=None, nfz=None,
                               derivSafety=FAULT_GEOM_DERIV_SAFETY,
                               slopeWarn=FAULT_GEOM_SLOPE_WARN,
                               offsetFail=FAULT_GEOM_OFFSET_FAIL,
                               verbose=True):
    """Validate a bFault_Rough_Geometry.txt against the fault grid of a case.

    Every argument is a plain number so the function can be used on a file with
    no case directory around it (the converter in scripts/convertFaultGeometry
    validates its own output this way).  validateFaultRoughGeometryForCase is
    the thin wrapper that pulls them off `par`.

    Raises FaultGeometryError listing EVERY problem found, so one case.setup
    tells the user the whole story rather than one item per re-run.
    Returns a dict of diagnostics (including a 'warnings' list) on success.
    """
    header = readFaultRoughGeometryHeader(fname)
    problems = []
    warnings = []

    # --- header node counts vs the case's fault grid ------------------------
    nnxExpect = _faultGeomNodeCount(fxmin, fxmax, dx, 'x', problems)
    nnzExpect = _faultGeomNodeCount(fzmin, fzmax, dz, 'z', problems)
    if nnxExpect is not None and header['nnx'] != nnxExpect:
        problems.append(
            f'header nnx = {header["nnx"]}, but the case fault grid has '
            f'{nnxExpect} nodes along strike '
            f'(round((fxmax {fxmax} - fxmin {fxmin})/dx {dx}) + 1). '
            f'The mesh indexes rough_geo with the case dx, so a mismatch '
            f'reads the wrong node or past the end of the array.')
    if nnzExpect is not None and header['nnz'] != nnzExpect:
        problems.append(
            f'header nnz = {header["nnz"]}, but the case fault grid has '
            f'{nnzExpect} nodes along dip '
            f'(round((fzmax {fzmax} - fzmin {fzmin})/dz {dz}) + 1).')
    # par.nfx/par.nfz are what the rest of the case (on_fault_vars_input.nc,
    # the stress loops) is built on; if the case file drifted from its own
    # fxmin/fxmax/dx that is a case bug worth naming here too.
    if nfx is not None and nnxExpect is not None and int(nfx) != nnxExpect:
        problems.append(
            f'par.nfx = {int(nfx)} disagrees with the fault grid implied by '
            f'par.fxmin/fxmax/dx ({nnxExpect}); the case file is inconsistent '
            f'with itself.')
    if nfz is not None and nnzExpect is not None and int(nfz) != nnzExpect:
        problems.append(
            f'par.nfz = {int(nfz)} disagrees with the fault grid implied by '
            f'par.fzmin/fzmax/dz ({nnzExpect}); the case file is inconsistent '
            f'with itself.')

    # --- spacing and origin -------------------------------------------------
    if not _faultGeomClose(header['dx'], dx):
        problems.append(
            f'header dx = {header["dx"]} m, but the case mesh has dx = {dx} m. '
            f'read_fault_rough_geometry derives the rough grid extent from the '
            f'header dx while insertFaultInterface indexes it with the mesh dx: '
            f'a surface sampled at a different spacing is silently stretched.')
    if not _faultGeomClose(header['fxmin'], fxmin):
        problems.append(
            f'header fxmin = {header["fxmin"]} m, but the case fault starts at '
            f'fxmin = {fxmin} m. The rough grid origin is the corner every '
            f'index is counted from.')
    if not _faultGeomClose(header['fzmin'], fzmin):
        problems.append(
            f'header fzmin = {header["fzmin"]} m, but the case fault bottom is '
            f'fzmin = {fzmin} m.')

    # --- provenance sidecar, when the file carries one ----------------------
    # A surface built elsewhere (another machine, another resolution, the
    # official converter) says what it is here. Disagreement with the file's
    # own header means the two were not written together -- the classic
    # "copied the wrong file in" symptom.
    prov = readFaultGeometryProvenance(fname)
    if prov:
        for key, fileVal in (('dx', header['dx']),
                             ('nnx', header['nnx']),
                             ('nnz', header['nnz']),
                             ('fxmin', header['fxmin']),
                             ('fzmin', header['fzmin'])):
            if key not in prov:
                continue
            try:
                provVal = float(prov[key])
            except ValueError:
                continue
            if abs(provVal - float(fileVal)) > FAULT_GEOM_COORD_TOL:
                problems.append(
                    f'the provenance sidecar '
                    f'{faultGeometryProvenancePath(fname)} says {key} = '
                    f'{prov[key]}, but the file header says {fileVal}. The '
                    f'geometry file and its provenance record were not written '
                    f'together -- one of them came from somewhere else. '
                    f'Regenerate the geometry (case.setup) or delete the stale '
                    f'sidecar.')

    # --- body: row count, shape, then the value-level checks ----------------
    # The header is reported first (above) so a wrong nnx reads as "wrong nnx",
    # not as the row count it implies.
    try:
        _, y, dydx, dydz = readFaultRoughGeometry(fname, header=header)
    except FaultGeometryError as exc:
        problems.append(str(exc))
        raise FaultGeometryError(
            f'{fname} does not match this case\'s fault grid '
            f'({len(problems)} problem(s)):\n  - ' + '\n  - '.join(problems))

    # --- finite values ------------------------------------------------------
    for label, arr in (('surface y', y), ('dy/dx', dydx), ('dy/dz', dydz)):
        bad = ~np.isfinite(arr)
        if bad.any():
            iz, ix = np.argwhere(bad)[0]
            problems.append(
                f'{label} column holds {int(bad.sum())} non-finite value(s) '
                f'(first at ix={ix}, iz={iz}, data row '
                f'{int(ix)*header["nnz"] + int(iz) + 3} of the file). '
                f'NaN/Inf propagates straight into the mesh coordinates.')

    finite = (np.isfinite(y).all() and np.isfinite(dydx).all()
              and np.isfinite(dydz).all())

    # --- surface stays inside the model domain ------------------------------
    if finite and ymin is not None and ymax is not None:
        if y.max() >= ymax or y.min() <= ymin:
            problems.append(
                f'surface y range [{y.min():.4g}, {y.max():.4g}] m is not '
                f'strictly inside the model domain y range [{ymin}, {ymax}] m. '
                f'insertFaultInterface maps off-fault nodes with '
                f'(ymax - peak)/ymax, which changes sign once |peak| reaches '
                f'the domain boundary and inverts the mesh. A surface this '
                f'large is usually a unit error (km supplied as m).')

    # --- derivative columns vs finite differences of the surface ------------
    if finite:
        _checkDerivativeColumn(y, dydx, dx, 1, 'dy/dx', problems, derivSafety)
        _checkDerivativeColumn(y, dydz, dz, 0, 'dy/dz', problems, derivSafety)

    # --- slope magnitude and per-cell offset --------------------------------
    maxSlope = 0.0
    maxOffset = 0.0
    if finite:
        maxSlope = float(max(np.abs(dydx).max(), np.abs(dydz).max()))

        # No silent stand-in for dy (rule 1). The threshold counts FAULT-NORMAL
        # cells, so substituting dx for an unknown dy silently changes the
        # units of the quantity being compared to 1.0, and the check stops
        # meaning what its message says.
        if dy is None:
            # Not a warning and not a silent dx substitution: the per-cell
            # offset is the ONLY hard limit on surface steepness, and a file
            # that has not been checked for element tangling must not be
            # accepted. Refuse instead of reporting a pass we did not earn.
            problems.append(
                'the fault-normal spacing dy was not supplied, so the '
                'per-cell fault-normal offset -- the only hard limit on '
                'surface steepness, and the check that catches element '
                'tangling -- cannot be evaluated. Pass dy (case.setup does, '
                'via par.dy).')
            dyCell = None
        else:
            dyCell = dy
            maxOffset = float(max(np.abs(dydx).max()*dx/dyCell,
                                  np.abs(dydz).max()*dz/dyCell))

        if dyCell is not None and maxOffset >= offsetFail:
            problems.append(
                f'the fault surface climbs by {maxOffset:.4g} of a '
                f'fault-normal cell (dy = {dyCell} m) between adjacent fault '
                f'nodes, which is >= {offsetFail}: max|dy/dx|*dx/dy = '
                f'{np.abs(dydx).max()*dx/dyCell:.4g}, max|dy/dz|*dz/dy = '
                f'{np.abs(dydz).max()*dz/dyCell:.4g}. At an offset of one cell '
                f'a column\'s fault node passes its neighbour\'s first '
                f'off-fault node layer and the elements inserted by '
                f'insertFaultInterface tangle.')

        if 0.0 < maxSlope < FAULT_GEOM_SLOPE_FLAT_WARN:
            warnings.append(
                f'maximum |dy/dx| or |dy/dz| is only {maxSlope:.4g}, below '
                f'{FAULT_GEOM_SLOPE_FLAT_WARN:g}: this surface is nearly '
                f'planar. If it is meant to be rough, check its UNITS -- a '
                f'surface supplied in km where EQdyna wants m is self-'
                f'consistent (y and both derivative columns scale together), '
                f'so no other check here can see it. The official TPV29 '
                f'surface has max |slope| = 0.65. If the fault really is '
                f'planar, prefer insertFaultType 1 or 2 over supplying a '
                f'flat type-3 file.')

        if maxSlope > slopeWarn:
            warnings.append(
                f'maximum |dy/dx| or |dy/dz| is {maxSlope:.4g} > {slopeWarn}: '
                f'the surface is rough enough to distort elements (the same '
                f'threshold scripts/generateFaultInterface warns at). Every '
                f'gated case that inserts a fault interface sits here '
                f'(test.drv.a6 0.38, test.tpv29 0.65, test.tpv10 0.577 = '
                f'cot 60 for its dipping plane), so this is a note, not a '
                f'defect. The hard limit is the per-cell offset, checked '
                f'above: the surface climbs {maxOffset:.4g} of a '
                f'fault-normal cell between adjacent fault nodes.')

    if problems:
        raise FaultGeometryError(
            f'{fname} does not match this case\'s fault grid '
            f'({len(problems)} problem(s)):\n  - ' + '\n  - '.join(problems))

    diagnostics = dict(fname=fname, nnx=header['nnx'], nnz=header['nnz'],
                       dx=header['dx'], fxmin=header['fxmin'],
                       fzmin=header['fzmin'],
                       yMin=float(y.min()), yMax=float(y.max()),
                       maxSlope=maxSlope, maxOffset=maxOffset,
                       provenance=prov, warnings=warnings)
    if verbose:
        print(f'validateFaultRoughGeometry: {fname} OK -- '
              f'nnx={header["nnx"]}, nnz={header["nnz"]}, dx={header["dx"]} m, '
              f'corner=({header["fxmin"]}, {header["fzmin"]}) m, '
              f'y in [{y.min():.1f}, {y.max():.1f}] m, '
              f'max|slope|={maxSlope:.4g}, '
              f'max per-cell offset={maxOffset:.4g} cell')
        if prov and prov.get('source'):
            print(f'validateFaultRoughGeometry: source -- {prov["source"]}'
                  + (f' (native spacing {prov["sourceDx"]} m)'
                     if prov.get('sourceDx') else ''))
        for w in warnings:
            print(f'validateFaultRoughGeometry: WARNING - {w}')
    return diagnostics


def validateFaultRoughGeometryForCase(par, fname=FAULT_ROUGH_GEOMETRY_FILE,
                                      verbose=True):
    """validateFaultRoughGeometry with the case's own parameters.

    Called by scripts/case.setup for every insertFaultType > 0 -- including the
    files case.setup just generated (types 1 and 2), so the generator is held
    to the same contract as a supplied file.

    par.faultGeometrySourceDx (defaultParameters.py) is the compset's
    declaration of the native sampling of the surface it ships; when it is set,
    par.dx is checked against it here too.  Compsets should also call
    requireFaultGeometryResolution from their own tooling so the failure comes
    before that tooling runs, not after.
    """
    requireFaultGeometryResolution(
        par.dx, getattr(par, 'faultGeometrySourceDx', None),
        getattr(par, 'faultGeometrySourceName',
                'the fault-geometry source this compset ships'))
    # Direct attribute access, not getattr(..., None) (rule 2). Every one of
    # these is defined in defaultParameters.py, so a missing one means the par
    # object is not a real case -- that must raise AttributeError here, not
    # silently disable the check that depends on it downstream. The `dy=None`
    # path in particular used to disable the per-cell offset check, the only
    # hard limit on surface steepness.
    return validateFaultRoughGeometry(
        fname,
        dx=par.dx, dz=par.dz, dy=par.dy,
        fxmin=par.fxmin, fxmax=par.fxmax,
        fzmin=par.fzmin, fzmax=par.fzmax,
        ymin=par.ymin, ymax=par.ymax,
        nfx=par.nfx, nfz=par.nfz,
        verbose=verbose)


def ensureFaultRoughGeometryForCase(par, fname=FAULT_ROUGH_GEOMETRY_FILE,
                                    verbose=True):
    """Make sure the case has a bFault_Rough_Geometry.txt that is right FOR
    THIS CASE, and return the validator diagnostics.

    Order of preference, and the reason for it:

      1. If the file is already there and already validates against this case,
         LEAVE IT EXACTLY AS IT IS.  This is what makes case.setup idempotent
         (run it twice, get the same bytes and the same verdict) and it is what
         stops a case author's deliberately placed surface from being silently
         clobbered -- the 50 m TPV29 run hit exactly that: importing
         user_defined_params used to overwrite the geometry file as a side
         effect, so a correct hand-built surface was replaced by whatever the
         compset's default path produced.  That is why compsets now hand
         case.setup a WRITER (par.faultGeometryWriter) instead of writing at
         import time.
      2. Otherwise, if the compset supplied par.faultGeometryWriter, call it,
         say out loud that the file was (re)generated and why, and validate the
         result -- which must pass.
      3. Otherwise validate whatever is there, so a case that ships a static
         file still gets the full check and an actionable failure.

    par.faultGeometryWriter is a callable taking the output path.
    """
    requireFaultGeometryResolution(
        par.dx, getattr(par, 'faultGeometrySourceDx', None),
        getattr(par, 'faultGeometrySourceName',
                'the fault-geometry source this compset ships'),
        availableDx=getattr(par, 'faultGeometrySourceAvailableDx', None))

    writer = getattr(par, 'faultGeometryWriter', None)
    if exists(fname):
        try:
            d = validateFaultRoughGeometryForCase(par, fname, verbose=verbose)
            if verbose:
                print(f'ensureFaultRoughGeometry: {fname} already matches this '
                      f'case; left untouched.')
            return d
        except FaultGeometryError as why:
            if writer is None:
                raise
            if verbose:
                print(f'ensureFaultRoughGeometry: regenerating {fname} -- the '
                      f'existing file does not match this case:\n{why}')
    elif writer is not None and verbose:
        print(f'ensureFaultRoughGeometry: generating {fname} for this case.')

    if writer is not None:
        writer(fname)
    return validateFaultRoughGeometryForCase(par, fname, verbose=verbose)


def writeFaultRoughGeometry(y, dydx, dydz, dx, fxmin, fzmin,
                            fname=FAULT_ROUGH_GEOMETRY_FILE, provenance=None):
    """Write arrays shaped (nnz, nnx), [iz, ix], in EQdyna's z-fastest order.

    `provenance`, when given, is a dict of free-form key/value strings written
    to the sidecar named by faultGeometryProvenancePath(fname).
    """
    y, dydx, dydz = np.asarray(y), np.asarray(dydx), np.asarray(dydz)
    if y.shape != dydx.shape or y.shape != dydz.shape:
        raise FaultGeometryError(
            f'writeFaultRoughGeometry: y {y.shape}, dy/dx {dydx.shape} and '
            f'dy/dz {dydz.shape} must all have the same (nnz, nnx) shape.')
    nnz, nnx = y.shape
    with open(fname, 'w') as f:
        f.write(f'{nnx}\t{nnz}\t0\n')
        f.write(f'{dx:.6f}\t{fxmin:.6f}\t{fzmin:.6f}\n')
        for ix in range(nnx):
            for iz in range(nnz):
                f.write(f'{y[iz, ix]:.7e}\t{dydx[iz, ix]:.7e}\t'
                        f'{dydz[iz, ix]:.7e}\n')
    if provenance is not None:
        fields = dict(provenance)
        # .12g, not g: the sidecar's numbers are cross-checked against the
        # file header to FAULT_GEOM_COORD_TOL (1 mm), and a dipping fault's
        # fzmin (= fzmin*sin(dip)) is an irrational-looking full-precision
        # float -- 6-significant-digit '%g' loses 2 cm of it and the check
        # then fires on its own rounding.
        fields.setdefault('dx', f'{dx:.12g}')
        fields.setdefault('nnx', str(nnx))
        fields.setdefault('nnz', str(nnz))
        fields.setdefault('fxmin', f'{fxmin:.12g}')
        fields.setdefault('fzmin', f'{fzmin:.12g}')
        writeFaultGeometryProvenance(fname, fields)
    return fname


# --- provenance sidecar ----------------------------------------------------
# bFault_Rough_Geometry.txt's first two rows are consumed POSITIONALLY by
# list-directed Fortran reads, so provenance cannot go in the file itself
# without risking that parse -- and adding bytes to a gated case's geometry
# file is exactly the sort of change that has to be proven inert.  It goes in a
# sidecar instead: same name plus '.provenance', "key: value" per line.  It is
# advisory -- EQdyna never reads it -- but when it IS present the validator
# cross-checks it against the file's own header, which is how a surface copied
# in from another machine or another resolution announces itself instead of
# being silently meshed.

FAULT_GEOM_PROVENANCE_SUFFIX = '.provenance'


def faultGeometryProvenancePath(fname=FAULT_ROUGH_GEOMETRY_FILE):
    return fname + FAULT_GEOM_PROVENANCE_SUFFIX


def writeFaultGeometryProvenance(fname=FAULT_ROUGH_GEOMETRY_FILE, fields=None):
    """Write the sidecar. Deterministic: no timestamps, so re-running a
    generator produces byte-identical output (idempotency)."""
    path = faultGeometryProvenancePath(fname)
    with open(path, 'w') as f:
        f.write(f'# provenance for {os.path.basename(fname)} -- advisory; '
                f'EQdyna does not read this file.\n')
        for k, v in (fields or {}).items():
            f.write(f'{k}: {v}\n')
    return path


def readFaultGeometryProvenance(fname=FAULT_ROUGH_GEOMETRY_FILE):
    """Return the sidecar's fields as a dict, or None when there is none."""
    path = faultGeometryProvenancePath(fname)
    if not exists(path):
        return None
    fields = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or ':' not in line:
                continue
            k, v = line.split(':', 1)
            fields[k.strip()] = v.strip()
    return fields
