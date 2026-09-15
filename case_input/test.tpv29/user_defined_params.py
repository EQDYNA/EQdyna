#! /usr/bin/env python3

# SCEC TPV29: right-lateral vertical strike-slip rough fault in an elastic
# half-space (strike.scec.org/cvws, TPV29_30_Description_v06, Jan 17 2015).
#
# Fault: 40 km x 20 km, reaches the surface, with the OFFICIAL randomly
# generated roughness (Hurst exponent 1) shipped in
# bFault_Rough_Geometry.tpv29.50m.txt.gz (exact 2:1 decimation of the official
# 25 m data file, at the 50 m SPEC resolution; see tpv29GeometryTools.py for
# the frame mapping and provenance).  Hypocenter 15 km from the left fault
# edge (x = -5 km), 10 km deep.
#
# Spec parameters encoded here:
#   material          rho 2670, Vs 3464, Vp 6000 (defaults in
#                     defaultParameters.py already match)
#   friction (SW)     mu_s 0.18, mu_d 0.12, d0 0.30 m
#   frictional        C0 = 0.4 MPa + 0.0002 MPa/m * (4000 m - depth) above
#   cohesion          4 km depth; 0.4 MPa below
#   nucleation        smoothed forced rupture, r_crit 4000 m, t0 0.5 s;
#                     par.tpv = 36 selects the forced-rupture-time formula in
#                     faulting.f90:swtwNucleation, which is exactly the
#                     TPV29 spec formula T(r) = (r + 0.081*r_crit*
#                     (1/(1-(r/r_crit)^2)-1)) / (0.7*Vs) with Vs fixed at
#                     3464 m/s (NUC_* constants in globalvar.f90)
#   initial stress    depth-dependent effective stress tensor (spec p.12-13)
#                     resolved onto the LOCAL rough-fault normal at every
#                     fault grid point (see loop below); gravity 9.8 exactly
#   run time          0 - 20 s
#   resolution        spec standard is 50 m (100 m acceptable). par.dx below
#                     is the fast-gate value; the full tier runs the 50 m spec.
#                     Any par.dx that is a multiple of the shipped 50 m source
#                     and divides 40/20 km works; anything else is refused by
#                     lib.requireFaultGeometryResolution with the reason.
#
# insertFaultType = 3: any value > 0 activates EQdyna's rough-fault machinery
# (read_fault_rough_geometry + insertFaultInterface + per-node normals).
# Values 1/2 would make scripts/generateFaultInterface (invoked at the end of
# case.setup) OVERWRITE bFault_Rough_Geometry.txt with a planar/synthetic
# fractal surface; with 3 that script is not called at all, and the official
# TPV29 surface is what EQdyna reads.
#
# This module does NOT write bFault_Rough_Geometry.txt. It computes the surface
# it needs in memory (for the hypocenter and the per-node stress resolution) and
# hands case.setup a WRITER (par.faultGeometryWriter) instead. Writing at import
# made every import of this file -- case.setup, plotRuptureDynamics, any tool --
# clobber whatever geometry file was on disk, which silently replaced a
# deliberately placed 50 m surface with the default one, and made case.setup
# non-idempotent. case.setup now keeps an already-correct file untouched and
# calls the writer only when the file is missing or does not match the case.

from defaultParameters import *
from math import *
from lib import *
import numpy as np
import tpv29GeometryTools as geoTools

par = parameters()

par.xmin, par.xmax = -35.0e3, 35.0e3
par.ymin, par.ymax = -30.0e3, 30.0e3
par.zmin, par.zmax = -35.0e3, 0.0e3

par.fxmin, par.fxmax = -20.0e3, 20.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -20.0e3, 0.0e3

par.dx = 500.   # fast gate; full tier overrides to the spec value
par.dy = par.dx
par.dz = par.dx

par.nmat = 1
par.vp, par.vs, par.rou = 6000., 3464., 2670.

par.term = 20.

par.C_elastic = 1
par.C_nuclea  = 1
par.insertFaultType = 3   # >0: rough fault; 3 (not 1/2): keep official geometry
par.friclaw = 1
par.tpv = 36              # selects the TPV29-spec forced-rupture-time formula
par.nucR = 4.0e3          # r_crit, m

par.dt = 0.5*par.dx/par.vp

# grav is fixed at exactly 9.8 m/s^2 in the spec
par.grav = 9.8
par.fric_sw_fs = 0.18
par.fric_sw_fd = 0.12
par.fric_sw_D0 = 0.30

par.nx, par.ny, par.nz = 2, 1, 2   # ny=1: keeps the fault plane off every MPI partition
par.HPC_ncpu  = par.nx*par.ny*par.nz
par.HPC_nnode = round(floor(par.HPC_ncpu/128)) + 1
par.HPC_queue = "normal"
par.HPC_time  = "10:00:00"
par.HPC_account = "EAR22013"
par.HPC_email = ""

# Fault grid
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)

# The compset ships the official surface at 50 m (the spec resolution), so
# par.dx must be an integer multiple of 50 m; declaring the source spacing here
# lets scripts/case.setup make the same check for any case that supplies its
# own geometry (lib.requireFaultGeometryResolution).
par.faultGeometrySourceDx   = geoTools.SHIPPED_DX
par.faultGeometrySourceName = geoTools.SHIPPED_NAME
par.faultGeometrySourceAvailableDx = [d for d, _ in geoTools.SHIPPED_SOURCES]

# Official rough geometry on the case grid ([iz, ix], iz=0 at fzmin). Needed
# here in memory for the hypocenter y and the per-node stress resolution below.
fY, fDydx, fDydz = geoTools.faultGridForCase(par.dx)

# How case.setup (re)writes bFault_Rough_Geometry.txt when it has to. Not
# called here: see the note at the top of this file.
par.faultGeometryWriter = lambda out: geoTools.writeBFault(
    fY, fDydx, fDydz, par.dx, out=out)

# Hypocenter: x = -5 km along strike, 10 km deep, ON the rough surface.
par.xsource, par.zsource = -5.0e3, -10.0e3
par.ysource = fY[round((par.zsource - par.fzmin)/par.dz),
                 round((par.xsource - par.fxmin)/par.dx)]

# --- TPV29 initial stresses: depth-dependent EFFECTIVE stress tensor -------
# (spec p.12-13; effective = total + Pf*I, Pf = 1000*9.8*depth hydrostatic)
#   sigma22_eff = -(2670-1000)*9.8*depth              (vertical, z_eq)
#   sigma11_eff = (Omega*b11 + (1-Omega))*sigma22_eff (fault-parallel, x_eq)
#   sigma33_eff = (Omega*b33 + (1-Omega))*sigma22_eff (fault-normal, y_eq)
#   sigma13     = Omega*b13*sigma22_eff               (x_eq-y_eq shear)
#   Omega = 1 (depth<=17 km), (22000-depth)/5000 (17-22 km), 0 (>22 km)
# Resolved per node onto the local rough-fault (n, s, d) frame constructed
# exactly as meshgen.f90:createMasterNode does from (dy/dx, dy/dz).
B11, B33, B13 = 1.025837, 0.974162, -0.158649

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    depth = abs(zcoor)
    # friction
    par.on_fault_vars[iz,ix,1] = par.fric_sw_fs
    # slip goes to zero at the fault border (left/right/bottom edges;
    # the free surface is NOT a border): make border nodes unbreakable
    if ix == 0 or ix == par.nfx-1 or iz == 0:
        par.on_fault_vars[iz,ix,1] = 10000.
    par.on_fault_vars[iz,ix,2] = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3] = par.fric_sw_D0
    # frictional cohesion C0, tapered in the uppermost 4 km
    if depth < 4000.:
        par.on_fault_vars[iz,ix,4] = 0.4e6 + 200.0*(4000. - depth)
    else:
        par.on_fault_vars[iz,ix,4] = 0.4e6
    # forced-rupture decay time t0
    par.on_fault_vars[iz,ix,5] = 0.5

    # initial effective tractions from the depth-dependent stress tensor;
    # at the free-surface node use a half-element depth (spec part 8)
    dEff = max(depth, par.dx/2.)
    omega = 1.0 if dEff <= 17000. else max(0.0, (22000. - dEff)/5000.)
    s22 = -(par.rou - 1000.)*par.grav*dEff  # effective vertical stress
    sxx = (omega*B11 + (1.-omega))*s22      # x_eq (fault-parallel)
    syy = (omega*B33 + (1.-omega))*s22      # y_eq (fault-normal)
    szz = s22                               # z_eq (vertical)
    sxy = omega*B13*s22                     # right-lateral driving shear

    pfx, pfz = fDydx[iz,ix], fDydz[iz,ix]
    nn = np.array([-pfx, 1.0, -pfz])/sqrt(pfx**2 + 1.0 + pfz**2)
    ss = np.array([1.0, pfx, 0.0])/sqrt(1.0 + pfx**2)
    dd = np.cross(ss, nn)
    S  = np.array([[sxx, sxy, 0.0],
                   [sxy, syy, 0.0],
                   [0.0, 0.0, szz]])
    t  = S @ nn
    par.on_fault_vars[iz,ix,7]  = nn @ t   # initial effective normal
    par.on_fault_vars[iz,ix,8]  = ss @ t   # initial strike shear
    par.on_fault_vars[iz,ix,49] = dd @ t   # initial dip shear

    par.on_fault_vars[iz,ix,46] = 0.0      # no creep/background slip rate

# --- Stations (spec parts 9-10) --------------------------------------------
# on-fault: (x, -depth) in km
par.st_coor_on_fault = [
    [-18.0,-15.6], [-15.0,-5.0], [-15.0,-12.0], [-11.0,-1.4], [-8.9,-10.1],
    [-5.0,0.0], [-5.0,-10.0], [-5.0,-16.0], [-4.2,-6.1], [0.0,-12.0],
    [4.3,-6.2], [4.6,-6.0], [5.0,0.0], [5.0,-12.0], [5.1,-5.7],
    [5.9,-4.7], [9.0,-15.3], [9.3,-15.3], [10.0,-5.0], [10.0,-11.0],
    [15.0,0.0], [15.0,-13.0], [16.7,-10.5], [17.0,-4.5]]
# off-fault: (x, y_eq, z_eq) in km; y_eq = TPV fault-normal coordinate
# (positive = far side), all stations at the surface
par.st_coor_off_fault = [
    [-20.0,-20.0,0.0], [0.0,-20.0,0.0], [20.0,-20.0,0.0],
    [-15.0,-3.0,0.0],  [0.0,-3.0,0.0],  [15.0,-3.0,0.0],
    [-15.0,3.0,0.0],   [0.0,3.0,0.0],   [15.0,3.0,0.0],
    [-20.0,20.0,0.0],  [0.0,20.0,0.0],  [20.0,20.0,0.0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
