#! /usr/bin/env python3

# SCEC TPV30: right-lateral vertical strike-slip rough fault in a
# NON-ASSOCIATIVE DRUCKER-PRAGER VISCOPLASTIC half-space
# (strike.scec.org/cvws, TPV29_30_Description_v06, Jan 17 2015).
#
# TPV30 is TPV29 with one change only.  Quoting the spec (p.11):
#   "The material properties are the only difference between the elastic
#    benchmark (TPV29), and the viscoplastic benchmark (TPV30)."
# Everything below that is not inside the "--- TPV30 off-fault plasticity"
# block is identical to case_input/test.tpv29 -- same official rough surface
# (md5 f6df5ececc9d95d0c6db43c9ba0c3c6e), same initial stress tensor, same
# slip-weakening friction, same nucleation, same stations, same 20 s run.
#
# Fault: 40 km x 20 km, reaches the surface, with the OFFICIAL randomly
# generated roughness (Hurst exponent 1) shipped in
# bFault_Rough_Geometry.tpv29.100m.txt (exact 4:1 decimation of the official
# 25 m data file; see tpv29GeometryTools.py for the frame mapping and
# provenance -- byte-identical copy of case_input/test.tpv29's tool and file,
# since TPV30 uses the identical surface per the spec quote above).
# Hypocenter 15 km from the left fault edge (x = -5 km), 10 km deep.
#
# Spec parameters encoded here:
#   material          rho 2670, Vs 3464, Vp 6000 (defaults in
#                     defaultParameters.py already match)
#   friction (SW)     mu_s 0.18, mu_d 0.12, d0 0.30 m
#   frictional        C0 = 0.4 MPa + 0.0002 MPa/m * (4000 m - depth) above
#   cohesion          4 km depth; 0.4 MPa below
#   nucleation        smoothed forced rupture, r_crit 4000 m, t0 0.5 s;
#                     par.tpv = 30 selects the shared forced-rupture-time
#                     formula in faulting.f90:swtwNucleation (Part 6, p.18,
#                     confirmed identical for TPV29 and TPV30 -- spec p.16:
#                     "Benchmarks TPV29 and TPV30 use linear slip-weakening
#                     friction" -- so TPV30 gets its OWN branch entry rather
#                     than impersonating TPV29/36 (rule 17 step 3; that
#                     impersonation pattern is exactly what hid TPV29's own
#                     nucleation gap for months)
#   initial stress    depth-dependent effective stress tensor (spec p.12-13)
#                     resolved onto the LOCAL rough-fault normal at every
#                     fault grid point (see loop below); gravity 9.8 exactly
#   off-fault         Drucker-Prager viscoplastic (spec Part 7, p.19-20):
#   plasticity        c = 1.18 MPa, nu = 0.1680, Tv = 0.05 s (real input
#                     slots as of item 24(b)/(c), 2026-09-17 -- see below),
#                     Pf hydrostatic (see the block below)
#   run time          0 - 20 s
#   resolution        spec standard is 50 m (100 m acceptable, spec p.15);
#                     this case is built at par.dx = 500 m as the FAST GATE
#                     (rule 17 step 4), matching test.tpv29's own gate.  The
#                     100 m file shipped here decimates to 500/200/100 m; the
#                     50 m spec file is NOT shipped yet (see README.md) --
#                     recorded, not run, per rule 17 step 5.
#
# --- Gaps in the OLD scratch/tpv30 draft, and their current status ----------
# (case_input_draft/test.tpv30/user_defined_params.py, written before item 24
# landed, documented six gaps G1-G6; status here, not carried forward stale):
#
#  G1 (Tv had no input slot) -- RESOLVED 2026-09-17 (item 24(b)):
#     par.viscoplasticRelaxTime below is a real input now.
#  G2 (off-fault deviatoric depth taper missing) -- RESOLVED 2026-09-17
#     (item 24(c)): par.devStrTaperDepthStart/End below are real inputs.
#  G3 (+7.3215 m element-centre depth offset in meshgen.f90) -- UNTOUCHED, per
#     owner decision 2026-09-16 (item 24(d)); LOAD-BEARING (c7c4f5f measured
#     that removing it moves the traction-check answer away from the correct
#     ratio), so this case does not attempt to work around it.
#  G4 (plastic-strain output window hardcoded to 5x2x8 km, sized for
#     test.drv.a6) -- RESOLVED 2026-09-17 (item 24(f)): par.plasticOutputHalfWidth
#     below is set for THIS case's 40x20 km fault instead of left at the
#     drv.a6-sized default.
#  G5 (b11 + b33 = 1.999999, not exactly 2, per the spec's own coefficients)
#     -- negligible (1e-6*|sv|, 164 Pa at 10 km depth), not a code gap, not
#     addressed (matches TPV29's own compset, which carries the same
#     coefficients unchanged).
#  G6 (on-fault initial traction allegedly HALF the spec value under
#     C_elastic=0) -- DOES NOT REPRODUCE. This was the BLOCKER the old draft
#     stopped on. Measured on test.drv.a6 (c7c4f5f, 2026-09-16): Tn ratio
#     1.0018, Ts ratio 0.9858 against setPlasticStress's own formula -- not
#     half. Root cause was item 26's `arn` MPI double-count, already fixed
#     2026-09-14; the old draft's G6 note predates that fix and is stale.
#
# insertFaultType = 3: any value > 0 activates EQdyna's rough-fault machinery
# (read_fault_rough_geometry + insertFaultInterface + per-node normals).
# Values 1/2 would make scripts/generateFaultInterface (invoked at the end of
# case.setup) OVERWRITE bFault_Rough_Geometry.txt with a planar/synthetic
# fractal surface; with 3 that script is not called at all, and the official
# TPV29/30 surface written below is what EQdyna reads.
#
# This module does NOT write bFault_Rough_Geometry.txt at import time. It
# computes the surface it needs in memory (hypocenter + per-node stress
# resolution) and hands case.setup a WRITER (par.faultGeometryWriter)
# instead -- same lazy pattern as test.tpv29 (writing at import made every
# import clobber a deliberately placed geometry file and made case.setup
# non-idempotent; see test.tpv29/user_defined_params.py's own note). The
# OLD scratch draft wrote at import time; that habit is NOT carried forward.

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

# --- TPV30 off-fault plasticity -------------------------------------------
# Spec p.11 / p.19 constitutive parameters:
#     cohesion c  = 1.18 MPa        -> par.coheplas
#     bulk friction nu = 0.1680     -> par.bulk   (phi = atan(nu) = 9.5366 deg)
#     relaxation time Tv = 0.05 s   -> par.viscoplasticRelaxTime (item 24(b);
#          resolveViscoplasticParams in scripts/lib.py writes it to bGlobal.txt
#          and readInputFiles.f90 reads it directly -- no more mesh-derived
#          2*dz/3464 stand-in for this case)
#     Pf = 1000*9.8*depth (hydrostatic) -> par.gamar = 0 (gamar is the
#          overpressurization coefficient; EQdyna carries the EFFECTIVE stress
#          tensor with eleporep = 0, and meshgen.f90:setPlasticStress builds it
#          as -(roumax - rhow*(gamar+1))*g*depth, so gamar = 0 gives exactly
#          the spec's sigma22 + Pf = -(2670-1000)*9.8*depth.)
# calcElemKU.f90's Drucker-Prager block is the spec return map verbatim:
#     yield  = ccosphi - sinphi*(strmea + porep)     = c cos(phi) - (sm+Pf) sin(phi)
#     rjust  = yield/taomax + (1-yield/taomax)*exp(-dt/tv)
#            = exp(-dt/Tv) + (1-exp(-dt/Tv)) Y/sqrt(J2)     [spec Part 7, Step 4]
# SCEC initial-condition method (spec TPV29_30_Description_v06, p.24; pathway
# item 97): this case uses METHOD 2, "explicit gravity and boundary
# tractions". C_elastic = 0 turns the gravity body force on and the stored
# stress is the TOTAL effective tensor built by meshgen.f90:setPlasticStress.
# test.tpv29 uses METHOD 1 (stress change, no gravity). The two are mutually
# exclusive per the spec. Pore pressure is implicit: it is folded into
# roumax - rhow*(gamar+1) and eleporep is hardcoded 0.0d0 (item 97(2)).
par.C_elastic = 0
par.coheplas  = 1.18e6    # plastic cohesion c, Pa (NOT the frictional C0 below)
par.bulk      = 0.1680    # bulk friction nu
par.gamar     = 0.0       # hydrostatic fluid pressure (no overpressurization)
par.roumax    = par.rou   # keep bGlobal's roumax tied to the material above
par.output_plastic = 1    # requires C_elastic == 0 (checkInputConsistency.f90:15)

# Item 24(b): SCEC spec p.19, Viscoplastic relaxation time Tv = 0.05 s exactly.
par.viscoplasticRelaxTime = 0.05

# Item 24(c): SCEC spec p.13, Omega(depth) taper on the off-fault deviatoric
# pre-stress, linear to zero between 17 and 22 km depth.
par.devStrTaperDepthStart = 17.0e3
par.devStrTaperDepthEnd   = 22.0e3

# Item 24(f): plastic-strain output window. The 5x2x8 km default (sized for
# test.drv.a6) would clip the vast majority of TPV30's 40x20 km fault -- set
# half-widths to cover the whole fault (|x|<=20km) plus margin, most of the
# y-domain, and down to (and slightly past) the deviatoric taper's 22 km
# floor, without reaching into the PML-adjacent domain edges (|x|<35km,
# |y|<30km, |z|<35km).
par.plasticOutputHalfWidth = (25.0e3, 20.0e3, 22.0e3)

# Off-fault initial stress tensor.  meshgen.f90:setPlasticStress parameterizes
# it as
#     szz = sv,  sxx = sv - dev*cos(2a),  syy = sv + dev*cos(2a),
#     sxy = dev*sin(2a),   sv = -(roumax - rhow*(gamar+1))*g*depth,  dev = |sv|*R
# The TPV29/30 spec (p.13) gives the same tensor as
#     sxx = b11*sv, syy = b33*sv, sxy = b13*sv,  szz = sv
#     b11 = 1.025837, b33 = 0.974162, b13 = -0.158649
# Matching term by term:  R*cos(2a) = b11 - 1,  R*sin(2a) = -b13, hence
par.str1ToFaultAngle     = 40.375111    # a, degrees
par.devStrToStrVertRatio = 0.160739092  # R
# Check (spec p.13): at 10 km depth this gives sqrt(J2) = 26.307 MPa and
# Y = 28.279 MPa, i.e. 93.0% of the yield stress -- exactly the spec's stated
# "initial deviatoric stress at hypocenter depth is 93% of the plastic yield
# stress", which independently confirms c, nu, gamar and the b-coefficients.

par.C_nuclea  = 1
par.insertFaultType = 3   # >0: rough fault; 3 (not 1/2): keep official geometry
par.friclaw = 1
par.tpv = 30              # this case IS TPV30; faulting.f90:swtwNucleation
                          # selects the spec's smoothed forced-rupture formula
                          # for TPV 29/30/36/37/201 (they share it -- Part 6)
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

# The compset ships the official surface at 100 m (see README.md for why the
# 50 m spec file is not shipped yet), so par.dx must be an integer multiple
# of 100 m for now; declaring the source spacing here lets scripts/case.setup
# make the same check (lib.requireFaultGeometryResolution).
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

# --- TPV29/30 initial stresses: depth-dependent EFFECTIVE stress tensor ----
# (spec p.12-13; effective = total + Pf*I, Pf = 1000*9.8*depth hydrostatic)
#   sigma22_eff = -(2670-1000)*9.8*depth              (vertical, z_eq)
#   sigma11_eff = (Omega*b11 + (1-Omega))*sigma22_eff (fault-parallel, x_eq)
#   sigma33_eff = (Omega*b33 + (1-Omega))*sigma22_eff (fault-normal, y_eq)
#   sigma13     = Omega*b13*sigma22_eff               (x_eq-y_eq shear)
#   Omega = 1 (depth<=17 km), (22000-depth)/5000 (17-22 km), 0 (>22 km)
# Resolved per node onto the local rough-fault (n, s, d) frame constructed
# exactly as meshgen.f90:createMasterNode does from (dy/dx, dy/dz). This is
# the ON-FAULT traction (fric() array); it is a SEPARATE quantity from the
# off-fault element stress par.devStrTaperDepthStart/End above tapers (that
# one feeds meshgen.f90:setPlasticStress, this one feeds the fric() slots
# below) -- both use the SAME Omega(depth) taper because the spec applies it
# to both, but they are two different code paths and both had to be checked.
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
