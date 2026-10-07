#! /usr/bin/env python3

# SCEC TPV26: right-lateral vertical strike-slip PLANAR fault in a linear
# elastic half-space (strike.scec.org/cvws, TPV26_27_Description_v13,
# Jan 9 2014).
#
# Fault: 40 km x 20 km, planar, vertical, reaches the surface. Hypocenter
# 15 km from the left edge, 10 km deep -- spec coordinate (x,y,z) =
# (-5000, 10000, 0), where the SPEC's y is depth and z is fault-normal; in
# EQdyna's own axes (x=strike, y=fault-normal, z=depth, negative down) that
# is par.xsource=-5.0e3, par.ysource=0.0, par.zsource=-10.0e3 -- the SAME
# hypocenter EQdyna already uses for test.tpv29/test.tpv30 (same fault
# dimensions, same geometry diagram, Part 2).
#
# This compset is a PLANAR analogue of test.tpv30's stress/nucleation
# machinery (no rough-fault geometry needed -- the geometry is a flat plane
# at y=0, same style as test.tpv8). See that case's user_defined_params.py
# for the derivation pattern this one follows.
#
# Spec parameters encoded here (TPV26_27_Description_v13):
#   material (Part 3, p.5)   rho=2670, Vs=3464, Vp=6000 (elastic)
#   friction (Part 3, p.8)   mu_s=0.18, mu_d=0.12, d0=0.30 m
#   frictional cohesion C0   0.40 MPa + 0.00072 MPa/m*(5000m-depth) if
#                            depth<=5000m, else 0.40 MPa (4.00 MPa at
#                            surface, p.8)
#   nucleation (Part 5,      smoothed forced rupture, r_crit=4000 m,
#     p.12)                  t0=0.5 s; par.tpv=26 selects the shared
#                            forced-rupture-time formula in
#                            faulting.f90:swtwNucleation (confirmed by
#                            direct spec read, row 150 investigation,
#                            2026-10-06 -- NOT assumed from TPV29/30's
#                            formula matching by coincidence: TPV26/27 Part
#                            5, p.12 is the identical r_crit=4000,
#                            0.081 taper coefficient, 0.7*Vs equation).
#   initial stress (Part 3,  depth-dependent effective stress tensor
#     p.6-7)                 (b11=0.926793, b33=1.073206, b13=-0.169029),
#                            Omega(depth) taper 15000-20000 m, gravity
#                            9.8 m/s^2 exactly. SCEC "Method 1" (stress
#                            change, spec Part 7): for the ELASTIC case the
#                            initial stress tensor need not be stored
#                            per-element at all -- only the shear/effective-
#                            normal stress RESOLVED ONTO THE FAULT is needed
#                            (spec Part 1 bullet 2), computed below the same
#                            way test.tpv8/test.tpv30 do it. No gravity body
#                            force, no off-fault stress tensor, no Method-2
#                            machinery required (confirmed: this is the row
#                            149/TPV12-13 gap's NEGATIVE case).
#   run time                 0 - 13.0 s (spec Part 3, p.9); the e2e harness
#                            overrides par.term to GATE_TERM_S regardless
#                            (testsys/matrix.py), so this value only matters
#                            for a by-hand full-resolution run.
#   resolution                spec standard: 100 m and 50 m (Part 3, p.9).
#                            This case is built at par.dx=500 m as the FAST
#                            GATE (rule 17 step 4); see full_specs.py for the
#                            recorded (not run) full-resolution tier.

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -35.0e3, 35.0e3
par.ymin, par.ymax = -30.0e3, 30.0e3
par.zmin, par.zmax = -35.0e3, 0.0e3

par.fxmin, par.fxmax = -20.0e3, 20.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -20.0e3, 0.0e3

par.dx = 500.   # fast gate; full tier overrides to the spec value (see full_specs.py)
par.dy = par.dx
par.dz = par.dx

par.nmat = 1
par.vp, par.vs, par.rou = 6000., 3464., 2670.

par.term = 13.0

par.C_nuclea = 1
par.friclaw = 1
par.tpv = 26
par.nucR = 4.0e3   # r_crit, m

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

# Hypocenter: x = -5 km along strike, 10 km deep, on the (planar) fault.
par.xsource, par.ysource, par.zsource = -5.0e3, 0.0, -10.0e3

# --- TPV26/27 shared initial stress tensor, resolved onto the fault -------
# (spec Part 3 p.6-7; effective = total + Pf*I, Pf=1000*9.8*depth hydrostatic)
#   sigma22_eff = -(2670-1000)*9.8*depth      (vertical / EQdyna z-depth)
#   sigma11_eff = (Omega*b11 + (1-Omega))*sigma22_eff   (fault-parallel, x)
#   sigma33_eff = (Omega*b33 + (1-Omega))*sigma22_eff   (fault-normal, y)
#   sigma13     = Omega*b13*sigma22_eff                 (x-y shear, right-lateral)
#   Omega = 1 (depth<=15 km), (20000-depth)/5000 (15-20 km), 0 (>20 km)
# Fault is planar at y=0 with normal=(0,1,0), strike=(1,0,0), dip=(0,0,-1):
# normal stress = sigma33_eff, strike shear = sigma13, dip shear = 0 exactly
# (spec: sigma23=sigma12=0, Part 3 p.7).
B11, B33, B13 = 0.926793, 1.073206, -0.169029

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    depth = abs(zcoor)
    par.on_fault_vars[iz,ix,1] = par.fric_sw_fs
    # slip goes to zero at the fault border (left/right/bottom edges; the
    # free surface -- iz at fzmax -- is NOT a border, spec Part 3 p.5)
    if ix == 0 or ix == par.nfx-1 or iz == 0:
        par.on_fault_vars[iz,ix,1] = 10000.
    par.on_fault_vars[iz,ix,2] = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3] = par.fric_sw_D0
    # frictional cohesion C0 (spec Part 3 p.8): 4.00 MPa at surface, tapered
    # linearly to 0.40 MPa by 5000 m depth, constant below.
    if depth < 5000.:
        par.on_fault_vars[iz,ix,4] = 0.40e6 + 0.00072e6*(5000. - depth)
    else:
        par.on_fault_vars[iz,ix,4] = 0.40e6
    # forced-rupture decay time t0
    par.on_fault_vars[iz,ix,5] = 0.5

    # initial effective tractions, at a half-element depth for the
    # free-surface node (spec Part 7, Method 2 remark; same convention as
    # test.tpv30)
    dEff = max(depth, par.dx/2.)
    omega = 1.0 if dEff <= 15000. else max(0.0, (20000. - dEff)/5000.)
    s22 = -(par.rou - 1000.)*par.grav*dEff  # effective vertical stress
    syy = (omega*B33 + (1.-omega))*s22      # fault-normal
    sxy = omega*B13*s22                     # right-lateral driving shear

    par.on_fault_vars[iz,ix,7]  = syy   # initial effective normal
    par.on_fault_vars[iz,ix,8]  = sxy   # initial strike shear
    par.on_fault_vars[iz,ix,49] = 0.0   # initial dip shear (sigma23=0 per spec)

    par.on_fault_vars[iz,ix,46] = 0.0   # no creep/background slip rate

# --- Stations (spec Parts 8-9, p.21-28) -------------------------------------
# on-fault: (x, -depth) in km, spec along-strike/down-dip coordinates
par.st_coor_on_fault = [
    [-5.0,0.0], [15.0,0.0], [-5.0,-5.0], [15.0,-5.0], [-15.0,-10.0],
    [-5.0,-10.0], [0.0,-10.0], [5.0,-10.0], [10.0,-10.0], [15.0,-10.0],
    [-5.0,-15.0], [15.0,-15.0]]
# off-fault: (x, y_eq, z_eq) in km; y_eq = TPV fault-normal coordinate
# (positive = far side), all stations at the surface
par.st_coor_off_fault = [
    [-5.0,3.0,0.0], [-5.0,-3.0,0.0], [5.0,3.0,0.0],
    [5.0,-3.0,0.0], [15.0,3.0,0.0], [15.0,-3.0,0.0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
