#! /usr/bin/env python3

# SCEC TPV35: Parkfield 2004 M6 VALIDATION benchmark (strike.scec.org/cvws,
# TPV35_Description_v05). Vertical right-lateral planar fault, 40 km x
# 15.5 km, x in [-30, 10] km, reaching the free surface, in a two-sided 1D
# layered elastic medium (a different 1D profile on each side of the fault),
# linear slip-weakening friction. Yield stress mu_s(x,z) and initial shear
# stress tau0(x,z) are SUPPLIED on a 100 m grid (Ma, Custodio, Archuleta &
# Liu 2008 inversion) and nucleation is built into them: around
# (x=0, depth 8.1 km) mu_s*sigma_n < tau0, so the fault fails spontaneously
# at t=0 with no artificial nucleation (par.C_nuclea = 0).
#
# Spec parameters encoded here (see tpv35Tools.py for the data files):
#   material          two 1D tables, near side (y<0) 8 layers + halfspace,
#                     far side (y>0) 9 layers + halfspace, min Vs 1100 m/s,
#                     halfspace Vp/Vs/rho 7300/4300/2800 (par.nmat/n2mat==5)
#   friction (SW)     mu_s from data; mu_d 0.30; d0 0.15 m; C0 0
#   initial stress    sigma_n 60 MPa constant; tau0 from data, pure
#                     right-lateral (positive strike shear); no gravity
#   borders           slip zero on the side/bottom edges: mu_s = 1000 there,
#                     the same device tpv8/tpv29 use
#   stations          7 on-fault at 8.1 km depth, x = -25..+5 km step 5;
#                     43 off-fault surface stations from the official list
#   run time          0 - 18 s (spec); par.term below is the gate term
#   resolution        spec 100 m (50 m optional). par.dx below is the gate
#                     value. Only 100 and 500 m both tile the fault and are
#                     integer multiples of the shipped grid; any other dx is
#                     refused by lib.requireFaultGeometryResolution.

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np
import tpv35Tools as tools

par = parameters()

par.xmin, par.xmax = -40.0e3, 20.0e3
par.ymin, par.ymax = -16.0e3, 16.0e3
par.zmin, par.zmax = -26.0e3, 0.0e3

par.fxmin, par.fxmax = -30.0e3, 10.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -15.5e3, 0.0e3

# Hypocenter (spec Part 3): x = 0, 8.1 km deep. Informational only here --
# nucleation comes from the data, C_nuclea = 0.
par.xsource, par.ysource, par.zsource = 0.0, 0.0, -8.1e3

par.dx = 500.
par.dy = par.dx
par.dz = par.dx

# Two-sided 1D velocity structure: n2mat == 5 columns
# [bottom depth, vp, vs, rho, side]; side -1 is y < fault plane (spec near
# side), +1 is y > fault plane (far side). Halfspace closed at 1e6 m.
par.mat = tools.materialTable(1.0e6)
par.nmat, par.n2mat = par.mat.shape
par.roumax = par.mat[:, 3].max()
par.vmaxPML = par.mat[:, 1].max()
par.vp, par.vs, par.rou = 7300., 4300., 2800.   # halfspace; not used for nmat>1

par.term = 5.

par.dt = 0.5*par.dx/par.vmaxPML
par.friclaw = 1
par.tpv = 35
par.C_nuclea = 0
par.faultStNormalStressSign = 'extension'   # spec: "Positive means extension."

par.nx = 2
par.ny = 1
par.nz = 2

# Fault grid
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)

# Declared so scripts/case.setup refuses a dx the shipped grid cannot give.
par.faultGeometrySourceDx = tools.SHIPPED_DX
par.faultGeometrySourceName = tools.SHIPPED_NAME
par.faultGeometrySourceAvailableDx = [tools.SHIPPED_DX]

# Official mu_s / tau0 on the case grid, exact decimation ([iz, ix], iz=0 at fzmin).
mu_s, tau0 = tools.faultGridForCase(par.dx, par.fx, par.fz)

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
par.fric_sw_fd = tools.MU_D
par.fric_sw_D0 = tools.D0
par.fric_cohesion = tools.COHESION
par.grav = 0.0
par.init_norm = -tools.NORMAL_STRESS
for ix, xcoor in enumerate(par.fx):
    for iz, zcoor in enumerate(par.fz):
        par.on_fault_vars[iz, ix, 1] = mu_s[iz, ix]
        if (abs(xcoor - par.fxmin) < 0.01 or abs(xcoor - par.fxmax) < 0.01
                or abs(zcoor - par.fzmin) < 0.01):
            par.on_fault_vars[iz, ix, 1] = 1000.      # no slip on the borders
        par.on_fault_vars[iz, ix, 2] = par.fric_sw_fd
        par.on_fault_vars[iz, ix, 3] = par.fric_sw_D0
        par.on_fault_vars[iz, ix, 4] = par.fric_cohesion
        par.on_fault_vars[iz, ix, 7] = par.init_norm   # negative compressive
        par.on_fault_vars[iz, ix, 8] = tau0[iz, ix]     # right-lateral positive

# Stations (km). On-fault: spec Part 3, depth 8.1 km, x = -25..+5 km step
# 5. setOnFaultStation (meshgen.f90) requires x AND z to EQUAL a fault node
# (it does not snap, unlike off-fault stations), so the depth is rounded to
# the nearest fault-grid plane for THIS dz: exactly 8.1 km at the 100 m spec
# resolution, 8.0 km at the 500 m gate (file names dp081 / dp080).
stDepthKm = round(8.1e3/par.dz)*par.dz/1.0e3
par.st_coor_on_fault = [[x, -stDepthKm] for x in (-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0)]
par.st_coor_off_fault = tools.offFaultStationsKm()
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
