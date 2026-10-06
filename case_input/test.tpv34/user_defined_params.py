#! /usr/bin/env python3

# SCEC TPV34: Imperial Fault, Model 1 (strike.scec.org/cvws,
# TPV34_Description_v10). Vertical right-lateral planar fault, 30 km x
# 15 km, x in [-15, 15] km, reaching the free surface, in the 3D CVM-H
# velocity structure (Imperial Valley), linear slip-weakening friction.
# Initial stresses scale with the local shear modulus, tau0 = 30 MPa *
# mu/mu0 and sigma0 = 60 MPa * mu/mu0, and nucleation is an ADDED shear
# stress (4.95 MPa * mu/mu0 in a 1.4 km radius, cosine taper to 2 km)
# around (x = 0, depth 7.5 km) at t = 0 -- so the fault fails spontaneously
# with no forced-rupture-time nucleation (par.C_nuclea = 0, no swtwNucleation
# TPV34 branch is needed or exists).
#
# Spec parameters encoded here (see tpv34Tools.py and README.md):
#   material          CVM-H sampled at the element centres of THIS box at
#                     the gate dx (tpv34_cvmh_grid_500m.txt.gz, extracted by
#                     extract_cvmh_grid.py per spec Part 2, min Vp/Vs clamp
#                     2984/1400 m/s); par.n2mat == 6 3D grid, nearest-cell
#   friction (SW)     mu_s 0.580; mu_d 0.450; d0 0.18 m; C0 = 0.000425 MPa/m
#                     * (2400 m - depth) above 2400 m depth
#   initial stress    tau0/sigma0 proportional to mu/mu0 at each fault node
#                     (mu = mean of the adjacent elements on both sides);
#                     pure right-lateral (positive strike shear); no gravity
#   borders           slip zero on the side/bottom edges: mu_s = 1000 there,
#                     the same device tpv8/tpv29/tpv35 use
#   stations          35 on-fault (5 along strike x 7 depths), 56 off-fault
#   run time          0 - 20 s (spec); par.term below is the gate term
#   resolution        spec: 25 - 50 m, submitter's choice (recorded as
#                     EXCLUDED in testsys/e2e/full_specs.py). par.dx below is
#                     the gate value and the only dx the shipped grid admits.

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np
import tpv34Tools as tools

par = parameters()

par.xmin, par.xmax = -30.0e3, 30.0e3
par.ymin, par.ymax = -16.0e3, 16.0e3
par.zmin, par.zmax = -26.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -15.0e3, 0.0e3

# Hypocentre (spec Part 3): x = 0, 7.5 km deep. Informational only here --
# nucleation is in the initial stress, C_nuclea = 0.
par.xsource, par.ysource, par.zsource = 0.0, 0.0, -7.5e3

par.dx = 500.
par.dy = par.dx
par.dz = par.dx

# 3D CVM-H material grid: n2mat == 6 columns [x, y, z, vp, vs, rho] on the
# element centres of the box above at par.dx (extract_cvmh_grid.py's
# default box MUST equal the box above).
par.mat = tools.materialGrid(par.dx)
par.nmat, par.n2mat = par.mat.shape
par.roumax = par.mat[:, 5].max()
par.vmaxPML = par.mat[:, 3].max()
par.vp, par.vs, par.rou = 2984., 1400., 2220.34   # clamp floor; not used for nmat>1

par.term = 5.

par.dt = 0.5*par.dx/par.vmaxPML
par.friclaw = 1
par.tpv = 34
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

# Declared so scripts/case.setup refuses a dx the shipped grid was not
# sampled at.
par.faultGeometrySourceDx = tools.SHIPPED_DX
par.faultGeometrySourceName = tools.SHIPPED_NAME
par.faultGeometrySourceAvailableDx = [tools.SHIPPED_DX]

# mu/mu0, tau0 (with nucleation), sigma0 and C0 on the fault grid [iz, ix].
mu_ratio, tau0, sigma0, C0 = tools.faultState(par.mat, par.dx, par.fx, par.fz)

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
par.fric_sw_fs = tools.MU_S
par.fric_sw_fd = tools.MU_D
par.fric_sw_D0 = tools.D0
par.grav = 0.0
for ix, xcoor in enumerate(par.fx):
    for iz, zcoor in enumerate(par.fz):
        par.on_fault_vars[iz, ix, 1] = par.fric_sw_fs
        if (abs(xcoor - par.fxmin) < 0.01 or abs(xcoor - par.fxmax) < 0.01
                or abs(zcoor - par.fzmin) < 0.01):
            par.on_fault_vars[iz, ix, 1] = 1000.      # no slip on the borders
        par.on_fault_vars[iz, ix, 2] = par.fric_sw_fd
        par.on_fault_vars[iz, ix, 3] = par.fric_sw_D0
        par.on_fault_vars[iz, ix, 4] = C0[iz, ix]
        par.on_fault_vars[iz, ix, 7] = -sigma0[iz, ix]  # negative compressive
        par.on_fault_vars[iz, ix, 8] = tau0[iz, ix]     # right-lateral positive

# Stations (km). On-fault depths are rounded to the nearest fault plane for
# THIS dz (2.4 km -> 2.5 km at the 500 m gate; setOnFaultStation needs an
# exact node). Off-fault stations are the spec's own coordinates (they snap).
par.st_coor_on_fault = tools.onFaultStationsKm(par.dz)
par.st_coor_off_fault = tools.offFaultStationsKm()
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
