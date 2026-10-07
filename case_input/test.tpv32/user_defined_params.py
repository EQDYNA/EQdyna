#! /usr/bin/env python3

# SCEC TPV32: right-lateral vertical strike-slip PLANAR fault in a linear
# elastic half-space with a CONTINUOUS, piecewise-linear 1D velocity
# structure (strike.scec.org/cvws, TPV31_32_Description_v03,
# scratch/specs/). Same fault geometry, stress tensor, nucleation, friction
# and station set as test.tpv31 -- see that case's header for the shared
# derivation. The ONLY difference from TPV31 is the velocity/density table
# (Part 2, p.5): TPV32's is continuous (no jump rows) and uses a different
# depth grid.
#
# Spec parameters NOT shared with TPV31:
#   material (Part 2, p.5)   CONTINUOUS piecewise-linear 1D velocity
#                            structure at depths
#                            [0,500,1000,1600,2400,3600,5000,9000,11000,
#                            15000] m. Reused via the same n2mat==4
#                            mechanism as test.tpv31 (no src/ changes).
#   run time / resolution    spec allows TPV32 node spacing in [25,50] m
#     (Part 2, p.9)          (mandatory 50 m for TPV31; TPV32 gets a range
#                            because the spec's own preliminary tests needed
#                            25 m near the low-velocity surface layer to get
#                            acceptable results). This case records 50 m (the
#                            upper, coarser end of that allowed range) as its
#                            full-resolution tier in full_specs.py -- a
#                            scheduling choice (rule 17 step 5), not a claim
#                            that 50 m is the spec's preferred value.

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -30.0e3, 30.0e3
par.ymin, par.ymax = -30.0e3, 30.0e3
par.zmin, par.zmax = -30.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -15.0e3, 0.0e3

par.dx = 500.   # fast gate; full tier overrides to the spec's 50 m (see full_specs.py)
par.dy = par.dx
par.dz = par.dx

# --- TPV32 continuous, piecewise-linear 1D velocity structure (spec Part 2,
# p.5, verbatim; no duplicate/jump rows -- unlike TPV31, truly continuous).
_depth = [0., 500., 1000., 1600., 2400., 3600., 5000., 9000., 11000., 15000.]
_vp    = [2200., 3000., 3600., 4400., 4800., 5250., 5500., 5750., 6100., 6300.]
_vs    = [1050., 1400., 1950., 2500., 2800., 3100., 3250., 3450., 3600., 3700.]
_rho   = [2200., 2450., 2550., 2600., 2600., 2620., 2650., 2720., 2750., 2900.]

_zDomainDepth = abs(par.zmin)   # 30000 m: material table must cover the WHOLE volume
_nLayers = round(_zDomainDepth/par.dz)
par.nmat = _nLayers + 1   # +1: catch-all row beyond the domain (meng2023a convention)
par.n2mat = 4
par.mat = np.zeros((par.nmat, par.n2mat))
for i in range(_nLayers):
    centerDepth = (i + 0.5)*par.dz
    bottomBound = (i + 1)*par.dz
    par.mat[i,0] = bottomBound
    par.mat[i,1] = np.interp(centerDepth, _depth, _vp)
    par.mat[i,2] = np.interp(centerDepth, _depth, _vs)
    par.mat[i,3] = np.interp(centerDepth, _depth, _rho)
par.mat[_nLayers,:] = [1.e6, _vp[-1], _vs[-1], _rho[-1]]   # catch-all, deepest value

par.vp, par.vs, par.rou = _vp[-1], _vs[-1], _rho[-1]   # reference (stable-timestep) values

par.term = 15.0

par.C_elastic = 1
par.C_nuclea = 0   # static nucleation patch below, not swtwNucleation (rule 17 step 3)
par.friclaw = 1
par.tpv = 32        # UNMATCHED in faulting.f90:swtwNucleation -> confirmed no-op

par.dt = 0.5*par.dx/_vp[-1]

par.grav = 9.8
par.fric_sw_fs = 0.580
par.fric_sw_fd = 0.450
par.fric_sw_D0 = 0.18

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

# Hypocenter: center of the fault, 7.5 km deep (spec Part 1, p.3).
par.xsource, par.ysource, par.zsource = 0.0, 0.0, -7.5e3

mu0 = 32.03812032e9   # spec Part 2 p.6-7: TPV5's reference shear modulus

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    depth = abs(zcoor)
    vs_local  = np.interp(depth, _depth, _vs)
    rho_local = np.interp(depth, _depth, _rho)
    mu_local  = rho_local*vs_local**2
    muRatio   = mu_local/mu0

    par.on_fault_vars[iz,ix,1] = par.fric_sw_fs
    # slip goes to zero at the fault border (left/right/bottom edges; the
    # free surface -- iz at fzmax -- is NOT a border, spec Part 1 p.3)
    if ix == 0 or ix == par.nfx-1 or iz == 0:
        par.on_fault_vars[iz,ix,1] = 10000.
    par.on_fault_vars[iz,ix,2] = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3] = par.fric_sw_D0
    # frictional cohesion C0 (spec Part 2, p.8)
    if depth <= 2400.:
        par.on_fault_vars[iz,ix,4] = 0.000425e6*(2400. - depth)
    else:
        par.on_fault_vars[iz,ix,4] = 0.0

    # initial effective tractions from the depth-dependent stress tensor
    # (spec Part 2 p.6-7): sigma33=sigma11=-60 MPa*(mu/mu0), sigma13=30
    # MPa*(mu/mu0); on_fault_vars normal=sigma33 (negative compressive),
    # strike shear=sigma13, dip shear=sigma23=0.
    syy = -60.0e6*muRatio       # initial effective normal (=sigma33, negative compressive)
    sxy = 30.0e6*muRatio        # initial strike shear (=sigma13)

    # nucleation shear stress (spec Part 2 p.7): r = 2D distance to
    # hypocenter IN THE FAULT PLANE (spec x, spec y=depth)
    r = sqrt(xcoor**2 + (depth - 7500.)**2)
    if r <= 1400.:
        tauNuke = 4.95e6*muRatio
    elif r <= 2000.:
        tauNuke = 2.475e6*(1. + cos(pi*(r - 1400.)/600.))*muRatio
    else:
        tauNuke = 0.0

    par.on_fault_vars[iz,ix,7]  = syy
    par.on_fault_vars[iz,ix,8]  = sxy + tauNuke   # regional + nucleation bump
    par.on_fault_vars[iz,ix,49] = 0.0             # initial dip shear (sigma23=0)

    par.on_fault_vars[iz,ix,46] = 0.0   # no creep/background slip rate

# --- Stations (spec Parts 5-6, p.12-21) -------------------------------------
# on-fault: (x, -depth) in km, spec along-strike/down-dip coordinates
# (verbatim from Part 5's 30-station table, p.12-13)
_onStrike = [0.0, 6.0, 12.0]
_onDepth  = [0.0, 0.2, 0.5, 1.0, 2.4, 3.0, 5.0, 7.5, 10.0, 12.0]
par.st_coor_on_fault = [[s, -d] for s in _onStrike for d in _onDepth]

# off-fault: (x, y_eq, z_eq) in km; y_eq = TPV fault-normal coordinate
# (positive = far side); verbatim from Part 6's 18-station table, p.18
# (6 "boreholes", each with 3 depths: 0, 0.5, 2.4 km)
_offBoreholes = [(0.0, -3.0), (12.0, -3.0), (0.0, -9.0),
                 (15.0, -9.0), (0.0, -15.0), (15.0, -15.0)]
_offDepth = [0.0, 0.5, 2.4]
par.st_coor_off_fault = [[s, off, -d] for s, off in _offBoreholes for d in _offDepth]

par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
