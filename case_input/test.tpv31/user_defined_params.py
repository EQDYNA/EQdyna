#! /usr/bin/env python3

# SCEC TPV31: right-lateral vertical strike-slip PLANAR fault in a linear
# elastic half-space with a DISCONTINUOUS 1D velocity structure
# (strike.scec.org/cvws, TPV31_32_Description_v03, scratch/specs/).
#
# Fault: 30 km along-strike x 15 km deep, planar, vertical, reaches the
# surface. Hypocenter at the center, 7.5 km deep -- spec (x,y,z) =
# (0, 7500, 0), where spec y=depth, z=fault-normal; EQdyna (x=strike,
# y=fault-normal, z=depth negative down): par.xsource=0.0, par.ysource=0.0,
# par.zsource=-7.5e3 (Part 1, p.3).
#
# This compset follows test.tpv26's pattern (planar fault, no rough-fault
# geometry, Method-2-style on_fault_vars stress specification, static
# nucleation patch instead of swtwNucleation): the spec's nucleation
# (Part 2 "Nucleation Shear Stress", p.7) is an additive, TIME-INDEPENDENT
# shear-stress bump near the hypocenter, not a forced-rupture-time formula,
# so par.tpv is left at an UNMATCHED value (31) and par.C_nuclea=0 -- see
# faulting.f90:swtwNucleation, which is a no-op for any unmatched par.tpv
# (confirmed by direct read, rule 17 step 3, 2026-10-06; same conclusion as
# test.tpv8/test.tpv10's static-patch cases).
#
# Spec parameters encoded here (TPV31_32_Description_v03):
#   material (Part 2, p.4)   DISCONTINUOUS 1D velocity structure, 3 jumps at
#                            depth 2400/5000/10000 m, piecewise-constant
#                            between jumps with LINEAR ramps in between (the
#                            table literally gives "2400-"/"2400+" duplicate
#                            rows -- a genuine piecewise-linear table with
#                            co-located jump points, reproduced exactly below
#                            via np.interp on the verbatim table, duplicate
#                            x-values included; no case grid point ever
#                            queries exactly at a jump depth since element/
#                            fault-node centers and nodes fall at multiples
#                            of dx=500 or (dx/2)-offset bin centers, none of
#                            which equal 2400/5000/10000).
#                            Reused via EQdyna's EXISTING n2mat==4 (1D
#                            layered, piecewise-constant-per-element)
#                            mechanism -- already live/gated via
#                            test.meng2023a -- by precomputing one row per
#                            par.dz mesh layer, each assigned the spec's
#                            exact np.interp-evaluated value at that layer's
#                            CENTER depth. No src/ changes needed.
#   initial stress (Part 2,  sigma22=0; sigma11=sigma33=(-60.00 MPa)(mu/mu0);
#     p.6-7)                 sigma13=(30.00 MPa)(mu/mu0); sigma23=sigma12=0.
#                            mu0=32.03812032 GPa (TPV5's reference modulus).
#                            Fault is the z=0 plane, normal=(0,0,1) in spec
#                            axes -> on_fault_vars normal=sigma33 (negative
#                            compressive, same convention as test.tpv26/29),
#                            shear=sigma13, dip shear=sigma23=0 (spec
#                            Stress Tensor Components table, p.6): this is
#                            SCEC's own "Method 2" (store only the on-fault
#                            traction, no off-fault stress tensor needed for
#                            a purely elastic case, Part 4) -- the SAME
#                            on_fault_vars[...,7]/[...,8] mechanism test.tpv8/
#                            test.tpv26/test.tpv29 already use.
#   nucleation shear stress  tau_nuke(r) = (4.95 MPa)(mu/mu0) if r<=1400 m;
#     (Part 2, p.7)          (2.475 MPa)(1+cos(pi*(r-1400)/600))(mu/mu0) if
#                            1400<=r<=2000 m; 0 otherwise. r = 2D distance to
#                            hypocenter in the fault plane. mu is evaluated
#                            at the LOCAL node depth (consistent with the
#                            spec's own worked check: at the hypocenter,
#                            depth=7500 m, total shear = 34.95 MPa*(mu/mu0),
#                            i.e. the SAME mu/mu0 factor for both the
#                            regional 30 MPa term and the 4.95 MPa nucleation
#                            term). ADDED directly to on_fault_vars[...,8],
#                            same style as test.tpv10's static patch.
#   friction (Part 2, p.8)   mu_s=0.580, mu_d=0.450, d0=0.18 m (linear
#                            slip-weakening, shared by TPV31 and TPV32);
#                            frictional cohesion C0 = 0.000425 MPa/m *
#                            (2400 m - depth) if depth<=2400 m, else 0
#                            (1.02 MPa at the surface, tapering to 0 by
#                            2400 m depth).
#   run time / resolution    0 - 15.0 s after nucleation (Part 2, p.9); spec
#                            standard resolution for TPV31 is 50 m node
#                            spacing (mandatory, not a range like TPV32).
#                            This case is built at par.dx=500 m as the FAST
#                            GATE (rule 17 step 4); see full_specs.py for the
#                            recorded (not run) 50 m full-resolution tier.
#   geometry (Part 1, p.3)   planar fault, no decimated/rough geometry file
#                            needed (rule 17 step 2: N/A for this case).

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

# --- TPV31 discontinuous 1D velocity structure (spec Part 2, p.4, verbatim) -
# depth (m), duplicate rows at the 3 jump depths reproduce the spec's exact
# "2400-"/"2400+" discontinuity-with-ramp table.
_depth = [0., 2400., 2400., 5000., 5000., 10000., 10000., 15000.]
_vp    = [4050., 4050., 4450., 5200., 5750., 5750., 6500., 6500.]
_vs    = [2250., 2250., 2550., 3050., 3450., 3450., 3800., 3800.]
_rho   = [2580., 2580., 2600., 2620., 2720., 2720., 3000., 3000.]

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
par.tpv = 31        # UNMATCHED in faulting.f90:swtwNucleation -> confirmed no-op

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
