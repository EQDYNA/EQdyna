#! /usr/bin/env python3

# SCEC TPV33: right-lateral vertical strike-slip PLANAR fault in a linear
# elastic half-space, centered in a 1.6 km thick LOW-VELOCITY FAULT ZONE
# (fault zone guided waves), different velocities on each outer side
# (strike.scec.org/cvws, TPV33_Description_v04, scratch/specs/).
#
# Fault: 16 km along-strike (-12000 m <= x <= 4000 m) x 10 km deep, planar,
# vertical, reaches the surface. Hypocenter 6 km from the left edge, 6 km
# deep -- spec (x,y,z) = (-6000, 6000, 0), where spec y=depth, z=fault-normal;
# EQdyna (x=strike, y=fault-normal, z=depth negative down), same mapping
# TPV34's own README documents (y_eq=z_tpv, z_eq=-y_tpv):
# par.xsource=-6.0e3, par.ysource=0.0, par.zsource=-6.0e3 (Part 1, p.3).
#
# This compset follows test.tpv31's pattern (planar fault, no rough-fault
# geometry, Method-2-style on_fault_vars stress specification, static
# nucleation patch instead of swtwNucleation): the spec's nucleation
# (Part 2 "Initial Stress", p.5) is an additive, TIME-INDEPENDENT shear-
# stress bump near the hypocenter, not a forced-rupture-time formula, so
# par.tpv is left at an UNMATCHED value (33) and par.C_nuclea=0 --
# faulting.f90:swtwNucleation is a no-op for any unmatched par.tpv (direct
# read, rule 17 step 3, confirmed 2026-10-07; same conclusion already
# reached for test.tpv8/test.tpv10/test.tpv31/test.tpv34).
#
# Rule 17 step 3 (board row 152): does the EXISTING n2mat==6 3D material
# grid (SCEC TPV34, meshgen.f90/meshgen.py setElementMaterial, nearest-cell,
# never interpolated) already cover a thin FAULT-PARALLEL low-velocity
# zone, or does it need a new material-grid mode? CONFIRMED: reuse, no
# src/ changes. TPV33's velocity structure (spec Part 2 p.4) depends ONLY
# on the spec's fault-normal coordinate z -- which maps directly onto
# EQdyna's y axis (no sign flip, same mapping TPV34's extract_cvmh_grid.py
# documents) -- and is CONSTANT along the spec's x (EQdyna x) and y/depth
# (EQdyna z). n2mat==6 rows are self-describing [x y z vp vs rho] grid
# cells read by NEAREST-CELL lookup per axis independently, so a grid that
# is uniform-spaced (= par.dx) along EQdyna y and carries just the TWO BOX
# EDGES along EQdyna x and z reproduces this exactly: every element's
# lookup along x/z always resolves to one of the two (identical-valued)
# edge samples, giving true x/z-independence, and the y lookup reproduces
# the spec's three bands at whatever precision par.dx (the gate's chosen
# resolution) gives -- same nearest-cell discretization principle as
# TPV34's dense 3-axis CVM-H grid, just exploiting that two of this
# case's three axes carry no material variation at all. (readInputFiles.f90
# buildMaterialGrid3D / checkInputConsistency.py build_material_grid3d's
# grid-covers-mesh invariant, PR #114/item 147, is satisfied the same way
# TPV34 satisfies it: the grid is built to exactly tile the declared mesh
# box at par.dx, not an independent/smaller source grid.)
#
# Spec parameters encoded here (TPV33_Description_v04):
#   material (Part 2, p.4)   rho=2670 kg/m^3 everywhere; Vs/Vp depend only
#                            on fault-normal distance z: 3248/5626 m/s if
#                            z<-800 m, 2165/3750 m/s (LVFZ) if -800<z<800 m,
#                            3464/6000 m/s if z>800 m. n2mat==6 grid above.
#   initial stress (Part 2,  sigma0(x,y)=60 MPa (constant, NOT mu/mu0-scaled
#     p.5)                   -- TPV33 has no material-scaling of stress,
#                            unlike TPV31/34/35); tau0(x,y) = 30 MPa*(1-Rtau)
#                            + tau_nuke(x,y), Rtau/Rx/Ry and tau_nuke exactly
#                            as the spec's own formulas (reproduced
#                            verbatim below, including the taper breakpoints
#                            -9800/1100 m along strike and 2300/8000 m depth).
#                            No gravity ("There is no gravity in the
#                            model.", Part 1 p.1) -- par.grav=0.0, same as
#                            test.tpv34.
#   nucleation shear stress  tau_nuke(r) = 3.150 MPa if r<=550 m;
#     (Part 2, p.5)          1.575 MPa*(1+cos(pi*(r-550)/250)) if
#                            550<=r<=800 m; 0 otherwise. r = 2D distance to
#                            hypocenter in the fault plane (spec x, depth).
#                            Worked check (spec p.6): total initial shear at
#                            the hypocenter = 33.15 MPa, yield = 33.00 MPa.
#   friction (Part 3, p.8)   mu_s=0.550, mu_d=0.450, d0=0.18 m; C0=0 (no
#                            frictional cohesion anywhere, unlike TPV31/34).
#   run time / resolution    0 - 13.0 s after nucleation (Part 2, p.7); spec
#                            resolution request for TPV33 is 12.5-25 m on
#                            the fault plane (mandatory range, narrower than
#                            TPV31/32's 25-50 m, because TPV33 requires
#                            twice the resolution of any earlier benchmark
#                            -- Part 2 p.7), 50 m through the low-velocity
#                            zone, 100 m outside it. This case is built at
#                            par.dx=400 m as the FAST GATE (rule 17 step 4);
#                            see full_specs.py for the recorded (not run)
#                            tier.
#   geometry (Part 1, p.3)   planar fault, no decimated/rough geometry file
#                            needed (rule 17 step 2: N/A for this case).
#
# Validation (rule 17 step 6): no independent evidence_tpv33_*.py script is
# shipped, the same documented choice already made for test.tpv26/test.tpv27/
# test.tpv31/test.tpv32 (no evidence_* script for any of those four either,
# see testsys/parity/README.md) -- this case's cross-backend fortran/jax
# agreement (CASE_BOUND/STATION_BOUND below) is the comparison it provides;
# an independent cross-code check against an EQdyna SCEC submission is left
# open the same way it is for those four cases, not silently claimed done.

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -26.0e3, 18.0e3
par.ymin, par.ymax = -22.0e3, 22.0e3
par.zmin, par.zmax = -20.0e3, 0.0e3

par.fxmin, par.fxmax = -12.0e3, 4.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -10.0e3, 0.0e3

par.dx = 400.   # fast gate; full tier recorded (not run) in full_specs.py
par.dy = par.dx
par.dz = par.dx

# --- n2mat==6 3D material grid (fault-normal low-velocity zone, spec Part 2
# p.4, verbatim thresholds +-800 m) -- see the rule-17-step-3 note above for
# why this reuses the EXISTING grid mechanism with no src/ change. Rows are
# constant along EQdyna x/z (just the two box edges) and uniformly spaced
# at par.dx along EQdyna y (the spec's fault-normal z), so every element's
# nearest-cell lookup along x/z always lands on one of the two
# identical-valued edge samples (true x/z independence) while the y lookup
# resolves the three bands at this gate's own par.dx.
_RHO = 2670.

def _vs_of_y(y):
    if y < -800.:
        return 3248.
    elif y > 800.:
        return 3464.
    else:
        return 2165.

def _vp_of_y(y):
    if y < -800.:
        return 5626.
    elif y > 800.:
        return 6000.
    else:
        return 3750.

_nYGrid = round((par.ymax - par.ymin)/par.dx)
_yCenters = par.ymin + par.dx*(np.arange(_nYGrid) + 0.5)
_xEdges = [par.xmin, par.xmax]
_zEdges = [par.zmin, par.zmax]
_rows = []
for _xv in _xEdges:
    for _yv in _yCenters:
        for _zv in _zEdges:
            _rows.append([_xv, _yv, _zv, _vp_of_y(_yv), _vs_of_y(_yv), _RHO])
par.mat = np.array(_rows)
par.nmat, par.n2mat = par.mat.shape
par.roumax = par.mat[:, 5].max()
par.vmaxPML = par.mat[:, 3].max()
par.vp, par.vs, par.rou = 6000., 3464., 2670.   # host-rock reference; not used for nmat>1

par.term = 13.0

par.C_elastic = 1
par.C_nuclea = 0   # static nucleation patch below, not swtwNucleation (rule 17 step 3)
par.friclaw = 1
par.tpv = 33        # UNMATCHED in faulting.f90:swtwNucleation -> confirmed no-op

par.dt = 0.5*par.dx/par.vmaxPML

par.grav = 0.0   # spec Part 1 p.1: there is no gravity in the model
par.fric_sw_fs = 0.550
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

# Hypocenter: 6 km from the left edge, 6 km deep (spec Part 1, p.3).
par.xsource, par.ysource, par.zsource = -6.0e3, 0.0, -6.0e3

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    depth = abs(zcoor)

    par.on_fault_vars[iz,ix,1] = par.fric_sw_fs
    # slip goes to zero at the fault border (left/right/bottom edges; the
    # free surface -- iz at fzmax -- is NOT a border, spec Part 1 p.3)
    if ix == 0 or ix == par.nfx-1 or iz == 0:
        par.on_fault_vars[iz,ix,1] = 10000.
    par.on_fault_vars[iz,ix,2] = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3] = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4] = 0.0   # frictional cohesion C0 = 0 (spec Part 3, p.8)

    # initial stress (spec Part 2, p.5): sigma0 constant, tau0 tapered with
    # an additive nucleation bump. Rx/Ry/Rtau/tau_nuke reproduced verbatim.
    sigma0 = 60.0e6
    if xcoor < -9800.:
        Rx = (-xcoor - 9800.)/10000.
    elif xcoor > 1100.:
        Rx = (xcoor - 1100.)/10000.
    else:
        Rx = 0.0
    if depth < 2300.:
        Ry = (-depth + 2300.)/10000.
    elif depth > 8000.:
        Ry = (depth - 8000.)/10000.
    else:
        Ry = 0.0
    Rtau = min(1.0, sqrt(Rx**2 + Ry**2))

    r = sqrt((xcoor + 6000.)**2 + (depth - 6000.)**2)
    if r <= 550.:
        tauNuke = 3.150e6
    elif r <= 800.:
        tauNuke = 1.575e6*(1. + cos(pi*(r - 550.)/250.))
    else:
        tauNuke = 0.0

    tau0 = 30.0e6*(1. - Rtau) + tauNuke

    par.on_fault_vars[iz,ix,7]  = -sigma0   # negative compressive
    par.on_fault_vars[iz,ix,8]  = tau0      # right-lateral positive
    par.on_fault_vars[iz,ix,49] = 0.0       # initial dip shear (sigma23=0)

    par.on_fault_vars[iz,ix,46] = 0.0   # no creep/background slip rate

# --- Stations (spec Parts 4-5, p.9-19) --------------------------------------
# on-fault: all 28 spec stations (Part 4, p.9-10), (along-strike km,
# down-dip km) -- verbatim 7 strikes x 4 depths.
_onStrike = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0]
_onDepth  = [2.0, 4.0, 6.0, 8.0]
par.st_coor_on_fault = [[s, -d] for s in _onStrike for d in _onDepth]

# off-fault: a VERIFIED SUBSET of the spec's 85 off-fault stations (Part 5,
# p.15-19) -- the four EARTH-SURFACE (depth 0) dense transects (strikes
# -6/4/8/12 km, each 9 stations at fault-normal offsets -1.6..1.6 km step
# 0.4 km, spec's own coordinates verbatim) plus the 12 "other stations
# surrounding the fault" (strikes -12/-4/4 km, offsets +-5/+-10 km, also
# surface, also verbatim). The spec's depth-6-km transects (at along-strike
# 0/4/8/12 km, offsets including an asymmetric +-0.1 km pair at the fault
# trace) are NOT transcribed here -- a documented gap (not invented data),
# left open the same way test.tpv26/27/31/32 leave rule 17 step 6 open; none
# of matrix.GATE_STATIONS' off-fault picks need them.
_offStrikesSurf = [-6.0, 4.0, 8.0, 12.0]
_offOffsets9 = [-1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6]
par.st_coor_off_fault = [[s, off, 0.0] for s in _offStrikesSurf for off in _offOffsets9]
_otherStrikes = [-12.0, -4.0, 4.0]
_otherOffsets = [-10.0, -5.0, 5.0, 10.0]
par.st_coor_off_fault += [[s, off, 0.0] for s in _otherStrikes for off in _otherOffsets]

par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
