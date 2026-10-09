#! /usr/bin/env python3

# SCEC TPV13: 60-degree dipping, planar, normal fault in a uniform
# non-associative Drucker-Prager PLASTIC half-space (strike.scec.org/cvws,
# TPV12_13_Description_v6.pdf, 2009-10-01). TPV13 is the plastic sibling of
# test.tpv12 (pathway_forward.md row 149); owner decision 2026-10-07 pairs
# TPV12/13 the same way TPV29/TPV30 are paired: TPV12 uses METHOD 1 (stress
# change, no gravity, C_elastic=1, per-node fault tractions set manually
# below). TPV13 uses METHOD 2 (explicit gravity + full off-fault stress
# tensor, C_elastic=0), as test.tpv30 does. Spec Part 2, p.4: "The choice of
# material properties is the only difference between TPV12 and TPV13" --
# geometry, friction, nucleation, on-fault traction formulas and stations
# below are IDENTICAL to test.tpv12 (byte-for-byte, except par.tpv/
# par.C_elastic/the plastic-material block).
#
# Owner clarification PR (board, 2026-10-07): the off-fault stress tensor
# needed here does NOT fit the shared setPlasticStress formula that
# test.tpv27/tpv30/drv.a6 use (that formula pins sxx+syy = 2*szz with a
# single deviatoric-rotation parameter -- a near-vertical strike-slip
# convention). TPV13's spec (Part 2, p.4-5) instead gives, for THIS fault's
# strike-along-x/dip-in-y-z geometry: sigma1 (vertical) = szz, sigma3
# (horizontal, normal to the fault trace) = syy, sigma2 = sxx =
# (sigma1+sigma3)/2 -- principal axes already aligned with x/y/z, no
# rotation, no shear term. meshgen.f90:setPlasticStress gets a new `TPV==13`
# branch (and its mirror in eqdyna3d.py's `build_solver_state`) rather than a
# new subroutine; every other case's formula (the `else` branch) is
# untouched.
#
# Geometry and nucleation reuse test.tpv10's established insertFaultType=1
# (planar dipping fault) + mod4dip machinery verbatim -- TPV10 is itself a
# 60-degree dipping normal fault at nearly the same material properties and
# fault footprint (par.vp,vs,rou below are IDENTICAL to test.tpv10's), and
# that mechanism is already gated fortran + python-jax, so this case adds
# NO new src/ code beyond the setPlasticStress branch above (rule 17 step 3).
#
# Fault: 30 km along-strike x 15 km down-dip (depth 12990.38 m), 60-degree
# dip, reaching the surface. Nucleation zone: 3 km x 3 km square centered
# 12 km down-dip (depth 10392.30 m), along-strike center.
#
# Spec parameters encoded here:
#   material       rho 2700, Vs 3300, Vp 5716 (Part 2, p.4)
#   friction (SW)  outside nucleation: mu_s 0.70; inside nucleation: mu_s
#                  0.54 (lower static friction IS the nucleation mechanism,
#                  per spec Part 1/2 -- "we specify a reduced static
#                  coefficient of friction within the nucleation zone");
#                  mu_d 0.10, d0 0.50 m, frictional cohesion c0 0.2 MPa
#                  everywhere (Part 2 "Friction Parameters and Nucleation").
#                  A node whose sub-fault cell straddles a nucleation-zone
#                  edge or corner gets the spec's own weighted-average mu_s
#                  (0.62 edge, 0.66 corner) rather than a hard inside/outside
#                  cutoff -- see MU_S_EDGE/MU_S_CORNER below.
#   initial stress depth-dependent principal-stress formulas resolved onto
#                  the planar 60-degree dipping fault (Part 2 "Initial Normal
#                  and Shear Stress on the Fault", p.5-6):
#                    (sigma_n - Pf) = 7390.01 Pa/m * down-dip distance,
#                                     tau = 0.549847 * (sigma_n - Pf)
#                    for down-dip distance < 13800 m; at/beyond 13800 m,
#                    (sigma_n - Pf) = 14427.98 Pa/m * down-dip distance,
#                    tau = 0 (stresses become isotropic, spec p.5). A node
#                    whose sub-fault cell straddles the 13800 m boundary gets
#                    the spec's required weighted average (area-overlap
#                    fraction, each side evaluated at its own portion's
#                    midpoint) rather than a hard cutoff -- see the
#                    DOWNDIP_SPLIT_M block below.
#   nucleation     C_nuclea = 0 (no artificial forced-rupture patch): the
#                  reduced static friction coefficient above, applied to the
#                  SAME depth-dependent shear/normal stress formula used
#                  everywhere else, is already what makes the initial shear
#                  exceed the (lower) yield stress inside the patch -- exactly
#                  the spec's own "What's New" mechanism, and exactly
#                  test.tpv10's established pattern (also C_nuclea=0).
#   strength       nodes at the along-strike edges (|x| = fxmax) and the
#   barrier        down-dip edge (z = fzmin, the deepest row) are given
#                  mu_s = 1000 (unbreakable) -- same border convention
#                  test.tpv10 already uses for this exact fault footprint.
#                  The free surface (z = 0 row) is NOT a border (spec: the
#                  fault reaches the surface and surface nodes may rupture).
#   run time       0 - 8 s after nucleation (spec: "you only need to run the
#                  model for 8 seconds" -- a powerful supershear rupture)
#   resolution     spec recommends 100 m node spacing; par.dx below (500 m)
#                  is the fast-gate value (rule 17 step 4)

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.dip = 60  # positive dipping angle tiles the fault to y+ (test.tpv10
              # convention); NOTE par.fzmin, par.ysource, par.zsource,
              # par.dz, par.dt are overwritten by mod4dip below.

par.xmin, par.xmax = -30.0e3, 30.0e3
par.ymin, par.ymax = -20.0e3, 40.0e3
par.zmin, par.zmax = -35.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fzmin, par.fzmax = -15.0e3, 0.0e3

par.xsource = 0.0
par.zsource = -12.0e3  # down-dip distance to the nucleation-patch center
ZSOURCE_DOWNDIP = abs(par.zsource)  # captured before mod4dip overwrites
                                    # par.zsource into a vertical coordinate

par.dx = 500.   # fast gate; spec recommends 100 m
par.dy = par.dx
par.dz = par.dx  # overwritten for the dipping fault by mod4dip below

par.nmat = 1
par.vp, par.vs, par.rou = 5716., 3300., 2700.

par.term = 8.   # spec: 0-8 s after nucleation

# --- TPV13 off-fault plasticity (spec Part 5, p.12) -----------------------
# Density/Vs/Vp are unchanged from TPV12 (par.vp,vs,rou above); cohesion and
# bulk friction below are the two parameters the spec adds.
par.C_elastic = 0
par.coheplas  = 5.0e6     # plastic cohesion c, Pa (spec Part 5)
par.bulk      = 0.85      # bulk friction nu (spec Part 5)
par.gamar     = 0.0       # hydrostatic fluid pressure (spec Part 2: "Fluid
                          # pressure Pf is hydrostatic with water table at
                          # the earth's surface"); meshgen.f90:setPlasticStress
                          # folds Pf into roumax - rhow*(gamar+1), gamar=0 ->
                          # Pf = rhow*grav*depth = 1000*9.8*depth exactly.
par.roumax    = par.rou   # keep bGlobal's roumax tied to the material above
par.output_plastic = 1    # requires C_elastic == 0 (checkInputConsistency.f90)

# TPV13's own return-map algorithm (spec Part 5, p.13-14) is an
# INSTANTANEOUS Drucker-Prager return map (no relaxation term at all) --
# unlike TPV30's Duvaut-Lions viscoplastic law (spec Part 7, Tv=0.05 s).
# EQdyna's shared calcElemKU.f90 Drucker-Prager block only implements the
# Duvaut-Lions form, rjust = exp(-dt/Tv) + (1-exp(-dt/Tv))*Y/sqrt(J2); in the
# limit Tv -> 0, exp(-dt/Tv) -> 0 and rjust -> Y/sqrt(J2) exactly, which IS
# TPV13's own instantaneous return map (spec p.14, Step 4: r =
# Y(trial)/sqrt(J2(trial))). At this case's gate dt (dt = 0.5*dx/vp, computed
# below), dt/Tv is far past IEEE double's exp() underflow-to-exactly-0.0
# threshold (~745), so this is not an approximation of the spec's algorithm
# at double precision -- it reproduces it bit-for-bit. 1e-5 s is chosen
# (not case.setup's pre-v5.9.0 mesh-derived 2*dz/3464 stand-in, which would
# not reach the same instantaneous limit at every dx) for exactly that
# reason, not tuned to any particular dx.
par.viscoplasticRelaxTime = 1.0e-5

# plasticOutputHalfWidth (library_output.f90's pstr.txt* output window): the
# 5x2x8 km default (sized for test.drv.a6) would clip this case's 30x15 km
# fault. Half-widths sized to cover the whole fault footprint (|x|<=15km)
# plus margin, the hanging-wall side the fault dips toward (domain
# ymin,ymax=-20e3,40e3), and down to (and past) the fault's 15 km down-dip
# extent (depth 12990.38 m), without reaching the PML-adjacent domain edges
# (|x|<30km, -20km<y<40km, |z|<35km).
par.plasticOutputHalfWidth = (20.0e3, 25.0e3, 20.0e3)

par.insertFaultType = 1  # insert planar dipping fault (test.tpv10's path)
par.friclaw = 1
par.tpv = 13
par.C_nuclea = 0  # no artificial nucleation; reduced static friction in the
                  # patch (below) is the whole nucleation mechanism
par.nucR = 1.5e3  # nucleation-patch half-width, m (3 km x 3 km square)

par.dz, par.fzmin, par.ysource, par.zsource, par.dt = mod4dip(par.dip,
    par.dx, par.fzmin, par.zsource, par.vp)
print(' ')
print('DIP: dz is    ', par.dz)
print('DIP: fzmin is ', par.fzmin)
print('DIP: ysource, zsource are', par.ysource, par.zsource)
print('DIP: dt is    ', par.dt)

par.nx, par.ny, par.nz = 2, 2, 1   # same partition as test.tpv10
par.HPC_ncpu  = par.nx*par.ny*par.nz
par.HPC_nnode = round(floor(par.HPC_ncpu/128)) + 1
par.HPC_queue = "normal"
par.HPC_time  = "00:10:00"
par.HPC_account = "EAR22013"
par.HPC_email = ""

# Fault grid
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)

par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))

par.fric_sw_fs_nuclea = 0.54   # reduced static friction, nucleation patch
par.fric_sw_fs        = 0.70   # elsewhere inside the rupture rectangle
par.fric_sw_fd        = 0.10
par.fric_sw_D0        = 0.50
par.fric_cohesion     = 0.2e6
par.grav              = 9.8    # exactly 9.8 m/s^2 per spec (Method 2,
                                # C_elastic=0: gravity body force is ON)

# down-dip-distance-dependent effective normal stress / shear-stress
# coefficients (spec Part 2, p.5-6)
SIGN_SHALLOW_PA_PER_M = 7390.01
SIGN_DEEP_PA_PER_M    = 14427.98
TAU_RATIO_SHALLOW     = 0.549847
DOWNDIP_SPLIT_M       = 13800.0

# Grid-alignment tolerance for the exact-equality classifications below. Both
# the nucleation-edge weighting (item 3) and the 13800 m stress boundary
# (item 4) rely on exact node-grid arithmetic (par.nucR and par.dx are exact
# multiples; mod4dip scales dz so the down-dip-distance node spacing is
# exactly par.dx -- see defaultParameters.mod4dip), so this is a float-noise
# guard, not a physical tolerance.
EPS_GRID = 1.0e-6 * par.dx

# Nucleation-zone edge/corner friction weighting (TPV12_13_Description_v6.pdf
# Part 2, "Friction Parameters and Nucleation", p.257-263 in the combined
# PDF): "A fault node represents a sub-fault with finite extent. If a
# sub-fault lies partly inside and partly outside the nucleation zone, the
# friction parameters should be a weighted average." The spec's own worked
# examples: a corner node gets mu_s = 0.75*0.70 + 0.25*0.54 = 0.66, an edge
# node (not a corner) gets mu_s = 0.50*0.70 + 0.50*0.54 = 0.62. Both are
# exact because par.nucR (1500 m) is an exact multiple of par.dx (500 m here,
# and of the spec's 100 m), so a node's own sub-fault cell (half-width
# par.dx/2 in both the along-strike and down-dip-distance directions) either
# lies fully inside, fully outside, straddles exactly one nucleation-zone
# edge (50/50 overlap), or straddles exactly one corner (75/25 overlap) --
# never a partial fraction other than those three.
MU_S_EDGE   = 0.5 *par.fric_sw_fs + 0.5 *par.fric_sw_fs_nuclea   # == 0.62
MU_S_CORNER = 0.75*par.fric_sw_fs + 0.25*par.fric_sw_fs_nuclea   # == 0.66

# 13800 m stress-boundary weighted average (same PDF, "Initial Normal and
# Shear Stress on the Fault", p.177-178/220-221): "A fault node represents a
# sub-fault with finite extent. If a sub-fault extends both above and below
# 13800 m down-dip, the initial normal and shear stresses should be a
# weighted average." The spec gives no explicit formula for this one (unlike
# the friction case above); the down-dip-distance node spacing is exactly
# par.dx (mod4dip), so each node's own down-dip-distance cell is
# [downDipDistance - par.dx/2, downDipDistance + par.dx/2]. For a node whose
# cell straddles 13800 m, this applies the same area-overlap-fraction
# principle as the friction weighting above: each side's stress is evaluated
# at the MIDPOINT of that side's own portion of the cell (exact for a linear
# formula), then the two values are combined weighted by each portion's
# fractional length of the cell.
HALF_DD = par.dx / 2.0

for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    downDipDistance = abs(zcoor)/sin(abs(par.dip)/180.*pi)

    # --- item 3: nucleation-zone friction, with edge/corner weighting -----
    dxn = abs(xcoor - par.xsource)
    ddn = abs(downDipDistance - ZSOURCE_DOWNDIP)
    x_inside  = dxn < par.nucR - EPS_GRID
    z_inside  = ddn < par.nucR - EPS_GRID
    x_on_edge = abs(dxn - par.nucR) < EPS_GRID
    z_on_edge = abs(ddn - par.nucR) < EPS_GRID

    if x_inside and z_inside:
        mu_s = par.fric_sw_fs_nuclea            # fully inside: 0.54
    elif x_on_edge and z_on_edge:
        mu_s = MU_S_CORNER                      # corner: 0.66
    elif (x_on_edge and z_inside) or (z_on_edge and x_inside):
        mu_s = MU_S_EDGE                        # edge (not corner): 0.62
    else:
        mu_s = par.fric_sw_fs                   # fully outside: 0.70
    par.on_fault_vars[iz,ix,1] = mu_s

    # strength barrier: along-strike edges and the deepest (down-dip) row;
    # the free surface (iz at zcoor==0) is NOT a border.
    if abs(abs(xcoor) - par.fxmax) < 0.01 or abs(zcoor - par.fzmin) < 0.01:
        par.on_fault_vars[iz,ix,1] = 1000.

    par.on_fault_vars[iz,ix,2] = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3] = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4] = par.fric_cohesion

    # --- item 4: 13800 m stress boundary, weighted average on straddle ----
    d_lo = downDipDistance - HALF_DD
    d_hi = downDipDistance + HALF_DD
    if d_hi <= DOWNDIP_SPLIT_M + EPS_GRID:
        sigEff = -SIGN_SHALLOW_PA_PER_M*downDipDistance
        tau    = -abs(TAU_RATIO_SHALLOW*sigEff)
    elif d_lo >= DOWNDIP_SPLIT_M - EPS_GRID:
        sigEff = -SIGN_DEEP_PA_PER_M*downDipDistance
        tau    = 0.0
    else:
        frac_shallow   = (DOWNDIP_SPLIT_M - d_lo) / par.dx
        frac_deep      = 1.0 - frac_shallow
        dd_shallow_mid = 0.5*(d_lo + DOWNDIP_SPLIT_M)
        dd_deep_mid    = 0.5*(DOWNDIP_SPLIT_M + d_hi)
        sigEff_shallow = -SIGN_SHALLOW_PA_PER_M*dd_shallow_mid
        tau_shallow    = -abs(TAU_RATIO_SHALLOW*sigEff_shallow)
        sigEff_deep    = -SIGN_DEEP_PA_PER_M*dd_deep_mid
        tau_deep       = 0.0
        sigEff = frac_shallow*sigEff_shallow + frac_deep*sigEff_deep
        tau    = frac_shallow*tau_shallow    + frac_deep*tau_deep
    par.on_fault_vars[iz,ix,7]  = sigEff  # initial effective normal stress
    par.on_fault_vars[iz,ix,8]  = 0.0     # no along-strike (horizontal) shear
    # positive initial dip stress moves the y+ (hanging) wall up; negative
    # moves it down -- negative drives normal-fault motion (test.tpv10
    # convention, same fault sense).
    par.on_fault_vars[iz,ix,49] = tau

    par.on_fault_vars[iz,ix,46] = 0.0     # no creep/background slip rate

# --- Stations (spec Parts 7-8) ----------------------------------------------
# on-fault: (x, z) in km, z = actual vertical coordinate (negative = depth),
# NOT down-dip distance -- same convention test.tpv10 already uses.
C = sin(par.dip/180.*pi)
par.st_coor_on_fault = [
    [0.0, 0.0], [4.5, 0.0], [12.0, 0.0],
    [0.0, -1.5*C], [0.0, -3.0*C], [0.0, -4.5*C], [0.0, -7.5*C],
    [4.5, -7.5*C], [12.0, -7.5*C], [0.0, -12.0*C]]
# off-fault: (x, y, z) in km
par.st_coor_off_fault = [
    [0, 1, 0], [0, -1, 0], [0, 2, 0], [0, -2, 0], [0, 3, 0], [0, -3, 0],
    [0, 0.5 + 0.3/tan(abs(par.dip)/180.*pi), -0.3],
    [0, -0.5 + 0.3/tan(abs(par.dip)/180.*pi), -0.3],
    [0, 1.0 + 0.3/tan(abs(par.dip)/180.*pi), -0.3],
    [0, -1.0 + 0.3/tan(abs(par.dip)/180.*pi), -0.3],
    [12, 3, 0], [12, -3, 0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
