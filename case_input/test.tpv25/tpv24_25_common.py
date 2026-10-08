#! /usr/bin/env python3
"""
Shared build logic for test.tpv24 and test.tpv25 (SCEC TPV24/25 "Branched
Fault" benchmarks, TPV24_25_Description_v07.pdf, fetched into
scratch/specs/). Row 153 checkpoint 2b.

STATUS: WIP skeleton, NOT yet gated. Geometry/stress/friction derivation
below is complete and self-checked (see derivation notes), but has NOT been
run through the solver yet. Do not freeze a reference or register in
testNameList.py/matrix.py off this file alone.

Geometry (Part 1/2, p.3-5): two vertical faults, both extruded along z
(reaching the surface, 15 km deep). Fault 1 (main): 28 km long, straight,
strike = code x-axis, at code y=0 -- same style-0 convention as every other
planar fault in this repo. Fault 2 (branch): 12 km long, meets fault 1 at a
junction point 12 km from fault 1's right (+x) end, trending 30 deg away
from fault 1's strike. The junction node is owned by fault 1 only (spec
p.5: "the branch fault ends at the junction point... slip goes to zero at
the junction, like any other fault border"); fault 2's own box EXCLUDES the
junction column (fxmin = +dx, matching testsys/parity/probe_branch_mesh.f90's
validated convention, which also established the non-isotropic-mesh
requirement below).

In code axes: junction at x=0. Fault 1 spans x in [-16000, 12000] (28 km,
12 km of it beyond the junction toward +x, matching "junction 12 km from
the main fault's end"). Fault 2 (style 2, faultDegenCode = 130 = 100+30)
follows y = -x*tan(30deg) for x>0 (probe_branch_mesh.f90's verified line
equation/winding, NOT the naive y=+x*tan(30) guess -- the sign is fixed by
library_degeneration.f90's wedge()/reorder() node permutation, unmodifiable,
see the probe's "v2: extrude-direction reversed" comment). Hypocenter is on
fault 1, 8 km from the junction (away from the branch, so rupture must
travel toward the junction before it can jump), 10 km deep: code (x,y,z) =
(-8000, 0, -10000).

MESH: style-2 wedge meshing needs dy = dx*tan(30deg) (NON-isotropic -- this
is NOT the style-1 dip convention of test.tpv36/37, which uses
dy=dx*cos(dip)/dz=dx*sin(dip); verified directly in
testsys/parity/probe_branch_mesh.f90, PR #143). dz = dx (planar fault,
reaches the surface, same as every style-0 fault in this repo).

BRANCH-LENGTH APPROXIMATION (document per coordinator instruction): the
branch's true along-strike length is 12000 m, so its true x-extent is
12000*cos(30deg) = 10392.3 m. At gate-coarse dx=1000 m that is not an
integer multiple of dx; we decimate (never interpolate, rule 17 step 2) to
the nearest achievable multiple, fxmax2 = 10*dx = 10000 m (branch x-extent
10.0 km, true along-strike length 10.0/cos(30) = 11.547 km -- a 3.8% short
branch, not the spec's 12 km). This also applies, smaller in relative terms,
at the 100 m/50 m submission resolutions (full_specs.py tier entries below
record this caveat without running those resolutions).

MATERIAL (Part "Material Properties", same table as TPV5): rho=2670 kg/m^3,
vs=3464 m/s, vp=6000 m/s.

STRESS (Part "Initial Stress Tensor", p.6-7): sigma11 = compressive vertical
stress = -(rho*g)*depth, rho*g = 2670*9.8 = 26166 Pa/m (g EXACTLY 9.8, spec's
own instruction, not 9.81/9.80665). Fluid pressure Pf = (1000*9.8)*depth =
9800 Pa/m, hydrostatic, water table at surface. SP == sigma11_effective =
sigma11+Pf = -(26166-9800)*depth = -16366.0*depth (negative=compressive,
matches code convention).
sigma22_eff (fault-parallel, == main-fault-strike direction) = b22*SP
sigma33_eff (fault-perpendicular, == main-fault-normal direction) = b33*SP
   -- "equals the negative of the total normal stress on the main fault"
      (spec p.6): this IS on_fault_vars[7] for fault 1 directly.
sigma23 (shear in the horizontal plane, + = right-lateral on MAIN fault)
   = b23*SP -- this IS on_fault_vars[8] for fault 1 directly.
b22/b33/b23 are the ONLY difference between TPV24 and TPV25 (spec p.6,
explicit), passed into build_params() below; all other formulas identical.

BRANCH RESOLUTION (derivation, not in the spec text verbatim -- the spec
gives the tensor in the MAIN FAULT's own (strike=x, normal=y) frame, which
IS the global (x,y) frame since fault 1's strike defines code-x; resolving
onto fault 2 is a standard in-plane stress-tensor rotation, independently
re-derived and sign-checked against the probe's y=-x*tan(phi) branch-line
equation (not the φ=+) before use):
  global sigma_xx = sigma22_eff = b22*SP, sigma_yy = sigma33_eff = b33*SP,
  sigma_xy = sigma23 = b23*SP (all at the SAME local depth -- fault 2 is
  vertical too, so "depth" is still just -z, independent of y).
  Fault 2's tangent is at angle theta=-phi from +x (probe's line is
  y=-x*tan(phi) for x>0, i.e. tangent direction (cos(phi),-sin(phi)) =
  (cos(theta),sin(theta)) with theta=-phi). Standard rotation of a 2D stress
  tensor onto a plane with tangent angle theta, normal n=(-sin(theta),
  cos(theta)):
    sigma_nn = sigma_xx*sin(theta)^2 + sigma_yy*cos(theta)^2
               - 2*sigma_xy*sin(theta)*cos(theta)
    sigma_nt = (sigma_yy-sigma_xx)*sin(theta)*cos(theta)
               + sigma_xy*cos(2*theta)
  With theta=-phi this reduces (sin(-phi)=-sin(phi), cos(-phi)=cos(phi),
  cos(-2phi)=cos(2phi)) to:
    sigma_nn = SP*[ b22*sin(phi)^2 + b33*cos(phi)^2 + 2*b23*sin(phi)*cos(phi) ]
    sigma_nt = SP*[ (b22-b33)*sin(phi)*cos(phi) + b23*cos(2*phi) ]
  Sanity check at phi=0 (theta=0): sigma_nn=SP*b33, sigma_nt=SP*b23 --
  reduces exactly to the main-fault formula above. Numerically (phi=30deg,
  sin*cos=0.4330127, cos(2phi)=0.5):
    TPV24 (b22=0.926793,b33=1.073206,b23=-0.169029):
       normal_coeff=0.8902098  shear_coeff=-0.1479345
    TPV25 (b22=1.119338,b33=0.880661,b23=0.138704):
       normal_coeff=1.0604563  shear_coeff=+0.1726990
  TPV24's branch normal_coeff < TPV25's (less compressively clamped at the
  same depth/SP) -- consistent with the spec's own "TPV24 releasing /
  TPV25 restraining" branch-sense narrative (narrative text, NOT itself a
  derivation input; this agreement is a sign-convention CHECK, not a
  substitute for the formula above).

FRICTION (Part "Friction Parameters and Nucleation", p.7-8): mu_s=0.18,
mu_d=0.12, d0=0.30 m. Cohesion C0: 3.00 MPa at the surface, 0.30 MPa at
depth>=4000 m, linear taper in between -- same functional form as
TPV22/23's cohesion (0.30 MPa floor, not 0.0): C0 = 300000 + 675*(4000-depth)
Pa for depth<=4000, else 300000 Pa (675 Pa/m = (3.00-0.30)e6/4000).
Nucleation: r_crit=4000 m, t0=0.50 s -- same T(r) smoothed forced-rupture
formula as TPV22/23/26/27/29/30/36/37/201 (faulting.f90/faulting.py
swtwNucleation, TPV in (...,24,25) branch added this row, commit fdd9ccc,
TPV24_25_Description_v07.pdf Part 5 p.13: same r_crit symbol/value, same
0.081 taper coefficient, same 0.7*Vs rupture speed, same Vs=3464 (matches
NUC_VS_FIXED), same t0=0.5s).

Resolution: spec requests submission at 100 m AND 50 m node spacing (p.8).
This gate build uses a COARSE dx=1000 m (minutes at 4 ranks, rule 17 step
4); the 100 m/50 m tiers are recorded in full_specs.py WITHOUT being run
(rule 17 step 5).

n-stress sign convention: par.faultStNormalStressSign = 'extension' (spec's
own sigma33 "negative=compression" convention matches defaultParameters.py's
default -- no override needed, same situation as TPV22/23).
"""
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

FAULT1_XRANGE = (-16000.0, 12000.0)   # main fault, 28 km, junction at x=0
FZMIN, FZMAX = -15000.0, 0.0          # both faults: 15 km deep, reach surface
BRANCH_ANGLE_DEG = 30.0
BRANCH_FXMIN = 1000.0   # excludes junction column (owned by fault 1), see docstring
BRANCH_FXMAX = 10000.0  # decimated from true 10392.3 m, see BRANCH-LENGTH note

MU_S, MU_D, D0 = 0.18, 0.12, 0.30
NUCR = 4000.0
TW_T0 = 0.50
RHOG = 2670.0 * 9.8        # 26166 Pa/m, g EXACTLY 9.8 (spec's own instruction)
PFG = 1000.0 * 9.8         # 9800 Pa/m
UNBREAKABLE_FS = 1000.0

_phi = BRANCH_ANGLE_DEG / 180.0 * pi
_sinphi, _cosphi = sin(_phi), cos(_phi)
BRANCH_TANPHI = tan(_phi)


def SP(depth_m):
    """sigma11_effective = sigma11 + Pf, depth_m positive-down."""
    return -(RHOG - PFG) * depth_m


def cohesion(depth_m):
    if depth_m <= 4000.0:
        return 300000.0 + 675.0 * (4000.0 - depth_m)
    return 300000.0


def build_on_fault_vars_main(fx, fz, nfx, nfz, b22, b33, b23):
    v = np.zeros((nfz, nfx, 100))
    for ix, xcoor in enumerate(fx):
        for iz, zcoor in enumerate(fz):
            depth = -zcoor
            sp = SP(depth)
            v[iz, ix, 1] = MU_S
            v[iz, ix, 2] = MU_D
            v[iz, ix, 3] = D0
            v[iz, ix, 4] = cohesion(depth)
            v[iz, ix, 5] = TW_T0
            v[iz, ix, 7] = b33 * sp   # fault-normal, already negative=compressive
            v[iz, ix, 8] = b23 * sp   # right-lateral shear on main fault
            on_left_edge = (ix == 0)
            on_right_edge = (ix == nfx - 1)
            on_bottom_edge = (iz == 0)
            if on_left_edge or on_right_edge or on_bottom_edge:
                v[iz, ix, 1] = UNBREAKABLE_FS
    return v


def build_on_fault_vars_branch(fx2, fz, nfx2, nfz, b22, b33, b23):
    normal_coeff = b22 * _sinphi**2 + b33 * _cosphi**2 + 2.0 * b23 * _sinphi * _cosphi
    shear_coeff = (b22 - b33) * _sinphi * _cosphi + b23 * cos(2.0 * _phi)
    v = np.zeros((nfz, nfx2, 100))
    for ix, xcoor in enumerate(fx2):
        for iz, zcoor in enumerate(fz):
            depth = -zcoor
            sp = SP(depth)
            v[iz, ix, 1] = MU_S
            v[iz, ix, 2] = MU_D
            v[iz, ix, 3] = D0
            v[iz, ix, 4] = cohesion(depth)
            v[iz, ix, 5] = TW_T0
            v[iz, ix, 7] = normal_coeff * sp
            v[iz, ix, 8] = shear_coeff * sp
            # true borders unbreakable (rule: PDF Part 3-style border
            # convention, matching every other TPV in this repo):
            # left edge is the junction-adjacent column (excluded from the
            # box already, so ix==0 here IS the first real branch node, one
            # dx away from the junction -- still a true border of fault 2,
            # per spec p.5 "slip goes to zero at the junction... like any
            # other border"), right edge, and bottom.
            on_left_edge = (ix == 0)
            on_right_edge = (ix == nfx2 - 1)
            on_bottom_edge = (iz == 0)
            if on_left_edge or on_right_edge or on_bottom_edge:
                v[iz, ix, 1] = UNBREAKABLE_FS
    return v


def build_params(tpv, b22, b33, b23):
    par = parameters()
    par.ntotft = 2

    par.xmin, par.xmax = -24.0e3, 20.0e3
    par.ymin, par.ymax = -10.0e3, 10.0e3
    par.zmin, par.zmax = -22.0e3, 0.0e3

    par.fxmin, par.fxmax = FAULT1_XRANGE
    par.fymin, par.fymax = 0.0, 0.0
    par.fzmin, par.fzmax = FZMIN, FZMAX

    par.faultgeom = [
        (FAULT1_XRANGE[0], FAULT1_XRANGE[1], 0.0, 0.0, FZMIN, FZMAX),
        (BRANCH_FXMIN, BRANCH_FXMAX, 0.0, 0.0, FZMIN, FZMAX),
    ]
    par.faultDegenCode = [0, int(100 + BRANCH_ANGLE_DEG)]  # style 0, style 2 @ 30 deg

    par.xsource, par.ysource, par.zsource = -8000.0, 0.0, -10000.0

    par.dx = 1000.0
    par.dz = par.dx
    par.dy = par.dx * BRANCH_TANPHI   # non-isotropic, see docstring

    par.nmat = 1
    par.vp, par.vs, par.rou = 6000.0, 3464.0, 2670.0

    par.term = 5.0       # GATE_TERM_S; spec's own 12.0s post-nucleation NOT
                          # used here -- see case docstring / PR body for the
                          # measured fraction of fault 2 ruptured at 5s.
    par.dt = 0.5 * par.dx / par.vp

    par.C_elastic = 1
    par.C_nuclea = 1
    par.insertFaultType = 0
    par.friclaw = 1
    par.tpv = tpv
    par.nucR = NUCR
    par.nucfault = 1

    par.nx, par.ny, par.nz = 2, 1, 2
    par.HPC_ncpu = par.nx * par.ny * par.nz
    par.HPC_nnode = round(floor(par.HPC_ncpu / 128)) + 1
    par.HPC_queue = "normal"
    par.HPC_time = "01:00:00"
    par.HPC_account = "EAR22013"
    par.HPC_email = ""

    par.nfx = round((FAULT1_XRANGE[1] - FAULT1_XRANGE[0]) / par.dx + 1)
    par.nfz = round((FZMAX - FZMIN) / par.dz + 1)
    par.fx = np.linspace(FAULT1_XRANGE[0], FAULT1_XRANGE[1], par.nfx)
    par.fz = np.linspace(FZMIN, FZMAX, par.nfz)
    par.on_fault_vars = build_on_fault_vars_main(par.fx, par.fz, par.nfx, par.nfz, b22, b33, b23)

    nfx2 = round((BRANCH_FXMAX - BRANCH_FXMIN) / par.dx + 1)
    fx2 = np.linspace(BRANCH_FXMIN, BRANCH_FXMAX, nfx2)
    fault2_vars = build_on_fault_vars_branch(fx2, par.fz, nfx2, par.nfz, b22, b33, b23)
    par.onFaultVarsPerFault = [par.on_fault_vars, fault2_vars]

    # On-fault stations: SCEC spec Part 7 (TPV24_25_Description_v07.pdf,
    # p.17-19), 8 stations on the main fault + 6 on the branch fault,
    # transcribed verbatim as (x_km, z_km) pairs. Along-strike x is measured
    # relative to the junction point (spec's own convention, p.17 Note),
    # which is exactly code x=0 here -- so the main-fault stations need no
    # conversion, and every one lands exactly on the dx=1000 m gate grid
    # (x in {-8,-2,9} km, z in {0,-5,-10} km, all multiples of dx/dz=1 km):
    #   faultst-080dp000, faultst-020dp000, faultst090dp000 (0 km down-dip)
    #   faultst-080dp050 (5 km down-dip)
    #   faultst-080dp100 (hypocenter), faultst-020dp100, faultst020dp100,
    #   faultst090dp100 (10 km down-dip)
    fault1_stations = [
        [-8.0, 0.0], [-2.0, 0.0], [9.0, 0.0],
        [-8.0, -5.0],
        [-8.0, -10.0], [-2.0, -10.0], [2.0, -10.0], [9.0, -10.0],
    ]
    # Branch-fault stations (spec p.19): location is given as distance ALONG
    # THE STRIKE OF THE BRANCH FAULT from the junction (L), but the gate
    # mesh's fault-2 box (probe_branch_mesh.f90-verified; see module
    # docstring) is parameterized by CODE-x, not along-strike arc length --
    # x_code = L*cos(30deg). setOnFaultStation (meshgen.f90,
    # report_dropped_onfault_st) matches a station's x AND z against a fault
    # node EXACTLY within tol, so a raw L*cos(30) value (e.g. 2.0 km ->
    # 1.732 km) would match no node and get silently dropped with no
    # faultst* file. Snapped to the nearest dx=1000 m fault-2 node instead
    # (spec's own explicit allowance, p.17/p.23: "you can move the station
    # to the nearest node"):
    #   L=1.0 km -> 0.866 km -> nearest node 1.0 km (== BRANCH_FXMIN, the
    #               near-junction node)
    #   L=2.0 km -> 1.732 km -> nearest node 2.0 km
    #   L=9.0 km -> 7.794 km -> nearest node 8.0 km
    # giving (snapped x_km, z_km):
    #   branchst020dp000 -> faultstft2_020dp000 (2.0, 0 km down-dip)
    #   branchst090dp000 -> faultstft2_080dp000 (2.0->8.0 snap, 0 km)
    #   branchst020dp050 -> faultstft2_020dp050 (2.0, 5 km down-dip)
    #   branchst010dp100 -> faultstft2_010dp100 (1.0, 10 km down-dip)
    #   branchst020dp100 -> faultstft2_020dp100 (2.0, 10 km down-dip)
    #   branchst090dp100 -> faultstft2_080dp100 (2.0->8.0 snap, 10 km)
    fault2_stations = [
        [2.0, 0.0], [8.0, 0.0],
        [2.0, -5.0],
        [1.0, -10.0], [2.0, -10.0], [8.0, -10.0],
    ]
    par.st_coor_on_fault = fault1_stations
    par.st_coor_on_fault_per_fault = [fault1_stations, fault2_stations]

    # Off-fault stations: SCEC spec Part 8 (p.23), 8 stations, all at the
    # earth's surface (0 km depth). Station name's first number is the
    # horizontal perpendicular offset from the MAIN fault (code y, since the
    # main fault is the code y=0 plane); positive = far side (spec's own
    # sign, carried directly into code y, no repo convention to reconcile --
    # this case has no prior off-fault sign history to match). [x_km, y_km,
    # z_km]:
    #   body030st-020dp000  / body-030st-020dp000  (x=-2.0, y=+-3.0)
    #   body030st020dp000   / body-006st020dp000   / body-042st020dp000 (x=2.0)
    #   body030st080dp000   / body-023st080dp000   / body-076st080dp000 (x=8.0... )
    # NOTE: spec's x for the x=8 group is written "st080" (8.0 km along
    # strike) -- same junction-relative convention as the on-fault table.
    par.st_coor_off_fault = [
        [-2.0, 3.0, 0.0], [-2.0, -3.0, 0.0],
        [2.0, 3.0, 0.0], [2.0, -0.6, 0.0], [2.0, -4.2, 0.0],
        [8.0, 3.0, 0.0], [8.0, -2.3, 0.0], [8.0, -7.6, 0.0],
    ]
    par.n_on_fault = len(par.st_coor_on_fault)
    par.n_off_fault = len(par.st_coor_off_fault)

    return par
