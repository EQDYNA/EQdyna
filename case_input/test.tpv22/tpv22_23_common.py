#! /usr/bin/env python3
"""
Shared build logic for test.tpv22 and test.tpv23 (SCEC TPV22/23 "Stepover
Benchmarks", TPV22_23_Description_v08.pdf, strike.scec.org/cvws/tpv22_23docs.html,
fetched into scratch/specs/tpv2223/).

Geometry (PDF Part 1/2, p.3-4): two vertical, planar, right-lateral
strike-slip faults, each 30 km along-strike x 20 km deep, reaching the
surface, with a 10 km along-strike overlap. In the spec's own (x,y,z) frame
(x: along-strike, y: depth positive down, z: horizontal, perpendicular to
strike -- the stepover direction):
  fault #1: -25000 <= x <= 5000,   0 <= y <= 20000,  z = 0
  fault #2:  -5000 <= x <= 25000,  0 <= y <= 20000,  z = -1600  (TPV22, right
             step, extensional) or z = +1000 (TPV23, left step, compressional)
hypocenter (x,y,z) = (-10000, 10000, 0), centered along-strike and at 10 km
depth in fault #1 (p.5).

EQdyna's own axis convention (this repo, e.g. case_input/test.multifault2) is
x=along-strike, y=horizontal fault-normal/stepover direction, z=depth
(negative down). So this module maps spec (x,y,z) -> code (x, z_spec, -y_spec)
i.e. code_y = spec_z (the stepover offset, unchanged sign and magnitude) and
code_z = -spec_y (depth, negative down). Material/friction/stress formulas
below are written directly in code axes.

GEOMETRY (restored, not invented -- mira/row17-tpv2223-rebased, row 17
widened audit): a prior iteration on this branch found checkInputConsistency.
f90's old row-17 guard (ERR_GEOM_MULTIFAULT_XZ_BAD) required EVERY fault to
share EXACTLY fault 1's x/z box -- the uniform x/z mesh belt used to be built
once, from fault 1's box alone (getLocalOneDimCoorArrAndSize, hardcoded
fltxyz(.,.,1)). TPV22/23's real geometry has fault #2 offset 20 km ALONG
STRIKE from fault #1 (only a 10 km overlap out of each fault's 30 km length),
so the two faults' TRUE x-extents differ -- exactly what that guard refused.
That prior iteration worked around the gap with a UNION-box scaffold (both
faults meshed over the same 50 km box, with the portion outside each fault's
own 30 km forced unbreakable) instead of fixing the mesh generator. The owner
caught this ("Why you invent new vars? ... the support was lost over time")
and the generator itself is now fixed (meshgen.f90's getLocalOneDimCoorArrAndSize
unions every fault's OWN box; checkInputConsistency.f90 no longer requires a
shared x/z extent) -- so each fault below gets its TRUE box, nothing shared,
no scaffold, no scaffold-stress fix, no per-resolution margin bump:
  fault #1: x in [-25000,  5000], z in [-20000, 0]  (code axes)
  fault #2: x in [ -5000, 25000], z in [-20000, 0]
A fault's TRUE borders (its own real left/right edge, and bottom [z=-20000])
get the unbreakable treatment, matching PDF Part 3, p.5: "Slip goes to zero
at the border of a fault ... a node which lies precisely on the border of a
fault should not be permitted to slip. This is a CHANGE from earlier
benchmarks" -- kept from the 2026-10-01 iteration (the one real fix from that
work, see build_on_fault_vars' inline comment). The TOP edge (z=0, the
surface trace) is deliberately NOT included: both independent references
(kaneko/SPECFEM3D, payne/EQdyna) show large early slip at the station
literally named "0 km down-dip" (onset ~5.35s, peak 2.5-2.8 m, both tpv22 and
tpv23) -- incompatible with that node being unbreakable.

Material (PDF p.5): rho=2670 kg/m^3, Vs=3464 m/s, Vp=6000 m/s (same as TPV5).
Stress (PDF p.5): sigma_ini = 60.00 MPa (normal, compressive, constant with
depth -- NOT depth-dependent, unlike TPV29). tau_ini = 29.38 MPa for depth <=
15000 m, 29.38 - 0.002938*(depth-15000) MPa for 15000 <= depth <= 20000 m
(the "soft bottom" shear-stress taper, p.5) -- this is a one-off linear
formula local to this benchmark, not the repo's devStrTaperDepthStart/End
viscoplastic (Omega) taper mechanism (scripts/lib.py resolveViscoplasticParams)
-- that mechanism exists for TPV29/30-style off-fault Drucker-Prager plastic
taper, unrelated: TPV22/23 is pure linear elastic (p.5), so it is not used here.
Friction (PDF p.6): mu_s=0.548, mu_d=0.373, d0=0.30 m, cohesion C0 =
0.0014 MPa/m * (5000-depth) for depth<=5000 m else 0.0 (7.0 MPa at the
surface, 0 at 5 km -- matches FRIC_SLOT_COHESION, globalvar.f90:25),
r_crit=3000 m, t0=0.50 s (forced-rupture decay time, FRIC_SLOT_TW_T0).
Nucleation (PDF p.6/10): smoothed forced rupture, same T(r) formula and same
numeric constants (Vs=3464, 0.081 taper coefficient, 0.7*Vs rupture speed) as
faulting.f90:swtwNucleation's TPV29/36/37/201 branch -- this IS that formula's
own source benchmark (see faulting.f90/faulting.py comments at the TPV==22/23
branch added for this mission).

n-stress sign convention (PDF Part 6, p.15): "Positive means extension" --
par.faultStNormalStressSign = 'extension', which is also defaultParameters.py's
default (no override needed).
"""
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

# Each fault's TRUE box, in CODE axes (x=strike, y=stepover/fault-normal,
# z=depth<=0). z is identical for both (20 km deep, reaching the surface);
# x differs -- that is the whole point of a stepover, and the restored mesh
# generator (meshgen.f90's getLocalOneDimCoorArrAndSize) now unions every
# fault's own box rather than requiring them to match.
FAULT1_XRANGE = (-25000.0, 5000.0)   # fault #1's TRUE along-strike extent
FAULT2_XRANGE = (-5000.0, 25000.0)   # fault #2's TRUE along-strike extent
FZMIN, FZMAX = -20000.0, 0.0         # both faults: 20 km deep, reaching the surface

MU_S, MU_D, D0 = 0.548, 0.373, 0.30
NUCR = 3000.0           # r_crit, m (PDF p.6)
TW_T0 = 0.50            # forced-rupture decay time, s (PDF p.6)
SIGMA_INI = 60.0e6      # normal stress, Pa, compressive (code convention: negative)
UNBREAKABLE_FS = 1000.0  # same sentinel test.tpv8/test.tpv29 use for domain-edge/
                         # inert nodes: fs this high never reaches failure.


def tau_ini(depth_m):
    """PDF p.5: 29.38 MPa for depth<=15000 m, soft-bottom taper to 14.69 MPa
    at depth=20000 m. depth_m is positive-down."""
    if depth_m <= 15000.0:
        return 29.38e6
    return 29.38e6 - 2938.0 * (depth_m - 15000.0)  # 0.002938 MPa/m = 2938 Pa/m


def cohesion(depth_m):
    """PDF p.6: C0 = 0.0014 MPa/m * (5000-depth) for depth<=5000 m, else 0."""
    if depth_m <= 5000.0:
        return 1400.0 * (5000.0 - depth_m)  # 0.0014 MPa/m = 1400 Pa/m
    return 0.0


def build_on_fault_vars(fx, fz, nfx, nfz):
    """One fault's on_fault_vars array, built on ITS OWN box (fx/fz span that
    fault's true extent -- see module docstring; no mesh-sharing scaffold, no
    "active" node distinction, every node here is real fault physics).
    Nucleation (swtwNucleation, TPV==22/23 branch) only acts on nodes where
    C_nuclea==1 AND this fault is par.nucfault -- no manual stress patch
    needed here (unlike the older "shear bump" style of test.tpv8/
    test.multifault2; TPV22/23 is forced-rupture-only, same as test.tpv29)."""
    v = np.zeros((nfz, nfx, 100))
    for ix, xcoor in enumerate(fx):
        for iz, zcoor in enumerate(fz):
            depth = -zcoor  # zcoor <= 0; depth positive-down
            v[iz, ix, 1] = MU_S
            v[iz, ix, 2] = MU_D
            v[iz, ix, 3] = D0
            v[iz, ix, 4] = cohesion(depth)
            v[iz, ix, 5] = TW_T0
            v[iz, ix, 7] = -SIGMA_INI               # init normal stress, negative=compressive
            v[iz, ix, 8] = tau_ini(depth)            # init strike (right-lateral) shear
            # True borders: this fault's own left/right edge (ix==0/nfx-1),
            # and bottom (iz==0, z=FZMIN). Per PDF Part 3 p.5 ("a node which
            # lies precisely on the border of a fault should not be permitted
            # to slip ... CHANGE from earlier benchmarks"). The TOP edge
            # (z=0, the surface trace) is deliberately EXCLUDED -- see module
            # docstring: both independent references show large early slip at
            # the "0 km down-dip" station, incompatible with that node being
            # unbreakable (measured directly: our own run showed exactly 0
            # slip there at term=15s with the top edge included, 2026-10-01
            # iteration log).
            on_left_edge = (ix == 0)
            on_right_edge = (ix == nfx - 1)
            on_bottom_edge = (iz == 0)
            if on_left_edge or on_right_edge or on_bottom_edge:
                v[iz, ix, 1] = UNBREAKABLE_FS
    return v


def build_params(tpv, fault2_z):
    """fault2_z is the stepover offset in CODE y (== spec z): -1600 (TPV22,
    extensional right step) or +1000 (TPV23, compressional left step), taken
    directly from PDF Parts 1/2 (p.3-4)."""
    par = parameters()
    par.ntotft = 2

    # Domain margins beyond the union footprint of both faults' TRUE x-extents
    # ([-25000,25000], ~7 km in x/z, matching test.tpv8's margin style; y
    # margin generous since the stepover offset is tiny, <=1.6 km, relative
    # to the domain). The mesh generator unions the faults' own boxes
    # internally (meshgen.f90); this domain footprint is unchanged by that --
    # it always had to cover both faults' full physical extent.
    par.xmin, par.xmax = -32.0e3, 32.0e3
    par.ymin, par.ymax = -10.0e3, 10.0e3
    par.zmin, par.zmax = -27.0e3, 0.0e3

    # ntotft==1 fallback fields (other scripts may still read these for
    # fault 1 specifically) -- fault 1's OWN true box, not a shared union.
    par.fxmin, par.fxmax = FAULT1_XRANGE
    par.fymin, par.fymax = 0.0, 0.0
    par.fzmin, par.fzmax = FZMIN, FZMAX

    # Each fault's TRUE box -- restored per-fault mesh extent (meshgen.f90's
    # getLocalOneDimCoorArrAndSize now unions every fault's own box; no
    # shared-box requirement left to work around).
    par.faultgeom = [
        (FAULT1_XRANGE[0], FAULT1_XRANGE[1], 0.0, 0.0, FZMIN, FZMAX),
        (FAULT2_XRANGE[0], FAULT2_XRANGE[1], fault2_z, fault2_z, FZMIN, FZMAX),
    ]

    # Hypocenter: spec (x,y,z)=(-10000,10000,0) -> code (x,y,z)=(-10000,0,-10000)
    par.xsource, par.ysource, par.zsource = -10000.0, 0.0, -10000.0

    # Resolution: ISOTROPIC dx=dy=dz (NOTES_tpv2223_iteration.md's own
    # resolution-scan evidence, iteration 5 -- owner-judged "good matches"):
    # 200 m for TPV22 (the extensional, harder-to-jump stepover: at coarser
    # 1000/400/400m fault #2 never ruptured past a 0.7%-of-nodes patch; 200m
    # isotropic is where it jumps, onset/peak-slip/peak-sr within 5-20% of
    # both independent references at all 3 gated stations), 250 m for TPV23
    # (jumps comfortably even gate-coarse; 250m isotropic sharpens onset/
    # peak-slip further). Both commensurate with the stepover offset
    # (checkInputConsistency's fault-y commensurability guard, now enforced
    # inside meshgen.f90's getLocalOneDimCoorArrAndSize/checkInputConsistency:
    # 1600/200=8, 1000/250=4) and with each fault's own x/z extent relative to
    # the union belt origin (every bound here is a multiple of 5000 m, itself
    # a multiple of both 200 and 250). No more per-resolution nuni_y margin
    # bump needed: the restored y-belt (union of every fault's own y) derives
    # its margin from the fault geometry itself, not a hand-set constant.
    par.dx = 200.0 if tpv == 22 else 250.0
    par.dz = par.dx
    par.dy = par.dx

    par.nmat = 1
    par.vp, par.vs, par.rou = 6000.0, 3464.0, 2670.0

    par.term = 15.0  # spec's own full duration (TPV22_23_Description_v08.pdf
    # p.7, "Run the model for times from 0.0 to 15.0 seconds"); GATE_TERM_S
    # (5.0) is too short even for fault #1's own far stations (onset ~5.37s
    # in both independent references) -- never mind fault #2 (~9-12s).

    par.dt = 0.5 * par.dx / par.vp  # isotropic dx=dy=dz; CFL on the common cell size

    par.C_elastic = 1
    par.C_nuclea = 1
    par.insertFaultType = 0   # planar fault, no rough-fault machinery
    par.friclaw = 1
    par.tpv = tpv             # faulting.f90/faulting.py TPV==22/23 branch
    par.nucR = NUCR
    par.nucfault = 1          # hypocenter is on fault #1 only (PDF p.5)

    par.nx, par.ny, par.nz = 2, 1, 2   # ny=1: keeps every fault y-plane off
                                        # every MPI partition boundary (both
                                        # faults' y get the full range on
                                        # every rank); 4 ranks via nx*nz.
    par.HPC_ncpu = par.nx * par.ny * par.nz
    par.HPC_nnode = round(floor(par.HPC_ncpu / 128)) + 1
    par.HPC_queue = "normal"
    par.HPC_time = "01:00:00"
    par.HPC_account = "EAR22013"
    par.HPC_email = ""

    # Fault 1's own box (ntotft==1 fallback fields other scripts/the netCDF
    # writer's fault-1 slot read directly -- scripts/case.setup's
    # netcdf_write_on_fault_vars, i==0 branch).
    par.nfx = round((FAULT1_XRANGE[1] - FAULT1_XRANGE[0]) / par.dx + 1)
    par.nfz = round((FZMAX - FZMIN) / par.dz + 1)
    par.fx = np.linspace(FAULT1_XRANGE[0], FAULT1_XRANGE[1], par.nfx)
    par.fz = np.linspace(FZMIN, FZMAX, par.nfz)
    par.on_fault_vars = build_on_fault_vars(par.fx, par.fz, par.nfx, par.nfz)

    # Fault 2's own box (independent along-strike extent, same z as fault 1).
    nfx2 = round((FAULT2_XRANGE[1] - FAULT2_XRANGE[0]) / par.dx + 1)
    fx2 = np.linspace(FAULT2_XRANGE[0], FAULT2_XRANGE[1], nfx2)
    fault2_vars = build_on_fault_vars(fx2, par.fz, nfx2, par.nfz)
    par.onFaultVarsPerFault = [par.on_fault_vars, fault2_vars]

    # On-fault stations, PDF Part 6 p.11 (7 on fault #1, 11 on fault #2),
    # code convention (x_km, z_km) with z = -depth_km (down-dip negative).
    fault1_stations = [
        [-10.0, 0.0], [0.0, 0.0], [-10.0, -5.0], [-10.0, -10.0],
        [-5.0, -10.0], [0.0, -10.0], [-10.0, -15.0],
    ]
    fault2_stations = [
        [0.0, 0.0], [5.0, 0.0], [20.0, 0.0], [4.0, -5.0], [5.0, -5.0],
        [0.0, -10.0], [5.0, -10.0], [6.5, -10.0], [10.0, -10.0], [20.0, -10.0],
        [5.0, -15.0],
    ]
    par.st_coor_on_fault = fault1_stations
    par.st_coor_on_fault_per_fault = [fault1_stations, fault2_stations]

    # Off-fault stations, PDF Part 7 p.18-19 (6 stations, all at the surface,
    # perpendicular distance from fault #1: + = far side = code +y).
    par.st_coor_off_fault = [
        [-5.0, 3.0, 0.0], [5.0, 3.0, 0.0], [15.0, 3.0, 0.0],
        [-5.0, -3.0, 0.0], [5.0, -3.0, 0.0], [15.0, -3.0, 0.0],
    ]
    par.n_on_fault = len(par.st_coor_on_fault)
    par.n_off_fault = len(par.st_coor_off_fault)

    return par
