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

MESH-SHARING CONSTRAINT (read, not guessed): checkInputConsistency.f90's row-17
multi-fault guards (ERR_GEOM_MULTIFAULT_XZ_BAD) require EVERY fault to share
EXACTLY fault 1's x/z box -- the uniform x/z mesh belt is built once, from
fault 1's box alone (checkInputConsistency.f90:69-81). TPV22/23's real
geometry has fault #2 offset 20 km ALONG STRIKE from fault #1 (only a 10 km
overlap out of each fault's 30 km length), so the two faults' TRUE x-extents
differ -- exactly the case this guard refuses outright.

Fix used here (not a physics shortcut): both faults are given the SAME shared
box, the UNION of their true along-strike extents, x in [-25000, 25000]
(50 km). Each fault's on_fault_vars then marks the portion of that shared box
OUTSIDE its own true 30 km extent as UNBREAKABLE (sw_fs = 1000, the same
"unbreakable border" technique already used in case_input/test.tpv29 and
test.tpv8 for domain-edge nodes) -- those extension nodes are mesh-only
scaffolding, never load-bearing physics, identical in spirit to the existing
border treatment, just extended over a wider inert region instead of a single
row. Each fault's TRUE borders (its own real left/right edge, and bottom
[z=-20000]) get the same unbreakable treatment, matching PDF Part 3, p.5:
"Slip goes to zero at the border of a fault ... a node which lies precisely
on the border of a fault should not be permitted to slip. This is a CHANGE
from earlier benchmarks". REVISED 2026-10-01 (mira/row17-tpv2223-rebased
iteration log): the TOP edge (z=0, the surface trace) is deliberately NOT
included, reversing this module's original reading. Both independent
references (kaneko/SPECFEM3D, payne/EQdyna) show large early slip at the
station literally named "0 km down-dip" (onset ~5.35s, peak 2.5-2.8 m, both
tpv22 and tpv23) -- incompatible with that node being unbreakable. See
build_on_fault_vars' inline comment for the full evidence and the remaining
ambiguity in the spec text.

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

# Shared geometry, in CODE axes (x=strike, y=stepover/fault-normal, z=depth<=0)
FXMIN, FXMAX = -25000.0, 25000.0   # UNION of fault1 [-25000,5000] and fault2 [-5000,25000]
FZMIN, FZMAX = -20000.0, 0.0       # both faults: 20 km deep, reaching the surface
FAULT1_XRANGE = (-25000.0, 5000.0)   # fault #1's TRUE along-strike extent
FAULT2_XRANGE = (-5000.0, 25000.0)   # fault #2's TRUE along-strike extent

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


def build_on_fault_vars(fx, fz, nfx, nfz, true_xrange, is_nucleation_fault):
    """One fault's on_fault_vars array, built on the SHARED box (fx/fz span
    the union [-25000,25000] x [-20000,0] for both faults -- see module
    docstring). true_xrange is this fault's own real along-strike extent;
    nodes outside it are UNBREAKABLE (mesh-sharing scaffolding, never ruptures).
    Nucleation (swtwNucleation, TPV==22/23 branch) only acts on nodes where
    C_nuclea==1 AND this fault is par.nucfault -- is_nucleation_fault is
    informational only here (no manual stress patch needed, unlike the older
    "shear bump" style of test.tpv8/test.multifault2; TPV22/23 is
    forced-rupture-only, same as test.tpv29)."""
    v = np.zeros((nfz, nfx, 100))
    xlo, xhi = true_xrange
    for ix, xcoor in enumerate(fx):
        active = (xcoor >= xlo - 0.5) and (xcoor <= xhi + 0.5)
        for iz, zcoor in enumerate(fz):
            depth = -zcoor  # zcoor <= 0; depth positive-down
            # ITERATION (mira/row17-tpv2223-rebased, 2026-10-01): scaffold
            # (mesh-sharing, outside this fault's true x-extent) nodes used to
            # fall through here with sigma_n=tau_ini=cohesion=0 (the
            # np.zeros() default) and ONLY the friction coefficient forced to
            # UNBREAKABLE_FS. That does not actually make them unbreakable:
            # shear strength = C0 + mu*max(0,sigma_n) = 0 + 1000*max(0,0) = 0
            # regardless of mu, so these nodes had ZERO frictional strength
            # and would slip freely under the slightest incident shear stress
            # -- found by inspection (not yet isolated by a measurement
            # against reference data; flagged as a found-but-unconfirmed
            # defect, see report). test.tpv29's own border treatment
            # (user_defined_params.py, sw_fs=10000 at ix==0/nfx-1/iz==0) does
            # NOT skip the physical stress assignment for border nodes --
            # every node gets the real in-situ stress state, and ONLY the
            # friction coefficient is forced absurdly high. Matched here: the
            # scaffold region now carries the SAME physical stress as an
            # active node at that (x,z); only friction is forced unbreakable.
            v[iz, ix, 1] = MU_S if active else UNBREAKABLE_FS
            v[iz, ix, 2] = MU_D
            v[iz, ix, 3] = D0
            v[iz, ix, 4] = cohesion(depth)
            v[iz, ix, 5] = TW_T0
            v[iz, ix, 7] = -SIGMA_INI               # init normal stress, negative=compressive
            v[iz, ix, 8] = tau_ini(depth)            # init strike (right-lateral) shear
            # True borders: left/right edge of THIS fault's own extent, and
            # bottom (z=fzmin). ITERATION (mira/row17-tpv2223-rebased,
            # 2026-10-01): top edge (z=0, the surface trace) REMOVED from the
            # border set -- both independent references (kaneko/SPECFEM3D,
            # payne/EQdyna) show large early slip (onset ~5.35s, peak
            # 2.5-2.8 m) at the station literally named "0 km down-dip" for
            # BOTH tpv22 and tpv23, which an unbreakable z=0 node cannot
            # produce (measured: our run showed exactly 0 slip there at
            # term=15s with top-edge included). Barall's own DayFD header
            # reports its "0 km down-dip" station's ACTUAL node at z=-800 m,
            # not 0 -- i.e. independent codes do not sample the literal
            # geometric border either; whether the spec's "a node which lies
            # precisely on the border... should not slip" is meant to apply to
            # the free surface trace at all is genuinely ambiguous from the
            # text alone, and the measured reference behavior is the
            # tie-breaker used here. Reverts the docstring's prior reading
            # (originally cited PDF p.5's "CHANGE from earlier benchmarks" as
            # covering top+bottom; evidence says bottom only matches, top does
            # not). Bottom and left/right are UNCHANGED (no reference evidence
            # against them).
            on_left_edge = abs(xcoor - xlo) < 0.5
            on_right_edge = abs(xcoor - xhi) < 0.5
            on_bottom_edge = abs(zcoor - FZMIN) < 0.5
            if on_left_edge or on_right_edge or on_bottom_edge:
                v[iz, ix, 1] = UNBREAKABLE_FS
    return v


def build_params(tpv, fault2_z):
    """fault2_z is the stepover offset in CODE y (== spec z): -1600 (TPV22,
    extensional right step) or +1000 (TPV23, compressional left step), taken
    directly from PDF Parts 1/2 (p.3-4)."""
    par = parameters()
    par.ntotft = 2

    # Domain margins beyond the fault union box (~7 km in x/z, matching
    # test.tpv8's margin style; y margin generous since the stepover offset
    # is tiny, <=1.6 km, relative to the domain).
    par.xmin, par.xmax = -32.0e3, 32.0e3
    par.ymin, par.ymax = -10.0e3, 10.0e3
    par.zmin, par.zmax = -27.0e3, 0.0e3

    # ntotft==1 fallback fields (other scripts may still read these for
    # fault 1); fault 1's box is also the shared union box (see docstring).
    par.fxmin, par.fxmax = FXMIN, FXMAX
    par.fymin, par.fymax = 0.0, 0.0
    par.fzmin, par.fzmax = FZMIN, FZMAX

    par.faultgeom = [
        (FXMIN, FXMAX, 0.0, 0.0, FZMIN, FZMAX),
        (FXMIN, FXMAX, fault2_z, fault2_z, FZMIN, FZMAX),
    ]

    # Hypocenter: spec (x,y,z)=(-10000,10000,0) -> code (x,y,z)=(-10000,0,-10000)
    par.xsource, par.ysource, par.zsource = -10000.0, 0.0, -10000.0

    # Gate-coarse resolution (minutes at 4 ranks -- NOT the spec's 100/50 m;
    # full_specs.py carries the spec resolution for later, scheduled runs).
    # dx/dz (strike/depth, in-plane fault resolution) independent of dy (the
    # cross-fault/volume-mesh spacing) -- same pattern test.tpv36/37 already
    # use (par.dy = par.dx*cos(dip) there). dy must be an integer multiple of
    # the stepover offset (checkInputConsistency.f90 ERR_GEOM_MULTIFAULT_Y_BAD):
    # 1600 m (TPV22) factors as 4*400; 1000 m (TPV23) factors as 2*500.
    par.dx = 1000.0
    par.dz = 1000.0
    par.dy = abs(fault2_z) / 4.0 if tpv == 22 else abs(fault2_z) / 2.0

    par.nmat = 1
    par.vp, par.vs, par.rou = 6000.0, 3464.0, 2670.0

    par.term = 15.0  # ITERATION (mira/row17-tpv2223-rebased, 2026-10-01): spec's
    # own full duration (TPV22_23_Description_v08.pdf p.7, "Run the model for
    # times from 0.0 to 15.0 seconds"); GATE_TERM_S (5.0) is too short even for
    # fault #1's own far stations (fault1st000dp000 onset ~5.37s in both
    # independent references) -- never mind fault #2 (~9-12s). NOT the
    # testsys gate value; this file is pre-gate physics investigation only
    # (see dispatch brief), term is restored before this becomes a gate case.
    # CFL: dt must be set from the SMALLEST element dimension, not dx --
    # test.tpv36/37 (also dx != dy != dz) use dz for exactly this reason when
    # dz is the small one. Here dy (400 m TPV22 / 500 m TPV23) is smaller
    # than dx=dz=1000 m; using 0.5*dx/vp (that is, ignoring dy) blew up the
    # run within the first few steps (IEEE_DIVIDE_BY_ZERO, slip ~1e37 m) --
    # measured directly, not inferred -- before this fix.
    par.dt = 0.5 * min(par.dx, par.dy, par.dz) / par.vp

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

    par.nfx = round((FXMAX - FXMIN) / par.dx + 1)
    par.nfz = round((FZMAX - FZMIN) / par.dz + 1)
    par.fx = np.linspace(FXMIN, FXMAX, par.nfx)
    par.fz = np.linspace(FZMIN, FZMAX, par.nfz)

    par.on_fault_vars = build_on_fault_vars(
        par.fx, par.fz, par.nfx, par.nfz, FAULT1_XRANGE, is_nucleation_fault=True)
    fault2_vars = build_on_fault_vars(
        par.fx, par.fz, par.nfx, par.nfz, FAULT2_XRANGE, is_nucleation_fault=False)
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
