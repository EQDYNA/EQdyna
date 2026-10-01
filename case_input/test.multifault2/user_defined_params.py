#! /usr/bin/env python3
"""
test.multifault2 -- Row 17 (multi-fault) plumbing case.

NOT a SCEC TPV benchmark and NOT tpv22/tpv23 (those are a separate, later,
reviewed mission -- pathway_forward.md item 17). This is a cheap two-fault
ROUTING/SMOKE case: two vertical, planar, PARALLEL strike-slip faults,
offset in y (the fault-normal direction) by 2000 m (4 * dy, landing exactly
on a mesh node line -- see checkInputConsistency.f90's new multi-fault
guards), same x/z extent (fault 1's uniform x/z mesh belt covers both).

Deliberately gives the two faults a DIFFERENT initial normal-stress
coefficient (7378 vs 5000 Pa/m) so a fault-2-reads-fault-1's-array bug (the
exact failure mode eqquasi's own multi-fault work found, and the one this
case exists to catch) shows up immediately in frt.txt's tnrm column, not
just in whether the run completes. Only fault 1 gets a forced-nucleation
patch (nucfault=1); fault 2 is a plumbing check (does it mesh, load, and
integrate correctly), not a rupture-physics case -- kept sub-critical on
purpose, short term (few tens of steps).

Based on case_input/test.tpv8 (nx=2, ny=2, nz=1 -> 4 ranks, dx=500).
"""

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.ntotft = 2

par.xmin, par.xmax = -22.0e3, 22.0e3
par.ymin, par.ymax = -10.0e3, 12.0e3
par.zmin, par.zmax = -22.0e3, 0.0e3

# fault 1's box also doubles as par.fxmin/fxmax/fymin/fymax/fzmin/fzmax, the
# ntotft==1 fallback fields other scripts (e.g. bFault_Rough_Geometry tools)
# may still read; keep it equal to faultgeom[0].
par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -15.0e3, 0.0e3

# Row 17: per-fault boxes (resolveFaultGeom, scripts/lib.py). Fault 2 is
# fault 1 translated by +2000 m in y (4 * dy = 4 * 500 m), same x/z extent
# -- checkInputConsistency.f90 refuses anything else for this release.
FAULT2_Y = 2.0e3
par.faultgeom = [
    (par.fxmin, par.fxmax, par.fymin, par.fymax, par.fzmin, par.fzmax),
    (par.fxmin, par.fxmax, FAULT2_Y,  FAULT2_Y,  par.fzmin, par.fzmax),
]

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -12.0e3

par.dx = 500.
par.dy = par.dx
par.dz = par.dx
par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700

# Short: this is a routing/plumbing smoke case, not a rupture-physics one.
par.term = 1.0

par.dt = 0.5*par.dx/par.vp
par.friclaw = 1
par.tpv = 8   # reuse tpv8's swtwNucleation branch; not claiming to BE tpv8

par.nucR = 1.5e3
par.nucfault = 1   # only fault 1 is nucleated

par.nx = 2
par.ny = 2
par.nz = 1

# Creating the fault interface -- fault 1's grid (also par.fx/par.fz/par.nfx/
# par.nfz/par.on_fault_vars, the ntotft==1 fallback names
# netcdf_write_on_fault_vars uses for fault 1 specifically).
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)

par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav = 9.8
par.fric_cohesion = 1.e6


def _build_on_fault_vars(fx, fz, nfx, nfz, norm_stress_coeff, nucleate):
    """One fault's on_fault_vars array (par's existing (nfz, nfx, 100)
    layout). norm_stress_coeff is Pa/m depth (negative*z -> compressive);
    fault 1 uses 7378 (test.tpv8's value), fault 2 uses 5000 -- deliberately
    DIFFERENT so an aliasing bug (fault 2 silently reading fault 1's array)
    is visible in the output, not just in whether the run completes."""
    v = np.zeros((nfz, nfx, 100))
    for ix, xcoor in enumerate(fx):
        for iz, zcoor in enumerate(fz):
            v[iz, ix, 1] = par.fric_sw_fs
            if abs(abs(xcoor) - fx[-1]) < 0.01 or abs(zcoor - fz[0]) < 0.01:
                v[iz, ix, 1] = 1000.
            v[iz, ix, 2] = par.fric_sw_fd
            v[iz, ix, 3] = par.fric_sw_D0
            v[iz, ix, 4] = par.fric_cohesion
            v[iz, ix, 7] = norm_stress_coeff * zcoor       # initial normal stress, Pa (negative, compressive)
            v[iz, ix, 8] = abs(0.55 * v[iz, ix, 7])        # initial shear stress
            if nucleate and abs(xcoor - par.xsource) <= 1.5e3 and abs(zcoor - par.zsource) <= 1.5e3:
                v[iz, ix, 8] = 1e6 + abs(1.005 * 0.76 * v[iz, ix, 7])
    return v


par.on_fault_vars = _build_on_fault_vars(par.fx, par.fz, par.nfx, par.nfz,
                                          norm_stress_coeff=7378., nucleate=True)
fault2_vars = _build_on_fault_vars(par.fx, par.fz, par.nfx, par.nfz,
                                    norm_stress_coeff=5000., nucleate=False)
# Row 17: per-fault on_fault_vars list (resolveOnFaultVarsPerFault,
# scripts/lib.py) -- fault 1's array here MUST be par.on_fault_vars itself
# (not a recomputed copy), so case.setup's ntotft==1 fallback path and this
# path write bit-identical fault-1 output.
par.onFaultVarsPerFault = [par.on_fault_vars, fault2_vars]

# Row 17: per-fault on-fault station lists (resolveOnFaultStationsPerFault).
# Fault 1 keeps test.tpv8's station set; fault 2 gets a smaller, distinct set
# so the two faults' station files cannot be confused for one another.
par.st_coor_on_fault = [[0.0, 0.0], [0.0, -4.5], [0.0, -7.5], [0.0, -12.0], [4.5, -7.5],
                        [12.0, -7.5], [4.5, 0.0], [12.0, 0.0]]
par.st_coor_on_fault_per_fault = [
    par.st_coor_on_fault,
    [[0.0, -7.5], [4.5, -7.5]],
]

# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0, 1, 0], [0, -1, 0], [0, 2, 0], [0, -2, 0], [0, 3, 0], [0, -3, 0],
                         [12, 6, 0], [12, -3, 0], [-12, 3, 0], [0, -0.5, -0.3], [0, 0.5, -0.3],
                         [0, -1, -0.3], [0, 1, -0.3], [12, -3, -12], [12, 3, -12]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
