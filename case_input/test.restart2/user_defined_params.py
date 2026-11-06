#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

# test.restart2 -- the mode = 2 (EQsimu restart) end-to-end gate, board row 16c.
# The fault state (stresses, state variable, slip rate, node velocities: fric
# slots 8, 49, 7, 47, 20, 23, 31-36) is read from the committed real
# fault.r.nc in this directory; see README.md for its provenance.

par = parameters()

par.mode = 2

par.xmin, par.xmax = -16.0e3, 16.0e3
par.ymin, par.ymax = -6.0e3, 6.0e3
par.zmin, par.zmax = -12.0e3, 0.0e3

# fault footprint = the fixture's grid: 81 x 17 nodes at dx = 250 m
par.fxmin, par.fxmax = -10.0e3, 10.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -4.0e3, 0.0e3

par.dx = 250.
par.dy = par.dx
par.dz = par.dx
par.nuni_y_plus, par.nuni_y_minus = 5, 5
par.enlarging_ratio = 1.1

par.C_elastic = 1
par.C_nuclea  = 0       # no artificial nucleation: the restart state drives the run
par.friclaw   = 3       # rate-and-state, aging law -- the law the source EQquasi run used
par.ntotft    = 1
par.insertFaultType = 0
par.tpv       = 1000
par.output_plastic = 0
par.outputGroundMotion = 0

par.term = 5.
par.dt   = 0.5*par.dx/par.vp

par.nx, par.ny, par.nz = 2, 1, 2

# BP5 rate-and-state parameters of the source EQquasi run (bp1001.qd.250)
par.fric_rsf_a      = 0.004
par.fric_rsf_b      = 0.03
par.fric_rsf_Dc     = 0.14
par.fric_rsf_deltaa = 0.036
par.fric_rsf_r0     = 0.6
par.fric_rsf_v0     = 1e-6
par.creep_slip_rate = 1e-9

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin, par.fxmax, par.nfx)
par.fz  = np.linspace(par.fzmin, par.fzmax, par.nfz)
par.on_fault_vars = np.zeros((par.nfz, par.nfx, 100))

# Per-node inputs. With mode = 2 the stresses, state variable and velocities
# (slots 7, 8, 20 here) are overwritten by the restart read; the values below
# are what a mode = 1 twin of this case starts from.
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    par.on_fault_vars[iz,ix,1]  = par.fric_sw_fs
    par.on_fault_vars[iz,ix,2]  = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]  = par.fric_sw_D0
    par.on_fault_vars[iz,ix,7]  = -25.0e6
    par.on_fault_vars[iz,ix,8]  = 13.25e6
    # a(x,z): the source run's BP5 map
    if abs(zcoor)>=18e3 or abs(xcoor)>=20e3 or abs(zcoor)<=2e3:
        par.on_fault_vars[iz,ix,9] = par.fric_rsf_a + par.fric_rsf_deltaa
    elif abs(zcoor)<=16e3 and abs(zcoor)>=4e3 and abs(xcoor)<=18e3:
        par.on_fault_vars[iz,ix,9] = par.fric_rsf_a
    else:
        tmp1 = (abs(abs(zcoor)-10e3) - 6e3)/2e3
        tmp2 = (abs(xcoor)-18e3)/2e3
        par.on_fault_vars[iz,ix,9] = par.fric_rsf_a + max(tmp1,tmp2)*par.fric_rsf_deltaa
    par.on_fault_vars[iz,ix,10] = par.fric_rsf_b
    par.on_fault_vars[iz,ix,11] = par.fric_rsf_Dc
    par.on_fault_vars[iz,ix,12] = par.fric_rsf_v0
    par.on_fault_vars[iz,ix,13] = par.fric_rsf_r0
    par.on_fault_vars[iz,ix,46] = par.creep_slip_rate
    par.on_fault_vars[iz,ix,20] = state_steady_state(par.on_fault_vars[iz,ix,9],
                                                par.on_fault_vars[iz,ix,10],
                                                par.on_fault_vars[iz,ix,11],
                                                par.on_fault_vars[iz,ix,12],
                                                par.on_fault_vars[iz,ix,13],
                                                par.on_fault_vars[iz,ix,8],
                                                par.on_fault_vars[iz,ix,7],
                                                par.on_fault_vars[iz,ix,46],
                                                par.friclaw)

par.st_coor_on_fault = [[0.0, -0.5], [0.0, -2.0], [5.0, -1.5]]
par.st_coor_off_fault = [[0, 2, 0], [5, -3, 0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
