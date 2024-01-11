#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -22.0e3, 22.0e3
par.ymin, par.ymax = -10.0e3, 12.0e3
par.zmin, par.zmax = -22.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -15.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -12.0e3

par.dx = 500.
par.dy = par.dx
par.dz = par.dx
par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700

par.term = 5.

par.dt = 0.5*par.dx/par.vp
par.friclaw = 1
par.tpv = 8

par.nucR = 1.5e3
            
# Creating the fault interface
par.nfx = int((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = int((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin,par.fxmax,par.nfx) # coordinates of fault grids along strike.
par.fz  = np.linspace(par.fzmin,par.fzmax,par.nfz) # coordinates of fault grids along dip.

# Create on_fault_vars array for on_fault varialbes.
par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
# functions are defined in lib.py under scripts/
# function lists:
# - shear_steady_state
# - state_steady_state
# - B1, defined in TPV104 and TPV105
# - B2 and B3, defined in TPV105
par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav = 9.8

for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,7]   = 7378.*zcoor     # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = abs(0.55*par.on_fault_vars[iz,ix,7])       # initial shear stress.
    if abs(xcoor-par.xsource)<=1.5e3 and abs(zcoor-par.zsource)<=1.5e3:
        par.on_fault_vars[iz,ix,8] = 1e6+abs(1.005*0.76*par.on_fault_vars[iz,ix,7])
    # tmp1  = B1(xcoor, 15.e3, 3.e3)
    # tmp2  = B1(-zcoor-7.5e3, 15.e3/2., 3.e3)
    # par.on_fault_vars[iz,ix,9]  = par.fric_rsf_a + (1. - tmp1*tmp2)*par.fric_rsf_deltaa
    # par.on_fault_vars[iz,ix,10] = par.fric_rsf_b # assign b in RSF 
    # par.on_fault_vars[iz,ix,11] = par.fric_rsf_Dc # assign Dc in RSF.
    # #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    # #  on_fault_vars[iz,ix,11] = minDc # a special Dc zone.
    # par.on_fault_vars[iz,ix,12] = par.fric_rsf_v0 # initial reference slip rate.
    # par.on_fault_vars[iz,ix,13] = par.fric_rsf_r0 # initial reference friction.
    
    # par.on_fault_vars[iz,ix,14] = par.fric_rsf_fw # 
    # par.on_fault_vars[iz,ix,15] = par.fric_rsf_vw  + par.fric_rsf_deltaavw0*(1. - tmp1*tmp2)  # 
    # par.on_fault_vars[iz,ix,16] = par.fric_tp_a_hy + par.fric_tp_deltaa_hy0*(1. - tmp1*tmp2)  #
    # par.on_fault_vars[iz,ix,17] = par.fric_tp_a_th
    # par.on_fault_vars[iz,ix,18] = par.fric_tp_rouc
    # par.on_fault_vars[iz,ix,19] = par.fric_tp_lambda
    # par.on_fault_vars[iz,ix,40] = par.fric_tp_h
    # par.on_fault_vars[iz,ix,41] = par.fric_tp_Tini
    # par.on_fault_vars[iz,ix,42] = par.fric_tp_pini
    
    # par.on_fault_vars[iz,ix,46] = par.creep_slip_rate # initial slip rates
    # #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    # #  on_fault_vars[iz,ix,46] = 0.03 # initial high slip rate patch.
    
    # par.on_fault_vars[iz,ix,20] = state_steady_state(par.on_fault_vars[iz,ix,9], 
                                                # par.on_fault_vars[iz,ix,10],
                                                # par.on_fault_vars[iz,ix,11],
                                                # par.on_fault_vars[iz,ix,12],
                                                # par.on_fault_vars[iz,ix,13],
                                                # par.on_fault_vars[iz,ix,8],
                                                # par.on_fault_vars[iz,ix,7],
                                                # par.on_fault_vars[iz,ix,46],
                                                # par.friclaw) # initial state var.

