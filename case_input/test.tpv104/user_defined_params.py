#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -38.0e3, 38.0e3
par.ymin, par.ymax = -20.0e3, 22.0e3
par.zmin, par.zmax = -38.0e3, 0.0e3

par.fxmin, par.fxmax = -18.0e3, 18.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -18.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -7.5e3

par.term = 5.
par.dx = 500.
par.dy = par.dx
par.dz = par.dx

par.dt = 0.5*par.dx/par.vp
par.friclaw = 4
par.tpv = 104

par.nucR = 3.0e3
par.nucdtau0 = 45.0e6

par.nx = 2
par.ny = 2
par.nz = 1
            
# Creating the fault interface
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dx + 1)
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
  
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs 
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,7]   = -120.0e6     # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = 40.0e6       # initial shear stress.
    
    tmp1  = B1(xcoor, 15.e3, 3.e3)
    tmp2  = B1(-zcoor-7.5e3, 15.e3/2., 3.e3)
    par.on_fault_vars[iz,ix,9]  = par.fric_rsf_a + (1. - tmp1*tmp2)*par.fric_rsf_deltaa
    par.on_fault_vars[iz,ix,10] = par.fric_rsf_b # assign b in RSF 
    par.on_fault_vars[iz,ix,11] = par.fric_rsf_Dc # assign Dc in RSF.
    #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    #  on_fault_vars[iz,ix,11] = minDc # a special Dc zone.
    par.on_fault_vars[iz,ix,12] = par.fric_rsf_v0 # initial reference slip rate.
    par.on_fault_vars[iz,ix,13] = par.fric_rsf_r0 # initial reference friction.
    
    par.on_fault_vars[iz,ix,14] = par.fric_rsf_fw # 
    par.on_fault_vars[iz,ix,15] = par.fric_rsf_vw  + par.fric_rsf_deltaavw0*(1. - tmp1*tmp2)  # 
    par.on_fault_vars[iz,ix,16] = par.fric_tp_a_hy + par.fric_tp_deltaa_hy0*(1. - tmp1*tmp2)  #
    par.on_fault_vars[iz,ix,17] = par.fric_tp_a_th
    par.on_fault_vars[iz,ix,18] = par.fric_tp_rouc
    par.on_fault_vars[iz,ix,19] = par.fric_tp_lambda
    par.on_fault_vars[iz,ix,40] = par.fric_tp_h
    par.on_fault_vars[iz,ix,41] = par.fric_tp_Tini
    par.on_fault_vars[iz,ix,42] = par.fric_tp_pini
    
    par.on_fault_vars[iz,ix,46] = par.creep_slip_rate # initial slip rates
    #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    #  on_fault_vars[iz,ix,46] = 0.03 # initial high slip rate patch.
    
    par.on_fault_vars[iz,ix,20] = state_steady_state(par.on_fault_vars[iz,ix,9], 
                                                par.on_fault_vars[iz,ix,10],
                                                par.on_fault_vars[iz,ix,11],
                                                par.on_fault_vars[iz,ix,12],
                                                par.on_fault_vars[iz,ix,13],
                                                par.on_fault_vars[iz,ix,8],
                                                par.on_fault_vars[iz,ix,7],
                                                par.on_fault_vars[iz,ix,46],
                                                par.friclaw) # initial state var.

