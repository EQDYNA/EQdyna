#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -30.0e3, 30.0e3
par.ymin, par.ymax = -28.0e3, 32.0e3
par.zmin, par.zmax = -34.0e3, 0.0e3

par.fxmin, par.fxmax = -6.0e3, 6.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -8.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -3.4e3

par.dx = 400.
par.nuni_y_plus=10
par.nuni_y_minus=10
par.enlarging_ratio = 1.025
par.nmat = 5
par.n2mat = 4
par.mat = np.zeros((par.nmat, par.n2mat))
# example here 1D velocity structure for the Cushing earthquake.
par.mat[0,:] = [1.19e3,  2.74e3, 1.45e3, 2.1e3] # top layer
par.mat[1,:] = [2.01e3,  5.75e3, 3.06e3, 2.4e3] # 2nd layer going downwards into the earth.
par.mat[2,:] = [4.94e3,  5.72e3, 3.4e3,  2.6e3] # 3rd
par.mat[3,:] = [10.94e3, 6.18e3, 3.62e3, 2.8e3] # 4th
par.mat[4,:] = [1.e6,    6.32e3, 3.67e3, 2.8e3] # rest



par.term = 5.
par.dt = 0.5*par.dx/6.32e3
par.friclaw = 2
par.tpv = 201

par.nucR = 25.e3
par.nucRuptVel = 1.5e3
par.nucdtau0 = -9999.
            
# Creating the fault interface
par.nfx = int((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = int((par.fzmax - par.fzmin)/par.dx + 1)
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
par.fric_sw_fs = 0.4
par.fric_sw_sd = 0.3
par.fric_sw_D0 = 0.3
par.fric_tw_t0 = 0.2
par.grav = 9.8

for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = 1000.
    par.on_fault_vars[iz,ix,2]   = 1000.
    if abs(xcoor)<2.25e3 and zcoor<-3.4e3+2.25e3 and zcoor>-3.4e3-2.25e3:
        tmp1 = 1.
        tmp2 = 1.
        if abs(xcoor)>=1.25e3:
            tmp1 = (2.25e3-abs(xcoor))/1.0e3
        if abs(zcoor--3.4e3)>=1.25e3:
            tmp2 = (2.25e3-abs(zcoor--3.4e3))/1.0e3
        par.on_fault_vars[iz,ix,1] = par.fric_sw_fs
        par.on_fault_vars[iz,ix,2] = par.fric_sw_fs - 0.1*tmp1*tmp2
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,5]   = par.fric_tw_t0
    par.on_fault_vars[iz,ix,7]   = -100.0e6     # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = 35.0e6       # initial shear stress.
   
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
    
    
print(dir(par))
