#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

# input parameters for drv.a6

par = parameters()

par.xmin, par.xmax = -40.0e3, 40.0e3
par.ymin, par.ymax = -31.0e3, 30.0e3
par.zmin, par.zmax = -40.0e3, 0.0e3

par.fxmin, par.fxmax = -25.0e3, 25.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -25.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = -12.0e3, 0.0, -12.0e3

par.C_elastic = 0
#par.C_nuclea  = 1
#par.C_degen   = 0
par.friclaw   = 4
par.rough_fault = 1
par.tpv       = 2802
par.enlarging_ratio = 1.
par.outputGroundMotion = 0
par.output_plastic     = 1

par.nmat, par.n2mat = 1,3
par.vp, par.vs, par.rou = 6.e3, 3.464e3, 2.67e3

par.term    = 5.
par.dx      = 500.
par.dt      = 0.5*par.dx/par.vp

par.nucR = 4.0e3
par.nucdtau0 = 50.0e6

par.nx, par.ny, par.nz = 2, 2, 1

par.fric_rsf_a = 0.01
par.fric_rsf_deltaa = 0.01
par.fric_rsf_b = 0.014
par.fric_rsf_Dc = 0.4
par.fric_rsf_r0 = 0.6
par.fric_rsf_v0 = 1e-6
par.fric_rsf_fw = 0.2
par.fric_rsf_vw = 0.1
par.fric_rsf_deltavw0 = 0.9

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
  
def getScaleCoeff(coorX,coorZ,wwx,wx,wwz1,wz1,wwz2,wz2):
    # wwx, symmetric about zero, where to start taper.
    # wx, length for tapering
    # wwz1, top edge to start tapering.
    # wwz2, bottom edge to start tapering.
    # wz1, top length for tapering.
    # wz2, bottom length for tapering.
    scaleCoeffX, scaleCoeffZ = 1., 1.
    if abs(coorX)<=wwx:
        scaleCoeffX = 1.
    elif abs(coorX)>wwx and abs(coorX)<(wwx+wx):
        scaleCoeffX = 1. - (abs(coorX)-wwx)/wx
    elif abs(coorX)>=(wwx+wx):
        scaleCoeffX = 0.
    
    if abs(coorZ)<=(wwz1-wz1) or abs(coorZ)>=(wwz2+wz2):
        scaleCoeffZ = 0.
    elif abs(coorZ)<=wwz1 and abs(coorZ)>=(wwz1-wz1):
        scaleCoeffZ = 1. - (wwz1 - abs(coorZ))/wz1
    elif abs(coorZ)<=(wwz2+wz2) and abs(coorZ)>=wwz2:
        scaleCoeffZ = 1. - (abs(coorZ) - wwz2)/wz2
    if abs(coorZ)>=wwz1 and abs(coorZ)<=wwz2:
        scaleCoeffZ = 1.
    return scaleCoeffX, scaleCoeffZ
    
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs 
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,7]   = -120.0e6     # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = 40.0e6       # initial shear stress.
    
    tmp1, tmp2 = getScaleCoeff(xcoor,zcoor,20.e3,2.e3,2.e3,2.e3,14.e3,1.e3)
    deltaa_tmp = par.fric_rsf_deltaa
    if abs(zcoor)>15.e3:
        deltaa_tmp = 2.*par.fric_rsf_deltaa
        
    par.on_fault_vars[iz,ix,9]  = par.fric_rsf_a + (1. - tmp1*tmp2)*deltaa_tmp
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
    
    
print(dir(par))
