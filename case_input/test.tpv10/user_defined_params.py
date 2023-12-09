#! /usr/bin/env python3

from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()
par.dip = -60 # negative dipping angle tiles the fault to y+

par.xmin, par.xmax = -22.0e3, 22.0e3
par.ymin, par.ymax = -10.0e3, 20.0e3
par.zmin, par.zmax = -22.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
#par.fymin, par.fymax = 0.0, -9999999999.
par.fzmin, par.fzmax = -15.0e3*sin(abs(par.dip)/180.*pi), 0.0e3

par.xsource = 0.0
par.ysource = 12.0e3*cos(abs(par.dip)/180.*pi)
par.zsource = -12.0e3*sin(abs(par.dip)/180.*pi)

par.dx   = 500.
par.dy = par.dx
par.dz = par.dx*sin(abs(par.dip)/180.*pi)

par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700

par.term = 5.

par.C_elastic = 1
par.insertFaultType = 1 # insert planar dipping fault
par.dt   = 0.5*par.dz/par.vp
par.friclaw = 1
par.tpv  = 10

par.C_nuclea = 0 # no artifical nulceation; instead, higher initial shear in a patch
par.nucR = 1.5e3
            
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
par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav       = 9.8
par.fric_cohesion = 0.2e6
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    downDipDistance = abs(zcoor)/sin(abs(par.dip)/180.*pi)
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4]   = par.fric_cohesion
    par.on_fault_vars[iz,ix,7]   = -7378.*downDipDistance # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = 0.
    par.on_fault_vars[iz,ix,49]  = abs(0.55*par.on_fault_vars[iz,ix,7])       # initial shear stress.
    if abs(xcoor-par.xsource)<=1.5e3 and abs(zcoor-par.zsource)<=1.5e3*cos(abs(par.dip)/180.*pi):
        par.on_fault_vars[iz,ix,49] = 0.2e6 + abs((0.76+0.0057)*par.on_fault_vars[iz,ix,7])
        #par.on_fault_vars[iz,ix,49] = -par.on_fault_vars[iz,ix,49]
    
print(dir(par))
