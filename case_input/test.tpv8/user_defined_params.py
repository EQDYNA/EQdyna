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

par.nx = 2
par.ny = 2
par.nz = 1
            
# Creating the fault interface
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
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
par.fric_cohesion = 1.e6 
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
  # assign a in RSF. a is a 2D distribution.
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4]   = par.fric_cohesion
    par.on_fault_vars[iz,ix,7]   = 7378.*zcoor     # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = abs(0.55*par.on_fault_vars[iz,ix,7])       # initial shear stress.
    if abs(xcoor-par.xsource)<=1.5e3 and abs(zcoor-par.zsource)<=1.5e3:
        par.on_fault_vars[iz,ix,8] = 1e6+abs(1.005*0.76*par.on_fault_vars[iz,ix,7])

par.st_coor_on_fault = [[0.0, 0.0], [0.0, -4.5], [0.0,-7.5], [0.0, -12.0], [4.5,-7.5], \
   [12.0, -7.5], [4.5, 0.0], [12.0, 0.0]]
   
# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,1,0], [0,-1,0], [0,2,0], [0,-2,0], [0,3,0], [0,-3,0], [12,6,0], \
                     [12,-3,0], [-12,3,0], [0,-0.5,-0.3], [0,0.5,-0.3], \
                     [0,-1,-0.3], [0,1,-0.3],[12,-3,-12], [12,3,-12]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)