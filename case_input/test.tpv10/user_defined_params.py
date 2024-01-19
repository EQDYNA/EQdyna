#! /usr/bin/env python3

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.dip = 60 # positive dipping angle tiles the fault to y+
# NOTE: for dipping faults, par.fzmin, par.ysource, par.zsource, par.dz, par.dt
# will be modified by mod4dip below. 
 
par.xmin, par.xmax = -30.0e3, 30.0e3
par.ymin, par.ymax = -20.0e3, 40.0e3
par.zmin, par.zmax = -35.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
#par.fymin, par.fymax = 0.0, -9999999999.
par.fzmin, par.fzmax = -15.0e3, 0.0e3

par.xsource = 0.0
#par.ysource = 0.0
par.zsource = -12.0e3 # down dip distance for dipping fault; negative downwards.

par.dx   = 500.
par.dy = par.dx
par.dz = par.dx # will be modified for dipping fault by mod4dip 

par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700
if par.nmat > 1:
    print('Layers for velocity structure is ', par.nmat)
    print('Maxmium and minimum S wave velocity is ', np.max(par.mat[:par.nmat,:],0)[2], np.min(par.mat[:par.nmat,:],0)[2])
    print('Maxmium and minimum P wave velocity is ', np.max(par.mat[:par.nmat,:],0)[1], np.min(par.mat[:par.nmat,:],0)[1])
    print('Maxmium and minimum Density is ', np.max(par.mat[:par.nmat,:],0)[3], np.min(par.mat[:par.nmat,:],0)[3])

par.term = 5.

par.C_elastic = 1
par.insertFaultType = 1 # insert planar dipping fault
#par.dt   = 0.5*par.dz/par.vp
par.friclaw = 1
par.tpv  = 10
par.C_nuclea = 0 # no artifical nulceation; instead, higher initial shear in a patch
par.nucR = 1.5e3

par.dz, par.fzmin, par.ysource, par.zsource, par.dt = mod4dip(par.dip,
    par.dx, par.fzmin, par.zsource, par.vp)
print(' ')
print('DIP: dz is    ',par.dz)
print('DIP: fzmin is ', par.fzmin)
print('DIP: ysource, zsource are' , par.ysource, par.zsource)
print('DIP: dt is    ', par.dt)

par.nx = 2
par.ny = 2
par.nz = 1
par.HPC_ncpu  = par.nx*par.ny*par.nz # Number of CPUs requested.
par.HPC_nnode = round(floor(par.HPC_ncpu/128)) + 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
par.HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
par.HPC_time  = "00:10:00" # WALLTIME, in hh:mm:ss format.
par.HPC_account = "EAR22013" # Project account to be charged SUs against.
par.HPC_email = ""#"dliu@ig.utexas.edu" # Email to receive job status.

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
    # positive initial dip stress to move y+ wall up.
    # negative initial dip stress to move y+ wall down. 
    par.on_fault_vars[iz,ix,49]  = -abs(0.55*par.on_fault_vars[iz,ix,7])       # initial shear stress.
    if abs(xcoor-par.xsource)<=1.5e3 and abs(zcoor-par.zsource)<=1.5e3*sin(abs(par.dip)/180.*pi):
        par.on_fault_vars[iz,ix,49] = -0.2e6 - abs((0.76+0.0057)*par.on_fault_vars[iz,ix,7])
        #par.on_fault_vars[iz,ix,49] = -par.on_fault_vars[iz,ix,49]
    
# (x,z) coordinate pairs for on-fault stations (in km).
C = sin(par.dip/180.*pi)
par.st_coor_on_fault = [[0.0, -1.5*C], [0.0, -3.0*C], [0.0,-4.5*C], [0.0,-7.5*C], 
    [0.0, 0.0], [4.5,0.0], [12.0, 0.0],
    [0.0, -12.0*C], [4.5,-7.5*C], [12.0, -7.5*C]]
print(par.st_coor_on_fault)
# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,1,0], [0,-1,0], [0,2,0], [0,-2,0], 
    [0,3,0], [0,-3,0], [12,3,0], [12,-3,0], [0, 0.5+0.3/tan(abs(par.dip)/180.*pi), -0.3],
    [0, -0.5+0.3/tan(abs(par.dip)/180.*pi), -0.3]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
