#! /usr/bin/env python3

from defaultParameters import *
from math import *
from lib import *
import numpy as np

par = parameters()

par.dip = 15 # positive dipping angle tiles the fault to y+
# NOTE: for dipping faults, par.fzmin, par.ysource, par.zsource, par.dz, par.dt
# will be modified by mod4dip below. 
faultWidth = 28.e3 # for a shallow dipping fault, convenient to discuss properties in the strike-dip map.

par.xmin, par.xmax = -20.0e3, 20.0e3
par.ymin, par.ymax = -10.0e3, faultWidth*cos(abs(par.dip)/180.*pi)+10.0e3
par.zmin, par.zmax = -15.0e3, 0.0e3

par.fxmin, par.fxmax = -15.0e3, 15.0e3
par.fymin, par.fymax = 0.0, faultWidth*cos(abs(par.dip)/180.*pi)
par.fzmin, par.fzmax = -faultWidth*sin(abs(par.dip)/180.*pi), 0.0e3

par.xsource = 0.0
par.ysource = 18.e3*cos(abs(par.dip)/180.*pi)
par.zsource = -18.e3*sin(abs(par.dip)/180.*pi) # down dip distance for dipping fault; negative downwards.

# Cell size on fault will be uniform along strike and dip.
par.dx   = 500.
par.dy = par.dx*cos(abs(par.dip)/180.*pi)
par.dz = par.dx*sin(abs(par.dip)/180.*pi) # will be modified for dipping fault by mod4dip 

par.nuni_y_plus = round(faultWidth/par.dx)+5
par.nuni_y_minus = 5

par.nmat = 1
par.vp, par.vs, par.rou = 6000, 3464, 2670
if par.nmat > 1:
    print('Layers for velocity structure is ', par.nmat)
    print('Maxmium and minimum S wave velocity is ', np.max(par.mat[:par.nmat,:],0)[2], np.min(par.mat[:par.nmat,:],0)[2])
    print('Maxmium and minimum P wave velocity is ', np.max(par.mat[:par.nmat,:],0)[1], np.min(par.mat[:par.nmat,:],0)[1])
    print('Maxmium and minimum Density is ', np.max(par.mat[:par.nmat,:],0)[3], np.min(par.mat[:par.nmat,:],0)[3])

par.term = 6.

par.C_elastic = 1
par.C_degen   = par.dip # degenrate hexahedrons to wedges
par.insertFaultType = 0 # don't insert planar dipping fault because of the low dipping angle
par.outputFinalSurfDisp = 1 # write out final surface disp.
par.dt = 0.5*par.dz/par.vp
par.friclaw = 1
par.tpv  = 37
par.C_nuclea = 1 # artifical nulceation; instead, higher initial shear in a patch
par.nucR = 4e3

print(' ')
print('DIP: dz is    ',par.dz)
print('DIP: fzmin is ', par.fzmin)
print('DIP: ysource, zsource are' , par.ysource, par.zsource)
print('DIP: dt is    ', par.dt)

par.nx = 2
par.ny = 2
par.nz = 2
par.HPC_ncpu  = par.nx*par.ny*par.nz # Number of CPUs requested.
par.HPC_nnode = round(floor(par.HPC_ncpu/128)) + 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
par.HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
par.HPC_time  = "02:00:00" # WALLTIME, in hh:mm:ss format.
par.HPC_account = "EAR22012" # Project account to be charged SUs against.
par.HPC_email = ""#"dliu@ig.utexas.edu" # Email to receive job status.

# Creating the fault interface
par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin,par.fxmax,par.nfx) # coordinates of fault grids along strike.
par.fz  = np.linspace(par.fzmin,par.fzmax,par.nfz) # coordinates of fault grids along dip.
print(par.nfx, par.nfz)
# Create on_fault_vars array for on_fault varialbes.
par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
# functions are defined in lib.py under scripts/
# function lists:
# - shear_steady_state
# - state_steady_state
# - B1, defined in TPV104 and TPV105
# - B2 and B3, defined in TPV105
par.fric_sw_fs = 0.575
par.fric_sw_fd = 0.450
par.fric_sw_D0 = 0.18
par.grav       = 9.8
par.fric_cohesion = 0.
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    downDipDistance = abs(zcoor)/sin(abs(par.dip)/180.*pi)
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.

    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,5]   = 0.5
    if downDipDistance<=8.0e3:
        par.on_fault_vars[iz,ix,4] = 0.001875e6*(8e3 - downDipDistance)

    par.on_fault_vars[iz,ix,7]   = -0.00424e6*downDipDistance # initial normal stress. Negative compressive.
    par.on_fault_vars[iz,ix,8]   = 0.0 # initial along strike shear
    # positive initial dip stress to move y+ wall up.
    # negative initial dip stress to move y+ wall down. 
    par.on_fault_vars[iz,ix,49]  = 0.00212e6*downDipDistance # initial along dip shear
    
# (x,z) coordinate pairs for on-fault stations (in km).
C = sin(par.dip/180.*pi)
par.st_coor_on_fault = [[0.,0.], [0., -1.*C], [0.0, -3.0*C], [0.0,-6.*C], [0.0,-12.*C], 
    [0.0, -18.*C], [0.0, -24.*C], [4., 0.], [4.,-3.*C], [4,-6*C], [4,-12*C], [4,-18*C],
    [8, 0], [8,-3*C], [8, -6*C], [8,-12*C], [8,-18*C],[12,0],[12,-3*C]]
print(par.st_coor_on_fault)
# (x,y,z) coordinates for off-fault stations (in km).
par.st_coor_off_fault = [[0,-9,0], [10,-9,0], [20,-9,0], [0,-3,0], 
    [0,-1,0], [10,-1,0], [0,1,0], [10,1,0], [0,3,0], [0,9,0],[10,9,0],
    [20,9,0], [0,15,0], [0,21,0], [0,27,0], [10,27,0], [20,27,0], [0,33,0],
    [0,39,0], [0,45,0], [10,45,0], [20,45,0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
