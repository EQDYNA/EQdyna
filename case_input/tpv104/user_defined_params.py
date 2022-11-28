#! /usr/bin/env python3

import numpy as np
from math import *

# model_domain (in meters)
xmin, xmax   = -42.0e3, 42.0e3
ymin, ymax   = -20.0e3, 22.0e3
zmin, zmax   = -42.0e3, 0.0e3

# fault geometry (in meters)
fxmin, fxmax = -18.0e3, 18.0e3
fymin, fymax = 0.0e3,   0.0e3    # for vertical strike-slip faults, we align faults along xz planes.
fzmin, fzmax = -18.0e3, 0.0e3 

dx           = 100.0e0 # cell size, spatial resolution
nuni_y_plus  = 70 # along the fault-normal dimension, the number of cells share the dx cell size.
nuni_y_minus = 70 
enlarging_ratio = 1.025e0 # along the fault-normal dimension (y), cell size will be enlarged at this ratio compoundly.

# Isotropic material propterty.
# Vp, Vs, Rou
vp, vs, rou = 6.0e3, 3.464e3, 2.67e3
init_norm = -25.0e6 # initial normal stress in Pa. Negative compressive.

# total simulation time and dt
term        = 15.
dt          = 0.008

# Controlling switches for EQquasi system
C_elastic   = 1 # elastic(1).
C_nuclea    = 1 # artificial nucleation (1), no (0). 
C_degen     = 0 # degenerate hexahedrals (1), no (0).
friclaw     = 3 # sw(1), tw(2), rsf_aging(3), rsf_slip(4), rsf_slip_srw(5).
ntotft      = 1 # number of total faults.
nucfault    = 1 # the fault id of nucleation fault. Should be no larger than ntotft
rough_fault = 0 # include rough fault yes(1) or not(0).
nt_out      = 20 # Every nt_out time steps, disp of the whole model and on-fault variables will be written out in netCDF format.
tpv          = 104 
# currently supported cases
# 104  (SCEC-TPV104)
# 105  (SCEC-TPV1053D)
# 2801 (DRV)
# 1001 (GM-cycle)

#################################
##### Frictional variables ######
#################################
# friclaw == 1, slip weakening
fric_sw_fs      = 0.18
fric_sw_fd      = 0.12
fric_sw_D0      = 0.3
# parameters needed for rsf.
# friclaw == 3/4, rsf with the aging law/slip law.
fric_rsf_a      = 0.01 
fric_rsf_b      = 0.014
fric_rsf_Dc     = 0.4
fric_rsf_deltaa = 0.01
fric_rsf_r0     = 0.6
fric_rsf_v0     = 1e-6
# additional parameters for friclaw == 5, rsf with the slip law and strong rate weakening.
fric_rsf_fw     = 0.2
fric_rsf_vw     = 0.1
fric_rsf_deltaavw0 = 0.9
# additional parameters for thermal pressurization. 
fric_tp_a_th    = 1.0e-6  # m^2/s
fric_tp_rouc    = 2.7d6   # J/(m^3K)
fric_tp_lambbda = 0.1d6   # paK^-1
fric_tp_h       = 0.02d0  # m
fric_tp_a_hy    = 4.0d-4  # m^2/s
fric_tp_deltaa_hy0 = 1.0d0# m^2/s
fric_tp_Tini    = 483.15  # K
fric_tp_pini    = 80.0d6  # Pa

#################################
#####   Initial stresses   ######
#################################

#init_norm       = -120.0e6

# Creating the fault interface
nfx     = int((fxmax - fxmin)/dx + 1)
nfz     = int((fzmax - fzmin)/dx + 1)
fx      = np.linspace(xmin,xmax,nfx) # coordinates of fault grids along strike.
fz      = np.linspace(zmin,zmax,nfz) # coordinates of fault grids along dip.

# Create on_fault_vars array for on_fault varialbes.
on_fault_vars = np.zeros((nfx,nfz,100))
def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) #+ rou*vs/2.0*slip_rate
  return res

def state_steady_state(a,b,d0,v0,r0,shear,norm,slip_rate):
  # calculate the state variable at steady state
	if friclaw == 3:
		tmp   = a*np.log(2.*sinh(abs(shear/norm/a))) - r0 - dlog(slip_rate/v0)
		state = d0/v0*np.exp(tmp/b)
	elif friclaw == 4 of friclaw == 5:
		state = a*np.log(2.*v0/slip_rate*sinh(abs(shear/norm/a)))
  return state
  
for ix, xcoor in enumerate(fx):
  for iz, zcoor in enumerate(fz):
  # assign a in RSF. a is a 2D distribution.
    on_fault_vars[ix,iz,1]   = fric_sw_fs 
    on_fault_vars[ix,iz,2]   = fric_sw_fd
    on_fault_vars[ix,iz,3]   = fric_sw_D0
    on_fault_vars[ix,iz,7]   = -120.0e6     # initial normal stress. Negative compressive.
    on_fault_vars[ix,iz,8]   = 40.0e6       # initial shear stress.
	
    if abs(xcoor)<=15e3:
		tmp1  = 1.0
	elif abs(xcoor)<18e3 .and. abs(xcoor)>15e3:
		tmp1  = 0.5*(1. + np.tanh(3e3/(abs(xcoor) - 18e3) + 3e3/(abs(xcoor) - 15e3)))
	else:
		tmp1 = 0.0
	if abs(zcoor - -7.5e3)<=7.5e3:
		tmp2 = 1.0
	elif abs(zcoor - -7.5e3)<10.5e3 and abs(zcoor - -7.5e3)>7.5e3:
		tmp2  = 0.5*(1. + np.tanh(3e3/(abs(zcoor -- 7.5e3) - 10.5e3) + 3e3/(abs(zcoor - -7.5e3) - 7.5e3)))
	else: 
		tmp2  = 0.0

    on_fault_vars[ix,iz,9]  = fric_rsf_a + (1. - tmp1*tmp2)*fric_rsf_deltaa
    on_fault_vars[ix,iz,10] = fric_rsf_b # assign b in RSF 
    on_fault_vars[ix,iz,11] = fric_rsf_Dc # assign Dc in RSF.
    #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    #  on_fault_vars[ix,iz,11] = minDc # a special Dc zone.
    on_fault_vars[ix,iz,12] = fric_rsf_v0 # initial reference slip rate.
    on_fault_vars[ix,iz,13] = fric_rsf_r0 # initial reference friction.
    
    on_fault_vars[ix,iz,14] = fric_rsf_fw # 
    on_fault_vars[ix,iz,15] = fric_rsf_vw  + fric_rsf_deltaavw0*(1. - tmp1*tmp2)  # 
	on_fault_vars[ix,iz,16] = fric_tp_a_hy + fric_tp_deltaa_hy0*(1. - tmp1*tmp2)  #
	
	on_fault_vars[ix,iz,46] = creep_slip_rate # initial slip rates
    #if (xcoor<=-18e3 and xcoor>=-30e3 and zcoor<=-4e3 and zcoor>=-16e3):
    #  on_fault_vars[ix,iz,46] = 0.03 # initial high slip rate patch.
	
    on_fault_vars[ix,iz,20] = state_steady_state(on_fault_vars[ix,iz,9], 
												on_fault_vars[ix,iz,10],
												on_fault_vars[ix,iz,11],
												on_fault_vars[ix,iz,12],
												on_fault_vars[ix,iz,13],
												on_fault_vars[ix,iz,8],
												on_fault_vars[ix,iz,7],
												on_fault_vars[ix,iz,46]) # initial state var.
	
    
###############################################
##### Domain boundaries for transferring ######
###############################################
xmin_trans, xmax_trans = -25e3, 25e3
zmin_trans = -25e3
ymin_trans, ymax_trans = -5e3, 5e3
dx_trans = 50 

####################################
##### HPC resource allocation ######
####################################
casename = str(tpv)
nx = 4
ny = 5
nz = 2

HPC_ncpu  = nx*ny*nz # Number of CPUs requested.
HPC_nnode = int(HPC_ncpu, 128) + 1 # Number of computing nodes. On LS6, one node has 128 CPUs.
HPC_queue = "normal" # q status. Depending on systems, job WALLTIME and Node requested.
HPC_time  = "02:00:00" # WALLTIME, in hh:mm:ss format.
HPC_account = "EAR22013" # Project account to be charged SUs against.
HPC_email = ""#"dliu@ig.utexas.edu" # Email to receive job status.

##############################################
##### Single station time series output ######
##############################################

# (x,z) coordinate pairs for on-fault stations (in km).
st_coor_on_fault = [[-36.0, 0.0], [-16.0,0.0], [0.0,0.0], [16.0,0.0], \
   [36.0,0.0], [-24.0,0.0], [-16.0,0.0], [0.0,-10.0], [16.0,-10.0], [0.0,-22.0]]
   
# (x,y,z) coordinates for off-fault stations (in km).
st_coor_off_fault = [[0,8,0], [0,8,-10], [0,16,0], [0,16,-10], [0,32,0], \
   [0,32,-10], [0,48,0], [16,8,0], [-16,8,0]]
n_on_fault  = len(st_coor_on_fault)
n_off_fault = len(st_coor_off_fault)



