#! /usr/bin/env python3
from math import *
import sys
import glob
import re
from os.path import exists

import numpy as np

# functions are defined in lib.py under scripts/
# function lists:
# - shear_steady_state
# - state_steady_state
# - B1, defined in TPV104 and TPV105
# - B2 and B3, defined in TPV105
# - loadFrtData, shared frt.txt* loader for plotRuptureDynamics/plotSlipAndRPT
# - tryint, alphanum_key, sort_nicely, generate_gif, seek_numbers_filename,
#   shared filename-sorting/gif helpers for plot_on_fault_vars

def shear_steady_state(a,b,v0,r0,load_rate,norm,slip_rate):
  # calculate shear stress at steady state
  res = -norm*a*asinh(slip_rate/2.0/v0*exp((r0+b*log(v0/load_rate))/a)) #+ rou*vs/2.0*slip_rate
  return res

def state_steady_state(a,b,d0,v0,r0,shear,norm,slip_rate, friclaw):
  # calculate the state variable at steady state
    if friclaw == 3:
        tmp   = a*log(2.*sinh(abs(shear/norm/a))) - r0 - log(slip_rate/v0)
        state = d0/v0*exp(tmp/b)
    elif friclaw == 4 or friclaw == 5:
        state = a*log(2.*v0/slip_rate*sinh(abs(shear/norm/a)))
    return state

# Mathematically smooth version of the boxcar function B1, B2, and B3
def B1(x,ww,w):
  if abs(x)<=ww:
    res = 1.0
  elif abs(x)>ww and abs(x)<ww+w: 
    res = 0.5*(1. + tanh(w/(abs(x)-ww-w) + w/(abs(x)-ww)))
  elif abs(x)>=ww+w:
    res = 0.0
  return res

def B2(y,ww,w):
  if y<0:
    print('z coordinates should be positive for B3')
    sys.exit()
  if y<w:
    res = 0.5*(1.+tanh(w/(w-y) - w/(y+1e-30)))
  elif y>=w and y<=ww:
    res = 1.0
  elif y>ww and y<ww+w: 
    res = 0.5*(1. + tanh(w/(y-ww-w) + w/(y-ww)))
  elif y>=ww+w:
    res = 0.0
  return res
def B3(y,ww,w):
  if y<0:
    print('z coordinates should be positive for B3')
    sys.exit()
  if y<=ww:
    res = 1.
  elif y>ww and y<ww+w:
    res = 0.5*(1. + tanh(w/(y-ww-w) + w/(y-ww)))
  elif y>=ww+w:
    res = 0.
  return res

def linear1(x,ww,w):
  if abs(x)<=ww:
    res = 1.0
  elif abs(x)>ww and abs(x)<ww+w:
    res = 1. - abs((abs(x)-ww))/w
  elif abs(x)>=ww+w:
    res = 0.0
  return res

def loadFrtData(par):
    """Load and grid the on-fault frt.txt* output written by EQdyna.

    Shared by plotRuptureDynamics and plotSlipAndRPT: loops over
    frt.txt{me} for me in range(nx*ny*nz), grids each row onto the
    (dip, strike) fault mesh, and accumulates seismic moment.

    frt.txt* file structure (1-indexed as in the on-fault output; see
    the FRIC_SLOT_* slot map in src/globalvar.f90 for the underlying
    per-node fault-variable layout):
      1-3,   coorx,y,z
      4,     rupture time
      5-9,   final slips,d,n, final sliprates,d.
      10,    peak slip rate.
      11,    final sliprate magnitude
      12-14, final tnrm, tstk, tdip
      15-20, vxm,vym,vzm, vxs,vys,vzs.
      21,    state variable
      22,    state var for normal stress variation (Shi and Day)

    Returns (xx, zz, rupt, rupt2d, fVarArr, magnitude):
      xx, zz    -- meshgrid of along-strike/along-dip coordinates, km
      rupt      -- (na*ma, 3) flat [xcoor, along-dip distance, rupture time]
      rupt2d    -- (ma, na, 100) gridded rupture-time/slip/stress panels
      fVarArr   -- (ma, na, 100) gridded fault-restart variables, passed
                   to generateNcRestart(faultVarArr)
      magnitude -- moment magnitude computed from summed slip*area*shearMod
    """
    nprocs = par.nx*par.ny*par.nz
    na     = round((par.fxmax-par.fxmin)/par.dx+1)
    ma     = round((par.fzmax-par.fzmin)/par.dz+1)
    rupt   = np.zeros((na*ma,3))
    rupt2d = np.zeros((ma,na,100))
    fVarArr= np.zeros((ma,na,100))

    [xx,zz] = np.meshgrid(par.fx,par.fz/sin(par.dip/180.*pi))
    xx = xx/1.e3
    zz = zz/1.e3#/sin(par.dip/180.*pi) # along dip distance
    moment = 0.

    for me in range(nprocs):
      fname = 'frt.txt' + str(me)
      if exists(fname):
        print('Post-processing ' + fname + ' ... ...')
        a = np.loadtxt(fname)
        n, m = a.shape
        for i in range(n):
            #!! use round() instead of int()!!
            ii = round((a[i,0] - par.fxmin)/par.dx)
            jj = round((a[i,2] - par.fzmin)/par.dz)

            rupt[jj*na+ii,0] = a[i,0]  # xcoor
            rupt[jj*na+ii,1] = -a[i,2]/sin(par.dip/180.*pi) # zcoor to along dip distance, reverse sign to positive numbers.
            rupt[jj*na+ii,2] = a[i,3]  # rupture time

            rupt2d[jj,ii,0]  = a[i,3]  # rupture time
            rupt2d[jj,ii,1]  = (a[i,4]**2 + a[i,5]**2)**0.5  # slip magnitude
            rupt2d[jj,ii,2]  = a[i,9]                        # peak slip rate
            rupt2d[jj,ii,3]  = a[i,10]                       # final slip rate
            shearMod = 3464**2*2800
            moment = moment + rupt2d[jj,ii,1]*par.dx*par.dx*shearMod
            rupt2d[jj,ii,4]  = a[i,12]/1.e6 # final shear stress
            rupt2d[jj,ii,5]  = a[i,11]/1.e6 # final normal stress
            rupt2d[jj,ii,6]  = a[i,13]/1.e6 # final dip shear
            rupt2d[jj,ii,7]  = a[i,4] # final slip s
            rupt2d[jj,ii,8]  = a[i,5] # final slip d

            #
            # fVarArr will be passed to the function generateNcRestart(faultVarArr):
            fVarArr[jj,ii,0]  = a[i,12] # shear_strike, Pa
            fVarArr[jj,ii,1]  = a[i,13] # shear_dip, Pa
            fVarArr[jj,ii,2]  = a[i,11] # effective_normal, Pa
            fVarArr[jj,ii,3]  = a[i,10] # slip_rate, m/s
            fVarArr[jj,ii,4]  = a[i,20] # state_variable
            fVarArr[jj,ii,5]  = a[i,21] # state_normal
            fVarArr[jj,ii,6]  = a[i,14] # vxm, m/s
            fVarArr[jj,ii,7]  = a[i,15] # vym
            fVarArr[jj,ii,8]  = a[i,16] # vzm
            fVarArr[jj,ii,9]  = a[i,17] # vxs
            fVarArr[jj,ii,10] = a[i,18] # vys
            fVarArr[jj,ii,11] = a[i,19] # vzs

    magnitude = 2/3*log10(moment*1.e7)-10.7

    return xx, zz, rupt, rupt2d, fVarArr, magnitude

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    return [tryint(c) for c in re.split('([0-9]+)', s)]

def sort_nicely(l):
    l.sort(key=alphanum_key, reverse=False)
    return l

def generate_gif():
  import imageio  # lazy: optional dependency, only needed for GIF generation
  filenames = glob.glob('.//*.png')
  filenames = sort_nicely(filenames)
  with imageio.get_writer('./on_fault_vars.gif', mode='I') as writer:
      for filename in filenames:
          image = imageio.v2.imread(filename)
          writer.append_data(image)

def seek_numbers_filename(x):
    return (x[6:10])
