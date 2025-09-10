#! /usr/bin/env python3
from math import *
from sys  import *

# functions are defined in lib.py under scripts/
# function lists:
# - shear_steady_state
# - state_steady_state
# - B1, defined in TPV104 and TPV105
# - B2 and B3, defined in TPV105

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
