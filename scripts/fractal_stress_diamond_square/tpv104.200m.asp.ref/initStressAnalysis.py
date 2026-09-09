from user_defined_params import par
from math import *
import numpy as np
import sys
np.set_printoptions(precision=2)
def calcInitStress(dev2Vert, str1ToFaultAngle, strike, dip):
    str1ToFaultAngle = str1ToFaultAngle/180.*pi
    strike = strike/180.*pi
    dip = dip/180.*pi
    
    vertStress = -1.
    #dev2Vert   = par.devStrToStrVertRatio
    #str1ToFaultAngle = par.str1ToFaultAngle
    #dip = par.dip

    H1Stress    = vertStress * (1. + dev2Vert*cos(2.*str1ToFaultAngle))
    H3Stress    = vertStress * (1. - dev2Vert*cos(2.*str1ToFaultAngle))
    shear     = abs(vertStress) * dev2Vert * sin(2.*str1ToFaultAngle)

    stressTensor = np.zeros((3,3))
    stressTensor[0,0] = H1Stress
    stressTensor[1,1] = H3Stress
    stressTensor[2,2] = vertStress
    stressTensor[0,1] = shear
    stressTensor[1,0] = shear
    
    un = np.zeros(3)
    us = np.zeros(3)
    ud = np.zeros(3)
    ud1 = np.zeros(3)
    
    un[0] = cos(strike)*sin(dip)
    un[1] = -sin(strike)*sin(dip)
    un[2] = cos(dip)
    us[0] = -sin(strike)
    us[1] = -cos(strike)
    us[2] = 0.
    ud[0] = cos(strike)*cos(dip)
    ud[1] = sin(strike)*cos(dip)
    ud[2] = sin(dip)
    ud1[0] = us[1]*un[2] - us[2]*un[1]
    ud1[1] = us[2]*un[0] - us[0]*un[2]
    ud1[2] = us[0]*un[1] - us[1]*un[0]
    T = np.dot(stressTensor, un)
    nT = np.dot(T,un)
    nS = np.dot(T,us)
    nD = np.dot(T,ud)
    
    print('Stress state tensor is ')
    print(stressTensor) 
    print('un,us,ud')
    print(un, us, ud, ud1)
    print('T ', T)
    print('nT ', "%.2f" % nT)
    print('nS ', "%.2f" % nS)
    print('nD ', "%.2f" % nD)
    
    print(' ')
    print('nS/nT = ', "%.2f" % nS/nT)
    print('nD/nT = ', "%.2f" % nD/nT)
    print('nShear/nT = ', "%.2f" % (nD**2+nS**2)**0.5/nT )
    
dev2Vert = float(sys.argv[1])
str1ToFaultAngle = float(sys.argv[2])
strike = float(sys.argv[3])
dip = float(sys.argv[4])
print(sys.argv)

calcInitStress(dev2Vert, str1ToFaultAngle, strike, dip)

