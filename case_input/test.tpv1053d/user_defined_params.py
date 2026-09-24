#! /usr/bin/env python3

from defaultParameters import parameters

par = parameters()

# Station n-stress sign (board row 22a): TPV105-3D's format spec
# (scratch/specs/TPV105_3D_formats, n-stress field, lines 114-115) says
# "Sign convention: Positive means compression." -- the opposite of the SCEC
# default in defaultParameters.py.
par.faultStNormalStressSign = 'compression'
par.term = 5.


