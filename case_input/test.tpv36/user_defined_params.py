#! /usr/bin/env python3
#
# test.tpv36 and test.tpv37 are identical except for par.tpv and the
# cohesion-buildup slope; both build via tpv36_37_common.buildParams().
from tpv36_37_common import buildParams

par = buildParams(tpv=36, cohesionSlope=0.0005e6)
