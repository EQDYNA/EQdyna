#! /usr/bin/env python3
"""
SCEC TPV22: right-step on right-lateral vertical strike-slip faults,
EXTENSIONAL step, 1.6 km stepover distance (TPV22_23_Description_v08.pdf,
Part 1, p.3; strike.scec.org/cvws/tpv22_23docs.html). Linear elastic.

Gate-coarse case (see tpv22_23_common.py for every parameter's provenance
and the mesh-sharing geometry fix this case needs). NOT the spec's 100/50 m
resolution -- full_specs.py carries that for a later, scheduled run.
"""
import tpv22_23_common as common

par = common.build_params(tpv=22, fault2_z=-1600.0)
