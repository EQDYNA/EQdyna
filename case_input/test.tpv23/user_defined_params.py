#! /usr/bin/env python3
"""
SCEC TPV23: left-step on right-lateral vertical strike-slip faults,
COMPRESSIONAL step, 1.0 km stepover distance (TPV22_23_Description_v08.pdf,
Part 2, p.4; strike.scec.org/cvws/tpv22_23docs.html). Linear elastic.

250 m isotropic resolution (see tpv22_23_common.py for every parameter's
provenance and the resolution-scan evidence this value is picked from --
NOTES_tpv2223_iteration.md iteration 5). Each fault meshed on its own TRUE
box (meshgen.f90's getLocalOneDimCoorArrAndSize unions every fault's own
box; no shared-box/scaffold workaround). NOT the spec's 100/50 m resolution
-- full_specs.py carries that for a later, scheduled run.
"""
import tpv22_23_common as common

par = common.build_params(tpv=23, fault2_z=1000.0)
