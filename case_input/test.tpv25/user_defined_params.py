#! /usr/bin/env python3
"""
SCEC TPV25: left-lateral main fault with a 30 deg RESTRAINING branch fault
(TPV24_25_Description_v07.pdf; see tpv24_25_common.py for the full
derivation -- geometry, stress resolution onto both faults, friction,
nucleation).

Gate-coarse dx=1000 m (NOT the spec's 100/50 m submission resolution --
full_specs.py carries those tiers, unrun). Branch length is decimated to
10.0 km (true 12.0 km) to land on a dx=1000 m mesh -- see common module's
BRANCH-LENGTH APPROXIMATION note.

b22/b33/b23 (Initial Stress Tensor Coefficients table, p.6/7) are the ONLY
difference from TPV24.

STATUS: WIP, not yet gated -- do not freeze a reference from this file
without a passing fortran + python-jax gate run first (rule 17 step 7).
"""
import tpv24_25_common as common

par = common.build_params(tpv=25, b22=1.119338, b33=0.880661, b23=0.138704)
