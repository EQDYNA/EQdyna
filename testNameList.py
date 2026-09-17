#! /usr/bin/env python3

nameList = ['test.drv.a6',   'test.tpv8', 'test.tpv10', 'test.tpv104',
            'test.tpv1053d', 'test.meng2023a', 'test.meng2023cb',
            'test.tpv29', 'test.tpv36', 'test.tpv37']
coreNumList = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# test.tpv30 (rough fault + Drucker-Prager viscoplasticity) is NOT registered
# here yet -- see case_input/test.tpv30/README.md "Gate status" and
# NOTES_tpv30_gate.md. The compset, geometry, and the swtwNucleation TPV==30
# branch (src/fortran/faulting.f90, src/python/eqdyna/faulting.py) are landed
# and verified correct at short duration (bit-exact vs Fortran through 24
# steps); a real, deterministic (numpy==jax, both != fortran) divergence
# appears by ~144 steps (t=6s of 20s) that a strict abs-max bound cannot
# pass and that is NOT yet root-caused. Registering the case here before that
# is closed would freeze a reference nobody has explained (rule 17 step 7:
# a new case is not added with a backend that fails).
