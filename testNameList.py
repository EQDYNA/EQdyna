#! /usr/bin/env python3

nameList = ['test.drv.a6',   'test.tpv8', 'test.tpv10', 'test.tpv104',
            'test.tpv1053d', 'test.meng2023a', 'test.meng2023cb',
            'test.tpv29', 'test.tpv36', 'test.tpv37', 'test.tpv30',
            'test.tpv22', 'test.tpv23', 'test.tpv35', 'test.tpv34']
coreNumList = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# test.tpv30 (rough fault + Drucker-Prager viscoplasticity) REGISTERED
# 2026-09-23 on the owner's gating decision, at the ONE 5 s term and dx=500 m,
# fortran + python-jax. The divergence that kept it out (numpy==jax, both !=
# fortran, by ~t=6 s) was root-caused and fixed in e1888e7 (a PML node never
# received its own elements' gravity); the 5 s gate catches that bug (with
# e1888e7 reverted: max|diff| 1.20e+08 vs bound 1e-10, session log KK).
#
# test.tpv22/test.tpv23 (SCEC TPV22/23 stepover benchmarks, two vertical
# planar faults each) REGISTERED (rule 17, item 17 section A, replacing
# test.multifault2 as the two-fault gate case -- see pathway_forward.md item
# 17): 4 ranks (par.nx,ny,nz=2,1,2), fortran + python-jax. These two run at
# their OWN gate term (15.0 s, testsys/matrix.py's CASE_TERM_OVERRIDE) rather
# than the everyday GATE_TERM_S (5.0 s) -- the SCEC spec
# (TPV22_23_Description_v08.pdf, Part 5) requires 15 s post-nucleation for
# fault #2 to rupture at all; see matrix.py's own comment for the measured
# wall-time cost of the exception.
