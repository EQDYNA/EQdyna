#! /usr/bin/env python3

nameList = ['test.drv.a6',   'test.tpv8', 'test.tpv10', 'test.tpv104',
            'test.tpv1053d', 'test.meng2023a', 'test.meng2023cb',
            'test.tpv29', 'test.tpv36', 'test.tpv37', 'test.tpv30']
coreNumList = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# test.tpv30 (rough fault + Drucker-Prager viscoplasticity) REGISTERED
# 2026-09-23 on the owner's gating decision, at the ONE 5 s term and dx=500 m,
# fortran + python-jax. The divergence that kept it out (numpy==jax, both !=
# fortran, by ~t=6 s) was root-caused and fixed in e1888e7 (a PML node never
# received its own elements' gravity); the 5 s gate catches that bug (with
# e1888e7 reverted: max|diff| 1.20e+08 vs bound 1e-10, session log KK).
