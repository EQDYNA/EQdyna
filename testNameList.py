#! /usr/bin/env python3

nameList = ['test.drv.a6',   'test.tpv8', 'test.tpv10', 'test.tpv104',
            'test.tpv1053d', 'test.meng2023a', 'test.meng2023cb',
            'test.tpv29', 'test.tpv36', 'test.tpv37', 'test.tpv30',
            'test.tpv22', 'test.tpv23', 'test.tpv35', 'test.tpv34',
            'test.tpv26', 'test.tpv27', 'test.tpv31', 'test.tpv32',
            'test.tpv33', 'test.tpv12']
coreNumList = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
# test.tpv12 (SCEC TPV12, 60-degree dipping planar normal fault, linear
# elastic) REGISTERED 2026-10-07 (board row 149, owner pairing decision):
# TPV12 is the METHOD 1 (C_elastic=1, stress change, no gravity) half of the
# TPV12/13 pair, built exactly as test.tpv29 is to test.tpv30's METHOD 2. 4
# ranks (par.nx,ny,nz=2,2,1, test.tpv10's own partition), fortran +
# python-jax, at the ONE 5 s GATE_TERM_S and dx=500 m. Confirmed (direct read
# of TPV12_13_Description_v6.pdf) NOT blocked: nucleation is reduced static
# friction in a patch (C_nuclea=0, test.tpv10's own mechanism, no forced-
# rupture code path at all), geometry is test.tpv10's already-gated
# insertFaultType=1 planar-dipping-fault + mod4dip machinery verbatim, and
# initial stress is per-node manual fault tractions (Method 1, no off-fault
# stress tensor, so setPlasticStress is never called for this case). TPV13
# (the Method 2 / Drucker-Prager sibling) is NOT added here -- see
# case_input/test.tpv12/README.md and pathway_forward.md row 149 for the
# still-open gap that blocks it.
# test.tpv33 (SCEC TPV33, planar vertical strike-slip, fault-parallel
# low-velocity zone) REGISTERED 2026-10-07 (board row 152): 4 ranks
# (par.nx,ny,nz=2,1,2), fortran + python-jax, at the ONE 5 s GATE_TERM_S and
# dx=400 m. Rule 17 step 3 (the question the board row posed): confirmed NOT
# blocked -- the EXISTING n2mat==6 3D nearest-cell material grid (TPV34, PR
# #93/#114) already covers a field that varies along one axis only (here,
# EQdyna y = the spec's fault-normal coordinate), no new material-grid
# architecture needed. See case_input/test.tpv33/README.md for the full
# evidence (file:line citations) and case_input/test.tpv33/user_defined_params.py
# for the grid construction. Nucleation is a static, time-independent
# additive shear-stress bump (swtwNucleation is a no-op for unmatched
# par.tpv, same as test.tpv8/test.tpv31/test.tpv34); initial stress is
# SCEC's own Method 2 (on-fault traction only), same mechanism as
# test.tpv26/29/31. Zero src/ changes.
# test.tpv31/test.tpv32 (SCEC TPV31/32, planar vertical strike-slip, 1D
# layered velocity structure -- discontinuous/continuous respectively)
# REGISTERED 2026-10-06 (board row 151): 4 ranks (par.nx,ny,nz=2,1,2),
# fortran + python-jax, at the ONE 5 s GATE_TERM_S and dx=500 m. Confirmed
# (TPV31_32_Description_v03, direct spec read) NOT blocked: nucleation is a
# static, time-independent additive shear-stress bump (swtwNucleation is a
# no-op for unmatched par.tpv, same as test.tpv8/test.tpv10), initial stress
# is SCEC's own Method 2 (on-fault traction only, no off-fault tensor) same
# mechanism as test.tpv26/29, and the 1D layered material reuses the
# EXISTING n2mat==4 mechanism (test.meng2023a) by precomputing one row per
# mesh z-layer from the spec's table. Zero src/ changes.
# test.tpv26/test.tpv27 (SCEC TPV26/27, planar vertical strike-slip, TPV27
# adding Drucker-Prager viscoplasticity) REGISTERED 2026-10-06 (board row
# 150): 4 ranks (par.nx,ny,nz=2,1,2), fortran + python-jax, at the ONE 5 s
# GATE_TERM_S and dx=500 m. Confirmed (TPV26_27_Description_v13, direct spec
# read) these use ORDINARY smoothed forced-rupture nucleation (Part 5,
# identical formula to TPV22/23/29/30/36/37/201) and SCEC Method 1/Method 2
# initial-stress machinery EQdyna already has from TPV29/30 -- NOT the
# row 149/TPV12-13 gravity-everywhere-with-C_elastic=1 gap.
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
