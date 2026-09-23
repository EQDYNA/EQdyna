#! /usr/bin/env python3
"""
Regression guard (4465c17's retry precondition) for `rsfNucleation`'s
TPV==2802 branch (drv.a6, C_nuclea=1) losing `Td` from its call.

Background: refactor round 2 (aca6979) dropped `Td` from
`rsfNucleation`'s signature (faulting.py:333, was
`def rsfNucleation(xp, finv, fric, Tn, Ts, Td, srS, srD, timeElapsed, nt)`)
and from its caller in `solveRSF` (faulting.py:477) while the TPV==2802
branch's body (faulting.py:351-363, mirroring faulting.f90's `elseif
(TPV == 2802)` arm at faulting.f90:367-379) still evaluates
`xp.sqrt(Ts ** 2 + Td ** 2)` -- NameError on step 1 of every drv.a6 python
cell. Nothing caught it: drv.a6's python cells are not in CI_CELLS, and the
landing check ran a case without RSF nucleation. Commit 4465c17 reverted
aca6979 and recorded the retry precondition this file satisfies: "a
unit/regression test that calls rsfNucleation with TPV=2802".

WHY DIRECT CALL, NOT A CALLER: `rsfNucleation`'s signature (ten plain
positional arrays/scalars, no MPI/case/mesh state) is not awkward to
construct -- it is the narrowest thing that reaches the branch, so this
calls it directly rather than through `solveRSF` (which would need un/us/ud,
arn, mr, massSlave/massMaster force-index arrays, and comps['fx_s'] etc.
that TPV==2802's own defect has nothing to do with).

ORACLE (rule 23 -- Fortran is the reference): `_oracle_tpv2802` below is a
scalar, per-node re-derivation of faulting.f90:344-384 written with Python's
`math` module -- a different library, a different (explicit per-node
if/else, not `xp.where`) control-flow shape than faulting.py -- and it is
never given faulting.rsfNucleation as an input or compared against a
captured "looks right" run of it. It is compared AGAINST
faulting.rsfNucleation's actual output at runtime.

Cheap (rule 9): three fault nodes, no case/mesh/build, pure Python; both
numpy and jax backends run in well under 50 ms combined.
"""
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PYSRC = os.path.join(ROOT, 'src', 'python')
sys.path.insert(0, PYSRC)

# --- three synthetic drv.a6-shaped fault nodes --------------------------
# Node 0, 1: inside the nucleation patch (radius < nucR) -> F > 0.
# Node 2: OUTSIDE the patch (radius > nucR) -> F == 0, dtau == 0, so Ts is
# unchanged for that node regardless of nt -- a cheap check that the F==0
# branch (faulting.f90:355-357's implicit "no if-body ran") is not silently
# turned into a nonzero perturbation.
NUC_RADIUS = [500.0, 1500.0, 4000.0]
NUC_R = 3000.0
NUC_T = 1.0
NUC_DTAU0_CASE = 2.5e6
TIME_ELAPSED = 0.4          # <= NUC_T, so G takes its exponential branch
TPV = 2802

TN = [-8.0e7, -6.5e7, -9.0e7]
TS = [3.0e6, 2.0e6, 4.0e6]
TD = [1.5e6, 1.0e6, 2.0e6]
SRS0 = [1.0e-9, 1.5e-9, 0.8e-9]
SRD0 = [0.5e-9, 0.7e-9, 0.4e-9]
VINI_X = [2.0e-10, 2.0e-10, 2.0e-10]
VINI_Z = [1.0e-10, 1.0e-10, 1.0e-10]
RSF_A = [0.010, 0.012, 0.011]
RSF_V0 = [1.0e-6, 1.0e-6, 1.0e-6]

# Pre-existing fric slots for the nt != 1 case: deliberately far from
# whatever the nt==1 formula would produce, so a mutation that re-derives
# state/theta_pc/nuc_dtau0 on EVERY step (not just nt==1) is caught.
PREV_STATE = [7.77, 8.88, 9.99]
PREV_THETA_PC = [1.11e7, 2.22e7, 3.33e7]
PREV_NUC_DTAU0 = [9.9e6, 8.8e6, 7.7e6]


def _oracle_tpv2802(nt):
    """Independent re-derivation of faulting.f90:344-384's TPV==2802 arm,
    scalar per node, using `math` -- not numpy, not faulting.py. Returns
    (Ts_out, state_out, theta_pc_out, nuc_dtau0_out), each a length-3 list.
    """
    n = len(TN)
    ts_out = [0.0] * n
    state_out = list(PREV_STATE)
    theta_pc_out = list(PREV_THETA_PC)
    nuc_dtau0_out = list(PREV_NUC_DTAU0)
    for i in range(n):
        # faulting.f90:354-357
        if NUC_RADIUS[i] < NUC_R:
            F = math.exp(NUC_RADIUS[i] ** 2 / (NUC_RADIUS[i] ** 2 - NUC_R ** 2))
        else:
            F = 0.0
        # faulting.f90:358-360
        if TIME_ELAPSED <= NUC_T:
            G = math.exp((TIME_ELAPSED - NUC_T) ** 2
                         / (TIME_ELAPSED * (TIME_ELAPSED - 2.0 * NUC_T)))
        else:
            G = 1.0
        # faulting.f90:369-377 -- ONLY on the first step.
        if nt == 1:
            nuc_dtau0_out[i] = NUC_DTAU0_CASE
            ttao = math.sqrt(TS[i] ** 2 + TD[i] ** 2)
            backSliprate = math.sqrt((SRS0[i] + VINI_X[i]) ** 2
                                      + (SRD0[i] + VINI_Z[i]) ** 2)
            state_out[i] = RSF_A[i] * math.log(
                2.0 * RSF_V0[i] / backSliprate
                * math.sinh(ttao / abs(TN[i]) / RSF_A[i]))
            theta_pc_out[i] = abs(TN[i])
        # faulting.f90:379,381
        dtau = nuc_dtau0_out[i] * F * G
        ts_out[i] = TS[i] + dtau
    return ts_out, state_out, theta_pc_out, nuc_dtau0_out


def _build_fric(xp, np_mod):
    NCOLS = 81  # gv.NUC_DTAU0 (80) + 1
    fric = np_mod.zeros((3, NCOLS))
    from eqdyna import globalvar as gv
    for i in range(3):
        fric[i, gv.VINI_X] = VINI_X[i]
        fric[i, gv.VINI_Z] = VINI_Z[i]
        fric[i, gv.RSF_A] = RSF_A[i]
        fric[i, gv.RSF_V0] = RSF_V0[i]
        fric[i, gv.STATE] = PREV_STATE[i]
        fric[i, gv.THETA_PC] = PREV_THETA_PC[i]
        fric[i, gv.NUC_DTAU0] = PREV_NUC_DTAU0[i]
    return xp.asarray(fric)


def _to_list(xp, a):
    import numpy as np
    if hasattr(a, 'tolist'):
        return np.asarray(a).tolist()
    return list(a)


def _check_nt(backend_name, nt, label):
    import numpy as np
    from eqdyna import backend as B
    from eqdyna import faulting
    from eqdyna import globalvar as gv

    xp = B.array_module(backend_name)
    fric = _build_fric(xp, np)
    Tn = xp.asarray(TN); Ts = xp.asarray(TS); Td = xp.asarray(TD)
    srS = xp.asarray(SRS0); srD = xp.asarray(SRD0)
    finv = {'nuc_radius': xp.asarray(NUC_RADIUS), 'nucR': NUC_R, 'nucT': NUC_T,
            'TPV': TPV, 'nucdtau0': NUC_DTAU0_CASE}

    fric_out, Ts_out = faulting.rsfNucleation(
        xp, finv, fric, Tn, Ts, Td, srS, srD, TIME_ELAPSED, nt)

    want_ts, want_state, want_theta_pc, want_dtau0 = _oracle_tpv2802(nt)
    got_ts = _to_list(xp, Ts_out)
    got_state = _to_list(xp, fric_out[:, gv.STATE])
    got_theta_pc = _to_list(xp, fric_out[:, gv.THETA_PC])
    got_dtau0 = _to_list(xp, fric_out[:, gv.NUC_DTAU0])

    def cmp(name, got, want):
        for i in range(3):
            if not math.isclose(got[i], want[i], rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(
                    '%s/%s node %d: got %r, want %r (oracle: faulting.f90:'
                    '344-384 TPV==2802 arm)' % (backend_name, name, i, got[i], want[i]))

    cmp('%s Ts' % label, got_ts, want_ts)
    cmp('%s fric[:,STATE]' % label, got_state, want_state)
    cmp('%s fric[:,THETA_PC]' % label, got_theta_pc, want_theta_pc)
    cmp('%s fric[:,NUC_DTAU0]' % label, got_dtau0, want_dtau0)
    print('  %s %s: Ts/STATE/THETA_PC/NUC_DTAU0 match the faulting.f90:344-384 oracle'
          % (backend_name, label))


def main():
    print('Regression guard: rsfNucleation TPV==2802 (4465c17 retry precondition)')
    fails = []
    for backend_name in ('numpy', 'jax'):
        for nt, label in ((1, 'nt=1 (first step)'), (2, 'nt=2 (not first step)')):
            try:
                _check_nt(backend_name, nt, label)
            except Exception as e:  # noqa: BLE001 -- report every failure, don't stop at the first
                fails.append('%s %s: %s: %s' % (backend_name, label, type(e).__name__, e))

    if fails:
        print('FAIL test_rsfNucleation_tpv2802_td')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_rsfNucleation_tpv2802_td (numpy and jax, nt==1 and nt!=1, '
          'against an independent faulting.f90:344-384 oracle)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
