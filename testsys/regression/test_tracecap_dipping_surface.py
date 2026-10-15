#! /usr/bin/env python3
"""
Regression guard (item 149, PR #212): the free-surface trace-corner cap on
tensile normal traction in solveSWTW (src/fortran/faulting.f90:199-215,
src/python/eqdyna/faulting.py's `on_trace` mask in build() + the cap in
solveSWTW). An audit of PR #212 found NO test exercised it (rule 14a: a
check must be able to come out both ways).

THE PHYSICS (see faulting.f90's comment at the cap, and TPV12_13_Description
p.10): at the free-surface trace of a non-vertical fault the traction-free
surface forces t_z = 0, i.e. Tn*un_z + Td*ud_z = 0, so a tensile trial Tn is
capped at |Td*ud_z/un_z|. Vertical faults (un_z=0) and buried nodes are
untouched.

Calls `faulting.build` (to pin the `on_trace` mask) and `faulting.solveSWTW`
(numpy backend) directly on a tiny hand-built 4-node state -- no case, no
binary, no MPI. Runs in well under 50 ms.

Four fault-node-pairs, each a (node id, cohesion-huge so the shear-strength
cap never fires, Ts=0) slip-weakening (friclaw=1) node with trial Td=-2e6 Pa
constant across every row:

  0: 60-degree-dipping fault, ON the surface trace (z=0), un_z=0.5,
     ud_z=0.8660254 (sin/cos of a 60-deg dip), TENSILE trial Tn=+5e6 Pa.
     Expected capped Tn = |Td*ud_z/un_z| = |-2e6*0.8660254/0.5| =
     3464101.6 Pa exactly (hand-computed) -- strictly less than the +5e6
     trial, so the cap must have fired.
  1: same dipping geometry as 0, COMPRESSIVE trial Tn=-5e6 Pa. Tn>0 is
     false, so the cap must NOT fire: output == trial, unchanged.
  2: same dipping geometry as 0, BURIED (z=-100 m, below -1e-3), tensile
     trial Tn=+5e6 Pa. Cap must NOT fire despite the tensile trial exceeding
     what the cap would otherwise impose: output == trial, unchanged.
  3: VERTICAL fault (un_z=0.0) at the surface (z=0), tensile trial
     Tn=+5e6 Pa. Cap must NOT fire (and must not divide by zero): output ==
     trial, unchanged and finite.

Together these four values are sensitive to every failure mode named in the
audit: deleting the cap (row 0 would read +5e6 instead of 3464101.6),
flipping its sign or dropping the abs() (row 0's capped value would come out
negative or with the wrong magsnitude), and firing where it must not (rows
2/3 would no longer equal their unmutated trial value).

A second, independent check: `assert_fortran_predicate_parity` greps
src/fortran/faulting.f90 for the three literal predicate terms the python
`on_trace` mask + cap must mirror (`> 0.0d0`, the `-1.0d-3` surface-trace
z-bound, the `1.0d-6` un_z floor) plus the `min(...abs(...` cap expression,
so a predicate edited out of ONE side alone is caught without running
either solver.

Exits non-zero on any failure.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

import numpy as np  # noqa: E402

from eqdyna import faulting, globalvar as gv  # noqa: E402

XP = np  # numpy backend, per the brief

DIP_UNZ = 0.5                 # cos(60 deg)
DIP_UDZ = 0.8660254037844387   # sin(60 deg)
TRIAL_TD = -2.0e6
EXPECT_CAP = abs(TRIAL_TD * DIP_UDZ / DIP_UNZ)  # 3464101.6151377548 Pa


def _build_state():
    """4 fault-node pairs: 0/1 share a 60-deg-dipping surface-trace geometry
    (tensile/compressive), 2 is the same geometry buried, 3 is a vertical
    fault at the surface. 8 mesh nodes (4 slave + 4 master), eq ids 0..23."""
    n = 4
    nsmp1 = np.arange(n)
    nsmp2 = np.arange(n) + n
    meshCoor = np.zeros((2 * n, 3))
    meshCoor[nsmp1, 2] = [0.0, 0.0, -100.0, 0.0]   # z of each slave node

    un = np.zeros((n, 3)); ud = np.zeros((n, 3)); us = np.zeros((n, 3))
    un[0] = un[1] = un[2] = [0.0, np.sqrt(1 - DIP_UNZ ** 2), DIP_UNZ]
    ud[0] = ud[1] = ud[2] = [0.0, -0.5, DIP_UDZ]
    un[3] = [1.0, 0.0, 0.0]          # vertical fault: un_z == 0
    ud[3] = [0.0, 0.0, 1.0]
    us[:] = [1.0, 0.0, 0.0]
    us[3] = [0.0, 1.0, 0.0]
    arn = np.full(n, 100.0)

    eq_ids = np.arange(2 * n * 3).reshape(2 * n, 3)
    fnms = np.full(2 * n, 1.0e6)

    S = dict(eq_ids=eq_ids, nsmp1=nsmp1, nsmp2=nsmp2, meshCoor=meshCoor,
             xsource=0.0, ysource=0.0, zsource=0.0, fnms=fnms,
             un=un, us=us, ud=ud, arn=arn, nftnd=n,
             friclaw=1, TPV=8, C_elastic=1.0, slipRateThres=1.0e-3)
    return S, eq_ids


def _fric_comps(n):
    """fric columns touched by solveSWTW's friclaw=1 path; COHESION huge so
    the shear-strength scale never engages (isolates the Tn cap). comps'
    Tn/Ts/Td are the TRIAL tractions solveSWTW receives, bypassing
    getNsdSlipSliprateTraction (not under test here)."""
    fric = np.zeros((n, gv.NUC_DTAU0 + 1))
    fric[:, gv.SW_FS] = 0.6
    fric[:, gv.SW_FD] = 0.5
    fric[:, gv.SW_D0] = 0.4       # nonzero: slip_weak divides slip/D0
    fric[:, gv.COHESION] = 1.0e12
    fric[:, gv.CUM_SLIP] = 0.0    # |slip| < 1e-10 -> slip_weak returns fs
    Tn = np.array([5.0e6, -5.0e6, 5.0e6, 5.0e6])
    Ts = np.zeros(n)
    Td = np.full(n, TRIAL_TD)
    comps = dict(Tn=Tn, Ts=Ts, Td=Td, srMag=np.zeros(n))
    return fric, comps


def check_cap_fires_and_abstains_correctly():
    S, eq_ids = _build_state()
    finv = faulting.build(S)
    want_trace = np.array([True, True, False, False])
    if not np.array_equal(finv['on_trace'], want_trace):
        raise AssertionError('build() on_trace mask = %r, want %r (surface-dipping, '
                             'surface-dipping, buried, vertical-surface)'
                             % (finv['on_trace'].tolist(), want_trace.tolist()))

    n = 4
    fric, comps = _fric_comps(n)
    force = np.zeros(2 * n * 3)
    fnft = np.full(n, 99999.0)
    fric_out, force_out, srMag = faulting.solveSWTW(
        XP, finv, fric, fnft, comps, force, timeElapsed=1.0, tr=None, friclaw=1)
    got = fric_out[:, gv.TRACT_NORM]

    want = np.array([EXPECT_CAP, -5.0e6, 5.0e6, 5.0e6])
    labels = ['(a) dipping surface tensile -- capped',
              '(b) dipping surface compressive -- unchanged',
              '(c) buried tensile -- unchanged',
              '(d) vertical-fault surface tensile -- unchanged, finite']
    for i, label in enumerate(labels):
        if not np.isfinite(got[i]):
            raise AssertionError('row %d %s: Tn = %r is not finite' % (i, label, got[i]))
        if abs(got[i] - want[i]) > 1e-6 * max(1.0, abs(want[i])):
            raise AssertionError('row %d %s: Tn = %.6f Pa, want %.6f Pa'
                                 % (i, label, got[i], want[i]))
    if not (got[0] < 5.0e6 - 1.0):
        raise AssertionError('row 0: capped Tn %.3f is not strictly below the +5e6 '
                             'trial -- the cap did not fire' % got[0])
    print('  PASS  on_trace mask == %r' % finv['on_trace'].tolist())
    for i, label in enumerate(labels):
        print('  PASS  %s: Tn = %.6f Pa' % (label, got[i]))


FORTRAN_PREDICATE_TERMS = [
    'nsdTractionVector(1) > 0.0d0',
    '> -1.0d-3',
    '> 1.0d-6',
    'min(nsdTractionVector(1), &',
    'abs(nsdTractionVector(3)*ud(3,iFaultNodePair,iFault)/un(3,iFaultNodePair,iFault))',
]


def check_fortran_predicate_parity():
    """Source-level parity: the Fortran block must still carry the same
    predicate (Tn>0, surface z>-1e-3, |un_z|>1e-6) and the same
    min(Tn, |Td*ud_z/un_z|) expression the python mask/cap mirrors. Grep-
    based on purpose: this is a guard against the two sides drifting apart,
    not a re-run of either solver."""
    path = os.path.join(ROOT, 'src', 'fortran', 'faulting.f90')
    with open(path) as f:
        text = f.read()
    missing = [t for t in FORTRAN_PREDICATE_TERMS if t not in text]
    if missing:
        raise AssertionError('src/fortran/faulting.f90 is missing %d trace-cap '
                             'predicate term(s) the python side mirrors: %r'
                             % (len(missing), missing))
    print('  PASS  src/fortran/faulting.f90 carries all %d trace-cap predicate terms'
          % len(FORTRAN_PREDICATE_TERMS))


def main():
    print('Regression guard: free-surface trace-corner tensile-Tn cap (item 149, PR #212)')
    failures = []
    for check in (check_cap_fires_and_abstains_correctly, check_fortran_predicate_parity):
        try:
            check()
        except AssertionError as e:
            failures.append('%s: %s' % (check.__name__, e))
            print('  FAIL  %s: %s' % (check.__name__, e))
    if failures:
        print('\nFAIL test_tracecap_dipping_surface (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_tracecap_dipping_surface')
    return 0


if __name__ == '__main__':
    sys.exit(main())
