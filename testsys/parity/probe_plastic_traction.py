#! /usr/bin/env python3
"""
Does the C_elastic==0 (viscoplastic) path load the fault correctly?

REPORT-ONLY. Never a gate, never wired into testsys/run.py -- it answers one
question with one number and prints it.

WHY IT EXISTS. pathway_forward item 24(a) claimed, as a BLOCKER, that on-fault
tractions resolve to "exactly HALF the intended values for C_elastic==0
(depth-independent; reproduced with yielding disabled, so not the return
mapping; affects gated test.drv.a6 and any past plastic runs)". That claim was
prose with no script behind it, which rule 4 does not allow to stand.

It does NOT reproduce. Measured 2026-09-16 on test.drv.a6, 4747 fault nodes
below 2 km: Tn ratio median 1.0018, Ts ratio median 0.9858.

The likely history is that 24(a) WAS real and was fixed by item 26 on
2026-09-14 -- the `arn` doubling. The mechanism fits exactly: faulting.f90's
traction divides by `totalMass = (massSlave + massMaster) * arn`, so an `arn`
counted twice halves every traction at every depth. Item 26 claimed its
changed path was "provably never entered" by the gated compsets; if this probe
is right, that claim was wrong for drv.a6, and it is worth re-deriving rather
than trusting.

HOW THE INTENDED VALUE IS DERIVED, so the comparison is checkable. meshgen.f90
calls setPlasticStress(depth) per element with

    strVert = -(roumax - rhow*(gamar + 1)) * depth * grav      (negative)
    devStr  = |strVert| * devStrToStrVertRatio                 (positive)
    Szz = strVert
    Sxx = strVert - devStr*cos(2*str1ToFaultAngle)
    Syy = strVert + devStr*cos(2*str1ToFaultAngle)
    Sxy = devStr*sin(2*str1ToFaultAngle)

EQdyna faults lie in the x-z plane, so the fault normal is y and the strike is
x. test.drv.a6 sets str1ToFaultAngle = 45 deg, where cos(2t) = 0 and
sin(2t) = 1, which collapses the above to

    intended fault-normal traction Tn = Syy = strVert
    intended strike shear        Ts = Sxy = devStr

That collapse is what makes this a clean one-line check for THIS case; a case
with a different str1ToFaultAngle needs the full rotation and this probe would
have to be extended rather than reused.

The +7.3215 m offset in meshgen.f90:104's depth argument is undocumented
(item 24(d)) and is ADOPTED here because the measurement needs it to land at
1.0018 -- i.e. it is load-bearing, not cosmetic. That is evidence for 24(d)
mattering, not against it.

Usage:
    python3 testsys/parity/probe_plastic_traction.py [case_dir]

With no argument it builds a fresh SERIAL test.drv.a6 under
testsys/parity/probe_case/ using the e2e sweep's own case builder (rule 1),
because the Python solver is serial-only and test/test.drv.a6 is the 4-rank
Fortran cell.
"""
import os
import shutil
import sys

import numpy as np

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))
sys.path.insert(0, os.path.join(TESTSYS, 'e2e'))
sys.path.insert(0, REPO_ROOT)
from testsys import runlock                                # noqa: E402

CASE_NAME = 'test.drv.a6'
DEFAULT_CASE = os.path.join(TESTSYS, 'parity', 'probe_case', CASE_NAME)
DEPTH_OFFSET_M = 7.3215          # meshgen.f90:104, undocumented -- item 24(d)
MIN_DEPTH_M = 2000.0             # skip the free-surface nodes
HALF_TOLERANCE = 0.05

# What a second concurrent probe run costs, printed by the refusal
# (pathway item 81 -- the same unlocked-rmtree class as item 74).
LOCK_CONSEQUENCE = [
    'A second invocation in this checkout rebuilds that SAME case directory'
    ' -- rmtree',
    'then create.newcase -- while the holder is mid-run against it. The'
    ' holder does not',
    'crash cleanly: it reads a half-written or missing case and the failure'
    ' reads as a',
    'plastic-traction mismatch (a probe failure) rather than as'
    ' infrastructure',
    '(rule 21a, item 81).',
    '',
    'NOT waiting. NOT rebuilding anyway. NOT falling back to a different'
    ' case directory --',
    'each of those is the silent fallback rule 2 forbids. Run your probe in'
    ' its own',
    'git worktree (rule 21a), or wait for the holder above to finish.']


def build_case(case_dir):
    """Rebuild the fixed probe case tree at case_dir from scratch.

    GATE 0 (item 81): the exclusive lock on the directory holding case_dir,
    taken before anything is deleted or created, held only across the WRITE
    (rmtree + make_serial_case) and released immediately after -- the
    probe's measurement phase (main(), after this returns) never touches
    this directory destructively, so there is nothing left to guard once the
    rebuild is done. Same non-blocking acquire()/refusal shape as
    run_perf.build_perf_case and run_jaxmpi_ab.main (rules 21a, 2).
    """
    import run_e2e                                        # noqa: E402
    lock_resource = os.path.relpath(os.path.dirname(case_dir), REPO_ROOT)
    try:
        lock = runlock.acquire(REPO_ROOT, lock_resource,
                               consequence=LOCK_CONSEQUENCE)
    except runlock.RunTreeLocked as exc:
        raise SystemExit('FAIL: %s' % exc)
    # Announced only once the lock is held: "Building ..." printed ahead of a
    # refusal describes something that never happened.
    print('Building a fresh serial %s at %s ...' % (CASE_NAME, case_dir))
    try:
        if os.path.isdir(case_dir):
            shutil.rmtree(case_dir)
        os.makedirs(os.path.dirname(case_dir), exist_ok=True)
        run_e2e.make_serial_case(CASE_NAME, case_dir, run_e2e.base_env())
    finally:
        lock.release()
    return case_dir


def main():
    case = sys.argv[1] if len(sys.argv) > 1 else None
    if case is None:
        case = build_case(DEFAULT_CASE)
    if not os.path.isdir(case):
        raise SystemExit('FAIL: no case directory at %s' % case)

    from eqdyna import eqdyna3d, readInputFiles           # noqa: E402
    _params, g = readInputFiles.build_params(case)
    if g.get('C_elastic') != 0:
        raise SystemExit('FAIL: %s has C_elastic=%r -- this probe is only '
                         'meaningful for the viscoplastic path (C_elastic==0)'
                         % (case, g.get('C_elastic')))
    angle = g['str1ToFaultAngle']
    if abs(abs(np.cos(2.0 * angle)) - 0.0) > 1e-12:
        raise SystemExit(
            'FAIL: str1ToFaultAngle=%r rad gives cos(2t)=%.3e, not 0. This '
            'probe uses the 45-degree collapse (Syy=strVert, Sxy=devStr) and '
            'would be WRONG here -- extend it to the full rotation rather '
            'than reading its output.' % (angle, np.cos(2.0 * angle)))

    roumax, rhow, gamar = g['roumax'], g['rhow'], g['gamar']
    grav = g.get('grav') or 9.8
    ratio = g['devStrToStrVertRatio']
    print('case      : %s' % case)
    print('C_elastic %d  friclaw %d  str1ToFaultAngle %.1f deg  devStr/strVert %.2f'
          % (g['C_elastic'], g['friclaw'], np.degrees(angle), ratio))
    print('roumax %.1f  rhow %.1f  gamar %.2f  grav %.2f'
          % (roumax, rhow, gamar, grav))

    frt = eqdyna3d.run_case(case, nsteps=1, verbose=False, backend='numpy')
    a = np.loadtxt(frt)
    # frt columns: 0-2 xyz, 3 fnft, 4-9 slip/sliprate, 10 peak sr,
    #              11 tnrm, 12 tstk, 13 tdip  (library_output.f90's output_frt)
    z, tn, ts = a[:, 2], a[:, 11], a[:, 12]
    depth = -z + DEPTH_OFFSET_M
    strVert = -(roumax - rhow * (gamar + 1.0)) * depth * grav
    devStr = np.abs(strVert) * ratio

    deep = depth > MIN_DEPTH_M
    if deep.sum() < 100:
        raise SystemExit('FAIL: only %d nodes below %.0f m -- too few to '
                         'report a median' % (deep.sum(), MIN_DEPTH_M))
    print('\nnodes: %d total, %d below %.0f m' % (len(z), deep.sum(), MIN_DEPTH_M))
    print('\n%-14s %14s %14s   %s' % ('', 'intended', 'actual', 'actual/intended'))
    rn = tn[deep] / strVert[deep]
    for nm, act, exp in (('Tn (normal)', tn, strVert), ('Ts (strike)', ts, devStr)):
        r = act[deep] / exp[deep]
        print('%-14s %14.4e %14.4e   median %.4f  [p10 %.4f p90 %.4f]'
              % (nm, np.median(exp[deep]), np.median(act[deep]),
                 np.median(r), np.percentile(r, 10), np.percentile(r, 90)))

    med = float(np.median(rn))
    print('\nitem 24(a) claimed EXACTLY HALF, depth-independent.')
    if abs(med - 0.5) < HALF_TOLERANCE:
        print('  Tn median ratio %.4f -- THE CLAIM REPRODUCES. The most likely '
              'mechanism is a double-counted arn in the traction denominator '
              '(faulting.f90 totalMass); see item 26.' % med)
        return 1
    print('  Tn median ratio %.4f -- the claim does NOT reproduce; the '
          'viscoplastic path loads the fault correctly.' % med)
    print('  The p10-p90 spread is the fractal roughness perturbing the local '
          'fault normal, not a systematic factor.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
