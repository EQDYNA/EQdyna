#! /usr/bin/env python3
"""
Regression guard: a WEDGE-DEGENERATE dipping fault may be split in y (rules 2, 10).

WHAT THIS PINS, and why it is the DIVIDE case.

MPI4arn's syncArnBoundary (src/fortran/meshgen.f90) adds a neighbour's arn
contribution back only when the fault has non-zero NOMINAL GRID extent in that
boundary's direction. The discriminator is fltxyz, and it is NOT the physical
dip -- the two ways checkIsOnFault selects fault nodes give opposite answers:

  insertFaultType > 0, C_degen == 0   (test.tpv10, dip 60)
      nodes chosen by `nodeCoor(2) == 0.0d0`, EXACT equality on the unblended
      grid y. The fault is ONE y-index plane however steeply it dips;
      insertFaultInterface displaces only the physical y. fltxyz y-extent is
      0.0/0.0, a y boundary COINCIDES with the fault, both ranks hold the
      COMPLETE tributary area -> DUPLICATE, no add-back. Guarded by
      test_fault_mpi_boundary_arn.py.

  C_degen > 3   (wedge degeneration; test.tpv36, test.tpv37, dip 15)
      nodes chosen by `|z + y*tan(C_degen)| < dx/100`. The fault SPANS a range
      of grid y, fltxyz y-extent is faultWidth*cos(dip), a y boundary CROSSES
      it, each rank holds a PARTIAL area -> DIVIDE, add back. THIS FILE.

checkFaultMPIAlignment (eqdyna3d.f90) treats the DIVIDE case as a rank-0
NOTICE, not a hard stop: a duplicated fault surface would have doubled arn
and halved every traction (a clean 0.5), and audited evidence (tpv36,
(npx,npy,npz)=(2,2,1) vs (2,1,2): tnrm/tstk/tdip ratio 1.000000, max|diff|
over all 22 canonical columns 1.0e-08 = output precision) confirms the
add-back is correct.

arn is a MESH-TIME quantity, computed once before any time stepping, so a
short run is a complete test of it: were arn doubled, the t=0 tractions would
already be halved. term is kept small deliberately.

Cheap-ish (rule 9): two 4-rank runs of a 1 s case. Minutes, not seconds --
it needs a real mesh and a real MPI decomposition, which is the whole point.
Skips with a LOUD notice, never a pass, if the binary is absent.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

CASE = 'test.tpv36'
TERM_S = 1.0
# ratio tolerance: a duplicate/divide error is a factor of 2, so anything near
# 1 settles it. This is deliberately loose on the RATIO and strict on the
# absolute agreement below.
RATIO_TOL = 0.01
ABS_TOL = 1.0e-6


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _binary():
    for c in (os.path.join(ROOT, 'bin', 'eqdyna'),
              os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        if os.path.exists(c):
            return c
    return None


def _run(tmp, tag, decomp):
    case_dir = os.path.join(tmp, tag)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE], env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed: %s' % (r.stderr or '')[-500:])
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'a') as f:
        f.write('\n# forced by test_dipping_fault_y_split.py\n'
                'par.nx, par.ny, par.nz = %d, %d, %d\npar.term = %r\n'
                % (decomp[0], decomp[1], decomp[2], TERM_S))
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed for %s: %s' % (tag, (r.stderr or '')[-500:]))
    nranks = decomp[0] * decomp[1] * decomp[2]
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    r = subprocess.run([mpirun, '-np', str(nranks), _binary()],
                       cwd=case_dir, env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(
            'eqdyna exited %d for %s at decomposition %r.\n%s'
            % (r.returncode, tag, decomp, (r.stdout + r.stderr)[-1500:]))
    return case_dir


def main():
    print('Regression guard: a wedge-degenerate dipping fault may be split in y')
    if _binary() is None:
        print('  SKIPPED (no bin/eqdyna or src/fortran/eqdyna; build it to enable).')
        print('  This is a SKIP, not a pass -- the DIVIDE branch is unverified in '
              'this run.')
        return 0
    from testsys import frt_canonical
    import numpy as np
    try:
        with tempfile.TemporaryDirectory() as tmp:
            d_split = _run(tmp, 'ysplit', (2, 2, 1))     # y IS split: DIVIDE path
            d_plain = _run(tmp, 'noysplit', (2, 1, 2))   # y not split: reference
            a = frt_canonical.canonical_from_case(d_plain)
            b = frt_canonical.canonical_from_case(d_split)
            if a.shape != b.shape:
                raise AssertionError('node counts differ: %r vs %r -- different '
                                     'discretisation, not a tolerance question'
                                     % (a.shape, b.shape))
            ac, bc = frt_canonical.align(a, b)
            worst = float(np.abs(ac - bc).max())
            bad = []
            for col, nm in ((11, 'tnrm'), (12, 'tstk'), (13, 'tdip')):
                m = np.abs(ac[:, col]) > 1.0
                if m.sum() < 50:
                    continue
                r = bc[m, col] / ac[m, col]
                med = float(np.median(r))
                if abs(med - 1.0) > RATIO_TOL:
                    bad.append('%s ratio median %.6f (n=%d) -- 0.5 would mean the '
                               'y-split DUPLICATED the fault surface and arn was '
                               'doubled' % (nm, med, m.sum()))
                else:
                    print('  %s ratio median %.6f over %d nodes' % (nm, med, m.sum()))
            if worst > ABS_TOL:
                bad.append('max |diff| over all 22 columns is %.3e, above %.1e'
                           % (worst, ABS_TOL))
            else:
                print('  max |diff| over all 22 columns %.3e' % worst)
            if bad:
                for b_ in bad:
                    print('  FAIL  %s' % b_)
                print('\nFAIL test_dipping_fault_y_split')
                return 1
    except AssertionError as e:
        print('  FAIL  %s' % e)
        print('\nFAIL test_dipping_fault_y_split')
        return 1
    print('\nSUCCESS test_dipping_fault_y_split '
          '(DIVIDE branch: y-split == unsplit)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
