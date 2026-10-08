#! /usr/bin/env python3
"""
Regression guard (row 153 checkpoint 2, three mesh bugs in the style-2
x-y branch mesh used by test.tpv24/test.tpv25): a style-2 fault (a second,
non-degenerating planar fault alongside a branch fault meshed with
faultDegenStyle==2) had THREE independent defects, each found by running
the real test.tpv24 case and reading the Fortran binary's own hard checks
(ERR_MESH_BAD_JACOBIAN=47, ERR_NUM_VELOCITY_NAN=63) rather than inferring
from a synthetic fixture:

  1. (59e590e) Fault 2's box used fymin=fymax=0.0 (a copy of fault 1's own
     planar y=0 box), which excludes every branch node (the branch line is
     strictly y<0 in its x-range) -- ZERO fault-2 nodes, no faultstft2_*.txt
     station file written at all.
  2. (dd01dcc) The style-2 wedge neworder permutation had the wrong z
     handedness for this port's actual z-grid convention -- NEGATIVE element
     Jacobian determinant (element 2553 on this fixture, det=-3.608e7),
     ERR_MESH_BAD_JACOBIAN=47.
  3. (7db4db3) replaceSlaveWithMasterNode's type-12/13 "always substitute"
     branch was applied across ALL faults, not scoped to the fault that
     produced the retag -- an ordinary fault-1 slave node near the branch
     got wrongly promoted to fault-2's master, orphaning the true fault-1
     slave node: zero mass -> NaN velocity, ERR_NUM_VELOCITY_NAN=63, at
     node (1000, ~0, -14000), time step 2.

Fixture: the real test.tpv24 case (case_input/test.tpv24), run as-is except
par.term is overridden to 3*par.dt (3 time steps -- enough for bug 3's NaN
to appear, same mesh/geometry as the gate case, no synthetic substitute).
4 ranks (par.nx,ny,nz=2,1,2, the case's own partition). Runs in ~1-2 s.

MUTATION TEST (both-ways, rule 14a) -- demonstrated by hand when this file
was written, against this EXACT fixture, by building bin/eqdyna at each
commit in a throwaway `git worktree` and rerunning:
  2be3cf7 (all three bugs present):        rc=47, element 2553, det=-36084391.82
  59e590e (bug 1 fixed, 2+3 present):      rc=47, element 2553, det=-36084391.82
  dd01dcc (bugs 1+2 fixed, 3 present):     rc=63, NaN at (1000, -2.5e-12, -14000), step 2
  7db4db3 (all three fixed):               rc=0, faultstft2_020dp000.txt present,
                                            680 frt rows across 4 ranks, no NaN.
This pins all three fixes in one place: reverting any ONE of them on this
fixture reproduces its own documented rc/location above (checked directly
for bugs 1-3's commits; a regression in any of the three files touched by
7db4db3/dd01dcc/59e590e is expected to reproduce one of these codes, not a
fourth, undocumented failure mode).

SKIPPED (scored as FAIL, not silently passed) without a built bin/eqdyna,
same contract as this tier's other Fortran-binary checks.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')
CASE_INPUT = os.path.join(REPO_ROOT, 'case_input', 'test.tpv24')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402


def _find_binary():
    for cand in (os.path.join(REPO_ROOT, 'bin', 'eqdyna'),
                 os.path.join(REPO_ROOT, 'src', 'fortran', 'eqdyna')):
        if os.path.exists(cand):
            return cand
    return None


def _build_tiny_case(d):
    """create.newcase-equivalent for test.tpv24, with par.term cut to 3
    steps -- cheap enough to run every day, same geometry as the gate case."""
    for f in os.listdir(CASE_INPUT):
        src = os.path.join(CASE_INPUT, f)
        if os.path.isfile(src):
            shutil.copy(src, d)
    for f in os.listdir(SCRIPTS):
        src = os.path.join(SCRIPTS, f)
        if os.path.isfile(src):
            shutil.copy(src, d)
    with open(os.path.join(d, 'user_defined_params.py'), 'a') as fh:
        fh.write('\npar.term = 3 * par.dt  # row153 ckpt2 mesh-bug regression: 3 steps is enough\n')
    env = dict(os.environ)
    env['PYTHONPATH'] = d + os.pathsep + env.get('PYTHONPATH', '')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=d, env=env,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, 'case.setup failed:\n%s\n%s' % (r.stdout, r.stderr)


def check_style2_branch_mesh():
    binary = _find_binary()
    assert binary is not None, ('no built bin/eqdyna (or src/fortran/eqdyna) '
                                '-- build with ./install-eqdyna.sh')
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    with tempfile.TemporaryDirectory() as d:
        _build_tiny_case(d)
        try:
            rc, out = mpirun_capture.run_rank_files(mpirun, binary, d, np=4, timeout=120)
        except subprocess.TimeoutExpired:
            raise AssertionError('test.tpv24 tiny fixture HUNG instead of completing')
        except RuntimeError as e:
            raise AssertionError(str(e))
        assert rc == 0, (
            'test.tpv24 tiny fixture (3 steps) exited %d, expected 0 -- a style-2 '
            'branch-mesh regression (bug 1/2/3 above) is the first suspect -- '
            'output:\n%s' % (rc, out[-2000:]))

        # bug 1 guard: fault-2 must have produced at least one on-fault
        # station file (zero fault-2 nodes silently drops every one).
        ft2_files = glob.glob(os.path.join(d, 'faultstft2_*.txt'))
        assert ft2_files, (
            'no faultstft2_*.txt station file written -- fault 2 has zero nodes '
            '(bug 1: fymin/fymax box regression)')

        # bugs 2/3 guard, redundant with rc==0 but checked explicitly so a
        # future change to error-code plumbing cannot silently defeat this
        # test by returning 0 on a real failure: frt rows exist and are
        # finite (a negative-Jacobian or zero-mass/NaN node, if it did not
        # abort for some other reason, would show up here as non-finite).
        frt_files = sorted(glob.glob(os.path.join(d, 'frt.txt*')))
        assert frt_files, 'no frt.txt<rank> files written at all'
        total_rows = 0
        for ff in frt_files:
            with open(ff) as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    total_rows += 1
                    vals = line.split()
                    assert all(_is_finite(v) for v in vals), (
                        'non-finite value in %s: %r' % (ff, line))
        assert total_rows > 0, 'frt.txt<rank> files are empty -- zero fault nodes total'


def _is_finite(tok):
    try:
        v = float(tok)
    except ValueError:
        return False
    return v == v and abs(v) != float('inf')


def main():
    print('Regression guard: row 153 checkpoint 2 -- style-2 branch mesh (3 bugs), test.tpv24 fixture')
    rc = 0
    try:
        check_style2_branch_mesh()
    except AssertionError as e:
        print('  FAIL  check_style2_branch_mesh: %s' % e)
        rc = 1
    else:
        print('  PASS  check_style2_branch_mesh')
    print('\n%s test_row153_style2_branch_mesh' % ('FAIL' if rc else 'SUCCESS'))
    return rc


if __name__ == '__main__':
    sys.exit(main())
