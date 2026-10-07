#! /usr/bin/env python3
"""
Regression guard (row 153 checkpoint 1 audit fix, victor-reyes finding):
before this PR, src/fortran/meshgen.f90's checkIsOnFault took neither its
C_degen==0 branch nor its C_degen>3.0d0 branch for 0 < C_degen <= 3 -- no
node was ever found on the fault, so the run was silently wrong (every
fault-node count would be 0). The checkpoint-1 refactor made this WORSE:
faultDegenStyle(i)/faultDegenAngle(i) (readfaultgeometry) now map that same
range to faultDegenStyle=0, which silently becomes a vertical planar fault
instead of the "no node is ever on the fault" non-behavior it had before.
src/python/eqdyna/meshgen.py's _check_is_on_fault_vec already raises
NotImplementedError for the identical range (the Python port came first);
this guard pins that BOTH languages now refuse it at INPUT TIME, before any
mesh work, with the same condition.

Refused in two places:
  1. src/fortran/readglobal (src/fortran/readInputFiles.f90), right after
     reading C_degen from bGlobal.txt -> abortRun(ERR_GEOM_DEGEN_UNSUPPORTED, ...)
  2. scripts/case.setup's create_model_input_file(), before bGlobal.txt is
     even written -> sys.exit(...) naming par.C_degen

Mutation test (both-ways, rule 14a) -- demonstrated by hand when this file
was written: commenting out the `if (C_degen > 0.0d0 .and. C_degen <=
3.0d0)` block in readInputFiles.f90 turned check_fortran_refuses_mid_range
RED (the run instead proceeded past the tiny bGlobal.txt and failed later
with ERR_INPUT_FILE_MISSING=21 on bFaultGeometry.txt, not 34); restoring it
turned the check green again. Commenting out case.setup's `if 0.0 <
par.C_degen <= 3.0` block turned check_case_setup_refuses_mid_range RED
(create_model_input_file ran to completion, writing bGlobal.txt with no
error); restoring it turned the check green again.

Cheap (rule 9): the Fortran check is one `mpirun -np 1` launch of the real
built binary against a 4-line bGlobal.txt in an empty tempdir -- it aborts
at readglobal, before any other input file is opened (same shape as
test_profile_env_strict.py's Fortran check). The case.setup check is one
subprocess launch of the real script against a minimal user_defined_params.py,
no mesh/case-input scaffolding needed since the refusal is the first thing
create_model_input_file does.

SKIPPED (scored as FAIL, not silently passed) without a built bin/eqdyna,
same contract as this tier's other Fortran-binary checks.
"""
import os
import subprocess
import sys
import tempfile

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
SRC_FORTRAN = os.path.join(REPO_ROOT, 'src', 'fortran')
SCRIPTS = os.path.join(REPO_ROOT, 'scripts')
ERRCODES = os.path.join(SRC_FORTRAN, 'errorCodes.f90')
CASE_SETUP = os.path.join(SCRIPTS, 'case.setup')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402  (item 95: rank-owned output past MPI_Abort)

MID_RANGE_VALUES = (1, 2, 3)  # C_degen is integer(kind=4) (globalvar.f90); (0, 3] contains only these three ints


def _find_binary():
    for cand in (os.path.join(REPO_ROOT, 'bin', 'eqdyna'),
                 os.path.join(SRC_FORTRAN, 'eqdyna')):
        if os.path.exists(cand):
            return cand
    return None


def _err_geom_degen_unsupported():
    for raw in open(ERRCODES, errors='replace'):
        s = raw.strip()
        if s.startswith('integer, parameter') and 'ERR_GEOM_DEGEN_UNSUPPORTED' in s:
            return int(s.split('=')[1].split('!')[0].strip())
    return None


def _minimal_bglobal(c_degen):
    # Only enough lines for readglobal to reach the C_degen read and the
    # refusal right after it -- the refusal fires before any later line is
    # read, so nothing past C_degen needs to be present or valid. C_degen is
    # integer(kind=4) (globalvar.f90) -- an int literal, not a float repr.
    return '0\n1\n0\n%d\n' % int(c_degen)


def check_fortran_refuses_mid_range():
    binary = _find_binary()
    assert binary is not None, ('no built bin/eqdyna (or src/fortran/eqdyna) '
                                '-- build with ./install-eqdyna.sh')
    want = _err_geom_degen_unsupported()
    assert want is not None, 'ERR_GEOM_DEGEN_UNSUPPORTED not found in %s' % ERRCODES
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    for c_degen in MID_RANGE_VALUES:
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'bGlobal.txt'), 'w') as f:
                f.write(_minimal_bglobal(c_degen))
            try:
                rc, out = mpirun_capture.run_rank_files(mpirun, binary, d, timeout=60)
            except subprocess.TimeoutExpired:
                raise AssertionError('bin/eqdyna C_degen=%r HUNG instead of aborting' % c_degen)
            except RuntimeError as e:
                raise AssertionError(str(e))
        assert rc == want, (
            'bin/eqdyna C_degen=%r exited %d, expected %d (ERR_GEOM_DEGEN_UNSUPPORTED) '
            '-- output:\n%s' % (c_degen, rc, want, out))
        assert 'C_degen' in out, 'FATAL block does not name C_degen for %r: %s' % (c_degen, out)


def check_fortran_accepts_boundary_values():
    """C_degen==0 and C_degen>3 (e.g. 3.1, the old tpv36/tpv37 value's
    neighbourhood) must NOT hit the new refusal -- they still fail later
    (missing bFaultGeometry.txt, ERR_INPUT_FILE_MISSING=21), proving the new
    check is scoped to (0, 3] and not accidentally widened."""
    binary = _find_binary()
    assert binary is not None, ('no built bin/eqdyna (or src/fortran/eqdyna) '
                                '-- build with ./install-eqdyna.sh')
    degen_refused = _err_geom_degen_unsupported()
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    for c_degen in (0, 15):  # 0 (planar) and 15 (tpv36/tpv37's real wedge-degeneration dip)
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'bGlobal.txt'), 'w') as f:
                f.write(_minimal_bglobal(c_degen))
            try:
                rc, out = mpirun_capture.run_rank_files(mpirun, binary, d, timeout=60)
            except subprocess.TimeoutExpired:
                raise AssertionError('bin/eqdyna C_degen=%r HUNG' % c_degen)
            except RuntimeError as e:
                raise AssertionError(str(e))
        assert rc != degen_refused, (
            'bin/eqdyna C_degen=%r (should be OUTSIDE the refused range) wrongly hit '
            'ERR_GEOM_DEGEN_UNSUPPORTED -- output:\n%s' % (c_degen, out))


_USER_PARAMS_TMPL = '''#! /usr/bin/env python3
from defaultParameters import parameters
par = parameters()
par.C_degen = %r
'''


def check_case_setup_refuses_mid_range():
    for c_degen in MID_RANGE_VALUES:
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'user_defined_params.py'), 'w') as f:
                f.write(_USER_PARAMS_TMPL % c_degen)
            env = dict(os.environ)
            env['PYTHONPATH'] = SCRIPTS + os.pathsep + env.get('PYTHONPATH', '')
            r = subprocess.run([sys.executable, CASE_SETUP], cwd=d, env=env,
                               capture_output=True, text=True, timeout=60)
        assert r.returncode != 0, (
            'case.setup with par.C_degen=%r (should be refused) exited 0 -- stdout:\n%s\nstderr:\n%s'
            % (c_degen, r.stdout, r.stderr))
        combined = r.stdout + r.stderr
        assert 'C_degen' in combined, (
            'case.setup refusal for par.C_degen=%r does not name C_degen -- stdout:\n%s\nstderr:\n%s'
            % (c_degen, r.stdout, r.stderr))


def main():
    checks = [check_fortran_refuses_mid_range,
              check_fortran_accepts_boundary_values,
              check_case_setup_refuses_mid_range]
    print('Regression guard: row 153 checkpoint 1 -- C_degen in (0,3] refused, both languages')
    rc = 0
    for c in checks:
        try:
            c()
        except AssertionError as e:
            print('  FAIL  %s: %s' % (c.__name__, e))
            rc = 1
        else:
            print('  PASS  %s' % c.__name__)
    print('\n%s test_row153_degen_range_refused' % ('FAIL' if rc else 'SUCCESS'))
    return rc


if __name__ == '__main__':
    sys.exit(main())
