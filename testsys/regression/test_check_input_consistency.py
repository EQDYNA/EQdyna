#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 23): the Python port must refuse
EXACTLY what src/fortran/checkInputConsistency.f90 refuses -- same
condition, same message text, same numbered exit code
(src/fortran/errorCodes.f90), at the same point in the run (before any
mesh or solver work).

Three checks are ported (checkInputConsistency.f90:7-19 <->
src/python/eqdyna/checkInputConsistency.py's `check`):

  1. C_elastic==0 and C_Q==1              -> ERR_CFG_Q_NEEDS_ELASTIC (11)
  2. C_Q==1 and rat>1                      -> ERR_CFG_Q_NEEDS_UNIFORM (12)
  3. output_plastic==1 and C_elastic!=0    -> ERR_CFG_PLASTIC_OUTPUT (13)

Checks 1/2 are UNREACHABLE through any real case: `C_Q` is not a
case-input field on either side (grepping src/fortran/readInputFiles.f90
and scripts/case.setup finds no reader/writer for it anywhere;
globalvar.f90:104 hardcodes it `= 0`, matching
docs/fortran_python_correspondence.md's note on calcQAttenuationCoeff.f90).
Building a case that sets C_Q would test a hook that does not exist in the
real system, so 1/2 are pinned as DIRECT unit calls against
checkInputConsistency.check() (forcing C_Q=1 via its keyword), not against
a case directory -- check_unreachable_q_checks below.

Check 3 IS reachable: scripts/case.setup writes both par.C_elastic and
par.output_plastic to bGlobal.txt unconditionally. It is pinned END TO
END: ONE minimal case, run through the REAL Fortran binary (bin/eqdyna,
built by ./install-eqdyna.sh) AND the REAL python port
(`python3 -m eqdyna ... --backend numpy`), asserting BOTH refuse with exit
code 13 and BOTH name the same condition. A second run of the identical
case with output_plastic=0 (valid) must have BOTH binaries run to
completion and write their frt output -- proving the guard only refuses,
it does not silently change physics for a valid config (the "frt
unchanged" requirement: this check gates NOTHING for a valid case, so its
frt is compared to nothing new here -- the full parity gate is
testsys/e2e/run_e2e.py's sweep, run separately).

Cheap (rule 9): a 5x5-node planar-fault case (dx=500m, 4 dt term for the
valid run so the python-numpy backend does not spend its time in the
solve loop), no mesh work needed for the refused runs at all (the check
fires before meshgen on both sides). Builds ONLY bin/eqdyna (already built
by the caller's ./install-eqdyna.sh); fails loudly, not silently, if it is
absent.

Mutation test (run by hand, not by this file, per the mission brief): comment
out the `output_plastic == 1 and C_elastic != 0` branch in
src/python/eqdyna/checkInputConsistency.py and re-run this script -- expect
check_plastic_output_refused_end_to_end to go RED naming exactly that
condition (the python side no longer refuses while the Fortran side still
does).
"""
import glob
import os
import shlex
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))

from eqdyna import checkInputConsistency as cic  # noqa: E402

MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE_NAME = 'test.tpv8'

PLASTIC_MSG = ('Plastic strains are only output for C_elastic=0. Set '
               'output_plastic=0 or C_elastic=0.')
ERR_CFG_PLASTIC_OUTPUT = 13

USER_PARAMS = '''#! /usr/bin/env python3
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -3.0e3, 3.0e3
par.ymin, par.ymax = -3.0e3, 3.1e3
par.zmin, par.zmax = -3.0e3, 0.0e3

par.fxmin, par.fxmax = -1.0e3, 1.0e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -2.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -500.0

par.dx = 500.
par.dy = par.dx
par.dz = par.dx
par.nmat = 1
par.vp, par.vs, par.rou = 5716, 3300, 2700

par.dt = 0.5*par.dx/par.vp
par.term = 4.*par.dt
par.friclaw = 1
par.tpv = 8
par.nucR = 300.

par.nx, par.ny, par.nz = 1, 1, 1

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin,par.fxmax,par.nfx)
par.fz  = np.linspace(par.fzmin,par.fzmax,par.nfz)

par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
par.fric_sw_fs = 0.76
par.fric_sw_fd = 0.448
par.fric_sw_D0 = 0.5
par.grav = 9.8
par.fric_cohesion = 1.e6
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    par.on_fault_vars[iz,ix,1]   = par.fric_sw_fs
    if abs(abs(xcoor) - par.fxmax)<0.01 or abs(zcoor-par.fzmin)<0.01:
        par.on_fault_vars[iz,ix,1] = 1000.
    par.on_fault_vars[iz,ix,2]   = par.fric_sw_fd
    par.on_fault_vars[iz,ix,3]   = par.fric_sw_D0
    par.on_fault_vars[iz,ix,4]   = par.fric_cohesion
    par.on_fault_vars[iz,ix,7]   = 7378.*zcoor
    par.on_fault_vars[iz,ix,8]   = abs(0.55*par.on_fault_vars[iz,ix,7])
    if abs(xcoor-par.xsource)<=300. and abs(zcoor-par.zsource)<=300.:
        par.on_fault_vars[iz,ix,8] = 1e6+abs(1.005*0.76*par.on_fault_vars[iz,ix,7])

par.st_coor_on_fault = [[0.0, 0.0]]
par.st_coor_off_fault = [[0,1,0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)

par.output_plastic = %d
'''


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(tmp, name, output_plastic):
    case_dir = os.path.join(tmp, 'case_%s' % name)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('create.newcase failed (%d):\n%s'
                           % (r.returncode, r.stdout[-2000:] + r.stderr[-2000:]))
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'w') as f:
        f.write(USER_PARAMS % output_plastic)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('case.setup failed for %s:\n%s' % (name, r.stdout + r.stderr))
    return case_dir


def _mpirun_rank_files(mpirun, binary, cwd, env=None, timeout=120):
    """(returncode, text) for `mpirun -np 1 <binary>` with the rank's stdout and
    stderr redirected to a rank-owned FILE (pathway item 95; the technique is
    test_stop_exit_status.py probe B, item 89). Open MPI discards
    already-flushed bytes it has not yet forwarded when MPI_Abort tears the
    job down, so text asserted from mpirun's piped capture past a refusal is
    a launcher property that can flake. The file cannot lose bytes the rank
    flushed. `exec` keeps the binary's own exit status; `$$` names the file.
    Raises when no rank file appears: then the binary never ran."""
    before = set(glob.glob(os.path.join(cwd, 'eqdyna-stdout.*.txt')))
    shell = 'exec %s > eqdyna-stdout.$$.txt 2>&1' % shlex.quote(binary)
    r = subprocess.run([mpirun, '-np', '1', 'sh', '-c', shell], cwd=cwd, env=env,
                       capture_output=True, text=True, timeout=timeout)
    files = sorted(set(glob.glob(os.path.join(cwd, 'eqdyna-stdout.*.txt'))) - before)
    if not files:
        raise RuntimeError('no rank-owned stdout file under %s -- %s never ran '
                           '(mpirun rc=%d): %s' % (cwd, binary, r.returncode,
                                                   (r.stdout + r.stderr)[-500:]))
    text = ''.join(open(f, errors='replace').read() for f in files)
    return r.returncode, text + r.stdout + r.stderr


def _run_fortran(case_dir):
    binary = None
    for cand in (os.path.join(ROOT, 'bin', 'eqdyna'),
                 os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        if os.path.exists(cand):
            binary = cand
            break
    if binary is None:
        return None
    return _mpirun_rank_files(MPIRUN, binary, case_dir)


def _run_python(case_dir, nsteps=4):
    r = subprocess.run([sys.executable, '-u', '-m', 'eqdyna', case_dir, str(nsteps),
                        '--backend', 'numpy'],
                       cwd=ROOT, env=_env(), capture_output=True, text=True, timeout=120)
    return r.returncode, r.stdout + r.stderr


def _fortran_truth():
    """Codes from errorCodes.f90 and reason strings from
    checkInputConsistency.f90, parsed from the Fortran SOURCE -- not from the
    Python module under test, so the guard cannot agree with itself."""
    import re
    codes = dict(re.findall(r'(ERR_CFG_\w+)\s*=\s*(\d+)',
                            open(os.path.join(ROOT, 'src', 'fortran', 'errorCodes.f90')).read()))
    src = open(os.path.join(ROOT, 'src', 'fortran', 'checkInputConsistency.f90')).read()
    calls = re.findall(r"abortRun\((ERR_CFG_\w+),\s*&\s*'([^']*)'\)", src)
    return [(name, int(codes[name]), msg) for name, msg in calls]


def check_unreachable_q_checks():
    truth = _fortran_truth()
    if len(truth) != 3:
        raise AssertionError('parsed %d abortRun calls from checkInputConsistency.f90, '
                             'expected 3: %r' % (len(truth), truth))
    fails = []
    cases = [(dict(C_elastic=0, output_plastic=0, rat=1.0, C_Q=1), truth[0]),
             (dict(C_elastic=1, output_plastic=0, rat=1.025, C_Q=1), truth[1]),
             (dict(C_elastic=1, output_plastic=1, rat=1.0, C_Q=0), truth[2])]
    for kw, (name, code, msg) in cases:
        try:
            cic.check(**kw)
            fails.append('%s: %r did not raise' % (name, kw))
        except cic.InputConsistencyError as e:
            if e.code != code:
                fails.append('%s: raised code %r, Fortran errorCodes.f90 says %d' % (name, e.code, code))
            if str(e) != msg:
                fails.append('%s: message %r != Fortran %r' % (name, str(e), msg))
    try:
        cic.check(C_elastic=1, output_plastic=0, rat=1.0, C_Q=1)
    except cic.InputConsistencyError as e:
        fails.append('check 2 raised at rat==1.0 (boundary; must only fire for rat>1): %r' % e)
    if fails:
        raise AssertionError('; '.join(fails))
    print('  PASS  all 3 checks: code and message equal to the Fortran source '
          '(errorCodes.f90 + checkInputConsistency.f90), by direct unit call')


def check_plastic_output_refused_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'plastic_bad', output_plastic=1)

        rc_f, out_f = _run_fortran(case_dir)
        if rc_f is None:
            raise AssertionError('bin/eqdyna (or src/fortran/eqdyna) not found; '
                                 'build it first with ./install-eqdyna.sh')
        if rc_f == 0:
            raise AssertionError('Fortran ACCEPTED output_plastic=1 with '
                                 'C_elastic=1 (exit 0) -- checkInputConsistency.f90 '
                                 'should have refused this')
        if rc_f != ERR_CFG_PLASTIC_OUTPUT:
            raise AssertionError('Fortran exited %d, expected %d (ERR_CFG_PLASTIC_OUTPUT)'
                                 % (rc_f, ERR_CFG_PLASTIC_OUTPUT))
        if PLASTIC_MSG not in out_f:
            raise AssertionError('Fortran refused (exit %d) but did not print the '
                                 'expected message:\n%s' % (rc_f, out_f[-1500:]))

        rc_p, out_p = _run_python(case_dir)
        if rc_p == 0:
            raise AssertionError('python port ACCEPTED output_plastic=1 with '
                                 'C_elastic=1 (exit 0) -- checkInputConsistency.py '
                                 'should have refused this')
        if rc_p != ERR_CFG_PLASTIC_OUTPUT:
            raise AssertionError('python port exited %d, expected %d '
                                 '(ERR_CFG_PLASTIC_OUTPUT):\n%s'
                                 % (rc_p, ERR_CFG_PLASTIC_OUTPUT, out_p[-1500:]))
        if PLASTIC_MSG not in out_p:
            raise AssertionError('python port refused (exit %d) but did not print '
                                 'the expected message:\n%s' % (rc_p, out_p[-1500:]))
        # The rest of Fortran's output, which the reason string alone would
        # not catch (PR #9 audit): the C_elastic echo and the FATAL block's
        # rank line, on BOTH sides.
        import re
        for tag, out in (('fortran', out_f), ('python', out_p)):
            if not re.search(r'Now, C_elastic =\s+1', out):
                raise AssertionError('%s did not echo "Now, C_elastic = 1" '
                                     '(checkInputConsistency.f90:16):\n%s' % (tag, out[-1500:]))
            if not re.search(r'rank\s+:\s+0', out):
                raise AssertionError('%s FATAL block has no "rank : 0" line '
                                     '(errorCodes.f90:156):\n%s' % (tag, out[-1500:]))
    print('  PASS  output_plastic=1/C_elastic=1 refused by BOTH binaries, '
          'exit %d, message matches' % ERR_CFG_PLASTIC_OUTPUT)


def check_valid_config_still_runs_both():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'valid', output_plastic=0)

        rc_f, out_f = _run_fortran(case_dir)
        if rc_f is None:
            raise AssertionError('bin/eqdyna (or src/fortran/eqdyna) not found; '
                                 'build it first with ./install-eqdyna.sh')
        if rc_f != 0:
            raise AssertionError('Fortran REJECTED a valid config (exit %d):\n%s'
                                 % (rc_f, out_f[-2000:]))
        frt_f = [f for f in os.listdir(case_dir) if f.startswith('frt.txt')]
        if not frt_f:
            raise AssertionError('Fortran run wrote no frt.txt* for a valid config')

        rc_p, out_p = _run_python(case_dir)
        if rc_p != 0:
            raise AssertionError('python port REJECTED a valid config (exit %d):\n%s'
                                 % (rc_p, out_p[-2000:]))
        if not os.path.isfile(os.path.join(case_dir, 'frt.txt0')):
            raise AssertionError('python port wrote no frt.txt0 for a valid config:\n%s'
                                 % out_p[-2000:])
    print('  PASS  valid config (output_plastic=0, C_elastic=1) runs to '
          'completion on both binaries and writes frt output')


def main():
    print('Regression guard: checkInputConsistency.f90 <-> checkInputConsistency.py')
    if not os.path.exists(os.path.join(ROOT, 'bin', 'eqdyna')) and \
       not os.path.exists(os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        print('FAIL test_check_input_consistency: bin/eqdyna (or src/fortran/eqdyna) '
              'not found -- build it with ./install-eqdyna.sh first (this guard does '
              'not skip silently)')
        return 1
    fails = []
    for c in (check_unreachable_q_checks,
              check_plastic_output_refused_end_to_end,
              check_valid_config_still_runs_both):
        try:
            c()
        except (AssertionError, RuntimeError) as e:
            fails.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if fails:
        print('\nFAIL test_check_input_consistency (%d check(s))' % len(fails))
        return 1
    print('\nSUCCESS test_check_input_consistency')
    return 0


if __name__ == '__main__':
    sys.exit(main())
