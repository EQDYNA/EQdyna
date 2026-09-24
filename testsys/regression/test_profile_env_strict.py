#! /usr/bin/env python3
"""
Regression guard: EQDYNA_PROFILE must be parsed STRICTLY, identically, on
both backends (lars-eriksson audit, 2026-09-23, item 3 follow-up): unset or
"" -> ON, "1" -> ON, "0" -> OFF, anything else -> hard failure naming the
variable, the bad value, and the accepted set.

Guards a real defect that shipped in 100a73b: Python's
`os.environ.get(OFF_ENV, '') != '0'` and Fortran's own startup parse
(`trim(envval) /= '0'`) both treated EQDYNA_PROFILE=off/false/2/<any typo>
as ON -- a silent fallback on bad input (rule 2). A mutation that reverts
either side back to "anything but '0' is ON" must turn the bogus-value
checks below RED, not just print a warning: `check_python_bogus_raises`
calls the real `profile_emit.enabled()` (no mock of the parse itself), and
`check_fortran_bogus_aborts_with_named_code` launches the real built
binary -- so a revert on either side is caught by running this file, not by
re-reading it.

Cheap (rule 9): the Python checks are in-process (~1 ms); the Fortran check
is one `mpirun -np 1` launch in an empty tempdir, which aborts at the
env-parse itself (eqdyna3d.f90, right after MPI_Init, before any input file
is opened) -- so it needs no case data and no full run, ~0.1-0.3 s.

SKIPPED (not silently PASS) without a built bin/eqdyna: per this tier's
existing contract (test_stop_exit_status.py's probe_real_binary), that
absence is itself scored as a FAIL, because the regression tier's declared
environment already requires the binary to be built.
"""
import glob
import os
import shlex
import subprocess
import sys
import tempfile

TESTSYS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTSYS)
SRC_FORTRAN = os.path.join(REPO_ROOT, 'src', 'fortran')
SRC_PYTHON = os.path.join(REPO_ROOT, 'src', 'python')
sys.path.insert(0, SRC_PYTHON)
from eqdyna import profile_emit  # noqa: E402

OFF_ENV = profile_emit.OFF_ENV
ACCEPTED_MARKERS = ('unset', '"1"', '"0"')
ERRCODES = os.path.join(SRC_FORTRAN, 'errorCodes.f90')


class _EnvCtx:
    """Set/unset OFF_ENV for one check, restore exactly afterward."""
    def __init__(self, value):
        self.value = value

    def __enter__(self):
        self._had = OFF_ENV in os.environ
        self._old = os.environ.get(OFF_ENV)
        if self.value is None:
            os.environ.pop(OFF_ENV, None)
        else:
            os.environ[OFF_ENV] = self.value
        return self

    def __exit__(self, *exc):
        if self._had:
            os.environ[OFF_ENV] = self._old
        else:
            os.environ.pop(OFF_ENV, None)
        return False


# --------------------------------------------------------------------------
# Python side: profile_emit.enabled()
# --------------------------------------------------------------------------
def check_python_unset_is_on():
    with _EnvCtx(None):
        assert profile_emit.enabled() is True, 'unset EQDYNA_PROFILE must be ON'


def check_python_empty_is_on():
    with _EnvCtx(''):
        assert profile_emit.enabled() is True, 'EQDYNA_PROFILE="" must be ON'


def check_python_one_is_on():
    with _EnvCtx('1'):
        assert profile_emit.enabled() is True, 'EQDYNA_PROFILE=1 must be ON'


def check_python_zero_is_off():
    with _EnvCtx('0'):
        assert profile_emit.enabled() is False, 'EQDYNA_PROFILE=0 must be OFF'


def check_python_bogus_raises():
    """THE mutation guard for the Python side: if enabled() regresses to
    `!= '0'`, every one of these returns True instead of raising, and the
    assertion below fires."""
    for bogus in ('off', 'false', 'FALSE', '2', 'yes', ' 0', '0 '):
        with _EnvCtx(bogus):
            try:
                profile_emit.enabled()
            except ValueError as exc:
                msg = str(exc)
                assert OFF_ENV in msg, (
                    'message does not name %s: %s' % (OFF_ENV, msg))
                assert bogus in msg, (
                    'message does not echo the bad value %r: %s' % (bogus, msg))
                for marker in ACCEPTED_MARKERS:
                    assert marker in msg, (
                        'message omits accepted value %s: %s' % (marker, msg))
            else:
                raise AssertionError(
                    'profile_emit.enabled() accepted EQDYNA_PROFILE=%r '
                    'instead of raising -- this is exactly the "anything '
                    'but 0 is ON" silent-fallback defect this guard exists '
                    'to catch' % bogus)


# --------------------------------------------------------------------------
# Fortran side: the real binary's startup parse (eqdyna3d.f90)
# --------------------------------------------------------------------------
def _find_binary():
    for cand in (os.path.join(REPO_ROOT, 'bin', 'eqdyna'),
                 os.path.join(SRC_FORTRAN, 'eqdyna')):
        if os.path.exists(cand):
            return cand
    return None


def _err_cfg_profile_env_invalid():
    for raw in open(ERRCODES, errors='replace'):
        s = raw.strip()
        if s.startswith('integer, parameter') and 'ERR_CFG_PROFILE_ENV_INVALID' in s:
            return int(s.split('=')[1].split('!')[0].strip())
    return None


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


def check_fortran_bogus_aborts_with_named_code():
    """Real bin/eqdyna, EQDYNA_PROFILE=bogus, an EMPTY tempdir -- legitimate
    because the parse runs right after MPI_Init, before readglobal opens any
    input file, so it aborts with ERR_CFG_PROFILE_ENV_INVALID rather than
    ERR_INPUT_FILE_MISSING. That distinction is the proof this is the
    strict-parse guard firing, not merely some unrelated refusal."""
    binary = _find_binary()
    # A missing binary is a FAILURE in this tier's declared environment (the
    # tier builds it first; a32c566 made the same call for the other guards).
    assert binary is not None, ('no built bin/eqdyna (or src/fortran/eqdyna) '
                                '-- build with ./install-eqdyna.sh')
    want = _err_cfg_profile_env_invalid()
    assert want is not None, (
        'ERR_CFG_PROFILE_ENV_INVALID not found in %s' % ERRCODES)
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    # 'bogus' plus the blank-padded values Python refuses: rule 23 requires the
    # SAME accepted set in both languages (PR #6 audit: Fortran's trim() let
    # "0 " through while profile_emit.enabled() raised).
    for bad in ('bogus', '0 ', ' 0', '1 '):
        env = dict(os.environ)
        env['EQDYNA_PROFILE'] = bad
        with tempfile.TemporaryDirectory() as d:
            try:
                rc, out = _mpirun_rank_files(mpirun, binary, d, env=env, timeout=60)
            except subprocess.TimeoutExpired:
                raise AssertionError(
                    'bin/eqdyna EQDYNA_PROFILE=%r HUNG instead of aborting' % bad)
            except RuntimeError as e:
                raise AssertionError(str(e))
        assert rc == want, (
            'bin/eqdyna EQDYNA_PROFILE=%r exited %d, expected %d '
            '(ERR_CFG_PROFILE_ENV_INVALID) -- output:\n%s'
            % (bad, rc, want, out))
        assert 'EQDYNA_PROFILE' in out, (
            'FATAL block does not name the variable for %r: %s' % (bad, out))
    return None


def main():
    print('Regression guard: EQDYNA_PROFILE strict parse (both backends)')
    py_checks = [check_python_unset_is_on, check_python_empty_is_on,
                 check_python_one_is_on, check_python_zero_is_off,
                 check_python_bogus_raises]
    rc = 0
    for c in py_checks:
        try:
            c()
            print('  PASS  %s' % c.__name__)
        except AssertionError as e:
            print('  FAIL  %s: %s' % (c.__name__, e))
            rc = 1

    try:
        note = check_fortran_bogus_aborts_with_named_code()
    except AssertionError as e:
        print('  FAIL  check_fortran_bogus_aborts_with_named_code: %s' % e)
        rc = 1
    else:
        if note is None:
            print('  PASS  check_fortran_bogus_aborts_with_named_code')
        else:
            # SKIPPED here is itself a failure: this tier's declared build
            # environment (testsys/run.py, CI) already requires a built
            # binary, so its absence is not an optional environment --
            # same contract as test_stop_exit_status.py's probe_real_binary.
            print('  FAIL  check_fortran_bogus_aborts_with_named_code: %s' % note)
            rc = 1

    print('\n%s test_profile_env_strict' % ('FAIL' if rc else 'SUCCESS'))
    return rc


if __name__ == '__main__':
    sys.exit(main())
