#! /usr/bin/env python3
"""
Regression guard: ntotft > 1 must WORK (not be refused) at the case.setup
level -- pathway_forward.md item 17's multi-fault campaign landed (Row 17).

HISTORY: this test used to pin the OPPOSITE -- that case.setup refused
ntotft > 1 outright, because its bStations.txt writer emitted one on-fault
station count where readInputFiles.f90:152 reads `(nonfs(i), i = 1, ntotft)`.
That writer/reader mismatch is fixed (scripts/lib.py's
resolveOnFaultStationsPerFault + case.setup's create_station_input_file now
write ntotft counts, grouped by fault); bFaultGeometry.txt similarly now
takes a per-fault box list (resolveFaultGeom / par.faultgeom) instead of
repeating one box for every fault. Leaving the old refusal test green after
removing the refusal would have been a silent-pass gate violation (the
refusal it asserted no longer exists), so it is replaced here with the
positive contract: ntotft=2 case.setup must SUCCEED and write files that
structurally match what the Fortran readers expect, and ntotft=1 must still
work identically (the fallback paths in resolveFaultGeom /
resolveOnFaultStationsPerFault / resolveOnFaultVarsPerFault are no-ops there).

This test is at the case.setup (file-writing) layer only -- cheap (rule 9):
two `create.newcase` + `case.setup` invocations, no build, no MPI, no
simulation. The full end-to-end 2-fault RUN (mesh, solve, output, and a
committed reference) is testsys/regression/test_multifault_two_fault_smoke.py,
against the dedicated case_input/test.multifault2 compset; this test uses
test.tpv8 with an appended par.ntotft=2 override (case.setup's own contract,
independent of any one compset) to keep the two tests orthogonal.

Exits non-zero on any failure.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_NAME = 'test.tpv8'


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(tmp, extra_params):
    """create.newcase + append extra_params + run the case's OWN case.setup.

    The case's copy is what a user runs, and create.newcase freezes a copy into
    the case directory -- so testing scripts/case.setup directly would test a
    file no user invokes.
    """
    case_dir = os.path.join(tmp, CASE_NAME)
    rc = subprocess.run([sys.executable,
                         os.path.join(ROOT, 'scripts', 'create.newcase'),
                         case_dir, CASE_NAME],
                        env=_env(), capture_output=True, text=True)
    if rc.returncode != 0:
        raise AssertionError('create.newcase failed (%d):\n%s'
                             % (rc.returncode, rc.stdout[-2000:] + rc.stderr[-2000:]))
    if extra_params:
        with open(os.path.join(case_dir, 'user_defined_params.py'), 'a') as f:
            f.write('\n' + extra_params + '\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=_env(), capture_output=True, text=True)
    return case_dir, r


TWO_FAULT_OVERRIDE = (
    "par.ntotft = 2\n"
    "par.nucfault = 1\n"
    "par.faultgeom = [(par.fxmin, par.fxmax, par.fymin, par.fymax, par.fzmin, par.fzmax),\n"
    "                 (par.fxmin, par.fxmax, 2000.0, 2000.0, par.fzmin, par.fzmax)]\n"
    "par.onFaultVarsPerFault = [par.on_fault_vars, par.on_fault_vars.copy()]\n"
    "par.st_coor_on_fault_per_fault = [par.st_coor_on_fault, [[0.0, -7.5]]]\n"
)


def check_two_faults_case_setup_succeeds():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir, r = _make_case(tmp, TWO_FAULT_OVERRIDE)
        if r.returncode != 0:
            raise AssertionError(
                'case.setup FAILED on par.ntotft=2 with faultgeom/onFaultVarsPerFault/'
                'st_coor_on_fault_per_fault all supplied (exit %d) -- this is the documented '
                'multi-fault contract (scripts/lib.py: resolveFaultGeom, '
                'resolveOnFaultVarsPerFault, resolveOnFaultStationsPerFault) and must work:\n%s'
                % (r.returncode, (r.stdout + r.stderr)[-3000:]))

        # bFaultGeometry.txt: two distinct boxes, second offset in y. Each
        # block is 5 lines (header + x/y/z bounds + the Row 153 checkpoint 2a
        # per-fault degeneration code), not 4.
        geom = open(os.path.join(case_dir, 'bFaultGeometry.txt')).read().split('\n')
        nonblank = [l for l in geom if l.strip()]
        if 'For fault No. 1' not in nonblank[0] or 'For fault No. 2' not in nonblank[5]:
            raise AssertionError('bFaultGeometry.txt does not have two "For fault No." blocks:\n%s'
                                 % geom)
        fault2_y = nonblank[7].split()
        if abs(float(fault2_y[0]) - 2000.0) > 1e-6 or abs(float(fault2_y[1]) - 2000.0) > 1e-6:
            raise AssertionError('bFaultGeometry.txt fault 2 y-line should be "2000.0 2000.0", '
                                 'got %r' % nonblank[6])

        # bStations.txt line 2 must carry ntotft=2 counts (readInputFiles.f90:152's contract).
        lines = open(os.path.join(case_dir, 'bStations.txt')).read().splitlines()
        counts = lines[1].split()
        if len(counts) != 2:
            raise AssertionError('bStations.txt line 2 should carry 2 on-fault station counts '
                                 '(ntotft=2), got %r' % lines[1])
        if int(counts[1]) != 1:
            raise AssertionError('bStations.txt fault 2 on-fault station count should be 1 '
                                 '(par.st_coor_on_fault_per_fault[1]), got %r' % counts[1])

        # on_fault_vars_input.nc: fault 2's variables are faultTag()-prefixed.
        import netCDF4
        with netCDF4.Dataset(os.path.join(case_dir, 'on_fault_vars_input.nc')) as ds:
            if 'sw_fs' not in ds.variables or 'ft2_sw_fs' not in ds.variables:
                raise AssertionError('on_fault_vars_input.nc should have both "sw_fs" (fault 1) '
                                     'and "ft2_sw_fs" (fault 2); got %r'
                                     % sorted(ds.variables.keys()))
        print('  PASS  ntotft=2 case.setup succeeds; bFaultGeometry.txt has 2 distinct boxes, '
              'bStations.txt has 2 on-fault counts, on_fault_vars_input.nc has both fault-1 and '
              'ft2_-prefixed variables')


def check_one_fault_still_works():
    """The multi-fault plumbing must not be satisfiable by breaking ntotft==1."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir, r = _make_case(tmp, 'par.ntotft = 1')
        if r.returncode != 0:
            raise AssertionError(
                'case.setup FAILED on the ordinary single-fault case '
                '(exit %d):\n%s' % (r.returncode, (r.stdout + r.stderr)[-2000:]))
        bst = os.path.join(case_dir, 'bStations.txt')
        lines = open(bst).read().splitlines()
        if len(lines) < 2 or len(lines[1].split()) != 1:
            raise AssertionError(
                'bStations.txt line 2 should carry exactly ntotft=1 on-fault '
                'station count; got %r' % (lines[1] if len(lines) > 1 else None))
        print('  PASS  ntotft=1 still works, bStations.txt line 2 has 1 count')


def check_ambiguous_multifault_still_refused():
    """ntotft>1 with NO per-fault geometry/stations/vars must still refuse --
    the scoped-down refusal (resolveFaultGeom /
    resolveOnFaultStationsPerFault / resolveOnFaultVarsPerFault), not the old
    blanket one. Guards against "fixed" by silently reusing fault 1's data
    for every fault."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir, r = _make_case(tmp, 'par.ntotft = 2\npar.nucfault = 1')
        if r.returncode == 0:
            raise AssertionError(
                'case.setup ACCEPTED par.ntotft=2 with NO par.faultgeom/onFaultVarsPerFault/'
                'st_coor_on_fault_per_fault -- it should refuse rather than silently reuse '
                'fault 1''s box/data for every fault (the exact aliasing bug multi-fault '
                'support exists to avoid)')
        out = r.stdout + r.stderr
        if 'par.faultgeom' not in out:
            raise AssertionError('refusal message does not name par.faultgeom -- not actionable:\n%s'
                                 % out[-1500:])
        print('  PASS  ntotft=2 with no per-fault geometry still refused, naming par.faultgeom')


def main():
    print('Regression guard: multi-fault case.setup contract (Row 17)')
    failures = []
    for c in (check_one_fault_still_works,
              check_two_faults_case_setup_succeeds,
              check_ambiguous_multifault_still_refused):
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_multifault_refused (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_multifault_refused')
    return 0


if __name__ == '__main__':
    sys.exit(main())
