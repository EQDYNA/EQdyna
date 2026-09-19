#! /usr/bin/env python3
"""
Regression guard: ntotft > 1 must be REFUSED with a named reason (rules 2, 10).

THE CAUSE, which is a WRITER/READER CONTRACT MISMATCH, not a parse bug:

    scripts/case.setup   writes ONE scalar (par.n_on_fault) on line 2
    readInputFiles.f90:152 reads ntotft integers, `(nonfs(i), i = 1, ntotft)`

They agree only at ntotft == 1. At ntotft == 2 the list-directed read runs off
the short line into the on-fault station coordinates and hits `0.0` where it
wants an integer -- a raw gfortran I/O error naming neither ntotft nor the
cause. `meshgen.f90:790`'s deliberate `abortRun(ERR_MESH_MULTIFAULT_MSNODE)`
(exit 45) for this case is UNREACHABLE through the normal path: the input
read above fails first.

Writing ntotft copies of the one count would silence the crash and produce
wrong station bookkeeping: `nonfs(i)` is the on-fault station count FOR FAULT
i, and `par.st_coor_on_fault` is a flat list with no per-fault assignment.
That assignment is a design decision (pathway_forward item 17).

WHAT THIS PINS. `case.setup` refuses ntotft > 1 with a non-zero exit and a
message that names ntotft, the file, the reader, and the item. It also pins
that ntotft == 1 still WORKS, so the guard cannot be "fixed" by refusing
everything.

Cheap (rule 9): two `case.setup` invocations, no build, no MPI, no simulation.
Exits non-zero on any failure.
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CASE_NAME = 'test.tpv8'

# The message must name each of these, or it is not actionable.
REQUIRED_IN_MESSAGE = ('ntotft', 'bStations.txt', 'readInputFiles.f90:152',
                       'item 17')


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
    return r


def check_two_faults_are_refused():
    with tempfile.TemporaryDirectory() as tmp:
        r = _make_case(tmp, 'par.ntotft = 2\npar.nucfault = 1')
        out = r.stdout + r.stderr
        if r.returncode == 0:
            raise AssertionError(
                'case.setup ACCEPTED par.ntotft = 2 (exit 0). It must refuse: '
                'bStations.txt would be written with one on-fault station '
                'count where readInputFiles.f90:152 reads ntotft of them, and '
                'the run dies later with an unactionable gfortran I/O error.')
        missing = [t for t in REQUIRED_IN_MESSAGE if t not in out]
        if missing:
            raise AssertionError(
                'case.setup refused ntotft=2 (exit %d) but its message does '
                'not name %r -- a refusal the reader cannot act on is only '
                'half a guard. Message:\n%s'
                % (r.returncode, missing, out[-1500:]))
        print('  PASS  ntotft=2 refused (exit %d), message names %s'
              % (r.returncode, ', '.join(REQUIRED_IN_MESSAGE)))


def check_one_fault_still_works():
    """The guard must not be satisfiable by refusing everything."""
    with tempfile.TemporaryDirectory() as tmp:
        r = _make_case(tmp, 'par.ntotft = 1')
        if r.returncode != 0:
            raise AssertionError(
                'case.setup FAILED on the ordinary single-fault case '
                '(exit %d) -- the multi-fault refusal is too broad:\n%s'
                % (r.returncode, (r.stdout + r.stderr)[-2000:]))
        bst = os.path.join(tmp, CASE_NAME, 'bStations.txt')
        if not os.path.isfile(bst):
            raise AssertionError('single-fault case.setup wrote no bStations.txt')
        lines = open(bst).read().splitlines()
        if len(lines) < 2 or len(lines[1].split()) != 1:
            raise AssertionError(
                'bStations.txt line 2 should carry exactly ntotft=1 on-fault '
                'station count; got %r' % (lines[1] if len(lines) > 1 else None))
        print('  PASS  ntotft=1 still works, bStations.txt line 2 has 1 count')


def check_the_fortran_reader_still_expects_ntotft():
    """If the reader is ever changed to read one value, this guard is moot and
    should be revisited rather than left asserting a stale contract."""
    src = os.path.join(ROOT, 'src', 'fortran', 'readInputFiles.f90')
    text = open(src, errors='replace').read()
    if 'nonfs(i), i = 1, ntotft' not in text.replace('  ', ' '):
        raise AssertionError(
            '%s no longer reads `(nonfs(i), i = 1, ntotft)`. The contract this '
            'guard pins has changed -- re-derive it instead of deleting this '
            'check.' % src)
    print('  PASS  readInputFiles.f90 still reads ntotft on-fault counts')


def main():
    print('Regression guard: ntotft > 1 must be refused with a named reason')
    failures = []
    for c in (check_the_fortran_reader_still_expects_ntotft,
              check_two_faults_are_refused,
              check_one_fault_still_works):
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
