#! /usr/bin/env python3
"""
Regression guard: test.multifault2, the row-17 two-fault case, on the
python-jax backend.

BLOCKER 3 (audited, 2026-10-01): the `60cb405` port mission's claim --
"test.multifault2 on python-jax matches the Fortran-verified reference to
near machine epsilon" -- had NOTHING in the repo testing it.
test.multifault2 is not registered in testNameList.py/matrix.py, and no
unit/regression test ran ntotft>1 through the python port at all before this
file. test_multifault_two_fault_smoke.py (the Fortran smoke test this module
sits beside) says outright, in its own docstring, that the python-jax
backend was ntotft==1-only "as of this writing" -- true when THAT docstring
was written, stale after `60cb405` actually ported ntotft>1 support
(src/python/eqdyna/readInputFiles.py, eqdyna3d.py, meshgen.py all handle
ntotft>1 now -- see readInputFiles.py:299-311, faultTag()).

THIS CLOSES THE GAP: runs case_input/test.multifault2 on python-jax (serial;
the standalone solver refuses any decomposition with npx*npy*npz > 1) and
compares its canonical frt output against the SAME existing, already
Fortran-verified committed reference
(test.reference.results/test.multifault2/frt.canonical.txt, 3782 rows) the
Fortran smoke test uses -- no new reference is created, so this does not
touch the reference-freeze rule (rule 7).

Mirrors, by IMPORTING (not duplicating), that module's per-column tolerance
structure (COL_BOUND / _per_column_bound_report: TOL_TIGHT=1e-9 for the
time/slip/slip-rate/velocity/state columns, TOL_PA=1e7 for the three
Pa-scale traction columns + theta_pc) and its two-faults-present-and-distinct
structural check (_check_two_faults_present_and_distinct), and the real
run_e2e.run_standalone (the actual `python3 -m eqdyna <case> --backend jax`
entry point every e2e python-jax cell goes through) to launch the solver --
rule 10a: behaviour, not a reimplementation.

par.term is left at case_input/test.multifault2's own COMMITTED value
(1.0 s), not overridden to matrix.GATE_TERM_S (5.0 s): this case is
deliberately NOT in testNameList.py/matrix.py (see that module's own
docstring), and the committed reference was generated at term=1.0s.
run_e2e.make_serial_case() is NOT reused here for exactly this reason -- it
also calls apply_term_override, which would silently compare this case
against the wrong physical endpoint. The serial-decomposition override is
duplicated here (3 lines) rather than pulling in that coupling.

MEASURED (this work, serial python-jax vs the committed reference): max|diff|
= 1.061e-08 Pa (column 13, tdip -- comfortably inside TOL_PA=1e7) and
5.811e-16 (column 4, tight group -- inside TOL_TIGHT=1e-9) -- the PR's parity
claim actually holds, not just asserted; this test is what makes a future
regression in the python port's multi-fault handling fail CI instead of
going unnoticed (there was nothing here to catch it before this file).

Not cheap (rule 9) -- builds a case and runs the actual jax solver -- but
bounded: dx=500m, term=1.0s (23 steps), serial, well under a minute.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

from testsys import frt_canonical  # noqa: E402
import run_e2e  # noqa: E402  (the real run_standalone -- not reimplemented)
import test_multifault_two_fault_smoke as fortran_smoke  # noqa: E402  (shared tolerance structure)

CASE_NAME = fortran_smoke.CASE_NAME
REFERENCE = fortran_smoke.REFERENCE


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _build_serial_case(tmp):
    """create.newcase + force a serial decomposition + case.setup.

    Deliberately NOT run_e2e.make_serial_case -- see module docstring: that
    helper also pins par.term to matrix.GATE_TERM_S (5.0s), which would
    compare this case's python-jax run against the committed reference at
    the WRONG term (the reference was generated at this case's own
    committed par.term=1.0s, and this case is deliberately not in
    matrix.py/testNameList.py, so matrix.GATE_TERM_S has no claim over it).
    """
    case_dir = os.path.join(tmp, 'case_jax')
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial by this test (the standalone '
                       'solver is serial-only); par.term is left untouched '
                       '-- see module docstring\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n')
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup failed (%d):\n%s' % (r.returncode, r.stdout + r.stderr))
    return case_dir, env


def _require_jax():
    try:
        import jax  # noqa: F401
    except ImportError as exc:
        raise AssertionError(
            'jax is not importable in this interpreter (%s) -- the '
            'python-jax backend cannot be exercised. Install it before '
            'running this guard.' % exc)


_CANONICAL = []  # one solver run per process, shared by both checks below


def _canonical_jax():
    if not _CANONICAL:
        _require_jax()
        with tempfile.TemporaryDirectory() as tmp:
            case_dir, env = _build_serial_case(tmp)
            run_e2e.run_standalone(case_dir, 'python-jax', device='cpu', env=env)
            _CANONICAL.append(frt_canonical.canonical_from_case(case_dir))
    return _CANONICAL[0]


def check_python_jax_matches_committed_reference():
    fortran_smoke._check_matches_committed_reference(_canonical_jax(), 'python-jax serial')


def check_python_jax_two_faults_present_and_distinct():
    fortran_smoke._check_two_faults_present_and_distinct(_canonical_jax(), 'python-jax serial')


def check_mutation_rupture_time_error_is_caught_by_tight_bound():
    """Same proof-of-load-bearing-ness as the Fortran smoke test's own
    check, against the SAME reference file -- confirms the imported
    tolerance structure is not vacuous for this module either."""
    fortran_smoke.check_mutation_rupture_time_error_is_caught_by_tight_bound()


def main():
    print('Regression guard: test.multifault2 two-fault case on python-jax (BLOCKER 3)')
    failures = []
    for c in (check_mutation_rupture_time_error_is_caught_by_tight_bound,
              check_python_jax_two_faults_present_and_distinct,
              check_python_jax_matches_committed_reference):
        try:
            c()
            print('  PASS  %s' % c.__name__)
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_multifault_two_fault_smoke_jax (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_multifault_two_fault_smoke_jax')
    return 0


if __name__ == '__main__':
    sys.exit(main())
