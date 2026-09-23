#! /usr/bin/env python3
"""
Regression guard: ONE term, everywhere (2026-09-23 owner decision,
superseding a same-day earlier two-term design) (rules 2, 6).

THE DEFECT SHAPE THIS GUARDS AGAINST, each pinned below:

  1. a `--term`-shaped flag reappearing on run_e2e.py's argparse, or
     `--release` disappearing from it -- reopening the fork this change
     closed.
  2. matrix.py regaining CASE_FULL_TERM_S, or compare.py regaining
     GATE_TERM_NAME / canonical_reference_name's two-way choice, or any
     public comparison entry point regaining a `term` parameter.
  3. a case ending up with two committed reference files on disk
     (frt.canonical.txt plus a retired frt.canonical.term5.txt) -- exactly
     what item 4 of the 2026-09-23 landing retired for
     test.tpv29/test.tpv36/test.tpv37.
  4. run.py's `release` tier drifting off `--release` (back onto a `--term`
     flag that no longer exists, or silently back to running nothing wider
     than the everyday sweep), or CI's smoke invocation drifting onto either
     flag.

Cheap (rule 9): imports, one workflow text parse, a handful of introspection/
signature checks, one throwaway-tempdir behavioural check, one monkeypatched
subprocess call. No solver, no MPI, no I/O beyond a few bytes in tempdirs.
Well under 1 s. Exits non-zero on any failure.
"""
import contextlib
import glob
import importlib.util
import inspect
import io
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')
RUN_PY = os.path.join(ROOT, 'testsys', 'run.py')
E2E_PY = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py')

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import compare, matrix  # noqa: E402


def _load_module(name, path):
    """A fresh module object from `path` (it is not a package) -- never the
    cached sys.modules entry, so each check gets an independent instance."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_matrix_has_no_second_term_table():
    if hasattr(matrix, 'CASE_FULL_TERM_S'):
        raise AssertionError(
            'matrix.CASE_FULL_TERM_S still exists -- the per-case "full" '
            'term table was supposed to be retired along with the --term '
            'flag (2026-09-23 one-term decision)')
    if not hasattr(matrix, 'GATE_TERM_S'):
        raise AssertionError(
            'matrix.GATE_TERM_S is gone -- this is the ONE term every cell '
            'runs at; it must still exist')
    print('  PASS  matrix.py carries GATE_TERM_S (%.1fs) and no '
          'CASE_FULL_TERM_S' % matrix.GATE_TERM_S)


def check_compare_has_no_second_reference_name():
    if hasattr(compare, 'GATE_TERM_NAME'):
        raise AssertionError(
            'compare.GATE_TERM_NAME still exists -- the gate-term reference '
            'name was supposed to be retired: there is only ONE reference '
            'file per case now (compare.CANONICAL_NAME)')
    if hasattr(compare, 'canonical_reference_name'):
        raise AssertionError(
            'compare.canonical_reference_name still exists -- its two-way '
            '"gate" vs "full" choice was supposed to be removed along with '
            'the term axis')
    checked = []
    for fn in (compare.load_reference, compare.load_coordinate_aligned,
              compare.compare_frt, compare.compare_cell, compare.drv_a6_gate):
        params = inspect.signature(fn).parameters
        if 'term' in params:
            raise AssertionError(
                'compare.%s still accepts a `term` parameter -- there is no '
                'term axis left to select' % fn.__name__)
        checked.append(fn.__name__)
    print('  PASS  compare.py carries no GATE_TERM_NAME, no '
          'canonical_reference_name, and no `term` parameter on %s'
          % ', '.join(checked))


def _canonical_variants(ref_root, case):
    """Every frt.canonical*.txt basename under ref_root/case, sorted -- the
    same glob the real check below walks against the committed reference
    tree; factored out so it can also be driven against a synthetic tmp tree
    (see the mutation check)."""
    ref_dir = os.path.join(ref_root, case)
    return sorted(os.path.basename(p)
                 for p in glob.glob(os.path.join(ref_dir, 'frt.canonical*.txt')))


def check_no_case_has_two_reference_files():
    bad = []
    for case in matrix.CASES:
        variants = _canonical_variants(compare.REFERENCE_ROOT, case)
        if variants != [compare.CANONICAL_NAME]:
            bad.append('%s: %r (want exactly [%r])'
                      % (case, variants, compare.CANONICAL_NAME))
    if bad:
        raise AssertionError(
            'case(s) with other than exactly one frt.canonical*.txt '
            'reference: %s' % '; '.join(bad))
    print('  PASS  every one of %d case(s) has exactly one reference file '
          '(%r) -- no second frt.canonical.term5.txt anywhere'
          % (len(matrix.CASES), compare.CANONICAL_NAME))


def check_mutation_a_planted_second_reference_file_is_caught():
    """_canonical_variants is what the real check above walks; prove it on a
    SYNTHETIC tmp reference tree with (a) exactly one canonical file --
    passes -- and (b) a planted second one, the retired term5 shape --
    fails -- without ever touching the real committed reference tree."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = os.path.join(tmp, 'test.probe')
        os.makedirs(case_dir)
        open(os.path.join(case_dir, compare.CANONICAL_NAME), 'w').close()
        variants_one = _canonical_variants(tmp, 'test.probe')
        if variants_one != [compare.CANONICAL_NAME]:
            raise AssertionError(
                'one committed-shaped reference file was not recognised as '
                'exactly one: %r' % variants_one)
        open(os.path.join(case_dir, 'frt.canonical.term5.txt'), 'w').close()
        variants_two = _canonical_variants(tmp, 'test.probe')
        if variants_two == [compare.CANONICAL_NAME]:
            raise AssertionError(
                'planting a second reference file (frt.canonical.term5.txt) '
                'was NOT detected -- the real check would not catch a case '
                'with two reference files')
    print('  PASS  a planted second reference file is detected (mutation-'
          'tested both ways: one file passes, a planted second one fails)')


def check_apply_term_override_takes_no_term_argument_and_is_unconditional():
    """run_e2e.apply_term_override(case_name, case_dir) -- no `term`
    parameter (signature check) -- and always appends GATE_TERM_S as the
    LAST par.term assignment, even over a pre-existing one (behavioural
    check against a throwaway tempdir; mutated in-line below)."""
    e2e = _load_module('run_e2e_for_term_axis_override', E2E_PY)
    params = list(inspect.signature(e2e.apply_term_override).parameters)
    if params != ['case_name', 'case_dir']:
        raise AssertionError(
            'apply_term_override signature is %r, want exactly '
            '(case_name, case_dir) -- no `term` parameter should exist' % params)
    with tempfile.TemporaryDirectory() as d:
        params_py = os.path.join(d, 'user_defined_params.py')
        with open(params_py, 'w') as f:
            f.write('par = object()\npar.term = 999.0  # simulated stale value\n')
        e2e.apply_term_override('test.probe', d)
        text = open(params_py).read()
    term_lines = [l for l in text.splitlines() if l.startswith('par.term')]
    if term_lines[-1] != ('par.term = %r' % matrix.GATE_TERM_S):
        raise AssertionError(
            'apply_term_override did not leave par.term = %r as the LAST '
            'assignment (case.setup reads the last one) -- got %r'
            % (matrix.GATE_TERM_S, term_lines))
    print('  PASS  apply_term_override(case_name, case_dir) takes no `term` '
          'argument and unconditionally appends par.term = %r as the LAST '
          'assignment, even over a pre-existing stale one'
          % matrix.GATE_TERM_S)


def check_run_e2e_help_lists_release_and_not_term():
    e2e = _load_module('run_e2e_for_term_axis_help', E2E_PY)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            e2e.main(['--help'])
    except SystemExit as exc:
        if exc.code != 0:
            raise AssertionError('run_e2e.py --help exited %r, want 0' % exc.code)
    else:
        raise AssertionError('run_e2e.py --help did not raise SystemExit')
    help_text = buf.getvalue()
    if '--release' not in help_text:
        raise AssertionError(
            'run_e2e.py --help does not mention --release -- the release '
            'selector this change added')
    if '--term' in help_text:
        raise AssertionError(
            'run_e2e.py --help still mentions --term -- the flag this '
            'change removed reappeared: %r' % help_text)
    print('  PASS  run_e2e.py --help lists --release and no --term')


def check_run_e2e_refuses_a_term_flag():
    e2e = _load_module('run_e2e_for_term_axis_refuse', E2E_PY)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stderr(buf):
            e2e.main(['--term', 'full'])
    except SystemExit as exc:
        if exc.code == 0:
            raise AssertionError(
                '`run_e2e.py --term full` was ACCEPTED (exit 0) -- the '
                '--term flag should not exist at all')
    else:
        raise AssertionError(
            '`run_e2e.py --term full` did not raise SystemExit -- argparse '
            'should refuse an unrecognised argument')
    print('  PASS  run_e2e.py refuses --term as an unrecognised argument '
          '(argparse exit-status message: %r)' % (buf.getvalue().strip()[-80:]))


def check_ci_smoke_invocation_passes_neither_term_nor_release():
    text = open(WORKFLOW, errors='replace').read()
    lines = [raw.split('#', 1)[0] for raw in text.splitlines()]
    ci_lines = [l for l in lines if 'run_e2e.py' in l and '--ci' in l]
    if not ci_lines:
        raise AssertionError('no `run_e2e.py --ci` invocation found in %s' % WORKFLOW)
    bad = [l.strip() for l in ci_lines if '--term' in l or '--release' in l]
    if bad:
        raise AssertionError(
            'CI smoke invocation(s) pass --term or --release, which would '
            'drift it off the portability-smoke selection: %r' % bad)
    print('  PASS  %d CI smoke invocation(s) pass neither --term nor --release'
          % len(ci_lines))


def _load_run_module():
    return _load_module('run_mod_for_term_axis', RUN_PY)


def check_release_tier_requests_the_release_selection():
    """run.py's `release` tier must construct a `run_e2e.py --release`
    command -- checked BEHAVIOURALLY (rule 10a): subprocess.call is
    monkeypatched to capture the argv this tier builds, and the real
    run_e2e.py is never invoked."""
    mod = _load_run_module()
    if 'release' not in mod.RUNNERS:
        raise AssertionError("run.py has no 'release' tier registered in RUNNERS")
    captured = []
    orig_call = mod.subprocess.call

    def fake_call(cmd, *a, **k):
        captured.append(list(cmd))
        return 0

    mod.subprocess.call = fake_call
    try:
        rc = mod.RUNNERS['release']()
    finally:
        mod.subprocess.call = orig_call
    if rc != 0:
        raise AssertionError('run_release() returned %r with a stubbed '
                             'subprocess.call (expected 0)' % rc)
    if len(captured) != 1:
        raise AssertionError(
            'expected the release tier to make exactly one subprocess call, '
            'got %d: %r' % (len(captured), captured))
    cmd = captured[0]
    if 'run_e2e.py' not in ' '.join(cmd):
        raise AssertionError('release tier did not invoke run_e2e.py: %r' % cmd)
    if '--release' not in cmd:
        raise AssertionError(
            'release tier must pass --release to run_e2e.py, got argv %r' % cmd)
    if '--term' in cmd:
        raise AssertionError(
            'release tier passed a --term flag, which no longer exists: %r' % cmd)
    if '--ci' in cmd:
        raise AssertionError(
            'release tier passed --ci, which would restrict it to the CI '
            'smoke cell list instead of the full runnable matrix: %r' % cmd)
    print('  PASS  release tier invokes %r' % cmd)


def main():
    print('Regression guard: ONE term everywhere (no --term flag, no second '
          'reference file, no per-case full-term table)')
    checks = [check_matrix_has_no_second_term_table,
              check_compare_has_no_second_reference_name,
              check_no_case_has_two_reference_files,
              check_mutation_a_planted_second_reference_file_is_caught,
              check_apply_term_override_takes_no_term_argument_and_is_unconditional,
              check_run_e2e_help_lists_release_and_not_term,
              check_run_e2e_refuses_a_term_flag,
              check_ci_smoke_invocation_passes_neither_term_nor_release,
              check_release_tier_requests_the_release_selection]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append(str(e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_term_axis (%d check(s))' % len(failures))
        return 1
    print('\nSUCCESS test_term_axis')
    return 0


if __name__ == '__main__':
    sys.exit(main())
