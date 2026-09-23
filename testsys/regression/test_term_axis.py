#! /usr/bin/env python3
"""
Regression guard: the TERM axis (2026-09-23 owner-approved test-methodology
change) fails closed, never silently substitutes, and CI/release stay on
their own sides of it (rules 2, 6).

THE DEFECT SHAPE THIS GUARDS AGAINST. `term` is an axis of a cell exactly
like `backend`: 'gate' runs every case at matrix.GATE_TERM_S regardless of
its own committed par.term; 'full' runs each case at its own par.term. A
case whose full term is not matrix.GATE_TERM_S (test.tpv29, test.tpv36,
test.tpv37) needs a SECOND reference file for the gate term
(compare.GATE_TERM_NAME). Three ways that could silently go wrong, each
pinned below:

  1. a gate-term cell for one of those three cases compares against the
     FULL-term reference instead of failing -- a 5 s run would look "wrong"
     against a 20 s/6 s reference for reasons that have nothing to do with a
     real regression, or worse, could be made to look right by chance.
  2. --term full stops selecting the one full-length reference every other
     guard and every committed frt.canonical.txt assumes.
  3. CI's smoke job (test.tpv8 only) drifts to running the full-term sweep,
     or the release tier drifts to running the cheap gate-term one -- either
     direction defeats the point of splitting them (portability smoke vs.
     physics-coverage release).

Cheap (rule 9): imports + one workflow text parse + one monkeypatched
subprocess call, well under 1 s. Exits non-zero on any failure.
"""
import importlib.util
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORKFLOW = os.path.join(ROOT, '.github', 'workflows', 'test.yml')
RUN_PY = os.path.join(ROOT, 'testsys', 'run.py')

if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
from testsys import compare, matrix  # noqa: E402


def _cases_needing_a_second_reference():
    """Cases whose full term != matrix.GATE_TERM_S -- the ones that MUST
    fail closed at term='gate' until compare.GATE_TERM_NAME is committed for
    them (rule 7: the owner generates it, this guard never does)."""
    return [c for c in matrix.CASES
            if matrix.CASE_FULL_TERM_S[c] != matrix.GATE_TERM_S]


def check_gate_term_fails_closed_without_a_term5_reference():
    needing = _cases_needing_a_second_reference()
    if not needing:
        raise AssertionError(
            'no case in matrix.CASES has a full term != matrix.GATE_TERM_S -- '
            'this guard has nothing to exercise, which means either the term '
            'axis or matrix.CASE_FULL_TERM_S regressed silently')
    checked = 0
    for case in needing:
        name = compare.canonical_reference_name(case, 'gate')
        if name != compare.GATE_TERM_NAME:
            raise AssertionError(
                '%s (full term %.1fs != gate term %.1fs) should select the '
                'gate-term reference name %r at term=gate, selected %r '
                'instead' % (case, matrix.CASE_FULL_TERM_S[case],
                            matrix.GATE_TERM_S, compare.GATE_TERM_NAME, name))
        ref_dir = os.path.join(compare.REFERENCE_ROOT, case)
        term5_path = os.path.join(ref_dir, compare.GATE_TERM_NAME)
        if os.path.isfile(term5_path):
            # The owner has generated this one already (rule 7) -- nothing to
            # prove FAILS CLOSED here anymore; skip rather than assert a
            # FileNotFoundError that would no longer be true.
            continue
        try:
            compare.reference_path(case, name)
        except FileNotFoundError as exc:
            if compare.GATE_TERM_NAME not in str(exc) or case not in str(exc):
                raise AssertionError(
                    '%s: FileNotFoundError did not name the missing file/case '
                    '(rule 2 -- a failure must say what it could not check): '
                    '%s' % (case, exc))
            checked += 1
        else:
            raise AssertionError(
                '%s: reference_path(%r) did not raise even though %s does '
                'not exist -- a gate-term cell for this case would silently '
                'compare against something other than a missing-file error '
                '(possibly the full-term reference), which is exactly the '
                'silent fallback rule 2 forbids' % (case, name, term5_path))
    print('  PASS  %d/%d case(s) needing a second reference fail closed '
          '(no silent fallback to the full-term reference): %s'
          % (checked, len(needing), ', '.join(needing)))


def check_full_term_selects_the_full_reference():
    for case in matrix.CASES:
        name = compare.canonical_reference_name(case, 'full')
        if name != compare.CANONICAL_NAME:
            raise AssertionError(
                '%s at term=full selected %r, expected the one full-length '
                'reference %r' % (case, name, compare.CANONICAL_NAME))
        # And that file must actually exist -- matrix.py's own import-time
        # consistency block already enforces this for every case, so this
        # is restating the guarantee at the call site this guard cares
        # about, not inventing a new one.
        compare.reference_path(case, name)
    print('  PASS  term=full selects %r for all %d case(s), and it exists'
          % (compare.CANONICAL_NAME, len(matrix.CASES)))


def check_gate_term_matches_full_term_reference_when_terms_are_equal():
    """A case whose own full term already equals matrix.GATE_TERM_S (e.g.
    test.tpv8) must reuse CANONICAL_NAME at term=gate -- no second file is
    needed or expected for it."""
    equal_cases = [c for c in matrix.CASES
                  if matrix.CASE_FULL_TERM_S[c] == matrix.GATE_TERM_S]
    if not equal_cases:
        raise AssertionError(
            'no case has full term == matrix.GATE_TERM_S -- test.tpv8 (CI\'s '
            'smoke cell) is expected to be one of these')
    for case in equal_cases:
        name = compare.canonical_reference_name(case, 'gate')
        if name != compare.CANONICAL_NAME:
            raise AssertionError(
                '%s: full term already equals GATE_TERM_S, so term=gate '
                'should reuse %r, selected %r instead'
                % (case, compare.CANONICAL_NAME, name))
    print('  PASS  %d case(s) with full term == GATE_TERM_S reuse %r at '
          'term=gate (no second reference needed): %s'
          % (len(equal_cases), compare.CANONICAL_NAME, ', '.join(equal_cases)))


def check_ci_smoke_invocations_never_request_full_term():
    text = open(WORKFLOW, errors='replace').read()
    lines = [raw.split('#', 1)[0] for raw in text.splitlines()]
    ci_lines = [l for l in lines if 'run_e2e.py' in l and '--ci' in l]
    if not ci_lines:
        raise AssertionError('no `run_e2e.py --ci` invocation found in %s' % WORKFLOW)
    bad = [l.strip() for l in ci_lines if '--term' in l]
    if bad:
        raise AssertionError(
            'CI smoke invocation(s) pass --term explicitly, which risks '
            'drifting to full: %r (they should rely on the gate default)'
            % bad)
    print('  PASS  %d CI smoke invocation(s) pass no --term (default gate)'
          % len(ci_lines))


def _load_run_module():
    spec = importlib.util.spec_from_file_location('run_mod_for_term_axis', RUN_PY)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def check_release_tier_requests_full_term():
    """run.py's `release` tier must construct a `run_e2e.py --term full`
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
    if '--term' not in cmd or 'full' not in cmd:
        raise AssertionError(
            'release tier must pass --term full to run_e2e.py, got argv %r'
            % cmd)
    if '--ci' in cmd:
        raise AssertionError(
            'release tier passed --ci, which would restrict it to the CI '
            'smoke cell list instead of the full runnable matrix: %r' % cmd)
    print('  PASS  release tier invokes %r' % cmd)


def check_full_term_table_matches_the_compsets():
    """matrix.CASE_FULL_TERM_S is bookkeeping ABOUT each compset's own
    par.term, which is what `--term full` actually runs -- so a compset whose
    par.term moves without the table moving would make `--term gate` pick
    the wrong reference file. Read the LAST `par.term =` assignment in
    case_input/<case>/user_defined_params.py, or, for a case that builds par
    elsewhere (tpv36/37: tpv36_37_common.buildParams), in the one other .py
    of that compset that assigns it. (wei-lin, 2026-09-23.)"""
    import glob, re
    pat = re.compile(r'^\s*par\.term\s*=\s*([0-9.eE+-]+)', re.M)
    bad = []
    for case, want in sorted(matrix.CASE_FULL_TERM_S.items()):
        d = os.path.join(ROOT, 'case_input', case)
        hits = pat.findall(open(os.path.join(d, 'user_defined_params.py')).read())
        if not hits:
            for f in sorted(glob.glob(os.path.join(d, '*.py'))):
                hits = hits or pat.findall(open(f).read())
        if not hits:
            bad.append('%s: no par.term assignment found under %s' % (case, d))
        elif float(hits[-1]) != want:
            bad.append('%s: compset par.term %s != matrix.CASE_FULL_TERM_S %s'
                       % (case, hits[-1], want))
    if bad:
        raise AssertionError('; '.join(bad))
    print('  PASS  CASE_FULL_TERM_S matches the compsets\' own par.term for '
          '%d case(s)' % len(matrix.CASE_FULL_TERM_S))


def main():
    print('Regression guard: the TERM axis fails closed and CI/release stay separated')
    checks = [check_full_term_table_matches_the_compsets,
              check_gate_term_fails_closed_without_a_term5_reference,
              check_full_term_selects_the_full_reference,
              check_gate_term_matches_full_term_reference_when_terms_are_equal,
              check_ci_smoke_invocations_never_request_full_term,
              check_release_tier_requests_full_term]
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
