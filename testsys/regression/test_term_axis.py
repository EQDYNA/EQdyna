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
    parameter (signature check) -- and always appends matrix.gate_term_for
    (case_name) as the LAST par.term assignment, even over a pre-existing
    one (behavioural check against a throwaway tempdir; mutated in-line
    below). Checked against 'test.probe', a name with NO entry in
    matrix.CASE_TERM_OVERRIDE, so this check is specifically about the
    DEFAULT path (matrix.GATE_TERM_S) -- the named exception for
    test.tpv22/test.tpv23 is its own check below."""
    e2e = _load_module('run_e2e_for_term_axis_override', E2E_PY)
    params = list(inspect.signature(e2e.apply_term_override).parameters)
    if params != ['case_name', 'case_dir']:
        raise AssertionError(
            'apply_term_override signature is %r, want exactly '
            '(case_name, case_dir) -- no `term` parameter should exist' % params)
    if 'test.probe' in matrix.CASE_TERM_OVERRIDE:
        raise AssertionError(
            "'test.probe' is in matrix.CASE_TERM_OVERRIDE -- this check needs "
            'a name with NO override entry to test the DEFAULT path')
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
          'assignment for an unregistered case name, even over a '
          'pre-existing stale one' % matrix.GATE_TERM_S)


def _assert_term_override_set_is_exactly(names, want):
    if set(names) != want:
        raise AssertionError('term-override case set is %r, want exactly %r'
                             % (set(names), want))


_WANT_TERM_OVERRIDE_CASES = frozenset({'test.tpv22', 'test.tpv23'})


def check_case_term_override_is_exactly_tpv22_tpv23():
    """matrix.CASE_TERM_OVERRIDE must name EXACTLY {test.tpv22, test.tpv23},
    each at 15.0s (TPV22_23_Description_v08.pdf Part 5's own full duration)
    -- never a third case, never a different number. This is the named,
    scoped exception rule 24/this module exists to bound: everything else
    in the sweep still runs at the one GATE_TERM_S."""
    _assert_term_override_set_is_exactly(matrix.CASE_TERM_OVERRIDE,
                                         _WANT_TERM_OVERRIDE_CASES)
    bad = {c: v for c, v in matrix.CASE_TERM_OVERRIDE.items() if v != 15.0}
    if bad:
        raise AssertionError(
            'matrix.CASE_TERM_OVERRIDE must carry 15.0s for every entry '
            '(the TPV22/23 spec duration), got %r' % bad)
    print('  PASS  matrix.CASE_TERM_OVERRIDE is exactly %s, each at 15.0s'
          % sorted(_WANT_TERM_OVERRIDE_CASES))


def check_mutation_a_third_case_in_term_override_is_caught():
    """Proves the exactness check above is load-bearing (rule 10a): a THIRD
    case name planted into a COPY of the real override set must make
    `_assert_term_override_set_is_exactly` raise -- never silently widen.
    Mutation-tested both ways, like check_mutation_a_planted_second_
    reference_file_is_caught above: the real set passes, the planted one
    fails."""
    _assert_term_override_set_is_exactly(matrix.CASE_TERM_OVERRIDE,
                                         _WANT_TERM_OVERRIDE_CASES)  # real set: passes
    mutated = dict(matrix.CASE_TERM_OVERRIDE)
    mutated['test.tpv8'] = 5.0  # a third case name, must be refused
    try:
        _assert_term_override_set_is_exactly(mutated, _WANT_TERM_OVERRIDE_CASES)
    except AssertionError:
        pass
    else:
        raise AssertionError(
            'planting a 3rd case (test.tpv8) into the term-override set was '
            'NOT detected -- the exactness guard is not load-bearing')
    print('  PASS  a planted 3rd case in the term-override set is detected '
          '(mutation-tested: real set passes, planted 3rd-case set fails)')


def check_apply_term_override_honours_the_named_exception():
    """The other half of the behavioural check above: for test.tpv22 (named
    in matrix.CASE_TERM_OVERRIDE), apply_term_override must append its OWN
    15.0s term, not matrix.GATE_TERM_S -- and test.tpv8 (NOT named) must
    still get matrix.GATE_TERM_S, in the same process, proving the
    exception is scoped by name and does not leak onto an ordinary case."""
    e2e = _load_module('run_e2e_for_term_axis_named_override', E2E_PY)
    for case_name, want_term in (('test.tpv22', matrix.CASE_TERM_OVERRIDE['test.tpv22']),
                                 ('test.tpv8', matrix.GATE_TERM_S)):
        with tempfile.TemporaryDirectory() as d:
            params_py = os.path.join(d, 'user_defined_params.py')
            with open(params_py, 'w') as f:
                f.write('par = object()\npar.term = 999.0  # simulated stale value\n')
            e2e.apply_term_override(case_name, d)
            text = open(params_py).read()
        term_lines = [l for l in text.splitlines() if l.startswith('par.term')]
        if term_lines[-1] != ('par.term = %r' % want_term):
            raise AssertionError(
                'apply_term_override(%r, ...) left %r as the last par.term '
                'assignment, want %r' % (case_name, term_lines[-1],
                                         'par.term = %r' % want_term))
    print('  PASS  apply_term_override honours the named exception '
          '(test.tpv22 -> %gs) without leaking it onto test.tpv8 (-> %gs)'
          % (matrix.CASE_TERM_OVERRIDE['test.tpv22'], matrix.GATE_TERM_S))


# --------------------------------------------------------------------------
# matrix.RELEASE_ONLY -- reintroduced 2026-10-02 (item 17 section B;
# PROJECT_RULES rule 24's own contingency clause: "If a future case is held
# out of the everyday run for cost, RELEASE_ONLY (or an equivalent) is
# reintroduced deliberately, not left implied"). This mirrors the
# CASE_TERM_OVERRIDE checks immediately above: an exactness guard, a
# mutation test proving it is load-bearing, and a behavioural check that the
# everyday selection excludes it while the release selection includes it.
# --------------------------------------------------------------------------
_WANT_RELEASE_ONLY_CASES = frozenset({'test.tpv22', 'test.tpv23'})


def _assert_release_only_set_is_exactly(names, want):
    if set(names) != want:
        raise AssertionError('RELEASE_ONLY case set is %r, want exactly %r'
                             % (set(names), want))


def check_release_only_is_exactly_tpv22_tpv23():
    """matrix.RELEASE_ONLY must name EXACTLY {test.tpv22, test.tpv23} --
    never a third case. A case held out of the everyday sweep for cost is a
    rare, deliberate, named exception (rule 24), not a growing list -- this
    is the mechanical fence around it."""
    _assert_release_only_set_is_exactly(matrix.RELEASE_ONLY,
                                        _WANT_RELEASE_ONLY_CASES)
    print('  PASS  matrix.RELEASE_ONLY is exactly %s'
          % sorted(_WANT_RELEASE_ONLY_CASES))


def check_mutation_a_third_case_in_release_only_is_caught():
    """Mutation-tested (rule 10a), same shape as the CASE_TERM_OVERRIDE
    mutation test above: a THIRD case name planted into a COPY of the real
    RELEASE_ONLY set must make `_assert_release_only_set_is_exactly` raise.
    Run live below; the actual red output is quoted in the dispatch report."""
    _assert_release_only_set_is_exactly(matrix.RELEASE_ONLY,
                                        _WANT_RELEASE_ONLY_CASES)  # real set: passes
    mutated = dict(matrix.RELEASE_ONLY)
    mutated['test.tpv36'] = 'planted -- must be refused'  # a third case, must be refused
    try:
        _assert_release_only_set_is_exactly(mutated, _WANT_RELEASE_ONLY_CASES)
    except AssertionError as exc:
        print('  (mutation correctly raised: %s)' % exc)
    else:
        raise AssertionError(
            'planting a 3rd case (test.tpv36) into RELEASE_ONLY was NOT '
            'detected -- the exactness guard is not load-bearing')
    print('  PASS  a planted 3rd case in RELEASE_ONLY is detected '
          '(mutation-tested: real set passes, planted 3rd-case set fails)')


def check_release_only_never_reads_pass_fail():
    """RELEASE_ONLY must be a static dict literal, never a function of a
    pass/fail result -- rule 24: "A cell may be held out of the everyday run
    for cost; it may never be held out because it fails." Checked two ways:
    (1) matrix.RELEASE_ONLY is a plain dict whose every value is a string
    reason naming COST as the basis (not e.g. a bool computed from a test
    result); (2) source-level, `testsys.compare` (the ONE comparison
    implementation, per this module's own docstring) is never imported by
    matrix.py at all -- a module that cannot even call compare.* cannot be
    deciding RELEASE_ONLY from a pass/fail read."""
    if not isinstance(matrix.RELEASE_ONLY, dict):
        raise AssertionError('matrix.RELEASE_ONLY is %r, not a plain dict -- '
                             'it must be static, named data' % type(matrix.RELEASE_ONLY))
    for case, reason in matrix.RELEASE_ONLY.items():
        if not isinstance(reason, str) or 'cost' not in reason.lower():
            raise AssertionError(
                '%s: RELEASE_ONLY reason does not name COST as the basis -- '
                'got %r' % (case, reason))
    import ast as _ast
    tree = _ast.parse(open(os.path.join(ROOT, 'testsys', 'matrix.py')).read())
    imported = set()
    for node in _ast.walk(tree):
        if isinstance(node, _ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, _ast.ImportFrom) and node.module:
            imported.add(node.module)
            imported.update('%s.%s' % (node.module, a.name) for a in node.names)
    bad = {m for m in imported if m == 'compare' or m.endswith('.compare')}
    if bad:
        raise AssertionError(
            'matrix.py IMPORTS %r -- a table module that can call the '
            'comparison implementation could decide RELEASE_ONLY from a '
            'pass/fail read; matrix.py must stay pure declared data' % bad)
    print('  PASS  matrix.RELEASE_ONLY is static named data, each reason '
          'citing COST, and matrix.py never imports testsys.compare')


def check_everyday_excludes_release_only_release_includes_it():
    """Behavioural, not just a dict inspection: run_e2e.select() with the
    DEFAULT args (no --release, no --cases/--backends) must not runnable any
    (case, backend) whose case is in RELEASE_ONLY, and must list them in
    `release_only`; the SAME module's select() with --release must runnable
    them, at matrix.CASE_BOUND[case]/matrix.gate_term_for(case) -- the SAME
    bound/term as every other case (no separate, looser release-only path)."""
    e2e = _load_module('run_e2e_for_release_only_selection', E2E_PY)
    import argparse
    everyday_args = argparse.Namespace(ci=False, cases=None, backends=None, release=False)
    release_args = argparse.Namespace(ci=False, cases=None, backends=None, release=True)
    ev_runnable, ev_unsupported, ev_label, ev_explicit, ev_release_only = e2e.select(everyday_args)
    rel_runnable, rel_unsupported, rel_label, rel_explicit, rel_release_only = e2e.select(release_args)
    ev_cases_in_runnable = {c for c, _ in ev_runnable if c in matrix.RELEASE_ONLY}
    if ev_cases_in_runnable:
        raise AssertionError(
            'the EVERYDAY selection (no --release) runnable %r -- these are '
            'matrix.RELEASE_ONLY cases and must be excluded by default'
            % sorted(ev_cases_in_runnable))
    ev_release_only_cases = {c for c, _, _ in ev_release_only}
    if ev_release_only_cases != set(matrix.RELEASE_ONLY):
        raise AssertionError(
            'the EVERYDAY selection\'s release_only list names cases %r, '
            'want exactly matrix.RELEASE_ONLY\'s keys %r'
            % (sorted(ev_release_only_cases), sorted(matrix.RELEASE_ONLY)))
    rel_cases_in_runnable = {c for c, _ in rel_runnable if c in matrix.RELEASE_ONLY}
    if rel_cases_in_runnable != set(matrix.RELEASE_ONLY):
        raise AssertionError(
            'the RELEASE selection (--release) runnable only %r of '
            'matrix.RELEASE_ONLY\'s cases %r -- --release must include every '
            'one of them' % (sorted(rel_cases_in_runnable), sorted(matrix.RELEASE_ONLY)))
    if rel_release_only:
        raise AssertionError(
            '--release selection still reports %r as held-back release_only '
            '-- --release must INCLUDE them, not merely acknowledge them'
            % rel_release_only)
    print('  PASS  everyday selection excludes %s (listed as release_only); '
          '--release includes them (runnable, same bound/term as every '
          'other case)' % sorted(matrix.RELEASE_ONLY))


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
    run_e2e.py is never invoked.

    Since item 3 (2026-09-30), `run_release()` first asks
    `_release_sweep_carry_forward()` whether committed evidence already
    covers this tree and, if so, skips the subprocess call entirely (by
    design -- see that function's docstring). This check is about the
    ARGV SHAPE of the sweep invocation, not the carry-forward decision
    (that has its own coverage), so it forces the sweep-needed path with
    a stub, regardless of whether real evidence happens to carry forward
    in whatever tree this runs in."""
    mod = _load_run_module()
    if 'release' not in mod.RUNNERS:
        raise AssertionError("run.py has no 'release' tier registered in RUNNERS")
    captured = []
    orig_call = mod.subprocess.call
    orig_carry_forward = mod._release_sweep_carry_forward

    def fake_call(cmd, *a, **k):
        captured.append(list(cmd))
        return 0

    def force_sweep_needed():
        return False, 'test stub: forcing the sweep-needed path'

    mod.subprocess.call = fake_call
    mod._release_sweep_carry_forward = force_sweep_needed
    try:
        rc = mod.RUNNERS['release']()
    finally:
        mod.subprocess.call = orig_call
        mod._release_sweep_carry_forward = orig_carry_forward
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


def check_release_tier_honours_carry_forward_skip():
    """The other half of item 3 (2026-09-30): when
    `_release_sweep_carry_forward()` says committed evidence already
    covers this tree, `run_release()` must return 0 and make ZERO
    subprocess calls -- not run the sweep anyway "to be safe". Forced
    here with a stub so it does not depend on whatever evidence this
    tree's own history happens to carry."""
    mod = _load_run_module()
    captured = []
    orig_call = mod.subprocess.call
    orig_carry_forward = mod._release_sweep_carry_forward

    def fake_call(cmd, *a, **k):
        captured.append(list(cmd))
        return 0

    def force_carry_forward():
        return True, 'test stub: forcing the carry-forward path'

    mod.subprocess.call = fake_call
    mod._release_sweep_carry_forward = force_carry_forward
    try:
        rc = mod.RUNNERS['release']()
    finally:
        mod.subprocess.call = orig_call
        mod._release_sweep_carry_forward = orig_carry_forward
    if rc != 0:
        raise AssertionError(
            'run_release() returned %r when carry-forward applies (expected 0)'
            % rc)
    if captured:
        raise AssertionError(
            'run_release() made %d subprocess call(s) when carry-forward '
            'applies -- it must skip the sweep entirely: %r'
            % (len(captured), captured))
    print('  PASS  release tier makes 0 subprocess calls when carry-forward '
          'applies')


def main():
    print('Regression guard: ONE term everywhere (no --term flag, no second '
          'reference file, no per-case full-term table)')
    checks = [check_matrix_has_no_second_term_table,
              check_compare_has_no_second_reference_name,
              check_no_case_has_two_reference_files,
              check_mutation_a_planted_second_reference_file_is_caught,
              check_apply_term_override_takes_no_term_argument_and_is_unconditional,
              check_case_term_override_is_exactly_tpv22_tpv23,
              check_mutation_a_third_case_in_term_override_is_caught,
              check_apply_term_override_honours_the_named_exception,
              check_release_only_is_exactly_tpv22_tpv23,
              check_mutation_a_third_case_in_release_only_is_caught,
              check_release_only_never_reads_pass_fail,
              check_everyday_excludes_release_only_release_includes_it,
              check_run_e2e_help_lists_release_and_not_term,
              check_run_e2e_refuses_a_term_flag,
              check_ci_smoke_invocation_passes_neither_term_nor_release,
              check_release_tier_requests_the_release_selection,
              check_release_tier_honours_carry_forward_skip]
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
