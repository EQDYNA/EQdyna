"""Unit tests for testsys/change_class.py, the ONE shared change classifier
(owner + victor-reyes design review, 2026-09-30).

Both ways (rule 14a): every tracked path classifies to exactly one bucket
and an unknown path defaults PHYSICS; the import-trace check actually goes
red if a sweep dependency were relabeled INTERNAL.
"""
import ast
import subprocess

from conftest import REPO_ROOT
from testsys import change_class


def _tracked_paths():
    out = subprocess.run(['git', 'ls-files'], cwd=REPO_ROOT,
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.splitlines() if p.strip()]


def test_totality_every_tracked_path_classifies_to_one_bucket():
    paths = _tracked_paths()
    assert len(paths) > 100, 'sanity: this should be the whole checkout'
    buckets = {change_class.PHYSICS, change_class.USER_FACING, change_class.INTERNAL}
    for p in paths:
        assert change_class.classify_path(p) in buckets, p


def test_unclassifiable_path_defaults_physics():
    for p in ('brand/new/unclassified/thing.xyz', 'a_new_root_file.txt',
             'testsys/a_future_file_nobody_has_named_yet.py'):
        assert change_class.classify_path(p) == change_class.PHYSICS, p


def test_representative_physics_paths():
    for p in ('src/fortran/driver.f90', 'src/python/eqdyna/fric.py',
             'case_input/test.tpv8/user_defined_params.py',
             'test.reference.results/test.tpv8/frt.canonical.txt',
             'testsys/matrix.py', 'testsys/compare.py',
             'testsys/e2e/run_e2e.py', 'testsys/e2e/run_e2e_full.py',
             'testNameList.py', 'install-eqdyna.sh',
             'src/fortran/makefile', 'scripts/case.setup',
             'scripts/lib.py', 'scripts/defaultParameters.py',
             'scripts/machines.py', 'scripts/src_hash.py',
             'scripts/create.newcase', 'testsys/common.py',
             'testsys/runlock.py', 'testsys/profile_record.py',
             'testsys/profile_schema.py', 'testsys/perf/ledger.py',
             'testsys/perf/run_numa_scaling.py', 'testsys/run.py',
             'testsys/change_class.py'):
        assert change_class.classify_path(p) == change_class.PHYSICS, p


def test_representative_user_facing_paths():
    for p in ('VERSION', 'README.md', 'Dockerfile', '.dockerignore',
             'ubuntu.env.sh', '.github/workflows/publish.yml',
             'docs/user/getting-started.md', 'docs/user/parameters.md'):
        assert change_class.classify_path(p) == change_class.USER_FACING, p


def test_representative_internal_paths():
    for p in ('pathway_forward.md', 'PROJECT_RULES.md', 'CLAUDE.md',
             'docs/notes/foo.md', 'docs/evidence/sweep-abc/summary.json',
             'docs/perf_ledger.jsonl', 'testsys/unit/test_matrix.py',
             'testsys/regression/test_dipping_fault_y_split.py',
             'testsys/regression/check_pretag_ci.py',
             'testsys/hooks/pre-commit', 'testsys/ci_shard.py',
             'testsys/ci_status.py', 'testsys/pr_policy.py',
             'testsys/check_board_separation.py',
             'testsys/check_git_config_no_credential.py',
             '.github/workflows/test.yml',
             'scripts/gmMain.m', 'scripts/calc_shear_mod.m',
             'scripts/figures/make_tpv29_tpv30_overlay.py'):
        assert change_class.classify_path(p) == change_class.INTERNAL, p


def test_check_star_is_internal_wherever_it_lives():
    """basename-based, not just testsys/regression/ -- a check_* file
    anywhere is a checker ABOUT the process, not part of what it checks."""
    assert change_class.classify_path('somewhere/else/check_something.py') == change_class.INTERNAL


# --------------------------------------------------------------------------
# is_output_change: rule 27's NARROW criterion, deliberately not the same
# as classify_path == PHYSICS (review item 6).
# --------------------------------------------------------------------------
def test_is_output_change_true_for_narrow_set_false_for_broader_physics():
    assert change_class.is_output_change('src/fortran/library_output.f90')
    assert change_class.is_output_change(
        'test.reference.results/test.tpv8/frt.canonical.txt')
    assert change_class.is_output_change('src/fortran/driver.f90')  # no diff given -> True
    # PHYSICS but NOT an output change on its own:
    assert not change_class.is_output_change('testsys/matrix.py')
    assert not change_class.is_output_change('testsys/e2e/run_e2e.py')
    assert not change_class.is_output_change('scripts/lib.py')


def test_is_output_change_comment_only_fortran_diff_is_false():
    diff = ('diff --git a/x b/x\n--- a/x\n+++ b/x\n@@ -1 +1 @@\n'
           '-! old note\n+! new note\n')
    assert not change_class.is_output_change('src/fortran/fric.f90', diff)


def test_is_output_change_openmp_sentinel_diff_is_true():
    diff = ('diff --git a/x b/x\n--- a/x\n+++ b/x\n@@ -1 +1 @@\n'
           '-! plain comment\n+!$omp parallel do\n')
    assert change_class.is_output_change('src/fortran/fric.f90', diff)


# --------------------------------------------------------------------------
# import-trace: everything the sweep actually imports must be PHYSICS.
# Both ways (rule 14a): the real set passes; a fixture that relabels
# ledger.py INTERNAL makes the SAME check fail.
# --------------------------------------------------------------------------
SWEEP_DEPENDENCY_PATHS = {
    'testsys.common': 'testsys/common.py',
    'testsys.runlock': 'testsys/runlock.py',
    'testsys.profile_record': 'testsys/profile_record.py',
    'testsys.profile_schema': 'testsys/profile_schema.py',
    'testsys.compare': 'testsys/compare.py',
    'testsys.matrix': 'testsys/matrix.py',
    'testsys.perf.ledger': 'testsys/perf/ledger.py',
    'testsys.perf.run_numa_scaling': 'testsys/perf/run_numa_scaling.py',
}


def _assert_all_physics(module_paths):
    bad = [(name, path) for name, path in module_paths.items()
          if change_class.classify_path(path) != change_class.PHYSICS]
    if bad:
        raise AssertionError('not classified PHYSICS: %r' % bad)


def _imported_module_names(py_path):
    """Every module name a `import X` / `from X import Y` statement in
    py_path names, at any depth in the source (not just top-level) --
    a real, if shallow (no transitive following), AST-based trace."""
    with open(py_path) as f:
        tree = ast.parse(f.read(), filename=py_path)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
            for alias in node.names:
                names.add('%s.%s' % (node.module, alias.name))
    return names


def test_run_e2e_actually_imports_the_named_sweep_dependencies():
    """Sanity that SWEEP_DEPENDENCY_PATHS is not a fictional list: run_e2e.py
    really does import at least profile_record, compare, and matrix (via
    `from testsys import ...`)."""
    import os
    names = _imported_module_names(
        os.path.join(REPO_ROOT, 'testsys', 'e2e', 'run_e2e.py'))
    assert 'profile_record' in names or any('profile_record' in n for n in names)
    assert any(n.endswith('compare') or n == 'compare' for n in names)
    assert any(n.endswith('matrix') or n == 'matrix' for n in names)


def test_sweep_dependencies_are_all_physics():
    _assert_all_physics(SWEEP_DEPENDENCY_PATHS)


def test_import_trace_check_goes_red_if_ledger_relabeled_internal(monkeypatch):
    real = change_class.classify_path
    monkeypatch.setattr(
        change_class, 'classify_path',
        lambda p: change_class.INTERNAL if p == 'testsys/perf/ledger.py' else real(p))
    try:
        _assert_all_physics(SWEEP_DEPENDENCY_PATHS)
    except AssertionError as exc:
        assert 'testsys/perf/ledger.py' in str(exc)
    else:
        raise AssertionError(
            'expected _assert_all_physics to fail once ledger.py is '
            'relabeled INTERNAL -- the import-trace check cannot see this '
            'class of regression')
