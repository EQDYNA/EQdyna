"""Unit tests for testsys/change_class.py, the ONE shared change classifier
(owner + victor-reyes design review, 2026-09-30).

Both ways (rule 14a): every tracked path classifies to exactly one bucket
and an unknown path defaults PHYSICS; the import-trace check actually goes
red if a sweep dependency were relabeled INTERNAL.
"""
import ast
import importlib
import importlib.util
import os
import subprocess
import sys

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


def test_totality_every_tracked_path_matches_an_explicit_rule():
    """victor-reyes audit of PR #62: classify_path's totality can never go
    red on its own (it ends in a bare `return PHYSICS`), so this asserts
    the STRONGER claim -- every tracked path matches an EXPLICIT rule
    (USER_FACING, INTERNAL, or a REVIEWED default-PHYSICS prefix/exact name
    in change_class.explicit_bucket), reporting by name any path that falls
    through to the bare, UNREVIEWED default. A brand-new top-level file
    outside every reviewed prefix makes this go red until it is reviewed
    and added to one of change_class.py's explicit lists."""
    paths = _tracked_paths()
    unreviewed = [p for p in paths if change_class.explicit_bucket(p) is None]
    assert not unreviewed, (
        'tracked path(s) fall through to the bare, unreviewed PHYSICS '
        'default -- review each and add it to change_class.py\'s '
        'USER_FACING/INTERNAL lists or DEFAULT_PHYSICS_REVIEWED_* set: %r'
        % unreviewed)


def test_explicit_bucket_mutation_self_check():
    """rule 10a: prove the totality test above actually exercises
    explicit_bucket's None case, not just a vacuous pass -- a path outside
    every list/prefix returns None, and one inside a reviewed prefix does not."""
    assert change_class.explicit_bucket('brand/new/unreviewed/thing.xyz') is None
    assert change_class.explicit_bucket('src/fortran/driver.f90') == change_class.PHYSICS
    assert change_class.explicit_bucket('README.md') == change_class.USER_FACING
    assert change_class.explicit_bucket('pathway_forward.md') == change_class.INTERNAL


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
             'docs/user/getting-started.md', 'docs/user/parameters.md',
             'Docker.guide.md'):
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
             'scripts/figures/make_tpv29_tpv30_overlay.py',
             'LICENSE', 'pastReleaseNotes.md', '.gitignore'):
        assert change_class.classify_path(p) == change_class.INTERNAL, p


def test_check_star_is_internal_only_under_testsys():
    """victor-reyes audit of PR #62, finding 4: scoped to testsys/ (root and
    any subdirectory, e.g. testsys/regression/) -- a check_* basename
    OUTSIDE testsys/ is NOT this repo's release-process tooling and must
    not be exempted from review just by that name."""
    assert change_class.classify_path('testsys/check_something.py') == change_class.INTERNAL
    assert change_class.classify_path(
        'testsys/regression/check_something_else.py') == change_class.INTERNAL
    assert change_class.classify_path('somewhere/else/check_something.py') == change_class.PHYSICS
    assert change_class.classify_path('scripts/check_something.py') == change_class.PHYSICS


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
# import-trace: everything the sweep actually imports must be PHYSICS. A
# REAL trace (victor-reyes audit of PR #62: the previous version was a
# hand-maintained dict, which could drift from what run_e2e.py/run.py
# actually import without this test noticing). Both ways (rule 14a): the
# real set passes; relabeling ledger.py INTERNAL, or a genuinely NEW
# INTERNAL import landing in run_e2e.py, makes the SAME mechanism fail.
# --------------------------------------------------------------------------
def _imported_module_names(py_path):
    """Every module name a `import X` / `from X import Y` statement in
    py_path names, at any depth in the source (not just top-level) --
    a real, if shallow (no transitive following), AST-based trace. Used
    only for scripts/case.setup below, which is a script (no .py suffix),
    not a package -- it cannot be `importlib.import_module`d the way
    testsys.run / testsys.e2e.run_e2e are."""
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


def _repo_relpath(abspath):
    return os.path.relpath(abspath, REPO_ROOT).replace(os.sep, '/')


def _new_repo_local_modules(before_names):
    """{module_name: repo-relative path} for every entry sys.modules gained
    since `before_names` whose __file__ lives inside REPO_ROOT -- stdlib and
    third-party imports (jax, numpy, mpi4py, ...) are excluded automatically
    because their __file__ lives elsewhere, never by name-listing them."""
    out = {}
    for name in set(sys.modules) - before_names:
        mod = sys.modules.get(name)
        f = getattr(mod, '__file__', None)
        if not f:
            continue
        f = os.path.abspath(f)
        if f.startswith(REPO_ROOT + os.sep):
            out[name] = _repo_relpath(f)
    return out


def _collect_repo_modules_transitively(mod, seen=None):
    """{module_name: repo-relative path}, walking `mod`'s own namespace for
    every attribute that is itself a module object with a __file__ inside
    REPO_ROOT, and recursing into each ONE FOUND -- so `run_e2e`'s `import
    profile_record` -> profile_record.py's own `import ledger` -> ledger.py
    is followed transitively. This reads each module's namespace directly
    (what `import X` / `from X import Y` actually BOUND as a name), never
    sys.modules-before/after diffing -- which breaks the moment this
    function (or an earlier test in the same process) has already imported
    the same module once, since a second import is served straight from
    the cache and adds NOTHING new to sys.modules to diff against. Stops
    (does not descend) the moment a node's own __file__ is NOT inside
    REPO_ROOT -- a stdlib/third-party module's internals are never walked,
    only its arrival at this repo's own module is ever recorded (it isn't,
    since it's excluded), so this cannot explode into numpy/jax internals."""
    if seen is None:
        seen = {}
    f = getattr(mod, '__file__', None)
    name = getattr(mod, '__name__', None)
    if not f or not name or name in seen:
        return seen
    f = os.path.abspath(f)
    if not f.startswith(REPO_ROOT + os.sep):
        return seen
    seen[name] = _repo_relpath(f)
    import types
    for val in vars(mod).values():
        if isinstance(val, types.ModuleType):
            _collect_repo_modules_transitively(val, seen)
    return seen


def discover_sweep_dependency_paths():
    """{module_name: repo-relative path} for every repo file the everyday
    sweep path (`testsys/run.py`, `testsys/e2e/run_e2e.py`) actually loads,
    via the three mechanisms this repo's sweep path uses:
      1. ordinary `import`, transitively -- see
         _collect_repo_modules_transitively above (so a future new `from
         testsys import X` in either module, or a new import inside a
         module they already import, is picked up with no test edit).
      2. testsys/run.py's own by-PATH loaders (`_load_machines`,
         `_load_src_hash` -- scripts/machines.py, scripts/src_hash.py):
         these do not register as an attribute of run.py's own namespace
         (the loaded module is a local inside the function, not a global),
         so they are CALLED for real and their return value's __file__ read.
      3. scripts/case.setup's own local (repo) imports (`from lib import
         ...`, `from machines import ...`): it is a script, not a package,
         so it is AST-parsed (see _imported_module_names) and each name
         resolvable to a real scripts/<name>.py file is included.
    """
    run_mod = importlib.import_module('testsys.run')
    run_e2e_mod = importlib.import_module('testsys.e2e.run_e2e')
    paths = {}
    for m in (run_mod, run_e2e_mod):
        _collect_repo_modules_transitively(m, paths)

    for loader_name, loader in (('scripts.machines (via run.py._load_machines)',
                                 run_mod._load_machines),
                                ('scripts.src_hash (via run.py._load_src_hash)',
                                 run_mod._load_src_hash)):
        loaded = loader()
        paths[loader_name] = _repo_relpath(os.path.abspath(loaded.__file__))

    case_setup = os.path.join(REPO_ROOT, 'scripts', 'case.setup')
    for name in _imported_module_names(case_setup):
        local = os.path.join(REPO_ROOT, 'scripts', name.split('.')[-1] + '.py')
        if os.path.isfile(local):
            paths['case.setup: ' + name] = _repo_relpath(local)
    return paths


def _assert_all_physics(module_paths):
    bad = [(name, path) for name, path in module_paths.items()
          if change_class.classify_path(path) != change_class.PHYSICS]
    if bad:
        raise AssertionError('not classified PHYSICS: %r' % bad)


def test_real_trace_finds_the_expected_shape_of_dependencies():
    """Sanity that the real trace is actually tracing something, not
    silently returning an empty/trivial result: it must find run.py's own
    file, run_e2e.py's own file, a plain `from testsys import ...` name
    (matrix.py), a by-path loader target (scripts/machines.py), and a
    case.setup local import (scripts/lib.py)."""
    paths = set(discover_sweep_dependency_paths().values())
    assert 'testsys/run.py' in paths
    assert 'testsys/e2e/run_e2e.py' in paths
    assert 'testsys/matrix.py' in paths
    assert 'scripts/machines.py' in paths
    assert 'scripts/lib.py' in paths
    assert len(paths) >= 8, 'suspiciously few real dependencies found: %r' % paths


def test_sweep_dependencies_are_all_physics():
    _assert_all_physics(discover_sweep_dependency_paths())


def test_import_trace_check_goes_red_if_ledger_relabeled_internal(monkeypatch):
    paths = discover_sweep_dependency_paths()
    assert 'testsys/perf/ledger.py' in paths.values(), (
        'sanity: the real trace no longer finds testsys/perf/ledger.py at '
        'all -- this test would pass for the wrong reason')
    real = change_class.classify_path
    monkeypatch.setattr(
        change_class, 'classify_path',
        lambda p: change_class.INTERNAL if p == 'testsys/perf/ledger.py' else real(p))
    try:
        _assert_all_physics(paths)
    except AssertionError as exc:
        assert 'testsys/perf/ledger.py' in str(exc)
    else:
        raise AssertionError(
            'expected _assert_all_physics to fail once ledger.py is '
            'relabeled INTERNAL -- the import-trace check cannot see this '
            'class of regression')


def test_import_trace_goes_red_on_a_real_added_internal_import(tmp_path):
    """victor-reyes's own named mutation: adding
    `from testsys.regression import check_release_due` to run_e2e.py must
    make the trace fail. Proved against a TEMP COPY of run_e2e.py's actual
    current source with that one line inserted -- the checked-in file is
    never touched -- loaded under a throwaway module name so it cannot
    clobber the real testsys.e2e.run_e2e already cached in sys.modules."""
    real_path = os.path.join(REPO_ROOT, 'testsys', 'e2e', 'run_e2e.py')
    with open(real_path) as f:
        source = f.read()
    marker = 'from testsys import common, compare, frt_canonical, matrix, runlock'
    assert marker in source, ('run_e2e.py\'s import line moved or changed -- '
                              'update this test\'s insertion point')
    mutated = source.replace(
        marker, marker + '\nfrom testsys.regression import check_release_due  '
                         '# MUTATION under test, never real', 1)
    mutated_file = tmp_path / 'run_e2e_mutated.py'
    mutated_file.write_text(mutated)

    spec = importlib.util.spec_from_file_location(
        'run_e2e_mutated_for_test', str(mutated_file))
    mod = importlib.util.module_from_spec(spec)
    before = set(sys.modules)
    try:
        spec.loader.exec_module(mod)
        added = _new_repo_local_modules(before)
        assert 'testsys/regression/check_release_due.py' in added.values(), (
            'the mutation\'s new import did not register in sys.modules '
            'with a repo-local __file__ -- this trace mechanism cannot see '
            'it: %r' % added)
        try:
            _assert_all_physics(added)
        except AssertionError as exc:
            assert 'check_release_due.py' in str(exc)
        else:
            raise AssertionError(
                'a real `from testsys.regression import check_release_due` '
                'added to run_e2e.py was NOT caught -- check_release_due.py '
                'classifies INTERNAL and must never be a real sweep import')
    finally:
        for name in set(sys.modules) - before:
            del sys.modules[name]
