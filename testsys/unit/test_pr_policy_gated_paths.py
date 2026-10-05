"""Unit tests for testsys/pr_policy.py's is_gated_path/touches_gated_paths
WIDENING (item 1, review CRITICAL, 2026-09-30): a path classified PHYSICS by
testsys/change_class.py must gate even outside the historical GATED_PREFIXES
(src/, testsys/, .github/) -- never narrowed, only unioned.

Both ways (rule 14a): a PHYSICS path outside GATED_PREFIXES now gates (it did
not before this change); a genuinely non-PHYSICS path outside GATED_PREFIXES
still does not gate, so the widening did not turn into "gate everything".
"""
from testsys import change_class, pr_policy


def test_physics_path_outside_gated_prefixes_now_gates():
    """case_input/, scripts/, testNameList.py, install-eqdyna.sh,
    test.reference.results/ are none of src/, testsys/, .github/ -- but each
    classifies PHYSICS, so each must gate per the review-1 union."""
    for p in ('case_input/test.tpv8/user_defined_params.py',
             'scripts/lib.py', 'scripts/machines.py', 'testNameList.py',
             'install-eqdyna.sh',
             'test.reference.results/test.tpv8/frt.canonical.txt'):
        assert change_class.classify_path(p) == change_class.PHYSICS, p
        assert pr_policy.is_gated_path(p), p
    assert pr_policy.touches_gated_paths(
        ['README.md', 'testNameList.py']) == ['testNameList.py']


def test_non_physics_path_outside_gated_prefixes_still_does_not_gate():
    """A genuinely INTERNAL or USER_FACING path outside GATED_PREFIXES must
    NOT gate -- the union only ADDS PHYSICS paths, it does not gate
    everything."""
    for p in ('README.md', 'VERSION', 'docs/notes/foo.md',
             'pathway_forward.md', 'docs/user/getting-started.md'):
        assert change_class.classify_path(p) != change_class.PHYSICS, p
        assert not pr_policy.is_gated_path(p), p
    assert pr_policy.touches_gated_paths(['README.md', 'docs/notes/foo.md']) == []


def test_historical_gated_prefixes_still_gate_even_when_internal():
    """testsys/unit/ and .github/ classify INTERNAL under change_class (no
    sweep needed) but must still gate under pr_policy: GATED_PREFIXES is a
    SEPARATE, unweakened answer to "does this need a PR" (module docstring)."""
    for p in ('testsys/unit/test_something.py', '.github/workflows/test.yml',
             'testsys/regression/check_pretag_ci.py'):
        assert change_class.classify_path(p) == change_class.INTERNAL, p
        assert pr_policy.is_gated_path(p), p


def _lane(files):
    """Mirror pr_policy.py's `pr-lane` CLI mode's own classification
    (main()'s `lane = 'full' if gated else 'fast'`) without needing a real
    git repo -- lane classification is a pure function of a changed-file
    set, see compute_pr_lane_files's docstring for how that set itself is
    obtained from a real PR diff (covered end-to-end, CLI included, by
    testsys/regression/test_pr_policy_guard.py's cases 10-13)."""
    return 'full' if pr_policy.touches_gated_paths(files) else 'fast'


def test_pr_lane_readme_only_is_fast():
    """A README-only diff -- the representative docs-only PR -- is LANE=fast."""
    assert _lane(['README.md']) == 'fast'


def test_pr_lane_src_only_is_full():
    """A diff touching src/ -- the representative code PR -- is LANE=full."""
    assert _lane(['src/fortran/driver.f90']) == 'full'


def test_pr_lane_mixed_diff_is_full():
    """A diff touching BOTH docs and a gated path is LANE=full -- one gated
    file is enough to pull the whole PR into the full lane, same as
    evaluate_commit_gate's mixed-commit case 4 in test_pr_policy_guard.py."""
    assert _lane(['README.md', 'testsys/pr_policy.py']) == 'full'
    assert _lane(['README.md', 'docs/notes/foo.md']) == 'fast'
