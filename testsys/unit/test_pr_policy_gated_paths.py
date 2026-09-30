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
