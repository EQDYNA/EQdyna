"""Unit tests for testsys/matrix.py -- the e2e sweep's support table.

The table is the thing that decides what a green sweep MEANS, so these tests
are about its honesty, not its contents: that every gated case declares a
bound, that no bound is looser than the repo's outer threshold, that an
unsupported cell is declared rather than absent, and that the coverage report
can never read as broader than the selection it describes.
"""
import pytest

from conftest import REPO_ROOT  # noqa: F401  (puts the repo root on sys.path)
from testsys import matrix


def test_every_gated_case_declares_a_bound_and_a_gate():
    # testNameList.py is the case list; matrix.py must cover it exactly.
    # (matrix.py also enforces this at import time -- this test states the
    # contract where a reader looks for it.)
    assert set(matrix.CASES) == set(matrix.CASE_BOUND) == set(matrix.GATE)


def test_no_case_bound_is_looser_than_the_outer_threshold():
    # Rule 5: one calibrated definition of pass. A per-case bound may be
    # tighter than THRESHOLD; it may never be looser.
    looser = {c: b for c, b in matrix.CASE_BOUND.items()
              if b is not None and b > matrix.THRESHOLD}
    assert looser == {}


def test_drv_a6_is_gated_on_a_flip_budget_not_a_scalar():
    # The bistable case is gated on bulk agreement PLUS an explicit flip
    # budget. A scalar max-abs bound cannot distinguish "a few hundred marginal
    # nodes flipped" from "everything drifted", so acquiring one here would be
    # a regression, and so would quietly widening the budget.
    assert matrix.CASE_BOUND['test.drv.a6'] is None
    assert matrix.GATE['test.drv.a6'] == 'flip-budget'
    assert matrix.DRV_A6['total_flip_bound'] == 450


def test_every_gated_case_has_a_committed_reference():
    # matrix.py enforces this at import; stated here where a reader looks for
    # the contract. A gated case with no reference cannot be compared, and
    # "could not compare" must never read as "passed".
    import os
    for c in matrix.CASES:
        ref = os.path.join(matrix.REPO_ROOT, 'test.reference.results', c,
                           'frt.canonical.txt')
        assert os.path.isfile(ref), '%s is gated but has no reference' % c


def test_default_selection_accounts_for_every_cell_in_the_table():
    runnable, unsupported = matrix.cells()
    assert len(runnable) + len(unsupported) == len(matrix.CASES) * len(matrix.BACKENDS)


def test_unsupported_table_is_empty_but_the_mechanism_still_works(monkeypatch):
    """UNSUPPORTED is EMPTY now -- friclaw=2 was the last gap, and closing it
    made test.meng2023a/test.meng2023cb runnable on both python backends for
    the first time. So this pins the MECHANISM against a synthetic entry
    rather than against whichever real gap happens to be open.

    Pinning a regression test to a real gap is how it dies the moment the gap
    is fixed: the previous version asserted `runnable == []` for
    test.meng2023a x python-jax and would now fail because that cell works.
    """
    assert matrix.UNSUPPORTED == {}, \
        'a newly declared-unsupported cell needs its reason reviewed here'

    monkeypatch.setitem(matrix.UNSUPPORTED, ('test.tpv8', 'python-jax'),
                        'synthetic reason for this test')
    runnable, unsupported = matrix.cells(cases=['test.tpv8'],
                                         backends=['python-jax'])
    assert runnable == []
    assert len(unsupported) == 1
    case, backend, reason = unsupported[0]
    assert (case, backend) == ('test.tpv8', 'python-jax')
    assert reason == 'synthetic reason for this test'


def test_every_gated_cell_runs_now_that_nothing_is_unsupported():
    """The sweep covers every cell of the gated table. Recorded as an
    assertion so losing coverage requires deleting a test, not just quietly
    editing a table."""
    runnable, unsupported = matrix.cells()
    assert unsupported == []
    assert len(runnable) == len(matrix.CASES) * len(matrix.BACKENDS)


def test_unsupported_reason_raises_for_a_supported_cell():
    # Asking why a supported cell was skipped is itself a bug: there is no
    # "skipped but fine" state in this table.
    with pytest.raises(KeyError):
        matrix.unsupported_reason('test.tpv8', 'fortran')


def test_unknown_case_or_backend_raises_rather_than_selecting_nothing():
    # A typo must not quietly produce an empty, trivially-green selection.
    with pytest.raises(ValueError):
        matrix.cells(cases=['test.tpv99'])
    with pytest.raises(ValueError):
        matrix.cells(backends=['fortran-ish'])


def test_coverage_report_states_both_the_selection_and_the_table_size():
    runnable, unsupported = matrix.cells(cases=['test.tpv8'], backends=['fortran'])
    text = '\n'.join(matrix.coverage_report(runnable, unsupported, 'one cell'))
    # The failure this guards: "SUCCESS (5/5 cases)" read as completeness
    # because the denominator was the list, not the table.
    assert '1 of %d cells' % (len(matrix.CASES) * len(matrix.BACKENDS)) in text
    assert 'NOT in this selection, %d cell(s)' % (
        len(matrix.CASES) * len(matrix.BACKENDS) - 1) in text


def test_coverage_report_names_every_unsupported_cell_with_its_reason():
    runnable, unsupported = matrix.cells()
    text = '\n'.join(matrix.coverage_report(runnable, unsupported, 'all'))
    for case, backend, reason in unsupported:
        assert case in text and backend in text
        assert reason.split('.')[0] in text


def test_ci_cells_are_all_in_the_table_and_within_the_measured_runner():
    for cell in matrix.CI_CELLS:
        assert cell[0] in matrix.CASES and cell[1] in matrix.BACKENDS
        rss = matrix.MEASURED_PEAK_RSS_GB.get(cell)
        if rss is not None:
            assert rss < matrix.CI_RUNNER_RAM_GB, (
                '%s x %s peaks at %.2f GB on a %.0f GB runner that is already '
                'holding the Fortran build' % (cell[0], cell[1], rss,
                                               matrix.CI_RUNNER_RAM_GB))


def test_every_backend_declares_what_artifacts_it_produces():
    assert set(matrix.ARTIFACTS) == set(matrix.BACKENDS)
    assert all(a for a in matrix.ARTIFACTS.values()), \
        'a backend that produces no compared artifact would run and gate nothing'
