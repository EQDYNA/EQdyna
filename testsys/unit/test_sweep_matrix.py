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


def test_unsupported_table_is_exactly_the_mpi_opt_outs(monkeypatch):
    """UNSUPPORTED is no longer empty: python-jax-mpi (an OPTIONAL execution
    mode of python-jax, not a fourth full backend every case must support --
    see matrix.py's _MPI_OPT_IN_REASON) has exactly one case opted in
    (matrix.PY_MPI_RANKS), so every OTHER case x python-jax-mpi is declared
    unsupported with that one recorded reason.

    This pins the CONTENT of that specific, real gap (so an accidental opt-in
    or opt-out is caught) and then pins the MECHANISM against a synthetic
    entry layered on top, same as the previous version of this test did
    against a since-closed friclaw=2 gap -- pinning a regression test to a
    real gap is how it dies the moment the gap is closed, so the synthetic
    half must survive that.
    """
    expected_unsupported_backend_pairs = {
        (c, 'python-jax-mpi') for c in matrix.CASES
        if c not in matrix.PY_MPI_RANKS
    }
    assert set(matrix.UNSUPPORTED) == expected_unsupported_backend_pairs
    assert len(matrix.UNSUPPORTED) == len(matrix.CASES) - len(matrix.PY_MPI_RANKS)
    # every one of those entries carries the SAME recorded policy reason --
    # a per-case ad-hoc reason here would mean the exclusion is a coverage
    # gap being rationalised case by case rather than one stated policy.
    reasons = set(matrix.UNSUPPORTED.values())
    assert len(reasons) == 1
    assert 'not opted into the optional MPI execution mode' in next(iter(reasons))
    # no case outside python-jax-mpi is unsupported: the three original
    # backends (fortran, python-numpy, python-jax) still cover every case.
    assert all(b == 'python-jax-mpi' for _c, b in matrix.UNSUPPORTED)

    monkeypatch.setitem(matrix.UNSUPPORTED, ('test.tpv8', 'python-jax'),
                        'synthetic reason for this test')
    runnable, unsupported = matrix.cells(cases=['test.tpv8'],
                                         backends=['python-jax'])
    assert runnable == []
    assert len(unsupported) == 1
    case, backend, reason = unsupported[0]
    assert (case, backend) == ('test.tpv8', 'python-jax')
    assert reason == 'synthetic reason for this test'


def test_every_cell_runs_except_the_declared_mpi_opt_outs():
    """The sweep covers every cell of the gated table except the python-jax-
    mpi cells that have not opted in. Recorded as an assertion so losing
    coverage on the three original backends requires deleting a test, not
    just quietly editing a table -- and so widening PY_MPI_RANKS is visible
    here as a runnable-count change rather than silent.
    """
    runnable, unsupported = matrix.cells()
    total = len(matrix.CASES) * len(matrix.BACKENDS)
    n_opted_in = len(matrix.PY_MPI_RANKS)
    n_opted_out = len(matrix.CASES) - n_opted_in
    assert len(unsupported) == n_opted_out
    assert len(runnable) == total - n_opted_out
    # nothing outside python-jax-mpi is ever declared unsupported
    assert all(b == 'python-jax-mpi' for _c, b, _r in unsupported)


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
