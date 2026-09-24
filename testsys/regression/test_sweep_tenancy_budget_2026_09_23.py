#! /usr/bin/env python3
"""
Regression guard: the 2026-09-23 speed-campaign tenancy fix (rules 2, 10).

BEFORE (docs/SESSION_LOG_2026-09-23_autopilot.md section NN): the sweep
billed every python cell at 1 core while a python-jax cell MEASURED 252-259%
CPU (this repo's own /usr/bin/time -v, cited on matrix.JAX_MEASURED_CORES;
corroborates the 249% figure section NN attributes to an earlier run), and
the default --jobs budget was a fixed `cores - 4` that ignored however many
cores a foreign job already held. Together those two under-counts let every
cell in the BEFORE sweep run at 2-3x its solo time under real contention.

Two independent things land together and each gets its own check here so one
regressing does not hide behind the other going green:

  1. cell_cost bills python-jax at its MEASURED core count
     (matrix.JAX_MEASURED_CORES, ceil'd -- currently 3) and fortran/
     python-jax-mpi at their real rank counts, never a guess. python-numpy
     LEFT the backend axis entirely 2026-09-23 (owner decision, same day):
     its old cost-1 branch (also a measured number -- eqdyna3d.
     _narrow_numpy_affinity pins it to one cpu) is gone with it, and
     cell_cost now RAISES for that backend name rather than silently
     costing it 1 -- see check_python_numpy_backend_is_rejected below, the
     mutation-tested guard against the axis widening back to include it.

  2. default_jobs_budget derives the default --jobs from cores measured FREE
     right now (testsys/perf/run_numa_scaling.cpu_busy_fractions, reused, not
     duplicated) minus a stated margin, never below the largest selected
     cell's own cost, and an explicit --jobs still wins outright with no
     measurement taken.

Cheap (rule 9): imports + plain-data assertions + direct calls against
synthetic cost/busy dicts. No solver, no MPI, no real /proc/stat sample (the
injected busy_fractions_fn stands in). Exits non-zero on any failure.
"""
import math
import os
import sys
import types

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

from testsys import matrix  # noqa: E402
import run_e2e  # noqa: E402


# --------------------------------------------------------------------------
# item 1: honest per-backend cell cost
# --------------------------------------------------------------------------
def check_fortran_cost_is_its_rank_count():
    case = 'test.tpv8'
    assert run_e2e.cell_cost(case, 'fortran') == matrix.FORTRAN_RANKS[case]
    print('  PASS  fortran cost == matrix.FORTRAN_RANKS[%r] (%d)'
          % (case, matrix.FORTRAN_RANKS[case]))


def check_jax_mpi_cost_is_its_rank_count():
    case = 'test.tpv8'  # the one case opted into python-jax-mpi
    assert run_e2e.cell_cost(case, 'python-jax-mpi') == matrix.PY_MPI_RANKS[case]
    print('  PASS  python-jax-mpi cost == matrix.PY_MPI_RANKS[%r] (%d)'
          % (case, matrix.PY_MPI_RANKS[case]))


def check_python_numpy_backend_is_rejected():
    """GUARD against the backend axis widening back to include python-numpy
    (2026-09-23 owner decision: 'I actually don't care numpy'). cell_cost
    must RAISE for it, not silently cost it 1 the way it did before the
    axis removed it -- a re-added numpy cell that still ran (billed for
    free) would be a worse regression than one that failed loudly."""
    try:
        run_e2e.cell_cost('test.tpv8', 'python-numpy')
    except ValueError as exc:
        print('  PASS  cell_cost raises ValueError for the retired '
              'python-numpy backend: %s' % exc)
        return
    raise AssertionError('cell_cost accepted "python-numpy" and returned a '
                         'cost -- the backend axis has silently widened back '
                         'to include it')


def check_python_jax_cost_is_measured_and_greater_than_one():
    """The defect this item exists to fix: python-jax billed at 1 undercounts
    a MEASURED ~2.5-core cell. The cost must equal ceil(matrix.
    JAX_MEASURED_CORES) exactly, and it must be strictly greater than the old
    hardcoded 1 -- a regression that reverted the branch to `return 1` must
    fail here, not merely differ from some unrelated number."""
    cost = run_e2e.cell_cost('test.tpv8', 'python-jax')
    expected = math.ceil(matrix.JAX_MEASURED_CORES)
    assert cost == expected, (
        'python-jax cell_cost is %r, expected ceil(matrix.JAX_MEASURED_CORES) '
        '= %r' % (cost, expected))
    assert cost > 1, (
        'python-jax cell_cost is %r -- this is the pre-fix undercounted '
        'value (1 core); a MEASURED jax cell costs more than one core '
        '(matrix.JAX_MEASURED_CORES = %.2f)' % (cost, matrix.JAX_MEASURED_CORES))
    print('  PASS  python-jax cost == ceil(matrix.JAX_MEASURED_CORES) == %d '
          '(> 1, the pre-fix value)' % cost)


def check_mutating_jax_measured_cores_changes_the_billed_cost():
    """Mutation, both ways: cell_cost must actually READ
    matrix.JAX_MEASURED_CORES rather than have the fix's number baked in as a
    second hardcoded literal. Set the constant to exactly 1.0 (the pre-fix
    value) and confirm cell_cost tracks it down to 1 -- proving the two are
    wired together, not merely coincidentally equal today."""
    real = matrix.JAX_MEASURED_CORES
    try:
        matrix.JAX_MEASURED_CORES = 1.0
        mutated_cost = run_e2e.cell_cost('test.tpv8', 'python-jax')
        assert mutated_cost == 1, (
            'setting matrix.JAX_MEASURED_CORES = 1.0 did not change '
            'cell_cost(..., "python-jax") to 1 -- cell_cost is not actually '
            'reading matrix.JAX_MEASURED_CORES (it is RED under this '
            'mutation, which is the point: this proves the wiring, it does '
            'not pass the real number back)')
    finally:
        matrix.JAX_MEASURED_CORES = real
    restored_cost = run_e2e.cell_cost('test.tpv8', 'python-jax')
    assert restored_cost == math.ceil(real)
    print('  PASS  cell_cost tracks matrix.JAX_MEASURED_CORES (mutated to '
          '1.0 -> cost 1; restored -> cost %d again)' % restored_cost)


def check_cell_cost_source_does_not_hardcode_python_jax_to_one():
    """Source-level tie, same shape as test_sweep_core_budget.py's own
    check: a future edit that special-cases python-jax back to a bare
    `return 1` (reverting the fix while leaving matrix.JAX_MEASURED_CORES in
    place, unread) must fail even though nothing above ran that exact code
    path with a hand-picked input the reverted branch happens to still
    satisfy."""
    path = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py')
    raw = open(path, errors='replace').read()
    src = '\n'.join(ln.split('#', 1)[0] for ln in raw.splitlines())
    assert "matrix.JAX_MEASURED_CORES" in src, (
        'run_e2e.py no longer references matrix.JAX_MEASURED_CORES anywhere '
        '-- the honest jax cost source may have been removed')
    # the exact branch: cell_cost's python-jax arm must read the constant.
    marker = "if backend == 'python-jax':"
    idx = src.find(marker)
    assert idx != -1, 'cell_cost no longer has a dedicated python-jax branch'
    branch = src[idx:idx + 200]
    assert 'JAX_MEASURED_CORES' in branch, (
        "cell_cost's python-jax branch no longer reads "
        'matrix.JAX_MEASURED_CORES: %r' % branch)
    print('  PASS  cell_cost\'s python-jax branch reads matrix.JAX_MEASURED_CORES')


# --------------------------------------------------------------------------
# item 2: tenancy-aware default budget
# --------------------------------------------------------------------------
def _fake_busy(fractions_by_cpu):
    """A busy_fractions_fn matching run_numa_scaling.cpu_busy_fractions's own
    signature (cpus, sample_s) -> {cpu: fraction}, injected so this guard
    never waits out a real >=2s /proc/stat sample and never depends on this
    box's own load."""
    def fn(cpus, sample_s=2.0):
        return {c: fractions_by_cpu[c] for c in cpus if c in fractions_by_cpu}
    return fn


def check_measured_free_cores_counts_correctly():
    """8 cpus, 3 busier than the ceiling -> 5 free."""
    fractions = {i: (0.9 if i < 3 else 0.05) for i in range(8)}
    free, busy, total, ceiling = run_e2e.measured_free_cores(
        total_cores=8, busy_ceiling=0.5, sample_s=0.0,
        busy_fractions_fn=_fake_busy(fractions))
    assert (free, busy, total, ceiling) == (5, 3, 8, 0.5), (
        (free, busy, total, ceiling))
    print('  PASS  measured_free_cores(8 cpus, 3 busy) -> free=5 busy=3 total=8')


def check_measured_free_cores_moves_when_the_busy_reading_moves():
    """The mutation this item exists to prove: inject a DIFFERENT busy
    reading over the SAME cpu set and confirm the free-core count changes
    with it -- a budget that is deaf to tenancy would report the same free
    count regardless of what the sample says."""
    quiet = {i: 0.05 for i in range(16)}
    busy = {i: 0.95 for i in range(16)}
    free_quiet, busy_quiet, _, _ = run_e2e.measured_free_cores(
        total_cores=16, busy_ceiling=0.5, sample_s=0.0,
        busy_fractions_fn=_fake_busy(quiet))
    free_busy, busy_busy, _, _ = run_e2e.measured_free_cores(
        total_cores=16, busy_ceiling=0.5, sample_s=0.0,
        busy_fractions_fn=_fake_busy(busy))
    assert free_quiet == 16 and busy_quiet == 0
    assert free_busy == 0 and busy_busy == 16
    assert free_quiet != free_busy, (
        'free-core count did not move between an all-idle and an all-busy '
        'injected reading -- the default budget is not tenancy-aware')
    print('  PASS  free-core count moves with the injected busy reading '
          '(all-idle free=%d, all-busy free=%d)' % (free_quiet, free_busy))


def check_measured_free_cores_raises_when_no_cpu_stat_can_be_read():
    """Rule 2: a machine this cannot measure (empty {} from the sampler,
    e.g. no /proc/stat) must RAISE, never silently assume idle or full."""
    try:
        run_e2e.measured_free_cores(total_cores=4, busy_ceiling=0.5,
                                    sample_s=0.0, busy_fractions_fn=_fake_busy({}))
    except RuntimeError as exc:
        print('  PASS  measured_free_cores raises when no cpu could be read: %s'
              % exc)
        return
    raise AssertionError('measured_free_cores did not raise on an empty '
                         'busy-fractions reading -- a silent fallback here '
                         'would put a guessed number in front of every cell')


def check_default_budget_uses_measured_free_cores_minus_margin():
    """default_jobs_budget's arithmetic, pinned against a monkeypatched
    measured_free_cores so this does not depend on this box's real load."""
    real_measure = run_e2e.measured_free_cores
    try:
        run_e2e.measured_free_cores = lambda: (40, 24, 64, 0.5)
        args = types.SimpleNamespace(jobs=None)
        cells = [('test.tpv8', 'python-jax')]  # cost 3, well under the budget
        budget = run_e2e.default_jobs_budget(cells, args)
        expected = 40 - run_e2e.JOBS_MARGIN_CORES
        assert budget == expected, (
            'default_jobs_budget returned %r, expected measured-free (40) - '
            'margin (%d) = %d' % (budget, run_e2e.JOBS_MARGIN_CORES, expected))
    finally:
        run_e2e.measured_free_cores = real_measure
    print('  PASS  default_jobs_budget == measured-free - margin (%d) when '
          'that exceeds the largest cell'
          % run_e2e.JOBS_MARGIN_CORES)


def check_default_budget_never_drops_below_the_largest_cell():
    """A box near-saturated (measured-free - margin goes to 0 or negative)
    must still budget at least the single most expensive selected cell, or
    that cell could never run at all."""
    real_measure = run_e2e.measured_free_cores
    try:
        run_e2e.measured_free_cores = lambda: (2, 62, 64, 0.5)  # free=2, margin 4 -> -2
        args = types.SimpleNamespace(jobs=None)
        cells = [('test.drv.a6', 'fortran')]  # cost == matrix.FORTRAN_RANKS
        largest = run_e2e.cell_cost(*cells[0])
        budget = run_e2e.default_jobs_budget(cells, args)
        assert budget == largest, (
            'default_jobs_budget returned %r on a near-saturated box, '
            'expected the largest selected cell\'s own cost %r (never below '
            'it, so every cell can still run)' % (budget, largest))
    finally:
        run_e2e.measured_free_cores = real_measure
    print('  PASS  default_jobs_budget clamps up to the largest cell (%d) '
          'when measured-free - margin would go non-positive' % largest)


def check_explicit_jobs_wins_with_no_measurement_taken():
    """--jobs, when given, must be used verbatim and must not even CALL
    measured_free_cores -- an explicit ask is answered exactly, the same
    discipline run_e2e.py already applies to an explicit --cases/--backends
    selection elsewhere in this file."""
    def _must_not_be_called():
        raise AssertionError('measured_free_cores was called even though '
                             '--jobs was given explicitly')
    real_measure = run_e2e.measured_free_cores
    try:
        run_e2e.measured_free_cores = _must_not_be_called
        args = types.SimpleNamespace(jobs=7)
        cells = [('test.tpv8', 'python-jax')]
        budget = run_e2e.default_jobs_budget(cells, args)
        assert budget == 7, 'explicit --jobs 7 was not honoured: got %r' % budget
    finally:
        run_e2e.measured_free_cores = real_measure
    print('  PASS  explicit --jobs 7 wins outright, no tenancy sample taken')


def check_run_e2e_does_not_use_the_old_fixed_cores_minus_4_default():
    """Source-level mutation guard: the exact pre-fix default,
    `max(1, (os.cpu_count() or 4) - 4)` assigned straight to `budget`, must
    be gone from main(). Its return would silently ignore tenancy again even
    if every function-level check above still passed against the (now
    unused) new functions."""
    path = os.path.join(ROOT, 'testsys', 'e2e', 'run_e2e.py')
    raw = open(path, errors='replace').read()
    src = '\n'.join(ln.split('#', 1)[0] for ln in raw.splitlines())
    bad = 'budget = args.jobs if args.jobs else max(1, (os.cpu_count() or 4) - 4)'
    assert bad not in src, (
        'run_e2e.py still assigns `budget` with the old fixed cores-4 '
        'formula -- the tenancy-aware default_jobs_budget this guard exists '
        'for has been bypassed')
    assert 'budget = default_jobs_budget(cells, args)' in src, (
        'run_e2e.py no longer assigns `budget` from default_jobs_budget(...) '
        '-- the tenancy-aware default may have been removed or bypassed')
    print('  PASS  run_e2e.py assigns budget from default_jobs_budget, not '
          'the old fixed cores-4 formula')


def main():
    print('Regression guard: 2026-09-23 sweep tenancy-aware budget + '
          'honest cell cost')
    checks = [
        check_fortran_cost_is_its_rank_count,
        check_jax_mpi_cost_is_its_rank_count,
        check_python_numpy_backend_is_rejected,
        check_python_jax_cost_is_measured_and_greater_than_one,
        check_mutating_jax_measured_cores_changes_the_billed_cost,
        check_cell_cost_source_does_not_hardcode_python_jax_to_one,
        check_measured_free_cores_counts_correctly,
        check_measured_free_cores_moves_when_the_busy_reading_moves,
        check_measured_free_cores_raises_when_no_cpu_stat_can_be_read,
        check_default_budget_uses_measured_free_cores_minus_margin,
        check_default_budget_never_drops_below_the_largest_cell,
        check_explicit_jobs_wins_with_no_measurement_taken,
        check_run_e2e_does_not_use_the_old_fixed_cores_minus_4_default,
    ]
    failures = []
    for c in checks:
        try:
            c()
        except AssertionError as e:
            failures.append('%s: %s' % (c.__name__, e))
            print('  FAIL  %s: %s' % (c.__name__, e))
    if failures:
        print('\nFAIL test_sweep_tenancy_budget_2026_09_23 (%d check(s))'
              % len(failures))
        return 1
    print('\nSUCCESS test_sweep_tenancy_budget_2026_09_23')
    return 0


if __name__ == '__main__':
    sys.exit(main())
