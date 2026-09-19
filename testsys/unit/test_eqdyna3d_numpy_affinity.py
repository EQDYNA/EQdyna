"""Unit cover for `eqdyna3d._narrow_numpy_affinity` (pathway item 33 Phase 3).

MEASURED, not assumed (see the function's own docstring and
`testsys/perf/run_scaling.py`): the numpy backend's hot kernels run on ONE
thread regardless of core count, and a caller-given affinity mask that spans
more than one NUMA node can make it MEASURABLY WORSE (single-threaded OS
migration between nodes stranding memory on the wrong one) -- 9443/9872/8900
ms/step flat across 1/2/4 cpus within one node, 19055 ms/step (~2x worse) at
16 cpus spanning two nodes, test.tpv104, this box. This is the regression
test for the fix, not an end-to-end timing assertion (timing is not
reproducible enough to gate on -- see PROJECT_RULES rule 5); it exercises the
AFFINITY-NARROWING behaviour itself, which is deterministic.

Not covered by the e2e sweep: every e2e cell runs single-process with
whatever affinity the test runner happens to have, and narrowing it to one
cpu changes wall time, never the solver's answer -- exactly the kind of
change the sweep cannot see (same reasoning as test_backend_jax.py's own
module docstring).
"""
import os

import pytest

from eqdyna import eqdyna3d as E

pytestmark = pytest.mark.skipif(
    not hasattr(os, 'sched_getaffinity'),
    reason='os.sched_getaffinity does not exist on this platform (e.g. '
           'macOS) -- _narrow_numpy_affinity is a documented no-op there.')


@pytest.fixture(autouse=True)
def _restore_affinity():
    """Every test in this file changes process-wide affinity; restore it so
    a test order dependency cannot leak into an unrelated later test (in this
    file or, since affinity is process-global, any other file run in the
    same pytest process)."""
    before = os.sched_getaffinity(0)
    yield
    os.sched_setaffinity(0, before)


def _cpu_and_a_sibling():
    """Two distinct real cpu ids present on this machine, for a mask test
    needs to actually be > 1 wide. Skips (does not fake a result) if this
    machine genuinely has only one cpu available to the process."""
    avail = sorted(os.sched_getaffinity(0))
    if len(avail) < 2:
        pytest.skip('fewer than 2 cpus available to this process on this '
                    'machine -- nothing to narrow')
    return avail[0], avail[1]


def test_narrows_a_wide_mask_to_exactly_one_cpu_chosen_by_pid():
    """CORRECTED (wei-lin, gate-axis-3 review): the first version always
    picked the LOWEST cpu in the mask, verified directly to make two
    concurrently-launched numpy processes both pin to cpu 0 -- fine alone,
    but it serialises testsys/e2e/run_e2e.py's own concurrent-cell sweep
    onto one physical core. The chosen cpu is now PID-dependent (still
    deterministic FOR a given process, still exactly one cpu -- the
    NUMA-migration fix is unchanged), so this asserts the actual formula
    rather than assuming 'lowest', which would make this test itself
    PID-order-dependent and occasionally wrong."""
    lo, hi = _cpu_and_a_sibling()
    os.sched_setaffinity(0, {lo, hi})
    os.environ.pop(E._NUMPY_WIDE_AFFINITY_ENV, None)
    E._narrow_numpy_affinity()
    ordered = sorted((lo, hi))
    expected = {ordered[os.getpid() % len(ordered)]}
    result = os.sched_getaffinity(0)
    assert result == expected
    assert len(result) == 1
    assert result.issubset({lo, hi})


def test_two_processes_narrow_to_different_cpus_when_pids_differ_in_parity():
    """The actual regression this correction fixes: two concurrent processes
    must NOT both land on the identical cpu just because they were given the
    identical starting mask. Simulated in-process by monkeypatching
    os.getpid rather than forking (deterministic, no process-launch cost),
    covering both parities of a 2-wide mask."""
    lo, hi = _cpu_and_a_sibling()
    chosen = []
    real_getpid = os.getpid
    for fake_pid in (2 * real_getpid(), 2 * real_getpid() + 1):
        os.sched_setaffinity(0, {lo, hi})
        os.environ.pop(E._NUMPY_WIDE_AFFINITY_ENV, None)
        os.getpid = lambda fp=fake_pid: fp
        try:
            E._narrow_numpy_affinity()
            chosen.append(next(iter(os.sched_getaffinity(0))))
        finally:
            os.getpid = real_getpid
    assert chosen[0] != chosen[1], (
        'two processes with different-parity pids and the same starting '
        'mask must land on different cpus, not the same one (%r)' % chosen)


def test_leaves_an_already_single_cpu_mask_alone():
    lo, _hi = _cpu_and_a_sibling()
    os.sched_setaffinity(0, {lo})
    os.environ.pop(E._NUMPY_WIDE_AFFINITY_ENV, None)
    E._narrow_numpy_affinity()
    assert os.sched_getaffinity(0) == {lo}


def test_override_env_var_disables_narrowing():
    lo, hi = _cpu_and_a_sibling()
    os.sched_setaffinity(0, {lo, hi})
    os.environ[E._NUMPY_WIDE_AFFINITY_ENV] = '1'
    try:
        E._narrow_numpy_affinity()
        assert os.sched_getaffinity(0) == {lo, hi}
    finally:
        del os.environ[E._NUMPY_WIDE_AFFINITY_ENV]


def test_run_case_only_narrows_for_the_numpy_backend(monkeypatch):
    """The jax backend must NOT be narrowed -- it measurably benefits from
    extra cores (unlike numpy; see the function's docstring). Patches
    `_narrow_numpy_affinity` itself rather than running a real case: this is
    a dispatch test (was the guard called for the right backend), not a
    physics test."""
    calls = []
    monkeypatch.setattr(E, '_narrow_numpy_affinity', lambda: calls.append(True))

    def fake_build_solver_state(case_dir):
        raise _StopEarly

    class _StopEarly(Exception):
        pass

    monkeypatch.setattr(E, 'build_solver_state', fake_build_solver_state)

    for backend in ('numpy', 'jax'):
        calls.clear()
        try:
            E.run_case('unused', nsteps=1, verbose=False, backend=backend)
        except _StopEarly:
            pass
        assert calls == ([True] if backend == 'numpy' else [])
