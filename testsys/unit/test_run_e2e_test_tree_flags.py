"""Unit tests for testsys/e2e/run_e2e.py's test/-tree coordination with the
calling testsys/run.py (owner, 2026-09-29; PR #56 audit fixing the first
pass at this):

  - EQDYNA_TEST_LOCK_HELD / EQDYNA_RUN_LOG_PATH env-var parsing (both ways).
  - acquire_test_lock's VERIFY-don't-trust behaviour (audit M2): a claimed
    holder is confirmed by a real acquire attempt, never trusted blindly.

TEST_ALREADY_ROTATED / EQDYNA_TEST_ALREADY_ROTATED are gone entirely (audit
M1): run.py no longer rotates test/ early (blocker B1), so there is nothing
left to mark "already rotated" -- run_e2e.py's own rotation is now
unconditional again, exactly as before PR #56.

Both ways (rule 14a): TEST_LOCK_HELD_CLAIMED reads true only on the exact
string '1'; acquire_test_lock takes its own lock when nothing really holds
one (claim or not) and defers when something genuinely does.
"""
import importlib.util
import os

import pytest

from conftest import REPO_ROOT

RUN_E2E_PATH = os.path.join(REPO_ROOT, 'testsys', 'e2e', 'run_e2e.py')

from testsys import runlock  # noqa: E402


@pytest.fixture(autouse=True)
def _clean_environ():
    saved = dict(os.environ)
    yield
    os.environ.clear()
    os.environ.update(saved)


def _load_run_e2e(env_overrides):
    """A FRESH module object each call -- TEST_LOCK_HELD_CLAIMED/RUN_LOG_PATH
    are read once, at import time, so re-importing the cached module would
    not observe an environment change made after the first load."""
    for k, v in env_overrides.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    spec = importlib.util.spec_from_file_location('run_e2e_under_test', RUN_E2E_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_lock_held_claimed_defaults_false_when_unset():
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': None, 'EQDYNA_RUN_LOG_PATH': None})
    assert mod.TEST_LOCK_HELD_CLAIMED is False
    assert mod.RUN_LOG_PATH is None


def test_lock_held_claimed_true_only_on_exact_string_one():
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': '1'})
    assert mod.TEST_LOCK_HELD_CLAIMED is True


def test_lock_held_claimed_false_on_any_other_value():
    """Guards against a truthy-string bug (e.g. `bool('0')` is True) --
    only the literal '1' means "the caller already holds it"."""
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': '0'})
    assert mod.TEST_LOCK_HELD_CLAIMED is False
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': 'true'})
    assert mod.TEST_LOCK_HELD_CLAIMED is False


def test_run_log_path_passes_through_verbatim():
    mod = _load_run_e2e({'EQDYNA_RUN_LOG_PATH': '/tmp/some/staged.log'})
    assert mod.RUN_LOG_PATH == '/tmp/some/staged.log'


# --------------------------------------------------------------------------
# acquire_test_lock (PR #56 audit M2): VERIFY a claimed holder, never trust
# it blindly.
# --------------------------------------------------------------------------
run_e2e = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': None, 'EQDYNA_RUN_LOG_PATH': None})


def test_acquire_test_lock_not_claimed_acquires_for_real(tmp_path):
    lock, refusal = run_e2e.acquire_test_lock(str(tmp_path), claimed_held=False)
    try:
        assert refusal is None
        assert lock is not None
        # A second, real attempt on the same tree now genuinely refuses.
        with pytest.raises(runlock.RunTreeLocked):
            runlock.acquire(str(tmp_path), 'test', announce=False)
    finally:
        lock.release()


def test_acquire_test_lock_not_claimed_refuses_when_already_locked(tmp_path):
    held = runlock.acquire(str(tmp_path), 'test', announce=False)
    try:
        lock, refusal = run_e2e.acquire_test_lock(str(tmp_path), claimed_held=False)
        assert lock is None
        assert refusal is not None and 'test' in refusal
    finally:
        held.release()


def test_acquire_test_lock_claimed_and_genuinely_held_does_not_reacquire(tmp_path):
    """The intended case: testsys/run.py really does hold the lock, and
    claimed_held=True is verified true -- acquire_test_lock takes no lock
    of its own (would deadlock/refuse if it tried) and reports no refusal."""
    held = runlock.acquire(str(tmp_path), 'test', announce=False)
    try:
        lock, refusal = run_e2e.acquire_test_lock(str(tmp_path), claimed_held=True)
        assert lock is None
        assert refusal is None
    finally:
        held.release()


def test_acquire_test_lock_claimed_but_not_actually_held_becomes_real_holder(tmp_path):
    """PR #56 audit M2's exact defect: EQDYNA_TEST_LOCK_HELD=1 is set (a
    stale or leaked env var) but NOBODY actually holds the lock. Blindly
    trusting it would run the sweep with zero protection against a real
    second invocation. acquire_test_lock must instead become the real
    holder -- verified here by confirming a SEPARATE acquire attempt on the
    same tree now genuinely refuses."""
    lock, refusal = run_e2e.acquire_test_lock(str(tmp_path), claimed_held=True)
    try:
        assert refusal is None
        assert lock is not None, (
            'a claimed-but-unverified holder must result in acquire_test_lock '
            'becoming the REAL holder, not silently running unprotected')
        with pytest.raises(runlock.RunTreeLocked):
            runlock.acquire(str(tmp_path), 'test', announce=False)
    finally:
        lock.release()
