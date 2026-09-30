"""Unit tests for testsys/e2e/run_e2e.py's test/-tree coordination with the
calling testsys/run.py (owner, 2026-09-29; PR #56 audit, two passes:
2026-09-29 and 2026-09-30):

  - EQDYNA_RUN_LOG_PATH env-var parsing.
  - acquire_test_lock ALWAYS takes the lock itself, unconditionally.

EQDYNA_TEST_LOCK_HELD / TEST_LOCK_HELD_CLAIMED are gone entirely (audit M2,
second finding, 2026-09-30): a first fix VERIFIED a claimed holder (a real
acquire attempt) rather than trusting it blindly, but that was still wrong
-- a failed acquire proves SOME process holds the lock, not that it is safe
to proceed unprotected riding on it, since that process could be an
unrelated concurrent sweep (rule 21a's exact incident). Since testsys/run.py
never legitimately holds this lock itself either (see
run.py._prepare_test_tree's docstring), there was no remaining legitimate
case for the flag, so the whole mechanism -- the env var, the
claimed_held parameter, and the two tests that exercised it
(test_acquire_test_lock_claimed_and_genuinely_held_does_not_reacquire and
test_acquire_test_lock_claimed_but_not_actually_held_becomes_real_holder) --
is deleted. TEST_ALREADY_ROTATED / EQDYNA_TEST_ALREADY_ROTATED were already
gone (audit M1): run.py no longer rotates test/ early (blocker B1), so
run_e2e.py's own rotation is unconditional again, exactly as before PR #56.

Both ways (rule 14a): acquire_test_lock takes its own lock when nothing
holds one and refuses when something genuinely does -- and it refuses even
when EQDYNA_TEST_LOCK_HELD=1 is set in its own environment, since that env
var is no longer read at all.
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
    """A FRESH module object each call -- RUN_LOG_PATH is read once, at
    import time, so re-importing the cached module would not observe an
    environment change made after the first load."""
    for k, v in env_overrides.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    spec = importlib.util.spec_from_file_location('run_e2e_under_test', RUN_E2E_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_run_log_path_defaults_none_when_unset():
    mod = _load_run_e2e({'EQDYNA_RUN_LOG_PATH': None})
    assert mod.RUN_LOG_PATH is None


def test_run_log_path_passes_through_verbatim():
    mod = _load_run_e2e({'EQDYNA_RUN_LOG_PATH': '/tmp/some/staged.log'})
    assert mod.RUN_LOG_PATH == '/tmp/some/staged.log'


def test_run_e2e_no_longer_defines_a_lock_held_flag_at_all():
    """The whole EQDYNA_TEST_LOCK_HELD mechanism is deleted, not merely
    unset -- a stray os.environ['EQDYNA_TEST_LOCK_HELD'] = '1' must have NO
    attribute on the module to read it back from."""
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': '1'})
    assert not hasattr(mod, 'TEST_LOCK_HELD_CLAIMED')


# --------------------------------------------------------------------------
# acquire_test_lock: always takes the lock itself, unconditionally.
# --------------------------------------------------------------------------
run_e2e = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': None, 'EQDYNA_RUN_LOG_PATH': None})


def test_acquire_test_lock_acquires_for_real(tmp_path):
    lock, refusal = run_e2e.acquire_test_lock(str(tmp_path))
    try:
        assert refusal is None
        assert lock is not None
        # A second, real attempt on the same tree now genuinely refuses.
        with pytest.raises(runlock.RunTreeLocked):
            runlock.acquire(str(tmp_path), 'test', announce=False)
    finally:
        lock.release()


def test_acquire_test_lock_refuses_when_already_locked(tmp_path):
    held = runlock.acquire(str(tmp_path), 'test', announce=False)
    try:
        lock, refusal = run_e2e.acquire_test_lock(str(tmp_path))
        assert lock is None
        assert refusal is not None and 'test' in refusal
    finally:
        held.release()


def test_acquire_test_lock_refuses_when_already_locked_even_with_stale_env_flag(
        tmp_path, monkeypatch):
    """PR #56 audit M2's second finding, made concrete: even if
    EQDYNA_TEST_LOCK_HELD=1 is (wrongly) present in THIS process's own
    environment when acquire_test_lock runs, it changes nothing --
    acquire_test_lock takes no env-var-driven shortcut of any kind any
    more, so a second invocation refuses exactly the same as it would with
    the flag absent."""
    monkeypatch.setenv('EQDYNA_TEST_LOCK_HELD', '1')
    held = runlock.acquire(str(tmp_path), 'test', announce=False)
    try:
        lock, refusal = run_e2e.acquire_test_lock(str(tmp_path))
        assert lock is None
        assert refusal is not None and 'test' in refusal
    finally:
        held.release()
