"""Unit test for testsys/e2e/run_e2e.py's EQDYNA_TEST_LOCK_HELD /
EQDYNA_TEST_ALREADY_ROTATED flags (owner, 2026-09-29): testsys/run.py sets
these when IT already took the test/ lock and rotation at its own start, so
run_e2e.py (invoked as its subprocess) must not repeat either -- but invoked
DIRECTLY (env vars unset), it must still lock and rotate exactly as before.

Both ways (rule 14a): set gives True, unset gives False -- covered for both
flags and both states, by reloading the module fresh under each environment
(the flags are read once, at module import time).
"""
import importlib.util
import os

from conftest import REPO_ROOT

RUN_E2E_PATH = os.path.join(REPO_ROOT, 'testsys', 'e2e', 'run_e2e.py')


def _load_run_e2e(env_overrides):
    """A FRESH module object each call -- these two flags are read once, at
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


def test_flags_default_false_when_unset():
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': None,
                         'EQDYNA_TEST_ALREADY_ROTATED': None})
    assert mod.TEST_LOCK_HELD is False
    assert mod.TEST_ALREADY_ROTATED is False


def test_flags_true_only_on_exact_string_one():
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': '1',
                         'EQDYNA_TEST_ALREADY_ROTATED': '1'})
    assert mod.TEST_LOCK_HELD is True
    assert mod.TEST_ALREADY_ROTATED is True


def test_flags_false_on_any_other_value():
    """Guards against a truthy-string bug (e.g. `bool('0')` is True) --
    only the literal '1' means "the caller already did this"."""
    mod = _load_run_e2e({'EQDYNA_TEST_LOCK_HELD': '0',
                         'EQDYNA_TEST_ALREADY_ROTATED': 'true'})
    assert mod.TEST_LOCK_HELD is False
    assert mod.TEST_ALREADY_ROTATED is False
