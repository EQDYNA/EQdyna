"""The shard_map decomposition was retired 2026-09-28 (PR #47). Its knobs
must be REFUSED, not silently ignored: a leftover EQDYNA_JAX_DEVICES would
otherwise run serial under a multi-device label."""
import os
import sys

import pytest

from conftest import REPO_ROOT

sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))
from eqdyna import backend as B  # noqa: E402


@pytest.mark.parametrize('key', ['EQDYNA_JAX_DEVICES', 'EQDYNA_SHARD_MODE',
                                 'EQDYNA_SHARD_SYNC'])
@pytest.mark.parametrize('name', ['numpy', 'jax'])
def test_retired_knob_is_refused(monkeypatch, key, name):
    monkeypatch.setenv(key, '4')
    with pytest.raises(RuntimeError, match='retired 2026-09-28'):
        B.array_module(name)


def test_clean_env_still_returns_numpy(monkeypatch):
    for k in list(os.environ):
        if k == 'EQDYNA_JAX_DEVICES' or k.startswith('EQDYNA_SHARD_'):
            monkeypatch.delenv(k)
    assert B.array_module('numpy') is B.np
