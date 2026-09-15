"""Unit cover for the two Py-only mechanisms port_jax.py wraps the time loop in:
`time_loop` (traced step count) and `enable_compilation_cache`.

Neither is covered by testsys/parity or testsys/accept. Those tiers gate the
NUMBERS the solver produces, and both mechanisms here are designed to leave
the numbers untouched -- so a regression in either is invisible to them:
`time_loop` running the wrong number of steps would look like a physics
change, and the cache switch degrading to uncached would just make every
process pay ~14 s of XLA compile again with nothing to attribute it to.
PROJECT_RULES rule 2 is the reason the second one must raise rather than
carry on.
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python', 'eqdyna'))

jax = pytest.importorskip('jax', reason='port_jax.py is the JAX backend; without '
                                        'jaxlib there is no time loop to test')
import jax.numpy as jnp  # noqa: E402
import port_jax  # noqa: E402


def _toy_step():
    """A multi-leaf carry whose update depends on the previous value, so a
    wrong trip count or a dropped carry leaf cannot pass by coincidence."""
    def step(carry, _):
        a, b, t = carry
        a = a * 1.5 + b
        b = b + jnp.sum(a)
        t = t + 1.0
        return (a, b, t), None
    return step


@pytest.mark.parametrize('nsteps', [0, 1, 2, 17])
def test_time_loop_matches_scan(nsteps):
    """time_loop must be indistinguishable from the lax.scan it replaced."""
    step = _toy_step()
    c0 = (jnp.arange(5.0), jnp.asarray(0.25), jnp.asarray(0.0))
    want = jax.jit(lambda c: jax.lax.scan(step, c, xs=None, length=nsteps)[0])(c0)
    got = port_jax.time_loop(step, c0, nsteps)
    assert len(got) == len(want)
    for g, w in zip(got, want):
        assert jnp.array_equal(g, w), 'time_loop diverged from lax.scan at nsteps=%d' % nsteps


def test_time_loop_step_count_is_traced_not_baked():
    """The whole point of time_loop: nsteps is a runtime value, not a shape.

    If it ever goes back to `lax.scan(..., length=nsteps)` this raises a
    ConcretizationTypeError, because a traced nsteps cannot be a scan length.
    That regression would silently reintroduce one full XLA compile per
    distinct step count (13.9 s on an A100 for test.tpv8 serial, measured).
    """
    step = _toy_step()
    c0 = (jnp.arange(5.0), jnp.asarray(0.25), jnp.asarray(0.0))
    outer = jax.jit(lambda n: port_jax.time_loop(step, c0, n))
    got = outer(6)
    want = jax.jit(lambda c: jax.lax.scan(step, c, xs=None, length=6)[0])(c0)
    for g, w in zip(got, want):
        assert jnp.array_equal(g, w)


@pytest.fixture(autouse=True)
def _reset_cache_switch():
    """enable_compilation_cache() is once-per-process by design; reset it so
    each test below actually exercises the code path rather than the memo."""
    saved = os.environ.get(port_jax._CACHE_ENV)
    yield
    port_jax._cache_enabled = False
    if saved is None:
        os.environ.pop(port_jax._CACHE_ENV, None)
    else:
        os.environ[port_jax._CACHE_ENV] = saved


def test_cache_dir_unusable_raises(tmp_path, monkeypatch):
    """Rule 2: a named cache directory that cannot be used is a hard failure.

    Falling back to uncached would be a silent ~11 s/process performance
    regression reported as a solver cost.
    """
    blocker = tmp_path / 'i-am-a-file'
    blocker.write_text('')
    monkeypatch.setenv(port_jax._CACHE_ENV, str(blocker / 'cache'))
    port_jax._cache_enabled = False
    with pytest.raises(OSError) as exc:
        port_jax.enable_compilation_cache()
    assert port_jax._CACHE_ENV in str(exc.value), \
        'the error must name the variable the operator has to fix'


def test_cache_off_is_honoured(monkeypatch):
    """`off` must leave JAX's cache dir alone -- it is how a cold compile is
    timed, so it must not quietly enable the cache anyway."""
    monkeypatch.setenv(port_jax._CACHE_ENV, 'off')
    port_jax._cache_enabled = False
    before = jax.config.jax_compilation_cache_dir
    port_jax.enable_compilation_cache()
    assert jax.config.jax_compilation_cache_dir == before


def test_cache_dir_is_configured_when_writable(tmp_path, monkeypatch):
    monkeypatch.setenv(port_jax._CACHE_ENV, str(tmp_path / 'jaxcache'))
    port_jax._cache_enabled = False
    port_jax.enable_compilation_cache()
    assert jax.config.jax_compilation_cache_dir == str(tmp_path / 'jaxcache')
    assert (tmp_path / 'jaxcache').is_dir()
    assert not (tmp_path / 'jaxcache' / '.eqdyna-write-probe').exists(), \
        'the write probe must be cleaned up, not left in the cache dir'
