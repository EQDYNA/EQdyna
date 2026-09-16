"""Unit cover for the three Py-only mechanisms backend.py wraps the time loop
in: `time_loop` (traced step count), `promote` (loop-invariant arrays passed
to jit as ARGUMENTS rather than closed-over HLO constants) and
`enable_compilation_cache`.

None is covered by the e2e sweep. That gates the
NUMBERS the solver produces, and all three mechanisms here are designed to
leave the numbers untouched -- so a regression in any of them is invisible to
it: `time_loop` running the wrong number of steps would look like a physics
change, the cache switch degrading to uncached would just make every process
pay ~14 s of XLA compile again with nothing to attribute it to, and
`promote` silently moving one more array into the promoted set would be a
peak-RSS win that quietly costs bit-identity (measured: promoting every array
moves test.tpv8's fric by 2.09e-07 -- far inside every accept bound, i.e.
exactly the kind of drift a tolerance gate cannot see). PROJECT_RULES rule 2
is the reason the cache switch must raise rather than carry on.
"""
import os
import sys

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src', 'python'))

jax = pytest.importorskip('jax', reason='these cover the JAX backend; without '
                                        'jaxlib there is no time loop to test')
jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp  # noqa: E402
from eqdyna import backend as B  # noqa: E402


def _toy_inv():
    """An `inv`-shaped dict covering all three classes promote sorts on: an
    integer index array, a float array named in _PROMOTED_FLOAT, a float array
    that is NOT, a list of integer arrays, and plain Python scalars."""
    return dict(
        idxIx=jnp.asarray([0, 2, 1, 2]),
        phi=jnp.asarray([1.5, 2.5, 3.5, 4.5]),        # in _PROMOTED_FLOAT
        det_w_p=jnp.asarray([0.25, 0.5, 0.75, 1.0]),  # NOT in _PROMOTED_FLOAT
        idxP3=[jnp.asarray([1, 0, 2, 1]), jnp.asarray([2, 2, 0, 0])],
        NEQ1=3, dt=0.125, C_elastic=1,
    )


def _toy_build_step():
    """Factory (not a pre-built closure) -- the shape time_loop now takes, so
    the step's `inv` lookups resolve against jit arguments. The carry update
    reads every class of inv entry and depends on the previous value, so a
    wrong trip count, a dropped carry leaf, or an inv entry lost by the
    dyn/sta split cannot pass by coincidence."""
    def build_step(inv):
        def step(carry, _):
            a, b, t = carry
            a = a * 1.5 + b + inv['phi'] * inv['dt'] + inv['det_w_p']
            a = a + jnp.zeros(inv['NEQ1']).at[inv['idxIx']].add(inv['phi'])[inv['idxP3'][0]]
            b = b + jnp.sum(a) * inv['C_elastic']
            t = t + 1.0
            return (a, b, t)
        return step
    return build_step


def _toy_carry():
    return (jnp.arange(4.0), jnp.asarray(0.25), jnp.asarray(0.0))


@pytest.mark.parametrize('nsteps', [0, 1, 2, 17])
def test_time_loop_matches_scan(nsteps):
    """time_loop must be indistinguishable from the lax.scan it replaced --
    including now that `inv` crosses the jit boundary as an argument."""
    build_step = _toy_build_step()
    inv = _toy_inv()
    c0 = _toy_carry()
    def scan_step(c, x):
        return build_step(inv)(c, x), None

    want = jax.jit(lambda c: jax.lax.scan(scan_step, c, xs=None, length=nsteps)[0])(c0)
    got = B.run_time_loop(jnp, build_step, inv, c0, nsteps)
    assert len(got) == len(want)
    for g, w in zip(got, want):
        assert jnp.array_equal(g, w), 'run_time_loop diverged from lax.scan at nsteps=%d' % nsteps


def test_time_loop_step_count_is_traced_not_baked():
    """The whole point of time_loop: nsteps is a runtime value, not a shape.

    If it ever goes back to `lax.scan(..., length=nsteps)` this raises a
    ConcretizationTypeError, because a traced nsteps cannot be a scan length.
    That regression would silently reintroduce one full XLA compile per
    distinct step count (13.9 s on an A100 for test.tpv8 serial, measured).
    """
    build_step = _toy_build_step()
    inv = _toy_inv()
    c0 = _toy_carry()
    got = B.run_time_loop(jnp, build_step, inv, c0, 6)

    def scan_step(c, x):
        return build_step(inv)(c, x), None

    want = jax.jit(lambda c: jax.lax.scan(scan_step, c, xs=None, length=6)[0])(c0)
    for g, w in zip(got, want):
        assert jnp.array_equal(g, w)


def test_time_loop_passes_the_1_based_step_number():
    """time_loop hands the step its 1-BASED step number, and for port_tp_jax
    that value IS physics, not bookkeeping: thermop writes the history column
    `nt-1` and evaluates its convolution kernel at age `(nt-j)*dt`
    (src/updateThermalPressurization.f90:22-33). It must be exactly the
    sequence `lax.scan(step, c, xs=jnp.arange(1, nsteps+1))` used to deliver
    before all three JAX ports were put on this one loop -- an off-by-one here
    would shift every past term's weight by one step and read as a physics
    change, not as a bug."""
    def build_step(inv):
        def step(carry, nt):
            seen, = carry
            return (seen.at[nt - 1].set(nt.astype(seen.dtype)),)
        return step

    n = 5
    got, = B.run_time_loop(jnp, build_step, _toy_inv(), (jnp.zeros(n),), n)
    assert jnp.array_equal(got, jnp.arange(1.0, n + 1)), \
        'the loop delivered %r, not the 1..nsteps sequence' % (got,)


def test_promote_promotes_int32_index_arrays_too():
    """to_device() hands every index array to the device as int32 now
    (~660 MB of addressing on test.tpv104, halved). promote classifies on
    dtype KIND, not width, so int32 must still land in `dyn` -- if a narrower
    index fell through to `static` it would become an HLO literal again and the
    footprint win would silently revert, with every gate still green."""
    dyn, sta = B.promote(jnp, dict(
        idxIx=jnp.asarray([0, 2, 1], dtype=jnp.int32),
        idxP3=[np.asarray([1, 0, 2], dtype=np.int32)],
        det_w_p=jnp.asarray([0.5, 1.0, 1.5]), dt=0.125))
    assert 'idxIx' in dyn and 'idxP3' in dyn, 'int32 index arrays must be promoted'
    assert 'det_w_p' in sta


def test_promote_is_a_partition():
    """Every `inv` key must land in exactly one side. A key dropped by the
    split surfaces as a loud KeyError in the kernel, but a key that MOVES from
    static to dynamic is the silent case: it changes what XLA may fold and can
    move the answer in the last bits with every gate still green."""
    inv = _toy_inv()
    dyn, sta = B.promote(jnp, inv)
    assert set(dyn) | set(sta) == set(inv)
    assert not (set(dyn) & set(sta))


def test_promote_promotes_only_the_verified_set():
    """The promoted set is an ALLOWLIST established by end-to-end digest
    comparison (see promote's docstring), not "every array". Integer index
    arrays and lists of them are promoted; a float array outside
    _PROMOTED_FLOAT stays a constant, precisely so its folding -- which
    test.tpv8/tpv104/tpv10/drv.a6 bit-identity was verified WITH -- is
    preserved."""
    inv = _toy_inv()
    dyn, sta = B.promote(jnp, inv)
    assert 'idxIx' in dyn and 'idxP3' in dyn, 'integer index arrays must be promoted'
    assert 'phi' in dyn, 'phi is in the verified _PROMOTED_FLOAT set'
    assert 'det_w_p' in sta, \
        'a float array outside _PROMOTED_FLOAT must stay a constant -- promoting it ' \
        'is a bit-identity change that needs its own end-to-end digest check first'
    assert 'det_w_p' not in B._PROMOTED_FLOAT


def test_promote_keeps_scalars_static_for_trace_time_branching():
    """The kernels branch on `inv['C_elastic']`/`inv['Ep']` with a plain
    Python `if` at trace time. A scalar promoted into the argument pytree
    would become a tracer and raise TracerBoolConversionError -- so scalars
    must stay static, and this pins that."""
    inv = _toy_inv()
    dyn, sta = B.promote(jnp, inv)
    for k in ('NEQ1', 'dt', 'C_elastic'):
        assert k in sta and not isinstance(sta[k], (jax.Array, np.ndarray))
    invd = {**sta, **dyn}
    assert bool(invd['C_elastic'] == 1), 'a static scalar must remain usable in a Python if'


def test_promote_accepts_numpy_and_jax_arrays_alike():
    """build() hands over jnp arrays, but the same dict is built from NumPy
    upstream; both must classify identically or the promoted set would depend
    on where the array came from."""
    dyn_j, sta_j = B.promote(jnp, dict(idxIx=jnp.asarray([0, 1]), phi=jnp.asarray([1.0]),
                                            det_w_p=jnp.asarray([2.0]), dt=0.5))
    dyn_n, sta_n = B.promote(jnp, dict(idxIx=np.asarray([0, 1]), phi=np.asarray([1.0]),
                                            det_w_p=np.asarray([2.0]), dt=0.5))
    assert set(dyn_j) == set(dyn_n) and set(sta_j) == set(sta_n)


@pytest.fixture(autouse=True)
def _reset_cache_switch():
    """enable_compilation_cache() is once-per-process by design; reset it so
    each test below actually exercises the code path rather than the memo."""
    saved = os.environ.get(B._CACHE_ENV)
    yield
    B._cache_enabled = False
    if saved is None:
        os.environ.pop(B._CACHE_ENV, None)
    else:
        os.environ[B._CACHE_ENV] = saved


def test_cache_dir_unusable_raises(tmp_path, monkeypatch):
    """Rule 2: a named cache directory that cannot be used is a hard failure.

    Falling back to uncached would be a silent ~11 s/process performance
    regression reported as a solver cost.
    """
    blocker = tmp_path / 'i-am-a-file'
    blocker.write_text('')
    monkeypatch.setenv(B._CACHE_ENV, str(blocker / 'cache'))
    B._cache_enabled = False
    with pytest.raises(OSError) as exc:
        B.enable_compilation_cache()
    assert B._CACHE_ENV in str(exc.value), \
        'the error must name the variable the operator has to fix'


def test_cache_off_is_honoured(monkeypatch):
    """`off` must leave JAX's cache dir alone -- it is how a cold compile is
    timed, so it must not quietly enable the cache anyway."""
    monkeypatch.setenv(B._CACHE_ENV, 'off')
    B._cache_enabled = False
    before = jax.config.jax_compilation_cache_dir
    B.enable_compilation_cache()
    assert jax.config.jax_compilation_cache_dir == before


def test_cache_dir_is_configured_when_writable(tmp_path, monkeypatch):
    monkeypatch.setenv(B._CACHE_ENV, str(tmp_path / 'jaxcache'))
    B._cache_enabled = False
    B.enable_compilation_cache()
    assert jax.config.jax_compilation_cache_dir == str(tmp_path / 'jaxcache')
    assert (tmp_path / 'jaxcache').is_dir()
    assert not (tmp_path / 'jaxcache' / '.eqdyna-write-probe').exists(), \
        'the write probe must be cleaned up, not left in the cache dir'


def test_numpy_time_loop_matches_the_same_contract():
    """The SAME step contract driven by the NUMPY branch of time_loop must
    produce the same trip count and the same 1-based nt sequence. This is what
    the unification actually claims -- one loop, two backends -- and nothing
    else in the suite compares the two loop DRIVERS directly."""
    seen = []

    def step(carry, nt):
        seen.append(nt)
        return carry + nt

    assert B.time_loop(np, step, 0, 5) == 15
    assert seen == [1, 2, 3, 4, 5]


def test_check_index_width_refuses_a_mesh_too_big_for_int32():
    """to_device narrows every index to int32. A mesh that exceeds 2**31-1
    would WRAP AROUND into a valid-looking but wrong index -- a silently
    corrupted scatter. It must raise, loudly, naming the numbers."""
    with pytest.raises(OverflowError) as exc:
        B.check_index_width(dict(NEQ1=3, N=3, E=2 ** 29))
    assert 'int32' in str(exc.value)
    B.check_index_width(dict(NEQ1=3, N=3, E=10))
