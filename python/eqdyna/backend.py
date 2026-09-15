"""The ONLY backend-aware module in the port.

Everything else -- driver.py, faulting.py, fric.py, assembleGlobalKU.py --
is written once against an array module `xp` (numpy or jax.numpy) and these
five helpers. That is the whole point: Fortran has one driver.f90 and one
faulting.f90, and every place this port grew a second copy is a place a fix
has to be made, and verified, twice. (It did: faulting.f90's swtwNucleation
forced-rupture branch was fixed in port.py and not port_jax.py, so
test.tpv29 passed on one backend and failed on the other at the identical
value.)

THE IN-PLACE CONTRACT
`scatter_add` takes the destination `out` and RETURNS it. On numpy it
mutates in place and hands back the same buffer; on jax it returns a new
array. Callers must use the return value and must not assume either. This
signature is what lets the unified kernel keep numpy's block-at-a-time
accumulation into one live `force` -- the shape that took test.tpv104's
peak RSS from 4.16 GB to 2.52 GB by never materialising the 148-entries-
per-element index/value function space.

WHAT THIS MODULE IS NOT
It answers exactly one question -- "how does THIS array library spell the
same operation" -- and nothing else. It knows no physics: there is no
hourglass mode, no fault node, no friction law in here. The two places where
the two backends genuinely compute a DIFFERENT answer are hourglass-kernel
decisions and live in assembleGlobalKU.py, next to the arithmetic they
change; an adapter that knew about hourglass modes would be the six port*
modules again, relocated.

One divergence is not under our control at all: XLA lowers a duplicate-index
scatter-add to atomics, so jax-gpu is already nondeterministic run-to-run
(measured 8.9e-08 on unmodified code). Any claim of bit-identity in this
port is a claim about CPU.
"""
import numpy as np


def is_jax(xp):
    """True when `xp` is jax.numpy. Checked by module name, not by import:
    importing jax to answer this would pay a jaxlib import for every numpy
    run (and fail on a numpy-only install)."""
    return xp.__name__ != 'numpy'


def scatter_add(xp, out, idx, val):
    """out[idx] += val with DUPLICATE indices accumulated. Returns `out`.

    numpy accumulates in index-array order; XLA's CPU scatter applies a
    duplicate index list in array order. Called once per block, in the same
    block order, the two therefore perform the same additions on the same
    partial sums -- which is load-bearing, because float addition is not
    associative and this port is gated on bit-identity.
    """
    if is_jax(xp):
        return out.at[idx].add(val)
    np.add.at(out, idx, val)
    return out


def setat(xp, a, idx, v):
    """a[idx] = v. Returns `a`."""
    if is_jax(xp):
        return a.at[idx].set(v)
    a[idx] = v
    return a


def addat(xp, a, idx, v):
    """a[idx] += v for UNIQUE indices (no duplicate-index accumulation
    intended). Separate from scatter_add so the two intents read
    differently at the call site."""
    if is_jax(xp):
        return a.at[idx].add(v)
    a[idx] += v
    return a


def time_loop(xp, step, carry, n):
    """Run `step(carry, nt)` for nt = 1..n and return the final carry.

    numpy: a plain Python for. jax: lax.fori_loop with a TRACED trip count,
    so one compiled executable serves every step count -- lax.scan's static
    `length` put nsteps into the cache key and made every distinct run
    length pay a full XLA compile.
    """
    if is_jax(xp):
        import jax
        return jax.lax.fori_loop(0, n, lambda i, c: step(c, i + 1), carry)
    for nt in range(1, n + 1):
        carry = step(carry, nt)
    return carry


def iadd(xp, a, b):
    """a + b, reusing `a`'s storage on numpy. Returns the sum.

    `a` MUST be a private temporary the caller owns -- this mutates it.

    Why it exists: `x = p*q; x = x + r; x = x + s` and `x = p*q; x += r;
    x += s` are the SAME left-associative sum on the same two roundings, so
    they are bit-identical; they differ only in whether numpy allocates a
    fresh block per `+`. The element-force blocks here are (E, 8) float64 --
    47 MB each on test.tpv104 -- and the kernel builds ~20 of them per step,
    so writing the unified kernel in pure functional style would hand numpy
    an allocation it did not previously make. On jax this is a plain `+`
    (arrays are immutable and XLA fuses the chain anyway).

    This is the helper that makes ONE kernel possible without the numpy
    column regressing from 2.52 GB back toward 4.16 GB on test.tpv104.
    """
    if is_jax(xp):
        return a + b
    a += b
    return a


# ---------------------------------------------------------------------------
# jax device/compile plumbing
#
# This section is why backend.py is ~200 lines and not the ~40 the shape of
# the four core helpers suggests. Every line of it is a MEASURED memory or
# compile-time result on test.tpv104, and none of it has anywhere else to
# live once the two kernels are one: it is exactly the backend-specific part
# that unification pushes down here instead of duplicating up there.
# ---------------------------------------------------------------------------

# Loop-invariant float arrays promoted from HLO constant to jit ARGUMENT.
# NOT applied to float arrays wholesale -- exactly this set was verified
# bit-identical end to end; see promote() for the counter-example that stops
# the list here.
_PROMOTED_FLOAT = ('phi', 'ss', 'dNx_i', 'dNy_i', 'dNz_i',
                   'dNx_p', 'dNy_p', 'dNz_p')

_I32_MAX = 2 ** 31 - 1
_CACHE_ENV = 'EQDYNA_JAX_CACHE_DIR'
_cache_enabled = False


def enable_compilation_cache():
    """Point JAX's on-disk compilation cache at a real directory.

    XLA backend compilation of this time loop costs 13.9 s (GPU) / 4.7 s
    (CPU) against 13.0 ms/step of actual compute, so a 120-step run spends
    90% of its "solve" inside the compiler. Warm, that drops to 6.5 / 2.8 s.

    NO silent fallback: a named directory that cannot be written raises.
    Quietly running uncached would report a compile time nobody could
    account for. EQDYNA_JAX_CACHE_DIR=off runs uncached deliberately.
    """
    global _cache_enabled
    if _cache_enabled:
        return
    import os
    import jax
    want = os.environ.get(_CACHE_ENV)
    if want == 'off':
        _cache_enabled = True
        return
    path = want or os.path.join(os.path.expanduser('~'), '.cache', 'eqdyna-jax')
    try:
        os.makedirs(path, exist_ok=True)
        probe = os.path.join(path, '.eqdyna-write-probe')
        with open(probe, 'w') as fh:
            fh.write('')
        os.remove(probe)
    except OSError as exc:
        raise OSError(
            'JAX persistent compilation cache directory %r is unusable (%s). '
            'Point %s at a writable directory, or set it to "off" to run '
            'uncached. No fallback: running uncached without saying so adds '
            '~14 s of XLA compile with nothing to attribute it to.'
            % (path, exc, _CACHE_ENV))
    jax.config.update('jax_compilation_cache_dir', path)
    jax.config.update('jax_persistent_cache_min_compile_time_secs', 0.5)
    jax.config.update('jax_persistent_cache_min_entry_size_bytes', 0)
    _cache_enabled = True


def to_device(xp, inv):
    """Move a host-built invariants dict onto the backend.

    numpy: returned unchanged -- it is already exactly what the kernel wants.

    jax: every array to device, with INDEX arrays narrowed to int32. On
    test.tpv104 the index arrays come to ~660 MB as int64 and are pure
    addressing -- every one is consumed by a gather or a scatter, never by
    float arithmetic -- so a narrower index type cannot change a single
    rounding. Verified end to end, digests unchanged.

    NO SILENT TRUNCATION: int32 addresses up to 2**31-1, and a mesh big
    enough to exceed that would wrap into a valid-looking but WRONG index --
    a silently corrupted scatter, the worst failure available here. Checked
    once, loudly, below.
    """
    if not is_jax(xp):
        return inv

    def conv(v):
        if isinstance(v, np.ndarray):
            if v.dtype.kind in 'iub':
                return xp.asarray(v, dtype=xp.int32) if v.dtype.kind != 'b' \
                    else xp.asarray(v)
            return xp.asarray(v)
        if isinstance(v, (list, tuple)) and len(v) and \
                all(isinstance(x, np.ndarray) for x in v):
            return [conv(x) for x in v]
        return v

    return {k: conv(v) for k, v in inv.items()}


def check_index_width(inv):
    """Refuse a mesh whose indices do not fit the int32 narrowing to_device
    applies on the jax path.

    int32 addresses up to 2**31-1. A mesh big enough to exceed that in
    equation or node count would WRAP AROUND into a valid-looking but wrong
    index -- a silently corrupted scatter, the worst failure available here.
    Checked once, loudly, against the mesh invariants, rather than trusted.
    """
    biggest = max(inv['NEQ1'], inv['N'], inv['E'] * 8)
    if biggest > _I32_MAX:
        raise OverflowError(
            'backend.check_index_width: this mesh needs more than int32 '
            'indices (max index magnitude %d > %d: NEQ+1=%d, N=%d, E*8=%d). '
            'Widen to_device\'s index narrowing for a mesh this size rather '
            'than letting the cast wrap around silently.'
            % (biggest, _I32_MAX, inv['NEQ1'], inv['N'], inv['E'] * 8))


def promote(xp, inv):
    """Split `inv` into (dynamic, static) for jit.

    dynamic -- passed as a jit ARGUMENT, so XLA sees a PARAMETER and the
    array stays exactly one device buffer: every integer array, plus the
    large float mesh arrays in _PROMOTED_FLOAT.
    static  -- closed over: every Python scalar the kernels branch on at
    TRACE time (C_elastic, Ep, friclaw, TPV, ...), plus the remaining small
    coefficient arrays.

    WHY. An array CLOSED OVER by a jitted function becomes an HLO *literal*,
    and a literal is paid for many times over: in the traced jaxpr's consts,
    again in the lowered MLIR module text, again when XLA parses and
    constant-folds it, and once more in the executable, which holds its own
    copy alongside the caller's still-live original. Measured on test.tpv104
    (inv is 1.03 GB of arrays), closing all of it over:
        module text 1200 MB; XLA compile RSS 3.25 -> 7.09 GB, high-water
        10.17 GB, 5.9 s; jax.live_arrays 1.15 -> 2.31 GB (the duplicate)
    while XLA's own memory_analysis reports arg 0.141 + out 0.141 + temp
    1.072 = 1.35 GB. The peak-RSS gap was the COMPILER holding copies of
    mesh constants, not the solver holding state. With this split:
        module text 1200 -> 171.5 MB; compile 10.18 GB/6.01 s -> 3.04 GB/
        2.87 s; peak RSS full run 9.93 -> 4.01 GB; 569.2 -> 530.8 ms/step.

    WHY THE LIST STOPS WHERE IT DOES -- the load-bearing part. Promoting an
    array removes XLA's opportunity to fold it into the expressions it feeds,
    and folding is not guaranteed association-preserving. Every promotion was
    therefore gated on sha256 over tobytes() of the WHOLE output at full step
    count, never allclose and never a spot check:
      * integers + _PROMOTED_FLOAT  -> BIT-IDENTICAL on tpv8/tpv104/tpv10/
        drv.a6. This is the landed configuration.
      * EVERY array promoted        -> NOT bit-identical (fric 2.09e-07).
        Inside this port's own nondeterminism floor, but a real change, and
        not taken, because the configuration above gets the memory anyway.
    A new array added to `inv` therefore defaults to STATIC. To promote it,
    add it here and re-run that digest check.
    """
    def isint(x):
        return hasattr(x, 'dtype') and x.dtype.kind in 'iu'

    dyn, sta = {}, {}
    for k, v in inv.items():
        if isint(v) or (k in _PROMOTED_FLOAT and hasattr(v, 'dtype')):
            dyn[k] = v
        elif isinstance(v, (list, tuple)) and len(v) and all(isint(x) for x in v):
            dyn[k] = list(v)
        else:
            sta[k] = v
    return dyn, sta


def run_time_loop(xp, build_step, inv, carry, n):
    """Drive `n` steps, with `inv` reaching the step the way each backend
    needs it. This is what driver.py calls; `time_loop` above is the
    primitive it is built on.

    Takes the step FACTORY rather than an already-built closure, so that on
    jax the step's `inv` lookups resolve against jit ARGUMENTS instead of
    closed-over constants -- see promote() for the measured reason. Building
    the step inside the jitted body is what makes that possible, and it is
    the single reason this helper exists rather than driver.py calling
    time_loop directly.
    """
    if not is_jax(xp):
        return time_loop(xp, build_step(inv), carry, n)

    import jax
    enable_compilation_cache()
    dyn, sta = promote(xp, inv)

    def body(dyn_arrays, c, steps):
        step = build_step({**sta, **dyn_arrays})
        return time_loop(xp, step, c, steps)

    out = jax.jit(body)(dyn, carry, n)
    jax.block_until_ready(out)
    return out


def array_module(name):
    """'numpy' | 'jax' -> the array module the whole port is written against.

    NO FALLBACK: 'jax' with jaxlib missing raises. A run labelled jax that
    silently WAS numpy makes every backend comparison meaningless, and a
    printed notice does not help -- it scrolls past in a log while the number
    it invalidates is what gets recorded.

    float64 is enabled BEFORE jax.numpy is imported. Silent float32 would
    fake a speedup and break parity against the double-precision Fortran;
    the ordering is the reason this import lives here and not at module top.
    """
    if name == 'numpy':
        return np
    if name != 'jax':
        raise ValueError("backend.array_module: name must be 'numpy' or 'jax' "
                         "(got %r)" % (name,))
    try:
        import jax
    except ImportError as exc:
        raise RuntimeError(
            "backend 'jax' was requested but jax is not importable (%s). "
            "Install jaxlib, or ask for 'numpy' explicitly. This does NOT "
            "fall back: a run reported as jax must have been jax." % exc)
    jax.config.update('jax_enable_x64', True)
    import jax.numpy as jnp
    return jnp


def mul_into(xp, out, a, b):
    """a * b, written into `out` on numpy; a plain `a * b` on jax.

    Same ufunc, same operands, same rounding -- `np.multiply(a, b, out=out)`
    and `out = a * b` differ only in WHERE the result lands, so this is
    bit-identical, not approximate.

    Why it exists: the element-force kernel builds ~30 blocks of shape (E, 8)
    per step and scatters each one immediately. Written functionally, numpy
    allocates and frees every one of them -- on test.tpv8 that is ~250 MB of
    malloc/munmap churn per step, and blocks that size go straight to mmap,
    so the cost is page faults rather than allocator bookkeeping. Measured:
    reusing one buffer per block shape is worth 3-5% of the step.

    `out` is ignored on jax, where arrays are immutable and XLA picks its own
    buffers; callers pass None there and must use the RETURN value either way.
    """
    if is_jax(xp):
        return a * b
    return np.multiply(a, b, out=out)


def store_into(xp, out, src):
    """Copy `src` into `out` on numpy; return `src` on jax.

    Needed where a running sum is SEEDED from a block that is itself a
    reused buffer: aliasing the two would let the next block overwrite the
    accumulator. Values are unchanged either way.
    """
    if is_jax(xp):
        return src
    np.copyto(out, src)
    return out
