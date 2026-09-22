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
_cache_path = None


def enable_compilation_cache(subdir=None):
    """Point JAX's on-disk compilation cache at a real directory.

    XLA backend compilation of this time loop costs 13.9 s (GPU) / 4.7 s
    (CPU) against 13.0 ms/step of actual compute, so a 120-step run spends
    90% of its "solve" inside the compiler. Warm, that drops to 6.5 / 2.8 s.

    NO silent fallback: a named directory that cannot be written raises.
    Quietly running uncached would report a compile time nobody could
    account for. EQDYNA_JAX_CACHE_DIR=off runs uncached deliberately.

    `subdir` GIVES THE CALLER ITS OWN CACHE DIRECTORY, and the MPI path uses
    it (driver.run_mpi passes 'rank<n>'). MEASURED REASON, not hygiene: with
    N ranks sharing one cache directory, `test.tpv8` x python-jax-mpi wedged
    intermittently -- three runs over three days lost -- and the same cell
    completed in 12.4 s with the cache off. JAX's persistent cache takes a
    per-key lock inside the cache directory, and the ranks contend on it
    whenever they compile at the same instant (which, after comm.Barrier(),
    is always) or whenever a foreign jax process on the box holds it.

    AND THERE IS NOTHING TO SHARE. Rank r's part_a/part_b are traced against
    rank r's OWN element counts (inv_l['Ei'], inv_l['Ep'], the halo length),
    so every rank compiles a DIFFERENT program and no two ranks could ever hit
    each other's entry. A shared directory buys exactly zero reuse and pays
    the whole lock contention, so splitting it loses nothing measurable. Reuse
    ACROSS runs -- the reuse that actually pays -- is preserved: rank r at the
    same rank count re-traces the same shapes and hits its own warm entry.
    """
    global _cache_enabled, _cache_path
    import os
    if _cache_enabled:
        # Loud on a second, DIFFERENT request. Returning silently would leave
        # the MPI path sharing the directory it asked not to share, which is
        # the wedge this argument exists to remove -- and it would look fixed.
        if subdir is not None and _cache_path != _resolved_cache_path(subdir):
            raise RuntimeError(
                'backend.enable_compilation_cache: the JAX compilation cache '
                'is already pointed at %r, but %r was requested. The first '
                'caller wins in JAX (jax.config is process-global and the '
                'cache is read at first compile), so honouring this silently '
                'would run under the wrong directory. Call this once, before '
                'the first compile.' % (_cache_path, _resolved_cache_path(subdir)))
        return
    import jax
    if os.environ.get(_CACHE_ENV) == 'off':
        _cache_enabled = True
        _cache_path = 'off'
        return
    path = _resolved_cache_path(subdir)
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
    _cache_path = path


def _resolved_cache_path(subdir=None):
    """The directory enable_compilation_cache would use. Separate so the
    already-enabled check above compares the same string the setter stores,
    rather than re-deriving it slightly differently."""
    import os
    want = os.environ.get(_CACHE_ENV)
    if want == 'off':
        return 'off'
    path = want or os.path.join(os.path.expanduser('~'), '.cache', 'eqdyna-jax')
    return os.path.join(path, subdir) if subdir else path


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


# ---------------------------------------------------------------------------
# EXPLICIT DOMAIN DECOMPOSITION (shard_map) -- what MPI does, in jax
#
# Fortran runs one rank per subdomain, assembles its own elements, and calls
# MPI4NodalQuant(nodalForceArr, 3) (driver.f90:27) to sum contributions at
# nodes shared between ranks. This section is that, with N host CPU devices
# standing in for N ranks:
#
#   ranks           -> jax.sharding.Mesh over N CPU devices
#   subdomain       -> a contiguous slab of the ELEMENT arrays, sharded on
#                      their leading axis inside shard_map
#   MPI4NodalQuant  -> nodal_sync() below, a collective at the same point in
#                      the step that driver.f90 calls MPI4NodalQuant
#
# The element kernel is UNCHANGED: assembleGlobalKU.py and faulting.py do not
# know this exists. What the decomposition needs from them is only that every
# element-axis array be sharded consistently, which is what _ELEM_GROUP
# declares, and that the nodal sum happen once per step, which is the one
# line added to driver.step.
#
# WHAT THIS DOES *NOT* DECOMPOSE, and it is the measured ceiling: the NODAL
# stages (velDispUpdate, faulting, the mass divide) stay REPLICATED -- every
# device computes all of them on the full nodal arrays. That is correct
# (every device holds the same post-collective force, so it computes the same
# nodal answer) and it is why the collective can be a plain psum. It also
# caps the speedup at Amdahl's law on the nodal fraction. See
# testsys/perf/run_shard_scaling.py for the measurement, and the session
# report for why the halo variant needs local equation renumbering to beat it.
# ---------------------------------------------------------------------------

SHARD_AXIS = 'd'          # the mesh axis name; 'ranks', spelled for jax
_DEVICES_ENV = 'EQDYNA_JAX_DEVICES'
_NOT_CONFIGURED = object()   # sentinel: _configure_host_devices not yet run;
                             # None is a real cached answer (env var unset),
                             # so it cannot double as "not yet computed".
_host_devices = _NOT_CONFIGURED   # set once by _configure_host_devices at import

# Which leading-axis count each element-axis array in assembleGlobalKU.build's
# dict is indexed by. EVERY key in that dict must appear here or in
# _REPLICATED: _split_sharded raises on an unclassified key rather than
# guessing, because guessing wrong means an array sharded that should be
# replicated (wrong answer, no error) or the reverse (shape error at best).
_ELEM_GROUP = {
    # interior elements (Ei)
    'lam_i': 'Ei', 'miu_i': 'Ei', 'dNx_i': 'Ei', 'dNy_i': 'Ei', 'dNz_i': 'Ei',
    'constk_w_i': 'Ei', 'conn_i': 'Ei', 'idxIx': 'Ei', 'idxIy': 'Ei',
    'idxIz': 'Ei', 'm_e_i': 'Ei', 'stress_i0': 'Ei',
    # PML elements (Ep)
    'lam_p': 'Ep', 'miu_p': 'Ep', 'dNx_p': 'Ep', 'dNy_p': 'Ep', 'dNz_p': 'Ep',
    'det_w_p': 'Ep', 'wx_p': 'Ep', 'wy_p': 'Ep', 'wz_p': 'Ep',
    'neg_det_w_p': 'Ep', 'conn_p': 'Ep', 'a1': 'Ep', 'a2': 'Ep', 'a3': 'Ep',
    'b1': 'Ep', 'b2': 'Ep', 'b3': 'Ep', 'idxP12': 'Ep', 'idxP3': 'Ep',
    'm_e_p': 'Ep', 'pml_init6': 'Ep',
    # all elements (E) -- hourglass control runs over every element
    'conn': 'E', 'phi': 'E', 'ss': 'E', 'idxH0': 'E', 'idxH1': 'E', 'idxH2': 'E',
}

# Nodal or scalar: identical on every device.
_REPLICATED = (
    'N', 'dt', 'rdampk', 'NEQ', 'NEQ1', 'int_nodes_idx', 'idx3_v',
    'pml_nodes_idx', 'idx12_v', 'a9', 'b9', 'Ei', 'Ep', 'E', 'C_elastic',
    'grav_const', 'ccosphi', 'sinphi', 'tv', 'shard_axis', 'shard_sync',
)

# The NODE-axis arrays velDispUpdate loops over. Sharding these as well is a
# MEASUREMENT MODE ONLY (EQDYNA_SHARD_MODE=element+nodal) and it computes the
# WRONG ANSWER: each device would then update only its own slice of the nodal
# arrays, so every device's velArr/dispArr is valid on its own nodes and stale
# elsewhere, while the next step's element kernel gathers velocities at the
# nodes of ITS elements -- which include nodes another device owns. Making it
# correct is exactly the halo problem (own + ghost node lists per device, and
# a boundary exchange instead of a full psum); this mode exists to measure
# what that would BUY before paying for it, and driver.run refuses to write
# output under it. Group lengths come from the arrays themselves, not from a
# second count that could drift.
_NODAL_GROUP = {'int_nodes_idx': 'Nint', 'idx3_v': 'Nint',
                'pml_nodes_idx': 'Npml', 'idx12_v': 'Npml',
                'a9': 'Npml', 'b9': 'Npml'}
_PAD_ONE_NODAL = ('b9',)
MODE_ENV = 'EQDYNA_SHARD_MODE'
SYNC_ENV = 'EQDYNA_SHARD_SYNC'
MODES = ('element', 'element+nodal')
SYNCS = ('psum', 'none')

# Index arrays that build() stores already .ravel()ed to (n*8,). They must be
# (n, 8) to be sharded on the element axis, and are ravelled back to (n_local*8,)
# INSIDE the shard_map body -- so the kernel sees exactly what it sees serially,
# and the scatter's within-device index order is unchanged.
_RAVELLED = ('idxIx', 'idxIy', 'idxIz', 'idxH0', 'idxH1', 'idxH2',
             'idxP12', 'idxP3')

# Padded-element fill. Zero everywhere makes a pad element contribute EXACTLY
# 0.0 (dN=0 -> every block is 0.0, and the scatter target is index 0, the sink
# that calcHourglassResist already scrubs), so padding cannot move a bit.
# EXCEPT the three PML divisors b1/b2/b3: calcPMLElemKU divides by them, and
# 0/0 is NaN, which would propagate through the scatter. They pad with 1.0.
_PAD_ONE = ('b1', 'b2', 'b3') + _PAD_ONE_NODAL


def shard_mode():
    """'element' (default, correct) or 'element+nodal' (timing only)."""
    import os
    m = os.environ.get(MODE_ENV, 'element')
    if m not in MODES:
        raise ValueError('%s=%r: must be one of %r' % (MODE_ENV, m, MODES))
    return m


def shard_sync():
    """'psum' (default, correct) or 'none' (collective removed; timing only).

    'none' exists to ATTRIBUTE the step cost: the difference between the two
    is the collective, measured rather than inferred. It computes the wrong
    answer -- each device keeps only its own elements' partial nodal sums --
    so it is refused for anything that writes results."""
    import os
    s = os.environ.get(SYNC_ENV, 'psum')
    if s not in SYNCS:
        raise ValueError('%s=%r: must be one of %r' % (SYNC_ENV, s, SYNCS))
    return s


def timing_only():
    """True when a knob is set that makes the ANSWER invalid but the TIMING
    valid. Callers that persist results must refuse, or mark, the output."""
    return shard_mode() != 'element' or shard_sync() != 'psum'


def jax_device_count():
    """How many CPU devices the decomposition was asked for. 1 = serial.

    NO FALLBACK: a non-integer or <1 value raises. Silently running serial
    after being asked for 16 devices would put a 1-core number in a scaling
    table under a 16-core label."""
    import os
    raw = os.environ.get(_DEVICES_ENV)
    if raw is None:
        return 1
    n = int(raw)   # ValueError on garbage, deliberately uncaught
    if n < 1:
        raise ValueError('%s=%r: device count must be >= 1' % (_DEVICES_ENV, raw))
    return n


def _configure_host_devices():
    """Turn EQDYNA_JAX_DEVICES=N into N host CPU devices, BEFORE jax imports.

    Returns None when EQDYNA_JAX_DEVICES is UNSET. An absent label is not a
    request for one device -- it makes no claim about device count at all,
    so array_module must not assert against jax.devices() in that case and
    jax is left to initialise on whatever platform it finds (GPU included).
    Returns the validated N (>=1) when the env var IS set; that is an
    explicit claim and stays checked strictly, including N==1.

    An explicit N pins the CPU platform (`JAX_PLATFORMS=cpu`) as well as the
    device count: `--xla_force_host_platform_device_count` only changes how
    many devices the CPU *backend* reports, so on a box with GPUs jax would
    otherwise still hand back the GPU devices and the count downstream would
    be comparing N host-CPU devices against however many GPUs are visible --
    the wrong comparison, and the reason the old assertion's wording ("N host
    CPU devices were requested") was a lie whenever nothing was requested and
    misleading even when something was: it named CPU devices while reporting
    whatever jax actually initialised. A pre-existing JAX_PLATFORMS pin to
    something other than 'cpu' is a second authority for the platform and is
    refused rather than silently overridden, the same stance already taken
    below for a pre-existing XLA_FLAGS.

    `--xla_force_host_platform_device_count` is read by XLA when the CPU
    backend is first INITIALISED (the first jax.devices()), not when jax is
    imported -- so `import jax` having happened already is harmless, but a
    prior jax.devices() is not: the flag becomes a silent no-op and the run
    reports N-core timings for a 1-device execution.

    That case is not detected here (there is no public "is the backend up"
    predicate, and guessing from sys.modules rejected legitimate callers --
    eqdyna3d imports jax before the solver runs). It is caught instead where
    it is exactly checkable: array_module compares jax.devices() against the
    count requested and raises. One check, at the only point where the answer
    is knowable.

    Called ONCE at the bottom of this module, i.e. at `import eqdyna.backend`,
    which is the earliest moment the package can act and is before any
    eqdyna module touches jax (eqdyna3d.Profile initialises the jax backend
    inside run_case, which is why doing this from array_module was too late).
    Idempotent: the second call returns the first answer instead of appending
    the flag twice."""
    import os
    global _host_devices
    if _host_devices is not _NOT_CONFIGURED:
        return _host_devices
    if _DEVICES_ENV not in os.environ:
        _host_devices = None
        return None
    n = jax_device_count()
    cur_platforms = os.environ.get('JAX_PLATFORMS', '')
    if cur_platforms and cur_platforms != 'cpu':
        raise RuntimeError(
            'JAX_PLATFORMS=%r is already set and %s=%d needs JAX_PLATFORMS=cpu '
            'to make the requested host CPU device count the platform jax '
            'actually initialises. Two authorities for the platform is one '
            'too many -- unset JAX_PLATFORMS or set it to cpu.'
            % (cur_platforms, _DEVICES_ENV, n))
    os.environ['JAX_PLATFORMS'] = 'cpu'
    if n == 1:
        _host_devices = 1
        return 1
    cur = os.environ.get('XLA_FLAGS', '')
    if 'xla_force_host_platform_device_count' in cur:
        raise RuntimeError(
            'XLA_FLAGS already sets xla_force_host_platform_device_count (%r) '
            'and %s=%d also asks for it. Two authorities for the device count '
            'is one too many -- pass it exactly one way.' % (cur, _DEVICES_ENV, n))
    os.environ['XLA_FLAGS'] = (cur + ' --xla_force_host_platform_device_count=%d' % n).strip()
    _host_devices = n
    return n


def nodal_sync(xp, inv, arr):
    """driver.f90:27's MPI4NodalQuant(nodalForceArr, 3).

    Serial (shard_axis None): identity, no copy, no collective -- the serial
    path is bit-for-bit what it was before this existed.

    Sharded: an all-reduce over the device mesh. Each device has assembled
    only ITS elements, so a node shared between subdomains holds a partial
    sum on each device that touches it; psum makes every device's nodal array
    the complete one, which is the postcondition every later stage
    (velDispUpdate, faulting, the mass divide) relies on.

    This is an O(NEQ) all-reduce where MPI moves O(boundary). Measured on
    test.tpv104 (NEQ+1 = 3.82e6 doubles = 30.5 MB): 14.8 ms at 4 devices,
    18.0 at 8, 19.4 at 16 -- see testsys/perf/run_shard_scaling.py."""
    axis = inv['shard_axis']
    if axis is None:
        return arr
    if inv['shard_sync'] == 'none':
        return arr                      # timing-only; see shard_sync()
    import jax
    return jax.lax.psum(arr, axis_name=axis)


def _pad_to(a, n, fill):
    """`a` extended along axis 0 to length n with `fill`."""
    if a.shape[0] == n:
        return a
    pad = [(0, n - a.shape[0])] + [(0, 0)] * (a.ndim - 1)
    return np.pad(np.asarray(a), pad, constant_values=fill)


def pad_counts(inv, ndev):
    """{group: padded length} -- each group's count rounded UP to a multiple
    of ndev, because shard_map requires the sharded axis to divide evenly."""
    raw = {'Ei': int(inv['Ei']), 'Ep': int(inv['Ep']), 'E': int(inv['E']),
           'Nint': int(inv['int_nodes_idx'].shape[0]),
           'Npml': int(inv['pml_nodes_idx'].shape[0])}
    return {g: -(-n // ndev) * ndev for g, n in raw.items()}


def _prep(key, v, ndev, counts, groups):
    """One inv value, host-side, ready to be sharded: element-axis arrays
    un-ravelled to (n, 8) where build() had flattened them, then padded."""
    n = counts[groups[key]]
    fill = 1.0 if key in _PAD_ONE else 0
    if key in _RAVELLED:
        v = np.asarray(v).reshape(-1, 8)
    return _pad_to(np.asarray(v), n, fill)


def _split_sharded(inv, ndev, counts, groups):
    """(dynamic, static, specs) for the sharded path.

    Unlike promote(), EVERY array is dynamic here -- an element-axis array
    closed over as an HLO literal would enter the manual region at its GLOBAL
    shape and could not be sharded at all. Raises on any key that is in
    neither table: a new entry in build()'s dict must be classified
    deliberately, not defaulted."""
    from jax.sharding import PartitionSpec as P
    unknown = [k for k in inv if k not in groups and k not in _REPLICATED]
    if unknown:
        raise KeyError(
            'backend._split_sharded: %r appear in assembleGlobalKU.build\'s '
            'invariants but are classified neither element-axis (_ELEM_GROUP) '
            'nor replicated (_REPLICATED). Classify them: sharding an array '
            'that should be replicated gives a wrong answer with no error.'
            % sorted(unknown))
    # `inv` has already been through to_device, so its values are jax Arrays,
    # not np.ndarray -- test for array-ness by dtype+ndim rather than by type,
    # or every array silently lands in `sta` (closed over at its GLOBAL shape)
    # and the kernel gets full-length operands against sharded carry entries.
    def isarr(x):
        return hasattr(x, 'dtype') and getattr(x, 'ndim', 0) >= 1

    dyn, sta, specs = {}, {}, {}
    for k, v in inv.items():
        elem = k in groups
        if isinstance(v, (list, tuple)) and len(v) and all(isarr(x) for x in v):
            dyn[k] = [_prep(k, x, ndev, counts, groups) for x in v] if elem else list(v)
            specs[k] = [P(SHARD_AXIS) if elem else P()] * len(v)
        elif isarr(v):
            dyn[k] = _prep(k, v, ndev, counts, groups) if elem else v
            specs[k] = P(SHARD_AXIS) if elem else P()
        else:
            sta[k] = v
    return dyn, sta, specs


def _restore_ravelled(d):
    """Undo _prep's (n, 8) reshape inside the body, so the kernel receives the
    flat index arrays it receives serially."""
    out = dict(d)
    for k in _RAVELLED:
        if k not in out:
            continue
        v = out[k]
        out[k] = [x.reshape(-1) for x in v] if isinstance(v, list) else v.reshape(-1)
    return out


def run_time_loop_sharded(xp, build_step, inv, carry, n, ndev, carry_shard):
    """run_time_loop's explicitly-decomposed twin. N devices, one element slab
    each, one nodal collective per step at driver.f90:27's position.

    `carry_shard` is a tuple parallel to `carry`, naming the element group of
    each entry that lives on the element axis ('Ei'/'Ep') and None for the
    nodal ones. driver.run states it next to carry0, where the shapes are;
    this module must not guess it from a shape (Ei == Ep is possible)."""
    import jax
    from jax.sharding import Mesh, NamedSharding, PartitionSpec as P
    if len(jax.devices()) != ndev:
        raise RuntimeError(
            'run_time_loop_sharded: asked for %d devices, jax has %d (%r). '
            'The device count must come from EQDYNA_JAX_DEVICES before jax is '
            'imported (see _configure_host_devices) -- a mismatch here means '
            'the flag did not take.' % (ndev, len(jax.devices()), jax.devices()))
    enable_compilation_cache()
    mesh = Mesh(np.array(jax.devices()), (SHARD_AXIS,))
    groups = dict(_ELEM_GROUP)
    if shard_mode() == 'element+nodal':
        groups.update(_NODAL_GROUP)
    counts = pad_counts(inv, ndev)
    dyn, sta, specs = _split_sharded(inv, ndev, counts, groups)
    sta['shard_axis'] = SHARD_AXIS
    sta['shard_sync'] = shard_sync()

    if len(carry_shard) != len(carry):
        raise ValueError('run_time_loop_sharded: carry_shard has %d entries for '
                         'a %d-entry carry' % (len(carry_shard), len(carry)))
    carry_p, cspecs = [], []
    for v, g in zip(carry, carry_shard):
        if g is None:
            carry_p.append(v); cspecs.append(P())
        else:
            carry_p.append(xp.asarray(_pad_to(np.asarray(v), counts[g], 0)))
            cspecs.append(P(SHARD_AXIS))
    carry_p = tuple(carry_p); cspecs = tuple(cspecs)

    def put(v, s):
        return jax.device_put(xp.asarray(v), NamedSharding(mesh, s))

    dyn = {k: ([put(x, s) for x, s in zip(v, specs[k])] if isinstance(v, list)
               else put(v, specs[k])) for k, v in dyn.items()}
    carry_p = tuple(put(v, s) for v, s in zip(carry_p, cspecs))

    def body(dyn_arrays, c, steps):
        step = build_step({**sta, **_restore_ravelled(dyn_arrays)})
        return time_loop(xp, step, c, steps)

    f = jax.jit(jax.shard_map(body, mesh=mesh,
                              in_specs=(specs, cspecs, P()), out_specs=cspecs,
                              check_vma=False))
    out = f(dyn, carry_p, xp.asarray(n))
    jax.block_until_ready(out)
    # Un-pad the element-axis carry entries, so callers see the real element
    # counts and a padded run is indistinguishable from an unpadded one.
    return tuple(v[:int(inv[g])] if g is not None else v
                 for v, g in zip(out, carry_shard))


def array_module(name):
    """'numpy' | 'jax' -> the array module the whole port is written against.

    NO FALLBACK: 'jax' with jaxlib missing raises. A run labelled jax that
    silently WAS numpy makes every backend comparison meaningless, and a
    printed notice does not help -- it scrolls past in a log while the number
    it invalidates is what gets recorded.

    float64 is enabled BEFORE jax.numpy is imported. Silent float32 would
    fake a speedup and break parity against the double-precision Fortran;
    the ordering is the reason this import lives here and not at module top.

    Device count is asserted ONLY when EQDYNA_JAX_DEVICES is explicitly set
    (_configure_host_devices returns None otherwise): an absent label makes
    no claim about device count, so jax is left to initialise on whatever
    platform/device count it finds -- GPUs included, and unconstrained.
    """
    if name == 'numpy':
        return np
    if name != 'jax':
        raise ValueError("backend.array_module: name must be 'numpy' or 'jax' "
                         "(got %r)" % (name,))
    n = _configure_host_devices()
    try:
        import jax
    except ImportError as exc:
        raise RuntimeError(
            "backend 'jax' was requested but jax is not importable (%s). "
            "Install jaxlib, or ask for 'numpy' explicitly. This does NOT "
            "fall back: a run reported as jax must have been jax." % exc)
    jax.config.update('jax_enable_x64', True)
    if n is not None and len(jax.devices()) != n:
        raise RuntimeError(
            '%s=%d host CPU devices were explicitly requested (JAX_PLATFORMS '
            'pinned to cpu) but jax initialised %d device(s) (%r) instead. '
            'Running on a different device count than the label says makes '
            'every scaling number wrong, so this does not proceed.'
            % (_DEVICES_ENV, n, len(jax.devices()), jax.devices()))
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


# Set --xla_force_host_platform_device_count at PACKAGE IMPORT, the earliest
# point available to us and before any eqdyna module initialises the jax
# backend. A no-op (not even an env write) unless EQDYNA_JAX_DEVICES is
# explicitly set, so the default (unset) serial/numpy/jax-on-whatever-jax-
# finds paths are untouched.
_configure_host_devices()
