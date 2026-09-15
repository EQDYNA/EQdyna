"""
JAX rewrite of the EQdyna tpv8 (friclaw=1) time loop -- phase 2 of the
feasibility spike. Mirrors python/eqdyna/port.py's NumPy `run()` exactly
(same formulas, same order, same known-bug reproduction for the PML
region cascade), but structured as one jitted, scanned step function
instead of a per-step Python-dispatched loop.

float64 is enabled FIRST, before any other jax import touches an array,
per the phase-2 requirement -- silent float32 would fake a speedup and
break parity against the double-precision Fortran/NumPy reference.
"""
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from eqdyna import nucleation
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from loading import load, region_damp  # noqa: E402  (must follow the x64 config line)
import kernels_jax

_CACHE_ENV = 'EQDYNA_JAX_CACHE_DIR'
_cache_enabled = False


def enable_compilation_cache():
    """Point JAX's persistent (on-disk) compilation cache at a real directory.

    Measured on test.tpv8 serial, 120 steps, jax 0.6.2, A100-SXM4-40GB pinned
    via CUDA_VISIBLE_DEVICES, 64-core EPYC 7532: XLA backend compilation of
    this time loop costs 13.9 s on the GPU and 4.7 s on the CPU, against
    13.0 ms/step (GPU) of actual compute. A 120-step run therefore spends 90%
    of its "solve" inside the compiler. With this cache warm, backend compile
    drops to 6.5 s (GPU) / 2.8 s (CPU) for every process after the first.

    Directory: $EQDYNA_JAX_CACHE_DIR, else ~/.cache/eqdyna-jax. Set
    EQDYNA_JAX_CACHE_DIR=off to run uncached (e.g. to time a cold compile).

    NO silent fallback (PROJECT_RULES rule 2): if a directory is named but
    cannot be created or written, this raises. Quietly running uncached would
    report a compile time nobody could account for -- exactly the confusion
    that made the 143 ms/step baseline unreadable in the first place.

    Cache entries are keyed by XLA on the lowered HLO, the backend and the
    jaxlib version, so a code change cannot be served a stale executable.
    """
    global _cache_enabled
    if _cache_enabled:
        return
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
            'Point %s at a writable directory, or set %s=off to run uncached. '
            '(No fallback: running uncached without saying so would add ~14 s '
            'of XLA compile to the solve phase with nothing to attribute it to.)'
            % (path, exc, _CACHE_ENV, _CACHE_ENV))
    jax.config.update('jax_compilation_cache_dir', path)
    jax.config.update('jax_persistent_cache_min_compile_time_secs', 0.5)
    jax.config.update('jax_persistent_cache_min_entry_size_bytes', 0)
    _cache_enabled = True


# Loop-invariant float arrays promoted from HLO constant to jit ARGUMENT by
# split_inv(). Promotion is NOT applied to float arrays wholesale: it is
# applied to exactly this set, because exactly this set was VERIFIED
# bit-identical end-to-end (see split_inv's docstring for the measurement and
# for the counter-example that stops the list here).
_PROMOTED_FLOAT = ('phi', 'ss', 'dNx_i', 'dNy_i', 'dNz_i', 'dNx_p', 'dNy_p', 'dNz_p')


def split_inv(inv):
    """Split build()'s `inv` into (dynamic, static) for jit.

    dynamic -- handed to `jax.jit` as an ARGUMENT, so XLA sees a PARAMETER and
    the array stays exactly one device buffer:
      * every INTEGER array (and list of integer arrays, e.g. idxP12): mesh
        connectivity, equation ids, the element-force scatter index. These
        feed gathers and a scatter, never float arithmetic, so promoting them
        cannot re-associate a floating-point expression.
      * the large float mesh arrays in _PROMOTED_FLOAT above.
    static -- stays closed over: every Python/NumPy scalar (N, NEQ1, dt,
    rdampk, Ei/Ep/E, C_elastic, grav_const, slipRateThres, ccosphi/sinphi/tv),
    which the kernels branch on at TRACE time (`if inv['C_elastic'] == 0:`,
    `if inv['Ep'] > 0:`), plus every remaining per-element/per-node
    coefficient array (lam/miu, constk_w_i, det_w_p, the PML damping ratios
    a1..b3/a9/b9, inv_mass_full, un/us/ud/arn, m_e_*, stress_i0, pml_init6).

    WHY promote at all (2026-09-15, peak-RSS work): an array CLOSED OVER by a
    jitted function becomes an HLO *literal*, and a literal is paid for many
    times over -- once in the traced jaxpr's consts, again rendered into the
    lowered MLIR module, again when XLA parses and constant-folds it, and once
    more in the executable, which then holds its OWN copy alongside the
    caller's still-live original. An isolated repro: 600 MB of closed-over
    arrays produced a 1200 MB module text and a 5.0 GB compile high-water,
    against 1.2 GB total when the same arrays were passed as arguments.
    Measured on test.tpv104 serial (E=735000, NEQ=3818583, jax-cpu; `inv` is
    1.03 GB of arrays), closing ALL of it over cost:
        module text    1200 MB
        jit lower      +810 MB resident
        XLA compile    RSS 3.25 -> 7.09 GB, high-water 10.17 GB, 5.9 s
        after execute  jax.live_arrays() 1.15 GB -> 2.31 GB (the duplicate)
    while XLA's own memory_analysis() for this loop reports only arg 0.141 +
    out 0.141 + temp 1.072 = 1.35 GB. The peak-RSS gap was the COMPILER
    holding copies of mesh constants, not the solver holding state.

    With this split, same case, same machine:
        module text          1200 MB  ->   171.5 MB
        XLA compile          10.18 GB high-water, 6.01 s -> 3.04 GB, 2.87 s
        peak RSS, full run    9.93 GB ->  4.01 GB   (/usr/bin/time -v, 120
                              steps, serial, incl. setup and frt write)
        solve                 90.4 s  ->  64.1 s
        ms/step, compile excl  569.2  ->  530.8     (-6.8%)
    It is a win on BOTH axes: nothing was traded for the footprint. Same
    change on test.tpv8 1.81 GB (was 3.43) / 229.0 ms/step (was 269.7),
    test.tpv10 3.98 GB (was 10.18), test.drv.a6 9.10 GB (was 23.77).

    WHY THE LIST STOPS WHERE IT DOES -- this is the load-bearing part.
    Promoting an array removes XLA's opportunity to fold it into the
    expressions it feeds (and, for the scatter index, to schedule the
    duplicate-index element-force accumulation the way it does when the index
    is a literal). Neither folding nor that schedule is guaranteed
    association-preserving, and float addition is not associative, so every
    promotion is a potential bit-identity change and was gated on a
    full-output end-to-end check -- sha256 over tobytes() of the WHOLE
    velArr/dispArr/force/fric/fnft, full step count, jax-cpu, never a spot
    check and never allclose:
      * integers + _PROMOTED_FLOAT -> BIT-IDENTICAL on every case this
        split serves: test.tpv8 (114 steps), test.tpv104, test.tpv10,
        test.drv.a6 (120/132/120 steps) -- all five arrays, ZERO differing
        elements, 0 fnft rupture-time flips. This is the landed configuration.
      * EVERY array promoted -> NOT bit-identical: velArr 1.53e-14, dispArr
        1.11e-15, force 3.36e-13, fric 2.09e-07 max abs, 0 fnft flips
        (test.tpv8). That is inside this port's own documented
        nondeterminism floor -- see kernels_jax.py's hourglass-fusion note,
        whose accepted fric figure is the same 2.09e-07 -- but it is a real
        change to the answer, and it is not taken, because the configuration
        above gets the memory without it.
      * everything except a1..b3/a9/b9/det_w_p -> also NOT bit-identical, so
        the responsible constant is one of the remaining small coefficient
        arrays, not the PML damping divisors; not narrowed further, because
        the bit-identical configuration above already meets the footprint
        target.
    A new array added to `inv` therefore defaults to STATIC (constant). To
    promote it, add it here and re-run that end-to-end digest check -- never
    on the strength of a spot check or an allclose.

    ALSO TRIED AND DROPPED (2026-09-15), recorded so it is not re-attempted
    blind: hoisting kernels_jax's per-step `jnp.concatenate` of the 18
    scatter-index blocks into build() as one host-side array. It is
    bit-identical on port_jax/port_rsf_jax, but it was MEMORY-NEUTRAL
    (test.tpv104 arg/temp unchanged, run peak 3.92 vs 3.97 GB -- the 18
    blocks and their concatenation are the same total bytes), and it is NOT
    bit-identical on port_tp_jax (see that module's run()). A change that
    buys no memory and costs bit-identity on one of five cases is not worth
    its diff.

    A scalar is never silently promoted to a traced array: it lands in
    `static`, and if a kernel then tries to branch on a traced value JAX
    raises TracerBoolConversionError -- loud, at the line, never a silent
    re-interpretation.
    """
    def _isint(x):
        return isinstance(x, (jax.Array, np.ndarray)) and x.dtype.kind in 'iu'

    dyn, sta = {}, {}
    for k, v in inv.items():
        if _isint(v) or (k in _PROMOTED_FLOAT and isinstance(v, (jax.Array, np.ndarray))):
            dyn[k] = v
        elif isinstance(v, (list, tuple)) and len(v) and all(_isint(x) for x in v):
            dyn[k] = list(v)
        else:
            sta[k] = v
    return dyn, sta


def time_loop(build_step, inv, carry0, nsteps):
    """Run `build_step(inv)`'s step `nsteps` times, with `nsteps` TRACED
    rather than baked in.

    `lax.scan` needs a static `length`, which put the step count into the
    compiled executable's cache key: every distinct step count paid a full
    XLA compile, and no on-disk cache entry could ever be reused across runs
    of different length. `lax.fori_loop` with a traced trip count lowers to
    the same XLA while-loop and measured identical per step -- 12.905 ms/step
    (scan) vs 12.965 ms/step (fori), GPU, tpv8 serial, 120 steps -- while one
    executable now serves every nsteps (verified: n=500 reused the n=120 cache
    entry, 7.15 s vs 13.9 s cold).

    Takes the step FACTORY (not an already-built closure) plus `inv`, so that
    the step's `inv` lookups resolve against jit ARGUMENTS rather than
    closed-over constants -- see split_inv() for the measured reason.

    The step is called as `step(carry, nt)` with nt the 1-BASED step number,
    the same value `lax.scan(step, c, xs=jnp.arange(1, nsteps+1))` used to
    deliver -- port_jax/port_rsf_jax ignore it, port_tp_jax needs it (it indexes
    the thermal-pressurization history column and evaluates the convolution
    kernel's age at `nt - j`). Passing it here is what lets ALL THREE JAX ports
    share this one loop, which is the point: the Fortran has exactly one time
    loop (driver.f90:10) with the friclaw branch INSIDE it (faulting.f90:17-18),
    and every place this port grew a second loop instead is a place a fix has to
    be made, and verified, twice.

    A caller whose carry contains an array shaped by nsteps (port_tp_jax.py's
    (nftnd, nsteps) history) still recompiles per step count -- the SHAPE is in
    the cache key even though the trip count is not -- and must pass a trip
    count no larger than that width, since `.at[:, nt-1].set()` would otherwise
    CLAMP silently. port_tp_jax allocates the history from the same `nsteps` it
    passes here, so the two cannot disagree; see the assertion there.
    """
    dyn, sta = split_inv(inv)

    def body(dyn_arrays, c, n):
        step = build_step({**sta, **dyn_arrays})
        return jax.lax.fori_loop(0, n, lambda i, cc: step(cc, i + 1)[0], c)

    return jax.jit(body)(dyn, carry0, nsteps)


def build(S):
    """Precompute every loop-invariant array ONCE (Python/NumPy side, static
    shapes), then hand them to a closure that lax.scan iterates. Nothing in
    the scanned step function does a nonzero()/boolean-mask reshape -- all
    such shape-determining operations happen here, outside jit, exactly like
    the NumPy port's pre-loop setup."""
    N = S['N']; dt = S['dt']; w = S['w']; rdampk = S['rdampk']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']

    is_int = elemType == 1
    E_int = np.nonzero(is_int)[0]
    E_pml = np.nonzero(elemType == 2)[0]

    int_nodes_idx = np.nonzero(ndof == 3)[0]
    idx3_v = eq_ids[int_nodes_idx, 0:3]

    pml_nodes_idx = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes_idx, 0], S['meshCoor'][pml_nodes_idx, 1],
                                 S['meshCoor'][pml_nodes_idx, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'])
    dampv_pml = np.zeros((pml_nodes_idx.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes_idx, :]
    a9 = 1.0 / dt - dampv_pml / 2.0; b9 = 1.0 / dt + dampv_pml / 2.0

    lam_i = mat[E_int, 3]; miu_i = mat[E_int, 4]
    dN_i = eleshp[E_int]
    dNx_i, dNy_i, dNz_i = dN_i[:, :, 0], dN_i[:, :, 1], dN_i[:, :, 2]
    constk_w_i = (-eledet[E_int]) * w
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel(); idxIy = eq_ids[conn_i, 1].ravel(); idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p, dNy_p, dNz_p = dN_p[:, :, 0], dN_p[:, :, 1], dN_p[:, :, 2]
    det_w_p = eledet[E_pml] * w
    conn_p = conn[E_pml]
    xc_p = S['meshCoor'][conn_p].mean(axis=1)
    d1p, d2p, d3p = region_damp(xc_p[:, 0], xc_p[:, 1], xc_p[:, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'])
    a1 = 1.0 / dt - d1p / 2.0; b1 = 1.0 / dt + d1p / 2.0
    a2 = 1.0 / dt - d2p / 2.0; b2 = 1.0 / dt + d2p / 2.0
    a3 = 1.0 / dt - d3p / 2.0; b3 = 1.0 / dt + d3p / 2.0
    nd_p = ndof[conn_p]; is12_p = nd_p == 12; is3_p = nd_p == 3
    idxP12 = [np.where(is12_p, eq_ids[conn_p, j], 0).ravel() for j in range(12)]
    idxP3 = [np.where(is3_p, eq_ids[conn_p, g], 0).ravel() for g in range(3)]

    nd_all = ndof[conn]
    slot0 = np.where(nd_all == 3, 0, 9); slot1 = np.where(nd_all == 3, 1, 10); slot2 = np.where(nd_all == 3, 2, 11)
    idxH0 = np.take_along_axis(eq_ids[conn], slot0[:, :, None], axis=2)[:, :, 0].ravel()
    idxH1 = np.take_along_axis(eq_ids[conn], slot1[:, :, None], axis=2)[:, :, 0].ravel()
    idxH2 = np.take_along_axis(eq_ids[conn], slot2[:, :, None], axis=2)[:, :, 0].ravel()

    nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    idxF_s = [eq_ids[nsmp1, d] for d in range(3)]
    idxF_m = [eq_ids[nsmp2, d] for d in range(3)]

    mass_full = np.concatenate(([1.0], S['nodalMassArr']))
    inv_mass_full = np.where(mass_full > 0, 1.0 / np.where(mass_full > 0, mass_full, 1.0), 0.0)

    # Milestone 10 (drv.a6, C_elastic==0) -- see kernels_numpy.py's top
    # docstring for the full derivation; EXACTLY 0.0-valued/no-op for every
    # C_elastic==1 caller (port_jax.py/port_tp_jax.py, and port_rsf_jax.py
    # for tpv104/tpv10).
    C_elastic_static = S['C_elastic']
    grav_const = ((1.0 - C_elastic_static) * S['grav']
                  * (S['roumax'] - (S['gamar'] + 1.0) * S['rhow']) / S['roumax'])
    m_e_all = mat[:, 2] * eledet
    m_e_i = m_e_all[E_int]
    m_e_p = m_e_all[E_pml]
    if C_elastic_static == 0:
        init_stress = S['init_stress']  # (E,6); raises (KeyError) loudly if missing.
        stress_i0 = init_stress[E_int].copy()
        pml_init6 = init_stress[E_pml].copy()
        ccosphi = S['ccosphi']; sinphi = S['sinphi']; tv = S['tv']
    else:
        stress_i0 = np.zeros((E_int.shape[0], 6))
        pml_init6 = np.zeros((E_pml.shape[0], 6))
        ccosphi = sinphi = tv = 0.0

    j = jnp.asarray

    # Every array below that is an INDEX (mesh connectivity, equation ids, the
    # element-force scatter targets, the fault node lists) goes to the device as
    # int32 rather than int64. There are a lot of them -- on test.tpv104
    # (735,000 elements, 763,537 nodes) conn/conn_i/conn_p + idxIx/idxIy/idxIz +
    # the twelve idxP12 + idxP3 + idxH0/1/2 + idx12_v/idx3_v come to ~660 MB as
    # int64 -- and they are pure addressing: every one of them is consumed by a
    # gather or a scatter, never by float arithmetic, so a narrower index type
    # cannot change a single rounding. Halving them removes ~330 MB from `inv`
    # (which split_inv then passes as jit ARGUMENTS, so this is 330 MB off both
    # the host arrays and the device buffers) and halves the bytes each gather/
    # scatter has to read to find its targets. Verified end-to-end, sha256 over
    # tobytes() of the full velArr/dispArr/force/fric/fnft at full step count:
    # digests UNCHANGED on test.tpv8 and test.tpv104, jax-cpu.
    #
    # NO SILENT TRUNCATION (rule 2/3): int32 addresses up to 2**31-1, and a mesh
    # big enough to exceed that in equation count or node count would wrap
    # around into a valid-looking but WRONG index -- a silently corrupted
    # scatter, the worst possible failure here. It is checked, once, loudly,
    # right here instead.
    _I32_MAX = 2 ** 31 - 1
    _biggest = max(NEQ + 1, N, conn.shape[0] * 8)
    if _biggest > _I32_MAX:
        raise OverflowError(
            'port_jax.build: this mesh needs more than int32 indices (max index '
            'magnitude %d > %d: NEQ+1=%d, N=%d, E*8=%d). The index arrays below are '
            'int32 to halve their footprint; widen them back to int64 for a mesh '
            'this size rather than letting the cast wrap around silently.'
            % (_biggest, _I32_MAX, NEQ + 1, N, conn.shape[0] * 8))

    def ji(x):
        """Index array -> int32 on device. Pure addressing, never arithmetic."""
        return jnp.asarray(x, dtype=jnp.int32)

    inv = dict(
        N=N, dt=dt, rdampk=rdampk, NEQ1=NEQ + 1,
        conn=ji(conn), phi=j(phi), ss=j(ss),
        int_nodes_idx=ji(int_nodes_idx), idx3_v=ji(idx3_v),
        pml_nodes_idx=ji(pml_nodes_idx), idx12_v=ji(idx12_v), a9=j(a9), b9=j(b9),
        lam_i=j(lam_i), miu_i=j(miu_i), dNx_i=j(dNx_i), dNy_i=j(dNy_i), dNz_i=j(dNz_i),
        constk_w_i=j(constk_w_i), conn_i=ji(conn_i), idxIx=ji(idxIx), idxIy=ji(idxIy),
        idxIz=ji(idxIz),
        lam_p=j(lam_p), miu_p=j(miu_p), dNx_p=j(dNx_p), dNy_p=j(dNy_p), dNz_p=j(dNz_p),
        det_w_p=j(det_w_p), a1=j(a1), b1=j(b1), a2=j(a2), b2=j(b2), a3=j(a3), b3=j(b3),
        idxP12=[ji(x) for x in idxP12], idxP3=[ji(x) for x in idxP3],
        idxH0=ji(idxH0), idxH1=ji(idxH1), idxH2=ji(idxH2),
        nsmp1=ji(nsmp1), nsmp2=ji(nsmp2), un=j(S['un']), us=j(S['us']), ud=j(S['ud']),
        arn=j(S['arn']),
        idxF_s=[ji(x) for x in idxF_s], idxF_m=[ji(x) for x in idxF_m],
        inv_mass_full=j(inv_mass_full), slipRateThres=S['slipRateThres'], C_elastic=S['C_elastic'],
        Ei=E_int.shape[0], Ep=E_pml.shape[0], E=conn.shape[0],
        grav_const=grav_const, m_e_i=j(m_e_i), m_e_p=j(m_e_p),
        stress_i0=j(stress_i0), pml_init6=j(pml_init6),
        ccosphi=ccosphi, sinphi=sinphi, tv=tv,
    )
    # swtwNucleation's forced-rupture time: geometry only, so it belongs
    # with the other loop invariants. None when this case does no forced
    # nucleation, which make_step branches on at TRACE time.
    if nucleation.enabled(S):
        import numpy as _np
        _r = nucleation.source_radius(_np, S['meshCoor'][S['nsmp1']],
                                      S['xsource'], S['ysource'], S['zsource'])
        inv['nuc_tr'] = jnp.asarray(nucleation.forced_rupture_time(
            _np, _r, S.get('TPV', 0), S.get('nucR', 0.0),
            S.get('nucRuptVel', 0.0)))
    else:
        inv['nuc_tr'] = None

    return inv


def make_step(inv):
    # Forced-rupture time is geometry-only, so it is computed ONCE here rather
    # than per step. None means this case does no forced nucleation.
    nuc_tr = inv.get('nuc_tr')
    dt = inv['dt']; rdampk = inv['rdampk']

    def step(carry, _):
        v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
        timeElapsed = timeElapsed + dt

        v1, velArr, dispArr, force, stress_i, s_p = kernels_jax.elastic_step(
            inv, v1, velArr, dispArr, force, stress_i, s_p, dt, rdampk)

        # ---- faulting (friclaw==1: solveSWTW; no Newton-Raphson -- that
        # branch (RSF, friclaw>=3) is not part of tpv8's executed path, so
        # nothing here resisted jitting) ----
        nsmp1, nsmp2 = inv['nsmp1'], inv['nsmp2']
        un, us, ud, arn = inv['un'], inv['us'], inv['ud'], inv['arn']
        idxF_s, idxF_m = inv['idxF_s'], inv['idxF_m']
        fx_s = force[idxF_s[0]]; fy_s = force[idxF_s[1]]; fz_s = force[idxF_s[2]]
        fx_m = force[idxF_m[0]]; fy_m = force[idxF_m[1]]; fz_m = force[idxF_m[2]]
        vx_s, vy_s, vz_s = velArr[nsmp1, 0], velArr[nsmp1, 1], velArr[nsmp1, 2]
        vx_m, vy_m, vz_m = velArr[nsmp2, 0], velArr[nsmp2, 1], velArr[nsmp2, 2]
        dx_s, dy_s, dz_s = dispArr[nsmp1, 0], dispArr[nsmp1, 1], dispArr[nsmp1, 2]
        dx_m, dy_m, dz_m = dispArr[nsmp2, 0], dispArr[nsmp2, 1], dispArr[nsmp2, 2]

        def rot(vx, vy, vz, dirvec):
            return vx * dirvec[:, 0] + vy * dirvec[:, 1] + vz * dirvec[:, 2]

        fN_s, fS_s, fD_s = rot(fx_s, fy_s, fz_s, un), rot(fx_s, fy_s, fz_s, us), rot(fx_s, fy_s, fz_s, ud)
        fN_m, fS_m, fD_m = rot(fx_m, fy_m, fz_m, un), rot(fx_m, fy_m, fz_m, us), rot(fx_m, fy_m, fz_m, ud)
        vN_s, vS_s, vD_s = rot(vx_s, vy_s, vz_s, un), rot(vx_s, vy_s, vz_s, us), rot(vx_s, vy_s, vz_s, ud)
        vN_m, vS_m, vD_m = rot(vx_m, vy_m, vz_m, un), rot(vx_m, vy_m, vz_m, us), rot(vx_m, vy_m, vz_m, ud)
        dN_s, dS_s, dD_s = rot(dx_s, dy_s, dz_s, un), rot(dx_s, dy_s, dz_s, us), rot(dx_s, dy_s, dz_s, ud)
        dN_m, dS_m, dD_m = rot(dx_m, dy_m, dz_m, un), rot(dx_m, dy_m, dz_m, us), rot(dx_m, dy_m, dz_m, ud)

        slipN = dN_m - dN_s; slipS = dS_m - dS_s; slipD = dD_m - dD_s
        srN = vN_m - vN_s; srS = vS_m - vS_s; srD = vD_m - vD_s
        srMag = jnp.sqrt(srN ** 2 + srS ** 2 + srD ** 2)

        fric = fric.at[:, 70].set(slipS); fric = fric.at[:, 71].set(slipD); fric = fric.at[:, 72].set(slipN)
        fric = fric.at[:, 73].set(srS); fric = fric.at[:, 74].set(srD)
        fric = fric.at[:, 75].set(jnp.maximum(fric[:, 75], srMag))
        fric = fric.at[:, 76].add(srMag * dt)

        massSlave = inv['fnms'][nsmp1] if False else inv['massSlave']
        massMaster = inv['massMaster']
        totalMass = (massSlave + massMaster) * arn
        C_elastic = inv['C_elastic']
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic
        Tmag = jnp.sqrt(Ts ** 2 + Td ** 2)

        fs = fric[:, 0]; fd = fric[:, 1]; D0 = fric[:, 2]
        slip = fric[:, 76]
        fricCoeff = jnp.where(jnp.abs(slip) < 1.0e-10, fs, fs - (fs - fd) * slip / D0)
        fricCoeff = jnp.where(slip >= D0, fd, fricCoeff)
        # swtwNucleation -- the SAME code path numpy uses (eqdyna.nucleation),
        # not a second copy. Absent entirely until now, which is why
        # test.tpv29 nucleated 0 of 3321 nodes: see nucleation.py's docstring.
        if nuc_tr is None:
            fricCoeff = jnp.minimum(fs, fricCoeff)
        else:
            fricCoeff = nucleation.apply(jnp, fricCoeff, fs, fd,
                                         fric[:, 4], nuc_tr, timeElapsed)

        effNorm = jnp.where((Tn + fric[:, 5]) > 0.0, 0.0, Tn + fric[:, 5])
        trialShear = fric[:, 3] - fricCoeff * effNorm
        over = Tmag > trialShear
        scale = jnp.where(over, trialShear / jnp.where(Tmag == 0, 1.0, Tmag), 1.0)
        Ts = Ts * scale; Td = Td * scale

        xN = Tn * un[:, 0] + Ts * us[:, 0] + Td * ud[:, 0]
        xNy = Tn * un[:, 1] + Ts * us[:, 1] + Td * ud[:, 1]
        xNz = Tn * un[:, 2] + Ts * us[:, 2] + Td * ud[:, 2]
        xTrac = jnp.stack([xN, xNy, xNz], axis=1) * arn[:, None]
        initN, initS, initD = fric[:, 6], fric[:, 7], fric[:, 48]
        xInit = jnp.stack([initN * un[:, 0] + initS * us[:, 0] + initD * ud[:, 0],
                            initN * un[:, 1] + initS * us[:, 1] + initD * ud[:, 1],
                            initN * un[:, 2] + initS * us[:, 2] + initD * ud[:, 2]], axis=1) * arn[:, None]
        delta_f = xTrac - xInit * C_elastic

        flt_idx = jnp.concatenate([idxF_s[0], idxF_m[0], idxF_s[1], idxF_m[1], idxF_s[2], idxF_m[2]])
        flt_val = jnp.concatenate([delta_f[:, 0], -delta_f[:, 0], delta_f[:, 1], -delta_f[:, 1],
                                    delta_f[:, 2], -delta_f[:, 2]])
        force = force.at[flt_idx].add(flt_val)
        force = force.at[0].set(0.0)

        fric = fric.at[:, 77].set(Tn); fric = fric.at[:, 78].set(Ts); fric = fric.at[:, 79].set(Td)

        need = fnft > 5000.0
        fnft = jnp.where(need & (srMag >= inv['slipRateThres']), timeElapsed, fnft)

        force = force * inv['inv_mass_full']

        new_carry = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
        return new_carry, None

    return step


def run(S, nsteps=None, verbose=True):
    nsteps = nsteps or S['nstep']
    enable_compilation_cache()
    inv = build(S)
    # a few arrays didn't fit the flat dict-literal above; attach them here.
    conn = S['conn']; elemType = S['elemType']
    E_pml = np.nonzero(elemType == 2)[0]
    inv['conn_p'] = jnp.asarray(conn[E_pml])
    inv['massSlave'] = jnp.asarray(S['fnms'][S['nsmp1']])
    inv['massMaster'] = jnp.asarray(S['fnms'][S['nsmp2']])

    N = S['N']; NEQ = S['NEQ']; nftnd = S['nftnd']
    v1 = jnp.zeros(NEQ + 1, dtype=jnp.float64)
    velArr = jnp.zeros((N, 3), dtype=jnp.float64)
    dispArr = jnp.zeros((N, 3), dtype=jnp.float64)
    force = jnp.zeros(NEQ + 1, dtype=jnp.float64)
    stress_i = inv['stress_i0']  # zeros for C_elastic==1; lithostatic pre-stress otherwise.
    s_p = jnp.zeros((inv['Ep'], 15), dtype=jnp.float64)
    fric = jnp.asarray(S['fric_init'].copy())
    fnft = jnp.full(nftnd, 99999.0, dtype=jnp.float64)
    timeElapsed = jnp.asarray(0.0, dtype=jnp.float64)

    carry0 = (v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed)
    carry = time_loop(make_step, inv, carry0, nsteps)
    jax.block_until_ready(carry)
    v1, velArr, dispArr, force, stress_i, s_p, fric, fnft, timeElapsed = carry
    return dict(velArr=np.asarray(velArr), dispArr=np.asarray(dispArr), fnft=np.asarray(fnft),
                fric=np.asarray(fric), force=np.asarray(force))
