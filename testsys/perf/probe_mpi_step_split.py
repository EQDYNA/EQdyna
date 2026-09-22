#! /usr/bin/env python3
"""Measure ONE rank's local per-step cost of driver.run_mpi's two jits, with
NO MPI at all, at an arbitrary (rank, nranks) decomposition. Report-only.

WHY THIS ISOLATES THE QUESTION. run_mpi's per-step cost splits into jitted
compute and host/MPI orchestration. MPI4NodalQuant.decompose needs no
communication to build a rank's local view (every rank holds the full mesh),
so a SINGLE process can build rank r of N's exact local inv/finv/carry and run
the exact same a_jit/b_jit. If that measures ~the 32-rank per-rank ms/step,
the plateau is inside the computation and MPI is not implicated; if it
measures ~1/N of the 1-rank number, the plateau is in the orchestration layer.

It also times three sub-blocks separately, at the same shapes, so a
rank-independent cost can be attributed to a code REGION and not to "the step":
    head = velDispUpdate + zero the force           (driver.py:40-83, :137)
    elem = assembleGlobalKU + calcHourglassResist   (the SCALABLE work)
    div  = the mass divide                          (driver.py:177)
These sub-jits are a MEASUREMENT ONLY -- they materialise state between blocks
that the real step keeps fused, so their sum is an upper bound on the real
step, not an identity. Reported as such and never summed into a verdict.

Usage:
    EQDYNAROOT=... PROBE_OUT=x.json \
        probe_mpi_step_split.py CASE_DIR RANK NRANKS REPS
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.environ['EQDYNAROOT']
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

from eqdyna import backend as B                      # noqa: E402
xp = B.array_module('jax')                            # x64 on, before jnp
import jax                                            # noqa: E402
from eqdyna import assembleGlobalKU as KU             # noqa: E402
from eqdyna import driver, eqdyna3d, faulting as FLT  # noqa: E402
from eqdyna import MPI4NodalQuant as MQ               # noqa: E402
from eqdyna import globalvar as gv                    # noqa: E402


def nbytes(tree):
    return int(sum(np.asarray(x).nbytes for x in jax.tree_util.tree_leaves(tree)
                   if hasattr(x, 'dtype')))


def timeit(fn, reps):
    """Median of `reps` blocked calls, warmed once first. Median not mean, so
    one page-fault outlier cannot set the figure."""
    jax.block_until_ready(fn(0))
    ts = []
    for i in range(reps):
        t0 = time.perf_counter()
        jax.block_until_ready(fn(i + 1))
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts)) * 1e3


def main():
    case_dir = sys.argv[1]
    rank, nranks, reps = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    t0 = time.perf_counter()
    S, _mesh = eqdyna3d.build_solver_state(case_dir)
    inv = KU.build(S)
    finv = FLT.build(S)
    finv['tr'] = (FLT.forced_rupture_time(np, finv)
                  if FLT.nucleation_enabled(finv) else None)
    loc = MQ.decompose(S, inv, finv, rank, nranks)
    inv_l, finv_l = loc['inv'], loc['finv']
    computed = loc['fault_computed_rows']
    build_s = time.perf_counter() - t0

    mass_h = np.concatenate(([1.0], S['nodalMassArr']))
    B.check_index_width(inv_l)
    inv_l = B.to_device(xp, inv_l)
    finv_l = B.to_device(xp, finv_l)
    mass = xp.asarray(mass_h)

    N = S['N']; NEQ = S['NEQ']; nftnd_l = int(finv_l['nftnd'])
    z = (lambda *a: xp.zeros(*a))
    # EXACTLY driver.run_mpi:384-389 -- note v1/velArr/dispArr/force take
    # their shapes from the GLOBAL S['N'] and S['NEQ'] on every rank.
    carry = (z(NEQ + 1), z((N, 3)), z((N, 3)), z(NEQ + 1),
             xp.asarray(inv_l['stress_i0']).copy(), z((inv_l['Ep'], 15)),
             xp.asarray(S['fric_init'][computed].copy()),
             xp.full(nftnd_l, gv.FNFT_SENTINEL),
             xp.asarray(0.0), z((nftnd_l, 0)), z((nftnd_l, 0)))

    B.enable_compilation_cache()
    scratch = KU.alloc_scratch(xp, inv_l)
    dyn, sta = B.promote(xp, inv_l)
    halo = xp.asarray(loc['halo_idx'])
    dt = inv_l['dt']; rdampk = inv_l['rdampk']

    names = ('v1', 'velArr', 'dispArr', 'force', 'stress_i', 's_p', 'fric',
             'fnft', 'timeElapsed', 'sliprate_hist', 'shear_hist')
    # GLOBAL-SIZED: shape comes from S['N'] / S['NEQ'], so it is the SAME at
    # 1 rank and at 32. That is the property a rank-independent cost has.
    glob = {'v1', 'velArr', 'dispArr', 'force'}
    cb = {n: int(np.asarray(a).nbytes) for n, a in zip(names, carry)}
    rep = dict(case=os.path.basename(case_dir), rank=rank, nranks=nranks,
               reps=reps, build_s=round(build_s, 2), N=int(N), NEQ=int(NEQ),
               Ei=int(inv_l['Ei']), Ep=int(inv_l['Ep']), E=int(inv_l['E']),
               Ei_global=int(inv['Ei']), Ep_global=int(inv['Ep']),
               nodes_int=int(np.asarray(inv_l['int_nodes_idx']).shape[0]),
               nodes_pml=int(np.asarray(inv_l['pml_nodes_idx']).shape[0]),
               nodes_int_global=int(np.asarray(inv['int_nodes_idx']).shape[0]),
               nodes_pml_global=int(np.asarray(inv['pml_nodes_idx']).shape[0]),
               nftnd_local=nftnd_l, nftnd_global=int(finv['nftnd']),
               halo_eqs=int(np.asarray(halo).shape[0]), carry_bytes=cb,
               carry_bytes_total=sum(cb.values()),
               carry_bytes_global_sized=sum(v for k, v in cb.items()
                                            if k in glob),
               mass_bytes=int(mass_h.nbytes), dyn_bytes=nbytes(dyn),
               sta_bytes=nbytes(sta), cpus=sorted(os.sched_getaffinity(0)))

    def a_body(dyn_arrays, c, nt, h):
        part_a, _ = driver.make_step_parts(xp, {**sta, **dyn_arrays}, finv_l,
                                           None, mass, scratch)
        c = part_a(c, nt)
        return c, c[driver.FORCE][h]

    def b_body(dyn_arrays, c, nt, h, delta):
        _, part_b = driver.make_step_parts(xp, {**sta, **dyn_arrays}, finv_l,
                                           None, mass, scratch)
        force = B.addat(xp, c[driver.FORCE], h, delta)
        c = c[:driver.FORCE] + (force,) + c[driver.FORCE + 1:]
        return part_b(c, nt)

    a_jit = jax.jit(a_body, donate_argnums=(1,))
    b_jit = jax.jit(b_body, donate_argnums=(1,))
    zero_delta = xp.zeros(np.asarray(halo).shape[0])

    # IS THE DONATION EFFECTIVE? A donated buffer is DELETED by the call that
    # consumed it. If donation silently failed the old buffers stay alive and
    # every step copies the whole carry. Asked of the buffers themselves, not
    # of a warning that may or may not have been emitted.
    t0 = time.perf_counter()
    carry, hv = a_jit(dyn, carry, 1, halo)
    jax.block_until_ready(hv)
    rep['compile_a_s'] = round(time.perf_counter() - t0, 2)
    pr = [carry[i] for i in (0, 1, 2, 3)]
    t0 = time.perf_counter()
    carry = b_jit(dyn, carry, 1, halo, zero_delta)
    jax.block_until_ready(carry)
    rep['compile_b_s'] = round(time.perf_counter() - t0, 2)
    rep['donated_deleted_by_b'] = {names[i]: bool(pr[k].is_deleted())
                                   for k, i in enumerate((0, 1, 2, 3))}
    pr = [carry[i] for i in (0, 1, 2, 3)]
    carry, hv = a_jit(dyn, carry, 2, halo)
    jax.block_until_ready(hv)
    rep['donated_deleted_by_a'] = {names[i]: bool(pr[k].is_deleted())
                                   for k, i in enumerate((0, 1, 2, 3))}
    carry = b_jit(dyn, carry, 2, halo, zero_delta)
    jax.block_until_ready(carry)

    # XLA's own accounting for the two executables. `temp` is buffer traffic
    # the step pays that the carry sizes do not explain.
    for nm, f, args in (('a', a_jit, (dyn, carry, 3, halo)),
                        ('b', b_jit, (dyn, carry, 3, halo, zero_delta))):
        an = f.lower(*args).compile().memory_analysis()
        rep['mem_' + nm] = dict(temp=int(an.temp_size_in_bytes),
                                arg=int(an.argument_size_in_bytes),
                                out=int(an.output_size_in_bytes),
                                alias=int(an.alias_size_in_bytes))

    st = {}
    cur = [carry]

    def run_a(i):
        c, h = a_jit(dyn, cur[0], i + 10, halo)
        cur[0] = c
        return h

    def run_b(i):
        cur[0] = b_jit(dyn, cur[0], i + 10, halo, zero_delta)
        return cur[0]

    def run_ab(i):
        c, h = a_jit(dyn, cur[0], i + 10, halo)
        cur[0] = b_jit(dyn, c, i + 10, halo, zero_delta)
        return cur[0]

    st['a_ms'] = timeit(run_a, reps)
    st['b_ms'] = timeit(run_b, reps)
    st['ab_ms'] = timeit(run_ab, reps)

    inv_f = {**sta, **dyn}
    v1, velArr, dispArr, force = cur[0][0], cur[0][1], cur[0][2], cur[0][3]

    @jax.jit
    def head(v1, velArr, dispArr, force):
        v1, velArr, dispArr = driver.velDispUpdate(xp, inv_f, v1, velArr,
                                                   dispArr, force, dt)
        return v1, velArr, dispArr, B.setat(xp, force, slice(None), 0.0)

    @jax.jit
    def elem(velArr, dispArr, force, stress_i, s_p):
        force, stress_i, s_p = KU.assembleGlobalKU(
            xp, inv_f, velArr, force, stress_i, s_p, dt, rdampk, scratch)
        return KU.calcHourglassResist(xp, inv_f, dispArr, velArr, force,
                                      rdampk), stress_i, s_p

    @jax.jit
    def div(force):
        force = B.setat(xp, force, slice(1, None), force[1:] / mass[1:])
        return B.setat(xp, force, 0, 0.0)

    st['head_ms'] = timeit(lambda i: head(v1, velArr, dispArr, force), reps)
    st['elem_ms'] = timeit(
        lambda i: elem(velArr, dispArr, force, cur[0][4], cur[0][5]), reps)
    st['div_ms'] = timeit(lambda i: div(force), reps)
    rep['timings_ms'] = {k: round(v, 3) for k, v in st.items()}

    out = os.environ['PROBE_OUT']
    json.dump(rep, open(out, 'w'), indent=1)
    print(json.dumps(rep, indent=1), flush=True)
    print('WROTE %s' % out, flush=True)


if __name__ == '__main__':
    main()
