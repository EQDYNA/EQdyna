#! /usr/bin/env python3
"""Strong-scaling measurement for the REAL-MPI jax path (one process per
rank, driver.run_mpi + MPI4NodalQuant.py) against Fortran MPI on the SAME
cpus in the SAME session. Report-only, never a gate.

WHY A SAME-SESSION FORTRAN COLUMN IS NOT OPTIONAL. The 14.39x Fortran figure
this campaign is measured against is from 2026-09-19 at v5.12.0, on that
day's tenancy. This box carries 10-25 foreign single-core jobs plus a
~15-core photogrammetry job, and the memory bandwidth they consume is shared
and not pinnable. A jax speedup divided by a stale Fortran denominator is not
a comparison, so this tool measures both, back to back, on the same cpu set.

CPU SELECTION IS BY MEASURED LOAD, NOT BY TOPOLOGY ORDER. Two things were
measured on this box and both are designed around here:
  * `mpirun --bind-to core --map-by core` always takes socket-0 cores 0..N-1,
    and some of those carry a tenant at 1.00 utilisation. A rank landing
    there ran 2.01x slower on identical work -- and in a halo-synchronised
    solver the SLOWEST rank sets the step, so blind binding converts a
    fluctuating average into a pinned straggler.
  * Per-NODE selection is not enough either: cpu 44 and cpu 46 are both NUMA
    node 5 and measured 7.94 vs 21.22 ms/step on identical work.
So placement is `mpirun --cpu-set <the k least-loaded cpus, measured now>
--bind-to cpu-list:ordered`, and the busy fraction of every chosen cpu is
recorded beside the number it produced.

THE STRICT CEILING CANNOT BE MET AT THIS TENANCY and pretending otherwise
just records SKIPs: two independent probes gave 0 of 64 cpus under 0.20,
while 0.45 gave 32 then 28. So this tool does not gate on a ceiling at all --
it takes the k least-loaded cpus it can find, REPORTS each one's busy
fraction, and refuses only if the chosen set's worst cpu exceeds
--max-busy (default 0.5), which is a statement about what was measured
rather than a filter that silently selects nothing.

PER-STEP BY DIFFERENCE, both engines, two step counts, fresh processes -- the
same technique and the same reason as run_scaling.py (`88f5227`): at 20-60
steps the fixed cost (mesh build, XLA compile, MPI init) is a large fraction
of the wall time, and subtracting it from one engine only makes the verdict a
property of the metric.

WHAT IT PRINTS FOR EVERY POINT, because a multi-rank number without these is
not interpretable on this box: rank count (counted from the ranks' own
stdout, not assumed from -np), per-rank element count, per-rank ms/step with
max/mean, EFFECTIVE_CORES per rank, the busy fraction of every cpu used, and
threads per rank.

Usage:
    python3 testsys/perf/run_mpi_scaling.py [--case test.tpv104]
        [--ranks 1,2,4,8,16] [--n-lo 20] [--n-hi 60] [--max-busy 0.5]
        [--skip-fortran] [--repeats 1]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('EQDYNAROOT') or os.path.dirname(os.path.dirname(TESTSYS))
PYTHON_PKG = os.path.join(ROOT, 'src', 'python')
OUT = os.path.join(ROOT, 'docs', 'perf_snapshots',
                   'mpi_scaling_%s_%s.json' % (time.strftime('%Y-%m-%d'),
                                          os.environ.get('EQDYNA_SNAPSHOT_TAG', 'tpv104')))
# PROJECT_RULES rule 19: dated and immutable, never a shared *_last.* file.
# A path whose name says "last" cites whatever ran most recently, so it is
# not evidence a board row can point at; two tools writing one such file
# silently discard each other.

sys.path.insert(0, TESTSYS)
import run_numa_scaling as numa      # noqa: E402
import run_scaling as rs             # noqa: E402

RANK_RE = re.compile(r'rank (\d+)/(\d+) wrote (\S+)\s+(.*)')


def least_loaded_cpus(nodes, k, exclude=()):
    """The k least-loaded cpus right now, with their busy fractions.

    Ordered by (busy fraction, node, cpu) so that among equally idle cpus the
    choice is still as NUMA-compact as the free set allows, and the tenant's
    hot cpus are avoided rather than taken because they sort first.

    `exclude` drops named cpus from the candidate set. MEASURED REASON, not a
    preference: /proc/stat read cpu 0 and cpu 1 at busy 0.00, and a 2-rank
    tpv104 run placed there took 1859 ms/step with rank 0 reporting
    EFFECTIVE_CORES 0.39, against 605 ms/step on cpus 8,9 and 591 on 8,16
    (20 steps each, same case, same binary, same session). A cpu that reads
    idle and then delivers 0.39 cores is not idle; it is stalled on something
    /proc/stat does not see, and in a halo-synchronised solver the slowest
    rank sets the step. Excluded cpus are recorded in the snapshot."""
    all_cpus = sorted(c for cs in nodes.values() for c in cs
                      if c not in set(exclude))
    busy = numa.cpu_busy_fractions(all_cpus)
    if not busy:
        raise SystemExit('FAIL: could not read per-cpu utilisation from '
                         '/proc/stat -- cannot choose a placement blind.')
    cpu2node = {c: n for n, cs in nodes.items() for c in cs}
    ranked = sorted((b, cpu2node[c], c) for c, b in busy.items() if b is not None)
    if len(ranked) < k:
        raise SystemExit('FAIL: read utilisation for %d cpus, need %d'
                         % (len(ranked), k))
    chosen = sorted(c for _b, _n, c in ranked[:k])
    return chosen, {c: busy[c] for c in chosen}


def jax_mpi_once(case_dir, nsteps, cpus, ranks, sync):
    """One mpirun of the jax MPI path. Returns (wall, [per-rank dicts])."""
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    env['EQDYNA_MPI_SYNC'] = sync
    env['PYTHONPATH'] = PYTHON_PKG + os.pathsep + env.get('PYTHONPATH', '')
    for key in ('XLA_FLAGS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                'EQDYNA_JAX_DEVICES'):
        env.pop(key, None)
    cmd = ['mpirun', '--cpu-set', ','.join(str(c) for c in cpus),
           '--bind-to', 'cpu-list:ordered', '-np', str(ranks),
           sys.executable, '-m', 'eqdyna', case_dir, str(nsteps),
           '--backend', 'jax', '--mpi']
    t0 = time.time()
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-2500:]); print(r.stderr[-2500:])
        return None, None
    per = []
    for line in r.stdout.splitlines():
        m = RANK_RE.match(line.strip())
        if m:
            d = dict(rank=int(m.group(1)), nranks=int(m.group(2)))
            for kv in m.group(4).split():
                k, _, v = kv.partition('=')
                d[k] = float(v) if re.match(r'^-?[\d.]+$', v) else v
            per.append(d)
    # EFFECTIVE_CORES comes from the solver's own line, matched separately so
    # a change to either print cannot silently drop it.
    eff = dict((int(a), float(b)) for a, b in
               re.findall(r'rank (\d+): .*EFFECTIVE_CORES ([\d.]+)', r.stdout))
    for d in per:
        if d['rank'] not in eff:
            raise RuntimeError('rank %d reported no EFFECTIVE_CORES -- a '
                               'per-rank number without it is not '
                               'interpretable on this box' % d['rank'])
        d['effective_cores'] = eff[d['rank']]
    # THE VACUOUS-GATE GUARD: a run that launched fewer ranks than asked, or
    # whose ranks did not all report, still produces a plausible wall time.
    if len(per) != ranks:
        raise RuntimeError(
            'asked for %d ranks, %d reported (%r). A rank count taken from '
            '-np rather than from the ranks themselves is how a multi-process '
            'measurement goes vacuous invisibly.'
            % (ranks, len(per), sorted(d['rank'] for d in per)))
    return wall, per


def per_step_jax_mpi(case_dir, cpus, ranks, n_lo, n_hi, sync):
    lo = jax_mpi_once(case_dir, n_lo, cpus, ranks, sync)
    hi = jax_mpi_once(case_dir, n_hi, cpus, ranks, sync)
    if lo[0] is None or hi[0] is None:
        return None
    ps = (hi[0] - lo[0]) / float(n_hi - n_lo)
    ms = [d['ms_per_step'] for d in hi[1]]
    return dict(ms_per_step=ps * 1e3, fixed_s=lo[0] - n_lo * ps,
                wall_lo_s=lo[0], wall_hi_s=hi[0],
                rank_ms=ms, rank_ms_max=max(ms), rank_ms_mean=sum(ms) / len(ms),
                straggler=max(ms) / (sum(ms) / len(ms)),
                sync=sync,
                mpi_ms=[d['mpi_ms_per_step'] for d in hi[1]],
                wait_ms=[d['wait_ms_per_step'] for d in hi[1]],
                eff=[d['effective_cores'] for d in hi[1]],
                Ei=[int(d['Ei']) for d in hi[1]],
                Ep=[int(d['Ep']) for d in hi[1]],
                halo=[int(d['halo_eqs']) for d in hi[1]],
                threads=[int(d['threads']) for d in hi[1]],
                cpus_allowed=[int(d['cpus_allowed']) for d in hi[1]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='test.tpv104')
    ap.add_argument('--ranks', default='1,2,4,8,16')
    ap.add_argument('--n-lo', type=int, default=20)
    ap.add_argument('--n-hi', type=int, default=60)
    ap.add_argument('--max-busy', type=float, default=0.5)
    ap.add_argument('--repeats', type=int, default=1)
    ap.add_argument('--syncs', default='halo,allreduce',
                    help='nodal-sync modes to measure, in order. `halo` moves '
                         'O(boundary) (MPI4NodalQuant\'s own pattern); '
                         '`allreduce` moves O(NEQ) with no ownership '
                         'bookkeeping. Both are correct; the difference is '
                         'the whole question.')
    ap.add_argument('--skip-fortran', action='store_true')
    ap.add_argument('--exclude-cpus', default='',
                    help='comma-separated cpus never to place a rank on. See '
                         'least_loaded_cpus: cpu 0 and cpu 1 read busy 0.00 '
                         'and then delivered 0.39 effective cores (1859 vs '
                         '605 ms/step on identical work), so excluding them '
                         'is a measurement, not a preference.')
    a = ap.parse_args()
    ranks = [int(x) for x in a.ranks.split(',') if x]
    syncs = [s for s in a.syncs.split(',') if s]
    for s in syncs:
        if s not in ('halo', 'allreduce'):
            raise SystemExit('FAIL: unknown sync mode %r' % s)
    for n in ranks:
        if n > 16:
            raise SystemExit('FAIL: %d ranks requested -- capped at 16.' % n)

    excl = [int(x) for x in a.exclude_cpus.split(',') if x.strip()]
    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology.')
    rs.CASE = a.case
    case_dir = rs.build_py_case(a.case)
    sha = rs.sh('git -C %s rev-parse --short HEAD' % ROOT).stdout.strip()
    work = tempfile.mkdtemp(prefix='mpiscaling.', dir=os.environ.get('TMPDIR', '/tmp'))
    dt = None
    if not a.skip_fortran:
        d0 = os.path.join(work, 'dtprobe')
        rs.make_case(d0, 1, 1, 1)
        _, dt = rs.read_term_dt(d0)
        shutil.rmtree(d0, ignore_errors=True)

    print('case %s  sha %s  host %s  n_lo/n_hi %d/%d  max-busy %.2f'
          % (a.case, sha, os.uname().nodename, a.n_lo, a.n_hi, a.max_busy),
          flush=True)
    rows = []
    for n in ranks:
        cpus, busy = least_loaded_cpus(nodes, n, exclude=excl)
        worst = max(busy.values())
        print('\n-- %d rank(s) -- cpus %s  busy %s  worst %.2f (max-busy %.2f), '
              'whole-box load %.1f'
              % (n, cpus, {c: round(b, 2) for c, b in busy.items()}, worst,
                 a.max_busy, os.getloadavg()[0]), flush=True)
        if worst > a.max_busy:
            print('  SKIPPED: the %d least-loaded cpus include one at %.2f > '
                  '%.2f. Not measured, not silently included.'
                  % (n, worst, a.max_busy), flush=True)
            rows.append(dict(ranks=n, skipped=True, cpus=cpus,
                             busy={str(c): b for c, b in busy.items()}))
            continue
        row = dict(ranks=n, cpus=cpus, busy={str(c): b for c, b in busy.items()},
                   loadavg=os.getloadavg(), n_lo=a.n_lo, n_hi=a.n_hi)
        for sync in syncs:
            best = None
            for _ in range(a.repeats):
                r = per_step_jax_mpi(case_dir, cpus, n, a.n_lo, a.n_hi, sync)
                if r is not None and (best is None
                                      or r['ms_per_step'] < best['ms_per_step']):
                    best = r
            if best is None:
                print('  jax-mpi/%-9s FAILED' % sync, flush=True)
                continue
            row['jax_' + sync] = best
            if sync == syncs[0]:
                row['jax'] = best
            print('  jax-mpi/%-9s %9.2f ms/step (by difference)  fixed %6.2fs  '
                  'wall %.1f/%.1fs' % (sync, best['ms_per_step'], best['fixed_s'],
                                       best['wall_lo_s'], best['wall_hi_s']),
                  flush=True)
            print('             per-rank ms/step %s' % [round(x, 2) for x in best['rank_ms']],
                  flush=True)
            print('             max/mean %.2fx  exchange %s ms/step  barrier '
                  'wait %s ms/step'
                  % (best['straggler'], [round(x, 3) for x in best['mpi_ms']],
                     [round(x, 2) for x in best['wait_ms']]), flush=True)
            print('             EFFECTIVE_CORES %s  threads %s  cpus_allowed %s'
                  % ([round(x, 2) for x in best['eff']], best['threads'],
                     best['cpus_allowed']), flush=True)
            print('             per-rank Ei %s  Ep %s  halo eqs %s'
                  % (best['Ei'], best['Ep'], best['halo']), flush=True)
        if not a.skip_fortran and n in rs.DECOMP:
            ps, fixed, lo, hi, _o1, _o2 = rs.per_step_fortran(
                work, n, 'leastloaded', cpus, nodes, dt, a.n_lo, a.n_hi)
            row['fortran'] = dict(ms_per_step=ps * 1e3, fixed_s=fixed,
                                  wall_lo_s=lo, wall_hi_s=hi,
                                  decomp=rs.DECOMP[n])
            print('  fortran    %9.2f ms/step (by difference)  fixed %6.2fs  '
                  'wall %.1f/%.1fs  decomp %s'
                  % (ps * 1e3, fixed, lo, hi, rs.DECOMP[n]), flush=True)
        rows.append(row)

    def first(key):
        return next((r[key]['ms_per_step'] for r in rows if r.get(key)), None)

    cols = ['fortran'] + ['jax_' + s for s in syncs]
    base = {c: first(c) for c in cols}
    hdr = '%-7s' % 'ranks'
    for c in cols:
        hdr += ' %13s %9s' % (c, 'speedup')
    print('\n' + hdr + ' %9s' % 'jax/fortran')
    for r in rows:
        line = '%-7d' % r['ranks']
        for c in cols:
            v = r.get(c, {}).get('ms_per_step')
            line += ' %13s %9s' % ('%.2f' % v if v else '-',
                                   '%.2fx' % (base[c] / v) if v and base[c] else '-')
        f = r.get('fortran', {}).get('ms_per_step')
        j = r.get('jax_' + syncs[0], {}).get('ms_per_step')
        line += ' %9s' % ('%.2fx' % (j / f) if f and j else '-')
        print(line)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(dict(case=a.case, sha=sha, host=os.uname().nodename,
                   date=time.strftime('%Y-%m-%d %H:%M'), n_lo=a.n_lo,
                   n_hi=a.n_hi, max_busy=a.max_busy, excluded_cpus=excl,
                   rows=rows),
              open(OUT, 'w'), indent=1)
    print('\nsaved %s' % OUT)
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
