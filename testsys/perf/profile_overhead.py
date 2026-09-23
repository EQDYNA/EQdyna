#! /usr/bin/env python3
"""Overhead gate for the per-rank profiler (EQDYNA_PROFILE=1): does turning
profiling ON cost measurably more than turning it off?

DO NOT RUN THIS FROM AN AGENT SESSION (owner instruction, 2026-09-23): it
launches the real solver, twice per step-count, pinned with numactl. It is
written and reviewed here, never executed here.

METHOD, mirroring this repo's one standard technique
(run_scaling.per_step_py / run_numa_scaling.per_step_and_fixed):
  1. Run the SAME case/backend at two step counts (n_lo, n_hi) with
     EQDYNA_PROFILE=1, and again with EQDYNA_PROFILE=0. Per-step cost is the
     DIFFERENCE (t(n_hi) - t(n_lo)) / (n_hi - n_lo) in each case -- this
     cancels the fixed cost (interpreter start, case load, JIT compile),
     exactly as every other per-step measurement in this repo does, and for
     the same reason: a raw wall-clock ratio would mix profiling overhead
     with compile/setup and drift when THOSE change, not when the profiler
     does.
  2. `overhead_ms = per_step(profile=1) - per_step(profile=0)`.
  3. A NOISE FLOOR is measured directly, not assumed: `--repeats` independent
     per_step(profile=0) measurements, same pinned placement; the floor is
     their max-min spread. This box's own run-to-run variance has been
     measured elsewhere in this repo to be large enough to manufacture a
     33% "regression" out of nothing (CLAUDE.md, "Measure, do not infer") --
     an overhead gate that skips measuring its OWN noise floor would repeat
     exactly that mistake, aimed at itself.
  4. FAIL if `overhead_ms > max(REL_OVERHEAD_BUDGET * per_step(profile=0),
     noise_floor_ms)`. Both terms are named, not folded into one constant,
     because they answer different questions: the relative budget is "is the
     profiler cheap enough to leave on by default", the noise floor is "is
     this even a measurement or box jitter".

PLACEMENT: `--cpus` is REQUIRED (no default, no auto-selection) -- an
overhead measurement on a placement nobody can name is not reproducible.
`--busy-ceiling` defaults to 0.1 (stricter than run_scaling's 0.45): this is
a two-run DIFFERENCE-of-a-difference, so it is more sensitive to transient
tenancy than a single per-step figure, and a busy box is not a probe of the
profiler's cost.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))          # testsys/perf/
TESTSYS = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(TESTSYS)
sys.path.insert(0, HERE)
sys.path.insert(0, TESTSYS)
import run_numa_scaling as numa   # noqa: E402
import run_scaling as rs          # noqa: E402  (build_py_case, numactl_prefix, PYTHON_PKG)

# Named, not folded together -- see module docstring point 4.
REL_OVERHEAD_BUDGET = 0.01   # 1% of the profile=0 per-step cost.

MPI_BACKENDS = ('fortran', 'python-jax-mpi')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')

# EXTENSION (2026-09-23, this landing): the MPI arms (fortran/python-jax-mpi
# at 4 ranks). CLAUDE.md's own convention -- "for MPI, difference each
# rank's own solve time, not the wrapper's wall clock" -- is honoured
# literally here: both arms read an EQDYNA_PROFILE-INDEPENDENT per-rank
# timer (fortran: compTime<rank>'s compTimeInSeconds(9), gated on
# writeCompTime which is hardcoded 1, never on EQDYNA_PROFILE; jax-mpi:
# driver.run_mpi's own unconditional 'rank R: T s,' stdout line -- see
# docs/run_profile.md and driver.py:670), differenced per rank between
# n_lo/n_hi, never the outer mpirun wall clock. This is a DELIBERATE
# departure from testsys/perf/run_scaling.py's own Fortran scaling gate
# (per_step_fortran), which differences the outer wall clock and is not
# wrong for ITS question (rank-count scaling) -- it is simply a different
# question from "did turning profiling on cost this run anything", which
# is what this module answers.
_JAXMPI_RANK_RE = re.compile(r'driver\.run_mpi rank (\d+): ([0-9.eE+-]+) s,')


def _mpi_pinned_cmd(ranks, cpus, node_map, argv):
    """One `mpirun --bind-to none -np <ranks> bash -c '...'` string, each
    rank pinned to its OWN cpu/numa node via $OMPI_COMM_WORLD_LOCAL_RANK --
    the identical technique and rationale as run_scaling.fortran_cmd's own
    docstring (OpenMPI's default --bind-to core binds each rank to its own
    choice of core BEFORE an outer numactl can move it, so it must be
    disabled and replaced by this explicit per-rank bind). `argv` is the
    program and its arguments as a plain list (never shell-embedded
    unescaped)."""
    import shlex
    cpu_list = ' '.join(str(c) for c in cpus[:ranks])
    cpu2node = rs.cpu_to_node(node_map)
    node_list = ' '.join(str(cpu2node[c]) for c in cpus[:ranks])
    prog = ' '.join(shlex.quote(a) for a in argv)
    inner = ('CPUS=(%s); NODES=(%s); exec numactl '
             '--physcpubind=${CPUS[$OMPI_COMM_WORLD_LOCAL_RANK]} '
             '--membind=${NODES[$OMPI_COMM_WORLD_LOCAL_RANK]} %s'
             % (cpu_list, node_list, prog))
    return "%s --bind-to none -np %d bash -c %s" % (MPIRUN, ranks,
                                                     shlex.quote(inner))


def build_fortran_term_case(case_name, dst, term, ranks):
    """create.newcase + case.setup for `case_name` at its NATURAL ranks-way
    decomposition (rs.DECOMP[ranks], never forced-serial: the fortran
    binary is the thing under test) and an explicit `par.term` override,
    via run_scaling.make_case -- reused, not reimplemented (rule 1).
    run_scaling.make_case reads its case name from the MODULE-LEVEL
    `rs.CASE` constant (its own scaling sweep hardcodes test.tpv104), so it
    is monkeypatched for the span of this one call and restored
    immediately after, success or failure."""
    nx, ny, nz = rs.DECOMP[ranks]
    old_case = rs.CASE
    rs.CASE = case_name
    try:
        rs.make_case(dst, nx, ny, nz, term=term)
    finally:
        rs.CASE = old_case
    return rs.read_term_dt(dst)   # (term_written, dt), both as case.setup wrote them


def fortran_rank_totals_s(case_dir, cpus, node_map, ranks, profile_flag):
    """One real mpirun of the shipped fortran binary in `case_dir` (whatever
    nstep its own bGlobal.txt encodes), EQDYNA_PROFILE forced to
    `profile_flag`. Returns {rank: compTimeInSeconds(9)} read back from
    each rank's own compTime<rank> file -- written by output_timeanalysis,
    gated on writeCompTime (hardcoded 1), independent of EQDYNA_PROFILE, so
    this per-rank timer exists on BOTH arms of the A/B."""
    binary = os.path.join(rs.ROOT, 'bin', 'eqdyna')
    if not os.path.exists(binary):
        raise RuntimeError('fortran overhead probe: %s missing -- build it '
                           'first (./install-eqdyna.sh -m ubuntu)' % binary)
    for name in os.listdir(case_dir):
        if name.startswith('compTime'):
            os.remove(os.path.join(case_dir, name))
    env = dict(os.environ)
    env['EQDYNA_PROFILE'] = profile_flag
    cmd = _mpi_pinned_cmd(ranks, cpus, node_map, [binary])
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True,
                       cwd=case_dir, env=env)
    if r.returncode != 0:
        raise RuntimeError(
            'fortran overhead probe (EQDYNA_PROFILE=%s) exited %d:\n%s\n%s'
            % (profile_flag, r.returncode, r.stdout[-1500:], r.stderr[-1500:]))
    totals = {}
    for rank in range(ranks):
        path = os.path.join(case_dir, 'compTime%d' % rank)
        if not os.path.exists(path):
            continue   # a rank owning 0 fault nodes still writes compTime;
                       # absence here means the launch started fewer real
                       # workers than `ranks`, which the caller must not hide.
        with open(path) as f:
            fields = f.read().split()
        totals[rank] = float(fields[8])   # compTimeInSeconds(9), 0-indexed 8
    if len(totals) != ranks:
        raise RuntimeError(
            'fortran overhead probe (EQDYNA_PROFILE=%s): found %d of %d '
            'compTime<rank> file(s) in %s -- a launch that started fewer '
            'real workers than %d ranks must not be able to produce a '
            'measurement' % (profile_flag, len(totals), ranks, case_dir, ranks))
    return totals


def jaxmpi_rank_totals_s(case_dir, cpus, node_map, ranks, nsteps, profile_flag):
    """One real mpirun of `python3 -u -m eqdyna <case_dir> <nsteps> --backend
    jax --mpi`, EQDYNA_PROFILE forced to `profile_flag`. Returns {rank:
    total_s} parsed from driver.run_mpi's own unconditional per-rank stdout
    line (`_JAXMPI_RANK_RE`) -- unconditional on EQDYNA_PROFILE exactly like
    the fortran arm's compTime file, so both arms of the A/B read a real
    per-rank number regardless of which side of the switch is being timed."""
    env = dict(os.environ)
    env['PYTHONPATH'] = rs.PYTHON_PKG
    env['PYTHONUNBUFFERED'] = '1'
    env['JAX_PLATFORMS'] = 'cpu'
    env['EQDYNA_PROFILE'] = profile_flag
    env['EQDYNA_MPI_SYNC'] = 'halo'
    for k in ('XLA_FLAGS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        env.pop(k, None)
    argv = [sys.executable, '-u', '-m', 'eqdyna', case_dir, str(nsteps),
           '--backend', 'jax', '--mpi']
    cmd = _mpi_pinned_cmd(ranks, cpus, node_map, argv)
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True,
                       cwd=REPO_ROOT, env=env)
    if r.returncode != 0:
        raise RuntimeError(
            'jax-mpi overhead probe (EQDYNA_PROFILE=%s, nsteps=%d) exited '
            '%d:\n%s\n%s' % (profile_flag, nsteps, r.returncode,
                             r.stdout[-1500:], r.stderr[-1500:]))
    totals = {int(m.group(1)): float(m.group(2))
             for m in _JAXMPI_RANK_RE.finditer(r.stdout)}
    if len(totals) != ranks:
        raise RuntimeError(
            'jax-mpi overhead probe (EQDYNA_PROFILE=%s, nsteps=%d): parsed '
            '%d of %d per-rank stdout line(s) -- a launch that started '
            'fewer real workers than %d ranks must not be able to produce '
            'a measurement. stdout tail:\n%s'
            % (profile_flag, nsteps, len(totals), ranks, ranks,
               r.stdout[-1500:]))
    return totals


def mpi_rank_totals(backend, case_dir_or_builder, cpus, node_map, ranks,
                    nsteps, profile_flag):
    if backend == 'fortran':
        return fortran_rank_totals_s(case_dir_or_builder, cpus, node_map,
                                     ranks, profile_flag)
    return jaxmpi_rank_totals_s(case_dir_or_builder, cpus, node_map, ranks,
                               nsteps, profile_flag)


def per_step_ms_mpi(backend, case_dir_lo, case_dir_hi, cpus, node_map, ranks,
                    n_lo, n_hi, profile_flag):
    """{rank: ms/step} by difference, per rank, between two REAL runs at
    n_lo/n_hi. For fortran, `case_dir_lo`/`case_dir_hi` are two case
    directories built at two different `par.term` values (nstep has no
    runtime override in the fortran binary); for jax-mpi they are the SAME
    directory (nsteps is a CLI argument, no rebuild needed) -- the caller
    passes the same path twice in that case."""
    lo = mpi_rank_totals(backend, case_dir_lo, cpus, node_map, ranks, n_lo,
                        profile_flag)
    hi = mpi_rank_totals(backend, case_dir_hi, cpus, node_map, ranks, n_hi,
                        profile_flag)
    dn = n_hi - n_lo
    out = {}
    for rank in sorted(set(lo) & set(hi)):
        ms = (hi[rank] - lo[rank]) * 1000.0 / dn
        if ms <= 0:
            raise SystemExit(
                'FAIL: rank %d non-positive per-step time (%.6g ms/step) for '
                '%s EQDYNA_PROFILE=%s (n_lo=%d t=%.4fs, n_hi=%d t=%.4fs) -- '
                'per CLAUDE.md convention this is a bug to raise on, not a '
                'number to report' % (rank, ms, backend, profile_flag, n_lo,
                                      lo[rank], n_hi, hi[rank]))
        out[rank] = ms
    return out


def timed_run(case_dir, nsteps, cpus, node_map, backend, profile_flag):
    """One numactl-pinned, fresh-process solve of `nsteps` steps, with
    EQDYNA_PROFILE forced to `profile_flag` ('0' or '1'). Mirrors
    run_scaling.time_one_py's subprocess/env discipline exactly (same reason:
    a fresh process per timing so neither jax's compile cache nor a stale
    thread pool carries over between the two arms being compared)."""
    # BUG FOUND RUNNING THIS FOR REAL (2026-09-23, this landing): `backend`
    # here is the matrix.py-style label ('python-numpy'/'python-jax'), but
    # eqdyna3d.run_case's own `backend` parameter takes its INTERNAL name
    # ('numpy'/'jax') -- see run_e2e.run_standalone's identical mapping.
    # Passing the label straight through raised ValueError from
    # _resolve_solver on the very first real invocation of this
    # never-before-run module. Fixed here (testsys/perf tooling, not src/).
    engine = {'python-numpy': 'numpy', 'python-jax': 'jax'}[backend]
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from eqdyna import eqdyna3d\n"
        "import time; t0 = time.time()\n"
        "eqdyna3d.run_case(%r, nsteps=%d, verbose=False, backend=%r)\n"
        "print('WALL', time.time() - t0)\n"
        % (rs.PYTHON_PKG, case_dir, nsteps, engine))
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    env['EQDYNA_PROFILE'] = profile_flag
    for k in ('XLA_FLAGS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        env.pop(k, None)
    cmd = rs.numactl_prefix(cpus, node_map) + [sys.executable, '-c', script]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            'overhead probe (EQDYNA_PROFILE=%s, nsteps=%d) exited %d:\n%s\n%s'
            % (profile_flag, nsteps, r.returncode, r.stdout[-1500:], r.stderr[-1500:]))
    for line in r.stdout.splitlines():
        if line.startswith('WALL'):
            return float(line.split()[1])
    raise RuntimeError('overhead probe (EQDYNA_PROFILE=%s, nsteps=%d) printed '
                       'no WALL line -- stdout:\n%s' % (profile_flag, nsteps,
                                                        r.stdout[-1500:]))


def per_step_ms(case_dir, cpus, node_map, backend, profile_flag, n_lo, n_hi):
    t_lo = timed_run(case_dir, n_lo, cpus, node_map, backend, profile_flag)
    t_hi = timed_run(case_dir, n_hi, cpus, node_map, backend, profile_flag)
    ps, _fixed = numa.per_step_and_fixed(t_lo, t_hi, n_lo, n_hi)
    if ps <= 0:
        raise SystemExit(
            'FAIL: non-positive per-step time (%.6g s) for EQDYNA_PROFILE=%s '
            '(n_lo=%d t=%.3fs, n_hi=%d t=%.3fs) -- per CLAUDE.md convention '
            'this is a bug to raise on, not a number to report'
            % (ps, profile_flag, n_lo, t_lo, n_hi, t_hi))
    return ps * 1000.0   # ms/step


def noise_floor_ms(case_dir, cpus, node_map, backend, n_lo, n_hi, repeats):
    """`repeats` independent profile=0 per-step measurements; the floor is
    their spread (max - min), in ms/step -- see module docstring point 3."""
    samples = [per_step_ms(case_dir, cpus, node_map, backend, '0', n_lo, n_hi)
               for _ in range(repeats)]
    return max(samples) - min(samples), samples


def _mean(xs):
    return sum(xs) / float(len(xs))


def noise_floor_ms_mpi(backend, case_dir_lo, case_dir_hi, cpus, node_map,
                       ranks, n_lo, n_hi, repeats):
    """`repeats` independent profile=0 per-rank-mean measurements; the floor
    is their spread (max - min), same construction as the serial
    `noise_floor_ms` and for the same reason."""
    samples = []
    for _ in range(repeats):
        per_rank = per_step_ms_mpi(backend, case_dir_lo, case_dir_hi, cpus,
                                   node_map, ranks, n_lo, n_hi, '0')
        samples.append(_mean(per_rank.values()))
    return max(samples) - min(samples), samples


def main_mpi(a, cpus, node_map):
    if a.backend == 'fortran':
        work = os.path.join(TESTSYS, 'overhead_case')
        os.makedirs(work, exist_ok=True)
        dst_probe = os.path.join(work, 'probe')
        _term, dt = build_fortran_term_case(a.case, dst_probe, 1.0, a.ranks)
        shutil.rmtree(dst_probe, ignore_errors=True)
        dst_lo = os.path.join(work, 'lo')
        dst_hi = os.path.join(work, 'hi')
        build_fortran_term_case(a.case, dst_lo, dt * a.n_lo, a.ranks)
        build_fortran_term_case(a.case, dst_hi, dt * a.n_hi, a.ranks)
        case_dir_lo, case_dir_hi = dst_lo, dst_hi
        cleanup = lambda: shutil.rmtree(work, ignore_errors=True)  # noqa: E731
    else:
        case_dir = rs.build_py_case(a.case)
        case_dir_lo = case_dir_hi = case_dir   # nsteps is a CLI arg; no rebuild
        cleanup = lambda: None  # noqa: E731

    try:
        off_per_rank = per_step_ms_mpi(a.backend, case_dir_lo, case_dir_hi,
                                       cpus, node_map, a.ranks, a.n_lo,
                                       a.n_hi, '0')
        on_per_rank = per_step_ms_mpi(a.backend, case_dir_lo, case_dir_hi,
                                      cpus, node_map, a.ranks, a.n_lo,
                                      a.n_hi, '1')
        floor_ms, floor_samples = noise_floor_ms_mpi(
            a.backend, case_dir_lo, case_dir_hi, cpus, node_map, a.ranks,
            a.n_lo, a.n_hi, a.repeats)
    finally:
        cleanup()

    off_mean = _mean(off_per_rank.values())
    on_mean = _mean(on_per_rank.values())
    overhead_ms = on_mean - off_mean
    budget_ms = REL_OVERHEAD_BUDGET * off_mean
    limit_ms = max(budget_ms, floor_ms)

    print('per-rank EQDYNA_PROFILE=0 (ms/step): %s'
         % {r: round(v, 4) for r, v in sorted(off_per_rank.items())})
    print('per-rank EQDYNA_PROFILE=1 (ms/step): %s'
         % {r: round(v, 4) for r, v in sorted(on_per_rank.items())})
    print('mean EQDYNA_PROFILE=0: %8.4f ms/step' % off_mean)
    print('mean EQDYNA_PROFILE=1: %8.4f ms/step' % on_mean)
    print('overhead (mean):       %8.4f ms/step' % overhead_ms)
    print('noise floor (%d repeats of profile=0, mean-per-rank): %.4f '
         'ms/step (samples: %s)' % (a.repeats, floor_ms,
                                     ['%.4f' % s for s in floor_samples]))
    print('budget: max(%.2f%% of profile=0 = %.4f ms/step, noise floor '
         '%.4f ms/step) = %.4f ms/step' % (REL_OVERHEAD_BUDGET * 100,
                                           budget_ms, floor_ms, limit_ms))

    if overhead_ms > limit_ms:
        print('FAIL profile_overhead: overhead %.4f ms/step exceeds budget '
             '%.4f ms/step' % (overhead_ms, limit_ms))
        return 1
    print('SUCCESS profile_overhead')
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--case', required=True)
    p.add_argument('--backend', required=True,
                   choices=('python-numpy', 'python-jax') + MPI_BACKENDS)
    p.add_argument('--cpus', required=True,
                   help='comma-separated cpu list, e.g. 6,7 -- required, no '
                        'default placement (module docstring). For an MPI '
                        'backend, the first --ranks entries are used, one '
                        'per rank.')
    p.add_argument('--ranks', type=int, default=4,
                   help='MPI ranks, fortran/python-jax-mpi only (default 4, '
                        'matching matrix.FORTRAN_RANKS/PY_MPI_RANKS for '
                        'test.tpv8). Ignored for the two serial backends.')
    p.add_argument('--n-lo', type=int, default=20)
    p.add_argument('--n-hi', type=int, default=60)
    p.add_argument('--repeats', type=int, default=3,
                   help='repeats of the profile=0 measurement, to establish '
                        'the noise floor (default 3)')
    p.add_argument('--busy-ceiling', type=float, default=0.1,
                   help='strict ceiling on the named cpus right now (default '
                        '0.1; run_scaling.py defaults to 0.45 for a single '
                        'per-step figure -- this is a difference-of-a-'
                        'difference and more sensitive to transient tenancy)')
    p.add_argument('--i-know-the-box-is-busy', action='store_true',
                   help='bypass the busy-ceiling refusal (explicit override, '
                        'never a silent default -- see numa.require_idle)')
    a = p.parse_args(argv)

    cpus = [int(c) for c in a.cpus.split(',') if c != '']
    if not cpus:
        raise SystemExit('FAIL: --cpus must name at least one cpu')

    node_map = numa.numa_topology()
    if not node_map:
        raise SystemExit('FAIL: numactl --hardware gave no topology; this '
                         'gate refuses to run unpinned')
    known = {c for cs in node_map.values() for c in cs}
    unknown = [c for c in cpus if c not in known]
    if unknown:
        raise SystemExit('FAIL: --cpus %s not in this box\'s topology %s'
                         % (unknown, sorted(known)))

    idle = numa.require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
    if idle is None:
        raise SystemExit(
            'FAIL: cpus %s are busier than --busy-ceiling=%.2f right now -- '
            'an overhead measurement on a busy box is not a measurement '
            '(pass --i-know-the-box-is-busy to override deliberately)'
            % (cpus, a.busy_ceiling))

    if a.backend in MPI_BACKENDS:
        if len(cpus) < a.ranks:
            raise SystemExit('FAIL: --cpus names %d cpu(s) but --ranks=%d '
                             'needs one each' % (len(cpus), a.ranks))
        return main_mpi(a, cpus, node_map)

    case_dir = rs.build_py_case(a.case)

    off_ps = per_step_ms(case_dir, cpus, node_map, a.backend, '0', a.n_lo, a.n_hi)
    on_ps = per_step_ms(case_dir, cpus, node_map, a.backend, '1', a.n_lo, a.n_hi)
    floor_ms, floor_samples = noise_floor_ms(case_dir, cpus, node_map, a.backend,
                                             a.n_lo, a.n_hi, a.repeats)

    overhead_ms = on_ps - off_ps
    budget_ms = REL_OVERHEAD_BUDGET * off_ps
    limit_ms = max(budget_ms, floor_ms)

    print('EQDYNA_PROFILE=0: %8.4f ms/step' % off_ps)
    print('EQDYNA_PROFILE=1: %8.4f ms/step' % on_ps)
    print('overhead:         %8.4f ms/step' % overhead_ms)
    print('noise floor (%d repeats of profile=0): %.4f ms/step (samples: %s)'
          % (a.repeats, floor_ms, ['%.4f' % s for s in floor_samples]))
    print('budget: max(%.2f%% of profile=0 = %.4f ms/step, noise floor '
          '%.4f ms/step) = %.4f ms/step' % (REL_OVERHEAD_BUDGET * 100,
                                            budget_ms, floor_ms, limit_ms))

    if overhead_ms > limit_ms:
        print('FAIL profile_overhead: overhead %.4f ms/step exceeds budget '
              '%.4f ms/step' % (overhead_ms, limit_ms))
        return 1
    print('SUCCESS profile_overhead')
    return 0


if __name__ == '__main__':
    sys.exit(main())
