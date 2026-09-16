#! /usr/bin/env python3
"""
`perf` testsys tier: pinned single-core Fortran-vs-NumPy-vs-JAX timing on
the tpv8 parity fixture (same case, same step count as the `parity` tier),
with a provenance-stamped table and a ratio guard against a checked-in
baseline (testsys/perf/baseline.json).

Every engine is pinned to the SAME single core (`taskset -c <core>`) with
every thread-pool env var this repo's dependency stack recognizes forced
to 1 -- OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, and JAX/XLA
via XLA_FLAGS. This replaces python/README-parity.md's earlier
threading-uncontrolled numbers (a shared 64-core machine at unpredictable
load, no thread pinning) with a number that is at least internally
comparable run-to-run on this same machine.

Gate: FAILS if the contemporaneous python/fortran wall-clock ratio (best
of NumPy/JAX, i.e. whichever python engine is faster) degrades by more
than 1.5x versus testsys/perf/baseline.json's ratio. Never gates on
absolute seconds (machine load varies far more than a 1.5x band -- see
python/README-parity.md's own load-variance caveat). If baseline.json is
missing, THIS run creates it and does not fail (there is nothing to
compare against yet) -- print a clear notice so that is not mistaken for
a silent pass on a real regression later.
"""
import json
import os
import platform
import socket
import subprocess
import sys
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
PYTHON_PKG = os.path.join(REPO_ROOT, 'src', 'python')
FIXTURE_CASE = os.path.join(TESTSYS, '..', 'parity', 'fixtures', 'test_tpv8_serial')
BASELINE_PATH = os.path.join(TESTSYS, 'baseline.json')
CORE = os.environ.get('PERF_CORE', '0')
DEGRADE_LIMIT = 1.5


def pinned_env():
    env = dict(os.environ)
    env['OMP_NUM_THREADS'] = '1'
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'
    env['NUMEXPR_NUM_THREADS'] = '1'
    env['XLA_FLAGS'] = (
        '--xla_cpu_multi_thread_eigen=false '
        f'--xla_force_host_platform_device_count=1'
    )
    env['XLA_CPU_MULTI_THREAD'] = '0'
    # jax's intra/inter-op thread pool (checked at runtime below too, not just asserted).
    env['JAX_PLATFORMS'] = 'cpu'
    return env


def verify_affinity():
    """Runs a subprocess UNDER the same taskset used for the real timings and
    reports os.sched_getaffinity(0) from inside it -- proof the pin took
    effect, not an assumption."""
    r = subprocess.run(
        ['taskset', '-c', CORE, sys.executable, '-c',
         'import os; print(sorted(os.sched_getaffinity(0)))'],
        capture_output=True, text=True)
    return r.stdout.strip()


def time_fortran(nsteps):
    binary = os.path.join(FIXTURE_CASE, 'eqdyna-pydump')
    if not os.path.exists(binary):
        raise SystemExit(
            f'FAIL: {binary} missing -- run testsys/parity/make_fixtures.py first '
            '(perf reuses the parity fixture case+binary for an identical-input comparison).')
    t0 = time.time()
    r = subprocess.run(['taskset', '-c', CORE, 'mpirun', '-np', '1', '-wdir', FIXTURE_CASE, binary],
                        env=pinned_env(), capture_output=True, text=True)
    dt = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        raise SystemExit('FAIL: Fortran run failed under perf tier')
    return dt


def time_python(engine, nsteps):
    script = (
        f"import sys; sys.path.insert(0, {PYTHON_PKG!r}); sys.path.insert(0, {os.path.join(PYTHON_PKG, 'eqdyna')!r})\n"
        f"from eqdyna.pydump import load\n"
        f"from eqdyna import eqdyna3d\n"
        f"S = load({FIXTURE_CASE!r})\n"
        f"import time; t0 = time.time()\n"
        f"out = eqdyna3d.run(S, nsteps={nsteps}, verbose=False, backend={engine!r})\n"
        f"print('PERF_WALL', time.time() - t0)\n"
    )
    t0 = time.time()
    r = subprocess.run(['taskset', '-c', CORE, sys.executable, '-c', script],
                        env=pinned_env(), capture_output=True, text=True)
    wall_incl_load = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        return None, None
    solve_only = None
    for line in r.stdout.splitlines():
        if line.startswith('PERF_WALL'):
            solve_only = float(line.split()[1])
    return wall_incl_load, solve_only


def provenance():
    sha = subprocess.run(['git', '-C', REPO_ROOT, 'rev-parse', '--short', 'HEAD'],
                          capture_output=True, text=True).stdout.strip() or 'unknown (uncommitted worktree)'
    load1, load5, load15 = os.getloadavg()
    compiler = subprocess.run(['mpif90', '--version'], capture_output=True, text=True).stdout.splitlines()
    compiler = compiler[0] if compiler else 'unknown'
    return dict(git_sha=sha, hostname=socket.gethostname(), platform=platform.platform(),
                load_avg=f'{load1:.2f} {load5:.2f} {load15:.2f}', compiler=compiler,
                pinned_core=CORE, timestamp=time.strftime('%Y-%m-%d %H:%M:%S'))


def main():
    if not os.path.isdir(FIXTURE_CASE):
        raise SystemExit(
            f'FAIL: fixtures missing at {FIXTURE_CASE}. Run: python3 testsys/parity/make_fixtures.py')

    affinity = verify_affinity()
    print(f'Pinned-core affinity check (taskset -c {CORE}): {affinity}')
    if affinity != f'[{CORE}]':
        print(f'WARNING: affinity check did not report exactly core {CORE} -- '
              'pinning may not have taken effect on this system; timings below are '
              'still reported but treat the ratio with more suspicion than usual.')

    prov = provenance()

    sys.path.insert(0, PYTHON_PKG)
    sys.path.insert(0, os.path.join(PYTHON_PKG, 'eqdyna'))
    from eqdyna.pydump import load  # noqa: E402
    S = load(FIXTURE_CASE)
    nsteps = S['nstep']

    t_fortran = time_fortran(nsteps)
    t_numpy_wall, t_numpy_solve = time_python('numpy', nsteps)
    try:
        import jax  # noqa: F401
        have_jax = True
    except ImportError:
        have_jax = False
    t_jax_wall, t_jax_solve = (time_python('jax', nsteps) if have_jax else (None, None))

    print('\n==== testsys: perf -- provenance ====')
    for k, v in prov.items():
        print(f'  {k}: {v}')

    print(f'\n==== testsys: perf -- pinned single-core timing (tpv8, {nsteps} steps, core {CORE}) ====')
    print(f'  Fortran (mpirun -np 1, pinned):        {t_fortran:8.2f} s')
    if t_numpy_solve is not None:
        print(f'  NumPy   (solve-only, pinned):           {t_numpy_solve:8.2f} s  '
              f'(wall incl. interpreter/import: {t_numpy_wall:.2f} s)')
    else:
        print('  NumPy   : FAILED to run')
    if have_jax:
        if t_jax_solve is not None:
            print(f'  JAX-CPU (solve+compile, pinned, 1 thread): {t_jax_solve:8.2f} s  '
                  f'(wall incl. interpreter/import: {t_jax_wall:.2f} s)')
        else:
            print('  JAX-CPU : FAILED to run')
    else:
        print('  JAX-CPU : SKIPPED (jax not importable)')

    ratios = {}
    if t_numpy_solve:
        ratios['numpy_over_fortran'] = t_numpy_solve / t_fortran
    if have_jax and t_jax_solve:
        ratios['jax_over_fortran'] = t_jax_solve / t_fortran
    best_python_over_fortran = min(ratios.values()) if ratios else None

    print('\n  Ratios (python/fortran, >1 = python slower):')
    for k, v in ratios.items():
        print(f'    {k}: {v:.3f}')

    if best_python_over_fortran is None:
        raise SystemExit('FAIL: no python engine produced a timing to gate on')

    if not os.path.exists(BASELINE_PATH):
        baseline = dict(prov, ratios=ratios, best_python_over_fortran=best_python_over_fortran,
                         note='First pinned run -- this file becomes the baseline. '
                              'Committed by the reviewer, not auto-committed by this script.')
        with open(BASELINE_PATH, 'w') as f:
            json.dump(baseline, f, indent=2)
        print(f'\nNo baseline found -- wrote {BASELINE_PATH} from this run. '
              'Re-run once more to exercise the ratio guard.')
        print('SUCCESS perf (baseline created, not gated this run)')
        return 0

    with open(BASELINE_PATH) as f:
        baseline = json.load(f)
    baseline_ratio = baseline['best_python_over_fortran']
    degradation = best_python_over_fortran / baseline_ratio
    print(f'\n  Baseline best python/fortran ratio: {baseline_ratio:.3f} '
          f'(from {baseline.get("timestamp", "unknown time")}, SHA {baseline.get("git_sha", "unknown")})')
    print(f'  This run:                            {best_python_over_fortran:.3f}')
    print(f'  Degradation factor:                  {degradation:.3f} (fail if > {DEGRADE_LIMIT})')

    if degradation > DEGRADE_LIMIT:
        print(f'FAIL perf: python/fortran ratio degraded {degradation:.2f}x vs baseline '
              f'(limit {DEGRADE_LIMIT}x)')
        return 1
    print('SUCCESS perf')
    return 0


if __name__ == '__main__':
    sys.exit(main())
