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


def timed_run(case_dir, nsteps, cpus, node_map, backend, profile_flag):
    """One numactl-pinned, fresh-process solve of `nsteps` steps, with
    EQDYNA_PROFILE forced to `profile_flag` ('0' or '1'). Mirrors
    run_scaling.time_one_py's subprocess/env discipline exactly (same reason:
    a fresh process per timing so neither jax's compile cache nor a stale
    thread pool carries over between the two arms being compared)."""
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from eqdyna import eqdyna3d\n"
        "import time; t0 = time.time()\n"
        "eqdyna3d.run_case(%r, nsteps=%d, verbose=False, backend=%r)\n"
        "print('WALL', time.time() - t0)\n"
        % (rs.PYTHON_PKG, case_dir, nsteps, backend))
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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--case', required=True)
    p.add_argument('--backend', required=True, choices=('python-numpy', 'python-jax'))
    p.add_argument('--cpus', required=True,
                   help='comma-separated cpu list, e.g. 6,7 -- required, no '
                        'default placement (module docstring)')
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
