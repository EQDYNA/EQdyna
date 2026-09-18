#! /usr/bin/env python3
"""
Item 40 re-check: pinned-single-core numpy-vs-jax per-step cost on
test.tpv29, independently re-run (the standing 4.0x figure was relayed via
another session's report and never re-run in this one).

Reuses `run_perf.py`'s pinned-single-core + per-step-by-difference machinery
(`pinned_env`, `build_perf_case`, `steady_state_per_step`, `verify_affinity`)
by repointing its PERF_CASE_NAME/PERF_CASE module globals at test.tpv29
instead of duplicating that logic (rule 1) -- run_perf.py's own
ceiling/pinning code is untouched, imported only.

Adds one thing run_perf.py does NOT have: a per-cpu busy-fraction readout
(imported from run_numa_scaling.cpu_busy_fractions, not reimplemented) on the
pinned core before each measurement, and a refusal at the same 0.2 default
ceiling every other tool in this sweep uses -- report-only, never overridden,
consistent with the session's ground rules for this dispatch.

n_lo/n_hi are the small (20, 60) step counts used elsewhere in this session's
tooling (run_scaling.py's own default), NOT tpv29's own full nstep (term=20s
at tpv29's dx would be thousands of steps -- far too slow pinned to 1 core).
Per-step-by-difference cancels the fixed cost regardless of which two step
counts are chosen; the case's own eventual physics horizon is irrelevant to a
per-step timing.

Usage:
    python3 testsys/perf/run_tpv29_pinned_compare.py [--core 0]
        [--n-lo 20] [--n-hi 60] [--busy-ceiling 0.2] [--i-know-the-box-is-busy]
"""
import argparse
import json
import os
import sys
import time

TESTSYS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TESTSYS)
import run_perf  # noqa: E402 (reused, not duplicated)
sys.path.insert(0, os.path.dirname(TESTSYS))
import run_numa_scaling as numa  # noqa: E402 (cpu_busy_fractions, not reimplemented)

OUT = os.path.join(TESTSYS, 'tpv29_pinned_last.json')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--core', default='0')
    ap.add_argument('--n-lo', type=int, default=20)
    ap.add_argument('--n-hi', type=int, default=60)
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    a = ap.parse_args()

    # Repoint run_perf's module globals at test.tpv29 instead of test.tpv8.
    run_perf.PERF_CASE_NAME = 'test.tpv29'
    run_perf.PERF_CASE = os.path.join(run_perf.TESTSYS, 'perf_case_tpv29', 'test.tpv29')
    run_perf.CORE = a.core

    print(f'Building a fresh {run_perf.PERF_CASE_NAME} case at {run_perf.PERF_CASE} ...')
    run_perf.build_perf_case()

    cpus = [int(a.core)]
    print(f'\nchecking cpu {cpus} busy fraction before measuring (ceiling {a.busy_ceiling}) ...')
    chk = numa.require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
    result = dict(case=run_perf.PERF_CASE_NAME, core=a.core, n_lo=a.n_lo, n_hi=a.n_hi,
                  busy_ceiling=a.busy_ceiling, overridden=a.i_know_the_box_is_busy,
                  timestamp=time.strftime('%Y-%m-%d %H:%M:%S'))
    if chk is None:
        result['status'] = 'SKIPPED (cpu busy, not overridden)'
        json.dump(result, open(OUT, 'w'), indent=2)
        print(f'\nSKIPPED: cpu {cpus} busier than ceiling -- see message above. Wrote {OUT}.')
        return 1
    load1, load5, load15, others, busy = chk
    result['load_avg'] = [load1, load5, load15]
    result['other_users'] = others
    result['cpu_busy'] = busy

    affinity = run_perf.verify_affinity()
    print(f'Pinned-core affinity check (taskset -c {a.core}): {affinity}')
    result['affinity_check'] = affinity

    print(f'\nmeasuring numpy pinned per-step ({a.n_lo} vs {a.n_hi} steps) ...')
    ps_numpy, fixed_numpy = run_perf.steady_state_per_step('numpy', a.n_lo, a.n_hi)
    print(f'\nre-checking cpu {cpus} busy fraction before the jax measurement ...')
    chk2 = numa.require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
    result['cpu_busy_before_jax'] = chk2[4] if chk2 else None
    if chk2 is None:
        result['status'] = 'PARTIAL (numpy measured, jax SKIPPED -- cpu went busy)'
        result['numpy_ms_per_step'] = None if ps_numpy is None else ps_numpy * 1e3
        json.dump(result, open(OUT, 'w'), indent=2)
        print(f'\nPARTIAL: cpu went busy before jax could be measured. Wrote {OUT}.')
        return 1
    print(f'\nmeasuring jax pinned per-step ({a.n_lo} vs {a.n_hi} steps) ...')
    ps_jax, fixed_jax = run_perf.steady_state_per_step('jax', a.n_lo, a.n_hi)

    if ps_numpy is None or ps_jax is None:
        result['status'] = 'FAILED'
        result['numpy_ms_per_step'] = None if ps_numpy is None else ps_numpy * 1e3
        result['jax_ms_per_step'] = None if ps_jax is None else ps_jax * 1e3
        json.dump(result, open(OUT, 'w'), indent=2)
        raise SystemExit(f'FAIL: one engine did not produce a timing (numpy={ps_numpy}, jax={ps_jax})')

    ratio = ps_numpy / ps_jax
    result.update(status='OK', numpy_ms_per_step=ps_numpy * 1e3, jax_ms_per_step=ps_jax * 1e3,
                  numpy_fixed_s=fixed_numpy, jax_fixed_s=fixed_jax,
                  numpy_over_jax=ratio)
    print(f'\n==== test.tpv29 pinned single-core (cpu {a.core}), per-step ({a.n_lo} vs {a.n_hi} steps) ====')
    print(f'  numpy: {ps_numpy*1e3:9.2f} ms/step  (fixed {fixed_numpy:.2f} s)')
    print(f'  jax:   {ps_jax*1e3:9.2f} ms/step  (fixed {fixed_jax:.2f} s)')
    print(f'  numpy/jax ratio: {ratio:.3f}x   (item 40\'s relayed figure: 4.0x)')
    json.dump(result, open(OUT, 'w'), indent=2)
    print(f'\nwrote {OUT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
