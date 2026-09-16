#! /usr/bin/env python3
"""
Regenerates, FROM SCRATCH, the two release-audit-flagged number sets behind
pathway_forward.md item 21 -- numbers previously cited only in committed
docstrings/commit messages, with no committed script to reproduce them
(PROJECT_RULES rule 4 "only fresh runs are evidence" / rule 6 "every
performance number carries its provenance"):

  (a) drv.a6 decomposition-chaos evidence behind Milestone 10's two-part
      acceptance criterion (testsys/matrix.py's DRV_A6 bounds; full
      derivation in pathway_forward.md item 18): the A/B/C existence-flip +
      timing-shift decomposition --
        A: fresh serial Fortran        vs committed 4-rank reference
        B: standalone Python (JAX)     vs that SAME fresh serial Fortran run
        C: standalone Python (JAX)     vs committed 4-rank reference
      -- plus the A/B flip-set overlap and the matched-arrival median
      |Delta fnft| in all three pairs.
  (b) (--gpu / --gpu-only) the gpu-tier timing comparison: standalone tpv8
      on GPU (JAX_PLATFORMS=cuda) vs CPU (JAX_PLATFORMS=cpu), run
      contemporaneously back-to-back on this box, with an explicit
      nvidia-smi contention caveat printed at the start.

REPORT-ONLY. This script is never a gate and is never wired into
testsys/run.py's unit/regression/e2e tiers -- it is invoked BY HAND,
exactly like testsys/parity/make_fixtures.py. It imports (does NOT
duplicate) make_serial_case / run_standalone from testsys/e2e/run_e2e.py
and load_coordinate_aligned / align_two_frt_files / flip_decomposition /
drv_a6_gate / abs_max_diff from testsys/compare.py, with the bounds from
testsys/matrix.py (rule 1) -- one implementation of "how do two frt runs
compare", every consumer. Those helpers used to live in
testsys/parity/test_standalone_acceptance.py, the separate `accept` tier
that has since been folded into the e2e sweep as a backend column.

THE NUMBERS BELOW WILL NOT MATCH 372/439/391 EXACTLY, and that is
expected, not a bug: matrix.py's DRV_A6 comments and pathway_forward.md
item 18 document that drv.a6's rupture-front arrival is genuinely
decomposition/reduction-order-sensitive (a PURE-Fortran, same-binary,
same-algorithm run flips ~7% of fault nodes' rupture arrival on a mere
MPI-decomposition change). The recorded 372/439/391 are ONE sample from
one session; matrix.DRV_A6['total_flip_bound']=450
is calibrated WITH HEADROOM above that sample, not an exact target this
script's re-run must reproduce. Do NOT edit DRV_A6_TOTAL_FLIP_BOUND (or
any other committed bound) to make a new sample from this script "match" --
if a fresh run overshoots 450, that is itself evidence worth reporting to
a human, not a reason to loosen the bound in the same commit (PROJECT_RULES
rule 5 / this repo's own anti-charter on that point).

COST / MACHINE COURTESY: drv.a6 serial Fortran build+run is ~5 minutes;
the standalone Python (JAX, CPU-pinned) leg is ~10 minutes; each --gpu
leg is ~1 minute. Check `uptime` FIRST. If 1-minute load average exceeds
~40 on a 64-core box, do NOT run this fresh -- inspect whatever
frt.txt0 / frt.txt0.fortran artifacts already exist on disk from a prior
invocation instead, and say so explicitly in your report rather than
adding load to a busy shared machine. Never run drv.a6's Fortran and
Python legs concurrently with each other or with any other tier
(run_parity.py / make_fixtures.py / the e2e sweep already document this
same one-at-a-time constraint).

Run:
    python3 testsys/parity/evidence_drv_a6_chaos.py            # part (a) only
    python3 testsys/parity/evidence_drv_a6_chaos.py --gpu       # (a) then (b)
    python3 testsys/parity/evidence_drv_a6_chaos.py --gpu-only  # (b) only

Writes a timestamped JSON snapshot to testsys/parity/evidence_output/
(gitignored -- a fresh chaotic sample is evidence-of-a-run, not golden
reference data; rule 4/rule 7) alongside the printed table.
"""
import argparse
import datetime
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time

import numpy as np

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
OUT_DIR = os.path.join(TESTSYS, 'evidence_output')

sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, 'testsys', 'e2e'))
from testsys import compare, matrix  # rule 1: reuse, don't duplicate
import run_e2e                       # the sweep's cell runners (same rule)

# CPU-pinned by default for part (a), exactly as the old accept tier was: a
# run labelled jax-on-cpu that silently landed on a contended GPU is a
# different measurement under the same name. Part (b) overrides per leg.
PLATFORM = os.environ.get('JAX_PLATFORMS', 'cpu')


# ---------------------------------------------------------------- provenance

def provenance():
    def first_line(cmd):
        r = subprocess.run(cmd, capture_output=True, text=True)
        out = (r.stdout or r.stderr).splitlines()
        return out[0] if out else '(unavailable)'

    sha = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=REPO_ROOT,
                          capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', 'status', '--porcelain'], cwd=REPO_ROOT,
                                 capture_output=True, text=True).stdout.strip())
    load1, load5, load15 = os.getloadavg()
    return dict(
        date_utc=datetime.datetime.utcnow().isoformat() + 'Z',
        host=socket.gethostname(), sha=sha, dirty=dirty,
        loadavg=f'{load1:.2f} {load5:.2f} {load15:.2f}', ncpu=os.cpu_count(),
        gfortran=first_line(['gfortran', '--version']),
        mpif90=first_line(['mpif90', '--version']),
    )


def print_provenance(prov, title):
    print(f'\n==== {title}: provenance ====')
    for k in ('date_utc', 'host', 'sha', 'dirty', 'loadavg', 'ncpu', 'gfortran', 'mpif90'):
        print(f'  {k:9s}: {prov[k]}')


# ---------------------------------------------------------------- part (a)

def build_serial_fortran_binary():
    """Clean rebuild of the DEFAULT src/eqdyna target (MACHINE=ubuntu) --
    matches drv_a6_gate's docstring's "a freshly-built src/eqdyna (default
    target, MACHINE=ubuntu)". Rebuilt every invocation (rule 4 -- a stale
    bin/eqdyna already caused one diagnosis detour this project, see
    test_standalone_acceptance.py's module docstring)."""
    src = os.path.join(REPO_ROOT, 'src', 'fortran')
    env = dict(os.environ, MACHINE='ubuntu')
    subprocess.run(['bash', '-c', 'rm -f *.o eqdyna'], cwd=src)
    t0 = time.time()
    r = subprocess.run(['make', 'MACHINE=ubuntu'], cwd=src, env=env,
                        capture_output=True, text=True)
    elapsed = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit(f'FAIL: make MACHINE=ubuntu exited {r.returncode}')
    binary = os.path.join(src, 'eqdyna')
    if not os.path.exists(binary):
        raise SystemExit('FAIL: src/eqdyna was not produced by the build')
    print(f'built src/eqdyna fresh in {elapsed:.1f}s')
    return binary


def run_fortran_serial(binary, case_dir):
    dst = os.path.join(case_dir, 'eqdyna')
    shutil.copy(binary, dst)
    t0 = time.time()
    r = subprocess.run(['mpirun', '-np', '1', dst], cwd=case_dir,
                        capture_output=True, text=True)
    elapsed = time.time() - t0
    if r.returncode != 0:
        print(r.stdout[-4000:]); print(r.stderr[-4000:])
        raise SystemExit(f'FAIL: serial Fortran eqdyna exited {r.returncode} on {case_dir}')
    frt = os.path.join(case_dir, 'frt.txt0')
    if not os.path.isfile(frt):
        raise SystemExit(f'FAIL: {frt} was not written by the serial Fortran run')
    print(f'serial Fortran test.drv.a6 run: {elapsed:.1f}s')
    return frt, elapsed


def _flip_overlap(mask_a, mask_b):
    return int(np.sum(mask_a & mask_b)), int(np.sum(mask_a | mask_b))


def regenerate_drv_a6_chaos():
    """Part (a): fresh serial Fortran build+run, fresh standalone-Python
    (JAX, CPU-pinned -- test_standalone_acceptance's import already sets
    JAX_PLATFORMS=cpu by default) run on the SAME case-input files, then
    A/B/C computed via THIS module's imported helpers only -- no
    duplicated comparison logic (rule 1)."""
    case_name = 'test.drv.a6'
    tmp = tempfile.mkdtemp(prefix='evidence_drv_a6_')
    try:
        case_dir = os.path.join(tmp, case_name)
        run_e2e.make_serial_case(case_name, case_dir, run_e2e.base_env())

        binary = build_serial_fortran_binary()
        fortran_frt, fortran_elapsed = run_fortran_serial(binary, case_dir)
        fortran_frt_saved = fortran_frt + '.fortran'
        shutil.copy(fortran_frt, fortran_frt_saved)  # preserve before python overwrites frt.txt0

        t0 = time.time()
        # overwrites frt.txt0 with the python output
        py_frt = run_e2e.run_standalone(case_dir, 'python-jax', device=PLATFORM)
        python_elapsed = time.time() - t0
        print(f'standalone Python (JAX, JAX_PLATFORMS={PLATFORM}) '
              f'test.drv.a6 run: {python_elapsed:.1f}s')

        # A: fresh serial Fortran vs committed reference.
        ok_a, diag_a = compare.drv_a6_gate(fortran_frt_saved)
        # C: standalone Python vs committed reference (the pair the sweep's
        # test.drv.a6 python cells gate on).
        ok_c, diag_c = compare.drv_a6_gate(py_frt)
        # B: standalone Python vs fresh serial Fortran -- NEITHER side is
        # the committed reference, so this needs align_two_frt_files, not
        # load_coordinate_aligned (which always loads test.reference.results/).
        ref_s_b, py_s_b = compare.align_two_frt_files(fortran_frt_saved, py_frt)
        diag_b = compare.flip_decomposition(ref_s_b, py_s_b)
        ok_b = bool(diag_b['ok_median'] and diag_b['ok_phys'] and diag_b['ok_flips'])

        overlap_ab, union_ab = _flip_overlap(diag_a['flip_mask'], diag_b['flip_mask'])

        def strip_mask(d):
            return {k: v for k, v in d.items() if k != 'flip_mask'}

        return dict(
            A=dict(label='fresh-serial-Fortran vs committed 4-rank reference',
                    ok=ok_a, **strip_mask(diag_a)),
            B=dict(label='standalone-python vs fresh-serial-Fortran',
                    ok=ok_b, **strip_mask(diag_b)),
            C=dict(label='standalone-python vs committed 4-rank reference',
                    ok=ok_c, **strip_mask(diag_c)),
            overlap_A_and_B=overlap_ab, union_A_or_B=union_ab,
            fortran_elapsed_s=fortran_elapsed, python_elapsed_s=python_elapsed,
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def print_drv_a6_table(r):
    print('\n==== part (a): drv.a6 decomposition-chaos evidence (FRESH SAMPLE) ====')
    print(f'  Fortran build+run: {r["fortran_elapsed_s"]:.1f}s   '
          f'Python (JAX/CPU) run: {r["python_elapsed_s"]:.1f}s')
    print(f'  {"":<3} {"total_flips":>12} {"/bound":>7} {"only_ref":>9} {"only_py":>8} '
          f'{"timing_shifts":>14} {"ruptured_both":>14} {"median|d fnft|":>15} {"phys_max_diff":>14}')
    for key in ('A', 'B', 'C'):
        d = r[key]
        print(f'  {key:<3} {d["total_flips"]:>12d} '
              f'{"/" + str(matrix.DRV_A6["total_flip_bound"]):>7} '
              f'{d["n_only_ref"]:>9d} {d["n_only_run"]:>8d} {d["n_timing_shifts"]:>14d} '
              f'{d["n_ruptured_both"]:>14d} {d["median_fnft_diff"]:>15.4f} {d["phys_max_diff"]:>14.4e}'
              f'   [{d["label"]}]')
    print(f'  A-intersect-B (flip-set overlap) = {r["overlap_A_and_B"]}; '
          f'A-union-B = {r["union_A_or_B"]}')
    print('  NOTE: these are CHAOTIC-SAMPLE quantities (see module docstring) -- the '
          'recorded values (A=372, B=439, C=391, overlap=212, union=599, median '
          '0.0417s in all three; pathway_forward.md item 18) are ONE sample; '
          'matrix.DRV_A6["total_flip_bound"]=450 is calibrated WITH HEADROOM above '
          'them, not an exact target this run is expected to reproduce.')
    for key in ('A', 'B', 'C'):
        d = r[key]
        balance = 'balanced' if d['n_only_ref'] == 0 and d['n_only_run'] == 0 else (
            'balanced' if 0.3 <= (d['n_only_ref'] + 1) / (d['n_only_run'] + 1) <= 3.3
            else 'ONE-SIDED (worth a second look)')
        print(f'  {key}: only_ref={d["n_only_ref"]} vs only_run={d["n_only_run"]} -> {balance}')


# ---------------------------------------------------------------- part (b)

def regenerate_gpu_timing():
    """Part (b): contemporaneous back-to-back standalone tpv8 timing on
    GPU (JAX_PLATFORMS=cuda) vs CPU (JAX_PLATFORMS=cpu), reusing
    make_serial_case/coordinate_aligned_diff from
    test_standalone_acceptance.py. Prints nvidia-smi utilization AT START
    as an explicit contention caveat: this is a shared box, and a timing
    ratio measured here is only as good as the other tenants' activity at
    that moment -- exactly why the accept tier's `run.py accept` pins
    JAX_PLATFORMS=cpu for its PASS/FAIL decision instead of depending on
    GPU wall-clock at all."""
    case_name = 'test.tpv8'

    nv = subprocess.run(
        ['nvidia-smi', '--query-gpu=index,name,utilization.gpu,memory.used,memory.total',
         '--format=csv'], capture_output=True, text=True)
    nvidia_smi_at_start = nv.stdout if nv.returncode == 0 else '(nvidia-smi unavailable)'
    print('\n==== part (b): nvidia-smi at start (CONTENTION CAVEAT) ====')
    print(nvidia_smi_at_start)

    def run_platform(platform):
        tmp = tempfile.mkdtemp(prefix=f'evidence_gpu_{platform}_')
        try:
            case_dir = os.path.join(tmp, case_name)
            run_e2e.make_serial_case(case_name, case_dir, run_e2e.base_env())
            t0 = time.time()
            frt = run_e2e.run_standalone(case_dir, 'python-jax', device=platform)
            elapsed = time.time() - t0
            ref_a, run_a = compare.load_coordinate_aligned(case_name, frt)
            max_abs_diff, ok, _ = compare.abs_max_diff(case_name, ref_a, run_a)
            return elapsed, max_abs_diff, ok
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    # GPU leg first (so any one-time device/warm-up cost lands on the GPU
    # number, not spilled onto CPU's), then CPU, back-to-back.
    gpu_elapsed, gpu_diff, gpu_ok = run_platform('cuda')
    cpu_elapsed, cpu_diff, cpu_ok = run_platform('cpu')

    return dict(gpu_elapsed_s=gpu_elapsed, cpu_elapsed_s=cpu_elapsed,
                gpu_max_abs_diff=gpu_diff, cpu_max_abs_diff=cpu_diff,
                gpu_ok=gpu_ok, cpu_ok=cpu_ok,
                nvidia_smi_at_start=nvidia_smi_at_start)


def print_gpu_table(r):
    print('\n==== part (b): gpu vs cpu timing (FRESH SAMPLE, CONTEMPORANEOUS) ====')
    print(f'  {"platform":<10} {"wall_s":>10} {"max_abs_diff":>14} {"pass_bound":>10}')
    bound = matrix.CASE_BOUND['test.tpv8']
    print(f'  {"gpu":<10} {r["gpu_elapsed_s"]:>10.2f} {r["gpu_max_abs_diff"]:>14.4e} '
          f'{"ok" if r["gpu_ok"] else "FAIL":>10} (bound {bound:e})')
    print(f'  {"cpu":<10} {r["cpu_elapsed_s"]:>10.2f} {r["cpu_max_abs_diff"]:>14.4e} '
          f'{"ok" if r["cpu_ok"] else "FAIL":>10} (bound {bound:e})')
    ratio = r['cpu_elapsed_s'] / r['gpu_elapsed_s'] if r['gpu_elapsed_s'] else float('nan')
    print(f'  cpu/gpu wall-time ratio = {ratio:.2f}x')
    print('  CONTENTION CAVEAT: this ratio was measured on a SHARED box; see the '
          'nvidia-smi snapshot above for other tenants\' utilization at the moment '
          'this ran. A low or negative apparent speedup here does not necessarily '
          'contradict a prior idle-box measurement, and a large one should not be '
          'over-trusted either.')


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gpu', action='store_true',
                     help='also run part (b) after part (a)')
    ap.add_argument('--gpu-only', action='store_true',
                     help='skip part (a) (drv.a6 chaos, ~15 min); run only part (b) (~2 min)')
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    out = dict(generated_utc=datetime.datetime.utcnow().isoformat() + 'Z')

    if not args.gpu_only:
        prov_a = provenance()
        print_provenance(prov_a, 'part (a): drv.a6 chaos')
        results_a = regenerate_drv_a6_chaos()
        print_drv_a6_table(results_a)
        out['part_a_drv_a6_chaos'] = dict(provenance=prov_a, **results_a)

    if args.gpu or args.gpu_only:
        prov_b = provenance()
        print_provenance(prov_b, 'part (b): gpu vs cpu timing')
        results_b = regenerate_gpu_timing()
        print_gpu_table(results_b)
        out['part_b_gpu_vs_cpu'] = dict(provenance=prov_b, **results_b)

    ts = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
    out_path = os.path.join(OUT_DIR, f'evidence_{ts}.json')
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2, default=str)
    print(f'\nWrote {out_path}')
    print('\nThis script is REPORT-ONLY: it always exits 0 (see module docstring) -- '
          'it is evidence, not a gate. A number here landing outside a committed '
          'bound is a prompt to go tell a human, not something this script itself '
          'fails on.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
