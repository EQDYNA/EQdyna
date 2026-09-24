#! /usr/bin/env python3
"""
THE test. One sweep, two nested loops, one comparison:

    for case in CASES:                       # testNameList.py, 8 cases
        for backend in BACKENDS:             # fortran, python-numpy, python-jax
            run(case, backend) -> canonical frt -> compare vs the ONE
                                  committed reference, at THAT CASE's bound

`backend` is an axis of the sweep, not a tier. The former `accept` tier is
gone: it ran the SAME solver against the SAME references with a SECOND
comparison implementation, over 5 listed cases while this tier gated 8, on the
jax backend only because it invoked the standalone with no --backend. Its
"SUCCESS accept (5/5 cases)" therefore reported five of five LISTED -- 5 of 16
case x backend combinations -- and read like completeness. Both of its halves
live here now: the run half below, the comparison half in testsys/compare.py.

WHAT A GREEN RUN MEANS, EXACTLY
It means every cell this run PRINTED as "will RUN" was run and gated. The
coverage block is printed before the first case starts and restated in the
summary, so the result can never be read as broader than it is. Cells outside
the selection are named, and cells the table declares unsupported are named
with their reason.

DECLARED-UNSUPPORTED IS NOT A SKIP
  - default sweep (no --cases/--backends): the table's unsupported cells are
    printed with their reasons and are NOT counted as coverage.
  - explicit selection naming an unsupported cell: hard FAILURE. The caller
    expected that cell to run; "I could not check this" and "this is fine"
    must not share an exit code (PROJECT_RULES rule 2).

Named gates, in order:
  1. testsys/regression/test_create_newcase.py -- cheap check first (rule 9).
  2. Fresh build via ./install-eqdyna.sh (only when a fortran cell is selected;
     the previous bin/eqdyna is removed first so a stale binary cannot paper
     over a broken build). EQDYNA_E2E_BIN overrides with a pre-built binary
     when a concurrent long job must not have bin/ disturbed -- it skips the
     build-and-install path, so a green run under it says nothing about that
     path (rule 16); it is printed loudly.
  3. Every selected cell: run, then compare.

Rule 8: the previous test/ tree is rotated to test.prev/ rather than deleted,
and every cell's output stays on disk under test/ -- including the python
cells, which the old accept tier ran in a tempfile directory and deleted, so a
failure destroyed its own evidence.

Rule 21a: that rotation acts on a FIXED path, so exactly one invocation per
checkout may hold it. testsys/runlock.py takes an exclusive flock on
$REPO_ROOT/test before anything here builds or rotates; a second concurrent
invocation REFUSES, names the holder's pid and start time, and exits non-zero
rather than renaming the first one's live tree out from under it (pathway item
70 -- the 2026-09-22 collision that killed a 1500.1 s cell and printed FAIL).

No timeouts: this is a shared box and dynamic-rupture runs are slow when it is
busy.
"""
import argparse
import datetime
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

E2E_DIR = os.path.dirname(os.path.abspath(__file__))
TESTSYS = os.path.dirname(E2E_DIR)
REPO_ROOT = os.path.dirname(TESTSYS)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from testsys import compare, frt_canonical, matrix, runlock  # noqa: E402
sys.path.insert(0, TESTSYS)  # testsys/ itself (profile_record.py lives there)
import profile_record  # noqa: E402  (append-only per-rank profile totals)

# Line-buffered stdout. Redirected to a file or through `tee`, Python block-
# buffers its OWN prints while subprocess children write straight to the fd --
# so the coverage block, the per-cell verdicts and the child output interleave
# wrongly, and none of this script's own lines survive a kill. A GitHub runner
# SIGTERM'd one of these runs (exit 143) and the log showed a silence and an
# exit code. The block that states what a run covered is the last thing that
# should be lost when a run dies.
sys.stdout.reconfigure(line_buffering=True)

MACHINE = os.environ.get('EQDYNA_TEST_MACHINE', 'ubuntu')
MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
BIN_OVERRIDE = os.environ.get('EQDYNA_E2E_BIN')


# GPU cells (--device cuda). XLA PREALLOCATES a fraction of the card at
# process start -- 0.75 by default -- so two jax cells landing on the same
# device cannot both have it and the second dies with RESOURCE_EXHAUSTED.
# That is what `test.tpv36` x python-jax did in the 2026-09-22 GPU sweep
# (pathway item 58); run ALONE the same cell passes at 7.06e-09, well inside
# its 1e-6 bound, so it was a harness/allocator interaction and never a
# parity failure. The sweep runs cells concurrently by design, so the fix is
# here, not in the case: one concurrent jax cell per VISIBLE DEVICE, each
# pinned to its own card, with the preallocation stated explicitly instead of
# inherited from XLA's default.
GPU_MEM_FRACTION = os.environ.get('EQDYNA_E2E_GPU_MEM_FRACTION', '0.75')


def _run(cmd, cwd, env):
    print('+ (%s) %s' % (os.path.basename(cwd), ' '.join(cmd)))
    return subprocess.call(cmd, cwd=cwd, env=env)


def visible_gpu_indices():
    """The CUDA device indices this sweep may use, as a list of ints.

    EQDYNA_E2E_GPUS='0,2' overrides; otherwise `nvidia-smi` is asked. Raises
    if neither yields a device -- a --device cuda sweep that quietly fell back
    to one slot, or to the CPU, would report a GPU measurement it never made
    (rule 2: "could not check" must not read as "passed").
    """
    want = os.environ.get('EQDYNA_E2E_GPUS')
    if want:
        idx = [int(t) for t in want.replace(',', ' ').split()]
        if not idx:
            raise RuntimeError('EQDYNA_E2E_GPUS=%r names no device' % want)
        return idx
    r = subprocess.run(['nvidia-smi', '--query-gpu=index', '--format=csv,noheader'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(
            'nvidia-smi failed (rc=%d): %s -- a --device cuda sweep cannot be '
            'scheduled without knowing how many cards exist. Set '
            'EQDYNA_E2E_GPUS to name them explicitly.'
            % (r.returncode, (r.stderr or '').strip()[:200]))
    idx = [int(line) for line in r.stdout.split() if line.strip().isdigit()]
    if not idx:
        raise RuntimeError('nvidia-smi reported no CUDA devices; --device cuda '
                           'has nothing to run on')
    return idx


class GpuSlots:
    """One concurrent jax cell per CUDA device, each pinned to its own card.

    Same all-or-nothing discipline as the core budget below, and for the same
    reason: a cell must take a whole device or wait. Nothing here caps MEMORY
    -- pinning is what makes the 0.75 preallocation safe, because the card is
    then exclusive for the cell's lifetime.
    """

    def __init__(self, devices):
        import threading
        self._cond = threading.Condition()
        self._free = list(devices)
        self.devices = list(devices)

    def acquire(self):
        with self._cond:
            while not self._free:
                self._cond.wait()
            return self._free.pop(0)

    def release(self, index):
        with self._cond:
            self._free.append(index)
            self._cond.notify_all()


def gpu_env(env, index):
    """The two variables a pinned GPU cell runs under. Separate from the
    launch so a test can assert them without starting a solver."""
    env = dict(env)
    env['CUDA_VISIBLE_DEVICES'] = str(index)
    # Explicit, not inherited: XLA's implicit 0.75 is exactly what made two
    # concurrent cells on one card fail, and an implicit number cannot be
    # read off a log when the next sweep misbehaves.
    env['XLA_PYTHON_CLIENT_MEM_FRACTION'] = str(GPU_MEM_FRACTION)
    return env


# --------------------------------------------------------------------------
# longest-first scheduling (item 2, 2026-09-23)
# --------------------------------------------------------------------------
LEDGER_PATH = os.path.join(REPO_ROOT, 'docs', 'perf_ledger.jsonl')


def load_ledger_wall_costs(ledger_path=LEDGER_PATH):
    """Latest measured wall_s per (case, backend) cell, read from
    docs/perf_ledger.jsonl's 'cell-wall-clock' rows (tool == 'run_e2e').

    This is the MEASURED cost source for start-order scheduling below --
    chosen over matrix.MEASURED_PEAK_RSS_GB's wall-time comments because the
    ledger is machine-readable, append-only (testsys/perf/ledger.py), and
    already the record of every past sweep's per-cell wall clock.

    Every sweep now runs at the same GATE_TERM_S (2026-09-23: one term,
    everywhere), so a ledger row's (case, backend) cost is comparable across
    every past sweep regardless of which selection (everyday/release/CI)
    produced it; this function takes whichever row is LATEST for that
    (case, backend). Scheduling needs a relative ORDER, not an exact
    duration -- it would not be an acceptable way to report a timing.

    Returns {} if the ledger does not exist yet (a fresh checkout has none);
    callers must treat a missing entry as UNMEASURED, never as zero/cheap
    (see schedule_order).
    """
    costs = {}
    if not os.path.isfile(ledger_path):
        return costs
    with open(ledger_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get('tool') != 'run_e2e' or row.get('metric') != 'cell-wall-clock':
                continue
            wall = row.get('wall_s')
            case, backend = row.get('case'), row.get('backend')
            if wall is None or case is None or backend is None:
                continue
            costs[(case, backend)] = wall  # later lines overwrite earlier: latest wins
    return costs


def schedule_order(cells, cost_estimates):
    """Longest-estimated-cost-first ordering of `cells`: a slow cell must
    never be the one that starts last, because the concurrent sweep's wall
    clock is bounded by whichever cell is STILL RUNNING, not by when it
    started (see the module docstring's WHY THIS MATTERS on cell concurrency).

    cost_estimates: {(case, backend): seconds}, e.g. load_ledger_wall_costs().
    A cell with NO entry is UNMEASURED; the stated rule for that case is the
    conservative one -- schedule it FIRST (float('inf')), never silently
    treated as cheap and landed last.

    Stable sort: cells with equal (or equally-unmeasured) cost keep their
    relative order from the input list, so the schedule is deterministic for
    an unchanged ledger and cell list.
    """
    return sorted(cells, key=lambda cb: -cost_estimates.get(cb, float('inf')))


def cell_cost(case, backend):
    """Core cost this cell BILLS against the concurrency budget below --
    module-level (not a closure inside main()) so a regression guard can call
    it directly and so mutating its inputs (matrix.JAX_MEASURED_CORES,
    matrix.FORTRAN_RANKS, matrix.PY_MPI_RANKS) is observable without
    re-running main().

    Every branch cites where its number came from -- MEASURED, not a guess
    (rule 2):
      fortran         -- its real rank count, matrix.FORTRAN_RANKS
                         (testNameList.coreNumList).
      python-jax-mpi  -- its real rank count, matrix.PY_MPI_RANKS: genuinely
                         N processes, not one.
      python-jax      -- matrix.JAX_MEASURED_CORES, ceil'd. XLA's CPU
                         backend threads its own ops even in a serial launch
                         with no explicit thread pinning; billing it at 1 (as
                         this function used to, unconditionally, before the
                         2026-09-23 speed campaign) undercounts a MEASURED
                         252% CPU cell by ~2.5x -- see matrix.py's own
                         citation on that constant for the measurement.
                         Rounded UP: under-billing a real cost is the defect
                         this exists to close.
      python-numpy    -- 1, and this is a MEASURED number too, not merely an
                         assumption: eqdyna3d._narrow_numpy_affinity pins the
                         standalone numpy backend to exactly one cpu
                         (flat 1/2/4-cpu scaling measured there, and ~2x
                         WORSE spanning NUMA nodes with a wider mask), and
                         `/usr/bin/time -v` on test.tpv8 x python-numpy this
                         session read "Percent of CPU this job got: 101%"
                         over the full 197.5s run -- one core, not several.
    """
    if backend == 'fortran':
        return max(1, matrix.FORTRAN_RANKS.get(case, 4))
    if backend == 'python-jax-mpi':
        return max(1, matrix.PY_MPI_RANKS[case])
    if backend == 'python-jax':
        return max(1, math.ceil(matrix.JAX_MEASURED_CORES))
    return 1


# Margin subtracted from measured-free cores before it becomes the default
# --jobs budget (see measured_free_cores/default_jobs_budget below): headroom
# for this process's own bookkeeping (the ThreadPoolExecutor, subprocess
# launch/reap) and for a foreign job to grow slightly between the tenancy
# sample and the first cell actually starting. Same value the sweep's OLD
# fixed "cores-4" default already reserved -- re-purposed onto a MEASURED
# free-core base instead of the machine's total core count, not re-derived.
JOBS_MARGIN_CORES = 4


def measured_free_cores(total_cores=None, busy_ceiling=None, sample_s=2.0,
                        busy_fractions_fn=None):
    """(free, busy, total): cores measured FREE at sweep start, by sampling
    EVERY cpu's own /proc/stat busy fraction over >=2s and counting how many
    read busier than `busy_ceiling`.

    Reuses testsys/perf/run_numa_scaling.cpu_busy_fractions -- the SAME
    instrument testsys/perf's own tools already use to pick idle cpus for a
    measurement -- rather than a second /proc/stat reader; the default
    `busy_ceiling` is testsys/perf/ledger.TENANCY_REFERENCE_CEILING (0.5),
    the same ceiling every perf-ledger row already records tenancy against,
    not a new number invented for this one caller.

    `busy_fractions_fn` is an injection point ONLY: a caller (this module's
    own default-budget wiring) leaves it None and gets the real /proc/stat
    sampler; a regression guard passes a synthetic {cpu: fraction} producer
    to prove the budget actually moves with what this function measures,
    without waiting out a real sample or depending on this box's own load.

    Raises if per-cpu stats could not be read for ANY requested cpu id
    (cpu_busy_fractions returns {} in that case, e.g. no /proc/stat at all --
    macOS) -- rule 2: "could not measure the tenancy" must not silently
    become "assume the box is idle" or "assume it is full"; either guess
    would put a wrong number in front of every cell this sweep launches. Pass
    --jobs explicitly on a machine this cannot measure.

    Returns (free, busy, total, busy_ceiling) -- the ceiling actually used
    travels back with the numbers so a caller can print it, per this
    project's own papercuts rule against a bare, unexplained verdict.
    """
    total = total_cores if total_cores is not None else (os.cpu_count() or 4)
    if busy_ceiling is None:
        sys.path.insert(0, os.path.join(TESTSYS, 'perf'))
        import ledger
        busy_ceiling = ledger.TENANCY_REFERENCE_CEILING
    if busy_fractions_fn is None:
        sys.path.insert(0, os.path.join(TESTSYS, 'perf'))
        import run_numa_scaling as numa
        busy_fractions_fn = numa.cpu_busy_fractions
    fractions = busy_fractions_fn(list(range(total)), sample_s=sample_s)
    if not fractions:
        raise RuntimeError(
            'measured_free_cores: no per-cpu /proc/stat reading for any of '
            '%d cpu id(s) -- the tenancy-aware default --jobs budget has no '
            'honest number to compute from on this machine. Pass --jobs '
            'explicitly.' % total)
    busy = sum(1 for f in fractions.values() if f > busy_ceiling)
    free = total - busy
    return free, busy, total, busy_ceiling


def default_jobs_budget(cells, args):
    """The core budget for this invocation: `--jobs` if given (wins outright,
    no measurement taken -- a caller who names a number gets exactly it,
    same discipline as an explicit --cases/--backends selection elsewhere in
    this file), otherwise a TENANCY-AWARE default derived from cores measured
    FREE right now, never below the single most expensive selected cell (a
    computed budget that could not even fit its own largest cell would make
    every run of this sweep hang or starve that cell, which is worse than
    merely being a slow default).

    Prints the line the budget was derived from before returning it --
    papercuts.md's own rule: never a bare number with no measured quantity
    beside it.
    """
    largest = max((cell_cost(*cb) for cb in cells), default=1)
    if args.jobs:
        print('e2e: core budget %d (explicit --jobs; no tenancy sample taken)'
              % args.jobs)
        return args.jobs
    free, busy, total, busy_ceiling = measured_free_cores()
    computed = free - JOBS_MARGIN_CORES
    budget = max(largest, computed)
    print('e2e: core budget %d = max(largest selected cell cost %d, '
          'measured-free %d - margin %d) -- %d of %d cpus measured busy '
          '(>%.0f%% utilised over a >=2.0s /proc/stat sample, '
          'busy_ceiling=%.2f)'
          % (budget, largest, free, JOBS_MARGIN_CORES, busy, total,
             busy_ceiling * 100.0, busy_ceiling))
    return budget


def base_env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([
        os.path.join(REPO_ROOT, 'bin'),
        os.path.join(REPO_ROOT, 'scripts'),
        env.get('PATH', ''),
    ])
    return env


# --------------------------------------------------------------------------
# the ONE term (2026-09-23 owner decision: no second term, no --term flag)
# --------------------------------------------------------------------------
def apply_term_override(case_name, case_dir):
    """Append a par.term override to the just-copied user_defined_params.py.
    Called from EVERY backend's case-build path -- make_serial_case
    (python-numpy/jax/jax-mpi) and run_fortran -- right after create.newcase
    and before case.setup runs, so every cell goes through the ONE case-build
    path rather than forking a second one.

    Appends `par.term = matrix.GATE_TERM_S`, unconditionally, for every case,
    every selection (everyday, release, CI) -- there is no per-case "full"
    term to fall back to. case_input/<case_name>/user_defined_params.py's own
    committed par.term (what a user gets from create.newcase, run by hand,
    outside the gate) is left on disk untouched above this override; this
    function only appends the LAST assignment, which is what case.setup
    reads.
    """
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# term override by testsys/e2e/run_e2e.py: '
                       'every gated cell runs at matrix.GATE_TERM_S '
                       'regardless of %r\'s own committed par.term.\n'
                       'par.term = %r\n' % (case_name, matrix.GATE_TERM_S))


# --------------------------------------------------------------------------
# cell runners
# --------------------------------------------------------------------------
def make_serial_case(case_name, case_dir, env):
    """create.newcase + the GATE_TERM_S override + force a SERIAL
    decomposition + case.setup.

    The standalone solver refuses npx*npy*npz > 1, and the cases spell their
    decomposition three different ways: three separate `par.nx = 2` lines
    (tpv8, tpv10, tpv104), one tuple `par.nx, par.ny, par.nz = 2, 2, 1`
    (drv.a6, tpv29), or not at all (tpv1053d, which inherits
    defaultParameters). APPENDING the override handles all three by
    construction -- it is the last assignment executed either way -- which
    matters because reading the decomposition instead of overriding it has
    twice left a case multi-rank and made the solver's refusal look like a
    coverage gap. par.HPC_ncpu keeps its original value; it only feeds the
    generated run.sh, which the standalone never reads.

    No Fortran binary is involved: case.setup writes the bFile/netCDF case
    inputs the standalone reads natively.
    """
    rc = subprocess.call([sys.executable,
                          os.path.join(REPO_ROOT, 'scripts', 'create.newcase'),
                          case_dir, case_name], env=env)
    if rc != 0:
        raise RuntimeError('create.newcase %s exited %d' % (case_name, rc))
    apply_term_override(case_name, case_dir)
    params = os.path.join(case_dir, 'user_defined_params.py')
    with open(params) as f:
        text = f.read().rstrip('\n')
    with open(params, 'w') as f:
        f.write(text + '\n\n# forced serial by testsys/e2e/run_e2e.py'
                       ' (the standalone solver is serial-only)\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n')
    rc = subprocess.call([sys.executable, 'case.setup'], cwd=case_dir, env=env)
    if rc != 0:
        raise RuntimeError('case.setup for %s exited %d' % (case_name, rc))


def run_standalone(case_dir, backend, device='cpu', env=None, gpu_index=None):
    """`python3 -m eqdyna <case_dir> --backend <numpy|jax>` -- the
    exact command a user would type, with the backend ALWAYS named.

    It used to be invoked with no --backend at all, which silently meant the
    jax default: ~1700 lines of NumPy port (port.py, port_rsf.py, port_tp.py)
    had never been exercised by any gate while the tier reported success.
    """
    solver = {'python-numpy': 'numpy', 'python-jax': 'jax'}[backend]
    env = dict(env or base_env())
    env['PYTHONPATH'] = os.path.join(REPO_ROOT, 'src', 'python')
    env['PYTHONUNBUFFERED'] = '1'
    # JAX_PLATFORMS pins the device: a run labelled jax-on-cpu that silently
    # landed on a contended GPU is a different measurement under the same name.
    env['JAX_PLATFORMS'] = device
    if device != 'cpu' and backend == 'python-jax':
        # REFUSES rather than defaults. A GPU cell launched without a slot is
        # the item-58 failure: it shares a card with whatever else the sweep
        # started and dies with RESOURCE_EXHAUSTED, which then reads as a
        # parity or port regression in the results table.
        if gpu_index is None:
            raise RuntimeError(
                'run_standalone: --device %s needs a GPU slot; none was '
                'reserved for %s. Concurrent cells must each pin their own '
                'card (see GpuSlots) -- launching without one is how item 58 '
                'turned an allocator collision into a "failed cell".'
                % (device, os.path.basename(case_dir)))
        env = gpu_env(env, gpu_index)
        print('+ (%s) pinned to CUDA device %d, '
              'XLA_PYTHON_CLIENT_MEM_FRACTION=%s'
              % (os.path.basename(case_dir), gpu_index, GPU_MEM_FRACTION))
    rc = subprocess.call([sys.executable, '-u', '-m', 'eqdyna',
                          case_dir, '--backend', solver],
                         cwd=REPO_ROOT, env=env)
    if rc != 0:
        raise RuntimeError('eqdyna.eqdyna3d %s --backend %s exited %d'
                           % (os.path.basename(case_dir), solver, rc))
    frt = os.path.join(case_dir, 'frt.txt0')
    if not os.path.isfile(frt):
        raise RuntimeError('%s was not written by eqdyna.eqdyna3d' % frt)
    return frt


def run_python_jax_mpi(case_name, case_dir, env):
    """`mpirun -np <ranks> python3 -m eqdyna <case_dir> --backend jax --mpi`
    -- the real-MPI execution mode of the python-jax backend (matrix.py's
    `python-jax-mpi` column). `ranks` is matrix.PY_MPI_RANKS[case_name]; the
    case must have opted in or this is not called (see run_cell).

    EQDYNA_MPI_SYNC=halo is PINNED here, not inherited from the calling
    environment and not left to the solver's own default: `allreduce` exists
    in the solver but is measurement-only (MPI does not promise a fixed
    reduction order), and a gate that silently ran whichever mode the shell
    happened to export would not be testing what it claims to.

    This is the vacuous-gate guard for this cell: a launch that silently
    started fewer workers than `ranks` still writes SOME frt.txt<rank> files
    and can still compare green against the reference if the surviving
    rank(s) happen to own every fault node once -- exactly the failure this
    cell exists to catch. Two checks, both against matrix.py DATA rather than
    an assumption baked into this function:
      1. the number of frt.txt<rank> files actually written must equal
         matrix.PY_MPI_EXPECTED_FRT_FILES[(case, ranks)] -- NOT `== ranks`,
         because a rank whose element slab never touches the fault legitimately
         writes none (see that dict's comment in matrix.py).
      2. the PRE-DEDUP total row count across those files -- before
         frt_canonical's dedup-by-coordinate collapses any double-ownership --
         must equal the committed reference's own (already-deduped) row
         count. Canonicalisation's dedupe assumes every physically-shared node
         is written once per owning rank and agrees byte-for-byte across
         owners (true for Fortran's face-shared boundary nodes); this port's
         partition instead assigns each fault node to exactly one rank, so a
         correct run's pre-dedup total must land exactly ON the reference
         count, not merely at or above it. A node owned twice (with identical
         values, so `align()`'s own duplicate-value check would not fire)
         would silently inflate this total past the reference count and be
         deduped away unnoticed by everything downstream -- this check is the
         only place that number is ever looked at.
    """
    ranks = matrix.PY_MPI_RANKS[case_name]
    env = dict(env)
    env['PYTHONPATH'] = os.path.join(REPO_ROOT, 'src', 'python')
    env['PYTHONUNBUFFERED'] = '1'
    env['EQDYNA_MPI_SYNC'] = 'halo'  # pinned -- see docstring; never inherited
    cmd = [MPIRUN, '-np', str(ranks), sys.executable, '-u', '-m', 'eqdyna',
           case_dir, '--backend', 'jax', '--mpi']
    rc = _run(cmd, REPO_ROOT, env)
    if rc != 0:
        raise RuntimeError('%s exited %d' % (' '.join(cmd), rc))

    frt_files = frt_canonical.frt_rank_files(case_dir)
    expected_files = matrix.PY_MPI_EXPECTED_FRT_FILES.get((case_name, ranks))
    if expected_files is None:
        raise RuntimeError(
            'no matrix.PY_MPI_EXPECTED_FRT_FILES entry for (%r, %d) -- this '
            'cell cannot be gated without a measured expected file count '
            '(rule 2: "could not check" must not read as "passed")'
            % (case_name, ranks))
    if len(frt_files) != expected_files:
        raise RuntimeError(
            'python-jax-mpi %s at %d ranks: observed %d frt.txt<rank> '
            'file(s) (%s), expected %d (matrix.PY_MPI_EXPECTED_FRT_FILES). A '
            'launch that started fewer real workers than %d ranks must not '
            'be able to compare green.'
            % (case_name, ranks, len(frt_files),
               ', '.join(os.path.basename(p) for p in frt_files),
               expected_files, ranks))

    pre_dedup_rows = frt_canonical.load_frt(frt_files).shape[0]
    ref_rows = compare.load_reference(case_name).shape[0]
    if pre_dedup_rows != ref_rows:
        raise RuntimeError(
            'python-jax-mpi %s at %d ranks: %d fault-node row(s) written '
            'across %d rank file(s) BEFORE dedup, expected exactly %d '
            '(the reference row count) -- this partition assigns every '
            'fault node to exactly one rank, so any node owned twice (even '
            'with identical values, which dedup would then silently '
            'collapse) must fail here.'
            % (case_name, ranks, pre_dedup_rows, len(frt_files), ref_rows))
    return case_dir


def run_fortran(case_name, case_dir, eqdyna_cmd, env):
    """create.newcase -> the GATE_TERM_S override -> case.setup -> mpirun ->
    plotRuptureDynamics.

    The term override (apply_term_override) is a direct python call sitting
    between the first and second subprocess steps below, not a third build
    path: it edits the SAME user_defined_params.py create.newcase just copied
    from case_input/<case_name>/, which case.setup then reads -- the identical
    mechanism make_serial_case uses for the python backends' term override.

    plotRuptureDynamics is what writes fault.dyna.r.nc, the second artifact
    this backend is compared on (matrix.ARTIFACTS); dropping it once already
    turned that comparison into a comparison of a stale file.
    """
    test_dir = os.path.dirname(case_dir)
    rc = _run(['create.newcase', os.path.basename(case_dir), case_name], test_dir, env)
    if rc != 0:
        raise RuntimeError('create.newcase %s exited %d' % (case_name, rc))
    apply_term_override(case_name, case_dir)
    steps = [
        (['./case.setup'], case_dir),
        ([MPIRUN, '-np', str(matrix.FORTRAN_RANKS[case_name]), eqdyna_cmd], case_dir),
        ([sys.executable, 'plotRuptureDynamics'], case_dir),
    ]
    for cmd, cwd in steps:
        rc = _run(cmd, cwd, env)
        if rc != 0:
            raise RuntimeError('`%s` exited %d' % (' '.join(cmd), rc))


def run_cell(case, backend, test_dir, eqdyna_cmd, env, device, gpu_slots=None):
    """Run one cell and return its run directory. Raises on any failure."""
    if backend == 'fortran':
        case_dir = os.path.join(test_dir, case)
        run_fortran(case, case_dir, eqdyna_cmd, env)
        return case_dir
    case_dir = os.path.join(test_dir, '%s.%s' % (case, backend))
    # make_serial_case is UNCHANGED for python-jax-mpi: the case setup (one
    # serial-decomposition set of bFile/netCDF inputs) is identical to the
    # other python backends. Only the launch differs -- driver.run_mpi does
    # its OWN Fortran-style domain decomposition of that same serial case
    # across MPI ranks; par.nx/ny/nz above is a different, unrelated
    # decomposition (the Fortran binary's, which never runs here).
    make_serial_case(case, case_dir, env)
    if backend == 'python-jax-mpi':
        run_python_jax_mpi(case, case_dir, env)
    elif device != 'cpu' and backend == 'python-jax':
        # Hold a device for exactly this cell's run, then give it back. The
        # cell blocks here rather than sharing a card -- item 58.
        if gpu_slots is None:
            raise RuntimeError(
                'run_cell: --device %s selected but no GpuSlots pool was '
                'passed; refusing to launch %s x %s onto an unreserved card'
                % (device, case, backend))
        index = gpu_slots.acquire()
        try:
            run_standalone(case_dir, backend, device=device, env=env,
                           gpu_index=index)
        finally:
            gpu_slots.release(index)
    else:
        run_standalone(case_dir, backend, device=device, env=env)
    return case_dir


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------
def select(args):
    """(runnable, declared_unsupported, label, explicit, release_only) for
    this invocation. release_only is non-empty only for the DEFAULT everyday
    selection (no --cases/--backends, no --ci, no --release) -- matrix.py's
    RELEASE_ONLY cells (2026-09-23 owner decision): SUPPORTED cells held out
    of the everyday sweep for wall-clock cost, restored by --release (the
    release tier -- everyday cells + RELEASE_ONLY, at the SAME GATE_TERM_S;
    this is a cell-selection split, never a term split) and by any EXPLICIT
    --cases/--backends ask, which names exactly what it wants and is
    answered exactly, not cost-filtered."""
    if args.ci:
        # --backends/--cases, WHEN COMBINED WITH --ci, filter matrix.CI_CELLS
        # itself rather than switching to matrix.cells() -- this is what lets
        # a CI matrix job ask for "the fortran slice of CI_CELLS" or "the
        # python slice of CI_CELLS" while matrix.CI_CELLS stays the ONE
        # declared, memory-measured source of truth (rule 6). The workflow
        # never spells out a cell list of its own; it only names an axis
        # subset, so widening CI_CELLS in matrix.py widens every job that
        # asks for 'all' on that axis with no workflow edit needed.
        wanted = list(matrix.CI_CELLS)
        if args.cases:
            want_cases = set(args.cases.split(','))
            wanted = [c for c in wanted if c[0] in want_cases]
        if args.backends:
            want_backends = set(args.backends.split(','))
            wanted = [c for c in wanted if c[1] in want_backends]
        runnable, unsupported = [], []
        for cell in wanted:
            if matrix.is_supported(*cell):
                runnable.append(cell)
            else:
                unsupported.append(cell + (matrix.unsupported_reason(*cell),))
        filt = ''
        if args.cases or args.backends:
            filt = (', filtered to cases=%s backends=%s'
                    % (args.cases or 'all', args.backends or 'all'))
        return (runnable, unsupported,
                'CI (declared cell list, chosen against a measured %.0f GB '
                'runner -- matrix.CI_CELLS%s)' % (matrix.CI_RUNNER_RAM_GB, filt),
                True, [])
    cases = args.cases.split(',') if args.cases else None
    backends = args.backends.split(',') if args.backends else None
    explicit = bool(cases or backends)
    if explicit:
        # A named ask is answered exactly, not cost-filtered: RELEASE_ONLY is
        # a default-selection policy, not a per-cell refusal (matrix.py's
        # cells() already includes these cells; that is unchanged here).
        runnable, unsupported = matrix.cells(cases, backends)
        release_only = []
        label = 'explicit: cases=%s backends=%s' % (args.cases or 'all',
                                                      args.backends or 'all')
    elif args.release:
        runnable, unsupported = matrix.cells(cases, backends)
        release_only = []
        label = 'default: every cell of the table (release, --release)'
    else:
        runnable, unsupported, release_only = matrix.everyday_cells(cases, backends)
        label = 'default: every cell of the table minus matrix.RELEASE_ONLY (everyday)'
    return runnable, unsupported, label, explicit, release_only


def memory_note(runnable):
    lines = ['measured peak RSS per cell (rule 6 -- the number travels with '
             'the decision):']
    for cell in runnable:
        rss = matrix.MEASURED_PEAK_RSS_GB.get(cell)
        lines.append('  %-16s %-13s %s'
                     % (cell[0], cell[1],
                        ('%.2f GB' % rss) if rss is not None
                        else 'not measured'))
    return lines


def _perf_meta(results, label, device, budget, sha, tree_dirty):
    """The snapshot dict for this sweep's per-cell wall clocks (owner policy
    2026-09-22: every gate run is a free perf data point). Platform per cell
    is what the launch PINS, never what was merely requested:
      - fortran: CPU-only solver, no GPU path exists.
      - python-numpy: host arrays only.
      - python-jax with --device cpu (the default): run_standalone exports
        JAX_PLATFORMS=cpu, which jax cannot override onto a GPU.
      - python-jax with --device cuda: the sweep records NO device evidence,
        so the honest platform is 'unknown' -- requested is not measured
        (the 812.78 ms/step incident is exactly a requested label recorded
        as a measurement).
      - python-jax-mpi: eqdyna3d's MPI branch forces the cpu platform when
        --device is left at 'auto' (run_e2e passes no --device) -- the very
        override that produced that incident is, here, the pin.

    `sha`/`tree_dirty` are the caller's values, captured once at sweep START
    (see `main`) -- never recomputed here at sweep END, which is the exact
    rule-24 shape of bug this field exists to avoid on its own axis."""
    cells = []
    for case, backend, ok, dt, _lines in results:
        if backend == 'fortran':
            ranks = matrix.FORTRAN_RANKS[case]
            plat, ev = 'cpu', 'fortran solver: CPU-only, src/fortran has no GPU path'
        elif backend == 'python-jax-mpi':
            ranks = matrix.PY_MPI_RANKS[case]
            plat, ev = 'cpu', ('eqdyna3d --mpi with --device left at "auto" '
                               'forces JAX_PLATFORMS=cpu (run_e2e passes no '
                               '--device)')
        elif backend == 'python-numpy':
            ranks = 1
            plat, ev = 'cpu', 'numpy backend: host arrays only, no GPU path'
        elif device == 'cpu':
            ranks = 1
            plat, ev = 'cpu', ('run_standalone pins JAX_PLATFORMS=cpu; jax '
                               'cannot land on a GPU under that pin')
        else:
            ranks = 1
            plat, ev = 'unknown', ('JAX_PLATFORMS=%s was requested but the '
                                   'sweep records no device evidence; '
                                   'requested is not measured' % device)
        cells.append(dict(case=case, backend=backend, ok=bool(ok),
                          seconds=dt, ranks=ranks, platform=plat,
                          platform_evidence=ev))
    return dict(tool='run_e2e', sha=sha, host=os.uname().nodename,
                date=time.strftime('%Y-%m-%d %H:%M'), label=label,
                device=device, jobs_budget=budget,
                tenancy_ceiling=ledger.TENANCY_REFERENCE_CEILING, cells=cells,
                tree_dirty=tree_dirty)


def _capture_perf(results, label, device, budget, sha, tree_dirty):
    """Append this sweep's per-cell wall clocks to the perf ledger. NEVER part
    of the verdict: any error, including import failure, degrades to a loud
    WARNING (ledger.capture_e2e_cells_or_warn does the same for errors past
    the import), and only cells that already PASSED produce rows -- a red
    cell can neither look green nor leave a timing behind."""
    try:
        sys.path.insert(0, os.path.join(TESTSYS, 'perf'))
        import ledger
        ledger.capture_e2e_cells_or_warn(
            _perf_meta(results, label, device, budget, sha, tree_dirty))
    except (Exception, SystemExit) as exc:          # noqa: BLE001
        print('WARNING: perf-ledger capture failed (%s: %s) -- the sweep '
              'verdict is unaffected.' % (type(exc).__name__, exc))


_MAX_DIFF_RE = re.compile(r'max\|diff\|=([0-9.eE+-]+)')


def capture_start_tree_state():
    """git status --porcelain captured at the very start of the sweep (main(),
    before Gate 0, before the build, before the test/ rotation) -- i.e. before
    this process writes anything at all.

    This is the ONLY correct place to answer "was the tree already dirty
    coming into this run". Captured at the END instead (the prior bug), the
    reading is contaminated by the sweep's OWN writes: docs/perf_ledger.jsonl
    (appended to every run), docs/perf_snapshots/e2e_cells_*.json (a new file
    every run) and docs/evidence/sweep-<sha>/ (a new directory every release
    run) are all tracked paths -- none is in .gitignore -- so `git status
    --porcelain` at the end of any real sweep is never empty. Proof:
    docs/evidence/sweep-679d8ad/summary.json, committed in f9c055e, reads
    tree_clean: False from a run whose `git status --porcelain` was 0 lines
    when it STARTED.

    Capturing here also makes moot any question of EXCLUDING the sweep's own
    output paths from the porcelain diff: snapshotted before those paths are
    touched, they cannot appear in this snapshot regardless of whether they
    are later written, so there is nothing to exclude.

    Returns (tree_clean, dirty_lines): dirty_lines is the raw porcelain
    output, one entry per path, empty if clean -- so a dirty-at-start run
    names its files instead of just reporting False."""
    status_r = subprocess.run(['git', '-C', REPO_ROOT, 'status', '--porcelain'],
                              capture_output=True, text=True)
    dirty_lines = [ln for ln in status_r.stdout.splitlines() if ln.strip()]
    return (len(dirty_lines) == 0), dirty_lines


def write_release_evidence(results, is_release, explicit, started_utc,
                           finished_utc, tree_clean, dirty_at_start):
    """docs/evidence/sweep-<shortsha>/summary.json -- written only for a
    default `--release` selection (no --cases/--backends), i.e. the release
    tier sweeping the full runnable matrix (everyday cells + RELEASE_ONLY).
    A filtered selection is not 'the full runnable matrix' and gets no
    evidence artifact under this name.

    `tree_clean`/`dirty_at_start` come from capture_start_tree_state(),
    called at the top of main() -- BEFORE this sweep wrote anything -- not
    recomputed here at the end (see that function's docstring for why the
    end-of-run reading was always False for a real sweep).

    A separate mission is writing the tag-time guard that READS this
    schema -- the field NAMES below are a contract: extend, never rename or
    remove (per the dispatch that added this function). `term` is one such
    field: as of the 2026-09-23 one-term change its value is the numeric
    GATE_TERM_S every cell in this sweep actually ran at, not a 'full'/'gate'
    mode string -- there is no second mode left to name. `dirty_at_start` is
    a new field (this fix): extending, not renaming, per that same contract.
    """
    if not is_release or explicit:
        return
    sha_r = subprocess.run(['git', '-C', REPO_ROOT, 'rev-parse', 'HEAD'],
                           capture_output=True, text=True)
    sha = sha_r.stdout.strip()
    if sha_r.returncode != 0 or len(sha) != 40:
        raise RuntimeError(
            'write_release_evidence: `git rev-parse HEAD` did not return a '
            'full 40-char sha (rc=%d, stdout=%r) -- refusing to write an '
            'evidence artifact with no sha' % (sha_r.returncode, sha))
    cells, n_success = [], 0
    for case, backend, ok, dt, lines in results:
        if ok:
            n_success += 1
        m = _MAX_DIFF_RE.search(lines[0]) if lines else None
        cells.append(dict(case=case, backend=backend,
                          verdict='SUCCESS' if ok else 'FAIL',
                          max_diff=float(m.group(1)) if m else None,
                          wall_s=dt))
    payload = dict(sha=sha, tree_clean=tree_clean, dirty_at_start=dirty_at_start,
                  term=matrix.GATE_TERM_S,
                  n_runnable=len(results), n_success=n_success, cells=cells,
                  started_utc=started_utc, finished_utc=finished_utc)
    out_dir = os.path.join(REPO_ROOT, 'docs', 'evidence', 'sweep-%s' % sha[:7])
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, 'summary.json')
    with open(path, 'w') as f:
        json.dump(payload, f, indent=2)
        f.write('\n')
    print('e2e: wrote release evidence %s (%d cells, %d success, sha=%s, '
          'tree_clean=%s)' % (path, len(cells), n_success, sha, tree_clean))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--cases', help='comma-separated subset of the case axis')
    ap.add_argument('--backends', help='comma-separated subset of the backend axis')
    ap.add_argument('--ci', action='store_true',
                    help='run matrix.CI_CELLS, the declared portability-smoke '
                         'cell list (2026-09-23: one case, all three backend '
                         'implementations, at matrix.GATE_TERM_S)')
    ap.add_argument('--release', action='store_true',
                    help='the RELEASE selection: every supported cell '
                         '(everyday cells + matrix.RELEASE_ONLY), used by '
                         '`run.py release`. Same GATE_TERM_S as every other '
                         'selection -- this flag only widens which CELLS run, '
                         'not the term they run at. Ignored (a named ask is '
                         'answered exactly) when --cases/--backends is given.')
    ap.add_argument('--jobs', type=int, default=None,
                    help='core budget for concurrent cells (default: TENANCY-'
                         'AWARE -- cores measured free right now, minus a '
                         'margin, never below the largest selected cell\'s '
                         'own cost; see default_jobs_budget). A fortran or '
                         'python-jax-mpi cell costs its rank count, a '
                         'python-jax cell its measured core count '
                         '(matrix.JAX_MEASURED_CORES, ceil\'d), a '
                         'python-numpy cell 1 (also measured). --jobs 1 runs '
                         'serially, which is what a 2-core CI runner should '
                         'use; the printed output is identical either way.')
    ap.add_argument('--device', default='cpu', choices=('cpu', 'cuda'),
                    help='JAX_PLATFORMS for the python-jax backend (default cpu)')
    args = ap.parse_args(argv)

    # Captured HERE -- first thing in main(), ahead of Gate 0, the build and
    # the test/ rotation -- so the reading is of the tree as the CALLER left
    # it, not as this sweep leaves it. See capture_start_tree_state()'s
    # docstring for why end-of-run capture (the prior bug) always read dirty.
    tree_clean_at_start, dirty_at_start = capture_start_tree_state()
    if not tree_clean_at_start:
        print('e2e: tree was NOT clean at sweep start (%d path(s)) -- a '
              '--release run from here will record tree_clean=False:'
              % len(dirty_at_start))
        for ln in dirty_at_start:
            print('  %s' % ln)

    runnable, unsupported, label, explicit, release_only = select(args)

    print('\n==== e2e sweep: coverage ====')
    print('term     : %gs (matrix.GATE_TERM_S, applied to every selected '
          'case regardless of its own committed par.term -- there is only '
          'one term)' % matrix.GATE_TERM_S)
    for line in matrix.coverage_report(runnable, unsupported, label, release_only):
        print(line)
    for line in memory_note(runnable):
        print(line)
    if args.device != 'cpu':
        print('device: JAX_PLATFORMS=%s for python-jax cells' % args.device)

    if explicit and unsupported:
        print('\ne2e: FAIL - the selection names %d cell(s) this table '
              'declares unsupported. A cell you asked for and could not be '
              'run is a failure, not a skip (rule 2):' % len(unsupported))
        for c, b, reason in unsupported:
            print('  - %s x %s: %s' % (c, b, reason))
        return 1
    if not runnable:
        print('\ne2e: FAIL - the selection runs zero cells; a sweep that '
              'compares nothing must never exit green')
        return 1

    # Gate 0 - ONE live invocation per run tree (rule 21a, pathway item 70).
    #
    # Taken HERE, ahead of Gate 1 and ahead of the build, rather than
    # immediately before the rotation further down: an invocation that is going
    # to be refused should be refused before it spends a Fortran build, and the
    # lock then covers bin/ for the same price. It is released by the kernel
    # when this process exits, however it exits (testsys/runlock.py explains
    # why that is flock and not an O_EXCL lockfile).
    try:
        lock = runlock.acquire(REPO_ROOT, 'test')
    except runlock.RunTreeLocked as exc:
        print('\ne2e: FAIL - %s' % exc)
        return 1
    print('e2e: holding %s (pid %d) - the test/ rotation below is serialised '
          'against every other invocation in this checkout' % (lock.path, lock.pid))

    # Gate 1 - cheap check before the expensive runs (rule 9).
    guard = os.path.join(REPO_ROOT, 'testsys', 'regression', 'test_create_newcase.py')
    if subprocess.call([sys.executable, guard]) != 0:
        print('e2e: FAIL - create.newcase guard failed; aborting before expensive runs')
        return 1

    # Gate 2 - fresh build, if and only if a fortran cell is selected.
    eqdyna_cmd = None
    if any(b == 'fortran' for _, b in runnable):
        if BIN_OVERRIDE:
            eqdyna_cmd = os.path.abspath(BIN_OVERRIDE)
            if not os.path.exists(eqdyna_cmd):
                print('e2e: FAIL - EQDYNA_E2E_BIN=%s does not exist' % BIN_OVERRIDE)
                return 1
            print('e2e: EQDYNA_E2E_BIN set - NOT rebuilding; this run says '
                  'nothing about the build-and-install path (rule 16). Using %s'
                  % eqdyna_cmd)
        else:
            bin_exe = os.path.join(REPO_ROOT, 'bin', 'eqdyna')
            if os.path.exists(bin_exe):
                os.remove(bin_exe)
            rc = subprocess.call(['./install-eqdyna.sh', '-m', MACHINE], cwd=REPO_ROOT)
            if rc != 0:
                print('e2e: FAIL - ./install-eqdyna.sh -m %s exited %d' % (MACHINE, rc))
                return 1
            if not os.path.exists(bin_exe):
                print('e2e: FAIL - build finished but bin/eqdyna does not exist')
                return 1
            eqdyna_cmd = 'eqdyna'
    else:
        print('e2e: no fortran cell in this selection - no Fortran build needed')

    # Rule 8 - preserve, never delete, the previous run's evidence.
    test_dir = os.path.join(REPO_ROOT, 'test')
    prev_dir = os.path.join(REPO_ROOT, 'test.prev')
    if os.path.isdir(test_dir):
        if os.path.isdir(prev_dir):
            shutil.rmtree(prev_dir)
        shutil.move(test_dir, prev_dir)
        print('e2e: preserved previous run as %s' % prev_dir)
    os.makedirs(test_dir)

    env = base_env()
    start = time.time()
    started_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

    # sha + tree_dirty for the profile-guard/collection calls below, captured
    # HERE at sweep START, before any cell runs -- the same rule-24 lesson
    # already applied to write_release_evidence's tree_clean: computed at the
    # END would read as clean if uncommitted work from DURING the sweep gets
    # committed before the sweep finishes, silently mislabelling a row that
    # was actually produced by a dirty tree. Computed once, not per cell: it
    # does not change mid-sweep on this axis either.
    sys.path.insert(0, os.path.join(TESTSYS, 'perf'))
    import ledger
    _sha_r = subprocess.run(['git', '-C', REPO_ROOT, 'rev-parse', '--short',
                             'HEAD'], capture_output=True, text=True)
    if _sha_r.returncode != 0 or not _sha_r.stdout.strip():
        raise RuntimeError('git rev-parse --short HEAD failed (rc=%d): %s'
                           % (_sha_r.returncode, _sha_r.stderr.strip()))
    sweep_sha = _sha_r.stdout.strip()
    sweep_tree_dirty = ledger.tree_dirty()

    # Cells are independent: each gets its own directory ('<case>.<backend>'),
    # reads the same read-only reference tree, and shares nothing else. So the
    # sweep runs them concurrently rather than one at a time.
    #
    # WHY THIS MATTERS: measured serially, the 20-cell sweep took 3496 s on a
    # 64-core box, and 67% of that was the python-numpy column alone
    # (2267 s over 6 cells; test.tpv29 964 s, test.drv.a6 553 s). Meanwhile at
    # most 4 cores were busy. Wall time is now bounded by the SLOWEST CELL, not
    # the sum.
    #
    # The budget is CORES, not memory: all 20 cells together peak at ~23.5 GB
    # (matrix.MEASURED_PEAK_RSS_GB) against 1 TB here. Per-cell cost is
    # cell_cost() above (module-level, MEASURED per backend -- see its own
    # docstring): a fortran or python-jax-mpi cell costs its real rank count,
    # a python-jax cell costs matrix.JAX_MEASURED_CORES (ceil'd, not 1 -- see
    # the 2026-09-23 speed-campaign citation there), a python-numpy cell costs
    # 1 (also measured, not assumed -- eqdyna3d._narrow_numpy_affinity pins it
    # to one cpu). --jobs sets the budget explicitly; left unset, the budget
    # is TENANCY-AWARE (default_jobs_budget: cores measured free right now,
    # minus a margin, never below the largest selected cell's own cost).
    # --jobs 1 restores the serial order exactly, which is what CI uses when
    # its runner has 2.
    #
    # Output is COLLECTED per cell and printed when that cell finishes, never
    # streamed, so concurrent cells cannot interleave their lines into an
    # unreadable log. The results table is re-sorted into table order
    # afterwards, so a parallel run and a serial run print identically.
    table_cells = [(c, b) for c in matrix.CASES for b in matrix.BACKENDS
                  if (c, b) in runnable]
    order = {cb: i for i, cb in enumerate(table_cells)}

    # Longest-first start order (item 2, 2026-09-23): reorders SUBMISSION
    # only -- `order` above (used to re-sort the printed results table) is
    # fixed to the table's own case x backend order regardless, so a
    # scheduled run and a --jobs 1 serial run still PRINT identically.
    ledger_costs = load_ledger_wall_costs()
    cells = schedule_order(table_cells, ledger_costs)
    print('e2e: start order (longest measured wall-clock first, '
          'docs/perf_ledger.jsonl latest cell-wall-clock row per cell; '
          'unmeasured cells scheduled first as the conservative choice):')
    for cb in cells:
        cost = ledger_costs.get(cb)
        print('  %-16s %-13s %s'
              % (cb[0], cb[1],
                 ('%.1fs (measured)' % cost) if cost is not None
                 else 'UNMEASURED -- scheduled first'))

    # One slot per visible card, created only when a GPU cell is actually
    # selected -- asking nvidia-smi on a CPU sweep would make a CPU-only box
    # fail for no reason.
    gpu_slots = None
    if args.device != 'cpu' and any(b == 'python-jax' for _, b in cells):
        gpu_slots = GpuSlots(visible_gpu_indices())
        print('e2e: %d CUDA device(s) %s -- at most %d concurrent python-jax '
              'cell(s), each pinned to its own card at '
              'XLA_PYTHON_CLIENT_MEM_FRACTION=%s (item 58)'
              % (len(gpu_slots.devices), gpu_slots.devices,
                 len(gpu_slots.devices), GPU_MEM_FRACTION))

    def run_one(cb):
        case, backend = cb
        t0 = time.time()
        try:
            case_dir = run_cell(case, backend, test_dir, eqdyna_cmd, env,
                                args.device, gpu_slots=gpu_slots)
            ok, lines = compare.compare_cell(case, backend, case_dir)
            if ok:
                # Mechanical profile guard + collection (owner mission,
                # 2026-09-23): every cell this sweep PASSES must also carry a
                # present, schema-valid, sum-checked profile.rank<r>.json per
                # rank (profile_schema.validate_run_dir, called inside
                # capture_run), appended to docs/run_profiles.jsonl.
                # capture_run RAISES on any problem -- missing/invalid
                # profile files included -- and there is deliberately no
                # `_or_warn` variant here (unlike the perf-ledger capture
                # below): a profile defect must be able to turn a physics-
                # passing cell into a failed one, not vanish into a warning.
                profile_record.capture_run(
                    case_dir, case=case, backend=backend,
                    ranks=cell_cost(case, backend),
                    term=str(matrix.GATE_TERM_S), sha=sweep_sha)
        except Exception as exc:                # noqa: BLE001 - reported, not swallowed
            ok, lines = False, ['%s: %s' % (type(exc).__name__, exc)]
        return (case, backend, ok, time.time() - t0, lines)

    budget = default_jobs_budget(cells, args)
    results = []
    if budget <= 1:
        for cb in cells:
            print('\n-- cell: %s x %s --' % cb)
            r = run_one(cb)
            results.append(r)
            for line in r[4]:
                print('   ' + line)
            print('%s %s x %s (%.1fs)'
                  % ('SUCCESS' if r[2] else 'FAIL', r[0], r[1], r[3]))
    else:
        import concurrent.futures as _cf
        import threading
        lock = threading.Lock()

        # ALL-OR-NOTHING core reservation. A cell costs more than one core
        # (a fortran cell costs its rank count), and the obvious spelling --
        # `for _ in range(cost): sem.acquire()` on a Semaphore -- DEADLOCKS,
        # because it takes units one at a time: with budget 6, two cells each
        # needing 4 can end up holding 3 apiece and both wait forever for a
        # fourth the other is holding. That is not hypothetical; it hung this
        # sweep for 80 minutes with 14 cells left, parent alive at 0.1% CPU
        # and no children. A Condition lets a thread take its whole cost
        # atomically or not at all, which cannot deadlock.
        cond = threading.Condition()
        free_cores = budget

        def reserve(cost):
            nonlocal free_cores
            with cond:
                while free_cores < cost:
                    cond.wait()
                free_cores -= cost

        def release(cost):
            nonlocal free_cores
            with cond:
                free_cores += cost
                cond.notify_all()

        def guarded(cb):
            # min(): a cell that costs more than the whole budget would
            # otherwise wait forever for cores that will never exist.
            cost = min(cell_cost(*cb), budget)
            reserve(cost)
            try:
                return run_one(cb)
            finally:
                release(cost)

        print('e2e: running %d cells concurrently, core budget %d '
              '(serial would be the sum of all cell times)'
              % (len(cells), budget))
        with _cf.ThreadPoolExecutor(max_workers=len(cells)) as pool:
            futs = {pool.submit(guarded, cb): cb for cb in cells}
            for fut in _cf.as_completed(futs):
                r = fut.result()
                with lock:
                    print('\n-- cell: %s x %s (%.1fs) --' % (r[0], r[1], r[3]))
                    for line in r[4]:
                        print('   ' + line)
                    print('%s %s x %s'
                          % ('SUCCESS' if r[2] else 'FAIL', r[0], r[1]))
                    results.append(r)
        results.sort(key=lambda r: order[(r[0], r[1])])
    elapsed = time.time() - start
    finished_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

    print('\n==== e2e sweep: results ====')
    print('%-16s %-13s %-8s %8s  %s' % ('case', 'backend', 'verdict', 'seconds',
                                        'measure'))
    for case, backend, ok, dt, lines in results:
        print('%-16s %-13s %-8s %8.1f  %s'
              % (case, backend, 'SUCCESS' if ok else 'FAIL', dt,
                 lines[0] if lines else ''))
    failed = [(c, b) for c, b, ok, _, _ in results if not ok]

    print('\n==== e2e sweep: SUMMARY ====')
    print('selection : %s' % label)
    print('ran       : %d of %d cells in the %d case x %d backend table '
          '(%d passed, %d failed)'
          % (len(results), len(matrix.CASES) * len(matrix.BACKENDS),
             len(matrix.CASES), len(matrix.BACKENDS),
             len(results) - len(failed), len(failed)))
    print('not gated : %d declared-unsupported cell(s): %s'
          % (len(unsupported),
             ', '.join('%s x %s' % (c, b) for c, b, _ in unsupported) or 'none'))
    print('release-only (not run in this everyday sweep): %d cell(s): %s'
          % (len(release_only),
             ', '.join('%s x %s' % (c, b) for c, b, _ in release_only) or 'none'))
    print('wall clock: %.1fs' % elapsed)
    # Every sweep is a free timing data point (owner policy 2026-09-22).
    # Placed BEFORE the verdict returns below but able to affect none of them:
    # _capture_perf swallows everything into a WARNING.
    _capture_perf(results, label, args.device, budget, sweep_sha,
                 sweep_tree_dirty)
    # Release-sweep evidence: only for a default (--release, no
    # --cases/--backends) run over the full runnable matrix -- see the
    # function's own docstring. Written before the pass/fail return below so
    # a failed release sweep still leaves its evidence on disk.
    write_release_evidence(results, args.release, explicit, started_utc,
                           finished_utc, tree_clean_at_start, dirty_at_start)
    if len(results) != len(runnable):
        print('e2e: FAIL - %d cell(s) were selected but %d produced a verdict; '
              'a cell that produced no verdict is a failure'
              % (len(runnable), len(results)))
        return 1
    if failed:
        print('e2e: FAIL - %d cell(s) failed: %s'
              % (len(failed), ', '.join('%s x %s' % cb for cb in failed)))
        return 1
    print('e2e: SUCCESS - %d/%d selected cells matched test.reference.results/ '
          '(this states its own scope; see the coverage block above)'
          % (len(results), len(runnable)))
    return 0


if __name__ == '__main__':
    if '--full' in sys.argv:
        # The full tier (SCEC spec resolution/duration, 16 ranks, report-only)
        # is a different gate shape, not a cell of this sweep -- see
        # run_e2e_full.py's module docstring.
        sys.path.insert(0, E2E_DIR)
        from run_e2e_full import main as main_full
        sys.exit(main_full())
    sys.exit(main())
