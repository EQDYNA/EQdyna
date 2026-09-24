#! /usr/bin/env python3
"""
Strong-scaling measurement (report-only; never a red/green gate).
pathway_forward item 33. Companion to `run_numa_scaling.py` (the JAX-only
locality tool fixed earlier this session, commit 03ba055/0b4684c) -- this
tool covers BOTH python backends (numpy and jax) and Fortran MPI, over the
full 1..32-core range, and follows the same pattern for the reasons below.

TWO DEFECTS FIXED HERE (found by reading this file, not assumed from a
prior report):

1. `PY_THREADS` stopped at 8 while `FORTRAN_RANKS` went to 32 -- the python
   side could not structurally answer whether either backend reaches 32
   cores. Now `PY_THREADS == FORTRAN_RANKS == [1,2,4,8,16,32]`, for BOTH
   numpy and jax. Never 64: a single process has no reason to span both
   sockets (owner constraint, same one `run_numa_scaling.py` enforces).

2. Every pin used bare `taskset -c <cores>` with `cores = range(n)`. Two
   bugs in one: (a) on this box (8 NUMA nodes x 8 cores, confirmed via
   `numactl --hardware`, not assumed) `range(16)` spans 2 nodes and
   `range(32)` spans 4, so a knee at 8->16 would be a NUMA-crossing artifact
   misread as a scaling limit; (b) `taskset` binds CPU only, not memory --
   first-touch allocation can land on any node, so a "pinned" run can still
   fault in memory remotely. Both python and Fortran paths now use
   `numactl --physcpubind=<cpus> --membind=<nodes>`, and every data point
   records exactly which cpus/nodes it used.

THE PIN IS ALSO THE FIX, NOT JUST THE FAIRNESS. Measured directly, on THIS
box, before writing this tool: OpenBLAS (scipy-openblas 0.3.29, this
project's numpy build) reads OPENBLAS_NUM_THREADS and probes its affinity
mask ONLY at library load (`import numpy`, transitively). Setting the env
var AFTER `import numpy` has ZERO effect (0.40s vs 0.36s default on a
4000x4000 matmul -- indistinguishable); numactl-restricting the process
to N cpus BEFORE python starts (which is what `taskset`/`numactl` always
do, since they exec a fresh process image) makes OpenBLAS auto-detect
exactly N threads with NO explicit OPENBLAS_NUM_THREADS needed at all
(0.80s at 4 cpus pinned+auto vs 0.81s at 4 cpus pinned+explicit=4 --
equal). So a bare `taskset` pin was already sizing BLAS's thread pool
correctly; explicit BLAS-thread env vars are not the fix here and this
tool does not set them. This REFRAMES, not confirms, pathway item 40's
"numpy anti-scales" hypothesis: that comparison was pinned-to-1-core vs
completely UNPINNED (sched_getaffinity sees all 64 cpus, so OpenBLAS
spawns up to 64 threads that then contend with this box's ~9 other
tenants for far fewer than 64 real cores -- a placement-lottery artifact
of an unfair baseline, not necessarily a property of the solver at a
given PINNED core count). There was, before this tool, no actual
pinned-N-core measurement for numpy at N>1 to check that against. This
tool's job is to produce that measurement honestly; see the session
report for what it found and whether a port-side fix followed from it.

Same reasoning for jax: `run_numa_scaling.py` already established (and
this tool follows) that XLA also auto-sizes its intra-op thread pool from
the process's cpu affinity -- do NOT force XLA_FLAGS/OMP_NUM_THREADS here;
let the pin speak for itself, exactly like the numa tool does.

TWO PLACEMENT POLICIES, so a knee can be attributed (same idea as
`run_numa_scaling.py`'s within-node/spread contrast, generalised to every
core count instead of one):
  COMPACT -- fill one node's cpus, then the next, etc. Minimum NUMA spread
             for a given core count.
  SPREAD  -- round-robin across nodes (one per node, then a second pass,
             ...). Maximum NUMA spread for the SAME core count. If compact
             scales and spread does not, the knee is locality; if both
             flatten together, it is the solver running out of parallel work.

NODE SELECTION IS NOW FREE-FIRST, NOT NODE0-FIRST (F4, 2026-09-17 remeasurement
of item 33). Both policies used to iterate `sorted(nodes)`, i.e. always start
at node 0. On this box node 0 (and neighbours) intermittently carry foreign
tenants -- a same-session sweep saw cpu 0/1/10/16-19/30/40-41 busy across four
attempts and skipped 12 of 24 configs, even though 60 of 64 cores were free
the whole time on nodes elsewhere. The tool was asking for the wrong 32, not
finding that the box lacked room. `free_node_map()` probes every cpu's busy
fraction ONCE per configuration (reusing `run_numa_scaling.cpu_busy_fractions`,
not reimplementing it) and restricts BOTH policies to nodes that are ENTIRELY
under the busy ceiling right now; `compact_cpus`/`spread_cpus` then run
unchanged over that filtered, still-node0-first-if-node0-is-free set. If the
free set cannot supply the requested core count, the configuration is SKIPPED
and recorded exactly as a busy-cpu skip is -- this widens WHERE the tool looks
for room, it does not loosen what it demands of that room (the ceiling itself,
`--busy-ceiling`, is untouched). When node 0 already has room this reduces to
the old behaviour exactly, since node 0 sorts first among the free nodes too.

PER-CPU BUSY CHECK BEFORE EVERY POINT, same discipline and same function
(`cpu_busy_fractions`/`require_idle`, imported from `run_numa_scaling.py`
rather than re-implemented -- one bug fixed once) as the locality tool:
a whole-box load average is unmeasurable-by-construction on a box with
long-running unrelated tenants (this one has ~9, load 6-9 all session),
so the ceiling is evaluated per-cpu, per-point, on exactly the cpus a
configuration is about to use. A busy point is SKIPPED and recorded, not
silently included and not silently dropped.

PER-STEP BY DIFFERENCE for the two python backends (`per_step_py`, same
technique as `run_numa_scaling.py`'s `per_step`): run n_lo and n_hi steps
in separate fresh processes (jax caches compiled code in-process, so a
second call in the same interpreter would pay no compile and cancel the
wrong term) and divide the wall-time difference by the step-count
difference, so fixed cost (interpreter start, case load, XLA compile)
cancels exactly. FORTRAN NOW USES THE SAME TECHNIQUE (`per_step_fortran`,
2026-09-19). It previously took one wall-clock run and reported
`wall / n_hi`, justified here as "no JIT to amortise, and its fixed cost is
genuinely small next to a multi-second MPI solve" -- true of a long
production run, false of the 20-60 step runs this tool actually does, where
mesh generation, input read, MPI init and netCDF open are a large fraction
of the wall time. Charging Fortran for its fixed cost while subtracting
jax's is a one-way bias in jax's favour and made the headline verdict a
property of the metric. Both engines now difference two runs.

Usage:
    python3 testsys/perf/run_scaling.py [--case test.tpv104]
        [--fortran-ranks 1,2,4,8] [--py-threads 1,2,4,8,16,32]
        [--policies compact,spread] [--n-lo 20] [--n-hi 60]
        [--busy-ceiling 0.2] [--i-know-the-box-is-busy] [--repeats 1]
        [--skip-fortran] [--backends numpy,jax]
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
# The checkout THIS FILE lives in, derived from __file__ and never from
# $EQDYNAROOT (item 77). `ROOT` above may point at a different checkout
# entirely; the case directory rebuilt below is built from TESTSYS, so the
# lock that guards it must be rooted the same way or it would lock one tree
# and destroy another.
REPO_ROOT = os.path.dirname(os.path.dirname(TESTSYS))
PYTHON_PKG = os.path.join(ROOT, 'src', 'python')
OUT = os.path.join(TESTSYS, 'scaling_last.json')

sys.path.insert(0, TESTSYS)
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.dirname(TESTSYS))  # testsys/ itself, for profile_record
import run_numa_scaling as numa  # noqa: E402  (numa_topology, cpu_busy_fractions, require_idle)
import perflib                   # noqa: E402  (acquire_case_lock, rebuild_serial_case)
import profile_record            # noqa: E402  (append-only per-rank profile totals)

# What a second concurrent invocation costs, printed by the refusal (item 77).
LOCK_CONSEQUENCE = [
    'A second invocation in this checkout rmtrees and rebuilds that SAME case'
    ' directory',
    'while the first is TIMING out of it. The first does not crash: it reports'
    ' seconds,',
    'and the collision arrives as a scaling knee or a speedup. FIVE tools call'
    ' this one',
    'builder -- run_scaling.py, run_mpi_scaling.py, run_jaxmpi_ab.py,'
    ' run_setup_probe.py',
    'and run_shard_scaling.py -- so the two colliding runs need not even be the'
    ' same tool',
    '(rule 21a, pathway item 77).',
    '',
    'NOT waiting. NOT rebuilding anyway. NOT falling back to a second case'
    ' directory --',
    'each of those is the silent fallback rule 2 forbids. Run your perf tool in'
    ' its own',
    'git worktree (rule 21a), or wait for the holder above to finish.']

# The lock, once per process. build_py_case has five callers and nothing stops
# one process from calling it twice; flock is per open file description, so a
# second acquire in THIS process would refuse against our own pid -- a false
# collision. The memo is not a fallback: the lock is genuinely held.
_case_lock = None

CASE = 'test.tpv104'  # tpv8 runs out of parallel work early (prior session notes); use a
                       # case with real per-step work -- tpv104 (friclaw 4, ~970k elements).
DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1), 8: (2, 2, 2),
          16: (4, 2, 2), 32: (4, 4, 2)}
FORTRAN_RANKS = [1, 2, 4, 8, 16, 32]
PY_THREADS = [1, 2, 4, 8, 16, 32]  # F1: was [1,2,4,8]; now matches FORTRAN_RANKS, both backends


def sh(cmd, **kw):
    r = subprocess.run(cmd, shell=True, text=True, capture_output=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f'FAIL ({r.returncode}): {cmd}\n{r.stderr[-800:]}')
    return r


# --------------------------------------------------------------------------
# placement: compact / spread cpu lists from `numa_topology()`'s {node: [cpu,...]}
# --------------------------------------------------------------------------
def compact_cpus(nodes, k):
    cpus = []
    for n in sorted(nodes):
        for c in nodes[n]:
            cpus.append(c)
            if len(cpus) == k:
                return cpus
    raise ValueError('compact_cpus: only %d cpus available, need %d' % (len(cpus), k))


def spread_cpus(nodes, k):
    node_ids = sorted(nodes)
    used = {n: 0 for n in node_ids}
    cpus = []
    while len(cpus) < k:
        progressed = False
        for n in node_ids:
            if used[n] < len(nodes[n]):
                cpus.append(nodes[n][used[n]])
                used[n] += 1
                progressed = True
                if len(cpus) == k:
                    break
        if not progressed:
            raise ValueError('spread_cpus: only %d cpus available, need %d' % (len(cpus), k))
    return cpus


def cpu_to_node(node_map):
    """{cpu: node}, inverted from `numa_topology()`'s {node: [cpu, ...]}.

    ONE COPY of the inversion. `nodes_of` (below), `fortran_cmd`'s per-rank
    `--membind` list and `run_mpi_scaling.least_loaded_cpus`'s NUMA-compact
    tie-break each built this dict themselves; two of them disagreeing would
    bind a rank's MEMORY to a node its CPU is not on, which is silent -- the
    run still completes and reports seconds.
    """
    return {c: n for n, cs in node_map.items() for c in cs}


def nodes_of(node_map, cpus):
    cpu2node = cpu_to_node(node_map)
    return sorted({cpu2node[c] for c in cpus})


def numactl_prefix(cpus, node_map):
    nodes = nodes_of(node_map, cpus)
    return ['numactl', '--physcpubind=%s' % ','.join(str(c) for c in cpus),
            '--membind=%s' % ','.join(str(n) for n in nodes)]


def free_node_map(all_nodes, busy_ceiling, override):
    """{node: cpus} restricted to the individual cpus currently under
    `busy_ceiling`, so `compact_cpus`/`spread_cpus` build placements out of
    room that actually exists right now instead of always starting at node 0
    (see module docstring, F4). Probes fresh on every call -- occupancy moves
    over the course of a multi-minute sweep, so a configuration built early
    from a stale free-set could target a cpu that has since gone busy (the
    per-cpu `require_idle` check downstream still catches that case, but
    asking for the wrong cpu in the first place is the defect being fixed).

    PER-CPU, NOT WHOLE-NODE (2026-09-18 fix, item 33 thread-scaling
    investigation). The original version required EVERY cpu on a node to be
    idle before offering ANY of that node's cpus -- correct for a box where
    interference clusters by node, wrong for THIS box: ~20-28 foreign
    single-core jobs with no cpu affinity of their own, so the Linux
    scheduler smears them across all 8 nodes and no node is ever seen fully
    idle in a snapshot even though ~40+ of 64 cores are free at any instant
    (measured: `free_node_map` returned `{}` on 6 consecutive samples over
    12s with `nproc` idle cores in the 40s the whole time). Requiring
    whole-node freedom in that regime means the tool can never select ANY
    placement, at ANY core count including k=1, and every run silently
    starves on SKIPPED rather than measuring -- worse than the coarse
    placement F4 fixed. Cpu-level filtering keeps `compact_cpus` filling
    node-by-node in cpu-id order exactly as before (so it stays AS compact
    as the actually-free cpus allow) but no longer demands more idle room
    than a configuration actually needs.

    `override` mirrors `--i-know-the-box-is-busy`: if the busy ceiling itself
    is being bypassed, there is nothing to filter FOR, so this reverts to the
    full topology in node-number order -- identical to pre-fix behaviour.

    Raises SystemExit if per-cpu utilisation could not be read AT ALL (same
    hard-failure discipline as `cpu_busy_fractions`/`require_idle`: a check
    that cannot evaluate must fail, not silently call every cpu free). An
    individual cpu whose utilisation could not be read is dropped from the
    free set (not assumed free, not used to invalidate cpus that WERE read)."""
    if override:
        return dict(all_nodes)
    all_cpus = sorted(c for cs in all_nodes.values() for c in cs)
    busy = numa.cpu_busy_fractions(all_cpus)
    if not busy:
        raise SystemExit(
            'FAIL: could not read per-cpu utilisation for cpus %s from '
            '/proc/stat -- cannot tell which cpus are free (rule 2: a '
            'check that cannot evaluate must fail).' % all_cpus)
    free = {}
    for n, cpus in all_nodes.items():
        idle_cpus = [c for c in cpus if busy.get(c) is not None and busy[c] <= busy_ceiling]
        if idle_cpus:
            free[n] = idle_cpus
    return free


def select_cpus(free_nodes, k, policy):
    """cpus for a config of size k, built from ONLY currently-free NUMA nodes
    (point 2/3 of the fix). Returns None -- not an exception -- if the free
    set cannot supply k cpus under the given policy; the caller records this
    exactly like a busy-cpu skip. This widens where the tool looks for room;
    it never lowers what it demands of that room."""
    fn = compact_cpus if policy == 'compact' else spread_cpus
    try:
        return fn(free_nodes, k)
    except ValueError:
        return None


# --------------------------------------------------------------------------
# case construction (Fortran: rank-count decomposition; python: forced serial)
# --------------------------------------------------------------------------
def make_case(dst, nx, ny, nz, term=None):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    sh(f'create.newcase {dst} {CASE}')
    p = os.path.join(dst, 'user_defined_params.py')
    s = open(p).read()
    s = re.sub(r'par\.nx\s*,\s*par\.ny\s*,\s*par\.nz\s*=.*',
               f'par.nx, par.ny, par.nz = {nx}, {ny}, {nz}', s)
    s = re.sub(r'^par\.nx = .*', f'par.nx = {nx}', s, flags=re.M)
    s = re.sub(r'^par\.ny = .*', f'par.ny = {ny}', s, flags=re.M)
    s = re.sub(r'^par\.nz = .*', f'par.nz = {nz}', s, flags=re.M)
    if term is not None:
        s = re.sub(r'^par\.term\s*=.*', f'par.term = {term!r}', s, flags=re.M)
    open(p, 'w').write(s)
    sh('./case.setup', cwd=dst)


def read_term_dt(dst):
    """(term, dt) exactly as case.setup wrote them to bGlobal.txt -- lines
    15/16 (0-indexed) per scripts/case.setup's fixed write order (term right
    after the 'nx ny nz' line and its blank separator; dt right after term).
    Fortran computes nstep = idnint(term/dt) at runtime (readInputFiles.f90);
    reading these back (rather than recomputing par.dt ourselves) means we
    use the SAME numbers the binary will use, not a second copy that could
    drift from case.setup's write order."""
    lines = open(os.path.join(dst, 'bGlobal.txt')).read().splitlines()
    return float(lines[15]), float(lines[16])


# --------------------------------------------------------------------------
# Fortran: rank-count decomposition, one wall-clock run per (n, policy)
# --------------------------------------------------------------------------
def fortran_cmd(n, cpus, node_map, binary):
    """mpirun with each of the n ranks numactl-pinned to its own cpu/node.
    A bare `mpirun --bind-to core` cannot be trusted to respect an outer
    numactl restriction (OpenMPI's own hwloc-based binder can reissue
    sched_setaffinity to any online cpu, which is not blocked by a plain
    affinity mask the way a cgroup would block it) -- so each rank pins
    itself explicitly via OMPI_COMM_WORLD_LOCAL_RANK, the correct and
    portable way to combine OpenMPI with numactl.

    `--bind-to none` IS REQUIRED (2026-09-18 fix, item 33 investigation):
    without it, OpenMPI's OWN default binding (--bind-to core, applied
    automatically whenever np>1) restricts each rank's cpu affinity mask
    to ITS OWN choice of core BEFORE the numactl call below ever runs. Our
    numactl then tries to move the rank to `cpus[rank]`, which is very
    often a DIFFERENT cpu than the one OpenMPI already bound it to --
    libnuma refuses with 'cpu argument N is out of range' (relative to the
    already-narrowed mask, not the system total) and the whole mpirun job
    aborts. Measured: `mpirun -np 2` with no `--bind-to` flag fails this
    way on rank 1 every time; adding `--bind-to none` (verified directly,
    not assumed) lets the explicit per-rank numactl below be the ONLY
    binding authority, exactly as this function's docstring already
    intended. n=1 never showed this bug (no second rank has to move to a
    non-default core), so `fortran np=1` cells passed silently while every
    n>1 fortran cell failed outright."""
    cpu_list = ' '.join(str(c) for c in cpus[:n])
    cpu2node = cpu_to_node(node_map)
    node_list = ' '.join(str(cpu2node[c]) for c in cpus[:n])
    inner = ('CPUS=(%s); NODES=(%s); exec numactl '
             '--physcpubind=${CPUS[$OMPI_COMM_WORLD_LOCAL_RANK]} '
             '--membind=${NODES[$OMPI_COMM_WORLD_LOCAL_RANK]} %s'
             % (cpu_list, node_list, binary))
    return "mpirun --bind-to none -np %d bash -c '%s'" % (n, inner)


def run_fortran(work, n, policy, cpus, node_map, term, nsteps):
    """`cpus` MUST be the already-selected, already-busy-checked cpu list
    from the caller (`select_cpus(free, n, policy)` + `measure_once`) --
    this function used to recompute cpus itself via
    `compact_cpus(node_map, n)` on the FULL, unfiltered topology, silently
    discarding the free-node-aware selection main() had just validated. On
    a box with scattered per-cpu (not per-node) foreign load that meant the
    busy check and the actual mpirun pin could target two different cpu
    sets -- the check passed on cpus the binary never used. Fixed by taking
    `cpus` as a parameter instead of a second, divergent computation."""
    nx, ny, nz = DECOMP[n]
    d = os.path.join(work, f'f{n}_{policy}_{nsteps}')
    make_case(d, nx, ny, nz, term=term)
    binary = os.path.join(ROOT, 'bin', 'eqdyna')
    if not os.path.exists(binary):
        raise SystemExit(f'FAIL: {binary} missing -- build it first '
                         '(./install-eqdyna.sh -m ubuntu).')
    t0 = time.time()
    r = sh(fortran_cmd(n, cpus, node_map, binary), cwd=d)
    wall = time.time() - t0
    # Profile-guard/collection wiring (item 3, testsys half, 2026-09-23):
    # captured HERE, before the rmtree two lines down destroys `d` -- this is
    # the only point in the fortran path where the run's own
    # profile.rank<r>.json files still exist on disk. Degraded to a WARNING,
    # never allowed to turn a scaling measurement red: this tool's verdict is
    # the ms/step figure, not profile coverage (unlike run_e2e's per-cell
    # hook, which deliberately CAN fail a cell -- see its own comment).
    try:
        _sha = sh(f'git -C {ROOT} rev-parse --short HEAD').stdout.strip()
        profile_record.capture_run(d, case=CASE, backend='fortran', ranks=n,
                                   term='perf-scaling-probe', sha=_sha,
            tree_dirty=profile_record.ledger.tree_dirty())
    except Exception as exc:                        # noqa: BLE001
        print('WARNING: profile-record capture failed for fortran n=%d (%s: '
             '%s) -- the scaling measurement above is unaffected.'
             % (n, type(exc).__name__, exc))
    shutil.rmtree(d, ignore_errors=True)
    return wall, r.stdout


def record_point(base, base_n, mode, n, ms_per_step, first_n):
    """Pure bookkeeping for the speedup baseline (item 91g), factored out of
    the callers' loops so it is directly unit-testable with no jax/subprocess.

    Mutates `base`/`base_n` (dicts keyed by mode) in place the first time
    `mode` is seen, and returns (speedup, baseline_n, note). `note` is a
    non-None warning string exactly when the baseline point being set is NOT
    `first_n` -- i.e. the intended n=1 (or whatever devices[0] is) baseline
    was skipped or failed and this mode's speedup column is grounded on a
    later point instead. The caller must not silently drop `note`."""
    note = None
    if mode not in base:
        base[mode] = ms_per_step
        base_n[mode] = n
        if n != first_n:
            note = ('baseline for mode %r is n=%d (n=%d, the requested '
                    'first point, was SKIPPED or FAILED) -- speedup values '
                    'for this mode are relative to n=%d, not n=%d'
                    % (mode, n, first_n, n, first_n))
    return base[mode] / ms_per_step, base_n[mode], note


def per_step_fortran(work, n, policy, cpus, node_map, dt, n_lo, n_hi):
    """PER-STEP BY DIFFERENCE for Fortran too (2026-09-19 fix, item 33
    headline measurement).

    This function did not exist: the Fortran branch took ONE wall-clock run
    of `n_hi` steps and reported `wall / n_hi` as ms/step, while the python
    branch used `per_step_py` and subtracted its fixed cost exactly. That is
    not a comparison -- it charges Fortran for mesh generation, input read,
    MPI setup and netCDF open, and charges jax for none of its equivalents.
    The bias runs ONE WAY, in jax's favour, and it is large at the short step
    counts this tool uses: measured at np=1/n_hi=15, `wall/n_hi` gave 1168
    ms/step against a by-difference value several times smaller. Any "jax
    beats Fortran" verdict taken off the old number was an artifact of the
    metric, not a property of either solver.

    The module docstring's original justification -- "a compiled binary with
    no JIT to amortise, and its fixed cost is genuinely small next to a
    multi-second MPI solve" -- is true only for a LONG run. At 15-60 steps
    the fixed cost is a large fraction of the wall time, so it must cancel
    the same way it cancels for python: two runs, separate processes, same
    pin, divide the difference.

    `term` is set to `dt * nsteps` because the solver computes
    `nstep = idnint(term/dt)` (readInputFiles.f90), so the requested step
    count is exact rather than approximate. Both step counts are returned
    with the figure -- item 40's 4.0x is permanently unreproducible for
    exactly the lack of them."""
    t_lo, out_lo = run_fortran(work, n, policy, cpus, node_map, dt * n_lo, n_lo)
    t_hi, out_hi = run_fortran(work, n, policy, cpus, node_map, dt * n_hi, n_hi)
    ps, fixed = numa.per_step_and_fixed(t_lo, t_hi, n_lo, n_hi)
    return ps, fixed, t_lo, t_hi, out_lo, out_hi


# --------------------------------------------------------------------------
# python: per-step by difference, numactl-pinned, fresh process each call
# --------------------------------------------------------------------------
def build_py_case(case_name):
    """create.newcase + forced serial + case.setup, under testsys/perf/.

    GATE 0 (item 77): the exclusive lock on the directory holding the case,
    taken before anything is imported, created or deleted. It lives HERE, in
    the function that does the rmtree, and not in main(): main() is the path
    this function is LEAST often reached by -- `run_mpi_scaling.py:408`,
    `run_jaxmpi_ab.py:169`, `run_setup_probe.py:142` and
    `run_shard_scaling.py:195` all call it directly as `rs.build_py_case` and
    never reach this module's main(). A lock in main() would guard one of five
    callers and look like it guarded all five.

    `perflib.acquire_case_lock` derives the resource from the directory
    actually about to be destroyed, so it follows the path rather than
    restating it -- and roots it in THIS checkout, never in $EQDYNAROOT.

    The lock is KEPT after this returns, and memoised: the five callers above
    time out of this directory afterwards, and flock is per open file
    description, so a second acquire from this same process would refuse
    against our own pid -- a false collision. The memo is not a fallback; the
    lock is genuinely held.
    """
    global _case_lock
    d = os.path.join(TESTSYS, 'scaling_case', case_name)
    if _case_lock is None:
        _case_lock = perflib.acquire_case_lock(d, LOCK_CONSEQUENCE)
    return perflib.rebuild_serial_case(case_name, d)


def time_one_py(case_dir, nsteps, cpus, node_map, backend):
    script = (
        "import sys; sys.path.insert(0, %r)\n"
        "from eqdyna import eqdyna3d\n"
        "import time; t0 = time.time()\n"
        "eqdyna3d.run_case(%r, nsteps=%d, verbose=False, backend=%r)\n"
        "print('WALL', time.time() - t0)\n" % (PYTHON_PKG, case_dir, nsteps, backend))
    env = dict(os.environ)
    env['JAX_PLATFORMS'] = 'cpu'
    # Do NOT force XLA_FLAGS/OMP_NUM_THREADS/OPENBLAS_NUM_THREADS: both XLA and
    # OpenBLAS auto-size their thread pools from the process's cpu affinity at
    # import/first-use time, which numactl already sets correctly below (see
    # module docstring for the measured evidence). Forcing them here would
    # hide whether the PORT itself pins correctly when called unpinned.
    env.pop('XLA_FLAGS', None)
    env.pop('OMP_NUM_THREADS', None)
    env.pop('OPENBLAS_NUM_THREADS', None)
    cmd = numactl_prefix(cpus, node_map) + [sys.executable, '-c', script]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-1500:]); print(r.stderr[-1500:])
        return None
    for line in r.stdout.splitlines():
        if line.startswith('WALL'):
            return float(line.split()[1])
    return None


def per_step_py(case_dir, cpus, node_map, backend, n_lo, n_hi):
    t_lo = time_one_py(case_dir, n_lo, cpus, node_map, backend)
    t_hi = time_one_py(case_dir, n_hi, cpus, node_map, backend)
    if t_lo is None or t_hi is None:
        return None, None, None, None
    # Profile-guard/collection wiring (item 3, testsys half, 2026-09-23):
    # `case_dir` is REUSED across calls (build_py_case, never rmtree'd
    # between them), so the profile.rank0.json resident right after the
    # n_hi call above is that run's own -- captured here, right after the
    # LARGER run, same convention as run_perf.py's single-nsteps capture.
    # Warn-only: this tool's verdict is the ms/step figure, not profile
    # coverage.
    try:
        _sha = sh(f'git -C {ROOT} rev-parse --short HEAD').stdout.strip()
        profile_record.capture_run(case_dir, case=CASE,
                                   backend='python-%s' % backend, ranks=1,
                                   term='perf-scaling-probe', sha=_sha,
            tree_dirty=profile_record.ledger.tree_dirty())
    except Exception as exc:                        # noqa: BLE001
        print('WARNING: profile-record capture failed for python-%s (%s: '
             '%s) -- the scaling measurement above is unaffected.'
             % (backend, type(exc).__name__, exc))
    ps, fixed = numa.per_step_and_fixed(t_lo, t_hi, n_lo, n_hi)
    return ps, fixed, t_lo, t_hi


def main():
    global CASE
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default=CASE)
    ap.add_argument('--fortran-ranks', default=','.join(str(x) for x in FORTRAN_RANKS))
    ap.add_argument('--py-threads', default=','.join(str(x) for x in PY_THREADS))
    ap.add_argument('--policies', default='compact,spread')
    ap.add_argument('--backends', default='numpy,jax')
    ap.add_argument('--n-lo', type=int, default=20)
    ap.add_argument('--n-hi', type=int, default=60)
    ap.add_argument('--repeats', type=int, default=1)
    ap.add_argument('--busy-ceiling', type=float, default=0.2)
    ap.add_argument('--i-know-the-box-is-busy', action='store_true')
    ap.add_argument('--skip-fortran', action='store_true')
    a = ap.parse_args()
    CASE = a.case
    fortran_ranks = [int(x) for x in a.fortran_ranks.split(',') if x]
    py_threads = [int(x) for x in a.py_threads.split(',') if x]
    policies = a.policies.split(',')
    backends = a.backends.split(',')

    for n in fortran_ranks + py_threads:
        if n > 32:
            raise SystemExit('FAIL: %d cores requested -- owner constraint is at most 32 '
                             '(one socket), never 64.' % n)

    nodes = numa.numa_topology()
    if not nodes:
        raise SystemExit('FAIL: numactl --hardware gave no topology; this experiment is '
                         'about NUMA placement and cannot run blind.')
    per_node = len(nodes[min(nodes)])
    print('topology: %d NUMA node(s) x %d cpu(s)' % (len(nodes), per_node))

    sha = sh(f'git -C {ROOT} rev-parse --short HEAD').stdout.strip()
    host = os.uname().nodename
    work = tempfile.mkdtemp(prefix='scaling.', dir=os.environ.get('TMPDIR', '/tmp'))
    rows = []
    skipped = []

    def measure_once(label, cpus):
        chk = numa.require_idle(cpus, a.busy_ceiling, a.i_know_the_box_is_busy)
        if chk is None:
            print(f'  {label:<28} SKIPPED (cpus busy, not overridden)')
            skipped.append(dict(label=label, cpus=cpus))
            return None
        return chk

    if not a.skip_fortran:
        print('\n-- Fortran MPI (%s) --' % CASE)
        d0 = os.path.join(work, 'dtprobe')
        make_case(d0, 1, 1, 1)
        _, dt = read_term_dt(d0)
        shutil.rmtree(d0, ignore_errors=True)
        base, base_n = {}, {}
        for policy in policies:
            for n in fortran_ranks:
                free = free_node_map(nodes, a.busy_ceiling, a.i_know_the_box_is_busy)
                label = f'fortran np={n:<3} {policy}'
                cpus = select_cpus(free, n, policy)
                if cpus is None:
                    print(f'  {label:<28} SKIPPED (no {n}-cpu {policy} placement among '
                          f'currently-free NUMA nodes {sorted(free)})')
                    skipped.append(dict(label=label, reason='no_free_node_placement',
                                        n=n, policy=policy, free_nodes=sorted(free)))
                    continue
                chk = measure_once(label, cpus)
                if chk is None:
                    continue
                best_ps, best_fixed, best_lo, best_hi = None, None, None, None
                for rep in range(a.repeats):
                    ps, fixed, t_lo, t_hi, _o1, _o2 = per_step_fortran(
                        work, n, policy, cpus, nodes, dt, a.n_lo, a.n_hi)
                    if best_ps is None or ps < best_ps:
                        best_ps, best_fixed, best_lo, best_hi = ps, fixed, t_lo, t_hi
                # Item 91g: the baseline point is labelled, never silently
                # swapped when the first n was skipped.
                speedup, bn, note = record_point(base, base_n, policy, n, best_ps,
                                                 fortran_ranks[0])
                if note:
                    print('  NOTE: ' + note, flush=True)
                used_nodes = nodes_of(nodes, cpus)
                row = dict(engine='fortran', n=n, policy=policy, ms_per_step=best_ps * 1e3,
                          speedup=speedup, baseline_n=bn, fixed_s=best_fixed, cpus=cpus, nodes=used_nodes,
                          n_lo=a.n_lo, n_hi=a.n_hi, wall_lo_s=best_lo, wall_hi_s=best_hi,
                          busy=chk, loadavg=os.getloadavg())
                rows.append(row)
                print(f'  {label:<28} {best_ps*1e3:9.2f} ms/step  speedup {speedup:5.2f}x (vs n={bn})  '
                      f'fixed {best_fixed:6.2f}s  n_lo/n_hi {a.n_lo}/{a.n_hi}  '
                      f'wall {best_lo:.1f}/{best_hi:.1f}s  cpus={cpus} nodes={used_nodes}',
                      flush=True)

    for backend in backends:
        py_backend = 'python-%s' % backend
        print('\n-- %s --' % py_backend)
        case_dir = build_py_case(CASE)
        base, base_n = {}, {}
        for policy in policies:
            for n in py_threads:
                free = free_node_map(nodes, a.busy_ceiling, a.i_know_the_box_is_busy)
                label = f'{py_backend} th={n:<3} {policy}'
                cpus = select_cpus(free, n, policy)
                if cpus is None:
                    print(f'  {label:<28} SKIPPED (no {n}-cpu {policy} placement among '
                          f'currently-free NUMA nodes {sorted(free)})')
                    skipped.append(dict(label=label, reason='no_free_node_placement',
                                        n=n, policy=policy, free_nodes=sorted(free)))
                    continue
                chk = measure_once(label, cpus)
                if chk is None:
                    continue
                best_ps, best_fixed, best_lo, best_hi = None, None, None, None
                for rep in range(a.repeats):
                    ps, fixed, t_lo, t_hi = per_step_py(case_dir, cpus, nodes, backend,
                                                        a.n_lo, a.n_hi)
                    if ps is None:
                        continue
                    if best_ps is None or ps < best_ps:
                        best_ps, best_fixed, best_lo, best_hi = ps, fixed, t_lo, t_hi
                if best_ps is None:
                    print(f'  {label:<28} FAILED')
                    continue
                speedup, bn, note = record_point(base, base_n, policy, n, best_ps,
                                                 py_threads[0])
                if note:
                    print('  NOTE: ' + note, flush=True)
                used_nodes = nodes_of(nodes, cpus)
                row = dict(engine=py_backend, n=n, policy=policy, ms_per_step=best_ps * 1e3,
                          speedup=speedup, baseline_n=bn, fixed_s=best_fixed, cpus=cpus, nodes=used_nodes,
                          n_lo=a.n_lo, n_hi=a.n_hi, wall_lo_s=best_lo, wall_hi_s=best_hi,
                          busy=chk, loadavg=os.getloadavg())
                rows.append(row)
                print(f'  {label:<28} {best_ps*1e3:9.2f} ms/step  speedup {speedup:5.2f}x (vs n={bn})  '
                      f'fixed {best_fixed:6.2f}s  n_lo/n_hi {a.n_lo}/{a.n_hi}  '
                      f'wall {best_lo:.1f}/{best_hi:.1f}s  cpus={cpus} nodes={used_nodes}',
                      flush=True)

    meta = dict(case=CASE, sha=sha, host=host, loadavg=os.getloadavg(),
               date=time.strftime('%Y-%m-%d %H:%M'),
               topology={str(k): v for k, v in nodes.items()},
               busy_ceiling=a.busy_ceiling, overridden=bool(a.i_know_the_box_is_busy),
               n_lo=a.n_lo, n_hi=a.n_hi, skipped=skipped, rows=rows)
    json.dump(meta, open(OUT, 'w'), indent=1)
    print(f'\nprovenance: {CASE} @ {sha} on {host}, loadavg {os.getloadavg()}; saved {OUT}')

    # Durable copy + ledger index (PROJECT_RULES rule 19 and its positive
    # counterpart). scaling_last.json stays what its name says -- the last
    # run -- but the run ALSO lands under a dated, immutable name, and every
    # measured point becomes one appended line in docs/perf_ledger.jsonl so
    # no number depends on someone transcribing it. Timestamp in the filename
    # so a same-day rerun cannot overwrite the earlier snapshot.
    import ledger
    snap_rel = os.path.join('docs', 'perf_snapshots', 'scaling_%s_%s.json'
                            % (time.strftime('%Y-%m-%d_%H%M%S'),
                               os.environ.get('EQDYNA_SNAPSHOT_TAG', CASE)))
    snap_abs = os.path.join(ROOT, snap_rel)
    os.makedirs(os.path.dirname(snap_abs), exist_ok=True)
    json.dump(meta, open(snap_abs, 'w'), indent=1)
    tenancy = ledger.box_tenancy(a.busy_ceiling)
    nledger = ledger.append_rows(
        ledger.rows_from_scaling_snapshot(meta, snap_rel, tenancy))
    print('snapshot %s; %d ledger row(s) appended to %s (box tenancy %d/%d '
          'cpus over %.2f)' % (snap_rel, nledger, ledger.LEDGER_RELPATH,
                               tenancy['busy'], tenancy['total'],
                               a.busy_ceiling))
    if skipped:
        print(f'{len(skipped)} configuration(s) SKIPPED as busy (not measured, not silently '
              f'dropped): {[s["label"] for s in skipped]}')
    shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    main()
