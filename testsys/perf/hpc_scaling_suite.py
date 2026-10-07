#! /usr/bin/env python3
"""
testsys/perf/hpc_scaling_suite.py -- board item 142 (owner, 2026-10-01; see
pathway_forward.md, grep '| 142 |' for the full row).

ONE SCRIPT, THREE MODES, NEVER MIXED IN ONE INVOCATION:

  --submit --machine M --account A --backend fortran|python-jax-mpi
           --test strong|weak|size --level medium|large|big [--max-nodes N]
      Builds a suite of scaling POINTS for ONE backend/test/level and either
      sbatch's them (a scheduler machine, scripts/machines.py) or runs them
      in place (a no-scheduler machine). Writes `manifest.json` under the
      suite directory describing every point's run directories, so --collect
      and --analyze never have to re-derive what was submitted.

  --collect --suite-dir DIR
      Meant to run ON THE HPC SIDE once --submit's jobs have finished.
      Validates every point's always-on per-rank `profile.rank<r>.json`
      files (testsys/profile_schema.py -- the SAME validator every other
      profiled run in this repo is held to) and packs manifest.json plus
      every validated profile file into one tarball.

  --analyze TGZ
      Runs LOCALLY. Unpacks the tarball and computes the metrics this
      module's `analyze()` documents below.

Backends run SEPARATELY: one --submit invocation is one backend, hence one
suite directory, one tarball, and one --analyze report. Never mix fortran
and python-jax-mpi points in the same suite -- `_refuse_mixed_backend` in
--collect enforces this mechanically rather than relying on nobody doing it
by hand.

TEST TYPES (owner spec, verbatim):
  strong  fixed TOTAL element count; cores double each step. Only
          par.nx/par.ny/par.nz change between points (the MPI decomposition
          that redistributes the SAME mesh across more ranks) -- the mesh
          itself (par.dx, hence element count) is held fixed across every
          point of a strong-scaling suite.
  weak    elements-PER-CORE held fixed; par.dx is refined as the rank grid
          grows, so total elements grows with core count.
  size    fixed core count; par.dx is swept to sweep total elements,
          looking for the elements-per-core sweet spot and bytes/element.
          This is the ONE test type this box must never actually run (see
          "What this session deliberately did not run" below) -- level
          'big' additionally refuses python-jax-mpi outright, Fortran only,
          enforced in `_validate_combo`, not just documented.

LEVELS: medium 10^7-10^8 total elements (sized for 1-16 LS6 nodes), large
10^9, big 10^10.

PER-STEP-BY-DIFFERENCE (CLAUDE.md "things that will bite you"): a single
wall-clock-over-step-count number is NOT how this repo measures per-step
cost -- fixed cost (interpreter start, MPI init, netCDF open, case load,
XLA compile where applicable) must be differenced out over TWO step counts,
and a non-positive result from that difference is a BUG to raise on, never a
number to report. Every point this suite submits is therefore actually TWO
runs, POINT_NSTEPS_LO and POINT_NSTEPS_HI, each recorded as its own run
directory in the manifest, so --analyze can call
`testsys.perf.run_numa_scaling.per_step_and_fixed(t_lo, t_hi, n_lo, n_hi)`
-- the ONE copy of this arithmetic already shared by run_scaling.py,
run_shard_scaling.py and run_numa_scaling.py's own `per_step` -- on each
RANK'S OWN `total_s` (never a wrapper wall clock: that is exactly the
>=4-rank pitfall this same CLAUDE.md section names, measured there as a
negative per-step figure). The same arithmetic is applied per-bucket
(`buckets_s`) to get a per-step BUCKET breakdown, not just a per-step total.

METRICS `analyze()` computes, all sourced from the per-rank profile files
reached through the manifest (nothing here invents a metric
`profile_emit.py` does not actually record):
  - per-step cost (per rank, per point)
  - bucket breakdown: element/fault/exchange/wait/io, the exact
    `profile_schema.BUCKET_KEYS` -- the buckets the emitter writes, no
    others invented here.
  - max/mean imbalance across ranks, per point
  - the 128->256-node step boundary, if the suite's points actually span
    it (reported, never fabricated when absent)
  - setup-and-io cost vs rank count (the "fixed" intercept
    `per_step_and_fixed` already separates out, tracked across points)
  - a parity check of every timed (hi) run against its case's registered
    e2e reference, reusing testsys/compare.py's `compare_cell` -- the ONE
    comparison in this repo, not a second one written here. A point whose
    mesh was refined away from the case's gate resolution (weak/size tests,
    and every strong-scaling point past the first) is NOT comparable to
    that reference; those points are DECLARED (not silently skipped), with
    the measured dx that makes them incomparable printed in the report.

WHAT THIS SESSION DELIBERATELY DID NOT RUN (2026-10-06, general-purpose-hpc-
scaling): the box this was built on measured swap already full (31Gi/31Gi)
in a prior session (pathway row 142's own prerequisite note), so no
`--submit --test size` sweep -- nor any other real multi-GB `--submit`
job -- was executed here. `--collect`/`--analyze` are exercised against (a)
the real, already-committed `testsys/regression/fixtures/profile_guard/`
profile files (a real test.tpv8 run, well within memory) and (b) small
synthetic profile.rank*.json fixtures built in-memory by this module's own
regression test, never a from-scratch "sweep dx until OOM" run.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))              # testsys/perf
TESTSYS = os.path.dirname(HERE)                                 # testsys/
REPO_ROOT = os.path.dirname(TESTSYS)

for _p in (HERE, TESTSYS, REPO_ROOT, os.path.join(REPO_ROOT, 'src', 'python')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_numa_scaling as numa          # noqa: E402  (per_step_and_fixed, NonPositiveMeasurement)
import profile_schema                    # noqa: E402  (validate_run_dir -- the one profile guard)
from testsys import matrix, compare      # noqa: E402  (the one e2e comparison)
from eqdyna.MPI4NodalQuant import DECOMP as GATED_DECOMP  # noqa: E402  (Fortran's own decomposition table)


def _load_machines():
    """scripts/machines.py, the ONE machine registry -- loaded by path
    (scripts/ is not a package), same mechanism as testsys/run.py's own
    `_load_machines`, so putting scripts/ on sys.path cannot shadow anything
    else."""
    import importlib.util
    path = os.path.join(REPO_ROOT, 'scripts', 'machines.py')
    spec = importlib.util.spec_from_file_location('eqdyna_machines', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


machines = _load_machines()

MANIFEST_NAME = 'manifest.json'
SUITE_SCHEMA_ID = 'eqdyna-hpc-scaling-suite/1'

BACKENDS = ('fortran', 'python-jax-mpi')
TEST_TYPES = ('strong', 'weak', 'size')
LEVELS = ('medium', 'large', 'big')

# Total-element TARGETS per level (owner spec: medium 10^7-10^8, large 10^9,
# big 10^10). One representative target per level, not a range, so every
# point this module builds has a single well-defined mesh to report --
# the full range is available via --target-elements for a caller who wants
# a different point in medium's two-decade span.
LEVEL_TARGET_ELEMENTS = {
    'medium': 3.0e7,
    'large': 1.0e9,
    'big': 1.0e10,
}

# Domain extents test.tpv8 (this suite's default base case) ships in
# scripts/defaultParameters.py (xmin/xmax, ymin/ymax, zmin/zmax). Not
# imported from there directly -- defaultParameters.py is a case-build-time
# script with its own print()/par-object side effects, not a side-effect-
# free constants module -- so these are a dated, commented COPY, same
# discipline as run_scaling.py's own CASE/term constants.
DOMAIN_X_M = 84.0e3   # xmax - xmin = 42e3 - (-42e3)
DOMAIN_Y_M = 42.0e3   # ymax - ymin = 22e3 - (-20e3)
DOMAIN_Z_M = 42.0e3   # zmax - zmin = 0 - (-42e3)
DOMAIN_VOLUME_M3 = DOMAIN_X_M * DOMAIN_Y_M * DOMAIN_Z_M

# The two step counts every point runs at, for per-step-by-difference. Kept
# small on purpose: this suite's job is to measure PER-STEP cost and bucket
# breakdown, not to run a case to physical completion -- see run_scaling.py's
# own n_lo/n_hi defaults for the same reasoning.
POINT_NSTEPS_LO = 5
POINT_NSTEPS_HI = 15

DEFAULT_CASE = 'test.tpv8'

# Fortran's own perf/e2e decomposition table (src/python/eqdyna/
# MPI4NodalQuant.DECOMP, imported above as GATED_DECOMP), reused so any rank
# count this suite shares with the gated backend uses the SAME (npx,npy,npz)
# split. Beyond this table (needed for the larger node counts 'large'/'big'
# submit), `_decompose` below falls back to a balanced factorization --
# documented at that function.
GATED_DECOMP = dict(GATED_DECOMP)


class SuiteError(RuntimeError):
    """Raised on any refused combination (big x python-jax-mpi, a mixed-
    backend suite directory, ...) -- a hard failure, never a silent skip."""


# --------------------------------------------------------------------------
# validation shared by --submit and --collect
# --------------------------------------------------------------------------
def _validate_combo(test, level, backend):
    if test not in TEST_TYPES:
        raise SuiteError('unknown --test %r (known: %s)' % (test, TEST_TYPES))
    if level not in LEVELS:
        raise SuiteError('unknown --level %r (known: %s)' % (level, LEVELS))
    if backend not in BACKENDS:
        raise SuiteError('unknown --backend %r (known: %s)' % (backend, BACKENDS))
    if level == 'big' and backend == 'python-jax-mpi':
        raise SuiteError(
            "level 'big' (10^10 elements) refuses --backend python-jax-mpi "
            "per the owner's spec (pathway_forward.md row 142): Fortran "
            "only at this level. Pass --backend fortran, or choose a "
            "smaller --level.")


# --------------------------------------------------------------------------
# decomposition
# --------------------------------------------------------------------------
def _factor_balanced(n):
    """(npx, npy, npz) for a rank count outside GATED_DECOMP: the most
    cube-like factorization of n into 3 integer factors (npx>=npy>=npz>=1).
    Deterministic and side-effect-free; used only for suite rank counts this
    repo's own gated decomposition table does not cover (>32)."""
    best = (n, 1, 1)
    best_spread = n  # max-min, smaller is more cube-like
    for a in range(1, int(round(n ** (1 / 3.0))) + 2):
        if n % a:
            continue
        rem = n // a
        for b in range(a, int(round(rem ** 0.5)) + 2):
            if rem % b:
                continue
            c = rem // b
            if b > c:
                continue
            spread = max(a, b, c) - min(a, b, c)
            if spread < best_spread:
                best_spread, best = spread, (c, b, a)
    npz, npy, npx = sorted(best)
    return npx, npy, npz


def _decompose(nranks):
    if nranks in GATED_DECOMP:
        return GATED_DECOMP[nranks]
    return _factor_balanced(nranks)


# --------------------------------------------------------------------------
# element-count <-> dx
# --------------------------------------------------------------------------
def dx_for_target_elements(target_elements):
    """Uniform cube spacing (m) whose element count over the default
    test.tpv8 domain volume is closest to `target_elements`. elements =
    (Lx/dx)*(Ly/dx)*(Lz/dx) = volume/dx**3 for a uniform mesh, so
    dx = (volume/elements)**(1/3)."""
    if target_elements <= 0:
        raise SuiteError('target_elements %r must be > 0' % (target_elements,))
    return (DOMAIN_VOLUME_M3 / float(target_elements)) ** (1.0 / 3.0)


def elements_for_dx(dx):
    if dx <= 0:
        raise SuiteError('dx %r must be > 0' % (dx,))
    return DOMAIN_VOLUME_M3 / (dx ** 3)


# --------------------------------------------------------------------------
# points: the suite's table of (ranks, dx) configurations
# --------------------------------------------------------------------------
def _node_counts(max_nodes, cores_per_node):
    """[1, 2, 4, ... ] node counts, doubling, capped at max_nodes (and never
    below 1)."""
    if cores_per_node is None:
        raise SuiteError(
            "machine's cores_per_node is not verified (scripts/machines.py) "
            "-- --submit refuses to guess a node's core count rather than "
            "silently assuming one (rule 2)")
    out, n = [], 1
    while n <= max_nodes:
        out.append(n)
        n *= 2
    if out[-1] != max_nodes:
        out.append(max_nodes)
    return out


def build_points(test, level, machine_cfg, max_nodes):
    """[(point_id, ranks, nnodes, dx, (npx,npy,npz)), ...] for one
    (test, level) combination. Pure function (no filesystem/subprocess), so
    it is directly unit-testable."""
    cores_per_node = machine_cfg['cores_per_node']
    target = LEVEL_TARGET_ELEMENTS[level]
    points = []

    if test == 'strong':
        # Fixed total elements; only the decomposition changes.
        dx = dx_for_target_elements(target)
        for nnodes in _node_counts(max_nodes, cores_per_node):
            ranks = nnodes * cores_per_node
            points.append(dict(point_id='strong_n%04d' % nnodes, ranks=ranks,
                                nnodes=nnodes, dx=dx,
                                decomp=_decompose(ranks)))
    elif test == 'weak':
        # Elements-per-core fixed: dx refines (shrinks) as ranks grow so
        # that target/ranks (elements per core at nnodes=1) stays constant.
        base_ranks = cores_per_node
        for nnodes in _node_counts(max_nodes, cores_per_node):
            ranks = nnodes * cores_per_node
            this_target = target * (ranks / float(base_ranks))
            dx = dx_for_target_elements(this_target)
            points.append(dict(point_id='weak_n%04d' % nnodes, ranks=ranks,
                                nnodes=nnodes, dx=dx,
                                decomp=_decompose(ranks)))
    else:  # size
        # Fixed core count (one node's worth, by default); dx swept across
        # a decade either side of the level's target so --analyze can look
        # for the elements-per-core sweet spot. The ACTUAL run of this sweep
        # "until memory runs out" is the step this session deliberately did
        # not execute (see module docstring) -- build_points only produces
        # the point TABLE, same as the other two test types.
        ranks = cores_per_node
        nnodes = 1
        for i, factor in enumerate((0.1, 0.3, 1.0, 3.0, 10.0)):
            dx = dx_for_target_elements(target * factor)
            points.append(dict(point_id='size_f%02d' % i, ranks=ranks,
                                nnodes=nnodes, dx=dx,
                                decomp=_decompose(ranks)))
    return points


# --------------------------------------------------------------------------
# case building (adapted from run_scaling.py's make_case/read_term_dt --
# same mechanism, generalized to also set par.dx)
# --------------------------------------------------------------------------
def _sh(cmd, cwd, env=None):
    rc = subprocess.call(cmd, shell=isinstance(cmd, str), cwd=cwd, env=env)
    if rc != 0:
        raise RuntimeError('`%s` (cwd=%s) exited %d' % (cmd, cwd, rc))


def build_point_case(dst, case_name, dx, decomp, nsteps, env=None):
    """create.newcase + par.dx/par.nx,ny,nz overrides + case.setup, twice:
    once to discover (term, dt) via bGlobal.txt (run_scaling.read_term_dt's
    exact mechanism -- the only place this repo already reads that file
    back), then again with par.term set to reach exactly `nsteps` steps at
    that case's own dt. Returns the case directory."""
    if os.path.exists(dst):
        shutil.rmtree(dst)
    npx, npy, npz = decomp
    _sh(['create.newcase', dst, case_name], cwd=os.path.dirname(dst), env=env)
    params = os.path.join(dst, 'user_defined_params.py')
    text = open(params).read()
    text = re.sub(r'^par\.dx\s*=.*', 'par.dx = %r' % (dx,), text, flags=re.M)
    if not re.search(r'^par\.dx\s*=', text, flags=re.M):
        text += '\npar.dx = %r\n' % (dx,)
    text += ('\n# hpc_scaling_suite.py point override\n'
             'par.nx, par.ny, par.nz = %d, %d, %d\n' % (npx, npy, npz))
    open(params, 'w').write(text)
    _sh(['./case.setup'], cwd=dst, env=env)

    term0, dt = (float(x) for x in
                 open(os.path.join(dst, 'bGlobal.txt')).read().splitlines()[15:17])
    term = nsteps * dt
    text = open(params).read()
    text += '\npar.term = %r  # hpc_scaling_suite.py: exactly %d steps at dt=%r\n' \
            % (term, nsteps, dt)
    open(params, 'w').write(text)
    _sh(['./case.setup'], cwd=dst, env=env)
    return dst


def run_point_fortran(case_dir, ranks, env=None):
    eqdyna_cmd = os.path.join(REPO_ROOT, 'bin', 'eqdyna')
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    _sh([mpirun, '-np', str(ranks), eqdyna_cmd], cwd=case_dir, env=env)


def run_point_python_jax_mpi(case_dir, ranks, env=None):
    mpirun = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
    _sh([mpirun, '-np', str(ranks), sys.executable, '-m', 'eqdyna',
         '--backend', 'jax', '--mpi'], cwd=case_dir, env=env)


# --------------------------------------------------------------------------
# --submit
# --------------------------------------------------------------------------
def submit_point(point, backend, case_name, machine, machine_cfg, account, out_dir, env):
    """Builds the point's two (lo/hi nsteps) cases and either sbatch's a
    job script (scheduler machine) or runs both in place (no-scheduler
    machine). Returns the point's manifest entry."""
    entry = dict(point)
    entry['backend'] = backend
    entry['runs'] = {}
    for tag, nsteps in (('lo', POINT_NSTEPS_LO), ('hi', POINT_NSTEPS_HI)):
        run_dir = os.path.join(out_dir, point['point_id'], tag)
        build_point_case(run_dir, case_name, point['dx'], point['decomp'], nsteps, env)
        entry['runs'][tag] = dict(run_dir=run_dir, nsteps=nsteps)

    if machine_cfg['scheduler'] == 'slurm':
        script_path = os.path.join(out_dir, point['point_id'], 'submit.sh')
        with open(script_path, 'w') as f:
            machines.write_slurm_header(
                f, jobname='hpcscale-%s' % point['point_id'], nnode=point['nnodes'],
                ncpu=point['ranks'], queue=machine_cfg['partition'] or 'normal',
                walltime=machine_cfg['walltime'] or '02:00:00', account=account, email='')
            mpirun = machine_cfg['mpirun']
            if mpirun is None:
                raise SuiteError(
                    "machine %r has no verified launcher (scripts/machines.py "
                    "mpirun=None) -- --submit refuses to guess one" % (machine,))
            for tag in ('lo', 'hi'):
                rd = entry['runs'][tag]['run_dir']
                if backend == 'fortran':
                    f.write('cd %s && %s -np %d %s\n'
                            % (rd, mpirun, point['ranks'],
                               os.path.join(REPO_ROOT, 'bin', 'eqdyna')))
                else:
                    f.write('cd %s && %s -np %d python3 -m eqdyna --backend jax --mpi\n'
                            % (rd, mpirun, point['ranks']))
        entry['submit_script'] = script_path
        out = subprocess.check_output(['sbatch', script_path], cwd=out_dir, env=env)
        entry['job_id'] = out.decode().strip()
    else:
        runner = run_point_fortran if backend == 'fortran' else run_point_python_jax_mpi
        for tag in ('lo', 'hi'):
            runner(entry['runs'][tag]['run_dir'], point['ranks'], env)
        entry['job_id'] = None
    return entry


def cmd_submit(args):
    _validate_combo(args.test, args.level, args.backend)
    machine_cfg = machines.machine(args.machine)
    max_nodes = args.max_nodes
    if max_nodes is None:
        max_nodes = {'medium': 16, 'large': 64, 'big': 256}[args.level]

    points = build_points(args.test, args.level, machine_cfg, max_nodes)
    out_dir = args.out_dir or os.path.join(
        REPO_ROOT, 'docs', 'evidence',
        'hpc-scaling-%s-%s-%s-%s' % (args.machine, args.backend, args.test, args.level))
    os.makedirs(out_dir, exist_ok=True)

    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([os.path.join(REPO_ROOT, 'bin'),
                                   os.path.join(REPO_ROOT, 'scripts'),
                                   env.get('PATH', '')])

    entries = [submit_point(p, args.backend, args.case, args.machine, machine_cfg,
                            args.account, out_dir, env) for p in points]

    manifest = dict(schema=SUITE_SCHEMA_ID, machine=args.machine, backend=args.backend,
                     test=args.test, level=args.level, case=args.case,
                     account=args.account, max_nodes=max_nodes,
                     nsteps_lo=POINT_NSTEPS_LO, nsteps_hi=POINT_NSTEPS_HI,
                     created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                     points=entries)
    with open(os.path.join(out_dir, MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print('hpc_scaling_suite: submitted %d point(s) to %s -- manifest at %s'
         % (len(entries), out_dir, os.path.join(out_dir, MANIFEST_NAME)))
    return 0


# --------------------------------------------------------------------------
# --collect
# --------------------------------------------------------------------------
def cmd_collect(args):
    suite_dir = args.suite_dir
    manifest_path = os.path.join(suite_dir, MANIFEST_NAME)
    if not os.path.isfile(manifest_path):
        raise SuiteError('no %s under %r -- was this suite directory produced '
                         'by --submit?' % (MANIFEST_NAME, suite_dir))
    with open(manifest_path) as f:
        manifest = json.load(f)
    if manifest.get('schema') != SUITE_SCHEMA_ID:
        raise SuiteError('manifest schema %r != %r' % (manifest.get('schema'), SUITE_SCHEMA_ID))

    backends_seen = {manifest['backend']}
    if len(backends_seen) != 1:
        raise SuiteError('a suite directory must carry exactly one backend; '
                         'saw %s -- backends run SEPARATELY' % (backends_seen,))

    validated = 0
    for point in manifest['points']:
        for tag, run in point['runs'].items():
            rows = profile_schema.validate_run_dir(run['run_dir'])  # raises on any problem
            if rows[0]['backend'] != manifest['backend']:
                raise SuiteError(
                    'point %s run %s: profile files report backend=%r but '
                    'manifest says backend=%r' % (point['point_id'], tag,
                                                  rows[0]['backend'], manifest['backend']))
            validated += len(rows)

    out_tgz = args.out or os.path.join(
        suite_dir, 'hpc_scaling_%s_%s_%s_%s.tgz'
        % (manifest['machine'], manifest['backend'], manifest['test'], manifest['level']))
    with tarfile.open(out_tgz, 'w:gz') as tar:
        tar.add(manifest_path, arcname=MANIFEST_NAME)
        for point in manifest['points']:
            for tag, run in point['runs'].items():
                for name in sorted(os.listdir(run['run_dir'])):
                    if profile_schema._RANK_FILE_RE.match(name):
                        arcname = os.path.join(point['point_id'], tag, name)
                        tar.add(os.path.join(run['run_dir'], name), arcname=arcname)
    print('hpc_scaling_suite: collected %d profile file(s) from %d point(s) into %s'
         % (validated, len(manifest['points']), out_tgz))
    return 0


# --------------------------------------------------------------------------
# --analyze
# --------------------------------------------------------------------------
def _load_point_profiles(tmpdir, manifest, point):
    """{'lo': [rows by rank], 'hi': [rows by rank]} for one manifest point,
    read back from the unpacked tarball."""
    out = {}
    for tag in ('lo', 'hi'):
        run_dir = os.path.join(tmpdir, point['point_id'], tag)
        out[tag] = profile_schema.validate_run_dir(run_dir)
    return out


# profile_emit.py's buckets split into two different KINDS of cost (same
# distinction CLAUDE.md draws for `setup`/`io` vs the step-loop buckets):
# PER_STEP_BUCKETS accumulate once per solver step and scale with nsteps, so
# per_step_and_fixed's differencing is the right arithmetic for them;
# FIXED_BUCKETS (setup, io) are a ONE-TIME cost each run pays regardless of
# how many steps it takes and do NOT scale with nsteps -- differencing them
# over nsteps is not "small", it is the WRONG QUESTION (lo and hi both pay
# roughly the same setup/io cost, so the difference is noise around zero and
# per_step_and_fixed correctly refuses it as non-positive/invalid). Fixed
# buckets are instead reported directly from the HI run (the longer, more
# representative measurement of the two).
PER_STEP_BUCKETS = ('element', 'fault', 'exchange', 'wait')
FIXED_BUCKETS = ('setup', 'io')
assert set(PER_STEP_BUCKETS) | set(FIXED_BUCKETS) == set(profile_schema.BUCKET_KEYS)


def per_rank_metrics(lo_rows, hi_rows):
    """Per-rank per-step-by-difference metrics for one point: total per-step
    cost, per-step-bucket breakdown (element/fault/exchange/wait), the fixed
    (non-per-step-scaling) intercept of the TOTAL, and the fixed buckets
    (setup/io) read directly off the hi run -- see PER_STEP_BUCKETS/
    FIXED_BUCKETS above for why those two groups are handled differently.
    All per-step arithmetic goes through run_numa_scaling.per_step_and_fixed,
    the ONE copy of it in this repo; it raises on any non-positive
    difference -- never reported as a number."""
    if len(lo_rows) != len(hi_rows):
        raise SuiteError('lo/hi profile rank counts differ: %d vs %d'
                         % (len(lo_rows), len(hi_rows)))
    per_rank = []
    for lo, hi in zip(lo_rows, hi_rows):
        if lo['rank'] != hi['rank']:
            raise SuiteError('lo/hi rank mismatch: %r vs %r' % (lo['rank'], hi['rank']))
        n_lo, n_hi = lo['nsteps'], hi['nsteps']
        ps_total, fixed_total = numa.per_step_and_fixed(lo['total_s'], hi['total_s'], n_lo, n_hi)
        buckets_per_step = {}
        for key in PER_STEP_BUCKETS:
            ps, _fixed = numa.per_step_and_fixed(
                lo['buckets_s'][key], hi['buckets_s'][key], n_lo, n_hi)
            buckets_per_step[key] = ps
        buckets_fixed = {key: hi['buckets_s'][key] for key in FIXED_BUCKETS}
        per_rank.append(dict(rank=lo['rank'], host=lo['host'],
                             per_step_s=ps_total, fixed_s=fixed_total,
                             buckets_per_step_s=buckets_per_step,
                             buckets_fixed_s=buckets_fixed))
    return per_rank


def imbalance(per_rank):
    vals = [r['per_step_s'] for r in per_rank]
    mean = sum(vals) / len(vals)
    return dict(max_s=max(vals), mean_s=mean, max_over_mean=max(vals) / mean)


def parity_for_point(manifest, point, tmpdir):
    """(status, lines): reuses testsys/compare.py's compare_cell against
    the HI run, but ONLY when this point's mesh is the case's unmodified
    gate resolution (dx == None sentinel meaning 'not overridden' is not
    tracked here, so the test is: is this the very FIRST point of a
    strong-scaling suite, at the level's own dx, for a case/backend the
    matrix table actually gates?). Every other point is DECLARED
    incomparable, with its own dx printed, per this repo's no-silent-skip
    rule -- never silently omitted."""
    case, backend = manifest['case'], manifest['backend']
    if not matrix.is_supported(case, backend):
        return 'declared', ['parity: %s x %s is not a gated cell in matrix.py '
                            '-- not applicable' % (case, backend)]
    if manifest['test'] != 'strong' or point['point_id'] != manifest['points'][0]['point_id']:
        return 'declared', [
            'parity: point %s ran at dx=%.3fm (%s test, level %s), which is '
            'NOT the case\'s registered gate resolution -- a mesh refined or '
            'coarsened away from the gate reference is not comparable to it; '
            'declared, not skipped' % (point['point_id'], point['dx'],
                                       manifest['test'], manifest['level'])]
    run_dir = os.path.join(tmpdir, point['point_id'], 'hi')
    try:
        ok, lines = compare.compare_cell(case, backend, run_dir)
    except Exception as e:
        # compare_cell raises (rather than returning ok=False) when the run
        # directory is missing an artifact it expected (e.g. no frt file at
        # all). That is itself a real finding -- surface it as a declared
        # failure with the exception text, never let it take the whole
        # --analyze command down for every other point in the suite.
        return 'fail', ['parity: compare_cell raised on %s x %s at %s: %r'
                        % (case, backend, run_dir, e)]
    return ('pass' if ok else 'fail'), lines


def analyze(tgz_path):
    """The report dict --analyze prints/returns. Pure function over an
    already-unpacked tarball's manifest+profiles (no argparse/CLI concerns),
    so it is directly unit-testable."""
    tmpdir = tempfile.mkdtemp(prefix='hpc_scaling_analyze_')
    try:
        with tarfile.open(tgz_path, 'r:gz') as tar:
            tar.extractall(tmpdir)
        with open(os.path.join(tmpdir, MANIFEST_NAME)) as f:
            manifest = json.load(f)

        report = dict(schema='eqdyna-hpc-scaling-report/1', machine=manifest['machine'],
                      backend=manifest['backend'], test=manifest['test'],
                      level=manifest['level'], case=manifest['case'], points=[])

        for point in manifest['points']:
            profiles = _load_point_profiles(tmpdir, manifest, point)
            per_rank = per_rank_metrics(profiles['lo'], profiles['hi'])
            imb = imbalance(per_rank)
            parity_status, parity_lines = parity_for_point(manifest, point, tmpdir)
            setup_io_fixed = sum(r['buckets_fixed_s']['setup'] + r['buckets_fixed_s']['io']
                                 for r in per_rank) / len(per_rank)
            report['points'].append(dict(
                point_id=point['point_id'], ranks=point['ranks'], nnodes=point['nnodes'],
                dx=point['dx'], decomp=point['decomp'],
                mean_per_step_s=sum(r['per_step_s'] for r in per_rank) / len(per_rank),
                max_per_step_s=imb['max_s'], max_over_mean=imb['max_over_mean'],
                setup_and_io_fixed_s=setup_io_fixed,
                bucket_breakdown_per_step_s={
                    key: sum(r['buckets_per_step_s'][key] for r in per_rank) / len(per_rank)
                    for key in PER_STEP_BUCKETS},
                parity_status=parity_status, parity_lines=parity_lines,
            ))

        # 128->256-node step boundary, if the suite's own points span it.
        by_nnodes = {p['nnodes']: p for p in report['points']}
        if 128 in by_nnodes and 256 in by_nnodes:
            a, b = by_nnodes[128]['mean_per_step_s'], by_nnodes[256]['mean_per_step_s']
            report['node_boundary_128_256'] = dict(
                per_step_s_128=a, per_step_s_256=b, ratio_256_over_128=b / a)
        else:
            report['node_boundary_128_256'] = None

        return report
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def cmd_analyze(args):
    report = analyze(args.tgz)
    text = json.dumps(report, indent=2, sort_keys=True)
    print(text)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(text + '\n')
    failed = [p['point_id'] for p in report['points'] if p['parity_status'] == 'fail']
    if failed:
        print('hpc_scaling_suite: FAIL -- parity mismatch at point(s): %s' % failed)
        return 1
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--submit', action='store_true')
    mode.add_argument('--collect', action='store_true')
    mode.add_argument('--analyze', metavar='TGZ')

    p.add_argument('--machine')
    p.add_argument('--account')
    p.add_argument('--backend', choices=BACKENDS)
    p.add_argument('--test', choices=TEST_TYPES)
    p.add_argument('--level', choices=LEVELS)
    p.add_argument('--max-nodes', type=int, default=None)
    p.add_argument('--case', default=DEFAULT_CASE)
    p.add_argument('--out-dir')
    p.add_argument('--suite-dir')
    p.add_argument('--out')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.submit:
            for name, val in (('--machine', args.machine), ('--account', args.account),
                              ('--backend', args.backend), ('--test', args.test),
                              ('--level', args.level)):
                if val is None:
                    raise SuiteError('--submit requires %s' % name)
            return cmd_submit(args)
        if args.collect:
            if not args.suite_dir:
                raise SuiteError('--collect requires --suite-dir')
            return cmd_collect(args)
        return cmd_analyze(argparse.Namespace(tgz=args.analyze, out=args.out))
    except SuiteError as exc:
        print('hpc_scaling_suite: FAIL -- %s' % exc)
        return 1


if __name__ == '__main__':
    sys.exit(main())
