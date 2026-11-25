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

  Row 142 additions (owner 2026-11-24): step pair settable (--nsteps-lo/-hi,
  default 100/200), --repeats N, --dry-run (writes the job script, builds
  nothing, submits nothing; prints the sbatch line with the placeholder
  <ALLOC>), --module/--mpirun/--walltime/--work-root/--bin for a machine whose
  registry entry is incomplete. The allocation is used only to build the
  sbatch command line: it is never written to the manifest or a job script.
  The manifest carries a provenance stamp; --collect records a parity verdict
  per timed run; --analyze reports overhead (MPI init + setup split + io)
  per rank count and refuses a profile without the setup split.

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
runs, the lo and hi step counts (default 100/200), each recorded as its own run
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
import math
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
SUITE_SCHEMA_ID = 'eqdyna-hpc-scaling-suite/2'

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
DEFAULT_NSTEPS_LO = 100   # owner 2026-11-24: 100/200 for the HPC levels; --nsteps-lo/-hi
DEFAULT_NSTEPS_HI = 200
ALLOC_PLACEHOLDER = '<ALLOC>'

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
# helpers: job id, provenance, job script
# --------------------------------------------------------------------------
def parse_job_id(stdout):
    """The job id out of sbatch's stdout. Vista's sbatch prints a multi-line
    banner before (or around) the id, so the whole output is NOT the id (the
    old code stored the banner). The id is the LAST line that is all digits;
    failing that (`--parsable` may print `<id>;<cluster>`, the classic reply
    is `Submitted batch job <id>`), the last all-digit token. Nothing found
    -> raise; the error text never quotes the output (a banner can carry the
    allocation name)."""
    text = stdout.decode() if isinstance(stdout, bytes) else stdout
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    for ln in reversed(lines):
        if ln.isdigit():
            return ln
    tokens = [t for ln in lines for t in ln.replace(';', ' ').split()]
    for t in reversed(tokens):
        if t.isdigit():
            return t
    raise SuiteError('no all-digit job id in sbatch output (%d non-empty line(s)); '
                     'output deliberately not echoed' % len(lines))


def _git(*args):
    return subprocess.check_output(['git', '-C', REPO_ROOT] + list(args)).decode().strip()


def provenance(machine, nsteps_lo, nsteps_hi, repeats, points):
    """The stamp every manifest carries (item 142 (f)). Deliberately has NO
    account/allocation field and no hostname: the manifest ends up in
    docs/evidence. A failing `git` raises -- an unstamped suite is not
    written."""
    return dict(
        git_sha=_git('rev-parse', 'HEAD'),
        git_dirty=bool(_git('status', '--porcelain')),
        machine=machine,
        loaded_modules=[m for m in os.environ.get('LOADEDMODULES', '').split(':') if m],
        nsteps_lo=nsteps_lo, nsteps_hi=nsteps_hi, repeats=repeats,
        rank_layout=[dict(point_id=p['point_id'], ranks=p['ranks'], nnodes=p['nnodes'],
                          decomp=list(p['decomp'])) for p in points],
        date_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))


_MODULE_RE = re.compile(r'^[A-Za-z0-9_.+/-]+$')
_WALLTIME_RE = re.compile(r'^\d{1,3}:[0-5]\d:[0-5]\d$')


def run_name(tag, repeat):
    return '%s_r%d' % (tag, repeat)


def job_commands(point, backend, mpirun, bin_path, work_root, repeats):
    """The run lines of a point's job script: lo then hi, for every repeat."""
    cmds = []
    for k in range(repeats):
        for tag in ('lo', 'hi'):
            rd = '%s/%s/%s' % (work_root.rstrip('/'), point['point_id'], run_name(tag, k))
            if backend == 'fortran':
                cmds.append('cd %s && %s -np %d %s' % (rd, mpirun, point['ranks'], bin_path))
            else:
                cmds.append('cd %s && %s -np %d python3 -m eqdyna --backend jax --mpi'
                            % (rd, mpirun, point['ranks']))
    return cmds


def write_job_script(path, *, jobname, nnodes, ranks, partition, walltime, modules, commands):
    """The sbatch script. It carries NO `#SBATCH -A`: the allocation is passed
    on the sbatch command line only (item 142 (i)), so it can never be written
    into a file that is later committed."""
    for m in modules:
        if not _MODULE_RE.match(m):
            raise SuiteError('module name %r is not a plain module name' % (m,))
    if not _WALLTIME_RE.match(walltime):
        raise SuiteError('walltime %r is not HH:MM:SS' % (walltime,))
    with open(path, 'w') as f:
        f.write('#! /bin/bash\n')
        f.write('#SBATCH -J %s\n#SBATCH -o a.eqdyna.log%%j\n#SBATCH -N %d\n#SBATCH -n %d\n'
                '#SBATCH -p %s\n#SBATCH -t %s\n' % (jobname, nnodes, ranks, partition, walltime))
        f.write('set -euo pipefail\n')
        for m in modules:
            f.write('module load %s\n' % m)
        for c in commands:
            f.write(c + '\n')


# --------------------------------------------------------------------------
# --submit
# --------------------------------------------------------------------------
def submit_point(point, backend, case_name, machine, machine_cfg, account, out_dir, env, opts):
    """Builds the point's lo/hi cases for every repeat and either sbatch's a
    job script (scheduler machine) or runs them in place (no-scheduler
    machine). `account` is used only to build the sbatch command line.
    Returns the point's manifest entry. `opts` carries nsteps_lo/nsteps_hi/
    repeats/modules/mpirun/bin_path/walltime/work_root/dry_run."""
    entry = dict(point)
    entry['backend'] = backend
    entry['runs'] = {'lo': [], 'hi': []}
    work_root = opts['work_root'] or out_dir
    for k in range(opts['repeats']):
        for tag, nsteps in (('lo', opts['nsteps_lo']), ('hi', opts['nsteps_hi'])):
            name = run_name(tag, k)
            run_dir = os.path.join(out_dir, point['point_id'], name)
            if not opts['dry_run']:
                build_point_case(run_dir, case_name, point['dx'], point['decomp'], nsteps, env)
            entry['runs'][tag].append(dict(name=name, run_dir=run_dir, nsteps=nsteps))

    if machine_cfg['scheduler'] == 'slurm':
        mpirun = opts['mpirun'] or machine_cfg['mpirun']
        if mpirun is None:
            raise SuiteError(
                "machine %r has no verified launcher (scripts/machines.py "
                "mpirun=None) -- pass --mpirun explicitly; --submit refuses to "
                "guess one" % (machine,))
        os.makedirs(os.path.join(out_dir, point['point_id']), exist_ok=True)
        script_path = os.path.join(out_dir, point['point_id'], 'submit.sh')
        write_job_script(
            script_path, jobname='hpcscale-%s' % point['point_id'], nnodes=point['nnodes'],
            ranks=point['ranks'], partition=machine_cfg['partition'] or 'normal',
            walltime=opts['walltime'] or machine_cfg['walltime'] or '02:00:00',
            modules=opts['modules'],
            commands=job_commands(point, backend, mpirun,
                                  opts['bin_path'] or os.path.join(REPO_ROOT, 'bin', 'eqdyna'),
                                  work_root, opts['repeats']))
        entry['submit_script'] = script_path
        if opts['dry_run']:
            entry['job_id'] = None
        else:
            out = subprocess.check_output(['sbatch', '--parsable', '-A', account, script_path],
                                          cwd=out_dir, env=env)
            entry['job_id'] = parse_job_id(out)
    else:
        if opts['dry_run']:
            raise SuiteError('--dry-run writes a job script; machine %r has no scheduler' % machine)
        runner = run_point_fortran if backend == 'fortran' else run_point_python_jax_mpi
        for tag in ('lo', 'hi'):
            for run in entry['runs'][tag]:
                runner(run['run_dir'], point['ranks'], env)
        entry['job_id'] = None
    return entry


def _check_steps(lo, hi, repeats):
    if not (isinstance(lo, int) and isinstance(hi, int) and 1 <= lo < hi):
        raise SuiteError('--nsteps-lo %r / --nsteps-hi %r must be ints with 1 <= lo < hi' % (lo, hi))
    if not (isinstance(repeats, int) and repeats >= 1):
        raise SuiteError('--repeats %r must be an int >= 1' % (repeats,))


def cmd_submit(args):
    _validate_combo(args.test, args.level, args.backend)
    _check_steps(args.nsteps_lo, args.nsteps_hi, args.repeats)
    if args.dry_run and args.account is not None:
        raise SuiteError('--dry-run takes no --account: the printed sbatch line carries the '
                         'placeholder %s and the real allocation is typed by hand' % ALLOC_PLACEHOLDER)
    machine_cfg = machines.machine(args.machine)
    max_nodes = args.max_nodes
    if max_nodes is None:
        max_nodes = {'medium': 16, 'large': 64, 'big': 256}[args.level]

    points = build_points(args.test, args.level, machine_cfg, max_nodes)
    out_dir = args.out_dir or os.path.join(
        REPO_ROOT, 'docs', 'evidence',
        'hpc-scaling-%s-%s-%s-%s' % (args.machine, args.backend, args.test, args.level))
    os.makedirs(out_dir, exist_ok=True)

    from testsys.common import mpi_child_env
    env = mpi_child_env(os.environ)
    env['EQDYNAROOT'] = REPO_ROOT
    env['PATH'] = os.pathsep.join([os.path.join(REPO_ROOT, 'bin'),
                                   os.path.join(REPO_ROOT, 'scripts'),
                                   env.get('PATH', '')])
    opts = dict(nsteps_lo=args.nsteps_lo, nsteps_hi=args.nsteps_hi, repeats=args.repeats,
                modules=args.module or [], mpirun=args.mpirun, bin_path=args.bin,
                walltime=args.walltime, work_root=args.work_root, dry_run=args.dry_run)

    entries = [submit_point(p, args.backend, args.case, args.machine, machine_cfg,
                            args.account, out_dir, env, opts) for p in points]
    if args.dry_run:
        for e in entries:
            print('hpc_scaling_suite: DRY RUN (nothing submitted) -- script %s\n'
                  '  sbatch -A %s %s' % (e['submit_script'], ALLOC_PLACEHOLDER, e['submit_script']))
        return 0

    manifest = dict(schema=SUITE_SCHEMA_ID, machine=args.machine, backend=args.backend,
                    test=args.test, level=args.level, case=args.case, max_nodes=max_nodes,
                    nsteps_lo=args.nsteps_lo, nsteps_hi=args.nsteps_hi, repeats=args.repeats,
                    created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                    provenance=provenance(args.machine, args.nsteps_lo, args.nsteps_hi,
                                          args.repeats, points),
                    points=entries)
    with open(os.path.join(out_dir, MANIFEST_NAME), 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print('hpc_scaling_suite: submitted %d point(s) to %s -- manifest at %s'
          % (len(entries), out_dir, os.path.join(out_dir, MANIFEST_NAME)))
    return 0


# --------------------------------------------------------------------------
# parity (recorded at --collect, where the run directories exist)
# --------------------------------------------------------------------------
def _case_gate_dx(case):
    """The dx the case's committed reference was frozen at: `par.dx` in
    case_input/<case>/user_defined_params.py. Raises if it states none."""
    path = os.path.join(REPO_ROOT, 'case_input', case, 'user_defined_params.py')
    with open(path) as f:
        m = re.search(r'^par\.dx\s*=\s*([0-9.eE+-]+)', f.read(), flags=re.M)
    if not m:
        raise SuiteError('%s states no par.dx; the gate resolution is unknown' % path)
    return float(m.group(1))


def parity_record(case, backend, run_dir, dx):
    """{'status': pass|fail|declared, 'lines': [...]} for one timed run.
    A timed run is comparable to the case's reference only if it ran at the
    reference's term (matrix.gate_term_for) AND mesh (the case's own par.dx);
    otherwise it is DECLARED with the measured mismatch, never skipped. A
    comparable run that compare_cell fails (or raises on) is `fail`."""
    mbackend = 'python-jax' if backend == 'python-jax-mpi' else backend
    if not matrix.is_supported(case, mbackend):
        return dict(status='declared', lines=['parity: %s x %s is not a gated cell in '
                                              'matrix.py' % (case, mbackend)])
    term, dt = (float(x) for x in
                open(os.path.join(run_dir, 'bGlobal.txt')).read().splitlines()[15:17])
    gate_term, gate_dx = matrix.gate_term_for(case), _case_gate_dx(case)
    reasons = []
    if not math.isclose(term, gate_term, rel_tol=1e-6):
        reasons.append('timed term %.6g s (dt %.6g) != reference term %.6g s' % (term, dt, gate_term))
    if not math.isclose(dx, gate_dx, rel_tol=1e-9):
        reasons.append('timed dx %.6g m != reference dx %.6g m' % (dx, gate_dx))
    if reasons:
        return dict(status='declared',
                    lines=['parity: not comparable to the case reference: ' + '; '.join(reasons)])
    try:
        ok, lines = compare.compare_cell(case, mbackend, run_dir)
    except Exception as e:  # compare_cell raises on a missing artifact: itself a finding
        return dict(status='fail', lines=['parity: compare_cell raised on %s x %s at %s: %r'
                                          % (case, mbackend, run_dir, e)])
    return dict(status='pass' if ok else 'fail', lines=list(lines))


# --------------------------------------------------------------------------
# --collect
# --------------------------------------------------------------------------
def _runs(point):
    """[(tag, run)] over every repeat of a manifest point (lo and hi lists)."""
    return [(tag, run) for tag in ('lo', 'hi') for run in point['runs'][tag]]


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
    if 'account' in manifest:
        raise SuiteError('manifest carries an `account` field; it must never be written')

    validated = 0
    for point in manifest['points']:
        for tag, run in _runs(point):
            rows = profile_schema.validate_run_dir(run['run_dir'])  # raises on any problem
            if rows[0]['backend'] != manifest['backend']:
                raise SuiteError(
                    'point %s run %s: profile files report backend=%r but '
                    'manifest says backend=%r' % (point['point_id'], run['name'],
                                                  rows[0]['backend'], manifest['backend']))
            validated += len(rows)
            if tag == 'hi':
                run['parity'] = parity_record(manifest['case'], manifest['backend'],
                                              run['run_dir'], point['dx'])
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2, sort_keys=True)

    out_tgz = args.out or os.path.join(
        suite_dir, 'hpc_scaling_%s_%s_%s_%s.tgz'
        % (manifest['machine'], manifest['backend'], manifest['test'], manifest['level']))
    with tarfile.open(out_tgz, 'w:gz') as tar:
        tar.add(manifest_path, arcname=MANIFEST_NAME)
        for point in manifest['points']:
            for tag, run in _runs(point):
                for name in sorted(os.listdir(run['run_dir'])):
                    if profile_schema._RANK_FILE_RE.match(name):
                        arcname = os.path.join(point['point_id'], run['name'], name)
                        tar.add(os.path.join(run['run_dir'], name), arcname=arcname)
    print('hpc_scaling_suite: collected %d profile file(s) from %d point(s) into %s'
          % (validated, len(manifest['points']), out_tgz))
    return 0


# --------------------------------------------------------------------------
# --analyze
# --------------------------------------------------------------------------
def _load_run_profiles(tmpdir, point, tag):
    """[rows by rank] for every repeat of `tag`, from the unpacked tarball."""
    return [profile_schema.validate_run_dir(os.path.join(tmpdir, point['point_id'], run['name']))
            for run in point['runs'][tag]]


# profile_emit.py's buckets split into two KINDS of cost. PER_STEP_BUCKETS
# accumulate once per solver step, so per_step_and_fixed's differencing is the
# right arithmetic for them; FIXED_BUCKETS (setup, io) are a ONE-TIME cost
# each run pays regardless of nsteps -- differencing them is the WRONG
# QUESTION (lo and hi pay the same, the difference is noise around zero and
# per_step_and_fixed correctly refuses it). Fixed buckets are reported
# directly from the HI run.
PER_STEP_BUCKETS = ('element', 'fault', 'exchange', 'wait')
FIXED_BUCKETS = ('setup', 'io')
assert set(PER_STEP_BUCKETS) | set(FIXED_BUCKETS) == set(profile_schema.BUCKET_KEYS)

# The `setup_detail_s` keys each backend's profile must carry (item 142 (b),
# (c)). Fortran times all five; python-jax-mpi fuses input read + meshing in
# build_solver_state and has no netCDF open, and its MPI is initialised at
# mpi4py import, outside any timer -- so it carries two keys and `mpi_init`
# is DECLARED not measured (overhead.init is null), never zero.
DETAIL_KEYS = {
    'fortran': ('mpi_init', 'input_read', 'meshing', 'netcdf_open', 'mass_assembly'),
    'python-jax-mpi': ('input_mesh', 'setup_exchange'),
}
DETAIL_SETUP_FLOOR_S = 0.05
DETAIL_SETUP_REL = 0.05


def _stats(vals):
    return dict(min=min(vals), mean=sum(vals) / len(vals), max=max(vals))


def per_rank_metrics(lo_rows, hi_rows):
    """Per-rank metrics for one (lo, hi) pair of runs. Per-step cost is the
    difference of each rank's OWN `loop_s` over the two step counts (never a
    wrapper clock, never total_s: CLAUDE.md "per-step-by-difference breaks at
    >=4 ranks"); run_numa_scaling.per_step_and_fixed raises on a non-positive
    result. Overhead = MPI init (if the backend times it) + setup + io of the
    hi run, with the setup split (`setup_detail_s`) carried per rank. A
    profile with no/incomplete `setup_detail_s` raises: it came from a solver
    that cannot give the split, and a report without it is not the report the
    owner asked for."""
    if len(lo_rows) != len(hi_rows):
        raise SuiteError('lo/hi profile rank counts differ: %d vs %d'
                         % (len(lo_rows), len(hi_rows)))
    per_rank = []
    for lo, hi in zip(lo_rows, hi_rows):
        if lo['rank'] != hi['rank']:
            raise SuiteError('lo/hi rank mismatch: %r vs %r' % (lo['rank'], hi['rank']))
        n_lo, n_hi = lo['nsteps'], hi['nsteps']
        ps_total, _ = numa.per_step_and_fixed(lo['loop_s'], hi['loop_s'], n_lo, n_hi)
        buckets_per_step = {}
        for key in PER_STEP_BUCKETS:
            ps, _fixed = numa.per_step_and_fixed(
                lo['buckets_s'][key], hi['buckets_s'][key], n_lo, n_hi)
            buckets_per_step[key] = ps
        buckets_fixed = {key: hi['buckets_s'][key] for key in FIXED_BUCKETS}

        backend = hi['backend']
        if backend not in DETAIL_KEYS:
            raise SuiteError('backend %r has no declared setup_detail_s keys' % (backend,))
        detail = hi.get('setup_detail_s')
        want = DETAIL_KEYS[backend]
        if not isinstance(detail, dict) or set(detail) != set(want):
            raise SuiteError(
                'rank %d (%s): setup_detail_s %r != required keys %s -- the profile '
                'comes from a solver build without the setup split'
                % (hi['rank'], backend, None if detail is None else sorted(detail), list(want)))
        in_setup = sum(v for k, v in detail.items() if k != 'mpi_init')
        setup = hi['buckets_s']['setup']
        slack = max(DETAIL_SETUP_REL * setup, DETAIL_SETUP_FLOOR_S)
        if not (in_setup <= setup + slack) or (backend == 'fortran' and not (in_setup >= setup - slack)):
            raise SuiteError(
                'rank %d (%s): setup_detail_s parts sum to %.6g s but buckets_s.setup is '
                '%.6g s (allowed slack %.3g s) -- the split does not account for setup'
                % (hi['rank'], backend, in_setup, setup, slack))
        init = detail.get('mpi_init')
        overhead = dict(init_s=init, setup_s=setup, io_s=hi['buckets_s']['io'],
                        total_s=(init or 0.0) + setup + hi['buckets_s']['io'])
        per_rank.append(dict(rank=lo['rank'], host=lo['host'], per_step_s=ps_total,
                             buckets_per_step_s=buckets_per_step,
                             buckets_fixed_s=buckets_fixed,
                             setup_detail_s=dict(detail), overhead=overhead))
    return per_rank


def imbalance(per_rank):
    vals = [r['per_step_s'] for r in per_rank]
    mean = sum(vals) / len(vals)
    return dict(max_s=max(vals), mean_s=mean, max_over_mean=max(vals) / mean)


def point_parity(point):
    """(status, lines) over every repeat's hi run, read from the manifest
    (recorded at --collect). A run with no record raises: an unchecked timed
    point is not reported as anything."""
    recs = []
    for run in point['runs']['hi']:
        if 'parity' not in run:
            raise SuiteError('point %s run %s carries no parity record -- re-run --collect'
                             % (point['point_id'], run['name']))
        recs.append(run['parity'])
    statuses = {r['status'] for r in recs}
    if not statuses <= {'pass', 'fail', 'declared'}:
        raise SuiteError('point %s: unknown parity status %s' % (point['point_id'], sorted(statuses)))
    status = 'fail' if 'fail' in statuses else ('declared' if 'declared' in statuses else 'pass')
    lines = []
    for r in recs:
        for ln in r['lines']:
            if ln not in lines:
                lines.append(ln)
    return status, lines


def _mean(vals):
    return sum(vals) / len(vals)


def analyze(tgz_path):
    """The report dict --analyze prints/returns. Pure function over the
    unpacked tarball's manifest+profiles, directly unit-testable."""
    tmpdir = tempfile.mkdtemp(prefix='hpc_scaling_analyze_')
    try:
        with tarfile.open(tgz_path, 'r:gz') as tar:
            tar.extractall(tmpdir)
        with open(os.path.join(tmpdir, MANIFEST_NAME)) as f:
            manifest = json.load(f)
        if manifest.get('schema') != SUITE_SCHEMA_ID:
            raise SuiteError('manifest schema %r != %r' % (manifest.get('schema'), SUITE_SCHEMA_ID))
        if 'provenance' not in manifest:
            raise SuiteError('manifest has no provenance stamp')

        report = dict(schema='eqdyna-hpc-scaling-report/2', machine=manifest['machine'],
                      backend=manifest['backend'], test=manifest['test'],
                      level=manifest['level'], case=manifest['case'],
                      provenance=manifest['provenance'], points=[])

        for point in manifest['points']:
            lo_runs = _load_run_profiles(tmpdir, point, 'lo')
            hi_runs = _load_run_profiles(tmpdir, point, 'hi')
            if len(lo_runs) != len(hi_runs):
                raise SuiteError('point %s: %d lo repeat(s) vs %d hi'
                                 % (point['point_id'], len(lo_runs), len(hi_runs)))
            reps = [per_rank_metrics(lo, hi) for lo, hi in zip(lo_runs, hi_runs)]
            imbs = [imbalance(pr) for pr in reps]
            rep_means = [i['mean_s'] for i in imbs]
            parity_status, parity_lines = point_parity(point)
            samples = [r['overhead'] for pr in reps for r in pr]
            detail_keys = list(reps[0][0]['setup_detail_s'])
            report['points'].append(dict(
                point_id=point['point_id'], ranks=point['ranks'], nnodes=point['nnodes'],
                dx=point['dx'], decomp=point['decomp'], repeats=len(reps),
                mean_per_step_s=_mean(rep_means),
                max_per_step_s=max(i['max_s'] for i in imbs),
                max_over_mean=_mean([i['max_over_mean'] for i in imbs]),
                repeat_per_step_s=rep_means,
                repeat_rel_range=(None if len(reps) == 1 else
                                  (max(rep_means) - min(rep_means)) / _mean(rep_means)),
                overhead_s=dict(
                    init=(None if samples[0]['init_s'] is None
                          else _stats([s['init_s'] for s in samples])),
                    setup=_stats([s['setup_s'] for s in samples]),
                    io=_stats([s['io_s'] for s in samples]),
                    total=_stats([s['total_s'] for s in samples])),
                setup_detail_mean_s={k: _mean([r['setup_detail_s'][k] for pr in reps for r in pr])
                                     for k in detail_keys},
                bucket_breakdown_per_step_s={
                    key: _mean([r['buckets_per_step_s'][key] for pr in reps for r in pr])
                    for key in PER_STEP_BUCKETS},
                parity_status=parity_status, parity_lines=parity_lines,
            ))

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
    p.add_argument('--nsteps-lo', type=int, default=DEFAULT_NSTEPS_LO)
    p.add_argument('--nsteps-hi', type=int, default=DEFAULT_NSTEPS_HI)
    p.add_argument('--repeats', type=int, default=1)
    p.add_argument('--dry-run', action='store_true',
                   help='write the job script(s) only; build no case, submit nothing')
    p.add_argument('--work-root', help='run-directory root as the JOB sees it (may contain $SCRATCH)')
    p.add_argument('--bin', help='solver binary as the JOB sees it (may contain $SCRATCH)')
    p.add_argument('--walltime', help='HH:MM:SS, overrides the machine default')
    p.add_argument('--module', action='append', help='`module load <name>` in the job (repeatable)')
    p.add_argument('--mpirun', help='launcher, overrides an unverified machine entry')
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.submit:
            required = [('--machine', args.machine), ('--backend', args.backend),
                        ('--test', args.test), ('--level', args.level)]
            if not args.dry_run:
                required.append(('--account', args.account))
            for name, val in required:
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
