#! /usr/bin/env python3
"""
Single entry point for EQdyna's tiered test system (PROJECT_RULES.md rule 3).

    python3 testsys/run.py unit          # fast pure-python unit tests (pytest, no MPI/Fortran)
    python3 testsys/run.py regression    # one guard per past incident (rule 10)
    python3 testsys/run.py e2e           # THE test: the everyday case x backend sweep vs test.reference.results/, at matrix.GATE_TERM_S (rule 7)
    python3 testsys/run.py e2e-ci        # the same sweep, restricted to matrix.CI_CELLS (2026-09-23: one case x 2 backends -- fortran, python-jax -- portability smoke)
    python3 testsys/run.py release       # THE test, over every supported cell (currently the same cells as `e2e` -- matrix.RELEASE_ONLY was retired 2026-09-23): the release gate, writes docs/evidence/sweep-<sha>/summary.json
    python3 testsys/run.py e2e-full      # SCEC cases at spec dx/term, 16 ranks, report-only (opt-in; hours)
    python3 testsys/run.py gpu           # the sweep's python-jax column on CUDA (one cell; needs a GPU)
    python3 testsys/run.py perf          # pinned single-core Fortran/NumPy/JAX timing, ratio-guarded
    python3 testsys/run.py profile-overhead   # EQDYNA_PROFILE on/off A/B, all 4 backends (needs EQDYNA_PROFILE_OVERHEAD_CPUS)
    python3 testsys/run.py readme        # the README-executes stranger-clone gate, FULL mode: clones this commit, runs README.md's fenced blocks for real (~2-4 min)
    python3 testsys/run.py all           # unit + regression + e2e, in that order (default; perf is opt-in, not in "all" -- it needs a Fortran build a fresh checkout does not have yet)

    python3 testsys/run.py all --machine ls6              # sets EQDYNA_TEST_MACHINE/EQDYNA_MPIRUN
                                                            # from scripts/machines.py and runs IN
                                                            # PLACE (e.g. inside an interactive
                                                            # allocation) -- ubuntu behaviour is
                                                            # unchanged when --machine is absent.
    python3 testsys/run.py all --machine ls6 --submit --account <alloc>
                                          # writes an sbatch script -- case.setup's shared SBATCH-header
                                          # writer (scripts/machines.write_slurm_header), body
                                          # `source install-eqdyna.sh -c ls6` + this same `run.py all
                                          # --machine ls6` -- submits it, and prints the job id and the
                                          # tarball name its own run.py invocation will produce. Refuses
                                          # on a machine with no scheduler (scripts/machines.py).

There is ONE test here -- e2e -- and backend is an axis of it, not a tier.
There is also ONE term (2026-09-23 owner decision, superseding a same-day
earlier two-term design): every cell of every selection -- `e2e`, `e2e-ci`,
`release` alike -- runs at matrix.GATE_TERM_S (5 s), applied through the ONE
existing override path regardless of a case's own committed par.term.
`release` is run_e2e.py's `--release` selector: it is the tier that gates a
release (PROJECT_RULES rule 15/16) and writes docs/evidence/sweep-<sha>/
summary.json (rule 24), currently over the SAME cell set as `e2e` --
matrix.RELEASE_ONLY, the mechanism that used to widen it, was retired
2026-09-23 (its only two occupants were both python-numpy cells, removed
along with the backend axis itself). CI runs e2e-ci ONLY; `release` is a
human- or conductor-scheduled tier, the same as the old `run.py all`'s e2e
slice was.

Two tiers that used to sit beside it are gone, for the same reason each time:
they were a SECOND way of asking a question e2e already answers, and a second
way of asking costs a second implementation to keep honest.
  * `accept` -- the same solver against the same references with its own
    comparison implementation over a shorter case list.
  * `parity` -- Python vs a Fortran serial oracle, per step, via the
    eqdyna-pydump build. It answered "at WHICH STEP does it diverge?", which
    is a debugging question, not a gating one; it root-caused the missing
    rsfNucleation branch and that work is done. It went stale the moment the
    Python tree was restructured -- its two test files were named for a
    `standalone/` package that no longer exists -- and keeping a stale tier
    green costs more than re-deriving a step dump the next time one is needed.
    Removed with it: make_fixtures.py, the pydump fixtures, and src/pydump*
    (the Fortran link-time seam they existed to drive).

Prints a per-test SUCCESS/FAIL line (from pytest or from each regression
script's own banner), a per-tier SUMMARY line, and exits non-zero if
anything in the requested scope failed. "The script ran" and "the script
passed" are always two different questions here (rule 3).
"""
import glob
import io
import os
import re
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)


def _load_by_path(name, relpath):
    """Load a scripts/ module by path, so putting scripts/ on sys.path
    cannot shadow any other module (same reasoning as the old
    _load_src_hash, generalised to the second caller: scripts/machines.py)."""
    import importlib.util
    path = os.path.join(REPO_ROOT, relpath)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_src_hash():
    """scripts/src_hash.py -- the ONE source-stamp implementation."""
    return _load_by_path('src_hash', os.path.join('scripts', 'src_hash.py'))


def _load_machines():
    """scripts/machines.py -- the ONE HPC-machine registry and the shared
    SBATCH-header writer case.setup also uses."""
    return _load_by_path('machines', os.path.join('scripts', 'machines.py'))


def build_submit_script(machines_mod, machine_name, tiers, account,
                         partition=None, walltime=None):
    """The sbatch script for `--submit`: case.setup's shared SBATCH header
    (scripts/machines.write_slurm_header) followed by a body that sources
    this machine's install-eqdyna.sh branch (the build/env source of truth --
    module loads and venv setup are NOT duplicated here) and re-invokes this
    exact `run.py <tiers> --machine <machine_name>`, then packs the run's log
    and the perf/profile rows it appended into one tarball. Replaces
    scripts/ls6_sweep.sbatch's hard-coded, LS6-only copy of this same shape.

    Raises ValueError -- never guesses -- when the machine has no scheduler,
    no account is given (this charges someone's allocation; no default -- rule
    2), or the registry is missing a value (e.g. grace's cores_per_node) that
    --partition/--time cannot supply."""
    m = machines_mod.machine(machine_name)
    if m['scheduler'] != 'slurm':
        raise ValueError(
            "machine %r has no scheduler (scripts/machines.py: scheduler=%r) "
            "-- --submit only applies to a slurm machine; drop --submit and "
            "run `testsys/run.py %s --machine %s` in place"
            % (machine_name, m['scheduler'], ' '.join(tiers), machine_name))
    if not account:
        raise ValueError(
            "--submit needs --account <allocation> -- an sbatch job charges "
            "an allocation and this refuses to guess one (PROJECT_RULES rule 2)")
    queue = partition or m['partition']
    if not queue:
        raise ValueError(
            "machine %r has no default partition in scripts/machines.py -- "
            "pass --partition explicitly" % machine_name)
    wt = walltime or m['walltime']
    if not wt:
        raise ValueError(
            "machine %r has no default walltime in scripts/machines.py -- "
            "pass --time explicitly" % machine_name)
    if m['cores_per_node'] is None:
        raise ValueError(
            "machine %r has no cores_per_node in scripts/machines.py yet -- "
            "fill it in before submitting a job there (rule 2: no guessed "
            "core count)" % machine_name)
    if not m['mpirun']:
        raise ValueError(
            "machine %r has no verified mpirun launcher in scripts/machines.py "
            "yet -- fill it in (see ls6's ibrun-vs-mpirun note) before "
            "submitting a job there (rule 2: no guessed launcher)" % machine_name)

    tier_str = ' '.join(tiers)
    buf = io.StringIO()
    machines_mod.write_slurm_header(
        buf, jobname='eqdyna-sweep', nnode=1, ncpu=m['cores_per_node'],
        queue=queue, walltime=wt, account=account, email='',
        output='eqdyna_sweep_%j.log')
    buf.write('\n')
    buf.write('cd "$SLURM_SUBMIT_DIR"\n')
    buf.write('source ./install-eqdyna.sh -c %s\n' % machine_name)
    buf.write('git rev-parse HEAD; module list 2>&1\n')
    # Refuse EARLY, before the (potentially hours-long) sweep, if the
    # sourced environment is missing a dependency the sweep needs -- this is
    # exactly what bit the owner on 2026-09-29: a bare `python3` with none of
    # jax/netCDF4/mpi4py, discovered only mid-regression via "No module named
    # netCDF4". Explicit check, not `set -e` (rule 2: no silent fallback).
    env_check_msg = (
        "eqdyna-sweep: environment check failed after "
        "source ./install-eqdyna.sh -c %s -- jax/netCDF4/mpi4py not all "
        "importable by this python3. Check that the venv activated (see "
        "install-eqdyna.sh own warning above if it did not) and EQDYNA_VENV."
        % machine_name)
    buf.write("python3 -c 'import jax, netCDF4, mpi4py' || {\n")
    buf.write('  echo "%s" >&2\n' % env_check_msg)
    buf.write('  exit 1\n')
    buf.write('}\n')
    buf.write('L0=$(wc -l < docs/perf_ledger.jsonl); '
              'P0=$(wc -l < docs/run_profiles.jsonl)\n')
    buf.write('python3 testsys/run.py %s --machine %s\n'
              % (tier_str, machine_name))
    buf.write('rc=$?\n')
    buf.write('echo "run.py %s --machine %s exit $rc"\n'
              % (tier_str, machine_name))
    buf.write('tail -n +$((L0 + 1)) docs/perf_ledger.jsonl '
              '> eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl\n')
    buf.write('tail -n +$((P0 + 1)) docs/run_profiles.jsonl '
              '> eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl\n')
    buf.write('tar czf eqdyna_sweep_$SLURM_JOB_ID.tgz '
              'eqdyna_sweep_$SLURM_JOB_ID.log '
              'eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl '
              'eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl\n')
    buf.write('exit $rc\n')
    return buf.getvalue()


def submit(machine_name, tiers, account, partition=None, walltime=None):
    """Write the submit script (under scratch/, gitignored -- it is a
    generated artifact, not source, and SLURM_SUBMIT_DIR is set from sbatch's
    OWN cwd, not the script's location, so it need not live at repo root),
    then `sbatch` it. Returns (returncode, stdout, stderr, script_path)."""
    machines_mod = _load_machines()
    script = build_submit_script(machines_mod, machine_name, tiers, account,
                                  partition, walltime)
    scratch = os.path.join(REPO_ROOT, 'scratch')
    os.makedirs(scratch, exist_ok=True)
    script_path = os.path.join(scratch, 'submit_eqdyna_sweep_%s.sh' % machine_name)
    with open(script_path, 'w') as f:
        f.write(script)
    os.chmod(script_path, 0o755)
    r = subprocess.run(['sbatch', script_path], cwd=REPO_ROOT,
                       capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr, script_path


def require_fresh_fortran_binary(binary=None, fsrc=None):
    """The regression tier drives bin/eqdyna directly; refuse (return the
    reason) unless it was built from THIS tree's src/fortran (owner-approved
    2026-09-24: a Sep-22 binary made test_profile_env_strict fail for no
    real reason, and a stale one could as easily pass on old code). A
    missing binary and a missing stamp refuse too. No bypass. Returns None
    when the binary is fresh."""
    src_hash = _load_src_hash()
    binary = binary or os.path.join(REPO_ROOT, 'bin', 'eqdyna')
    ok, msg = src_hash.check_binary(binary, fsrc or src_hash.FSRC)
    print(('testsys: ' if ok else 'testsys: REFUSED - ') + msg)
    return None if ok else msg
TIERS = ('unit', 'regression', 'e2e')

# Line-buffered: redirected to a log (CI pipes this through `tee`), Python
# block-buffers this script's own prints while the children write straight to
# the fd, so the tier banners and the SUMMARY would land out of order, and
# would vanish entirely if the job were killed.
sys.stdout.reconfigure(line_buffering=True)


def run_unit():
    d = os.path.join(TESTSYS, 'unit')
    print('\n==== testsys: unit ====')
    return subprocess.call([sys.executable, '-m', 'pytest', '-v', d], cwd=REPO_ROOT)


def run_regression():
    """Each regression/test_*.py is a standalone script with its own
    SUCCESS/FAIL banner and sys.exit (matching the shape of the incident
    it guards) rather than pytest functions -- run each one and gate on
    its exit code."""
    print('\n==== testsys: regression ====')
    scripts = sorted(glob.glob(os.path.join(TESTSYS, 'regression', 'test_*.py')))
    if not scripts:
        print('regression: FAIL - no regression scripts found (misconfigured testsys/)')
        return 1
    overall = 0
    for script in scripts:
        name = os.path.basename(script)
        print(f'-- {name} --')
        rc = subprocess.call([sys.executable, script], cwd=REPO_ROOT)
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {name} (exit {rc})')
        overall = overall or rc
    return overall


def _e2e(*args):
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py')]
                           + list(args), cwd=REPO_ROOT)


def run_e2e():
    """THE sweep, EVERYDAY selection (run_e2e.py's default, no flags): every
    supported case x backend cell, at matrix.GATE_TERM_S regardless of its
    own committed par.term. (matrix.RELEASE_ONLY, which used to hold cells
    back from this selection for cost, was retired 2026-09-23.)"""
    print('\n==== testsys: e2e ====')
    return _e2e()


def run_e2e_ci():
    """The sweep restricted to matrix.CI_CELLS -- 2026-09-23 (owner-approved
    test-methodology change): one case (test.tpv8), both remaining backend
    implementations (fortran, python-jax), at the gate term. Its job is
    portability (clean checkout, fresh deps, a different MPI from this box's),
    not physics coverage -- that job belongs to `e2e` and `release` now. This
    is a declared selection, not a magic env var: the run prints which cells
    it covered and which it did not."""
    print('\n==== testsys: e2e-ci ====')
    return _e2e('--ci')


def run_release():
    """THE sweep, RELEASE selection (run_e2e.py --release): every supported
    cell, at the SAME matrix.GATE_TERM_S as every other selection. This is the
    release gate (PROJECT_RULES rule 15/16) -- the tier `run.py all`'s e2e
    slice used to be before CI stopped running the wider sweep. Currently
    selects the SAME cells as `e2e` (matrix.RELEASE_ONLY, the mechanism that
    used to widen this selection beyond the everyday one, was retired
    2026-09-23 along with the python-numpy backend axis -- its only two
    occupants were both python-numpy cells); kept as its own flag because it
    is also the deliberate pre-tag/evidence invocation, not merely a cell
    filter. Run over the full default selection (no --cases/--backends), it
    also writes docs/evidence/sweep-<shortsha>/summary.json
    (run_e2e.write_release_evidence).

    `run.py unit regression release` in ONE invocation is SAFE: checked
    2026-09-23 (rule-24 tree_clean fix) that no regression script under
    testsys/regression/ writes to a real tracked repo path -- every script
    that touches docs/perf_ledger.jsonl-shaped data does so inside its own
    tempfile.mkdtemp() sandbox or reads the committed ledger read-only
    (test_perf_ledger.py); test_sweep_speed_2026_09_23.py imports run_e2e
    but calls schedule_order() directly, never main(), so it runs no sweep.
    run_e2e.py's tree_clean is captured at the very start of ITS OWN
    process (run_e2e.capture_start_tree_state, called first thing in
    main()), so even if a future regression script broke this invariant by
    writing to a real tracked path, `release`'s reading would still be
    correct -- it would just, correctly, read as dirty. This tier is not
    required to run standalone; nothing before it in the same invocation
    is known to dirty the tree, and if that ever changes the release
    evidence will say so rather than silently passing."""
    print('\n==== testsys: release ====')
    return _e2e('--release')


def run_gpu():
    """The sweep's python-jax column on CUDA. A selection of the one sweep,
    not a tier of its own: JAX_PLATFORMS=cuda, same cases, same references,
    same per-case bounds. Fails loudly if there is no CUDA device -- "no GPU
    here" and "the GPU column passed" must not share an exit code (rule 2)."""
    print('\n==== testsys: gpu ====')
    try:
        import jax
    except Exception as exc:                        # noqa: BLE001
        print('FAIL gpu: jax is not importable (%s) -- this selection cannot '
              'be gated on this machine' % exc)
        return 1
    devs = [d for d in jax.devices() if 'cuda' in str(d).lower() or 'gpu' in str(d).lower()]
    if not devs:
        print('FAIL gpu: no CUDA device/plugin visible to jax (pip install '
              '"jax[cuda12]" on a GPU box). Reporting this as a failure, not a '
              'skip: a green line here would claim GPU coverage that did not '
              'happen.')
        return 1
    return _e2e('--backends', 'python-jax', '--cases', 'test.tpv8',
                '--device', 'cuda')


def run_readme():
    """The stranger-clone gate, FULL mode (testsys/regression/
    test_readme_executes.py, same file the seconds-scale regression tier
    runs in its FAST mode -- one parser/engine, two depths). Sets
    EQDYNA_README_GATE=full so that script clones THIS commit into a temp
    dir and actually executes README.md's Requirements/Install/Quick
    start/Python-solver fenced blocks in a bare env -i-equivalent
    environment, then checks the documented output files exist and are
    non-empty. ~2-4 minutes -- opt-in, never swept into `all`, for the same
    reason e2e-full and perf are not (rule 9: a fresh checkout's `run.py
    regression` must stay seconds-scale)."""
    print('\n==== testsys: readme ====')
    env = dict(os.environ)
    env['EQDYNA_README_GATE'] = 'full'
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'regression', 'test_readme_executes.py')],
        cwd=REPO_ROOT, env=env)


def run_e2e_full():
    print('\n==== testsys: e2e-full ====')
    return subprocess.call(
        [sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py'), '--full'],
        cwd=REPO_ROOT)


def run_scaling():
    print('\n==== testsys: scaling ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_scaling.py')], cwd=REPO_ROOT)


def run_perf():
    print('\n==== testsys: perf ====')
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'perf', 'run_perf.py')],
                            cwd=REPO_ROOT)


def run_profile_overhead():
    """Zero-perf-cost gate for the always-on per-rank profiler
    (EQDYNA_PROFILE=1 vs 0): testsys/perf/profile_overhead.py, once per
    backend (python-numpy, python-jax, fortran, python-jax-mpi), NEVER
    concurrently -- each subprocess.call blocks until the previous one has
    fully exited, one measurement in flight at a time on this box.

    Opt-in, same shape as `perf`/`scaling`: requires `EQDYNA_PROFILE_OVERHEAD_CPUS`
    (comma-separated cpu list, e.g. "16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31")
    naming the cpus to pin every arm to -- refused, not defaulted, if unset:
    an overhead measurement on a placement nobody named is not reproducible
    (profile_overhead.py's own module docstring). `EQDYNA_PROFILE_OVERHEAD_CASE`
    defaults to test.tpv8 (the case this gate's owner measured against).
    `EQDYNA_PROFILE_OVERHEAD_RANKS` defaults to 4 for the two MPI backends."""
    print('\n==== testsys: profile-overhead ====')
    cpus = os.environ.get('EQDYNA_PROFILE_OVERHEAD_CPUS')
    if not cpus:
        print('profile-overhead: FAIL - EQDYNA_PROFILE_OVERHEAD_CPUS is not '
              'set. This gate refuses to run unpinned; export the cpu list '
              'this box reserves for it (see this function\'s docstring) and '
              're-run.')
        return 1
    case = os.environ.get('EQDYNA_PROFILE_OVERHEAD_CASE', 'test.tpv8')
    ranks = os.environ.get('EQDYNA_PROFILE_OVERHEAD_RANKS', '4')
    overall = 0
    for backend in ('python-numpy', 'python-jax', 'fortran', 'python-jax-mpi'):
        print(f'-- profile-overhead: {backend} --')
        cmd = [sys.executable, os.path.join(TESTSYS, 'perf', 'profile_overhead.py'),
              '--case', case, '--backend', backend, '--cpus', cpus]
        if backend in ('fortran', 'python-jax-mpi'):
            cmd += ['--ranks', ranks]
        rc = subprocess.call(cmd, cwd=REPO_ROOT)
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} profile-overhead {backend} (exit {rc})')
        overall = overall or rc
    return overall


RUNNERS = {'unit': run_unit, 'regression': run_regression, 'e2e': run_e2e,
           'e2e-ci': run_e2e_ci, 'release': run_release, 'perf': run_perf,
           'gpu': run_gpu, 'scaling': run_scaling, 'e2e-full': run_e2e_full,
           'profile-overhead': run_profile_overhead, 'readme': run_readme}
# 'all' stays unit+regression+e2e only (TIERS below, `e2e` the everyday
# selection) -- perf requires a Fortran build and a generated baseline that a
# fresh checkout does not have; it is opt-in, invoked by name, not swept into
# 'all'. e2e-ci and gpu are SELECTIONS of e2e, not tiers of their own;
# `release` is the same-cells (2026-09-23: RELEASE_ONLY retired), same-term
# selection of the same one sweep, kept as its own flag for the pre-tag
# evidence artifact it writes, not for a second implementation. 'all'
# runs the everyday sweep; CI runs e2e-ci (2026-09-23: a one-case portability
# smoke, not physics coverage); `release` is the human-/conductor-scheduled
# release gate, whose narrower or wider coverage the sweep itself prints
# either way.
# e2e-full additionally needs EQDYNA_FULL_LAUNCH=yes-hours (see
# testsys/e2e/run_e2e_full.py) -- spec-resolution SCEC runs are hours long
# and user-scheduled, never automatic.
OPTIONAL_TIERS = ('e2e-ci', 'release', 'perf', 'gpu', 'scaling', 'e2e-full',
                  'profile-overhead', 'readme')
# `readme` is the README-executes stranger-clone gate, FULL mode (~2-4 min).
# Requesting `release` also runs it (mission: "include it in run.py
# release") -- a release gate that never re-proves the README a new user
# would actually follow is not the release gate, it is the release gate
# minus the one artifact every new user reads first. See RELEASE_EXPANDS.
RELEASE_EXPANDS = ('release', 'readme')


def parse_argv(argv):
    """Split argv into (tiers, machine, submit, account, partition, walltime).
    Manual, not argparse: the existing contract is a bare list of tier names,
    and this only adds a handful of `--flag value` pairs ahead of or after
    them -- argparse's positional/optional interleaving would change how an
    unknown tier is reported (tested by callers today)."""
    tiers = []
    machine = None
    submit_flag = False
    account = None
    partition = None
    walltime = None
    flags_with_value = {
        '--machine': 'machine', '--account': 'account',
        '--partition': 'partition', '--time': 'walltime',
    }
    i = 1
    while i < len(argv):
        a = argv[i]
        if a in flags_with_value:
            if i + 1 >= len(argv):
                raise SystemExit('testsys/run.py: %s needs a value' % a)
            value = argv[i + 1]
            if a == '--machine':
                machine = value
            elif a == '--account':
                account = value
            elif a == '--partition':
                partition = value
            elif a == '--time':
                walltime = value
            i += 2
        elif a == '--submit':
            submit_flag = True
            i += 1
        else:
            tiers.append(a)
            i += 1
    return tiers, machine, submit_flag, account, partition, walltime


def main(argv):
    # Several tiers may be named in one invocation, in the order given --
    # that is how CI asks for "unit regression e2e-ci" without a second
    # entry point and without an env var deciding coverage behind its back.
    tiers, machine_name, submit_flag, account, partition, walltime = parse_argv(argv)
    requested = tiers or ['all']
    unknown = [t for t in requested if t not in TIERS + OPTIONAL_TIERS + ('all',)]
    if unknown:
        print(f'unknown tier(s) {unknown}')
        print(f'usage: python3 testsys/run.py [{"|".join(TIERS + OPTIONAL_TIERS)}|all] '
              f'[--machine <m>] [--submit --account <alloc> [--partition <p>] [--time <hh:mm:ss>]]')
        return 2

    if submit_flag:
        if not machine_name:
            print('testsys/run.py: --submit requires --machine <m>')
            return 2
        try:
            rc, out, err, script_path = submit(machine_name, requested, account,
                                                partition, walltime)
        except ValueError as exc:
            print('testsys/run.py: %s' % exc)
            return 2
        sys.stdout.write(out)
        if rc != 0:
            print('testsys/run.py: sbatch failed (exit %d): %s' % (rc, err.strip()))
            return 1
        jobid_match = re.search(r'(\d+)\s*$', out.strip())
        jobid = jobid_match.group(1) if jobid_match else '<jobid>'
        print('testsys/run.py: submitted %s (job %s); results will land in '
              'eqdyna_sweep_%s.tgz (see %s)' % (script_path, jobid, jobid, script_path))
        return 0

    if machine_name:
        machines_mod = _load_machines()
        try:
            m = machines_mod.machine(machine_name)
        except ValueError as exc:
            print('testsys/run.py: %s' % exc)
            return 2
        if not m['mpirun']:
            print('testsys/run.py: machine %r has no verified mpirun launcher '
                  'in scripts/machines.py yet -- refusing to guess one '
                  '(rule 2); fill it in first' % machine_name)
            return 2
        os.environ['EQDYNA_TEST_MACHINE'] = machine_name
        os.environ['EQDYNA_MPIRUN'] = m['mpirun']
        print('testsys/run.py: --machine %s -> EQDYNA_TEST_MACHINE=%s EQDYNA_MPIRUN=%s'
              % (machine_name, machine_name, m['mpirun']))

    selected = []
    for tier in requested:
        if tier == 'all':
            expansion = TIERS
        elif tier == 'release':
            expansion = RELEASE_EXPANDS
        else:
            expansion = (tier,)
        for t in expansion:
            if t not in selected:
                selected.append(t)

    if 'regression' in selected and require_fresh_fortran_binary() is not None:
        return 1

    results = {tier: RUNNERS[tier]() for tier in selected}

    print('\n==== testsys: SUMMARY ====')
    print('tiers run: %s' % ', '.join(selected))
    overall = 0
    for tier in selected:
        rc = results[tier]
        print(f'{"SUCCESS" if rc == 0 else "FAIL"} {tier} (exit {rc})')
        overall = overall or rc

    return 1 if overall else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
