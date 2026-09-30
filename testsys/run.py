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

LOGS LIVE WITH THE RUN, the same place on every machine (owner, 2026-09-29):
every invocation's own full console output -- including every subprocess's
output, via a real `tee` -- is ALSO written under `<repo>/test/`. Rotation
(`test/` -> `test.prev/`, rule 8) happens exactly where it always did --
inside testsys/e2e/run_e2e.py, immediately before the e2e cells actually
run, once every one of ITS OWN refusal gates (create.newcase, the build, a
stale binary stamp) has passed (PR #56 audit, blocker B1: hoisting the
rotation any earlier meant a refused or Ctrl-C'd run -- a stale binary, a
bad build, a SECOND run.py's own require_fresh_fortran_binary check -- still
destroyed the PREVIOUS real sweep's evidence). A selection touching the e2e
run tree (`e2e`, `e2e-ci`, `release`, `gpu` -- never more than one of them in
the same invocation, which would rotate twice and destroy the first result)
stages its console log OUTSIDE test/ (`scratch/run.<pid>.log`); the
exclusive test/ lock stays entirely inside run_e2e.py's own Gate 0, taken
unconditionally every time it runs (a same-day follow-up finding: holding
it any earlier, from run.py itself, collided with regression scripts that
invoke run_e2e.py directly against this same repo as part of their OWN
testing). run_e2e.py moves the staged log into the freshly-rotated
`test/run.log` once it actually rotates. A selection that
never touches that tree (e.g. bare `unit`/`regression`) never rotates and
appends straight to `test/run.<tiers>.log` instead. `--submit`'s job keeps
its SLURM log and packed perf/profile deltas in `scratch/` (never inside the
rotating tree, where a jobid-named file would be destroyed after two more
rotations) and packs them via an EXIT trap so a walltime kill or an
env-check failure still packs whatever ran.
"""
import glob
import io
import json
import os
import re
import shutil
import subprocess
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
# run.py never acquires testsys.runlock itself (PR #56 audit follow-up: an
# earlier version of this file held the test/ lock for its own whole
# invocation, which collided with regression scripts -- e.g.
# test_src_stamp.py -- that invoke testsys/e2e/run_e2e.py directly against
# this SAME real repo_root as part of the 'regression' tier itself. Locking
# is left entirely to run_e2e.py's own Gate 0, taken unconditionally, every
# time it runs, exactly as before this whole feature existed). Nothing here
# needs testsys.runlock at all any more.


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


# The modules every RUNNER below needs from whatever python3 is currently
# active -- checked in TWO places, both citing this same tuple so they
# cannot drift apart: here, in-process, for `--machine <m>` run IN PLACE
# (main()'s machine_name branch); and as a `python3 -c 'import ...'` line in
# the generated --submit job body (build_submit_script), right after
# `source ./install-eqdyna.sh -c <m>` and before the sweep. Owner hit this
# twice on 2026-09-29 running in place with the venv not activated: first
# "No module named netCDF4", then (after xarray was still missing from the
# checked set) "No module named xarray".
REQUIRED_MODULES = ('jax', 'netCDF4', 'xarray', 'mpi4py')


def missing_required_modules():
    """Which of REQUIRED_MODULES this SAME python3 (sys.executable -- every
    RUNNER below invokes it via subprocess) cannot import. Empty means all
    present."""
    missing = []
    for mod in REQUIRED_MODULES:
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    return missing


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
    import_list = ', '.join(REQUIRED_MODULES)
    buf = io.StringIO()
    # The SLURM stdout/stderr capture (`#SBATCH -o`), and everything this
    # job packs alongside it (perf/profile-ledger deltas, the tarball), all
    # live in scratch/ (gitignored) -- NEVER inside test/ (PR #56 audit m2).
    # Each is named by $SLURM_JOB_ID, so it would otherwise accumulate
    # inside the ROTATING test/ tree across many separate job submissions;
    # rule 8 keeps only ONE level of history (test.prev/), so the second
    # rotation after any given job would delete that job's own artifacts
    # for good. Only test/run.log (or test/run.<tiers>.log) -- the SAME
    # filename every run, by design meant to track "this run vs the last
    # one" -- belongs inside the rotating tree, and testsys/run.py's own
    # inner invocation below already handles that by itself.
    machines_mod.write_slurm_header(
        buf, jobname='eqdyna-sweep', nnode=1, ncpu=m['cores_per_node'],
        queue=queue, walltime=wt, account=account, email='',
        output='scratch/eqdyna_sweep_%j.log')
    buf.write('\n')
    buf.write('cd "$SLURM_SUBMIT_DIR"\n')
    buf.write('source ./install-eqdyna.sh -c %s\n' % machine_name)
    buf.write('git rev-parse HEAD; module list 2>&1\n')
    # Registered BEFORE the env-check below (PR #56 audit m1): a walltime
    # kill (SLURM sends SIGTERM before SIGKILL) or the env-check's own
    # `exit 1` must still pack whatever ran, not lose it because the script
    # died before reaching what used to be unconditional packing lines at
    # the very bottom. `pack` is defensive about what actually exists --
    # L0/P0 are unset if the env-check failed before they were computed.
    buf.write('pack() {\n')
    buf.write('  rc=$?\n')
    buf.write('  files="eqdyna_sweep_$SLURM_JOB_ID.log"\n')
    buf.write('  if [ -n "$L0" ] && [ -n "$P0" ]; then\n')
    buf.write('    tail -n +$((L0 + 1)) docs/perf_ledger.jsonl '
              '> "scratch/eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl"\n')
    buf.write('    tail -n +$((P0 + 1)) docs/run_profiles.jsonl '
              '> "scratch/eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl"\n')
    buf.write('    files="$files eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl '
              'eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl"\n')
    buf.write('  fi\n')
    buf.write('  ( cd scratch && tar czf "eqdyna_sweep_$SLURM_JOB_ID.tgz" $files )\n')
    buf.write('  echo "eqdyna-sweep: packed scratch/eqdyna_sweep_'
              '$SLURM_JOB_ID.tgz ($files) -- kept outside test/, which '
              'rotates; this run'"'"'s own console output is in <repo>/test/ '
              '(test/run.log or test/run.<tiers>.log), rotated with the run"\n')
    buf.write('  exit $rc\n')
    buf.write('}\n')
    buf.write('trap pack EXIT\n')
    # Refuse EARLY, before the (potentially hours-long) sweep, if the
    # sourced environment is missing a dependency the sweep needs -- this is
    # exactly what bit the owner on 2026-09-29: a bare `python3` with none of
    # jax/netCDF4/mpi4py (and, found the same day, xarray) importable,
    # discovered only mid-regression as "No module named netCDF4"/"xarray".
    # Explicit check, not `set -e` (rule 2: no silent fallback).
    env_check_msg = (
        "eqdyna-sweep: environment check failed after "
        "source ./install-eqdyna.sh -c %s -- %s not all importable by this "
        "python3. Check that the venv activated (see install-eqdyna.sh's "
        "own warning above if it did not) and EQDYNA_VENV."
        % (machine_name, import_list))
    buf.write("python3 -c 'import %s' || {\n" % import_list)
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
    its exit code.

    EXCLUDES testsys.regression_sweep_exclusions.EXCLUDED_FROM_SWEEP
    (2026-09-30, owner requirement (c)): test_release_complete.py evaluates
    its assertions against a tag that already exists -- meaningful only at
    the moment of auditing a just-cut release, never at an arbitrary later
    commit. Sweeping it in here unconditionally meant it ran on EVERY
    commit and PR (and, via this exact command, inside
    .github/workflows/publish.yml's build-and-push gate too), and once a
    release's post-tag board/Release steps land after its tag (rule 8: a
    pushed tag's tree is immutable), two of its checks fail FOREVER for
    that tag's own CI run -- which then poisoned every commit after it.
    It stays a real, directly-runnable test_*.py file; it is simply not
    part of THIS sweep. Paired with check_pretag_ci.py (rule 15/25) as the
    explicit release-time check instead."""
    from testsys import regression_sweep_exclusions
    print('\n==== testsys: regression ====')
    scripts = sorted(
        p for p in glob.glob(os.path.join(TESTSYS, 'regression', 'test_*.py'))
        if os.path.basename(p) not in regression_sweep_exclusions.EXCLUDED_FROM_SWEEP)
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


# Env vars for the e2e subprocess ONLY (currently just EQDYNA_RUN_LOG_PATH)
# -- set by main() right after _prepare_test_tree() returns, read by _e2e()
# when it builds that ONE subprocess's env. Never written into os.environ:
# every OTHER subprocess in the same run.py invocation (unit's pytest, every
# regression/test_*.py script -- some of which, e.g. test_src_stamp.py,
# themselves launch run_e2e.py "directly" to test ITS OWN standalone
# behaviour) would otherwise inherit it for no reason (PR #56 audit M2).
# An EARLIER version of this dict also carried EQDYNA_TEST_LOCK_HELD, a
# "the caller already holds the lock" claim -- removed entirely (M2, second
# finding, 2026-09-30): run.py never legitimately holds this lock itself
# (see _prepare_test_tree's docstring), so the flag could only ever be
# leaked or hand-set, and even a "verified" claim was unsafe (see
# run_e2e.acquire_test_lock's docstring). run_e2e.py now always takes the
# lock itself, unconditionally.
_E2E_EXTRA_ENV = {}


def _e2e(*args):
    env = dict(os.environ)
    env.update(_E2E_EXTRA_ENV)
    return subprocess.call([sys.executable, os.path.join(TESTSYS, 'e2e', 'run_e2e.py')]
                           + list(args), cwd=REPO_ROOT, env=env)


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


def _release_sweep_carry_forward():
    """(skip, message) -- item 3 (owner, 2026-09-30, "fewer full sweeps and
    fewer releases"): does the most recently committed
    docs/evidence/sweep-*/summary.json already cover this tree's
    release-physics content (testsys.content_key.compute_release_physics,
    over testsys.change_class.is_release_physics_path's closed, owner-named
    set)? If so, skip=True -- no release-physics path has changed since
    that sweep, and it can be carried forward instead of re-run.

    Fails toward running the full sweep (skip=False) on ANY doubt: no
    evidence found, an unreadable summary.json, evidence with no
    release_physics_key recorded (written before this field existed), or a
    git/HEAD resolution error -- rule 2, same posture as needs_sweep.py and
    check_pretag_ci.py's own error handling.

    Advisory only for THIS function's caller (run_release, deciding whether
    to spend the minutes-to-hours re-running run_e2e.py): the actual GATE
    enforced at tag time is testsys/regression/check_pretag_ci.py's
    evaluate_sweep_candidate, which independently re-derives the same
    answer against whatever sha is about to be tagged. A bug here can only
    cost an unnecessary sweep, never let a stale one through -- that
    guarantee lives entirely in check_pretag_ci.py, unchanged by this
    function's verdict.

    Imports testsys.content_key lazily, on purpose (test_src_stamp.py check
    B runs a sandboxed `testsys/run.py` with only run.py itself and no
    sibling testsys/*.py present, to prove the stale-binary refusal fires
    before anything heavier loads; a module-level import here would break
    that gate with an ImportError instead of the REFUSED it expects, for a
    tier -- `release` -- that sandboxed run never even selects)."""
    from testsys import content_key
    found = sorted(glob.glob(os.path.join(REPO_ROOT, 'docs', 'evidence',
                                          'sweep-*', 'summary.json')))
    if not found:
        return False, ('no docs/evidence/sweep-*/summary.json found -- '
                       'running the full sweep')
    head_r = subprocess.run(['git', '-C', REPO_ROOT, 'rev-parse', 'HEAD'],
                            capture_output=True, text=True)
    head_sha = head_r.stdout.strip()
    if head_r.returncode != 0 or len(head_sha) != 40:
        return False, 'could not resolve HEAD -- running the full sweep'
    try:
        head_key = content_key.compute_release_physics(REPO_ROOT, head_sha)
    except Exception as exc:  # pragma: no cover - defensive, git itself failing
        return False, ('could not compute this tree\'s release_physics_key '
                       '(%s) -- running the full sweep' % exc)
    reasons = []
    for p in found:
        try:
            with open(p) as fh:
                data = json.load(fh)
        except (OSError, ValueError) as exc:
            reasons.append('%s: unreadable (%s)' % (p, exc))
            continue
        swept_key = data.get('release_physics_key')
        swept_sha = data.get('sha')
        if swept_key is None:
            reasons.append('%s: no release_physics_key recorded (older evidence)' % p)
            continue
        if swept_key == head_key:
            return True, ('%s (swept sha %s) already covers this tree\'s '
                         'release-physics content (release_physics_key %s) '
                         '-- carrying that evidence forward, not re-running '
                         'the sweep' % (p, swept_sha, head_key))
        reasons.append('%s: release_physics_key %s != this tree\'s %s'
                       % (p, swept_key, head_key))
    return False, ('no committed sweep evidence has a matching '
                   'release_physics_key -- running the full sweep:\n  %s'
                   % '\n  '.join(reasons))


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

    Item 3 (owner, 2026-09-30): before running anything, checks
    `_release_sweep_carry_forward()` -- if the most recent committed
    release evidence already covers this tree's release-physics content, it
    SKIPS the sweep and reports the carried-forward evidence instead of
    invoking run_e2e.py again. Any doubt runs the full sweep (rule 2); see
    that function's own docstring for why this is safe even if it is wrong.

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
    skip, msg = _release_sweep_carry_forward()
    print('testsys: release: %s' % msg)
    if skip:
        return 0
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

# Tiers whose RUNNER ends up invoking testsys/e2e/run_e2e.py (via _e2e
# above), i.e. touches the shared, rule-8-rotated `test/` run tree.
# `e2e-full` is excluded: it rotates `test.full/`, a SEPARATE tree, and is
# opt-in/hours-long, never part of an everyday selection.
TOUCHES_TEST_TREE = frozenset({'e2e', 'e2e-ci', 'release', 'gpu'})


class _Tee:
    """Duplicate this process's stdout/stderr -- including every
    subprocess.call'd child's inherited fd, which a Python-level
    sys.stdout wrapper cannot reach -- to `path` as well as the terminal,
    via a real `tee` child. Same place on every machine (owner, 2026-09-29):
    ubuntu, Docker, and ls6 all have coreutils' `tee`."""

    def __init__(self, path, append=False):
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        self._saved_out = os.dup(1)
        self._saved_err = os.dup(2)
        self._proc = subprocess.Popen(['tee', '-a', path] if append
                                      else ['tee', path], stdin=subprocess.PIPE)
        os.dup2(self._proc.stdin.fileno(), 1)
        os.dup2(self._proc.stdin.fileno(), 2)

    def close(self):
        # flush()/write() below can raise BrokenPipeError if the tee child
        # already died (m3, PR #56 audit); restore the real fds FIRST, in a
        # nested try, so a dead tee never leaves this process's own
        # stdout/stderr pointed at a closed pipe.
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except (BrokenPipeError, OSError):
            pass
        os.dup2(self._saved_out, 1)
        os.dup2(self._saved_err, 2)
        os.close(self._saved_out)
        os.close(self._saved_err)
        try:
            self._proc.stdin.close()
        except (BrokenPipeError, OSError):
            pass
        rc = self._proc.wait()
        if rc != 0:
            print('testsys: WARNING - the tee child that mirrors console '
                  'output to the log exited %d -- the log may be truncated '
                  'or missing lines' % rc, file=sys.stderr)


def _tiers_slug(selected):
    return '-'.join(selected)


def _prepare_test_tree(selected):
    """Where this invocation's own console log goes, and the env vars (if
    any) to pass to the e2e subprocess -- rule 8 (preserve, never delete)
    applied to `test/`.

    NEITHER ROTATION NOR THE LOCK IS TAKEN HERE (PR #56 audit, blocker B1
    and a same-day follow-up finding). The FIRST attempt at this fix hoisted
    BOTH to run.py's own start for an e2e-family selection -- rotation was
    always wrong that early (a refused or Ctrl-C'd invocation would still
    rotate the PREVIOUS real sweep's evidence into test.prev/, rule 8's one
    level of history, and the NEXT run's rotation would then delete it for
    good), and holding the lock that early turned out to be wrong TOO: it
    collided with regression scripts that invoke testsys/e2e/run_e2e.py
    directly, against this SAME real repo_root, as part of the 'regression'
    tier ITSELF (test_src_stamp.py's stale-EQDYNA_E2E_BIN check, among
    others) -- their child run_e2e.py inherited no claim, tried a REAL
    acquire, and correctly refused, but for the WRONG reason (the ancestor
    run.py's own held lock, not the stale binary the test was checking for).
    Before PR #56 this collision could not happen: nothing held the lock
    until run_e2e.py's OWN invocation, which runs LAST (the `e2e` tier,
    after `unit`/`regression`). Locking, like rotation, is now left exactly
    where it always was: inside run_e2e.py's own Gate 0, taken
    unconditionally, every time it runs -- whether invoked as run.py's `e2e`
    subprocess or by a regression script testing it directly.

    What IS done here, for a selection touching TOUCHES_TEST_TREE: the
    console log is staged OUTSIDE test/, at scratch/run.<pid>.log, so this
    process never touches test/ or test.prev/ at all, and the staged path
    is returned in e2e_env (EQDYNA_RUN_LOG_PATH) for the CALLER to pass
    ONLY to the e2e subprocess's own environment (see _E2E_EXTRA_ENV /
    _e2e), never into os.environ, which every OTHER subprocess in this same
    invocation would otherwise inherit for no reason. run_e2e.py moves that
    staged log into test/run.log once its OWN gates pass and it rotates.

    A selection that never touches TOUCHES_TEST_TREE never rotates, and
    appends to test/run.<tiers>.log directly (no staging needed: there is
    no rotation to protect against; PR #56 audit m6 notes this path is
    unlocked -- two concurrent bare `unit`-only invocations can interleave
    lines in that shared append-mode file -- accepted as a follow-up, not
    fixed here).

    Returns (log_path, append, e2e_env).
    """
    test_dir = os.path.join(REPO_ROOT, 'test')
    if not any(t in TOUCHES_TEST_TREE for t in selected):
        os.makedirs(test_dir, exist_ok=True)
        return os.path.join(test_dir, 'run.%s.log' % _tiers_slug(selected)), True, {}

    scratch = os.path.join(REPO_ROOT, 'scratch')
    os.makedirs(scratch, exist_ok=True)
    staged_log = os.path.join(scratch, 'run.%d.log' % os.getpid())
    return staged_log, False, {'EQDYNA_RUN_LOG_PATH': staged_log}


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
        print('testsys/run.py: submitted %s (job %s); once the job '
              'completes, this run\'s own console output lands in '
              '<repo>/test/ (test/run.log or test/run.<tiers>.log, rotated '
              'with the run), and scratch/eqdyna_sweep_%s.tgz holds the '
              'SLURM log plus perf/profile deltas (kept outside test/ so '
              'rotation never deletes it -- packed via an EXIT trap even on '
              'a walltime kill or an env-check failure)'
              % (script_path, jobid, jobid))
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
        # Same early refusal as the --submit job body, for the run-IN-PLACE
        # path (e.g. inside an interactive allocation): the owner hit this
        # today running --machine ls6 without --submit and without having
        # activated the venv -- "No module named xarray" mid-regression.
        missing = missing_required_modules()
        if missing:
            print('testsys/run.py: REFUSED - %s not importable by %s -- the '
                  'venv for %s is not active (see install-eqdyna.sh -c %s\'s '
                  'own warning if it printed one) or EQDYNA_VENV points '
                  'elsewhere. Activate it (`source install-eqdyna.sh -c %s`), '
                  'or use --submit so the sbatch job does it for you.'
                  % (', '.join(missing), sys.executable, machine_name,
                     machine_name, machine_name))
            return 2

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

    # PR #56 audit m4: two e2e-family tiers in ONE invocation would each
    # rotate test/ -> test.prev/ in turn, the second one destroying the
    # first's fresh results (rule 8 keeps only one level of history).
    touching = [t for t in selected if t in TOUCHES_TEST_TREE]
    if len(touching) > 1:
        print('testsys/run.py: selection runs more than one e2e-family tier '
              'in a single invocation (%s) -- each rotates test/ again, '
              'destroying the previous one\'s results within this SAME run. '
              'Run them separately (two invocations, or two worktrees).'
              % ', '.join(touching))
        return 2

    log_path, append, e2e_env = _prepare_test_tree(selected)

    global _E2E_EXTRA_ENV
    _E2E_EXTRA_ENV = e2e_env

    tee = _Tee(log_path, append=append)
    try:
        print('testsys: console output also written to %s' % log_path)

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
    finally:
        tee.close()


if __name__ == '__main__':
    sys.exit(main(sys.argv))
