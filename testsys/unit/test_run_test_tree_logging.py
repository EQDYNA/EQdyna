"""Unit tests for testsys/run.py's test/ logging and rotation (owner,
2026-09-29: "logs must be saved with the runs, the SAME place on every
machine"; PR #56 audit fixing the first pass at this).

Both ways (rule 14a): an e2e-family selection stages its log and takes the
lock without rotating; a bare unit/regression one appends without rotating;
a locked tree refuses; a REFUSED run (blocker B1) leaves test/ and
test.prev/ completely untouched; two e2e-family tiers in one invocation
(m4) refuse rather than double-rotating.

M3 isolation: every test here either patches run.REPO_ROOT to a tmp_path
before anything can touch a real path, or stubs out the one function that
would (_prepare_test_tree), AND the module-level `_clean_environ` fixture
below restores os.environ around every test -- main()'s --machine handling
mutates it directly (os.environ['EQDYNA_TEST_MACHINE'] = ...), which
monkeypatch cannot undo on its own since it never made that change itself.
"""
import importlib.util
import os

import pytest

from conftest import REPO_ROOT

_spec = importlib.util.spec_from_file_location(
    'run_under_test', os.path.join(REPO_ROOT, 'testsys', 'run.py'))
run = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run)

from testsys import runlock  # noqa: E402


@pytest.fixture(autouse=True)
def _clean_environ():
    """Restores os.environ after EVERY test in this file, regardless of
    what run.main()/_prepare_test_tree mutated directly (PR #56 audit M3:
    EQDYNA_TEST_MACHINE=ls6, EQDYNA_MPIRUN etc. were leaking into later
    tests in the same pytest process)."""
    saved = dict(os.environ)
    yield
    os.environ.clear()
    os.environ.update(saved)


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    """A throwaway repo_root: run._prepare_test_tree/_Tee/runlock all take a
    repo_root argument or read run.REPO_ROOT, so pointing that at tmp_path
    isolates every test here from this checkout's REAL test/ tree."""
    monkeypatch.setattr(run, 'REPO_ROOT', str(tmp_path))
    return tmp_path


def test_touches_test_tree_names_every_e2e_family_tier():
    """Declarative check of the set _prepare_test_tree gates rotation on."""
    assert run.TOUCHES_TEST_TREE == frozenset({'e2e', 'e2e-ci', 'release', 'gpu'})
    assert 'e2e-full' not in run.TOUCHES_TEST_TREE  # separate tree, opt-in


def test_e2e_family_selection_stages_log_without_rotating_or_locking(sandbox):
    """PR #56 audit blocker B1 (plus a same-day follow-up finding): run.py
    itself must NEVER rotate test/ -> test.prev/, and -- the follow-up --
    must never take the test/ lock either, only stage a log outside it.
    Rotation AND locking happen lazily, inside run_e2e.py, once ITS OWN
    gates pass (early lock-holding here collided with regression scripts,
    e.g. test_src_stamp.py, that invoke run_e2e.py directly against this
    same repo_root as part of the 'regression' tier itself)."""
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'previous_evidence.txt').write_text('old run')

    log_path, append, e2e_env = run._prepare_test_tree(['unit', 'regression', 'e2e'])

    assert append is False
    # Staged OUTSIDE test/, under scratch/ -- never inside the tree that
    # might still need rotating.
    assert log_path == str(sandbox / 'scratch' / ('run.%d.log' % os.getpid()))
    assert e2e_env == {'EQDYNA_RUN_LOG_PATH': log_path}
    # Nothing in os.environ -- this goes to the e2e subprocess's env only
    # (PR #56 audit M2).
    assert 'EQDYNA_RUN_LOG_PATH' not in os.environ
    # test/ is UNTOUCHED: no rotation happened here.
    assert (test_dir / 'previous_evidence.txt').read_text() == 'old run'
    assert not (sandbox / 'test.prev').exists()
    # No lock taken either -- a real acquire on the same tree still succeeds.
    from testsys import runlock
    lock = runlock.acquire(str(sandbox), 'test', announce=False)
    lock.release()


def test_non_e2e_selection_never_rotates_and_appends(sandbox):
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'evidence.txt').write_text('kept')

    log_path, append, e2e_env = run._prepare_test_tree(['unit', 'regression'])

    assert append is True
    assert e2e_env == {}
    assert log_path == str(test_dir / 'run.unit-regression.log')
    assert not (sandbox / 'test.prev').exists()
    assert (test_dir / 'evidence.txt').read_text() == 'kept'


def test_e2e_family_selection_does_not_refuse_when_tree_already_locked(sandbox):
    """_prepare_test_tree itself never touches the lock at all any more, so
    it must succeed (staging a log) even while something else genuinely
    holds test/'s lock -- the REFUSAL now belongs entirely to run_e2e.py's
    own Gate 0 (testsys/unit/test_run_e2e_test_tree_flags.py's
    acquire_test_lock tests, and testsys/regression/test_e2e_run_tree_lock.py
    for the real end-to-end refusal)."""
    held = runlock.acquire(str(sandbox), 'test')
    try:
        log_path, append, e2e_env = run._prepare_test_tree(['e2e'])
        assert e2e_env == {'EQDYNA_RUN_LOG_PATH': log_path}
    finally:
        held.release()


def test_multiple_e2e_family_tiers_refuse_without_touching_test_tree(sandbox, monkeypatch):
    """PR #56 audit m4: `e2e` and `gpu` in ONE invocation would each rotate
    test/ in turn, the second one destroying the first's fresh results."""
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'evidence.txt').write_text('kept')
    monkeypatch.setitem(run.RUNNERS, 'e2e',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'RUNNERS must not run when the multi-e2e-tier check refuses')))
    monkeypatch.setitem(run.RUNNERS, 'gpu',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'RUNNERS must not run when the multi-e2e-tier check refuses')))

    rc = run.main(['run.py', 'e2e', 'gpu'])

    assert rc == 2
    assert (test_dir / 'evidence.txt').read_text() == 'kept'
    assert not (sandbox / 'test.prev').exists()
    assert not (sandbox / '.test.lock').exists()  # never even reached the lock


def test_b1_refused_run_leaves_test_and_test_prev_untouched(sandbox, monkeypatch):
    """PR #56 audit blocker B1, the exact regression the audit asked for:
    a run refused by require_fresh_fortran_binary (standing in for a stale
    binary) must leave test/ and test.prev/ completely untouched -- this
    goes RED on cecc09e, where _prepare_test_tree rotated test/ BEFORE this
    check ever ran, for any selection ('all') that included e2e."""
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'previous_evidence.txt').write_text('do not touch')

    monkeypatch.setattr(run, 'require_fresh_fortran_binary',
                        lambda *a, **k: 'stale binary (forced for this test)')
    monkeypatch.setitem(run.RUNNERS, 'unit',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'RUNNERS must not run at all once regression refuses first')))
    monkeypatch.setitem(run.RUNNERS, 'e2e',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'e2e must never run when the fresh-binary check already refused')))

    rc = run.main(['run.py', 'all'])

    assert rc == 1
    assert (test_dir / 'previous_evidence.txt').read_text() == 'do not touch'
    assert not (sandbox / 'test.prev').exists()


def test_tee_writes_own_and_subprocess_output_to_file(sandbox):
    """Both writers land in the file: a raw write to fd 1 (what this
    process's own C-level/os.write output looks like, and exactly what a
    subprocess's inherited fd is) and a real subprocess's stdout -- checked
    via os.write rather than print()/sys.stdout, which pytest's own default
    capture reassigns to a SEPARATE Python-level object decoupled from the
    numeric fd this test manipulates (verified: print() output vanished from
    the tee'd file under pytest's default capture even though _Tee dup2's fd
    1 correctly -- an artifact of the test harness, not of _Tee)."""
    log_path = str(sandbox / 'run.log')
    tee = run._Tee(log_path)
    try:
        os.write(1, b'own fd-level line\n')
        run.subprocess.call(['echo', 'child process line'])
    finally:
        tee.close()
    with open(log_path) as f:
        content = f.read()
    assert 'own fd-level line' in content
    assert 'child process line' in content


def test_tee_close_warns_but_does_not_raise_when_tee_child_already_dead(sandbox, capsys):
    """PR #56 audit m3: a dead tee child must not crash close() -- BrokenPipe
    on the flush/write is swallowed, and a non-zero tee exit is reported."""
    log_path = str(sandbox / 'run.log')
    tee = run._Tee(log_path)
    tee._proc.kill()
    tee._proc.wait()
    tee.close()  # must not raise
    err = capsys.readouterr().err
    assert 'WARNING' in err and 'tee' in err


def test_missing_required_modules_detects_unimportable_name(monkeypatch):
    monkeypatch.setattr(run, 'REQUIRED_MODULES', ('os', 'sys'))
    assert run.missing_required_modules() == []

    monkeypatch.setattr(run, 'REQUIRED_MODULES',
                        ('os', 'this_module_does_not_exist_xyz'))
    assert run.missing_required_modules() == ['this_module_does_not_exist_xyz']


def test_main_machine_flag_refuses_when_modules_missing(monkeypatch):
    monkeypatch.setattr(run, 'REQUIRED_MODULES',
                        ('this_module_does_not_exist_xyz',))
    monkeypatch.setitem(run.RUNNERS, 'unit',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'RUNNERS must not run when the module check refuses')))
    rc = run.main(['run.py', 'unit', '--machine', 'ls6'])
    assert rc == 2


def test_main_machine_flag_proceeds_when_modules_present(monkeypatch, tmp_path):
    """_load_machines() needs the REAL REPO_ROOT (it loads scripts/
    machines.py by path), so REPO_ROOT itself is not sandboxed here --
    _prepare_test_tree is stubbed instead (PR #56 audit M3), which is what
    would otherwise touch this checkout's real test/."""
    log_path = str(tmp_path / 'test' / 'run.unit.log')
    monkeypatch.setattr(run, 'REQUIRED_MODULES', ('os', 'sys'))
    monkeypatch.setattr(run, '_prepare_test_tree',
                        lambda selected: (log_path, True, {}))
    monkeypatch.setitem(run.RUNNERS, 'unit', lambda: 0)
    rc = run.main(['run.py', 'unit', '--machine', 'ls6'])
    assert rc == 0
    assert os.path.isfile(log_path)


def test_build_submit_script_env_check_includes_xarray_before_sweep():
    machines = run._load_machines()
    script = run.build_submit_script(machines, 'ls6', ['unit'], account='A')
    assert "import jax, netCDF4, xarray, mpi4py" in script
    check_pos = script.find("import jax, netCDF4, xarray, mpi4py")
    sweep_pos = script.find('python3 testsys/run.py unit --machine ls6')
    assert check_pos != -1 and sweep_pos != -1
    assert check_pos < sweep_pos


def test_build_submit_script_packing_stays_in_scratch_and_traps_exit():
    """PR #56 audit m1 (EXIT trap so a walltime kill or env-check failure
    still packs) and m2 (nothing job-scoped lives inside the rotating
    test/ tree, where two more rotations would delete it)."""
    machines = run._load_machines()
    script = run.build_submit_script(machines, 'ls6', ['unit'], account='A')
    lines = script.splitlines()
    assert '#SBATCH -o scratch/eqdyna_sweep_%j.log' in lines
    assert 'trap pack EXIT' in lines
    # Nothing job-scoped (log, ledger delta, profile delta, tarball) ever
    # targets test/ -- it all stays in scratch/, which is never rotated.
    assert 'test/eqdyna_sweep' not in script
    assert 'scratch/eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl' in script
    assert 'scratch/eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl' in script
    assert 'tar czf "eqdyna_sweep_$SLURM_JOB_ID.tgz"' in script
    # The trap is registered before the env-check's own `exit 1`, so a
    # failure there still triggers packing.
    trap_pos = script.find('trap pack EXIT')
    env_check_exit_pos = script.find("import jax, netCDF4, xarray, mpi4py' || {")
    assert 0 <= trap_pos < env_check_exit_pos
