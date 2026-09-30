"""Unit tests for testsys/run.py's test/ logging and rotation-hoisting
(owner, 2026-09-29): "logs must be saved with the runs, the SAME place on
every machine" -- every invocation's console output lands under
`<repo>/test/`, and a selection touching the e2e run tree gets its rotation
hoisted to run.py's own start, ahead of unit/regression.

Both ways (rule 14a): an e2e-family selection rotates and a bare
unit/regression one appends without rotating; a locked tree refuses.
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


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    """A throwaway repo_root: run._prepare_test_tree/_Tee/runlock all take a
    repo_root argument or read run.REPO_ROOT, so pointing that at tmp_path
    isolates every test here from this checkout's REAL test/ tree."""
    monkeypatch.setattr(run, 'REPO_ROOT', str(tmp_path))
    monkeypatch.delenv('EQDYNA_TEST_ALREADY_ROTATED', raising=False)
    monkeypatch.delenv('EQDYNA_TEST_LOCK_HELD', raising=False)
    return tmp_path


def test_touches_test_tree_names_every_e2e_family_tier():
    """Declarative check of the set _prepare_test_tree gates rotation on --
    the behavioural exercise below covers one representative member ('e2e');
    re-running it per member would re-acquire the SAME lockfile in one
    process, which flock refuses by design (see the refusal test)."""
    assert run.TOUCHES_TEST_TREE == frozenset({'e2e', 'e2e-ci', 'release', 'gpu'})
    assert 'e2e-full' not in run.TOUCHES_TEST_TREE  # separate tree, opt-in


def test_e2e_family_selection_rotates_and_marks_env(sandbox):
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'previous_evidence.txt').write_text('old run')

    log_path, append = run._prepare_test_tree(['unit', 'regression', 'e2e'])

    assert append is False
    assert log_path == str(test_dir / 'run.log')
    assert os.environ.get('EQDYNA_TEST_ALREADY_ROTATED') == '1'
    assert os.environ.get('EQDYNA_TEST_LOCK_HELD') == '1'
    prev = sandbox / 'test.prev'
    assert (prev / 'previous_evidence.txt').read_text() == 'old run'
    assert test_dir.is_dir() and not (test_dir / 'previous_evidence.txt').exists()


def test_non_e2e_selection_never_rotates_and_appends(sandbox):
    test_dir = sandbox / 'test'
    test_dir.mkdir()
    (test_dir / 'evidence.txt').write_text('kept')

    log_path, append = run._prepare_test_tree(['unit', 'regression'])

    assert append is True
    assert log_path == str(test_dir / 'run.unit-regression.log')
    assert not (sandbox / 'test.prev').exists()
    assert (test_dir / 'evidence.txt').read_text() == 'kept'
    assert 'EQDYNA_TEST_ALREADY_ROTATED' not in os.environ
    assert 'EQDYNA_TEST_LOCK_HELD' not in os.environ


def test_e2e_family_selection_refuses_when_tree_already_locked(sandbox):
    held = runlock.acquire(str(sandbox), 'test')
    try:
        with pytest.raises(run.TestTreeLocked):
            run._prepare_test_tree(['e2e'])
    finally:
        held.release()


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
    """Uses the REAL REPO_ROOT (unlike the sandbox-based tests above) so
    _load_machines() finds the real scripts/machines.py; _prepare_test_tree
    is stubbed instead, since this test is about the module-check gate, not
    test/ rotation -- covered separately above."""
    monkeypatch.setattr(run, 'REQUIRED_MODULES', ('os', 'sys'))
    monkeypatch.setattr(run, '_prepare_test_tree',
                        lambda selected: (str(tmp_path / 'run.log'), True))
    monkeypatch.setitem(run.RUNNERS, 'unit', lambda: 0)
    rc = run.main(['run.py', 'unit', '--machine', 'ls6'])
    assert rc == 0


def test_build_submit_script_env_check_includes_xarray_before_sweep():
    machines = run._load_machines()
    script = run.build_submit_script(machines, 'ls6', ['unit'], account='A')
    assert "import jax, netCDF4, xarray, mpi4py" in script
    check_pos = script.find("import jax, netCDF4, xarray, mpi4py")
    sweep_pos = script.find('python3 testsys/run.py unit --machine ls6')
    assert check_pos != -1 and sweep_pos != -1
    assert check_pos < sweep_pos


def test_build_submit_script_output_and_packing_target_test_dir():
    machines = run._load_machines()
    script = run.build_submit_script(machines, 'ls6', ['unit'], account='A')
    lines = script.splitlines()
    assert '#SBATCH -o scratch/eqdyna_sweep_%j.log' in lines
    # nothing writes eqdyna_sweep_* directly at repo root any more
    assert '> eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl' not in script
    assert '> "test/eqdyna_sweep_$SLURM_JOB_ID.ledger.jsonl"' in script
    assert '> "test/eqdyna_sweep_$SLURM_JOB_ID.profiles.jsonl"' in script
    assert 'mv "scratch/eqdyna_sweep_$SLURM_JOB_ID.log" "test/eqdyna_sweep_$SLURM_JOB_ID.log"' in script
    assert 'tar czf "test/eqdyna_sweep_$SLURM_JOB_ID.tgz" -C test' in script
