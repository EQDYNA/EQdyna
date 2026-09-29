"""Unit tests for testsys/run.py's `--machine`/`--submit` support (owner
request 2026-09-29): the generated sbatch job script, and the CLI's refusal
to guess a scheduler, an account, or a core count. `sbatch` itself is
stubbed -- these tests never touch a real scheduler.

Both ways (rule 14a): a slurm machine with an account builds a script and a
non-slurm machine (or a slurm one with no account) refuses.
"""
import importlib.util
import io
import os

from conftest import REPO_ROOT

_spec = importlib.util.spec_from_file_location(
    'run_under_test', os.path.join(REPO_ROOT, 'testsys', 'run.py'))
run = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run)

machines = run._load_machines()


def test_parse_argv_extracts_flags_and_leaves_tiers():
    tiers, machine, submit_flag, account, partition, walltime = run.parse_argv(
        ['run.py', 'unit', 'regression', '--machine', 'ls6', '--submit',
         '--account', 'ACCT-9', '--partition', 'dev', '--time', '00:30:00'])
    assert tiers == ['unit', 'regression']
    assert machine == 'ls6'
    assert submit_flag is True
    assert account == 'ACCT-9'
    assert partition == 'dev'
    assert walltime == '00:30:00'


def test_parse_argv_no_machine_leaves_everything_none():
    tiers, machine, submit_flag, account, partition, walltime = run.parse_argv(
        ['run.py', 'all'])
    assert tiers == ['all']
    assert machine is None and submit_flag is False
    assert account is None and partition is None and walltime is None


def test_build_submit_script_header_and_body_for_ls6():
    script = run.build_submit_script(machines, 'ls6', ['unit', 'regression'],
                                     account='ACCT-1')
    lines = script.splitlines()
    assert lines[0] == '#! /bin/bash'
    assert '#SBATCH -J eqdyna-sweep' in lines
    assert '#SBATCH -o eqdyna_sweep_%j.log' in lines
    assert '#SBATCH -N 1' in lines
    assert '#SBATCH -n 128' in lines          # ls6 cores_per_node
    assert '#SBATCH -p normal' in lines       # ls6 default partition
    assert '#SBATCH -A ACCT-1' in lines
    assert 'source ./install-eqdyna.sh -c ls6' in script
    assert 'python3 testsys/run.py unit regression --machine ls6' in script
    assert 'exit $rc' in script
    # Never duplicates a module load / venv step -- install-eqdyna.sh -c ls6
    # is the one place that happens.
    assert 'module load' not in script
    assert 'pip install' not in script


def test_build_submit_script_partition_and_time_overrides():
    script = run.build_submit_script(machines, 'ls6', ['e2e-ci'],
                                     account='A', partition='dev',
                                     walltime='00:15:00')
    lines = script.splitlines()
    assert '#SBATCH -p dev' in lines
    assert '#SBATCH -t 00:15:00' in lines


def test_build_submit_script_refuses_on_non_slurm_machine():
    for name in ('ubuntu', 'macos'):
        try:
            run.build_submit_script(machines, name, ['unit'], account='A')
            assert False, 'expected ValueError for %s' % name
        except ValueError as exc:
            assert 'no scheduler' in str(exc)


def test_build_submit_script_refuses_without_account():
    try:
        run.build_submit_script(machines, 'ls6', ['unit'], account=None)
        assert False, 'expected ValueError'
    except ValueError as exc:
        assert '--account' in str(exc)


def test_build_submit_script_refuses_when_cores_per_node_unknown():
    """grace's registry entry deliberately carries cores_per_node=None
    (not yet verified) -- submit must refuse rather than guess (rule 2)."""
    try:
        run.build_submit_script(machines, 'grace', ['unit'], account='A',
                                partition='sn4')
        assert False, 'expected ValueError'
    except ValueError as exc:
        assert 'cores_per_node' in str(exc)


def test_main_submit_writes_script_and_calls_sbatch(monkeypatch):
    calls = []

    class FakeCompleted:
        returncode = 0
        stdout = 'Submitted batch job 424242\n'
        stderr = ''

    def fake_run(cmd, cwd, capture_output, text):
        calls.append(cmd)
        return FakeCompleted()

    monkeypatch.setattr(run.subprocess, 'run', fake_run)
    # RUNNERS must never be invoked on the --submit path: this dispatches a
    # job, it does not run tiers in this process.
    monkeypatch.setitem(run.RUNNERS, 'unit',
                        lambda: (_ for _ in ()).throw(AssertionError(
                            'RUNNERS must not run locally under --submit')))

    script_path = os.path.join(REPO_ROOT, 'scratch',
                               'submit_eqdyna_sweep_ls6.sh')
    try:
        rc = run.main(['run.py', 'unit', '--machine', 'ls6', '--submit',
                      '--account', 'ACCT-7'])
        assert rc == 0
        assert len(calls) == 1
        assert calls[0][0] == 'sbatch'
        assert calls[0][1] == script_path
        assert os.path.exists(script_path)
        with open(script_path) as f:
            content = f.read()
        assert 'python3 testsys/run.py unit --machine ls6' in content
        assert '#SBATCH -A ACCT-7' in content
    finally:
        if os.path.exists(script_path):
            os.remove(script_path)


def test_main_submit_on_non_slurm_machine_refuses_without_calling_sbatch(monkeypatch):
    called = []
    monkeypatch.setattr(run.subprocess, 'run',
                        lambda *a, **k: called.append(1))
    rc = run.main(['run.py', 'unit', '--machine', 'ubuntu', '--submit',
                  '--account', 'ACCT-7'])
    assert rc == 2
    assert not called


def test_main_machine_flag_sets_env_without_submit(monkeypatch):
    monkeypatch.delenv('EQDYNA_TEST_MACHINE', raising=False)
    monkeypatch.delenv('EQDYNA_MPIRUN', raising=False)
    seen_env = {}

    def fake_unit():
        seen_env['EQDYNA_TEST_MACHINE'] = os.environ.get('EQDYNA_TEST_MACHINE')
        seen_env['EQDYNA_MPIRUN'] = os.environ.get('EQDYNA_MPIRUN')
        return 0

    monkeypatch.setitem(run.RUNNERS, 'unit', fake_unit)
    rc = run.main(['run.py', 'unit', '--machine', 'ls6'])
    assert rc == 0
    assert seen_env == {'EQDYNA_TEST_MACHINE': 'ls6', 'EQDYNA_MPIRUN': 'ibrun'}


def test_main_without_machine_leaves_env_untouched(monkeypatch):
    monkeypatch.delenv('EQDYNA_TEST_MACHINE', raising=False)
    monkeypatch.delenv('EQDYNA_MPIRUN', raising=False)
    monkeypatch.setitem(run.RUNNERS, 'unit', lambda: 0)
    rc = run.main(['run.py', 'unit'])
    assert rc == 0
    assert 'EQDYNA_TEST_MACHINE' not in os.environ
    assert 'EQDYNA_MPIRUN' not in os.environ
