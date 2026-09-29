"""Unit tests for scripts/machines.py -- the one HPC-machine registry and the
SBATCH-header writer it shares with scripts/case.setup (owner request
2026-09-29: "a refactored way to add [HPC] support", "you shouldn't reinvent
everything").

Both ways (rule 14a): a known machine resolves and an unknown one raises; a
slurm machine reads as slurm and a local one reads as no-scheduler.
"""
import io

import machines


def test_known_machines_present_with_required_shape():
    for name in ('ls6', 'grace', 'ubuntu', 'macos'):
        m = machines.machine(name)
        for key in ('scheduler', 'partition', 'cores_per_node', 'walltime',
                    'mpirun', 'notes'):
            assert key in m, (name, key)
        assert m['scheduler'] in (None, 'slurm'), (name, m['scheduler'])
        assert m['mpirun'], name


def test_unknown_machine_raises_naming_known_machines():
    try:
        machines.machine('nonexistent-hpc')
        assert False, 'expected ValueError'
    except ValueError as exc:
        msg = str(exc)
        assert 'nonexistent-hpc' in msg
        for name in machines.MACHINES:
            assert name in msg


def test_ls6_and_grace_are_slurm_with_ibrun_launcher():
    for name in ('ls6', 'grace'):
        m = machines.machine(name)
        assert m['scheduler'] == 'slurm', name
        # docs/user/troubleshooting.md: "on ls6 and grace the launcher is
        # ibrun -> srun".
        assert m['mpirun'] == 'ibrun', name


def test_ubuntu_and_macos_have_no_scheduler_and_use_mpirun():
    for name in ('ubuntu', 'macos'):
        m = machines.machine(name)
        assert m['scheduler'] is None, name
        assert m['mpirun'] == 'mpirun', name


def test_write_slurm_header_emits_every_required_directive():
    buf = io.StringIO()
    machines.write_slurm_header(
        buf, jobname='j', nnode=2, ncpu=8, queue='normal', walltime='01:00:00',
        account='ACCT-1', email='a@b.c')
    out = buf.getvalue()
    lines = out.splitlines()
    assert lines[0] == '#! /bin/bash'
    assert '#SBATCH -J j' in lines
    assert '#SBATCH -N 2' in lines
    assert '#SBATCH -n 8' in lines
    assert '#SBATCH -p normal' in lines
    assert '#SBATCH -t 01:00:00' in lines
    assert '#SBATCH -A ACCT-1' in lines
    assert '#SBATCH --mail-user=a@b.c' in lines


def test_write_slurm_header_output_defaults_and_overrides():
    buf = io.StringIO()
    machines.write_slurm_header(buf, jobname='j', nnode=1, ncpu=1, queue='q',
                                walltime='00:01:00', account='A', email='')
    assert '#SBATCH -o a.eqdyna.log%j' in buf.getvalue().splitlines()

    buf2 = io.StringIO()
    machines.write_slurm_header(buf2, jobname='j', nnode=1, ncpu=1, queue='q',
                                walltime='00:01:00', account='A', email='',
                                output='sweep_%j.log')
    assert '#SBATCH -o sweep_%j.log' in buf2.getvalue().splitlines()


def test_case_setup_header_matches_shared_writer_for_equivalent_params():
    """case.setup's `_write_slurm_header(f)` (which reads from `par`) must
    produce the SAME header machines.write_slurm_header produces for the
    equivalent explicit arguments -- a behavioural check (rule 10a), not a
    pin on case.setup's source text, so either side may be refactored freely
    as long as the two stay in lockstep."""
    class FakePar:
        casename = 'tpv-test'
        HPC_nnode = 3
        HPC_ncpu = 384
        HPC_queue = 'normal'
        HPC_time = '00:20:00'
        HPC_account = 'EAR22013'
        HPC_email = ''

    direct = io.StringIO()
    machines.write_slurm_header(
        direct, jobname=FakePar.casename, nnode=FakePar.HPC_nnode,
        ncpu=FakePar.HPC_ncpu, queue=FakePar.HPC_queue,
        walltime=FakePar.HPC_time, account=FakePar.HPC_account,
        email=FakePar.HPC_email)

    import importlib.util
    from importlib.machinery import SourceFileLoader
    import os
    case_setup_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        'scripts', 'case.setup')
    # case.setup has no .py suffix, so spec_from_file_location cannot infer a
    # loader from the extension -- give it one explicitly (it IS python source).
    loader = SourceFileLoader('case_setup_under_test', case_setup_path)
    spec = importlib.util.spec_from_file_location('case_setup_under_test',
                                                   case_setup_path, loader=loader)
    # case.setup imports `from user_defined_params import par` at module
    # level, which this fixture does not have -- stub it out via sys.modules
    # before exec, since we only need `_write_slurm_header`.
    import sys
    import types
    fake_udp = types.ModuleType('user_defined_params')
    fake_udp.par = FakePar
    sys.modules['user_defined_params'] = fake_udp
    try:
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        via_case_setup = io.StringIO()
        mod._write_slurm_header(via_case_setup)
    finally:
        del sys.modules['user_defined_params']

    assert via_case_setup.getvalue() == direct.getvalue()
