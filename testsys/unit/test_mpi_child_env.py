"""Every testsys mpirun child gets HWLOC_COMPONENTS excluded and no DISPLAY.

hwloc GL/OpenCL/CUDA/NVML discovery hangs mpirun at 0% CPU on some boxes.
The helper is `testsys.common.mpi_child_env`; every launch point below must
route its child env through it. Both directions are guarded: the property holds
on the real launcher env, and a launcher that stops calling the helper (or a new
unlisted mpirun launcher) fails the source checks.
"""
import os
import re

import pytest

from testsys import common

REPO = common.REPO_ROOT
EXCL = '-gl,-opencl,-cuda,-nvml'

# file that builds an env for an mpirun child -> number of env builders using the helper
LAUNCHERS = {
    'testsys/e2e/run_e2e.py': 1,            # base_env(): fortran + python-jax-mpi cells
    'testsys/e2e/run_e2e_full.py': 1,
    'testsys/perf/run_perf.py': 1,
    'testsys/perf/hpc_scaling_suite.py': 1,
    'testsys/perf/run_mpi_scaling.py': 1,
    'testsys/perf/run_scaling.py': 1,
    'testsys/perf/profile_overhead.py': 2,
}


def test_default_excludes_hwloc_components_and_drops_display():
    out = common.mpi_child_env({'DISPLAY': ':0', 'PATH': '/x'})
    assert (out.get('HWLOC_COMPONENTS'), 'DISPLAY' in out, out['PATH']) == (EXCL, False, '/x')


def test_user_set_hwloc_components_is_never_overridden():
    assert common.mpi_child_env({'HWLOC_COMPONENTS': 'x86'})['HWLOC_COMPONENTS'] == 'x86'


def test_caller_environment_is_not_mutated(monkeypatch):
    monkeypatch.setenv('DISPLAY', ':7')
    monkeypatch.delenv('HWLOC_COMPONENTS', raising=False)
    common.mpi_child_env(os.environ)
    assert (os.environ.get('DISPLAY'), 'HWLOC_COMPONENTS' in os.environ) == (':7', False)


def test_absent_display_is_not_an_error():
    assert 'DISPLAY' not in common.mpi_child_env({})


def test_run_e2e_base_env_applies_it(monkeypatch):
    monkeypatch.setenv('DISPLAY', ':7')
    monkeypatch.delenv('HWLOC_COMPONENTS', raising=False)
    from testsys.e2e import run_e2e
    env = run_e2e.base_env()
    assert (env.get('HWLOC_COMPONENTS'), 'DISPLAY' in env) == (EXCL, False)


def test_run_e2e_base_env_respects_user_value(monkeypatch):
    monkeypatch.setenv('HWLOC_COMPONENTS', 'x86')
    from testsys.e2e import run_e2e
    assert run_e2e.base_env()['HWLOC_COMPONENTS'] == 'x86'


@pytest.mark.parametrize('rel,n', sorted(LAUNCHERS.items()))
def test_launcher_routes_child_env_through_helper(rel, n):
    src = open(os.path.join(REPO, rel)).read()
    calls = len(re.findall(r'^\s*(?:env\s*=\s*|return\s+)(?:common\.)?mpi_child_env\(', src, re.M))
    assert calls == n, '%s: %d mpi_child_env env builders, expected %d' % (rel, calls, n)


def test_every_mpirun_launching_testsys_file_is_listed():
    """A new file that launches mpirun must be added to LAUNCHERS, not skipped."""
    pat = re.compile(r"""(?:MPIRUN|mpirun|'mpirun'|"mpirun")\s*,\s*['"]-|mpirun --bind-to none -np|%s --bind-to none""")
    found = set()
    for sub in ('e2e', 'perf'):
        d = os.path.join(REPO, 'testsys', sub)
        for f in os.listdir(d):
            if f.endswith('.py') and pat.search(open(os.path.join(d, f)).read()):
                found.add('testsys/%s/%s' % (sub, f))
    assert found == set(LAUNCHERS), 'launchers not in LAUNCHERS: %s; stale: %s' % (
        sorted(found - set(LAUNCHERS)), sorted(set(LAUNCHERS) - found))


# files that pass `env=mpi_child_env(...)` inline on their own mpirun subprocess call
INLINE_LAUNCHERS = ('testsys/regression/mpirun_capture.py',
                    'testsys/regression/test_dipping_fault_y_split.py',
                    'testsys/regression/test_equilibrium_dump.py',
                    'testsys/regression/test_fault_mpi_boundary_arn.py',
                    'testsys/regression/test_fault_time_series_output.py')

# regression/parity launchers that STILL bypass the helper (debt; fixing one means
# moving it to INLINE_LAUNCHERS, which this test forces: a fixed file left here is stale)
KNOWN_BYPASS = (
    'testsys/parity/evidence_c_degen_port.py',
    'testsys/parity/evidence_drv_a6_chaos.py',
    'testsys/parity/evidence_wedge_kernel.py',
    'testsys/regression/row165_locked_variant.py',
    'testsys/regression/test_mesh_screening.py',
    'testsys/regression/test_offfault_station_depth_selection.py',
    'testsys/regression/test_pml_element_type_grid_centre.py',
    'testsys/regression/test_row120_mpi_station_output.py',
    'testsys/regression/test_row127_station_ownership.py',
    'testsys/regression/test_stop_exit_status.py',
    'testsys/regression/test_unknown_cli_arg.py',
)
_CALL = re.compile(r"subprocess\.run\(\s*\[[^\]]*\b(?:mpirun|MPIRUN|mpirun_bin)\b")


def _bypass_scan():
    """regression/parity files with a subprocess.run mpirun call and no inline helper."""
    raw, ok = set(), set()
    for sub in ('regression', 'parity'):
        d = os.path.join(REPO, 'testsys', sub)
        for f in os.listdir(d):
            if not f.endswith('.py'):
                continue
            src = open(os.path.join(d, f)).read()
            rel = 'testsys/%s/%s' % (sub, f)
            if _CALL.search(src) or f == 'mpirun_capture.py':
                (ok if re.search(r'env=mpi_child_env\(', src) else raw).add(rel)
    return raw, ok


@pytest.mark.parametrize('rel', INLINE_LAUNCHERS)
def test_inline_launcher_passes_helper_env_to_its_mpirun_call(rel):
    src = open(os.path.join(REPO, rel)).read()
    assert re.search(r'env=mpi_child_env\(', src), '%s launches mpirun without mpi_child_env' % rel


def test_no_unlisted_regression_or_parity_launcher_bypasses_helper():
    raw, _ = _bypass_scan()
    assert raw - set(KNOWN_BYPASS) == set(), 'new mpirun launcher bypasses mpi_child_env: %s' % sorted(raw - set(KNOWN_BYPASS))


def test_known_bypass_list_has_no_stale_entries():
    raw, _ = _bypass_scan()
    stale = set(KNOWN_BYPASS) - raw
    assert stale == set(), 'fixed (or vanished) launchers still in KNOWN_BYPASS: %s' % sorted(stale)


def _fake_mpirun(tmp_path):
    p = tmp_path / 'fakempirun'
    p.write_text('#!/bin/sh\n'
                 'echo "HWLOC_COMPONENTS=${HWLOC_COMPONENTS-<unset>}"\n'
                 'echo "DISPLAY=${DISPLAY-<unset>}"\n'
                 'shift 2\n'      # drop `-np N`
                 'exec "$@"\n')
    p.chmod(0o755)
    return str(p)


@pytest.mark.parametrize('caller_env', ['inherit', 'explicit'])
def test_run_rank_files_child_gets_exclusion_and_no_display(tmp_path, monkeypatch, caller_env):
    from testsys.regression import mpirun_capture
    monkeypatch.setenv('DISPLAY', ':7')
    monkeypatch.delenv('HWLOC_COMPONENTS', raising=False)
    env = None if caller_env == 'inherit' else dict(os.environ)
    work = tmp_path / 'w'
    work.mkdir()
    rc, text = mpirun_capture.run_rank_files(_fake_mpirun(tmp_path), 'true', str(work), env=env)
    assert (rc, 'HWLOC_COMPONENTS=%s' % EXCL in text, 'DISPLAY=<unset>' in text) == (0, True, True), text
