"""scripts/case.submit: refuses a batch.hpc with no project account, and
otherwise hands it to sbatch and returns sbatch's status. sbatch is a stub
on PATH that records its arguments, so nothing is ever really submitted."""
import os
import subprocess
import sys

from conftest import REPO_ROOT

SUBMIT = os.path.join(REPO_ROOT, 'scripts', 'case.submit')


def _run(tmp_path, account, sbatch_rc=0):
    (tmp_path / 'batch.hpc').write_text(
        '#!/bin/bash\n#SBATCH -J t\n#SBATCH -A %s\n  ibrun eqdyna\n' % account)
    bindir = tmp_path / 'bin'; bindir.mkdir(exist_ok=True)
    log = tmp_path / 'sbatch.called'
    (bindir / 'sbatch').write_text('#!/bin/sh\necho "$@" > %s\nexit %d\n' % (log, sbatch_rc))
    os.chmod(bindir / 'sbatch', 0o755)
    env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ['PATH'])
    r = subprocess.run([sys.executable, SUBMIT], cwd=tmp_path, env=env,
                       capture_output=True, text=True)
    return r, log


def test_empty_account_is_refused_before_sbatch(tmp_path):
    r, log = _run(tmp_path, '')
    assert r.returncode == 1 and 'no project account' in r.stderr
    assert not log.exists()


def test_account_set_calls_sbatch_and_returns_its_status(tmp_path):
    r, log = _run(tmp_path, 'EAR24033', sbatch_rc=0)
    assert r.returncode == 0 and log.read_text().strip() == 'batch.hpc'
    r, _ = _run(tmp_path, 'EAR24033', sbatch_rc=3)
    assert r.returncode == 3


def test_missing_batch_hpc_is_refused(tmp_path):
    r = subprocess.run([sys.executable, SUBMIT], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 1 and 'run ./case.setup first' in r.stderr
