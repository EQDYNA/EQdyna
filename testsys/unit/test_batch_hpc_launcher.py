"""batch.hpc (the LS6 per-case SLURM script case.setup writes) launches the
solver with the launcher scripts/machines.py records for ls6, never with an
`ibrun -np N` line (2026-10-01: TACC's ibrun takes -n/-o, not -np; the
testsys --submit sweep already moved to `mpirun -np N`, measured on LS6
2026-09-29, while case.setup kept writing the invalid ibrun call)."""
import os
import re
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from machines import machine  # noqa: E402


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'scripts'), env.get('PATH', '')])
    return env


def test_batch_hpc_uses_the_registry_launcher_with_the_run_sh_rank_count(tmp_path):
    case_dir = str(tmp_path / 'case')
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, 'test.tpv8'], env=_env(), capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=_env(),
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]

    run_sh = open(os.path.join(case_dir, 'run.sh')).read()
    batch = open(os.path.join(case_dir, 'batch.hpc')).read()
    ranks = re.search(r'mpirun -np (\d+) eqdyna', run_sh).group(1)

    launcher = machine('ls6')['mpirun']
    assert '%s -np %s eqdyna' % (launcher, ranks) in batch, batch
    assert not re.search(r'\bibrun\b.*-np\b', batch), batch
