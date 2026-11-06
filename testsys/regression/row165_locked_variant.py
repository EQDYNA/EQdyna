#! /usr/bin/env python3
"""Row 165 helpers: build a synthetic locked-node variant of the test.restart2
fixture, run it on fortran (4 ranks) and python-jax/numpy (serial), and compare
the per-node frt output.  Imported by test_row165_locked_variant_parity.py and
by the scratch investigation scripts.

The committed fixture is a crop ABOVE the source's locked band (README.md of
case_input/test.restart2).  The variant overwrites the first LOCKED_ROWS dip rows
(z = -4 km ... ) with a locked state: slip_rate = LOCKED_V and the six node
velocities (vxm..vzs) = 0, leaving the stresses and state_variable as stored.
"""
import os
import shutil
import subprocess
import sys

import netCDF4
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, 'testsys'))
import frt_canonical  # noqa: E402

CASE = 'test.restart2'
FIXTURE = os.path.join(REPO, 'case_input', CASE, 'fault.r.nc')
VEL_VARS = ('vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs')


def make_variant(dst, rows=range(0, 4), locked_v=1e-30, zero_vel=True):
    """Copy the committed fixture to dst, locking the given dip rows."""
    shutil.copyfile(FIXTURE, dst)
    ds = netCDF4.Dataset(dst, 'r+')
    try:
        for r in rows:
            ds['slip_rate'][r, :] = locked_v
            if zero_vel:
                for v in VEL_VARS:
                    ds[v][r, :] = 0.0
    finally:
        ds.close()


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = REPO
    env['PATH'] = os.pathsep.join([os.path.join(REPO, 'bin'), os.path.join(REPO, 'scripts'), env['PATH']])
    env['OMP_NUM_THREADS'] = '2'
    return env


def stage(root, name, fixture, term, extra_params=''):
    """create.newcase + term override + fixture swap; returns case_dir."""
    env = _env()
    os.makedirs(root, exist_ok=True)
    case_dir = os.path.join(root, name)
    if os.path.isdir(case_dir):
        raise RuntimeError('refusing to reuse %s' % case_dir)
    subprocess.run(['create.newcase', name, CASE], cwd=root, env=env, check=True,
                   stdout=subprocess.DEVNULL)
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'a') as f:
        f.write('\npar.term = %r\n%s\n' % (term, extra_params))
    shutil.copyfile(fixture, os.path.join(case_dir, 'fault.r.nc'))
    return case_dir


def run_fortran(case_dir, ranks=4, timeout=600, bin_env=None):
    env = _env()
    env.update(bin_env or {})
    subprocess.run(['./case.setup'], cwd=case_dir, env=env, check=True, stdout=subprocess.DEVNULL, timeout=timeout)
    subprocess.run(['mpirun', '-np', str(ranks), 'eqdyna'], cwd=case_dir, env=env, check=True,
                   stdout=open(os.path.join(case_dir, 'run.out'), 'w'), stderr=subprocess.STDOUT, timeout=timeout)
    return frt_canonical.canonical_from_case(case_dir)


def run_python(case_dir, backend='jax', timeout=900, extra_env=None):
    env = _env()
    env['PYTHONPATH'] = os.path.join(REPO, 'src', 'python')
    env['JAX_PLATFORMS'] = 'cpu'
    env.update(extra_env or {})
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'a') as f:
        f.write('\npar.nx = 1\npar.ny = 1\npar.nz = 1\n')
    subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env, check=True, stdout=subprocess.DEVNULL, timeout=timeout)
    subprocess.run([sys.executable, '-u', '-m', 'eqdyna', case_dir, '--backend', backend, '--device', 'cpu'],
                   cwd=REPO, env=env, check=True, stdout=open(os.path.join(case_dir, 'run.out'), 'w'),
                   stderr=subprocess.STDOUT, timeout=timeout)
    return frt_canonical.canonical_from_case(case_dir)


def diff_report(a, b, label):
    """Per-column max|diff| of two canonical frt arrays, with the node of the worst row."""
    if a.shape != b.shape:
        raise RuntimeError('%s: shapes %s vs %s' % (label, a.shape, b.shape))
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise RuntimeError('%s: non-finite values (nan/inf) in output' % label)
    d = np.abs(a - b)
    out = []
    for c in range(3, a.shape[1]):
        i = int(d[:, c].argmax())
        out.append((c, d[i, c], i, a[i, 0:3], a[i, c], b[i, c], int((d[:, c] > 0).sum())))
    return out
