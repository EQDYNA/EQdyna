import numpy as np
import subprocess

import os
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

def load_old(case, commit='9cac0fd'):
    path = 'test.reference.results/%s/frt.canonical.txt' % case
    r = subprocess.run(['git', 'show', '%s^:%s' % (commit, path)], cwd=REPO,
                        capture_output=True, text=True, check=True)
    import io
    return np.loadtxt(io.StringIO(r.stdout))

def load_new(case):
    return np.loadtxt('%s/test.reference.results/%s/frt.canonical.txt' % (REPO, case))

SENTINEL = 1.0e4
FNFT_COL = 3

for case in ['test.tpv13', 'test.tpv27', 'test.tpv30', 'test.drv.a6']:
    old = load_old(case)
    new = load_new(case)
    assert old.shape == new.shape, case
    assert np.max(np.abs(old[:, :3] - new[:, :3])) == 0.0, 'coords not aligned: %s' % case
    z = old[:, 2]
    for label, mask in (('deep', z < -2000.0), ('shallow', z >= -2000.0)):
        ro = old[mask, FNFT_COL] < SENTINEL
        rn = new[mask, FNFT_COL] < SENTINEL
        frac_old = ro.mean()
        frac_new = rn.mean()
        both = ro & rn
        dt = np.abs(old[mask, FNFT_COL][both] - new[mask, FNFT_COL][both])
        dt_med = np.median(dt) if dt.size else float('nan')
        dt_p90 = np.percentile(dt, 90) if dt.size else float('nan')
        print('%-14s %-8s n=%5d ruptured_old=%d ruptured_new=%d frac_old=%.4f frac_new=%.4f '
              'Delta_frac=%+.4f  dt_med=%.4f dt_p90=%.4f'
              % (case, label, mask.sum(), ro.sum(), rn.sum(), frac_old, frac_new,
                 frac_new - frac_old, dt_med, dt_p90))
