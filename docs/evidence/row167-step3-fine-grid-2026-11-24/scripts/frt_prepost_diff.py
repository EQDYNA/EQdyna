"""Pre/post final-stress and rupture-time diff of two run dirs (frt.txt* concatenated, deduped on rounded x,y,z).
usage: python3 -I frt_prepost_diff.py <pre_run_dir> <post_run_dir>
Columns (library_output.f90): 0-2 xyz, 3 rupture time (99999 = not ruptured), 4-6 final slips, 11-13 final tnrm/tstk/tdip."""
import glob, sys
import numpy as np
def load(d):
    rows = []
    for f in glob.glob(d + '/frt.txt*'):
        rows += [[float(x) for x in l.split()] for l in open(f) if l.strip()]
    a = np.array(rows); a[:, :3] = np.round(a[:, :3], 1)
    a = a[np.lexsort((a[:, 2], a[:, 1], a[:, 0]))]
    _, i = np.unique(a[:, :3], axis=0, return_index=True)
    return a[np.sort(i)]
a, b = load(sys.argv[1]), load(sys.argv[2])
assert a.shape == b.shape and np.array_equal(a[:, :3], b[:, :3]), 'node sets differ'
d = np.abs(a - b)
print('nodes', a.shape[0], 'max diff by col:', ' '.join('%d:%.3g' % (i, d[:, i].max()) for i in (3, 4, 5, 6, 11, 12, 13)))
m = d[:, 11] > 1e3
print('rows with |dtnrm|>1 kPa:', int(m.sum()), 'z range', (a[m, 2].min(), a[m, 2].max()) if m.any() else None)
