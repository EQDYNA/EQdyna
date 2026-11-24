"""Row 167 step 3: map the elements rewired by PR #281 onto the changed reference rows.
usage: map_rows.py <pre.npz> <post.npz> <frt_old> <frt_new> [thr_Pa]
pre/post.npz from dump_mesh.py run at the tree before / after #281 (same case).
Prints: rewired element count, meshCoor equality, fault rows touched by their nodes (E),
rows whose stress cols 11-13 changed by > thr (R), R-E and E-R. R-E must be 0."""
import sys
import numpy as np
pre, post, fo, fn = sys.argv[1:5]
thr = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-3
a, b = np.load(pre), np.load(post)
X, ca, cb = a['meshCoor'], a['conn'], b['conn']
ch = np.where((ca != cb).any(1))[0]
print('rewired elements', len(ch), '| meshCoor equal', np.array_equal(a['meshCoor'], b['meshCoor']),
      '| elem_type equal', np.array_equal(a['elem_type'], b['elem_type']), '| nsmp equal', np.array_equal(a['nsmp'], b['nsmp']))
old, new = np.loadtxt(fo), np.loadtxt(fn)
assert old.shape == new.shape and np.array_equal(old[:, :3], new[:, :3]), 'frt grid differs'
d = np.abs(new - old)
print('frt rows', len(old), '| cols with any change', np.where(d.any(0))[0].tolist(), '| max stress diff', d[:, 11:14].max())
key = lambda p: tuple(np.round(p, 0).astype(int))
rk = {key(r[:3]): i for i, r in enumerate(old)}
nodes = np.unique(np.concatenate([ca[ch].ravel(), cb[ch].ravel()]))
E = {rk[key(X[n])] for n in nodes if key(X[n]) in rk}
R = set(np.where(d[:, 11:14].max(1) > thr)[0].tolist())
print('E (fault rows touched)', len(E), '| R (stress changed >%g)' % thr, len(R), '| R-E', len(R - E), '| E-R', len(E - R))
assert len(ch) > 0 and len(R) > 0 and not (R - E), 'a changed row lies outside the rewired elements'
