import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from eqdyna.port import load
from eqdyna import port_jax as pj

case_dir = sys.argv[1]
nsteps = int(sys.argv[2]) if len(sys.argv) > 2 else None
golden_path = sys.argv[3] if len(sys.argv) > 3 else None

t0 = time.time()
S = load(case_dir)
t1 = time.time()
print('load time', t1 - t0, 's; N=', S['N'], 'E=', S['E'], 'nftnd=', S['nftnd'])

out = pj.run(S, nsteps=nsteps, verbose=False)
t2 = time.time()
print('jax compile+solve time', t2 - t1, 's for', nsteps or S['nstep'], 'steps')

fric = out['fric']
nsmp1 = S['nsmp1']
coords = S['meshCoor'][nsmp1]
cols = np.column_stack([
    coords[:, 0], coords[:, 1], coords[:, 2],
    out['fnft'],
    fric[:, 70], fric[:, 71], fric[:, 72], fric[:, 73], fric[:, 74], fric[:, 75],
    fric[:, 46],
    fric[:, 77], fric[:, 78], fric[:, 79],
    fric[:, 30], fric[:, 31], fric[:, 32],
    fric[:, 33], fric[:, 34], fric[:, 35],
    fric[:, 19], fric[:, 22],
])
py_out_path = os.path.join(case_dir, 'frt_jax.txt')
np.savetxt(py_out_path, cols, fmt='%18.7e')
print('wrote', py_out_path)

if golden_path:
    golden = np.loadtxt(golden_path)
    labels = ['x', 'y', 'z', 'fnft(rupt.time)', 'slipS', 'slipD', 'slipN', 'srS', 'srD',
              'peakSr', 'finalSr(47,unused)', 'Tn(78)', 'Ts(79)', 'Td(80)',
              'vxm(31)', 'vym(32)', 'vzm(33)', 'vxs(34)', 'vys(35)', 'vzs(36)',
              'rsfstate(20,unused)', 'thetaPc(23,unused)']
    print(f"{'col':22s} {'max_abs_diff':>16s} {'max_rel_diff':>16s}")
    max_abs_overall = 0.0
    for i, lab in enumerate(labels):
        d = np.abs(cols[:, i] - golden[:, i])
        denom = np.maximum(np.abs(golden[:, i]), 1e-300)
        rel = d / denom
        max_abs_overall = max(max_abs_overall, d.max())
        print(f"{lab:22s} {d.max():16.6e} {rel.max():16.6e}")
    print('overall max abs diff:', max_abs_overall)
