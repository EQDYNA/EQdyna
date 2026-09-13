import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'eqdyna'))
from eqdyna.port import load

engine = sys.argv[1]  # 'numpy' or 'jax'
case_dir = sys.argv[2]
nsteps = int(sys.argv[3]) if len(sys.argv) > 3 else None
golden_path = sys.argv[4] if len(sys.argv) > 4 else None

if engine == 'numpy':
    from eqdyna import port_tp as mod
elif engine == 'jax':
    from eqdyna import port_tp_jax as mod
else:
    raise ValueError(engine)

t0 = time.time()
S = load(case_dir)
t1 = time.time()
print('load time', t1 - t0, 's; N=', S['N'], 'E=', S['E'], 'nftnd=', S['nftnd'])

with np.errstate(divide='ignore', invalid='ignore'):
    out = mod.run(S, nsteps=nsteps, verbose=False)
t2 = time.time()
print(engine, 'solve time', t2 - t1, 's for', nsteps or S['nstep'], 'steps')

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
out_path = os.path.join(case_dir, f'frt_{engine}.txt')
np.savetxt(out_path, cols, fmt='%18.7e')
print('wrote', out_path)

if golden_path:
    golden = np.loadtxt(golden_path)
    labels = ['x', 'y', 'z', 'fnft(rupt.time)', 'slipS', 'slipD', 'slipN', 'srS', 'srD',
              'peakSr', 'finalSr(47)', 'Tn(78)', 'Ts(79)', 'Td(80)',
              'vxm(31)', 'vym(32)', 'vzm(33)', 'vxs(34)', 'vys(35)', 'vzs(36)',
              'rsfstate(20)', 'thetaPc(23)']
    print(f"{'col':22s} {'max_abs_diff':>16s} {'max_rel_diff':>16s}")
    max_abs_overall = 0.0
    for i, lab in enumerate(labels):
        d = np.abs(cols[:, i] - golden[:, i])
        denom = np.maximum(np.abs(golden[:, i]), 1e-300)
        rel = d / denom
        max_abs_overall = max(max_abs_overall, d.max())
        print(f"{lab:22s} {d.max():16.6e} {rel.max():16.6e}")
    print('overall max abs diff:', max_abs_overall)
