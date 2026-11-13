#! /usr/bin/env python3
"""Probe 2 / guard-rerun analysis. usage: dip_ana.py <sweep_dir> <out_dir> [guard_dir]
sweep_dir/d<dip>_dx<dx>/faultst<STN>.txt (columns t, h-slip, h-rate, h-shear, v-slip, v-rate, v-shear, n-stress; n negative = compression).
Writes dip_sweep.csv, fig_probe2_nstress_vs_dip.png; with guard_dir also guard_tpv12.csv (dp000 n-stress at 4 s)."""
import csv, glob, os, re, sys
import numpy as np
sweep, out = sys.argv[1], sys.argv[2]
guard = sys.argv[3] if len(sys.argv) > 3 else None
THRESH = 1.0      # MPa, fixed before the first run (README / notes)
def fl(x):
    if 'E' not in x and re.search(r'\d[-+]\d{3}$', x):
        x = x[:-4] + 'E' + x[-4:]
    return float(x)
def load(p):
    rows = []
    for l in open(p):
        s = l.strip()
        if s and not s.startswith('#') and not s.startswith('t '):
            rows.append([fl(x) for x in s.split()])
    return np.array(rows)
def at(a, te):
    i = np.searchsorted(a[:, 0], te, side='right') - 1
    return i
res = {}
rows = []
for d in sorted(glob.glob(sweep + '/d*_dx*/')):
    m = re.search(r'/d(\d+)_dx(\d+)/$', d)
    if not m or not os.path.exists(d + 'DONE'):
        continue
    dip, dx = int(m.group(1)), int(m.group(2))
    a = load(d + 'faultst000dp000.txt'); b = load(d + 'faultst000dp010.txt'); c = load(d + 'faultst000dp030.txt')
    sl = lambda x: np.hypot(x[:, 1], x[:, 4])
    sr = np.hypot(a[:, 2], a[:, 5]); k = np.argmax(sr > 1e-3) if sr.max() > 1e-3 else -1
    i5, ie = at(a, 5.0), len(a) - 1
    r = dict(dip=dip, dx=dx, nsamples=len(a), t_end=a[-1, 0], arrival_trace_s=(a[k, 0] if k >= 0 else float('nan')),
             n_at_5s=a[i5, 7], n_at_end=a[ie, 7], n_max=a[:, 7].max(), n_min=a[:, 7].min(),
             slip_trace_end=sl(a)[ie], slip_1km_end=sl(b)[ie], slip_3km_end=sl(c)[ie],
             ratio_trace_to_1km=sl(a)[ie] / max(sl(b)[ie], 1e-12), n_1km_max=b[:, 7].max())
    rows.append(r); res[(dip, dx)] = r
with open(out + '/dip_sweep.csv', 'w', newline='') as f:
    keys = list(rows[0]); w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow({k: ('%.6g' % v if isinstance(v, float) else v) for k, v in r.items()})
print('sweep rows', len(rows))
dxs = sorted({r['dx'] for r in rows}); dips = sorted({r['dip'] for r in rows})
for dip in dips:
    rr = [res[(dip, dx)] for dx in dxs if (dip, dx) in res]
    if len(rr) < 2:
        print('dip %d: only %d grid(s) done' % (dip, len(rr))); continue
    both = all(r['n_max'] > THRESH for r in rr); grow = rr[-1]['n_max'] > rr[0]['n_max']
    print('dip %d: n_max %s MPa (dx %s) -> TENSION APPEARS = %s (both>%.1f: %s, finer larger: %s)' %
          (dip, [round(r['n_max'], 4) for r in rr], dxs, both and grow, THRESH, both, grow))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for dx in dxs:
    pts = [(r['dip'], r['n_max'], r['n_at_end'], r['ratio_trace_to_1km']) for r in rows if r['dx'] == dx]
    pts.sort()
    ax[0].plot([p[0] for p in pts], [p[1] for p in pts], 'o-', label='max over time, dx=%d m' % dx)
    ax[0].plot([p[0] for p in pts], [p[2] for p in pts], 's--', label='at term, dx=%d m' % dx)
    ax[1].plot([p[0] for p in pts], [p[3] for p in pts], 'o-', label='dx=%d m' % dx)
ax[0].axhline(THRESH, color='r', ls=':', label='+%.1f MPa criterion' % THRESH); ax[0].axhline(0, color='k', lw=0.5)
ax[0].set_xlabel('dip (deg)'); ax[0].set_ylabel('trace (dp000) n-stress (MPa, + = tension)'); ax[0].legend(fontsize=7)
ax[0].set_title('test.tpv36, dip changed only: trace n-stress')
ax[1].set_xlabel('dip (deg)'); ax[1].set_ylabel('slip(trace) / slip(1 km down dip) at term'); ax[1].legend(fontsize=7); ax[1].set_title('trace slip loss')
fig.tight_layout(); fig.savefig(out + '/fig_probe2_nstress_vs_dip.png', dpi=130)
if guard:
    g = []
    for d in sorted(glob.glob(guard + '/gw_dx*/')):
        if not os.path.exists(d + 'DONE'): continue
        dx = int(re.search(r'dx(\d+)', d).group(1)); a = load(d + 'faultst000dp000.txt'); i = at(a, 4.0)
        g.append([dx, len(a), a[i, 0], a[i, 7], a[:, 7].max(), np.hypot(a[i, 1], a[i, 4])])
    with open(out + '/guard_tpv12.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['dx_m', 'nsamples', 't_eval_s', 'n_at_4s_MPa', 'n_max_MPa', 'slip_trace_at_4s_m'])
        for r in g: w.writerow(['%g' % x for x in r])
    print('guard rows', g)
