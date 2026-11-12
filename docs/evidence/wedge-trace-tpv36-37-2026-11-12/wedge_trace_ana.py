#! /usr/bin/env python3
"""Trace-node (faultst*dp000) behaviour of the wedge-degenerated TPV36/37 runs vs EQdyna's own 50 m archive.
usage: wedge_trace_ana.py RUNS_DIR ARCHIVE_DIR REFERENCE_DIR OUT_DIR T_EVAL
RUNS_DIR/<t36|t37>_dx<dx>/ are runs; ARCHIVE_DIR/<tpv36|tpv37>/eqdyna-v5.3.3-<50m|100m>-2024-nstress-corrected are read-only
(n-stress negative = compression, same sign convention as the runs); REFERENCE_DIR = test.reference.results (the gate reference,
500 m, 5 s). Station columns: t, h-slip, h-rate, h-shear, v-slip, v-rate, v-shear, n-stress."""
import csv, glob, os, re, sys
import numpy as np
runs, arch, ref, out, teval = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])
os.makedirs(out, exist_ok=True)

def fl(x):
    # Fortran E format drops the letter for 3-digit exponents: 0.14-104 means 0.14E-104
    if 'E' not in x and re.search(r'\d[-+]\d{3}$', x):
        x = x[:-4] + 'E' + x[-4:]
    return float(x)

def load(p):
    rows = []
    for l in open(p):
        s = l.strip()
        if not s or s.startswith('#') or s.startswith('t '):
            continue
        rows.append([fl(x) for x in s.split()])
    return np.array(rows)

STN = ['000dp000', '000dp010', '000dp030', '000dp060', '000dp120', '000dp180', '040dp000', '080dp000', '120dp000']
def metrics(a, te):
    sr = np.hypot(a[:, 2], a[:, 5])
    sl = np.hypot(a[:, 1], a[:, 4])
    k = np.argmax(sr > 1e-3) if sr.max() > 1e-3 else -1
    i = np.searchsorted(a[:, 0], te, side='right') - 1
    m = a[:i + 1]
    return dict(arrival_s=(a[k, 0] if k >= 0 else float('nan')), slip_te=sl[i], n_te=a[i, 7],
                n_min=m[:, 7].min(), n_max=m[:, 7].max(), tmax_avail=a[-1, 0])

rows = []
srcs = []
for d in sorted([g for g in glob.glob(runs + '/t3[67]_dx*/') if re.search(r'_dx\d+/$', g)]):
    nm = os.path.basename(d.rstrip('/'))
    if not os.path.exists(d + 'DONE'):
        continue
    case = 'tpv' + nm[1:3]
    dx = float(re.search(r'dx(\d+)', nm).group(1))
    srcs.append((case, 'wedge run', dx, d + 'faultst'))
for case in ('tpv36', 'tpv37'):
    for tag, dx in (('100m', 100.), ('50m', 50.)):
        d = '%s/%s/eqdyna-v5.3.3-%s-2024-nstress-corrected/' % (arch, case, tag)
        if os.path.isdir(d):
            srcs.append((case, 'archive', dx, d + 'faultst'))
data = {}
for case, kind, dx, pre in srcs:
    for s in STN:
        p = pre + s + ('.txt' if kind == 'wedge run' else '')
        if not os.path.exists(p):
            continue
        a = load(p)
        data[(case, kind, dx, s)] = a
        m = metrics(a, teval)
        rows.append([case, kind, dx, s] + [m[k] for k in ('arrival_s', 'slip_te', 'n_te', 'n_min', 'n_max', 'tmax_avail')])
with open(out + '/trace_summary.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['case', 'source', 'dx_m', 'station', 'arrival_s', 'slip_at_t_eval_m', 'n_at_t_eval_MPa', 'n_min_MPa', 'n_max_MPa', 'tmax_s'])
    for r in rows:
        w.writerow([r[0], r[1], '%g' % r[2], r[3]] + ['%.5g' % x for x in r[4:]])
print('summary rows', len(rows), 'teval', teval)

# control: 500 m run vs gate reference, t <= 5 s, on the reference stations
ctl = []
for case in ('tpv36', 'tpv37'):
    for s in ('000dp180', '000dp010', '120dp030'):
        k = (case, 'wedge run', 500., s)
        pr = '%s/%s/stations/faultst%s.txt' % (ref, 'test.' + case, s)
        if k not in data or not os.path.exists(pr):
            continue
        a, b = data[k], load(pr)
        t = b[b[:, 0] <= 5.0, 0]
        for col, nmx in ((1, 'h-slip'), (4, 'v-slip'), (7, 'n-stress')):
            ya = np.interp(t, a[:, 0], a[:, col]); yb = np.interp(t, b[:, 0], b[:, col])
            ctl.append([case, s, nmx, np.abs(ya - yb).max() / max(np.abs(yb).max(), 1e-6), np.abs(yb).max(), len(t)])
with open(out + '/control_vs_gate_reference.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['case', 'station', 'quantity', 'max_abs_diff_over_ref_max', 'ref_max_abs', 'n_samples_le_5s'])
    for r in ctl:
        w.writerow([r[0], r[1], r[2], '%.3e' % r[3], '%.5g' % r[4], r[5]])
print('control rows', len(ctl))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
for case in ('tpv36', 'tpv37'):
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    for (c, kind, dx, s), a in sorted(data.items(), key=lambda kv: (kv[0][1], -kv[0][2])):
        if c != case:
            continue
        lab = '%s %g m' % ('wedge run' if kind == 'wedge run' else 'archive', dx)
        ls = '-' if kind == 'wedge run' else '--'
        if s == '000dp000':
            ax[0].plot(a[:, 0], a[:, 7], ls, label=lab)
        if s == '000dp010':
            ax[1].plot(a[:, 0], a[:, 7], ls, label=lab)
        if s == '000dp000':
            ax[2].plot(a[:, 0], np.hypot(a[:, 1], a[:, 4]), ls, label=lab)
    for x, t in zip(ax, ('%s trace n-stress (dp000)' % case, 'dp010 n-stress', 'dp000 slip')):
        x.set_title(t); x.set_xlabel('t (s)'); x.set_xlim(0, teval); x.grid(alpha=.3)
    ax[0].set_ylabel('MPa (neg = compression)'); ax[2].set_ylabel('m'); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig('%s/fig_%s_trace.png' % (out, case), dpi=130); plt.close(fig)
