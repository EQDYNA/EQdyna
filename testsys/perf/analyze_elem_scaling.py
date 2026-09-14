#! /usr/bin/env python3
"""
Fit and plot seconds-per-timestep vs elements-per-rank for EQdyna.

Reads elem_scaling_last.json (written by run_elem_scaling.py) and answers one
question: is per-step cost a straight line through ~zero in elements/rank,
independent of how those elements arrive?

Outputs a table, an ordinary-least-squares fit (slope, intercept, R^2,
residuals), per-rank-count sub-fits (the memory-bandwidth test), and
elem_scaling.png.

Usage: python3 testsys/perf/analyze_elem_scaling.py [--json f.json] [--png f.png]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def ols(x, y):
    """y = a + b x. Returns slope, intercept, R^2, residuals, stderr(slope)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    b, a = np.polyfit(x, y, 1)
    fit = a + b * x
    res = y - fit
    ss_res = float((res ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot else float('nan')
    se_b = float(np.sqrt(ss_res / (n - 2) / ((x - x.mean()) ** 2).sum())) if n > 2 else float('nan')
    return b, a, r2, res, se_b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json', default=os.path.join(HERE, 'elem_scaling_last.json'))
    ap.add_argument('--png', default=os.path.join(HERE, 'elem_scaling.png'))
    a = ap.parse_args()

    d = json.load(open(a.json))
    rows = sorted(d['rows'], key=lambda r: r['elem_max'])
    prov = d['provenance']

    x = np.array([r['elem_max'] for r in rows], float)
    y = np.array([r['sec_per_step'] for r in rows], float)
    ranks = np.array([r['ranks'] for r in rows])
    dxs = np.array([r['dx'] for r in rows])

    print(f"provenance: {prov['case']} @ {prov['sha']}"
          f"{' +dirty-src' if prov.get('dirty') else ''} on {prov['host']} "
          f"({prov['cpu']}, {prov['cores']} cores), {prov['fc']}, {prov['mpi']}, "
          f"netCDF {prov['netcdf']}, {prov['date']}")
    print()
    hdr = (f"{'dx':>5} {'ranks':>5} {'decomp':>8} {'nstep':>6} {'elem/rank':>11} "
           f"{'imbal':>6} {'loop s':>9} {'ms/step':>9} {'ms/step':>9} {'ms/step':>9} "
           f"{'ms/step':>9} {'ns/elem':>8} {'load':>5} {'cont':>5}")
    print(hdr)
    print(f"{'':>5} {'':>5} {'':>8} {'':>6} {'(max)':>11} {'':>6} {'(max)':>9} "
          f"{'LOOP':>9} {'kernel':>9} {'halo':>9} {'wall':>9} {'/step':>8} "
          f"{'pre':>5} {'end':>5}")
    print('-' * len(hdr))
    for r in rows:
        imb = r['elem_max'] / r['elem_mean'] - 1
        print(f"{r['dx']:5.0f} {r['ranks']:5d} {'x'.join(map(str, r['decomp'])):>8} "
              f"{r['nstep']:6d} {r['elem_max']:11,} {imb:+5.1%} "
              f"{r['t_loop_max']:9.2f} {r['sec_per_step']*1e3:9.3f} "
              f"{r['sec_per_step_kernel']*1e3:9.3f} {r['sec_per_step_mpi']*1e3:9.3f} "
              f"{r['sec_per_step_wall']*1e3:9.3f} "
              f"{r['sec_per_step']/r['elem_max']*1e9:8.1f} "
              f"{r['load_before']:5.1f} {r.get('contended', 0):5d}")
    print('LOOP = comp(9)-comp(1)-comp(8) per step (the measurand).  kernel = '
          'comp(3..6)/nstep (velDisp+KU+hourglass+faulting, pure per-element work).')
    print('halo = MPICommTime/nstep (sendrecv + mpi_barrier; includes 2 setup '
          'exchanges out of nstep+2).  cont = other eqdyna procs at run end.')

    b, a0, r2, res, se = ols(x, y)
    print()
    print('Global OLS  time/step = a + b * (elements per rank)')
    print(f'  slope  b = {b*1e9:.3f} ns per element per step  (+/- {se*1e9:.3f} 1-sigma)')
    print(f'  interc a = {a0*1e3:.3f} ms per step')
    print(f'  R^2      = {r2:.5f}')
    print(f'  intercept as a fraction of the smallest measured time/step: '
          f'{a0/y.min():.1%}')
    print('  residuals (measured - fit), ms/step, and % of measured:')
    for r, rr in zip(rows, res):
        print(f'    dx={r["dx"]:>4.0f} np={r["ranks"]:>3d}  {rr*1e3:+9.3f}  '
              f'{rr/r["sec_per_step"]:+7.1%}')
    # proportional (zero-intercept) fit, the hypothesis in its strict form
    b0 = float((x * y).sum() / (x * x).sum())
    res0 = y - b0 * x
    r2_0 = 1 - float((res0 ** 2).sum()) / float(((y - y.mean()) ** 2).sum())
    print()
    print(f'Strict hypothesis (forced through origin): b = {b0*1e9:.3f} ns/elem/step, '
          f'R^2 = {r2_0:.5f}')
    print('  worst proportional residual: '
          f'{max(abs(res0 / y)):.1%}')

    # An unweighted OLS over a 1000x dynamic range is dominated by the largest
    # points: it will report R^2 ~ 1 and a small intercept no matter how badly
    # the small-per-rank points behave. Two estimators that are sensitive to
    # relative error are therefore reported alongside.
    w = 1.0 / y ** 2                       # relative-error weighting
    sw, swx, swy = w.sum(), (w * x).sum(), (w * y).sum()
    swxx, swxy = (w * x * x).sum(), (w * x * y).sum()
    det = sw * swxx - swx ** 2
    bw = (sw * swxy - swx * swy) / det
    aw = (swy - bw * swx) / sw
    relres = (y - aw - bw * x) / y
    print()
    print('Relative-error (1/y^2 weighted) fit -- the estimator that can actually '
          'see the small-per-rank points:')
    print(f'  slope  b = {bw*1e9:.3f} ns per element per step')
    print(f'  interc a = {aw*1e3:.3f} ms per step '
          f'(= {aw/y.min():.1%} of the fastest configuration\'s time/step)')
    print(f'  rms relative residual = {np.sqrt((relres**2).mean()):.2%}, '
          f'max |relative residual| = {np.abs(relres).max():.2%}')

    lx, ly = np.log(x), np.log(y)
    pexp, lc = np.polyfit(lx, ly, 1)
    lfit = lc + pexp * lx
    lr2 = 1 - ((ly - lfit) ** 2).sum() / ((ly - ly.mean()) ** 2).sum()
    se_p = float(np.sqrt(((ly - lfit) ** 2).sum() / (len(x) - 2)
                         / ((lx - lx.mean()) ** 2).sum()))
    print()
    print('Power-law fit  time/step = C * (elements per rank)^p  '
          '(p = 1 exactly iff strictly proportional):')
    print(f'  p = {pexp:.4f} +/- {se_p:.4f} (1-sigma), C = {np.exp(lc):.4g}, '
          f'R^2(log) = {lr2:.5f}')
    print('  per-point deviation from the power law:')
    for r, dv in zip(rows, np.exp(ly - lfit) - 1):
        print(f'    dx={r["dx"]:>4.0f} np={r["ranks"]:>3d}  {dv:+7.1%}')

    print()
    print('Per-rank-count sub-fits (bandwidth test -- slope should be constant '
          'if only element count matters):')
    print(f"  {'ranks':>5} {'n':>3} {'slope ns/elem/step':>20} {'intercept ms':>13} {'R^2':>8}")
    slopes = {}
    for n in sorted(set(ranks)):
        m = ranks == n
        if m.sum() >= 2:
            if m.sum() == 2:
                bb = float(np.diff(y[m]) / np.diff(x[m]))
                aa = float(y[m][0] - bb * x[m][0])
                rr2 = float('nan')
            else:
                bb, aa, rr2, _, _ = ols(x[m], y[m])
            slopes[int(n)] = bb
            print(f'  {n:5d} {m.sum():3d} {bb*1e9:20.3f} {aa*1e3:13.3f} {rr2:8.5f}')
    if len(slopes) > 1:
        lo, hi = min(slopes.values()), max(slopes.values())
        print(f'  slope spread across rank counts: {hi/lo:.2f}x '
              f'({lo*1e9:.1f} -> {hi*1e9:.1f} ns/elem/step)')

    print()
    print('MATCHED SETS -- same elements/rank (within +/-6%) reached at different '
          'rank counts and resolutions. This is the hypothesis stated as a '
          'prediction: every member of a set must cost the same per step.')
    used = set()
    for i in range(len(rows)):
        if i in used:
            continue
        grp = [j for j in range(len(rows))
               if abs(x[j] / x[i] - 1) <= 0.06 and rows[j]['ranks'] != rows[i]['ranks']]
        if not grp:
            continue
        grp = sorted(set([i] + grp), key=lambda j: rows[j]['ranks'])
        used.update(grp)
        base = grp[0]
        print(f'  set @ ~{x[base]:,.0f} elem/rank:')
        for j in grp:
            print(f'    dx={rows[j]["dx"]:>4.0f} np={rows[j]["ranks"]:>3d} '
                  f'({"x".join(map(str, rows[j]["decomp"]))})  '
                  f'{x[j]:>9,.0f} elem/rank  {y[j]*1e3:8.3f} ms/step  '
                  f'{y[j]/x[j]*1e9:7.1f} ns/elem/step  '
                  f'-> x{y[j]/y[base]:.2f} vs the 1st row '
                  f'(elements x{x[j]/x[base]:.3f})')

    print()
    print('All near-size pairs across different rank counts:')
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if rows[i]['ranks'] == rows[j]['ranks']:
                continue
            ratio = x[j] / x[i]
            if 0.62 <= ratio <= 1.62:
                print(f'  dx={rows[i]["dx"]:>4.0f}/np={rows[i]["ranks"]:<3d} '
                      f'{x[i]:>9,.0f} el @ {y[i]*1e3:7.3f} ms   vs   '
                      f'dx={rows[j]["dx"]:>4.0f}/np={rows[j]["ranks"]:<3d} '
                      f'{x[j]:>9,.0f} el @ {y[j]*1e3:7.3f} ms   '
                      f'| elem x{ratio:.2f}, time x{y[j]/y[i]:.2f}, '
                      f'cost/elem x{(y[j]/x[j])/(y[i]/x[i]):.2f}')

    # ---------------- plot ----------------
    colors = {1: '#1b6ca8', 4: '#2e9e5b', 16: '#e08214', 48: '#b2182b'}
    marks = {500: 'o', 250: 's', 125: '^'}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.9))

    xs = np.logspace(np.log10(x.min() * 0.6), np.log10(x.max() * 1.6), 200)
    ax1.plot(xs, a0 + b * xs, '-', color='0.35', lw=1.3,
             label=f'OLS: {b*1e9:.2f} ns/elem + {a0*1e3:.2f} ms')
    ax1.plot(xs, b0 * xs, '--', color='0.55', lw=1.2,
             label=f'strict proportional: {b0*1e9:.2f} ns/elem')
    for r, xi, yi in zip(rows, x, y):
        ax1.plot(xi, yi, marks[r['dx']], color=colors[r['ranks']], ms=8,
                 mec='k', mew=0.6, zorder=3)
    ax1.set_xscale('log'); ax1.set_yscale('log')
    ax1.set_xlabel('elements per rank (max over ranks)')
    ax1.set_ylabel('wall seconds per timestep')
    ax1.set_title('Per-step cost vs elements per rank')
    ax1.grid(True, which='both', alpha=0.25, lw=0.5)
    ax1.legend(loc='upper left', fontsize=8, framealpha=0.9)

    # connect matched-size sets: if only elements/rank mattered, each connector
    # would be flat.
    seen = set()
    for i in range(len(rows)):
        if i in seen:
            continue
        grp = sorted({i} | {j for j in range(len(rows))
                            if abs(x[j] / x[i] - 1) <= 0.06
                            and rows[j]['ranks'] != rows[i]['ranks']})
        if len(grp) > 1:
            seen.update(grp)
            ax2.plot([x[j] for j in grp], [y[j] / x[j] * 1e9 for j in grp],
                     '-', color='0.55', lw=1.0, zorder=2)
    for r, xi, yi in zip(rows, x, y):
        ax2.plot(xi, yi / xi * 1e9, marks[r['dx']], color=colors[r['ranks']],
                 ms=8, mec='k', mew=0.6, zorder=3)
    ax2.axhline(b * 1e9, color='0.35', lw=1.3,
                label=f'OLS slope {b*1e9:.2f} ns/elem/step')
    ax2.set_xscale('log')
    ax2.set_xlabel('elements per rank (max over ranks)')
    ax2.set_ylabel('ns per element per timestep')
    ax2.set_title('Normalised cost: flat iff strictly proportional')
    ax2.grid(True, which='both', alpha=0.25, lw=0.5)
    ax2.margins(x=0.08, y=0.16)

    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls='', marker='o', ms=8, mec='k', mew=0.6,
                color=colors[n], label=f'{n} rank' + ('s' if n > 1 else ''))
         for n in sorted(colors) if n in set(ranks.tolist())]
    h += [Line2D([], [], ls='', marker=marks[v], ms=8, color='0.7', mec='k',
                 mew=0.6, label=f'dx = {v} m') for v in sorted(set(dxs.tolist()), reverse=True)]
    ax2.legend(handles=h + [Line2D([], [], color='0.35', lw=1.3,
                                   label=f'OLS slope {b*1e9:.1f}')],
               loc='lower left', fontsize=8, ncol=2, framealpha=0.9,
               bbox_to_anchor=(0.0, 0.0))

    fig.suptitle(f"EQdyna {prov['case']} @ {prov['sha']} -- {prov['host']}, "
                 f"{prov['cpu'].strip()}, {prov['fc'].split('(')[0].strip()}, "
                 f"{prov['mpi']}, {prov['date'][:10]}", fontsize=8.5, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(a.png, dpi=150)
    print(f'\nwrote {a.png}')


if __name__ == '__main__':
    main()
