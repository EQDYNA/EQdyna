#!/usr/bin/env python3
"""Row 149 force budget figures, from the CSVs written by budget_ana.py.
usage: budget_figs.py <csv_dir> <fig_dir>"""
import sys, os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = [('t12', 'TPV12, 60 deg dip (elastic)', 'd'), ('t12v90', 'TPV12 with dip = 90 (stand-in)', 'd'),
         ('t8', 'TPV8, vertical strike-slip (registered)', 's'), ('t13', 'TPV13, 60 deg dip (plastic)', 'd')]
DXS = [(200, '-'), (150, '--')]
LINES = [('BtS', 'tab:blue', 'stress divergence B^T sigma'), ('hg', 'tab:red', 'hourglass (displacement)'),
         ('damp', 'tab:gray', 'damping (stiffness + hourglass)'), ('grav', 'tab:green', 'gravity'),
         ('app', 'k', 'fault traction applied')]


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def series(rows, run, row, key):
    sel = [r for r in rows if r['run'] == run and r['row'] == str(row)]
    return [float(r['t_s']) for r in sel], [float(r[key]) for r in sel]


def part(rows, run, row, side, name, d):
    if name == 'damp':
        t, a = series(rows, run, row, 't%s_kdamp_%s' % (side, d))
        _, b = series(rows, run, row, 't%s_hgdamp_%s' % (side, d))
        _, c = series(rows, run, row, 't%s_mdamp_%s' % (side, d))
        return t, [x + y + z for x, y, z in zip(a, b, c)]
    if name == 'app':
        return series(rows, run, row, 't%sapp_%s' % (side, d))
    return series(rows, run, row, 't%s_%s_%s' % (side, name, d))


def fig_timeseries(rows, fig_dir, slipdir):
    fig, axs = plt.subplots(len(CASES), 2, figsize=(11, 12), sharex=True, sharey='row')
    for i, (c, title, sd) in enumerate(CASES):
        d = sd if slipdir else 'n'
        for j, (side, sname) in enumerate((('M', 'master side (y+; hanging wall when dipping)'),
                                           ('S', 'slave side (y-; footwall when dipping)'))):
            ax = axs[i][j]
            for dx, ls in DXS:
                run = '%s_dx%d' % (c, dx)
                for name, col, lab in LINES:
                    t, y = part(rows, run, 1, side, name, d)
                    if t:
                        ax.plot(t, y, ls, color=col, lw=1.1, label='%s, %d m' % (lab, dx))
            ax.axhline(0, color='0.7', lw=0.5)
            ax.set_title('%s\n%s' % (title, sname), fontsize=9)
            if j == 0:
                ax.set_ylabel('force / tributary area (MPa)')
            if i == len(CASES) - 1:
                ax.set_xlabel('time (s)')
    axs[0][1].legend(fontsize=7, ncol=2, loc='upper left')
    what = 'slip direction (dip; strike for TPV8)' if slipdir else 'fault-normal direction'
    fig.suptitle('Trace node (x = 0, free surface): force parts per side, %s; solid 200 m, dashed 150 m' % what, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(fig_dir, 'fig_timeseries_%s.png' % ('slipdir' if slipdir else 'normal'))
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print('wrote', out)


def fig_vs_dx(at, fig_dir):
    fig, axs = plt.subplots(1, len(CASES), figsize=(13, 4.6), sharey=True)
    for i, (c, title, sd) in enumerate(CASES):
        ax = axs[i]
        sel = [r for r in at if r['run'].split('_dx')[0] == c]
        for name, col, lab in LINES[:3]:
            for side, mk, sl in (('M', 'o', 'y+ side'), ('S', 's', 'y- side')):
                pts = sorted((float(r['dx_m']), float(r['t%s_%s_n' % (side, name)])) for r in sel
                             if r['row'] == '1' and name != 'damp')
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], '-' + mk, color=col, ms=5,
                            mfc='none' if side == 'S' else col, label='%s, %s' % (lab, sl))
        for row, ls in (('1', '-'), ('2', '--'), ('3', ':')):
            pts = sorted((float(r['dx_m']), float(r['Tcode_n'])) for r in sel if r['row'] == row)
            if pts:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], ls + 'k', marker='x', label='total n-traction, node row ' + row)
        ax.axhline(0, color='0.7', lw=0.5)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('dx (m)')
        ax.set_xlim(130, 220)
    axs[0].set_ylabel('MPa at t = 4 s')
    axs[0].legend(fontsize=6.5, loc='lower left')
    fig.suptitle('Fault-normal force parts at the trace node (per side) and total normal traction of the top three node rows, at 4 s vs dx', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(fig_dir, 'fig_parts_vs_dx.png')
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print('wrote', out)


def fig_surface_violation(es, fig_dir):
    """Top-row element stress at 4 s: fault-normal stress and the free-surface traction components it should not carry."""
    fig, axs = plt.subplots(1, len(CASES), figsize=(13, 4.2), sharey=True)
    for i, (c, title, sd) in enumerate(CASES):
        ax = axs[i]
        sel = [r for r in es if r['run'].split('_dx')[0] == c and r['row'] == '1' and float(r['t_s']) > 1.0
               and float(r['x_centroid_m']) > 0]
        for key, col in (('sig_nn', 'tab:blue'), ('sig_zn', 'tab:red'), ('sig_zz', 'tab:purple'), ('sig_zd', 'tab:orange')):
            for side, mk in (('1', 'o'), ('-1', 's')):
                pts = sorted((float(r['dx_m']), float(r[key])) for r in sel if r['side'] == side)
                if pts:
                    ax.plot([p[0] for p in pts], [p[1] for p in pts], '-' + mk, color=col, ms=5,
                            mfc='none' if side == '-1' else col, label='%s, %s side' % (key, 'y+' if side == '1' else 'y-'))
        ax.axhline(0, color='0.7', lw=0.5)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('dx (m)')
        ax.set_xlim(130, 220)
    axs[0].set_ylabel('one-point element stress (MPa) at 4 s')
    axs[0].legend(fontsize=6.5)
    fig.suptitle('Top-row elements beside the trace: sig_nn and the free-surface traction components (sig_zn, sig_zd, sig_zz; zero in the continuum at z = 0)', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(fig_dir, 'fig_surface_violation.png')
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print('wrote', out)


if __name__ == '__main__':
    csv_dir, fig_dir = sys.argv[1:3]
    os.makedirs(fig_dir, exist_ok=True)
    ts = read(os.path.join(csv_dir, 'timeseries.csv'))
    fig_timeseries(ts, fig_dir, False)
    fig_timeseries(ts, fig_dir, True)
    fig_vs_dx(read(os.path.join(csv_dir, 'at4s.csv')), fig_dir)
    fig_surface_violation(read(os.path.join(csv_dir, 'elem_stress.csv')), fig_dir)
