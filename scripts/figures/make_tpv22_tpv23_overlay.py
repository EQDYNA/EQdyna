#! /usr/bin/env python3
"""TPV22/TPV23 rupture-time overlay: our FORTRAN run vs two independent SCEC
cvws submissions, Kaneko (SPECFEM3D, 100 m, independent method) and Barall
(FaultMod, 100 m, independent code/group). python-jax is explicitly out of
scope for this mission (dispatch brief) and is not plotted here.

Adapted from make_tpv29_tpv30_overlay.py (that script's own TPV29/30 geometry
is single-fault; TPV22/23 has two faults per case, separated here by the
fault-local y/stepover coordinate: y=0 is fault #1, y=fault2_z (nonzero) is
fault #2). 2x2 grid: rows = fault #1 / fault #2, columns = tpv22 / tpv23.

Barall resolution picked: FaultMod 100 m (`barall-faultmod-100m-2013`), to
match Kaneko's own 100 m -- a same-resolution, cross-method/cross-code pair.
Barall's 50 m and DayFD variants are archived (scec_archive/) but not
overlaid here, to keep the panel legible; they are available for a follow-up
figure if needed.

WHY THE CAPTION IS COMPUTED, NOT WRITTEN -- same discipline as
make_tpv29_tpv30_overlay.py: every number below is measured from the arrays
being plotted, at plot time, never hand-entered.

Regenerate:
    python3 scripts/figures/make_tpv22_tpv23_overlay.py \\
        --results '<repo>/test/test.%s'
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NEVER = 999.0   # ours: 99999 s sentinel; kaneko/barall: 1e9/1e10 s sentinel

FAULT2_Z = {'tpv22': -1600.0, 'tpv23': 1000.0}   # code-y stepover offset


def load_ours(path, fault2_z):
    """frt.canonical.txt: col0 x(m), col1 y(m, fault-selector), col2 z(m,
    depth<=0), col3 rupture time(s). Returns {1: (x,z,t) km/km/s, 2: (...)}."""
    a = np.loadtxt(path)
    out = {}
    for ift, ysel in ((1, 0.0), (2, fault2_z)):
        m = np.abs(a[:, 1] - ysel) < 0.5
        out[ift] = (a[m, 0] / 1e3, a[m, 2] / 1e3,
                    np.where(a[m, 3] >= NEVER, np.nan, a[m, 3]))
    return out


def load_cplot(path):
    """SCEC TPV22/23 cplot_N: along-strike(m), down-dip(m, >=0), rupture
    time(s); sentinel 1e9/1e10. '#' comments plus one unparseable 'j k t'
    header line, same quirk as TPV29/30's cplot."""
    rows = []
    for line in open(path, errors='replace'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        p = line.split()
        try:
            rows.append([float(p[0]), float(p[1]), float(p[2])])
        except (ValueError, IndexError):
            continue
    b = np.array(rows)
    return b[:, 0] / 1e3, -b[:, 1] / 1e3, np.where(b[:, 2] >= NEVER, np.nan, b[:, 2])


def to_grid(x, y, t):
    xs, ys = np.unique(x), np.unique(y)
    g = np.full((ys.size, xs.size), np.nan)
    g[np.searchsorted(ys, y), np.searchsorted(xs, x)] = t
    return xs, ys, g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results', default=os.path.join(ROOT, 'test', 'test.%s'),
                    help='printf-style template for the completed run dir, '
                         '%%s = tpv22|tpv23 (default: test/test.<case>)')
    ap.add_argument('--archive', default=os.path.join(ROOT, 'scec_archive'),
                    help='READ-ONLY cross-code archive (symlinked into this '
                         'worktree to the main checkout -- see dispatch brief)')
    ap.add_argument('--kaneko', default=os.path.join(
        ROOT, 'scratch', 'specs', 'tpv2223', 'scec_cross_code'),
        help='kaneko/payne reference dir (per-case subdir)')
    ap.add_argument('--out', default=os.path.join(
        ROOT, 'docs', 'figures', 'tpv22_23', 'tpv22_tpv23_cplot_overlay.png'))
    args = ap.parse_args()

    levels = np.arange(0, 15.1, 1.0)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), dpi=200,
                             sharex='col', sharey='row')
    facts = {}

    for col, case in enumerate(('tpv22', 'tpv23')):
        rdir = args.results % case if '%s' in args.results else os.path.join(
            args.results, 'test.' + case)
        ours = load_ours(os.path.join(rdir, 'frt.canonical.txt'),
                         FAULT2_Z[case])
        barall_dir = os.path.join(args.archive, case,
                                  'barall-faultmod-100m-2013')
        kaneko_dir = os.path.join(args.kaneko, case, 'kaneko')

        for row, ift in enumerate((1, 2)):
            ax = axes[row][col]
            xo, zo, To = to_grid(*ours[ift])
            ax.contour(xo, zo, To, levels=levels, colors='C3', linewidths=1.2,
                      zorder=3)
            f = dict(rupt=int(np.isfinite(To).sum()), tot=int(To.size))

            bp = os.path.join(barall_dir, 'cplot_%d' % ift)
            if os.path.isfile(bp):
                xb, zb, Tb = to_grid(*load_cplot(bp))
                ax.contour(xb, zb, Tb, levels=levels, colors='0.4',
                          linewidths=1.3, linestyles='--', zorder=2)
                f['barall_rupt'] = int(np.isfinite(Tb).sum())
                f['barall_tot'] = int(Tb.size)
            kp = os.path.join(kaneko_dir, 'cplot_%d.txt' % ift)
            if os.path.isfile(kp):
                xk, zk, Tk = to_grid(*load_cplot(kp))
                ax.contour(xk, zk, Tk, levels=levels, colors='C0',
                          linewidths=1.1, linestyles=':', zorder=1)
                f['kaneko_rupt'] = int(np.isfinite(Tk).sum())
                f['kaneko_tot'] = int(Tk.size)

            if ift == 1:
                ax.plot(-10, -10, 'k*', ms=10, zorder=4)  # hypocentre
            ax.set_title('%s, fault #%d  [ours %d/%d ruptured]'
                         % (case, ift, f['rupt'], f['tot']), fontsize=8.5)
            if col == 0:
                ax.set_ylabel('Depth (km)')
            if row == 1:
                ax.set_xlabel('Along strike (km)')
            facts[(case, ift)] = f

    axes[0][1].legend(
        handles=[Line2D([], [], color='C3', lw=1.2, label='EQdyna today (fortran)'),
                Line2D([], [], color='0.4', lw=1.3, ls='--',
                        label='Barall, FaultMod, 100 m (2013)'),
                Line2D([], [], color='C0', lw=1.1, ls=':',
                        label='Kaneko, SPECFEM3D, 100 m (2013)')],
        fontsize=7.5, loc='lower right', framealpha=0.9)

    lines = ['Every count below is measured from the arrays plotted above, at plot time.']
    for case in ('tpv22', 'tpv23'):
        for ift in (1, 2):
            f = facts[(case, ift)]
            lines.append(
                '%s fault#%d: ours %d/%d ruptured; barall %s/%s; kaneko %s/%s'
                % (case, ift, f['rupt'], f['tot'],
                   f.get('barall_rupt', '-'), f.get('barall_tot', '-'),
                   f.get('kaneko_rupt', '-'), f.get('kaneko_tot', '-')))
    lines.append('EQdyna mesh (200 m isotropic, tpv22 / 250 m isotropic, tpv23; '
                 'restored per-fault mesh extent, each fault on its own true '
                 'box) vs 100 m references -- a cross-code check at a finer, '
                 'not resolution-matched, gate. python-jax not plotted (out of scope).')
    fig.text(0.01, -0.02, '\n'.join(lines), fontsize=6.8, va='top', family='monospace')
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, bbox_inches='tight')
    print('wrote %s' % args.out)
    for l in lines:
        print(l)


if __name__ == '__main__':
    main()
