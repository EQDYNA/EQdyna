#! /usr/bin/env python3
"""Single-case, single-run-dir rupture-time (cplot) overlay for a TPV22/TPV23
RESOLUTION DIAGNOSTIC run (mira/row17-tpv2223-rebased, isotropic-resolution
investigation, 2026-10-01/02) -- e.g. test/test.tpv22.r400 (400 m isotropic)
or test/test.tpv23.r500 (500 m isotropic), run dirs that do not share a
common --results %-template with the OTHER case's own (different) resolution
suffix, so make_tpv22_tpv23_overlay.py's combined 2x2 figure (which processes
BOTH cases in one call) cannot be reused as-is for a single off-pattern dir.

Reuses, does not reimplement, make_tpv22_tpv23_overlay.py's own
load_ours/load_cplot/to_grid/FAULT2_Z (imported, not copied). 1x2 panel
(fault #1, fault #2) for ONE case/run-dir, same three curves (ours, Barall
FaultMod 100m, Kaneko SPECFEM3D 100m), same computed-at-plot-time caption
discipline as the original script.

Regenerate (example):
    python3 scripts/figures/make_tpv2223_resolution_cplot.py \\
        --case tpv22 --run-dir test/test.tpv22.r400 \\
        --out docs/figures/tpv22_23/tpv22_r400_cplot_overlay.png
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_tpv22_tpv23_overlay import load_ours, load_cplot, to_grid, FAULT2_Z  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--case', required=True, choices=('tpv22', 'tpv23'))
    ap.add_argument('--run-dir', required=True)
    ap.add_argument('--archive', default=os.path.join(ROOT, 'scec_archive'))
    ap.add_argument('--kaneko', default=os.path.join(
        ROOT, 'scratch', 'specs', 'tpv2223', 'scec_cross_code'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--label', default=None,
                    help='extra string for the figure title (e.g. "400 m isotropic")')
    args = ap.parse_args()

    case = args.case
    ours = load_ours(os.path.join(args.run_dir, 'frt.canonical.txt'), FAULT2_Z[case])
    barall_dir = os.path.join(args.archive, case, 'barall-faultmod-100m-2013')
    kaneko_dir = os.path.join(args.kaneko, case, 'kaneko')

    levels = np.arange(0, 15.1, 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=200, sharey=True)
    facts = {}
    for col, ift in enumerate((1, 2)):
        ax = axes[col]
        xo, zo, To = to_grid(*ours[ift])
        ax.contour(xo, zo, To, levels=levels, colors='C3', linewidths=1.3, zorder=3)
        f = dict(rupt=int(np.isfinite(To).sum()), tot=int(To.size))

        bp = os.path.join(barall_dir, 'cplot_%d' % ift)
        if os.path.isfile(bp):
            xb, zb, Tb = to_grid(*load_cplot(bp))
            ax.contour(xb, zb, Tb, levels=levels, colors='0.4', linewidths=1.3,
                      linestyles='--', zorder=2)
            f['barall_rupt'], f['barall_tot'] = int(np.isfinite(Tb).sum()), int(Tb.size)
        kp = os.path.join(kaneko_dir, 'cplot_%d.txt' % ift)
        if os.path.isfile(kp):
            xk, zk, Tk = to_grid(*load_cplot(kp))
            ax.contour(xk, zk, Tk, levels=levels, colors='C0', linewidths=1.1,
                      linestyles=':', zorder=1)
            f['kaneko_rupt'], f['kaneko_tot'] = int(np.isfinite(Tk).sum()), int(Tk.size)
        if ift == 1:
            ax.plot(-10, -10, 'k*', ms=10, zorder=4)
        ax.set_title('%s, fault #%d  [ours %d/%d ruptured]' % (case, ift, f['rupt'], f['tot']),
                     fontsize=9)
        ax.set_xlabel('Along strike (km)')
        if col == 0:
            ax.set_ylabel('Depth (km)')
        facts[ift] = f

    axes[1].legend(
        handles=[Line2D([], [], color='C3', lw=1.3, label='EQdyna today (fortran)'),
                Line2D([], [], color='0.4', lw=1.3, ls='--', label='Barall, FaultMod, 100 m'),
                Line2D([], [], color='C0', lw=1.1, ls=':', label='Kaneko, SPECFEM3D, 100 m')],
        fontsize=7, loc='lower right', framealpha=0.9)

    lab = (' (%s)' % args.label) if args.label else ''
    fig.suptitle('%s resolution diagnostic%s -- run dir %s' % (case, lab, args.run_dir),
                 fontsize=9, x=0.01, ha='left')
    lines = ['Every count below is measured from the arrays plotted above, at plot time.']
    for ift in (1, 2):
        f = facts[ift]
        lines.append('fault#%d: ours %d/%d ruptured; barall %s/%s; kaneko %s/%s'
                     % (ift, f['rupt'], f['tot'], f.get('barall_rupt', '-'), f.get('barall_tot', '-'),
                        f.get('kaneko_rupt', '-'), f.get('kaneko_tot', '-')))
    fig.text(0.01, -0.02, '\n'.join(lines), fontsize=7, va='top', family='monospace')
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, bbox_inches='tight')
    print('wrote %s' % args.out)
    for l in lines:
        print(l)


if __name__ == '__main__':
    main()
