#! /usr/bin/env python3
"""TPV12 rupture-time overlay: our FORTRAN run vs one independent SCEC cvws
submission, Michael Barall (FaultMod, 100 m, independent code/group).
python-jax is out of scope for this figure (the cplot/station reference is
a single-code cross-check, same role as make_tpv22_tpv23_overlay.py's barall
panel); see that script's own docstring for the general pattern this one
follows.

Barall's submission ships at 100 m / 1600 steps (8 s), the TPV12 spec's own
resolution. This repo's GATE-TIER run is coarser (see the run dir's own
par.dx/par.term) -- the figure title records which run directory produced
the overlaid contour, so a coarse-vs-spec-resolution figure is never
mistaken for a spec-resolution validation.

Regenerate:
    python3 scripts/figures/make_tpv12_overlay.py --run-dir test/test.tpv12
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NEVER = 999.0   # ours: 99999 s sentinel; barall: 1e9/1e10 s sentinel
DIP_DEG = 60.0  # test.tpv12's par.dip; down-dip distance = |z| / sin(dip)


def load_ours(frt_path):
    """frt.canonical.txt: col0 x(m, along-strike), col2 z(m, depth<=0),
    col3 rupture time(s). Converts to (along-strike km, down-dip km)."""
    a = np.loadtxt(frt_path)
    x = a[:, 0] / 1e3
    downdip = np.abs(a[:, 2]) / np.sin(np.deg2rad(DIP_DEG)) / 1e3
    t = np.where(a[:, 3] >= NEVER, np.nan, a[:, 3])
    return x, downdip, t


def load_cplot(path):
    """Barall's cplot: along-strike(m), down-dip(m), rupture time(s);
    sentinel 1e9/1e10. '#' comments plus one unparseable 'j k t' header
    line, same quirk as TPV22/23's/TPV29/30's cplot."""
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
    return b[:, 0] / 1e3, b[:, 1] / 1e3, np.where(b[:, 2] >= NEVER, np.nan, b[:, 2])


def to_grid(x, y, t):
    xs, ys = np.unique(x), np.unique(y)
    g = np.full((ys.size, xs.size), np.nan)
    g[np.searchsorted(ys, y), np.searchsorted(xs, x)] = t
    return xs, ys, g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--run-dir', required=True,
                    help='completed test.tpv12 case directory (any '
                         'resolution/term -- a resolution mismatch against '
                         'the 100 m archive is reported in the title, not '
                         'hidden)')
    _SHARED = os.path.expanduser(os.path.join(
        '~', 'shared_dataset', 'scec_cvws.tpv1213', 'raw'))
    ap.add_argument('--archive', default=_SHARED,
                    help='READ-ONLY cross-code archive (barall submission)')
    ap.add_argument('--out', default=os.path.join(
        ROOT, 'docs', 'figures', 'tpv12', 'tpv12_cplot_overlay.png'))
    ap.add_argument('--gate-resolution', action='store_true',
                    help='label the title as the everyday GATE-resolution '
                         'run rather than a spec-resolution run (set this '
                         'until the owner supplies a spec-resolution run '
                         'directory)')
    args = ap.parse_args()

    xo, zo, To = load_ours(os.path.join(args.run_dir, 'frt.canonical.txt'))
    xg, zg, Tg = to_grid(xo, zo, To)

    barall_path = os.path.join(args.archive, 'tpv12',
                               'barall-faultmod-100m-2009', 'cplot')
    xb, zb, Tb = load_cplot(barall_path)
    xbg, zbg, Tbg = to_grid(xb, zb, Tb)

    levels = np.arange(0, 8.1, 0.5)
    fig, ax = plt.subplots(figsize=(7.5, 5.0), dpi=200)
    c1 = ax.contour(xg, zg, Tg, levels=levels, colors='C3', linewidths=1.3)
    c2 = ax.contour(xbg, zbg, Tbg, levels=levels, colors='0.3',
                    linewidths=1.1, linestyles='--')
    ax.set_xlabel('along strike (km)')
    ax.set_ylabel('down-dip distance (km)')
    ax.invert_yaxis()

    n_ours = int(np.isfinite(To).sum())
    n_barall = int(np.isfinite(Tb).sum())
    res_tag = ('GATE-resolution run (coarser than the 100 m spec)'
               if args.gate_resolution else 'run directory as supplied')
    ax.set_title('TPV12 rupture time: EQdyna (solid) vs Barall/FaultMod '
                 '100 m (dashed)\n%s -- ruptured nodes: ours %d, barall %d'
                 % (res_tag, n_ours, n_barall))
    h1 = plt.Line2D([0], [0], color='C3', lw=1.3, label='EQdyna')
    h2 = plt.Line2D([0], [0], color='0.3', lw=1.1, ls='--',
                    label='Barall/FaultMod 100 m')
    ax.legend(handles=[h1, h2], loc='upper right', fontsize=8)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out)
    print('wrote %s (ours ruptured=%d/%d, barall ruptured=%d/%d)'
         % (args.out, n_ours, To.size, n_barall, Tb.size))


if __name__ == '__main__':
    main()
