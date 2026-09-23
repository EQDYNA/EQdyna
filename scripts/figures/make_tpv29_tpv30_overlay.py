#! /usr/bin/env python3
"""TPV29 vs TPV30 rupture-time overlay, ours against the 2015 SCEC submissions.

WHY BOTH CASES ARE ON ONE FIGURE. Alone, TPV30's panel invites the reading
"500 m is too coarse" -- the spec asks for 50 m and accepts 100 m, and our
gate reference is 500 m, so resolution is the obvious suspect. TPV29 is the
control that kills it: SAME dx, same mesh, same rough geometry, same 2015
comparison, and its front is clean. Whatever breaks TPV30 is not the mesh.

WHAT THE FIGURE SHOWS. TPV30's rupture-time field is not monotonic outward --
it carries dozens of isolated closed contours disconnected from the front.
That is the signature of nodes rupturing spontaneously rather than being
reached by a propagating front. 304 of 3321 TPV30 fault nodes (9.2%) sit
above C0 + mus*(-sigma_n) at t=0, and each seeds its own island.

WHY THE RECORDED NUMBERS MISSED IT. The prior comparison reported 99.2%
rupture-extent overlap and a 6.3% median slip difference, both reassuring.
Extent cannot see this: a spuriously-ruptured node still counts as ruptured.
No recorded metric measured the front's COHERENCE, which is the thing that
is wrong.

This is a bug report, not a validation figure. It is Fortran-vs-Fortran
across ~a decade, so it certifies neither run against the benchmark.

Regenerate:  python3 scripts/figures/make_tpv29_tpv30_overlay.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, 'docs', 'figures', 'tpv30', 'tpv29_vs_tpv30_cplot_overlay.png')
NEVER = 999.0   # both writers use a large sentinel; ours 99999, the 2015 files 1000


def load_ours(path):
    """frt.canonical.txt: col0 x(m), col2 z(m), col3 rupture time (s)."""
    a = np.loadtxt(path)
    return a[:, 0] / 1e3, a[:, 2] / 1e3, np.where(a[:, 3] >= NEVER, np.nan, a[:, 3])


def load_archive(path):
    """SCEC cplot: along-strike(m), down-dip(m), rupture time(s), '#' comments
    plus one unparseable 'j k t' header line."""
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        p = line.split()
        try:
            rows.append([float(p[0]), float(p[1]), float(p[2])])
        except (ValueError, IndexError):
            continue
    b = np.array(rows)
    return b[:, 0] / 1e3, -np.abs(b[:, 1]) / 1e3, np.where(b[:, 2] >= NEVER, np.nan, b[:, 2])


def to_grid(x, y, t):
    xs, ys = np.unique(x), np.unique(y)
    g = np.full((ys.size, xs.size), np.nan)
    g[np.searchsorted(ys, y), np.searchsorted(xs, x)] = t
    return xs, ys, g


def main():
    cases = [('test.tpv29', 'tpv29'), ('test.tpv30', 'tpv30')]
    levels = np.arange(0, 13.1, 1.0)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.6), dpi=200, sharex=True)

    for ax, (case, arc) in zip(axes, cases):
        xo, yo, To = to_grid(*load_ours(
            os.path.join(ROOT, 'test.reference.results', case, 'frt.canonical.txt')))
        xp, yp, Tp = to_grid(*load_archive(
            os.path.join(ROOT, 'scec_archive', arc, 'eqdyna-v3.1-100m-2015', 'cplot')))
        ax.contour(xp, yp, Tp, levels=levels, colors='0.5', linewidths=1.4, linestyles='--')
        ax.contour(xo, yo, To, levels=levels, colors='C3', linewidths=1.1)
        ax.plot(-5, -10, 'k*', ms=9)
        ax.set_ylabel('Depth (km)')
        n = int(np.isfinite(To).sum())
        ax.set_title('%s  --  ours 500 m (red) vs 2015 100 m (grey dashed), 1 s contours '
                     '[%d/%d nodes ruptured]' % (case, n, To.size), fontsize=8.5)

    axes[-1].set_xlabel('Along strike (km)')
    axes[0].legend(handles=[Line2D([], [], color='C3', lw=1.1, label='EQdyna today, 500 m'),
                            Line2D([], [], color='0.5', lw=1.4, ls='--',
                                   label='EQdyna v3.1, 100 m (2015 SCEC submission)')],
                   fontsize=7, loc='lower right', framealpha=0.9)
    fig.text(0.01, -0.015,
             'TPV29 is the CONTROL: same dx, same mesh, same rough geometry -- its front is clean, so TPV30\'s\n'
             'scattered closed contours are NOT a resolution artifact. They are nodes rupturing spontaneously:\n'
             '304 of 3321 TPV30 nodes (9.2%) are above C0 + mus*(-sigma_n) at t=0. Extent overlap (99.2%) and\n'
             'median slip difference (6.3%) cannot see this -- a spuriously-ruptured node still counts as ruptured.\n'
             'Fortran-vs-Fortran across ~a decade: certifies neither run against the benchmark.',
             fontsize=6.5, va='top', family='monospace')
    fig.tight_layout()
    fig.savefig(OUT, bbox_inches='tight')
    print('wrote %s' % OUT)


if __name__ == '__main__':
    main()
