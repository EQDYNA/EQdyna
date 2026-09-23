#! /usr/bin/env python3
"""TPV29 vs TPV30 rupture-time overlay, ours against the 2015 SCEC submissions.

WHY THE CAPTION IS COMPUTED, NOT WRITTEN. The previous version of this file
carried a hand-written caption ("304 of 3321 TPV30 nodes (9.2%) are above
C0 + mus*(-sigma_n) at t=0", "TPV29 is the CONTROL ... its front is clean").
Those numbers described the 500 m gate run on master. The moment the geometry
fix landed and the resolution moved to 100 m they described nothing, while
still reading as a finding. Every number in the caption below is now MEASURED
from the very arrays being plotted, at plot time, by the functions in this
file. A caption cannot go stale if it cannot be written by hand.

WHAT IS PLOTTED. One panel per case: our rupture-time contours (red) over the
2015 EQdyna v3.1 SCEC submission at the same or finer resolution (grey
dashed), 1 s contours, hypocentre starred.

THE THREE MEASURED QUANTITIES, and why each is here:

  max tau/strength at t=0, and the count of nodes at or above strength.
      The initial on-fault traction resolved on the mesh-facet normal in the
      assigned un/us/ud frame -- testsys/regression/
      test_rough_fault_normal_consistency.py:_tauOverStrength, imported here
      rather than re-implemented so the figure and the gate cannot drift. A
      node at or above 1.0 fails at t=0 without a front reaching it.

  SEEDS: a ruptured node strictly earlier than every ruptured 8-neighbour.
      This is the coherence metric the older comparison lacked. Rupture-extent
      overlap cannot see spontaneous failure -- a spuriously ruptured node
      still counts as ruptured -- and neither can a median slip difference.
      8-connectivity is not a free choice: it reproduces the 500 m numbers
      this figure is read against (28 for TPV30, 2 for TPV29) where
      4-connectivity gives 31.

  median |dt| vs the 2015 submission, nearest node, ours -> archive.
      Direction pinned the same way, by reproducing the 500 m table (TPV29
      0.1443 s; the reverse direction gives 0.1577 s).

WHAT THIS FIGURE DOES NOT ESTABLISH. It is Fortran-vs-Fortran across roughly a
decade of the same code family. It certifies neither run against the
benchmark; only an independent cross-code comparison does that.

Regenerate:
    python3 scripts/figures/make_tpv29_tpv30_overlay.py                  # dx=100 runs
    python3 scripts/figures/make_tpv29_tpv30_overlay.py --results REFS   # 500 m refs
"""
import argparse
import importlib.util
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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


def seed_count(T):
    """Ruptured nodes strictly earlier than every ruptured 8-neighbour."""
    R = np.isfinite(T)
    nz, nx = T.shape
    n = 0
    for j in range(nz):
        for i in range(nx):
            if not R[j, i]:
                continue
            nb = [T[j+dj, i+di]
                  for dj, di in ((1,0),(-1,0),(0,1),(0,-1),
                                 (1,1),(1,-1),(-1,1),(-1,-1))
                  if 0 <= j+dj < nz and 0 <= i+di < nx and R[j+dj, i+di]]
            if nb and T[j, i] < min(nb):
                n += 1
    return n


def median_dt(ours, arc):
    """For each of OUR ruptured nodes, |dt| to the nearest ruptured archive
    node; median over our nodes. `ours`/`arc` are (x, y, t) triples in km/s."""
    xo, yo, To = ours
    xa, ya, Ta = arc
    m = np.isfinite(To); xo, yo, To = xo[m], yo[m], To[m]
    m = np.isfinite(Ta); xa, ya, Ta = xa[m], ya[m], Ta[m]
    P = np.stack([xa, ya], 1)
    d = np.empty(xo.size)
    for k in range(xo.size):
        i = np.argmin((P[:, 0]-xo[k])**2 + (P[:, 1]-yo[k])**2)
        d[k] = abs(Ta[i] - To[k])
    return float(np.median(d))


def tau_over_strength(case_dir, dx):
    """(max, count>=1, n) of t=0 tau/strength, via the REGRESSION module's own
    function so the figure and the gate can never disagree."""
    rp = os.path.join(ROOT, 'testsys', 'regression',
                      'test_rough_fault_normal_consistency.py')
    spec = importlib.util.spec_from_file_location('rfn_for_figure', rp)
    rfn = importlib.util.module_from_spec(spec)
    src = open(rp).read().replace("if __name__ == '__main__':\n    main()", "")
    exec(compile(src, rp, 'exec'), rfn.__dict__)
    g = rfn.loadGeoTools(case_dir)
    with rfn.inDir(case_dir):
        y, dydx, dydz = g.faultGridForCase(dx)
    fX, fZ = g.meshConsistentDerivatives(y, dx)
    r = rfn._tauOverStrength(y, dydx, dydz, fX, fZ, dx)[1:-1, 1:-1]
    return float(r.max()), int((r >= 1.0).sum()), int(r.size)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results', default=os.path.join(ROOT, 'runs', 'dx100.%s'),
                    help='printf-style template for the result dir, %%s = tpv29|tpv30 '
                         '(default: the dx=100 runs)')
    ap.add_argument('--dx', type=float, default=100.0, help='our resolution, m')
    ap.add_argument('--archive-res', default='100m',
                    help='which 2015 submission to overlay (100m|50m|25m)')
    ap.add_argument('--archive', default=os.path.join(ROOT, 'scec_archive'),
                    help='the 2015 SCEC submissions, READ-ONLY. Not tracked in '
                         'git, so a worktree does not have its own copy; point '
                         'this at the main checkout. Never symlinked into a '
                         'run directory -- a mis-set flag in a run that writes '
                         'through a symlink would corrupt the oracle.')
    ap.add_argument('--out', default=os.path.join(
        ROOT, 'docs', 'figures', 'tpv30', 'tpv29_vs_tpv30_cplot_overlay.png'))
    args = ap.parse_args()

    levels = np.arange(0, 13.1, 1.0)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.6), dpi=200, sharex=True)
    facts = {}

    for ax, case in zip(axes, ('tpv29', 'tpv30')):
        rdir = args.results % case if '%s' in args.results else os.path.join(
            args.results, 'test.' + case)
        ours = load_ours(os.path.join(rdir, 'frt.canonical.txt'))
        xo, yo, To = to_grid(*ours)
        if not os.path.isdir(args.archive):
            raise FileNotFoundError(
                'the 2015 SCEC submissions are not at %s. They are untracked, '
                'so a git worktree has none; pass --archive pointing at the '
                'main checkout. Refusing to plot ours alone -- a one-curve '
                '"overlay" is not a comparison.' % args.archive)
        arcdir = os.path.join(args.archive, case,
                              'eqdyna-v3.1-%s-2015' % args.archive_res)
        arc = load_archive(os.path.join(arcdir, 'cplot'))
        xp, yp, Tp = to_grid(*arc)
        ax.contour(xp, yp, Tp, levels=levels, colors='0.5', linewidths=1.4, linestyles='--')
        ax.contour(xo, yo, To, levels=levels, colors='C3', linewidths=1.1)
        ax.plot(-5, -10, 'k*', ms=9)
        ax.set_ylabel('Depth (km)')

        mx, nat, ntot = tau_over_strength(
            os.path.join(ROOT, 'case_input', 'test.' + case), args.dx)
        f = dict(rupt=int(np.isfinite(To).sum()), tot=int(To.size),
                 seeds=seed_count(To), taumax=mx, atstr=nat, ntau=ntot,
                 mdt=median_dt(ours, arc))
        for res in ('100m', '50m'):
            p = os.path.join(args.archive, case,
                             'eqdyna-v3.1-%s-2015' % res, 'cplot')
            f['mdt_' + res] = median_dt(ours, load_archive(p)) if os.path.isfile(p) else None
        facts[case] = f

        ax.set_title('%s -- ours %g m (red) vs 2015 %s (grey dashed), 1 s contours '
                     '[%d/%d ruptured, %d seeds]'
                     % (case, args.dx, args.archive_res, f['rupt'], f['tot'], f['seeds']),
                     fontsize=8.5)

    axes[-1].set_xlabel('Along strike (km)')
    axes[0].legend(handles=[Line2D([], [], color='C3', lw=1.1,
                                   label='EQdyna today, %g m' % args.dx),
                            Line2D([], [], color='0.5', lw=1.4, ls='--',
                                   label='EQdyna v3.1, %s (2015 SCEC submission)'
                                         % args.archive_res)],
                   fontsize=7, loc='lower right', framealpha=0.9)

    def line(c):
        f = facts[c]
        s = ('%s  t=0 max tau/strength %.5f, %d of %d nodes at or above strength;  '
             '%d of %d ruptured;  %d seeds;  median |dt| vs 2015 %s'
             % (c.upper(), f['taumax'], f['atstr'], f['ntau'], f['rupt'], f['tot'],
                f['seeds'], ', '.join('%s %.4f s' % (r, f['mdt_' + r])
                                      for r in ('100m', '50m')
                                      if f.get('mdt_' + r) is not None)))
        return s

    fig.text(0.01, -0.015,
             'Every number below is measured from the arrays plotted above, at plot time.\n'
             + line('tpv29') + '\n' + line('tpv30') + '\n'
             'SEED = a ruptured node strictly earlier than all 8 of its ruptured neighbours: the count of\n'
             'fronts that started on their own rather than arriving. Rupture-extent overlap and median slip\n'
             'difference cannot see this -- a spuriously ruptured node still counts as ruptured.\n'
             'Fortran-vs-Fortran across ~a decade: certifies neither run against the benchmark.',
             fontsize=6.5, va='top', family='monospace')
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, bbox_inches='tight')
    print('wrote %s' % args.out)
    for c in ('tpv29', 'tpv30'):
        print(line(c))


if __name__ == '__main__':
    main()
