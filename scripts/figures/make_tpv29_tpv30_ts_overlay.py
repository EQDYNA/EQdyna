#! /usr/bin/env python3
"""TPV29/TPV30 station time-series overlays, ours against the 2015 SCEC submissions.

Companion to make_tpv29_tpv30_overlay.py (which does the rupture-time cplot).
This does the two time-series families the SCEC benchmark defines:

  on-fault  (faultst<strike>dp<depth>)  slip, slip rate, shear stress
  off-fault (body<normal>st<strike>dp<depth>)  displacement, velocity

WHAT IS AND IS NOT ESTABLISHED. Fortran-vs-Fortran across roughly a decade of
the same code family: a cross-version CONSISTENCY check. It is not an
independent validation and no panel here is an accuracy claim.

NO INTERPOLATION ANYWHERE, in either axis:

  SPACE. At dx = 100 m our mesh nodes coincide with the archive's. Every SCEC
  station name encodes its position in units of 100 m (faultst-042dp061 is
  x = -4200 m, depth = 6100 m), so every station is exactly on a node of both
  meshes. `exact_stations` re-derives this from the NAME and any station that
  is not a multiple of the mesh spacing is EXCLUDED and listed, never snapped
  to a neighbour. At the 500 m gate this rule excluded 11 of 24; at 100 m it
  should exclude none, and the figure states the count it actually got.

  TIME. Ours steps at dt = 0.5*dx/vp = 8.3333 ms, the 2015 run at 8.0 ms, so
  the sample times genuinely differ. Each series is therefore drawn against
  ITS OWN time axis. Resampling one onto the other would invent values and
  would be indistinguishable, by eye, from a physical phase difference --
  which is exactly the thing these panels exist to show.

THE FORTRAN EXPONENT DEFECT is handled on READ, not by rewriting the run's
output. `library_output.f90` writes e15.7 and loses the 'E' when the exponent
needs four characters ("0.1341312-114"). scripts/correctSCECStOutputFormat.py
fixes that destructively, in place; this reader repairs the same pattern in
memory so the run directory stays exactly as the solver left it.

Regenerate:
    python3 scripts/figures/make_tpv29_tpv30_ts_overlay.py \
        --results '<root>/runs/dx100.%s' --dx 100 \
        --archive <main-checkout>/scec_archive
"""
import argparse
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_EXP = re.compile(r'(\d\.\d+)([-+]\d\d\d+)')     # 0.1341312-114 -> 0.1341312E-114

# --- station selection -----------------------------------------------------
# Chosen to SPAN the fault, not to flatter it. Fault: both along-strike
# directions from the x = -5 km hypocentre (backward -15, forward +15), the
# free surface and depth (0 -> 13 km), and the hypocentre itself.
FAULT_ST = ['faultst-150dp120',   # far backward along strike, deep
            'faultst-050dp000',   # hypocentre strike, free surface
            'faultst-050dp100',   # THE HYPOCENTRE (x=-5 km, 10 km deep)
            'faultst050dp000',    # forward, free surface
            'faultst100dp050',    # forward, mid-depth
            'faultst150dp130']    # far forward along strike, deep
# Body: near/far NORMAL to the fault (3 vs 20 km) and both sides (+/-3 km),
# plus along-strike contrast (0 vs +15 km) to catch forward directivity.
BODY_ST = ['body-030st000dp000',  # near side, 3 km off, mid-strike
           'body030st000dp000',   # FAR side, 3 km off -- the across-fault pair
           'body-030st150dp000',  # near side, 3 km off, forward directivity
           'body-200st000dp000']  # far field, 20 km off

FAULT_COLS = [(1, 'h-slip (m)'), (2, 'h-slip rate (m/s)'), (3, 'h-shear stress (MPa)')]
BODY_COLS = [(1, 'h-disp (m)'), (5, 'n-disp (m)'), (2, 'h-vel (m/s)'), (6, 'n-vel (m/s)')]


def read_ts(path):
    """Numeric block of a SCEC-format station file, repairing the lost-'E'
    exponent in memory. Raises if the file has no numeric rows -- an empty
    panel must not be silently drawn as a flat line."""
    rows = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        if line.lstrip()[0].isalpha():        # the ' t h-slip ...' name line
            continue
        try:
            rows.append([float(v) for v in _EXP.sub(r'\1E\2', line).split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError('%s: no numeric rows' % path)
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def station_coords(name):
    """(x, y, depth) in metres, re-derived from the SCEC station NAME.
    faultst-042dp061 -> (-4200, None, 6100); body-030st150dp000 -> (15000, -3000, 0)."""
    m = re.match(r'^faultst(-?\d+)dp(-?\d+)$', name)
    if m:
        return int(m.group(1))*100, None, int(m.group(2))*100
    m = re.match(r'^body(-?\d+)st(-?\d+)dp(-?\d+)$', name)
    if m:
        return int(m.group(2))*100, int(m.group(1))*100, int(m.group(3))*100
    raise ValueError('unrecognised SCEC station name %r' % name)


def exact_stations(names, dx):
    """Split into (exact, excluded): a station is EXACT when every coordinate
    it names is an integer multiple of the mesh spacing, so both meshes have a
    node there. Excluded stations are returned, never snapped."""
    exact, excluded = [], []
    for n in names:
        c = [v for v in station_coords(n) if v is not None]
        (exact if all(abs(v) % dx == 0 for v in c) else excluded).append(n)
    return exact, excluded


def panel(ax, ours, arc, col, label, first):
    if arc is not None:
        ax.plot(arc[:, 0], arc[:, col], color='0.5', lw=1.4, ls='--', zorder=1)
    ax.plot(ours[:, 0], ours[:, col], color='C3', lw=1.0, zorder=2)
    ax.set_ylabel(label, fontsize=6.5)
    ax.tick_params(labelsize=6)
    ax.grid(alpha=0.25, lw=0.4)


def make_family(case, rdir, adir, names, cols, dx, out, kind):
    ex, excl = exact_stations(names, dx)
    nrow, ncol = len(ex), len(cols)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6*ncol, 1.35*nrow), dpi=190,
                             squeeze=False, sharex=True)
    stats = []
    for i, name in enumerate(ex):
        ours = read_ts(os.path.join(rdir, name + '.txt'))
        ap = os.path.join(adir, name)
        arc = read_ts(ap) if os.path.isfile(ap) else None
        for j, (col, label) in enumerate(cols):
            ax = axes[i][j]
            panel(ax, ours, arc, col, label, i == 0)
            if j == 0:
                x, y, d = station_coords(name)
                loc = ('x %+.1f km, %.1f km deep' % (x/1e3, d/1e3) if y is None
                       else 'x %+.1f km, y %+.1f km' % (x/1e3, y/1e3))
                ax.text(0.02, 0.92, '%s\n%s' % (name, loc), transform=ax.transAxes,
                        fontsize=5.6, va='top', family='monospace')
        if arc is not None:
            # peak-amplitude ratio on the primary column, ours vs 2015
            po, pa = np.abs(ours[:, cols[0][0]]).max(), np.abs(arc[:, cols[0][0]]).max()
            stats.append((name, po, pa, po/pa if pa else np.nan))
    for j in range(ncol):
        axes[-1][j].set_xlabel('Time (s)', fontsize=7)

    # Title left-aligned and legend hard right on the SAME band, so a long
    # case name can never run under the legend box (it did when both were
    # centred/upper-right).
    fig.suptitle('%s -- %s stations, dx = %g m' % (case, kind, dx),
                 fontsize=9, x=0.01, ha='left')
    fig.legend(handles=[Line2D([], [], color='C3', lw=1.0, label='EQdyna today, %g m' % dx),
                        Line2D([], [], color='0.5', lw=1.4, ls='--',
                               label='EQdyna v3.1, 100 m (2015 SCEC submission)')],
               fontsize=6.5, loc='upper right', bbox_to_anchor=(0.995, 1.0),
               framealpha=0.9)

    lines = ['Every number here is measured from the arrays plotted above, at plot time.',
             '%d of %d %s stations are EXACT nodes of both meshes at dx=%g m (all coordinates '
             'multiples of the spacing); %d excluded, never snapped: %s'
             % (len(ex), len(names), kind, dx, len(excl), ', '.join(excl) if excl else 'none'),
             'peak |%s| ours / 2015: ' % cols[0][1]
             + ', '.join('%s %.3f' % (n.replace('faultst', 'ft').replace('body', 'b'), r)
                         for n, _, _, r in stats),
             'Sample times differ (ours dt=%.4f s, 2015 dt=0.0080 s); each series is drawn on '
             'its own axis and neither is resampled.' % (0.5*dx/6000.),
             'Fortran-vs-Fortran across ~a decade: a cross-version consistency check, not an '
             'independent validation and not an accuracy claim.']
    fig.text(0.01, 0.005, '\n'.join(lines), fontsize=6.0, va='top', family='monospace')
    fig.tight_layout(rect=[0, 0.055, 1, 0.965])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print('wrote %s  (%d stations, %d excluded)' % (out, len(ex), len(excl)))
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results', default=os.path.join(ROOT, 'runs', 'dx100.%s'))
    ap.add_argument('--dx', type=float, default=100.0)
    ap.add_argument('--archive', default=os.path.join(ROOT, 'scec_archive'),
                    help='2015 submissions, READ-ONLY; untracked, so pass the '
                         'main checkout when running from a worktree')
    ap.add_argument('--outdir', default=os.path.join(ROOT, 'docs', 'figures', 'tpv30'))
    args = ap.parse_args()
    if not os.path.isdir(args.archive):
        raise FileNotFoundError(
            '2015 submissions not at %s -- pass --archive. Refusing to draw '
            'ours alone; a one-curve "overlay" is not a comparison.' % args.archive)

    for case in ('tpv29', 'tpv30'):
        rdir = args.results % case if '%s' in args.results else args.results
        adir = os.path.join(args.archive, case, 'eqdyna-v3.1-100m-2015')
        make_family(case, rdir, adir, FAULT_ST, FAULT_COLS, args.dx,
                    os.path.join(args.outdir, '%s_ts_fault_overlay.png' % case), 'on-fault')
        make_family(case, rdir, adir, BODY_ST, BODY_COLS, args.dx,
                    os.path.join(args.outdir, '%s_ts_body_overlay.png' % case), 'off-fault body')


if __name__ == '__main__':
    main()
