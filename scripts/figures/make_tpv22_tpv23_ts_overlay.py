#! /usr/bin/env python3
"""TPV22/TPV23 on-fault (+ one off-fault) station time-series overlay: our
FORTRAN run vs Kaneko (SPECFEM3D, independent) and Barall (FaultMod,
independent) SCEC cvws submissions. python-jax out of scope (dispatch brief).

Adapted from make_tpv29_tpv30_ts_overlay.py. Differences forced by TPV22/23's
geometry: two faults, so the 3 stations this mission's own evidence script
(testsys/parity/evidence_tpv22_23_scec_comparison.py) already compares are
used here too, rather than re-deriving a station set from mesh spacing --
fault1st000dp000 (fault #1), fault2st000dp000 / fault2st050dp050 (fault #2,
named WITHOUT the repo's own 'ft2_' file-name prefix in the archives).

ALL 7 non-time columns are plotted (h-slip/-rate/-stress, v-slip/-rate/-stress,
n-stress) -- the SCEC on-fault format is identical across all three sources
(8 columns, time series in e15.7/e15.7), so there is no reason to drop the
v-/n- columns the way the TPV29/30 script had to size its own station list
around mesh-exactness.

OFF-FAULT: Kaneko's archive here (~/shared_dataset/scec_cvws.tpv2223/) was
fetched for exactly the 3 on-fault stations this mission needed (see that
module's docstring) and carries no body* files, so no Kaneko off-fault
overlay is possible. Barall DID serve body* files, and one of them --
body030st050dp000 (x=+5 km, y=+3 km i.e. far side, surface) -- matches one of
this case's own off-fault stations exactly (case_input/test.tpv22/
tpv22_23_common.py's st_coor_off_fault includes [5.0, 3.0, 0.0]). That pair is
plotted; no fabricated comparison is made where neither source has a match.

Regenerate:
    python3 scripts/figures/make_tpv22_tpv23_ts_overlay.py
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
_EXP = re.compile(r'(\d\.\d+)([-+]\d\d\d+)')     # lost-'E' exponent repair, belt-and-braces

# (archive file name, our file name, human label)
ON_FAULT = [
    ('fault1st000dp000', 'faultst000dp000.txt', 'fault #1, 0 km strike, 0 km down-dip'),
    ('fault2st000dp000', 'faultstft2_000dp000.txt', 'fault #2, 0 km strike, 0 km down-dip'),
    ('fault2st050dp050', 'faultstft2_050dp050.txt', 'fault #2, 5 km strike, 5 km down-dip'),
]
ON_FAULT_COLS = [(1, 'h-slip (m)'), (2, 'h-slip-rate (m/s)'), (3, 'h-shear-str (MPa)'),
                (4, 'v-slip (m)'), (5, 'v-slip-rate (m/s)'), (6, 'v-shear-str (MPa)'),
                (7, 'n-stress (MPa)')]

OFF_FAULT = [('body030st050dp000', 'body030st050dp000.txt',
             'off-fault, x=+5 km, y=+3 km (far side, surface)')]
OFF_FAULT_COLS = [(1, 'h-disp (m)'), (2, 'h-vel (m/s)'),
                  (5, 'n-disp (m)'), (6, 'n-vel (m/s)')]


def read_ts(path):
    rows = []
    for line in open(path, errors='replace'):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        if s[0].isalpha():
            continue
        try:
            rows.append([float(v) for v in _EXP.sub(r'\1E\2', s).split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError('%s: no numeric rows' % path)
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def panel(ax, series, cols_idx):
    for (label, arr, color, ls, lw, z) in series:
        if arr is not None and arr.shape[1] > cols_idx:
            ax.plot(arr[:, 0], arr[:, cols_idx], color=color, lw=lw, ls=ls, zorder=z)


def make_family(case, rdir, kaneko_dir, barall_dir, stations, cols, out, kind):
    nrow, ncol = len(stations), len(cols)
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.5*ncol, 1.6*nrow), dpi=190,
                             squeeze=False, sharex=True)
    stats = []
    for i, (arc_name, our_name, label) in enumerate(stations):
        ours = read_ts(os.path.join(rdir, our_name))
        kp = os.path.join(kaneko_dir, arc_name + ('.txt' if kind == 'on-fault' else ''))
        bp = os.path.join(barall_dir, arc_name)
        kan = read_ts(kp) if os.path.isfile(kp) else None
        bar = read_ts(bp) if os.path.isfile(bp) else None
        for j, (col, clabel) in enumerate(cols):
            ax = axes[i][j]
            series = [('ours', ours, 'C3', '-', 1.1, 3),
                     ('barall', bar, '0.4', '--', 1.3, 2),
                     ('kaneko', kan, 'C0', ':', 1.1, 1)]
            panel(ax, series, col)
            ax.set_ylabel(clabel, fontsize=6.3)
            ax.tick_params(labelsize=6)
            ax.grid(alpha=0.25, lw=0.4)
            if j == 0:
                ax.text(0.02, 0.92, label, transform=ax.transAxes,
                       fontsize=5.8, va='top', family='monospace')
        prim = cols[0][0]
        po = float(np.abs(ours[:, prim]).max())
        pk = float(np.abs(kan[:, prim]).max()) if kan is not None else None
        pb = float(np.abs(bar[:, prim]).max()) if bar is not None else None
        stats.append((arc_name, po, pk, pb))
    for j in range(ncol):
        axes[-1][j].set_xlabel('Time (s)', fontsize=7)

    fig.suptitle('%s -- %s stations (fortran vs independent SCEC cvws submissions)'
                % (case, kind), fontsize=9, x=0.01, ha='left')
    fig.legend(handles=[Line2D([], [], color='C3', lw=1.1, label='EQdyna today (fortran)'),
                        Line2D([], [], color='0.4', lw=1.3, ls='--', label='Barall, FaultMod, 100 m'),
                        Line2D([], [], color='C0', lw=1.1, ls=':', label='Kaneko, SPECFEM3D, 100 m')],
               fontsize=6.5, loc='upper right', bbox_to_anchor=(0.995, 1.0), framealpha=0.9)

    def fmt(v):
        return ('%.3f' % v) if v is not None else 'n/a'
    lines = ['Every number here is measured from the arrays plotted above, at plot time.',
             'peak |%s|  ours / kaneko / barall:  ' % cols[0][1]
             + '; '.join('%s %s/%s/%s' % (n, fmt(po), fmt(pk), fmt(pb))
                        for n, po, pk, pb in stats),
             'EQdyna mesh (200 m isotropic, tpv22 / 250 m isotropic, tpv23; 15 s '
             'term; restored per-fault mesh extent, each fault on its own true '
             'box) vs 100 m/15 s independent references -- a cross-code check at '
             'a finer, not resolution-matched, gate. python-jax not plotted (out of scope).']
    fig.text(0.01, 0.005, '\n'.join(lines), fontsize=6.2, va='top', family='monospace')
    fig.tight_layout(rect=[0, 0.07, 1, 0.94])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print('wrote %s' % out)
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--results', default=os.path.join(ROOT, 'test', 'test.%s'))
    # barall* and kaneko/payne both live, read-only, in one place now:
    # ~/shared_dataset/scec_cvws.tpv2223/raw/<case>/<label>/ (item 17
    # section C, 2026-10-02) -- see that bundle's MANIFEST.json/README row.
    _SHARED = os.path.expanduser(os.path.join(
        '~', 'shared_dataset', 'scec_cvws.tpv2223', 'raw'))
    ap.add_argument('--kaneko', default=_SHARED)
    ap.add_argument('--archive', default=_SHARED)
    ap.add_argument('--outdir', default=os.path.join(ROOT, 'docs', 'figures', 'tpv22_23'))
    args = ap.parse_args()

    for case in ('tpv22', 'tpv23'):
        rdir = args.results % case if '%s' in args.results else args.results
        kaneko_dir = os.path.join(args.kaneko, case, 'kaneko')
        barall_dir = os.path.join(args.archive, case, 'barall-faultmod-100m-2013')
        make_family(case, rdir, kaneko_dir, barall_dir, ON_FAULT, ON_FAULT_COLS,
                   os.path.join(args.outdir, '%s_ts_onfault_overlay.png' % case),
                   'on-fault')
        make_family(case, rdir, kaneko_dir, barall_dir, OFF_FAULT, OFF_FAULT_COLS,
                   os.path.join(args.outdir, '%s_ts_offfault_overlay.png' % case),
                   'off-fault')


if __name__ == '__main__':
    main()
