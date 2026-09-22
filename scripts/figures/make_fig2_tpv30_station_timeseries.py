#! /usr/bin/env python3
"""Figure 2 -- TPV30 on-fault station time series: EQdyna's committed 500 m
Fortran gate reference against the owner's own 2015 EQdyna v3.1 100 m SCEC
submission.  Slip, slip rate and shear stress, one panel per quantity per
station.

Run:  python3 scripts/figures/make_fig2_tpv30_station_timeseries.py
Out:  docs/figures/tpv30/fig2_tpv30_station_timeseries.png

Data paths, print geometry, framing, the station-selection rule and the shared
cached y-ranges all live in tpv30_figlib.py -- read its module docstring
first.  Same framing as fig1: a CROSS-VERSION CONSISTENCY CHECK, not a
validation and not an accuracy claim.

The four stations are EXACTLY coincident between the two runs (both
coordinates are multiples of 500 m), so no nearest-node substitution and no
interpolation is involved anywhere in this figure -- unlike
evidence_tpv30_scec_comparison.py's METRIC 3, which deliberately accepts a
nearest-node offset of up to 283 m in order to use all 24 spec stations.
"""
import argparse
import os

import numpy as np

import tpv30_figlib as fl

OUT_PNG = os.path.join(fl.OUT_DIR, 'fig2_tpv30_station_timeseries.png')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--recompute-scales', action='store_true')
    args = ap.parse_args()

    plt = fl.apply_style()
    sc = fl.shared_scales(recompute=args.recompute_scales)

    nrow, ncol = len(fl.STATIONS), len(fl.TS_PANELS)
    fig, axes = plt.subplots(nrow, ncol, figsize=(fl.CANVAS_WIDTH_IN, 6.1),
                             sharex=True, layout='constrained')

    final = []
    for i, st in enumerate(fl.STATIONS):
        ours, arch = fl.load_station_pair(st)
        s_m, d_m = fl.cmp29.station_coords_m(st)
        for j, (key, label) in enumerate(fl.TS_PANELS):
            ax = axes[i, j]
            ax.plot(arch['t'], arch[key], color='0.45', lw=0.9)
            ax.plot(ours['t'], ours[key], color='#b2182b', lw=0.8, ls='--')
            ax.set_xlim(*sc['ts_t_s'])
            ax.set_ylim(*sc['ts_' + key])
            ax.locator_params(axis='y', nbins=4)
            ax.locator_params(axis='x', nbins=5)
            ax.axhline(0.0, color='0.8', lw=0.4, zorder=0)
            if i == 0:
                ax.set_title(label, loc='center')
        # Row identity only on the left column (rule 6: collapse repeated
        # labels); the quantity and its unit live in the column header.
        axes[i, 0].set_ylabel(f'({s_m/1e3:+.0f}, {d_m/1e3:.0f}) km')
        final.append((st, ours['h_slip'][-1], arch['h_slip'][-1]))

    # One shared x-label, on the bottom-centre panel.  fig.supxlabel is NOT
    # used: the layout engine places it against the figure edge, below the
    # rect footnote() reserves, and it lands on top of the note.
    axes[-1, 1].set_xlabel('Time since nucleation (s)')
    fig.suptitle('TPV30 on-fault stations: EQdyna today (500 m) vs '
                 'EQdyna v3.1 (100 m, 2015)', fontsize=fl.PT['title'] * fl.K)

    LEGEND = ([plt.Line2D([], [], color='0.45', lw=0.9),
               plt.Line2D([], [], color='#b2182b', lw=0.8, ls='--')],
              ['2015, 100 m', 'ours, 500 m'])

    note = (
        'Fortran-vs-Fortran cross-version consistency check across ~a decade '
        'of code evolution -- NOT an independent validation and NOT an '
        'accuracy claim.\n'
        'Different node spacings (ours 500 m, 2015 100 m) and different time '
        'steps (0.0417 s vs 0.0080 s); curves are plotted at their own '
        'samples, nothing is resampled or interpolated.\n'
        'Stations: of the 24 spec stations, 13 fall exactly on the 500 m '
        'grid; one is shown per quadrant of the fault plane about the '
        'hypocentre.  Both coordinates coincide exactly -- zero offset.\n'
        'Row labels are (along-strike, down-dip) station coordinates. Every '
        'panel shows the ALONG-STRIKE (horizontal) component. Rows share x; '
        'each column shares one y-range taken from the union of both runs '
        'over all four stations (docs/figures/tpv30/shared_scales.json).\n'
        'The 2015 100 m submission ranks 14/14 against the portal group '
        'median (its own PROVENANCE.md), so agreement here certifies neither '
        'run against the benchmark.')
    ylegend = fl.footnote(fig, note, width_chars=104, extra_in=0.24)
    fig.legend(*LEGEND, loc='lower center', bbox_to_anchor=(0.5, ylegend),
               ncol=2, handlelength=1.8, columnspacing=1.6,
               handletextpad=0.5, borderaxespad=0.0)

    os.makedirs(fl.OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=fl.DPI)
    print(f'wrote {OUT_PNG}')
    print(fl.print_geometry_note())
    print('effective printed pt: labels %.1f, ticks %.1f, titles %.1f, '
          'legend %.1f, note %.1f'
          % tuple(v * fl.K for v in (fl.PT['label'], fl.PT['tick'],
                                     fl.PT['title'], fl.PT['legend'],
                                     fl.PT['annot'])))
    print('final along-strike slip (m): station / ours / 2015 / %diff')
    for st, a, b in final:
        print('  %-18s %8.3f %8.3f %8.1f%%'
              % (st, a, b, abs(a - b) / abs(b) * 100.0 if b else float('nan')))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
