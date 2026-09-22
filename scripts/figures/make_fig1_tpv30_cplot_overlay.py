#! /usr/bin/env python3
"""Figure 1 -- TPV30 rupture-time (cplot) overlay: EQdyna's committed 500 m
Fortran gate reference over the owner's own 2015 EQdyna v3.1 100 m SCEC
submission, on one set of axes.

Run:  python3 scripts/figures/make_fig1_tpv30_cplot_overlay.py
Out:  docs/figures/tpv30/fig1_tpv30_cplot_overlay.png

Data paths, print geometry, framing and the shared cached scales all live in
tpv30_figlib.py -- read its module docstring first; in particular this is a
CROSS-VERSION CONSISTENCY CHECK, not a validation and not an accuracy claim,
and the two runs are at different node spacings with nothing interpolated.

CONTOUR INTERVAL, stated because it was not free to choose: the TPV29/30
specification (TPV29_30_Description_v06, Part 11, p.38) defines the cplot FILE
FORMAT only -- three columns, strike / down-dip / rupture time -- and
prescribes no plot interval; the archived cplot is raw nodal data and carries
none either.  So the figure shows BOTH: 0.5 s filled bands (the interval SCEC
contour plots conventionally use) carry the 2015 field, and the two line sets
that must be read against each other are drawn at 2.0 s, which is what stays
legible at 146 mm (at 1.0 s the coarse run's front closes into unreadable
blobs where it decelerates).  Panel (b) is node-by-node and does not depend on either
choice.
"""
import argparse
import os

import numpy as np

import tpv30_figlib as fl
from testsys.parity import evidence_tpv30_scec_comparison as cmp30

OUT_PNG = os.path.join(fl.OUT_DIR, 'fig1_tpv30_cplot_overlay.png')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--recompute-scales', action='store_true')
    args = ap.parse_args()

    plt = fl.apply_style()
    sc = fl.shared_scales(recompute=args.recompute_scales)
    p, ref, x15, dip15, T15 = fl.load_fields()
    d = fl.colocated_difference(ref, x15, dip15, T15)

    t_ours = fl.masked_times(ref['rupt'], cmp30.SENTINEL_CURRENT)
    t_2015 = fl.masked_times(T15, cmp30.SENTINEL_ARCHIVE)
    xo, yo = ref['x'] / 1e3, ref['dip'] / 1e3          # km, dip positive down
    xa, ya = x15 / 1e3, dip15 / 1e3
    hyp = (p['xsource'] / 1e3, -p['zsource'] / 1e3)

    tlo, thi = sc['rupture_time_s']
    band = sc['rupture_time_band_interval_s']
    line = sc['rupture_time_contour_interval_s']
    dlo, dhi = sc['dt_s']

    fig, (axa, axb) = plt.subplots(
        2, 1, figsize=(fl.CANVAS_WIDTH_IN, 6.5), layout='constrained')

    def attach_cbar(ax, mappable, label, ticks, extend='neither'):
        """Colour bar tied to the AXES bbox, not to the subplot slot: these
        panels are aspect='equal' and therefore much shorter than their slot,
        and a constrained-layout `fig.colorbar(ax=...)` would stretch the bar
        over the whole slot height (it did, on the first render)."""
        cax = ax.inset_axes([1.022, 0.0, 0.020, 1.0])
        cb = fig.colorbar(mappable, cax=cax, extend=extend)
        cb.set_label(label, fontsize=fl.PT['label'] * fl.K)
        cb.set_ticks(ticks)
        cb.ax.tick_params(labelsize=fl.PT['tick'] * fl.K, width=0.5, length=2)
        cb.outline.set_linewidth(0.5)
        return cb

    # ---- (a) overlay -------------------------------------------------------
    cf = axa.contourf(xa, ya, t_2015, levels=np.arange(tlo, thi + band, band),
                      cmap='viridis', extend='neither')
    lv = np.arange(line, thi, line)
    axa.contour(xa, ya, t_2015, levels=lv, colors='white',
                linewidths=0.55, linestyles='dashed')
    cnow = axa.contour(xo, yo, t_ours, levels=lv, colors='black',
                       linewidths=0.75)
    axa.clabel(cnow, lv[::2], fmt='%.0f', fontsize=fl.PT['annot'] * fl.K,
               inline=True, inline_spacing=2)
    axa.plot(*hyp, marker='*', ms=7, mfc='red', mec='k', mew=0.4, ls='none')

    attach_cbar(axa, cf, 'Rupture time (s)', fl.endpoint_ticks(tlo, thi))

    axa.set_title(f'(a)  Rupture-time contours every {line:.0f} s '
                  f'(bands: 2015, {band} s)', loc='left')

    # ---- (b) node-by-node difference --------------------------------------
    pm = axb.pcolormesh(xo, yo, d['dt'], cmap='RdBu_r', vmin=dlo, vmax=dhi,
                        shading='nearest', rasterized=True)
    ii, jj = np.where(d['only_ours'])
    axb.plot(xo[jj], yo[ii], ls='none', marker='^', ms=2.2, mfc='none',
             mec='k', mew=0.45)
    ii, jj = np.where(d['only_2015'])
    axb.plot(xo[jj], yo[ii], ls='none', marker='v', ms=2.2, mfc='none',
             mec='0.25', mew=0.45)
    axb.plot(*hyp, marker='*', ms=7, mfc='red', mec='k', mew=0.4, ls='none')

    attach_cbar(axb, pm, r'$\Delta t$ = ours $-$ 2015 (s)',
                fl.endpoint_ticks(dlo, dhi, mid=0.0), extend='both')

    axb.set_title(f'(b)  $\\Delta t$ at the {sc["n_both"]} nodes ruptured in '
                  f'both (blue: ours earlier)', loc='left')

    # ONE figure-level legend, placed outside the axes.  A legend box inside
    # either panel necessarily hides data (panel b is filled edge to edge), so
    # the space is reserved instead of stolen -- rule 7, layout before
    # anything else.
    LEGEND = ([plt.Line2D([], [], color='k', lw=0.75),
                plt.Line2D([], [], color='0.35', lw=0.55, ls='--'),
                plt.Line2D([], [], color='red', marker='*', ls='none', ms=5,
                           mec='k', mew=0.4),
                plt.Line2D([], [], ls='none', marker='^', ms=3.4, mfc='none',
                           mec='k', mew=0.45),
                plt.Line2D([], [], ls='none', marker='v', ms=3.4, mfc='none',
                           mec='0.25', mew=0.45)],
               ['ours, 500 m', '2015, 100 m', 'hypocentre',
                f'ruptured only in ours ({sc["n_only_ours"]})',
                f'only in 2015 ({sc["n_only_2015"]})'])

    for ax in (axa, axb):
        ax.set_xlim(xo.min(), xo.max())
        ax.set_ylim(yo.max(), 0.0)          # depth positive down, 0 at top
        ax.set_aspect('equal')
        ax.set_ylabel('Distance down-dip (km)')
        ax.locator_params(axis='x', nbins=9)
        ax.locator_params(axis='y', nbins=5)
    axa.tick_params(labelbottom=False)
    axb.set_xlabel('Distance along strike (km)')

    fig.suptitle('TPV30 rupture time: EQdyna today (500 m) vs EQdyna v3.1 '
                 '(100 m, 2015)', fontsize=fl.PT['title'] * fl.K)

    note = (
        f'Fortran-vs-Fortran cross-version consistency check across ~a decade '
        f'of code evolution -- NOT an independent validation and NOT an '
        f'accuracy claim.\n'
        f'Different node spacings: ours {ref["dx"]:.0f} m '
        f'({ref["rupt"].shape[1]}$\\times${ref["rupt"].shape[0]} fault nodes), '
        f'2015 {float(x15[1]-x15[0]):.0f} m ({x15.size}$\\times${dip15.size}). '
        f'Nothing is interpolated: each field is contoured on its own grid, '
        f'and all {d["n_nodes"]} of our nodes fall exactly on 2015 nodes '
        f'(500 divides 100), so (b) is a direct node-to-node difference.\n'
        f'Extent overlap {sc["extent_overlap_pct"]:.1f}% of 2015-ruptured '
        f'nodes; median |dt| {sc["dt_median_s"]:.2f} s, median signed dt '
        f'{np.nanmedian(d["dt"]):+.2f} s. Colour limit $\\pm${dhi:.1f} s '
        f'covers 90% of them, not the {sc["dt_max_s"]:.1f} s maximum; '
        f'{sc["dt_saturated_pct"]:.0f}% saturate (arrows).\n'
        f'The 2015 100 m submission ranks 14/14 against the portal group '
        f'median (its own PROVENANCE.md), so agreement here certifies neither '
        f'run against the benchmark.')
    ylegend = fl.footnote(fig, note, width_chars=104, extra_in=0.24)
    fig.legend(*LEGEND, loc='lower center', bbox_to_anchor=(0.5, ylegend),
               ncol=5, handlelength=1.4, columnspacing=1.0,
               handletextpad=0.4, borderaxespad=0.0)

    os.makedirs(fl.OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=fl.DPI)
    print(f'wrote {OUT_PNG}')
    print(fl.print_geometry_note())
    print('effective printed pt: labels %.1f, ticks %.1f, titles %.1f, '
          'legend %.1f, note %.1f'
          % tuple(v * fl.K for v in (fl.PT['label'], fl.PT['tick'],
                                     fl.PT['title'], fl.PT['legend'],
                                     fl.PT['annot'])))
    print('extent overlap %.1f%% | dt median %.3f s | only-ours %d | '
          'only-2015 %d'
          % (sc['extent_overlap_pct'],
             float(np.nanmedian(np.abs(d['dt']))),
             sc['n_only_ours'], sc['n_only_2015']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
