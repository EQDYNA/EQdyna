#! /usr/bin/env python3
"""
scec_compare.py -- THE one SCEC/USGS TPV overlay tool: N result sets (our
runs/references, the scec_archive EQdyna submissions, other modellers' SCEC
CVWS submissions) for ANY benchmark, one plot mode per invocation.

    scec_compare.py --plot cplot      --models P1 P2 [...]   rupture-time contours
    scec_compare.py --plot res-series --models REF P2 [...]  REF vs each, one panel each
    scec_compare.py --plot ts-fault   --models P1 P2 [...]   on-fault time series
    scec_compare.py --plot ts-body    --models P1 P2 [...]   off-fault time series

    --fault N     multi-fault / branch cases (TPV22/23, TPV24/25): fault 1 or 2
                  (ours: frt split geometrically, faultstft2_* stations;
                  CVWS: cplot_N[.txt] / cplot_main|cplot_branch,
                  fault<N>st* / branchst*)

Each PATH's layout is sniffed (scec_readers.RUPTURE_READERS and the station
layouts), never declared. Styles/labels are assigned from the ORDER of
--models, so a model keeps its colour and dash pattern across every mode.

Layers: scec_readers.py (one registered reader per format, all normalised to
(strike, down-dip, t) and station-series records), scec_plots.py (one
function per mode), this file (argument parsing only). The format hazards
handled are listed in scec_readers.py's docstring. NEVER modifies its inputs.

PRINT GEOMETRY: Elsevier full width 190 mm (--print-width-mm) at k = 1.0, so
every font size is the printed point size; 400 dpi -> 2992 px across.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scec_readers import Model, header_element_size_m          # noqa: E402
import scec_plots                                                # noqa: E402


def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--plot', required=True,
                    choices=('cplot', 'res-series', 'ts-fault', 'ts-body'))
    ap.add_argument('--models', nargs='+', required=True,
                    help='ordered model paths (our run/reference dirs, archive '
                         'or CVWS submission dirs, or a single rupture file)')
    ap.add_argument('--labels', nargs='*', default=None,
                    help='optional per-model legend labels, same order')
    ap.add_argument('--fault', type=int, default=None,
                    help='fault number on a multi-fault / branch case')
    ap.add_argument('--component', nargs='+', default=['h'], choices=('h', 'v', 'n'),
                    help='h = along-strike, v = down-dip, n = fault-normal; '
                         'several give component-major columns')
    ap.add_argument('--max-stations', type=int, default=4)
    ap.add_argument('--stations', nargs='+', default=None,
                    help='plot exactly these stations (any layout prefix)')
    ap.add_argument('--sentinel', type=float, default=None,
                    help='override the never-ruptured fill value (default: '
                         'detected per file)')
    ap.add_argument('--contour-step', type=float, default=None,
                    help='contour interval, s (default: a nice step from the data)')
    ap.add_argument('--hypo', nargs=2, type=float, default=None,
                    metavar=('STRIKE_KM', 'DOWNDIP_KM'), help='star the hypocentre')
    ap.add_argument('--seeds', action='store_true',
                    help='caption the per-model seed count (8-neighbour)')
    ap.add_argument('--tau-strength', default=None, metavar='CASE_DIR',
                    help='caption t=0 tau/strength of a rough-fault case_input dir')
    ap.add_argument('--note', default=None,
                    help='one caption line appended verbatim (e.g. the run term)')
    ap.add_argument('--out', default=None)
    ap.add_argument('--print-width-mm', type=float, default=190.0,
                    help='printed width: Elsevier 90/140/190 mm, AGU 146 mm')
    ap.add_argument('--dpi', type=int, default=400)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.labels and len(args.labels) != len(args.models):
        raise SystemExit('--labels must have the same length as --models')
    args.width_in = args.print_width_mm / 25.4          # k = 1.0 by construction
    # DejaVu Sans averages ~3.95 pt per character; wrap from the PRINTED width
    args.wrap_chars = max(40, int(args.width_in * 72.0 / 3.95 * 0.96))
    if args.out is None:
        args.out = f'scec_compare_{args.plot.replace("-", "_")}.png'

    models = [Model(p, i, (args.labels[i] if args.labels else None),
                    args.sentinel, args.fault)
              for i, p in enumerate(args.models)]

    print('=== models (style assigned by position in --models) ===')
    for m in models:
        print(f'  [{m.index}] {m.kind:8s} {m.path}  reader={m.reader}')
        print(f'        label="{m.label()}"  colour={m.color}  '
              f'linestyle={m.linestyle}')
        if m.rupture_file and os.path.isfile(m.rupture_file):
            es = header_element_size_m(m.rupture_file)
            if es is not None:
                print(f'        header element_size = {es:g} m (reported only; '
                      f'spacing used in figures is measured from coordinates)')

    if args.plot == 'cplot':
        scec_plots.plot_cplot(models, args)
    elif args.plot == 'res-series':
        scec_plots.plot_res_series(models, args)
    else:
        scec_plots.plot_ts(models, args.plot.split('-')[1], args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
