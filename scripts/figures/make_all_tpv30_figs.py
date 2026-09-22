#! /usr/bin/env python3
"""Driver: regenerate the complete TPV30 SCEC-comparison figure set.

    python3 scripts/figures/make_all_tpv30_figs.py [--recompute-scales]

Recomputes the shared scales FIRST (so both figures are drawn against the same
cached ranges -- rule 5: after changing a shared range, every figure that
shares it is regenerated), then runs each make_figN_*.py in order.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = ('make_fig1_tpv30_cplot_overlay.py',
        'make_fig2_tpv30_station_timeseries.py')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--recompute-scales', action='store_true',
                    help='rebuild docs/figures/tpv30/shared_scales.json from '
                         'the union of both datasets before drawing')
    args = ap.parse_args()

    if args.recompute_scales:
        sys.path.insert(0, HERE)
        import tpv30_figlib as fl
        fl.shared_scales(recompute=True)
        print(f'rebuilt {fl.SCALES_JSON}')

    rc = 0
    for f in FIGS:
        print(f'\n--- {f} ---')
        rc |= subprocess.call([sys.executable, os.path.join(HERE, f)])
    return rc


if __name__ == '__main__':
    raise SystemExit(main())
