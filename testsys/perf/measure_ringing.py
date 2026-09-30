#! /usr/bin/env python3
"""
measure_ringing.py -- post-front on-fault stress ringing, and what a damping
change does to the rupture (pathway_forward.md row 34).

WHAT IT MEASURES

For each on-fault station file it is given (default: the case's on-fault
matrix.GATE_STATIONS), on ONE shear-stress column (default `v-shear-stress`,
column 7, the down-dip shear the tpv36/tpv37 zigzags live in):

  t_arr   first time the station's slip-rate magnitude
          sqrt(h-slip-rate^2 + v-slip-rate^2) exceeds ARRIVAL_SLIPRATE
          (1e-3 m/s, the SCEC rupture-time convention). A station that never
          reaches it did not rupture; its ringing is reported as
          NOT-RUPTURED, never as 0.
  window  [t_arr + LEAD, min(t_arr + LEAD + LENGTH, t_end)], LEAD = 0.5 s,
          LENGTH = 2.0 s. LEAD excludes the stress drop itself: at dx = 500 m
          the grid-scale period is 1/(Vs/2dx) = 0.29 s and the drop at these
          stations completes inside ~0.3 s of arrival, so 0.5 s is past both.
          A window shorter than MIN_WINDOW (1.0 s) is refused, not measured
          (rule 2) -- the station ruptured too late in the record.
  detrend cubic polynomial least-squares fit over the window, subtracted.
          The residual keeps everything above roughly 1/LENGTH ~ 0.5 Hz,
          i.e. the 1.4-3.5 Hz mesh-scale band row 34 is about, and drops the
          slow post-drop stress evolution that continues while slip continues.
  ringing RMS and peak-to-peak (MPa) of the residual, and its dominant
          frequency (Hz) from the residual's spectrum.

The same three numbers are computed for the reference (committed 500 m gate
station) and for the archive (published 50 m SCEC station) when given, so the
reader sees what "flat" measures on this metric: the 50 m archive is the
noise floor a 500 m run would have to reach.

Also, per station, max |run - other| (MPa) over the whole record after
linear resampling of `other` onto the run's time axis.

RUPTURE (when --run and --case are given): the run's frt.txt* are
canonicalised with testsys.frt_canonical and aligned to the case's committed
frt.canonical.txt. Reported: nodes ruptured in run vs ref (a flip count),
max and mean |rupture-time change| over nodes ruptured in both, and the
fault-wide peak slip rate (frt column 10, FRIC_SLOT_SLIPRATE_MAX) for run
and ref with its ratio. A damping that lowers ringing while moving any of
these has not fixed anything -- that is the row's kill criterion.

BOTH OUTCOMES (rule 14a): with --max-p2p MPa the exit code is 1 if ANY
measurable station's ringing peak-to-peak exceeds the bound (UNFIXED) and 0
if none does (FIXED). A station that has not ruptured by the end of the
record (at the 5 s gate term two of tpv37's three on-fault gate stations have
not) is printed as NOT-RUPTURED and does not vote; if NO station is
measurable the exit code is 2 -- "could not check" is a third answer, not a
pass (rule 2). Without --max-p2p the tool only reports (exit 0).

--stations overrides the station list with files under --run (or under the
case's committed reference dir when there is no --run): use it to reach
on-fault stations that DO rupture inside the term but are not gate stations
(tpv36/37: faultst000dp120, faultst040dp180, faultst080dp180,
faultst000dp240). A station with no committed reference file says so on its
own line; nothing is skipped silently.

USAGE

  # one run directory against the committed reference and the 50 m archive
  python3 testsys/perf/measure_ringing.py --case test.tpv37 \\
      --run test/test.tpv37 \\
      --archive scec_archive/tpv37/eqdyna-v5.3.3-50m-2024-nstress-corrected

  # the committed reference alone, no run: how much does the GATE ring?
  python3 testsys/perf/measure_ringing.py --case test.tpv37 --max-p2p 0.2

  # arbitrary station files
  python3 testsys/perf/measure_ringing.py --station a/faultst000dp180.txt

The archive is READ-ONLY and irreplaceable (scec_archive/*/PROVENANCE.md);
this tool opens it for reading only.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO_ROOT)
from testsys import compare, frt_canonical, matrix  # noqa: E402

ARRIVAL_SLIPRATE = 1e-3   # m/s, SCEC rupture-time convention
LEAD = 0.5                # s after arrival before the window opens
LENGTH = 2.0              # s, window length
MIN_WINDOW = 1.0          # s, shortest window this tool will measure
DETREND_DEGREE = 3
DEFAULT_COLUMN = 'v-shear-stress'
FRT_COL_RUPT = 3          # frt column 4: rupture time (99999 = never)
FRT_COL_PEAK_SLIPRATE = 9  # frt column 10: FRIC_SLOT_SLIPRATE_MAX
FRT_NEVER = 9.0e4


class NotMeasurable(Exception):
    """The station cannot support the metric; the reason is the message."""


def load_station(path):
    names, data = compare.read_station_file(path)
    cols = {n: data[:, i] for i, n in enumerate(names)}
    for need in ('t', 'h-slip', 'v-slip', 'h-slip-rate', 'v-slip-rate'):
        if need not in cols:
            raise ValueError('%s: no %r column (have %s)' % (path, need, names))
    return cols


def arrival_time(cols):
    sr = np.hypot(cols['h-slip-rate'], cols['v-slip-rate'])
    hit = np.nonzero(sr > ARRIVAL_SLIPRATE)[0]
    if hit.size == 0:
        raise NotMeasurable('NOT-RUPTURED (slip rate never > %g m/s)' % ARRIVAL_SLIPRATE)
    final_slip = float(np.hypot(cols['h-slip'][-1], cols['v-slip'][-1]))
    return float(cols['t'][hit[0]]), float(sr.max()), final_slip


def ringing(cols, column):
    """(t_arr, peak_sliprate, final_slip, t0, t1, rms, p2p, f_dom) for one station."""
    t = cols['t']
    x = cols[column]
    t_arr, peak_sr, final_slip = arrival_time(cols)
    t0 = t_arr + LEAD
    t1 = min(t0 + LENGTH, float(t[-1]))
    if t1 - t0 < MIN_WINDOW:
        raise NotMeasurable('window [%.2f, %.2f] s shorter than %.1f s (arrived at %.2f s, '
                            'record ends %.2f s)' % (t0, t1, MIN_WINDOW, t_arr, t[-1]))
    sel = (t >= t0) & (t <= t1)
    tw, xw = t[sel], x[sel]
    if tw.size < 8:
        raise NotMeasurable('only %d samples in window' % tw.size)
    tc = tw - tw.mean()
    coef = np.polyfit(tc, xw, DETREND_DEGREE)
    res = xw - np.polyval(coef, tc)
    rms = float(np.sqrt(np.mean(res ** 2)))
    p2p = float(res.max() - res.min())
    dt = float(np.median(np.diff(tw)))
    spec = np.abs(np.fft.rfft(res * np.hanning(res.size)))
    freqs = np.fft.rfftfreq(res.size, dt)
    spec[0] = 0.0
    f_dom = float(freqs[int(np.argmax(spec))])
    return t_arr, peak_sr, final_slip, t0, t1, rms, p2p, f_dom


def max_abs_diff_resampled(run_cols, other_cols, column):
    t = run_cols['t']
    other = np.interp(t, other_cols['t'], other_cols[column])
    lo, hi = other_cols['t'][0], other_cols['t'][-1]
    inside = (t >= lo) & (t <= hi)
    if not inside.any():
        raise NotMeasurable('no time overlap with comparison series')
    return float(np.max(np.abs(run_cols[column][inside] - other[inside])))


def on_fault_gate_stations(case):
    st = matrix.GATE_STATIONS[case]
    if isinstance(st, dict):
        names = st.get('on_fault') or st.get('on') or st.get('fault')
        if names is None:
            names = [n for v in st.values() for n in v if str(n).startswith('faultst')]
    else:
        names = [n for n in st if str(n).startswith('faultst')]
    return [n if n.endswith('.txt') else n + '.txt' for n in names]


def archive_path(archive_dir, station_file):
    base = station_file[:-4] if station_file.endswith('.txt') else station_file
    for cand in (os.path.join(archive_dir, base), os.path.join(archive_dir, base + '.txt')):
        if os.path.isfile(cand):
            return cand
    raise FileNotFoundError('%s: no %s in archive' % (archive_dir, base))


def rupture_block(case, run_dir):
    ref = compare.load_reference(case)
    run = frt_canonical.canonical_from_case(run_dir)
    try:
        ref_a, run_a = frt_canonical.align(ref, run)
    except ValueError as e:
        # A different mesh (e.g. a dx=250 resolution control against the
        # dx=500 reference) has no node-by-node rupture comparison. Say so
        # explicitly; the caller prints it and the station numbers stand
        # on their own (rule 2: not a skip, a stated non-result).
        return {'not_comparable': 'run has %d fault nodes, reference %d: %s'
                % (run.shape[0], ref.shape[0], e)}
    rr, rn = ref_a[:, FRT_COL_RUPT], run_a[:, FRT_COL_RUPT]
    ref_rupt, run_rupt = rr < FRT_NEVER, rn < FRT_NEVER
    both = ref_rupt & run_rupt
    flips = int(np.sum(ref_rupt != run_rupt))
    d = np.abs(rn[both] - rr[both]) if both.any() else np.array([0.0])
    pk_ref = float(ref_a[:, FRT_COL_PEAK_SLIPRATE].max())
    pk_run = float(run_a[:, FRT_COL_PEAK_SLIPRATE].max())
    return {'nodes': int(ref_a.shape[0]), 'ref_ruptured': int(ref_rupt.sum()),
            'run_ruptured': int(run_rupt.sum()), 'flips': flips,
            'drt_max': float(d.max()), 'drt_mean': float(d.mean()),
            'peak_sr_ref': pk_ref, 'peak_sr_run': pk_run,
            'peak_sr_ratio': pk_run / pk_ref if pk_ref > 0 else float('nan')}


def fmt_ring(r):
    t_arr, peak_sr, final_slip, t0, t1, rms, p2p, f = r
    return ('t_arr=%.3f s  win=[%.2f,%.2f]  rms=%.4f MPa  p2p=%.4f MPa  f=%.2f Hz  '
            'peakSR=%.3f m/s  slip(t_end)=%.3f m' % (t_arr, t0, t1, rms, p2p, f, peak_sr, final_slip))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--case', help='gated case name, e.g. test.tpv37 (selects reference + gate stations)')
    ap.add_argument('--run', help='run directory holding station files and frt.txt*')
    ap.add_argument('--station', action='append', default=[], help='explicit station file(s)')
    ap.add_argument('--stations', help='comma-separated station names under --run, replacing the gate list')
    ap.add_argument('--archive', help='directory of published (e.g. 50 m) station files, read-only')
    ap.add_argument('--column', default=DEFAULT_COLUMN, help='stress column name (default %s)' % DEFAULT_COLUMN)
    ap.add_argument('--max-p2p', type=float, default=None,
                    help='exit 1 if any measured station ringing peak-to-peak (MPa) exceeds this')
    ap.add_argument('--label', default='run')
    args = ap.parse_args(argv)

    if not (args.case or args.station):
        ap.error('give --case and/or --station')

    # station list: explicit files, else the case's on-fault gate stations
    if args.station:
        targets = [(os.path.basename(p), p) for p in args.station]
    else:
        if args.stations:
            names = [n if n.endswith('.txt') else n + '.txt' for n in args.stations.split(',')]
        else:
            names = on_fault_gate_stations(args.case)
        if args.run:
            targets = [(n, os.path.join(args.run, n)) for n in names]
        else:
            targets = [(n, compare.station_reference_path(args.case, n)) for n in names]
            args.label = 'reference'

    print('measure_ringing: column=%s lead=%.2f s length=%.2f s detrend=deg%d arrival>%g m/s'
          % (args.column, LEAD, LENGTH, DETREND_DEGREE, ARRIVAL_SLIPRATE))
    worst_p2p, measured = 0.0, 0
    for name, path in targets:
        print('\n[%s] %s' % (name, path))
        run_cols = load_station(path)
        try:
            r = ringing(run_cols, args.column)
            worst_p2p = max(worst_p2p, r[6])
            measured += 1
            print('  %-10s %s' % (args.label, fmt_ring(r)))
        except NotMeasurable as e:
            print('  %-10s NOT MEASURABLE: %s' % (args.label, e))
        if args.case and args.run:
            try:
                ref_p = compare.station_reference_path(args.case, name)
            except FileNotFoundError:
                ref_p = None
                print('  %-10s no committed reference station file for %s (not a gate station)' % ('reference', name))
            if ref_p is not None:
                ref_cols = load_station(ref_p)
                try:
                    print('  %-10s %s' % ('reference', fmt_ring(ringing(ref_cols, args.column))))
                except NotMeasurable as e:
                    print('  %-10s NOT MEASURABLE: %s' % ('reference', e))
                print('  max|%s - reference| = %.4f MPa' % (args.label, max_abs_diff_resampled(run_cols, ref_cols, args.column)))
        if args.archive:
            arc_cols = load_station(archive_path(args.archive, name))
            try:
                print('  %-10s %s' % ('archive', fmt_ring(ringing(arc_cols, args.column))))
            except NotMeasurable as e:
                print('  %-10s NOT MEASURABLE: %s' % ('archive', e))
            print('  max|%s - archive|   = %.4f MPa' % (args.label, max_abs_diff_resampled(run_cols, arc_cols, args.column)))

    if args.case and args.run:
        b = rupture_block(args.case, args.run)
        if 'not_comparable' in b:
            print('\n[rupture] %s vs committed frt.canonical.txt: NOT COMPARABLE -- %s'
                  % (args.label, b['not_comparable']))
        else:
            print('\n[rupture] %s vs committed frt.canonical.txt (%d nodes)' % (args.label, b['nodes']))
            print('  ruptured: ref %d  run %d  flips %d' % (b['ref_ruptured'], b['run_ruptured'], b['flips']))
            print('  |d rupture time|: max %.4f s  mean %.4f s' % (b['drt_max'], b['drt_mean']))
            print('  fault peak slip rate: ref %.4f  run %.4f m/s  ratio %.4f'
                  % (b['peak_sr_ref'], b['peak_sr_run'], b['peak_sr_ratio']))

    if args.max_p2p is not None:
        if measured == 0:
            print('\nVERDICT: no measurable station (%d listed) -> COULD NOT CHECK' % len(targets))
            return 2
        bad = worst_p2p > args.max_p2p
        print('\nVERDICT: worst p2p %.4f MPa over %d measurable station(s) vs bound %.4f MPa -> %s'
              % (worst_p2p, measured, args.max_p2p, 'UNFIXED' if bad else 'FIXED'))
        return 1 if bad else 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
