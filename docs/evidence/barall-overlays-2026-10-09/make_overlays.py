#! /usr/bin/env python3
"""Rebuild this directory's FaultMod CVWS overlays and the verdict numbers.

    python3 make_overlays.py figures RESULTS OUTDIR [case ...]
    python3 make_overlays.py numbers RESULTS [case ...]

RESULTS/<case>/ holds one fortran run's station files plus frt.canonical.txt
(python3 -m testsys.frt_canonical <run_dir>). The FaultMod submissions are read,
read-only, from the shared-dataset store ($SCEC_CVWS, default
~/shared_dataset). Figures come from scripts/figures/scec_compare.py only.

`numbers` is the verdict table in README.md. It is computed on EXACTLY-shared
node coordinates and station names, with no spatial interpolation. One
exception: the TPV24/25 branch (fault 2) shares no node, so it is matched
to the nearest FaultMod node within 111 m, measured in (along-branch from the
run's own junction, down-dip) coordinates. Station arrival is the first time
|h-slip-rate| > 1 mm/s on the fault, or |h-velocity| > 1 cm/s off the fault.
"""
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
FIG = os.path.join(ROOT, 'scripts', 'figures')
sys.path.insert(0, FIG)
SD = os.environ.get('SCEC_CVWS', os.path.expanduser('~/shared_dataset'))

FAULTMOD = {
    'tpv12': 'scec_cvws.tpv1213/raw/tpv12/barall-faultmod-100m-2009',
    'tpv13': 'scec_cvws.tpv1213/raw/tpv13/barall-faultmod-100m-2009',
    'tpv26': 'scec_cvws.tpv2627/raw/tpv26/barall-faultmod-100m-2013',
    'tpv27': 'scec_cvws.tpv2627/raw/tpv27/barall-faultmod-100m-2013',
    'tpv31': 'scec_cvws.tpv3132/raw/tpv31/barall-faultmod-50m-2014',
    'tpv32': 'scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014',
    'tpv33': 'scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015',
    'tpv24': 'scec_cvws.tpv2425/raw/tpv24/barall-faultmod-100m-2013',
    'tpv25': 'scec_cvws.tpv2425/raw/tpv25/barall-faultmod-100m-2013',
}
DX = {'tpv33': 400, 'tpv24': 250}            # every other case: 500 m
TERM = {'tpv24': 12.0, 'tpv25': 12.0}        # every other case: the 5 s gate term
JUNCTION_X = {'tpv24': 750, 'tpv25': 500}    # measured from the runs' frt


def faults(c):
    return [1, 2] if c in ('tpv24', 'tpv25') else [None]


def note(c, f):
    dx, term = DX.get(c, 500), TERM.get(c, 5.0)
    s = (f'EQdyna fortran, gate config dx {dx} m, run to {term:g} s '
         + ('(5 s gate term; spec term is longer)' if term == 5.0
            else '(branch-feature run, not a gated case yet)'))
    if f == 2:
        s += (f'; NOTE this run puts the branch junction at x = {JUNCTION_X[c]} m, '
              'not the spec x = 0 (WIP case: branch fxmin fixed at 1000 m while '
              'the wedge line starts at fxmin - dx), and the branch is decimated '
              "to a 10 km x-extent; along-strike is measured from this run's junction")
    return s


def figures(results, out, cases):
    os.makedirs(out, exist_ok=True)
    for c in cases:
        for f in faults(c):
            for mode in ('cplot', 'ts-fault', 'ts-body'):
                if f == 2 and mode != 'cplot':
                    continue            # no branch station is comparable; see README
                tag = c + (f'_fault{f}' if f else '')
                cmd = [sys.executable, os.path.join(FIG, 'scec_compare.py'), '--plot', mode]
                if mode != 'cplot':
                    cmd += ['--component', 'h', 'v', 'n']
                if f:
                    cmd += ['--fault', str(f)]
                cmd += ['--models', os.path.join(results, c), os.path.join(SD, FAULTMOD[c]),
                        '--dpi', '200', '--note', note(c, f),
                        '--out', os.path.join(out, f'{tag}_{mode.replace("-", "_")}.png')]
                r = subprocess.run(cmd, capture_output=True, text=True)
                print(c, f, mode, 'rc', r.returncode)
                if r.returncode:
                    print(r.stderr[-400:])


def _arrival(a, thr):
    i = np.nonzero(np.abs(a[:, 2]) > thr)[0]
    return a[i[0], 0] if i.size else np.nan


def _dt_stats(to, tb, term):
    rb = np.isfinite(tb) & (tb < term - 0.25)        # FaultMod arrivals inside our term
    ro = np.isfinite(to)
    both = rb & ro
    early = ro & ~(np.isfinite(tb) & (tb < term + 0.25))
    dt = (to - tb)[both]
    s = (f'FaultMod ruptured {rb.sum()}, ours also {both.sum()} '
         f'({100.0 * both.sum() / max(1, rb.sum()):.1f}%), ours ruptured where FaultMod '
         f'is not by term+0.25 s: {early.sum()}')
    if dt.size:
        s += (f'; |dt| median {np.median(abs(dt)):.3f} p90 {np.percentile(abs(dt), 90):.3f} '
              f'max {abs(dt).max():.3f} s, mean(ours-FaultMod) {dt.mean():+.3f} s')
    return s


def numbers(results, cases):
    import scec_readers as R
    from scipy.spatial import cKDTree
    for c in cases:
        term = TERM.get(c, 5.0)
        for f in faults(c):
            o = R.Model(os.path.join(results, c), 0, fault=f)
            b = R.Model(os.path.join(SD, FAULTMOD[c]), 1, fault=f)
            fo, fb = o.load_rupture(), b.load_rupture()
            tag = c + (f' fault {f}' if f else '')
            if f == 2:
                B = np.column_stack([fb['strike'], fb['downdip']])
                d, i = cKDTree(B).query(np.column_stack([fo['strike'], fo['downdip']]))
                print(f'## {tag} (nearest FaultMod node, max {d.max():.0f} m): '
                      + _dt_stats(fo['t'], fb['t'][i], term))
                continue
            kb = {(round(x, 1), round(y, 1)): t for x, y, t in zip(fb['strike'], fb['downdip'], fb['t'])}
            keys = [(round(x, 1), round(y, 1)) for x, y in zip(fo['strike'], fo['downdip'])]
            sh = np.array([k in kb for k in keys])
            to = fo['t'][sh]
            tb = np.array([kb[k] for k, s in zip(keys, sh) if s])
            print(f'## {tag} (shared nodes {sh.sum()}): ' + _dt_stats(to, tb, term))
            if f == 1:
                xs = fo['strike'][sh]
                for lab, sel in (('x<0', xs < 0), ('x>=0', xs >= 0)):
                    print(f'   main fault {lab}: ' + _dt_stats(to[sel], tb[sel], term))
            dx = fo['dx_m']
            for which, thr in (('fault', 1e-3), ('body', 1e-2)):
                po, pb = o.station_files(which), b.station_files(which)
                for k in sorted(set(po) & set(pb)):
                    if any(abs(v) % dx > 1e-6 and abs(abs(v) % dx - dx) > 1e-6 for v in k):
                        continue        # off our node grid: not comparable
                    ao, _ = R.load_station(po[k], which)
                    ab, _ = R.load_station(pb[k], which)
                    ao, abt = ao[ao[:, 0] <= term], ab[ab[:, 0] <= term]
                    line = f'   {R.station_name(k, which):22s} arrival {_arrival(ao, thr):6.2f} vs {_arrival(ab, thr):6.2f} s'
                    if which == 'fault':
                        line += (f'  slip@term {ao[-1, 1]:+.3f} vs {np.interp(ao[-1, 0], ab[:, 0], ab[:, 1]):+.3f} m'
                                 f'  peak rate {abs(ao[:, 2]).max():.3g} vs {abs(abt[:, 2]).max():.3g} m/s')
                    else:
                        line += f'  peak vel {abs(ao[:, 2]).max():.3g} vs {abs(abt[:, 2]).max():.3g} m/s'
                    print(line)


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('figures', 'numbers'):
        raise SystemExit(__doc__)
    mode, results = sys.argv[1], os.path.abspath(sys.argv[2])
    if mode == 'figures':
        figures(results, os.path.abspath(sys.argv[3]), sys.argv[4:] or list(FAULTMOD))
    else:
        numbers(results, sys.argv[3:] or list(FAULTMOD))
