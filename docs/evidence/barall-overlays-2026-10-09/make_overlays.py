#! /usr/bin/env python3
"""Rebuild this directory's FaultMod CVWS overlays and the verdict numbers.

    python3 make_overlays.py figures RESULTS OUTDIR [case ...]
    python3 make_overlays.py numbers RESULTS [case ...]
    python3 make_overlays.py resolution RUNS OUTDIR [--no-figures] [case ...]

RUNS/<case>_dx<dx>/ (resolution mode) are runs of one case at several dx,
each run to its spec term; see resolution/README.md. Resolution mode times
station arrivals on the VECTOR rate (hypot of h and v slip rate on the fault,
max |velocity| over three components off it); `numbers` keeps h-only.

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


def _dt_metrics(to, tb, term):
    """The verdict numbers as a dict. Our arrivals after `term` are masked,
    so a run longer than the window is compared on the same window as one
    stopped at it (a run stopped at term has none, so the gate-term numbers
    are unchanged)."""
    rb = np.isfinite(tb) & (tb < term - 0.25)        # FaultMod arrivals inside our term
    ro = np.isfinite(to) & (to <= term)
    both = rb & ro
    early = ro & ~(np.isfinite(tb) & (tb < term + 0.25))
    dt = (to - tb)[both]
    m = dict(faultmod_ruptured=int(rb.sum()), both=int(both.sum()),
             ruptured_pct=100.0 * both.sum() / max(1, rb.sum()), ours_only=int(early.sum()))
    if dt.size:
        m.update(median=float(np.median(abs(dt))), p90=float(np.percentile(abs(dt), 90)),
                 max=float(abs(dt).max()), mean=float(dt.mean()))
    return m


def _dt_stats(to, tb, term):
    m = _dt_metrics(to, tb, term)
    s = (f'FaultMod ruptured {m["faultmod_ruptured"]}, ours also {m["both"]} '
         f'({m["ruptured_pct"]:.1f}%), ours ruptured where FaultMod '
         f'is not by term+0.25 s: {m["ours_only"]}')
    if 'median' in m:
        s += (f'; |dt| median {m["median"]:.3f} p90 {m["p90"]:.3f} '
              f'max {m["max"]:.3f} s, mean(ours-FaultMod) {m["mean"]:+.3f} s')
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


# ---------------------------------------------------------------------------
# resolution study (resolution/README.md): the same numbers, per dx, per window
# ---------------------------------------------------------------------------
SPEC_TERM = {'tpv12': 8.0, 'tpv13': 8.0, 'tpv26': 13.0, 'tpv27': 13.0,
             'tpv31': 15.0, 'tpv32': 15.0, 'tpv33': 13.0}
GATE_DX = {c: DX.get(c, 500) for c in SPEC_TERM}


def _arrival_vec(a, which, thr):
    """Arrival on the vector rate: |slip rate| (h, v) on the fault, max |v|
    over the three components off it. The gate-term `numbers` mode reads the
    horizontal rate only, which is blind on a dip-slip fault (TPV12/13)."""
    r = np.hypot(a[:, 2], a[:, 5]) if which == 'fault' else np.abs(a[:, [2, 4, 6]]).max(1)
    i = np.nonzero(r > thr)[0]
    return a[i[0], 0] if i.size else np.nan


def _on_grid(k, dx):
    return all(abs(v) % dx < 1e-6 or abs(abs(v) % dx - dx) < 1e-6 for v in k)


def case_metrics(run_dir, c, term, grid_dx):
    """Verdict numbers for one run of case `c` on the window [0, term].
    Stations are restricted to the `grid_dx` node grid (the gate dx), so every
    dx of one case is scored on the same station set."""
    import scec_readers as R
    o = R.Model(run_dir, 0)
    b = R.Model(os.path.join(SD, FAULTMOD[c]), 1)
    fo, fb = o.load_rupture(), b.load_rupture()
    kb = {(round(x, 1), round(y, 1)): t for x, y, t in zip(fb['strike'], fb['downdip'], fb['t'])}
    keys = [(round(x, 1), round(y, 1)) for x, y in zip(fo['strike'], fo['downdip'])]
    sh = np.array([k in kb for k in keys])
    tb = np.array([kb[k] for k, s in zip(keys, sh) if s])
    m = dict(case=c, run=os.path.basename(run_dir), dx=float(fo['dx_m']), term=term,
             shared_nodes=int(sh.sum()), **_dt_metrics(fo['t'][sh], tb, term))
    st, one_only = [], []
    for which, thr in (('fault', 1e-3), ('body', 1e-2)):
        po, pb = o.station_files(which), b.station_files(which)
        for k in sorted(set(po) & set(pb)):
            if not (_on_grid(k, grid_dx) and _on_grid(k, fo['dx_m'])):
                continue
            ao, _ = R.load_station(po[k], which)
            ab, _ = R.load_station(pb[k], which)
            ao, abt = ao[ao[:, 0] <= term], ab[ab[:, 0] <= term]
            r = dict(name=R.station_name(k, which), kind=which,
                     t_ours=_arrival_vec(ao, which, thr), t_fm=_arrival_vec(abt, which, thr))
            if which == 'fault':
                sb = np.hypot(ab[:, 1], ab[:, 4])
                r.update(slip_ours=float(np.hypot(ao[-1, 1], ao[-1, 4])),
                         slip_fm=float(np.interp(ao[-1, 0], ab[:, 0], sb)),
                         tauh0_ours=float(ao[0, 3]), tauh0_fm=float(ab[0, 3]),
                         tauv0_ours=float(ao[0, 6]), tauv0_fm=float(ab[0, 6]),
                         sn0_ours=float(ao[0, 7]), sn0_fm=float(ab[0, 7]))
            else:
                r.update(pv_ours=float(abs(ao[:, [2, 4, 6]]).max()), pv_fm=float(abs(abt[:, [2, 4, 6]]).max()))
            fo_ok, fb_ok = np.isfinite(r['t_ours']), np.isfinite(r['t_fm'])
            # one code only: one arrives before term - 0.25 s, the other not by term
            if (fo_ok and not fb_ok and r['t_ours'] < term - 0.25) or \
               (fb_ok and not fo_ok and r['t_fm'] < term - 0.25):
                one_only.append(r['name'])
            st.append(r)
    m['stations'] = st
    m['one_code_only'] = one_only
    fa = [r for r in st if r['kind'] == 'fault' and np.isfinite(r['t_ours']) and np.isfinite(r['t_fm'])]
    bo = [r for r in st if r['kind'] == 'body' and r['pv_fm'] > 1e-3]
    m['fault_station_dt_median'] = float(np.median([r['t_ours'] - r['t_fm'] for r in fa])) if fa else None
    m['fault_slip_ratio_median'] = float(np.median([r['slip_ours'] / r['slip_fm'] for r in fa
                                                     if abs(r['slip_fm']) > 1e-3])) if fa else None
    m['body_pv_ratio_median'] = float(np.median([r['pv_ours'] / r['pv_fm'] for r in bo])) if bo else None
    m['match'] = bool(m['ruptured_pct'] >= 95.0 and m.get('median', 9) <= 0.2
                      and m.get('p90', 9) <= 0.5 and not one_only)
    return m


def _runs(runs, c):
    """RUNS/<case>_dx<dx>/ directories of case c, coarsest first."""
    out = []
    for d in os.listdir(runs):
        if d.startswith(c + '_dx') and os.path.isfile(os.path.join(runs, d, 'frt.canonical.txt')):
            rest = d[len(c) + 3:]
            if rest.replace('.', '').isdigit():
                out.append((float(rest), os.path.join(runs, d)))
    return sorted(out, reverse=True)


def _f(v, fmt):
    return '—' if v is None or (isinstance(v, float) and not np.isfinite(v)) else format(v, fmt)


def _nan_to_none(v):
    """NaN (no arrival) -> JSON null, so metrics.json is strict JSON."""
    if isinstance(v, dict):
        return {k: _nan_to_none(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_nan_to_none(x) for x in v]
    if isinstance(v, (float, np.floating)) and not np.isfinite(v):
        return None
    return v


def resolution(runs, out, cases, do_figures=True):
    """Per case and per window (5 s gate term, spec term): metrics for every
    dx found under RUNS, a markdown table, metrics.json, and one overlay set
    per case with every dx plus FaultMod."""
    import json
    os.makedirs(out, exist_ok=True)
    allm, md = [], []
    for c in cases:
        rs = _runs(runs, c)
        if not rs:
            print(c, 'no runs under', runs)
            continue
        md += [f'### {c}', '',
               '| window | dx (m) | ruptured | median \\|dt\\| | p90 | mean | max | ours-only nodes '
               '| one-code-only stations | fault-st dt median | slip ratio | body pv ratio | MATCH |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        for term in sorted({5.0, SPEC_TERM[c]}):
            for dx, d in rs:
                m = case_metrics(d, c, term, GATE_DX[c])
                allm.append(m)
                md.append(f'| {term:g} s | {dx:g} | {m["ruptured_pct"]:.1f}% ({m["both"]}/{m["faultmod_ruptured"]}) '
                          f'| {_f(m.get("median"), ".3f")} | {_f(m.get("p90"), ".3f")} '
                          f'| {_f(m.get("mean"), "+.3f")} | {_f(m.get("max"), ".2f")} | {m["ours_only"]} '
                          f'| {len(m["one_code_only"])} {" ".join(m["one_code_only"])} '
                          f'| {_f(m["fault_station_dt_median"], "+.2f")} | {_f(m["fault_slip_ratio_median"], ".2f")} '
                          f'| {_f(m["body_pv_ratio_median"], ".2f")} | {"yes" if m["match"] else "no"} |')
                print(md[-1], flush=True)
        md.append('')
        if do_figures:
            spec = SPEC_TERM[c]
            for mode in ('cplot', 'ts-fault', 'ts-body'):
                cmd = [sys.executable, os.path.join(FIG, 'scec_compare.py'), '--plot', mode]
                if mode != 'cplot':
                    cmd += ['--component', 'h', 'v', 'n']
                cmd += ['--models'] + [d for _, d in rs] + [os.path.join(SD, FAULTMOD[c])]
                cmd += ['--labels'] + ['EQdyna'] * len(rs) + ['FaultMod']   # the tool appends ', dx N m'
                cmd += ['--dpi', '200', '--note',
                        f'EQdyna fortran, 4 ranks, scratch runs at dx '
                        f'{", ".join("%g" % dx for dx, _ in rs)} m, each run to the spec term '
                        f'{spec:g} s (not gated cells)',
                        '--out', os.path.join(out, f'{c}_{mode.replace("-", "_")}_res.png')]
                r = subprocess.run(cmd, capture_output=True, text=True)
                print(c, mode, 'rc', r.returncode, flush=True)
                if r.returncode:
                    print(r.stderr[-400:])
    with open(os.path.join(out, 'metrics.json'), 'w') as fh:
        json.dump(_nan_to_none(allm), fh, indent=1, default=float)
    with open(os.path.join(out, 'table.md'), 'w') as fh:
        fh.write('\n'.join(md) + '\n')


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('figures', 'numbers', 'resolution'):
        raise SystemExit(__doc__)
    mode, results = sys.argv[1], os.path.abspath(sys.argv[2])
    if mode == 'figures':
        figures(results, os.path.abspath(sys.argv[3]), sys.argv[4:] or list(FAULTMOD))
    elif mode == 'resolution':
        rest = [a for a in sys.argv[4:] if a != '--no-figures']
        resolution(results, os.path.abspath(sys.argv[3]), rest or list(SPEC_TERM),
                   '--no-figures' not in sys.argv)
    else:
        numbers(results, sys.argv[3:] or list(FAULTMOD))
