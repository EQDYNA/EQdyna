#! /usr/bin/env python3
"""
scec_plots.py -- the PLOTTING LAYER of the SCEC overlay tooling: one
function per plot mode, all fed by scec_readers.Model records.

    plot_cplot       rupture-time contours, N models on one axes
    plot_res_series  rupture-time contours, model 0 against each other model
                     in its own panel (a resolution or code series), shared
                     contour levels
    plot_ts          station time series, on-fault ('fault') or off-fault
                     ('body'), stations in rows and quantities in columns

EVERY NUMBER IN EVERY CAPTION IS COMPUTED FROM THE ARRAYS BEING PLOTTED, at
plot time. Nothing here is hand-typed.

NO INTERPOLATION ACROSS RESOLUTIONS, EVER. Each model is contoured on its own
grid; quantitative comparison uses EXACT coordinate matches only, and only
stations that land exactly on every model's node grid are plotted -- the
excluded ones are counted and explained. Nearest-node statistics (the retired
one-off TPV29/30 overlay's median |dt|) are therefore NOT offered.

PRINT GEOMETRY: laid out at k = canvas_width / print_width = 1.0, so every
font size is the printed point size: labels 9, ticks 7.5, titles 9, legend
7.5, caption 7.
"""
import importlib.util
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scec_readers import FAULT_COLS, BODY_COLS, station_name

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CAP_PT = 7.0        # caption point size AT PRINT SIZE (k = 1)
CAP_LEAD = 1.32     # line spacing multiple

FRAMING_EQDYNA = ('CROSS-VERSION CONSISTENCY CHECK between EQdyna results -- '
                  'NOT an independent validation, NOT an accuracy claim: '
                  'agreement certifies no run against the benchmark.')
FRAMING_CROSS_CODE = ('CROSS-CODE COMPARISON against independent SCEC CVWS '
                      'submissions -- the benchmark has no analytic answer; '
                      'agreement means agreement between codes, not accuracy.')


def framing(models):
    if any(m.is_independent for m in models):
        return FRAMING_CROSS_CODE
    return FRAMING_EQDYNA


# ================================================================== numbers

def as_grid(strike, downdip, t):
    """Reshape a node list to a lattice WITHOUT interpolating. Returns
    (X, Y, T) if the nodes form a complete rectangular lattice, else None so
    the caller falls back to a Delaunay contour on the nodes as given."""
    ux = np.unique(np.round(strike, 3))
    uy = np.unique(np.round(downdip, 3))
    if ux.size * uy.size != strike.size:
        return None
    T = np.full((uy.size, ux.size), np.nan)
    ix = np.searchsorted(ux, np.round(strike, 3))
    iy = np.searchsorted(uy, np.round(downdip, 3))
    T[iy, ix] = t
    if np.isnan(T).all():
        return None
    return ux, uy, T


def nice_step(tmax, target=9):
    raw = tmax / max(target, 1)
    if raw <= 0 or not np.isfinite(raw):
        return 1.0
    mag = 10.0 ** np.floor(np.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if raw <= m * mag:
            return m * mag
    return 10 * mag


def exact_match_stats(a, b):
    """Compare two rupture-time fields on the EXACT intersection of their node
    coordinates (keys rounded to 0.1 m; no snapping, no nearest-neighbour)."""
    ka = {(round(float(x), 1), round(float(d), 1)): i
          for i, (x, d) in enumerate(zip(a['strike'], a['downdip']))}
    kb = {(round(float(x), 1), round(float(d), 1)): i
          for i, (x, d) in enumerate(zip(b['strike'], b['downdip']))}
    common = ka.keys() & kb.keys()
    if not common:
        return dict(n_common=0, n_both=0, n_only_a=0, n_only_b=0,
                    median_dt=None, max_dt=None, overlap=None)
    ia = np.fromiter((ka[k] for k in common), int, len(common))
    ib = np.fromiter((kb[k] for k in common), int, len(common))
    ta, tb = a['t'][ia], b['t'][ib]
    ra, rb = np.isfinite(ta), np.isfinite(tb)
    both = ra & rb
    dt = np.abs(ta[both] - tb[both])
    denom = int(both.sum() + (rb & ~ra).sum())
    return dict(n_common=len(common), n_both=int(both.sum()),
                n_only_a=int((ra & ~rb).sum()), n_only_b=int((rb & ~ra).sum()),
                median_dt=float(np.median(dt)) if dt.size else None,
                max_dt=float(np.max(dt)) if dt.size else None,
                overlap=(float(both.sum()) / denom) if denom else None)


def ruptured_area_km2(strike, downdip, t, dx):
    """(ruptured node count) x (node cell area), each model at its OWN
    spacing -- the one definition that means the same at 25 m and 500 m."""
    n = int(np.isfinite(t).sum())
    return n * dx * dx / 1.0e6, n


def seed_count(T):
    """Ruptured nodes strictly earlier than every ruptured 8-neighbour: fronts
    that started on their own rather than arriving (retired TPV29/30 one-off;
    8-connectivity reproduces its 500 m numbers)."""
    R = np.isfinite(T)
    Tp = np.pad(np.where(R, T, np.inf), 1, constant_values=np.inf)
    nz, nx = T.shape
    nb = np.full(T.shape, np.inf)
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dj or di:
                nb = np.minimum(nb, Tp[1 + dj:1 + dj + nz, 1 + di:1 + di + nx])
    return int((R & np.isfinite(nb) & (T < nb)).sum())


def tau_over_strength(case_dir, dx):
    """(max, count >= 1, n) of t=0 tau/strength via the REGRESSION module's
    own function, so the figure and the gate can never disagree."""
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


# =============================================================== rupture

def _field_rows(models):
    fields, rows = [], []
    for m in models:
        f = m.load_rupture()
        fields.append(f)
        area, nrup = ruptured_area_km2(f['strike'], f['downdip'], f['t'], f['dx_m'])
        g = as_grid(f['strike'], f['downdip'], f['t'])
        rows.append(dict(label=m.label(f['dx_m']), short=m.label(), n=f['t'].size,
                         nrup=nrup, area=area, dx=f['dx_m'], grid=g,
                         tmax=float(np.nanmax(f['t'])) if nrup else float('nan'),
                         sent=f['sentinel_how'], rule=f['dipwise'], src=f['source'],
                         seeds=seed_count(g[2]) if g is not None else None))
    return fields, rows


def _levels(rows, args):
    tmax = max((r['tmax'] for r in rows if np.isfinite(r['tmax'])), default=1.0)
    step = args.contour_step or nice_step(tmax)
    return step, np.arange(step, tmax + step, step)    # SHARED across all models


def _contour(ax, m, f, g, levels, label_every=None):
    st = m.style()
    kw = dict(levels=levels, colors=[st['color']], linestyles=[st['linestyle']],
              linewidths=st['linewidth'])
    if g is None:
        cs = ax.tricontour(f['strike'] / 1e3, f['downdip'] / 1e3, f['t'], **kw)
    else:
        X, Y, T = g
        cs = ax.contour(X / 1e3, Y / 1e3, T, **kw)
    if label_every is not None:
        ax.clabel(cs, levels[::label_every], fmt='%g', fontsize=6.0, inline=True)


def _decorate_fault_axes(ax, xs, ds, ax_h, args, ylabel=True):
    ax.set_xlabel('Along-strike distance (km)', fontsize=9, labelpad=1.5)
    if ylabel:
        ax.set_ylabel('Down-dip distance (km)', fontsize=9, labelpad=1.5)
    ax.tick_params(labelsize=7.5, pad=1.5)
    ax.set_xlim(xs.min() / 1e3, xs.max() / 1e3)
    ax.set_ylim(ds.max() / 1e3, ds.min() / 1e3)      # down-dip positive DOWN
    ax.locator_params(axis='x', nbins=9)
    # tick density from the PRINTED height of the axes, not a fixed count
    ax.locator_params(axis='y', nbins=max(4, int(round(ax_h * 2.5))))
    ax.grid(alpha=0.18, linewidth=0.4)
    if args.hypo:
        ax.plot(args.hypo[0], args.hypo[1], 'k*', ms=8, zorder=5)


def _rupture_captions(models, rows, pairs, step, args):
    cap = [framing(models),
           f'Rupture-time contours, {step:g} s interval, labelled every '
           f'{2 * step:g} s on {rows[0]["label"]}. Each model is contoured on '
           f'its OWN node grid; nothing is interpolated or resampled. Node '
           f'spacings: ' + ', '.join(f'{r["dx"]:g} m' for r in rows) + '.',
           'Ruptured nodes / ruptured area / latest arrival, in legend order: '
           + '; '.join(f'{r["nrup"]}/{r["n"]}, {r["area"]:.0f} km2, '
                       f'{r["tmax"]:.2f} s' for r in rows)
           + '. Never-ruptured fill detected per file, legend order: '
           + ', '.join(r['sent'].split(' ')[0] for r in rows) + ' s.']
    vs = []
    for m, st in pairs:
        if st['n_common'] == 0:
            vs.append(f'{m.label()}: NO node coordinate shared exactly, no '
                      f'pointwise statistic reported (interpolating would '
                      f'invent one)')
        else:
            vs.append(f'{m.label()}: {st["n_common"]} nodes, '
                      f'{st["median_dt"]:.3f} / {st["max_dt"]:.2f} s, '
                      f'{100 * st["overlap"]:.1f}%')
    cap.append(f'Against {rows[0]["label"]} on EXACTLY-shared node '
               f'coordinates only -- median |dt| / max |dt| / rupture-extent '
               f'overlap: ' + '; '.join(vs) + '.')
    if args.seeds:
        cap.append('Seeds (a ruptured node strictly earlier than every ruptured '
                   '8-neighbour), legend order: ' + ', '.join(
                       str(r['seeds']) if r['seeds'] is not None else 'n/a (no lattice)'
                       for r in rows) + '.')
    if args.tau_strength:
        mx, nat, ntot = tau_over_strength(args.tau_strength, rows[0]['dx'])
        cap.append(f't = 0 max tau/strength {mx:.5f}; {nat} of {ntot} nodes at '
                   f'or above strength (testsys/regression/'
                   f'test_rough_fault_normal_consistency.py, {rows[0]["dx"]:g} m).')
    for m in models:
        for n in m.notes:
            cap.append(f'{m.label()}: {n}.')
    if args.note:
        cap.append(args.note)
    return cap


def _rupture_report(rows, pairs):
    report = ['=== per-model summary (by-product; the figure is the result) ===']
    for r in rows:
        report.append(f'  {r["label"]}')
        report.append(f'      source={r["src"]}  {r["nrup"]}/{r["n"]} nodes '
                      f'ruptured  area={r["area"]:.1f} km2  latest={r["tmax"]:.3f} s')
        report.append(f'      never-ruptured fill: {r["sent"]}   {r["rule"]}')
    report.append('=== pairwise, EXACT coordinate matches only (no snapping, '
                  'no interpolation) ===')
    for m, st in pairs:
        report.append(f'  {m.label()} vs {rows[0]["short"]}: '
                      f'common={st["n_common"]} both_ruptured={st["n_both"]} '
                      f'only_ref={st["n_only_a"]} only_other={st["n_only_b"]} '
                      + (f'median|dt|={st["median_dt"]:.3f}s '
                         f'max|dt|={st["max_dt"]:.3f}s '
                         f'overlap={100 * st["overlap"]:.1f}%'
                         if st['n_common'] else ''))
    return report


def _legend(target, models, rows, **kw):
    handles = [plt.Line2D([], [], **m.style(), label=r['label'])
               for m, r in zip(models, rows)]
    target.legend(handles=handles, fontsize=7.5, ncol=min(len(models), 2),
                  frameon=False, handlelength=3.0, columnspacing=1.4,
                  borderpad=0.2, labelspacing=0.3, **kw)


def plot_cplot(models, args):
    fields, rows = _field_rows(models)
    step, levels = _levels(rows, args)
    pairs = [(m, exact_match_stats(fields[0], f))
             for m, f in zip(models[1:], fields[1:])]
    text, cap_in = caption_block(_rupture_captions(models, rows, pairs, step, args), args)

    # layout in INCHES: the axes get true 1:1 geometry
    xs = np.concatenate([f['strike'] for f in fields])
    ds = np.concatenate([f['downdip'] for f in fields])
    span_x = (xs.max() - xs.min()) / 1e3
    span_d = (ds.max() - ds.min()) / 1e3
    W = args.width_in
    L, Rm = 0.62, 0.10
    ax_w = W - L - Rm
    ax_h = ax_w * (span_d / max(span_x, 1e-9))
    top_in = 0.30 + 0.20 * int(np.ceil(len(models) / 2))
    bot_in = 0.42 + cap_in
    H = top_in + ax_h + bot_in
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L / W, bot_in / H, ax_w / W, ax_h / H])
    for m, f, r in zip(models, fields, rows):
        _contour(ax, m, f, r['grid'], levels, 2 if m.index == 0 else None)
    _decorate_fault_axes(ax, xs, ds, ax_h, args)
    _legend(ax, models, rows, loc='lower center', bbox_to_anchor=(0.5, 1.005))
    finish(fig, text, _rupture_report(rows, pairs), args)
    return rows


def plot_res_series(models, args):
    """Model 0 is the reference; every other model gets its own panel with
    the reference under it, all on the SAME contour levels."""
    if len(models) < 2:
        raise SystemExit('--plot res-series needs a reference and >= 1 other model')
    fields, rows = _field_rows(models)
    step, levels = _levels(rows, args)
    pairs = [(m, exact_match_stats(fields[0], f))
             for m, f in zip(models[1:], fields[1:])]
    text, cap_in = caption_block(_rupture_captions(models, rows, pairs, step, args), args)

    xs = np.concatenate([f['strike'] for f in fields])
    ds = np.concatenate([f['downdip'] for f in fields])
    span_x = (xs.max() - xs.min()) / 1e3
    span_d = (ds.max() - ds.min()) / 1e3
    npan = len(models) - 1
    W = args.width_in
    L, Rm, VS = 0.62, 0.10, 0.62
    ax_w = W - L - Rm
    ax_h = ax_w * (span_d / max(span_x, 1e-9))
    top_in = 0.30 + 0.20 * int(np.ceil(len(models) / 2))
    bot_in = cap_in + 0.10
    H = top_in + npan * ax_h + npan * VS + bot_in
    fig = plt.figure(figsize=(W, H))
    for k, (m, f, r, (_, st)) in enumerate(zip(models[1:], fields[1:], rows[1:], pairs)):
        y0 = bot_in + (npan - 1 - k) * (ax_h + VS) + 0.42
        ax = fig.add_axes([L / W, y0 / H, ax_w / W, ax_h / H])
        _contour(ax, models[0], fields[0], rows[0]['grid'], levels, 2)
        _contour(ax, m, f, r['grid'], levels)
        _decorate_fault_axes(ax, xs, ds, ax_h, args)
        ax.set_title(r['short'] + (f'  [{st["n_common"]} shared nodes, median '
                                   f'|dt| {st["median_dt"]:.3f} s]'
                                   if st['n_common'] else '  [no shared node]'),
                     fontsize=8.5, pad=2)
    _legend(fig, models, rows, loc='upper center',
            bbox_to_anchor=(0.5, 1.0 - 0.04 / H))
    finish(fig, text, _rupture_report(rows, pairs), args)
    return rows


# =========================================================== time series

def station_table(models, which, args):
    """Intersect the discovered station sets and keep ONLY stations that land
    exactly on every model's node grid. Excluded ones are counted with the
    reason -- never snapped, never interpolated."""
    per = [m.station_files(which) for m in models]
    spacing = []
    for m in models:
        try:
            spacing.append(m.load_rupture()['dx_m'])
        except SystemExit:
            spacing.append(float('nan'))
    union = set().union(*per) if per else set()
    keep, excl = [], []
    for key in sorted(union):
        missing = [models[i].label() for i, d in enumerate(per) if key not in d]
        if missing:
            excl.append((key, f'absent from {len(missing)} of {len(models)} models'))
            continue
        offgrid = []
        for i, dx in enumerate(spacing):
            if not np.isfinite(dx):
                continue
            if any(abs(c) % dx > 1e-6 and abs(abs(c) % dx - dx) > 1e-6 for c in key):
                offgrid.append(f'{models[i].label()} (dx {dx:g} m)')
        if offgrid:
            excl.append((key, 'coordinates not on the node grid of ' + ', '.join(offgrid)))
            continue
        keep.append(key)
    return per, keep, excl, spacing


def _ts_columns(which, comps):
    """(component, quantity) per column, component-major; a quantity with no
    column for that component (h/v body has no n... fault n is shear only)
    is skipped."""
    quants = FAULT_COLS if which == 'fault' else BODY_COLS
    order = ['slip', 'sliprate', 'shear'] if which == 'fault' else ['disp', 'vel']
    return [(c, q) for c in comps for q in order if quants[q].get(c) is not None]


def _col_title(which, comp, q):
    spec = (FAULT_COLS if which == 'fault' else BODY_COLS)[q]
    if which == 'fault' and comp == 'n':
        return f'n-normal stress ({spec["unit"]})'
    return f'{comp}-{spec["name"]} ({spec["unit"]})'


def plot_ts(models, which, args):
    from scec_readers import load_station, parse_station_request
    per, keep, excl, spacing = station_table(models, which, args)
    if not keep:
        raise SystemExit(
            f'No {which} station is common to all models AND on every '
            f'model\'s node grid. Discovered per model: '
            + '; '.join(f'{m.label()}={len(d)}' for m, d in zip(models, per))
            + f'. {len(excl)} candidate(s) excluded -- first reasons: '
            + '; '.join(f'{station_name(k, which)}: {r}' for k, r in excl[:3]))

    if args.stations:
        shown = []
        for s in args.stations:
            k = parse_station_request(s, which)
            if k not in keep:
                why = dict(excl).get(k, 'not discovered in any model')
                raise SystemExit(f'--stations {s}: not comparable -- {why}')
            shown.append(k)
    else:
        # deterministic, spread evenly through the sorted station list
        n_show = min(args.max_stations, len(keep))
        idx = np.unique(np.linspace(0, len(keep) - 1, n_show).round().astype(int))
        shown = [keep[i] for i in idx]

    quants = FAULT_COLS if which == 'fault' else BODY_COLS
    comps = args.component
    cols = _ts_columns(which, comps)
    if not cols:
        raise SystemExit(f'--component {" ".join(comps)}: no {which}-station column')
    ncol, nrow = len(cols), len(shown)

    # layout in INCHES (fix the layout, never shrink the font)
    W = args.width_in
    L, Rm, PANEL_H, HS_IN, WS_IN = 0.62, 0.10, 1.05, 0.22, 0.52
    legend_rows = int(np.ceil(len(models) / 2))
    top_in = 0.26 + 0.20 * legend_rows + 0.22
    axes_h = nrow * PANEL_H + (nrow - 1) * HS_IN
    panel_w = (W - L - Rm - (ncol - 1) * WS_IN) / ncol
    fig = plt.figure(figsize=(W, top_in + axes_h + 0.40))
    axes = fig.subplots(nrow, ncol, squeeze=False, sharex=True)

    notes, tmax_all, peaks = [], 0.0, {cq: [] for cq in cols}
    for r, key in enumerate(shown):
        for m in models:
            a, note = load_station(per[m.index][key], which)
            if note:
                notes.append(note)
            tmax_all = max(tmax_all, float(a[:, 0].max()))
            for c, (comp, q) in enumerate(cols):
                col = quants[q][comp]
                if col >= a.shape[1]:
                    continue
                axes[r][c].plot(a[:, 0], a[:, col], **m.style())
                peaks[(comp, q)].append(float(np.max(np.abs(a[:, col]))))
        axes[r][0].set_ylabel(station_name(key, which), fontsize=8, labelpad=1.5)
    for c, (comp, q) in enumerate(cols):
        axes[0][c].set_title(_col_title(which, comp, q), fontsize=9.0, pad=3)
        lo = min(axes[r][c].get_ylim()[0] for r in range(nrow))
        hi = max(axes[r][c].get_ylim()[1] for r in range(nrow))
        for r in range(nrow):                      # shared y per COLUMN
            axes[r][c].set_ylim(lo, hi)
            axes[r][c].tick_params(labelsize=7.5, pad=1.5)
            axes[r][c].grid(alpha=0.18, linewidth=0.4)
            axes[r][c].locator_params(axis='y', nbins=4)
    for a in axes[-1]:
        a.set_xlim(0, tmax_all)
        a.locator_params(axis='x', nbins=6)

    handles = [plt.Line2D([], [], **m.style(), label=m.label(spacing[m.index]))
               for m in models]
    fig.legend(handles=handles, fontsize=7.5, loc='upper center',
               ncol=min(len(models), 2), frameon=False, handlelength=3.0,
               columnspacing=1.4, borderpad=0.2, labelspacing=0.3)

    kind = 'on-fault' if which == 'fault' else 'off-fault (body)'
    reasons = {}
    for k, r in excl:
        reasons[r] = reasons.get(r, 0) + 1
    comp_s = '/'.join(comps)
    cap = [framing(models),
           f'{kind} stations, {comp_s}-component; rows are stations, columns '
           f'quantities, y-range shared down each column. Node spacings: '
           + ', '.join(f'{m.label()}' for m in models) + '.',
           f'{len(shown)} of {len(keep)} comparable stations shown ('
           + ('named by --stations' if args.stations else
              'evenly sampled from the sorted list')
           + f'); {len(excl)} of {len(keep) + len(excl)} discovered stations '
           f'excluded because their coordinates are not shared exactly -- '
           f'none snapped, none interpolated.',
           (f'Peak |{comps[0]}| over the panels shown: ' + '; '.join(
               f'{quants[q]["name"]} {max(peaks[(c, q)]):.4g} {quants[q]["unit"]}'
               for c, q in cols if peaks[(c, q)]) if len(comps) == 1 else
            'Peak |value| over the panels shown: ' + '; '.join(
               f'{_col_title(which, c, q).split(" (")[0]} {max(peaks[(c, q)]):.4g} '
               f'{quants[q]["unit"]}' for c, q in cols if peaks[(c, q)]))
           + f'; traces span 0-{tmax_all:.2f} s.']
    if notes:
        cap.append('Header/data column-count mismatch (data trusted): '
                   + '; '.join(sorted(set(notes))[:2]) + '.')
    for m in models:
        for n in m.notes:
            cap.append(f'{m.label()}: {n}.')
    if args.note:
        cap.append(args.note)
    # now that the caption is known, give the figure exactly the height it
    # needs and re-place the axes -- no overlap, no font below 7 pt
    text, cap_in = caption_block(cap, args)
    bot_in = cap_in + 0.40
    H = top_in + axes_h + bot_in
    fig.set_size_inches(W, H)
    fig.subplots_adjust(left=L / W, right=1.0 - Rm / W,
                        top=1.0 - top_in / H, bottom=bot_in / H,
                        hspace=HS_IN / PANEL_H, wspace=WS_IN / panel_w)
    fig.legends[0].set_bbox_to_anchor((0.5, 1.0 - 0.04 / H),
                                      transform=fig.transFigure)
    fig.text(L / W + (1.0 - (L + Rm) / W) / 2.0, (cap_in + 0.03) / H,
             'Time (s)', fontsize=9, ha='center', va='bottom')

    report = [f'=== {kind} stations: {len(keep)} comparable, {len(excl)} '
              f'excluded (by-product; the figure is the result) ===']
    for m, d in zip(models, per):
        report.append(f'  discovered in {m.label()}: {len(d)}')
    report.append('  plotted: ' + ', '.join(station_name(k, which) for k in shown))
    for r, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        report.append(f'  excluded x{n}: {r}')
    for n in sorted(set(notes)):
        report.append('  format note: ' + n)
    for m in models:
        for n in m.notes:
            report.append(f'  {m.label()}: {n}')
    finish(fig, text, report, args)
    return shown, keep, excl


# ================================================================ output

def caption_block(lines, args):
    """Wrap the computed caption to the PRINT width; return it with the
    vertical space it needs (inches) so the caller reserves that space BEFORE
    placing axes: the caption never overlaps, the font never drops below 7 pt."""
    text = '\n'.join(_wrap(l, args.wrap_chars) for l in lines)
    n = text.count('\n') + 1
    return text, (n * CAP_PT * CAP_LEAD + 8.0) / 72.0


def finish(fig, text, report, args):
    fig.text(0.006, 0.006, text, fontsize=CAP_PT, va='bottom', ha='left',
             linespacing=CAP_LEAD)
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    for line in report:
        print(line)
    w_in, h_in = fig.get_size_inches()
    print(f'\nwrote {args.out}')
    print(f'  canvas {w_in:.2f} x {h_in:.2f} in at {args.dpi} dpi = '
          f'{w_in * args.dpi:.0f} x {h_in * args.dpi:.0f} px; print width '
          f'{args.print_width_mm:g} mm -> k = canvas/print = 1.000, so the '
          f'set font sizes ARE the printed point sizes '
          f'(labels 9, ticks 7.5, titles 9, legend 7.5, caption {CAP_PT:g}); '
          f'{w_in * args.dpi / (args.print_width_mm / 25.4):.0f} dpi at print width.')
    print('\n--- caption as placed on the figure (every number computed from '
          'the plotted arrays) ---')
    print(text)


def _wrap(s, n):
    out, line = [], ''
    for w in s.split():
        if len(line) + len(w) + 1 > n:
            out.append(line)
            line = w
        else:
            line = (line + ' ' + w).strip()
    out.append(line)
    return '\n    '.join(out)
