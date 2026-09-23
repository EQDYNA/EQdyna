#! /usr/bin/env python3
"""
scec_compare.py -- overlay N EQdyna results (our runs and/or SCEC archive
submissions) for ANY SCEC/USGS TPV benchmark, one plot mode per invocation.

    scec_compare.py --plot cplot    --models P1 P2 [P3 ...]   rupture-time contours
    scec_compare.py --plot ts-fault --models P1 P2 [P3 ...]   on-fault time series
    scec_compare.py --plot ts-body  --models P1 P2 [P3 ...]   off-fault time series

Each PATH is either one of OUR run/reference directories or an archive
submission directory; the script sniffs which by layout, it is never declared.
Styles/labels are assigned deterministically from the ORDER of --models, so a
given model keeps the same colour, dash pattern and label across all three
plot modes and re-running with the same arguments reproduces the same figure.

WHAT THIS IS, AND IT GOES ON EVERY FIGURE
-----------------------------------------
Comparing an EQdyna run against one of the scec_archive submissions is a
CROSS-VERSION CONSISTENCY CHECK against one modeller's own earlier submission.
It is NOT an independent validation and NOT an accuracy claim: agreement
certifies neither run against the benchmark. The archived submissions were
themselves produced by EQdyna (v3.1 ... v5.3.3) by the same author.

EVERY NUMBER IN EVERY CAPTION IS COMPUTED FROM THE ARRAYS BEING PLOTTED, at
plot time. Nothing here is hand-typed. (A previous one-off shipped a
hand-written caption that was stale on arrival.)

NEVER MODIFIES scec_archive/ -- it is opened read-only, always.

FORMAT HAZARDS THIS HANDLES EXPLICITLY (each one has cost time before)
----------------------------------------------------------------------
* SENTINELS DIFFER BY WRITER AND BY ERA. Never-ruptured is 1000.0 s in the
  2015-2020 archives (tpv29/30/34/35/104/105-3d) but 99999.0 s in the 2024
  ones (tpv36/37) and in our own writer. A `> 1e4` mask lets 1000 through and
  it contours as a real rupture front. detect_sentinel() finds it per file
  from the data (a repeated maximum >= 10x the 99.9th percentile of the rest)
  and every figure states the value it used. --sentinel overrides.
* THE cplot HEADER ENDS IN AN UNPARSEABLE ` j  k  t` LINE after the `#`
  block, so np.loadtxt(..., comments='#') alone chokes. read_numeric() finds
  the first all-numeric line by scanning only the header and hands loadtxt a
  skiprows.
* COLUMN LAYOUT DIFFERS. Ours (frt.canonical.txt): col0 x(m), col1 y(m),
  col2 z(m, negative down), col3 rupture time. Archive cplot: along-strike(m),
  down-dip(m, POSITIVE DOWN), rupture time(s). Our SCECRuptureTime.txt is
  already in the archive's convention. All three are normalised to
  (strike_m, downdip_m positive-down, t_s) on read.
* DOWN-DIP IS NOT DEPTH ON A DIPPING FAULT. From frt.canonical.txt the
  down-dip distance is derived as hypot(y, z) when the node cloud shows the
  fault is dipping (|corr(y,z)| > 0.95) and as -z when it is vertical; rough
  vertical faults (tpv29/30) have y that is roughness noise, uncorrelated
  with z, and using hypot there would corrupt every coordinate by ~0.4%.
  The rule used is decided from the data and printed.
* NO INTERPOLATION ACROSS RESOLUTIONS, EVER. Each model is contoured on its
  OWN grid (rule 17 step 2: a rough fault surface is one fixed sampling of one
  random realisation). Quantitative comparison uses EXACT coordinate matches
  only -- the intersection of the node coordinate sets, no snapping, no
  nearest-neighbour. Stations likewise: only stations whose coordinates land
  exactly on every model's node grid are plotted, and the excluded ones are
  counted and explained.
* A STATION HEADER CAN LIE ABOUT ITS OWN COLUMN COUNT -- true of ARCHIVED
  files, which this script still reads. Before pathway item 67 (fixed
  2026-09-23), today's writer declared "Time series in 11 columns" and wrote
  8 whenever friclaw < 3 (src/fortran/library_output.f90:66 vs :97); the
  writer now declares the branch's true column count. The 2015 archive
  predates that fix and still declares 8 while writing 8 (so it never lied),
  but older EQdyna-written references on disk (e.g.
  scripts/fractal_stress_diamond_square/tpv104.200m.asp.ref/fort.51, an
  archived reference left as-is) still carry the old, wrong header. The DATA
  is trusted, and any mismatch is reported rather than silently coped with --
  that coping logic stays for exactly this reason.
* AN ARCHIVE DIRECTORY CAN CARRY A STALE TEMPLATE HEADER. tpv30's 100 m cplot
  still says `problem = LVFZ, date = 6-2-2012` while its own faultst files
  correctly say TPV30 / Dunyu Liu / 2015. cplot headers are therefore NEVER
  used as provenance; labels come from the path (self-describing in this
  archive: `eqdyna-v3.1-100m-2015`) and node spacing is MEASURED from the
  coordinates rather than read from any header.

PRINT GEOMETRY
--------------
Figures are laid out for Elsevier full width, 190 mm (--print-width-mm), at
k = canvas_width / print_width = 1.0, so every font size below is literally
the printed point size: axis labels 9 pt, ticks 7.5 pt, panel titles 9.5 pt,
legend 8 pt, caption 7 pt. 400 dpi -> 2992 px across at 190 mm (>= 300 dpi).
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- constants

# Deterministic, colour-blind-safe. Index into these by POSITION IN --models,
# never by dict iteration or mtime, so the same command reproduces the same
# figure and a model keeps its identity across cplot / ts-fault / ts-body.
COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7',
          '#56B4E9', '#E69F00', '#000000', '#F0E442']
LINESTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1)),
              (0, (1, 1)), (0, (4, 1, 1, 1, 1, 1))]

FRAMING = ('CROSS-VERSION CONSISTENCY CHECK between EQdyna results -- NOT an '
           'independent validation, NOT an accuracy claim: agreement '
           'certifies no run against the benchmark.')

# SCEC station-name convention: distances in units of 100 m.
STATION_UNIT_M = 100.0
RE_FAULTST = re.compile(r'^faultst(-?\d+)dp(-?\d+)(?:\.txt)?$')
RE_BODYST = re.compile(r'^body(-?\d+)st(-?\d+)dp(-?\d+)(?:\.txt)?$')

# Column layout of the SCEC station time series. Identical in the 2015 archive
# and in today's writer (library_output.f90 output_onfault_st/output_offfault_st)
# for the first 8 / 7 columns; RSF runs append psi/T/p to the fault files.
FAULT_COLS = {
    'slip':       dict(h=1, v=4, unit='m',    name='slip'),
    'sliprate':   dict(h=2, v=5, unit='m/s',  name='slip rate'),
    'shear':      dict(h=3, v=6, unit='MPa',  name='shear stress'),
}
FAULT_DECLARED_MIN = 8
BODY_COLS = {
    'disp': dict(h=1, v=3, n=5, unit='m',   name='displacement'),
    'vel':  dict(h=2, v=4, n=6, unit='m/s', name='velocity'),
}
BODY_DECLARED_MIN = 7

CAP_PT = 7.0        # caption point size AT PRINT SIZE (k = 1)
CAP_LEAD = 1.32     # line spacing multiple

RE_NCOL_DECL = re.compile(r'Time series in\s+(\d+)\s+columns', re.I)
RE_ELEMSIZE = re.compile(r'element_size\s*=\s*([0-9.eEdD+-]+)')


# ------------------------------------------------------------------- I/O


def read_numeric(path, maxcols=None):
    """Load the numeric block of a SCEC-style text file.

    Both the archive and our writer put an unparseable field-name line AFTER
    the '#' comment block (` j  k  t` in cplot/SCECRuptureTime.txt,
    ` t h-slip h-slip-rate ...` in the station files), which defeats
    np.loadtxt(comments='#'). Scan only the header for the first line whose
    tokens all parse as floats, then let loadtxt do the bulk in C -- the 25 m
    tpv30 cplot is 1.28 M rows and a pure-python parse of it is not free.
    """
    skip = 0
    with open(path, 'r') as fh:
        for i, line in enumerate(fh):
            if i > 400:
                raise SystemExit(f'{path}: no numeric data in first 400 lines')
            s = line.strip()
            if not s or s.startswith('#'):
                skip = i + 1
                continue
            try:
                [float(t) for t in s.split()]
            except ValueError:
                skip = i + 1
                continue
            skip = i
            break
    a = np.loadtxt(path, comments='#', skiprows=skip, ndmin=2)
    if maxcols is not None and a.shape[1] > maxcols:
        a = a[:, :maxcols]
    return a


def header_text(path, nlines=40):
    out = []
    with open(path, 'r') as fh:
        for i, line in enumerate(fh):
            if i >= nlines:
                break
            out.append(line.rstrip('\n'))
    return out


def header_declared_ncol(path):
    for line in header_text(path):
        m = RE_NCOL_DECL.search(line)
        if m:
            return int(m.group(1))
    return None


def header_element_size_m(path):
    """Parse element_size from a header. Returned for REPORTING ONLY -- every
    spacing used in a figure or a metric is measured from the coordinates
    instead, because archive headers are demonstrably stale (see module
    docstring, tpv30's `problem = LVFZ` cplot)."""
    for line in header_text(path):
        m = RE_ELEMSIZE.search(line)
        if m:
            try:
                return float(m.group(1).replace('D', 'E').replace('d', 'e'))
            except ValueError:
                return None
    return None


def detect_sentinel(t, override=None):
    """Find the 'never ruptured' fill value from the data itself.

    1000.0 s in the 2015-2020 archives, 99999.0 s in the 2024 ones and in our
    own writer. Criterion: the maximum is repeated AND is at least 10x the
    99.9th percentile of everything below it -- which is true of both
    conventions by orders of magnitude and false for a genuine latest arrival.
    Returns (sentinel_or_None, how_string).
    """
    if override is not None:
        return float(override), f'{override:g} s (--sentinel)'
    t = np.asarray(t, float)
    finite = t[np.isfinite(t)]
    if finite.size == 0:
        return None, 'no finite values'
    vmax = float(finite.max())
    n_at = int((finite == vmax).sum())
    rest = finite[finite < vmax]
    if n_at >= 2 and (rest.size == 0 or vmax >= 10.0 * np.percentile(rest, 99.9)):
        return vmax, f'{vmax:g} s (detected, {n_at} nodes)'
    return None, 'none detected (every node ruptured)'


# --------------------------------------------------------------- the model


class Model:
    """One result set: either one of OUR run/reference directories or one
    archive submission directory. Which it is, is sniffed, never declared."""

    def __init__(self, path, index, label=None, sentinel=None):
        self.path = os.path.abspath(path)
        self.index = index
        self.sentinel_override = sentinel
        self.notes = []          # format problems worth reporting, not hiding
        self._dx_cache = None    # MEASURED spacing, filled by load_rupture()
        self.color = COLORS[index % len(COLORS)]
        self.linestyle = LINESTYLES[index % len(LINESTYLES)]
        self.linewidth = float(np.clip(1.6 - 0.18 * index, 0.7, 1.6))
        self._sniff()
        self.user_label = label

    # -- discovery -------------------------------------------------------

    def _sniff(self):
        p = self.path
        if os.path.isfile(p):
            self.dir = os.path.dirname(p)
            self.rupture_file = p
            self.kind = 'eqdyna' if 'frt' in os.path.basename(p) else 'archive'
        elif os.path.isdir(p):
            self.dir = p
            cands = [('archive', os.path.join(p, 'cplot')),
                     ('eqdyna', os.path.join(p, 'frt.canonical.txt')),
                     ('eqdyna', os.path.join(p, 'SCECRuptureTime.txt'))]
            self.kind, self.rupture_file = None, None
            for kind, f in cands:
                if os.path.isfile(f):
                    self.kind, self.rupture_file = kind, f
                    break
            if self.kind is None:
                # our raw per-rank output, not yet canonicalised
                if glob.glob(os.path.join(p, 'frt.txt*')):
                    self.kind, self.rupture_file = 'eqdyna', 'frt.txt*'
                elif glob.glob(os.path.join(p, 'faultst*')):
                    self.kind = 'archive' if not glob.glob(
                        os.path.join(p, 'faultst*.txt')) else 'eqdyna'
                else:
                    raise SystemExit(
                        f'{p}: not recognisable as a model directory -- no '
                        f'cplot, frt.canonical.txt, SCECRuptureTime.txt, '
                        f'frt.txt* or faultst* found.')
        else:
            raise SystemExit(f'{p}: no such file or directory')

    def station_files(self, which):
        """Discover stations by PARSING names, never by enumerating a list.
        Archive files have no extension, ours end in .txt; both are globbed.
        Names the writer mangled (`faultst***dp075.txt`, the pre-item-22
        i3.3 overflow on negative strike) are counted and reported, not
        silently dropped."""
        rx = RE_FAULTST if which == 'fault' else RE_BODYST
        pat = 'faultst*' if which == 'fault' else 'body*'
        found, bad = {}, []
        for f in sorted(glob.glob(os.path.join(self.dir, pat))):
            base = os.path.basename(f)
            if base.endswith('.zip') or os.path.isdir(f):
                continue
            m = rx.match(base)
            if not m:
                bad.append(base)
                continue
            g = [int(v) * STATION_UNIT_M for v in m.groups()]
            key = tuple(g)            # (strike, dip) or (offset, strike, dip)
            found[key] = f
        if bad:
            self.notes.append(
                f'{len(bad)} {which}-station file(s) with unparseable names '
                f'skipped (e.g. {bad[0]}) -- pre-fix i3.3 strike overflow')
        return found

    # -- rupture-time field ---------------------------------------------

    def load_rupture(self):
        """-> dict(strike, downdip, t, sentinel_how, dx_m, grid) with down-dip
        POSITIVE DOWN in metres, for every input flavour."""
        if self.rupture_file == 'frt.txt*':
            files = sorted(glob.glob(os.path.join(self.dir, 'frt.txt*')))
            a = np.vstack([read_numeric(f) for f in files])
            src = f'{len(files)} x frt.txt<rank>'
            xs, dd, t, rule = self._frt_to_strike_dip(a)
        elif os.path.basename(self.rupture_file).startswith('frt'):
            a = read_numeric(self.rupture_file)
            if a.shape[1] < 4:
                raise SystemExit(f'{self.rupture_file}: {a.shape[1]} columns, '
                                 f'need >= 4 (x, y, z, rupture time)')
            src = os.path.basename(self.rupture_file)
            xs, dd, t, rule = self._frt_to_strike_dip(a)
        else:
            a = read_numeric(self.rupture_file, maxcols=3)
            if a.shape[1] != 3:
                raise SystemExit(f'{self.rupture_file}: {a.shape[1]} columns, '
                                 f'expected 3 (strike, down-dip, t)')
            src = os.path.basename(self.rupture_file)
            xs, dd, t = a[:, 0], a[:, 1], a[:, 2]
            rule = 'file is already (strike, down-dip positive-down, t)'
        sent, how = detect_sentinel(t, self.sentinel_override)
        tm = t.astype(float).copy()
        if sent is not None:
            tm[t >= sent] = np.nan
        self._dx_cache = self._measure_dx(xs, dd)
        return dict(strike=xs, downdip=dd, t=tm, t_raw=t, sentinel=sent,
                    sentinel_how=how, dipwise=rule, source=src,
                    dx_m=self._dx_cache)

    def _frt_to_strike_dip(self, a):
        """frt.canonical.txt / frt.txt<rank>: col0 x, col1 y, col2 z (negative
        down), col3 rupture time. Down-dip distance is -z on a vertical fault
        and hypot(y, z) on a dipping one -- decided FROM THE DATA, because a
        rough vertical fault (tpv29/30) has y = roughness, uncorrelated with
        z, and hypot would silently bias every coordinate off the node grid."""
        x, y, z, t = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
        if np.std(y) > 1.0 and np.std(z) > 1.0:
            r = float(np.corrcoef(y, z)[0, 1])
        else:
            r = 0.0
        if abs(r) > 0.95:
            return x, np.hypot(y, z), t, f'dipping fault (corr(y,z)={r:+.3f}): down-dip = hypot(y,z)'
        return x, -z, t, f'vertical fault (corr(y,z)={r:+.3f}): down-dip = -z'

    @staticmethod
    def _measure_dx(xs, dd):
        """Node spacing MEASURED from the coordinates (headers are stale)."""
        out = []
        for v in (xs, dd):
            u = np.unique(np.round(np.asarray(v, float), 3))
            if u.size > 1:
                out.append(float(np.median(np.diff(u))))
        return float(np.median(out)) if out else float('nan')

    # -- labelling -------------------------------------------------------

    def label(self, dx_m=None):
        if self.user_label:
            base = self.user_label
        else:
            base = os.path.basename(self.dir) or self.dir
            if self.kind == 'archive':
                # Archive dir names are self-describing: eqdyna-v3.1-100m-2015
                m = re.match(r'eqdyna3?d?-(v[\d.]+)-(\d+)m-(\d{4})$', base)
                if m:
                    bench = os.path.basename(os.path.dirname(self.dir))
                    # the dir's OWN declared resolution goes in the label, so
                    # three submissions of one benchmark never share a legend
                    # entry; the MEASURED spacing is appended beside it below
                    # and a disagreement between the two is then visible.
                    base = (f'{bench} SCEC {m.group(1)} '
                            f'{m.group(2)} m ({m.group(3)})')
        if dx_m is None:
            dx_m = self._dx_cache
        if dx_m is not None and np.isfinite(dx_m):
            base = f'{base}, dx {dx_m:g} m'
        return base

    def style(self):
        return dict(color=self.color, linestyle=self.linestyle,
                    linewidth=self.linewidth)


# ---------------------------------------------------------------- helpers


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
    coordinates. No snapping, no nearest-neighbour, no interpolation: keys are
    coordinates rounded to 0.1 m and a node that is not shared is simply not
    compared. Returns the by-product summary numbers."""
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
    """Ruptured area = (ruptured node count) x (node cell area), each model at
    its OWN spacing. Deliberately not a corner rule: this is a by-product
    number, and a cell count is the one definition that means the same thing
    at 25 m and at 500 m."""
    n = int(np.isfinite(t).sum())
    return n * dx * dx / 1.0e6, n


# ------------------------------------------------------------- plot: cplot


def plot_cplot(models, args):
    fields, rows = [], []
    for m in models:
        f = m.load_rupture()
        fields.append(f)
        area, nrup = ruptured_area_km2(f['strike'], f['downdip'], f['t'], f['dx_m'])
        rows.append(dict(label=m.label(f['dx_m']), short=m.label(), n=f['t'].size,
                         nrup=nrup, area=area, dx=f['dx_m'],
                         tmax=float(np.nanmax(f['t'])) if nrup else float('nan'),
                         sent=f['sentinel_how'], rule=f['dipwise'], src=f['source']))

    tmax = max((r['tmax'] for r in rows if np.isfinite(r['tmax'])), default=1.0)
    step = nice_step(tmax)
    levels = np.arange(step, tmax + step, step)      # SHARED across all models

    pairs = []
    for m, f in zip(models[1:], fields[1:]):
        pairs.append((m, exact_match_stats(fields[0], f)))

    # ---- captions: short on the figure, full table to stdout -------------
    cap = [FRAMING,
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
    text, cap_in = caption_block(cap, args)

    # ---- layout in INCHES: the axes get true 1:1 geometry ----------------
    xs = np.concatenate([f['strike'] for f in fields])
    ds = np.concatenate([f['downdip'] for f in fields])
    span_x = (xs.max() - xs.min()) / 1e3
    span_d = (ds.max() - ds.min()) / 1e3
    W = args.width_in
    L, Rm = 0.62, 0.10                       # y-label gutter, right margin
    ax_w = W - L - Rm
    ax_h = ax_w * (span_d / max(span_x, 1e-9))
    legend_rows = int(np.ceil(len(models) / 2))
    top_in = 0.30 + 0.20 * legend_rows       # title + legend block
    bot_in = 0.42 + cap_in                   # x-label + caption
    H = top_in + ax_h + bot_in
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([L / W, bot_in / H, ax_w / W, ax_h / H])

    for m, f in zip(models, fields):
        st = m.style()
        g = as_grid(f['strike'], f['downdip'], f['t'])
        if g is None:
            cs = ax.tricontour(f['strike'] / 1e3, f['downdip'] / 1e3, f['t'],
                               levels=levels, colors=[st['color']],
                               linestyles=[st['linestyle']],
                               linewidths=st['linewidth'])
        else:
            X, Y, T = g
            cs = ax.contour(X / 1e3, Y / 1e3, T, levels=levels,
                            colors=[st['color']], linestyles=[st['linestyle']],
                            linewidths=st['linewidth'])
        if m.index == 0:
            ax.clabel(cs, levels[::2], fmt='%g', fontsize=6.0, inline=True)

    ax.set_xlabel('Along-strike distance (km)', fontsize=9, labelpad=1.5)
    ax.set_ylabel('Down-dip distance (km)', fontsize=9, labelpad=1.5)
    ax.tick_params(labelsize=7.5, pad=1.5)
    ax.set_xlim(xs.min() / 1e3, xs.max() / 1e3)
    ax.set_ylim(ds.max() / 1e3, ds.min() / 1e3)      # down-dip positive DOWN
    ax.locator_params(axis='x', nbins=9)
    # tick density from the PRINTED height of the axes, not a fixed count:
    # a 28 km down-dip fault (tpv36) gets a near-square panel and needs more
    # y ticks than a 2:1 one (tpv30)
    ax.locator_params(axis='y', nbins=max(4, int(round(ax_h * 2.5))))
    ax.grid(alpha=0.18, linewidth=0.4)

    handles = [plt.Line2D([], [], **m.style(), label=r['label'])
               for m, r in zip(models, rows)]
    ax.legend(handles=handles, fontsize=7.5, loc='lower center',
              bbox_to_anchor=(0.5, 1.005), ncol=min(len(models), 2),
              frameon=False, handlelength=3.0, columnspacing=1.4,
              borderpad=0.2, labelspacing=0.3)

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
    finish(fig, text, report, args)
    return rows


# ------------------------------------------------------- plot: time series


def station_table(models, which, args):
    """Intersect the discovered station sets and keep ONLY stations that land
    exactly on every model's node grid. Excluded ones are counted with the
    reason -- never snapped, never interpolated."""
    per = [m.station_files(which) for m in models]
    spacing = []
    for m in models:
        try:
            f = m.load_rupture()
            spacing.append(f['dx_m'])
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
            # station coordinates are exact multiples of dx or the station is
            # not a node of that model's mesh
            if any(abs(c) % dx > 1e-6 and abs(abs(c) % dx - dx) > 1e-6 for c in key):
                offgrid.append(f'{models[i].label()} (dx {dx:g} m)')
        if offgrid:
            excl.append((key, 'coordinates not on the node grid of ' + ', '.join(offgrid)))
            continue
        keep.append(key)
    return per, keep, excl, spacing


def _sccode(v):
    """SCEC filename convention: zero-padded magnitude, '-' only if negative
    (no '+'), so a reconstructed name is the real filename."""
    n = int(round(v / STATION_UNIT_M))
    return f'-{abs(n):03d}' if n < 0 else f'{n:03d}'


def station_name(key, which):
    if which == 'fault':
        s, d = key
        return f'faultst{_sccode(s)}dp{_sccode(d)}'
    b, s, d = key
    return f'body{_sccode(b)}st{_sccode(s)}dp{_sccode(d)}'


def load_station(path, which):
    a = read_numeric(path)
    decl = header_declared_ncol(path)
    need = FAULT_DECLARED_MIN if which == 'fault' else BODY_DECLARED_MIN
    note = None
    if decl is not None and decl != a.shape[1]:
        note = (f'{os.path.basename(path)} header declares {decl} columns, '
                f'file has {a.shape[1]} -- data trusted')
    if a.shape[1] < need:
        raise SystemExit(f'{path}: {a.shape[1]} columns, need >= {need}')
    return a, note


def plot_ts(models, which, args):
    per, keep, excl, spacing = station_table(models, which, args)
    if not keep:
        raise SystemExit(
            f'No {which} station is common to all models AND on every '
            f'model\'s node grid. Discovered per model: '
            + '; '.join(f'{m.label()}={len(d)}' for m, d in zip(models, per))
            + f'. {len(excl)} candidate(s) excluded -- first reasons: '
            + '; '.join(f'{station_name(k, which)}: {r}' for k, r in excl[:3]))

    # deterministic, spread evenly through the sorted station list
    n_show = min(args.max_stations, len(keep))
    idx = np.unique(np.linspace(0, len(keep) - 1, n_show).round().astype(int))
    shown = [keep[i] for i in idx]

    quants = FAULT_COLS if which == 'fault' else BODY_COLS
    order = ['slip', 'sliprate', 'shear'] if which == 'fault' else ['disp', 'vel']
    comp = args.component
    ncol, nrow = len(order), len(shown)

    # ---- layout in INCHES (rule 7: fix the layout, never shrink the font).
    # Space for the legend, the column titles, the shared x-label and the
    # computed caption is RESERVED before the axes are placed, so nothing can
    # land on top of anything else at any --print-width-mm.
    W = args.width_in
    L, Rm, PANEL_H, HS_IN, WS_IN = 0.62, 0.10, 1.05, 0.22, 0.52
    legend_rows = int(np.ceil(len(models) / 2))
    top_in = 0.26 + 0.20 * legend_rows + 0.22      # legend + column titles
    axes_h = nrow * PANEL_H + (nrow - 1) * HS_IN
    panel_w = (W - L - Rm - (ncol - 1) * WS_IN) / ncol

    cap_probe = ['x'] * 6
    fig = plt.figure(figsize=(W, 1.0))              # resized below
    axes = None

    def build(cap_in):
        bot = cap_in + 0.40                          # caption + shared x-label
        H = top_in + axes_h + bot
        f = plt.figure(figsize=(W, H))
        ax = f.subplots(nrow, ncol, squeeze=False, sharex=True)
        f.subplots_adjust(left=L / W, right=1.0 - Rm / W,
                          top=1.0 - top_in / H, bottom=bot / H,
                          hspace=HS_IN / PANEL_H, wspace=WS_IN / panel_w)
        return f, ax, H, bot

    plt.close(fig)
    fig, axes, H, bot_in = build(0.0)               # provisional; rebuilt below

    notes, tmax_all, peaks = [], 0.0, {q: [] for q in order}
    for r, key in enumerate(shown):
        for m in models:
            a, note = load_station(per[m.index][key], which)
            if note:
                notes.append(note)
            tmax_all = max(tmax_all, float(a[:, 0].max()))
            for c, q in enumerate(order):
                spec = quants[q]
                col = spec.get(comp)
                if col is None or col >= a.shape[1]:
                    continue
                axes[r][c].plot(a[:, 0], a[:, col], **m.style())
                peaks[q].append(float(np.max(np.abs(a[:, col]))))
        axes[r][0].set_ylabel(station_name(key, which), fontsize=8,
                              labelpad=1.5)
    for c, q in enumerate(order):
        spec = quants[q]
        axes[0][c].set_title(f'{comp}-{spec["name"]} ({spec["unit"]})',
                             fontsize=9.0, pad=3)
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
               columnspacing=1.4, borderpad=0.2, labelspacing=0.3,
               bbox_to_anchor=(0.5, 1.0 - 0.04 / H))

    kind = 'on-fault' if which == 'fault' else 'off-fault (body)'
    reasons = {}
    for k, r in excl:
        reasons[r] = reasons.get(r, 0) + 1
    cap = [FRAMING,
           f'{kind} stations, {comp}-component; rows are stations, columns '
           f'quantities, y-range shared down each column. Node spacings: '
           + ', '.join(f'{m.label()}' for m in models) + '.',
           f'{len(shown)} of {len(keep)} comparable stations shown (evenly '
           f'sampled from the sorted list); {len(excl)} of '
           f'{len(keep) + len(excl)} discovered stations excluded because '
           f'their coordinates are not shared exactly -- none snapped, none '
           f'interpolated.',
           'Peak |' + comp + '| over the panels shown: ' + '; '.join(
               f'{quants[q]["name"]} {max(peaks[q]):.4g} {quants[q]["unit"]}'
               for q in order if peaks[q])
           + f'; traces span 0-{tmax_all:.2f} s.']
    if notes:
        cap.append('Header/data column-count mismatch (data trusted): '
                   + '; '.join(sorted(set(notes))[:2]) + '.')
    # now that the caption is known, give the figure exactly the height it
    # needs and re-place the axes -- no overlap, no font below 7 pt
    text, cap_in = caption_block(cap, args)
    bot_in = cap_in + 0.40
    H_new = top_in + axes_h + bot_in
    fig.set_size_inches(W, H_new)
    fig.subplots_adjust(left=L / W, right=1.0 - Rm / W,
                        top=1.0 - top_in / H_new, bottom=bot_in / H_new,
                        hspace=HS_IN / PANEL_H, wspace=WS_IN / panel_w)
    fig.legends[0].set_bbox_to_anchor((0.5, 1.0 - 0.04 / H_new),
                                      transform=fig.transFigure)
    fig.text(L / W + (1.0 - (L + Rm) / W) / 2.0, (cap_in + 0.03) / H_new,
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


# ------------------------------------------------------------------ output


def caption_block(lines, args):
    """Wrap the computed caption to the PRINT width and return the text plus
    the vertical space it needs, in inches, so the caller can reserve that
    space BEFORE creating the figure -- rule 7, layout before font shrinking:
    the caption never overlaps an axis and the font never drops below 7 pt."""
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


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--plot', required=True,
                    choices=('cplot', 'ts-fault', 'ts-body'))
    ap.add_argument('--models', nargs='+', required=True,
                    help='ordered list of model paths (our run/reference dirs '
                         'or scec_archive submission dirs); style and label '
                         'are assigned by position, so the same command '
                         'always yields the same figure')
    ap.add_argument('--labels', nargs='*', default=None,
                    help='optional per-model legend labels, same order')
    ap.add_argument('--component', default='h', choices=('h', 'v', 'n'),
                    help='h = along-strike, v = down-dip, n = fault-normal '
                         '(n is body-station only)')
    ap.add_argument('--max-stations', type=int, default=4)
    ap.add_argument('--sentinel', type=float, default=None,
                    help='override the never-ruptured fill value (default: '
                         'detected per file; 1000 s in 2015-2020 archives, '
                         '99999 s in 2024 ones and in our writer)')
    ap.add_argument('--out', default=None)
    ap.add_argument('--print-width-mm', type=float, default=190.0,
                    help='printed width: Elsevier 90/140/190 mm, AGU 146 mm')
    ap.add_argument('--dpi', type=int, default=400)
    args = ap.parse_args(argv)

    if args.labels and len(args.labels) != len(args.models):
        raise SystemExit('--labels must have the same length as --models')
    args.width_in = args.print_width_mm / 25.4          # k = 1.0 by construction
    # DejaVu Sans averages ~3.95 pt per character; wrap the caption from the
    # PRINTED width so it cannot run off the canvas at any --print-width-mm.
    args.wrap_chars = max(40, int(args.width_in * 72.0 / 3.95 * 0.96))
    if args.out is None:
        args.out = f'scec_compare_{args.plot.replace("-", "_")}.png'

    models = [Model(p, i, (args.labels[i] if args.labels else None), args.sentinel)
              for i, p in enumerate(args.models)]

    print('=== models (style assigned by position in --models) ===')
    for m in models:
        print(f'  [{m.index}] {m.kind:8s} {m.path}')
        print(f'        label="{m.label()}"  colour={m.color}  '
              f'linestyle={m.linestyle}')
        es = None
        for probe in ('cplot', 'frt.canonical.txt'):
            f = os.path.join(m.dir, probe)
            if os.path.isfile(f):
                es = header_element_size_m(f)
        if es is not None:
            print(f'        header element_size = {es:g} m (reported only; '
                  f'spacing used in figures is measured from coordinates)')

    if args.plot == 'cplot':
        plot_cplot(models, args)
    else:
        plot_ts(models, args.plot.split('-')[1], args)
    return 0


if __name__ == '__main__':
    sys.exit(main())
