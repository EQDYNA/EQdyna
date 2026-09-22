#! /usr/bin/env python3
"""Shared data paths, loaders, print-geometry and shared colour/axis scales for
the TPV30 SCEC comparison figures (fig1 = cplot overlay, fig2 = station time
series).

WHAT THE FIGURES CLAIM, AND WHAT THEY DO NOT
--------------------------------------------
Both figures are EQdyna-Fortran vs EQdyna-Fortran, roughly a decade apart:
EQdyna's own committed dx = 500 m TPV30 gate reference against the OWNER'S OWN
2015 EQdyna v3.1 dx = 100 m submission to the SCEC/USGS cvws portal.  That is
a CROSS-VERSION CONSISTENCY CHECK -- "recognizably the same physics" -- and it
is NEITHER an independent validation NOR an accuracy claim.  The wording here
deliberately matches NOTES_tpv30_gate.md's own framing; see that file's
"Rule 17 step 6" section.

Two things a reader must be told, and which the figures themselves print:

  * The two runs are at DIFFERENT node spacings (500 m vs 100 m).  Nothing is
    interpolated anywhere in these figures.  Each rupture-time field is
    contoured on its OWN grid, and every node-to-node and station-to-station
    comparison uses EXACT coordinate coincidence only (dx = 500 divides
    dx = 100, so every 500 m node has an exact 100 m twin).  The rough fault
    surface is one fixed random realisation and must not be resampled
    (PROJECT_RULES rule 17 step 2; scripts/lib.requireFaultGeometryResolution
    refuses finer-than-source dx for the same reason).
  * The 2015 100 m submission's own standing in the portal is what it is --
    14/14, farthest from the group median among 14 TPV30 submissions, per
    scec_archive/tpv30/eqdyna-v3.1-100m-2015/PROVENANCE.md, improving
    monotonically 100 m (14/14) -> 50 m (12/14) -> 25 m (8/14).  Agreement
    with it therefore certifies NEITHER run against the benchmark.

The separate, open, unresolved python-numpy/python-jax vs Fortran divergence
on this same case (NOTES_tpv30_gate.md "STOP -- finding, not a landing";
pathway item 19(b)) is NOT touched, illustrated or implied by these figures.
test.tpv30 is deliberately unregistered in testNameList.py / testsys/matrix.py
and these figures do not change that.

PRINT GEOMETRY (rule 1: settled before any styling)
---------------------------------------------------
Printed width  : 146 mm  = 5.7480 in   (AGU text width, \\includegraphics at
                                        full text width, fraction 1.0)
Canvas width   : 5.7480 in             -> k = canvas / print = 1.000
Font multiplier: x k = x 1.000, i.e. every rcParam pt below IS the printed pt.
Raster         : 400 dpi -> 2299 px across at 146 mm (>= 300 dpi, rule 3).

DATA PATHS (declared here and nowhere else)
-------------------------------------------
  ARCHIVE_DIR  read-only 2015 submission.  scec_archive/ is gitignored and
               lives only in the primary checkout, so it is resolved through
               $EQDYNA_SCEC_ARCHIVE if set, else the primary checkout's path,
               else this tree's own scec_archive/.  NEVER written to.
  OURS_TS_DIR  the 13 on-fault station time series from the very run that
               produced test.reference.results/test.tpv30/frt.canonical.txt
               (verified bit-identical, max|diff| = 0.0 over 3321 x 22 values;
               see that directory's PROVENANCE.md).  Committed, because the
               frozen reference itself carries no time series -- only a final
               -state frt.canonical.txt and a time-independent fault.dyna.r.nc
               -- so fig2 would otherwise not be regenerable.
  REF_FRT      test.reference.results/test.tpv30/frt.canonical.txt, read via
               the committed comparison script's own loader (rule 1: reuse).
"""
import json
import os
import re
import sys

import numpy as np

FIGDIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(FIGDIR))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Reuse, do not reimplement: the committed report-only comparison scripts
# already own the cplot reader, the frt.canonical gridder, the case-parameter
# reader and the station-name -> metres mapping.
from testsys.parity import evidence_tpv30_scec_comparison as cmp30  # noqa: E402
from testsys.parity import evidence_tpv29_scec_comparison as cmp29  # noqa: E402

_PRIMARY_CHECKOUT = '/home/utig5/dliu/EQdyna'
ARCHIVE_DIR = (os.environ.get('EQDYNA_SCEC_ARCHIVE')
               or next((p for p in (
                   os.path.join(REPO_ROOT, 'scec_archive', 'tpv30',
                                'eqdyna-v3.1-100m-2015'),
                   os.path.join(_PRIMARY_CHECKOUT, 'scec_archive', 'tpv30',
                                'eqdyna-v3.1-100m-2015'))
                   if os.path.isdir(p)), ''))
OUT_DIR = os.path.join(REPO_ROOT, 'docs', 'figures', 'tpv30')
OURS_TS_DIR = os.path.join(OUT_DIR, 'data', 'eqdyna-500m-fortran-gate')
SCALES_JSON = os.path.join(OUT_DIR, 'shared_scales.json')

ARCHIVE_LABEL = 'EQdyna v3.1, 100 m, 2015 SCEC submission'
OURS_LABEL = 'EQdyna today, 500 m, committed gate reference'

# The four on-fault stations fig2 shows.  SELECTION RULE, applied before
# looking at any result: of the 24 stations the TPV29/30 spec asks for, our
# 500 m run resolves 13 at EXACTLY the archive's own coordinates (the other 11
# spec stations are not on 500 m multiples and are excluded rather than
# snapped -- zero offset, no interpolation, nothing nearest-neighbour).  From
# those 13, one is taken in each quadrant of the fault plane about the
# hypocenter: backward / forward along strike x near-surface / deep.  That is
# what "spans the fault" means here; no station was swapped after plotting.
STATIONS = ('faultst-150dp120', 'faultst-050dp000',
            'faultst050dp120', 'faultst150dp000')

# 8-column SCEC on-fault time-series layout, identical in the 2015 archive and
# in today's writer (the archive header says 8 columns; ours says "11 columns"
# and then writes 8 -- a real header bug in today's writer, reported, not
# worked around silently: the DATA is 8 columns in both and is read as such).
TS_COLS = ('t', 'h_slip', 'h_slip_rate', 'h_shear', 'v_slip', 'v_slip_rate',
           'v_shear', 'n_stress')
# Which column each fig2 panel column shows, with its printed axis label.
# Short headers: at 146 mm each of the three columns is ~1.5 in wide, so
# "Along-strike slip rate (m s-1)" simply does not fit -- the string is
# shortened and the "along-strike" qualifier stated once in the figure note
# (rule 7: fewer strings, not smaller ones).
TS_PANELS = (('h_slip', 'Slip (m)'),
             ('h_slip_rate', r'Slip rate (m s$^{-1}$)'),
             ('h_shear', 'Shear stress (MPa)'))

# ------------------------------------------------------------ print geometry
PRINT_WIDTH_MM = 146.0
PRINT_WIDTH_IN = PRINT_WIDTH_MM / 25.4
CANVAS_WIDTH_IN = PRINT_WIDTH_IN
K = CANVAS_WIDTH_IN / PRINT_WIDTH_IN
DPI = 400
PT = dict(label=9.5, tick=7.5, title=9.5, legend=7.5, annot=7.0)


def apply_style():
    """One serif family, every size multiplied by k (= 1.0 here), so the
    rcParam pt values ARE the printed pt values."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['DejaVu Serif'],
        'mathtext.fontset': 'dejavuserif',
        'axes.labelsize': PT['label'] * K,
        'axes.titlesize': PT['title'] * K,
        'xtick.labelsize': PT['tick'] * K,
        'ytick.labelsize': PT['tick'] * K,
        'legend.fontsize': PT['legend'] * K,
        'font.size': PT['annot'] * K,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'lines.linewidth': 0.8,
        'legend.frameon': False,
        'savefig.dpi': DPI,
        'figure.dpi': DPI,
    })
    return plt


def print_geometry_note():
    return (f'Printed width {PRINT_WIDTH_MM:.0f} mm; canvas '
            f'{CANVAS_WIDTH_IN:.3f} in, k = {K:.3f}; {DPI} dpi '
            f'({round(CANVAS_WIDTH_IN * DPI)} px across).')


# -------------------------------------------------------------------- fields
def load_fields():
    """Both rupture-time fields, each on its own native grid, plus the case's
    own parameters.  Nothing is regridded or interpolated."""
    if not ARCHIVE_DIR or not os.path.isfile(os.path.join(ARCHIVE_DIR, 'cplot')):
        raise SystemExit(
            f'no readable 2015 archive cplot at {ARCHIVE_DIR!r}.  scec_archive/ '
            f'is gitignored; set $EQDYNA_SCEC_ARCHIVE to the '
            f'eqdyna-v3.1-100m-2015 directory.')
    p = cmp30._case_params(cmp30.CASE)
    ref = cmp30.load_reference(cmp30.REFERENCE_FILE, p)
    x15, dip15, T15 = cmp29.load_rupture_grid(os.path.join(ARCHIVE_DIR, 'cplot'))
    # Same smoke check the comparison script does: stop rather than draw a
    # misaligned overlay.
    if not (float(x15.min()) == float(ref['x'].min()) == p['fxmin']
            and float(x15.max()) == float(ref['x'].max()) == p['fxmax']
            and float(dip15.min()) == 0.0 and float(dip15.max()) == -p['fzmin']):
        raise SystemExit('strike/dip extents disagree between the 500 m '
                         'reference and the 100 m archive grid -- refusing to '
                         'draw a misaligned overlay.')
    return p, ref, x15, dip15, T15


def masked_times(T, sentinel):
    """Rupture time with never-ruptured nodes set to NaN, so contour() and
    pcolormesh() leave them blank instead of drawing through them."""
    out = np.array(T, dtype=float)
    out[~(out < sentinel)] = np.nan
    return out


def colocated_difference(ref, x15, dip15, T15):
    """Delta t = ours - 2015 at the 3321 EXACTLY coincident nodes (dx = 500
    divides dx = 100), plus the two extent-disagreement node sets.  Raises if
    any 500 m node has no exact 100 m twin -- a grid bug must not look like
    noise."""
    ix = np.rint((ref['x'] - x15[0]) / (x15[1] - x15[0])).astype(int)
    iz = np.rint((ref['dip'] - dip15[0]) / (dip15[1] - dip15[0])).astype(int)
    if (np.abs(x15[ix] - ref['x']).max() > 1e-6
            or np.abs(dip15[iz] - ref['dip']).max() > 1e-6):
        raise SystemExit('500 m nodes do not land exactly on 100 m archive '
                         'nodes -- coordinate mapping is wrong.')
    T15_at = T15[np.ix_(iz, ix)]
    r15 = T15_at < cmp30.SENTINEL_ARCHIVE
    rnow = ref['rupt'] < cmp30.SENTINEL_CURRENT
    dt = np.where(r15 & rnow, ref['rupt'] - T15_at, np.nan)
    return dict(dt=dt, both=r15 & rnow, only_ours=rnow & ~r15,
                only_2015=r15 & ~rnow, n_nodes=int(r15.size))


# ------------------------------------------------------------- time series --
_EXP_FIX = re.compile(r'(\d\.\d+)([-+]\d{2,})')


def load_timeseries(path):
    """One SCEC on-fault station file -> dict of 1-D arrays.

    Applies scripts/correctSCECStOutputFormat.py's own repair in memory
    (Fortran E15.7 drops the 'E' when the exponent needs three digits), so no
    line is ever silently dropped as unparseable -- which is what a
    "skip lines that do not parse" reader would do to exactly the smallest,
    most numerous values."""
    rows = []
    with open(path, errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = _EXP_FIX.sub(r'\1E\2', s).split()
            try:
                vals = [float(v) for v in parts]
            except ValueError:
                continue          # the one 't h-slip ...' field-name row
            if len(vals) == len(TS_COLS):
                rows.append(vals)
    a = np.asarray(rows)
    if a.ndim != 2 or a.shape[1] != len(TS_COLS):
        raise SystemExit(f'{path}: could not read {len(TS_COLS)} numeric '
                         f'columns (got {a.shape})')
    return {c: a[:, i] for i, c in enumerate(TS_COLS)}


def load_station_pair(name):
    """(ours, 2015) for one station, both at the SAME physical point --
    membership of STATIONS already guarantees exact coincidence."""
    ours_p = os.path.join(OURS_TS_DIR, name + '.txt')
    arch_p = os.path.join(ARCHIVE_DIR, name)
    for p in (ours_p, arch_p):
        if not os.path.isfile(p):
            raise SystemExit(f'{p}: missing station time series')
    return load_timeseries(ours_p), load_timeseries(arch_p)


# ---------------------------------------------------- shared, cached scales --
def shared_scales(recompute=False):
    """Colour and axis ranges computed ONCE from the UNION of both datasets and
    cached to shared_scales.json, so every panel meant to be compared carries
    the identical range and so reruns are stable (rule 5).  Delete the JSON or
    pass --recompute-scales to rebuild it."""
    if os.path.isfile(SCALES_JSON) and not recompute:
        with open(SCALES_JSON) as f:
            return json.load(f)

    p, ref, x15, dip15, T15 = load_fields()
    d = colocated_difference(ref, x15, dip15, T15)
    t_ours = masked_times(ref['rupt'], cmp30.SENTINEL_CURRENT)
    t_2015 = masked_times(T15, cmp30.SENTINEL_ARCHIVE)
    tmax = float(np.ceil(max(np.nanmax(t_ours), np.nanmax(t_2015)) * 2) / 2)

    adt = np.abs(d['dt'][np.isfinite(d['dt'])])
    # Symmetric limit = smallest 0.5 s multiple covering 90% of the
    # co-ruptured nodes.  NOT the maximum (9.6 s) and not the 98th percentile
    # (6.8 s): the distribution is heavy-tailed, and on a +-7 s bar the 1.08 s
    # median -- the number the figure exists to let a reader see -- is
    # invisible.  The saturated fraction is stated on the figure and the bar
    # is drawn with both extend arrows, so nothing is hidden by the choice.
    dtlim = float(np.ceil(np.percentile(adt, 90) * 2) / 2)

    sc = dict(rupture_time_s=[0.0, tmax],
              rupture_time_contour_interval_s=2.0,
              rupture_time_band_interval_s=0.5,
              dt_s=[-dtlim, dtlim],
              dt_median_s=float(np.median(adt)),
              dt_p90_s=float(np.percentile(adt, 90)),
              dt_max_s=float(adt.max()),
              dt_saturated_pct=float(100.0 * (adt > dtlim).mean()),
              n_both=int(d['both'].sum()),
              n_only_ours=int(d['only_ours'].sum()),
              n_only_2015=int(d['only_2015'].sum()),
              extent_overlap_pct=100.0 * d['both'].sum()
              / max(1, d['both'].sum() + d['only_2015'].sum()))

    for key, _label in TS_PANELS:
        lo, hi = np.inf, -np.inf
        for st in STATIONS:
            a, b = load_station_pair(st)
            lo = min(lo, a[key].min(), b[key].min())
            hi = max(hi, a[key].max(), b[key].max())
        pad = 0.05 * (hi - lo)
        sc['ts_' + key] = [float(lo - pad), float(hi + pad)]
    tend = min(max(load_station_pair(STATIONS[0])[0]['t']),
               max(load_station_pair(STATIONS[0])[1]['t']))
    sc['ts_t_s'] = [0.0, float(np.ceil(tend))]

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(SCALES_JSON, 'w') as f:
        json.dump(sc, f, indent=2, sort_keys=True)
    return sc


def footnote(fig, text, width_chars=116, extra_in=0.0):
    """Place the framing note INSIDE the canvas and reserve space for it.

    Rule 1 is about the SAVED width, not the requested figsize: an unwrapped
    fig.text overflows the canvas and bbox_inches='tight' then saves a wider
    PNG, which silently multiplies k and shrinks every printed font.  (That
    happened here on the first render -- 4051 px = 10.1 in wide instead of
    5.748 in, i.e. k = 1.76 and 9.5 pt labels printing at 5.4 pt.)  So the
    text is hard-wrapped to the canvas, constrained_layout's rect is shrunk to
    make room, and the caller saves WITHOUT bbox_inches."""
    import textwrap
    wrapped = '\n'.join(textwrap.fill(p, width_chars)
                        for p in text.strip().split('\n'))
    nlines = wrapped.count('\n') + 1
    h_in = nlines * PT['annot'] * K * 1.35 / 72.0 + 0.07
    frac = (h_in + extra_in) / fig.get_figheight()
    fig.get_layout_engine().set(rect=(0.0, frac, 1.0, 1.0 - frac))
    fig.text(0.004, 0.004, wrapped, fontsize=PT['annot'] * K, va='bottom',
             ha='left', linespacing=1.35)
    # Figure-fraction y of the top of the note block: anything the caller
    # reserved `extra_in` for (a figure legend, typically) anchors here.
    return h_in / fig.get_figheight()


def endpoint_ticks(lo, hi, mid=None):
    """Colour bars are ALWAYS ticked at both endpoints plus the midpoint (or
    zero for a diverging scale) -- rule 4."""
    m = 0.5 * (lo + hi) if mid is None else mid
    return [lo, m, hi]
