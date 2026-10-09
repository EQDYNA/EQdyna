#! /usr/bin/env python3
"""
scec_readers.py -- the READER LAYER of the SCEC overlay tooling.

Every on-disk layout a SCEC/USGS TPV comparison meets is one REGISTERED
reader below, and every reader normalises to the same two records:

    rupture field   dict(strike, downdip, t, ...) -- metres, down-dip
                    POSITIVE DOWN, t in s with never-ruptured masked to NaN
    station series  (array, note) -- the numeric block of one SCEC station
                    file, columns exactly as the SCEC format defines them

Adding a format is adding a registry entry (RUPTURE_READERS /
FAULT_STATION_LAYOUTS / BODY_STATION_LAYOUTS), never a new script. Which
layout a directory is, is SNIFFED from what is on disk, never declared.

Layouts registered today (each has a fixture in
testsys/unit/test_scec_compare.py that fails if its reader breaks):

  ours        frt.canonical.txt, frt.txt<rank>, SCECRuptureTime.txt;
              faultst*.txt / faultstft<N>_*.txt / body*.txt
  archive     cplot + faultst* / body* (no extension), 2015-2024 eras: the
              never-ruptured sentinel is 1000 s (2015-2020) or 99999 s (2024)
  CVWS        the same plus multi-fault cplot_<N>[.txt] / fault<N>st*
  modellers   (FaultMod, DayFD, SPECFEM3D, ...; 1e9/1e10 s sentinels) and
              the branch layout cplot_main / cplot_branch / branchst*

FORMAT HAZARDS HANDLED EXPLICITLY (each one has cost time before)
-----------------------------------------------------------------
* SENTINELS DIFFER BY WRITER AND BY ERA. 1000.0 s in the 2015-2020 archives,
  99999.0 s in the 2024 ones and in our writer, 1e9/1e10 s in the CVWS
  modellers'. A `> 1e4` mask lets 1000 through and it contours as a real
  rupture front. detect_sentinel() finds it per file from the data.
* THE cplot HEADER ENDS IN AN UNPARSEABLE ` j  k  t` LINE (and FaultMod's
  files carry a second `#` block AFTER it). read_numeric() scans the header.
* COLUMN LAYOUT DIFFERS. Ours (frt.canonical.txt): x, y, z (negative down),
  t. Archive cplot: along-strike, down-dip (POSITIVE DOWN), t.
* DOWN-DIP IS NOT DEPTH ON A DIPPING FAULT: hypot(y, z) when the node cloud
  is dipping (|corr(y,z)| > 0.95), -z when vertical (rough faults have y =
  roughness noise, uncorrelated with z). Decided from the data and printed.
* MULTIPLE FAULTS IN ONE frt. frt.txt<rank> is never split per fault (only
  station files are tagged ft<N>_). With --fault N: fault 1 is the y = 0
  plane (code convention); fault >= 2 is every other node, its along-strike
  coordinate x for a fault parallel to fault 1 (TPV22/23 stepover) and the
  distance from the junction with y = 0 for a branch (TPV24/25).
* OUR ON-FAULT STATION dp ON A DIPPING FAULT CAN BE DEPTH. The writer divides
  by sin(fltxyz(2,4,j)) (library_output.f90 output_onfault_st), which is 90
  deg whenever C_degen == 0 (TPV12/13), so `dp013` there is 1.3 km DEPTH,
  1.5 km down-dip. Decided from the data per model (see
  Model.station_files) and printed; correct writers (TPV36/37) are untouched.
* A STATION HEADER CAN LIE ABOUT ITS OWN COLUMN COUNT (pre-item-67 writer):
  the DATA is trusted and the mismatch reported.
* THE FORTRAN e15.7 EXPONENT LOSES ITS 'E' at 3-digit exponents
  ("0.1341312-114"); repaired IN MEMORY on read, never in the run dir.
* cplot HEADERS ARE NEVER PROVENANCE (tpv30's 100 m cplot still says
  `problem = LVFZ`): labels come from the path, spacing is MEASURED.

NEVER MODIFIES its inputs -- every file is opened read-only.
"""
import glob
import os
import re

import numpy as np

# SCEC station-name convention: distances in units of 100 m.
STATION_UNIT_M = 100.0

# Column layout of the SCEC station time series (identical in every archive
# and in today's writer, library_output.f90 output_onfault_st/offfault_st).
FAULT_COLS = {
    'slip':     dict(h=1, v=4, unit='m',   name='slip'),
    'sliprate': dict(h=2, v=5, unit='m/s', name='slip rate'),
    'shear':    dict(h=3, v=6, n=7, unit='MPa', name='shear stress'),
}
FAULT_DECLARED_MIN = 8
BODY_COLS = {
    'disp': dict(h=1, v=3, n=5, unit='m',   name='displacement'),
    'vel':  dict(h=2, v=4, n=6, unit='m/s', name='velocity'),
}
BODY_DECLARED_MIN = 7

RE_NCOL_DECL = re.compile(r'Time series in\s+(\d+)\s+columns', re.I)
RE_ELEMSIZE = re.compile(r'element_size\s*=\s*([0-9.eEdD+-]+)')
RE_LOST_E = re.compile(r'(\d\.\d+)([-+]\d\d\d+)')     # 0.1341312-114


# =================================================================== text I/O

def _floats(s):
    """Parse a whitespace line as floats, repairing a lost-'E' exponent."""
    try:
        return [float(t) for t in s.split()]
    except ValueError:
        return [float(t) for t in RE_LOST_E.sub(r'\1E\2', s).split()]


def read_numeric(path, maxcols=None):
    """Numeric block of a SCEC-style text file.

    Both the archive and our writer put an unparseable field-name line AFTER
    the '#' block, which defeats np.loadtxt(comments='#'). Scan only the
    header for the first all-numeric line, then let loadtxt do the bulk in C
    (the 25 m tpv30 cplot is 1.28 M rows). A lost-'E' exponent anywhere makes
    loadtxt raise; only then is the file re-parsed line by line with the
    repair, so every file that parsed before parses identically now.
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
                _floats(s)
            except ValueError:
                skip = i + 1
                continue
            skip = i
            break
    try:
        a = np.loadtxt(path, comments='#', skiprows=skip, ndmin=2)
    except ValueError:
        rows = []
        with open(path, 'r') as fh:
            for i, line in enumerate(fh):
                s = line.strip()
                if i < skip or not s or s.startswith('#'):
                    continue
                rows.append(_floats(s))
        a = np.array(rows, float, ndmin=2)
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
    """REPORTING ONLY -- spacing used anywhere is measured from coordinates
    (archive headers are demonstrably stale)."""
    for line in header_text(path):
        m = RE_ELEMSIZE.search(line)
        if m:
            try:
                return float(m.group(1).replace('D', 'E').replace('d', 'e'))
            except ValueError:
                return None
    return None


def detect_sentinel(t, override=None):
    """The 'never ruptured' fill value, from the data: the maximum is
    repeated AND >= 10x the 99.9th percentile of everything below it --
    true of 1000 / 99999 / 1e9 / 1e10 by orders of magnitude, false for a
    genuine latest arrival. -> (sentinel_or_None, how_string)."""
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


def measure_dx(xs, dd):
    """Node spacing MEASURED from the coordinates (headers are stale)."""
    out = []
    for v in (xs, dd):
        u = np.unique(np.round(np.asarray(v, float), 3))
        if u.size > 1:
            out.append(float(np.median(np.diff(u))))
    return float(np.median(out)) if out else float('nan')


# ===================================================== rupture-field readers

def frt_to_strike_dip(a):
    """frt columns x, y, z (negative down), t -> (strike, downdip, t, rule).
    Down-dip is -z on a vertical fault and hypot(y, z) on a dipping one --
    decided FROM THE DATA (a rough vertical fault has y = roughness)."""
    x, y, z, t = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    if np.std(y) > 1.0 and np.std(z) > 1.0:
        r = float(np.corrcoef(y, z)[0, 1])
    else:
        r = 0.0
    if abs(r) > 0.95:
        return x, np.hypot(y, z), t, f'dipping fault (corr(y,z)={r:+.3f}): down-dip = hypot(y,z)'
    return x, -z, t, f'vertical fault (corr(y,z)={r:+.3f}): down-dip = -z'


FAULT1_Y_TOL_M = 0.5   # fault 1 is the code y = 0 plane


def select_fault(a, fault):
    """Rows of an frt array belonging to fault `fault` and the along-strike
    coordinate to use for them. frt.txt<rank> is never split per fault, so
    the split is geometric: fault 1 = the y = 0 plane, fault >= 2 = the rest.
    -> (rows, strike_m, rule)."""
    on1 = np.abs(a[:, 1]) < FAULT1_Y_TOL_M
    if fault == 1:
        return a[on1], a[on1, 0], 'fault 1 = nodes with |y| < 0.5 m'
    rest = a[~on1]
    if fault != 2:
        raise SystemExit(f'--fault {fault}: an frt split is defined for fault '
                         f'1 (y = 0) and fault 2 (every other node) only')
    if rest.shape[0] == 0:
        raise SystemExit('--fault 2: every frt node is on y = 0 -- single fault')
    x, y = rest[:, 0], rest[:, 1]
    if np.ptp(y) < 1.0:
        return rest, x, (f'fault 2 = nodes off y = 0, parallel at y = '
                         f'{np.median(y):g} m: along-strike = x')
    slope, icept = np.polyfit(x, y, 1)
    xj = -icept / slope
    s = np.hypot(x - xj, y)
    return rest, s, (f'fault 2 = nodes off y = 0, a BRANCH at '
                     f'{np.degrees(np.arctan(abs(slope))):.2f} deg meeting '
                     f'y = 0 at x = {xj:.1f} m: along-strike = distance '
                     f'from that junction')


def _read_frt_array(a, fault):
    if a.shape[1] < 4:
        raise SystemExit(f'frt: {a.shape[1]} columns, need >= 4 (x, y, z, t)')
    if fault is None:
        return frt_to_strike_dip(a)
    rows, strike, rule = select_fault(a, fault)
    _, dd, t, dip_rule = frt_to_strike_dip(rows)
    return strike, dd, t, rule + '; ' + dip_rule


def _probe_cplot(d, fault):
    if fault in (None, 1):
        names = ['cplot']
    else:
        names = []
    if fault is not None:
        names += [f'cplot_{fault}', f'cplot_{fault}.txt']
        names += {1: ['cplot_main'], 2: ['cplot_branch']}.get(fault, [])
    for n in names:
        if os.path.isfile(os.path.join(d, n)):
            return os.path.join(d, n)
    return None


def _read_cplot(path, fault):
    a = read_numeric(path, maxcols=3)
    if a.shape[1] != 3:
        raise SystemExit(f'{path}: {a.shape[1]} columns, expected 3 '
                         f'(strike, down-dip, t)')
    # some CVWS writers (SPECFEM3D) carry a signed down-dip; the archive
    # convention is positive down
    return a[:, 0], np.abs(a[:, 1]), a[:, 2], 'file is already (strike, down-dip positive-down, t)'


def _probe_file(name):
    def probe(d, fault):
        p = os.path.join(d, name)
        return p if os.path.isfile(p) else None
    return probe


def _read_frt_file(path, fault):
    return _read_frt_array(read_numeric(path), fault)


def _probe_frt_ranks(d, fault):
    return os.path.join(d, 'frt.txt*') if glob.glob(os.path.join(d, 'frt.txt*')) else None


def _read_frt_ranks(pattern, fault):
    files = sorted(glob.glob(pattern))
    return _read_frt_array(np.vstack([read_numeric(f) for f in files]), fault)


def _probe_scec_rt(d, fault):
    if fault not in (None, 1):
        return None
    p = os.path.join(d, 'SCECRuptureTime.txt')
    return p if os.path.isfile(p) else None


# (name, kind, probe(dir, fault) -> path|None, read(path, fault) -> 4-tuple)
# PROBE ORDER IS THE PRE-REFACTOR ORDER (cplot, frt.canonical.txt,
# SCECRuptureTime.txt, frt.txt*), so a directory that held several of these
# resolves to the same file it always did.
RUPTURE_READERS = [
    ('cvws-cplot',          'archive', _probe_cplot,                     _read_cplot),
    ('eqdyna-frt-canonical', 'eqdyna', _probe_file('frt.canonical.txt'), _read_frt_file),
    ('eqdyna-scec-rt',       'eqdyna', _probe_scec_rt,                   _read_cplot),
    ('eqdyna-frt-ranks',     'eqdyna', _probe_frt_ranks,                 _read_frt_ranks),
]


def read_rupture_file(path, fault=None):
    """A single rupture file given by path: frt* is ours, anything else is a
    3-column (strike, down-dip, t) cplot-style file."""
    if os.path.basename(path).startswith('frt'):
        return _read_frt_file(path, fault)
    return _read_cplot(path, fault)


# ===================================================== station-file layouts

# Each layout: (name, glob patterns, compiled name regex). The regex's groups
# are the station coordinates in 100 m units: (strike, dip) for a fault,
# (offset, strike, dip) for a body station.
def fault_station_layouts(fault):
    if fault in (None, 1):
        return [('scec-faultst', ['faultst*'],
                 re.compile(r'^(?:faultst|fault1st)(-?\d+)dp(-?\d+)(?:\.txt)?$'),
                 'faultstft'),
                ('cvws-fault1st', ['fault1st*'],
                 re.compile(r'^fault1st(-?\d+)dp(-?\d+)(?:\.txt)?$'), None)]
    out = [('eqdyna-ft-tag', [f'faultstft{fault}_*'],
            re.compile(rf'^faultstft{fault}_(-?\d+)dp(-?\d+)(?:\.txt)?$'), None),
           ('cvws-faultNst', [f'fault{fault}st*'],
            re.compile(rf'^fault{fault}st(-?\d+)dp(-?\d+)(?:\.txt)?$'), None)]
    if fault == 2:
        out.append(('cvws-branchst', ['branchst*'],
                    re.compile(r'^branchst(-?\d+)dp(-?\d+)(?:\.txt)?$'), None))
    return out


BODY_STATION_LAYOUTS = [
    ('scec-body', ['body*'],
     re.compile(r'^body(-?\d+)st(-?\d+)dp(-?\d+)(?:\.txt)?$'), None),
]


def discover_stations(d, which, fault=None):
    """-> (found {coord-tuple-in-m: path}, bad [unparseable basenames]).
    Stations are discovered by PARSING names, never from a list."""
    layouts = fault_station_layouts(fault) if which == 'fault' else BODY_STATION_LAYOUTS
    found, bad, seen = {}, [], set()
    for _, pats, rx, exclude_prefix in layouts:
        for pat in pats:
            for f in sorted(glob.glob(os.path.join(d, pat))):
                base = os.path.basename(f)
                if f in seen or base.endswith('.zip') or os.path.isdir(f):
                    continue
                if exclude_prefix and base.startswith(exclude_prefix):
                    continue          # another fault's tagged file
                seen.add(f)
                m = rx.match(base)
                if not m:
                    bad.append(base)
                    continue
                found[tuple(int(v) * STATION_UNIT_M for v in m.groups())] = f
    return found, bad


def load_station(path, which):
    """-> (array, note). The data is trusted over a header that misstates
    its own column count; the mismatch is returned, not hidden."""
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


def _sccode(v):
    """SCEC filename convention: zero-padded magnitude, '-' only if negative."""
    n = int(round(v / STATION_UNIT_M))
    return f'-{abs(n):03d}' if n < 0 else f'{n:03d}'


def station_name(key, which):
    if which == 'fault':
        s, d = key
        return f'faultst{_sccode(s)}dp{_sccode(d)}'
    b, s, d = key
    return f'body{_sccode(b)}st{_sccode(s)}dp{_sccode(d)}'


def parse_station_request(name, which):
    """A --stations entry -> coordinate key. Any layout's prefix is accepted
    (faultst / fault2st / branchst / faultstft2_ / body), only the numbers
    matter."""
    rx = (r'st(?:ft\d+_)?(-?\d+)dp(-?\d+)(?:\.txt)?$' if which == 'fault'
          else r'^body(-?\d+)st(-?\d+)dp(-?\d+)(?:\.txt)?$')
    m = re.search(rx, os.path.basename(name))
    if not m:
        raise SystemExit(f'--stations {name!r}: not a SCEC {which}-station name')
    return tuple(int(v) * STATION_UNIT_M for v in m.groups())


# ================================================================ the model

COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7',
          '#56B4E9', '#E69F00', '#000000', '#F0E442']
LINESTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1)),
              (0, (1, 1)), (0, (4, 1, 1, 1, 1, 1))]

RE_EQDYNA_ARCHIVE = re.compile(r'eqdyna3?d?-(v[\d.]+)-(\d+)m-(\d{4})$')
# CVWS modeller dirs as laid out in ~/shared_dataset/scec_cvws.*:
# <user>[.<n>]-<code>-<res>m-<year>, or a bare <user> dir.
RE_CVWS_ARCHIVE = re.compile(r'^([a-z]+)(?:\.\d+)?-([a-z0-9]+)-([\d.]+)m-(\d{4})$')
CODE_NAMES = {'faultmod': 'FaultMod', 'dayfd': 'DayFD', 'specfem3d': 'SPECFEM3D'}


class Model:
    """One result set: one of OUR run/reference directories, or one archive /
    CVWS submission directory. Which it is, is sniffed, never declared."""

    def __init__(self, path, index, label=None, sentinel=None, fault=None):
        self.path = os.path.abspath(path)
        self.index = index
        self.sentinel_override = sentinel
        self.fault = fault
        self.notes = []          # format problems worth reporting, not hiding
        self._dx_cache = None    # MEASURED spacing, filled by load_rupture()
        self._field = None
        self.color = COLORS[index % len(COLORS)]
        self.linestyle = LINESTYLES[index % len(LINESTYLES)]
        self.linewidth = float(np.clip(1.6 - 0.18 * index, 0.7, 1.6))
        self._sniff()
        self.user_label = label

    # -- discovery -------------------------------------------------------

    def _sniff(self):
        p = self.path
        self.reader = None
        if os.path.isfile(p):
            self.dir = os.path.dirname(p)
            self.rupture_file = p
            self.kind = 'eqdyna' if 'frt' in os.path.basename(p) else 'archive'
            self.reader = 'file'
            return
        if not os.path.isdir(p):
            raise SystemExit(f'{p}: no such file or directory')
        self.dir = p
        self.kind, self.rupture_file = None, None
        for name, kind, probe, _ in RUPTURE_READERS:
            f = probe(p, self.fault)
            if f:
                self.kind, self.rupture_file, self.reader = kind, f, name
                return
        if glob.glob(os.path.join(p, 'faultst*')):
            self.kind = 'archive' if not glob.glob(
                os.path.join(p, 'faultst*.txt')) else 'eqdyna'
            return
        if glob.glob(os.path.join(p, 'fault[0-9]st*')) or glob.glob(
                os.path.join(p, 'branchst*')):
            self.kind = 'archive'          # CVWS per-fault station layout
            return
        multi = sorted(os.path.basename(f) for f in glob.glob(os.path.join(p, 'cplot_*')))
        hint = (f' It holds per-fault files ({", ".join(multi)}): pass --fault N.'
                if multi and self.fault is None else '')
        raise SystemExit(
            f'{p}: not recognisable as a model directory -- no cplot, '
            f'frt.canonical.txt, SCECRuptureTime.txt, frt.txt* or faultst* '
            f'found{" for fault %s" % self.fault if self.fault else ""}.{hint}')

    @property
    def is_independent(self):
        """True for a CVWS submission by a code other than EQdyna. Decided
        from the directory name (headers are never provenance)."""
        if self.kind != 'archive':
            return False
        return not os.path.basename(self.dir).lower().startswith('eqdyna')

    def station_files(self, which):
        found, bad = discover_stations(self.dir, which, self.fault)
        if bad:
            why = (' -- pre-fix i3.3 strike overflow' if any('*' in b for b in bad)
                   else '')
            self.notes.append(f'{len(bad)} {which}-station file(s) with '
                              f'unparseable names skipped (e.g. {bad[0]}){why}')
        if which == 'fault' and self.kind == 'eqdyna':
            found = self._depth_named_dp(found)
        return found

    def _depth_named_dp(self, found):
        """Our writer names a dipping-fault station by DEPTH whenever its
        dip angle is the C_degen == 0 default of 90 deg (see module
        docstring). Decided from the data: only when the field is dipping,
        the names do NOT all land on the down-dip node grid read as down-dip,
        and DO all land on it read as depth / sin(dip)."""
        try:
            f = self.load_rupture()
        except SystemExit:
            return found
        sin_dip, ddx = f['sin_dip'], f['dx_dip_m']
        if sin_dip is None or not found or not np.isfinite(ddx):
            return found

        def on(v):
            r = abs(v) % ddx
            return r < 1e-6 or abs(r - ddx) < 1e-6
        if all(on(k[1]) for k in found):
            return found
        conv = {(k[0], round(k[1] / sin_dip / STATION_UNIT_M) * STATION_UNIT_M): v
                for k, v in found.items()}
        if not all(on(k[1]) for k in conv):
            return found
        self.notes.append(
            f'on-fault station dp read as DEPTH and converted to down-dip '
            f'(/ sin dip = {sin_dip:.4f}, rounded to 100 m): the names do not '
            f'land on the {ddx:g} m down-dip grid as down-dip, and all do as '
            f'depth -- writer bug, library_output.f90 output_onfault_st')
        return conv

    # -- rupture-time field ---------------------------------------------

    def load_rupture(self):
        """-> dict(strike, downdip, t, ...) with down-dip POSITIVE DOWN in
        metres, for every registered layout. Cached: a model is read once."""
        if self._field is not None:
            return self._field
        if self.rupture_file is None:
            raise SystemExit(f'{self.dir}: no rupture-time file for this model')
        if self.reader == 'file':
            xs, dd, t, rule = read_rupture_file(self.rupture_file, self.fault)
            src = os.path.basename(self.rupture_file)
        else:
            read = {n: r for n, _, _, r in RUPTURE_READERS}[self.reader]
            xs, dd, t, rule = read(self.rupture_file, self.fault)
            src = (f'{len(glob.glob(self.rupture_file))} x frt.txt<rank>'
                   if self.reader == 'eqdyna-frt-ranks'
                   else os.path.basename(self.rupture_file))
        sent, how = detect_sentinel(t, self.sentinel_override)
        tm = t.astype(float).copy()
        if sent is not None:
            tm[t >= sent] = np.nan
        self._dx_cache = measure_dx(xs, dd)
        u = np.unique(np.round(np.asarray(dd, float), 3))
        dx_dip = float(np.median(np.diff(u))) if u.size > 1 else float('nan')
        sin_dip = None
        if 'dipping' in rule and self.kind == 'eqdyna':
            # sin(dip) = depth / down-dip, from the same node cloud
            z = self._frt_depths()
            ok = dd > 1.0
            if z is not None and z.size == dd.size and ok.any():
                sin_dip = float(np.median(z[ok] / dd[ok]))
        self._field = dict(strike=xs, downdip=dd, t=tm, t_raw=t, sentinel=sent,
                           sentinel_how=how, dipwise=rule, source=src,
                           dx_m=self._dx_cache, dx_dip_m=dx_dip,
                           sin_dip=sin_dip)
        return self._field

    def _frt_depths(self):
        try:
            if self.reader == 'eqdyna-frt-ranks':
                a = np.vstack([read_numeric(f) for f in sorted(glob.glob(self.rupture_file))])
            else:
                a = read_numeric(self.rupture_file)
            if self.fault is not None:
                a, _, _ = select_fault(a, self.fault)
            return np.abs(a[:, 2])
        except (SystemExit, ValueError, IndexError):
            return None

    # -- labelling -------------------------------------------------------

    def label(self, dx_m=None):
        if self.user_label:
            base = self.user_label
        else:
            base = os.path.basename(self.dir) or self.dir
            bench = os.path.basename(os.path.dirname(self.dir))
            if self.kind == 'archive':
                m = RE_EQDYNA_ARCHIVE.match(base)
                c = RE_CVWS_ARCHIVE.match(base)
                if m:
                    # the dir's OWN declared resolution goes in the label, so
                    # three submissions of one benchmark never share a legend
                    # entry; the MEASURED spacing is appended beside it below
                    base = (f'{bench} SCEC {m.group(1)} '
                            f'{m.group(2)} m ({m.group(3)})')
                elif c:
                    who, code, res, yr = c.groups()
                    base = (f'{bench} CVWS {who.capitalize()} ({yr}), '
                            f'{CODE_NAMES.get(code, code)} {res} m')
                elif self.is_independent:
                    base = f'{bench} CVWS {base}'
        if self.fault is not None:
            base = f'{base} [fault {self.fault}]'
        if dx_m is None:
            dx_m = self._dx_cache
        if dx_m is not None and np.isfinite(dx_m):
            base = f'{base}, dx {dx_m:g} m'
        return base

    def style(self):
        return dict(color=self.color, linestyle=self.linestyle,
                    linewidth=self.linewidth)
