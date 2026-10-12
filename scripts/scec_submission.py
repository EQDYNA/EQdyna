#! /usr/bin/env python3
"""
scec_submission.py -- package one EQdyna run directory as a SCEC/USGS CVWS
submission zip, in the layout of `scec_archive/<tpv>/<version-dx-year>/`.

    scec_submission.py RUN_DIR --tpv 26 --out tpv26_eqdyna_100m.zip
        [--commit SHA] [--date YYYY/MM/DD] [--label eqdyna-v5.25.0-100m-2026]

WHAT GOES IN THE ZIP (one flat directory, `<tpv>/<label>/`):

  faultst<s>dp<d>         one per on-fault station (fault 1; `fault<N>st...`
                          for fault N >= 2 from our `faultstft<N>_` files)
  body<n>st<s>dp<d>       one per off-fault station
  cplot                   rupture time (j = along strike, k = down-dip
                          POSITIVE DOWN, t), `cplot_<N>` per fault when the
                          run has more than one
  PROVENANCE.md           code, version, commit, source stamp, dx, dt, and a
                          sha256 per packaged file. No paths, no hostnames.

No `.txt` extension: the SCEC upload form and the archive name files bare.

DATA ARE NEVER RE-FORMATTED. Every numeric row of a station file is copied
byte-for-byte from the run; the one repair is the Fortran e15.7 lost-'E'
exponent (`0.1341312-114` -> `0.1341312E-114`), which the CVWS server cannot
parse. Only the HEADER is rewritten, to the spec's field list (TPV26_27
Description v13 Part 8/9/10; the same list in every TPV spec since TPV29):
problem, author, date, code, code_version, commit, source_stamp,
element_size, time_step, num_time_steps, location, the column descriptions,
and the field-name line. cplot rows come from frt.canonical.txt (or the
deduplicated frt.txt<rank> set) through scec_readers.frt_to_strike_dip; a
never-ruptured node is written as 1.0E+09, the value the spec names
(EQdyna's frt carries 99999).

AUTHOR IS A PLACEHOLDER. The run's own headers carry a fictitious
`Project=San-Ti / Author=Sophon` pair; both are dropped and `author=` is
written as `<FILL IN>` for the submitter to set before upload. The script
never writes a personal name.

REFUSES rather than guesses when: the station headers disagree on
element_size or time_step, the measured dt disagrees with the header's,
the run's own banner names a version other than VERSION, the frt node
spacing disagrees with element_size, or a station name cannot be parsed.

NEVER MODIFIES the run directory.
"""
import argparse
import datetime
import glob
import hashlib
import os
import re
import subprocess
import sys
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, 'figures'))
import scec_readers as R      # noqa: E402

CODE_NAME = 'EQdyna'
AUTHOR_PLACEHOLDER = '<FILL IN>'
NEVER_RUPTURED = 1.0e9           # spec Part 10: "use the value 1.0E+09"

RE_FAULT = re.compile(r'^faultst(?:ft(\d+)_)?(-?\d+)dp(-?\d+)\.txt$')
RE_BODY = re.compile(r'^body(-?\d+)st(-?\d+)dp(-?\d+)\.txt$')
RE_BANNER_VERSION = re.compile(r'Welcome to EQdyna\s+([0-9][0-9A-Za-z.\-+]*)')
RE_BANNER_STAMP = re.compile(r'\(src\s+([0-9a-f]+)\)')
RE_DT = re.compile(r'time_step\s*=\s*([0-9.eEdD+-]+)')


# ------------------------------------------------------------- source files

def split_station(path):
    """-> (header_lines, field_line, data_lines). The field-name line is the
    first non-'#' line that does not parse as numbers; data lines are raw
    text (lost-'E' repaired), never re-formatted."""
    header, field, data = [], None, []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip('\n')
            if not s.strip():
                continue
            if field is None and data == []:
                if s.lstrip().startswith('#'):
                    header.append(s.strip())
                    continue
                try:
                    R._floats(s)
                except ValueError:
                    field = s.strip()
                    continue
            if s.lstrip().startswith('#'):
                continue
            fixed = R.RE_LOST_E.sub(r'\1E\2', s)
            data.append(fixed)
    if field is None:
        raise SystemExit(f'{path}: no field-name line (t h-slip ...)')
    if not data:
        raise SystemExit(f'{path}: no data rows')
    return header, field, data


def find_stations(run_dir):
    """Our station files live in the run dir, or in `stations/` for a
    committed reference set. -> list of (kind, fault, coords, path)."""
    dirs = [run_dir, os.path.join(run_dir, 'stations')]
    out, bad = [], []
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, 'faultst*.txt'))
                        + glob.glob(os.path.join(d, 'body*.txt'))):
            b = os.path.basename(f)
            m = RE_FAULT.match(b)
            if m:
                ft = int(m.group(1)) if m.group(1) else 1
                out.append(('fault', ft, (m.group(2), m.group(3)), f))
                continue
            m = RE_BODY.match(b)
            if m:
                out.append(('body', None, m.groups(), f))
                continue
            bad.append(b)
        if out:
            break
    if bad:
        raise SystemExit(f'unparseable station names: {bad}')
    if not out:
        raise SystemExit(f'{run_dir}: no faultst*.txt / body*.txt found')
    return out


def submission_name(kind, fault, coords):
    if kind == 'body':
        return 'body%sst%sdp%s' % coords
    pre = 'faultst' if fault == 1 else f'fault{fault}st'
    return '%s%sdp%s' % ((pre,) + coords)


def load_frt(run_dir):
    p = os.path.join(run_dir, 'frt.canonical.txt')
    if os.path.isfile(p):
        return R.read_numeric(p), 'frt.canonical.txt'
    ranks = sorted(glob.glob(os.path.join(run_dir, 'frt.txt*')))
    if not ranks:
        raise SystemExit(f'{run_dir}: no frt.canonical.txt or frt.txt<rank>')
    a = np.vstack([R.read_numeric(f) for f in ranks])
    # a node on a subdomain boundary is written by every rank that holds it
    _, keep = np.unique(np.round(a[:, :3], 3), axis=0, return_index=True)
    return a[np.sort(keep)], f'{len(ranks)} x frt.txt<rank> (deduplicated)'


# ------------------------------------------------------------- provenance

def repo_version():
    with open(os.path.join(REPO, 'VERSION')) as fh:
        return fh.read().strip()


def run_banner(run_dir):
    """(version, source stamp) from the run's own log banner, or (None, None)."""
    for name in ('run.log', 'eqdyna.log'):
        p = os.path.join(run_dir, name)
        if os.path.isfile(p):
            with open(p, errors='replace') as fh:
                txt = fh.read(20000)
            v = RE_BANNER_VERSION.search(txt)
            s = RE_BANNER_STAMP.search(txt)
            return (v.group(1) if v else None), (s.group(1) if s else None)
    return None, None


def git_head(path):
    try:
        return subprocess.check_output(
            ['git', '-C', path, 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL, text=True).strip()
    except (subprocess.CalledProcessError, OSError):
        return None


def one_value(values, what):
    u = sorted(set(values))
    if len(u) != 1:
        raise SystemExit(f'station headers disagree on {what}: {u}')
    return u[0]


# ------------------------------------------------------------- writers

def station_text(meta, hdr, field, data, kind):
    location = [h for h in hdr if re.match(r'#\s*location', h, re.I)]
    columns = [h for h in hdr if re.match(r'#\s*Column\s*#', h, re.I)]
    if kind == 'body':
        columns = _dedupe_body_columns(columns)
    ncol = len(data[0].split())
    lines = [f'# problem={meta["problem"]}',
             f'# author={AUTHOR_PLACEHOLDER}',
             f'# date={meta["date"]}',
             f'# code={CODE_NAME}',
             f'# code_version={meta["version"]}',
             f'# commit={meta["commit"]}']
    if meta.get('stamp'):
        lines.append(f'# source_stamp={meta["stamp"]}')
    lines += [f'# element_size={meta["dx"]:g} m',
              f'# time_step={meta["dt"]:.12e}',
              f'# num_time_steps={len(data)}']
    lines += [h.replace('# location =', '# location=') for h in location]
    lines.append(f'# Time series in {ncol} columns')
    lines += columns
    lines += ['#', '# The line below lists the names of the data fields:', field]
    return '\n'.join(lines + data) + '\n'


def _dedupe_body_columns(cols):
    """The archive-era body header listed `Column #3` twice; today's writer
    does not, but a header is never trusted to be right twice."""
    seen, out = set(), []
    for c in cols:
        k = re.match(r'#\s*Column\s*#\s*(\d+)', c).group(1)
        if k not in seen:
            seen.add(k)
            out.append(c)
    return out


def cplot_text(meta, strike, downdip, t, fault, nft):
    lines = [f'# problem={meta["problem"]}',
             f'# author={AUTHOR_PLACEHOLDER}',
             f'# date={meta["date"]}',
             f'# code={CODE_NAME}',
             f'# code_version={meta["version"]}',
             f'# commit={meta["commit"]}']
    if meta.get('stamp'):
        lines.append(f'# source_stamp={meta["stamp"]}')
    lines.append(f'# element_size={meta["dx"]:g} m')
    if nft > 1:
        lines.append(f'# fault={fault}')
    lines += ['# Column #1 = Distance along strike (m)',
              '# Column #2 = Distance down-dip (m), positive underground',
              '# Column #3 = Rupture time (s); 1.0E+09 = never ruptured',
              '#', 'j k t']
    rows = ['%15.7E %15.7E %15.7E' % r for r in zip(strike, downdip, t)]
    return '\n'.join(lines + rows) + '\n'


# ------------------------------------------------------------- main

def package(run_dir, tpv, out, commit=None, date=None, label=None):
    run_dir = os.path.abspath(run_dir)
    stations = find_stations(run_dir)
    parsed = [(k, ft, c, p, split_station(p)) for k, ft, c, p in stations]

    # dx and dt: one value across every station header, dt also MEASURED
    dxs = [R.header_element_size_m(p) for *_, p, _ in parsed]
    if any(v is None for v in dxs):
        raise SystemExit('a station header has no element_size')
    dx = one_value([round(v, 6) for v in dxs], 'element_size')
    dt_meas = []
    for *_, p, (hdr, _, data) in parsed:
        tcol = np.array([float(r.split()[0]) for r in data])
        if tcol.size > 1:
            d = np.diff(tcol)
            if np.ptp(d) > 1e-6 * np.median(d):
                raise SystemExit(f'{p}: unequal time steps')
            dt_meas.append(float(np.median(d)))
        hdt = [float(m.group(1)) for h in hdr for m in [RE_DT.search(h)] if m]
        if hdt and dt_meas and abs(hdt[0] - dt_meas[-1]) > 0.51e-4:
            raise SystemExit(f'{p}: header time_step {hdt[0]} != measured '
                             f'{dt_meas[-1]:.6g}')
    dt = float(np.median(dt_meas))
    if np.ptp(dt_meas) > 1e-9 * dt:
        raise SystemExit(f'stations disagree on dt: {sorted(set(dt_meas))}')

    version = repo_version()
    run_version, stamp = run_banner(run_dir)
    if run_version is not None and run_version != version:
        raise SystemExit(f'run banner says EQdyna {run_version}, VERSION says '
                         f'{version}: package from a checkout at the run\'s version')
    commit = commit or git_head(run_dir)
    if not commit:
        raise SystemExit('no commit SHA: pass --commit')

    # rupture time, per fault
    frt, frt_src = load_frt(run_dir)
    faults = sorted({ft for k, ft, *_ in stations if k == 'fault'}) or [1]
    nft = len(faults)
    cplots = []
    for ft in faults:
        strike, dd, t, rule = R._read_frt_array(frt, ft if nft > 1 else None)
        sent, _ = R.detect_sentinel(t)
        t = np.asarray(t, float).copy()
        if sent is not None:
            t[t >= sent] = NEVER_RUPTURED
        measured = R.measure_dx(strike, dd)
        if abs(measured - dx) > 1e-3 * dx:
            raise SystemExit(f'frt node spacing {measured:g} m != station '
                             f'element_size {dx:g} m')
        cplots.append((ft, strike, dd, t, rule, int((t >= NEVER_RUPTURED).sum())))

    date = date or datetime.date.today().strftime('%Y/%m/%d')
    year = date.split('/')[0]
    label = label or f'eqdyna-v{version}-{dx:g}m-{year}'
    meta = dict(problem=f'TPV{tpv}', date=date, version=version,
                commit=commit, stamp=stamp, dx=dx, dt=dt)

    members = []
    for kind, ft, coords, path, (hdr, field, data) in parsed:
        members.append((submission_name(kind, ft, coords),
                        station_text(meta, hdr, field, data, kind), path))
    for ft, strike, dd, t, _, _ in cplots:
        name = 'cplot' if nft == 1 else f'cplot_{ft}'
        members.append((name, cplot_text(meta, strike, dd, t, ft, nft), None))
    names = [m[0] for m in members]
    if len(set(names)) != len(names):
        raise SystemExit('two source files map to one submission name')

    prov = [f'# TPV{tpv} -- {CODE_NAME} v{version}, {dx:g} m, packaged {date}', '',
            f'- code: {CODE_NAME}', f'- code_version: {version} (VERSION)',
            f'- commit: {commit}',
            f'- source_stamp: {stamp or "(no run banner)"} (scripts/src_hash.py, from the run banner)',
            f'- element_size: {dx:g} m', f'- time_step: {dt:.12e} s',
            f'- rupture time from: {frt_src}',
            f'- author: {AUTHOR_PLACEHOLDER} -- set in every file header before upload',
            f'- files: {sum(1 for n in names if n.startswith("fault"))} on-fault, '
            f'{sum(1 for n in names if n.startswith("body"))} off-fault, '
            f'{len(cplots)} cplot', '',
            'sha256 of each packaged file:', '']
    prov += [f'    {hashlib.sha256(txt.encode()).hexdigest()}  {n}'
             for n, txt, _ in members]
    members.append(('PROVENANCE.md', '\n'.join(prov) + '\n', None))

    zdt = tuple(int(v) for v in date.split('/')) + (0, 0, 0)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or '.', exist_ok=True)
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, txt, _ in sorted(members):
            zi = zipfile.ZipInfo(f'tpv{tpv}/{label}/{name}', date_time=zdt)
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, txt)
    return dict(out=out, label=label, members=names, dx=dx, dt=dt,
                version=version, commit=commit, stamp=stamp, frt=frt_src,
                unruptured=[c[5] for c in cplots], rules=[c[4] for c in cplots])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run_dir')
    ap.add_argument('--tpv', required=True, help='benchmark number, e.g. 26')
    ap.add_argument('--out', required=True, help='zip path to write')
    ap.add_argument('--commit', default=None,
                    help='commit SHA of the build (default: git HEAD of RUN_DIR)')
    ap.add_argument('--date', default=None, help='YYYY/MM/DD (default today)')
    ap.add_argument('--label', default=None,
                    help='directory name in the zip (default eqdyna-v<ver>-<dx>m-<year>)')
    a = ap.parse_args(argv)
    r = package(a.run_dir, str(a.tpv).lower().lstrip('tpv'), a.out,
                a.commit, a.date, a.label)
    n_f = sum(1 for n in r['members'] if n.startswith('fault'))
    n_b = sum(1 for n in r['members'] if n.startswith('body'))
    print(f'wrote {r["out"]}: {len(r["members"]) + 1} files '
          f'({n_f} on-fault, {n_b} off-fault, '
          f'{len(r["members"]) - n_f - n_b} cplot, PROVENANCE.md) under {r["label"]}/')
    print(f'  EQdyna {r["version"]}  commit {r["commit"]}  stamp {r["stamp"]}')
    print(f'  dx {r["dx"]:g} m  dt {r["dt"]:.12e} s  rupture from {r["frt"]}; '
          f'never-ruptured nodes {r["unruptured"]}; {"; ".join(r["rules"])}')
    print(f'  author field is "{AUTHOR_PLACEHOLDER}" in every header: fill it before upload')
    return 0


if __name__ == '__main__':
    sys.exit(main())
