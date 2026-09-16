#! /usr/bin/env python3
"""
Lay the fetched raw files out in scec_archive/<benchmark>/<ver>-<res>-<year>/
per scec_archive/README.md, and write one PROVENANCE.md per directory.

Version / resolution / date come from the FILE HEADERS, never from guesses.
A benchmark whose cplot carries no header (tpv104, tpv105-3d) takes them from
that submission's own time-series headers, which do carry them.
"""
import html, json, os, re, shutil, sys, hashlib, glob

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, 'raw')
ARCHIVE = os.path.abspath(os.path.join(HERE, '..', '..', 'scec_archive'))
PAGES = os.path.join(HERE, 'pages')
FETCH_DATE = '2026-09-14'
CHECKSUMS = os.path.join(HERE, 'CHECKSUMS.sha256')
CGI = 'https://strike.scec.org/cvws/cgi-bin/cvws.cgi'

HDR = re.compile(r'^#\s*([A-Za-z_]+)\s*=\s*(.*?)\s*$')


def headers(path):
    """Header key/value pairs; tolerant of the one-space indent some
    EQdyna writers put in front of every comment line."""
    out = {}
    with open(path, errors='replace') as f:
        for n, line in enumerate(f):
            s = line.lstrip()
            if not s.startswith('#'):
                if out or n > 30:
                    break
                continue
            m = HDR.match(s.rstrip('\n'))
            if m:
                out.setdefault(m.group(1).lower(), m.group(2))
    return out


def _cand(path, fname):
    h = headers(path)
    if not h:
        return None
    es = h.get('element_size', '')
    m = re.search(r'([\d.]+)', es)
    return dict(code=h.get('code', ''), version=h.get('code_version', ''),
                element_size=float(m.group(1)) if m else None,
                element_size_raw=es, date=h.get('date', ''),
                author=h.get('author', ''),
                problem=h.get('problem', '') or h.get('project', ''),
                first_line=open(path, errors='replace').readline().rstrip('\n'),
                src=fname, time_step=h.get('time_step', ''),
                num_time_steps=h.get('num_time_steps', ''))


def label_res(label):
    """Resolution (m) named in the public submission label, if any."""
    m = re.search(r'(\d+(?:\.\d+)?)\s*m\b', label)
    return float(m.group(1)) if m else None


def label_ver(label):
    m = re.search(r'[Vv]?(\d+\.\d+(?:\.\d+)?)', label.replace('3DMPI', ''))
    return m.group(1) if m else ''


def meta(bm, user, label=''):
    """Metadata from the FILE HEADERS.

    cplot and the time series are written by different parts of the code and
    occasionally disagree; when they do, the one whose element_size matches
    the resolution named in the public submission label wins, and the
    disagreement is reported so PROVENANCE.md can record it.
    """
    d = os.path.join(RAW, bm, user)
    cands = {}
    cp = os.path.join(d, 'cplot')
    if os.path.exists(cp):
        cands['cplot'] = _cand(cp, 'cplot')
    for f in sorted(x for x in os.listdir(d) if x != 'cplot'):
        c = _cand(os.path.join(d, f), f)
        if c:
            cands['timeseries'] = c
            break
    cands = {k: v for k, v in cands.items() if v}
    if not cands:
        return None, []
    want = label_res(label)
    order = ['cplot', 'timeseries']
    if want is not None:
        order.sort(key=lambda k: (k not in cands or
                                  cands[k]['element_size'] != want))
    pick = next(k for k in order if k in cands)
    m = dict(cands[pick])
    m['picked_from'] = pick
    # cplot carries a single hand-entered date that is sometimes a leftover
    # from another run (tpv30/dliu still says "LVFZ / Benchun Duan / 6-2-2012").
    # Every time-series file stamps its own write time, so the year comes from
    # there whenever a time series exists.
    m['date_for_year'] = (cands['timeseries']['date'] if 'timeseries' in cands
                          and cands['timeseries']['date'] else m['date'])
    m['cplot_date'] = cands.get('cplot', {}).get('date', '')
    m['ts_date'] = cands.get('timeseries', {}).get('date', '')
    # disagreements worth recording
    diffs = []
    if len(cands) == 2:
        a, b = cands['cplot'], cands['timeseries']
        for k in ('code', 'version', 'element_size', 'date', 'author',
                  'problem'):
            if (a[k] or b[k]) and a[k] != b[k]:
                diffs.append((k, a[k], b[k]))
    if not m['version'] and label:
        m['version'] = label_ver(label)
        m['version_from'] = 'public submission label (no code_version in any '
        'file header)'
    return m, diffs


def dx_from_cplot(path):
    """Fallback resolution: smallest positive along-strike step in cplot."""
    xs = []
    for line in open(path, errors='replace'):
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) >= 3:
            try:
                xs.append(float(p[0]))
            except ValueError:
                continue
        if len(xs) > 500:
            break
    d = sorted({round(abs(b - a), 3) for a, b in zip(xs, xs[1:]) if b != a})
    return d[0] if d else None


def year_of(datestr):
    """(year, ISO date) from the several date spellings the headers use:
    ' 2/27/2015 11:17:30', '2016/01/30', '6-2-2012', '11/ 6/2024 ...'."""
    s = datestr.strip()
    m = re.match(r'(\d{4})[/-](\d{1,2})[/-](\d{1,2})', s)
    if m:
        return m.group(1), (f'{m.group(1)}-{int(m.group(2)):02d}-'
                            f'{int(m.group(3)):02d}')
    m = re.match(r'\s*(\d{1,2})\s*[/-]\s*(\d{1,2})\s*[/-]\s*(\d{4})', s)
    if m:
        return m.group(3), (f'{m.group(3)}-{int(m.group(1)):02d}-'
                            f'{int(m.group(2)):02d}')
    m = re.search(r'(\d{4})', s)
    return (m.group(1), s) if m else (None, s)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def verify_checksums(checksums_path=CHECKSUMS, archive=ARCHIVE):
    """Check every file scec_archive/ is supposed to hold against
    CHECKSUMS.sha256 (standard `sha256sum` format: '<hex>  ./relpath').
    scec_archive/ itself is gitignored (484 MB, an owner-retained asset, not
    reproducible from this repo per its own README) -- this 80 KB manifest is
    the tracked, auditable stand-in: it proves the on-disk archive has not
    silently drifted, without putting the archive itself in git (the mistake
    f2c9851 made via an unmatchable .gitignore inline-comment pattern).

    Returns (ok, mismatched, missing, extra) counts/lists; prints a report;
    does not raise (rule 2 -- the CALLER decides the exit code, same as every
    other entry point in this file)."""
    if not os.path.isfile(checksums_path):
        print(f'FAIL: no checksum manifest at {checksums_path}', file=sys.stderr)
        return None
    expected = {}
    with open(checksums_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            digest, relpath = line.split('  ', 1)
            expected[relpath[2:] if relpath.startswith('./') else relpath] = digest

    ok, mismatched, missing = [], [], []
    for relpath, digest in sorted(expected.items()):
        full = os.path.join(archive, relpath)
        if not os.path.isfile(full):
            missing.append(relpath)
            continue
        actual = sha(full)
        (ok if actual == digest else mismatched).append(relpath)

    on_disk = set()
    for root, _, files in os.walk(archive):
        for fn in files:
            on_disk.add(os.path.relpath(os.path.join(root, fn), archive))
    extra = sorted(on_disk - set(expected))

    print(f'CHECKSUMS.sha256: {len(expected)} entries, archive has '
          f'{len(on_disk)} files')
    print(f'  OK        : {len(ok)}')
    print(f'  MISMATCH  : {len(mismatched)}')
    for r in mismatched:
        print(f'    {r}')
    print(f'  MISSING   : {len(missing)}  (listed in manifest, not on disk)')
    for r in missing:
        print(f'    {r}')
    print(f'  EXTRA     : {len(extra)}  (on disk, not in manifest)')
    for r in extra:
        print(f'    {r}')
    return ok, mismatched, missing, extra


# ---------------------------------------------------------------- metrics
def table_rows(doc):
    doc = re.sub(r'(?is)<(script|style).*?</\1>', '', doc)
    out = []
    for tr in re.findall(r'(?is)<tr[^>]*>(.*?)</tr>', doc):
        cs = [re.sub(r'\s+', ' ', html.unescape(re.sub(r'<[^>]+>', ' ', c))).strip()
              for c in re.findall(r'(?is)<t[dh][^>]*>(.*?)</t[dh]>', tr)]
        if cs:
            out.append(cs)
    return out


ID = re.compile(r'^\((\d+)\)\s+(\S+)$')


def rms_matrix(path):
    """(names, {name: {other: ms}}) from a metric_cvv1_<bm>_ac_0.html page."""
    rows = table_rows(open(path, errors='replace').read())
    names, legend, data = [], {}, {}
    for r in rows:
        m = ID.match(r[0])
        if not m:
            continue
        n = m.group(2)
        vals = r[1:]
        if all(v == '' or re.fullmatch(r'-?[\d.]+', v) for v in vals) and \
                any(re.fullmatch(r'-?[\d.]+', v) for v in vals):
            if n not in data:
                names.append(n)
                data[n] = vals
        else:
            legend.setdefault(n, r[1])
    out = {}
    for n in names:
        # the diagonal cell is PRESENT but empty, so vals lines up with names
        out[n] = {o: v for o, v in zip(names, data[n]) if v}
    # symmetry check -- the published matrix is symmetric
    for a in names:
        for b, v in out[a].items():
            assert out[b][a] == v, f'asymmetric {a}/{b}: {v} vs {out[b][a]}'
    return names, legend, out


# ---------------------------------------------------------------- surfdef
def surfdef_table(bm):
    p = os.path.join(PAGES, f'{bm}_surfdef.html')
    if not os.path.exists(p):
        return None
    out = {}
    for r in table_rows(open(p, errors='replace').read()):
        if len(r) >= 3 and re.fullmatch(r'[\d.]+', r[1] or 'x'):
            out[r[0]] = (r[1], r[2])
    return out


# ------------------------------------------------- served-vs-uploaded text
#
# The CGI prepends ONE space to lines that start with '#' AT COLUMN 0 and
# leaves every other line alone.  fetch_all_dliu.unrender() reverses that by
# stripping one space from every line that starts with ' #' -- which is right
# only when the upload's comment lines really were at column 0.
#
# EQdyna's own Fortran writer (src/library_output.f90) emits comment lines
# with list-directed WRITE(51,*), which always prepends one blank, and writes
# the time_step / num_time_steps lines under explicit formats:
#     on-fault   '( a14,...)'  '# time_step ='      13 chars -> 1 pad blank
#     off-fault  '( a14,...)'  '# time_step='       12 chars -> 2 pad blanks
#     on-fault   '( a19,...)'  '# num_time_steps =' 18 chars -> 1 pad blank
#     off-fault  '( a19,...)'  '# num_time_steps='  17 chars -> 2 pad blanks
# So a submission written by that code MUST show a two-space '  # time_step='
# line in its off-fault files -- and it does, served exactly two spaces wide,
# proving the CGI prepended nothing to it.  Its sibling comment lines are
# therefore already one space in, and unrender() over-strips them.
#
# Detection is per submission, on that unambiguous two-space ruler.
RULER = re.compile(r'^  #\s*(time_step|num_time_steps)=')


def indenting_writer(srcdir):
    """True when this submission's time series came from the indenting writer."""
    for f in sorted(os.listdir(srcdir)):
        if f == 'cplot':
            continue
        with open(os.path.join(srcdir, f), errors='replace') as fh:
            for i, line in enumerate(fh):
                if i > 25:
                    break
                if RULER.match(line):
                    return True
    return False


def restore(path):
    """Put back the one space unrender() should not have taken."""
    txt = open(path, errors='replace').read()
    if not any(l.startswith('#') for l in txt.split('\n')):
        return False
    out = '\n'.join(' ' + l if l.startswith('#') else l
                    for l in txt.split('\n'))
    open(path, 'w').write(out)
    return True


# ---------------------------------------------------------------- layout
def res_tag(size):
    if size is None:
        return 'unkres'
    return f'{int(round(size))}m' if abs(size - round(size)) < 1e-6 \
        else f'{size:g}m'


def dirname(m, dxfall):
    code = (m['code'] or 'eqdyna').lower().replace(' ', '')
    ver = 'v' + m['version'].strip()
    size = m['element_size'] if m['element_size'] else dxfall
    yr = year_of(m.get('date_for_year') or m['date'])[0] or 'unkyear'
    return f'{code}-{ver}-{res_tag(size)}-{yr}'


def main():
    manifest = json.load(open(os.path.join(HERE, 'manifest.json')))
    idx = []
    for e in manifest:
        bm, user = e['benchmark'], e['user']
        src = os.path.join(RAW, bm, user)
        m, diffs = meta(bm, user, e['label'])
        cplot = os.path.join(src, 'cplot')
        dxf = dx_from_cplot(cplot) if os.path.exists(cplot) else None
        if m is None:
            print(f'!! {bm}/{user}: no headers at all', file=sys.stderr)
            continue
        if m['element_size'] is None:
            m['element_size'] = dxf
        dn = dirname(m, dxf)
        dst = os.path.join(ARCHIVE, bm, dn)
        os.makedirs(dst, exist_ok=True)
        files = sorted(os.listdir(src))
        indents = indenting_writer(src)
        sizes, preexisting = {}, {}
        for f in files:
            s, d = os.path.join(src, f), os.path.join(dst, f)
            if indents and f != 'cplot':
                restore(s)          # undo the over-eager unrender, in place
            if os.path.exists(d):
                # rule 7: archival records are read-only -- never overwrite.
                # Compare instead and record the verdict.
                preexisting[f] = ('identical' if sha(s) == sha(d)
                                  else 'DIFFERS -- kept the existing file')
            else:
                shutil.copy2(s, d)
            sizes[f] = os.path.getsize(d)
        idx.append(dict(benchmark=bm, user=user, label=e['label'], dir=dn,
                        meta=m, files=files, sizes=sizes,
                        preexisting=preexisting, header_diffs=diffs,
                        indenting_writer=indents,
                        failed=e.get('failed', {}), dx_from_cplot=dxf,
                        dst=dst))
    # the TPV36/37 surface-deformation files are not served by the CGI; they
    # were e-mailed and are published as zips linked from tpv3*_surfdef.html.
    # They carry the portal user id of the *_arc listing, so they attach to
    # the directory whose (benchmark, user) matches that id.
    for bm, user, zf in [
            ('tpv36', 'dliu.3',
             'TPV36_EQdyna_v5_3_3_50m_dliu_SCECFinalSurfDisp.zip'),
            ('tpv37_arc', 'dliu.3',
             'TPV37_EQdyna_v5_3_3_50m_dliu_SCECFinalSurfDisp.zip')]:
        src = os.path.join(HERE, 'surfdef', zf)
        for e in idx:
            if e['benchmark'] == bm and e['user'] == user:
                d = os.path.join(e['dst'], zf)
                if not os.path.exists(d):
                    shutil.copy2(src, d)
                e['surfdef_zip'] = zf
                e['files'].append(zf)
                e['sizes'][zf] = os.path.getsize(d)
    json.dump(idx, open(os.path.join(HERE, 'layout.json'), 'w'), indent=1,
              default=str)
    print(f'{len(idx)} submissions laid out under {ARCHIVE}')


if __name__ == '__main__':
    if '--verify' in sys.argv:
        result = verify_checksums()
        sys.exit(1 if result is None or any(result[1:]) else 0)
    main()
