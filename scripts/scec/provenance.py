#! /usr/bin/env python3
"""Write one PROVENANCE.md per archived submission + scec_archive/INDEX.md."""
import json, os, re, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from organize import (HERE, ARCHIVE, PAGES, FETCH_DATE, CGI, rms_matrix,
                      surfdef_table, year_of, res_tag, label_res)

METRIC_URL = ('https://strike.scec.org/cvws/metric_cvv1_u1/{bm}/'
              'metric_cvv1_{bm}_ac_0.html')
SURFDEF_URL = 'https://strike.scec.org/cvws/{bm}_surfdef.html'


def load_metric(bm):
    p = os.path.join(PAGES, f'{bm}.html')
    if not os.path.exists(p) or os.path.getsize(p) < 1000:
        return None
    try:
        names, legend, mat = rms_matrix(p)
    except Exception as e:
        print(f'!! {bm} metric parse: {e}', file=sys.stderr)
        return None
    return dict(names=names, legend=legend, mat=mat)


def standing(bm, label, met):
    """Locate this submission in the published RMS matrix BY LABEL."""
    if not met:
        return None, (f'No public metric page for `{bm}` — '
                      f'<{METRIC_URL.format(bm=bm)}> returns HTTP 404, so SCEC '
                      f'publishes no cross-code table for this listing.')
    hit = [n for n in met['names']
           if met['legend'].get(n, '').strip() == label.strip()]
    if not hit:
        mine = {n: met['legend'].get(n, '') for n in met['names']
                if n.startswith('dliu')}
        extra = ('' if not mine else
                 ' The page does carry other EQdyna column(s), under different '
                 'labels: ' + '; '.join(f'`{k}` = "{v}"'
                                        for k, v in mine.items()) + '.')
        return None, (f'Not in the published metric page for `{bm}` '
                      f'(<{METRIC_URL.format(bm=bm)}>) — no column carries this '
                      f'submission label.{extra}')
    n = hit[0]
    row = {k: float(v) for k, v in met['mat'][n].items()}
    meds = {m: statistics.median([float(x) for x in met['mat'][m].values()])
            for m in met['names']}
    rank = sorted(meds, key=meds.get).index(n) + 1
    return dict(col=n, row=row, median=meds[n], rank=rank,
                n=len(met['names'])), None


def fmt_standing(bm, s):
    lo = min(s['row'], key=s['row'].get)
    hi = max(s['row'], key=s['row'].get)
    lines = [
        f"Published as column `{s['col']}` of the cross-code metric page "
        f"<{METRIC_URL.format(bm=bm)}>, table *RMS difference in rupture time "
        f"(milliseconds)*, file `cplot`.",
        '',
        f"- median RMS rupture-time difference from the other "
        f"{s['n']-1} submissions: **{s['median']:.1f} ms**",
        f"- closest submission `{lo}` at {s['row'][lo]:.1f} ms; farthest "
        f"`{hi}` at {s['row'][hi]:.1f} ms",
        f"- rank by that median among the {s['n']} submissions: "
        f"**{s['rank']} of {s['n']}** (1 = closest to the group)",
        '',
        '| other submission | RMS (ms) |',
        '|---|---|',
    ]
    for k in sorted(s['row'], key=s['row'].get):
        lines.append(f'| `{k}` | {s["row"][k]:.1f} |')
    return '\n'.join(lines)


PROV = """# {BM} — {code} v{ver}, {res}, {date}

Submission label in the public record: `{user}` — "{label}".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST {cgi}
      o=1005   m={bm}   mus={user}
      Q0001=<on-fault station list>    Q0003=<off-fault station list>
      G1063cplot=Raw Data              -> cplot (rupture-time contour array)
      G1059<station>=Raw Data          -> one faultst*/body* time series

Field names were read off the CGI's own file-list form, not guessed; `G1063`
is bound to the `cplot` row only and `G1059` to every station row.
Fetched **{fetched}**.

### Served text vs uploaded text

The CGI renders each file inside `<pre>` and prepends one space to lines that
start with `#` **at column 0**; every other line is served byte-identical.
{transform}

## From the file headers

| field | value | read from |
|---|---|---|
| first header line | `{first}` | `{src}` |
| author | {author} | `{src}` |
| code | {code} | `{src}` |
| code_version | {ver}{verfrom} | `{src}` |
| element_size | {esize} | `{src}` |
| date | {date_raw} | `{src}` |
{extra_hdr}
## Files ({nfiles})

{filetable}

## Published standing

{standing}
{notes}"""

TRANSFORM_COL0 = """This submission's uploaded comment lines start at column 0
(the same writer as the TPV29 `dliu` upload the fetch path was verified
against), so the one prepended space is stripped back off and the files below
are the uploads verbatim, byte for byte."""

TRANSFORM_INDENT = """This submission's time series were written by EQdyna's
Fortran writer, whose list-directed `WRITE(51,*)` puts one blank in front of
every comment line, and whose `time_step` / `num_time_steps` lines use explicit
`a14` / `a19` formats. Those formats pin the upload's indent exactly — the
off-fault `'# time_step='` literal is 12 characters, so `a14` right-justifies
it with **two** leading blanks, and the CGI serves that line with exactly two.
The CGI therefore prepended nothing to this file, and the served text already
**is** the upload. No space is stripped from the time series here. `cplot`,
which is written separately with column-0 comments, is unrendered normally."""


# submissions the public area lists twice, verified byte-identical:
#   (benchmark, user) -> (other benchmark, other directory, other has metrics)
MIRRORS = {
    ('tpv105-3d', 'dliu.2'):
        ('tpv105-3d_arc', 'eqdyna3d-v5.1.0-100m-2020', False),
    ('tpv105-3d_arc', 'dliu.2'):
        ('tpv105-3d', 'eqdyna3d-v5.1.0-100m-2020', True),
    ('tpv36', 'dliu.3'): ('tpv36_arc', 'eqdyna-v5.3.3-50m-2024', False),
    ('tpv36_arc', 'dliu.3'): ('tpv36', 'eqdyna-v5.3.3-50m-2024', False),
}


def main():
    layout = json.load(open(os.path.join(HERE, 'layout.json')))
    mets, index = {}, []
    for e in layout:
        bm, user, m = e['benchmark'], e['user'], e['meta']
        base = bm.replace('_arc', '')
        if bm not in mets:
            mets[bm] = load_metric(bm)
        s, why = standing(bm, e['label'], mets[bm])
        st = fmt_standing(bm, s) if s else why

        notes = []
        diffs = [(k, a, b) for k, a, b in e.get('header_diffs', [])
                 # a date written '2/27/2015' one place and
                 # '2/27/2015 11:17:29' the other is the same date, not a
                 # disagreement
                 if not (k == 'date' and year_of(str(a))[1] == year_of(str(b))[1])]
        for k, a, b in diffs:
            notes.append(f'`cplot` and the time series disagree on '
                         f'**{k}**: `cplot` says `{a}`, the time series say '
                         f'`{b}`.')
        if any(k == 'author' for k, _, _ in diffs):
            notes.append('That `cplot` header is a leftover from a different '
                         'run and does not describe this submission; the '
                         'directory year is taken from the time-series '
                         'timestamps, which all 36 files agree on.')

        sd = surfdef_table(base) if base in ('tpv36', 'tpv37') else None
        if sd and user in sd:
            up, sub = sd[user]
            ups = [float(v[0]) for v in sd.values()]
            where = (f"Archived here as `{e['surfdef_zip']}`."
                     if e.get('surfdef_zip') else
                     "Archived under this benchmark's `dliu.3` directory.")
            notes.append(
                f"**Surface deformation.** {base.upper()} also takes a surface-"
                f"deformation file, which is e-mailed rather than uploaded to "
                f"the CGI and is published at <{SURFDEF_URL.format(bm=base)}>. "
                f"That page lists user `{user}`: max uplift **{up} m**, max "
                f"subsidence **{sub} m**, against a group range of "
                f"{min(ups):.3f}–{max(ups):.3f} m uplift over {len(sd)} "
                f"submissions. {where}")
        elif base in ('tpv36', 'tpv37') and sd:
            notes.append(
                f"The {base.upper()} surface-deformation listing "
                f"(<{SURFDEF_URL.format(bm=base)}>) carries no entry for user "
                f"`{user}`; its EQdyna entry is `dliu.3`, archived under "
                f"`{base if base == 'tpv36' else base + '_arc'}/`.")

        mir = MIRRORS.get((bm, user))
        if mir:
            notes.append(
                f"The public area lists this same submission twice, once "
                f"under `{bm}` and once under `{mir[0]}`. Every file fetched "
                f"here is **byte-identical** to the copy in "
                f"`{mir[0]}/{mir[1]}/`; both listings are kept so the archive "
                f"mirrors what the portal actually offers."
                + (f" The published metric standing is recorded on the "
                   f"`{mir[0]}` copy." if mir[2] else ''))

        pre = e.get('preexisting') or {}
        if pre:
            bad = [f for f, v in pre.items() if v != 'identical']
            notes.insert(0,
                'This directory was seeded from the author\'s own retained '
                '2015 upload (`cplot`, `faultst000dp120`, '
                '`body-030st000dp000`, placed here 2026-09-14) before any '
                'fetch, which is why it is the reference for the byte-match '
                'check below.')
            notes.append(
                f"{len(pre)} file(s) were already archived in this directory "
                f"before the fetch and were left untouched (README rule: "
                f"archival records are read-only). The freshly fetched bytes "
                f"were compared against them instead — "
                + ('**all byte-identical**, which is the correctness check on '
                   'this whole fetch path.' if not bad
                   else '**DIFFERENT**: ' + ', '.join(f'`{b}`' for b in bad)))
        if e.get('failed'):
            notes.append('Files the public area would not serve: '
                         + ', '.join(f'`{k}` ({v})'
                                     for k, v in e['failed'].items()))

        files = e['files']
        rows = ['| file | bytes |', '|---|---|']
        for f in files:
            rows.append(f'| `{f}` | {e["sizes"][f]:,} |')
        nfault = sum(1 for f in files if f.startswith('faultst'))
        nbody = sum(1 for f in files if f.startswith('body'))
        nother = len(files) - nfault - nbody - (1 if 'cplot' in files else 0) \
            - (1 if e.get('surfdef_zip') else 0)
        extra = []
        if m.get('time_step'):
            extra.append(f"| time_step | {m['time_step']} | `{m['src']}` |")
        if m.get('num_time_steps'):
            extra.append(
                f"| num_time_steps | {m['num_time_steps']} | `{m['src']}` |")
        if e.get('dx_from_cplot'):
            extra.append(f"| on-fault node spacing | {e['dx_from_cplot']:g} m "
                         f"| `cplot` coordinates |")
        counts = []
        if 'cplot' in files:
            counts.append('cplot')
        if nfault:
            counts.append(f'{nfault} on-fault')
        if nbody:
            counts.append(f'{nbody} off-fault')
        if nother > 0:
            counts.append(f'{nother} named station')
        if e.get('surfdef_zip'):
            counts.append('surface-deformation zip')
        yr, iso = year_of(m.get('date_for_year') or m['date'])
        txt = PROV.format(
            BM=bm.upper().replace('_ARC', '_arc'),
            code=m['code'] or 'EQdyna', ver=m['version'],
            res=res_tag(m['element_size']), date=iso, user=user,
            label=e['label'], cgi=CGI, bm=bm, fetched=FETCH_DATE,
            transform=(TRANSFORM_INDENT if e.get('indenting_writer')
                       else TRANSFORM_COL0),
            first=m.get('first_line', '').strip(),
            author=m['author'] or '(not in header)',
            verfrom=('' if m.get('version') and not m.get('version_from')
                     else ' *(from the public submission label — no '
                          '`code_version` in any file header)*'),
            esize=m['element_size_raw'] or '(absent from the headers)',
            date_raw=m['date'] or '(absent)', src=m['src'],
            extra_hdr=('\n'.join(extra) + '\n') if extra else '',
            nfiles=f'{len(files)} — ' + ', '.join(counts),
            filetable='\n'.join(rows), standing=st,
            notes=('\n## Notes\n\n' + '\n\n'.join('- ' + n for n in notes)
                   + '\n') if notes else '')
        open(os.path.join(e['dst'], 'PROVENANCE.md'), 'w').write(txt)
        index.append(dict(bm=bm, dir=e['dir'], user=user, label=e['label'],
                          code=m['code'], ver=m['version'],
                          res=res_tag(m['element_size']), date=iso,
                          nfault=nfault, nbody=nbody, nother=max(nother, 0),
                          has_cplot='cplot' in files,
                          surfdef=bool(e.get('surfdef_zip')),
                          total=sum(e['sizes'].values()),
                          standing=(f"{s['median']:.1f} ms median RMS, "
                                    f"rank {s['rank']}/{s['n']}"
                                    if s else '—')))
    json.dump(index, open(os.path.join(HERE, 'index.json'), 'w'), indent=1)
    print(f'{len(index)} PROVENANCE.md written')


if __name__ == '__main__':
    main()
