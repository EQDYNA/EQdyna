#! /usr/bin/env python3
"""Write scec_archive/INDEX.md from index.json + enumeration.json."""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from organize import HERE, ARCHIVE, FETCH_DATE

idx = json.load(open(os.path.join(HERE, 'index.json')))
idx.sort(key=lambda e: (e['bm'].replace('_arc', 'z_arc'), e['res']))

L = []
w = L.append
w('# INDEX — EQdyna submissions to the SCEC/USGS code-verification project')
w('')
w(f'Everything archived here, fetched from the public SCEC cvws portal on '
  f'**{FETCH_DATE}**. Layout and rules: `README.md`. '
  f'Re-fetch tooling: `scratch/scec_archive/`.')
w('')
w('## What is archived')
w('')
w('| benchmark | directory | portal id | submission label | code | version | '
  'resolution | date | cplot | on-fault | off-fault | surf-def | size | '
  'published standing (RMS rupture time) |')
w('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for e in idx:
    other = e['nother']
    on = str(e['nfault']) + (f' + {other} named' if other else '')
    w(f"| `{e['bm']}` | `{e['dir']}` | `{e['user']}` | {e['label']} | "
      f"{e['code']} | {e['ver']} | {e['res']} | {e['date']} | "
      f"{'yes' if e['has_cplot'] else 'no'} | {on} | {e['nbody']} | "
      f"{'yes' if e['surfdef'] else '—'} | {e['total']/1e6:.0f} MB | "
      f"{e['standing']} |")
w('')
w(f"**{len(idx)} submissions across 11 benchmark listings: "
  f"{sum(e['nfault']+e['nbody']+e['nother'] for e in idx)} time series, "
  f"{sum(1 for e in idx if e['has_cplot'])} cplot arrays, "
  f"{sum(1 for e in idx if e['surfdef'])} surface-deformation archives, "
  f"{sum(e['total'] for e in idx)/1e6:.0f} MB** (about half that on disk — "
  f"the filesystem compresses it).")
w('')
w('"Published standing" is the median RMS rupture-time difference from every '
  'other submission in that benchmark\'s published `cplot` metric table, with '
  'the rank by that median (1 = closest to the group). A dash means SCEC '
  'publishes no metric column for that submission — see the reasons below.')
w('')
w('## File types')
w('')
w('| type | what | how it was obtained |')
w('|---|---|---|')
w('| `cplot` | rupture-time contour array (x, z, t) | CGI `G1063cplot=Raw '
  'Data` |')
w('| `faultst*` | on-fault time series, 8 columns | CGI '
  '`G1059<station>=Raw Data` |')
w('| `body*` | off-fault time series, 7 columns | CGI '
  '`G1059<station>=Raw Data` |')
w('| named stations (`4064_donna`, ...) | TPV35 uses named receivers instead '
  'of the `body*` grid | CGI `G1059<station>=Raw Data` |')
w('| `TPV3*_..._SCECFinalSurfDisp.zip` | final surface displacement field, '
  'TPV36/37 only | direct download from `tpv3*_surfdef.html` (e-mailed '
  'submission, never served by the CGI) |')
w('')
w("These are the only data the portal exposes. The CGI's own file-list form "
  "was scraped for every `<input>`/`<button>`: besides the `Raw Data` buttons "
  "above it offers only `Select`, `Graph` and `X-Y Plot`, which return "
  "rendered PNGs, not files.")
w('')
w('## Enumeration')
w('')
w('The benchmark list came from the CGI itself (step 1, `G0012=Go -->`): '
  '**56 benchmarks**. Each one\'s user list was then requested '
  '(`G1045<benchmark>= Select `) and searched for `dliu` and its suffixed '
  'variants. Eleven listings carry one:')
w('')
enum = json.load(open(os.path.join(HERE, 'enumeration.json')))
w('| benchmark listing | dliu submissions offered |')
w('|---|---|')
for bm in sorted(enum, key=lambda b: b.replace('_arc', 'z_arc')):
    w(f"| `{bm}` | " + ', '.join(f'`{u}`' for u in enum[bm]['users']) + ' |')
w('')
w('Every other benchmark listing was checked and carries no `dliu*`. The '
  'portal does hold EQdyna results under other ids — `duan`/`duan.2`/`duan.3` '
  '(Benchun Duan, TPV5-13, 16, 17, 24-27, 205, 210), `luo`/`luo.2`/`luo.3` '
  '(Bin Luo, TPV28, 31, 32, 102) and `payne`/`payne.2` (Ryan Payne, TPV22, '
  '23) — but those are other authors\' submissions and are not archived here. '
  '`liu` on TPV8 is Yi Liu, a boundary-integral code, not this owner.')
w('')
w('## Known to exist but not fetchable')
w('')
w('| submission | evidence it exists | why it is not here |')
w('|---|---|---|')
w('| `tpv30` / `dliu.4` — "Dunyu Liu - Finite Element - EQdyna viscous '
  'hg-100m" | column (6) of the published TPV30 metric page, 155.8 ms from '
  '`barall` | the public area does not list it; `mus=dliu.4` returns '
  '*Forbidden Operation* |')
w('| `tpv36` / `dliu.2` — "EQdyna.v5.3.3.50m.dliu" | column (2) of the '
  'published TPV36 metric page, 82.9 ms from `barall` | the public area lists '
  'only `dliu.3` (the large-domain run); `mus=dliu.2` returns *Forbidden '
  'Operation*. The metric page and the public area are snapshots of '
  'different submission sets |')
w('')
w('## Portal quirks worth knowing')
w('')
w('- **The metric page and the public area disagree about who is listed.** '
  'For TPV36 the scored EQdyna column is `dliu.2`, but only `dliu.3` is '
  'served; for TPV37 the scored `dliu` is served but the companion `dliu.2` '
  '(100 m) has no column. Matching is therefore done on the *submission '
  'label*, never on the user id.')
w('- **`*_arc` listings.** `tpv36_arc` and `tpv105-3d_arc` serve files '
  'byte-identical to their non-`_arc` namesakes, so both are archived and '
  'each PROVENANCE says so. `tpv37_arc`/`dliu.3` is **not** a mirror — it is '
  'a distinct large-domain 50 m run, and it is the id the TPV37 '
  'surface-deformation listing refers to.')
w('- **`tpv30`/`dliu` ships a stale `cplot` header** reading '
  '`problem = LVFZ`, `author = Benchun Duan`, `date = 6-2-2012`. Its 36 time '
  'series all say `Dunyu Liu`, 2015-02-27, 100 m, and the rupture-time array '
  'is on the TPV30 fault; the header text is a leftover.')
w('- **`tpv34`/`dliu.2` headers disagree on resolution**: `cplot` says '
  '`element_size=25 m`, code_version 4.0 (matching the portal label '
  '`EQdyna3d_v4.0_3DMPI_Ada_25m`), while the time-series headers still carry '
  'the previous run\'s `50 m x 50 m on fault` and version 3.2.3. The '
  'time step, 0.0020 s over 10001 steps, matches 25 m.')
w('- **TPV36/37 headers carry no `code_version`** (they say '
  '`Project=San-Ti`, `Author=Sophon`). The version in the directory name is '
  'the one in the portal submission label, `EQdyna.v5.3.3`; element size and '
  'date do come from the headers.')
w('')
w('## Verification')
w('')
w('The fetch path was checked against ground truth before anything else was '
  'trusted: `tpv29/eqdyna-v3.1-100m-2015/` already held the author\'s own '
  'retained 2015 upload of `cplot`, `faultst000dp120` and '
  '`body-030st000dp000`. The freshly fetched copies are **byte-identical** to '
  'all three (SHA-256 match; `cplot` 2,982,531 B, `faultst000dp120` 348,217 '
  'B, `body-030st000dp000` 303,279 B). Those three files were left in place, '
  'not overwritten.')
w('')
w('One transform separates the served text from the upload: the CGI prepends '
  'a single space to lines starting with `#` **at column 0**. Submissions '
  'whose writer already indents its comment lines (every EQdyna run from '
  'v4.1 on — the off-fault `  # time_step=` line is served exactly two spaces '
  'wide, which is what the Fortran `a14` format produces unaided) are stored '
  'without stripping anything. Each PROVENANCE.md states which rule applied.')
w('')

open(os.path.join(ARCHIVE, 'INDEX.md'), 'w').write('\n'.join(L))
print(f'INDEX.md: {len(L)} lines')
