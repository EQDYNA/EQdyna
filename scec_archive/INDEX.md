# INDEX — EQdyna submissions to the SCEC/USGS code-verification project

Everything archived here, fetched from the public SCEC cvws portal on **2026-09-14**. Layout and rules: `README.md`. Re-fetch tooling: `scratch/scec_archive/`.

## What is archived

| benchmark | directory | portal id | submission label | code | version | resolution | date | cplot | on-fault | off-fault | surf-def | size | published standing (RMS rupture time) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `tpv104` | `eqdyna3d-v4.1-100m-2016` | `dliu` | Dunyu Liu - Finite Element - EQdyna3d_V4.1 | EQdyna3d | 4.1 | 100m | 2016-10-06 | yes | 9 | 6 | — | 6 MB | — |
| `tpv105-3d` | `eqdyna3d-v5.1.0-100m-2020` | `dliu.2` | D.Liu, B. Luo - FEM - EQdyna_3D_5.1.0-100m | EQdyna3D | 5.1.0 | 100m | 2020-09-11 | yes | 13 | 6 | — | 10 MB | 44.6 ms median RMS, rank 3/5 |
| `tpv105-3d` | `eqdyna3d-v5.1.0-200m-2020` | `dliu` | D.Liu, B. Luo - FEM - EQdyna_3D_5.1.0-200m | EQdyna3D | 5.1.0 | 200m | 2020-09-11 | yes | 13 | 6 | — | 4 MB | — |
| `tpv105-3d_arc` | `eqdyna3d-v5.1.0-100m-2020` | `dliu.2` | D.Liu, B. Luo - FEM - EQdyna_3D_5.1.0-100m | EQdyna3D | 5.1.0 | 100m | 2020-09-11 | yes | 13 | 6 | — | 10 MB | — |
| `tpv29` | `eqdyna-v3.1-100m-2015` | `dliu` | Dunyu Liu - Finite Element - EQdyna-100m | EQdyna | 3.1 | 100m | 2015-02-27 | yes | 24 | 12 | — | 15 MB | 87.2 ms median RMS, rank 15/15 |
| `tpv29` | `eqdyna-v3.1-50m-2015` | `dliu.2` | Dunyu Liu - Finite Element - EQdyna-50m | EQdyna | 3.1 | 50m | 2015-02-27 | yes | 24 | 12 | — | 36 MB | 33.0 ms median RMS, rank 12/15 |
| `tpv30` | `eqdyna-v3.1-100m-2015` | `dliu` | Dunyu Liu - Finite Element - EQdyna-100m | EQdyna | 3.1 | 100m | 2015-02-27 | yes | 24 | 12 | — | 15 MB | 287.7 ms median RMS, rank 14/14 |
| `tpv30` | `eqdyna-v3.1-25m-2015` | `dliu.3` | Dunyu Liu - Finite Element - EQdyna - 25m | EQdyna | 3.1 | 25m | 2015-03-06 | yes | 24 | 12 | — | 83 MB | 55.1 ms median RMS, rank 8/14 |
| `tpv30` | `eqdyna-v3.1-50m-2015` | `dliu.2` | Dunyu Liu - Finite Element - EQdyna-50m | EQdyna | 3.1 | 50m | 2015-02-27 | yes | 24 | 12 | — | 36 MB | 120.6 ms median RMS, rank 12/14 |
| `tpv34` | `eqdyna3d-v4.0-25m-2016` | `dliu.2` | Dunyu Liu_EQdyna3d_v4.0_3DMPI_Ada_25m | EQdyna3d | 4.0 | 25m | 2016-03-06 | yes | 35 | 56 | — | 147 MB | — |
| `tpv34` | `eqdyna3d-v3.2.3-50m-2016` | `dliu` | Dunyu Liu - Finite Element - EQdyna -50m | EQdyna3d | 3.2.3 | 50m | 2016-02-08 | yes | 35 | 56 | — | 67 MB | 32.4 ms median RMS, rank 13/14 |
| `tpv35` | `eqdyna3d-v4.1.3-100m-2017` | `dliu` | Dunyu Liu - Finite Element - EQdyna3Dv4.1.3 - 100m | EQdyna3D | 4.1.3 | 100m | 2017-01-17 | yes | 7 + 42 named | 0 | — | 16 MB | — |
| `tpv36` | `eqdyna-v5.3.3-50m-2024` | `dliu.3` | EQdyna.v5.3.3.50m.dliu.largeDomain | EQdyna | 5.3.3 | 50m | 2024-11-06 | yes | 19 | 22 | yes | 202 MB | — |
| `tpv36_arc` | `eqdyna-v5.3.3-50m-2024` | `dliu.3` | EQdyna.v5.3.3.50m.dliu.largeDomain | EQdyna | 5.3.3 | 50m | 2024-11-06 | yes | 19 | 22 | — | 190 MB | — |
| `tpv37` | `eqdyna-v5.3.3-100m-2024` | `dliu.2` | EQdyna.v5.3.3.100m.dliu | EQdyna | 5.3.3 | 100m | 2024-09-16 | yes | 4 | 4 | — | 14 MB | — |
| `tpv37` | `eqdyna-v5.3.3-50m-2024` | `dliu` | EQdyna.v5.3.3.50m.dliu | EQdyna | 5.3.3 | 50m | 2024-09-16 | yes | 19 | 19 | — | 113 MB | 60.9 ms median RMS, rank 7/8 |
| `tpv37_arc` | `eqdyna-v5.3.3-50m-2024` | `dliu.3` | EQdyna.v5.3.3.50m.dliu.largeDomain | EQdyna | 5.3.3 | 50m | 2024-11-06 | yes | 22 | 22 | yes | 253 MB | — |

**17 submissions across 11 benchmark listings: 655 time series, 17 cplot arrays, 2 surface-deformation archives, 1216 MB** (about half that on disk — the filesystem compresses it).

"Published standing" is the median RMS rupture-time difference from every other submission in that benchmark's published `cplot` metric table, with the rank by that median (1 = closest to the group). A dash means SCEC publishes no metric column for that submission — see the reasons below.

## File types

| type | what | how it was obtained |
|---|---|---|
| `cplot` | rupture-time contour array (x, z, t) | CGI `G1063cplot=Raw Data` |
| `faultst*` | on-fault time series, 8 columns | CGI `G1059<station>=Raw Data` |
| `body*` | off-fault time series, 7 columns | CGI `G1059<station>=Raw Data` |
| named stations (`4064_donna`, ...) | TPV35 uses named receivers instead of the `body*` grid | CGI `G1059<station>=Raw Data` |
| `TPV3*_..._SCECFinalSurfDisp.zip` | final surface displacement field, TPV36/37 only | direct download from `tpv3*_surfdef.html` (e-mailed submission, never served by the CGI) |

These are the only data the portal exposes. The CGI's own file-list form was scraped for every `<input>`/`<button>`: besides the `Raw Data` buttons above it offers only `Select`, `Graph` and `X-Y Plot`, which return rendered PNGs, not files.

## Enumeration

The benchmark list came from the CGI itself (step 1, `G0012=Go -->`): **56 benchmarks**. Each one's user list was then requested (`G1045<benchmark>= Select `) and searched for `dliu` and its suffixed variants. Eleven listings carry one:

| benchmark listing | dliu submissions offered |
|---|---|
| `tpv104` | `dliu` |
| `tpv105-3d` | `dliu`, `dliu.2` |
| `tpv105-3d_arc` | `dliu.2` |
| `tpv29` | `dliu`, `dliu.2` |
| `tpv30` | `dliu`, `dliu.2`, `dliu.3` |
| `tpv34` | `dliu`, `dliu.2` |
| `tpv35` | `dliu` |
| `tpv36` | `dliu.3` |
| `tpv36_arc` | `dliu.3` |
| `tpv37` | `dliu`, `dliu.2` |
| `tpv37_arc` | `dliu.3` |

Every other benchmark listing was checked and carries no `dliu*`. The portal does hold EQdyna results under other ids — `duan`/`duan.2`/`duan.3` (Benchun Duan, TPV5-13, 16, 17, 24-27, 205, 210), `luo`/`luo.2`/`luo.3` (Bin Luo, TPV28, 31, 32, 102) and `payne`/`payne.2` (Ryan Payne, TPV22, 23) — but those are other authors' submissions and are not archived here. `liu` on TPV8 is Yi Liu, a boundary-integral code, not this owner.

## Known to exist but not fetchable

| submission | evidence it exists | why it is not here |
|---|---|---|
| `tpv30` / `dliu.4` — "Dunyu Liu - Finite Element - EQdyna viscous hg-100m" | column (6) of the published TPV30 metric page, 155.8 ms from `barall` | the public area does not list it; `mus=dliu.4` returns *Forbidden Operation* |
| `tpv36` / `dliu.2` — "EQdyna.v5.3.3.50m.dliu" | column (2) of the published TPV36 metric page, 82.9 ms from `barall` | the public area lists only `dliu.3` (the large-domain run); `mus=dliu.2` returns *Forbidden Operation*. The metric page and the public area are snapshots of different submission sets |

## Portal quirks worth knowing

- **The metric page and the public area disagree about who is listed.** For TPV36 the scored EQdyna column is `dliu.2`, but only `dliu.3` is served; for TPV37 the scored `dliu` is served but the companion `dliu.2` (100 m) has no column. Matching is therefore done on the *submission label*, never on the user id.
- **`*_arc` listings.** `tpv36_arc` and `tpv105-3d_arc` serve files byte-identical to their non-`_arc` namesakes, so both are archived and each PROVENANCE says so. `tpv37_arc`/`dliu.3` is **not** a mirror — it is a distinct large-domain 50 m run, and it is the id the TPV37 surface-deformation listing refers to.
- **`tpv30`/`dliu` ships a stale `cplot` header** reading `problem = LVFZ`, `author = Benchun Duan`, `date = 6-2-2012`. Its 36 time series all say `Dunyu Liu`, 2015-02-27, 100 m, and the rupture-time array is on the TPV30 fault; the header text is a leftover.
- **`tpv34`/`dliu.2` headers disagree on resolution**: `cplot` says `element_size=25 m`, code_version 4.0 (matching the portal label `EQdyna3d_v4.0_3DMPI_Ada_25m`), while the time-series headers still carry the previous run's `50 m x 50 m on fault` and version 3.2.3. The time step, 0.0020 s over 10001 steps, matches 25 m.
- **TPV36/37 headers carry no `code_version`** (they say `Project=San-Ti`, `Author=Sophon`). The version in the directory name is the one in the portal submission label, `EQdyna.v5.3.3`; element size and date do come from the headers.

## Verification

The fetch path was checked against ground truth before anything else was trusted: `tpv29/eqdyna-v3.1-100m-2015/` already held the author's own retained 2015 upload of `cplot`, `faultst000dp120` and `body-030st000dp000`. The freshly fetched copies are **byte-identical** to all three (SHA-256 match; `cplot` 2,982,531 B, `faultst000dp120` 348,217 B, `body-030st000dp000` 303,279 B). Those three files were left in place, not overwritten.

One transform separates the served text from the upload: the CGI prepends a single space to lines starting with `#` **at column 0**. Submissions whose writer already indents its comment lines (every EQdyna run from v4.1 on — the off-fault `  # time_step=` line is served exactly two spaces wide, which is what the Fortran `a14` format produces unaided) are stored without stripping anything. Each PROVENANCE.md states which rule applied.
