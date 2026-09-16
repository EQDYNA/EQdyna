# TPV36_arc — EQdyna v5.3.3, 50m, 2024-11-06

Submission label in the public record: `dliu.3` — "EQdyna.v5.3.3.50m.dliu.largeDomain".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv36_arc   mus=dliu.3
      Q0001=<on-fault station list>    Q0003=<off-fault station list>
      G1063cplot=Raw Data              -> cplot (rupture-time contour array)
      G1059<station>=Raw Data          -> one faultst*/body* time series

Field names were read off the CGI's own file-list form, not guessed; `G1063`
is bound to the `cplot` row only and `G1059` to every station row.
Fetched **2026-09-14**.

### Served text vs uploaded text

The CGI renders each file inside `<pre>` and prepends one space to lines that
start with `#` **at column 0**; every other line is served byte-identical.
This submission's time series were written by EQdyna's
Fortran writer, whose list-directed `WRITE(51,*)` puts one blank in front of
every comment line, and whose `time_step` / `num_time_steps` lines use explicit
`a14` / `a19` formats. Those formats pin the upload's indent exactly — the
off-fault `'# time_step='` literal is 12 characters, so `a14` right-justifies
it with **two** leading blanks, and the CGI serves that line with exactly two.
The CGI therefore prepended nothing to this file, and the served text already
**is** the upload. No space is stripped from the time series here. `cplot`,
which is written separately with column-0 comments, is unrendered normally.

## From the file headers

| field | value | read from |
|---|---|---|
| first header line | `# Project=San-Ti` | `body-010st000dp000` |
| author | Sophon | `body-010st000dp000` |
| code | EQdyna | `body-010st000dp000` |
| code_version | 5.3.3 *(from the public submission label — no `code_version` in any file header)* | `body-010st000dp000` |
| element_size | 50.0000000000000 | `body-010st000dp000` |
| date | 11/ 6/2024 11:22:48 | `body-010st000dp000` |
| time_step | 0.0006  s | `body-010st000dp000` |
| num_time_steps | 34001 | `body-010st000dp000` |
| on-fault node spacing | 50 m | `cplot` coordinates |

## Files (42 — cplot, 19 on-fault, 22 off-fault)

| file | bytes |
|---|---|
| `body-010st000dp000` | 4,015,252 |
| `body-010st100dp000` | 4,016,721 |
| `body-030st000dp000` | 4,015,317 |
| `body-090st000dp000` | 4,015,419 |
| `body-090st100dp000` | 4,017,027 |
| `body-090st200dp000` | 4,017,655 |
| `body010st000dp000` | 4,015,121 |
| `body010st100dp000` | 4,016,539 |
| `body030st000dp000` | 4,014,968 |
| `body090st000dp000` | 4,014,590 |
| `body090st100dp000` | 4,015,375 |
| `body090st200dp000` | 4,016,162 |
| `body150st000dp000` | 4,014,474 |
| `body210st000dp000` | 4,014,606 |
| `body270st000dp000` | 4,014,941 |
| `body270st100dp000` | 4,016,018 |
| `body270st200dp000` | 4,016,827 |
| `body330st000dp000` | 4,015,176 |
| `body390st000dp000` | 4,015,266 |
| `body450st000dp000` | 4,015,323 |
| `body450st100dp000` | 4,016,505 |
| `body450st200dp000` | 4,017,220 |
| `cplot` | 15,172,252 |
| `faultst000dp000` | 4,559,547 |
| `faultst000dp010` | 4,557,175 |
| `faultst000dp030` | 4,557,129 |
| `faultst000dp060` | 4,557,053 |
| `faultst000dp120` | 4,556,847 |
| `faultst000dp180` | 4,556,781 |
| `faultst000dp240` | 4,556,851 |
| `faultst040dp000` | 4,559,995 |
| `faultst040dp030` | 4,557,194 |
| `faultst040dp060` | 4,557,128 |
| `faultst040dp120` | 4,556,902 |
| `faultst040dp180` | 4,556,781 |
| `faultst080dp000` | 4,560,961 |
| `faultst080dp030` | 4,557,324 |
| `faultst080dp060` | 4,557,235 |
| `faultst080dp120` | 4,556,994 |
| `faultst080dp180` | 4,556,826 |
| `faultst120dp000` | 4,561,259 |
| `faultst120dp030` | 4,557,396 |

## Published standing

No public metric page for `tpv36_arc` — <https://strike.scec.org/cvws/metric_cvv1_u1/tpv36_arc/metric_cvv1_tpv36_arc_ac_0.html> returns HTTP 404, so SCEC publishes no cross-code table for this listing.

## Notes

- **Surface deformation.** TPV36 also takes a surface-deformation file, which is e-mailed rather than uploaded to the CGI and is published at <https://strike.scec.org/cvws/tpv36_surfdef.html>. That page lists user `dliu.3`: max uplift **1.031 m**, max subsidence **0.462 m**, against a group range of 0.986–1.069 m uplift over 8 submissions. Archived under this benchmark's `dliu.3` directory.

- The public area lists this same submission twice, once under `tpv36_arc` and once under `tpv36`. Every file fetched here is **byte-identical** to the copy in `tpv36/eqdyna-v5.3.3-50m-2024/`; both listings are kept so the archive mirrors what the portal actually offers.
