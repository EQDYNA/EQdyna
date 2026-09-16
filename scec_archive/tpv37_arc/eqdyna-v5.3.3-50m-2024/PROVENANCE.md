# TPV37_arc — EQdyna v5.3.3, 50m, 2024-11-06

Submission label in the public record: `dliu.3` — "EQdyna.v5.3.3.50m.dliu.largeDomain".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv37_arc   mus=dliu.3
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
| date | 11/ 6/2024 22:39: 7 | `body-010st000dp000` |
| time_step | 0.0005  s | `body-010st000dp000` |
| num_time_steps | 40801 | `body-010st000dp000` |
| on-fault node spacing | 50 m | `cplot` coordinates |

## Files (46 — cplot, 22 on-fault, 22 off-fault, surface-deformation zip)

| file | bytes |
|---|---|
| `body-010st000dp000` | 4,818,234 |
| `body-010st100dp000` | 4,820,017 |
| `body-030st000dp000` | 4,818,310 |
| `body-090st000dp000` | 4,818,432 |
| `body-090st100dp000` | 4,820,354 |
| `body-090st200dp000` | 4,821,101 |
| `body010st000dp000` | 4,818,082 |
| `body010st100dp000` | 4,819,792 |
| `body030st000dp000` | 4,817,910 |
| `body090st000dp000` | 4,817,507 |
| `body090st100dp000` | 4,818,397 |
| `body090st200dp000` | 4,819,314 |
| `body150st000dp000` | 4,817,408 |
| `body210st000dp000` | 4,817,525 |
| `body270st000dp000` | 4,817,898 |
| `body270st100dp000` | 4,819,165 |
| `body270st200dp000` | 4,820,108 |
| `body330st000dp000` | 4,818,166 |
| `body390st000dp000` | 4,818,269 |
| `body450st000dp000` | 4,818,338 |
| `body450st100dp000` | 4,819,717 |
| `body450st200dp000` | 4,820,565 |
| `cplot` | 15,172,252 |
| `faultst000dp000` | 5,471,305 |
| `faultst000dp010` | 5,468,468 |
| `faultst000dp030` | 5,468,414 |
| `faultst000dp060` | 5,468,327 |
| `faultst000dp120` | 5,468,086 |
| `faultst000dp180` | 5,467,981 |
| `faultst000dp240` | 5,468,085 |
| `faultst040dp000` | 5,471,905 |
| `faultst040dp030` | 5,468,494 |
| `faultst040dp060` | 5,468,417 |
| `faultst040dp120` | 5,468,147 |
| `faultst040dp180` | 5,467,981 |
| `faultst080dp000` | 5,473,176 |
| `faultst080dp030` | 5,468,645 |
| `faultst080dp060` | 5,468,538 |
| `faultst080dp120` | 5,468,267 |
| `faultst080dp180` | 5,468,065 |
| `faultst120dp000` | 5,473,285 |
| `faultst120dp030` | 5,468,725 |
| `faultst120dp060` | 5,468,634 |
| `faultst120dp120` | 5,468,367 |
| `faultst120dp180` | 5,468,172 |
| `TPV37_EQdyna_v5_3_3_50m_dliu_SCECFinalSurfDisp.zip` | 11,906,185 |

## Published standing

No public metric page for `tpv37_arc` — <https://strike.scec.org/cvws/metric_cvv1_u1/tpv37_arc/metric_cvv1_tpv37_arc_ac_0.html> returns HTTP 404, so SCEC publishes no cross-code table for this listing.

## Notes

- **Surface deformation.** TPV37 also takes a surface-deformation file, which is e-mailed rather than uploaded to the CGI and is published at <https://strike.scec.org/cvws/tpv37_surfdef.html>. That page lists user `dliu.3`: max uplift **0.928 m**, max subsidence **0.425 m**, against a group range of 0.888–0.977 m uplift over 8 submissions. Archived here as `TPV37_EQdyna_v5_3_3_50m_dliu_SCECFinalSurfDisp.zip`.
