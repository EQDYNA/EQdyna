# TPV37 — EQdyna v5.3.3, 100m, 2024-09-16

Submission label in the public record: `dliu.2` — "EQdyna.v5.3.3.100m.dliu".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv37   mus=dliu.2
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
| element_size | 100.000000000000 | `body-010st000dp000` |
| date | 9/16/2024 15:54:43 | `body-010st000dp000` |
| time_step | 0.0022  s | `body-010st000dp000` |
| num_time_steps | 10200 | `body-010st000dp000` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (9 — cplot, 4 on-fault, 4 off-fault)

| file | bytes |
|---|---|
| `body-010st000dp000` | 1,204,443 |
| `body-010st100dp000` | 1,204,752 |
| `body-030st000dp000` | 1,204,446 |
| `body-090st000dp000` | 1,204,473 |
| `cplot` | 3,806,152 |
| `faultst000dp000` | 1,367,751 |
| `faultst000dp010` | 1,367,485 |
| `faultst000dp030` | 1,367,471 |
| `faultst000dp060` | 1,367,454 |

## Published standing

Not in the published metric page for `tpv37` (<https://strike.scec.org/cvws/metric_cvv1_u1/tpv37/metric_cvv1_tpv37_ac_0.html>) — no column carries this submission label. The page does carry other EQdyna column(s), under different labels: `dliu` = "EQdyna.v5.3.3.50m.dliu".

## Notes

- The TPV37 surface-deformation listing (<https://strike.scec.org/cvws/tpv37_surfdef.html>) carries no entry for user `dliu.2`; its EQdyna entry is `dliu.3`, archived under `tpv37_arc/`.
