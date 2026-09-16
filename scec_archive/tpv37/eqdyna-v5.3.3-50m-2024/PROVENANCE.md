# TPV37 — EQdyna v5.3.3, 50m, 2024-09-16

Submission label in the public record: `dliu` — "EQdyna.v5.3.3.50m.dliu".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv37   mus=dliu
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
| date | 9/16/2024 18:54:16 | `body-010st000dp000` |
| time_step | 0.0011  s | `body-010st000dp000` |
| num_time_steps | 20400 | `body-010st000dp000` |
| on-fault node spacing | 50 m | `cplot` coordinates |

## Files (39 — cplot, 19 on-fault, 19 off-fault)

| file | bytes |
|---|---|
| `body-010st000dp000` | 2,408,967 |
| `body-010st100dp000` | 2,409,731 |
| `body-030st000dp000` | 2,408,987 |
| `body-090st000dp000` | 2,409,033 |
| `body-090st100dp000` | 2,409,810 |
| `body-090st200dp000` | 2,409,947 |
| `body010st000dp000` | 2,408,894 |
| `body010st100dp000` | 2,409,638 |
| `body030st000dp000` | 2,408,766 |
| `body090st000dp000` | 2,408,335 |
| `body090st100dp000` | 2,409,019 |
| `body090st200dp000` | 2,409,597 |
| `body150st000dp000` | 2,408,205 |
| `body210st000dp000` | 2,408,291 |
| `body270st000dp000` | 2,408,601 |
| `body270st100dp000` | 2,409,488 |
| `body270st200dp000` | 2,410,036 |
| `body330st000dp000` | 2,408,880 |
| `body390st000dp000` | 2,409,004 |
| `cplot` | 15,172,252 |
| `faultst000dp000` | 2,735,570 |
| `faultst000dp010` | 2,734,430 |
| `faultst000dp030` | 2,734,405 |
| `faultst000dp060` | 2,734,347 |
| `faultst000dp120` | 2,734,255 |
| `faultst000dp180` | 2,734,247 |
| `faultst000dp240` | 2,734,256 |
| `faultst040dp000` | 2,735,880 |
| `faultst040dp030` | 2,734,443 |
| `faultst040dp060` | 2,734,383 |
| `faultst040dp120` | 2,734,277 |
| `faultst040dp180` | 2,734,247 |
| `faultst080dp000` | 2,736,312 |
| `faultst080dp030` | 2,734,517 |
| `faultst080dp060` | 2,734,465 |
| `faultst080dp120` | 2,734,314 |
| `faultst080dp180` | 2,734,247 |
| `faultst120dp000` | 2,736,503 |
| `faultst120dp030` | 2,734,566 |

## Published standing

Published as column `dliu` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv37/metric_cvv1_tpv37_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 7 submissions: **60.9 ms**
- closest submission `li` at 39.1 ms; farthest `barall` at 80.4 ms
- rank by that median among the 8 submissions: **7 of 8** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `li` | 39.1 |
| `wzhang` | 48.9 |
| `kutschera.2` | 54.5 |
| `yang` | 60.9 |
| `wang` | 67.7 |
| `ma` | 68.7 |
| `barall` | 80.4 |

## Notes

- The TPV37 surface-deformation listing (<https://strike.scec.org/cvws/tpv37_surfdef.html>) carries no entry for user `dliu`; its EQdyna entry is `dliu.3`, archived under `tpv37_arc/`.
