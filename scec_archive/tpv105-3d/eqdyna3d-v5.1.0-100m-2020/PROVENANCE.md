# TPV105-3D — EQdyna3D v5.1.0, 100m, 2020-09-11

Submission label in the public record: `dliu.2` — "D.Liu, B. Luo - FEM - EQdyna_3D_5.1.0-100m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv105-3d   mus=dliu.2
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
| first header line | `# SCECTPV105-3D` | `body-060st-120dp000` |
| author | D.Liu & B.Luo | `body-060st-120dp000` |
| code | EQdyna3D | `body-060st-120dp000` |
| code_version | 5.1.0 | `body-060st-120dp000` |
| element_size | 100.000000000000 | `body-060st-120dp000` |
| date | 9/11/2020 15:56:33 | `body-060st-120dp000` |
| time_step | 0.0080  s | `body-060st-120dp000` |
| num_time_steps | 1876 | `body-060st-120dp000` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (20 — cplot, 13 on-fault, 6 off-fault)

| file | bytes |
|---|---|
| `body-060st-120dp000` | 210,721 |
| `body-060st120dp000` | 210,721 |
| `body-090st000dp000` | 210,770 |
| `body060st-120dp000` | 210,721 |
| `body060st120dp000` | 210,721 |
| `body090st000dp000` | 210,770 |
| `cplot` | 4,172,149 |
| `faultst-090dp075` | 323,338 |
| `faultst-120dp030` | 323,338 |
| `faultst-120dp120` | 323,339 |
| `faultst-150dp075` | 323,338 |
| `faultst-180dp075` | 323,338 |
| `faultst000dp030` | 323,338 |
| `faultst000dp075` | 323,338 |
| `faultst000dp120` | 323,338 |
| `faultst090dp075` | 323,371 |
| `faultst120dp030` | 323,338 |
| `faultst120dp120` | 323,339 |
| `faultst150dp075` | 323,338 |
| `faultst180dp075` | 323,338 |

## Published standing

Published as column `dliu.2` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv105-3d/metric_cvv1_tpv105-3d_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 4 submissions: **44.6 ms**
- closest submission `barall.2` at 34.7 ms; farthest `vyas` at 75.2 ms
- rank by that median among the 5 submissions: **3 of 5** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `barall.2` | 34.7 |
| `wang.2` | 35.5 |
| `chen.2` | 53.8 |
| `vyas` | 75.2 |

## Notes

- The public area lists this same submission twice, once under `tpv105-3d` and once under `tpv105-3d_arc`. Every file fetched here is **byte-identical** to the copy in `tpv105-3d_arc/eqdyna3d-v5.1.0-100m-2020/`; both listings are kept so the archive mirrors what the portal actually offers.
