# TPV104 — EQdyna3d v4.1, 100m, 2016-10-06

Submission label in the public record: `dliu` — "Dunyu Liu - Finite Element - EQdyna3d_V4.1".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv104   mus=dliu
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
| first header line | `# SCECTPV104` | `body-060st-120dp000` |
| author | D.Liu | `body-060st-120dp000` |
| code | EQdyna3d | `body-060st-120dp000` |
| code_version | 4.1 | `body-060st-120dp000` |
| element_size | 100.000000000000 | `body-060st-120dp000` |
| date | 10/ 6/2016 20:28: 8 | `body-060st-120dp000` |
| time_step | 0.0080  s | `body-060st-120dp000` |
| num_time_steps | 1501 | `body-060st-120dp000` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (16 — cplot, 9 on-fault, 6 off-fault)

| file | bytes |
|---|---|
| `body-060st-120dp000` | 182,218 |
| `body-060st120dp000` | 182,218 |
| `body-090st000dp000` | 182,218 |
| `body060st-120dp000` | 182,218 |
| `body060st120dp000` | 182,218 |
| `body090st000dp000` | 182,218 |
| `cplot` | 2,765,550 |
| `faultst-090dp075` | 236,222 |
| `faultst-120dp030` | 236,222 |
| `faultst-120dp120` | 236,222 |
| `faultst000dp030` | 236,222 |
| `faultst000dp075` | 236,222 |
| `faultst000dp120` | 236,222 |
| `faultst090dp075` | 236,222 |
| `faultst120dp030` | 236,222 |
| `faultst120dp120` | 236,222 |

## Published standing

Not in the published metric page for `tpv104` (<https://strike.scec.org/cvws/metric_cvv1_u1/tpv104/metric_cvv1_tpv104_ac_0.html>) — no column carries this submission label.
