# TPV34 — EQdyna3d v3.2.3, 50m, 2016-02-08

Submission label in the public record: `dliu` — "Dunyu Liu - Finite Element - EQdyna -50m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv34   mus=dliu
      Q0001=<on-fault station list>    Q0003=<off-fault station list>
      G1063cplot=Raw Data              -> cplot (rupture-time contour array)
      G1059<station>=Raw Data          -> one faultst*/body* time series

Field names were read off the CGI's own file-list form, not guessed; `G1063`
is bound to the `cplot` row only and `G1059` to every station row.
Fetched **2026-09-14**.

### Served text vs uploaded text

The CGI renders each file inside `<pre>` and prepends one space to lines that
start with `#` **at column 0**; every other line is served byte-identical.
This submission's uploaded comment lines start at column 0
(the same writer as the TPV29 `dliu` upload the fetch path was verified
against), so the one prepended space is stripped back off and the files below
are the uploads verbatim, byte for byte.

## From the file headers

| field | value | read from |
|---|---|---|
| first header line | `# This is the file header:` | `cplot` |
| author | D.Liu | `cplot` |
| code | EQdyna3d | `cplot` |
| code_version | 3.2.3 | `cplot` |
| element_size | 50 m | `cplot` |
| date | 2016/01/30 | `cplot` |
| on-fault node spacing | 50 m | `cplot` coordinates |

## Files (92 — cplot, 35 on-fault, 56 off-fault)

| file | bytes |
|---|---|
| `body-030st-100dp000` | 605,701 |
| `body-030st-100dp024` | 605,701 |
| `body-030st-200dp000` | 605,701 |
| `body-030st-200dp024` | 605,701 |
| `body-030st000dp000` | 605,701 |
| `body-030st000dp024` | 605,701 |
| `body-030st100dp000` | 605,701 |
| `body-030st100dp024` | 605,701 |
| `body-030st200dp000` | 605,701 |
| `body-030st200dp024` | 605,701 |
| `body-090st-100dp000` | 605,701 |
| `body-090st-100dp024` | 605,701 |
| `body-090st-200dp000` | 605,701 |
| `body-090st-200dp024` | 605,701 |
| `body-090st000dp000` | 605,701 |
| `body-090st000dp024` | 605,701 |
| `body-090st100dp000` | 605,701 |
| `body-090st100dp024` | 605,701 |
| `body-090st200dp000` | 605,701 |
| `body-090st200dp024` | 605,701 |
| `body-150st-150dp000` | 605,701 |
| `body-150st-150dp024` | 605,701 |
| `body-150st000dp000` | 605,701 |
| `body-150st000dp024` | 605,701 |
| `body-150st150dp000` | 605,701 |
| `body-150st150dp024` | 605,701 |
| `body000st-200dp000` | 605,701 |
| `body000st-200dp024` | 605,701 |
| `body000st200dp000` | 605,701 |
| `body000st200dp024` | 605,701 |
| `body030st-100dp000` | 605,701 |
| `body030st-100dp024` | 605,701 |
| `body030st-200dp000` | 605,701 |
| `body030st-200dp024` | 605,701 |
| `body030st000dp000` | 605,701 |
| `body030st000dp024` | 605,701 |
| `body030st100dp000` | 605,701 |
| `body030st100dp024` | 605,701 |
| `body030st200dp000` | 605,701 |
| `body030st200dp024` | 605,701 |
| `body090st-100dp000` | 605,701 |
| `body090st-100dp024` | 605,701 |
| `body090st-200dp000` | 605,701 |
| `body090st-200dp024` | 605,701 |
| `body090st000dp000` | 605,701 |
| `body090st000dp024` | 605,701 |
| `body090st100dp000` | 605,701 |
| `body090st100dp024` | 605,701 |
| `body090st200dp000` | 605,701 |
| `body090st200dp024` | 605,701 |
| `body150st-150dp000` | 605,701 |
| `body150st-150dp024` | 605,701 |
| `body150st000dp000` | 605,701 |
| `body150st000dp024` | 605,701 |
| `body150st150dp000` | 605,701 |
| `body150st150dp024` | 605,701 |
| `cplot` | 8,502,906 |
| `faultst-060dp000` | 695,660 |
| `faultst-060dp010` | 695,660 |
| `faultst-060dp024` | 695,660 |
| `faultst-060dp050` | 695,660 |
| `faultst-060dp075` | 695,660 |
| `faultst-060dp100` | 695,660 |
| `faultst-060dp120` | 695,660 |
| `faultst-120dp000` | 695,660 |
| `faultst-120dp010` | 695,660 |
| `faultst-120dp024` | 695,660 |
| `faultst-120dp050` | 695,660 |
| `faultst-120dp075` | 695,660 |
| `faultst-120dp100` | 695,660 |
| `faultst-120dp120` | 695,660 |
| `faultst000dp000` | 695,660 |
| `faultst000dp010` | 695,660 |
| `faultst000dp024` | 695,660 |
| `faultst000dp050` | 695,660 |
| `faultst000dp075` | 695,660 |
| `faultst000dp100` | 695,660 |
| `faultst000dp120` | 695,660 |
| `faultst060dp000` | 695,660 |
| `faultst060dp010` | 695,660 |
| `faultst060dp024` | 695,660 |
| `faultst060dp050` | 695,660 |
| `faultst060dp075` | 695,660 |
| `faultst060dp100` | 695,660 |
| `faultst060dp120` | 695,660 |
| `faultst120dp000` | 695,660 |
| `faultst120dp010` | 695,660 |
| `faultst120dp024` | 695,660 |
| `faultst120dp050` | 695,660 |
| `faultst120dp075` | 695,660 |
| `faultst120dp100` | 695,660 |
| `faultst120dp120` | 695,660 |

## Published standing

Published as column `dliu` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv34/metric_cvv1_tpv34_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 13 submissions: **32.4 ms**
- closest submission `bydlon.2` at 12.2 ms; farthest `roten` at 44.6 ms
- rank by that median among the 14 submissions: **13 of 14** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `bydlon.2` | 12.2 |
| `chen.2` | 15.4 |
| `bai` | 19.8 |
| `kaneko.2` | 23.2 |
| `barall` | 23.6 |
| `barall.2` | 23.6 |
| `ma` | 32.4 |
| `ma.2` | 32.5 |
| `chen` | 33.3 |
| `bydlon` | 35.9 |
| `daub` | 38.2 |
| `kaneko` | 42.3 |
| `roten` | 44.6 |

## Notes

- `cplot` and the time series disagree on **date**: `cplot` says `2016/01/30`, the time series say `2/ 8/2016 19:12:23`.
