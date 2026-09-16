# TPV34 — EQdyna3d v4.0, 25m, 2016-03-06

Submission label in the public record: `dliu.2` — "Dunyu Liu_EQdyna3d_v4.0_3DMPI_Ada_25m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv34   mus=dliu.2
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
| code_version | 4.0 | `cplot` |
| element_size | 25 m | `cplot` |
| date | 2016/03/07 | `cplot` |
| on-fault node spacing | 25 m | `cplot` coordinates |

## Files (92 — cplot, 35 on-fault, 56 off-fault)

| file | bytes |
|---|---|
| `body-030st-100dp000` | 1,210,701 |
| `body-030st-100dp024` | 1,210,701 |
| `body-030st-200dp000` | 1,210,701 |
| `body-030st-200dp024` | 1,210,701 |
| `body-030st000dp000` | 1,210,701 |
| `body-030st000dp024` | 1,210,701 |
| `body-030st100dp000` | 1,210,701 |
| `body-030st100dp024` | 1,210,701 |
| `body-030st200dp000` | 1,210,701 |
| `body-030st200dp024` | 1,210,701 |
| `body-090st-100dp000` | 1,210,701 |
| `body-090st-100dp024` | 1,210,701 |
| `body-090st-200dp000` | 1,210,701 |
| `body-090st-200dp024` | 1,210,701 |
| `body-090st000dp000` | 1,210,701 |
| `body-090st000dp024` | 1,210,701 |
| `body-090st100dp000` | 1,210,701 |
| `body-090st100dp024` | 1,210,701 |
| `body-090st200dp000` | 1,210,701 |
| `body-090st200dp024` | 1,210,701 |
| `body-150st-150dp000` | 1,210,701 |
| `body-150st-150dp024` | 1,210,701 |
| `body-150st000dp000` | 1,210,701 |
| `body-150st000dp024` | 1,210,701 |
| `body-150st150dp000` | 1,210,701 |
| `body-150st150dp024` | 1,210,701 |
| `body000st-200dp000` | 1,210,701 |
| `body000st-200dp024` | 1,210,701 |
| `body000st200dp000` | 1,210,701 |
| `body000st200dp024` | 1,210,701 |
| `body030st-100dp000` | 1,210,701 |
| `body030st-100dp024` | 1,210,701 |
| `body030st-200dp000` | 1,210,701 |
| `body030st-200dp024` | 1,210,701 |
| `body030st000dp000` | 1,210,701 |
| `body030st000dp024` | 1,210,701 |
| `body030st100dp000` | 1,210,701 |
| `body030st100dp024` | 1,210,701 |
| `body030st200dp000` | 1,210,701 |
| `body030st200dp024` | 1,210,701 |
| `body090st-100dp000` | 1,210,701 |
| `body090st-100dp024` | 1,210,701 |
| `body090st-200dp000` | 1,210,701 |
| `body090st-200dp024` | 1,210,701 |
| `body090st000dp000` | 1,210,701 |
| `body090st000dp024` | 1,210,701 |
| `body090st100dp000` | 1,210,701 |
| `body090st100dp024` | 1,210,701 |
| `body090st200dp000` | 1,210,701 |
| `body090st200dp024` | 1,210,701 |
| `body150st-150dp000` | 1,210,701 |
| `body150st-150dp024` | 1,210,701 |
| `body150st000dp000` | 1,210,701 |
| `body150st000dp024` | 1,210,701 |
| `body150st150dp000` | 1,210,701 |
| `body150st150dp024` | 1,210,701 |
| `cplot` | 30,343,079 |
| `faultst-060dp000` | 1,390,660 |
| `faultst-060dp010` | 1,390,660 |
| `faultst-060dp024` | 1,390,660 |
| `faultst-060dp050` | 1,390,660 |
| `faultst-060dp075` | 1,390,660 |
| `faultst-060dp100` | 1,390,660 |
| `faultst-060dp120` | 1,390,660 |
| `faultst-120dp000` | 1,390,660 |
| `faultst-120dp010` | 1,390,660 |
| `faultst-120dp024` | 1,390,660 |
| `faultst-120dp050` | 1,390,660 |
| `faultst-120dp075` | 1,390,660 |
| `faultst-120dp100` | 1,390,660 |
| `faultst-120dp120` | 1,390,660 |
| `faultst000dp000` | 1,390,660 |
| `faultst000dp010` | 1,390,660 |
| `faultst000dp024` | 1,390,660 |
| `faultst000dp050` | 1,390,660 |
| `faultst000dp075` | 1,390,660 |
| `faultst000dp100` | 1,390,660 |
| `faultst000dp120` | 1,390,660 |
| `faultst060dp000` | 1,390,660 |
| `faultst060dp010` | 1,390,660 |
| `faultst060dp024` | 1,390,660 |
| `faultst060dp050` | 1,390,660 |
| `faultst060dp075` | 1,390,660 |
| `faultst060dp100` | 1,390,660 |
| `faultst060dp120` | 1,390,660 |
| `faultst120dp000` | 1,390,660 |
| `faultst120dp010` | 1,390,660 |
| `faultst120dp024` | 1,390,660 |
| `faultst120dp050` | 1,390,660 |
| `faultst120dp075` | 1,390,660 |
| `faultst120dp100` | 1,390,660 |
| `faultst120dp120` | 1,390,660 |

## Published standing

Not in the published metric page for `tpv34` (<https://strike.scec.org/cvws/metric_cvv1_u1/tpv34/metric_cvv1_tpv34_ac_0.html>) — no column carries this submission label. The page does carry other EQdyna column(s), under different labels: `dliu` = "Dunyu Liu - Finite Element - EQdyna -50m".

## Notes

- `cplot` and the time series disagree on **version**: `cplot` says `4.0`, the time series say `3.2.3`.

- `cplot` and the time series disagree on **element_size**: `cplot` says `25.0`, the time series say `50.0`.

- `cplot` and the time series disagree on **date**: `cplot` says `2016/03/07`, the time series say `3/ 6/2016 20:10: 4`.
