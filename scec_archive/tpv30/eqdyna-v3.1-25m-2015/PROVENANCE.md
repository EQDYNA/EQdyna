# TPV30 — EQdyna v3.1, 25m, 2015-03-06

Submission label in the public record: `dliu.3` — "Dunyu Liu - Finite Element - EQdyna - 25m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv30   mus=dliu.3
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
| first header line | `# problem = TPV30` | `cplot` |
| author | Dunyu Liu | `cplot` |
| code | EQdyna | `cplot` |
| code_version | 3.1 | `cplot` |
| element_size | 25 m x 25 m on fault | `cplot` |
| date | 3-6-2015 | `cplot` |
| on-fault node spacing | 25 m | `cplot` coordinates |

## Files (37 — cplot, 24 on-fault, 12 off-fault)

| file | bytes |
|---|---|
| `body-030st-150dp000` | 908,279 |
| `body-030st000dp000` | 908,277 |
| `body-030st150dp000` | 908,278 |
| `body-200st-200dp000` | 908,280 |
| `body-200st000dp000` | 908,278 |
| `body-200st200dp000` | 908,279 |
| `body030st-150dp000` | 908,278 |
| `body030st000dp000` | 908,276 |
| `body030st150dp000` | 908,277 |
| `body200st-200dp000` | 908,279 |
| `body200st000dp000` | 908,277 |
| `body200st200dp000` | 908,280 |
| `cplot` | 47,449,138 |
| `faultst-042dp061` | 1,043,217 |
| `faultst-050dp000` | 1,043,215 |
| `faultst-050dp100` | 1,043,216 |
| `faultst-050dp160` | 1,043,214 |
| `faultst-089dp101` | 1,043,218 |
| `faultst-110dp014` | 1,043,218 |
| `faultst-150dp050` | 1,043,214 |
| `faultst-150dp120` | 1,043,217 |
| `faultst-180dp156` | 1,043,219 |
| `faultst000dp120` | 1,043,215 |
| `faultst043dp062` | 1,043,216 |
| `faultst046dp060` | 1,043,214 |
| `faultst050dp000` | 1,043,214 |
| `faultst050dp120` | 1,043,213 |
| `faultst051dp057` | 1,043,216 |
| `faultst059dp047` | 1,043,216 |
| `faultst090dp153` | 1,043,215 |
| `faultst093dp153` | 1,043,217 |
| `faultst100dp050` | 1,043,215 |
| `faultst100dp110` | 1,043,214 |
| `faultst150dp000` | 1,043,213 |
| `faultst150dp130` | 1,043,214 |
| `faultst167dp105` | 1,043,218 |
| `faultst170dp045` | 1,043,215 |

## Published standing

Published as column `dliu.3` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv30/metric_cvv1_tpv30_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 13 submissions: **55.1 ms**
- closest submission `kozdon.2` at 38.3 ms; farthest `dliu` at 245.3 ms
- rank by that median among the 14 submissions: **8 of 14** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `kozdon.2` | 38.3 |
| `duru.2` | 40.2 |
| `shi.2` | 43.4 |
| `kozdon` | 45.0 |
| `shi` | 45.4 |
| `barall.2` | 46.8 |
| `barall` | 55.1 |
| `dliu.2` | 76.1 |
| `ma.2` | 82.3 |
| `duru` | 92.5 |
| `ma` | 111.8 |
| `dliu.4` | 205.9 |
| `dliu` | 245.3 |

## Notes

- `cplot` and the time series disagree on **problem**: `cplot` says `TPV30`, the time series say ``.
