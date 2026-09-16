# TPV29 — EQdyna v3.1, 100m, 2015-02-27

Submission label in the public record: `dliu` — "Dunyu Liu - Finite Element - EQdyna-100m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv29   mus=dliu
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
| first header line | `# TPV29` | `cplot` |
| author | Dunyu Liu | `cplot` |
| code | EQdyna | `cplot` |
| code_version | 3.1 | `cplot` |
| element_size | 100 m x 100 m on fault | `cplot` |
| date | 2/27/2015 | `cplot` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (37 — cplot, 24 on-fault, 12 off-fault)

| file | bytes |
|---|---|
| `body-030st-150dp000` | 303,281 |
| `body-030st000dp000` | 303,279 |
| `body-030st150dp000` | 303,280 |
| `body-200st-200dp000` | 303,282 |
| `body-200st000dp000` | 303,280 |
| `body-200st200dp000` | 303,281 |
| `body030st-150dp000` | 303,280 |
| `body030st000dp000` | 303,278 |
| `body030st150dp000` | 303,279 |
| `body200st-200dp000` | 303,281 |
| `body200st000dp000` | 303,279 |
| `body200st200dp000` | 303,282 |
| `cplot` | 2,982,531 |
| `faultst-042dp061` | 348,219 |
| `faultst-050dp000` | 348,217 |
| `faultst-050dp100` | 348,218 |
| `faultst-050dp160` | 348,216 |
| `faultst-089dp101` | 348,220 |
| `faultst-110dp014` | 348,220 |
| `faultst-150dp050` | 348,216 |
| `faultst-150dp120` | 348,219 |
| `faultst-180dp156` | 348,221 |
| `faultst000dp120` | 348,217 |
| `faultst043dp062` | 348,218 |
| `faultst046dp060` | 348,216 |
| `faultst050dp000` | 348,216 |
| `faultst050dp120` | 348,215 |
| `faultst051dp057` | 348,218 |
| `faultst059dp047` | 348,218 |
| `faultst090dp153` | 348,217 |
| `faultst093dp153` | 348,219 |
| `faultst100dp050` | 348,217 |
| `faultst100dp110` | 348,216 |
| `faultst150dp000` | 348,215 |
| `faultst150dp130` | 348,216 |
| `faultst167dp105` | 348,220 |
| `faultst170dp045` | 348,217 |

## Published standing

Published as column `dliu` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv29/metric_cvv1_tpv29_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 14 submissions: **87.2 ms**
- closest submission `chen` at 46.0 ms; farthest `ma` at 111.6 ms
- rank by that median among the 15 submissions: **15 of 15** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `chen` | 46.0 |
| `dliu.2` | 63.8 |
| `barall` | 66.3 |
| `barall.2` | 78.6 |
| `bai` | 86.3 |
| `kozdon.2` | 86.4 |
| `kozdon` | 86.8 |
| `duru.2` | 87.6 |
| `shi.2` | 88.8 |
| `duru` | 89.7 |
| `shi` | 90.4 |
| `gabriel` | 91.5 |
| `ma.2` | 100.0 |
| `ma` | 111.6 |

## Notes

- This directory was seeded from the author's own retained 2015 upload (`cplot`, `faultst000dp120`, `body-030st000dp000`, placed here 2026-09-14) before any fetch, which is why it is the reference for the byte-match check below.

- 3 file(s) were already archived in this directory before the fetch and were left untouched (README rule: archival records are read-only). The freshly fetched bytes were compared against them instead — **all byte-identical**, which is the correctness check on this whole fetch path.
