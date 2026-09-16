# TPV30 — EQdyna v3.1, 50m, 2015-02-27

Submission label in the public record: `dliu.2` — "Dunyu Liu - Finite Element - EQdyna-50m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv30   mus=dliu.2
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
| first header line | `# TPV30` | `cplot` |
| author | Dunyu Liu | `cplot` |
| code | EQdyna | `cplot` |
| code_version | 3.1 | `cplot` |
| element_size | 50 m x 50 m on fault | `cplot` |
| date | 2/27/2015 | `cplot` |
| on-fault node spacing | 50 m | `cplot` coordinates |

## Files (37 — cplot, 24 on-fault, 12 off-fault)

| file | bytes |
|---|---|
| `body-030st-150dp000` | 605,779 |
| `body-030st000dp000` | 605,777 |
| `body-030st150dp000` | 605,778 |
| `body-200st-200dp000` | 605,780 |
| `body-200st000dp000` | 605,778 |
| `body-200st200dp000` | 605,779 |
| `body030st-150dp000` | 605,778 |
| `body030st000dp000` | 605,776 |
| `body030st150dp000` | 605,777 |
| `body200st-200dp000` | 605,779 |
| `body200st000dp000` | 605,777 |
| `body200st200dp000` | 605,780 |
| `cplot` | 11,884,729 |
| `faultst-042dp061` | 695,717 |
| `faultst-050dp000` | 695,715 |
| `faultst-050dp100` | 695,716 |
| `faultst-050dp160` | 695,714 |
| `faultst-089dp101` | 695,718 |
| `faultst-110dp014` | 695,718 |
| `faultst-150dp050` | 695,714 |
| `faultst-150dp120` | 695,717 |
| `faultst-180dp156` | 695,719 |
| `faultst000dp120` | 695,715 |
| `faultst043dp062` | 695,716 |
| `faultst046dp060` | 695,714 |
| `faultst050dp000` | 695,714 |
| `faultst050dp120` | 695,713 |
| `faultst051dp057` | 695,716 |
| `faultst059dp047` | 695,716 |
| `faultst090dp153` | 695,715 |
| `faultst093dp153` | 695,717 |
| `faultst100dp050` | 695,715 |
| `faultst100dp110` | 695,714 |
| `faultst150dp000` | 695,713 |
| `faultst150dp130` | 695,714 |
| `faultst167dp105` | 695,718 |
| `faultst170dp045` | 695,715 |

## Published standing

Published as column `dliu.2` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv30/metric_cvv1_tpv30_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 13 submissions: **120.6 ms**
- closest submission `dliu.3` at 76.1 ms; farthest `dliu.4` at 274.9 ms
- rank by that median among the 14 submissions: **12 of 14** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `dliu.3` | 76.1 |
| `kozdon.2` | 111.4 |
| `duru.2` | 113.8 |
| `shi.2` | 117.7 |
| `kozdon` | 117.9 |
| `shi` | 119.6 |
| `barall.2` | 120.6 |
| `barall` | 128.2 |
| `ma` | 146.5 |
| `ma.2` | 155.0 |
| `duru` | 164.9 |
| `dliu` | 171.7 |
| `dliu.4` | 274.9 |
