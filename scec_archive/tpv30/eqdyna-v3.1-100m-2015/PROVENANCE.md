# TPV30 — EQdyna v3.1, 100m, 2015-02-27

Submission label in the public record: `dliu` — "Dunyu Liu - Finite Element - EQdyna-100m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv30   mus=dliu
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
| first header line | `# problem = LVFZ` | `cplot` |
| author | Benchun Duan | `cplot` |
| code | EQdyna | `cplot` |
| code_version | 3.1 | `cplot` |
| element_size | 100 m x 100 m on fault | `cplot` |
| date | 6-2-2012 | `cplot` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (37 — cplot, 24 on-fault, 12 off-fault)

| file | bytes |
|---|---|
| `body-030st-150dp000` | 303,280 |
| `body-030st000dp000` | 303,278 |
| `body-030st150dp000` | 303,279 |
| `body-200st-200dp000` | 303,281 |
| `body-200st000dp000` | 303,279 |
| `body-200st200dp000` | 303,280 |
| `body030st-150dp000` | 303,279 |
| `body030st000dp000` | 303,277 |
| `body030st150dp000` | 303,278 |
| `body200st-200dp000` | 303,280 |
| `body200st000dp000` | 303,278 |
| `body200st200dp000` | 303,281 |
| `cplot` | 2,982,542 |
| `faultst-042dp061` | 348,219 |
| `faultst-050dp000` | 348,217 |
| `faultst-050dp100` | 348,218 |
| `faultst-050dp160` | 348,216 |
| `faultst-089dp101` | 348,219 |
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

Published as column `dliu` of the cross-code metric page <https://strike.scec.org/cvws/metric_cvv1_u1/tpv30/metric_cvv1_tpv30_ac_0.html>, table *RMS difference in rupture time (milliseconds)*, file `cplot`.

- median RMS rupture-time difference from the other 13 submissions: **287.7 ms**
- closest submission `dliu.2` at 171.7 ms; farthest `ma.2` at 560.2 ms
- rank by that median among the 14 submissions: **14 of 14** (1 = closest to the group)

| other submission | RMS (ms) |
|---|---|
| `dliu.2` | 171.7 |
| `dliu.3` | 245.3 |
| `kozdon.2` | 278.2 |
| `duru.2` | 281.2 |
| `kozdon` | 284.1 |
| `shi.2` | 285.6 |
| `shi` | 287.7 |
| `barall.2` | 288.0 |
| `barall` | 295.4 |
| `duru` | 325.2 |
| `ma` | 360.8 |
| `dliu.4` | 429.8 |
| `ma.2` | 560.2 |

## Notes

- `cplot` and the time series disagree on **date**: `cplot` says `6-2-2012`, the time series say `2/27/2015  1:38:49`.

- `cplot` and the time series disagree on **author**: `cplot` says `Benchun Duan`, the time series say `Dunyu Liu`.

- `cplot` and the time series disagree on **problem**: `cplot` says `LVFZ`, the time series say ``.

- That `cplot` header is a leftover from a different run and does not describe this submission; the directory year is taken from the time-series timestamps, which all 36 files agree on.


Note: the cplot header carries a stale template from a prior LVFZ run
(author = Benchun Duan, EQdyna's original author and a group member; date
2012) while all 36 time series are this run's (Dunyu Liu, 2015-02-27, 100 m).
The data is the TPV30 run; only the header text is stale, and it originates
within the same group.
