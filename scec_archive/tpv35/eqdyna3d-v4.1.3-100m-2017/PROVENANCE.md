# TPV35 — EQdyna3D v4.1.3, 100m, 2017-01-17

Submission label in the public record: `dliu` — "Dunyu Liu - Finite Element - EQdyna3Dv4.1.3 - 100m".

## Source

Public Area of the SCEC/USGS Code Verification Web Site comparison CGI. No
credentials and no cookies; all state travels in the POST body.

    POST https://strike.scec.org/cvws/cgi-bin/cvws.cgi
      o=1005   m=tpv35   mus=dliu
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
| first header line | `# This is the file header:` | `cplot` |
| author | D.Liu | `cplot` |
| code | EQdyna3D | `cplot` |
| code_version | 4.1.3 | `cplot` |
| element_size | 100.0 m | `cplot` |
| date | 2017/01/17 | `cplot` |
| on-fault node spacing | 100 m | `cplot` coordinates |

## Files (50 — cplot, 7 on-fault, 42 named station)

| file | bytes |
|---|---|
| `4064_donna` | 272,969 |
| `4065_eades` | 272,969 |
| `4066_froel` | 272,969 |
| `4067_gold` | 272,969 |
| `4069_jack` | 272,969 |
| `4070_joaqu` | 272,969 |
| `4071_middl` | 272,969 |
| `4072_redh` | 272,969 |
| `4074_viney` | 272,969 |
| `4097_scn` | 272,969 |
| `4098_c01` | 272,969 |
| `4099_tm2` | 272,969 |
| `4100_c02` | 272,969 |
| `4101_tm3` | 272,969 |
| `4102_c03` | 272,969 |
| `4103_c04` | 272,969 |
| `4104_c4a` | 272,969 |
| `4107_cow` | 272,969 |
| `4108_coh` | 272,969 |
| `4109_z04` | 272,969 |
| `4110_z06` | 272,969 |
| `4111_z07` | 272,969 |
| `4112_z08` | 272,969 |
| `4113_z09` | 272,969 |
| `4114_z11` | 272,969 |
| `4115_prk` | 272,969 |
| `4117_z15` | 272,969 |
| `4118_pg1` | 272,969 |
| `4119_gh2` | 272,969 |
| `4121_gh3` | 272,969 |
| `4122_pg3` | 272,969 |
| `4124_pg5` | 272,969 |
| `4126_sc1` | 272,969 |
| `4127_sc2` | 272,969 |
| `4128_sc3` | 272,969 |
| `4129_36510` | 272,969 |
| `4131_vc1` | 272,969 |
| `4132_pgd` | 272,969 |
| `4133_vc2` | 272,969 |
| `4135_vc4` | 272,969 |
| `4136_vc5` | 272,969 |
| `8486_nphob` | 272,969 |
| `cplot` | 2,705,030 |
| `faultst-050dp081` | 313,419 |
| `faultst-100dp081` | 313,419 |
| `faultst-150dp081` | 313,419 |
| `faultst-200dp081` | 313,419 |
| `faultst-250dp081` | 313,419 |
| `faultst000dp081` | 313,419 |
| `faultst050dp081` | 313,419 |

## Published standing

No public metric page for `tpv35` — <https://strike.scec.org/cvws/metric_cvv1_u1/tpv35/metric_cvv1_tpv35_ac_0.html> returns HTTP 404, so SCEC publishes no cross-code table for this listing.

## Notes

- `cplot` and the time series disagree on **code**: `cplot` says `EQdyna3D`, the time series say `EQdyna3d`.

- `cplot` and the time series disagree on **problem**: `cplot` says `TPV35`, the time series say ``.
