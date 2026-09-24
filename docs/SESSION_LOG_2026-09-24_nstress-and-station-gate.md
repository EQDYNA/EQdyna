# Session log 2026-09-24: n-stress sign (22a) and station gate (wei-lin)

Started 12:15 CDT at master `5cc0999`. Grant as stated by the coordinator: src/, testsys/ and .github/ changes go as serial PRs, each with CI green plus a victor-reyes audit, then squash-merge. Board, docs and evidence push direct. No tags this campaign. Board writer this session: the conductor. Rule 21e ("a single board line is done by the conductor directly") and rule 21c (one writer) were followed; the persona default of routing the board to zofia-kaminska was not, and that conflict is named here.

## M1: normal-stress sign per SCEC spec (board row 22a), 12:15-13:20

### PRs
| PR | merge SHA | what | open to merged | audit |
|---|---|---|---|---|
| #20 | `782463f` | per-case n-stress sign; `nsign` gate; meng station list moved onto the 400 m grid | 17:35:13Z to 17:54:08Z, 18.9 min | MERGE. 4 advisories: 2 fixed in the PR (unparseable name fails; read_bglobal refusals pinned), 2 routed to rows 116 and 117 |
| #21 | open | phantom all-zero `faultst000dp000.txt` from a rank with no on-fault station | -- | in progress |

Evidence pushed direct: `bbb1f6a` (the 23-cell sweep, 311 s at load ~15, plus the meng re-gate).

### Spec conventions (fetched into scratch/specs/ today)
TPV36_37_Description_v12, uploadTPV10_11_v3, uploadTPV103 and uploadTPV89 were fetched today; Signconvention3d was already there.
- Extension-positive: TPV29/30, TPV10/11, TPV36/37.
- Compression-positive: TPV103/104 and TPV105-3D.
- No n-stress field requested: TPV8, which takes the SCEC default, extension.

### Contradictions with what was assumed
- **TPV104's spec says compression-positive, not extension.** The brief assumed only TPV105-3D was compression-positive. uploadTPV103 (2008-10-18) says "Positive means compression". test.tpv104 now declares it.
- **Owner's published data (row 117):** TPV36/37 (2024) read +25.44 at dp060 against an extension-positive spec, so that submission is flipped. TPV104 (2016) reads -120 against a compression-positive spec, so it is also flipped. scec_archive was not touched.
- **meng2023a/cb had no real on-fault station.** The default list is on a 500 m grid; at dx=400, 12 of 13 stations were dropped silently. The one file written was a phantom (below).
- **Phantom station (PR #21).** Before the fix, `eqdyna3d.f90:268` forced `numOfOnFaultStCount=1` on a rank that matched no on-fault station. `output_onfault_st` then wrote an all-zero `faultst000dp000.txt` through `xonfs(:,0,0)`, out of bounds. That gave tpv104, tpv1053d and drv.a6 14 files for 13 stations. Where a real (0,0) station exists, the result depends on write order.
- **tpv29/30 write 13 of 24 on-fault stations** at 500 m, and say nothing about the rest (row 116).

### Gate evidence
- Real mutation: tpv8 re-run with `faultStNormalStressSign='compression'`. The `nsign` gate FAILED on 5 of 5 buried stations; the unmutated run passed.
- Fresh values: tpv29 dp120 -182.0745 MPa (it was +182.07).

### Board
- 22a moved to D, CLOSED.
- 114 updated: in flight, mira-volkov.
- 116 new: on-fault drops are silent.
- 117 new: published-submission sign, blocked on the owner.

### In flight at M1 close
- mira-volkov: row 114, the port of station output to python (`mira/row114-station-output`).
- iris-vermeulen: the station/nc gate (`iris/station-gate`, `iris/station-refs`).
