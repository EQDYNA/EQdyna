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
- **meng2023a/cb had no real on-fault station.** The default list is on a 500 m grid; at dx=400, all 13 stations were dropped silently. The one file written was a phantom (below).
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

## M2: station time series and jax nc in the gate (13:20-17:10)

### Interruption
At about 14:0x a session limit (429) killed both specialists. mira-volkov (row 114) had pushed 3 commits but not her parity sweep. iris-vermeulen (the gate) had nothing committed; her compare.py/matrix.py diff and staged refs existed only in her worktree. The limit reset at 14:50. Following the persona rule "after a limit hit, stop dispatching and work directly", the conductor salvaged both and finished the work. Only the victor-reyes per-PR audits were dispatched after that.

### PRs (serial; CI green and a victor-reyes audit on each)
| PR | merge SHA | what | open to merged | audit |
|---|---|---|---|---|
| #22 | `9444d9a` | 116: on-fault drops are loud | 14.2 min | MERGE; 3 lows fixed |
| #23 | `a338ee9` | 114: the port writes station files | 33.9 min | BLOCK (silent drops; parity JSON uncommitted), fixed, MERGE |
| #24 | `1204ce4` | station series and jax nc in the gate | 28.7 min | BLOCK (a NaN passed), fixed, MERGE |

Direct pushes:
- Station references: 11 rule-7 commits ending `b997d06`, reason "owner gate design".
- Parity evidence: `6439a37`.

### Stations and thresholds
These are per case in `testsys/matrix.py` (GATE_STATIONS, STATION_BOUND) and in board row 119. The error normalization divides each station's worst difference by that quantity's case-level maximum, clamped at a floor of 1e-6. Each bound is the worst observed fortran-vs-jax error times about 121.5, rounded up to a bound already in use.

### Contradictions and self-corrections
- mira's notes explained Fortran tpv104's extra `faultst000dp000` as a rank-boundary duplicate. PR #21 found the real cause: an empty rank's station count was forced to 1.
- iris measured "the port writes zero station files" for meng/tpv36/tpv37. That came from her pre-22a build and pre-wedge-fix branch. On the merged port every case matches Fortran's file set.
- M1's log said meng dropped "12 of 13" stations. It was 13 of 13; the one file it wrote was the phantom (corrected in PR #22).
- The first everyday sweep of the gate died with SIGSEGV (exit -11). jax cells now open netCDF4/HDF5 from concurrent threads, which is not thread-safe here. `_NC_LOCK` serializes the opens (papercut logged).
- The owner asked for the jax nc "at the case bound". The one existing nc comparison (THRESHOLD 1e-3) was kept under rule 5, and the measured diffs went to row 121 for an owner decision.

### Sweep wall time (`python3 testsys/run.py e2e`, 23 cells)
| state | wall | 1-min load |
|---|---|---|
| 5cc0999 (frt-only) | 311 s | ~15 |
| 9444d9a | 312 s | 27 |
| 1204ce4 tree (final gate) | 329 s | 15-36 |

### Board
- Closed: 114, 116, 119.
- New: 120 (jax-mpi writes no stations) and 121 (nc tolerance, owner decision).
