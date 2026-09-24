# Session log 2026-09-24: autopilot (wei-lin)

Budget: 24 h from 07:30 CDT. Grant: serial PRs squash-merged to master (rule 25), with board, docs and rule text pushed direct. No tags were cut this milestone.

## M1: queue items 94, 96, 22a, 109, 78, 90, 95, 91, 105, 98 worked; 106, 93, 72, 97 and the source stamp added (07:30-10:20)

### PRs, serial, all CI green plus a victor-reyes audit
| PR | merge SHA | item | open→merged (min) | audit |
|---|---|---|---|---|
| #10 | `6831295` | 94: dropped off-fault stations named | 4.6 | MERGE |
| #11 | `74d02ee` | 109: python cell stdout/stderr kept | 19.9 | BLOCK→fixed→MERGE |
| #12 | `6c52408` | source-hash stamp on bin/eqdyna (owner-approved) | 10.5 | MERGE, advisory fixed |
| #13 | `5de359a` | 78 every workflow green + 106 nc/cell-set | 6.4 | MERGE |
| #14 | `faf7909` | 95: FATAL text from rank files | 12.3 | BLOCK→fixed→MERGE |
| #15 | `21e91bd` | 91: six perf defects (iris-vermeulen) | 10.9 | BLOCK→fixed→MERGE |
| #16 | `6f3d609` | 93: off-fault header/filename | 11.6 | MERGE, advisory fixed |
| #17 | `dde54b9` | 91g residual in run_scaling (row 115) | 5.6 | MERGE |

### Findings that changed what we know
- **94 is wider than tpv8.** Four gated cases drop off-fault stations: tpv8 4/15, tpv36 3/22, tpv37 3/22, tpv10 2/10. The cause is `setSurfaceStation`, which matches depth exactly while x and y snap to the nearest node. The drops are now loud. Refuse or snap is an owner decision (row 94, section B).
- **96 is not a defect.** `faulting.f90`'s `dtau` is never written; every assignment is to rsfNucleation's own local.
- **22a is verified WRONG for TPV29/30.** The writer outputs n-stress with compression positive (+182.07 MPa at tpv29 dp120). The spec says "Positive means extension", and the 2015 archive reads -181.44. TPV105-3D, whose spec says compression positive, is correct. The fix is a Fortran output change and the owner's call.
- **The Python port writes no station files at all** (new row 114).
- **78:** v5.16.0's failed publish run was tag-triggered, so no pre-tag check could see it. The fix therefore has two halves, a pre-tag check and a post-tag check.

### Contradictions and self-corrections
- Row 95's first classification (by grep for `abort` in the test file) missed `test_fault_geometry_guard`, whose abort is in the Fortran. The audit caught it and the guard was converted in the same PR.
- Row 91a's ledger wiring (iris) would have filed shard points as run_scaling rows. The audit caught it and it was removed before merge; row 91 stays open for it.
- Two board commands could never go red: row 105's `git log master..branch` after a squash-merge, and row 90's `grep -c` counting comments. Both are now in rule 14a (zofia, `f826a87`).

### Board and rules
zofia's pass (`8a92e32` board, `f826a87` rules) was salvaged from her worktree after a 429 killed her before the push; the fast suite passed before the push. This pass (conductor) closed 72, 93, 97 and 115, and updated 91.

### Interruptions
A session limit (429) hit around 09:1x; it reset at 09:50. Nothing was lost: every branch was pushed or committed, and zofia's commits were recovered from her worktree.
