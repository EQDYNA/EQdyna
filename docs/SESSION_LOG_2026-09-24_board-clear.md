# Session log 2026-09-24 -- "the 7, clear the board" (wei-lin)

Budget as read: close board rows 91, 120, 66, 100, 92, 109, 34 on commands that ran, or restate each so its own
command decides it. Owner rulings relayed mid-session added 94, 56, 76, 87, 44, 112, 73, 49, 63, root cleanup,
docs rule, README rewrite, docs site. Grant: merge + minor/patch tags on master of this repo per rule 25
(squash PRs for src/testsys/.github; docs/board direct). No tag, publish or Pages flip this session without the owner.

## Milestone 1 -- master red fixed; 66 and 100 decided; board carries the queue (22:40)

- **Master red (rule 3a).** CI run 36082084976 on b1c6d95 failed `test.tpv8 x python-jax` on profile_schema's 5%
  sum check (0.803 s of 15.398 s); physics green. The brief said CI was green, which was stale. Root cause, measured:
  (1) the 'write stations' phase (#23) had no bucket; (2) `Profile.phase` started its clock after the entry sync,
  whose first call imports jax (3.45 s here). PR #27 -> f186860: unaccounted 1.1245 s / 35.79 s (3.1%) ->
  0.0008 s / 24.91 s. Guard mutation-checked. The owner session landed #28 (8425260, `SUM_FLOOR_S=2.0`) in parallel.
  Contradiction recorded: the floor admits 13% on the 15 s cell although the attributed gap is 0.003%, and keep/drop
  is owner-held (row 122). Lesson compounded as rule 5b.
- **Row 66 closed.** tpv29, 32 ranks, (4,2,4) vs (8,1,4): max|diff| 9.708830e-15 over 3321 nodes (bound 1e-10); both
  runs vs reference <= 1.24e-14; the y split sat on the fault (7224 vs 3698 raw rows). The row's own command was
  unrunnable (run_e2e fixes ranks at `FORTRAN_RANKS`), so the case build was replicated by hand.
- **Row 100 recorded.** Serial tpv104 peak RSS: fortran 946632 kB, numpy 2398380 kB, jax 3848884 kB
  (was 0.95 / 4.16 / 10.41 GB). Load 33-36/64, so the wall clock is not comparable.
- **Board** f8332d1, **rules** 012150c (26, 1a, 5b), both by zofia-kaminska. Rule 15 "notes lead the README" vs
  rule 26 "3-line pointer" is reconcilable: `test_release_complete.py:178` accepts a pointer line.
- **Collision observed.** The owner session edited the main checkout mid-session (profile_schema floor), against
  rule 21b, and was invisible to `git worktree list` as an agent. Nothing reverted; it landed as #28.
- **In flight:** mira-volkov (94), anya-petrov (125 README), victor-reyes (#29 audit). Branches ready, serial
  behind #29: wei/row44-common-loud, wei/row56-no-device-auto, wei/row87-push-filter (live-checked: its push started
  no run), iris/perf-round2 (91/76/92-probe/112), wei/user-docs-style (lands with 125).
