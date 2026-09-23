
## run 2026-09-22 20:46  case test.tpv104  merged sha b183a1e  n_lo/n_hi 20/60  min-eff 0.95
box at start: 10/64 cpus busy over 0.50 -> 54 free; loadavg 10.79

-- 4 ranks -- cpus [2, 3, 4, 5] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 10.79
  arm base   302.25 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [312.3, 312.2, 312.3, 312.3]
           exchange [3.337, 1.974, 4.072, 1.674] ms/step  barrier wait [9.3, 11.28, 1.49, 6.55] ms/step  compile 0.6s  point wall 77s
  arm b      303.16 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [305.5, 305.6, 305.6, 305.6]
           exchange [12.886, 8.829, 8.44, 4.577] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 0.1s  point wall 72s
  arm a      298.83 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [308.2, 308.3, 308.2, 308.3]
           exchange [1.946, 2.091, 1.998, 4.124] ms/step  barrier wait [8.05, 6.44, 6.23, 1.6] ms/step  compile 0.6s  point wall 72s
  arm ab     279.59 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [293.1, 293.2, 293.2, 293.2]
           exchange [11.7, 5.816, 5.402, 5.449] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 71s
  VERDICT 4 ranks arm b    0.997x vs base  (302.25 -> 303.16 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm a    1.011x vs base  (302.25 -> 298.83 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm ab   1.081x vs base  (302.25 -> 279.59 ms/step; min eff base 1.00 arm 1.00)

-- 8 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 9] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 15.47
  arm base   186.29 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [203.3, 203.3, 203.4, 203.4, 203.4, 203.4, 203.4, 203.4]
           exchange [0.74, 1.288, 3.267, 3.428, 1.259, 2.412, 3.128, 2.198] ms/step  barrier wait [49.3, 48.02, 89.72, 92.48, 2.9, 3.96, 4.09, 1.77] ms/step  compile 1.0s  point wall 75s
  arm b      181.19 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [196.7, 197.4, 197.6, 197.6, 197.6, 197.7, 197.7, 197.7]
           exchange [62.453, 57.471, 86.567, 90.726, 4.583, 3.817, 3.973, 4.727] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 1.0s  point wall 79s
  arm a      183.58 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [200.8, 200.8, 200.8, 200.8, 200.8, 200.9, 200.8, 200.9]
           exchange [1.159, 1.483, 1.475, 0.849, 2.166, 1.469, 1.736, 2.491] ms/step  barrier wait [95.65, 84.98, 52.75, 53.28, 3.95, 2.62, 6.26, 1.52] ms/step  compile 1.0s  point wall 76s
  arm ab     176.30 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [191.1, 192.1, 192.7, 192.9, 192.8, 192.9, 193.0, 193.0]
           exchange [100.324, 55.464, 95.776, 52.17, 9.499, 6.83, 4.51, 2.715] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 1.0s  point wall 71s
  VERDICT 8 ranks arm b    1.028x vs base  (186.29 -> 181.19 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm a    1.015x vs base  (186.29 -> 183.58 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm ab   1.057x vs base  (186.29 -> 176.30 ms/step; min eff base 1.00 arm 1.00)

-- 12 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0, 10: 0.0, 11: 0.0, 12: 0.0, 13: 0.0} worst 0.00; box 9/64 busy -> 55 free; loadavg 17.99
  arm base   119.45 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [134.5, 134.5, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6]
           exchange [0.865, 1.294, 4.292, 4.083, 3.812, 2.158, 3.723, 3.728, 3.524, 2.09, 3.102, 3.273] ms/step  barrier wait [33.35, 34.56, 30.88, 7.26, 31.3, 6.39, 4.84, 6.11, 4.81, 3.36, 1.48, 1.74] ms/step  compile 0.9s  point wall 145s
  arm b      116.05 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [129.3, 130.0, 130.1, 130.1, 130.1, 130.1, 130.1, 130.1, 130.1, 130.1, 130.1, 130.1]
           exchange [35.726, 34.336, 9.275, 31.255, 8.986, 30.169, 4.879, 5.679, 4.47, 4.018, 8.271, 6.207] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 147s
  arm a      125.69 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [139.2, 139.2, 139.2, 139.2, 139.2, 139.3, 139.2, 139.3, 139.3, 139.3, 139.3, 139.3]
           exchange [1.45, 1.232, 1.14, 0.927, 1.231, 1.782, 1.624, 1.7, 2.189, 1.58, 1.238, 1.855] ms/step  barrier wait [33.39, 29.42, 29.93, 33.92, 14.09, 4.98, 8.62, 2.93, 3.46, 4.76, 9.13, 8.62] ms/step  compile 0.8s  point wall 141s
  arm ab     116.50 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [127.4, 128.1, 128.6, 129.1, 129.8, 129.6, 130.3, 130.4, 130.4, 130.4, 130.4, 130.4]
           exchange [46.389, 41.477, 35.413, 30.352, 34.555, 34.746, 32.294, 8.444, 32.081, 3.752, 3.894, 2.047] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 143s
  VERDICT 12 ranks arm b    1.029x vs base  (119.45 -> 116.05 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm a    0.950x vs base  (119.45 -> 125.69 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm ab   1.025x vs base  (119.45 -> 116.50 ms/step; min eff base 1.00 arm 1.00)

snapshot /home/utig5/dliu/EQdyna/.claude/worktrees/agent-ab9486378831c8e09/docs/perf_snapshots/jaxmpi_ab_2026-09-22_204652.json
12 ledger row(s) appended to docs/perf_ledger.jsonl

## run 2026-09-22 21:06  case test.tpv104  merged sha b183a1e  n_lo/n_hi 20/60  min-eff 0.95
box at start: 10/64 cpus busy over 0.50 -> 54 free; loadavg 18.33

-- 4 ranks -- cpus [2, 4, 5, 6] busy {2: 0.0, 4: 0.0, 5: 0.0, 6: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 18.33

## STEP 0 -- STALE-MERGE VERDICT (before any measurement)

Branch `mira/jaxmpi-step-2026-09-22` (0ab3fa7) branched at 596f004; master is
44afd4a. Merged master into it in this worktree -> `mira/jaxmpi-merged-2026-09-22`.

MERGE WAS CLEAN on both source files. The only conflict was
`docs/perf_ledger.jsonl`, which is append-only JSONL: resolved by keeping BOTH
sides (153 lines, every line re-parsed as JSON). `src/python/eqdyna/driver.py`
auto-merged; `backend.py` came through from master untouched.

DIFF VERDICT, `git show master:src/python/eqdyna/driver.py` vs the merged copy:
ONLY the halo/barrier changes appear. Three hunks, all of them the branch's:
  - the per-step `comm.Barrier()` put under `if prof:`
  - `barrier_in_step=bool(prof)` added to the report
  - `barrier_in_step=True` added to the profiled report
NO REVERTED LINES. `B.enable_compilation_cache(subdir='rank%d' % rank)` (item
61, the per-rank XLA cache-directory fix) is present at driver.py:448, and
`backend.py` is byte-identical to master (0 diff lines). The A/B harness
re-asserts the cache-fix string before it will stage any arm.

## METHOD

Four arms, a 2x2 over the two independent changes, each arm a PAIR OF FILE
VERSIONS taken from git and copied into src/python/eqdyna:
  base = master driver + master MPI4NodalQuant  (Sendrecv, unconditional barrier)
  b    = merged driver + master MPI4NodalQuant  (barrier under prof)
  a    = master driver + merged MPI4NodalQuant  (non-blocking halo)
  ab   = merged both                            (the branch)
Staging verified by sha1 and by three content predicates per arm before launch.

(b) IS INVISIBLE UNDER THE STEP PROFILE -- with `prof` on, arm b runs the
barrier exactly as base does. So the metric is the perf tier's, not the step
profile: per-step by DIFFERENCE over 20 vs 60 steps, over each rank's OWN
solve clock (max over ranks), compile/setup subtracted as the common term.
`--warmup` discards one n_lo run first so both timed runs hit a warm XLA cache.
No total wall clock is reported anywhere.

ONE CPU SET PER RANK COUNT, chosen once (least-loaded, cpus 0 and 1 excluded by
measurement) and reused by all four arms, so placement cannot masquerade as the
effect. Every point carries each rank's EFFECTIVE_CORES measured AFTER LAUNCH
from the rank's own stdout; a point with any rank below 0.95 is REJECTED.

## PASS 1 (2026-09-22 20:46, arm order base,b,a,ab)

box at start 10/64 cpus busy over 0.50 -> 54 free, loadavg 10.79.
EVERY ONE OF THE 12 POINTS MEASURED min eff 1.00 AND WAS ACCEPTED -- the
window really was open; nothing was rejected.

 ranks   cpus            base ms/step    b        a        ab
   4     2-5               302.25      0.997x   1.011x   1.081x
   8     2-9               186.29      1.028x   1.015x   1.057x
  12     2-13              119.45      1.029x   0.950x   1.025x

The 0.48x-at-4-ranks shape DID NOT REPRODUCE: arm a at 4 ranks is 1.011x.
Neither did ~2.2x at 8/12: arm a is 1.015x and 0.950x there.
ab is not a+b at any rank count (4 ranks: 1.081 vs 1.011*0.997=1.008), which
says run-to-run spread is the same size as the effect. Replicate passes 2 and 3
launched with SHUFFLED arm order to break the arm/position confound.

Direct barrier cost, read off arm base where t_wait is the real unconditional
barrier: 4 ranks [9.3, 11.28, 1.49, 6.55] ms/step of a 302 ms step -- mean 7.2
ms, 2.4%. That is the CEILING on what (b) can buy, and removing the barrier
does not delete the imbalance it was measuring, it only moves it into the
exchange (arm b exchange rose to [12.9, 8.8, 8.4, 4.6] from base's [3.3, 2.0,
4.1, 1.7]).
  arm ab     283.32 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [301.8, 301.8, 301.9, 301.9]
           exchange [47.621, 64.985, 2.401, 1.451] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 1.1s  point wall 72s
  arm a      294.47 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [308.1, 308.1, 308.3, 308.3]
           exchange [1.078, 0.803, 1.154, 1.414] ms/step  barrier wait [63.39, 48.29, 1.79, 0.99] ms/step  compile 0.8s  point wall 71s
  arm b      282.49 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [296.3, 296.3, 296.4, 296.5]
           exchange [41.054, 57.8, 2.47, 2.548] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 70s
  arm base   286.81 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [298.4, 298.5, 298.6, 298.6]
           exchange [0.913, 1.477, 0.81, 1.315] ms/step  barrier wait [58.03, 39.48, 2.34, 0.61] ms/step  compile 0.7s  point wall 70s
  VERDICT 4 ranks arm ab   1.012x vs base  (286.81 -> 283.32 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm a    0.974x vs base  (286.81 -> 294.47 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm b    1.015x vs base  (286.81 -> 282.49 ms/step; min eff base 1.00 arm 1.00)

-- 8 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 9] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 14.92
  arm ab     181.59 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [194.7, 195.2, 196.0, 196.1, 196.1, 196.2, 196.3, 196.3]
           exchange [103.012, 60.879, 58.076, 95.389, 4.953, 2.632, 3.799, 3.82] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.9s  point wall 70s
  arm a      191.64 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [202.5, 202.5, 202.5, 202.5, 202.5, 202.5, 202.5, 202.5]
           exchange [0.901, 1.737, 1.428, 1.132, 1.788, 2.027, 2.169, 2.12] ms/step  barrier wait [97.22, 90.85, 49.44, 52.77, 3.94, 2.14, 2.29, 2.73] ms/step  compile 0.7s  point wall 71s
  arm b      180.69 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [195.5, 196.3, 196.5, 196.5, 196.5, 196.5, 196.5, 196.5]
           exchange [59.609, 53.873, 84.659, 85.291, 5.344, 3.746, 4.221, 5.041] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.9s  point wall 70s
  arm base   187.11 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [204.5, 204.5, 204.6, 204.6, 204.6, 204.6, 204.6, 204.6]
           exchange [1.094, 1.592, 3.379, 3.211, 3.115, 2.256, 1.805, 2.647] ms/step  barrier wait [51.93, 50.69, 94.82, 89.22, 3.18, 2.33, 3.5, 2.96] ms/step  compile 1.0s  point wall 76s
  VERDICT 8 ranks arm ab   1.030x vs base  (187.11 -> 181.59 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm a    0.976x vs base  (187.11 -> 191.64 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm b    1.036x vs base  (187.11 -> 180.69 ms/step; min eff base 1.00 arm 1.00)

-- 12 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0, 10: 0.0, 11: 0.0, 12: 0.0, 14: 0.0} worst 0.00; box 11/64 busy -> 53 free; loadavg 19.79
  arm ab     112.57 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [125.5, 126.2, 126.7, 127.2, 127.9, 128.4, 127.8, 128.5, 128.6, 128.6, 128.6, 128.6]
           exchange [43.824, 41.158, 33.448, 32.419, 34.626, 33.055, 31.696, 10.788, 4.287, 30.14, 2.278, 2.997] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 1.0s  point wall 141s
  arm a      124.06 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1, 137.1]
           exchange [0.931, 1.352, 1.725, 1.324, 1.009, 1.456, 1.689, 1.23, 1.712, 1.14, 2.161, 1.519] ms/step  barrier wait [32.18, 31.83, 5.67, 18.26, 46.97, 44.48, 2.88, 15.39, 2.23, 14.83, 1.7, 14.53] ms/step  compile 0.8s  point wall 155s
  arm b      121.10 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [133.0, 133.5, 133.6, 133.7, 133.7, 133.7, 133.7, 133.7, 133.7, 133.7, 133.7, 133.7]
           exchange [35.804, 34.18, 8.343, 18.094, 44.82, 46.074, 4.684, 4.29, 5.065, 14.414, 16.254, 15.983] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 158s
  arm base   121.52 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [137.0, 137.0, 137.0, 137.0, 137.0, 137.0, 137.2, 137.0, 137.0, 137.1, 137.1, 137.0]
           exchange [3.878, 0.844, 1.27, 1.463, 3.902, 3.706, 1.857, 3.502, 3.171, 2.681, 2.94, 3.48] ms/step  barrier wait [16.1, 31.89, 31.48, 6.1, 39.16, 44.86, 2.69, 12.37, 13.35, 1.62, 1.53, 12.16] ms/step  compile 0.9s  point wall 161s
  VERDICT 12 ranks arm ab   1.080x vs base  (121.52 -> 112.57 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm a    0.980x vs base  (121.52 -> 124.06 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm b    1.003x vs base  (121.52 -> 121.10 ms/step; min eff base 1.00 arm 1.00)

snapshot /home/utig5/dliu/EQdyna/.claude/worktrees/agent-ab9486378831c8e09/docs/perf_snapshots/jaxmpi_ab_2026-09-22_210655.json
12 ledger row(s) appended to docs/perf_ledger.jsonl

## run 2026-09-22 21:26  case test.tpv104  merged sha b183a1e  n_lo/n_hi 20/60  min-eff 0.95
box at start: 10/64 cpus busy over 0.50 -> 54 free; loadavg 21.31

-- 4 ranks -- cpus [2, 3, 4, 5] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 21.31
  arm a      252.52 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [278.3, 278.3, 278.3, 278.3]
           exchange [1.474, 2.553, 3.196, 2.077] ms/step  barrier wait [5.36, 5.19, 2.96, 3.82] ms/step  compile 1.5s  point wall 71s
  arm ab     251.82 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [280.0, 279.9, 279.9, 279.9]
           exchange [14.832, 4.769, 4.273, 11.242] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 1.7s  point wall 70s
  arm base   294.29 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [307.5, 307.6, 307.6, 307.6]
           exchange [1.561, 4.514, 2.69, 4.933] ms/step  barrier wait [8.7, 6.6, 5.7, 3.27] ms/step  compile 0.8s  point wall 72s
  arm b      293.89 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0]   rank ms [306.7, 306.7, 306.7, 306.8]
           exchange [17.033, 8.208, 8.462, 4.704] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 72s
  VERDICT 4 ranks arm a    1.165x vs base  (294.29 -> 252.52 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm ab   1.169x vs base  (294.29 -> 251.82 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 4 ranks arm b    1.001x vs base  (294.29 -> 293.89 ms/step; min eff base 1.00 arm 1.00)

-- 8 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 10] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 10: 0.0} worst 0.00; box 10/64 busy -> 54 free; loadavg 14.56
  arm a      185.63 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [201.5, 201.5, 201.5, 201.5, 201.6, 201.6, 201.6, 201.6]
           exchange [0.964, 1.076, 1.691, 1.673, 1.774, 1.951, 2.457, 2.054] ms/step  barrier wait [100.74, 50.85, 94.93, 49.31, 3.2, 2.88, 2.48, 2.88] ms/step  compile 1.0s  point wall 71s
  arm ab     177.09 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [191.8, 192.7, 193.2, 193.5, 193.3, 193.5, 193.5, 193.5]
           exchange [99.584, 57.747, 95.709, 54.334, 3.135, 2.614, 8.696, 2.742] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 1.0s  point wall 71s
  arm base   188.68 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [204.2, 204.3, 204.3, 204.3, 204.4, 204.3, 204.4, 204.4]
           exchange [0.845, 1.683, 3.783, 3.544, 3.511, 1.627, 2.805, 2.602] ms/step  barrier wait [53.88, 51.39, 102.89, 97.0, 3.14, 3.64, 2.44, 1.98] ms/step  compile 0.9s  point wall 72s
  arm b      182.94 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [195.6, 196.4, 196.5, 196.5, 196.5, 196.6, 196.6, 196.6]
           exchange [64.711, 58.943, 98.526, 92.135, 3.814, 4.484, 3.908, 3.847] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.8s  point wall 70s
  VERDICT 8 ranks arm a    1.016x vs base  (188.68 -> 185.63 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm ab   1.065x vs base  (188.68 -> 177.09 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 8 ranks arm b    1.031x vs base  (188.68 -> 182.94 ms/step; min eff base 1.00 arm 1.00)

-- 12 ranks -- cpus [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0, 10: 0.0, 11: 0.0, 12: 0.0, 13: 0.0} worst 0.00; box 9/64 busy -> 55 free; loadavg 19.85
  arm a      122.47 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [136.4, 136.4, 136.5, 136.4, 136.5, 136.5, 136.5, 136.5, 136.5, 136.5, 136.5, 136.5]
           exchange [1.262, 0.934, 1.77, 1.291, 1.538, 0.829, 1.225, 1.447, 1.252, 1.644, 1.682, 2.144] ms/step  barrier wait [30.72, 30.5, 6.0, 12.61, 39.02, 34.03, 9.46, 10.04, 9.49, 2.3, 2.48, 1.89] ms/step  compile 0.8s  point wall 139s
  arm ab     113.51 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [126.2, 126.9, 127.4, 128.4, 127.9, 128.9, 128.9, 128.4, 129.0, 129.0, 129.1, 129.1]
           exchange [40.161, 39.541, 33.931, 34.792, 28.444, 32.256, 8.041, 29.765, 25.259, 6.453, 2.87, 2.407] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 0.9s  point wall 138s
  arm base   119.63 ms/step   min eff 1.00 ACCEPTED
           eff [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [134.5, 134.5, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6, 134.6]
           exchange [0.947, 0.771, 3.952, 3.426, 1.66, 3.739, 1.616, 3.236, 3.608, 2.962, 2.914, 3.206] ms/step  barrier wait [30.82, 30.69, 30.49, 31.42, 7.42, 9.85, 3.91, 9.43, 6.47, 1.46, 0.98, 7.55] ms/step  compile 0.9s  point wall 131s
  arm b      113.09 ms/step   min eff 0.99 ACCEPTED
           eff [0.99, 0.99, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]   rank ms [126.9, 127.5, 127.9, 128.0, 128.9, 128.4, 128.9, 128.9, 128.9, 128.9, 128.9, 128.9]
           exchange [33.932, 32.382, 28.35, 19.374, 28.114, 19.19, 28.292, 6.608, 4.356, 2.626, 2.946, 17.113] ms/step  barrier wait [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ms/step  compile 1.0s  point wall 146s
  VERDICT 12 ranks arm a    0.977x vs base  (119.63 -> 122.47 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm ab   1.054x vs base  (119.63 -> 113.51 ms/step; min eff base 1.00 arm 1.00)
  VERDICT 12 ranks arm b    1.058x vs base  (119.63 -> 113.09 ms/step; min eff base 1.00 arm 0.99)

snapshot /home/utig5/dliu/EQdyna/.claude/worktrees/agent-ab9486378831c8e09/docs/perf_snapshots/jaxmpi_ab_2026-09-22_212647.json
12 ledger row(s) appended to docs/perf_ledger.jsonl

## run 2026-09-22 21:47  case test.tpv104  merged sha b183a1e  n_lo/n_hi 20/60  min-eff 0.95
box at start: 9/64 cpus busy over 0.50 -> 55 free; loadavg 11.79

-- 16 ranks -- cpus [2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19] busy {2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 8: 0.0, 9: 0.0, 11: 0.0, 12: 0.0, 13: 0.0, 14: 0.0, 15: 0.0, 16: 0.0, 17: 0.0, 18: 0.0, 19: 0.0} worst 0.00; box 9/64 busy -> 55 free; loadavg 11.79

## AGGREGATE OVER THREE REPETITIONS AT 4 / 8 / 12 RANKS (test.tpv104)

Arm order was SHUFFLED between repetitions (base,b,a,ab | ab,a,b,base |
a,ab,base,b) so that "which arm ran first" cannot masquerade as the effect.
Every speedup below is PAIRED: divided by the base arm of the SAME repetition
on the SAME cpu set.

ALL 36 POINTS MEASURED min EFFECTIVE_CORES >= 0.99. Nothing was rejected.
The window was genuinely open: box 9-11 of 64 cpus busy over 0.50 throughout.

NOISE FLOOR, measured -- base-arm spread across the three repetitions:
  4 ranks  5.2%   8 ranks  1.3%   12 ranks  1.7%
At 4 ranks the noise is LARGER than anything either change does. That alone
disqualifies 4 ranks as a place to decide this.

PAIRED SPEEDUP vs base (mean +- sd over three repetitions):
  ranks   b (barrier)        a (non-blocking)    ab (both)
    4     1.005 +- 0.008     1.050 +- 0.083      1.087 +- 0.064
    8     1.032 +- 0.003     1.003 +- 0.019      1.051 +- 0.015
   12     1.030 +- 0.022     0.969 +- 0.013      1.053 +- 0.022

THE 0.48x-AT-4-RANKS SHAPE DID NOT REPRODUCE. Arm a at 4 ranks measured
1.011, 0.974, 1.165 over three repetitions -- mean 1.050, sd 0.083, straddling
1.0 with a spread that the base arm alone also shows (5.2%). There is no
4-rank penalty here to special-case, and equally no 4-rank benefit.
NEITHER DID ~2.2x AT 8 AND 12. Arm a there is 1.003 and 0.969.

(b) IS THE ONE THAT STANDS. At 8 ranks it is 1.028 / 1.036 / 1.031 -- three
independent repetitions inside 0.8% of each other, a +3.2% effect against a
1.3% noise floor. At 12 ranks +3.0%. At 4 ranks it is 1.005, i.e. neutral and
not a loss, which is what the gated cell needs.

(a) DOES NOT STAND ALONE. 1.003 at 8 ranks is nothing, and at 12 ranks it is
0.969 -- BELOW 1.0 in all three repetitions (0.950, 0.980, 0.977) against a
1.7% noise floor. On its own evidence the non-blocking exchange is a small
REGRESSION at 12 ranks, not a speedup.

INTERACTION ab/(a*b) = 1.032 / 1.016 / 1.055 at 4 / 8 / 12: the pair is
SUPERADDITIVE, so ab (1.051, 1.053) beats b alone (1.032, 1.030) by about 2%
-- around the noise floor. The mechanism the timing columns support: with the
unconditional barrier in place every rank enters the exchange already
synchronised, so there is no jitter for Irecv/Isend to overlap; remove the
barrier and there is. Note the attribution columns invert between arms --
base/a put 7-40 ms/step in `wait` and 1-3 ms in `exchange`, while b/ab report
`wait` 0.00 by construction and 8-41 ms in `exchange`. The exchange column in
b/ab is therefore mostly WAITING FOR THE SLOWEST NEIGHBOUR, not transfer
cost, exactly as the branch's own `barrier_in_step` note warns.

Extending to 16 and 32 ranks next -- the author's wave argument was measured
at 32, and 12 is the largest rank count above.
  arm base   552.02 ms/step   min eff 0.25 REJECTED (< 0.95)
           eff [0.98, 0.97, 0.99, 0.99, 0.99, 0.98, 0.99, 1.0, 0.99, 0.98, 0.99, 0.94, 0.25, 0.26, 0.26, 0.26]   rank ms [414.5, 414.4, 414.5, 414.8, 414.6, 414.8, 414.5, 414.8, 414.8, 414.6, 414.6, 414.6, 414.6, 415.1, 415.0, 414.8]
           exchange [5.466, 6.629, 7.836, 19.151, 11.474, 16.021, 9.18, 14.103, 14.972, 8.782, 10.276, 9.817, 5.784, 56.475, 53.369, 12.381] ms/step  barrier wait [310.8, 302.61, 297.44, 314.87, 295.09, 313.67, 304.22, 290.82, 294.6, 278.42, 275.77, 280.02, 28.72, 49.51, 64.31, 64.09] ms/step  compile -8.2s  point wall 481s

## 16 AND 32 RANKS: WINDOW CLOSED, POINT REJECTED, SWEEP STOPPED

Attempted 16 and 32 ranks to test the author's wave argument at the rank count
it was argued at (32). The 16-rank base point was REJECTED by the filter:

  min EFFECTIVE_CORES 0.25. Ranks 0-11 measured 0.94-1.00; ranks 12-15, on
  cpus 16,17,18,19, measured 0.25 / 0.26 / 0.26 / 0.26 -- and all four of
  those cpus had read /proc/stat busy 0.00 immediately before launch. That is
  the same "reads idle, delivers a quarter of a core" pathology already on
  record for cpus 0 and 1, on four more cpus. The point also returned
  compile -8.2 s, i.e. the difference was invalid as well as contended.

WHOSE LOAD. Not a foreign tenant: a SIBLING SESSION OF OURS started an e2e
sweep at 21:50 (three `python3 -m eqdyna .../test/test.{meng2023a,tpv36,tpv10}
.python-numpy --backend numpy` out of /home/utig5/dliu/EQdyna/test/), taking
the box from loadavg ~11 to ~35. I stopped MY 16/32 sweep rather than compete
with it, and killed only processes I had started.

THE 4 / 8 / 12 DATA IS UNAFFECTED: all three repetitions completed by ~21:47,
before those jobs began, and their base-arm spread (1.3% at 8 ranks, 1.7% at
12) confirms the box was quiet while they ran.

SO THE WAVE ARGUMENT IS UNTESTED ABOVE 12 RANKS by this campaign. What is
tested is that it does not pay at 4, 8 or 12.

## PARITY

`testsys/e2e/run_e2e.py --cases test.tpv8 --backends python-jax-mpi` on the
merged tree: SUCCESS, max|diff| = 1.220703e-10 against a bound of 1.0e-08,
1891 fault nodes compared, 15.8 s, peak RSS 5.34 GB. No reference and no bound
was touched. The cell did NOT wedge -- which is the positive confirmation that
the per-rank XLA cache fix survived the merge and is compatible with both of
the branch's changes.

## RECOMMENDATION

LAND (b), the barrier. It is first a correctness-of-the-metric fix -- an
unconditional per-step global rendezvous was running in PRODUCTION purely to
serve a profiling attribution -- and it measures +3.2% at 8 ranks (sd 0.3%
over three repetitions, against a 1.3% noise floor), +3.0% at 12, and 1.005
at 4, i.e. neutral where the gate lives. Nothing about it needs a rank-count
branch.

DO NOT LAND (a) ON ITS OWN EVIDENCE. Non-blocking exchange measured 1.003 at
8 ranks and 0.969 at 12 -- below 1.0 in all three repetitions at 12 -- and at
4 ranks its spread (sd 0.083) exceeds any effect. It pays only in combination
with (b), and that combination beats (b) alone by ~2%, which is at the noise
floor. The case for taking it anyway is not performance: Irecvs posted before
any Isend cannot deadlock for ANY neighbour graph, whereas ascending-order
Sendrecv is deadlock-free only because the slab decomposition happens to give
two neighbours. That is a robustness argument, and it is the owner's to weigh.

THE 0.48x-AT-4-RANKS CAVEAT DOES NOT REPRODUCE. Three repetitions at 4 ranks
gave arm a at 1.011, 0.974, 1.165. There is no 4-rank cliff in this data and
therefore nothing to special-case. Equally, the ~2.2x at 8 and 12 does not
reproduce. The branch author's decision not to claim a speedup was correct;
the correction this measurement makes is that it was correct for a different
reason -- not because 4 ranks is pathological, but because (a) does nothing at
any rank count measured.
