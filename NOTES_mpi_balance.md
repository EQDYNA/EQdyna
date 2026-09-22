# jax-MPI element-balance campaign (item 43 follow-on)

Base SHA: 09da166 (worktree HEAD == origin/master; re-sync checkout was a no-op).
Worktree: /home/utig5/dliu/EQdyna/.claude/worktrees/agent-a40ab3dc8bab919fc

## M1 (2026-09-22) BEFORE: the partition table, reproduced

`testsys/perf/scaling_case/test.tpv104` built via run_scaling.build_py_case;
S built by eqdyna3d.build_solver_state; cuts taken from MQ._cuts on the real
elemType array. No solver run needed -- the partition is pure index algebra.

    E=735000  Ei=516096  Ep=218904  nftnd=2701  N=763537 NEQ=3818583
    PML element indices: [0, 734999]; INTERIOR element indices: [31956, 703493]
      -> elements 0..31955 are ALL PML (31956 of them), and 703494..734999 too
         (31506).  These are the two x-face PML slabs, contiguous at the ends
         of the mesh generator's nested loop order.

    nranks=4   Ei [104112, 153936, 153216, 104832]              spread 1.479x, 0 zeros
    nranks=8   Ei [27153, 76959, 76959, 76977, 76608, 76608, 76833, 27999]
                                                                spread 2.835x, 0 zeros
    nranks=16  Ei [0, 27153, 38052, ... , 27999, 0]              2 zeros
    nranks=32  Ei [0, 0, 8064, 19089, ..., 8124, 0, 0]           4 zeros
               E  [12217, 12217, 17593, 24942, ... , 17633, 12217, 12216]
               E spread 2.086x (25478/12216)

Matches the brief's measured defect (4 zero-Ei ranks of 32, others ~19k) and
the 4xA100 confirmation (104112/104832/153216/153936, 1.479x) exactly.

## M1a ARITHMETIC RESULT -- contiguity is what forces the zero-Ei rank

With a CONTIGUOUS partition, rank 0's interval is [0, e1). Every element below
31956 is PML, so Ei(rank 0) > 0 requires e1 > 31956, i.e. rank 0 must hold at
least 31956 PML elements. Equal-weight cutting gives rank 0 a weight share of
total/32, so Ei(rank 0) > 0 at 32 ranks requires

    31956*w < (516096 + 218904*w)/32
    803688*w < 516096
    w < 0.642

i.e. a PML element would have to cost LESS than 0.64 of an interior element.
calcPMLElemKU computes 15 split stress components and 12 force blocks against
calcElemKU's 6 and 3, so w < 0.642 is not a physically available value.

CONCLUSION: "no rank owns zero interior elements" and "small work spread" are
NOT simultaneously satisfiable by a contiguous partition at 32 ranks for any
honest PML weight. The two goals are in conflict, and which one matters is an
empirical question: a zero-Ei rank is only a defect if it is UNDER-LOADED.
That is what M2 measures.

## M2 HYPOTHESIS to test
PML_WEIGHT=3.0 (MPI4NodalQuant.py:85) is a guess, documented as a guess. If
the true cost ratio b/a (PML/interior, per element, in the jax kernel) is far
below 3.0, the all-PML ranks 0/1/30/31 are under-loaded and every other rank
waits for them -- which is the 105 ms of 172.6 ms barrier wait. Fit a and b
from the 32 per-rank (Ei, Ep, own-work) triples that one 32-rank run already
prints, then recut with the measured ratio.

## M2 (2026-09-22) the cost ratio, and a null result

Source: docs/perf_snapshots/mpi_scaling_2026-09-22_012512_tpv104_1to32_cacheoff.json
(sha 3004ab9, the pre-existing sweep). own_work = rank_ms - barrier_wait -
exchange - compile/nsteps.

32 ranks, two classes:
  4 all-PML ranks   Ei=0      Ep=12217   own 36.5 ms/step   wait 95-105 ms
  26 mixed ranks    Ei~19227  Ep~5808    own 114.8 ms/step  wait 4-40 ms
The 4 largest barrier waits belong EXACTLY to the 4 zero-Ei ranks. They are
UNDER-loaded, so PML_WEIGHT=3.0 over-weights PML.
  -> b = 2.99e-3 ms/PML elem, a = 5.07e-3 ms/interior elem, b/a = 0.59
  -> no-intercept lstsq over all 32 (Ei,Ep,own): 0.563
  -> same fit at 16 ranks: 1.067 (2-class 1.163);  at 8 ranks: 1.079
The drift with rank count is real: a no-intercept ELEMENT model absorbs each
rank's nodal work into the element coefficients (a rises 1.77e-3 -> 2.46e-3
-> 5.09e-3 from 8 to 32 ranks), so the ratio is stated as a RANGE 0.56-1.16,
"about one, never three".

### NULL RESULT: the cut algorithm was never the problem
A minimax-optimal contiguous split (bisection on the maximum part weight,
greedy feasibility, pad to exactly nranks by halving the heaviest part) was
implemented and compared against the existing quantile cuts on the real
elemType array at 4/8/16/32 ranks and w in {3.0, 1.0, 0.56}. IDENTICAL to
three decimals in every one of the 12 cells. Dropped -- not carried as code.

### Weight scan on the real mesh (predicted elem-work, a/b above)
        nranks=32:  w=3.0  4 zero-Ei, work spread 3.350x
                    w=1.0  2 zero-Ei, work spread 1.601x
                    w=0.56 0 zero-Ei, work spread 1.004x
Threshold arithmetic: Ei>0 on every rank at 32 needs w < 0.642 (M1a). So at
32 ranks the choice is DISCONTINUOUS: below 0.642 rank 0 swallows the whole
31956-element leading PML block (33087 elements, 1.49x the mean element
count, but work-balanced at ratio 0.59); at or above it, rank 0 gets a strict
subset of that block and Ei=0. There is no contiguous partition with both a
small ELEMENT spread and Ei>0 everywhere at 32 ranks. Work balance is the
criterion; Ei>0 is the canary for "w exceeds the true ratio".

## M3 BEFORE numbers, my own hands, sha 14f0b7f (base + cache fix only)
docs/perf_snapshots/mpi_scaling_2026-09-22_035613_BEFORE_w3.0.json
test.tpv104, halo sync, per-step by difference over 40/160 steps, --warmup,
cpus chosen by measured load, EFFECTIVE_CORES 1.00 on every rank.

  ranks    ms/step   speedup   max barrier wait / step   zero-Ei ranks
    1       750.30     1.00x        0.0 ms                  0
    2       479.43     1.56x        2.5                     0
    4       327.10     2.29x       77.2                     0
    8       172.63     4.35x       15.8                     0
   16       171.64     4.37x      100.1                     2
   32       134.12     5.59x       75.8                     4

16 ranks (171.64) is NO BETTER than 8 (172.63) -- and 16 is exactly where the
first zero-Ei ranks appear. Independent corroboration of M2.

XLA compile with the per-rank cache dir + warmup: 0.56-1.62 s at every rank
count, against 4.43-4.93 s in the cache-off snapshot. See M5.
