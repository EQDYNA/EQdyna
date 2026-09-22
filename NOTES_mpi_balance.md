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
