# Session log — 2026-09-22 — rough-fault normal, TPV29/TPV30

Conductor: wei-lin. Branch `wei/rough-fault-normal-2026-09-22`, worktree
`/home/utig5/dliu/EQdyna.wt-wei-roughnormal`, branched from master `0e83e76`
(v5.14.0 tagged, VERSION 5.15.0).

**STATUS: NOT MERGED. Held at the gate for an owner decision.** The defect is
real and the fix is verified faithful, but the measured outcome is split — it
improves TPV30 by 2.5x and degrades TPV29 by 1.9x on the same independent
oracle, and the owner's stated visual acceptance test is only partly met.

---

## The defect

`case_input/test.tpv29/tpv29GeometryTools.py` decimated the official SCEC
surface correctly (`y[::step]`, rule 17 step 2) and then decimated the
supplied DERIVATIVE columns the same way, at both call sites
(`decimateToDx`, `convertOfficial25m`). A finite difference is a property of a
grid at a spacing, not a value of a surface. At the 500 m gate the normals
described a surface sampled 20x finer than the mesh they were attached to.

Re-derived independently, not inherited. SCEC's own `dF/dx`, `dF/dy` columns
ARE `np.gradient(F, 25)`: median relative difference **0.0215% / 0.0304%**,
corr 1.000000, max|diff| 2.4e-04 over all 1601x801 nodes.

## Where the audit's reproducer and mine first disagreed

My first reproducer computed `tau/strength` as the analytic `S·n` resolution
and got **0.897 at dx=500 with either derivative set** — no defect visible.
The audit's number was 1.16995. The analytic path is self-consistent for ANY
normal, so it cannot see this bug; it is the wrong instrument (rule 4b). The
defect only appears when the traction is taken on the MESH FACET normal and
resolved in the ASSIGNED frame — the C_elastic=0 reconstruction. With that
instrument I reproduce the audit's structure:

| dx | normals | max tau/strength | at or above 1 (3081 interior) |
|---|---|---|---|
| 500 | as-shipped (subsampled) | 1.15461 | 103 |
| 500 | mesh-consistent | **0.89717** | **0** |
| 100 | as-shipped | 0.90474 | 0 |
| 100 | mesh-consistent | 0.89731 | 0 |

Audit reported 1.16995 / 165 and 0.90473; my numbers differ by ~1% and in the
count, most likely `edge_order` or the exact facet-normal definition. Same
verdict, different third digit — recorded rather than reconciled.

Misalignment, the coordinator's metric, median
`|atan(assigned dy/dx) - atan(centered difference of decimated y)|`:

| | dx=100 | dx=200 | dx=500 |
|---|---|---|---|
| shipped file, before | 0.00047 | 0.00183 | 0.01012 rad |
| shipped file, after | 0.00000 | 0.00137 | 0.00967 rad |
| **pipeline output, after** | 0 | 0 | **0 exactly** |

The middle row is the metric applied to the FILE — re-subsampling the new
file's derivatives, an operation the pipeline no longer performs. What the
pipeline produces is exact: `faultGridForCase(500)` returns derivatives
identical to `np.gradient` of its own surface, max|diff| 0.000e+00, asserted
by the new guard.

## What landed on the branch

| commit | what |
|---|---|
| `702b56b` | tool fix (`meshConsistentDerivatives`, both call sites), the `checkDerivativesAgainstOfficial` self-check, regression guard + committed 61x61 SCEC fixture, `.gitignore` |
| `51d37a0` | shipped 100 m / 50 m surfaces regenerated; TPV30's byte-identical copy; README/root-README provenance corrected (rule 11) |
| `bbe6f1a` | `test.reference.results/test.tpv29/` regenerated (rule 7, own commit) |
| `642f119` | `test.reference.results/test.tpv30/` regenerated (rule 7, own commit) |

Surface column bit-identical in both shipped files (max|diff| 0.000e+00);
only the two derivative columns moved (100 m: 3.89e-03 / 3.42e-03; 50 m:
9.58e-04 / 8.94e-04).

Gate run: `python3 testsys/run.py unit regression` SUCCESS, including the new
guard's 8 checks. `run.py all` NOT run — see blockers.

## The outcome, measured against the owner's own 2015 SCEC submissions

Nearest-node match onto `scec_archive/<case>/eqdyna-v3.1-{100m,50m}-2015/cplot`.

| case | | vs 100 m | vs 50 m | mean bias vs 100 m |
|---|---|---|---|---|
| tpv29 | old (subsampled) | 0.0773 (p90 0.2513) | 0.0940 | -0.0082 s |
| tpv29 | new (mesh-consistent) | **0.1443** (p90 0.4073) | 0.1640 | **-0.1593 s** |
| tpv30 | old (subsampled) | 1.0773 (p90 3.2465) | 1.2205 | -1.4596 s |
| tpv30 | new (mesh-consistent) | **0.4358** (p90 1.0851) | 0.5585 | **-0.5165 s** |

median |dt| in seconds. Oracle self-consistency on the same node set (archive
100 m vs archive 50 m): tpv29 **0.0160 s**, tpv30 0.1000 s.

**TPV30 improves 2.5x (3.0x at p90).** That is the case the mechanism is
about: C_elastic=0, the reconstruction path, gain 6.30 into shear.

**TPV29 degrades 1.9x, and that is not noise** — the oracle's own resolution
sensitivity is 0.0160 s, so a 0.0773 -> 0.1443 s move is five times outside
it. TPV29 is C_elastic=1 and never reconstructs traction from nodal forces,
so it is structurally blind to the leak and gets none of the benefit.

Mechanism offered as a hypothesis, NOT settled (rule 4a — no convergence
sweep was run): at 500 m the subsampled normals retained 25 m roughness
variability that acted as an accidental sub-grid model, slowing the front;
mesh-consistent normals describe a genuinely smoother fault, which ruptures
faster. The new mean bias is -0.159 s (ours arrives earlier), consistent with
that direction. The discriminating experiment is a dx=250 run of test.tpv29
both ways — affordable at reduced `par.term`, not run here (item 42: another
`run.py all` was in flight all session).

## The visual acceptance test — NOT met

Regenerated `tpv29_vs_tpv30_cplot_overlay.png` from the fixed runs. TPV30's
isolated closed contours are REDUCED, not gone. Counting the seeds directly
(a ruptured node strictly earlier than every ruptured neighbour):

| case | old | new |
|---|---|---|
| test.tpv30 | 48 | **28** |
| test.tpv29 | 2 | 2 |

42% fewer, still plainly visible in the figure. The misaligned normal was a
real contributor and is not the whole cause. 28 is where the next
investigation starts.

## Corrections to the record

- Pathway item 24(a)'s `test.drv.a6` line — "p10-p90 spread 0.93-1.08,
  fractal roughness perturbing the local fault normal" — is **not** this
  defect. Scope verified by grep: exactly two compsets read this geometry
  path (`test.tpv29`, `test.tpv30`). `test.drv.a6` uses
  `insertFaultType=2`/`seedId=1` (the fractal `generateFaultInterface` path)
  and never reads the SCEC file or its derivative columns; same for
  `liu2020.fdc.rough.250` and `bp1001.fdc.rough.250`. Whether the fractal
  generator has its own version of this problem is **open and unexamined**.
  None of those three references was touched.
- The original brief's five-reference scope was corrected to two by the owner
  mid-session; only `test.tpv29` and `test.tpv30` moved.

## Missing coverage, for iris-vermeulen

Nothing asserts that `un`/`us`/`ud` are consistent with the mesh facets they
bound. `readInputFiles.f90:341-345` checks only `|dy/dx|*dx/dy < 1`, a
mesh-VALIDITY bound, which a wrong-spacing derivative passes comfortably —
that is why this survived. The new guard covers the TPV29/30 geometry TOOL;
it does not cover the solver-side frame, nor the fractal generator's output.
Also `testsys/parity/probe_plastic_traction.py:72` hard-codes
`DEPTH_OFFSET_M = 7.3215`, so it cannot detect the separate +7.3215 m offset
issue by construction.

## Blockers before this can merge

1. **Owner decision on the TPV29 trade** — accept a 1.9x degradation on a
   gated benchmark's independent validation in exchange for removing a real
   defect and a 2.5x TPV30 improvement. A third option exists and is a
   design decision I did not take: use the mesh-consistent normal for the
   FE-facet projection while keeping the fine-scale normal for the prescribed
   initial traction (a Fortran+Python change, not a geometry-tool change).
2. **`python3 testsys/run.py all`** on the branch — not run; another
   `run.py all` (PID 3917167, worktree `wei-item61`) was in flight the whole
   session and item 42 forbids stacking. Expect the other 8 cases bit-identical
   (no other compset reads this geometry), but that is reasoning, not a run.
3. Board rows for zofia, with the evidence above and today's date.

## Housekeeping

- Isolation was broken mid-session: the first three commits were made in the
  MAIN CHECKOUT rather than a worktree, caught by the coordinator. Recovered
  without loss — committed, `git worktree add`, main checkout returned to
  `master` at `0e83e76` clean. Recorded because it is exactly the rule this
  role exists to enforce.
- Run outputs left on disk as evidence (rule 8): `test/test.tpv29`,
  `test/test.tpv30` in the main checkout (gitignored), and the scratchpad
  reproducers.
- One perf snapshot from the single-cell diagnostic run
  (`e2e_cells_2026-09-22_153023_3999634.json`) was moved to the scratchpad
  rather than committed: it records a cell that FAILED by design and does not
  belong in the ledger.
