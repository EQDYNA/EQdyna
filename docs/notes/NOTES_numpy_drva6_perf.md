# NOTES_numpy_drva6_perf.md -- test.drv.a6 python-numpy perf: measurement
method that beats box noise, one landed candidate, one profile, two
STOP-and-report candidates

Base SHA: d470482 (origin/master at session start). Branch:
`origin/mira/numpy-drva6-perf-2026-09-23`.

## Case confirmed real, full-scale (not a toy grid)

`case_input/test.drv.a6/user_defined_params.py`: `par.friclaw = 4` (RSF),
`par.C_elastic = 0` (viscoplastic, so `_drucker_prager` runs every step),
`par.output_plastic = 1`, `par.fric_tp_*` keys present at
`on_fault_vars[:,:,16:19]`/`[:,:,40:42]` -- TP inputs are set even though
friclaw 4 (not 5) never reads `updateThermalPressurization` (confirmed by
profile below: 0 calls to that module). Built via
`eqdyna3d.build_solver_state` (the port's own build path, through
`testsys/e2e/run_e2e.make_serial_case`, imported read-only -- testsys/ itself
untouched): **E=1,982,128** interior+PML elements, N=2,037,036 nodes,
NEQ=9,049,842 equations, nftnd=5,151 fault nodes, full committed `par.term`
gives nstep=120. This is the single biggest case of the four named in the
mission brief (bigger than tpv36's E=947,820 profiled previously).

## Method (task step 1): in-process, interleaved, pinned, real state

Whole-run wall-clock is unusable on this box (previous mission's own
conclusion, reconfirmed here -- see "Whole-run wall time is still noise"
below). Built instead: a `microbench_hg.py` script, in-process, single Python
process, one `numactl --cpunodebind=4,5 --membind=4,5` pin for the whole run:

1. `eqdyna3d.build_solver_state(case_dir)` on the REAL drv.a6 case -- the
   port's own build path, not a synthetic mesh.
2. `driver.run(S, nsteps=5, xp=np)` -- 5 REAL solver steps, to get actual
   (not zero-initialized) `dispArr`/`velArr` for the kernel call (rupture
   has not nucleated yet at step 5, but the arrays are real floating-point
   state from real physics, not degenerate zeros a multiply could special-case).
3. Rebuild `inv = KU.build(S)` (a pure function of S, safe to call twice) and
   `scratch = KU.alloc_scratch(np, inv)` fresh, matching what `driver.run`
   itself does before its loop.
4. Two closures, `run_old()` (`calcHourglassResist(..., scratch=None)`,
   the pre-existing allocate-fresh path) and `run_new()`
   (`scratch=scratch`, the buffer-reuse path), each zeroing `force` first
   (`B.scatter_add` mutates in place on numpy, so a stale `force` would let
   reps 2..N add onto reps 1..N-1's result).
5. Interleaved A,B,A,B,... x 12 reps, `time.perf_counter()` per call,
   report median/min/max/spread, and drop the result if `|gain| < spread`.

## Step 1 result: the buffer-reuse candidate from the 38af040 branch

Cherry-picked (files, not the branch -- that branch's other 16 files
diverged too far from current master's testsys/CI restructuring to apply
cleanly, and only `assembleGlobalKU.py`/`driver.py` are in this mission's
scope anyway): `alloc_scratch` grows a `stage_hg=(E,8)` buffer; on numpy's
non-fused branch (numpy is always non-fused; `_fuse_modes` returns
`B.is_jax(xp)`), `-(phi*r)` becomes `np.multiply(...,out=buf);
np.negative(...,out=buf)`, reusing ONE buffer across all 12 (4 modes x 3
dims) blocks/step instead of allocating fresh. jax is untouched
(`scratch` is always `None` there -- `alloc_scratch` returns `None` when
`B.is_jax(xp)`, so `buf_hg` is always `None` and jax takes the `fuse`
branch, which this diff does not touch at all).

**Bit-identity, single call, drv.a6, numpy:**
```
bit-identical OLD vs NEW single call: True
sha256 old 85eb1a06813f34887096f7cf988929cd145652c81cfffaa13f1c6eb75398d327
sha256 new 85eb1a06813f34887096f7cf988929cd145652c81cfffaa13f1c6eb75398d327
```

**Kernel-level timing, 12 reps each, interleaved, `calcHourglassResist` in
isolation (drv.a6, E=1,982,128, numactl-pinned):**

| variant | median | min | max | std | spread (max-min) |
|---|---|---|---|---|---|
| OLD (`scratch=None`) | 2.591851 s | 2.478938 s | 2.631708 s | 0.046537 s | 0.152770 s |
| NEW (`scratch=buf`)  | 2.438566 s | 2.363583 s | 2.489929 s | 0.039774 s | 0.126347 s |

Median gain = 0.153285 s = **5.91% of OLD's median**, and 0.153 s is bigger
than EITHER variant's own max-min spread (0.153, 0.127) -- **outside the
spread, real**, not the previous mission's INCONCLUSIVE verdict on the same
change measured by whole-run wall clock.

Per-step context: the same 5-real-step `driver.run` that produced the warm
state logged 4811.9 ms/step; `calcHourglassResist` alone (isolated call,
OLD) is 2591.9 ms of that -- consistent with the profile below, where it is
the single largest self-time sink. 0.153 s / 4.812 s ~= **3.2% of a whole
step**, in the same ballpark as this file's own precedent (`stage`/`stage_p`
in `assembleGlobalKU`, "measured 3-5% of the step").

**End-to-end frt parity, drv.a6, numpy, 30 real steps (`par.term=1.25`,
scratch case, NOT the committed reference -- this is a same-code-same-input
A/B, not the release parity gate):**

| run | sha256(frt.txt0) | lines |
|---|---|---|
| patched (working tree with only these 2 hunks applied) | `2e0f2b277e42137fcc741e285b287e968e4b4861fb7004e132deabca8c1d8007` | 5151 |
| baseline (temporary `git stash` of the 2 hunks -> exact master `assembleGlobalKU.py`/`driver.py`, same case dir, fresh run, then restored via `stash apply`+`stash drop`) | `2e0f2b277e42137fcc741e285b287e968e4b4861fb7004e132deabca8c1d8007` | 5151 |

**Bit-identical, 5151 nonzero lines both sides** (not an empty-vs-empty
false green -- papercuts.md's own rule).

**python-jax, same drv.a6 short case, `JAX_PLATFORMS=cpu` (never gate
bit-identity on GPU):**

| run | sha256(frt.txt0) |
|---|---|
| patched | `0ba080a750432043724280f12a99d55edbcf67907e10898760cd036c7d352c1c` |
| baseline (stashed back to master, same procedure as above) | `0ba080a750432043724280f12a99d55edbcf67907e10898760cd036c7d352c1c` |

Bit-identical (expected: the new branch is reached only when `scratch is not
None`, which is never true for jax by construction of `alloc_scratch`).

**Whole-run wall time is still noise, confirmed again on drv.a6** (do not
use it as evidence): the SAME two 30-step runs above logged 221.5 s
(patched) vs 138.2 s (baseline) wall clock -- i.e. whole-run timing said the
*slower* variant was faster, the opposite of the kernel-level, interleaved,
median-of-12 result above. This is exactly why step 1's method exists.

**Verdict: LANDED.** Buffer-reuse in `calcHourglassResist` is bit-identical
on both backends and has a real, spread-exceeding, kernel-level gain
(5.91%, ~3.2%/step) on the biggest case in the sweep.

## Step 2: drv.a6 python-numpy profile (30 real steps, cProfile, patched
tree -- the fresh uninstrumented runs above are the parity oracle, never
this profiled run, per the prior mission's own finding that cProfile
perturbs this port's floats)

`E=1,982,128`, `nftnd=5,151`, 331.3 s wall / 30 steps, 14,443,773 calls:

| ncalls | self (tottime) | cum | file:line | note |
|---|---|---|---|---|
| 30 | 58.972 s | 99.850 s | `assembleGlobalKU.py:497 calcHourglassResist` | the step-1 candidate, already landed |
| 30 | 45.433 s | 141.581 s | `assembleGlobalKU.py:262 assembleGlobalKU` | self time, excl. `_pml`/`_c`/`_drucker_prager` below |
| 540 | 40.149 s | 50.974 s | `assembleGlobalKU.py:86 _c` | 18 calls/step (9 interior + 9 PML); intermediate `(dN*v)` product reallocated every call |
| 124 | 36.186 s | 36.186 s | `numpy c_einsum` (via `_contract_modes`, line 83) | numpy-only spelling, already the memory-frugal one per its own docstring |
| 900 | 19.283 s | 19.283 s | `{numpy.ufunc 'at'}` (`backend.scatter_add`'s `np.add.at`) | 30 calls/step |
| 612 | 11.818 s | 11.818 s | `{numpy.ufunc 'reduce'}` | `_c`'s `.sum(axis=1)` + others |
| 30 | 11.260 s | 16.756 s | `driver.py:40 velDispUpdate` | |
| 30 | 10.388 s | 26.986 s | `assembleGlobalKU.py:373 _pml` | |
| 30 | 9.442 s | 10.467 s | `assembleGlobalKU.py:336 _drucker_prager` | **drv.a6-specific** (C_elastic==0); sqrt/maximum/where/exp, each reallocating |
| 1 | 9.511 s | 12.343 s | `assembleGlobalKU.py:95 build` | ONE-TIME, not per-step |
| 1 | 4.010 s | 24.928 s | `meshgen.py:323 build_node_coordinates` | ONE-TIME |
| (rest) | -- | -- | `faulting.py`/`fric.py`/`updateThermalPressurization.py`/`library_output.py` | **combined < 0.5 s of 331 s (0.15%)** -- friclaw 4 RSF + TP-input-but-unused physics is NOT a per-step cost at this E; `nftnd=5151 << E=1,982,128` |

This generalizes the previous mission's tpv36 finding (`calcHourglassResist`
dominant, `_c`/einsum/`add.at` next) to drv.a6's RSF+TP-input+plasticity
case: **the friction-law/TP machinery is not why drv.a6 is slow; E is**, and
the per-step kernel cost is the same handful of O(E) array ops on every
case, drv.a6 included. `_drucker_prager` (C_elastic==0's viscoplastic
return-map) is the one genuinely drv.a6-specific line in the profile, and it
is a distant 8th at 2.9% of a step, not a hidden dominant cost.

## Candidates acted on vs STOP-and-report

**Acted on:** `calcHourglassResist` buffer reuse (above). Landed.

**STOP-and-report, not attempted:**

- **`np.add.at` -> `bincount`/segment-sum** (900 calls, 19.283 s self,
  ~5.8% of a step): unchanged from the tpv36 finding -- this changes
  summation ORDER (`backend.scatter_add`'s own docstring: numpy accumulates
  in index-array order, a bincount-style group-then-sum does not), and this
  port is gated on bit-identity, not magnitude. Per the brief's own
  instruction, a STOP, not an attempt.
- **`_c`'s intermediate-product buffer reuse** (540 calls, 40.149 s self,
  ~12% of a step -- second largest sink after the landed change):
  same TECHNIQUE as the landed change (`np.multiply(dN, v, out=buf)` then
  `buf.sum(axis=1)`, bit-identical to `(dN*v).sum(axis=1)` -- a plain
  multiply is a plain multiply whether written into a fresh or reused
  buffer), but NOT a STOP by the brief's own gate (no summation reorder, no
  new dependency). Not attempted this session: `_c` is called from TWO
  call-sites with DIFFERENT shapes ((Ei,8) from `assembleGlobalKU`, (Ep,8)
  from `_pml`), and safely sharing one buffer across both without aliasing
  a value that is still live (unlike `calcHourglassResist`'s single call-site,
  single shape) needs its own scalar-first check before being trusted at
  bit-identity, which this session's remaining budget did not spend on an
  unverified rewrite. Estimated payoff, from `_c`'s OWN allocation share of
  its 40.149 s (the sum-reduction itself, not the allocation, likely
  dominates the 74 ms/call average): plausibly 5-15% of `_c`'s tottime, i.e.
  roughly **0.5-2.0 s of a 30-step profiled run, ~0.6-1.2% of a step** --
  smaller than the landed change and not free to get right; flagged for a
  future session with its own microbenchmark-first pass, not silently
  skipped.
- **`_drucker_prager` buffer reuse** (30 calls, 9.442 s self, ~2.9% of a
  step, drv.a6-specific): candidate of the same general shape (sqrt/
  maximum/where/exp each reallocate an `(Ei,)` array), lower priority than
  `_c` (smaller absolute sink, more distinct ops to convert to `out=` forms,
  `xp.where` in particular has no direct `out=` equivalent and would need a
  `np.copyto`-based rewrite to stay bit-identical -- more surface for a
  subtle divergence than a single multiply-negate pair). Not attempted.
- **einsum (`_contract_modes`)**: NOT a candidate. Its own docstring already
  states the memory-frugal numpy spelling was chosen and a matmul-based
  alternative was tried and DIVERGED (FMA contraction, caught at step 5 on
  a prior change). No safe alternate spelling identified; correctly
  untouched.

## Reproduction

```
export EQDYNAROOT=$(pwd) PYTHONPATH=$(pwd)/src/python
# drv.a6 scratch case, full committed term (E=1,982,128, nstep=120):
python3 -c "
import sys; sys.path.insert(0,'testsys'); sys.path.insert(0,'testsys/e2e'); sys.path.insert(0,'.')
import run_e2e
run_e2e.make_serial_case('test.drv.a6', '<dir>', run_e2e.base_env(), term='full')"
numactl --cpunodebind=4,5 --membind=4,5 python3 microbench_hg.py <dir> 5 12
# short-term (par.term=1.25, ~30 steps) case + profile:
numactl --cpunodebind=4,5 --membind=4,5 python3 -m cProfile -o out.prof -m eqdyna <short-dir> --backend numpy
```
`microbench_hg.py` itself is scratch (not committed under `src/python`; its
logic is reproduced verbatim in this file's "Method" section above for
anyone who needs to rebuild it).
