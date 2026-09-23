# NOTES_numpy_wedge_perf.md -- python-numpy wedge-path (test.tpv36/37) profile

Base SHA: dc92983daa975ff52204f25e2f5da3a597539dfc (origin/master at start of session)

## Setup
Scratch case: `scripts/create.newcase <dir> test.tpv36`, forced serial
(`par.nx=par.ny=par.nz=1`, required -- standalone python backend refuses
npx*npy*npz>1), `par.term` overridden to 0.2 (19 steps) and 0.4 (37 steps)
for short cProfile/timing runs. `case.setup` reports "Estimated total cells
... 316800" (Fortran cell estimate for the WHOLE domain including ghost/
non-fault padding); the actual mesh element count read back from
`build_solver_state` is **E=947,820** (Ei=685,236 interior, Ep=262,584 PML,
6,720 wedge (elemType 11/12, 0.7% of E) -- confirms the wedge branch itself
is NOT the extra per-step cost: wedge elements are folded into E_int at
BUILD time (assembleGlobalKU.build's `(elemType==1)|(elemType>10)` mask) and
run through the exact same per-step kernel as every other interior element,
with no extra per-step branch).

Both backends profiled with `python3 -m cProfile -o ... -m eqdyna <case> --backend {numpy,jax}`,
pinned `numactl --cpunodebind=4,5 --membind=4,5`, `JAX_PLATFORMS=cpu`.

## Profile verdict: numpy's cost is real Python-level compute; jax's is invisible to cProfile

**numpy, 19 steps, cProfile-instrumented, 149.09 s total:**

| ncalls | self (tottime) | cum | file:line |
|---|---|---|---|
| 19 | 57.609 s | 69.735 s | assembleGlobalKU.py:496 `calcHourglassResist` |
| 342 | 21.018 s | 24.422 s | assembleGlobalKU.py:86 `_c` |
| 19 | 20.893 s | 57.721 s | assembleGlobalKU.py:261 `assembleGlobalKU` |
| 81 | 11.742 s | 11.742 s | numpy `c_einsum` (via `_contract_modes`'s numpy branch, line 83) |
| 589 |  5.826 s |  5.826 s | `{method 'at' of numpy.ufunc}` (via backend.py:56 `scatter_add`'s `np.add.at`) |
| 19 |  4.787 s | 11.219 s | assembleGlobalKU.py:372 `_pml` |
|  1 |  4.393 s |  5.588 s | assembleGlobalKU.py:95 `build` (one-time, not per-step) |

**jax, 19 steps, cProfile-instrumented, 21.01 s total:** the ENTIRE 19-step
solve is one `jax.jit(body)` call (backend.py:347); cProfile sees exactly
one frame, `jax/_src/api.py:3106 block_until_ready`, 6.428 s self/cum --
XLA compiles the whole step loop into one executable and cProfile cannot
see inside it. This is a structural fact about jax, not a profiling gap on
our side: **cProfile cannot be used to attribute jax's per-step cost to a
line**, only to confirm that (a) it is one opaque call and (b) the
non-loop, host-side setup (`build`, `compute_hourglass`,
`compute_element_shape`, `assemble_mass`, `build_elements` -- all numpy,
run once, same on both backends) is a small, comparable fraction of the
total on both columns.

## Verdict on the candidates named in the task

- **Wedge-specific python loops** (candidate: "degenerate-element kernels"):
  REFUTED as a per-step cost. The only genuinely wedge-specific code (a
  per-node closed-form mass formula, a masked local-derivative swap) runs
  ONCE at build time in assembleGlobalMass.py (`_contm_wedge_node_mass`,
  `compute_element_shape`'s `is_wedge` mask), confirmed above at 6,720 of
  947,820 elements. The per-step kernel (assembleGlobalKU.py) has no wedge
  branch at all.
- **Per-element python loops jax vectorises**: none found in the per-step
  path (driver.step -> part_a/part_b -> assembleGlobalKU/calcHourglassResist/
  velDispUpdate/faulting). Every remaining `for e in range(E)` in the repo
  (assembleGlobalMass.py:228,390) is build-time only.
- **np.add.at (scatter-add)**: CONFIRMED slow in absolute terms (5.826 s
  self time / 589 calls, ~9.9 ms/call over the profiled run) but is NOT the
  dominant term (10.4% of numpy's ~55.9 s of self time in named functions).
  NOT changed: bincount/segment-sum alternatives group-then-sum in a
  DIFFERENT order than add.at's index-array order (backend.py:45-57's own
  docstring states this is load-bearing for bit-identity), so per the
  brief's own instruction this candidate is a STOP, not an attempt --
  flagged here, not silently tried and reverted.
- **Repeated allocation per step**: CONFIRMED and the one candidate acted
  on (see below). `calcHourglassResist` builds a fresh (E,8) = 947,820 x 8 x
  8B = 60.7 MB negated-product block 12 times per step (4 modes x 3 dims)
  with no buffer reuse, unlike assembleGlobalKU's own force blocks (which
  already reuse `scratch['stage']`/`scratch['stage_p']`, per this file's own
  `mul_into` docstring, "measured 3-5% of the step" for that kernel).
- **dtype promotion**: not observed; both backends run float64 throughout
  (see backend.array_module's `jax_enable_x64` and to_device's int-only
  narrowing).

## Change made: buffer reuse in calcHourglassResist (numpy only)

`src/python/eqdyna/assembleGlobalKU.py`:
- `alloc_scratch` now also allocates `stage_hg = np.empty((inv['E'], 8))`
  (numpy only; `None` on jax, unchanged).
- `calcHourglassResist` takes a new trailing `scratch=None` parameter
  (backward compatible: every existing caller -- driver.py,
  evidence_wedge_kernel.py, probe_mpi_step_split.py -- passes exactly 6
  positional args and gets the unchanged functional path when `scratch` is
  not supplied).
- On the numpy, non-fused branch only (fuse is always False for numpy --
  `_fuse_modes` returns `B.is_jax(xp)`), `-(phi[:, m, :] * r[d][:, None])`
  becomes `np.multiply(..., out=buf); np.negative(buf, out=buf)`, reusing
  ONE buffer across all 12 blocks/step instead of allocating fresh each
  time. Same two operations (a multiply, then a sign-bit flip), same
  roundings, just not reallocated -- exactly the precedent this file's own
  `mul_into`/`iadd` docstrings already establish and gate bit-identity on.

`src/python/eqdyna/driver.py`: one-line change, the ALREADY-IN-SCOPE
`scratch` (closed over by `make_step_parts`, already threaded to
`assembleGlobalKU`) is now also passed to the `calcHourglassResist` call at
line 140. This is the only touch to driver.py; the sink itself lives in
assembleGlobalKU.py.

## Bit-identity evidence

Isolated synthetic check (random (E,4,8) phi against random r, E=5000):
functional `-(a*b)` vs `multiply(out=)+negate(out=)` -- exact
`np.array_equal` on all 12 blocks, byte for byte.

End-to-end, test.tpv36 scratch case, numpy backend, UNINSTRUMENTED (no
cProfile -- see papercut below):

| run | steps | sha256(frt.txt0) |
|---|---|---|
| baseline (unmodified, run 1) | 19 | e4bb3b4f1f8e391810365bda2e679ed65e15b08867eb2f0a7b7502861a5b6c47 |
| baseline (unmodified, run 2) | 19 | e4bb3b4f... (identical) |
| patched  (this change, run 1)| 19 | e4bb3b4f... (identical) |
| baseline (unmodified) | 37 | 2706da4366c07fe058156d686f66ea7909528d6b474bde8eeffcd5d7cb1fd8f7 |
| patched (this change) | 37 | 2706da43... (identical) |

**Byte-identical at both step counts.** `wc -l` on both frt.txt0 files: 3477
lines (19-step case), nonzero and matching -- not an empty-vs-empty false
green.

### Papercut hit during this check (added to ~/code/papercuts.md)
The FIRST "divergence" observed was spurious: comparing a **cProfile-
instrumented** numpy run's frt.txt0 (sha 3ebf7cc0...) against an
**uninstrumented** run of the SAME unmodified code gave DIFFERENT hashes,
while two uninstrumented runs of unmodified code (and two uninstrumented
runs of the patched code) all agreed with each other. cProfile's
instrumentation itself perturbs this run's floating-point result (most
likely by changing thread-pool/BLAS scheduling inside `np.einsum`'s
threaded path -- not yet root-caused further, and not this change's
concern). **A cProfile-captured frt.txt0 must never be used as a parity
oracle**; only uninstrumented runs were compared above.

## Timing: INCONCLUSIVE on this box right now, and said so rather than reported

Attempted a before/after ms/step measurement, term=0.2 (19 steps) and
term=0.4 (37 steps), pinned `numactl --cpunodebind=4,5 --membind=4,5`,
interleaved baseline/patched runs. Wall times for IDENTICAL code, same
term, back to back:

- unmodified, 19 steps: 50.587 s, 42.429 s, 145.256 s (2.9x spread)
- patched, 19 steps: 41.707 s, 83.205 s (2.0x spread)
- unmodified, 37 steps: 205.261 s, 220.792 s
- patched, 37 steps: 136.469 s

`/proc/loadavg` read 35-48 throughout (consistent with the stated shared-box
load). A 3-5% effect (the size measured for the identical technique already
in assembleGlobalKU) is not extractable from swings this size with this few
samples -- reporting a "before/after %" from these numbers would be
reporting noise, not the change. Per CLAUDE.md ("Measure, do not infer"):
this is flagged INCONCLUSIVE rather than rounded to a number that looks
clean. The change is retained because (a) it is proven bit-identical above,
(b) it is structurally the same, already-verified-safe technique as the
adjacent `mul_into`/`stage`/`stage_p` buffers in the same file, and (c) it
strictly reduces peak allocation churn (one persisted (E,8) buffer instead
of 12 fresh (E,8) temporaries/step, ~728 MB/step of malloc/munmap avoided
at E=947,820) regardless of whether the wall-clock win is currently
measurable through this box's contention.

## Why this does not close the 4.3x gap, and what would

The wedge cells (tpv36/37) are not slow because of wedge-specific code --
they are slow because E is large (947,820, by far the biggest of the three
cases named in the brief) and the numpy path executes the per-step kernel
as ~20 separate Python-level array operations, each a full pass over
memory with its own allocation, while jax's XLA JIT (backend.py:347
`run_time_loop`) compiles the ENTIRE step (all of assembleGlobalKU,
calcHourglassResist, velDispUpdate, faulting) into ONE fused executable
with no per-op allocation and no Python dispatch overhead. That gap scales
with E, which is exactly why the wedge cases (largest E) show the largest
gap and why this is "numpy generally" (tpv29) at smaller scale, per the
brief's own framing.

Closing more of it without an ask-first escalation is not available under
this session's constraints:
- **np.add.at -> bincount/segment-sum**: STOPPED, not attempted, per the
  brief's own instruction -- changes summation order, and this port is
  gated on bit-identity, not magnitude.
- **True operator fusion in numpy** (numexpr, numba, Cython) would need a
  new build-system dependency -- an explicit stop-and-ask trigger in this
  agent's charter ("Optimization requires a C extension. Ask before
  introducing a build-system dependency"), not taken unilaterally here.

## Reproduction
```
scripts/create.newcase <dir> test.tpv36
# append to <dir>/user_defined_params.py: par.term = 0.2 (or 0.4); par.nx=par.ny=par.nz=1
cd <dir> && python3 case.setup
numactl --cpunodebind=4,5 --membind=4,5 python3 -m cProfile -o out.prof -m eqdyna <dir> --backend numpy
numactl --cpunodebind=4,5 --membind=4,5 python3 -m cProfile -o out.prof -m eqdyna <dir> --backend jax
```
