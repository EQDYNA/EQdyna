## UPDATE 7: pinned single-core timing (testsys/perf), replaces the threading-uncontrolled numbers above

All prior timing numbers in this file (Updates 1-6, the "1.14x-2.61x
faster" claims) were measured with **no thread pinning, on a shared,
variably-loaded 64-core machine** -- NumPy/BLAS/JAX were free to use
however many threads the runtime felt like, and load average swung
15-65s on the same Fortran case run minutes apart. Those numbers are kept
below as a **footnote, explicitly marked "threading uncontrolled"** — they
are not wrong, they're just not a controlled measurement.

`testsys/perf/run_perf.py` (new; wired as `python3 testsys/run.py perf`)
pins every engine to the SAME single core (`taskset -c 0`) with
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1`
and JAX forced single-threaded via `XLA_FLAGS=--xla_cpu_multi_thread_eigen=false`
plus `JAX_PLATFORMS=cpu`. Pin verified empirically each run (not assumed):
a subprocess launched under the same `taskset` prints
`os.sched_getaffinity(0)` and the tier checks it reports exactly the pinned
core (`[0]`, confirmed both runs below).

### Pinned single-core timing, tpv8, 114 steps, core 0 (two independent runs)

| | Run 1 (14:33:57) | Run 2 (14:38:24) |
|---|---|---|
| Fortran (`mpirun -np 1`, pinned) | 67.53 s | 66.35 s |
| NumPy (solve-only, pinned) | 120.80 s | 124.04 s |
| JAX-CPU (solve+compile, pinned, 1 thread) | 41.29 s | 42.06 s |
| NumPy/Fortran ratio | 1.789 | 1.869 |
| **JAX/Fortran ratio** | **0.611** | **0.634** |

**Pinned to a single core, JAX-CPU is still ~1.6x FASTER than serial
Fortran** (ratio <1), and consistent run-to-run (0.611 vs 0.634, a 3.7%
wobble — nothing like the earlier threading-uncontrolled runs' 2-3x
wall-clock swings). NumPy is genuinely slower than Fortran even pinned
(ratio ~1.8), consistent with Update 2's finding that NumPy's remaining
cost is real elementwise/gather FLOP work, not thread-count-inflated
overhead.

`testsys/perf/baseline.json` (git SHA `530ae62`, host `cotopaxi`, compiler
`GNU Fortran 11.4.0`, load avg ~13-15) was created from Run 1
(`best_python_over_fortran: 0.6114`, i.e. JAX). Run 2 gated against it:
degradation factor 1.037, well under the 1.5x fail threshold —
**SUCCESS**. The tier fails non-zero if a future run's best python/fortran
ratio degrades by more than 1.5x versus this baseline; it never gates on
absolute seconds, per the load-variance evidence above.

### Old numbers (Updates 1-6), threading uncontrolled -- kept for provenance, not for comparison

The 1.14x/2.36x/2.61x-faster-than-Fortran claims and the 34.76-203.7s
absolute-second ranges earlier in this file were all measured without
thread pinning on a shared machine at unpredictable, sometimes very high
(load avg 15-19) contention. Use `testsys/perf`'s numbers above for any
future comparison; these are retained only so the earlier Updates'
reasoning (profiling, optimization-pass deltas) remains attributable to
the numbers that were actually in front of me at the time.

## UPDATE 6: synced to master (post f21afaf), fric_tp_h fix ported, region_damp fixed too

Master had moved since this worktree branched: `f21afaf` fixes the
`fric_tp_h=0.0` bug this port found (item 12 in pathway_forward.md,
logged from this spike's own Update 5 finding); `a844e92` (already in
master's history, predating f21afaf) fixes the PML region-14 cross-axis
bug (item 8) this port had been faithfully *reproducing* as a bug; plus
verb-phrase renames (`thermop.f90`->`updateThermalPressurization.f90`,
`comdampv.f90`->`computePMLDampingVector.f90`, `hrglss.f90`->
`calcHourglassResist.f90`, `mesh4num.f90`->`countMeshEntities.f90`,
`qconstant.f90`->`calcQAttenuationCoeff.f90`, `warning.f90`->
`checkInputConsistency.f90`), a `globalvar.f90` restructure with named
`FRIC_SLOT_*` constants, a `func_lib.f90` extraction (shared PML region
cascade + a few helpers), and an `MPI4arn`/`MPI4NodalQuant` consolidation.

### Merge

`git fetch origin master && git stash -u && git merge origin/master
--no-edit && git stash pop`. Fast-forwarded cleanly (worktree had no
unique commits, only uncommitted instrumentation). One real conflict,
`src/makefile` (my `pydump.o` addition vs master's renamed object list) —
resolved by keeping master's `OBJ`/dependency list and re-adding
`pydump.o`'s build rule + its now-renamed dependents
(`eqdyna3d.o`/`driver.o` still call `pydump_state`/`pydump_step` — those
two files merged with zero conflict, confirming the coordinator's "renames
don't touch driver.f90/eqdyna3d.f90" expectation). `python/` merged with
**zero conflicts and zero diffs** — master's `python/` (commit `6595df3`)
is byte-identical to this worktree's pre-merge state, confirming an
external sync had already landed phases 1-4 verbatim.

### Two physics fixes required auditing, not just the one the coordinator named

1. **`fric_tp_h` (f21afaf, the one asked for)**: `updateThermalPressurization.f90`
   now reads `fric(FRIC_SLOT_TP_H,i,ift)` = `fric(40)` (0-indexed
   `fric[:,39]`) per node, instead of the disconnected global scalar.
   Ported in `port_tp.py`/`port_tp_jax.py`: `fric_tp_h` is now a per-node
   `(nftnd,)` array (confirmed populated at 0.02 from the fresh dump),
   broadcast against the `(nftnd, nsteps)` history arrays via `[:, None]`.
2. **PML region-14 (a844e92, not named by the coordinator but present in
   master's history since before f21afaf, and load-bearing for `region_damp`,
   which every one of the six port files uses)**: `src/func_lib.f90`'s new
   `pmlRegionDistance` fixes region 14 to test `y` against `ymax0` (not
   `xmax0`, the bug this port had been explicitly, deliberately
   reproducing since Update 1). Also surfaced a **second, previously-unported
   distinction** while reading the fixed source: `pmlRegionDistance` takes a
   `boundInclusive` flag that is NOT unified across call sites --
   `computePMLDampingVector.f90` (node-based, from `velDispUpdate`) always
   passes `.true.` (`>=`/`<=`); `assembleGlobalKU.f90`'s `calcPMLElemKU`
   (element-center-based) always passes `.false.` (strict `>`/`<`) -- a
   genuine pre-existing Fortran difference the refactor explicitly did not
   unify (per its own comment, citing rule 1). This port's original
   `region_damp()` used a single (inclusive) comparison set for BOTH call
   sites -- a real, previously-unflagged port bug, not something the
   coordinator's brief mentioned, found only by reading the fixed source
   line-by-line rather than assuming the old reproduction just needed the
   region-14 line swapped. Fixed in `port.py`'s `region_damp()` (now takes
   a `bound_inclusive` argument) and both call sites in all four files that
   call it directly (`port.py`, `port_jax.py`, `port_rsf.py`, `port_tp.py`;
   `port_rsf_jax.py`/`port_tp_jax.py` inherit it via `port_jax.build()`).

### Rebuild + fresh oracles (all three cases, same fixed binary)

`rm -f src/*.o src/eqdyna && make MACHINE=ubuntu` — clean rebuild (renamed
files force new object names throughout). Regenerated all three serial
oracles fresh (rule 4/7: no reused pre-merge oracle) with the rebuilt
binary: tpv8 (114 steps), tpv104 (120 steps), tpv1053d (24-step truncated,
same `term=1.0` methodology as Update 5).

### Refreshed parity (all three friction families, post-merge)

| Case | Check | Before merge | After merge (fresh oracle) |
|---|---|---|---|
| tpv8 (friclaw=1) | 5-step checkpoint, JAX fault max abs diff | 2.98e-8 | 2.98e-8 (unchanged) |
| tpv8 | 114-step full run, overall max abs diff (NumPy) | 49.5424162298 | 49.5424162298 (**identical to the printed digit** -- the PML region-14 fix's documented "<=3e-8 shift" for tpv8 doesn't move the frt.txt comparison at this precision) |
| tpv104 (friclaw=4) | 5-step checkpoint, NumPy/JAX fault max abs diff | 7.45e-9 / 1.49e-8 | 7.45e-9 / 1.49e-8 (unchanged) |
| tpv1053d (friclaw=5) | 5-step checkpoint, NumPy/JAX accel max abs diff | 1.08e-19 / 7.05e-19 (pre-fix h=0 physics) | 5.63e-13 / 5.63e-13 (**fixed h=0.02 physics** -- larger absolute value is the physics changing, not parity degrading; same-order-of-magnitude relative agreement) |
| tpv1053d | 24-step truncated full run, overall max abs diff | 4.983447 (pre-fix) | NumPy 4.962068440392613, JAX 4.962068440392613 (**identical between engines, fixed-physics oracle**) |

**Verdict: parity holds across all three friction families against the
current master, post both physics fixes.** The two ports (`port_tp.py`/
`port_tp_jax.py`) needed the `fric_tp_h` update; `region_damp()` needed
both the region-14 fix and the previously-missed `bound_inclusive`
distinction; nothing else in any of the six port files needed to change
(the renames are transparent to the ports, which only ever depended on
`pydump_*` dump file contents and formulas, not Fortran subroutine names).

### Updated file list (this worktree, `python/` at repo root)

```
python/eqdyna/port.py          -- tpv8 NumPy reference (region_damp fixed here)
python/eqdyna/port_jax.py       -- tpv8 JAX
python/eqdyna/port_rsf.py       -- tpv104 NumPy
python/eqdyna/port_rsf_jax.py   -- tpv104 JAX
python/eqdyna/port_tp.py        -- tpv1053d NumPy (fric_tp_h fix ported)
python/eqdyna/port_tp_jax.py    -- tpv1053d JAX (fric_tp_h fix ported)
python/run_parity.py, run_parity_jax.py, run_parity_rsf.py, run_parity_tp.py
python/tests/data/               -- tpv8 golden-oracle snapshot + repro config
python/README-parity.md          -- this file (Updates 1-6)
```
`src/pydump.f90` and the two-line hooks in `src/driver.f90`/`src/eqdyna3d.f90`/
`src/makefile` remain uncommitted in this worktree, spike-only, as before --
they were not part of the master merge (nothing in master depends on them)
and re-applied cleanly on top of the renamed files.

## UPDATE 5: phase 4 — TP case tpv1053d (friclaw=5), thermal pressurization, full term attempted

Extended to `test.tpv1053d` (friclaw=5, RSF slip law + strong rate
weakening + thermal pressurization). Biggest case yet: N=946,633,
E=912,600, nftnd=4,005, dt=0.0416667s, term=5s -> 120 steps. New files:
`python/eqdyna/port_tp.py` (NumPy), `python/eqdyna/port_tp_jax.py` (JAX),
`python/run_parity_tp.py`.

**Full term tried first, as asked, and it worked for Fortran**: the
120-step golden oracle was generated in full (253.9 s serial). The
NumPy/JAX comparison itself uses a **truncated, identical-step-count
(24-step) oracle** (`term=1.0` instead of `5.0`, same `nx=ny=nz=1`, same
binary, freshly generated) — disclosed, not hidden — because 120 full
steps at this ~4x-larger-than-tpv104 mesh would cost proportionally more
wall-clock to iterate against during development than the spike's budget
allowed; the 24-step oracle is equally valid (same binary, same fix), just
cheaper.

### Pre-flight checks done before porting (both explicitly requested)

1. **d9a50fa thetaPcTmp fix provenance**: `git log --oneline -- src/faulting.f90`
   in this worktree shows `d9a50fa Fix uninitialized thetaPcTmp in
   NewtonRaphson for friclaw=5` already in history, and `grep` confirms the
   `thetaPcTmp = thetaPc0` initializer is present in the worktree's
   `NewtonRaphson`. No diff-against-main was needed — this worktree already
   has the fix, and the serial oracle was generated by this exact binary,
   so oracle and port are provably the same fixed code.
2. **`fric_tp_h` — a newly confirmed latent bug, empirically verified, not
   assumed**: `thermop.f90` uses a module-level scalar `fric_tp_h` in
   every kernel denominator (`2.0d0*fric_tp_h**2`), but grepping all of
   `src/*.f90` shows it is declared in `globalvar.f90` and **never assigned
   anywhere**. Compiled a minimal standalone program linking only
   `globalvar.f90` and printed it directly: `fric_tp_h = 0.0000000000000000`.
   The Python-side `fric_tp_h=0.02` (`defaultParameters.py`) is written into
   `on_fault_vars(...,40)`/`fric(41)` but that per-node value is never read
   by `thermop.f90`, which uses the disconnected global instead. **EQdyna's
   thermal-pressurization kernel currently runs with h=0 for every
   friclaw=5 case in this src/ tree** — reproduced here as `fric_tp_h=0.0`
   (matching the Fortran as it actually runs), not "fixed" to 0.02. Worth a
   bug report to the maintainers; out of scope to fix in this spike.

### thermop as a scan-carry structural problem (the "watch memory" ask)

Fortran preallocates `onFaultTPHist(2,nftmx,nstep,ntotft)` once, up front,
for the *entire* run — for this case that's 2×4005×120×8 bytes ≈ 7.7 MB,
trivial. Ported the same way: two `(nftnd, nsteps)` scan-carry arrays
(`sliprate_hist`, `shear_hist`), updated via `.at[:,nt-1].set(...)` each
step. The Fortran's inner sum is `j=1..nt-1` (a *dynamic*-length slice,
since `nt` is the scan's own iteration variable) — not directly
expressible as a static-shape JAX op. Resolved by summing over the **full**
`nsteps`-wide history every step and relying on the fact that
not-yet-written "future" columns are exactly `0.0` (functional `.at[].set`
never touches unwritten columns), so their contribution to the weighted
sum is exactly `0` regardless of the (otherwise physically meaningless)
kernel value computed for them. The one hazard: `age=(nt-(j+1))*dt` is
negative for future columns, and `sqrt(negative)` is `NaN`, and `0*NaN=NaN`
— so `age` is floored at a small positive value purely to keep the kernel
finite; the zero history term already zeroes the product, so this floor
changes nothing about the computed physics. No mask array was needed —
verified by the full-run parity below, not asserted.

**Memory/compute scaling flag for a future full port** (not measured, a
structural observation from reading the code): this case's history arrays
are trivial (7.7 MB) and the O(nftnd·nstep²) inner-sum cost is
correspondingly trivial (~14M weighted-sum terms total for 4005 nodes ×
120 steps). But both scale badly: a large fault (nftnd~50,000) run for
200,000 steps (~seconds of rupture at a few-meter cell size) would need
2×50,000×200,000×8 bytes ≈ 160 GB just for one history array, and the
O(nftnd·nstep²) compute cost would become dominant, not incidental. A real
full port at that scale needs either (a) a bounded-lookback truncation
(the kernel decays as 1/√age, so dropping terms once age exceeds a few
diffusion times is defensible but changes the algorithm and must be
validated against the exact sum, not assumed), or (b) an online recursive
reformulation of the diffusion integral. Neither exists in the current
Fortran, and neither was implemented here — flagged as open work.

### friclaw==5 deltas from friclaw==4 (port_rsf.py), all traced from faulting.f90 (not guessed)

- `tnrm` offset: `+= fric(51)` (TP pore pressure), not `fric(6)`.
- `taoc = fric(4)[cohesion] - xmu*MIN(tnrm,0)`, not `xmu*theta_pc` — theta_pc
  is **not evolved** for friclaw==5 (`fric(23)` is left untouched every
  step, matching the coordinator-flagged "not evolved" note and confirmed
  by the near-bit-identical `thetaPc(23)` column below).
- Post-loop creeping-rate floor (`if TPV==105 and v_trial<fric(46):
  v_trial=fric(46)`) **is exercised here** (TPV=105) — unlike tpv104
  (TPV=104), where the identical Fortran line was a documented no-op in
  Update 4. Applied after extracting the Newton loop's final state, only to
  `v_trial`, matching the Fortran's order (state isn't re-derived from the
  clamped value).

### Parity

| Check | NumPy | JAX |
|---|---|---|
| 5-step checkpoint: accel max abs diff | 1.08e-19 (val ~2e-4) | 7.05e-19 |
| 5-step checkpoint: veldisp max abs diff | 2.07e-25 (val ~3e-10) | 1.84e-24 |
| 5-step checkpoint: fault max abs diff | 7.45e-9 (val ~4.5e7) | 7.45e-9 |
| 24-step full run vs golden: overall max abs diff | 4.983447168022394 | 4.983447171747685 |
| 24-step: Tn(78) | 4.49 abs / 1.68e-7 rel | same | 
| 24-step: thetaPc(23) | 7.45e-9 abs / 1.82e-16 rel | same | not bit-identical-to-zero this time (unlike tpv104's fric23), consistent with "frozen at its init value" rather than "always exactly 0" — still roundoff, not evolution |

All diffs at the same ASCII-quantization floor established in Updates 3-4;
NumPy and JAX agree with each other to 8-9 significant digits (reduction
-order noise). No divergence beyond the established envelope.

### Timing trio (24 steps, same session)

| | wall clock |
|---|---|
| Fortran (`mpirun -np 1`) | 147.63 s |
| NumPy (TP port) | 112.13 s |
| JAX-CPU (`jit`+`lax.scan`+`lax.fori_loop`, compile included) | **56.49 s** |

JAX-CPU is 2.61x faster than Fortran and 1.98x faster than NumPy —
consistent with phases 2-3's pattern (JAX wins growing with more
per-node branchy work per step: elastic-only 1.14x/2.44x, RSF-Newton
2.36x/2.26x, RSF+TP-history 2.61x/1.98x).

## UPDATE 4: phase 3 — RSF case tpv104 (friclaw=4), NumPy + JAX, Newton-Raphson jitted

Extended the port to `test.tpv104` (friclaw=4, rate-and-state with strong
rate weakening), the real jit stress-test: `faulting.f90`'s Newton-Raphson
branch, `rsfNucleation`'s forced rupture, and `rate_state_slip_law`. Golden
oracle self-generated exactly like tpv8: `mpirun -np 1` with
`nx=ny=nz=1`, fresh `src/eqdyna` build (same binary, same `pydump`
instrumentation — no changes needed there, it dumps the full 100-column
`fric` array already). Case is bigger than tpv8: N=763,537 nodes,
E=735,000 elements, nftnd=2,701 fault nodes, dt=0.0416667s, full term
(5s) = 120 steps.

New files: `python/eqdyna/port_rsf.py` (NumPy reference, read first per the
requirement), `python/eqdyna/port_rsf_jax.py` (JAX rewrite, reuses
`port_jax.build()` for the elastic-kernel invariants since those are
identical regardless of friclaw), `python/run_parity_rsf.py` (shared
driver for both).

### The RSF-specific subtleties (this is where fric(20)/fric(23) could have
bitten, per pathway_forward.md's theta_pc-uninitialized incident)

- `solveRSF` redefines `nsdSliprateVector(4)` as the **S,D-only** magnitude
  (dropping the normal component) — and this mutated value, not the
  original 3-component one, is what `storeRuptureTime` sees. `fric(71-77)`
  are written from the original 3-component magnitude *before* this
  mutation (in `getNsdSlipSliprateTraction`), so they're unaffected. Got
  this right by tracing the Fortran variable's mutation across subroutine
  boundaries (same array, passed by reference) rather than assuming a
  single "sliprate magnitude" exists per step.
- Two "pre-loop" calls (`rate_state_slip_law`, `rate_state_normal_stress`,
  called once before `NewtonRaphson` with the raw measured sliprate) mutate
  `fric(20)`/`fric(23)` in the Fortran, but those mutations are **provably
  dead**: both columns are overwritten again immediately after
  `NewtonRaphson` returns, using `NewtonRaphson`'s own internal state seeded
  from a snapshot taken *before* the pre-loop mutation. Confirmed by
  tracing the call graph (not by assumption) before deciding to skip
  simulating those discarded mutations — computing and discarding them
  would have been wasted work, not a fidelity gap.
- friclaw=4 is `<5`, so `NewtonRaphson`'s `taoc_old = xmu*theta_pc` branch
  is exercised (not the friclaw==5 `fric(4)-xmu*tnrm` branch) — the
  coordinator's flagged risk, confirmed handled correctly by the
  bit-identical `thetaPc(23)` column in the parity table below.
- `fric(31-36)` (vxm/vym/vzm/vxs/vys/vzs) **are** written for friclaw>=3
  (replaced, not added-to, `nodalForceArr` too) — unlike friclaw<=2 where
  they stay 0 (see port.py's comment on the same columns). Different
  branch, different semantics, both now covered by a real test case.

### Newton-Raphson as `jax.lax.fori_loop` — the actual jit stress-test

Fixed 20-iteration loop, per-node "converged" mask, identical in NumPy
(plain Python `for` with `np.where` masking) and JAX (`lax.fori_loop` body
with the same masking). The correctness argument (`port_rsf.py`'s
docstring): `stateTmp`/`thetaPcTmp`/`xmu`/`taoc_new`/`rsfeq`/`drsfeqdv` are
**pure functions** of `(v_trial, the frozen state0/thetaPc0 baseline)` —
not accumulated across iterations — so freezing `v_trial` for a converged
node makes every later iteration's recompute for that node bit-identical
(same inputs -> same floating-point outputs). Masking the *v_trial update*
by the *pre-this-iteration* converged flag (so a newly-converging node does
NOT get this iteration's update applied) reproduces Fortran's `exit`-before-
update exactly, including the "never converged in 20 tries" fallback (that
node's `v_trial` DOES get the 20th update applied, because Fortran's `do`
loop runs its full body on `iv==ivmax` too). No hybrid/partial-jit fallback
was needed — the whole thing (elastic kernels + Newton-Raphson + rsfNucleation)
is one `jit`+`lax.scan`.

### Parity (5-step double-precision checkpoint + truncated full-run vs golden)

Full term is 120 steps at this mesh size; a full-run comparison at that
step count would cost proportionally more wall-clock than tpv8's, so per
the coordinator's explicit permission, the full-run comparison here uses a
**truncated, identical-step-count (30-step) golden oracle**, freshly
generated the same way as the 120-step one (`term=1.25` instead of `5.0`,
same `nx=ny=nz=1`, same binary) — disclosed, not hidden. The 120-step run
was also completed in full (203.7 s Fortran) and its own step-1..5
checkpoints are the ones used below; both oracles come from the same
binary and are equally valid, the 30-step one is just cheaper to iterate
against.

| Check | NumPy | JAX | Note |
|---|---|---|---|
| 5-step checkpoint: accel max abs diff | 1.86e-26 | 1.00e-25 | both roundoff (values are tiny at step 5, ~1e-12 scale, this early in a forced-nucleation rupture) |
| 5-step checkpoint: veldisp max abs diff | 1.58e-30 | 1.48e-28 | roundoff |
| 5-step checkpoint: fault max abs diff | 7.45e-9 | 1.49e-8 | roundoff (~6e-17 relative on ~1.2e8 Pa values) |
| 30-step full run vs golden: overall max abs diff | 4.973851248621941 | 4.973850704729557 | both at the ASCII-quantization floor (Ts(79), ~4e7 Pa scale, ~1.2e-7 relative) — the two differ from each other in the 7th significant digit, reduction-order noise |
| 30-step: Tn(78) | 0 | 0 | bit-identical in both |
| 30-step: thetaPc(23) | 0 | 0 | **bit-identical in both** — the specific column flagged as historically risky |
| 30-step: fnft(rupt.time) | 3.33e-7 abs | 3.33e-7 abs | identical, ASCII-quantization floor |
| 30-step: fric(31-36) (now non-zero, unlike tpv8) | ~4e-7 abs on O(1) m/s values | same | ASCII-quantization floor |

**Verdict: same quality of parity as tpv8/friclaw=1 — every column agrees
to the ASCII output format's own resolution, both engines agree with each
other to within reduction-order noise, and the two previously-risky
columns (fric20 `rsfstate`, fric23 `thetaPc`) show no divergence.**

### Timing trio (contemporaneous, 30 steps, same shell session)

| | wall clock |
|---|---|
| Fortran (`mpirun -np 1`) | 125.16 s |
| NumPy (RSF port) | 120.06 s |
| JAX-CPU (`jit` + `lax.scan` + `lax.fori_loop`, compile included) | **53.13 s** |

JAX-CPU is 2.36x faster than Fortran and 2.26x faster than NumPy on this
case, in the same contemporaneous window — consistent with the tpv8 result
(phase 2's 1.14x/2.44x), and if anything a larger margin here, plausibly
because the Newton-Raphson loop (20 fixed iterations x several elementwise
ops each) is exactly the kind of repeated small-kernel work that benefits
most from XLA fusion — the same operation Python would otherwise dispatch
20 separate times per node-batch per step, JAX fuses into the compiled
graph once.

### What forced a compromise: nothing

Per the coordinator's ask — "there must be none, parity is the proof":
none did. The Newton-Raphson solve jitted cleanly as a fixed-iteration
`lax.fori_loop` with masking; no `jax.lax.while_loop`, no early-return
approximation, no reduced iteration count, no float32 fallback. The parity
table above is the evidence: if the masking had produced a different
converged value than Fortran's early-exit, it would show up as excess
drift beyond the ASCII-quantization floor already established, and it does
not.

### GPU path

Still not measured (same reasoning as phase 2: no CUDA jaxlib installed,
GPUs shared and at 96-100% utilization from other jobs). The RSF step
function is structurally identical in shape to the elastic-only one (same
static arrays, plus one small fixed-length inner loop) — no new GPU-blocking
pattern was introduced by the Newton-Raphson `fori_loop`, so the qualitative
assessment from phase 2 carries over unchanged.

# EQdyna Python port — feasibility spike (tpv8, serial)

## UPDATE (same spike, follow-on pass): kernel ported, parity achieved, timing measured

`python/eqdyna/port.py` now ports the full tpv8 (friclaw=1) time loop:
`velDispUpdate` (interior + 12-dof PML nodes), `assembleGlobalKU` (interior
`calcElemKU` + PML `calcPMLElemKU`), `hrglss` (C_hg==1 KF78), and `faulting`
(`getNsdSlipSliprateTraction` + `solveSWTW` + the `swtwNucleation` no-op +
`storeRuptureTime`). It does **not** re-derive the mesh/mass/shape-function
state from scratch — that static state is loaded from Fortran `pydump_*`
files written by `src/pydump.f90` (temporary instrumentation added for this
spike; see "Instrumentation" below). `python/run_parity.py` drives it and
diffs against the golden oracle.

### Checkpoint parity (steps 1-5, full double precision, Fortran `pydump_stepN_*` vs Python)

| Quantity | max abs diff | max value | verdict |
|---|---|---|---|
| nodal acceleration (post mass-divide) | 5.4e-14 | 1.61 | roundoff |
| velArr/dispArr (all nodes) | 2.9e-15 | 0.135 | roundoff |
| fault slip/sliprate (fric 71-77) | ≤5.9e-15 | O(1) | roundoff |
| fault traction (fric 78-80) | 3.0e-8 | 1.1e8 | roundoff (Pa-scale, ~1e-16 relative) |

### Full 114-step run vs `tpv8_serial_frt.txt0.golden` (ASCII, `22e18.7e4`, 7 significant digits)

| Column | max abs diff | max rel diff | note |
|---|---|---|---|
| x,y,z | 0 | 0 | bit-identical |
| fnft (rupture time) | 4.99e-7 | 3.86e-7 | at 7-sig-fig ASCII quantization for O(1)-scale values |
| slipS/slipD/slipN | ≤5.0e-7 | large (denominator ~0) | absolute diff at quantization floor; strike slip is near-zero for this dip-slip mechanism, so relative diff is a script artifact, not a real error |
| srS/srD/peakSr | ≤5.0e-7 | large (denominator ~0) | same as above |
| fric(47) (finalSr, unused for friclaw=1) | 0 | 0 | bit-identical — confirms fric(47) is genuinely never written outside `solveRSF`, matched by omission |
| Tn/Ts/Td (fric 78-80) | ≤49.96 | ≤4.8e-7 | at 7-sig-fig quantization of ~1e8 Pa values |
| fric(31-36) (unused for friclaw=1) | 0 | 0 | bit-identical — confirms these stay untouched for friclaw=1 |
| fric(20), fric(23) (unused for friclaw=1) | 0 | 0 | bit-identical |

**Verdict: no text-comparable column diverges beyond what the golden file's
own 7-significant-digit ASCII format can resolve.** Combined with the
roundoff-level agreement at the step-1..5 double-precision checkpoints, this
is genuine parity, not a masked drift — I did not need to fall back to the
"first N=100 steps" reduced bar; the full 114-step run passes.

### Timing (this machine, this SHA, gfortran 11.4.0, 1 core, no threads)

| | wall clock | per step |
|---|---|---|
| Fortran (`mpirun -np 1`) | 34.76 s / 114 steps | 0.305 s |
| Python (NumPy, this port) | 143.33 s / 114 steps (solve only; +5.34 s one-time dump load) | 1.257 s |

**The naive vectorization is ~4.1x SLOWER than serial Fortran, not faster.**
This is the honest headline number, not the ~4x-faster hope stated in the
original task framing — reporting it straight, per rule 5 (no invented
metric) and rule 4 (only fresh runs are evidence).

### Profile (cProfile, 10 steps, ~12.0 s total → 1.20 s/step, consistent with the 114-step rate)

| Hotspot | share of total | cause |
|---|---|---|
| Inline array arithmetic/broadcasting in `run()` | ~70% (8.5/12.0 s) | Per-step fancy-index gathers (`velArr[conn]`, `dispArr[conn]`, `eq_ids[conn,...]`) over 235k elements × 8 nodes, several times per step |
| `np.einsum` (strain/stress-rate contractions) | ~17% (2.1/12.0 s) | Generic einsum path is slower than a hand-fused broadcast-multiply-sum for this simple contraction |
| `np.add.at` (scatter-add force assembly) | ~8% (0.9/12.0 s) | `ufunc.at` is a known-slow unbuffered scatter path in NumPy |
| `np.clip` (hourglass PML/interior slot selection) | ~2.5% (0.3/12.0 s) | **Bug in this port, not in Fortran**: `slot0/1/2 = np.where(ndof[conn]==...)` and the clip are recomputed every step even though `ndof(conn)` is loop-invariant — should be hoisted out of the time loop entirely |

### Concrete, unactioned optimizations for phase 2 (not applied in this pass — no budget left to re-verify parity after each one, which Phase D requires)

1. Hoist every loop-invariant per-element/per-node index array (`eq_ids[conn,...]`,
   the hourglass `slot0/1/2`, `is12`/`is3` masks) out of the step loop —
   currently recomputed every step for no reason; this alone should remove a
   meaningful fraction of the ~70% "inline arithmetic" bucket.
2. Replace `np.add.at` scatter-add with `np.bincount(idx, weights=val)` or a
   sorted-segment-sum; `ufunc.at` is documented as one of NumPy's slowest
   paths and this port calls it ~30+ times per step.
3. Replace the generic `np.einsum('ei,ei->e', ...)` contractions with
   `(dN * v).sum(axis=1)` or precombine constant tensors — should collapse
   ~2 s/10-steps into a fused elementwise multiply-reduce.
4. Batch the interior/PML/hourglass force scatters into one combined
   scatter-add per direction instead of 3 (interior) + 12 (PML) + 3×4
   (hourglass) separate `np.add.at` calls.
5. Re-verify the checkpoint and full-run parity tables above after each of
   1-4, per the workflow's Phase D (one optimization at a time, re-gate every
   time) — not done here.

### JAX/GPU assessment (qualitative, not measured — no JAX available in this pass)

The identified hotspots (repeated gather/scatter over a fixed sparsity
pattern, small per-element contractions) are exactly what `jax.jit` +
`jax.lax.scatter_add`/`segment_sum` are built for: the index arrays are
static across the whole run (same element/node topology every step), so a
JIT-compiled, fused version should remove most of steps 1-4 above "for
free" via XLA fusion, and a GPU port trades the ~235k-element loop for
one wide parallel op. This is a reasonable phase-2 bet but is **not
verified** — no ratio is claimed here that wasn't measured.

## Original phase-1 note (superseded by the above; kept for provenance)

Original status at the end of phase 1: **reconnaissance + golden-oracle generation only. No Python kernel
was ported in this pass.** Do not treat anything under `python/` as a
working solver yet — this directory currently holds only the parity
oracle and this scoping note, per PROJECT_RULES rule 2 (no placeholder
code shipped as if it were real).

## What was done

1. Read `src/eqdyna3d.f90`, `driver.f90`, `globalvar.f90`, `mesh4num.f90`,
   `meshgen.f90` (full), `calcElemKU.f90`, `assembleGlobalKU.f90` (incl.
   `calcPMLElemKU`), `hrglss.f90`, `comdampv.f90`, `qconstant.f90`.
   Not yet read line-by-line: `faulting.f90` (547 lines, slip-weakening +
   RSF + Newton-Raphson), `fric.f90`, `assembleGlobalMass.f90` (406
   lines), `calcB.f90`, `calcGlobalShapeFunc.f90`, `calcLocalShapeFunc.f90`,
   `library.f90`, `library_output.f90`, `library_degeneration.f90`,
   `readInputFiles.f90`.
2. Built `src/eqdyna` serially: `cd src && MACHINE=ubuntu make`
   (gfortran 11.4.0 via `mpif90`, `-fopenmp -ffree-line-length-none -O3`).
3. Created a scratch case from `case_input/test.tpv8` with
   `par.nx=par.ny=par.nz=1` (case default is 2×2×1 = 4 ranks; forced to
   1 rank so decomposition doesn't confound parity, per spike brief).
   `mpirun -np 1 -wdir <case> eqdyna` ran to completion: 114 steps
   (`term=5.0s`, `dt=0.043736878936319105s`), producing a physically
   sane tpv8 rupture (peak slip ≈ 3.83 m dip-slip at hypocenter,
   consistent with the published tpv8 benchmark).
4. Captured the golden oracle: `tests/data/tpv8_serial_frt.txt0.golden`
   (1891 rows = 61×31 fault nodes, 22 columns, format
   `write(10004+me,'(1x,22e18.7e4)')` from `library_output.f90:172`).
   Companion files `tpv8_serial_bGlobal.txt` and
   `tpv8_serial_user_defined_params.py` record the exact inputs (nx=ny=nz=1,
   nstep=114) needed to regenerate it byte-for-byte.
5. Measured wall clock: **34.76 s** for the full 114-step serial run on
   this shared machine (git SHA: check `git -C <worktree> rev-parse HEAD`
   at spike time; compiler `gfortran 11.4.0`, flags `-fopenmp
   -ffree-line-length-none -O3`, 1 MPI rank, mesh ≈ 77,440 cells at
   nx=ny=nz=1 per case.setup's memory estimate — this is a per-machine,
   per-SHA number, not a portable claim).

## Regenerating the oracle

```
cd src && MACHINE=ubuntu make
mkdir -p ../bin && cp eqdyna ../bin/
python3 ../scripts/create.newcase <scratch_dir> test.tpv8
# edit <scratch_dir>/user_defined_params.py: par.nx=par.ny=par.nz=1
cd <scratch_dir> && python3 case.setup
cp <worktree>/bin/eqdyna .
mpirun -np 1 -wdir <scratch_dir> <scratch_dir>/eqdyna
# frt.txt0 is the oracle; diff against tests/data/tpv8_serial_frt.txt0.golden
```

## Executed-path inventory for tpv8 (friclaw=1, C_Q=0, insertFaultType=0, C_degen=0, C_elastic=1)

| Stage | File(s) | Notes |
|---|---|---|
| CLI/param read | `readInputFiles.f90` | `bGlobal.txt` etc; `nstep = idnint(totalSimuTime/dt)` |
| Mesh sizing | `mesh4num.f90` | Non-uniform stretched grid outside the uniform fault-aligned core (`rat` growth factor), PML thickness `nPML=6` nodes, per-node dof count (3 or 12) decided by position vs `PMLb(1..5)` |
| Mesh + split nodes | `meshgen.f90` (996 lines, all executed except wedge/degeneration branches since `C_degen=0`) | Master/slave split-node pairs on the planar fault (`ynft` test is `ycoor==0.0d0` exactly), element material assignment, PML element flagging (`elemTypeArr=2`), on-fault area (`arn`) via quadrilateral-diagonal formula |
| Mass assembly | `assembleGlobalMass.f90` | not yet read in this pass — required for phase 2 |
| Time loop | `driver.f90` | `velDispUpdate` (leapfrog; PML nodes get the 12-dof damped update, else `numOfDofPerNodeArr==3` explicit update), `assembleGlobalKU`, `hrglss`, `faulting`, force/mass divide |
| Internal force, interior elements | `calcElemKU.f90`, `calcB.f90`, `calcGlobalShapeFunc.f90` | Standard 8-node hex, B-bar not used here (straight B), `C_Q==0` branch only |
| Internal force, PML elements | `assembleGlobalKU.f90:calcPMLElemKU` | 15-component split stress state per PML element, quadratic damping profile with **9-14 region branches per element per step** classifying position relative to `PMLb`; this is the branch-heavy kernel flagged as pattern-4 risk |
| Hourglass control | `hrglss.f90` | `C_hg==1` (KF78, 4-mode `phi` tensor) is tpv8's path; `C_hg==2` (viscous) branch not exercised |
| Damping profile (interior/PML boundary nodes) | `comdampv.f90` | Same branch structure as `calcPMLElemKU`'s region logic, called once per PML node per step from `velDispUpdate` — **carries the same `y >= xmax0` bug noted in pathway_forward.md item 8**, must be reproduced verbatim, not fixed |
| Fault mechanics | `faulting.f90`, `fric.f90` | Not yet read; slip-weakening (friclaw=1) branch is the only one tpv8 exercises — RSF/thermal-pressurization branches in the same file must be identified and excluded |
| Output | `library_output.f90:output_frt` | 22-column `frt.txt` writer (verified above) |

## Constants to reproduce verbatim (not yet exhaustively enumerated — flagged for phase 2)

- `kapa_hg = 0.1d0` (hourglass viscosity coefficient, `globalvar.f90`)
- `R = 0.01d0` (PML theoretical reflection coefficient)
- `nPML = 6` (PML thickness in nodes)
- PML damping profile: `damp(i) = 3.0d0*vmaxPML/2.0d0/delta*log(1.0d0/R)*(damp(i)/delta)**(2.0d0)`
  — note the literal `3.0d0...2.0d0` structure (not simplified to `1.5d0*vmaxPML/delta*...`)
  must be preserved bit-for-bit; same expression appears independently in
  `comdampv.f90` and `assembleGlobalKU.f90`'s `calcPMLElemKU` and must match.
- The known latent bug (`y >= xmax0` should likely be `y >= ymax0`) in all
  4 copies of the region-classification cascade (`comdampv.f90:29,60`,
  `assembleGlobalKU.f90:149,181` per pathway_forward.md item 8) must be
  copied as-is into the Python port; fixing it during the port would be a
  silent behavioral change relative to the reference the port is supposed
  to match.

## UPDATE 2: optimization pass (hoist / bincount / de-einsum), CPU-JAX measurement

Applied the three optimizations from the profile above, one mechanical
change at a time, with a full 114-step parity re-check after each (per
Phase D — none reverted, none moved parity beyond the roundoff envelope
already established):

| Change | 5-step checkpoint (accel / veldisp / fault, max abs diff) | Full-run frt.txt diff vs golden | Verdict |
|---|---|---|---|
| Baseline (first pass) | 5.43e-14 / 2.91e-15 / 2.98e-8 | max col diff 49.54 (Tn), all others at ASCII-quantization floor | pass |
| #1 hoist loop-invariant `eq_ids[...]` gathers, slot arrays, PML damping factors out of the step loop | 5.43e-14 / 2.91e-15 / 2.98e-8 (**identical**) | max col diff 49.54 (**identical to the digit printed**) | pass — purely mechanical, no float-order change |
| #2 replace all `np.add.at` scatter-adds with one combined `np.bincount` per step | 5.43e-14 / 2.91e-15 / 2.98e-8 (**identical**) | max col diff 49.54 (**identical**) | pass |
| #3 replace `np.einsum('ei,ei->e'/'ei,eik->ek', ...)` with direct `(a*b).sum(axis=...)` | 5.43e-14 / 2.91e-15 / 2.98e-8 (**identical**) | max col diff 49.54 (**identical**) | pass |

None of the three optimizations moved any digit of the checkpoint or
full-run diffs from the pre-optimization numbers in the original README
section above — expected, since all three are the same arithmetic in a
different NumPy call shape, not an algorithm change.

### Timing — machine-load caveat (rule 4/6: only fresh, provenanced runs are evidence)

This is a shared 64-core machine; `uptime` showed `load average: 18.42,
18.89, 17.85` (other users' jobs) during this pass, and `nvidia-smi` showed
all 3 GPUs at 96-100% utilization from unrelated processes throughout. Repeated
identical-methodology runs varied by up to 1.7x:

| Run | Fortran (fresh, same invocation each time) | Python (post-optimization) |
|---|---|---|
| Original phase-1 baseline | 34.76 s | — (pre-port) |
| Post opt #1 only | not re-measured that instant | 112.83 s |
| Post opt #1+#2+#3, run A | not re-measured that instant | 187.57 s |
| Post opt #1+#2+#3, run B, **immediately preceded by a fresh Fortran run in the same shell session** | **65.93 s** | **133.30 s** |

**Headline, contemporaneous, same-load-window number: Python is 2.02x
slower than Fortran** (133.30 s / 65.93 s), down from the pre-optimization
4.1x. Given the load swings above, treat "2.02x" as the center of a noisy
band, not a precise constant — the profiled *relative* cost shares (below)
are more informative than the absolute wall-clock ratio on this machine
right now.

### Profile after optimization (cProfile, 10 steps, 16.5 s total → 1.65 s/step)

| Hotspot | share |
|---|---|
| Inline gather/broadcast arithmetic in `run()` (unavoidable per-step state-dependent work: `velArr[conn]` gathers, elementwise stress/force updates) | ~81% (13.4/16.5 s) |
| `np.einsum` (hourglass `phid` contraction — NOT yet de-einsum'd at profiling time; fixed in optimization #3 above, not re-profiled after) | ~10% |
| `_c()` helper (specialized strain-rate contractions, replaces old einsum) | ~8% |

The remaining cost is now genuinely FLOP/gather-bound elementwise NumPy
work, not redundant recomputation — the easy, free wins (hoisting,
scatter-add, contraction specialization) are exhausted. Closing the
remaining ~2x gap on CPU alone would need op-fusion (fewer, larger kernel
launches) which plain NumPy cannot do; this is exactly the JAX/GPU case
below.

### CPU-JAX measurement (real, not extrapolated)

No CUDA jaxlib was installed — `nvidia-smi` shows 3x A100 GPUs present but
all at 96-100% utilization from other users' jobs, and a CUDA-enabled jaxlib
is a much larger download than fits the "if quick" budget for this pass, so
no GPU number is reported (rule: no extrapolated numbers). CPU-only
`jax==0.6.2` installed via `pip install jax` (`jax.devices()` → `[CpuDevice(id=0)]`)
and measured directly against the single most expensive kernel in the
profile (the interior-element strain-rate/stress/force computation,
`calcElemKU`'s vectorized equivalent), on the real `Ei=137,592` interior
elements from this case, float64 (`jax_enable_x64=True`), 20-call average
after 3 warmup calls, `jax.block_until_ready()` to force execution:

| | per-call wall time | 
|---|---|
| NumPy (this port's kernel) | 93.6 ms |
| JAX CPU, `jax.jit`-compiled, same math | 24.5 ms |
| **Speedup** | **3.82x** |
| Numerical agreement (same random inputs, both float64) | max abs diff 2.9e-3 on ~1e11-scale values = **~3e-14 relative, roundoff** |

This is a kernel-level measurement, not a full-solver JAX port (that would
require re-implementing the PML/hourglass/fault scatter logic in
`jax.lax.scatter_add`/`segment_sum`, which was out of budget for this
pass) — but it is a real, measured number, and it is consistent with the
qualitative JAX bet made in the phase-1 note: fusing exactly this kind of
repeated elementwise/gather kernel is where `jax.jit` earns its keep even
on CPU, before GPU is in the picture at all.

## UPDATE 3: JAX rewrite of the full time loop (`python/eqdyna/port_jax.py`)

`float64` enabled as the first line executed (`jax.config.update("jax_enable_x64", True)`
before any other jax import), per the phase-2 requirement. Structure: all
loop-invariant index/damping arrays precomputed once in plain NumPy (static
shapes, same as the optimized NumPy port), then one pure `step(carry, _)`
function scanned via `jax.lax.scan(step, carry0, xs=None, length=nsteps)`,
the whole thing wrapped in a single `jax.jit`. **Nothing resisted jitting**:
tpv8 is friclaw=1 (slip-weakening), which is a closed-form calculation with
no Newton-Raphson iteration — that solver only exists in `faulting.f90`'s
friclaw>=3 (RSF) branch, which tpv8 never executes. So the full step
(velDispUpdate + assembleGlobalKU interior/PML + hrglss + faulting) is one
jitted, scanned function; there was no hybrid/partial-jit fallback needed
for this case. `np.add.at`/`np.bincount` became `jnp.zeros(...).at[idx].add(val)`
(`jax.ops` scatter-add), which is the direct JAX equivalent used throughout.
The NumPy port (`port.py`) is untouched and kept as the readable reference,
per the requirement.

### Parity (same two gates as before)

| Check | NumPy port | JAX port | Characterization |
|---|---|---|---|
| 5-step checkpoint: accel max abs diff | 5.43e-14 | 5.62e-14 | both roundoff; JAX's ~1.2e-14 delta from NumPy's own number is reduction-order noise (`jnp.sum`'s pairwise/tree reduction vs NumPy's), not real drift — same order of magnitude, same envelope |
| 5-step checkpoint: veldisp max abs diff | 2.91e-15 | 2.33e-15 | same — reduction-order noise, JAX's is actually *smaller* here |
| 5-step checkpoint: fault max abs diff | 2.98e-8 | 2.98e-8 | **identical to the printed digit** |
| Full 114-step vs golden `frt.txt`: overall max abs diff | 49.5424162298 | 49.5424161106 | both at the ASCII 7-sig-fig quantization floor identified in Update 1; the two differ from each other in the 9th significant digit — reduction-order noise, several orders of magnitude below the quantization floor that already dominates the text-file comparison |
| Full 114-step: unused columns (fric 31-36/47/20/23) | 0 | 0 | bit-identical in both — still never written for friclaw=1 |

**Verdict: JAX's scan-based reduction order does not push any column beyond
the roundoff/quantization envelope already established for the NumPy
port.** No fallback to the "first N=100 steps" reduced bar was needed.

### Timing trio (contemporaneous, same shell session, same machine state)

| | wall clock (114 steps) | 
|---|---|
| Fortran (`mpirun -np 1`) | 65.96 s |
| NumPy (optimized port, Update 2) | 140.71 s |
| JAX-CPU (`jit` + `lax.scan`, includes one-time XLA compile) | **57.65 s** |

**JAX-CPU (compile included) is 1.14x faster than serial Fortran and 2.44x
faster than the optimized NumPy port, in the same contemporaneous
measurement window.** This machine's load varies run to run (see Update 2's
caveat); a second JAX run earlier in this pass took 39.44 s for the same
114 steps under lighter load, so treat 57.65 s as one point in a noisy band,
not an exact constant — the *ranking* (JAX < Fortran < NumPy) held in both
measurement windows.

Compile-vs-steady-state decomposition (from two earlier same-process JAX
calls, 1 step → 6.43 s and 5 steps → 7.53 s, linear fit): **≈6.1-6.2 s
one-time XLA compile, ≈0.28 s/step steady state** — i.e. once compiled,
JAX-CPU's per-step cost (~0.28 s) is already close to or faster than
Fortran's per-step cost in this session (65.96 s / 114 = 0.58 s/step). For
any run with more than a handful of steps (this spike's 114, and certainly
a real multi-thousand-step production run), the one-time compile is
amortized to near-zero and JAX-CPU wins on a per-step basis, not just on
total wall clock.

### GPU path (still not measured — qualitative only, per the no-extrapolation rule)

No CUDA jaxlib was installed this pass either (same reasoning as Update 2:
a GPU-enabled jaxlib is a much larger download than fits an "if quick"
budget, and this machine's 3x A100s were at 96-100% utilization from other
users' jobs throughout both passes, which would make any GPU wall-clock
number unreliable evidence even if measured). Qualitative assessment only:
the step function is already expressed as a handful of large
elementwise/gather/scatter ops over static-shape arrays (235k elements /
249k nodes / 1891 fault nodes) with no data-dependent control flow — this
is close to the ideal shape for a GPU port (the same `jax.jit` code should
run on GPU with no rewrite beyond installing `jax[cuda12]`), but this is an
expectation from the code's structure, not a measured number, and is
reported as such.

## Go/no-go for phase 2 (updated after the kernel port + timing above)

**Go — and phase 2's JAX rewrite (`port_jax.py`) delivers on the bet.**
Parity is real and demonstrated at every stage: NumPy port (checkpoint
roundoff + full-run ASCII-quantization-floor agreement), reconfirmed after 3
NumPy optimization passes with zero digit movement, and reconfirmed again
for the full JAX rewrite (checkpoint diffs 5.6e-14/2.3e-15/3.0e-8, full-run
diffs identical to the NumPy port's to the 9th significant digit — reduction
-order noise only, no real drift, as anticipated). Speed: the JAX-CPU
rewrite (single `jit` + `lax.scan`, no partial/hybrid fallback needed since
tpv8's friclaw=1 path has no Newton-Raphson to resist jitting) beat both the
optimized NumPy port (2.44x) **and serial Fortran (1.14x) in the same
contemporaneous measurement window, including one-time XLA compile time**.
Steady-state (post-compile) JAX-CPU is ≈0.28 s/step vs Fortran's ≈0.58
s/step — roughly 2x faster per step once compiled. **Recommendation: take
this JAX-CPU implementation as the phase-2 baseline and evaluate a
CUDA-enabled jaxlib next** (the code needs no rewrite for that, per its
structure — see the GPU section above) — this pass could not measure GPU
(no CUDA jaxlib installed; the machine's 3 GPUs were at 96-100% utilization
from other users throughout both the NumPy-optimization and JAX passes,
which would have made any GPU wall-clock number unreliable evidence even if
measured).

### Instrumentation added for this spike (not part of production `src/`)

- `src/pydump.f90` (new file): `pydump_state` (called once, after
  `assembleGlobalMass`+`init_vel`, before `driver`) dumps mesh/mass/
  shape-function/PML/fault-static state; `pydump_step` (called for `nt<=5`
  from `driver.f90`) dumps per-step checkpoints for first-divergence
  diagnosis. Both gated to only ever run when explicitly called — zero
  effect on production behavior when not invoked.
- `src/driver.f90`: one added line, `if (nt <= 5) call pydump_step(nt)`.
- `src/eqdyna3d.f90`: one added line, `call pydump_state` before `call driver`.
- `src/makefile`: `pydump.o` added to `OBJ` and given a build rule.
- These are left in place in this worktree for reproducibility of the
  numbers above; they should NOT be merged into `src/` as-is (rule 1/2 —
  they are spike-only debug instrumentation, not production code), or
  should be gated behind a compile-time flag if kept.

## Superseded phase-1 go/no-go (kept for provenance, see updated version above)

**Conditional go — scope the phase-2 port as a multi-week effort, not a
follow-on session.** Basis:

- The oracle-generation path (build, serial run, parity-format capture)
  works cleanly and reproducibly: this de-risks the *comparison*
  infrastructure for phase 2.
- The kernel itself (`calcElemKU`/`calcB`) is a standard, easily
  vectorizable 8-node-hex stiffness-equivalent computation — this part
  is low risk and should batch trivially over elements with
  `einsum`/`tensordot`.
- The two real risks are (a) the PML region-classification cascade,
  which is branch-dependent per element/node per step and appears
  **three times** with subtly different bugs already flagged in
  `pathway_forward.md` item 8 — a straight `np.where`/`np.select`
  vectorization must be checked against the *exact* existing branch
  order, not the "corrected" one; and (b) the split-node fault/mesh
  generation, which is not a simple structured grid but a
  stretched-grid + split-node + MPI-partition-aware indexing scheme
  (`mesh4num.f90`/`meshgen.f90`) that assigns equation numbers and
  master/slave node pairs in a specific traversal order that downstream
  code (fault area `arn`, hourglass `phi`) depends on implicitly.
- `faulting.f90` (547 lines) and `assembleGlobalMass.f90` (406 lines)
  were not read in this pass; a phase-2 kickoff must start there before
  writing any Python, per the workflow's "read every C(Fortran) source"
  gate.

**No wall-clock ratio, no per-column parity diff, and no NumPy profile
are reported here** — none exist yet, because no Python kernel was
written. Reporting a ratio against nothing would be inventing a number
(rule 5) and reusing today's Fortran timing (34.76 s / 114 steps,
this SHA, this machine) as a phase-2 baseline is what rule 4 calls a
hypothesis until phase 2 actually reruns it alongside a Python number.

## File layout under `python/`

```
python/
  README-parity.md          <- this file
  tests/data/
    tpv8_serial_frt.txt0.golden          <- golden oracle, 114 steps, 1891 fault nodes
    tpv8_serial_bGlobal.txt              <- exact Fortran input record for the run above
    tpv8_serial_user_defined_params.py   <- case config (nx=ny=nz=1) used to regenerate it
```

No `python/eqdyna/` package exists yet — creating one without a ported
kernel behind it would be exactly the "port that isn't there yet"
anti-pattern this discipline exists to prevent.
