# Item 32 — drv.a6 mesh-refinement experiment (dx=500 -> dx=250)

Worktree: /home/utig5/dliu/EQdyna/.claude/worktrees/agent-a062e5e0475453665
Branch: item32-dx250-refinement off master 23971a8
Scratchpad: /tmp/claude-16759/-home-utig5-dliu-EQdyna/e618789a-38a5-4235-9ca3-c8a3a0fe9ba6/scratchpad

## Question
Is drv.a6's ~400-node marginal rupture-arrival population a discretisation
artifact (FRACTION shrinks at dx=250) or genuine bistability at this roughness
and friction (fraction persists)?

## Hypotheses and kill criteria (stated before any run)
- dx=500 baseline: fault nodes 101x51 = 5151, flips (fortran-serial vs
  committed 4-rank ref) expected 372 (7.22%). If fresh reproduction is not 372,
  STOP and report — baseline problem, per brief.
- dx=250: fault nodes 201x101 = 20301 (x3.94).
  - Constant COUNT (~370, fraction ~1.8%) => discretisation artifact.
  - Constant FRACTION (~7.2%, count ~1470) => genuine bistability.
  - In-between (fraction ~3-6%) => discriminates neither; say so.

## Design decision: SAME realisation, not a re-seeded one (checked scoping)
scripts/generateFaultInterface synthesizes the fractal surface on an
lx = par.nfx grid from np.random.seed(seedId=1, hardcoded at call site);
regenerating at nfx=201 draws a completely different random field — that would
confound dx with realisation. Instead: synthesize the dx=500 101x101 scaled
periodic field exactly as the shipped generator does, then Fourier
(zero-pad, 101->202, odd length so no Nyquist ambiguity) interpolate it to the
250 m grid. The synthesis is a fixed trigonometric polynomial, so this
evaluates the SAME surface at the new nodes exactly — no invented roughness
(consistent with the spirit of rule 17 step 2: nothing is interpolated in the
sampled-data sense; the underlying band-limited function is known in closed
form). Coincident nodes must match the dx=500 file to write precision
(fmt='%f', 1e-6 m); this is checked and printed before any run.
The deterministic dip term (nz-1-j)*dx*cos(dip), dip=90, is identical at
coincident nodes by construction ((100-2j)*250 = (50-j)*500).

Perturbation pair probing the marginal population: fortran-serial vs
fortran-4-rank (summation order only), same as the dx=500 A-pair. Parent's
reasoning accepted: the population is defined by roundoff sensitivity, and the
measured dx=500 pairs (372-450 across five very different perturbations) show
the count measures population size, not perturbation identity, so one cheap
pair suffices. Residual risk noted: if the dx=250 pair count sits near a
decision boundary, which perturbation was used could matter at the +-20% level
seen at dx=500 (372 vs 450); conclusions will only be drawn at the 2x level.

## Flip definition (identical at both dx, from testsys/compare.py +
matrix.DRV_A6, hands off): canonical frt col 3 = fnft; ruptured = fnft < 1e4;
flip = ruptured-in-one-only OR (ruptured both AND |dfnft| > 1.0 s);
median|dfnft| over all both-ruptured nodes. matrix.DRV_A6 and the committed
reference are NOT touched.

## Cost estimate (pre-launch, to be refined against fresh dx=500 timing)
Elements: dx=500 160x122x80 = 1.56M; dx=250 320x244x160 = 12.5M (x8).
Steps: dt = 0.5*dx/6000 -> 0.04167 s (120 steps) vs 0.02083 s (240 steps): x2.
Work x16. Sweep history puts drv.a6 fortran 4-rank at O(1-2 min) => dx=250
4-rank ~20-40 min, serial ~1.5-2.5 h. Under the 6 h ceiling. Memory: box has
1007 GB total, 591 GB available; fortran serial dx=500 RSS will be measured,
x8 expected at dx=250 (estimate_HPC_resource heuristic ~1.54 KB/cell ->
~19 GB serial at dx=250). No constraint.

## Constraint log
- 2026-09-21 20:51 CDT: parent sweep PID 779367 still running; no solver runs
  launched yet. Setup only.

## Run log
- 2026-09-21 20:55 CDT: built bin/eqdyna via `./install-eqdyna.sh -m ubuntu`
  (worktree, branch item32-dx250-refinement, base 23971a8). Build log:
  scratchpad/item32_build.log.
- 2026-09-21 21:00 CDT: three scratch cases prepped via
  scratch/item32/prep_cases.sh (create.newcase test.drv.a6 + edits +
  case.setup, all inside the worktree):
  - dx500_serial: stock geometry (seedId=1), header 101x51 dx=500, 5153 rows,
    forced par.nx=ny=nz=1.
  - dx250_serial, dx250_np4: par.dx 500->250; case-local generateFaultInterface
    replaced by scratch/item32/generateFaultInterface.dx250 (Fourier-refined
    SAME realisation). Header 201x101 dx=250, 20303 rows. Both files
    md5-identical (78bb154286b0447b37f607850ef906ae). Validator passed
    (ensureFaultRoughGeometry: "already matches this case; left untouched"
    on idempotent re-check).
  - Geometry checks, measured: Fourier refinement coincident-node max diff
    7.958e-13 m pre-write; WRITTEN files agree at all 5151 coincident nodes
    to 0.0 m (fmt=%f). Peak roughness 532.61 m (coarse) vs 538.37 m (fine
    grid sampling same surface between coarse nodes) -- same amplitude scale.
- 2026-09-21 21:02 CDT: background watcher launched
  (scratch/item32/run_dx500_serial.sh): waits for owner sweep PID 779367,
  then `/usr/bin/time -v mpirun -np 1 bin/eqdyna` in dx500_serial.
  No solver run started while the sweep is alive.
- 2026-09-21 21:42 CDT: owner sweep PID 779367 ended; watcher started
  dx500_serial automatically.
- 2026-09-21 21:47 CDT: dx500_serial DONE. `mpirun -np 1 bin/eqdyna`,
  wall 5:12.30, peak RSS 2.29 GB (time -v), frt.txt0 2,050,098 bytes.
- 2026-09-21 21:48 CDT: dx=500 REPRODUCTION (flips.py, committed flip
  definition): committed 4-rank ref vs fresh fortran-serial =
  **372 flips / 5151 nodes = 7.222%** (only-ref 153, only-run 176,
  timing-shifts 43, median|dfnft| 0.0417 s, phys_max 2.122e7).
  EXACTLY matches the recorded 372 (item 18/21 A=372) and the recorded
  2.12e7 max diff. Baseline is sound.
  Note for interpretation: only ~1250 of 5151 nodes rupture at all
  (ruptured-both 1075), so 372 flips is ~30% of the RUPTURED extent.
  dx=250 fractions will be reported both per-total-node and per-ruptured.
- 2026-09-21 21:50 CDT: refined dx=250 cost estimate from measured dx=500
  serial (312 s, 2.29 GB): x8 elements, x2 steps -> serial ~85 min,
  ~19 GB; 4-rank ~25-35 min. Under the 6 h ceiling; box has 591 GB free.
  Launched scratch/item32/run_dx250.sh in background: dx250_np4
  (mpirun -np 4) then dx250_serial (mpirun -np 1), sequential.
- 2026-09-21 22:06 CDT: dx250_np4 DONE. wall 18:24.13, peak RSS 4.07 GB
  (whole mpirun tree, time -v). Writes frt.txt1 + frt.txt3 only (the two
  y-side ranks own the fault; same as 4-rank behaviour generally -- the
  canonical artifact is decomposition-independent). CONTENDED measurement:
  a 32-rank jax-MPI scaling job (another session's, tpv104) is running on
  this box -- wall times are upper bounds.
  Canonical check: 20301 rows exactly as predicted (201x101); 5533 ruptured
  (27.25%) vs ~24.3% ruptured extent at dx=500 -- rupture character
  comparable across resolutions.
- 2026-09-21 22:06-22:47 CDT: dx250_serial first attempt KILLED at ~40 min
  together with the harness background task (external stop, not a solver
  error; logs preserved as run.log.killed1/time.log.killed1).
- 2026-09-21 22:48 CDT: dx250_serial relaunched DETACHED
  (setsid nohup /usr/bin/time -v mpirun -np 1 bin/eqdyna, pid 1162639,
  scratch/item32/run_dx250_serial_detached.sh). Expect ~85+ min (contended).

---

## RESULT — recovered and measured 2026-09-22 by wei-lin (conductor #4)

The agent that designed and launched this experiment was rate-limited after the
last entry above. **Both dx=250 solver runs had in fact COMPLETED** and sat
unanalysed in its worktree (`scratch/item32/`), while the P1 item-32 board row
still read "dx=250 run in flight 2026-09-21, result pending". `time.log` for
`dx250_serial` ends `Exit status: 0`, written 2026-09-22 00:04; no process of
this experiment was still alive. Recovered by running the experiment's OWN
`flips.py` (which imports the committed `testsys.compare.flip_decomposition`
and `matrix.DRV_A6` verbatim — no new flip definition was introduced).

**The instrument was verified on a known answer first**, before the new number
was quoted:

```
dx500 serial vs committed 4-rank ref (baseline re-derivation)
  nodes 5151  ruptured-both 1075  only-A 176  only-B 153  timing-shifts(>1.0s) 43
  TOTAL FLIPS 372   FRACTION 7.222%   median|dfnft| 0.0417 s   phys_max 2.122e+07
```

Exactly the recorded item 18/21 `A=372` and the recorded 2.12e7 — to every
digit printed.

**The experiment:**

```
dx250 4-rank vs dx250 serial
  nodes 20301  ruptured-both 4554  only-A 979  only-B 966  timing-shifts(>1.0s) 397
  TOTAL FLIPS 2342  FRACTION 11.536%  median|dfnft| 0.0208 s  phys_max 4.143e+07
```

### Against the pre-registered kill criteria (stated above, before any run)

| hypothesis | predicted at dx=250 | measured |
|---|---|---|
| discretisation artifact | constant COUNT ~370, fraction ~1.8% | count **2342** (6.30x), fraction **11.536%** |
| genuine bistability | constant FRACTION ~7.2%, count ~1470 | fraction **grew 1.60x** |

**The discretisation-artifact hypothesis is refuted, and not marginally** — the
flip count grew 6.30x against a node count that grew 3.94x, i.e. the marginal
population grew *super-proportionally* to resolution. The bistability branch
predicted the right direction but understated the magnitude; refinement made
the marginal population a LARGER fraction of the fault, not a constant one. Per
ruptured-both nodes the same way: 372/1075 = 34.6% at dx=500 vs 2342/4554 =
51.4% at dx=250.

`median|dfnft|` is **0.0208 s at dx=250 against 0.0417 s at dx=500 — exactly
half, which is exactly one time step at half the cell size.** The
roundoff-decided-arrest signature item 32 identified is unchanged in character;
only its extent grew.

### The caveat, stated rather than buried

The two comparisons are not perfectly matched: the dx=500 pair is a fresh
serial run against a committed 4-rank reference frozen at an earlier commit,
while the dx=250 pair is a serial and a 4-rank run from the SAME binary in the
same session. That asymmetry can only INFLATE the dx=500 number (it carries
extra sources of difference the dx=250 pair does not), so it cannot manufacture
the observed growth — the direction of the finding is safe against it. A
same-session dx=500 4-rank run would tighten the ratio and has not been done.

Also contended: `dx250_np4` was measured while another session's 32-rank
jax-MPI job was on the box (noted above). That affects its WALL TIME only; flip
counts are not a timing metric.

### Status: not closed here

This is a measurement, not a ruling. Whether it closes the P1 item-32 numerics
question, and whether a bistable-at-this-roughness fault should continue to be
gated by a flip BUDGET at all, is the owner's call. The raw artifacts
(`frt.txt0` serial 8.1 MB, `frt.txt1`/`frt.txt3` 4-rank, and both `time.log`s)
are NOT committed; they remain in worktree `agent-a062e5e0475453665`, which is
deliberately NOT reaped for that reason (rule 8).
