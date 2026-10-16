# TPV30 initial on-fault stress -- root-cause finding (2026-09-17)

Model override note: dispatched on `claude-sonnet-5` after Fable-5 hit a
session-level 429 on the first attempt at this mission (matches this repo's
own prior precedent for that failure mode, per the coordinator's dispatch
message).

Worktree: agent-a1d538c13522ab1b1. Base: master `8d5a5df` (resynced at Step 0,
confirmed `git log -1` matches and `git status` clean before starting).

This is a SEPARATE finding from `NOTES_tpv30_gate.md` (the numpy/jax-vs-Fortran
~30%-at-t=20s divergence). That finding is untouched here. This one is about
whether the C_elastic=0 (plastic) on-fault initial traction is even the
value the spec/compset intends, on ANY backend.

## Step 1 -- controlled experiment (dx=500 both cases, ISOLATES resolution)

Both `case_input/test.tpv29/` and `case_input/test.tpv30/` already default to
`par.dx = 500` (confirmed by reading, not assumed). Station lists are
byte-identical, including `[0.0, -12.0]` (faultst000dp120) at the same index.
Ran BOTH fresh (fresh `./install-eqdyna.sh -m ubuntu` build off master
8d5a5df, `create.newcase` + `case.setup` + `mpirun -np 4 bin/eqdyna`, no
harness) in `scratch/tpv30_initstress/test.tpv29` and `.../test.tpv30`.

`case.setup` for both printed the SAME geometry (nnx=81, nnz=41, dx=500,
corner (-20000,-20000), y in [-1718.3,1590.1], max|slope|=0.6513) confirming
one shared mesh.

**Result** (`faultst000dp120.txt`, first output row, t=0.0417 s):

(Sign note, 2026-09-24, board row 22a: these n-stress values were written
by the pre-22a writer, compression-positive. From 22a on, TPV29/30 station
files are extension-positive per their spec, so the same stress reads
negative -- e.g. -181.504 MPa.)

| case  | h-shear-stress (MPa) | n-stress (MPa) |
|---|---|---|
| tpv29 (C_elastic=1) | 27.7913 | 181.504 |
| tpv30 (C_elastic=0) | 33.1126 | 181.823 |

Same dx=500 mesh, same station, ~19% shear divergence -- **matches the
coordinator's dx=200-vs-100m numbers in direction and rough magnitude, at
IDENTICAL resolution.** This falsifies the "resolution artifact" explanation:
the divergence is present with zero resolution confound. CONFIRMED per
Step 1's own decision rule.

Checked across all 13 on-fault stations tpv29 ships (`compare_stations.py`,
kept as scratch, gitignored): shear ratio (tpv30/tpv29) ranges **0.88 to
1.33** station to station -- NOT a constant multiplicative factor, NOT a
uniform sign. Normal stress agrees to <1% at every station. This asymmetry
(shear volatile, normal stable) is itself a clue used in Step 2.

## Step 2 -- mechanism, located and confirmed with a printed value

`src/fortran/faulting.f90:58-134` (`getNsdSlipSliprateTraction`) and
`:136-184` (`solveSWTW`):

```fortran
nsdTractionVector(1) = (... dynamic FE term ...) + nsdInitTractionVector(1)*C_elastic   ! line 119
nsdTractionVector(2) = (... dynamic FE term ...) + nsdInitTractionVector(2)*C_elastic   ! line 123
nsdTractionVector(3) = (... dynamic FE term ...) + nsdInitTractionVector(3)*C_elastic   ! line 127
...
nodalForceArr(...) = nodalForceArr(...) + xyzTractionVector(j) - xyzInitTractionVector(j)*C_elastic  ! line 179
nodalForceArr(...) = nodalForceArr(...) - xyzTractionVector(j) + xyzInitTractionVector(j)*C_elastic  ! line 180
```

`nsdInitTractionVector` is read directly (line 71-73) from `fric(FRIC_SLOT_INIT_NORM/
INIT_STRIKE_SHEAR/INIT_DIP_SHEAR, ...)` -- the ANALYTIC, locally-rough-normal-
resolved traction that BOTH `test.tpv29/user_defined_params.py` and
`test.tpv30/user_defined_params.py` compute with byte-identical Python code
(same `S = [[sxx,sxy,0],[sxy,syy,0],[0,0,szz]]`, same `t = S @ nn`, same local
`nn,ss,dd` from the rough-fault slopes `fDydx,fDydz`).

**Printed proof** (`python3 -c` reading each case's own
`on_fault_vars_input.nc`, `init_strike_shear` variable, grid index
`ix=40,iz=16` = the `(x=0,z=-12km)` node):

```
tpv29 fric-array init strike shear at (ix=40,iz=16): 27.79126670088692 MPa
tpv30 fric-array init strike shear at (ix=40,iz=16): 27.79126670088692 MPa
max abs diff over whole 41x81 array: 0.0
```

**Both compsets carry the IDENTICAL, spec-correct analytic traction, bit for
bit.** tpv29's actual FE run at t=0.0417s reports 27.7913 MPa -- matching its
own analytic value to 6 sig figs (the `*C_elastic=1` multiplier applies it
directly, dynamic term ~0 since off-fault elements start at zero stress when
`C_elastic==1`: `meshgen.f90:104` only calls `setPlasticStress` `if
(C_elastic==0)`). tpv30's actual FE run reports 33.1126 MPa -- **NOT its own
27.79 MPa analytic value**, because `*C_elastic=0` zeroes that whole term.
What tpv30 actually reports instead is purely the FE-assembled internal-force
reconstruction (`assembleGlobalKU`/`calcElemKU`) of the volumetric
Drucker-Prager pre-stress field set element-by-element by
`meshgen.f90:952-979 setPlasticStress`, which uses a single GLOBAL
`str1ToFaultAngle`/`devStrToStrVertRatio`-based tensor (same formula, same
numbers as the analytic on-fault one) but applied UNIFORMLY per element in
GLOBAL x,y,z, never locally rotated onto each rough-fault node's own normal --
unlike the fric()-array computation, which explicitly re-projects via
`nn,ss,dd` built from the LOCAL slope at that exact node.

**This explains the per-station 0.88-1.33 ratio spread directly**: normal
stress (`nn @ t`, dominated by the `syy`/`szz`-like terms, `nn ~= [0,1,0]` to
first order) is nearly slope-insensitive to first order and matches to <1%.
Shear stress (`ss @ t`, `dd @ t`) is exactly the term that flips sign/mixes
components fastest as the local slope perturbs `ss,dd` away from the
coordinate axes -- so it is the one an FE reconstruction from a
globally-oriented tensor gets wrong by an amount that tracks each node's own
roughness, not a constant scalar bug (ruling out a simple sign/factor/
double-count explanation).

**Root cause, stated plainly**: `C_elastic=0` (every plastic case) silently
discards the analytically correct, spec-conformant on-fault initial traction
that the case's own Python setup computed and shipped in
`on_fault_vars_input.nc`, and substitutes a traction reconstructed purely from
FE divergence of a volumetric stress field that is NOT locally rotated to the
rough fault surface. The published TPV30 (SCEC v3.1, matches tpv29's own
27.79-28.18 MPa band, per `NOTES_tpv30_gate.md`'s earlier numbers) sides with
the DISCARDED analytic value, not with the FE-reconstructed one -- **the
plastic path is the wrong one, as the brief's own hypothesis suspected,
confirmed rather than assumed.**

## Step 3 -- blast radius: test.drv.a6 IS affected, and worse than assumed

`case_input/test.drv.a6/user_defined_params.py`: `C_elastic=0` (line 22),
`insertFaultType=2` (fractal rough fault, `seedId=1`, max|slope| 0.38 per
`NOTES_tpv30_gate.md`'s own prior measurement) -- same code path, smaller
roughness than TPV30's 0.65 but nonzero. `getNsdSlipSliprateTraction` is
friclaw-agnostic (called before branching to `solveSWTW` vs `solveRSF`), so
the `C_elastic` zeroing applies to drv.a6 (friclaw=4, RSF) exactly as it does
to tpv29/30 (friclaw=1, slip-weakening).

drv.a6's own compset **states** an on-fault initial condition explicitly:
`user_defined_params.py:134-135`:
```python
par.on_fault_vars[iz,ix,7] = -120.0e6   # initial normal stress
par.on_fault_vars[iz,ix,8] = 40.0e6     # initial shear stress
```
uniform, not depth-dependent. Read the ALREADY-SHIPPED, gated reference
`test.reference.results/test.drv.a6/faultst000dp075.txt` (frozen from a real
past run, not re-run here): first output row, t=0.0417s:
**h-shear-stress = 24.68 MPa, n-stress = 76.89 MPa.**

**Neither number is anywhere near the compset's own stated 40 MPa / -120
MPa.** This is the SAME mechanism: `C_elastic=0` zeroes the explicit
`-120e6`/`40e6` fric()-array values drv.a6's own file sets, and substitutes
whatever `setPlasticStress`'s `str1ToFaultAngle=45`, `devStrToStrVertRatio=
0.33` DEPTH-DEPENDENT volumetric tensor reconstructs via FE at that fault
depth (7.5 km, not the hypocenter's 12 km) -- which is visibly depth-varying
and nothing like a uniform 40 MPa. This corroborates (rather than
contradicts) the earlier `c7c4f5f` "Tn ratio 1.0018, Ts ratio 0.9858" note in
`test.tpv30/user_defined_params.py`'s own G6 section -- re-reading that note
now: **it was checking FE-reconstruction fidelity AGAINST setPlasticStress's
OWN formula (a numerical-accuracy self-check), not against the compset's
stated -120/40 MPa initial-condition intent.** That earlier check never
established that drv.a6 starts near 40 MPa shear; it only established that
the (wrong, discarding) FE-reconstruction path is internally consistent with
itself to ~1-2%. Both statements are true and are about different things.

**Scope conclusion**: this is NOT tpv30-specific. `C_elastic=0` is shared code
(`faulting.f90`'s `getNsdSlipSliprateTraction`/`solveSWTW`, exercised
identically regardless of friclaw) used by every plastic case shipped today
(`test.drv.a6`, gated, flip-budget reference; `test.tpv30`, unregistered).
**Per the mission brief's own stop-trigger** ("if the divergence touches
shared code used by every C_elastic=0 case... STOP and report"): this
qualifies. `test.drv.a6`'s currently-passing, flip-budget-gated reference was
frozen from a run where the fault's actual initial shear/normal traction is
NOT the value its own compset states -- the reference is *internally
consistent* (a valid regression baseline: "does the code still do what it did
before") but says nothing about whether drv.a6's documented physics intent is
what actually runs.

## What was NOT done (deliberately, per the brief's constraints)

- **No fix attempted.** The two candidate directions -- (a) always apply the
  fric()-array analytic on-fault traction regardless of `C_elastic` (removing
  the `*C_elastic` multiplier at `faulting.f90:119,123,127,179,180`), or (b)
  keep discarding it but make `setPlasticStress`'s volumetric tensor locally
  rotated to the rough-fault normal so the FE reconstruction matches -- are
  both real physics/design decisions with consequences for EVERY currently-
  gated C_elastic=0 case, not a one-line correctness fix. Direction (a) would
  change drv.a6's frozen reference (its flip-budget gate would almost
  certainly need re-measurement) and tpv30's already-frozen candidate
  reference in `test.reference.results/test.tpv30/`. Per the mission's own
  constraint ("if it changes a REFERENCE-backed case's results, STOP and
  report rather than update the reference yourself"), this is exactly that
  case -- not attempted.
- **No reference regenerated, no bound moved, no test registered/changed.**
  `testNameList.py` and `testsys/matrix.py` are untouched (tpv30 was already
  un-registered per `NOTES_tpv30_gate.md`'s prior STOP; left that way).
- **Did not re-run `test.drv.a6` fresh** -- read its already-committed,
  gated `test.reference.results/test.drv.a6/faultst000dp075.txt` instead
  (a real, frozen run's output, not inferred/guessed). If a fresher
  cross-check is wanted, re-running it is cheap (par.term=5s) but was not
  needed to establish the finding.
- **Did not touch the numpy/jax-vs-Fortran divergence** in
  `NOTES_tpv30_gate.md`. That finding's own t=1s bit-exact
  fortran-vs-numpy check is *consistent* with (not contradicted by) this one:
  both backends evidently implement the SAME `C_elastic`-zeroing mechanism
  (they agree with each other on the wrong t=0 value), so that divergence's
  root cause is downstream of t=0, separate from this finding.

## Scratch artifacts (gitignored, kept as evidence, not committed)

- `scratch/tpv30_initstress/test.tpv29/`, `.../test.tpv30/` -- full fresh
  case dirs (create.newcase + case.setup + a real 20s, 4-rank run each),
  built from master 8d5a5df.
- `scratch/tpv30_initstress/compare_stations.py` -- the 13-station shear/
  normal comparison script (diagnostic only, not promoted to testsys/).

## Handoff

This needs an owner decision on which of the two directions above (or a
third) is the intended physics for C_elastic=0, followed by: (1) implementing
it once in Fortran, matched in the Python port (per this repo's own
CLAUDE.md: "change it in BOTH or say plainly which one you changed and why"),
(2) re-measuring `test.drv.a6`'s flip-budget gate from scratch (its current
gate is a valid regression baseline for the WRONG behavior; a fix will
require deliberately re-freezing it, a reviewed act per rule 7, not a side
effect), and (3) re-attempting TPV30 promotion only after both C_elastic
paths are known to produce the SAME analytically-intended on-fault traction --
at which point the SEPARATE numpy/jax-vs-Fortran divergence in
`NOTES_tpv30_gate.md` still needs to be closed before gating.
