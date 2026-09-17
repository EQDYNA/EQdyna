# CLAUDE.md

Guidance for Claude Code working in this repository.

NOTE: the `CLAUDE.md` one directory up (`/home/utig5/dliu/CLAUDE.md`) describes
**EQdyna.2Dcycle**, a different project. This file is the one for EQdyna (3D).

`PROJECT_RULES.md` is authoritative — 17 rules, and they are enforced by
`testsys/regression/`, not just written down. Read it before changing a gate,
a bound, or a reference. `pathway_forward.md` is the live status board: open
items there have a `Command` column that must actually run today.

## What this is

A 3D finite-element dynamic-rupture code, in **two implementations of the same
solver**:

```
src/fortran/   26 .f90 + makefile     — the production solver, MPI
src/python/    eqdyna/ — 14 .py       — a serial port, numpy and jax backends
```

Every Python module is named after its Fortran counterpart and is meant to be
read beside it: `faulting.py` ↔ `faulting.f90`, `fric.py` ↔ `fric.f90`,
`driver.py` ↔ `driver.f90`, and so on. Two files have no counterpart by design:
`backend.py` (the numpy/jax array-module adapter) and `__main__.py`.
`assembleGlobalKU.py` covers four Fortran files (`assembleGlobalKU`,
`calcElemKU`, `calcHourglassResist`, `calcElemMass`) because they are one fused
loop in the port.

When you change physics, change it in BOTH or say plainly which one you
changed and why. The port exists so a fix can be verified twice.

## Build and test

```
./install-eqdyna.sh -m ubuntu        # ubuntu/ls6/macos; builds src/fortran, installs bin/eqdyna
export EQDYNAROOT=$(pwd); PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH

python3 testsys/run.py all           # unit + regression + the sweep (24 cells, ~1100 s)
python3 testsys/run.py unit regression       # seconds — run this constantly
python3 testsys/e2e/run_e2e.py --cases test.tpv8 --backends fortran   # one cell
```

CI covers 19 of 24 cells (`testsys/matrix.py`'s `CI_CELLS`, widened
2026-09-16 once the matrix split removed the shared-runner memory ceiling
that used to cap it at 10) split across parallel jobs — `build`,
`unit-regression`, `e2e-ci-fortran-a`/`-b`, `e2e-ci-python-cheap`/`-meng`/
`-tpv29` — not one invocation; see `.github/workflows/test.yml` for the exact
per-job commands. `run.py all` is the wider LOCAL gate; it is not a
reproduction of CI.
Rule 16 was itself wrong about this until 2026-09-16.

## There is ONE test

The e2e sweep, and backend is an AXIS of it, not a tier:

```
for case in testNameList.nameList:        # 8 cases
    for backend in (fortran, python-numpy, python-jax):
        run → canonical frt → compare against ONE committed reference
              at THAT CASE's bound
```

`testsys/matrix.py` is the table (which cells exist, the one bound per case).
`testsys/compare.py` is the only comparison. A cell is SUPPORTED or DECLARED
UNSUPPORTED with a recorded reason — there is no third state and no skip.

References are one `frt.canonical.txt` per case: deduped on rounded (x,y,z)
then lexsorted, so a result is a statement about the PHYSICS and not about the
decomposition. That is what lets a 4-rank Fortran run, a serial Fortran run and
a serial Python run all compare to the same artifact.

Regenerate with `python3 -m testsys.frt_canonical <case_dir>` and commit it as
its own reviewed change (rule 7). Never regenerate a reference to make a cell
go green.

## Adding or reviving a TPV benchmark

**PROJECT_RULES.md rule 17 is the full recipe.** The short version, and the
reason each step exists:

1. **Fetch the official spec first** into `scratch/specs/`, cite it by name and
   part. If the spec states no number, do not invent one — carry an EXCLUDED
   entry in `full_specs.py` saying so.
2. **Decimate geometry, never interpolate.** A fault surface is a fixed
   sampling of one random realisation. `lib.requireFaultGeometryResolution`
   refuses any dx finer than, or not a multiple of, the shipped source.
3. **Check the CODE has that TPV's branch.** This is the step that cost the
   most: `swtwNucleation` branched on `TPV in {201,36,37}` and omitted 29 — the
   case the formula comes FROM — so `test.tpv29` declared `par.tpv = 36` to
   reach its own physics, and that misdeclaration hid a port gap that nucleated
   **0 of 3321** fault nodes. If a case must impersonate another TPV to work,
   that is a code bug, not a config trick.
4. **Gate coarse (minutes at 4 ranks), freeze ONE reference**, and register the
   case in BOTH `testNameList.py` and `matrix.py` (`CASE_BOUND` and `GATE`).
5. **Record the spec-resolution tier in `full_specs.py` without running it** —
   running it is a scheduling decision.
6. **Validate against something independent, as a committed SCRIPT.** TPV29's
   cross-code numbers are prose-only with a gitignored baseline; that is
   pathway item 28 and it is the part of the TPV29 work done wrong.
7. **Run all three backends, and all three must pass.** A new case is 3 cells,
   not 1. **Supporting a TPV means supporting it on every backend** -- a case is
   never added with its Python columns declared UNSUPPORTED to be filled in
   later. If the port lacks a feature the case needs, port the feature first.
   `UNSUPPORTED` records a gap that already exists; it is not a runway for new
   ones.

## Things that will bite you

- **`checkIsOnFault` tests `nodeCoor(2)==0.0d0` by exact float equality.** The
  rough-fault y-blend is deliberately kept OUT of `nodeCoor` (it goes to
  `ycoort` → `meshCoor`). Blending in place — the obvious simplification —
  breaks every rough-fault node. See pathway item 11.
- **`ntotft > 1` is refused in `case.setup`.** Multi-fault input is
  unimplemented: `bStations.txt` line 2 must carry ntotft station counts and
  the writer emits one. Items 7/9/10 sit behind that. See item 17.
- **jax-GPU is nondeterministic run to run** — XLA lowers the duplicate-index
  scatter-add to atomics (measured 8.9e-08 between identical runs). jax-CPU is
  exactly reproducible. No gate may assume bit-identity on GPU.
- **A jitted function that closes over an array pays it as an HLO literal**
  (~7×). Pass arrays as ARGUMENTS. And never build `jax.jit` inside a run
  function — it recompiles every call.
- **`driver.f90:30` DIVIDES** by the mass. The port used `force * inv_mass`
  for friclaw 4/5 and that reciprocal-multiply was a real divergence worth 38
  rupture-arrival flips on drv.a6.
- **`writeCompTime` is initialised to 0 and read from nowhere**, so the
  per-stage timing output is unreachable without editing source. A 511× error
  in `compTimeInSeconds(2)` survived in it undetected.
- **The perf tier gates on PER-STEP cost**, fixed compile/setup subtracted by
  difference over two step counts. A total-wall-clock number mixes in XLA
  compile (15% of a 114-step jax run, 2% of a numpy one) and drifts without the
  solver changing.

## Measure, do not infer

Four claims in this repo turned out not to hold when finally measured: the
"exactly HALF" traction blocker (item 24a — actual ratio 1.0018), a 33% JAX
regression (an artifact of compile time in the metric), a friclaw-4 coverage
gap (friclaw 4 is drv.a6 and tpv104, the two largest cases), and a
compile-fraction figure taken in the wrong environment. Two nearly became
fixes. Prefer a measurement over a plausible mechanism, and when you report a
number, say what produced it.
