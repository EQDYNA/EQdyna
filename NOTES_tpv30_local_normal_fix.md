# TPV30 / C_elastic=0 initial-traction fix — Step 1 result: premise falsified (2026-09-17)

Model override note: dispatched on `claude-sonnet-5` per the coordinator's
explicit instruction (Fable-5 hit a session-level 429 on the first attempt,
matching this repo's own repeated precedent for that failure mode).

Worktree: `agent-a848b868d19f5d172`. Resynced at Step 0: `git log -1` ==
`2b8a9a1` (matches the brief's "around 2b8a9a1 or later"), `git status`
clean before and after this investigation — **no fix was written, no file
was changed**, per the outcome below.

This is a follow-on to `NOTES_tpv30_initial_stress.md` (this session's own
prior finding, read first, in full, per the brief). That note is not
contradicted on its headline number (tpv30 33.05-33.11 MPa vs tpv29/spec
~27.8 MPa at `faultst000dp120`, dx=500 both) — it is contradicted on the
**mechanism** it offered for the "Option B" direction, which is what this
note is about.

## Verdict

**Step 1 falsifies Option B's premise. Stopping here, per the mission's own
explicit branch ("if the Step 1 trace falsifies option B's premise, stop
there and report — that is a complete outcome"). No fix was implemented.**

The *literal* Step-1 question ("does the C_elastic=0 element stress reach
the fault node via the split-node/adjacent-element path") is **YES** —
traced and cited below. But the *reason* Option B was proposed — "the
element tensor gets resolved onto the NOMINAL (flat) fault plane instead of
the LOCAL (rough) normal" — does **not** hold on inspection: there is no
nominal-plane resolution anywhere in this pipeline to swap for a local one.
Both sides of the comparison (the FE-reconstructed traction, and the
analytic `fric()`-array traction) already use the identical local rough
normal, sourced from the identical file, applied to an identical stress
tensor. `setPlasticStress` has nothing in it that needs to "resolve onto
the local normal" — it never resolves onto any normal at all; the one and
only place a normal is applied is downstream, in `faulting.f90`, and it is
already the local one, for both terms.

## Step 1 trace, with evidence

**(a) The element stress DOES reach the fault node (literal ask, confirmed).**

`src/fortran/faulting.f90:75-127` (`getNsdSlipSliprateTraction`):
- `massSlave`/`massMaster` come from `fnms` (assembled nodal mass).
- `xyzNodalQuant(k,j,1) = nodalForceArr(...)` (line 81) reads the GLOBAL x,y,z
  internal force at the slave/master node.
- `nodalForceArr` is populated by `assembleGlobalKU.f90:26-39`, which calls
  `calcElemKU` with `stressArr(stressCompIndexArr(nel)+1:...+12)` as the
  element's CURRENT stress state.
- For `C_elastic==0` elements, that `stressArr` state is what
  `meshgen.f90:104` set at mesh-build time: `call setPlasticStress(depth,
  elemCount)` — `meshgen.f90:952-979`.
- So at t≈0, before any dynamics, `nodalForceArr` at a fault node is exactly
  the FE-assembled internal force implied by `setPlasticStress`'s pre-stress
  field. **Confirmed: this is the "dynamic FE term" in
  `nsdTractionVector(1/2/3)` (`faulting.f90:117-127`) at C_elastic=0.**

**(b) Both sides of the comparison already use the LOCAL (rough) normal —
this is the part that falsifies Option B's stated mechanism.**

- `nsdNodalQuant(1,k,j)` etc. (`faulting.f90:89-97`) dot the GLOBAL x,y,z
  quantity with `un(:,iFaultNodePair,iFault)` / `us(...)` / `ud(...)` — this
  is the ONLY place a normal/strike/dip frame is applied to the FE
  reconstruction, and it is applied identically whether `C_elastic` is 0 or 1.
- `un`/`us`/`ud` for `insertFaultType>0` (tpv30 sets `par.insertFaultType=3`,
  drv.a6 sets 2 — both branches take this code path, there is no further
  case on the value) are built at `meshgen.f90:892-905` from **per-node**
  `pfx`, `pfz` — NOT a fixed/nominal value. `pfx`/`pfz` are set fresh on
  every (ix,iy,iz) node by `insertFaultInterface` (`func_lib.f90:80-122`),
  which reads `rough_geo(2,...)`/`rough_geo(3,...)` — i.e. the `dy/dx`,
  `dy/dz` columns of `bFault_Rough_Geometry.txt`
  (`readInputFiles.f90:310-331`, the file this session's brief pointed at).
- `nsdInitTractionVector` (the analytic term, `faulting.f90:71-73`) reads
  `fric(FRIC_SLOT_INIT_NORM/STRIKE_SHEAR/DIP_SHEAR,...)`, which
  `case_input/test.tpv30/user_defined_params.py:273-283` computed as
  `t = S @ nn`, `nn/ss/dd` built from `fDydx[iz,ix]`, `fDydz[iz,ix]` — i.e.
  the SAME `bFault_Rough_Geometry.txt`, written by the SAME tool
  (`tpv29GeometryTools.py`) that supplies the Fortran-side `rough_geo` array
  (checked: `writeBFault`'s "z-fastest" write order,
  `for ix: for iz: write y[iz,ix],dydx[iz,ix],dydz[iz,ix]`, at
  `tpv29GeometryTools.py:161-163`, matches `readInputFiles.f90`'s
  `nnz*(ixx-1)+izz` read index exactly — no off-by-one, no transpose).
- **So `un`/`us`/`ud` (Fortran, FE path) and `nn`/`ss`/`dd` (Python, analytic
  path) are the same local rough-surface vectors, from the same source
  file, for the same node.** There is no "nominal" normal anywhere in this
  pipeline to have been used by mistake.

**(c) The element tensor itself is algebraically IDENTICAL to the analytic
tensor — also not a target for Option B.**

`meshgen.f90:964-976` (`setPlasticStress`, GLOBAL x,y,z frame, only
xx/yy/zz/xy populated, `xz=yz=0`):
```
strVert = -(roumax - rhow*(gamar+1))*depth*grav          ! zz
devStr  = |strVert| * devStrToStrVertRatio * taper(depth)
Sxx = strVert - devStr*cos(2*str1ToFaultAngle)
Syy = strVert + devStr*cos(2*str1ToFaultAngle)
Sxy = devStr*sin(2*str1ToFaultAngle)
```
`case_input/test.tpv30/user_defined_params.py:266-271` (same frame, labeled
x_eq/y_eq/z_eq, same taper `omega`):
```
s22 = -(par.rou - 1000.)*par.grav*dEff                    ! = strVert
sxx = (omega*B11 + (1-omega))*s22
syy = (omega*B33 + (1-omega))*s22
sxy = omega*B13*s22
```
With `R*cos(2a) = B11-1`, `R*sin(2a) = -B13` (the compset's own comment,
`user_defined_params.py:167`), substitution gives `Sxx=strVert*B11`,
`Syy=strVert*(2-B11)=strVert*B33` (since `B11+B33=1.999999≈2`, item G5, a
164 Pa @ 10 km rounding artifact, negligible), `Sxy=strVert*B13` — **matching
the Python `sxx/syy/sxy` term for term at `omega=1`.** Confirmed by hand,
using this compset's own constants (`str1ToFaultAngle=40.375111`,
`devStrToStrVertRatio=0.160739092`, `B11=1.025837`, `B33=0.974162`,
`B13=-0.158649`).

`calcElemKU.f90:48-59` confirms the component ordering used in (b)/(c):
`stress(1)=xx, (2)=yy, (3)=zz, (4)=yz, (5)=xz, (6)=xy` (from the strainrate
build), so `setPlasticStress`'s `component+1/+2/+3/+6` really are
`xx/yy/zz/xy` as assumed above, and `calcElemKU.f90` (grepped in full) has
**no reference to `un`/`us`/`ud`/`str1ToFaultAngle`/any local frame at all**
— the element-level stress update is pure global-Cartesian, exactly as
`setPlasticStress` populates it. There is no second, hidden "nominal
resolution" step at the element level either.

## What this leaves unexplained (and is NOT this note's job to fix)

If both the tensor and the projection are identical and local on both
sides, the measured 12-33% station-to-station shear ratio spread (this
session's own fresh `NOTES_tpv30_initial_stress.md` numbers, not
re-measured in this dispatch — flagged as **inherited**, not fresh, below)
is not "nominal vs local normal." The likely mechanism, offered as a
hypothesis and explicitly **not implemented**:

- Reconstructing a nodal traction from an element-uniform stress tensor via
  `∫ B^T σ dV` on a geometrically DISTORTED (rough) hex element, then
  mass-splitting across a slave/master node pair, is not algebraically
  equal to a pointwise `σ·n` evaluation — it is an FE approximation whose
  error should scale with element distortion (i.e. with local slope), and
  should shrink with `dx`. This was NOT measured with a resolution sweep in
  this dispatch (out of scope for a Step-1 source trace) — the cheapest
  next diagnostic, if this is pursued further, is a dx=500/250/100 sweep at
  the same station, checking whether the FE/analytic ratio converges to 1.0.
- A quadratic-vs-bilinear-form sensitivity argument explains the observed
  asymmetry (normal stress `n^T S n` stable to <1%, shear `s^T S n` volatile
  12-33%) independent of any bug: a stationary quadratic form is
  insensitive to first-order errors in `n`; a bilinear form between two
  orthogonal directions is not.
- **Independent, pre-existing, already-committed corroboration**:
  `testsys/parity/probe_plastic_traction.py` (commit `c7c4f5f`, this
  session, already landed) measured the same FE-reconstruction-vs-intended
  gap for `test.drv.a6` (max|slope| 0.38, much less rough than tpv30's
  0.65) and found `Tn` median ratio 1.0018, `Ts` median ratio 0.9858, p10-p90
  spread 0.93-1.08, attributed explicitly ("not a systematic factor") to
  "the fractal roughness perturbing the local fault normal." **Important
  caveat when comparing the two**: that probe's "intended" reference is the
  NOMINAL flat-fault value (`Syy=strVert`, `Sxy=devStr`, valid only because
  drv.a6's `str1ToFaultAngle=45°` collapses `cos(2a)=0`) — NOT the
  locally-rotated analytic value. tpv30's own comparison in
  `NOTES_tpv30_initial_stress.md`, by contrast, was against the ALREADY
  locally-rotated `fric()`-array value (both sides local) and still showed a
  much wider 12-33% spread — consistent with roughness-scaled discretization
  error (0.65 vs 0.38 slope), but this is a hypothesis, not confirmed by a
  fresh resolution sweep in this dispatch.

## Consequence for Option A (still correctly rejected, now with harder evidence)

Because the FE-reconstructed term at C_elastic=0 already produces something
close in magnitude and sign to the intended local traction (drv.a6's own
probe: ratios near 1.0), Option A ("drop the `*C_elastic` multiplier,"
i.e. additively re-apply `nsdInitTractionVector` on top of the existing FE
term) would come close to **doubling** the on-fault traction for every
already-shipped plastic case, not correcting a zeroed-out value. This
corroborates, with a concrete number, the owner's original reasoning for
rejecting Option A — it is not merely plausible, it is measured (via the
pre-existing probe) to be a real risk.

## Why Step 2/3/4 were not attempted

Per the mission's own explicit instruction: "If the element stress never
actually reaches the fault nodes through that path: the premise for B is
wrong — STOP, do not write a fix, report back... If the Step 1 trace
falsifies option B's premise, stop there and report — that is a complete
outcome per this project's own established pattern this session." The
literal reachability question passed, but the specific mechanistic premise
that would make "make setPlasticStress resolve onto the local normal" a
well-defined, locatable code change did not survive tracing — there is no
nominal-resolution code path in `setPlasticStress`, in `calcElemKU`, or in
the local-frame projection in `faulting.f90` to redirect from nominal to
local; all three already are local. Implementing "Option B" as literally
specified would mean writing an unmotivated change with no identified
defect behind it — exactly what the cardinal rules forbid ("a plausible-
looking scheme that is wrong is worse than an admitted gap").

**Not done, deliberately:** no fix written in Fortran or Python; no test run
against tpv29 post-fix (there is no post-fix); no drv.a6 before/after
comparison run (there is no "after"); no full-sweep re-verification (no
change was made to verify); `test.drv.a6`'s reference untouched;
`test.tpv30`'s registration untouched; no compset-declaration-inertness
flag applies (nothing changed the meaning of `on_fault_vars[...,7/8]`).

## What a reviewer would attack

- The "12-33% is FE-discretization, roughness-scaled" explanation is a
  hypothesis corroborated by one independent data point (drv.a6's probe,
  smaller roughness, smaller — but non-zero — spread) and a
  quadratic/bilinear sensitivity argument, **not** a resolution sweep run in
  this dispatch. It should not be treated as settled.
- I did not re-run tpv29/tpv30 fresh in this dispatch to re-confirm the
  27.79/33.11 MPa numbers myself — those are **inherited** from this
  session's own immediately-prior, already-fresh
  `NOTES_tpv30_initial_stress.md` run, not independently reproduced here.
  What IS fresh in this dispatch is the full source/algebra trace in Step 1
  above (un/us/ud construction, rough_geo indexing match, tensor-term
  algebra, calcElemKU component ordering) — all read directly from the
  worktree at `2b8a9a1`.
- If a resolution sweep were run and did NOT converge toward ratio 1.0, that
  would falsify the discretization-error hypothesis too, and the mechanism
  would still be open. This note does not claim to have found the true root
  cause — only that Option B, as specified, has no valid target.
- The owner may reasonably want a THIRD direction evaluated (e.g. a
  node-wise correction/calibration of the C_elastic=0 static traction toward
  the analytic value, rather than either "always add" (A) or "fix
  setPlasticStress's resolution" (B, no target)) — that is a numerical-
  method design decision (Ingrid Lindqvist's or Rafael Santos's territory
  per this repo's own routing), not something to invent unilaterally here.
