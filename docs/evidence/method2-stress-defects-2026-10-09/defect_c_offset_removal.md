# Defect C: `setPlasticStress` depth-offset removal -- measured impact

Board PR #164 (2026-10-09). Owner decision (relayed, superseding item 24(d)'s
2026-09-16 "do what fortran does ... not to be derived or changed" ruling):
"it should be universal right?" / "we are setting stresses at the center of
each cell" -- the undocumented `+ 7.3215d0` shift in `setPlasticStress`'s
depth argument (`src/fortran/meshgen.f90:201`, mirrored verbatim in
`src/python/eqdyna/meshgen.py`) is removed universally, in both backends, for
every `C_elastic==0` case. This is an intended physics change, not a
tolerance-chase: the pre-fix depth argument put the true free surface
(z=0 m) 7.3215 m below where `setPlasticStress` thought it was, leaving a
nonzero `sigma_zz` at the real free surface -- an unbalanced traction, since
the free surface carries no boundary stress -- which made every Method-2
case sink from step 1.

Mechanism, confirmed by measurement (rule 4a: outcome measured with the cause
present and then removed, not inferred from code reading alone): the uniform
additive offset does not corrupt the lithostatic stress GRADIENT (which
`assembleGlobalKU.f90:21`'s gravity body force correctly balances against),
but it does leave a nonzero traction at the TRUE surface specifically.

## Cases affected

`grep -l "par.C_elastic = 0" case_input/*/user_defined_params.py` -- exactly
5 cases, matching the owner's expectation:

```
test.tpv13       (registered, gated fortran + python-jax, FaultMod benchmark)
test.tpv27       (registered, gated fortran + python-jax, FaultMod benchmark)
test.tpv30       (registered, gated fortran + python-jax, FaultMod benchmark)
test.drv.a6      (registered, gated fortran + python-jax [flip-budget],
                  NO benchmark -- self-comparison only, per owner instruction)
test.drv.a6.v2   (NOT registered in testNameList.py/matrix.py; experimental;
                  NO benchmark -- self-comparison only)
```

All other registered cases (18 of 22 in `testNameList.nameList`) use
`C_elastic==1` (Method 1), never call `setPlasticStress`, and read exactly 0
for the removed term either way -- confirmed unaffected below (two of those
18, `test.tpv22`/`test.tpv23`, are release-only and not run in the everyday
sweep this section reports).

## Step-1 surface residual acceleration (EQDYNA_DUMP_EQUIL=1, nt==1)

**Re-measured 2026-10-09 by a committed script** (`measure_surface_residual.py`,
this directory; raw output in `surface_residual_output.txt`, same directory),
not restated from an earlier ad-hoc run -- the earlier table below this one
labelled a SIGNED mean as `|a_z|`; this correction reports both, separately,
for every case and both before/after states. The script reverts the single
Defect C call-site line in `src/fortran/meshgen.f90`, rebuilds `bin/eqdyna`,
measures, then restores the committed source and rebuilds again, verified
clean via `git status` before it exits -- it never leaves the tree dirty.

Free-surface (z>-1e-6 m) non-fault (faultTag==0), non-PML (dof<=3) nodes,
serial fortran backend, gate dx (500 m):

| case          | mean&#124;a_z&#124; before | mean&#124;a_z&#124; after | signed mean before | signed mean after |
|---------------|------------------:|----------:|------------------:|----------:|
| test.tpv13    | 0.2066            | 0.0206    | 0.2066            | -0.0020   |
| test.tpv27    | 0.1795            | 0.0000    | 0.1795            | -0.0000   |
| test.tpv30    | 0.1795            | 0.0078    | 0.1795            | -0.0000   |
| test.drv.a6   | 0.1581            | 0.0481    | 0.1581            | 0.0481    |

Before the fix every node sinks the same direction, so the absolute and
signed means agree exactly (within print precision). After the fix
tpv13/tpv27/tpv30 have an absolute mean close to an order of magnitude
smaller AND a signed mean that is itself near zero (and slightly negative,
not the old uniform positive sink) -- the residual no longer has a
systematic one-sided direction, it is now symmetric FEM
constant-stress-element discretization noise. drv.a6's absolute and signed
means still agree after the fix (0.0481 both) because its remaining
residual (investigated below) is itself one-sided, unlike tpv13/27/30's.

`test.drv.a6.v2` (not registered in `testNameList.py`/`matrix.py`) was NOT
re-measured in this correction pass -- it is not part of this PR's gated
cell set and the earlier "measured directly, same magnitude" claim for it
is left as previously reported, unverified in this pass.

**drv.a6/drv.a6.v2's nonzero after-fix remainder (0.0481 m/s^2), investigated:**

- Hypothesis tested: fault-roughness-driven mesh distortion near the fault.
  REFUTED by measurement -- binning the after-fix residual by distance from
  the fault trace (|y|, dof<=3 nodes only) gives a FLAT mean|a_z| from the
  fault out to the far field (|y|<500 m: 0.0495; |y| in [500,8000) m:
  0.0482-0.0482; |y|>8000 m: 0.0481) -- no decay with distance, so this is
  not a fault-roughness artifact.
- Alternative explanation: `test.drv.a6`/`test.drv.a6.v2` build a
  DEPTH-VARYING density/velocity layer stack (`case_input/test.drv.a6/
  user_defined_params.py:51-68`, and the case comment there: "sets neither
  par.gamar nor par.roumax"), unlike tpv13/27/30's single uniform density.
  The gravity body force (`assembleGlobalKU.f90:21`) uses one global
  `roumax`, not the local per-element density. For a layered model this
  could leave a uniform, depth-independent residual at the surface layer
  regardless of fault roughness -- consistent with what was measured. **This
  explanation is UNVERIFIED** -- plausible from code reading and consistent
  with the flat-vs-distance measurement, but not itself tested (no planar
  variant of drv.a6 was built, and no direct per-element gravity-vs-density
  check was run). NOT fixed in this PR: it is a separate, pre-existing
  defect, independent of Defect C (present, smaller, even after this fix;
  would also have been present at the old offset, just harder to see under
  the larger offset-driven residual), out of this PR's scope, and neither
  drv.a6 nor drv.a6.v2 has a SCEC/FaultMod benchmark to gate it against.

## Rupture-time diff vs committed (pre-fix) reference, gate dx (500 m), both backends

Matched by rounded (x,y,z); `dt` computed only over nodes ruptured in BOTH
runs (sentinel 99999 = never ruptured excluded from the dt sample, but
included in the ruptured-fraction count). drv.a6's fortran row is the
OLD-vs-NEW committed `frt.canonical.txt` directly (reproducible via
`rupture_fraction_table.py`, this directory; output in
`rupture_fraction_output.txt`) -- the new reference was itself regenerated
from a fresh fortran run (commit `3338a8e`), so this is exactly what "the
fortran run vs the pre-fix reference" measures; a prior pass of this table
had arithmetic drift in that cell (-0.0392/+0.0002), corrected here. The
python-jax row is a separate measurement (an actual jax run vs the old
reference) and is not re-derived in this pass.

| case | backend | deep (z<-2km) dt_med / p90 (s) | deep Δruptured-frac | shallow dt_med / p90 (s) | shallow Δruptured-frac |
|---|---|---:|---:|---:|---:|
| tpv13 | fortran | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0000 | +0.0000 |
| tpv13 | python-jax | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0000 | +0.0000 |
| tpv27 | fortran | 0.0000 / 0.0000 | -0.0010 | 0.0000 / 0.0417 | +0.0000 |
| tpv27 | python-jax | 0.0000 / 0.0000 | -0.0010 | 0.0000 / 0.0417 | +0.0000 |
| tpv30 | fortran | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0417 | +0.0000 |
| tpv30 | python-jax | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0417 | +0.0000 |
| drv.a6 (self-cmp only, no benchmark) | fortran | 0.0417 / 0.2917 | +0.0000 | 0.2500 / 1.1667 | -0.0376 |
| drv.a6 (self-cmp only, no benchmark) | python-jax | 0.0417 / 0.2125 | -0.0015 | 0.2500 / 1.0250 | -0.0297 |

**Owner stop-gates (deep median |dt| > 0.02 s, or ruptured-fraction change >
1%), applied to tpv13/27/30 (the benchmark-gated cases): NOT triggered, on
either backend.** Deep medians are exactly 0.0000 s and fraction changes are
<=0.1% everywhere. This confirms, by measurement, the owner's expectation
that the 0.12 MPa offset (7.3 m x hydrostatic gradient) does essentially
nothing at seismogenic depth against tens of MPa of lithostatic stress.
drv.a6's self-comparison numbers are reported for completeness only, per
instruction -- NOT held to this stop-gate and NOT a benchmark verdict.

## Shallow-station (<2km) slip-rate / final-slip ratio, after/before, fortran, gate dx

| case | station | peak slip-rate ratio | final-slip ratio |
|---|---|---:|---:|
| tpv13 | faultst000dp013 (on-fault) | ~1.001-1.003 | ~1.001-1.003 |
| tpv13 | body*dp000 (off-fault, surface) | ~0.973-1.003 | ~0.973-1.003 |
| tpv27 | body*dp000 (off-fault, surface) | ~0.973-1.003 | ~0.973-1.003 |
| tpv30 | body*dp000 (off-fault, surface) | ~0.973-1.003 | ~0.973-1.003 |

All within a few percent, consistent with a small, physically-expected,
surface-localized correction.

## Elastic (Method-1, C_elastic==1) cases: byte-identity

**Case/cell counts, corrected:** 22 cases are registered in
`testNameList.nameList`. 4 are Method-2 (`C_elastic==0`: tpv13/tpv27/tpv30/
drv.a6, listed above); the other 18 are Method-1 (`C_elastic==1`), never
call `setPlasticStress`, and are structurally unaffected by this change. Of
those 18, `test.tpv22`/`test.tpv23` are `matrix.RELEASE_ONLY` (owner
decision 2026-10-02) and were NOT run in this everyday sweep -- they are not
part of the byte-identity claim below, which covers the remaining 16.

Confirmed via a full `python3 testsys/run.py all`
sweep on the fixed tree (wall clock 4428.4 s; `docs/perf_snapshots/
e2e_cells_2026-10-09_180840_3768237.json`): of 66 gated/declared cells, the
run.py-all e2e tier reported `e2e: FAIL - 7 cell(s) failed`, and ALL 7 are the
4 Method-2 cases (`test.tpv13` x {fortran,python-jax}, `test.tpv27` x
{fortran,python-jax}, `test.tpv30` x {fortran,python-jax}, `test.drv.a6` x
fortran) against the now-stale (pre-fix) reference -- expected, closed by the
reference-refreeze commit below.

**Byte-identity claim, downgraded to what the committed snapshot actually
records.** The snapshot (`e2e_cells_2026-10-09_180840_3768237.json`) records
only pass/fail (`"ok": true/false`) per cell, not a per-cell max&#124;diff&#124;
-- the earlier wording here ("`max|diff|=0.000000e+00`... 16 of 16") implied a
diff table that was never committed. What the snapshot actually gives,
read directly (`test.drv.a6`/`test.tpv13`/`test.tpv27`/`test.tpv30` x
fortran: `ok=false`, expected, stale reference at sweep time; every other
fortran cell -- `test.tpv8`, `test.tpv10`, `test.tpv104`, `test.tpv1053d`,
`test.meng2023a`, `test.meng2023cb`, `test.tpv29`, `test.tpv36`,
`test.tpv37`, `test.tpv35`, `test.tpv34`, `test.tpv26`, `test.tpv31`,
`test.tpv32`, `test.tpv33`, `test.tpv12`: `ok=true`): **16 of 16 elastic
fortran cells PASS** against their committed reference (`abs-max` gate,
`testsys/compare.py`) -- pass/fail, not a re-asserted byte-identity number.
`run.py all` is `python3 testsys/run.py all` run with no other changes on
the fixed tree (`run.py` itself only writes a
`docs/evidence/sweep-<sha>/summary.json` for `release`, not `all`).
`test.tpv27 x python-jax` additionally hit one unrelated,
non-reproduced infra flake on this run (`plotRuptureDynamics` exited 2 with
"file not found" against its own absolute path, immediately after an
identical invocation for `test.tpv31.python-jax` succeeded one line above it
in the log) -- re-run in isolation immediately after, it passed clean
(`max|diff|=1.63e-07` vs bound `1.0e-04`); not attributed to this fix.

After the reference-refreeze commit (below), `python3 testsys/e2e/run_e2e.py
--cases test.tpv13,test.tpv27,test.tpv30,test.drv.a6 --backends fortran` and
`--cases test.tpv13,test.tpv27,test.tpv30 --backends python-jax` both report
all cells SUCCESS (see that commit message for the per-cell numbers).

## 250 m resolution

Not run for this PR (owner decision 2026-10-09: "skip the 250 m sweep ...
the offset itself does not depend on dx, and at gate dx the change is zero
at output resolution").

## `test.drv.a6` x `python-jax`, re-checked against the re-frozen reference

`python3 testsys/e2e/run_e2e.py --cases test.drv.a6 --backends python-jax`
(2026-10-09, post audit-fix pass): **SUCCESS.** `nodes=5151
ruptured-in-both=1071 (matched-arrival=1018)`, `flips=326/450 [ok]
(existence=273: ref-only 138, run-only 135; timing-shift>1s=53)`,
`median|dfnft|=0.0417s/0.1s [ok]`, `phys_max=3.705779e+07/1.1e+09 [ok]` --
all three of `test.drv.a6`'s flip-budget gate components pass. `nc` and
`station` are DECLARED UNSUPPORTED for this cell (chaotic rupture-arrival
bistability, matrix.py, pre-existing and unrelated to this fix).

## Plastic-traction probe, re-measured

`testsys/parity/probe_plastic_traction.py` (report-only, never a gate) was
re-run 2026-10-09 against the current (post-removal, `DEPTH_OFFSET_M = 0.0`)
tree: **Tn median ratio 1.0018** (Ts median ratio 0.9859), unchanged from the
figure measured under the pre-removal offset. The offset is not load-bearing
for this ratio either way -- the probe's own 45-degree `str1ToFaultAngle`
collapse (`Tn = Syy = strVert`, depending on depth but not on the additive
constant's sign relative to the traction ROTATION) does not amplify or
cancel a uniform depth shift the way the step-1 vertical-equilibrium residual
above does.

## Both-ways regression test

`testsys/regression/test_defect_c_depth_offset_removed.py` (registered in
`testsys/ci_shard.py` shard 1): asserts the live depth expression in both
`src/fortran/meshgen.f90` and `src/python/eqdyna/meshgen.py` carries no
additive offset, and that the SAME checker rejects an in-memory fixture with
the historical offset re-inserted into either file (rule 14a both-ways).
Generalized 2026-10-09: the Python-side check (`PYTHON_DEPTH_RE`) now
matches ANY trailing `+ <number>` on the depth expression, not only the
literal `7.3215`, and the test carries a second both-ways fixture that
reinserts a DIFFERENT constant (`+ 3.0`, never historically present) to
prove the generalization actually catches a new offset, not just a verbatim
re-insertion of the old one. The Fortran-side check already caught any
additive constant (it compares the whole call-site argument against the
single expected expression, not just a `7.3215` substring), so only the
Python side needed the regex widened.

## Reference re-freeze sha256 (authoritative; corrects commit `3338a8e`)

Commit `3338a8e`'s message abbreviates each sha256 with a `...` elision for
readability, and one of those abbreviations is a copy/paste artifact (the
`test.tpv27` "new" `frt.canonical.txt` line shows a suffix that is actually
`test.tpv13`'s). Per the owner's instruction, `3338a8e`'s own message is NOT
rewritten (no history edit). The FULL, independently re-verified
(`sha256sum`, 2026-10-09, against the files as currently committed) values
are:

| case / file                          | sha256 (new, verified) |
|---------------------------------------|----------------------------------------------------------------|
| `test.tpv13` `fault.dyna.r.nc`        | `88398fb034ea8cbe533d63c95d31d74b22458119b775cae8bd715accc086fd7e` |
| `test.tpv13` `frt.canonical.txt`      | `3c31b2514a6369736f685ac8e75d30780d0eedadff875184eefa6703d2f12dc7` |
| `test.tpv27` `fault.dyna.r.nc`        | `75bcd2070f00d95ed5b5d2590945bbeed3d2fe12719d806f9ab00afdfc5588cb` |
| `test.tpv27` `frt.canonical.txt`      | `bce624bec7c0b381b056508be5c92ceb3765fa72527a2c5523e68cfc8fb6e097` |
| `test.tpv30` `fault.dyna.r.nc`        | `4e88cec91a6400cf9bd1dceb1dec8449da8375517609d3751319a5e7dbcef3a2` |
| `test.tpv30` `frt.canonical.txt`      | `e024a8a3f0b98c9f5a1a2737c45a75dce1498530477661436bd95602f6528625` |
| `test.drv.a6` `fault.dyna.r.nc`       | `6bd7408ae0ec67495595aecdef7c7e259d342d122ca14ce7f805fee8882bcb20` |
| `test.drv.a6` `frt.canonical.txt`     | `923268fc397e6ebb7e54b5cc1155d957ab77768a48994c16b20e552bdced8fa4` |

These 8 are the hashes of the files as they sit in this worktree right now
(`test.reference.results/<case>/{fault.dyna.r.nc,frt.canonical.txt}`); this
table, not `3338a8e`'s elided prose, is the citeable source for the PR body.

## Scripts committed with this evidence

- `measure_surface_residual.py` / `surface_residual_output.txt` -- the
  step-1 surface-residual table above.
- `rupture_fraction_table.py` / `rupture_fraction_output.txt` -- the
  ruptured-fraction/dt table's drv.a6-fortran correction above.
