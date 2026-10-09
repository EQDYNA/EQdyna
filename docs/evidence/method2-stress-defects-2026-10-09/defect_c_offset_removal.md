# Defect C: `setPlasticStress` depth-offset removal -- measured impact

Board PR #164 (2026-10-09). Owner decision (relayed, superseding item 24(d)'s
2026-09-16 "do what fortran does ... not to be derived or changed" ruling):
"it should be universal right?" / "we are setting stresses at the center of
each cell" -- the undocumented `+ 7.3215d0` shift in `setPlasticStress`'s
depth argument (`src/fortran/meshgen.f90:181`, mirrored verbatim in
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

All other registered cases (17 of 22 in `testNameList.nameList`) use
`C_elastic==1` (Method 1), never call `setPlasticStress`, and read exactly 0
for the removed term either way -- confirmed unaffected below.

## Step-1 surface residual acceleration (EQDYNA_DUMP_EQUIL=1, nt==1)

Mean |a_z| (m/s^2) over free-surface (z>-1e-6 m) non-fault, non-PML
(dof<=3) nodes, fortran backend, gate dx (500 m):

| case          | before (+7.3215) | after (0) | notes |
|---------------|------------------:|----------:|-------|
| test.tpv13    | 0.1616            | -0.0014   | ~100x drop |
| test.tpv27    | 0.1427            | ~0.0000   | fully removed (planar fault, uniform density) |
| test.tpv30    | 0.1427            | -0.0000   | fully removed (planar fault, uniform density) |
| test.drv.a6   | 0.1341            | 0.0408    | nonzero remainder -- see below |
| test.drv.a6.v2| 0.1341            | 0.0408    | measured directly (not inferred from drv.a6); same magnitude |

**drv.a6/drv.a6.v2's nonzero after-fix remainder (0.0408 m/s^2), investigated:**

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
included in the ruptured-fraction count).

| case | backend | deep (z<-2km) dt_med / p90 (s) | deep Δruptured-frac | shallow dt_med / p90 (s) | shallow Δruptured-frac |
|---|---|---:|---:|---:|---:|
| tpv13 | fortran | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0000 | +0.0000 |
| tpv13 | python-jax | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0000 | +0.0000 |
| tpv27 | fortran | 0.0000 / 0.0000 | -0.0010 | 0.0000 / 0.0417 | +0.0000 |
| tpv27 | python-jax | 0.0000 / 0.0000 | -0.0010 | 0.0000 / 0.0417 | +0.0000 |
| tpv30 | fortran | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0417 | +0.0000 |
| tpv30 | python-jax | 0.0000 / 0.0000 | +0.0000 | 0.0000 / 0.0417 | +0.0000 |
| drv.a6 (self-cmp only, no benchmark) | fortran | 0.0417 / 0.2917 | +0.0002 | 0.2500 / 1.1667 | -0.0392 |
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

17 of the 22 registered cases in `testNameList.nameList` never call
`setPlasticStress` (guarded by `if (C_elastic == 0)`) and are structurally
unaffected by this change. Confirmed via a full `python3 testsys/run.py all`
sweep on the fixed tree (wall clock 4428.4 s; `docs/perf_snapshots/
e2e_cells_2026-10-09_180840_3768237.json`): of 66 gated/declared cells, the
run.py-all e2e tier reported `e2e: FAIL - 7 cell(s) failed`, and ALL 7 are the
4 Method-2 cases (`test.tpv13` x {fortran,python-jax}, `test.tpv27` x
{fortran,python-jax}, `test.tpv30` x {fortran,python-jax}, `test.drv.a6` x
fortran) against the now-stale (pre-fix) reference -- expected, closed by the
reference-refreeze commit below. Every OTHER fortran cell in that same sweep
(`test.tpv8`, `test.tpv10`, `test.tpv104`, `test.tpv1053d`, `test.meng2023a`,
`test.meng2023cb`, `test.tpv29`, `test.tpv36`, `test.tpv37`, `test.tpv35`,
`test.tpv34`, `test.tpv26`, `test.tpv31`, `test.tpv32`, `test.tpv33`,
`test.tpv12`) reported `max|diff|=0.000000e+00` against its committed
reference -- exact byte identity, 16 of 16 elastic fortran cells, `run.py
all` is `python3 testsys/run.py all` run with no other changes on the fixed
tree (`run.py` itself only writes a `docs/evidence/sweep-<sha>/summary.json`
for `release`, not `all`; the raw per-cell table above is this run's own
evidence). `test.tpv27 x python-jax` additionally hit one unrelated,
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

## Both-ways regression test

`testsys/regression/test_defect_c_depth_offset_removed.py` (registered in
`testsys/ci_shard.py` shard 1): asserts the live depth expression in both
`src/fortran/meshgen.f90` and `src/python/eqdyna/meshgen.py` carries no
`7.3215` offset, and that the SAME checker rejects an in-memory fixture with
the historical offset re-inserted into either file (rule 14a both-ways).
