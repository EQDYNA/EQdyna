# TPV37 gating notes (mira-volkov, worktree agent-a19d02377f7a99d91)

## Step 0 -- pre-flight (rule 17 step 3)
Confirmed the Fortran TPV37 branch already exists:
  src/fortran/faulting.f90:405 `if (TPV == 201 .or. TPV==36 .or. TPV==37 .or. TPV==29)`
Confirmed the Python port already handles it:
  src/python/eqdyna/faulting.py:109 `if TPV in (29, 36, 37, 201):`
Diffed case_input/test.tpv36 vs case_input/test.tpv37/user_defined_params.py:
  only par.tpv (36->37) and on_fault_vars[iz,ix,4] taper coefficient
  (0.0005e6 -> 0.001875e6) differ. Everything else (dip=15, dx=500, term=6,
  C_degen=par.dip, mesh layout, nx=2 (4 ranks)) is identical to tpv36.

## Step 1 -- coarse gate, fortran, 4 ranks (rule 17 step 4)
Ran outside the matrix framework (bootstrap -- matrix.py's import-time check
needs a reference to exist before a case can be added to CASES):
  create.newcase <scratch>/test.tpv37 test.tpv37
  cd <scratch>/test.tpv37 && python3 case.setup   # generated run.sh (mpirun -np 4)
  python3 clean.py
  mpirun -np 4 <worktree>/bin/eqdyna              # ~4 min wall, 4 ranks
  python3 plotRuptureDynamics
Output: frt.txt0..3, fault.dyna.r.nc (349723 bytes -- byte-size-identical to
tpv36's own reference, expected since only tpv-branch and initial stress
differ, mesh geometry is identical), cRuptureDynamics.png.
Sanity check (visual, cRuptureDynamics.png): rupture nucleates centrally,
propagates outward with elliptical rupture-time contours, Mag 6.78, peak
slip rate up to 5.6 m/s, dip-slip dominated (Slip d up to 1.75 m vs Slip s
~0.02 m) -- consistent with TPV37's dipping-fault geometry. Sane.

## Step 2 -- freeze reference (rule 17 step 4)
  python3 -m testsys.frt_canonical <scratch>/test.tpv37
  -> 4 rank files -> 3477 canonical rows (IDENTICAL row count to tpv36's
     frt.canonical.txt -- expected, same mesh/decomposition).
Copied into test.reference.results/test.tpv37/:
  frt.canonical.txt, fault.dyna.r.nc, cRuptureDynamics.png
(matches the file set tpv36's reference dir carries).

## Step 3 -- registration (rule 17 step 4)
testNameList.py: appended 'test.tpv37' to nameList, 4 to coreNumList.
testsys/matrix.py: CASE_BOUND['test.tpv37']=1e-6, GATE['test.tpv37']='abs-max',
  following tpv36's exact precedent/comment style. Provisional bound (same
  convention tpv36 used before its own fresh numpy/jax numbers were in hand);
  updated with actual measured values after step 5 below.

## Step 1b -- fortran e2e cell, via run_e2e.py (bootstrap confirmation)
  python3 testsys/e2e/run_e2e.py --cases test.tpv37 --backends fortran
  -> SUCCESS, max|diff|=0.000000e+00 bound=1.0e-06, 3477 fault nodes compared,
     nc: SUCCESS (fault.dyna.r.nc vs threshold=1e-3), 220.0s wall.
  Bit-exact, as expected (Fortran side is completely unchanged by this work).

## Step 4 -- independent validation (rule 17 step 6)
Ran `testsys/parity/evidence_tpv36_37_vs_published.py` (already wired for
both cases) against a FRESH local run of both test.tpv36 and test.tpv37
(create.newcase+case.setup+4-rank mpirun+plotRuptureDynamics, same scratch
tree), pointed at the owner's published SCEC submission
(`scec_archive/tpv36|tpv37/eqdyna-v5.3.3-50m-2024`, copied read-only from the
shared checkout's untracked, gitignored cache -- same bytes, not regenerated).
tpv36 fresh local run's frt.canonical.txt diffed BYTE-IDENTICAL against the
committed test.reference.results/test.tpv36/frt.canonical.txt (sanity check
that the "fresh" side is the same physics as what's already gated).

RESULT (both cases): 19 common stations. The factor-of-2 "arn doubling" check
-- INITIAL down-dip shear (a pure INPUT, resolution-independent) -- passes
cleanly: ratio local(500m)/published(50m) = median 1.0000, min 1.0000, max
1.0000 (n=15 stations with nonzero initial dd-shear), for BOTH tpv36 and
tpv37. No unit/factor bug in the loading construction for either case.

## Step 5 -- all three backends (rule 17 step 7)
  python3 testsys/e2e/run_e2e.py --cases test.tpv37 --backends fortran
    -> SUCCESS, max|diff|=0.000000e+00, bound 1e-6, 220.0s, nc SUCCESS.
  python3 testsys/e2e/run_e2e.py --cases test.tpv37 --backends python-numpy,python-jax
    -> test.tpv37 x python-jax   SUCCESS  271.7s   max|diff|=5.256410e-09
    -> test.tpv37 x python-numpy SUCCESS 1269.5s   max|diff|=6.929140e-09
    (both ran CONCURRENTLY -- numpy's wall time is inflated by CPU
    contention with jax sharing the box, not a real 4.7x numpy/jax gap;
    the important number is max|diff|, both far inside the 1e-6 bound.)
  Updated CASE_BOUND['test.tpv37'] comment in testsys/matrix.py with these
  fresh measured numbers (was provisional 1e-6-by-convention before this).
  ALL THREE BACKENDS PASS. No UNSUPPORTED cells, no feature gaps -- tpv37
  needed no port work beyond what tpv36's own C_degen/wedge landing already
  provided (confirmed at the top of this file: both faulting.f90 and
  faulting.py already branch on TPV==37 explicitly).

## Housekeeping alongside registration
Widening test.tpv37's fortran cell into matrix.CI_CELLS (automatic, same
mechanism as tpv36) needed the same follow-up tpv36's own re-gating commit
made: added test.tpv37 to `.github/workflows/test.yml`'s e2e-ci-fortran-a
job (Group A) and updated the "20 of 27"/"27 cells" cell-count comments in
CLAUDE.md and testsys/matrix.py to 21 of 30 / 30 cells.

## Step 6 -- unit + regression tier
  python3 testsys/run.py unit regression
  -> SUCCESS unit (exit 0), SUCCESS regression (exit 0), exit code 0.
  (test_ci_workflow_coverage.py is part of the regression tier and the tier
  reported SUCCESS overall -- run.py fails non-zero on any regression-test
  failure, so this covers the CI_CELLS/workflow-coverage check even though
  the individual per-file line for it fell outside the captured tail.)

## Step 7 -- prove tpv36 unmoved (all 3 backends, fresh)
  python3 testsys/e2e/run_e2e.py --cases test.tpv36 --backends fortran,python-numpy,python-jax
    -> test.tpv36 x fortran      SUCCESS   201.2s  max|diff|=0.0 (bit-exact)
    -> test.tpv36 x python-jax   SUCCESS   234.1s  max|diff|=5.867293e-09
    -> test.tpv36 x python-numpy SUCCESS  1772.0s  max|diff|=7.894064e-09
  Both python numbers are EXACT reproductions of the numbers already
  committed in testsys/matrix.py's CASE_BOUND['test.tpv36'] comment and in
  pathway_forward.md's item 19(a) closure entry (5.867293e-09 /
  7.894064e-09) -- proves tpv36's own gate did not move at all from this
  registration, only additions were made to testNameList.py/matrix.py.
  (numpy wall time inflated by CPU contention with the concurrently-running
  jax cell, same effect observed for tpv37 in step 5 -- not a regression.)

## DONE -- summary
All rule-17 steps complete. tpv37 registered on all 3 backends, bit-exact
on fortran, roundoff-level agreement on numpy/jax, independent validation
passed cleanly, tpv36 proven unmoved, unit+regression green.

Final down-dip slip and peak dd-shear differ station-by-station between the
500 m coarse run and the 50 m published run (e.g. tpv37 faultst000dp060:
pub -1.047 m vs loc -0.290 m final dd-slip) -- expected, not a defect: this
is the SAME coarse-vs-fine resolution gap already documented in tpv36's own
CASE_BOUND comment ("the published tpv36/37 standings rank a run this coarse
last and improve monotonically with resolution"). The independent-validation
step is passing the check it exists to catch (a doubling/halving in the
initial loading), not claiming SCEC-standings-level accuracy at 500 m.
