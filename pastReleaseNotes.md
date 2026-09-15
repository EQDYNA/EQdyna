# Past release notes\

# News in 2026
* 20260914 v5.6.0 release notes
  * New - supplied fault geometry (insertFaultType=3) is validated before use: grid, spacing, origin, row count, NaN/Inf, derivative columns and per-cell element-tangling offset, in Python at case.setup and again in Fortran at read time.
  * New - scripts/convertFaultGeometry resamples a supplied (x, z, y) surface onto a case's fault grid and validates its own output.
  * New - test.tpv29 ships the official surface at 50 m as well as 100 m, so the full-resolution tier runs from a clean checkout with no download.
  * Fix - refusing a run now actually fails: bare `stop` and `stop 'message'` both exit 0 under gfortran, and 13 sites used one of them, including the two mesh-alignment gates. All fatal paths go through src/errorCodes.f90 with named codes (1-125) and MPI_Abort, so a bad run ends instead of hanging other ranks.
  * Fix - the full-resolution tier raised KeyError before running any case; its TPV29 entry used the wrong key names.
  * Fix - geometry validator missed a globally rescaled surface (units error) and local corruption up to 3000x tolerance; both now caught or warned.
  * Change - no silent fallbacks (rule 2): par.dy is required rather than standing in as dx, and a fault grid too small to check is refused rather than checked weakly.
  * Docs - README carries a generated exit-code table that cannot drift from the source; measured TPV29 speeds at 100 m and 50 m on 48 ranks.
  * Full details: the v5.6.0 GitHub Release, git tag message, and pathway_forward.md.
* 20260914 v5.5.0 release notes
  * New - test.tpv29 (SCEC TPV29, official 25 m rough-fault geometry) gated as the 8th benchmark case: fast tier at dx=500 m/4 ranks (2,1,2)/~3 min, full-tier spec entry at dx=50 m/term=20 s; frozen reference in test.reference.results/test.tpv29.
  * Fix - fault-on-MPI-boundary bug: arn (fault nodal area) was double-counted whenever an MPI partition boundary coincided with the fault plane (symmetric y-domains), halving every on-fault traction term; fixed, with a hard-stop guard (checkFaultMPIAlignment) and a decomposition-invariance regression test.
  * Add - insertFaultType=3 (case-supplied fault geometry, e.g. TPV29's official surface) documented as a first-class mode; case.setup no longer invokes the geometry generator for it.
  * Fix - memory_estimate now reports cells/rank, cells total, ranks, and both memory figures separately (previously multiplied rank 0's count by the rank count, so it looked rank-invariant).
  * Full details: the v5.5.0 GitHub Release, git tag message, and pathway_forward.md.
* 20260914 v5.4.0 release notes
  * New - python/eqdyna/standalone: a fully standalone Python EQdyna (zero Fortran involved), now covering all five gated cases (tpv8, tpv104, tpv1053d, tpv10, drv.a6) with a committed acceptance tier (`testsys/run.py accept`) — four at roundoff-level agreement, drv.a6 under a documented chaos-aware criterion for its rupture-arrest bistability.
  * New - Drucker-Prager viscoplasticity ported to the standalone solver (test.drv.a6).
  * New - optional GPU execution (JAX/CUDA), verified via `testsys/run.py gpu`.
  * Add - testsys/ parity, accept, gpu, perf, and scaling tiers.
  * Fix - PML boundary tests standardized to inclusive bounds; mesh-time checkPMLAlignment guard added.
  * Fix - retired dead fric slot 6 (unused surface-pore-pressure reads).
  * Refactor - root cleanup: netCDF4 replaces xarray in check.test.py; orphaned testNameListWhole.py removed; shared loading/kernels modules factored out of the Python port.
  * Fix - scripts/lib.py: lazy imageio import, guarded optional dependencies.
  * Full details: the v5.4.0 GitHub Release, git tag message, and pathway_forward.md.

* 20260913 v5.3.6 release notes
  * Bug - PML region-14 corner damping used an x-bound for the y coordinate; fixed in all 4 sites. test.drv.a6 shear traction moved up to ~29 MPa; reference regenerated; guarded by a regression test.
  * Bug - thermal pressurization ran with shear-zone half-width h=0 since 2022 (fric_tp_h never assigned); now reads per-node fric(40)=0.02 m per the TPV105-3D spec. Verified against a clean-room rebuild of the SCEC-verified v5.2.0. tpv1053d reference regenerated.
  * Refactor - output-neutral cleanup, bit-gated: misc/ deleted; tabs to spaces; dead code removed; requireInputFile(); 6 files/subroutines renamed to verb-phrase names; globalvar.f90 restructured with a named FRIC_SLOT_* fric() slot map; func_lib.f90 shared helpers; MPI sendrecv consolidation; Python plotting dedup (~2,600 net lines removed).
  * New - python/: vectorized Python EQdyna (NumPy + JAX) covering slip-weakening, rate-and-state, and thermal-pressurization friction; parity and timing documented in python/README-parity.md (JAX-CPU 1.1-2.6x faster than serial Fortran).
  * New - tiered testing system testsys/ (unit, regression, e2e); CI gates on real exit codes.

* 20260909 v5.3.5 release notes
  * Add - `testsys/`, a tiered test system (unit/regression/e2e) with a single entry point `testsys/run.py [unit|regression|e2e|all]`; `testAll.py` is now a thin wrapper delegating to `testsys/e2e/run_e2e.py` (PROJECT_RULES.md rule 3). Gate: `python3 testsys/run.py unit` (33 tests, pure Python, no MPI/Fortran, <1s) then `regression` (2 guards, one per past incident, rule 10) then `e2e` (full create.newcase -> case.setup -> mpirun -> plotRuptureDynamics -> check.test.py pipeline against test.reference.results/, rule 7).
  * Fix - `scripts/lib.py`: `B2`/`B3`'s negative-`y` guard called `sys.exit()`, but the file did `from sys import *`, which does not bind the name `sys` -- the guard raised `NameError` instead of exiting. Changed to `import sys`. Covered by new unit tests `test_B2_negative_y_exits_cleanly_not_with_nameerror` / `test_B3_...` in `testsys/unit/test_lib.py`.
  * Fix - `check.test.py`'s `compare_nc_files` set `metadata_equal = f1.identical(f2)`, which requires bit-exact data in addition to matching attrs -- conflating the function's one calibrated 1e-3 threshold (rule 5) with an uncalibrated bit-exact check. This turned ordinary parallel-MPI floating-point non-determinism into a false `FAIL metadata` (observed: `test.tpv10/fault.dyna.r.nc`'s `shear_strike`, off by ~2.9e-8 between reruns, well under threshold). Changed the metadata check to compare variable sets and attrs only. Covered by new unit tests in `testsys/unit/test_check_comparisons.py` (within-threshold-but-not-bit-exact -> SUCCESS; differing attrs -> FAIL).
  * Fix - `testsys/e2e/run_e2e.py` dropped the `plotRuptureDynamics` step from the old `testAll.py`'s per-case sequence, so `fault.dyna.r.nc` (one of `check.test.py`'s `fileNameList` comparisons) was never regenerated. Restored: `create.newcase` -> `case.setup` -> `mpirun` -> `plotRuptureDynamics` -> `check.test.py`.
  * Fix - `install-eqdyna.sh` ran `chmod -R 755 scripts` on every build (PROJECT_RULES.md rule 13), flipping the mode bit on all ~60+ tracked files under `scripts/` -- including data files (`*.m`, `*.txt`, `*.mat`) never meant to be executable -- as a side effect of running the e2e test tier. Replaced with `chmod 755` of only the specific entry-point scripts (`case.setup`, `clean.py`, `create.newcase`, `generateFaultInterface`, `plotRuptureDynamics`, `plotSlipAndRPT`), matching their git-tracked exec bits. README's Installation section updated to match (`chmod 755 install-eqdyna.sh`, dropping the recursive `scripts` chmod); resolves pathway_forward.md item 6.
  * Add - `install-eqdyna.sh` to the exec-bit checks in `testsys/regression/test_create_newcase.py` (same guard family as `create.newcase`'s exec bit): it is invoked as `./install-eqdyna.sh` by `testsys/e2e/run_e2e.py`'s fresh-rebuild step, and a lost exec bit there fails a full e2e run the same way.
  * Add - `test.prev/` to `.gitignore` (byproduct of the e2e tier's rule-8 evidence-preservation step; was showing untracked in `git status`, rule 12).
  * Gate: build green; `testsys/run.py unit` 33/33 SUCCESS; `regression` 2/2 SUCCESS; `e2e` 5/5 cases, `check.test.py` 15/15 comparisons SUCCESS. Verified on ubuntu 22.04, gfortran 11.4.0, Open MPI 4.1.1, at commit c4b13e5 + this release's changes.

* 20260909 v5.3.4 release notes
  * Fix - uninitialized `thetaPcTmp` in `NewtonRaphson` (src/faulting.f90) for friclaw==5: the copy-back `thetaPc0 = thetaPcTmp` ran unconditionally but `thetaPcTmp` was only assigned inside the friclaw<5 branch, so fric(23)/frt.txt col 22 held a deterministic uninitialized stack value every step. Physics for friclaw==5 is unaffected (theta_pc is not used); frt.txt cols 1-21 are bit-identical before and after. test.tpv1053d golden reference regenerated from the fixed binary and doubles as the regression test.
  * Fix - src/makefile: FC on ubuntu changed from mpif90.mpich to mpif90; mpif90.mpich does not exist on the target system and the ubuntu build was broken.
  * Fix - restored the executable bit on scripts/create.newcase and scripts/generateFaultInterface. Without it, create.newcase was PATH-shadowed by an unrelated project's script of the same name and the test suite could not run from a clean checkout.
  * New - output_src_evol subroutine (src/library_output.f90) writes per-fault-node final slip rate to binary src_evol files, for on-fault state visualization/AI use.
  * Change - output_gm and output_src_evol now run every time step (previously every 10th step via mod(nt,10)==1) in src/driver.f90.
  * Add - PROJECT_RULES.md, a 14-rule project rule book, and pathway_forward.md (formerly docs/PROJECT_STATUS.md), a living status board of open issues, re-check intervals, and a tasks-done log.
  * Add - .gitignore for bin/, __pycache__/, *.mod, *.pyc, scratch/, src/eqdyna, test/.
  * Update - scripts/plotRuptureDynamics: fixed duplicate subplot-axis variable names and added axis labels.
  * Update - scripts/clean.py: also purge src_evol and src* output files.
  * Update - README.md: reworded collaboration and benchmark-performance sections, added DOI links to references.
  * Known issue - src/netcdf_io.f90:76-105 hardcodes fault index 1 instead of the loop variable in on-fault netcdf-input assignment, so multi-fault (ntotft>=2) on-fault netcdf input is not applied correctly. No current test case exercises ntotft>=2. Tracked in pathway_forward.md item 7.

# News in 2024
* 20241006 v5.3.3 release notes:
  * New - verification against new SCEC/USGS Spontaneous Rupture Code Verification benchmarks [TPV36&37](https://strike.scec.org/cvws/tpv36_37docs.html) for 15 deg shallow dipping thrusting. 
  * New - exclusive model setup python script user_defined_params.py for TPV36&37 are under /case_input/test.tpv36 and /case_input/test.tpv37, respectively. 
  * Performance: 512 cores are used for 50 m resolution TPV36 on Lonestar6 at TACC using 4 hours and 40 minutes. 
  * New - previous feature of degeneration of hexahedrons (Hughes, 2000) for complex fault geometry is incorporated in the new EQdyna architecture with TPV36&37.
  * New - autotesting workflow is added on GitHub for developers. 
  * New - supporting MacOS (M3 chip tested).
  * Refctor - rename file and subroutine names for clarification.
  * Reference: Hughes, 2000, The Finite Element Method: Linear Static and Dynamic Finite Element Analysis, Dover Publications.  

* 20240221 v5.3.2 release notes
  * New - new link for EQdyna. https://github.com/EQDYNA/EQdyna.git
  * New - new organization EQDYNA is created. 
  * Bug - dz for dipping fault. Verified against TPV10.  
  * Bug - avoid int() in Python and FORTRAN, use round() or nint() instead.
  * Add - MATLAB scripts for GM postprocessing. 
  * Add - str1ToFaultAngle and devStrToStrVertRatio for assigning stresses for plastic models.
  * Add - case test.drv.a6.v2 
  * Refactor - move adjustable parameters out of EQdyna. 
  * Refactor - rename file and function names to reflect their intents, for easy search. 
  * Refactor - driver for seismic wave propagation + faulting only. 
  * Refactor - MPI communications for nodal quantities, driver.f90, depreciate PMLwhg.f90 and contm.f90, refactor qdct3.f90, rename qdct3 to ku, offFaultStationSCEC. 
  * Change - Positive dip angles for faults tilting to y+.
  * Change - Use empirical estimate for memory usage. 
  * Update - ubuntu.env.sh

# News in 2023
* 20231209 v5.3.1 release notes
  * Python utility generateFaultInterface is created to generate fractal and dipping fault interface;
  * update the test system to compare the numerical residuals between test and reference results; 
  * support a new test test.tpv10, a 60 deg dipping fault;
  * support a new test test.drv.a6, a fractal fault plastic model for ground motion application;
  * refactoring and bug fixes. 

* 20231130 v5.3.0 release notes
  * refactor faulting.f90 and meshgen.f90;
  * introduce a system for quick testing by running testAll.py; 
  * test system now supports tpv8, tpv104, tpv1053d, tpv1053d.6c, meng2023a, meng2023cb;
  * update input parameter system with a single defaultParameters.py and customized user_defined_params.py;

* 20230327 
  * *```EQdyna```* works on Ubuntu now. 
  * install-eqdyna.sh can support multiple systems. 
  * test-all.sh will create and run four pre-defined cases with 4 CPUs with a few minutes. 
* 20230321 *```eqdyna.docker```* is published via [dunyuliu/eqdyna.docker](https://hub.docker.com/repository/docker/dunyuliu/eqdyna.docker/general) 

# News in 2022
* A set of python utlities that make *```EQdyna```* much easier to use.
* Thermal pressurization implemented and benchmarked against TPV105-3D.