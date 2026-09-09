# Past release notes\

# News in 2026
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