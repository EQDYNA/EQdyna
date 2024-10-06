# Past release notes\

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