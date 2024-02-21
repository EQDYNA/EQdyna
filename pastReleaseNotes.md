# Past release notes

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