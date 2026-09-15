# News in 2026
* 20260915 v5.7.1 release notes
  * Fix - v5.7.0's CI went red: the Python/JAX backend inside e2e needs 10.4 GB for test.tpv104 and a GitHub runner has 7 GB, so the job was SIGTERM'd (exit 143) with no output. CI now gates test.tpv8 (measured 3.59 GB), and e2e prints its case list so a green check states its own coverage instead of implying five.
  * Fix - the e2e python child runs unbuffered; under the previous buffering a resource kill produced zero output, so the log showed seven minutes of silence and a bare exit code.
  * Change - release workflow (rule 15): the tag and GitHub Release now come AFTER CI is green on the pushed commit, not before. v5.7.0 was tagged on a local green that could not model the runner's memory; a released tag pointing at a red commit makes the Releases page the authoritative wrong answer.
  * Note - all five accept cases run and pass locally. The exclusion is the runner's memory, not the cases. Python peak RSS on test.tpv104 is 10.4 GB (jax) / 4.2 GB (numpy) against Fortran's 0.95 GB; reducing that is tracked in pathway_forward.md.
  * For past release notes, please refer to pastReleaseNotes.md.





# Introduction to *```EQdyna```*

*```EQdyna```* is a parallel finite element software to simulate earthquake spontaneous dynamic rupture, seismic wave propagation and high frequency deterministic ground motions. It has a focus to simulate earthquakes on geometrically complex fault systems with other heterogeneities in mind. The core is written in FORTRAN90 with a set of python utilities to setup cases and allocate HPC resources.  <br/>

*```EQdyna```* is highly efficient and scalable due to its adoptions of built-in mesh generation, explicit time integration, under-integrated hexahedral elements, degenerated wedge elements for complex fault geometries, hourglass controls, perfectly matched layer and 3D parallelization with MPI. <br/> 

Other features include 
* 3D velocity structure to simulate basin effects.
* Frequency independent Q by coarsed-grained memory scheme for seismic attenuation. 
* Various frictional constitutions such as slip-weakening and various forms of rate- and state-friction.
* Normal stress evolutions by an additional state variable.
* Drucker-prager off-fault viscoplasticity.
* Dynamic relaxation for earthquake cycle applications, etc.

*```EQdyna```* has been extensively verified against benchmark problems from [SCEC/USGS Spontaneous Rupture Code Verification Project](https://strike.scec.org/cvws/).

*```EQdyna```* is also part of the fully dynamic earthquake cycle simulator *```EQsimu```* [(*Liu et al., 2020, GJI*)](https://www.researchgate.net/publication/346814142_EQsimu_a_3-D_finite_element_dynamic_earthquake_simulator_for_multicycle_dynamics_of_geometrically_complex_faults_governed_by_rate-_and_state-dependent_friction).

# Verification & Benchmarks

EQdyna is verified against SCEC/USGS Spontaneous Rupture Code Verification
benchmarks (https://strike.scec.org/cvws/). The fast tier gates every
commit at coarse resolution against frozen references
(`python3 testsys/run.py e2e`); the full tier reproduces each gated
benchmark at its official spec resolution and duration, report-only, on
16 ranks (`python3 testsys/run.py e2e-full`).

The fast tier is a single sweep over `case x backend`: each case below is run
on the Fortran solver (MPI) and on the standalone Python solver (NumPy and
JAX), and every cell is compared against the same committed reference for that
case, at that case's one tolerance. The run prints which cells it covered and
which it did not, so a green result states its own scope. CI runs the subset
that fits a 7 GB runner (`python3 testsys/run.py unit regression e2e-ci`);
`python3 testsys/run.py all` runs the whole sweep.

| case | physics | SCEC benchmark | fast tier (dx/cores/~time) | full-tier spec (dx/term) |
|---|---|---|---|---|
| [test.tpv8](case_input/test.tpv8/README.md) | strike-slip, slip-weakening | [TPV8](https://strike.scec.org/cvws/tpv89docs.html) | 500m/4/~48s | 100m/15s |
| [test.tpv10](case_input/test.tpv10/README.md) | dipping normal fault | [TPV10](https://strike.scec.org/cvws/tpv10_11docs.html) | 500m/4/~48s | 100m/15s |
| [test.tpv104](case_input/test.tpv104/README.md) | strike-slip, rate-and-state | [TPV104](https://strike.scec.org/cvws/tpv103_104docs.html) | 500m/4/~48s | 50m/12s |
| [test.tpv1053d](case_input/test.tpv1053d/README.md) | RSF + thermal pressurization | [TPV105-3D](https://strike.scec.org/cvws/tpv105_3D_docs.html) | 500m/4/~48s | excluded (no spec dx) |
| [test.drv.a6](case_input/test.drv.a6/README.md) | fractal fault + plasticity | internal | 500m/4/~48s | excluded (no published spec) |
| [test.meng2023a](case_input/test.meng2023a/README.md) | layered velocity structure | internal | 400m/4/~48s | excluded (no published spec) |
| [test.meng2023cb](case_input/test.meng2023cb/README.md) | layered velocity, multi-patch | internal | 400m/4/~48s | excluded (no published spec) |
| [test.tpv29](case_input/test.tpv29/README.md) | strike-slip, fractal rough fault | [TPV29](https://strike.scec.org/cvws/tpv29_30docs.html) | 500m/4/~3 min | 50m/20s |

# Environment
*Optional (Python solver on GPU)*: `pip install "jax[cuda12]"` — the standalone Python solver (`python -m eqdyna.standalone`) then runs on NVIDIA GPUs; verify with `python3 testsys/run.py gpu`.

*```EQdyna```* requires <br/>
  - FORTRAN compiler (gfortran/intel FORTRAN)
  - MPI (mpich/intel MPI)
  - netCDF (libnetcdf libnetcdff)

Pre-staging and post-processing require Python packages <br/>
  - Python3
  - numpy>=1.20 (as of 20240201, matplotlib requires numpy>=1.20)
  - matplotlib
  - xarray
  - netCDF4
```
bash ubuntu.env.sh
```
will install the required packages through apt-get and pip on Ubuntu 22. <br/>

# Installation
```
git clone https://github.com/EQDYNA/EQdyna.git
cd EQdyna
chmod 755 install-eqdyna.sh
./install-eqdyna.sh -m ubuntu # ubuntu/ls6/macos
export EQDYNAROOT=$(pwd)
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
python3 testsys/run.py all # the sweep: 8 cases x 3 backends (fortran, python-numpy,
                           # python-jax) against one canonical reference each.
                           # Cells run concurrently; ~730 s on a 64-core box.
                           # It prints which cells it ran, so a pass states its scope.
```
For bash, please insert the following lines in .bashrc
```
export EQDYNAROOT=/path/to/EQdynaRootDirectory
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
```

# Quick Start Guide
Three steps are needed to run a new case. <br/>
```
create.newcase $caseDirectoryName $predefinedCompset
cd $caseDirectoryName
./case.setup
bash run.sh # or ./case.submit to submit batch job on LS6.
```
Replace $caseDirectoryName with the directory name you want to create. <br/>
Replace $predefinedCompset with one of the following supported compsets <br/>
* test.drv.a6, for determinisitc ground motion with fractal fault and plasticity
* test.tpv8
* test.tpv10
* test.tpv36
* test.tpv37
* test.tpv104
* test.tpv1053d
* test.tpv29 <br/>
[TPV+number is the naming convention of [SCEC/USGS Spontaneous Rupture Code Verification Excercise](https://strike.scec.org/cvws/).] <br/>

For a customized case, please choose the most relevant predefined compset and modify ```user_defined_param.py``` accordingly. <br/>

# Benchmark computational performance and resource
* TPV36 & 37: 4.7 hours for 50 m resolution using 512 CPUs on Lonestar6 at TACC. <br/>
* TPV104: 0.4 hours for 15-sec simulation (1875 time steps) using 40 CPUs on Lonestar6.  <br/>

# Exit codes
<!-- BEGIN EXIT CODES (generated from src/errorCodes.f90; do not edit by hand) -->

When a run is refused or fails, EQdyna prints a `FATAL` block naming the code
and the reason, and exits with that code. Codes are kept in 1-125 so the number
in the source is the number the shell reports (a status above 255 wraps).

How faithfully the number survives depends on the launcher. `mpirun` (Open MPI,
hydra) reports it directly. `srun` reports the **maximum** status across tasks,
so a straggler killed while `MPI_Abort` tears the job down yields 137 or 143 and
masks the code; `sbatch` reports the wrapper script's status unless the script
ends with `exit $?`. On ls6 and grace the launcher is `ibrun` -> `srun`. So treat
a **non-zero status** as the reliable signal and the printed `FATAL` block as
authoritative; the specific number is advisory under `srun`.

**Generic** (1-9)

| code | name | meaning |
|-----:|------|---------|
| 1 | `ERR_GENERIC` | unclassified fatal error |

**Configuration and parameters** (11-19)

| code | name | meaning |
|-----:|------|---------|
| 11 | `ERR_CFG_Q_NEEDS_ELASTIC` | C_Q=1 requires C_elastic=1 |
| 12 | `ERR_CFG_Q_NEEDS_UNIFORM` | C_Q=1 requires rat=1.0 (uniform elements) |
| 13 | `ERR_CFG_PLASTIC_OUTPUT` | output_plastic=1 requires C_elastic=0 |

**Input files** (21-29)

| code | name | meaning |
|-----:|------|---------|
| 21 | `ERR_INPUT_FILE_MISSING` | a required FE_*.txt / data file is absent |

**Fault geometry** (31-39)

| code | name | meaning |
|-----:|------|---------|
| 31 | `ERR_GEOM_ROUGH_INVALID` | bFault_Rough_Geometry.txt does not match this mesh |

**Mesh generation and element quality** (41-49)

| code | name | meaning |
|-----:|------|---------|
| 41 | `ERR_MESH_STRESS_ARR_SMALL` | sizeOfStressDofIndexArr exceeds 5*sizeOfEqNumIndexArr |
| 42 | `ERR_MESH_COUNT_MISMATCH` | meshgen's node/element/equation tallies disagree |
| 43 | `ERR_MESH_EQNUM_MISMATCH` | eqNumIndexArrLocTag /= sizeOfEqNumIndexArr |
| 44 | `ERR_MESH_FAULT_MISMATCH` | nftnd0 /= nftnd (meshgen vs countMeshEntities) |
| 45 | `ERR_MESH_MULTIFAULT_MSNODE` | master-node construction cannot handle ntotft>1 |
| 46 | `ERR_MESH_BAD_WEDGE` | degenerate wedge built with unequal node ids |
| 47 | `ERR_MESH_BAD_JACOBIAN` | non-positive Jacobian determinant (inverted element) |
| 48 | `ERR_MESH_MATERIAL_UNSET` | an element has no material property assigned |

**MPI and domain decomposition** (51-59)

| code | name | meaning |
|-----:|------|---------|
| 51 | `ERR_MPI_FAULT_ALIGNMENT` | a rank boundary in y coincides with the fault plane |

**Numerics and runtime state** (61-69)

| code | name | meaning |
|-----:|------|---------|
| 61 | `ERR_NUM_PML_ALIGNMENT` | element centre lies exactly on a PML bound |
| 62 | `ERR_NUM_PML_DAMPING` | negative PML damping vector component |
| 63 | `ERR_NUM_VELOCITY_NAN` | NaN velocity during time stepping |
| 64 | `ERR_NUM_NEGATIVE_DEPTH` | negative depth passed to a B-function |

**External libraries** (71-79)

| code | name | meaning |
|-----:|------|---------|
| 71 | `ERR_NETCDF` | a NetCDF call returned an error |

Exit status 0 means the run completed. Every fatal path goes through
`abortRun` in `src/errorCodes.f90`, which calls `MPI_Abort` so the whole job
ends instead of one rank stopping while the others block in a collective.

<!-- END EXIT CODES -->

# Collaboration
We try our best to make *```EQdyna```* easy to use but doing great science is our priority. We welcome collaborations, comments and suggestions. Please reach out to Drs Benchun Duan (bduan@tamu.edu) and Dunyu Liu (dliu@ig.utexas.edu).

# Reference
1. [Duan and Oglesby (2006)]( https://doi.org/10.1029/2005JB004138). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults, J. Geophys. Res., 111.
2. [Duan (2012)](https://doi.org/10.1029/2011JB009124). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res., 117.
3. [Luo and Duan (2018)](https://doi.org/10.1029/2017JB015320). Dynamics of non-planar thrust faults governed by various friction laws, J. Geophys. Res. Solid Earth, 123.
4. [Liu and Duan (2018)](https://doi.org/10.1785/0120170374). Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure, Bull. Seismol. Soc. Am., 108(4).
