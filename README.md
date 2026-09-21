# News in 2026
* 20260921 v5.14.0 release notes
  * New - **the jax backend now runs under real MPI**: one process per rank, elements split across ranks, halo exchange at `driver.f90:27`'s position -- `mpirun -np N python3 -m eqdyna <case_dir> --backend jax --mpi`. Measured FASTER than Fortran at every rank count measured on `test.tpv104` (per-step by difference, Fortran measured in the same session): 1 rank 795.05 vs 1034.50 ms/step, 2 ranks 499.49 vs 504.83, 4 ranks 192.60 vs 250.17 -- 1.30x / 1.01x / 1.30x in the gated `halo` sync mode. Scaling tracks Fortran's curve (3.8x at 4 ranks against Fortran's same-session 4.14x), where the alternative `shard_map` route flattened at 2.76x on 16 devices.
  * New - **the e2e sweep gains a fourth backend-axis value, `python-jax-mpi`**, per-case opt-in via `matrix.PY_MPI_RANKS` with `test.tpv8` at 4 ranks as the first entry. It compares through the existing canonicalisation against the SAME committed reference at the case's existing bound, and it HARD-FAILS unless it counted the expected number of per-rank `frt.txt<rank>` files -- a launch that silently starts one rank would otherwise compare green.
  * Known limitation - **8 and 16 ranks were NOT measured and the ratio bar is UNANSWERED.** A foreign tenant held the box at loadavg 94-100 with 63 of 64 cpus busy; the 8-rank attempt under that load recorded EFFECTIVE_CORES 0.22-0.60 per rank and is kept labelled as a record of the tenancy, not of the solver. The ~12-14x-at-16-ranks question is not claimed from a 4-rank point; the re-run command and its precondition (~16 cpus under 0.3 busy) are on the board, `pathway_forward.md` item 43.
  * Note - **`halo` is the gated sync mode; `allreduce` is a measurement mode only**, because MPI does not promise reduction order and the gate must not compare a reduction-order-dependent result.
  * Note - no physics changed and no reference was regenerated. Also contains today's earlier master commits: `e83c02e` (rule-book drift), `732331c` (board header and four stale claims) and `e9c3fa2` (the dropped-bump fix and its two-directional release guard).
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
16 ranks. These are multi-hour runs, so the tier REFUSES to start unless you
opt in explicitly:

```
EQDYNA_FULL_LAUNCH=yes-hours python3 testsys/run.py e2e-full
```

Without the variable it prints what it would run and exits non-zero, so the
tier can never be triggered by accident (for example by `run.py all`).

Each case reaches its spec resolution by DECIMATING a shipped fault surface --
every value is an official one, no interpolation -- so a clean checkout needs
no downloads. test.tpv29 ships its surface at 100 m (3.5 MB) and 50 m (14 MB);
a request for a dx finer than the shipped source, or one that is not a
multiple of it, is refused rather than interpolated (`scripts/lib.py`'s
`requireFaultGeometryResolution`). 25 m needs the official SCEC file and
`tpv29GeometryTools.convertOfficial25m`.

The fast tier is a single sweep over `case x backend`: each case below is run
on the Fortran solver (MPI) and on the standalone Python solver (NumPy and
JAX), and every cell is compared against the same committed reference for that
case, at that case's one tolerance. The run prints which cells it covered and
which it did not, so a green result states its own scope. CI runs the subset
this table declares memory-safe (`testsys/matrix.py`'s `CI_CELLS` /
`MEASURED_PEAK_RSS_GB`, 25 of 30 cells as of 2026-09-21 -- `len(matrix.CI_CELLS)`
is the number, and `testsys/matrix.py`'s own closing comment states it), split
across parallel jobs so each cell group gets its own 7 GB runner rather than
sharing one, and so no case queues behind another that does not need to
(`.github/workflows/test.yml`: `build`, then `unit-regression`
[`testsys/run.py unit regression`]; `e2e-ci-fortran-a`
[`--backends fortran --cases test.tpv29,test.tpv1053d,test.tpv104,test.tpv8,test.tpv37`]
and `e2e-ci-fortran-b`
[`--backends fortran --cases test.drv.a6,test.meng2023a,test.meng2023cb,test.tpv10,test.tpv36`]
-- 5 of the 10 cases each;
`e2e-ci-python-cheap` [tpv8 both backends + tpv10/tpv104/tpv1053d x jax +
tpv36/tpv37 both backends];
`e2e-ci-python-meng` [meng2023a/meng2023cb both backends + tpv29 x jax]; and
`e2e-ci-python-tpv29` [tpv29 x python-numpy alone -- 1711 s measured locally,
an order of magnitude past every other python cell, isolated so its
irreducible cost doesn't queue behind or ahead of anything else] -- all in
parallel after `build`) -- the same total coverage as one local invocation of
`python3 testsys/run.py unit regression e2e-ci`. `python3 testsys/run.py all`
runs the whole sweep. Still excluded from CI: `test.drv.a6` on either python
backend (9.57 GB measured on jax; numpy never measured) and
`test.tpv10`/`test.tpv104`/`test.tpv1053d` on python-numpy (never measured --
their jax columns are the ones now in CI).

| case | physics | SCEC benchmark | fast tier (dx / ranks / fortran s / numpy s / jax s) | full-tier spec (dx/term) |
|---|---|---|---|---|
| [test.tpv8](case_input/test.tpv8/README.md) | strike-slip, slip-weakening | [TPV8](https://strike.scec.org/cvws/tpv89docs.html) | 500m/4/13.3/74.3/20.4 | 100m/15s |
| [test.tpv10](case_input/test.tpv10/README.md) | dipping normal fault | [TPV10](https://strike.scec.org/cvws/tpv10_11docs.html) | 500m/4/39.0/331.3/98.5 | 100m/15s |
| [test.tpv104](case_input/test.tpv104/README.md) | strike-slip, rate-and-state | [TPV104](https://strike.scec.org/cvws/tpv103_104docs.html) | 500m/4/33.6/309.5/86.1 | 50m/12s |
| [test.tpv1053d](case_input/test.tpv1053d/README.md) | RSF + thermal pressurization | [TPV105-3D](https://strike.scec.org/cvws/tpv105_3D_docs.html) | 500m/4/51.0/368.3/98.1 | excluded (no spec dx) |
| [test.drv.a6](case_input/test.drv.a6/README.md) | fractal fault + plasticity | internal | 500m/4/96.5/687.6/209.5 | excluded (no published spec) |
| [test.meng2023a](case_input/test.meng2023a/README.md) | layered velocity structure | internal | 400m/4/49.7/400.1/99.0 | excluded (no published spec) |
| [test.meng2023cb](case_input/test.meng2023cb/README.md) | layered velocity, multi-patch | internal | 400m/4/50.7/346.3/113.9 | excluded (no published spec) |
| [test.tpv29](case_input/test.tpv29/README.md) | strike-slip, fractal rough fault | [TPV29](https://strike.scec.org/cvws/tpv29_30docs.html) | 500m/4/148.3/1108.5/229.4 | 50m/20s |

`test.tpv36` and `test.tpv37` are the ninth and tenth gated cases (dipping
thrust, wedge degeneration; `abs-max` gate at 1e-6, registered 2026-09-17 at
v5.10.0) and are absent from the table above only because their per-backend
fast-tier wall times have never been recorded as a sweep measurement -- an
unrecorded number is left unrecorded rather than estimated (rule 6).
`testsys/matrix.py` is the authoritative cell list: 10 cases, `CASE_BOUND` and
`GATE` per case.

# Environment
*Optional (Python solver on GPU)*: `pip install "jax[cuda12]"` — the Python solver (`python3 -m eqdyna <case_dir> --backend jax`) then runs on NVIDIA GPUs; verify with `python3 testsys/run.py gpu`.

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
python3 testsys/run.py all # unit + regression + the sweep: 10 cases x 3 backends
                           # (fortran, python-numpy, python-jax) against one
                           # canonical reference each, 30 cells.
                           # Cells run concurrently. The three most recent
                           # recorded 30-cell sweeps on this 64-core box took
                           # 2263.3 s, 2171.3 s and 2123.2 s wall clock.
                           # It prints which cells it ran, so a pass states its scope.
```
For bash, please insert the following lines in .bashrc
```
export EQDYNAROOT=/path/to/EQdynaRootDirectory
PATH=$EQDYNAROOT/bin:$EQDYNAROOT/scripts:$PATH
```

# Performance

Every number here is measured, and names the command that produced it. If you
cannot reproduce one, that is a bug in this section -- report it.

**Per-step solver cost**, `test.tpv8`, pinned to a single core, fixed
setup/compile subtracted by difference over two step counts
(`python3 testsys/run.py perf`, 2026-09-16, load 4.35):

| engine | ms/step | vs fortran |
|---|---|---|
| fortran | 305.8 | 1.00 |
| jax-cpu | 191.3 | 0.63 (faster) |
| numpy | 529.5 | 1.73 (slower) |

Gated on the PER-STEP ratio, not total wall clock: at 114 steps XLA compile is
15% of a jax run and 2% of a numpy one, so a total-time metric drifts when the
compiler changes and the solver does not. A baseline recording a different
metric is refused rather than compared against.

**Peak RSS** (`/usr/bin/time -v`, one case at a time; `testsys/matrix.py`'s
`MEASURED_PEAK_RSS_GB` is the authoritative record and the reason CI runs only
the cells that fit a 7 GB runner):

| case | numpy | jax |
|---|---|---|
| test.tpv8 | 1.45 GB | 1.91 GB |
| test.tpv10 | -- | 3.23 GB |
| test.tpv104 | 2.29 GB | 3.42 GB |
| test.tpv1053d | -- | 3.95 GB |
| test.drv.a6 | -- | 9.57 GB |

**Full sweep**: 30 cells. The three most recent recorded sweeps on this
64-core box, each a release gate run of `python3 testsys/run.py all`, took
**2263.3 s** (v5.11.0's gate, `docs/SESSION_LOG_2026-09-17_v5.10.0-followon.md:689`),
**2171.3 s** (2026-09-19, reused as v5.12.0's gate) and **2123.2 s** (v5.10.0's
gate). Wall clock varies with this box's foreground load, which other tenants
supply -- see `pathway_forward.md` item 42 -- so treat the spread, not any one
figure, as the number. The older "24 cells in ~1100 s" claim here was a
24-cell-era figure carrying no provenance and has been dropped rather than
rescaled. Cells are independent and each uses about one core, so sweep
throughput comes from running many at once, not from scaling one cell.

**Core scaling of a single jax-CPU run is NOT currently characterised.** One
ratio is measured on current code -- 193 ms/step on 1 core vs 13.8 ms/step
unpinned, 14x. The older `2.00x / 3.96x / 7.48x on 2/4/8 cores` curve predates
the solver restructure and has not been reproduced; it also has a knee at
exactly 8 cores, which is this box's NUMA node size (8 nodes x 8 cores, 2
sockets x 32), so it may be measuring memory locality rather than the solver.
`testsys/perf/run_numa_scaling.py` separates the two and REFUSES to run on a
loaded box -- a bandwidth measurement taken under other users' memory traffic
measures them. Tracked as pathway_forward item 33; do not cite a scaling
figure until it lands.


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
<!-- BEGIN EXIT CODES (generated from src/fortran/errorCodes.f90; do not edit by hand) -->

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
| 22 | `ERR_INPUT_FILE_STALE` | an input file is present but predates this binary's input format |

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
| 51 | `ERR_MPI_FAULT_ALIGNMENT` | a rank boundary in y coincides with the fault plane (not raised as of v5.8.2; downgraded to a NOTICE, see syncArnBoundary) |
| 52 | `ERR_MPI_BAD_NEIGHBOR` | point-to-point exchange with a rank outside 0..npx*npy*npz-1 |

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
`abortRun` in `src/fortran/errorCodes.f90`, which calls `MPI_Abort` so the whole job
ends instead of one rank stopping while the others block in a collective.

<!-- END EXIT CODES -->

# Collaboration
We try our best to make *```EQdyna```* easy to use but doing great science is our priority. We welcome collaborations, comments and suggestions. Please reach out to Drs Benchun Duan (bduan@tamu.edu) and Dunyu Liu (dliu@ig.utexas.edu).

# Reference
1. [Duan and Oglesby (2006)]( https://doi.org/10.1029/2005JB004138). Heterogeneous fault stresses from previous earthquakes and the effect on dynamics of parallel strike-slip faults, J. Geophys. Res., 111.
2. [Duan (2012)](https://doi.org/10.1029/2011JB009124). Dynamic rupture of the 2011 Mw 9.0 Tohoku-Oki earthquake: Roles of a possible subducting seamount, J. Geophys. Res., 117.
3. [Luo and Duan (2018)](https://doi.org/10.1029/2017JB015320). Dynamics of non-planar thrust faults governed by various friction laws, J. Geophys. Res. Solid Earth, 123.
4. [Liu and Duan (2018)](https://doi.org/10.1785/0120170374). Scenario Earthquake and Ground‐Motion Simulations in North China: Effects of Heterogeneous Fault Stress and 3D Basin Structure, Bull. Seismol. Soc. Am., 108(4).
