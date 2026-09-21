# Past release notes\

# News in 2026
* 20260921 v5.13.1 release notes
  * Fix - **v5.13.0 shipped with no version bump at all, so its tagged tree declares itself 5.12.0 at runtime.** No commit in `v5.12.0..v5.13.0` touched `VERSION` (`git log --oneline v5.12.0..v5.13.0 -- VERSION` is empty), and the runtime banner and README's News block were stale to match, so a v5.13.0 binary prints `Welcome to EQdyna 5.12.0` and misattributes every run to the previous release. The tag is annotated and pushed and is **not** being re-pointed (rule 8): v5.13.0's tree stays permanently self-inconsistent, and that is recorded here and in `pastReleaseNotes.md` rather than hidden. v5.13.1 is the next tag whose tree names its own version again. v5.13.0 also never got a GitHub Release object (`gh release view v5.13.0` -> `release not found`), so the Releases page skips from v5.12.0 to v5.13.1.
  * Fix - **`test_release_complete.py` could not see that failure, by construction -- the more important half of this release.** It keyed every check off `VERSION` and SKIPPED entirely when no tag matched it, so the only direction it tested was "VERSION ahead of the tags". The direction that failed was the mirror image: VERSION *behind* an already-released tag. At `b014178` the guard was green throughout, re-verifying v5.12.0 -- a release that was already complete -- and never looked at v5.13.0. New check `check_version_not_behind_newest_tag` compares `VERSION` against the highest semver tag reachable from HEAD (`git tag --merged HEAD`, so an unmerged maintenance tag is not a false positive) and runs **before and independently of** the skip. A dropped bump now fails the cheap regression tier at the first commit after the tag instead of surviving to the next release.
  * Docs - **README's sweep counts corrected against the code, not against other docs.** The sweep is 10 cases x 3 backends = **30 cells** (`nameList` in `testNameList.py` at the repo root, `len(matrix.CASE_BOUND)`), not the 8 x 3 = 24 README claimed; CI gates **25 of 30** (`len(matrix.CI_CELLS)`), not 19 of 24 -- `testsys/matrix.py`'s own closing comment already said 25 of 30 while README said otherwise. The e2e-ci job descriptions now match `.github/workflows/test.yml`: the two fortran groups pass **5 cases each of 10** (`test.yml:184`, `:228`, both named in full), and the cheap python group also runs `test.tpv36`/`test.tpv37` on both python backends (`test.yml:275`), which README omitted entirely. `test.tpv36`/`test.tpv37` are also now named beneath the benchmark table, which lists 8 of the 10 cases.
  * Docs - **the "~1100 s" full-sweep figure was a 24-cell-era number with no provenance (rule 6), and is replaced by measurements that name their source.** The three most recent recorded 30-cell sweeps on this box: 2263.3 s (v5.11.0's gate, `docs/SESSION_LOG_2026-09-17_v5.10.0-followon.md:689`), 2171.3 s (2026-09-19, reused as v5.12.0's gate), 2123.2 s (v5.10.0's gate). The spread, not any one figure, is the number -- wall clock on this shared box tracks other tenants' load (`pathway_forward.md` item 42). Two unmeasured timings (`test.tpv36`/`test.tpv37` per-backend fast-tier seconds) are recorded as unmeasured rather than estimated.
  * Note - **v5.13.0's own contents are written up retroactively in `pastReleaseNotes.md`**, since that tag shipped without notes: refactor rounds 4, 5 and 6 (round 5 a provable no-op -- AST-identical after docstring strip; round 6 byte-identical across all three batch-script writers); `run_scaling.py`'s Fortran metric-bias fix; the tpv36/tpv37 symlink-coupling regression guard plus `PROJECT_RULES.md` rule 18; and the `publish.yml` fix for v5.12.0's Docker image never publishing. No physics changed and no reference was regenerated in v5.13.0 or in v5.13.1.





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
* 20260919 v5.13.0 release notes -- WRITTEN RETROACTIVELY AT v5.13.1 (tag dated 2026-09-19)
  * This tag shipped with no release notes and no version bump. `VERSION`, the Fortran runtime banner (`src/fortran/eqdyna3d.f90:17`) and README's News block all still read 5.12.0 on the v5.13.0 tagged tree (`git log --oneline v5.12.0..v5.13.0 -- VERSION` is empty), so a v5.13.0 binary prints `Welcome to EQdyna 5.12.0`. The tag is annotated and pushed; it is NOT being re-pointed (rule 8), so that tree stays permanently self-inconsistent and this entry is the record of it. There is also no GitHub Release object for it (`gh release view v5.13.0` -> `release not found`), so the Releases page skips from v5.12.0 to v5.13.1.
  * What it actually contained: refactor rounds 4, 5 and 6 (round 5 a provable no-op -- AST-identical after docstring strip; round 6 byte-identical output across all three batch-script writers); `testsys/perf/run_scaling.py`'s Fortran metric-bias fix; the tpv36/tpv37 symlink-coupling regression guard plus `PROJECT_RULES.md` rule 18; the `publish.yml` fix for v5.12.0's Docker image never having published; and a start on the 8-case/24-cell doc-drift correction, which did not reach README's own copies of those counts -- those are fixed in v5.13.1. No physics changed and no reference was regenerated.
* 20260919 v5.12.0 release notes
  * Fix - **three latent `ntotft>=2` (multi-fault) bugs in the Fortran solver, guarded by new regression tests** (pathway items 7/9/10). Confirmed true no-ops for the `ntotft==1` case every currently-gated benchmark runs. Still not exercisable end-to-end: multi-fault input is refused loudly at `case.setup` (pathway item 17, deferred by owner decision), so these fixes have no gated case to prove them end-to-end yet.
  * Fix - **CI ran the same commit twice** on any push to a branch with an open PR (`push` and `pull_request` firing independently). Grouped the workflow's concurrency key by head commit so the second run cancels the first instead of both completing (pathway item 35).
  * Change - **JAX/NumPy core-scaling tooling (item 33): NUMA placement fixed twice more.** `run_scaling.py` now probes per-cpu occupancy fresh per configuration and places work on genuinely free NUMA nodes instead of always starting at node 0, and MPI ranks are pinned via a per-rank `numactl` rankfile rather than relying on `mpirun --bind-to` alone. Cleanest measurement to date: JAX compact placement peaks at 2.52x (16 cores), 1.98x (32 cores); NumPy remains flat (no multi-core benefit) across the same range. The standing "why does JAX plateau past 8 cores" question (scatter-add vs memory-bandwidth) stayed genuinely open after a follow-up microbenchmark's own re-verification found its first result did not reproduce across repeat runs on this shared box -- recorded as NOT SETTLED rather than closed on a single pass.
  * Refactor - **three-round duplication-trim pass across code comments and session-log docs, no behaviour change.** Round 1 trimmed comment/prose bloat in source; round 2 audited every regression-test guard for a deletable duplicate and found none (each ties to a dated incident); round 3 consolidated restated detail across three older session logs into pointers at `pathway_forward.md`'s own item rows. Net -356 lines across the three rounds (18 files touched total). Each round gated on a fresh full 30/30 `python3 testsys/run.py all` before landing.
  * Known issue - **`test.tpv30` is not gated: a fault-local equilibrium defect in the `C_elastic=0` (viscoplastic) path, root cause still under investigation.** Its initial on-fault shear at `faultst000dp120`, dx=500 m, reads 33.05 MPa and relaxes 4.3 MPa toward 27.79 MPa over the first ~2 s of simulated time. The owner's own three published TPV30 submissions (100/50/25 m) all read a stable, resolution-independent ~27.8-27.9 MPa with NO comparable relaxation -- ruling out a coarse-mesh discretization artifact (which would grow with element size) as the explanation. See `pathway_forward.md` items 39/19(b) for the full, still-evolving investigation record.
  * Known issue - **grid-scale on-fault stress ringing at 500 m is still not damped** -- the one alternative hourglass control this code ships (`C_hg=2`) over-damps the rupture itself instead of selectively suppressing the mesh-scale mode (unchanged from v5.11.1). See `pathway_forward.md` item 34.
  * CORRECTION - **the Docker image did NOT actually publish for v5.12.0.** The tag's publish run (`publish.yml`) failed its own pre-push gate: `test_stress_i0_carry_aliasing.py` (item 41's guard) exercises both backends, and the image deliberately ships without jax (`Dockerfile:25-27`) -- neither change was wrong on its own, they simply had not met until this tag. `test.yml`'s main test tier is unaffected and stayed green throughout; only the container did not ship. Fix (`9d1977d`: add jax to `publish.yml`'s CI-only gate overlay) verified via `workflow_dispatch` (build+gate, no push) without burning a second tag; it rides the next tag. Same pattern as v5.8.3/v5.8.4 -- a release whose code is fine but whose image needs a follow-up tag to actually publish.
* 20260917 v5.11.1 release notes
  * Fix - **latent aliasing bug in the Python port's time-stepping carry state (`driver.py:182`).** Under the NumPy backend, `xp.asarray(inv['stress_i0'])` returned the SAME array object rather than a copy, aliasing the per-step carry slot to the dict entry built once at solver-state construction. Harmless today (nothing else re-read that value), but latent against the carry's own deliberate in-place mutation every step. Fixed with an explicit `.copy()`, verified as a byte-identical no-op by a fresh full sweep: unit 249 SUCCESS, regression SUCCESS, e2e **30/30 cells SUCCESS** (2275.0s), all max|diff| figures unchanged from v5.11.0. New regression guard: `testsys/regression/test_stress_i0_carry_aliasing.py`.
  * Known issue - **`test.tpv30` is not gated: a fault-local equilibrium defect in the `C_elastic=0` (viscoplastic) path, root cause still under investigation.** Its initial on-fault shear at `faultst000dp120`, dx=500 m, reads 33.05 MPa and relaxes 4.3 MPa toward 27.79 MPa over the first ~2 s of simulated time. The owner's own three published TPV30 submissions (100/50/25 m) all read a stable, resolution-independent ~27.8-27.9 MPa with NO comparable relaxation -- ruling out a coarse-mesh discretization artifact (which would grow with element size) as the explanation. See `pathway_forward.md` items 39/19(b) for the full, still-evolving investigation record.
  * Known issue - **grid-scale on-fault stress ringing at 500 m is not damped by the one alternative hourglass control this code ships (`C_hg=2`)** -- it was measured this release and found to over-damp the rupture itself (some stations arrest entirely) rather than selectively suppress the mesh-scale mode. See `pathway_forward.md` item 34.
* 20260917 v5.11.0 release notes
  * New - **CI coverage widened 10 -> 25 of 30 cells.** `test.tpv36`/`test.tpv37`'s python-numpy and python-jax columns are now gated in CI (previously fortran-only), admitted on measured peak RSS (2.91-4.34 GB against the 7 GB runner) -- closing exactly the kind of backend-coverage gap that let a real wedge-degeneration bug through last release.
  * Change - **JAX/NumPy core-scaling tooling fixed and measured (item 33).** `testsys/perf/run_scaling.py` now actually tests up to 32 cores (was capped at 8) with NUMA-aware placement (`numactl`, replacing a bare `taskset` that could span multiple memory nodes). JAX measured at 2.35x speedup on 32 cores. NumPy does not meaningfully scale with more cores -- root-caused to its BLAS backend deciding its thread count once, at import time -- and a landed fix bounds the worst-case regression from bad core placement without making NumPy multi-core (a real fix needs a kernel rewrite, out of scope here).
  * Change - CI's python-numpy cells are now thread-pinned, removing a placement-lottery effect that previously made the same cell's wall-clock vary several-fold run to run.
  * Known issue - **`test.tpv30` is not gated: a fault-local equilibrium defect in the `C_elastic=0` (viscoplastic) path, root cause still under investigation.** Its initial on-fault shear at `faultst000dp120`, dx=500 m, reads 33.05 MPa and relaxes 4.3 MPa toward 27.79 MPa over the first ~2 s of simulated time. The owner's own three published TPV30 submissions (100/50/25 m) all read a stable, resolution-independent ~27.8-27.9 MPa with NO comparable relaxation -- ruling out a coarse-mesh discretization artifact (which would grow with element size) as the explanation. Everything upstream of the fault nodes is independently verified correct (the stress tensor to 1e-6, the mesh morph to 5e-4 m, the FE force assembly itself, all shared geometric/inertial quantities bit-identical or exact to the elastic case at all 1921 common fault nodes) -- the defect is isolated to a force imbalance specifically at the fault nodes in this code path. Two architectural fixes were investigated and both retracted before landing once this resolution-independence evidence came in; neither is the answer. `test.drv.a6` shares the same `C_elastic=0` path and its reference will very likely need regeneration once this is fixed (a separate, reviewed commit, rule 7). See `pathway_forward.md` items 39/19(b) for the full, still-evolving investigation record.
* 20260917 v5.10.0 release notes
  * New - **TPV36 and TPV37 are both now gated on all three backends** (fortran, python-numpy, python-jax) -- the concrete payoff of the "all backends or it does not land" policy. TPV37 needed no new port work (both `faulting.f90` and `faulting.py` already carried its branch from TPV36's own landing) and lands bit-exact on fortran, 6.929140e-09 / 5.256410e-09 on numpy/jax against a 1e-6 bound. Worth restating why this policy matters: gating TPV36 on every backend is what exposed a real silent-wrong-physics bug -- the Python port's wedge-degeneration interior-force dispatch excluded elemType>10, so wedge elements carried mass and hourglass resistance but contributed ZERO interior force. A fortran-only gate would never have caught it; the all-backends gate did (see `pathway_forward.md` item 19(a)).
  * Change - **`bGlobal.txt` input-file contract change.** Four previously-hardcoded viscoplastic values are now real case inputs (item 24 b/c/e/f): the nucleation relaxation time Tv, the off-fault deviatoric depth taper, the case's own shear modulus, and the plastic-strain output window -- each defaulting to exactly the constant the solver used to compile in, so an unconfigured case is unchanged. `case.setup` now writes a new trailing block in `bGlobal.txt` unconditionally, and the solver refuses a `bGlobal.txt` written by an older `case.setup` loudly (new error code `ERR_INPUT_FILE_STALE=22`), rather than silently defaulting. If you have your own `user_defined_params.py` from before this release, re-run `case.setup` before your next run.
* 20260916 v5.9.0 release notes
  * New - **`test.tpv36` gated on all three backends** (fortran bit-exact, python-numpy 7.894064e-09, python-jax 5.867293e-09, bound 1e-6): the first TPV benchmark added since the all-backends policy, and it lands WITH its Python columns rather than promising them later. The Python port gains C_degen wedge-degeneration DYNAMICS (mesh support shipped earlier this cycle) -- `calcGlobalShapeFunc.f90:22-28`'s shape-function merge and `calcElemKU.f90`'s wedge-element interior-force dispatch, both verified directly against the real, unmodified Fortran kernels by a new jax-aware parity probe (`testsys/parity/evidence_wedge_kernel.py`).
  * Fix - a real regression caught before it shipped: the first attempt at this port passed every check (including a from-scratch kernel probe) but broke `test.tpv8 x python-jax` 10x over its bound in the full sweep -- reverted within the hour, root-caused to a genuine NumPy floating-point association difference (`np.einsum`'s BLAS-batched-matmul dispatch vs its per-element form, ULP-level, on a 3-term contraction) that only the nonlinear rupture dynamics amplified past the tightest bound in the suite. Fixed by restoring the exact pre-port einsum form for every non-wedge element; wedge elements alone take the new per-element path.
  * New - CI widened from 10 to 20 of 27 cells: matrixing the workflow into parallel jobs (v5.8.6) removed the shared-runner memory ceiling that had excluded `tpv10`/`tpv104`/`tpv1053d` (jax), `meng2023a`/`meng2023cb`/`tpv29` (both backends), all newly measured and added to `matrix.MEASURED_PEAK_RSS_GB`. `test.drv.a6 x python-jax` (9.57 GB) is the one still-excluded cell.
* 20260916 v5.8.7 release notes
  * Fix - v5.8.6's `verify-published-image` CI job reported a false FAIL (a failed git network round-trip inside a container run on a separate job/runner was treated as a confirmed lightweight tag). Fixed to report UNVERIFIED instead, matching an existing pattern elsewhere in the same file. v5.8.6's image had, in fact, already published successfully.
* 20260916 v5.8.6 release notes
  * Fix - the Docker image had never actually built the solver: `ubuntu:22.04`'s `mpich` package depends on `libgfortran5` (the runtime library) but not `gfortran` (the compiler `mpif90` shells out to for every build), so `make` failed with `gfortran: command not found` inside the Dockerfile's own build step. `gfortran` added explicitly to both `Dockerfile` and `install-eqdyna.sh`'s ubuntu branch, which had the identical latent gap. A second issue found right behind it: `--no-install-recommends` gave the image a STRICTER dependency closure than `install-eqdyna.sh -m ubuntu` itself uses, silently dropping `mpif.h`; the flag is now dropped.
  * New - `.github/workflows/publish.yml` gained `workflow_dispatch` (debug iterations no longer cost a release); CI parallelized 660s -> 303s (-54%, coverage unchanged); two release-guard bugs fixed (annotated-tag remote-peel fallback, `GH_TOKEN` for the Release check).
  * CORRECTION: the image DID actually publish for this tag (`docker push` succeeded) -- what failed was a false FAIL in the SEPARATE verify-published-image job, itself a bug in the annotated-tag fallback added above (it treated a failed network round-trip the same as a confirmed lightweight tag). Fixed in v5.8.7.
* 20260916 v5.8.5 release notes
  * Fix - v5.8.4's Dockerfile fix addressed the failing log line, not the pattern: `.github/workflows/publish.yml` had two more copies of the same `--break-system-packages` flag (the CI-only pytest/scipy overlay in both the `build-and-push` and `verify-published-image` jobs). Both fixed; a repo-wide grep confirmed no more occurrences (`install-eqdyna.sh`'s own use of the flag is in the macOS/Homebrew branch, which does need it -- checked, not a matching bug). CORRECTION: the image still did not publish for this tag either -- a missing `gfortran` dependency (found next) blocked the Docker build itself one layer further in; fixed in v5.8.6.
* 20260916 v5.8.4 release notes
  * Fix - the Dockerfile added in v5.8.3 could not actually build: `ubuntu:22.04`'s stock pip3 (22.0.2) does not support the `--break-system-packages` flag (added in pip 23.0.1), so `RUN pip3 install --break-system-packages ...` failed with exit 2 on the very first real build v5.8.3's tag push triggered. Flag removed (22.04's system Python has no PEP 668 externally-managed-environment marker, so none is needed). CORRECTION: the image still did not publish for v5.8.4 either -- `publish.yml` had two more copies of the identical flag one step further down the pipeline; fixed in v5.8.5, the first tag to actually publish an image.
* 20260916 v5.8.3 release notes
  * Fix - item 23 (Drucker-Prager viscoplastic kernel, `calcElemKU.f90:127-161`, friclaw 4) had zero isolated test; `testsys/regression/test_drucker_prager_kernel.py` now compiles the real, unmodified Fortran routine and diffs it against the Python port on 4 synthetic stress states, no MPI/case/mesh required.
  * Fix - item 28's TPV29 cross-code comparison is now a committed, reproducible script (`testsys/parity/evidence_tpv29_scec_comparison.py`), not prose. Running it found a real error in the previously recorded claim: Mw 7.45 has no reproducible provenance and is physically impossible for this run (Kanamori relation puts it at 4.7x the run's actual seismic moment); corrected to 7.034, derivation left inline.
  * New - `src/python/eqdyna/meshgen.py` ports `C_degen` wedge-degeneration MESH generation (the fault-adjacent hex-to-2-wedge split), verified an exact element/node-count match against a fresh Fortran build on test.tpv36. Its DYNAMICS remain explicitly refused (not silently wrong) -- test.tpv36/test.tpv37 stay ungated pending that further port.
  * New - `Dockerfile` + `.github/workflows/publish.yml` added, gating a tagged release's image on the repo's own test tiers before push. CORRECTION: the image did NOT actually publish for v5.8.3 itself -- `ubuntu:22.04`'s stock pip3 does not support a flag the Dockerfile used; the gate this release added caught its own Dockerfile's latent bug on the very first real build. Fixed in v5.8.4, which is the first tag to actually publish an image.
  * New - `scripts/scec/` tracks a sha256 manifest (`CHECKSUMS.sha256`) of the 693-file SCEC archive plus an `organize.py --verify` mode, so the archive's integrity is checkable without re-fetching 484 MB.
  * Housekeeping - `git gc` reclaimed loose-object bloat (307 MiB loose -> 1.85 MiB loose; 283.57 MiB stays packed, since the 693 blobs a prior `.gitignore` bug tracked are still reachable through history -- a real fix needs an owner-approved history rewrite, not done here).
* 20260916 v5.8.2 release notes
  * Fix - test.tpv36/test.tpv37 (dipping, wedge-degenerate faults) run under MPI again. `checkFaultMPIAlignment` was correctly refusing a y-split of the fault (exit 51); the `MPI_ERR_RANK` first reported was the MPI library's secondary noise from partners that had already `MPI_Abort`'ed, not the real cause. The divide-vs-duplicate `arn` accounting for this branch is now audited (same case with/without a y-split, 3416 fault nodes: `tnrm`/`tstk`/`tdip` ratio 1.000000, max |diff| over all 22 canonical columns 1.0e-08) and was already correct, so the refusal is now a rank-0 NOTICE carrying that evidence, pinned by new regression guard `test_dipping_fault_y_split.py`.
  * Fix - a data-integrity bug found auditing this release: a trailing inline comment on `.gitignore`'s `scec_archive/` line made the whole line, comment text included, the literal (unmatchable) pattern, so `scec_archive/` was never actually ignored. 693 files / ~1.16 GB of fetched SCEC submissions were tracked into git by v5.8.2's own predecessor commit despite every doc in the repo saying that data stays untracked. Pattern fixed and the files untracked going forward in this release; the blobs remain in git history pending a deliberate, separately-approved history rewrite.
  * Fix - `requireValidNeighbor` (new error code 52, `ERR_MPI_BAD_NEIGHBOR`) on both `mpi_sendrecv` call sites, so an out-of-range MPI rank names its call site, direction, side and decomposition instead of leaving the reader to guess.
  * Fix - corrected in-code explanation of what discriminates the DUPLICATE vs DIVIDE fault-boundary branches: the NOMINAL GRID y-extent in `fltxyz`, not the physical dip angle -- a 60-degree inserted fault (test.tpv10) behaves like a vertical one. Also corrected a docstring claim that the DIVIDE case was "not producible by any case in this codebase today" (that survey covered only the 8 gated cases).
  * New - interim `C_degen != 0` refusal in the Python port (`src/python/eqdyna/eqdyna3d.py`): the port was documented as `C_degen==0` scope in four docstrings with nothing enforcing it; it now refuses rather than silently mis-treating a degenerate wedge as a hex and the fault as a y=0 plane. test.tpv36/test.tpv37 remain UNGATED (not added to `testNameList.py` or `testsys/matrix.py`) until wedge degeneration is actually ported.
  * Fix - moved the SCEC re-fetch TOOLING (`scripts/scec/`: fetch_all_dliu.py, organize.py, make_index.py, provenance.py, pure stdlib) out of gitignored `scratch/` and into the tracked tree, so TPV29's cross-code validation baseline is regenerable from a clean clone even though the 484 MB it fetches correctly stays untracked.
  * Fix - `testsys/regression/test_release_complete.py`: a released version must be released everywhere. v5.8.0 and v5.8.1 both skipped the `pathway_forward.md` Tasks-done row and `gh release create`; this guard checks the runtime banner, README's leading News block, a dated Tasks-done row, an annotated tag, and (best-effort, network) the pushed tag and GitHub Release, for whatever version `VERSION` currently names.
  * Change - `PROJECT_RULES.md` rule 17 step 7 and `CLAUDE.md`: supporting a TPV means supporting it on every backend -- a case is never added with its Python columns declared UNSUPPORTED to be filled in later.
* 20260916 v5.8.1 release notes
  * Fix - per-stage timers shared ONE global `startTimeStamp`, written by eight sites across six files, so any nesting silently truncated the outer measurement. `compTimeInSeconds(2)` was understated **511x** (0.0001 s where the true cost is 0.0511 s), not the "<=1% bias, harmless" the issue recorded. Each site now uses a local. Output bit-identical.
  * Fix - `case.setup` REFUSES `par.ntotft > 1` with a named reason. It previously produced a `bStations.txt` whose line 2 carried one on-fault station count where `readInputFiles.f90:152` reads `ntotft` of them, and the run died in list-directed input with `Bad integer for item 2` and exit 2. That crash also made the deliberate multi-fault stop at `meshgen.f90:790` unreachable.
  * Fix - `test_stop_exit_status.py`'s real-binary probe -- the only check that a refused run EXITS rather than hangs under mpirun -- had been silently SKIPPING on every ordinary tree, because it looked for `src/eqdyna` while `install-eqdyna.sh` moves the binary to `bin/`. Live again, and running in CI for the first time.
  * Fix - `test_ci_dependencies.py` never scanned `src/python`, the solver package, despite its own docstring claiming every third-party import is checked against CI's pip line.
  * Fix - `testsys/perf/run_scaling.py` set PYTHONPATH to a directory that does not exist.
  * Fix - README's performance numbers were stale in three ways: embedded v5.7.1 release notes quoting 10.4 GB for tpv104 (now 3.42), a case table reporting `~48s` for cases that range 13.3-1108.5 s, and no statement of what produced any figure. There is now a Performance section where every number names its command, and it says plainly that single-run core scaling is NOT characterised.
  * Fix - PROJECT_RULES rule 16 told the reader to reproduce CI's entry point and then quoted a command CI does not run. Six of seventeen rules cited deleted files.
  * New - rule 17 and `CLAUDE.md`: the TPV revival workflow, written from what TPV29 cost. The load-bearing step is checking the CODE has that TPV's branch before trusting a run.
  * Closed by measurement, not by fixing - item 24(a)'s "tractions are exactly HALF for C_elastic==0" does NOT reproduce (Tn ratio 1.0018 over 4747 nodes); item 11's counting-vs-allocating divergence is structurally impossible AND already aborts; item 32's additive flip model is refuted -- drv.a6's flips are a fixed marginal population of ~400 roundoff-decided nodes, not a sum of error sources.
* 20260916 v5.8.0 release notes
  * New - `src/fortran/` + `src/python/`: the two implementations of the same solver side by side. 26 .f90, and 14 .py of which every one but `backend.py` (the numpy/jax adapter) and `__main__.py` maps 1:1 onto a Fortran file.
  * New - friclaw 2 (time-weakening) ported: test.meng2023a and test.meng2023cb run on the Python backends for the first time. The sweep is 24 cells with ZERO declared-unsupported.
  * Fix - TPV 29 was missing from `swtwNucleation`'s smoothed forced-rupture list, so test.tpv29 had been declaring `par.tpv = 36` to reach its own spec formula. That hid a port gap which nucleated 0 of 3321 fault nodes against a reference of 2974: 7.10e7 -> 9.9e-15.
  * Fix - the port ended each step with `force * inv_mass` where `driver.f90:30` DIVIDES. friclaw 4/5 moved toward the reference; drv.a6 x python-numpy 461 -> 423 rupture-arrival flips under its UNTOUCHED 450 budget.
  * Fix - the parallel sweep could deadlock (multi-core cells reserved units one at a time). Observed: 80 minutes, 14 of 24 cells unrun, no error. Guarded.
  * Change - ONE test: the `accept` and `parity` tiers are gone, along with the pydump Fortran instrumentation they needed. ~5,700 lines net removed.
  * Perf - solver 1,497 -> 842 lines; tpv104 numpy peak RSS 2.52 -> 2.29 GB. The perf tier now gates on PER-STEP cost with compile subtracted, after a total-wall-clock metric reported a 33% JAX regression that did not exist.
* 20260915 v5.7.1 release notes
  * Fix - v5.7.0's CI went red: the Python/JAX backend inside e2e needs 10.4 GB for test.tpv104 and a GitHub runner has 7 GB, so the job was SIGTERM'd (exit 143) with no output. CI now gates test.tpv8 (measured 3.59 GB at the time; 1.91 GB jax / 1.45 GB numpy as of v5.8.0), and e2e prints its case list so a green check states its own coverage instead of implying five.
  * Fix - the e2e python child runs unbuffered; under the previous buffering a resource kill produced zero output, so the log showed seven minutes of silence and a bare exit code.
  * Change - release workflow (rule 15): the tag and GitHub Release now come AFTER CI is green on the pushed commit, not before. v5.7.0 was tagged on a local green that could not model the runner's memory; a released tag pointing at a red commit makes the Releases page the authoritative wrong answer.
  * Note - all five accept cases run and pass locally. The exclusion is the runner's memory, not the cases. Python peak RSS on test.tpv104 is 10.4 GB (jax) / 4.2 GB (numpy) against Fortran's 0.95 GB; reducing that is tracked in pathway_forward.md. [SUPERSEDED: the jit closure fix in v5.8.0 took tpv104 jax to 3.42 GB and the solver restructure took numpy to 2.29 GB. See README's Performance section for current figures.]
* 20260916 v5.8.0 release notes
  * New - `src/fortran/` and `src/python/`: the two implementations of the same solver now sit side by side under one parent. 26 .f90 + makefile, and 14 .py of which every one but `backend.py` (the numpy/jax adapter) and `__main__.py` maps 1:1 onto a Fortran file.
  * New - the Python solver mirrors the Fortran design. The three per-friclaw solver copies (port/port_rsf/port_tp, each duplicated again for JAX) collapse into ONE driver whose friclaw dispatch sits inside `solveSWTW`, exactly as `faulting.f90:17-18` does it. Solver 1,497 -> 842 lines; tpv104 peak RSS 2.52 -> 2.29 GB.
  * New - friclaw 2 (time-weakening) ported. `test.meng2023a` and `test.meng2023cb` run on the Python backends for the first time (4.63e-09 / 3.15e-09 / 5.05e-09 / 4.40e-09 against a 1e-3 bound). The sweep is 24 cells with ZERO declared-unsupported.
  * Fix - TPV 29 was missing from `swtwNucleation`'s smoothed forced-rupture list `{201, 36, 37}` in the Fortran, so `case_input/test.tpv29` had been declaring `par.tpv = 36` to reach its own spec formula. That misdeclaration hid a real gap: the port implemented only the degenerate tr=1e9 branch, nucleated 0 of 3321 fault nodes against a reference of 2974, and failed at max|diff| 7.10e7 on BOTH backends at the identical value. Now 2974 of 2974, zero flips, 9.9e-15.
  * Fix - the port deviated from `driver.f90:30`. `port_rsf.py` and `port_tp.py` ended each step with `force * inv_mass`; Fortran DIVIDES, and a reciprocal-multiply is not the same rounding. The unified driver divides for every friclaw, so friclaw 1 is bit-identical to the old port and friclaw 4/5 move TOWARD the reference: `test.drv.a6 x python-numpy` 461 -> 423 flips against its UNTOUCHED 450 budget, `test.tpv104` still green at 4.15e-07 against 1e-5. The case that had been failing was never a tolerance problem.
  * Fix - the parallel sweep could DEADLOCK. Multi-core cells were reserved with `sem.acquire()` in a loop, taking units one at a time, so two 4-rank cells could each hold 3 of a 6-core budget and wait forever. Observed: 80 minutes, 14 of 24 cells unrun, parent at 0.1% CPU, no error and no timeout. Replaced with an all-or-nothing Condition allocator; guarded by `testsys/regression/test_sweep_core_budget.py`, which was verified to reproduce the deadlock before being trusted (a first version passed on the broken code).
  * Fix - `test_stop_exit_status.py`'s real-binary probe -- the only check that a refused run EXITS rather than hangs under mpirun -- looked for `src/eqdyna` while `install-eqdyna.sh` moves the binary to `bin/`. It had been printing SKIPPED on every ordinary tree. Live again, and it runs in CI for the first time: `mpirun -np 2` on an empty case exits 21, no hang.
  * Fix - `test_ci_dependencies.py` never scanned `src/python`, the solver package, despite its own docstring claiming every third-party import is checked against CI's pip line. Now scans all 14 modules.
  * Fix - `testsys/perf/run_scaling.py` set PYTHONPATH to a directory that does not exist; the tier would import-fail on invocation.
  * Fix - PROJECT_RULES rule 16 told the reader to reproduce CI's entry point verbatim and then quoted `testsys/run.py all`, while CI runs `unit regression e2e-ci` (10 of 24 cells). The rule was wrong in exactly the way it exists to prevent. Six of sixteen rules cited deleted files and now name what enforces them today; rule 8 was entirely about `testAll.py`'s unconditional `rm -rf test`, and the behaviour is now the opposite.
  * Change - references are one `frt.canonical.txt` per case instead of per-rank `frt.txt*`. A result is now a statement about the physics rather than about the decomposition: dedupe on rounded (x,y,z), then lexsort, so a 4-rank Fortran run, a serial Fortran run and a serial Python run all compare to the same artifact at the same bound. Verified lossless before the per-rank files were replaced -- 358 duplicate groups across eight cases, agreeing to 0.000e+00 in every one of the 22 columns -- and `canonicalize()` now REFUSES to dedupe rows that disagree.
  * Change - ONE test. The `accept` and `parity` tiers are gone: `accept` was the same solver against the same references with a second comparison implementation over a shorter case list, reporting "5/5 cases" for 5 of 16 cells; `parity` answered a debugging question, went stale when the Python tree was restructured, and needed the pydump Fortran instrumentation to exist. Removed with them: `src/pydump.f90`, `src/pydump_noop.f90`, the PYDUMP=1 makefile seam and both call sites, `check.test.py`, `testAll.py`, an orphaned perf elem-scaling harness nothing invoked, and `testsys/parity/fixtures/` -- frozen copies of `scripts/` that had drifted to a 187-line `lib.py` against the live 949 and were kept lint-clean by an explicit carve-out. ~5,700 lines net removed.
  * Perf - the perf tier gates on PER-STEP cost with the fixed setup/compile subtracted, by difference over two step counts. It had gated on total wall clock and reported a 1.329x degradation that did not exist: compile is 15% of a 114-step JAX run and 2% of a NumPy one, so the mixture moved for JAX while NumPy/Fortran stayed identical to four digits. Per-step jax/fortran is 0.626 against the old baseline's 0.611 -- flat. Fresh figures, pinned core: fortran 305.8 ms/step, jax-cpu 191.3, numpy 529.5. A baseline recording a different metric is now REFUSED rather than compared against.
  * Note - `test.drv.a6` at the reference's own 4-rank decomposition is BIT-EXACT (0 flips). The sensitivity is entirely decomposition-driven: a serial Fortran run against that same 4-rank reference still gives 329 arrival flips of 5151, so the Python columns' 415-423 is roughly 329 inherited plus ~90 of actual port difference. The diagnostic experiment is a serial Fortran reference, not a finer mesh.
  * Known gap - `solveRSF` divides by `v_trial` unguarded, so a node with zero background creep and zero slip rate gives 0/0 -> NaN and Newton-Raphson diverges silently. Present IDENTICALLY in the Fortran (`solveRSF:257-258`): a latent bug in both languages, not a port divergence, and unreachable by any gated case. Logged, not fixed in this release.
  * Known gap - friclaw 3 (RSF ageing law) is not covered by any gated case. friclaw 1, 2, 4 and 5 are.
* 20260915 v5.7.0 release notes
  * Perf - Python setup 14.0 s -> 1.46 s (9.6x), byte-identical: mass assembly, element build, equation numbering and coordinate/fault builders vectorized with order-preserving reductions, each keeping its scalar loop as a named oracle.
  * Perf - JAX solver 2.1x end-to-end: the jit was rebuilt inside run(), so every call paid a full ~14 s XLA compile (GPU utilisation measured at 1.6%). Now a traced trip count plus a persistent compilation cache; hourglass modes pre-summed before the scatter (-16% exec).
  * Perf - NumPy solver 1.90x (1042 -> 550 ms/step), bit-identical: duplicate contractions eliminated, loop-invariant products hoisted, the two 300 MB/step concatenates removed.
  * New - e2e is backend-parameterized: the same gated cases run through Fortran and through the Python/JAX standalone, against the same references. jax added to CI.
  * Note - the JAX GPU path is nondeterministic run-to-run (XLA lowers the duplicate-index scatter-add to atomics); measured 8.9e-08 between identical runs of unmodified code. jax-cpu is exactly reproducible. Bit-identity is not available on GPU and no gate should assume it.
  * Known gap - the accept tier has never exercised the NumPy backend; it runs the default (jax). Forcing --backend numpy gives 4/5, drv.a6 at 461 flips against a 450 bound. Pre-existing, unresolved.
* 20260915 v5.6.2 release notes
  * Fix - scipy was used by scripts/convertFaultGeometry but never added to CI's pip install, so the unit tier passed locally and failed in CI. scipy added; its absence now raises a message naming it and pointing at the exact-decimation path, which needs no scipy.
  * New - testsys/regression/test_ci_dependencies.py: every third-party import under scripts/ and testsys/ must be installed in CI, or be lazy AND guarded. Reads the package list out of the workflow file so the two cannot drift. It immediately found a second instance (imageio), now guarded.
  * Change - python/eqdyna/standalone: `--backend jax` no longer falls back to numpy, and a new `--device auto|cpu|gpu` pins the JAX platform. Neither falls back (rule 2): a run reported as jax-on-gpu must have been that. Each run now prints the device it actually used.
* 20260915 v5.6.1 release notes
  * Fix - a bare `make` built nothing, so `install-eqdyna.sh` produced no binary and every tier needing one failed. v5.6.0 fixed the `clean` target, which had accidentally been the thing making bare `make` build (`clean:eqdyna`, no space, made clean depend on eqdyna). Default goal is now pinned explicitly.
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