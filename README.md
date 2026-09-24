# News in 2026
* 20260923 v5.16.2 release notes
  * Fix - **a live `GITHUB_TOKEN` was being baked into every image `publish.yml` published.** `actions/checkout`'s `persist-credentials: true` default writes the token into `.git/config`; its cleanup is a POST step that runs only after `docker build` has already run `COPY . /opt/eqdyna` (`.dockerignore` deliberately keeps `.git` for two history-reading guards), so a real, if short-lived, credential landed in a layer of every published image. Fixed with `persist-credentials: false` plus `testsys/check_git_config_no_credential.py` run as the first line of both in-image gate blocks -- it refuses (exit 2) rather than passes when there is nothing to scan. The mechanism predates this release (unchanged since checkout v2.0.0); this is the release that closes it. **Not remediated: every image already published, `v5.16.1` included, still carries a now-expired token in a layer** -- deleting and re-pushing those is an outward-facing publish decision held open as pathway item 86, not attempted here.
  * New - **`publish.yml`'s release-path premises are now pinned, not assumed.** `testsys/regression/test_publish_image_fetch_depth.py` grew from 5 to 13 checks: `packages: write` on the push job, the `if: github.event_name == 'push'` gate on both the push step and the verify job, the two in-image gate blocks pinned as EQUIVALENT to each other (not merely each present), and the shell-level `if` gating the `:latest` tag (no YAML `if:` covers it). `actions/checkout` bumped `@v2` -> `@v4` (`@v2` is node12, unsupported since April 2022) with `fetch-depth: 0` still explicit.
  * Fix - **the station-file header column count now follows the friclaw branch instead of a fixed number.** 8 columns for the common case, 11 when friclaw's extra fields apply; new guard `testsys/regression/test_station_header_column_count.py` (item 67; `test/test.tpv8/faultst000dp075.txt` now reads 8 columns / 8 names / 8 data fields, was 11/8/8 before).
  * New - **`runlock` guards `probe_plastic_traction.build_case()`** against the same concurrent-run collision rule 21a already closed for the e2e harness (item 81); guard checks 17 -> 18.
  * Note - **`v5.16.1`'s missing GitHub Release object was backfilled and the regression tier re-verified green on master** (item 83; rules 15c and 3a added for the tag-push/Release-create ordering that caused it). Rule 14a added (a board row's evidence command must be able to read both ways); rule 15d added (the board-only commit and the VERSION+README commit are two commits, not one, per rule 21c's hook); a stale hardcoded rule count in `CLAUDE.md` corrected to point at the book's own index instead of a number.
  * Note - patch release: no e2e SWEEP was run, one cell only (`test.tpv8` x fortran, `max|diff|=3.051760e-11` against `1.0e-08`, 1891 fault nodes -- item 67 touched `src/fortran` only). No perf claim. No reference regenerated.
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
every surface value is an official one, no interpolation -- so a clean
checkout needs no downloads. The fault-normal DERIVATIVE columns that travel
with that surface are recomputed at the target spacing rather than
decimated, because a finite difference belongs to a grid at a spacing and
not to the surface (see `case_input/test.tpv29/README.md`). test.tpv29 ships its surface at 100 m (3.5 MB) and 50 m (14 MB);
a request for a dx finer than the shipped source, or one that is not a
multiple of it, is refused rather than interpolated (`scripts/lib.py`'s
`requireFaultGeometryResolution`). 25 m needs the official SCEC file and
`tpv29GeometryTools.convertOfficial25m`.

The fast tier is a single sweep over `case x backend`: each case below is run
on the Fortran solver (MPI) and on the standalone Python solver (NumPy and
JAX), and every cell is compared against the same committed reference for that
case, at that case's one tolerance, at the ONE gate term
(`testsys/matrix.py`'s `GATE_TERM_S`, 5 s) -- there is no second term.
`python3 testsys/run.py e2e` (and `run.py all`) run the everyday cell
selection at that term; `python3 testsys/run.py release` runs the SAME term
over the SAME cell selection -- as of 2026-09-23 `matrix.RELEASE_ONLY` is
DELETED, not empty: numpy dropped from every gate that day (owner decision,
"I actually don't care numpy... I will use Jax anyway"; PR #4, `91afba4`),
taking its only member (`test.tpv36`/`test.tpv37` x python-numpy) with it --
and is the tier that gates a release, writing
`docs/evidence/sweep-<shortsha>/summary.json`. A `--term` flag and a per-case
"full term" existed for less than a day (2026-09-23), each case whose full
term differed from 5 s carrying a second reference file
(`frt.canonical.term5.txt`); the owner retired both the same day (rule 7,
rule 24; `pathway_forward.md` item 102), so every gated case now has exactly
ONE committed reference (`frt.canonical.txt`, plus `fault.dyna.r.nc` for
`test.tpv36`/`test.tpv37`), gated at 5 s regardless of that case's own
`par.term`. **What this gives up**: `test.tpv29` (par.term 20 s),
`test.tpv36`/`test.tpv37` (par.term 6 s) are not run to their full physical
term by any gate -- the owner's own accounting puts 47-63% of the fault
rupturing after the 5 s cutoff across those three cases, and none of that
late-time rupture is gated anywhere. `testsys/e2e/full_specs.py` still
records each case's official spec resolution and duration; those runs stay
defined and unrun (the separate, opt-in, report-only `e2e-full` tier below,
not a release gate). The TPV29/TPV30 cross-code overlay against the 2015 SCEC
submissions (`case_input/test.tpv29/README.md`) remains the only check on
physics past 5 s, and it is manual. The owner accepted this tradeoff
explicitly on 2026-09-23.

The run prints which cells it covered and which it did not, so a green result
states its own scope. **CI no longer runs this sweep for physics coverage.**
CI's one e2e job is a portability SMOKE test -- `testsys/matrix.py`'s
`CI_CELLS` is exactly `test.tpv8` x `{fortran, python-jax}` at
the gate term (numpy dropped the same day it dropped from every other gate,
2026-09-23) (`.github/workflows/test.yml`: `build`, `unit-regression`
(2026-09-23: SHARDED 3 ways -- `testsys/ci_shard.py run 1`/`2`/`3` in
parallel, statically partitioning every `testsys/regression/test_*.py` plus
the `pytest testsys/unit` tier; `testsys/ci_shard.py verify` and
`testsys/regression/test_ci_shard_coverage.py` guard that the 3 shards'
union always covers every script and never double-runs one -- together the
same coverage as one local `testsys/run.py unit regression`, which stays the
unsharded local command), then `e2e-ci-smoke`
[`testsys/run.py e2e-ci`, i.e. `run_e2e.py --ci`] in parallel after `build`) --
the same total coverage as one local invocation of
`python3 testsys/run.py unit regression e2e-ci`. Its purpose is a clean
checkout with freshly `pip install`ed dependencies and (for the fortran cell)
the CI runner's own mpich, distinct from this dev box's Open MPI 4.1.1 -- not
a claim about the other 9 cases or the other two backends' physics, which the
non-CI `e2e` and `release` tiers cover instead.

Wall-clock figures below predate the ONE-term simplification (2026-09-23) and
were measured running each case at its own committed `par.term` -- a wall
clock `run.py release` no longer produces, since release now runs the SAME
5 s `GATE_TERM_S` as `run.py e2e`/`run.py all` over the SAME cells (`matrix.RELEASE_ONLY`
is deleted, not term-differentiated; see above and
`pathway_forward.md` items 102 and 104). `test.tpv29`'s gate-term cells are therefore
faster than the 148.3/1108.5/229.4 s shown below; that faster number has not
been measured yet and is not estimated here (rule 6).

| case | physics | SCEC benchmark | full-term wall clock, historical (dx / ranks / fortran s / numpy s / jax s) | full-tier spec (dx/term) |
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
v5.10.0), and `test.tpv30` is the eleventh (fortran + python-jax, `abs-max`
gate at 1e-10, registered 2026-09-23, PR #5 `567e723` -- see
`pathway_forward.md` item 19(b)). All three are absent from the table above
only because their per-backend fast-tier wall times have never been recorded
as a sweep measurement -- an unrecorded number is left unrecorded rather than
estimated (rule 6).
`testsys/matrix.py` is the authoritative cell list: 11 cases, `CASE_BOUND` and
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
python3 testsys/run.py all # unit + regression + the sweep: 11 cases x
                           # {fortran, python-jax} against one canonical
                           # reference each, plus test.tpv8 x python-jax-mpi
                           # = 23 cells (numpy dropped from every gate and
                           # test.tpv30 registered, both 2026-09-23;
                           # `--backend numpy` still runs by hand, ungated).
                           # Cells run concurrently; the 30-cell/2263.3-
                           # 2123.2s figures below predate both changes and
                           # are historical.
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
`MEASURED_PEAK_RSS_GB` is the authoritative measured record, kept for the
local `run.py all`/`run.py release` sweeps to reason about which cells fit
which box. It stopped being the reason CI's cell list is what it is on
2026-09-23 (`ef7196c`): CI no longer covers physics at all, it runs one fixed
portability smoke set, `matrix.CI_CELLS` -- `test.tpv8` x {fortran, python-jax}
(numpy dropped from this set 2026-09-23 along with every other gate) --
chosen for portability, not memory):

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
| 14 | `ERR_CFG_PROFILE_ENV_INVALID` | EQDYNA_PROFILE is set to a value other than unset, "", "1" or "0" |

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
