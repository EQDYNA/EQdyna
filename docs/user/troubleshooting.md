# Troubleshooting

## Exit codes

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
| 15 | `ERR_CFG_NSTRESS_SIGN_INVALID` | bGlobal.txt's station n-stress sign is neither +1 nor -1 |

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
