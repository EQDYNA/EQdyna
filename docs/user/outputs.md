# Output Files

A case directory fills up with several kinds of output while `run.sh` (or
a direct `python3 -m eqdyna` invocation) runs. This page describes each
file family and its columns; see Benchmarks for how these files are
compared against a reference result.

## On-fault station time series: `faultst<strike>dp<depth>.txt`

One file per requested on-fault station (`st_coor_on_fault` in
`user_defined_params.py`), for every MPI rank that owns that station. Each
file opens with commented header lines (starting with `#`) naming the
project, the element size, the time step, and the meaning of every column,
followed by one row per time step. With slip-weakening or time-weakening
friction the file has 8 data columns: time, along-strike slip, slip rate
and shear stress, down-dip slip, slip rate and shear stress, and normal
stress. A rate-and-state friction law adds 3 more columns: the state
variable, temperature, and pore pressure.

### Normal-stress sign convention

The normal-stress column's sign follows whichever convention that case's
own SCEC benchmark specification states -- some benchmarks report positive
as extension, others as compression, and the two disagree between
benchmarks. Set `faultStNormalStressSign` in `user_defined_params.py` to
`'extension'` or `'compression'` to match your case's own spec; `case.setup`
refuses to run with any other value, and both the Fortran and Python
solvers write that column with the sign this setting selects.

## Off-fault station time series: `body<y>st<x>dp<z>.txt`

One file per requested off-fault station (`st_coor_off_fault`). The file
name encodes the coordinates that were requested; the header's own
location line instead reports the mesh node the station was actually
matched to, which can differ slightly from the request on a coarse mesh.

## Fault restart/rupture file: `frt.txt<rank>`

One file per MPI rank that owns at least one fault node -- a rank whose
subdomain never touches the fault writes none. Each row holds one fault
node's coordinates, rupture time, final slip, slip rate, traction, and
friction-law state. `plotRuptureDynamics` and `plotSlipAndRPT` read these
files back in to produce the rupture-dynamics plots on the Benchmarks page.

## Field output: `fault.dyna.r.nc`

A NetCDF file of the whole fault plane's slip, slip rate, and traction,
sampled every `nt_out` time steps, for visualization and for resuming a
later stage of a multicycle run.

## Whole-fault time series: `src_ts<rank>`

Set `outputFaultTimeSeries = 1` in `user_defined_params.py` (default `0`,
off) to record every time-varying on-fault quantity at every fault node on
every time step -- the full history behind the station files and
`frt.txt`, so a later analysis does not need a rerun. Fortran solver only;
the Python solver refuses the setting. Each MPI rank that owns fault nodes
writes two files:

* `src_ts<rank>` -- binary, 64-bit floats in the machine's native byte order
  (little-endian on x86 and ARM), no record markers. One record per time
  step; inside a record, all fields of node 1, then all fields of node 2, and
  so on. Record `k` (counting from 0) is at time `(k + 1) * dt`.
* `src_ts_hdr<rank>` -- text. Lines starting with `#` give `nnode`, `nfield`,
  `dt`, `nstep`, `friclaw`, and one `# field <index> <name> <unit> <meaning>`
  line per field; the remaining rows are the `x y z` coordinates (m) of each
  node, in the binary's node order.

Fields, in file order:

| # | name | unit | meaning |
|---|------|------|---------|
| 0 | `sliprate` | m/s | slip-rate magnitude (the value in `src_evol`) |
| 1 | `sliprate_strike` | m/s | slip rate along strike |
| 2 | `sliprate_dip` | m/s | slip rate along dip |
| 3 | `slip_strike` | m | slip along strike |
| 4 | `slip_dip` | m | slip along dip |
| 5 | `slip_normal` | m | fault-normal relative displacement (opening) |
| 6 | `slip_path` | m | path-integrated slip, running sum of slip rate times `dt` |
| 7 | `sliprate_peak` | m/s | running maximum slip-rate magnitude |
| 8 | `traction_normal` | Pa | normal traction, negative in compression |
| 9 | `traction_strike` | Pa | shear traction along strike |
| 10 | `traction_dip` | Pa | shear traction along dip |
| 11 | `traction_shear` | Pa | shear-traction magnitude |

Rate-and-state friction (`friclaw` 3, 4, 5) adds `rsf_state` (the state
variable: theta in s for the ageing law, psi, dimensionless, for the slip
law), `theta_pc` (Pa) and `theta_pc_dot` (Pa/s), the normal-stress memory
variable and its rate, and the particle velocities of the two split nodes,
`vel_master_x/y/z` and `vel_slave_x/y/z` (m/s). Thermal pressurization
(`friclaw = 5`) also adds `tp_pore_change` (Pa, pore-pressure change from
its initial value) and `tp_temperature` (K). Constant inputs such as the
friction parameters and initial stresses are not repeated per step.

The tractions are the solver's values in Pa; the station files print the
same quantities in MPa with the case's sign convention. With rate-and-state
friction and a nonzero background slip rate, the station files add that
background rate to slip and slip rate; these files do not.

The size is 8 bytes per field per node per step: a 12-field slip-weakening
run with 10,000 fault nodes and 3,000 steps writes about 2.9 GB. Read the
files with `scripts/readFaultTimeSeries.py`:

```python
from readFaultTimeSeries import read_fault_ts
ts = read_fault_ts('src_ts_hdr0')      # the header names its binary
ts['slip_strike'].shape                # (nstep, nnode)
ts['xyz'], ts['t']                     # node coordinates, record times
```

or `python3 scripts/readFaultTimeSeries.py <case_dir>` for a per-file
summary. Every other output is unchanged by this setting.

## Multiple faults

A case with `par.ntotft > 1` (see [Parameters](parameters.md#multiple-faults);
`test.tpv22`/`test.tpv23` are the worked examples) tags every PER-FAULT output
filename with `ft<N>_` for fault `N >= 2`, and leaves fault 1's names
untouched -- a single-fault case's files are byte-for-byte unaffected:

* `faultst<tag><strike>dp<depth>.txt` -- on-fault stations, `tag` empty for
  fault 1, `ft2_` for fault 2 (e.g. `faultstft2_000dp100.txt`).
* `src_evol<tag><rank>` -- slip-rate magnitude at every fault node, every
  time step (64-bit floats, node order as in `src_ts`), written when
  `outputGroundMotion` is set; same `tag` convention, per rank.
* `src_ts<tag><rank>` and `src_ts_hdr<tag><rank>` -- the whole-fault time
  series above, same `tag` convention, per rank.
* `fault.dyna.r.nc` is the exception to the tags: ONE file for all faults, in
  EQquasi's restart layout -- untagged variables with dims
  `(nid_fault, dip, strike)`, each fault at the `(0, 0)` corner of its own
  slice, zero-padded to the largest fault. A fully dynamic cycle (`mode = 2`)
  reads EQquasi's `fault.r.nc` in the same layout.
* `SCECRuptureTime<suffix>.txt` and `cRuptureDynamics<suffix>.png`, written by
  `plotRuptureDynamics` -- `suffix` empty for fault 1, `_ft2` for fault 2
  (e.g. `SCECRuptureTime_ft2.txt`).

`frt.txt<rank>` is the one exception: it is NOT split per fault. One file per
MPI rank still covers every fault that rank owns nodes on, same as the
single-fault case (`library_output.f90`'s `output_frt`); the per-rank/per-node
row data already carries each row's own fault id.

## Other outputs

* `compTime<rank>` -- a per-rank wall-clock breakdown by simulation phase.
* `pstr.txt<rank>` -- off-fault plastic strain, written only when `output_plastic` is set.
* `gm<rank>` -- ground-motion velocity time series, written only when `outputGroundMotion` is set.
* `finalSurfDisp.txt<rank>` -- final surface displacement, written only when `outputFinalSurfDisp` is set.
* `cRuptureDynamics.png` -- the rupture-dynamics summary plot `run.sh` produces at the end of every run.
