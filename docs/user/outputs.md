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

## Multiple faults

A case with `par.ntotft > 1` (see [Parameters](parameters.md#multiple-faults);
`test.tpv22`/`test.tpv23` are the worked examples) tags every PER-FAULT output
filename with `ft<N>_` for fault `N >= 2`, and leaves fault 1's names
untouched -- a single-fault case's files are byte-for-byte unaffected:

* `faultst<tag><strike>dp<depth>.txt` -- on-fault stations, `tag` empty for
  fault 1, `ft2_` for fault 2 (e.g. `faultstft2_000dp100.txt`).
* `src_evol<tag><rank>` -- source-time-function restart data, same `tag`
  convention, per rank.
* `fault.dyna.r.nc` is the exception to the tags: ONE file for all faults, in
  EQquasi's restart layout -- untagged variables with dims
  `(nid_fault, dip, strike)`, each fault at the `(0, 0)` corner of its own
  slice, zero-padded to the largest fault. A fully dynamic cycle (`mode = 2`)
  reads EQquasi's `fault.r.nc` in the same layout.
* `SCECRuptureTime<suffix>.txt` and `cRuptureDynamics<suffix>.png`, written by
  `plotRuptureDynamics` -- same `suffix` convention.

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
