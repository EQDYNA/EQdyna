# Whole-fault time-series output (`outputFaultTimeSeries`): verification

Built with `install-eqdyna.sh -m ubuntu` (gfortran, Open MPI), 4 MPI ranks,
`OMP_NUM_THREADS=1`. `base` = the binary built from origin/master 42ce254;
`off` / `on` = this branch with the switch 0 / 1.

* `runall.sh CASE TERM base off on` -- create, set up, run, plotRuptureDynamics.
* `cmpdirs.py A B` -- byte comparison of every output except `src_ts*`, logs,
  timing/profile files and `bGlobal.txt` (which gains the switch line); a file
  differing only in its `# date =` header line is counted separately.
* `verify_ts.py CASE_DIR` -- checks `src_ts*` against every independent output.

## Results

`bitidentity.out`: test.tpv8 (outputGroundMotion on), test.tpv1053d (TP),
test.tpv104 (RSF slip law), test.tpv22 (two faults, outputGroundMotion on):
base vs off vs on -- no differing file; `src_evol*` (including
`src_evolft2_*`) byte-identical.

`verify.out` -- every check within half a unit of the last printed digit of
the text output it is compared with (`ok` = True):

* `src_evol<tag><rank>` == field 0 `sliprate`, bit for bit (tpv8, tpv22 both faults).
* last record == `frt.txt*` and `fault.dyna.r.nc`, every shared column.
* station files `faultst*` (both faults in tpv22): time, slip, slip rate,
  tractions, state, temperature, pore pressure. RSF/TP stations add the
  background creep rate to slip and slip rate; the check adds it back.
* slip_strike/dip vs cumulative sum of sliprate_strike/dip * dt: max 8.9e-15 m;
  slip_path vs cumsum(sliprate) * dt (slip weakening): 0.
* size: 8 bytes per field per node per step on every rank and fault.

A 7 s test.tpv22 run (420 steps, 15504 nodes per fault) gave the same verdicts
(frt/nc/station all within rounding, slip_strike integral 8.9e-15). Fault 2
does not nucleate in this case within 7 s (slip rate ~1e-15 m/s), so its
time-varying content is the tractions, which match the fault-2 station files.

Wall time on vs off: tpv8 10.14 s vs 10.04 s; tpv22 132.5 s on vs 133.8 s base.
