# Row 142: Vista one-node repeat, fortran, test.tpv8 strong medium -- 2026-11-25

Evidence only. Question: turn the single-run Vista smoke (row 142, `428b030`, 17.49 ms/step by difference, no repeat, no element
count) into a repeated measurement with the overhead split and the element count. Command (run on the Vista login node in the
suite's venv; the allocation was passed on the `--account` command line and is in no file here):

```
python3 testsys/perf/hpc_scaling_suite.py --submit --machine vista --account <ALLOC> --backend fortran --test strong --level medium \
  --max-nodes 1 --repeats 5 --out-dir $SCRATCH/eqdyna-scaling/rep1 --work-root "$SCRATCH/eqdyna-scaling/rep1" \
  --bin "$SCRATCH/eqdyna-vista-cur/bin/eqdyna" --module netcdf --walltime 00:30:00
python3 testsys/perf/hpc_scaling_suite.py --collect --suite-dir $SCRATCH/eqdyna-scaling/rep1 --out <tgz>
python3 testsys/perf/hpc_scaling_suite.py --analyze <tgz>      # here
```

Build: master `dcf50e3` (provenance `git_sha`, `git_dirty` false), `./install-eqdyna.sh -m vista` on the login node, aarch64 ELF, rc 0.
Job: one node, partition `gg`, 144 ranks (decomp 6x6x4), COMPLETED, 80 s wall (`sacct`); 5 repeats x {100, 200} steps + the parity point.
Files: `report.txt` (the `--analyze` output, paths scrubbed to `$SCRATCH`), `manifest.json` (scrubbed). The raw tarball (1444 profile
files, 174 KB) stays on Vista scratch; its profiles carry the run directory, which names the user.

## Result (MEASURED: the same script re-run with repeats on one node; no independent oracle for per-step cost, so not VERIFIED)

| quantity | value |
|---|---|
| mesh elements (the solver's grid lines, manifest `elements`) | 2,874,510 (19,962 per rank) |
| dx | 170.4545 m |
| per-step by difference, mean over 144 ranks' own solve times | 19.136 ms (max-rank 19.625 ms, max/mean 1.00003) |
| per-step, 5 repeats | 17.826, 19.625, 19.170, 19.465, 19.595 ms; relative range 9.4 % (first repeat lowest) |
| first smoke (single run, `428b030`) | 17.49 ms: inside the repeat range (17.83 lowest), 8.6 % below the repeat mean |
| cost per element-step per core | 0.959 microseconds (19.136e-3 s x 144 / 2,874,510) |
| bucket split per step (mean) | element 15.08 ms, wait 2.78 ms, exchange 0.91 ms, fault 0.051 ms |
| overhead per rank, mean (min-max) | total 1.784 s (1.567-3.280); MPI init 1.636 s (1.428-1.853); setup 0.125 s (0.103-0.139); io 0.023 s (0-1.547) |
| setup detail, mean | input read 0.038 s, mass assembly 0.015 s, meshing 0.011 s, netCDF open 0.061 s, MPI init 1.636 s |

The fixed cost is dominated by MPI init (92 % of the mean overhead; 1.64 s equals 86 steps of solve). `exchange` is below the noise floor on
65 of 720 rank-repeat samples (reported by `--analyze`, not dropped).

Parity point: FAIL, unchanged and not touched: frt max|diff| 1.0e-07 against the 1.0e-08 bound at row 940 column 13 (the one entry
diagnosed 2026-11-25 as 1.3e-14 of its column scale, cross-architecture round-off); `fault.dyna.r.nc` missing (the suite never runs
plotRuptureDynamics); every station comparison ok (worst e 1.6e-10 against 1.0e-07). The parity scheme is the owner's; no bound or
reference was changed.

Gap found: `--analyze` does not print `elements` in its report (the manifest carries it): the number above is read from `manifest.json`.
