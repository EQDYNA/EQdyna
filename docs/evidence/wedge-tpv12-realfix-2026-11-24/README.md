# TPV12 60-degree wedge mesh, trace node, with the real fix (row 167 step 2) -- 2026-11-24

Evidence only. Question (owner): does the real fix (PR #281, geometric side test) remove the trace-node tension on the TPV12
60-degree wedge case (row 149 wedge test, `par.C_degen = 60`)? Runs: master `a22b9dd` (src identical to `97f9c62`), 4 ranks,
term 4.5 s, cases built exactly as `docs/evidence/row149-wedge-mesh-2026-11-11/` (`run_case.sh`, `set_case.py`, params file),
OMP/BLAS threads 1, `HWLOC_COMPONENTS=-gl,-opencl,-cuda,-nvml`. Column 8 of `faultst000dp000` is n-stress, + = tension.

```
faultst000dp000_dx500.txt: samples 206 t_end 4.505  n_max +5.0989 MPa at t=2.580  n(4s) -11.1134  n_min -11.1668  samples>0.1MPa 27  >1MPa 18  
faultst000dp000_dx200.txt: samples 514 t_end 4.496  n_max +9.8515 MPa at t=2.414  n(4s) -5.7244  n_min -7.6966  samples>0.1MPa 65  >1MPa 40  
```

Verdict: NO. The fix removes the grid-divergent steady tension (as built, PR #269 archive: +30.06 / +84.97 MPa at 4 s at 500 / 200 m;
now -11.11 / -5.72 MPa) but a transient tensile spike remains: +5.10 MPa at 2.58 s (500 m) and +9.85 MPa at 2.41 s (200 m), above
1 MPa for 18 and 40 samples. It grows with refinement (5.10 -> 9.85), cause NOT explained. The earlier scratch guard
(`wedge-trace-audit-dipsweep-2026-11-13`) gave +5.0989 / +9.8515 and -11.0436 / -5.7244 (500 / 200 m): the real fix reproduces the
guard's n_max to 5 figures at both grids and n(4 s) to 5 figures at 200 m (500 m differs by 0.07 MPa), so the guard was the same fix
in effect. Same-script re-run, MEASURED, not VERIFIED (no independent oracle at this resolution). Not measured: finer than 200 m,
why tpv36 is clean at every dip.
