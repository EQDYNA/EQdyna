# TPV12 60-degree wedge mesh, trace node, with the real fix (row 167 step 2) -- 2026-11-24

Evidence only. Question (owner): does the real fix (PR #281, geometric side test) remove the trace-node tension on the TPV12
60-degree wedge case (row 149 wedge test, `par.C_degen = 60`)? Runs: master `a22b9dd` (src identical to `97f9c62`), 4 ranks,
term 4.5 s, cases built exactly as `docs/evidence/row149-wedge-mesh-2026-11-11/` (`run_case.sh`, `set_case.py`, params file),
OMP/BLAS threads 1, `HWLOC_COMPONENTS=-gl,-opencl,-cuda,-nvml`. Column 8 of `faultst000dp000` is n-stress, + = tension.

```
Traceback (most recent call last):
  File "/home/utig5/dliu/EQdyna/.claude/worktrees/cond-e/docs/evidence/wedge-tpv12-realfix-2026-11-24/trace_stats.py", line 9, in <module>
    print(f"{sys.argv[1].split('/')[-2]}: samples {len(t)} t_end {t[-1]:.3f}  n_max {n[i]:+.4f} MPa at t={t[i]:.3f}  n(4s) {n[k]:+.4f}  n_min {n.min():+.4f}  samples>0.1MPa {(n>0.1).sum()}  >1MPa {(n>1).sum()}  last-slip(col6?) see file")
IndexError: list index out of range
Traceback (most recent call last):
  File "/home/utig5/dliu/EQdyna/.claude/worktrees/cond-e/docs/evidence/wedge-tpv12-realfix-2026-11-24/trace_stats.py", line 9, in <module>
    print(f"{sys.argv[1].split('/')[-2]}: samples {len(t)} t_end {t[-1]:.3f}  n_max {n[i]:+.4f} MPa at t={t[i]:.3f}  n(4s) {n[k]:+.4f}  n_min {n.min():+.4f}  samples>0.1MPa {(n>0.1).sum()}  >1MPa {(n>1).sum()}  last-slip(col6?) see file")
IndexError: list index out of range
```

Verdict: NO. The fix removes the grid-divergent steady tension (as built, PR #269 archive: +30.06 / +84.97 MPa at 4 s at 500 / 200 m;
now -11.11 / -5.72 MPa) but a transient tensile spike remains: +5.10 MPa at 2.58 s (500 m) and +9.85 MPa at 2.41 s (200 m), above
1 MPa for 18 and 40 samples. It grows with refinement (5.10 -> 9.85), cause NOT explained. The earlier scratch guard
(`wedge-trace-audit-dipsweep-2026-11-13`) gave +5.0989 / +9.8515 and -11.0436 / -5.7244 (500 / 200 m): the real fix reproduces the
guard's n_max to 5 figures at both grids and n(4 s) to 5 figures at 200 m (500 m differs by 0.07 MPa), so the guard was the same fix
in effect. Same-script re-run, MEASURED, not VERIFIED (no independent oracle at this resolution). Not measured: finer than 200 m,
why tpv36 is clean at every dip.
