# Wedge-mesh trace-node audit and dip sweep (row 149, P4 paused) -- 2026-11-13

**SCRATCH DIAGNOSTIC. No source, case, test, reference or gate changed. The "guard" below is a one-line scratch edit used only to test a hypothesis (`scratch_guard.diff`); it is NOT a fix, was never audited, and every number here is unaudited (single runs, no independent oracle).** Status: gated, not validated.

## Verdict

1. **Probe 1 (attachment audit, wedge mesh, dx=500 m): the first fault node pair below the trace has a footwall brick wired to the hanging-wall (master) id, at both 15 and 60 deg.** Counts (`audit_summary.txt`): as built, master-attached elements on the footwall side 2 of 22 at 60 deg and 2 of 22 at 15 deg, slave-attached on the hanging-wall side 0 of 20 at both; with the scratch guard 0 of 20 / 0 of 22 at both dips (42 distinct node-element pairs read each). Signature: nodal mass of pair 1 is split x0.625 slave / x1.375 master (as built) against x1.000 / x1.000 (guard). Figure: `fig_probe1_attachment.png` (60 deg as built, 15 deg as built, 60 deg guard; red edge = the footwall brick on a master id).
   Reading of the cause (`replaceSlaveWithMasterNode`'s vertical-fault "one cell beside the plane" clause firing on a wedge mesh) is the predecessor's inference from the guard's effect; the code line was not re-audited this session.
2. **On the TPV12 wedge at 60 deg the guard removes most of the trace tension, not all** (`guard_tpv12.csv`, dp000 n-stress, + = tension, evaluated at 4 s, 206 and 514 samples):

   | dx | unguarded at 4 s (PR #269 archive) | guard at 4 s | guard max over time | trace slip at 4 s, unguarded -> guard (insertion mesh ref) |
   |---|---|---|---|---|
   | 500 m | +30.0603 MPa | -11.0436 MPa | +5.0989 MPa | 1.9187 -> 7.042 m (6.527) |
   | 200 m | +84.9702 MPa | -5.7244 MPa | +9.8515 MPa | 1.6545 -> 7.995 m (7.640) |

   The residual transient spike (peak at 2.58 s / 2.41 s) exceeds +0.1 MPa for 0.57 s (500 m) and 0.65 s (200 m), +1 MPa for 0.37 s / 0.39 s (88 / 182 samples positive, 18 / 40 above +1 MPa). Guard also restores trace slip to the insertion-mesh level (insertion values from `row149-wedge-mesh-2026-11-11/wedge_vs_insertion.csv`).
3. **Probe 2 (test.tpv36, dip changed only): no trace tension at 30, 45, 60 deg** (`dip_sweep.csv`, 6 runs; the pre-written criterion was n_max > +1 MPa on both grids with the finer larger: met at no dip). n_max: 30 deg 4.5e-5 (500 m) / 1.7e-6 (250 m); 45 deg 1.3e-8 (500 m); 60 deg 4.2e-6 (500 m) / 6.9e-7 (250 m) MPa. n at term is compressive, -7.4 to -17.1 MPa. 15 deg 500 m reproduces the earlier evidence (n_max +0.0192, slip 0.0892 m, arrival 9.059 s). Trace/1 km-down-dip slip ratio 0.17-0.30 (figure `fig_probe2_nstress_vs_dip.png`). The same mis-wiring exists at 15 deg (probe 1) with no tension, so the mis-wiring alone does not give tension on tpv36.
4. **Extra probes (outside the owner's two, flagged):** X1, test.tpv36 at 60 deg, 500 m, along-dip shear sign flipped (normal sense): trace n_max +5.517 MPa unguarded (slip 0.0117 m at 13 s) and +4.662 MPa guarded (slip 0.2324 m): tension in the normal sense with and without the guard, so "mis-wiring bites only in the normal sense" is refuted (run dirs `x60n`, `xg60n`; no CSV written for these, numbers read from the raw station files and the scratch notes). X2 (15 deg, normal sense) is inconclusive: n_max +0.206 MPa and trace slip 4.5e-17 m, it did not rupture.

## Not run / differences from the brief's plan
- **45 deg at 250 m: not run** (stopped unfinished). **15 deg at 250 m: not rerun** (stopped at 3.8 of 13 s); the value used is the archived `wedge-trace-tpv36-37-2026-11-12/trace_summary.csv` row (n_max +0.0023198, slip 0.12716 m, arrival 7.371 s, 16 ranks). So 45 deg and 15 deg have one grid each in `dip_sweep.csv` and the script's tension test is evaluated at 30 and 60 deg only.
- **Amendment 1 (written before any fine-grid run):** 45 deg at dt factor 0.5 diverged (1e69 at 12.99 s), so 45 deg uses factor 0.25. Control: 60 deg 500 m at 0.25 gives n_max +6.8e-6 MPa, n at term -10.12 MPa, slip 0.1083 m, against 4.2e-6, -10.07, 0.1049 at 0.5: insensitive to dt. Criteria were unchanged.
- Forced extra change: model bottom held 7.753 km below the fault bottom (28 km fault does not fit above -15 km at 45+ deg). A first 15 deg run was spoiled by editing a running script and rerun.
- Spike duration in item 2 was re-measured here from the raw station files; unguarded values are quoted from PR #269's archive, not rerun. The earlier claim "n_max ~1e-8 to 1e-5" is slightly low: largest is 4.5e-5.

## Unexplained
- Why tpv36 is clean at every dip although it carries the same mis-wiring, while TPV12 wedge at 60 deg has +30 MPa (as built) and, with the guard, a +5 to +10 MPa transient. X1 shows the sense of slip is not the difference; what is (TPV12 stress law, surface strength, dy=dx/2 cells, ...) is not isolated. PR #274 already excluded surface strength.
- The residual guarded spike.
- Whether the guard is a correct treatment of the vertical-fault clause on wedge meshes (it only tests the hypothesis).

## Files
`audit_summary.txt`, `*_attach.csv`, `*_pairs.csv` (a=as built, g=guard; 15/60 deg), `trace_audit_ana.py`, `fig_probe1.py`; `dip_ana.py` (regenerates `dip_sweep.csv`, `guard_tpv12.csv`, fig_probe2; regenerated this session, identical to predecessor's), `run_dip.sh`, `set_dip_case.py`; `scratch_guard.diff` (clean vs guard tree, -U0, srcStamp line dropped; includes the probe-1 dump routine `diagAudit` and the one-line guard). Raw runs and the notes (`NOTES_audit.md`, with the kill criteria written before the runs) stay in the untracked scratch area.
