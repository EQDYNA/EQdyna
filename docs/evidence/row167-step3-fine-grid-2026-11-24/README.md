# Row 167 step 3: fine-grid pre/post evidence for the TPV36/37 reference change (PR #284)

Pre = master `d720df6` (mis-wired wedge cells), post = PR head `3628339` (geometric side test, PR #281).
Same case, same decomposition (16 ranks), same params; only the binary differs (md5 in `meta/binaries.txt`).
Cross-version consistency, not validation against an independent code.

## dx 250 m (par.dx 250, fitted to 249.999 m by the belt rule)

- TPV36 (term 13 s): rupture time median |dt| 0 s, max 0.72 s; 13868 vs 13888 ruptured nodes; station series differ at 0-3 km depth, converge below (0.003 at dp120/dp180).
- TPV37 (term 9 s): rupture time and slip IDENTICAL (`meta/frt_diff_tpv37_dx250.txt`, column 3 max diff 0), but final stress columns differ on 209 fault nodes at 0-7.2 km depth, up to 6.9 MPa in tnrm. So the two binaries DO produce different results; TPV37 simply does not rupture the affected shallow nodes within its term, so a rupture-time/slip overlay cannot see the defect.
- Why TPV37 looked unchanged: the pre-fix mesh carries the defect at every dx (`meta/prefix_mesh_screen.txt`, screen-only build at ecb8ccb, pre-fix meshgen): sides_wrong_master = 120 (gate dx), 240 (250 m), 300 (200 m), identical for tpv36 and tpv37 (the mesh does not depend on the case beyond geometry, which both share).
- Gate dx reference diff (d720df6 vs 3628339): TPV36 and TPV37 each change 99 rows in columns 11-13 (final tnrm/tstk/tdip) at depth 0-7.2 km; rupture time (col 3), slips and rates unchanged to 1e-15.
- Figures: `dx250/tpv3?_dx250_cplot.png` (rupture-time contours pre vs post), `dx250/tpv3?_dx250_ts_own.png` (same-name station series, no interpolation, `scripts/stdiff.py`), `dx250/tpv3?_dx250_ts-fault.png`.
- `scec_compare.py --plot ts-fault` keeps only stations exactly on the node grid (its documented rule, `scec_plots.py:station_table`: never snapped); at dx 249.999 m 18 of 19 stations are 0.016 m off-grid and excluded, hence `ts_own` (same file name in both runs). The figure tool was not changed.

Runs: `scripts/run_case.sh`, `scripts/set_case.py`, `scripts/queue2.sh`; run metadata `meta/*.run.meta`.

## Independent four-cell re-run on the PR content head (conductor, 2026-11-24 21:52 CST)

`meta/four_cell_rerun_3628339.txt`: `run_e2e.py --cases test.tpv36,test.tpv37 --backends fortran,python-jax` on a clean tree at `3628339` (binary md5 in `meta/binaries.txt`, post): 4 of 4 cells SUCCESS against the regenerated references, exit 0, 271.6 s. Evidence-only commits since `3628339` change nothing under src/ or testsys/.
