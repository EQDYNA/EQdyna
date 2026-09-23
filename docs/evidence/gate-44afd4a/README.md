# Gate evidence for `44afd4a` (rough-fault normals recomputed at the mesh dx)

**WHAT:** the full-suite run that gated the merge of
`wei/rough-fault-normal-2026-09-22`. Run by the conductor on the MERGED tree
(branch + `origin/master` 83ce3c2, merge commit `64091b7`), not inherited from
the branch — the branch's own transcript is not a substitute for this run.

**FILES**
- `run_all_sweep_2026-09-22.log.gz` — `python3 testsys/run.py all`, 2,292,974
  bytes raw, md5 `5b9ac9352605714506b6c4917418281e`, round-trip verified.
  Built by `./install-eqdyna.sh -m ubuntu` in
  `/home/utig5/dliu/EQdyna.wt-wei-roughnormal` at `64091b7`. Exit 0.
  unit SUCCESS, regression SUCCESS, e2e SUCCESS; 31 of 40 cells ran, 31
  passed, 0 failed, 9 declared-unsupported (`python-jax-mpi` column),
  2706.1 s.

**THE CELLS THIS COMMIT IS ABOUT**

    test.tpv29  fortran       max|diff| = 0.000000e+00   bound 1.0e-10
    test.tpv29  python-numpy  max|diff| = 8.733900e-15   bound 1.0e-10
    test.tpv29  python-jax    max|diff| = 7.981116e-14   bound 1.0e-10

The Fortran cell is bit-exact against the reference regenerated in `bbe6f1a`,
i.e. a fresh run on the merged tree reproduces that reference exactly.

**THE CELLS THAT MUST NOT HAVE MOVED, and did not**

    test.drv.a6  fortran / python-numpy / python-jax   all SUCCESS

`test.drv.a6`, `liu2020.fdc.rough.250` and `bp1001.fdc.rough.250` build their
roughness with the fractal generator (`insertFaultType=2`) and never enter
`tpv29GeometryTools.py`. Their references and bounds are untouched by this
commit.

**THE BOUND ON THE CLAIM.** `test.tpv30` is UNREGISTERED (pathway 19(b)), so
its regenerated reference (`642f119`) is exercised by NO cell in this sweep.
What the sweep does cover for tpv30 is the geometry underneath it, via two
regression checks:

- `shipped_surfaces_mesh_consistent` — over all 3 shipped surfaces the stored
  `dy/dx`, `dy/dz` columns match `np.gradient` of their own `y` at their own
  header `dx` to max|diff| `1.454e-14` (worst:
  `test.tpv29/bFault_Rough_Geometry.tpv29.50m.txt`), bound `1.0e-04`.
- `tpv30_copy_is_identical` — `case_input/test.tpv30/tpv29GeometryTools.py` is
  byte-identical to `test.tpv29`'s, so those checks speak for both compsets.

Nobody should later read this gate as having validated tpv30's reference. It
did not, and cannot until the case is registered.

**VERIFY**

    gunzip -c run_all_sweep_2026-09-22.log.gz | md5sum   # 5b9ac9352605714506b6c4917418281e
    gunzip -c run_all_sweep_2026-09-22.log.gz | tail -20
