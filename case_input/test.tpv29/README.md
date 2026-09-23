# test.tpv29

SCEC TPV29: a vertical right-lateral strike-slip fault on a self-similar
**fractal rough surface** (Hurst 1), 40 x 20 km, reaching the free surface, in
a homogeneous elastic halfspace with slip-weakening friction. SCEC supplies the
exact surface at 25 m sampling so every code runs the same roughness.

## Fault geometry: source, shipped files, regeneration

**Authoritative source.** `tpv29_tpv30_geometry_25m_data.txt`, from
https://strike.scec.org/cvws/tpv29_30docs.html ->
`download/tpv29_tpv30_geometry_25m_data.zip` (66 MB unzipped). Not in the repo.

**Shipped files.** The *surface* column of both is an *exact decimation* of
that 25 m data -- every `y` value is an official one, nothing is interpolated
-- so a clean checkout runs both tiers with no download:

| file | spacing | used for |
|---|---|---|
| `bFault_Rough_Geometry.tpv29.50m.txt` | 50 m | full tier (the spec resolution) |
| `bFault_Rough_Geometry.tpv29.100m.txt` | 100 m | 100 m runs and every coarser gate/trial dx |

**The two derivative columns are RECOMPUTED, not decimated, and that is not a
detail.** `y` is a sampling of a surface, so every n-th value is still an
official value. `dy/dx` and `dy/dz` are not values of that surface: a
finite-difference derivative is a property of a GRID AT A SPACING, and
subsampling one leaves a slope describing a surface sampled n times finer
than the mesh it is attached to. The official columns themselves are
`np.gradient(F, 25)` -- verified 2026-09-22 against SCEC's own file, median
relative difference 0.0215% (dF/dx) / 0.0304% (dF/dy), corr 1.000000 over all
1601x801 nodes -- which is exactly why they are specific to the 25 m grid. So
every derivative this compset emits is `np.gradient(y, dx, edge_order=2)` at
the TARGET spacing: the same operation SCEC performed, at the mesh's own dx
(`tpv29GeometryTools.meshConsistentDerivatives`). `convertOfficial25m`
refuses to run unless recomputing at 25 m reproduces SCEC's supplied columns,
and `testsys/regression/test_rough_fault_normal_consistency.py` holds the
whole chain. Until 2026-09-22 both call sites subsampled the derivative
columns instead; at the 500 m gate that attached 25 m normals to a 500 m
mesh.

**Regenerating.** `tpv29GeometryTools.py` is the converter, and it validates
whatever it writes (`lib.validateFaultRoughGeometry`) before reporting success:

```bash
python3 tpv29GeometryTools.py                 # bFault_Rough_Geometry.txt for par.dx
python3 tpv29GeometryTools.py --dx 100        # ... for an explicit dx
# rebuild a shipped file (or any new resolution) from the official download:
python3 tpv29GeometryTools.py --official tpv29_tpv30_geometry_25m_data.txt \
        --dx 50 --out bFault_Rough_Geometry.tpv29.50m.txt
```

**Verifying the provenance chain.** Re-run the last command against the
official download and `cmp` the result with the shipped file -- they are
byte-identical (and it will print the derivative-faithfulness percentages
above before writing anything). Failing that, `lib.validateFaultRoughGeometry` checks any
file's header, row count, finite values and derivative-vs-finite-difference
consistency on its own, and every generated file carries a
`.provenance` sidecar naming its source and native spacing.

**Resolution rule.** A supplied surface can only be *coarsened* exactly;
refining it would interpolate roughness the benchmark does not have. So
`par.dx` must be a multiple of the finest shipped spacing (50 m) that divides
the 40 x 20 km fault. The compset declares that spacing as
`par.faultGeometrySourceDx`, and `lib.requireFaultGeometryResolution` refuses
anything else up front, naming the requested dx, the source spacing and what
would satisfy it. (Before v5.6.0 only the 100 m surface shipped, which made the
50 m full tier impossible to set up. Decimations of the 50 m file to
100/200/500 m are bit-identical to decimations of the 100 m file, so no gated
result moved.)

**`case.setup` does not clobber your file.** The geometry is written by
`par.faultGeometryWriter`, which `case.setup` calls only when
`bFault_Rough_Geometry.txt` is missing or does not match the case. A correct
file you placed yourself is kept byte for byte, and re-running `case.setup` is
idempotent.

**Convention for other supplied-geometry compsets** (TPV30 uses this identical
surface): ship the resolutions your tiers need, ship or reference the converter
that derives them from the authoritative source, set
`par.faultGeometrySourceDx` / `par.faultGeometryWriter`, and document the
source and the regeneration command here.


![reference result at fast-tier resolution](cRuptureDynamics.png)

| tier | dx (m) | term (s) | ranks (decomp) | ~cells | wall time | reference |
|---|---|---|---|---|---|---|
| fast (gate) | 500 | 20 | 4 (2,1,2) | 0.88 M | ~3 min | frozen golden, `test.reference.results/test.tpv29/` |
| full (spec) | 50 | 20 | see note | 119 M | ~0.7 h on 1024 ranks (est., HPC) | report-only; geometry from the shipped 50 m file |

**`ny=1` is deliberate**: the fault plane sits at y=0, and a partition boundary
landing on it used to halve on-fault tractions (fixed; a hard-stop guard now
refuses that configuration). Keep the fault off rank boundaries — `ny=1`, or
`par.ymax /= -par.ymin`.

Coarse-resolution caveat: at 500 m the fractal surface is a 20:1 decimation of
the official data, i.e. a smoother fault. It is a valid regression baseline,
not the benchmark. Validation at matched resolution: the current code at 100 m
agrees with EQdyna's published 2015 100 m submission to **0.006 s median**
rupture time (500 m: 0.076 s); slip agrees between the 100 m and 500 m runs to
~4% (2.90 vs 2.78 m max).

Spec: 50 m preferred / 100 m acceptable, 0-20 s
(https://strike.scec.org/cvws/tpv29_30docs.html). Cross-code results:
https://strike.scec.org/cvws/metric_cvv1_u1/tpv29/metric_cvv1_tpv29_ac_0.html
Past submissions archived in `scec_archive/tpv29/`.

Provenance: reference frozen from EQdyna at the fault-MPI fix, cotopaxi,
gfortran 11.4.0 / Open MPI 4.1.1, 2026-09-14.
