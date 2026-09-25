# test.tpv30 python-vs-Fortran divergence — root cause, 2026-09-23

Base: origin/master `7484b48`. Worktree `agent-ac002d35f6721ffd3`.
Predecessor record: `NOTES_tpv30_gate.md` (2026-09-17), which localised the
divergence only to "between t=1 s (bit-exact) and t=6 s (widespread)".

## Root cause

`src/python/eqdyna/assembleGlobalKU.py`, `_pml()` — the gravity body force of
a **PML element** was scattered only into the 3-dof-node group sum
(`groups[2]` → `idxP3[2]`) and never into the twelve per-slot blocks
(`idxP12`) that a **12-dof PML node** receives.

Reference, `src/fortran/assembleGlobalKU.f90`:

```
:19-20   al(1:ned,1:nen) = rdampm*velArr(...)                 ! rdampm == 0.0d0 always
         al(3,1:nen)     = al(3,1:nen) + (1-C_elastic)*grav*(roumax-(gamar+1)*rhow)/roumax
:23      call calcElemMass(elemass(1:nee,nel), al, elresf)    ! elresf(3,j) = -grav_const*elmass
:44          efPML((i-1)*12+9+j) = elresf((i-1)*3+j)          ! SEEDS slots 10,11,12
:48      call calcPMLElemKU(...)                              ! f(..+10..12) -= det*w*(b*s0)
:51-58       if (ndof(node)==12) then
                 do j = 1, 12                                 ! <-- ALL TWELVE, incl. slot 12
                     nodalForceArr(eq) = nodalForceArr(eq) + efPML((i-1)*12+j)
:59-65       elseif (ndof(node)==3) then                      ! the group sums the port had
```

So Fortran gives a 12-dof PML node the gravity term through slot 12; the port
gave it to nobody. `grav_const` carries a `(1 - C_elastic)` factor and is
EXACTLY 0.0 for every elastic case, which is why nine of ten gated cases never
saw it.

## How it was found (method)

1. **Bisected in TIME, per step, not by re-running.** Added an env-gated dump
   to BOTH sides at the same point in the step (after `faulting`, after the
   `driver.f90:30` mass divide): per-fault-node tractions + cumulative slip
   every step, and the full interior-element stress field at steps 1-4.
   Fortran: `dumpFaultTract` in `src/fortran/driver.f90` (diagnostic only,
   `EQDYNA_DUMP_FAULT=1`; REMOVED before the commit). Python:
   `scratch_hunt/pydump.py`, monkeypatching `backend.time_loop`.
2. Case: `case_input/test.tpv30` at dx=500, **forced serial** (`par.nx=ny=nz=1`,
   `mpirun -np 1`) on both sides, ONE `case.setup`, the python dir a byte copy
   of the Fortran dir, so the two runs read identical input bytes and share a
   fault-node ordering (verified: fault-node coords max|diff| = 0.0, element
   centroids max|diff| = 0.0, 3321 fault nodes and 635040 interior elements on
   both sides).
3. Per-step relative difference of the interior-element stress sum of squares:

   | step | rel. diff | max fault-traction diff |
   |---|---|---|
   | 1 | 1.08e-12 | 2.4e-07 Pa |
   | 2 | 8.18e-13 | 3.6e-07 Pa |
   | **3** | **1.03e-05** | 3.6e-07 Pa |
   | 4 | 4.81e-05 | 4.8e-07 Pa |
   | 31 | 2.39e-02 | 1.3e-04 Pa |
   | 42 | 4.30e-02 | 2.4e+00 Pa (81 nodes > 1 Pa) |
   | 144 (t=6 s) | 2.7e-01 | 3.4e+07 Pa (2259+ nodes) |

   A seven-order jump at step 3 off a 1e-12 noise floor: an ONSET, not growth.
4. Full element-field diff at step 3 named the elements: 35376 of them, all on
   the interior/PML interface (the entire bottom interior layer, 10080
   elements, plus a 408-element perimeter ring on every other z-level). Worst
   element centroid `(-31872.67, 304.29, -35233.37)`:

   ```
   Fortran  sxx=syy=szz=-5.7674909688e+08   (isotropic: below the 22 km
                                              deviatoric taper, devStr == 0)
   port     sxx=syy=  -5.7678196298e+08     szz=-5.7684768362e+08
   ```

   The port's excess is `dszz = 3 * dsxx`, `dsxx == dsyy` — exactly
   `(lam+2mu)/lam` for this material (lam = mu = 3.2e10), i.e. a spurious
   **vertical** strain rate. Gravity acts in z.
5. Timing corroborates the chain: a wrong force on a 12-dof PML node at step 1
   gives a wrong PML-node velocity at step 2, which gives a wrong PML element
   force on the 3-dof interface nodes at step 2, which first reaches an
   INTERIOR element's stress at step 3. Step 3 is the earliest step at which
   this defect can be visible in the interior field, and it is the step it
   appears.

## The fix

`src/python/eqdyna/assembleGlobalKU.py`, `_pml()`: seed the f12 block with
`-(m_e_p * grav_const)` before it is scattered to `idxP12[11]`, mirroring
`assembleGlobalKU.f90:44`, and drop the separate `groups[2]` addition (the
group now picks it up from the block, which is also the Fortran's own
association: `((f7+f8)+f9) + (gravity + t)`).

**The PORT was changed, not the Fortran.** The Fortran is the oracle for this
campaign, and its behaviour here is also the physically correct one: a PML node
carrying mass must carry its own body force.

Verification, same instrumented 4-step run, same inputs:

```
                    before fix        after fix
nt=1   max|dstress|  0.000000e+00     0.000000e+00
nt=2   max|dstress|  1.192093e-07     1.192093e-07
nt=3   max|dstress|  9.858674e+04     2.384186e-07   (35376 -> 0 elements > 1 Pa)
nt=4   max|dstress|  4.498748e+05     2.384186e-07   (69448 -> 0 elements > 1 Pa)
```

2.4e-07 Pa on a 5.8e+08 Pa field is 4e-16 relative — the roundoff floor.

## Gate evidence (all run in this worktree, after the fix)

**test.tpv30, full 20 s, against the ALREADY-COMMITTED
`test.reference.results/test.tpv30/frt.canonical.txt` (NOT regenerated).**
Fortran run fresh at 4 ranks in the same campaign as a check that the
committed reference is still the oracle; python cells serial, 3321 fault
nodes each:

```
fortran-4rank   max|diff| = 0.000000e+00     (reference reproduces exactly)
python-numpy    max|diff| = 5.184742e-14     was 4.035089e+08
python-jax      max|diff| = 6.300516e-14     was 4.035089e+08
```

Per-column: columns 3 (rupture time) and 11/12/13 (normal/strike/dip
traction) are EXACTLY 0.0 on both python cells. The residual is entirely in
columns 4-9 (slip/slip-rate) at 1e-14..1e-16 absolute, on nodes whose values
are themselves ~1e-13 — the E18.7E4 print floor, not physics.

**Registered cells (`testsys/e2e/run_e2e.py`), bounds as committed:**

```
test.tpv8     python-numpy  SUCCESS  max|diff|=3.861189e-10  bound=1.0e-08
test.tpv8     python-jax    SUCCESS  max|diff|=1.983643e-10  bound=1.0e-08
test.tpv104   python-numpy  SUCCESS  max|diff|=4.152048e-07  bound=1.0e-05
test.drv.a6   python-numpy  SUCCESS  flips=332/450  median|dfnft|=0.0417s/0.1s  phys_max=2.892923e+07/1.1e+09
test.drv.a6   python-jax    SUCCESS  flips=377/450  median|dfnft|=0.0417s/0.1s  phys_max=1.410980e+07/1.1e+09
python3 testsys/run.py unit regression  -> SUCCESS unit, SUCCESS regression
```

**test.drv.a6 before/after, same box, same session, same committed reference
and the same unchanged 450 flip budget** (measured by reverting only
`assembleGlobalKU.py` and re-running):

```
                 flips before -> after      phys_max before -> after
python-numpy       415  ->  332             2.894505e+07 -> 2.892923e+07
python-jax         423  ->  377             1.691466e+07 -> 1.410980e+07
```

Both pass either way, so this is not a gate change — but the fix moves both
cells TOWARD the reference, and python-numpy's 332 is now essentially at the
329 flips that `matrix.py`'s own DRV_A6 comment records as the floor imposed
by DECOMPOSITION ALONE (serial Fortran vs the committed 4-rank reference).
No bound and no reference was touched.

## Blast radius

`grav_const == 0.0` exactly whenever `C_elastic == 1`, and `x + (-0.0) == x`
for every x, so the change is bit-for-bit identity on every elastic case. It
can only reach python cells whose case sets `C_elastic = 0`.

## Still NOT done by this note

`test.tpv30` remains UNREGISTERED in `testNameList.py` and `testsys/matrix.py`.
No reference was regenerated. No bound was calibrated.

