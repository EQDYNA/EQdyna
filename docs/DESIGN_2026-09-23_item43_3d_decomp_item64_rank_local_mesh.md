# DESIGN — item 43 (3D decomposition) + item 64 Stages 1-3 (rank-local mesh)

Date: 2026-09-23. Branch base: `2e1624d`. Status: **DESIGN ONLY, nothing implemented.**
Author: dliu (agent). Everything below was read out of the tree at `2e1624d`; every
number is labelled measured-elsewhere or not-yet-measured. No run was made for this
document.

---

## 0. What is being decided

Two changes land together:

1. **3D split** — an n-rank jax-MPI run decomposes as `(npx,npy,npz)` matching
   Fortran (`testsys/perf/run_scaling.py:168` `DECOMP`), replacing the 1D contiguous
   element cut at `src/python/eqdyna/MPI4NodalQuant.py:244`.
2. **Rank-local mesh** — each rank generates only its own nodes/elements with
   rank-local numbering (Fortran's design: `src/fortran/meshgen.f90:65-67`,
   `src/fortran/countMeshEntities.f90:20-70`), replacing "build the proven serial mesh
   on every rank, then restrict".

Motivation (measured elsewhere, not re-derived here): the 1D slab tracks Fortran to 4
ranks and turns over between 4 and 8 (`docs/SESSION_LOG_2026-09-23_autopilot.md` §EE) —
halo holds ~50 ms/step at both 8 and 16 ranks while compute halves, which is what a
slab's non-shrinking surface does. Expected return from rank-local (item 64 Stage 0,
settled): **≤ 7.92 s of setup per run and ~51 GB of peak memory at 32 ranks, with NO
movement in the differenced per-step number.** Section 6 argues both halves of that are
optimistic and says what to measure first.

---

## 1. THE GATE QUESTION: does frt canonicalisation stay honest under rank-local numbering?

### 1.1 Answer

**Yes, and neither `testsys/frt_canonical.py` nor `testsys/compare.py` changes — but frt
alone stops being a sufficient gate, and that is the finding.**

Read the path, not the memory of it:

- `src/python/eqdyna/eqdyna3d.py:479` calls
  `library_output.write_frt(path, mesh['meshCoor'], mesh['nsmp'][rows], fnft_1idx, fric_1idx)`.
  `library_output.py:106-122` builds each row from `build_frt_rows(meshCoor, nsmp, ...)`:
  **columns 0-2 are physical coordinates fetched through `nsmp`.** A node number appears
  only as the ADDRESS used to fetch a coordinate; it is never written.
- `testsys/frt_canonical.py:131` keys the dedupe on `np.round(arr[:, :3], 6)` and
  `:150` lexsorts on `(x,y,z)`. No column of the canonical artifact is a function of
  node numbering, equation numbering, rank count, or element order.
- `testsys/compare.py:87` → `frt_canonical.canonical_from_case(run_dir)` → `:79`
  `glob('frt.txt*')`. The gate is already rank-count agnostic; that is why a 4-rank
  Fortran run (which ALREADY uses rank-local numbering — this is not new territory,
  Fortran has never had global numbering) compares to the same artifact as serial Python.

So the canonicalisation was never a statement about numbering. It is a statement about
**position**. Rank-local numbering does not weaken it.

### 1.2 What rank-local numbering actually gives up, stated precisely

The property being lost is NOT in the comparison, it is in the construction. Today
(`MPI4NodalQuant.py:35-52`, `:159-187`) every index array is the gated serial array under
an *injective relabelling*, so a partition bug can only produce a **loud** failure:
`_relabel` raises on any index the rank does not hold. Under rank-local generation,
`conn`, `eq_ids`, `nsmp`, `meshCoor` are **generated** locally, and three new failure
classes become reachable:

| class | what the frt gate does with it | verdict |
|---|---|---|
| **A. wrong coordinate fetched for a right node** (nsmp ↔ meshCoor misregistration) | `frt_canonical.align` (`:172-198`) gates coordinate agreement at `coord_tol=1e-9` m; a real misregistration is one grid spacing, i.e. metres — nine orders above the floor. Also `align` raises on any node-count mismatch. | **CAUGHT, loudly** |
| **B. right coordinate, wrong physics AT a fault node** (bad local eq number on or near the fault) | abs-max gate at the case bound. | **CAUGHT** |
| **C. right coordinate, wrong physics AWAY from the fault** (bad numbering in a rank's far interior or PML) | only observable at the fault, and only after the error propagates there within the gated step count. The gated cases are deliberately coarse and short (rule 17 step 4, "minutes at 4 ranks"). A defect in a far corner of a PML slab plausibly never reaches the fault in the gated window. | **NOT RELIABLY CAUGHT** |

Class C is the honest hole and it is the reason this section is the gate on the plan.
It is not fixed by changing the comparison — the comparison is right. It is fixed by
**adding gates that observe the mesh directly**, which are cheap because the serial mesh
still exists and is still gated:

- **U1 (unit, ms).** The 1D partition function (port of
  `meshgen.f90:543-550,588-596`) against a hand-written table of
  `(global_size, np, rank) -> (local_size, offset)` for every `(n, np)` the gated cases
  use. Pure integer arithmetic, no mesh.
- **U2 (unit, ms).** Each rank's 1D line is a **bitwise slice** of the global line:
  `np.array_equal(local_line, global_line[off:off+n_loc])`. See §3.2 — this is what
  protects `align`'s 1e-9.
- **U3 (unit, seconds — the replacement for the lost construction property).** For
  `test.tpv8` at `(2,2,1)` and `(2,2,2)`: build the serial mesh once, build each rank's
  mesh, and assert **bitwise** that
  `meshCoor_local == meshCoor_serial[:, g_of_l]`, `g_of_l[conn_local] == conn_serial[my_elems]`,
  `g_of_l[nsmp_local] == nsmp_serial[my_fault_rows]`, and that `eq_ids_local` is an
  injective relabelling of `eq_ids_serial[g_of_l]` preserving the sink at 0. `g_of_l` is
  the analytic `(ix,iy,iz)` map of §3.3 — i.e. the identity is CHECKED, where today it is
  constructed. This restores the "cannot introduce a numbering bug" guarantee at the
  rank counts where a serial mesh is affordable, and converts it into evidence.
- **R1 (runtime, every MPI run, O(1) reductions).** Conservation across the partition,
  which is what U3 cannot cover at 32 ranks: allreduce of uniquely-owned node count and
  element count against `N_global`/`E_global`; allreduce of uniquely-owned lumped mass
  against the serial total (Fortran's own reason for `MPI4arn`); allreduce of `arn` over
  fault nodes. These are physics-level invariants (mass conservation, partition of
  unity) and they go red on a duplicated or dropped entity that produced no frt symptom.

### 1.3 Two traps in the canonicalisation that this change walks into

- **Keep single-owner fault writing.** `MPI4NodalQuant.decompose:309-328` assigns each
  fault node to exactly one rank, so `canonicalize`'s duplicate branch
  (`frt_canonical.py:136-148`) never fires for python-mpi. Fortran's references DO carry
  duplicates (1.0-3.6%, docstring `:26`) because Fortran lets every touching rank write.
  Adopting Fortran's behaviour here would activate a check that requires the duplicated
  rows to agree at **exactly 0.0 spread in all 22 columns** — but a shared fault node's
  `arn` is a cross-rank SUM (`MPI4arn`, `meshgen.f90:212`) whose addition order differs
  per rank, so the last bits could differ and `canonicalize` would **raise**, correctly,
  and block the sweep. Do not switch to duplicate writing as part of this change.
- **`driver.run_mpi:382-388`'s owned-count allreduce becomes MORE load-bearing.** Today
  ownership is derived from the global `conn`, so it is consistent by construction. Under
  rank-local it is derived analytically from index boxes (§3.4) and can be wrong in both
  directions. Keep the check, and add its mirror: assert no fault node is owned TWICE
  (allreduce of an ownership count per global fault row, == 1 everywhere), not just that
  the total matches.

### 1.4 One consequence to state before anyone runs it

`test.tpv8` at 4 ranks is currently a 1D element cut; under `DECOMP` it becomes `(2,2,1)`
— a **different partition**, hence a different summation order at boundary equations,
hence **different last bits**. It must still pass at the case's bound. It will not be
bit-identical to today's 4-rank output, and the committed reference does not move
(rule 7). Anyone expecting bit-identity here is expecting the wrong thing; the
bit-identity claim in `MPI4NodalQuant._relabel`'s docstring (`:166-167`) is about the
relabelling at a FIXED partition, not across partitions.

---

## 2. The 3D split

### 2.1 Where `npx/npy/npz` come from

Not from a new table. `readInputFiles.build_params` already reads `npx/npy/npz` from
`bGlobal.txt`; `eqdyna3d.py:221-223` currently REFUSES anything but `(1,1,1)`. The design:

- `build_solver_state(case_dir, decomp=(1,1,1), rank=0)` accepts the decomposition;
  the `(1,1,1)` refusal narrows to "the serial entry point requires `(1,1,1)`".
- `run_case_mpi` reads `(npx,npy,npz)` from the case input and **refuses** unless
  `npx*npy*npz == comm.Get_size()`. No inference of a decomposition from the rank count
  inside the solver — a silent choice of partition is the class of thing that makes a
  measurement unattributable.
- The perf harness and the e2e cell set `npx/npy/npz` in `bGlobal.txt` from the same
  `DECOMP` table Fortran uses (`run_scaling.py:168`). Agreement with Fortran is then
  **structural — both backends read the same file** — not a coincidence maintained by
  two tables. That is what makes the Fortran-vs-python ms/step comparison like-for-like.

### 2.2 Rank → (mex,mey,mez)

Verbatim from `calcXyzMPIId` (`meshgen.f90:472-481`), z fastest:

```
mex = rank // (npy*npz)
mey = (rank - mex*npy*npz) // npz
mez =  rank - mex*npy*npz - mey*npz
```

This must match exactly, or python rank *r* holds a different subdomain than Fortran rank
*r* and any per-rank number (element count, halo size, ms/step) compares the wrong pair.
U1 covers it.

### 2.3 The 1D partition, verbatim

From `meshgen.f90:543-550` and `:588-596`, per dimension, with `n = global_size`,
`P = npx|npy|npz`:

```
per   = (n + P - 1) // P
resid = (n + P - 1) - per*P
n_loc = per if m < (P - resid) else per + 1
off   = (per-1)*m            if m <= (P - resid)      # 0-based offset into the global line
      = (per-1)*m + (m - P + resid)   otherwise
```

Note the **overlap of one node plane** between neighbours (stride `per-1`, size `per`):
node planes are shared, elements are not. That is the invariant the whole halo rests on —
each rank creates elements only for `ix,iy,iz >= 2` of its own slab
(`countMeshEntities.f90:59`, `meshgen.f90:86`), so the element sets partition exactly
while the shared plane carries partial nodal sums. Port the arithmetic literally,
including the `<` vs `<=` asymmetry between `:546` and `:588` (they are genuinely
different comparisons in the Fortran — reproduce, do not tidy).

---

## 3. Rank-local mesh generation

### 3.1 What changes, file by file

- `src/python/eqdyna/meshgen.py`
  - `one_dim_coor_array` (`:75`) gains `(m, P)` and returns the **slice** plus the offset.
    Its serial-only disclaimer at `:65-70` is retired.
  - `build_grid_lines` (`:153`) gains `decomp` and `(mex,mey,mez)`, returns the three
    local lines and their global offsets. `PMLb` and `model_bound` keep coming from the
    **global** arrays (`:136-149`) — they are model-wide constants and must not become
    per-rank, or `setNumDof`'s PML test (`meshgen.f90:606-608`) reclassifies nodes at
    subdomain boundaries.
  - `build_node_coordinates` (`:323`), `build_elements` (`:410`),
    `build_equation_numbers` (`:898`), `build_fault_geometry` (`:1140`),
    `build_station_matching` (`:1283`) are already written as functions of
    `(xline, yline, zline)`. **They need no structural change** — handed a local line they
    produce a local mesh with local numbering, exactly as the Fortran's identical loop
    nest does. This is the reason the change is affordable: the port's mesh builders are
    already parameterised the way Fortran's are.
  - Genuinely new: the boundary-condition test. `countMeshEntities.f90:34-35` fixes a node
    when it sits on a **model** boundary, compared against `modelBoundCoor` — which stays
    global. A local-line-derived `min/max` would fix every subdomain face and weld the
    model shut. This is the single most dangerous line in the change; U3 catches it
    (`eq_ids` would stop being an injective relabelling) and so does R1 (equation count).
- `src/python/eqdyna/MPI4NodalQuant.py`
  - `decompose`'s element cut (`:241-273`), touched-set restriction (`:279-286`),
    neighbour discovery loops (`:296-308`, `:312-316`) and the whole global→local
    relabelling block (`:330-368`) **are deleted**. They exist to recover locally what the
    rank-local builder now produces directly, and each of the two `for s in range(nranks)`
    loops is O(nranks × N_global) on the global `conn` — impossible once `conn` is local
    and unnecessary once the split is structured.
  - What remains and grows: the analytic neighbour/halo plan (§4), fault ownership (§3.4),
    `report`, `sync_mode`, `step_profile`, `exchange`.
  - `PML_WEIGHT` (`:85`) and `_cuts` (`:88`) **go away entirely**. The 3D split is
    geometric, not work-weighted; a structured split cannot honour an element-weight
    balance. Item 60 already measured that the weighted recut bought 1.4% at 32 ranks —
    i.e. nothing outside run-to-run variation — so this removes a knob that was shown not
    to matter. Say so in the commit rather than letting it look like a regression.
- `src/python/eqdyna/eqdyna3d.py` — `build_solver_state` signature (§2.1); `run_case_mpi`
  passes the decomposition and writes `frt.txt<rank>` from the **local** `meshCoor`/`nsmp`
  (already local-indexed after the change; `:479` needs no other edit).
- `src/python/eqdyna/driver.py` — `run_mpi`'s global mass check (`:393-399`) must be
  re-expressed as a rank-local check plus an allreduce of the bad count, or it reinstates
  a global array on every rank and gives back part of the memory this change is for.
  `mass_l` (`:404`) then comes straight from the local build instead of being gathered
  through `loc['local_eqs']`.

### 3.2 The coordinate-identity constraint (non-negotiable)

`getLocalOneDimCoorArrAndSize` allocates the **global** 1D array, fills it, and then slices
(`meshgen.f90:541,553-568,588-596`). Fortran's mesh is rank-local in its 3D arrays but its
coordinates come from a globally-accumulated line. The port must do the same.

Why it is not optional: the line is built by cumulative geometric stretching
(`arr[i] = arr[i-1] + grid*rat^k`, `meshgen.py:124-133`). Restarting that accumulation from
a rank-local origin changes the result in the last bits of a cumulative product; at
|x| ~ 1e4-1e5 m a relative 1e-12 is 1e-8..1e-7 m — **one to two orders above
`align`'s `coord_tol=1e-9`**, so the e2e cell would fail with "canonicalised rows do not
describe the same nodes" and nothing would say the cause was a re-accumulated line. Slice,
never re-accumulate. U2 gates it bitwise.

Cost of keeping the global lines: three arrays of O(n^{1/3}) doubles (`nx+ny+nz`, order
1e3 entries) on every rank. Negligible against the 1.65 GB this change is targeting.

### 3.3 Global↔local identity, where it is still needed

Only three places need it, and all three are satisfied by the analytic map — no stored
global tables:

```
g_node(ix,iy,iz) = f(off_x+ix, off_y+iy, off_z+iz)   # f = the serial numbering's own nest
```

1. **U3** (test-only): to compare against the serial mesh.
2. **Fault ownership** (§3.4): expressed in global `(ifx,ifz)` fault-plane indices, which
   are O(nftnd) ints — thousands, not millions.
3. **Halo pairing** (§4.2): replaced by a coordinate handshake, not by global ids.

Explicitly dying: `decompose`'s `halo_eq_global` (`:224-226`, "the only form in which two
ranks can agree on a shared equation"). Under rank-local numbering there is no global
equation id. Its replacement is §4.2 — and it must be a REPLACEMENT, not a deletion: the
invariant it enforced is the one thing standing between a mispaired halo and a silently
wrong answer.

### 3.4 Fault-node ownership under a structured split

Today: "lowest rank whose element slab touches both `nsmp1` and `nsmp2`"
(`:309-324`), computed from the global `conn`. Under rank-local: a fault node has global
fault-plane indices `(ifx, ifz)`; the set of ranks whose index boxes contain it follows
from the three 1D partition tables (which every rank can hold — they are O(P)); the owner
is the lowest such `rank` by `calcXyzMPIId` ordering. Every rank computes the same answer
from the same table with no communication, as today. Gated by §1.3's two allreduces
(sum == nftnd, and per-node ownership count == 1).

**A 3D split changes which ranks own zero fault nodes.** With `(2,2,1)` for `test.tpv8` at
4 ranks, the fault is the `y=0` plane and the split is in x and y, so the ranks on the
far-y side plausibly own none — where the 1D element cut gave all four ranks fault nodes
(matrix.py:98-103 records the measurement: 132/829/806/124, 4 of 4 files).
`matrix.PY_MPI_EXPECTED_FRT_FILES[('test.tpv8',4)] = 4` is DATA and will likely have to
become 2. **Re-measure it and change it deliberately; do not predict it in the table.**
This is the item-43 lesson restated: the empty-fault-rank shape needed 4 ranks to appear
under a slab, and under a 3D split it appears at different ranks — so the smoke test must
run at `(2,2,1)` AND `(2,2,2)`, not at one of them.

---

## 4. The halo under a 3D split

### 4.1 6 neighbours, not 26 — and why that is correct

Read `assembleGlobalMass.f90:94-168`: Fortran loops `do ixyz = 1,3` and exchanges only the
two FACE neighbours in each direction, sending the whole face plane
(`abc(ixyz) = numxyz(j)*numxyz(k)` nodes) plus the fault boundary terms
(`addFaultBoundaryTerm`, `:172`). Edge- and corner-shared nodes (4-way and 8-way sharing)
are never exchanged directly: the x-phase result is written back into `quantArray` before
the y-phase reads it, so an edge contribution reaches its other sharers by **relay**.
6 exchanges, not 26, and the correctness argument is the ORDERING, not the neighbour set.

Two options for the port:

- **(a) Match Fortran: three sequential phases, 6 Sendrecv.** Message count fixed at 6
  regardless of rank count; volume = the three face areas. Requires three device→host→device
  round trips per step instead of one, because phase k+1's send buffer depends on phase k's
  received values. On the current numbers (halo stage ~50 ms/step at 8-16 ranks under the
  1D slab) tripling the host round trips is a real risk, not a detail.
- **(b) One round over the full sharing set (up to 26 neighbours), keeping the current
  `exchange` shape.** Each rank sends its partial at each shared equation to every rank
  sharing it and adds everything received; one round, one host round trip, correct provided
  the sharing sets are exact. Volume is slightly higher (edges/corners sent to more peers),
  message count higher, host round trips 1 instead of 3.

**Recommendation: (b)**, because it preserves `exchange`'s current contract
(`MPI4NodalQuant.py:482-502`) — one gather, one host hop, one add — and because the halo
stage's measured cost is dominated by something other than transport (0.116 ms of a 64.74
ms step at 16 ranks, `:493-495`). But (a) is the fallback if the 26-way pairing wave (§4.3)
measures badly, and the choice must be MADE BY MEASUREMENT at Stage 4, not now.

Either way the sharing sets are **analytic** under a structured split (the shared set with
the +x neighbour is this rank's last x-plane, etc.), which is what makes them affordable
once the global `conn` is gone.

### 4.2 Replacing the `halo_eq_global` symmetry invariant

Once per run, at setup, per neighbour: exchange the `(x,y,z)` coordinates and per-node dof
counts of the shared plane, in the order the exchange will use, and assert **equal length
and bitwise equality**. Cost O(surface) once; it converts "both sides enumerate the plane
in the same order" from a claim into a checked precondition. Without it a transposed
`(iy,iz)` loop on one side produces a silently wrong sum — the `conn`-free analogue of
exactly the bug `_relabel` currently makes impossible.

### 4.3 Does a 3D neighbour graph break the "no Isend/Irecv" decision?

The decline is recorded as measured at 4, 8, 12, 16 and 32 ranks under the slab
(`:493-495`). **It does not transfer, and it must be re-measured.** Two separate points:

- **Deadlock safety still holds.** The brief's worry — that ascending `Sendrecv` is safe
  only because a slab gives two neighbours — is not the binding reason. The ordering
  argument is general: with every rank walking its neighbour list in ascending rank order,
  take the smallest still-incomplete pair `(a,b)`, `a<b`, under the `(min,max)` ordering.
  All of `a`'s pairs with partners `<b` have the same min and a smaller max, so they are
  complete; hence `a` is posting to `b`. All of `b`'s pairs with partners `<a` have min
  `<a`, so they are complete; and `b` prefers `a` over any partner in `(a,b)`; hence `b` is
  posting to `a`. The pair matches and progress is made. By induction the exchange
  completes on ANY neighbour graph. So 26 neighbours do not introduce a deadlock.
- **Performance does not follow.** That same argument shows the exchange can SERIALISE
  into a wave — a rank waits for a busy neighbour before its later pairs can match. With 2
  neighbours the wave is 2 deep; with 26 it can be much deeper, and 26 sequential Sendrecvs
  each pay a latency. The 0.116 ms figure has no predictive value here. **Stage 4 must
  re-run the halo microbenchmark at (2,2,2), (4,2,2), (4,4,2) before the non-blocking
  question is called either way**, and the answer must be attributed (transport vs wait),
  because the ~50 ms halo residual is currently unattributed and a wave would look exactly
  like it.

---

## 5. Staging — where the green gates are

Each stage is independently gated, and the first gate is at the end of **Stage 1**, not at
the end of the plan.

| stage | what lands | independent gate | parity gate |
|---|---|---|---|
| **1. Partition arithmetic + sliced lines** | `one_dim_coor_array`/`build_grid_lines` take `(m,P)`; nothing else calls them with `P>1` yet | **U1 + U2** (`testsys/run.py unit`, ms) | serial path byte-unchanged → the whole existing 30-cell sweep is the parity gate, and it must be **bit-identical**, not within bound (no solver line changed) |
| **2. Rank-local mesh build, no MPI** | `build_solver_state(decomp, rank)` builds one rank's mesh; single-process, loop over ranks in one test process | **U3** — bitwise identity against the serial mesh at `(2,2,1)` and `(2,2,2)` on `test.tpv8` | same: serial sweep unchanged; U3 IS the numbering parity gate, and it is the one that replaces the lost construction guarantee (§1.2) |
| **3. Memory/time measurement, before any solve** | nothing lands; measure per-rank setup peak RSS and wall time at 1/4/8/32 under rank-local vs global | **the kill measurement of §6.1** — if the mesh-proportional fraction of the 1.65 GB is small, stop here and report | none (no solver path) |
| **4. Analytic halo + 3D `decompose`** | neighbour/sharing sets, coordinate handshake (§4.2), ownership (§3.4); `exchange` unchanged or (a)/(b) chosen by measurement | R1 reductions + the §4.2 handshake, on a 2-step run at `(2,2,1)` and `(2,2,2)` | **e2e `test.tpv8` python-jax-mpi at 4 ranks, at the case bound, against the unchanged committed reference** — the existing cell, now on a 3D split |
| **5. Rank counts and the empty-fault shape** | opt-in table + expected-file counts re-measured | re-measured `PY_MPI_EXPECTED_FRT_FILES`, committed as data with the measurement quoted | same cell at 8 `(2,2,2)` — the rank count at which the zero-fault-node rank moves |
| **6. Scaling** | `run_scaling.py` reads the same `DECOMP` for both backends | per-rank solve-time differencing only (never wrapper wall clock — the ≥4-rank trap) | the Fortran-vs-python table of §0, re-taken at 1-32 |

Stage 3 is a stage on purpose: it is the cheapest point at which the whole thing can be
abandoned having spent two days instead of two weeks.

### 5.1 Interaction with the in-flight TPV30 port-divergence hunt

Stages 1-3 touch `meshgen.py`, `eqdyna3d.build_solver_state`, and test files only —
**no `faulting.py`, no `driver.py`, no `fric.py`.** They rebase cleanly over that fix.

Stage 4 edits `driver.run_mpi` at `:366-404` (mass check, `mass_l`, `decompose` call site)
and consumes `faulting.build`'s dict through `_restrict_fault`, which keys on leading-axis
length `n == nftnd` (`MPI4NodalQuant.py:406-422`). A TPV30 fix that **adds a per-fault-node
array** to `faulting.build` is absorbed silently and correctly; one that changes a fault
array's **leading axis** (e.g. per-node → per-node-per-something) breaks the restriction
silently in the old code and must be re-checked once the fix lands. That is the one
coupling worth naming; everything else is textual.

---

## 6. What would make this NOT worth building

### 6.1 The memory return may be roughly half what is quoted — measure first (Stage 3)

The settled figure is **1.65 GB per rank, 52.8 GB summed at 32**. That is the peak of
`build_solver_state`, and only its **mesh-proportional** part is recoverable. Not
mesh-proportional: the input-file arrays, `on_fault_vars_input.nc`, station tables,
material tables, and whatever the interpreter and numpy themselves hold. Nothing in the
tree tells me that split, and I did not measure it for this document. If, say, 40% of the
1.65 GB is rank-invariant, the 32-rank return is ~32 GB, not ~51 GB, and the per-rank floor
(which is what actually caps rank count on a smaller machine) falls to ~0.7 GB, not ~0.05
GB. **Stage 3 measures that split before Stage 4 is written.** Kill criterion: if the
mesh-proportional fraction is below ~50%, the memory prize is not worth the numbering risk
of §1.2 class C and the plan stops with the 3D split only (which does not require
rank-local mesh — §6.4).

### 6.2 The time return is bounded by 7.92 s, not 10 s, and probably well under

The 10 s figure comes from the serial 4.42 s; the honest ceiling is the **concurrent**
7.92 s, since that is what a run actually pays at 32-way. The recoverable part is again
only the mesh-proportional fraction, and the concurrency slowdown (1.79x) is itself
evidence of memory-bandwidth contention that shrinks as the per-rank build shrinks — which
helps, but is second-order. Against a `test.tpv104` run of hundreds of steps at ~130
ms/step, a few seconds of setup is not a reason to do this. **The memory is the only prize;
if §6.1 comes back small, there is no prize.**

### 6.3 The owner's explicit stop condition: global arrays becoming untenable

Checked, and **not currently triggered.** After the change the surviving global-extent
objects are: the three 1D coordinate lines (O(1e3) doubles), the 1D partition tables
(O(P)), and the fault-plane ownership indices (O(nftnd) ~ 1e3-5e3). All are
O(n^{1/3}) or O(fault area), none is O(N). The two O(nranks × N_global) loops in
`decompose` (`:296-308`, `:312-316`) disappear rather than needing a global array.

The condition WOULD trigger on a case with a **per-node 3D input** (a 3D material volume,
a per-node initial-stress field) that cannot be evaluated pointwise from `(x,y,z)`. No
gated case has one today. If one appears mid-implementation: **report and stop; do not
trade away the remaining global identity silently.**

### 6.4 The two changes are separable, and that is a live option

The 3D split does **not** require rank-local mesh: it can be done on the current
"build serial, restrict" design by replacing the contiguous element cut with a structured
`(ix,iy,iz)`-box selection, keeping global numbering, keeping the construction guarantee
of §1.2 intact, and keeping `halo_eq_global`. That version buys the whole halo-surface
argument (the thing the §0 measurement actually motivates) and none of the memory. If
Stage 3 kills the memory return, **this is what should land instead** — and it is a
smaller, safer change than the combined one. The combined landing is the owner's call and
I am not reopening it; I am recording that the fallback exists and is cheap, so that a bad
Stage-3 number has somewhere to go other than "abandon item 43".

### 6.5 The residual that three refutations have already survived

Per-step cost on this path sits 103-126 ms above T1/N at every rank count ≥ 4, and
transport, placement and element balance have each been measured and each failed to
explain it (`MPI4NodalQuant.step_profile`, `:458-479`; item 60). **A 3D split is a fourth
candidate explanation for the same residual**, and it has the same shape as the three that
failed: a real defect (a slab's surface does not shrink) whose repair may not move the
outcome curve. Rule 4a applies — falsify against the OUTCOME curve. The honest prediction
to write down BEFORE Stage 6 runs: if 3D is the cause, the 8- and 16-rank ratios should
fall from 1.56x/1.58x toward the ~1.07x of 2-4 ranks. **If they land near 1.4x, the halo
surface was not the binding constraint and the residual is still unattributed** — that is a
reportable negative result, not a reason to add a fifth mechanism.

---

## 7. Open questions a reviewer would attack

1. **Class C (§1.2) is mitigated, not closed.** U3 covers the rank counts where a serial
   mesh fits; R1 covers conservation everywhere. Neither proves correct numbering at 32
   ranks on `test.tpv104`. The remaining argument is that Fortran has lived with exactly
   this exposure since 2006 and the references it produced are the ones everything is gated
   against — which is an appeal to precedent, and a reviewer is entitled to say so.
2. **Choosing (a) vs (b) in §4.1 is deferred to a measurement that has not been made.** If
   (b) measures badly and (a) triples the host round trips, there may be no good option and
   the halo stage stays where it is.
3. **`PML_WEIGHT` removal is justified by item 60's null result**, i.e. by a measurement
   that the knob did not matter — not by a measurement that a geometric split is balanced.
   A structured split on a mesh with a thick PML on five faces will NOT be work-balanced,
   and the imbalance shows up as halo wait time, which is exactly the signal §6.5 is trying
   to attribute. `report` must keep printing per-rank `(Ei, Ep)` so the imbalance stays
   visible in the measurement instead of hiding inside it.
4. **No number in this document was measured by me.** Everything quoted is from the
   session log, `matrix.py`'s comments, or module docstrings, and is labelled as such.
