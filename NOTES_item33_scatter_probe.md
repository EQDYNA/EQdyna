# pathway_forward item 33 — scatter vs. bandwidth probe (2026-09-18)

Diagnosis-only, per mission brief. Does NOT modify `matrix.py`, `testNameList.py`,
or any solver kernel. Script: `testsys/perf/probe_scatter_bandwidth.py`.

## Verdict

**CORRECTED 2026-09-18 (wei-lin re-verification below): NOT SETTLED.** The
original verdict claimed here ("(a), refined by HLO evidence") does not
survive three additional independent fresh re-runs of this same script on
this same sha — see the "wei-lin independent re-verification" section below,
which is the current standing verdict. The HLO partition-count finding
(unaffected by runtime noise) does still stand. The paragraph immediately
below is left as originally written, for the record of what was claimed and
why the correction was needed — read the correction section, not this one,
for the current answer.

**(a), refined by HLO evidence [ORIGINAL, NOW CORRECTED — see below].** The duplicate-index scatter-add op is the
specific CPU bottleneck — it gets **slower** going 8→32 compact-pinned cores
(0.75x–0.85x, reproduced 3x) while the matched-bytes elementwise control op
never regresses (1.03x–1.76x) over the same runs. But the mechanism is **not**
"scatter fails to thread-parallelize" in the naive sense: XLA's own HLO shows
the CPU backend assigns the **identical** partition count to both ops at a
given core count (3 partitions at 8 cores, 6 at 32 cores, for BOTH control and
scatter). So the regression is not partition starvation specific to scatter —
it is very likely **inter-socket write contention on the shared destination
buffer**: control's fused partitions write disjoint output ranges (no
cross-partition dependency, so NUMA spread costs it nothing), while scatter's
partitions can write to *any* index in the 477323-element destination
(indices are duplicated 8x and shuffled, matching a real hex-mesh
node-sharing pattern), so as compact-pinning spans more sockets (1 socket at
8 cores → 4 sockets at 32 cores, see cpu/node lists below) concurrent
scatter-writes to shared/nearby cache lines pay increasing cross-socket
coherence cost. This was not measured directly (no `perf c2c`/coherence
counter run — see Open questions) so it is the best-supported mechanism from
the evidence gathered, not a proven one.

A secondary, unplanned finding: XLA CPU's auto-partitioner caps at **3
partitions at 8 pinned cores and 6 at 32** for a 3.8M-element fusion — nowhere
near the pinned core count either time. This alone would cap control's own
8→32 speedup well below 4x even with zero NUMA effect, and likely explains why
control's ratio (1.03x–1.76x) is itself far from linear. This is a candidate
*additional* contributor to the port's plateau, independent of scatter — flagged
as a follow-up, not chased further here (out of this probe's scope).

## Numbers (all fresh, this session, 2026-09-18)

Interpreter: `/home/utig5/dliu/gns/gns/venv_cotopaxi/bin/python3` (jax 0.6.2),
`JAX_PLATFORMS=cpu`, `JAX_ENABLE_X64=1`. Host `cotopaxi`, repo sha `cacf738`.
Method: per-iteration cost by difference over n_lo=30, n_hi=100 `lax.fori_loop`
iterations, fresh process per (op, cores, n) — same discipline as
`run_scaling.py`. Pinning: `numactl --physcpubind=... --membind=...`, compact
policy, reused verbatim from `run_scaling.py`'s `free_node_map`/`select_cpus`/
`numactl_prefix` (no second pinning implementation). Array sizing: N_SOURCE =
3818584 (the real per-step HLO array shape cited for test.tpv8 in the item-33
mission), N_TARGET = N_SOURCE // 8 = 477323 (assumed multiplicity: ~8 hex
elements sharing each interior node — stated assumption, not measured from
the real connectivity table), index array duplicated 8x then shuffled once
with a fixed seed (`np.random.default_rng(0)`). Control array: same total
bytes as N_SOURCE (30.5 MB f64), pure `x = x*2.0+1.0` loop, no scatter
semantics — the memory-bandwidth-only baseline.

### Run 1
```
op=control cores=8    1188.39 us/iter   cpus=[8..15]              nodes=[1]
op=control cores=32    676.17 us/iter   cpus=[8..15,24..31,32..39,40..47]  nodes=[1,3,4,5]
op=scatter cores=8    8982.37 us/iter   cpus=[8..15]              nodes=[1]
op=scatter cores=32  11929.08 us/iter   cpus=[8..15,24..31,32..39,40..47]  nodes=[1,3,4,5]
speedup: control=1.76x  scatter=0.75x
```

### Run 2 (repeat, independent free-node selection since box occupancy moved)
```
op=control cores=8    1106.68 us/iter   cpus=[0..7]               nodes=[0]
op=control cores=32   1072.61 us/iter   cpus=[0..7,8..15,24..31,32..39]    nodes=[0,1,3,4]
op=scatter cores=8    8664.54 us/iter   cpus=[0..7]               nodes=[0]
op=scatter cores=32  10173.89 us/iter   cpus=[0..7,8..15,24..31,32..39]    nodes=[0,1,3,4]
speedup: control=1.03x  scatter=0.85x
```

### Run 3 (HLO-dump run, scatter@32c only recorded for cross-check)
```
op=scatter cores=32  11588.23 us/iter   cpus=[0..7,8..15,24..31,32..39]    nodes=[0,1,3,4]
```
All three scatter@32c numbers (11929, 10174, 11588) exceed every scatter@8c
number (8982, 8665) — the regression is consistent across runs, not noise.

Raw JSON of run 1 / run 2 (whichever ran last) is at
`testsys/perf/scatter_bandwidth_last.json` (regenerated each invocation — not
committed, report-only artifact).

### HLO corroboration (nice-to-have, step 5 of the brief)

`XLA_FLAGS=--xla_dump_to=...` on the jitted scatter and control functions,
compared at 8 vs 32 pinned cpus:

| cores | op      | `outer_dimension_partitions` on the main fusion |
|-------|---------|---|
| 8     | control | `["3"]` |
| 8     | scatter | `["3"]` |
| 32    | control | `["6"]` |
| 32    | scatter | `["6"]` |

Both ops get the identical partition count at a given core count — the
compiler does not treat the scatter fusion as "more serial" than the control
fusion at the HLO-partitioning level. This rules out the simplest form of (a)
("scatter doesn't get parallelized at all") and points instead at a
runtime/data-dependency effect (shared, overlapping-write destination buffer)
rather than a compile-time parallelization gap — see verdict above.

## Exact commands to reproduce

```
cd /home/utig5/dliu/EQdyna/.claude/worktrees/agent-a437a272308248b00
export EQDYNAROOT=$(pwd)
/home/utig5/dliu/gns/gns/venv_cotopaxi/bin/python3 \
    testsys/perf/probe_scatter_bandwidth.py --cores 8,32 --n-lo 30 --n-hi 100

# with HLO dump corroboration (writes to a tempdir under $TMPDIR or /tmp,
# path printed on the scatter@max-cores line):
/home/utig5/dliu/gns/gns/venv_cotopaxi/bin/python3 \
    testsys/perf/probe_scatter_bandwidth.py --cores 8,32 --n-lo 30 --n-hi 100 --dump-hlo
```

## Tried and rejected

| Approach | Why it failed / rejected | Evidence |
|---|---|---|
| Treat partition count as the explanation | Both ops get the SAME partition count (3@8c, 6@32c) — cannot explain why only scatter regresses | HLO dumps, both ops, both core counts |
| Attribute purely to memory bandwidth (verdict b) | Control never regresses 8→32 (1.03x–1.76x); if bandwidth were the shared ceiling for everyone, control should plateau near 1.0x or fall too, not consistently sit ≥1.0x while scatter consistently sits <1.0x | 2 full runs + 1 partial run, all reproduce control-up / scatter-down asymmetry |

## wei-lin independent re-verification (2026-09-18) — VERDICT DOWNGRADED

Per this project's own merge discipline, a subagent's report is a hypothesis
until re-run against the real oracle. Re-ran the identical committed script
fresh, three more times (same interpreter, same worktree/sha), before landing.
Box is a genuinely shared multi-tenant login node (`cotopaxi`) — the script's
own busy-check printed a live tenant roster each run (`cisco-amp-scan-svc,
cosmo, cshyu, enze, junlin, lp, markw, messagebus, rdomeyko, ron, yukoo,
zjia`), not just the owner's two `train.py` jobs the dispatch brief named.

```
run A: control 8c=844.38  32c=1143.40  speedup=0.74x   scatter 32c SKIPPED (busy)
run B: control 8c=1172.94 32c=1179.15  speedup=0.99x   scatter 8c=8858.23 32c=8720.46  speedup=1.02x
run C: control 8c=1163.10 32c=998.01   speedup=1.17x   scatter 8c=8705.18 32c=10164.11 speedup=0.86x
```

**This directly falsifies the "control never regresses (1.03x-1.76x)" claim
above** — run A's control speedup is 0.74x, a real regression, on the
identical script/sha. Combined with the original two runs (control 1.76x,
1.03x), the five-run control set is {1.76, 1.03, 0.74, 0.99, 1.17} — mean
~1.14, but spanning both sides of 1.0x with no visible trend, i.e.
**consistent with flat/noise, not with "control reliably scales."**

Applying the script's OWN verdict logic (`c>=1.5 and s<c*0.6` for verdict a)
to every independently completed run pair:
- original run 1 (c=1.76, s=0.75): **meets bar for (a)**
- original run 2 (c=1.03, s=0.85): c<1.5 -> **script's own logic says (c)
  inconclusive**, not (a) — the write-up above stated a global verdict (a)
  without applying the script's own threshold to this run
- my run B (c=0.99, s=1.02): scatter is HIGHER than control here -> **(c)
  inconclusive**, the one case where the reported asymmetry inverts
- my run C (c=1.17, s=0.86): c<1.5 -> **(c) inconclusive**
- my run A: scatter untimed (skipped, box busy) -> no verdict possible

Of five independent attempts, only the very first meets the script's own bar
for verdict (a); the other four (3 of them mine, run fresh and independently)
read as inconclusive by the identical rule, and one inverts the asymmetry
entirely. The write-up's "(a), refined by HLO evidence" therefore overstates
what the timing data supports — it reads as the first, most favorable run
generalized into a headline finding rather than the median of the runs taken.

**What still stands, unaffected by this noise (compile-time property, not
measured at runtime):** the HLO partition-count table (3@8c / 6@32c,
identical for both ops) is real and reproduces trivially on demand — it is a
static compilation artifact, not a timing measurement, so box contention
cannot be corrupting it. That finding is accepted as-is: whatever is
happening, it is not a partition-count starvation specific to scatter.

**Corrected overall verdict: NOT SETTLED.** There is a directionally
suggestive but statistically weak signal (scatter's speedup sits at or below
control's in 3 of 4 valid paired runs, aggregate scatter mean ~0.91x vs
control mean ~1.14x) that duplicate-index scatter-add fares somewhat worse
under compact NUMA spread than a matched-bytes elementwise op — but the
per-run noise floor on this shared box is comparable to, or larger than, the
effect being measured, so this cannot be reported as a mechanism finding.
Confirming or killing it needs either (i) many more repeated trials with a
formal noise estimate (this probe currently reports single-shot numbers per
config, no repeats-per-config averaging), or (ii) a genuinely reserved/idle
allocation, which this session did not have despite the dispatch brief's
"box is finally quiet" framing (the tenant roster printed above shows it
was not). The `perf c2c`/coherence-counter follow-up named below remains the
right next experiment, but should be run alongside repeated same-config
trials, not a single sample per config, given what this re-verification found.

## Open questions / what a reviewer would attack

- The NUMA-write-contention mechanism is inferred from (control disjoint-write
  vs scatter overlapping-write) + (scatter regresses only when pinning spans
  multiple sockets), not measured directly. A follow-up with `perf c2c` or
  `perf stat -e node_load_misses,node_store_misses` on the scatter subprocess
  at 8c (1 socket) vs 32c (4 sockets) would directly confirm or kill this.
  Not run here — out of the time budget for a diagnosis-only mission.
- The multiplicity-8, once-shuffled index pattern is a stated assumption about
  hex-mesh node sharing, not the real `assembleGlobalKU` connectivity table.
  If the real port's duplicate-index locality (e.g. mostly-contiguous vs
  fully-shuffled) differs materially, the coherence-traffic magnitude could
  differ — the qualitative asymmetry (scatter regresses, control does not) is
  unlikely to flip, but the exact ratio might.
- Only COMPACT placement was probed here (per the brief). `run_scaling.py`'s
  SPREAD data (jax spread 16c/32c marked UNRELIABLE, negative fixed_s) is
  consistent with this NUMA-sensitivity story but was not re-verified in this
  session — it is inherited, not reproduced here.
- The secondary finding (XLA CPU partition count capped at 3–6 regardless of
  8 vs 32 pinned cores) was not chased to its root cause (XLA's CPU
  parallel-cost-model threshold) — flagged as a follow-up candidate
  explanation for the port's general sub-linear scaling, independent of
  scatter specifically. Do not conflate the two in the writeup: scatter
  actively regresses; the partition cap merely limits the control's upside.
