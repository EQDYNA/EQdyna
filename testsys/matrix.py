#! /usr/bin/env python3
"""
The e2e sweep's TABLE: which (case, backend) cells exist, which are supported,
and the ONE bound each case is gated at.

There is one test in this repo -- the e2e sweep -- and it is two nested loops:

    for case in CASES:
        for backend in BACKENDS:
            run(case, backend) -> canonical frt -> compare against the ONE
                                  committed reference for that case, at THAT
                                  CASE's bound

`backend` is an axis of the sweep, not a tier: it is a column of this table,
not a separate runner. Adding a GPU backend is appending to BACKENDS.

WHY THE TABLE IS DATA, NOT CONTROL FLOW
A cell is either SUPPORTED -- it runs, and it is gated -- or it is DECLARED
UNSUPPORTED here with a recorded reason. There is no third state, no skip, and
no silent absence. PROJECT_RULES rule 2 requires "I could not check this" and
"this is fine" to produce different exit codes, so asking the sweep for a cell
this table declares unsupported is a hard FAILURE (see cells()).

WHY THE BOUND BELONGS TO THE CASE, NOT THE BACKEND
A per-backend bound makes the backends incomparable: "jax passed and numpy
failed" would not be a statement about the solver, it would be a statement
about which number each column was allowed. One bound per case, applied
identically to every backend, is what makes a row of this sweep quantitative.

The bounds below are UNCHANGED from the accept tier they came from (rule 5:
one calibrated definition of pass; and no bound is ever moved to make a cell
go green). Each is 15-120x the roundoff observed for that case on the jax
backend, and every one is at or below THRESHOLD, the repo's single outer
sanity bound.
"""
import os
import sys

TESTSYS = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTSYS)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from testNameList import nameList as _NAME_LIST, coreNumList as _CORE_NUM_LIST

# The case list is testNameList.py's, not a second copy of it (rule 1). A case
# added there with no entry below is a hard import-time failure, not a silently
# ungated case -- see the consistency block at the bottom of this module.
CASES = tuple(_NAME_LIST)
FORTRAN_RANKS = dict(zip(_NAME_LIST, _CORE_NUM_LIST))

BACKENDS = ('fortran', 'python-jax', 'python-jax-mpi')
# python-numpy LEFT the backend axis entirely (2026-09-23 owner decision:
# "I actually don't care numpy... I will use Jax anyway"). It is not a
# UNSUPPORTED cell and not a third cell state -- it is not a column of this
# table at all, the same way there is no 'python-torch' column. The numpy
# CODE (src/python/eqdyna, `--backend numpy`) is untouched and still runnable
# by hand; only the gate stopped exercising it. See git history for the
# retired numpy cells (bounds, RSS measurements, RELEASE_ONLY entries) --
# they are evidence of what the axis USED to cover, not of what it covers now.

# Which artifacts each backend produces, and therefore what gets compared.
#   frt -- the canonical fault-node table (testsys/frt_canonical.py). Every
#          backend produces it; it is what makes the sweep ONE comparison.
#   nc  -- fault.dyna.r.nc, written by scripts/plotRuptureDynamics from the frt
#          output. DELIBERATELY still compared, for the fortran column only:
#          it is strictly lossy relative to frt (12 resampled dip x strike
#          variables vs 18 physics columns at native node resolution), so it
#          adds no physics coverage -- but it is the ONLY thing that exercises
#          plotRuptureDynamics, post-processing that would otherwise be gated
#          by nothing. It is backend DATA here, so a green run states that the
#          python columns did not cover it.
#   nsign -- the SIGN of on-fault station column 8 (n-stress) at the first
#          recorded step, at every buried on-fault station, against the case's
#          SCEC spec convention (NSTRESS_CONVENTION below; board row 22a).
#          Fortran only until the port writes station files (board row 114).
ARTIFACTS = {
    'fortran': ('frt', 'nc', 'nsign'),
    'python-jax': ('frt',),
    'python-jax-mpi': ('frt',),
}

# NSTRESS_CONVENTION -- per case, the sign convention its SCEC spec states for
# the on-fault station n-stress column, and where it says so (board row 22a;
# spec texts fetched into scratch/specs/, 2026-09-24). The specs DISAGREE, so
# this is per case, never a constant: until row 22a library_output.f90 negated
# column 8 for every case (compression-positive everywhere), which was right
# only for TPV104 and TPV105-3D.
#
# This table is INDEPENDENT of each case's own par.faultStNormalStressSign
# (case_input/<case>/, default in scripts/defaultParameters.py):
# testsys/regression/test_nstress_sign_convention.py checks the declaration
# against this table, and compare.nstress_sign_gate checks the written files
# against it. A case with no SCEC spec (or a spec that asks for no n-stress
# field) takes the SCEC default, 'extension', and says so.
NSTRESS_CONVENTION = {
    'test.tpv8': ('extension', 'uploadTPV89 requests no n-stress field; SCEC '
                  'default (Signconvention3d: normal slip > 0 for extension)'),
    'test.tpv10': ('extension', 'uploadTPV10_11_v3, n-stress: "Positive means '
                   'extension."'),
    'test.tpv104': ('compression', 'uploadTPV103 ("Uploading Data for TPV103 '
                    'and TPV104", 2008-10-18), n-stress: "Positive means '
                    'compression."'),
    'test.tpv1053d': ('compression', 'TPV105_3D_formats:114-115, n-stress: '
                      '"Positive means compression."'),
    'test.tpv29': ('extension', 'TPV29_30_Description_v06:1174-1177, n-stress: '
                   '"Positive means extension."'),
    'test.tpv30': ('extension', 'TPV29_30_Description_v06:1174-1177 (TPV29 and '
                   'TPV30 share the format), "Positive means extension."'),
    'test.tpv36': ('extension', 'TPV36_37_Description_v12 (2024-07-15), '
                   'n-stress: "Positive means extension."'),
    'test.tpv37': ('extension', 'TPV36_37_Description_v12 (2024-07-15), '
                   'n-stress: "Positive means extension."'),
    'test.meng2023a': ('extension', 'not a SCEC benchmark; SCEC default'),
    'test.meng2023cb': ('extension', 'not a SCEC benchmark; SCEC default'),
    'test.drv.a6': ('extension', 'not a SCEC benchmark; SCEC default'),
}

# python-jax-mpi is a FOURTH value on the `backend` axis: real MPI (one
# process per rank, driver.run_mpi + MPI4NodalQuant.py owning the Fortran-
# style domain decomposition, jax owning the local element kernel), launched
# under mpirun rather than as a single serial process. It is NOT a flag on
# python-jax, NOT a new tier and NOT a new comparison -- it produces the same
# frt artifact, canonicalised and compared against the SAME one reference at
# the SAME CASE_BOUND every other backend uses (no per-backend bound, see the
# module docstring's WHY above).
#
# PY_MPI_RANKS: which cases opt into this backend, and at how many ranks. A
# case absent here has NO python-jax-mpi cell -- it is DECLARED UNSUPPORTED
# below, not silently missing. Opting more cases in is the owner's suite-cost
# call (see the UNSUPPORTED reason string), not a policy requirement.
PY_MPI_RANKS = {
    'test.tpv8': 4,
}

# PY_MPI_EXPECTED_FRT_FILES: the number of `frt.txt<rank>` files a cell must
# produce, keyed by (case, ranks) -- DATA, not `== ranks` assumed in code.
# A rank writes a frt file only if it OWNS a fault node, i.e. is the lowest
# rank holding it (MPI4NodalQuant.owned_mask): a box that never meets the
# fault writes none, and neither does a box whose fault nodes all sit on a
# plane it shares with a lower rank. So this number CAN be less than the rank
# count, and a launch that silently started fewer workers than asked must not
# be indistinguishable from one that started the expected count.
#
# Which ranks write depends on the DECOMPOSITION, and it changed with item 64
# (2026-09-24): the python-jax-mpi split is now Fortran's (npx,npy,npz) =
# MPI4NodalQuant.DECOMP[ranks], rank-local boxes, replacing the 1D element
# slab under which all 4 tpv8 ranks owned fault nodes (132/829/806/124).
# Re-measured, not predicted: `mpirun -np 4 python3 -m eqdyna <serial tpv8
# case> --backend jax --mpi` on the item-64 branch logged decomposition
# (2, 2, 1), fault computed 961/0/961/0 and owned 961/0/930/0 -- the y split
# of tpv8's asymmetric y line (-10/+12 km) falls off the fault plane, so the
# two mey=1 boxes never meet the fault -- and wrote frt.txt0 and frt.txt2,
# 1891 rows pre-dedup (= the reference's row count).
PY_MPI_EXPECTED_FRT_FILES = {
    ('test.tpv8', 4): 2,
}

THRESHOLD = 1e-3  # PROJECT_RULES rule 5's one outer sanity bound.

# --- the ONE term (2026-09-23 owner decision, superseding the two-term design
# a same-day earlier change had landed) -------------------------------------
# Every gated case in this sweep runs at GATE_TERM_S, period: there is no
# second term, no `--term` flag, and no per-case "full" term to fall back to.
# run_e2e.py applies it through the ONE existing override path
# (apply_term_override), unconditionally, for every cell it runs -- the
# everyday sweep, the release sweep (the same cells, plus evidence) and
# CI's smoke selection all run the same 5 s.
#
# case_input/<case>/user_defined_params.py's own committed `par.term` (what
# `create.newcase` hands a user who builds the case by hand, outside the
# gate) is DELIBERATELY left alone by this change and by run_e2e.py -- it is
# the compset's own default for standalone use, not a second source of truth
# the gate reads or reconciles against. The gate always overrides it.
GATE_TERM_S = 5.0

# One bound per case. None means "this case is not gated on a scalar" -- see
# GATE/DRV_A6 below. Every case in CASES must appear here.
CASE_BOUND = {
    'test.tpv8': 1e-8,        # observed 8.23e-11 (python-jax vs reference)
    'test.tpv104': 1e-5,      # observed 3.39e-07
    'test.tpv1053d': 1e-4,    # observed 6.65e-06
    'test.tpv10': 1e-6,       # observed 3.02e-08 (dipping fault, insertFaultType==1)
    # TIGHTENED 2026-09-17 (item 37 -- same defect item 19a caught and fixed
    # for tpv36's bound, found here by inspection of this file's own
    # commented history a few lines below, not by a fresh claim): the
    # THRESHOLD these three carried was justified by "no python run has ever
    # been measured for these", which stopped being true once haruto's
    # CI_CELLS widening pass measured all three (`/usr/bin/time -v`,
    # full-length, one at a time, this dev box, 2026-09-16) AND confirmed
    # each one passes its cell via testsys.compare.compare_cell on that same
    # run. Observed: meng2023a numpy/jax 4.628301e-09/3.229082e-09,
    # meng2023cb numpy/jax 5.054473e-09/4.395842e-09, tpv29 numpy/jax
    # 9.876230e-15/1.276548e-13. Bounds set the same way tpv36's was: the
    # worst observation times test.tpv8's own headroom ratio
    # (1e-8 / 8.23e-11 = ~121.5x), then rounded to the nearest bound already
    # in use at that order of magnitude rather than an arbitrary new number.
    'test.meng2023a': 1e-6,    # worst 4.63e-09 * ~121.5 = 5.6e-7; rounds to
                               # the same 1e-6 already used for tpv10/tpv36.
    'test.meng2023cb': 1e-6,   # worst 5.05e-09 * ~121.5 = 6.1e-7; same bound.
    'test.tpv30': 1e-10,       # registered 2026-09-23 (owner gating decision), 5 s / 500 m:
                               # observed python-jax 1.909216e-14 vs the 5 s Fortran
                               # reference (51b7649), twice (15:34 at 95d4220, 22:27 at
                               # 51b7649); * ~121.5 = 2.3e-12, rounded UP to tpv29's
                               # 1e-10: two orders above that, ~5200x headroom (tpv29
                               # carries ~780x), kept because it is the near-epsilon
                               # bound already in use. Sensitivity: with
                               # e1888e7 (PML-node gravity) reverted the 5 s cell reads
                               # 1.202240e+08 -- caught by 18 orders of magnitude.
    'test.tpv29': 1e-10,       # worst 1.28e-13 * ~121.5 = 1.6e-11; rounded UP
                               # to 1e-10 (headroom ~780x, not the smaller
                               # ~121x -- these are near-machine-epsilon
                               # diffs and 1e-11 would leave no margin for
                               # ordinary cross-run noise at that scale).
    'test.drv.a6': None,      # chaotically bistable -- flip-budget gate, below.
    # C_degen>3 (wedge-degeneration), dip=15 (tpv36's par.dip), dx=500,
    # par.term=6 -- a COARSE regression gate (rule 17 step 4), NOT a claim
    # about physical accuracy at SCEC-standings resolution (the published
    # tpv36/37 standings rank a run this coarse last and improve
    # monotonically with resolution). Measured AFTER item 19(a)'s jax-only
    # regression fix (assembleGlobalMass.py's eleshp einsum dispatch), on a
    # reference frozen fresh for this commit -- the merge's original
    # measurement predates that fix and is not reused: observed python-jax
    # 5.867293e-09, python-numpy 7.894064e-09 (both near-zero peak-slip-rate
    # components on marginally-ruptured nodes -- roundoff, not physics).
    # 1e-6 is ~127x the worst observation, the same headroom ratio
    # test.tpv8 carries (1e-8 / 8.23e-11).
    'test.tpv36': 1e-6,
    # test.tpv37: identical wedge-degeneration machinery as tpv36 (same
    # dip=15, dx=500, par.term=6, C_degen=par.dip) -- only par.tpv=37 and one
    # initial-stress-taper coefficient differ from tpv36's compset (confirmed
    # by diff: 0.001875e6 vs tpv36's 0.0005e6 in on_fault_vars[iz,ix,4]).
    # Same coarse-regression-gate reasoning as tpv36 (rule 17 step 4, NOT a
    # SCEC-standings accuracy claim). Fortran cell bit-exact against a
    # freshly-frozen reference (max|diff|=0.0). Measured fresh, this commit:
    # python-jax 5.256410e-09, python-numpy 6.929140e-09 (both near-zero
    # peak-slip-rate-type components on marginally-ruptured nodes -- roundoff,
    # not physics, same shape as tpv36's own observation). 1e-6 is ~144x the
    # worst observation, in the same headroom family as tpv8/tpv36.
    'test.tpv37': 1e-6,
}
# test.tpv30 was held out of the gate from 2026-09-17 to 2026-09-23 by a real
# divergence (numpy==jax, both != Fortran by up to 4.0e8 Pa at t=20 s). It was
# root-caused and fixed in e1888e7: a PML node never received its own
# elements' gravity. It is registered above, at 5 s, on the owner's decision.

GATE = {
    'test.tpv8': 'abs-max',
    'test.tpv104': 'abs-max',
    'test.tpv1053d': 'abs-max',
    'test.tpv10': 'abs-max',
    'test.meng2023a': 'abs-max',
    'test.meng2023cb': 'abs-max',
    'test.tpv30': 'abs-max',
    'test.tpv29': 'abs-max',
    'test.tpv36': 'abs-max',
    'test.tpv37': 'abs-max',
    # test.drv.a6 (C_elastic==0 viscoplastic, friclaw==4, fractal-rough, long
    # duration) has genuinely bistable rupture arrivals, so a scalar max-abs
    # bound cannot distinguish "a few hundred marginal nodes flipped" from
    # "everything drifted a little". Gated on bulk agreement PLUS an explicit
    # flip budget instead.
    'test.drv.a6': 'flip-budget',
}

# test.drv.a6's flip-budget gate. These numbers are measured, not chosen;
# testsys/parity/evidence_drv_a6_chaos.py regenerates them and they must not be
# changed without re-running it fresh (rule 4).
DRV_A6 = {
    'fnft_col': 3,               # canonical frt column layout: x, y, z, fnft, ...
    'rupture_sentinel': 1.0e4,   # fnft sentinel is 99999.0 in Fortran
                                 # (eqdyna3d.f90:139) and in the port; below
                                 # 1e4 means "this node ruptured".
    'timing_shift_s': 1.0,       # a ruptured-in-both node counts as a flip only
                                 # past this |delta fnft|.
    'median_fnft_bound': 0.1,    # seconds; measured 0.0417 s in all three of
                                 # A/B/C -- 10x margin.
    'phys_max_bound': 1.1e9,     # 30x the measured 3.996e7 ceiling over the 18
                                 # non-coordinate/non-fnft columns, restricted
                                 # to matched-arrival nodes (those columns carry
                                 # ~1e8 Pa stress fields).
    'total_flip_bound': 450,     # A=372 (Fortran-only decomposition floor),
                                 # B=439 (python vs serial fortran). 450 is
                                 # minimal headroom above the larger measured
                                 # number, fixed BEFORE C (391) was measured.
                                 # A LATER, independent Fortran-vs-Fortran run
                                 # (serial vs the committed 4-rank reference,
                                 # 2026-09-15) measured 329 flips, max diff
                                 # 2.12e7, p50 4.03e-03, p99 5.04e+05, and a
                                 # median arrival difference of 0.0417 s =
                                 # exactly one time step. Recorded, not
                                 # substituted: A and this run are different
                                 # experiments and 450 stays where it was set.
                                 # That 329 is the floor imposed by
                                 # decomposition ALONE is what withdrew the
                                 # case -- python-jax 397 and python-numpy 461
                                 # straddle a budget whose own floor is 329.
}

# DECLARED UNSUPPORTED CELLS. A reason, verified in the source, not a guess.
# Absence from this table means "supported"; presence means the sweep refuses
# to pretend it covered the cell.
_MPI_OPT_IN_REASON = (
    'not opted into the optional MPI execution mode (no PY_MPI_RANKS entry). '
    'PROJECT_RULES rule 17 step 7 requires a CASE to be supported on every '
    'BACKEND IMPLEMENTATION; python-jax-mpi is an optional execution MODE of '
    'the python-jax backend, whose own cell is supported and gated here. '
    'Opting more cases in is a suite-cost decision for the owner, not a '
    'policy requirement -- it is not a 40-cell obligation.'
)
UNSUPPORTED = {
    # Every case except test.tpv8 has not opted into python-jax-mpi
    # (PY_MPI_RANKS above). This is NOT the same finding as the empty block
    # this table used to be for the three original backends -- those three
    # ARE covered for every case; python-jax-mpi is an optional execution mode
    # of python-jax and only test.tpv8 has been gated on it so far.
    (c, 'python-jax-mpi'): _MPI_OPT_IN_REASON
    for c in _NAME_LIST if c not in PY_MPI_RANKS
}

# RELEASE_ONLY -- RETIRED 2026-09-23 (owner decision, same session as the
# numpy-axis removal above). Its only two occupants ever were
# ('test.tpv36', 'python-numpy') and ('test.tpv37', 'python-numpy') -- a
# suite-COST flag that existed solely to move an expensive numpy cell out of
# the everyday sweep. With python-numpy gone from BACKENDS entirely, both
# occupants are gone with it (a stale entry here would fail the import-time
# consistency check the same way a stale MEASURED_PEAK_RSS_GB key would), so
# the flag itself is deleted rather than kept empty as dead state. There is
# now no cell held out of the everyday sweep for cost: `run.py e2e` and
# `run.py release` select the identical cell set (see
# run_e2e.py's `select()`, which keeps --release as the rule-24
# evidence-writing mode -- docs/evidence/sweep-<sha>/summary.json -- not as a
# wider cell selection anymore).

# MEASURED_PEAK_RSS_GB stays as measured evidence -- it is what the RELEASE
# tier's evidence artifact and the local `run.py all`/`run.py release` sweeps
# use to reason about which cells fit which box. It is no longer what decides
# CI_CELLS (see the redefinition below, 2026-09-23): CI stopped covering
# physics with this table on 2026-09-23 (owner-approved test-methodology
# change) and now runs exactly one portability smoke cell set. The numbers
# below remain the measured, provenanced record they always were; exceeding a
# GitHub ubuntu-22.04 runner's 7 GB produces a SIGTERM with no output (exit
# 143), a resource kill and not a test result.
MEASURED_PEAK_RSS_GB = {
    # Full-length runs, one case at a time on the 64-core development box,
    # measured with `/usr/bin/time -v` (exact peak) unless noted otherwise.
    # python-numpy rows retired 2026-09-23 with the backend axis itself
    # (owner decision) -- a stale key here would fail this module's own
    # import-time consistency check (every key must be in CASES x BACKENDS).
    # The historical numbers (tpv8 1.45 GB/62.5s, meng2023a 2.545 GB/313.4s,
    # meng2023cb 2.543 GB/419.7s, tpv29 2.705 GB/1711.0s, tpv36 2.91 GB/
    # 1415.9s, tpv37 2.91 GB/95.4s) are in git history, not reproduced here.
    ('test.tpv8', 'python-jax'): 1.91,     # 30.7 s wall
    ('test.tpv104', 'python-jax'): 3.42,
    ('test.tpv10', 'python-jax'): 3.23,
    ('test.tpv1053d', 'python-jax'): 3.95,
    # Sampled every 20 s (not /usr/bin/time -v): a LOWER BOUND on true peak.
    # Either way this cell cannot run on a 7 GB runner.
    ('test.drv.a6', 'python-jax'): 9.57,
    # Each of these also PASSED its bound on this same run
    # (testsys.compare.compare_cell), not just produced a number:
    # meng2023a jax 3.229082e-09, meng2023cb jax 4.395842e-09, tpv29 jax
    # 1.276548e-13, all against THRESHOLD=1e-3.
    ('test.meng2023a', 'python-jax'): 4.025,      # 70.7 s wall
    ('test.meng2023cb', 'python-jax'): 4.020,     # 65.1 s wall
    ('test.tpv29', 'python-jax'): 4.052,     # 198.0 s wall -- by far the
                                              # most expensive cell in the
                                              # table.
    # Peak RSS is a SETUP-time property here, not a full-run one: a
    # truncated par.term=0.3 (28 steps) run peaked within 0.09% of the
    # full-length (par.term=6, 556 steps) run it stands in for -- jax needs
    # at least one full post-compile step; 28 is comfortably enough.
    ('test.tpv36', 'python-jax'): 4.34,     # 38.6 s wall (truncated, validated method)
    ('test.tpv37', 'python-jax'): 4.20,     # 30.1 s wall (truncated, validated method)
    # python-jax-mpi, test.tpv8, 4 ranks: SUM of per-rank peak RSS.
    # RE-MEASURED 2026-09-24 for item 64 (rank-local boxes, (2,2,1)): each
    # rank wrapped in `/usr/bin/time -f %M` under mpirun (the process's own
    # getrusage maxrss, not a sampled poll): 805 MB max rank, 3.10 GB summed,
    # full 5 s gate term. The SAME method on the pre-item-64 1D
    # build-then-restrict code, same session: 1193 MB max, 4.65 GB summed.
    # (The 5.34 GB this row carried before was a 5 s-interval RSS poll of the
    # older code on 2026-09-21 -- a different instrument, kept in git history.)
    # NOT added to CI_CELLS.
    ('test.tpv8', 'python-jax-mpi'): 3.10,
}
CI_RUNNER_RAM_GB = 7.0

# JAX_MEASURED_CORES -- billed CORE cost of a python-jax cell for
# testsys/e2e/run_e2e.py's concurrency budgeting (cell_cost), not a physics
# bound and not per-case: XLA's CPU backend threads its own ops even though
# nothing in the launch path (run_standalone) requests any thread count, so a
# single serial `python3 -m eqdyna ... --backend jax` process is NOT one
# core's worth of load, and billing it as 1 (the sweep's behaviour before the
# 2026-09-23 speed campaign) undercounts every concurrent jax cell by
# ~2.5x -- exactly what let the BEFORE sweep run every cell at 2-3x its solo
# time (docs/SESSION_LOG_2026-09-23_autopilot.md section NN: tpv8 numpy
# 209.4s in-sweep vs 66s solo, tpv29 numpy 648.7 vs ~312, box load 35-58).
#
# MEASURED directly this session (not merely re-cited): `/usr/bin/time -v`
# around `python3 -m eqdyna <serial test.tpv8 case> --backend jax`, gate term
# (par.term=5.0), box at ~33/64 cpus busy (uptime load 36) --
# "Percent of CPU this job got: 252%" over a 53.7s wall run. This
# corroborates, rather than merely repeats, the 249% figure section NN
# attributes to an earlier run ("iris-X") under similar tenancy.
#
# Rounded UP (ceil, applied where this is consumed) -- under-billing a real
# cost is the defect this constant exists to close; a fractional core is
# never truncated down to fewer cores than were actually observed busy.
JAX_MEASURED_CORES = 2.52

# CI_CELLS -- REDEFINED 2026-09-23 (owner-approved test-methodology change).
#
# CI no longer runs the e2e sweep for physics coverage; that job is the local/
# release tiers' now (`run.py e2e`, `run.py release`, both at GATE_TERM_S --
# there is only one term). CI's remaining e2e job is a SMOKE TEST: one case,
# both remaining backend IMPLEMENTATIONS (fortran, python-jax), at
# GATE_TERM_S (5 s), whose purpose is portability -- a clean checkout,
# fresh-installed dependencies, and (for the
# fortran cell) the RUNNER's own mpich rather than this box's Open MPI 4.1.1
# -- not physics regression. test.tpv8 is the smallest gated case, and its one
# committed frt.canonical.txt is simply THE reference -- there is no second
# file for any cell, gate or release.
#
# Physics coverage across the full case x backend table is now the job of the
# WIDER tiers: `run.py e2e` (every case, run locally/on demand) and
# `run.py release` (every case, the release gate, PROJECT_RULES rule 15/16)
# -- both at the same GATE_TERM_S and, since RELEASE_ONLY's retirement
# (2026-09-23), the same cell SET. Neither of those is CI; CI's old
# 25-of-30-cell memory-driven selection over MEASURED_PEAK_RSS_GB is retired
# along with the jobs that ran it.
#
# BACKEND AXIS (2026-09-23, same owner decision that removed python-numpy
# from BACKENDS): the smoke set is test.tpv8 x every remaining backend
# IMPLEMENTATION -- fortran, python-jax. python-jax-mpi stays the per-case
# opt-in execution MODE it always was (see PY_MPI_RANKS above) and is not a
# third CI_CELLS row.
CI_CELLS = (
    ('test.tpv8', 'fortran'),
    ('test.tpv8', 'python-jax'),
)


def is_supported(case, backend):
    return (case, backend) not in UNSUPPORTED



def unsupported_reason(case, backend):
    """The recorded reason a cell is not gated. Raises for a supported cell:
    asking why a supported cell was skipped is itself a bug (rule 2)."""
    if is_supported(case, backend):
        raise KeyError('%s x %s is supported -- it has no unsupported reason'
                       % (case, backend))
    return UNSUPPORTED[(case, backend)]


def gate_description(case):
    if GATE[case] == 'abs-max':
        return 'max|diff| <= %.1e' % CASE_BOUND[case]
    return ('flips <= %d, median|dfnft| <= %.2gs, phys <= %.1e'
            % (DRV_A6['total_flip_bound'], DRV_A6['median_fnft_bound'],
               DRV_A6['phys_max_bound']))


def cells(cases=None, backends=None):
    """Return (runnable, declared_unsupported) for the requested selection.

    `runnable` is the list of (case, backend) the sweep will execute;
    `declared_unsupported` is the list of (case, backend, reason) INSIDE the
    requested selection that this table refuses to run.

    The caller decides what an unsupported cell in the selection means, and
    the sweep's answer is: it fails (see run_e2e.py). A default sweep
    (cases=None, backends=None) selects every cell of the table, so the
    unsupported ones are always PRINTED; an explicit selection that names an
    unsupported cell is a hard failure, because the caller asked for coverage
    that does not exist.
    """
    sel_cases = tuple(cases) if cases else CASES
    sel_backends = tuple(backends) if backends else BACKENDS
    for c in sel_cases:
        if c not in CASES:
            raise ValueError('unknown case %r -- known cases: %s' % (c, ', '.join(CASES)))
    for b in sel_backends:
        if b not in BACKENDS:
            raise ValueError('unknown backend %r -- known backends: %s'
                             % (b, ', '.join(BACKENDS)))
    runnable, unsupported = [], []
    for c in sel_cases:
        for b in sel_backends:
            if is_supported(c, b):
                runnable.append((c, b))
            else:
                unsupported.append((c, b, UNSUPPORTED[(c, b)]))
    return runnable, unsupported


def coverage_report(runnable, declared_unsupported, selection_label):
    """The lines a run prints BEFORE it starts, so its own output states
    exactly what it is about to cover. Requirement zero of this sweep: a green
    result must never be readable as broader than it is."""
    total = len(CASES) * len(BACKENDS)
    selected = len(runnable) + len(declared_unsupported)
    lines = [
        'selection: %s' % selection_label,
        'coverage : %d of %d cells in the full %d case x %d backend table'
        % (selected, total, len(CASES), len(BACKENDS)),
        'cases    (%d): %s' % (len(CASES), ', '.join(CASES)),
        'backends (%d): %s' % (len(BACKENDS), ', '.join(BACKENDS)),
        'will RUN %d cell(s):' % len(runnable),
    ]
    for c, b in runnable:
        lines.append('  %-16s %-13s artifacts=%-7s gate: %s'
                     % (c, b, '+'.join(ARTIFACTS[b]), gate_description(c)))
    lines.append('declared UNSUPPORTED, %d cell(s) -- declared, not absent:'
                 % len(declared_unsupported))
    for c, b, reason in declared_unsupported:
        lines.append('  %-16s %-13s %s' % (c, b, reason))
    chosen = set(runnable) | set((c, b) for c, b, _ in declared_unsupported)
    not_selected = [(c, b) for c in CASES for b in BACKENDS if (c, b) not in chosen]
    lines.append('NOT in this selection, %d cell(s)%s'
                 % (len(not_selected),
                    (': ' + ', '.join('%s x %s' % cb for cb in not_selected))
                    if not_selected else ''))
    return lines


# --- consistency, enforced at import time -------------------------------------
# A case added to testNameList.py with no entry here would otherwise run ungated
# or crash mid-sweep. Fail at import, naming the case.
# Every gated case must have a committed reference. Absent = broken checkout,
# and a checkout that cannot compare must not be able to report a pass.
for _c in CASES:
    _ref = os.path.join(REPO_ROOT, 'test.reference.results', _c,
                        'frt.canonical.txt')
    if not os.path.isfile(_ref):
        raise RuntimeError(
            '%s is gated but %s is missing -- restore it from git. A gated case '
            'with no reference cannot be compared, and "could not compare" must '
            'never read as "passed" (rule 2).' % (_c, _ref))

_missing_bound = [c for c in CASES if c not in CASE_BOUND]
_missing_gate = [c for c in CASES if c not in GATE]
_extra = [c for c in set(CASE_BOUND) | set(GATE) if c not in CASES]
if _missing_bound or _missing_gate or _extra:
    raise RuntimeError(
        'testsys/matrix.py is out of step with testNameList.py: '
        'missing CASE_BOUND for %r, missing GATE for %r, entries for unknown '
        'cases %r. Every gated case declares its bound and its gate here.'
        % (_missing_bound, _missing_gate, _extra))
for _c, _g in GATE.items():
    if _g not in ('abs-max', 'flip-budget'):
        raise RuntimeError('%s has unknown gate %r' % (_c, _g))
    if _g == 'abs-max' and CASE_BOUND[_c] is None:
        raise RuntimeError('%s is gated abs-max but has no bound' % _c)
    if _g == 'flip-budget' and CASE_BOUND[_c] is not None:
        raise RuntimeError('%s is gated flip-budget but also carries a scalar '
                           'bound -- one gate per case' % _c)
if set(NSTRESS_CONVENTION) != set(CASES) or any(
        v[0] not in ('extension', 'compression') for v in NSTRESS_CONVENTION.values()):
    raise RuntimeError(
        'testsys/matrix.py NSTRESS_CONVENTION must name every gated case, each '
        "'extension' or 'compression' with its spec citation: missing %r, "
        'unknown %r' % (sorted(set(CASES) - set(NSTRESS_CONVENTION)),
                        sorted(set(NSTRESS_CONVENTION) - set(CASES))))
for (_c, _b) in list(UNSUPPORTED) + list(CI_CELLS):
    if _c not in CASES or _b not in BACKENDS:
        raise RuntimeError('table entry for unknown or ungated cell %r x %r'
                           % (_c, _b))
for (_c, _b) in list(MEASURED_PEAK_RSS_GB):
    if _c not in CASES or _b not in BACKENDS:
        raise RuntimeError('measurement for unknown cell %r x %r' % (_c, _b))
# RELEASE_ONLY's own consistency check retired with the flag itself
# (2026-09-23) -- there is nothing left to validate.
