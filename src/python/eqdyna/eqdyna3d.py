"""eqdyna3d.py <- src/eqdyna3d.f90. The main entry point.

Builds the solver state `S` ENTIRELY from the case inputs -- bGlobal.txt,
bModelGeometry.txt, bFaultGeometry.txt, bMaterial.txt,
on_fault_vars_input.nc -- via readInputFiles.py, then meshgen.py and
assembleGlobalMass.py, hands it to driver.py, and writes frt.txt0 through
library_output.py. No Fortran is in the loop: there are ZERO pydump_*.txt
reads here (pydump.py exists only for run_parity.py's per-step diagnostic).

`backend` is the one argument eqdyna3d.f90 does not have. It selects numpy
or jax.numpy and is threaded down to a SINGLE driver/faulting/fric/
assembleGlobalKU implementation; there is no per-backend and no per-friclaw
solver module to dispatch between.

SCOPE, enforced by loud refusals in build_solver_state rather than by silent
partial runs: ntotft==1, serial (npx==npy==npz==1) OR -- python-jax-mpi only --
one rank's box of an MPI4NodalQuant.DECOMP decomposition, C_degen==0 (planar) OR
C_degen>3 (dipping, wedge-degeneration -- meshgen.py's build_elements/
build_fault_geometry port library_degeneration.f90's wedge()/reorder(); see
those functions' docstrings and testsys/parity/evidence_c_degen_port.py).
C_degen>3's MESH is fully ported and verified against Fortran on
test.tpv36. Its DYNAMICS are ALSO ported: assembleGlobalMass.py's
compute_element_shape/assemble_mass apply calcGlobalShapeFunc.f90:22-28's
elemTypeArr==11/12 shape-function merge, and assembleGlobalKU.py's element
dispatch routes elemType>10 through the interior kernel
(assembleGlobalKU.f90:25) -- see both modules' docstrings for what changed
and why compute_hourglass/calcElemKU/calcHourglassResist/calcElemMass
needed no changes of their own. Verified directly against Fortran (not
end-to-end self-consistency alone) by testsys/parity/evidence_wedge_kernel.py
on synthetic hex/wedge/shallow-dip element geometries; see that script's
docstring for the exact tolerance and what was checked. friclaw 1-5 are all
implemented.

S-dict provenance, field by field:
  N, E, NEQ, nen, ned            -- meshgen.py M1-M3 (nen=8, ned=3 are the
                                    globalvar.f90 hex-element/dof-per-node
                                    constants, never case input)
  dt, w, rdampk, rdampm, kapa_hg,
  R, nPML, vmaxPML, PMLb, grav,
  C_elastic, roumax, rhow, gamar,
  slipRateThres, xsource, ysource,
  zsource, nucR, nucRuptVel,
  nucdtau0, nucT, TPV, C_nuclea,
  nucfault, friclaw, nstep       -- readInputFiles.read_bglobal + the
                                    globalvar.f90 PARAMETER constants
  meshCoor, conn, elemType, mat  -- meshgen.py M1/M2, converted from this
                                    port's 1-indexed convention to 0-indexed
  eledet, eleshp, ss, phi        -- assembleGlobalMass.py M7.5
  nodalMassArr, fnms             -- assembleGlobalMass.py M6.5
  ndof, eq_ids                   -- meshgen.py M3
  nsmp1, nsmp2, un, us, ud, arn  -- meshgen.py M4
  fric_init                      -- readInputFiles.read_on_fault_vars
  v1_init                        -- assembleGlobalMass.init_vel. NOTE
                                    (pre-existing): driver.run zero-inits v1
                                    and never reads this. For mode==1 with
                                    fric(31:36) identically 0.0 that is a
                                    verified no-op; a case with mode==2 or
                                    nonzero fric(31:36) WOULD need it wired
                                    in. Flagged, not silently worked around.

C_elastic==0 (test.drv.a6) additionally needs ccosphi/sinphi/tv
(readInputFiles.f90's readmaterial) and init_stress (meshgen.f90:104's
setPlasticStress lithostatic pre-stress). Both are computed here, once,
before the loop, exactly as the Fortran does; assembleGlobalKU.build gates
their USE on C_elastic==0 and raises loudly if they are missing rather than
silently zero-initialising a plastic run.
"""
import argparse
import contextlib
import os
import sys
import time

import numpy as np

from . import assembleGlobalMass, checkInputConsistency, driver, func_lib, library_output, meshgen, readInputFiles
from . import backend as _backend
from . import profile_emit as _profile_emit

# The friction laws this port implements. EVERY one of them is served by the
# SAME code -- eqdyna/{driver,faulting,fric,assembleGlobalKU,backend}.py --
# with the friclaw dispatch inside faulting.py exactly where faulting.f90:21-22
# puts it, and with the backend as an argument rather than a second module.
# build_solver_state is therefore friclaw-agnostic: there is no per-friclaw
# and no per-backend solver module left to dispatch between (the two tables
# of three port*.py modules that used to stand here are deleted).
SUPPORTED_FRICLAW = (1, 2, 3, 4, 5)

DEFAULT_BACKEND = 'jax'

_NUMPY_WIDE_AFFINITY_ENV = 'EQDYNA_NUMPY_ALLOW_WIDE_AFFINITY'


def _narrow_numpy_affinity():
    """numpy's hot kernels (calcHourglassResist/assembleGlobalKU's elementwise
    ops, np.add.at scatter-adds -- 85%+ of a step by the pathway item 40/
    testsys/perf profile) run on ONE thread regardless of core count: OpenBLAS
    threading only engages BLAS-dispatched calls (matmul/dot), which is a
    small fraction of this port's per-step cost. MEASURED on this box
    (testsys/perf/run_scaling.py, test.tpv104, numactl-pinned, per-step by
    difference): 1 cpu 9443 ms/step, 2 cpus (same NUMA node) 9872 ms/step,
    4 cpus (same node) 8900 ms/step -- flat within noise, no speedup from
    extra cores. 16 cpus SPANNING TWO NUMA NODES: 19055 ms/step, ~2x WORSE.
    Mechanism: a single-threaded process given an affinity mask wider than
    one NUMA node can be migrated by the OS scheduler between cpus on
    DIFFERENT nodes over the run, stranding its memory on whichever node it
    happened to first-touch -- a real regression, not noise (see
    testsys/perf/run_scaling.py's docstring for the full measurement).

    So there is no configuration that makes more cores help the numpy
    backend, and a wide/unpinned affinity mask (e.g. a caller's plain
    `python3 -m eqdyna` with no taskset/numactl at all, which sees every cpu
    on the box) can make it WORSE by construction. This narrows the numpy
    backend to exactly the FIRST cpu in whatever affinity mask it was given
    -- never wider than the caller's own choice, just never spanning more
    than it can use.

    Does nothing for jax (which does show real, if limited, benefit from
    extra cores -- see the same script's `python-jax` rows) and does nothing
    where `os.sched_getaffinity` does not exist (macOS has no such syscall;
    this is an optimisation, not a correctness check, so absence of the API
    is a silent no-op here by design, not a rule-2 violation).

    `EQDYNA_NUMPY_ALLOW_WIDE_AFFINITY=1` disables this, for anyone
    deliberately re-measuring raw multi-core numpy behaviour (e.g. this
    tool's own scaling curve, after a future kernel rewrite that might
    finally give numpy something to parallelise).

    CORRECTION (wei-lin, gate-axis-3 review, before this landed): the first
    version of this always chose `min(current)` -- the SAME lowest-numbered
    cpu for every process. Verified directly (two `python3 -m eqdyna`
    invocations launched concurrently, `/proc/<pid>/status`'s
    `Cpus_allowed_list`): both narrowed to cpu 0 exactly. That is fine for
    ONE process at a time (this function's original, tested use case, and
    still true for it), but `testsys/e2e/run_e2e.py --jobs N`'s own
    concurrent-cell sweep (this repo's wider LOCAL gate, `run.py all`) then
    forces every simultaneous numpy cell onto the identical physical core
    while 63 others sit idle -- reproduced directly: two concurrent
    tpv8-scale cells both measured `Cpus_allowed_list: 0`. Picking a cpu
    DETERMINISTICALLY BY PID from within the given mask instead keeps the
    single-process behaviour (still narrows to exactly one cpu, so the
    NUMA-migration fix above is unchanged) while spreading concurrent
    processes across different cpus without any inter-process coordination.
    """
    if os.environ.get(_NUMPY_WIDE_AFFINITY_ENV) == '1':
        return
    getter = getattr(os, 'sched_getaffinity', None)
    setter = getattr(os, 'sched_setaffinity', None)
    if getter is None or setter is None:
        return
    current = getter(0)
    if len(current) > 1:
        ordered = sorted(current)
        chosen = ordered[os.getpid() % len(ordered)]
        setter(0, {chosen})


def _resolve_solver(friclaw, backend):
    """Returns a run(S, nsteps, verbose) callable for `friclaw` under
    `backend` -- one solver for every friclaw, so this only VALIDATES the
    pair and binds the backend; there is no module to pick between.

    NO FALLBACK (PROJECT_RULES rule 2). backend='jax' with jaxlib missing is a
    hard failure, not a quiet demotion to NumPy.

    This used to print a NOTICE and return the NumPy solver. That is fatal the
    moment backends are compared: a run labelled `jax` would BE numpy, and a
    "JAX is no faster than NumPy" result would be NumPy measured twice. A
    printed notice does not help either -- it scrolls past in a log while the
    number it invalidates is what gets recorded.
    """
    if friclaw not in SUPPORTED_FRICLAW:
        raise NotImplementedError('_resolve_solver: friclaw=%d is not implemented '
                                   '(implemented: %r)' % (friclaw, list(SUPPORTED_FRICLAW)))
    if backend not in ('numpy', 'jax'):
        raise ValueError("_resolve_solver: backend must be 'jax' or 'numpy' (got %r)" % backend)
    return lambda S, nsteps=None, verbose=True: run(
        S, nsteps=nsteps, verbose=verbose, backend=backend)


def active_device(backend):
    """What actually ran -- recorded, never inferred from the request.

    A backend matrix that prints the REQUESTED device is worthless: JAX_PLATFORMS
    can be overridden, a GPU can be busy or invisible, and a run asking for gpu
    can legitimately land on cpu. Every timing row must carry what actually ran.

    `backend` is required: reporting jax.devices()[0] for a NUMPY run named a
    GPU that the run never touched -- observed, and exactly the kind of label
    that would put a bogus device into a perf table.
    """
    if backend == 'numpy':
        return 'host CPU (numpy; jax not used)'
    import jax
    d = jax.devices()[0]
    return '%s:%d (%s)' % (d.platform, d.id, getattr(d, 'device_kind', '?'))


def report_dropped_stations(xonfs, x4nds, anonfs, off_matches, meshCoor,
                             off_z_valid, tol):
    """Port of report_dropped_onfault_st / report_dropped_offfault_st
    (library_output.f90; board rows 116 and 94). On-fault stations: unchanged
    -- name, on stdout, every requested station that matched no fault node and
    so gets no file (neither snaps nor refuses; out of this mission's scope).
    ntotft == 1 is the only case case.setup allows, so on-fault stations are
    all fault 1.

    Off-fault (row 94, owner ruling 2026-09-24): depth now snaps to the
    nearest node (build_station_matching), so every station whose (x,y) is
    inside the mesh matches SOME node -- possibly not the one requested. This
    reports that difference the same way report_dropped_offfault_st
    (library_output.f90) does:
      - matched at a different node than requested: a SNAP, named "snapped
        off-fault station" (and "...would otherwise be a dropped off-fault
        station", so the line stays equally loud and equally grep-able as
        the true-drop case below).
      - matched on no node at all: a true DROP, "dropped off-fault station",
        unchanged.

    SNAP is gated on the DEPTH difference alone (|actual z - requested z| >
    `tol`), not the full 3-axis distance -- x and y already snapped to the
    nearest node before this fix (silently, forever) and that is unchanged
    and out of scope; gating broadly would flood this NOTICE with stations
    that were never at risk of being dropped (measured: it would have
    reported 19 for test.tpv36/test.tpv37, 0 of them depth-caused, and 6 for
    test.tpv10 instead of the 2 this fix actually recovers). Still prints the
    full (x,y,z) requested vs actual and full 3-axis distance for each
    reported station, since a depth-snapped station can shift in x/y too
    (test.tpv10 stations 9/10). `tol` is the caller's `params['tol']` --
    finding 5 (row 94 audit, 2026-09-25): this used to hardcode a second
    `1.0e-5` literal here instead of reading the one place this port already
    reproduces Fortran globalvar.f90:209's `tol = 1.0d-5` PARAMETER
    (readInputFiles.build_params).

    Finding 6 (row 94 audit, 2026-09-25): a true DROP used to say "match no
    grid node (outside the mesh)" unconditionally -- not always true (the
    Fortran side's tpv8 station 11 is dropped by a pre-existing y-partition-
    boundary gap in setSurfaceStation, not by being outside the mesh; the
    serial Python port carries the identical un-gated iy==0/iy==ny-1 case in
    principle, just unmasked at different domain sizes than that MPI split).
    `off_z_valid[i-1]` (build_station_matching's return, mirroring Fortran's
    x4ndsZValidPersist) is the one cause this function HAS checked, so a
    drop's message says either "requested depth is outside the physical mesh
    band" (a CONFIRMED cause) or that the cause is not further diagnosed here
    (depth was in-band; x, y, or the y-boundary gap above are candidates, but
    none is checked by this function, so none is named).

    Coordinates arrive in metres. `meshCoor` is build_node_coordinates'
    1-indexed (row 0 unused) array; `nc` (off_matches' second element) is
    already that same 1-indexed node id."""
    # Fortran's order: off-fault first (checkOffFaultStationCoverage,
    # eqdyna3d.f90:126), then on-fault (checkOnFaultStationCoverage).
    off_matched = {sc: nc for sc, nc in off_matches}
    off_dropped = [i for i in range(1, x4nds.shape[1] + 1) if i not in off_matched]
    off_snapped = []
    for i in range(1, x4nds.shape[1] + 1):
        if i in off_matched:
            actual = meshCoor[off_matched[i]]
            requested = x4nds[:, i - 1]
            if abs(float(actual[2]) - float(requested[2])) > tol:
                dist = float(np.linalg.norm(actual - requested))
                off_snapped.append((i, requested, actual, dist))
    if off_snapped:
        print(' NOTICE: %d of %d requested off-fault stations do not sit exactly '
              'on a grid z-plane' % (len(off_snapped), x4nds.shape[1]))
        print('   (setSurfaceStation, meshgen.f90: depth now snaps to the nearest node '
              'instead of requiring an exact grid-plane match; x and y unchanged, they '
              'already snapped to the nearest node)')
        for i, requested, actual, dist in off_snapped:
            print('   snapped off-fault station %d (would otherwise be a dropped '
                  'off-fault station): requested x,y,z =%10.3f%10.3f%10.3f km, '
                  'actual x,y,z =%10.3f%10.3f%10.3f km, distance =%8.3f km'
                  % (i, requested[0] / 1000.0, requested[1] / 1000.0, requested[2] / 1000.0,
                     actual[0] / 1000.0, actual[1] / 1000.0, actual[2] / 1000.0,
                     dist / 1000.0))
    if off_dropped:
        print(' WARNING: %d of %d requested off-fault stations match no grid '
              'node and get NO body* file'
              % (len(off_dropped), x4nds.shape[1]))
        for i in off_dropped:
            if off_z_valid[i - 1]:
                cause = ('cause not checked here -- requested depth is within the physical '
                          'mesh band; the miss may be x or y outside the mesh, or a known '
                          'y-partition-boundary gap in setSurfaceStation')
            else:
                cause = ('checked cause: requested depth is outside the physical, non-PML '
                          'mesh band')
            print('   dropped off-fault station %d at x,y,z =%10.3f%10.3f%10.3f km (%s)'
                  % (i, x4nds[0, i - 1] / 1000.0, x4nds[1, i - 1] / 1000.0,
                     x4nds[2, i - 1] / 1000.0, cause))
    on_matched = {sc for fs, sc, ift in anonfs}
    on_dropped = [i for i in range(1, xonfs.shape[1] + 1) if i not in on_matched]
    if on_dropped:
        print(' WARNING: %d of %d requested on-fault stations match no fault '
              'node and get NO faultst* file' % (len(on_dropped), xonfs.shape[1]))
        print('   (setOnFaultStation, meshgen.f90: along-strike x and depth z '
              'must both equal a fault node within tol)')
        for i in on_dropped:
            print('   dropped on-fault station %d (fault 1) at x,z =%10.3f%10.3f km'
                  % (i, xonfs[0, i - 1] / 1000.0, xonfs[1, i - 1] / 1000.0))
    return on_dropped, off_dropped


def build_solver_state(case_dir, part=None):
    """Builds the full `S` dict eqdyna/driver.py's `run()` expects,
    reading ONLY case-input files (bGlobal.txt/bModelGeometry.txt/
    bFaultGeometry.txt/bMaterial.txt/on_fault_vars_input.nc) via
    readInputFiles.py, then meshgen.py/assembleGlobalMass.py -- no pydump_*.txt
    anywhere. Returns (S, mesh) where mesh is a dict of the raw 1-indexed
    builder outputs (meshCoor, nsmp, eq_nums, ...) library_output needs later
    (S itself is converted to loading.load()'s 0-indexed convention).

    `part` (an MPI4NodalQuant.Partition) builds ONE RANK'S BOX instead, with
    rank-local numbering, the way meshgen.f90 does: the same builders handed
    the rank's local line slices plus the GLOBAL model bounds. No global mesh
    is built. This function does NO communication, so `nodalMassArr`, `fnms`
    and `arn` come back as this rank's PARTIAL sums -- exactly what Fortran
    holds after its element loop and before MPI4arn / MPI4NodalQuant --
    and MPI4NodalQuant.setup_exchange completes them. `mesh` then also
    carries what that needs (n_local, flt_lists, the 1-indexed partials) and
    the two global censuses the run checks its ownership against. part=None
    is the serial path, and its arithmetic is unchanged.
    """
    params, g = readInputFiles.build_params(case_dir)
    # checkInputConsistency.f90:7-19 <-> checkInputConsistency.check, called
    # at the SAME point Fortran calls it: every bFile has been read, no mesh
    # or solver work has started (eqdyna3d.f90:104). Raises
    # InputConsistencyError (caught in main()'s _abort) rather than
    # NotImplementedError, since these are configuration refusals with a
    # numbered exit code to match, not scope gaps in this port.
    checkInputConsistency.check(g['C_elastic'], g['output_plastic'], params['rat'])
    if g['ntotft'] != 1:
        raise NotImplementedError('build_solver_state: only ntotft==1 is supported (got %d)'
                                   % g['ntotft'])
    if g['friclaw'] not in SUPPORTED_FRICLAW:
        raise NotImplementedError('build_solver_state: friclaw=%d is not implemented '
                                   '(implemented: %r)' % (g['friclaw'], list(SUPPORTED_FRICLAW)))
    case_decomp = (g['npx'], g['npy'], g['npz'])
    if part is None and case_decomp != (1, 1, 1):
        raise NotImplementedError('build_solver_state: only serial (npx=npy=npz=1) is supported '
                                   '(got %r)' % (case_decomp,))
    if part is not None and case_decomp not in ((1, 1, 1), part.dims):
        # The python MPI path takes its split from MPI4NodalQuant.DECOMP (the
        # table the Fortran is run with). A case that names a DIFFERENT split
        # in bGlobal.txt is refused rather than silently re-split.
        raise NotImplementedError(
            'build_solver_state: bGlobal.txt names npx,npy,npz=%r but this '
            'run is decomposed %r (MPI4NodalQuant.DECOMP[%d])'
            % (case_decomp, part.dims, part.nranks))
    # C_degen: meshgen.py's build_elements/build_fault_geometry (and
    # readInputFiles.py's fltxyz(2,4,i) derivation) now port BOTH C_degen==0
    # (planar fault, tpv8/tpv104/tpv10/drv.a6) and C_degen>3
    # (wedge-degeneration, tpv36/tpv37 -- library_degeneration.f90's wedge()/
    # reorder()). 0<C_degen<=3 stays refused: checkIsOnFault itself takes
    # neither if/elseif branch there (isOnFault stays 0 unconditionally), a
    # degenerate Fortran behavior this port does not silently mimic.
    if not (params['C_degen'] == 0.0 or params['C_degen'] > 3.0):
        raise NotImplementedError('build_solver_state: only C_degen==0 or C_degen>3 is '
                                   'supported (got %r)' % params['C_degen'])
    # C_degen>3 + C_elastic==0 (plastic): NOT supported. meshgen.f90:104 only
    # calls setPlasticStress once per (ix,iy,iz) grid point, using whichever
    # `elemCount` is current AFTER wedge() has (possibly) split it into two
    # sub-elements -- so the type-11 (below-fault) sub-element's lithostatic
    # depth is simply never written by the Fortran for that slot, a quirk
    # this port has not verified end to end. Every currently supported
    # C_degen>3 case (tpv36/tpv37) has C_elastic==1, so elem_depth/
    # init_stress are never read regardless (see build_elements' docstring)
    # -- refusing the untested combination explicitly rather than silently
    # feeding it a depth value with no Fortran-verified meaning.
    if params['C_degen'] > 3.0 and g['C_elastic'] == 0:
        raise NotImplementedError(
            'build_solver_state: C_degen>3 (wedge-degenerate elements) combined with '
            'C_elastic==0 (plastic) is not supported -- meshgen.f90:104\'s '
            'setPlasticStress is never called for the type-11 wedge sub-element '
            '(it runs once per grid point, against whichever elemCount wedge() left '
            'current -- the type-12 slot), a depth-assignment quirk this port has not '
            'verified; C_elastic==1 (tpv36/tpv37) is unaffected since elem_depth/'
            'init_stress are never read then.')
    # (The insertFaultType>0 x friclaw==5 refusal that stood here is GONE, and
    # not by relaxing it: faulting.f90:205-212's min_norm/max_norm clamp was
    # ported in port_rsf.py and missing from port_tp.py, so the combination was
    # genuinely unimplemented. There is now one solveRSF, the clamp is in it,
    # and every friclaw reaches the same code -- so there is nothing left to
    # refuse.)

    material = readInputFiles.read_bmaterial(
        os.path.join(case_dir, 'bMaterial.txt'), g['nmat'], g['n2mat'])

    # Row 114 -- station output. bStations.txt is read unconditionally (every
    # case.setup-generated case dir has one, scripts/case.setup:111); MATCHING
    # (meshgen.build_station_matching, which needs the GLOBAL grid lines) runs
    # only on the serial path -- see run_case_mpi's loud warning for why the
    # python-jax-mpi path leaves st_on_idx/st_off_idx empty instead.
    xonfs, x4nds = readInputFiles.read_bstations(os.path.join(case_dir, 'bStations.txt'))
    st_on_total, st_off_total = xonfs.shape[1], x4nds.shape[1]

    xline, yline, zline, pmlb, bounds = meshgen.build_grid_lines(params)
    if part is None:
        anonfs, off_matches, off_z_valid = meshgen.build_station_matching(
            xline, yline, zline, params, xonfs, x4nds,
            pmlb['zmin0'], bounds[2][1])
        st_on_idx = np.array([fs - 1 for fs, sc, ift in anonfs], dtype=np.int64)
        st_on_strike_m = np.array([xonfs[0, sc - 1] for fs, sc, ift in anonfs])
        st_on_depth_m = np.array([xonfs[1, sc - 1] for fs, sc, ift in anonfs])
        st_off_idx = np.array([nc - 1 for sc, nc in off_matches], dtype=np.int64)
        st_off_x_m = np.array([x4nds[0, sc - 1] for sc, nc in off_matches])
        st_off_y_m = np.array([x4nds[1, sc - 1] for sc, nc in off_matches])
        st_off_z_m = np.array([x4nds[2, sc - 1] for sc, nc in off_matches])
        # Row 94: the snap report (and the header stamp write_offfault_stations
        # uses) needs the ACTUAL matched node location, which needs meshCoor --
        # not built yet at this point. report_dropped_stations is therefore
        # called further down, once meshCoor exists (see there).
    else:
        st_on_idx = np.zeros(0, dtype=np.int64)
        st_on_strike_m = st_on_depth_m = np.zeros(0)
        st_off_idx = np.zeros(0, dtype=np.int64)
        st_off_x_m = st_off_y_m = st_off_z_m = np.zeros(0)
    # fltxyz(2,4,1) (readInputFiles.f90:139-143), the fault dip angle in
    # radians used by output_onfault_st's down-dip-distance conversion
    # (library_output.f90:51/68) -- same formula as meshgen.py's
    # build_fault_geometry `fdip`, recomputed here (pure function of
    # params['C_degen']) rather than plumbed through that builder's return.
    fault_dip_rad = (params['C_degen'] * np.pi / 180.0 if params['C_degen'] > 3.0
                      else 90.0 * np.pi / 180.0)

    model_bound = None
    if part is not None:
        # The global lines are O(nx+ny+nz) and kept only for the two
        # censuses; every builder below sees the rank's SLICES. PMLb and the
        # model bounds stay global (meshgen.f90:32-37), or setNumDof and the
        # fixed-boundary test would reclassify nodes on subdomain faces.
        glines = (xline, yline, zline)
        (xline, yline, zline), offsets = part.slice_lines(xline, yline, zline)
        model_bound = bounds
    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(
        xline, yline, zline, params, model_bound=model_bound)
    if part is None:
        # Row 94: the ACTUAL matched node's (x,y,z) per off-fault station,
        # meshCoor being 1-indexed with row 0 unused (build_node_coordinates)
        # and `nc` (off_matches' second element) already that same 1-indexed
        # node id -- read it straight, no offset. A station this run did not
        # match at all (x or y outside the mesh) never appears in off_matches
        # and so never appears here either; report_dropped_stations reports
        # that case from off_matched/off_dropped directly, not from this array.
        st_off_actual_m = (meshCoor[np.array([nc for sc, nc in off_matches], dtype=np.int64)]
                            if off_matches else np.zeros((0, 3)))
        report_dropped_stations(xonfs, x4nds, anonfs, off_matches, meshCoor,
                                 off_z_valid, params['tol'])
    else:
        st_off_actual_m = np.zeros((0, 3))
    conn, elem_type, mat, elem_depth = meshgen.build_elements(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    num_dof, eq_start, eq_nums, total_eqs = meshgen.build_equation_numbers(
        xline, yline, zline, params, pmlb, model_bound=model_bound)
    un, us, ud, arn = meshgen.build_fault_geometry(
        xline, yline, zline, params, nsmp, model_bound=model_bound)

    fric = readInputFiles.read_on_fault_vars(
        os.path.join(case_dir, 'on_fault_vars_input.nc'), params['fxmin'], params['fzmin'],
        params['dx'], params['dz'], meshCoor, nsmp)

    # calcGlobalShapeFunc.f90:22-28 (called unconditionally, for EVERY
    # element, from assembleGlobalMass.f90:35) special-cases elemTypeArr
    # 11/12 by merging shape-function rows 3+4 and 7+8 (Hughes p.125's
    # standard hex-to-wedge collapse fix-up). PORTED: assembleGlobalMass.py's
    # compute_element_shape now takes `elem_type` and applies this merge to
    # the derivative rows before the Jacobian is built (exactly where the
    # Fortran applies it); `assemble_mass` takes `elem_type` too, for the
    # SEPARATE row-4 (shape-function VALUE) merge contm consumes
    # (`_contm_wedge_node_mass`). `compute_hourglass` needed no change: it
    # only consumes compute_element_shape's `eleshp`/`xs` output and has no
    # elemTypeArr branch of its own in the Fortran either (confirmed by
    # reading assembleGlobalMass.f90:330-377 directly).
    # assembleGlobalKU.py's element dispatch (`E_int`) was ALSO fixed to
    # route elemType>10 through the interior kernel, matching
    # assembleGlobalKU.f90:25's `elemTypeArr(nel)==1 .or. elemTypeArr(nel)>10`
    # -- without that fix wedge elements would silently contribute zero
    # interior force while still contributing mass/hourglass, a second,
    # independent gap found by reading the full call chain, not just the
    # mass-lumping formula the original refusal here named.
    # Verified directly against Fortran (real Fortran subroutine calls, not
    # a re-implementation) by testsys/parity/evidence_wedge_kernel.py, on a
    # non-degenerate hex (regression guard), a degenerate wedge, and a
    # 15-degree shallow-dip wedge (tpv36's par.dip) -- see that script's
    # docstring for the exact tolerance and what was/was not checked.

    xl = meshCoor[conn]
    det, eleshp3, xs = assembleGlobalMass.compute_element_shape(xl, elem_type)
    ss, phi48 = assembleGlobalMass.compute_hourglass(xl, xs, mat, eleshp3)
    nodalMassArr, fnms = assembleGlobalMass.assemble_mass(
        conn, mat, det, elem_type, num_dof, eq_start, eq_nums, total_eqs,
        meshCoor.shape[0] - 1)
    v1 = assembleGlobalMass.init_vel(nsmp, eq_nums, fric, total_eqs)

    N = meshCoor.shape[0] - 1  # drop the unused row-0
    E = conn.shape[0]

    # Milestone 10 (drv.a6, C_elastic==0): readInputFiles.f90's readmaterial
    # computes ccosphi/sinphi AFTER bMaterial.txt is read but BEFORE the time
    # loop -- reproduced here, verbatim formula.
    # tv is NOT derived here any more: it is read from bGlobal.txt, exactly as
    # readglobal now reads it (pathway item 24(b)). scripts/case.setup writes
    # the pre-v5.9.0 value 2*dz/3464 for any case that does not set
    # par.viscoplasticRelaxTime, so an unchanged case is unchanged here too.
    # g['bulk']/g['coheplas'] are read unconditionally by
    # readInputFiles.read_bglobal regardless of C_elastic (same bGlobal.txt
    # line for every case) -- always computed here too, harmless for
    # C_elastic==1 cases since assembleGlobalKU gates its USE on
    # C_elastic==0 (see those modules' build()).
    ccosphi = g['coheplas'] * np.cos(np.arctan(g['bulk']))
    sinphi = np.sin(np.arctan(g['bulk']))
    tv = g['tv']

    # meshgen.f90:104's setPlasticStress (called for EVERY element, both
    # interior and PML, only when C_elastic==0): lithostatic per-element
    # pre-stress, Voigt order [xx,yy,zz,yz,xz,xy] (calcB.f90's b(4,*)/
    # b(5,*)/b(6,*) confirm 4=yz,5=xz,6=xy) -- ALWAYS computed (cheap,
    # harmless when unused) so assembleGlobalKU.build can
    # raise loudly (not silently default) if C_elastic==0 but this key is
    # somehow missing, rather than silently zero-initializing plastic runs.
    strVert = -(g['roumax'] - g['rhow'] * (g['gamar'] + 1.0)) * elem_depth * 9.8
    # func_lib.dev_str_depth_taper is SCEC TPV29/30's Omega(depth), the port of
    # func_lib.f90's devStrDepthTaper (pathway item 24(c)). It is exactly 1.0
    # when the case does not configure a taper, so this is bit-for-bit the
    # pre-v5.9.0 `np.abs(strVert) * ratio` for every such case.
    devStr = np.abs(strVert) * g['devStrToStrVertRatio'] * func_lib.dev_str_depth_taper(
        elem_depth, g['devStrTaperDepthStart'], g['devStrTaperDepthEnd'])
    theta2 = 2.0 * g['str1ToFaultAngle']
    init_stress = np.zeros((E, 6))
    init_stress[:, 0] = strVert - devStr * np.cos(theta2)  # xx
    init_stress[:, 1] = strVert + devStr * np.cos(theta2)  # yy
    init_stress[:, 2] = strVert  # zz
    init_stress[:, 5] = devStr * np.sin(theta2)  # xy

    # ---- convert to loading.load()'s 0-indexed convention ----
    # -1 fixed-boundary sentinel -> sink 0; see meshgen.pack_eq_ids (a
    # vectorized, all-integer repack of the per-node loop this replaced).
    ndof0, eq_ids = meshgen.pack_eq_ids(num_dof, eq_nums, N, ncols=12)

    PMLb = np.array([pmlb['xmax0'], pmlb['xmin0'], pmlb['ymax0'], pmlb['ymin0'], pmlb['zmin0'],
                      pmlb['maxdx'], pmlb['maxdy'], pmlb['maxdz']])

    S = dict(
        N=N, E=E, NEQ=total_eqs, nen=8, ned=3, nftnd=nftnd, ntotft=1,
        nstep=g['nstep'], dt=g['dt'], w=assembleGlobalMass._W, rdampk=g['rdampk'], rdampm=0.0,
        kapa_hg=0.1, R=params['R'], nPML=params['nPML'], vmaxPML=g['vmaxPML'], PMLb=PMLb,
        grav=9.8, C_elastic=g['C_elastic'], roumax=g['roumax'], rhow=g['rhow'],
        gamar=g['gamar'], slipRateThres=g['slipRateThres'], xsource=g['xsource'],
        ysource=g['ysource'], zsource=g['zsource'], nucR=g['nucR'], nucRuptVel=g['nucRuptVel'],
        nucdtau0=g['nucdtau0'], nucT=g['nucT'], TPV=g['TPV'], C_nuclea=g['C_nuclea'],
        nucfault=g['nucfault'], friclaw=g['friclaw'], insertFaultType=g['insertFaultType'],
        meshCoor=meshCoor[1:], conn=conn - 1, elemType=elem_type, mat=mat,
        eledet=det, eleshp=np.transpose(eleshp3, (0, 2, 1)), ss=ss,
        phi=np.transpose(phi48, (0, 2, 1)),
        nodalMassArr=nodalMassArr[1:], fnms=fnms[1:], v1_init=v1[1:],
        ndof=ndof0, eq_ids=eq_ids,
        nsmp1=nsmp[:, 0] - 1, nsmp2=nsmp[:, 1] - 1,
        un=un[1:], us=us[1:], ud=ud[1:], arn=arn[1:], fric_init=fric[1:, 1:101],
        ccosphi=ccosphi, sinphi=sinphi, tv=tv, init_stress=init_stress,
        # Row 114 -- station output.
        dx=params['dx'], fault_dip_rad=fault_dip_rad,
        # Column 8 (n-stress) sign: the case's SCEC spec convention, +1
        # extension-positive / -1 compression-positive, read from bGlobal.txt
        # exactly as readInputFiles.f90 reads it (board row 22a, PR #20).
        nStressOutSign=float(g['nStressOutSign']),
        st_on_idx=st_on_idx, st_on_strike_m=st_on_strike_m, st_on_depth_m=st_on_depth_m,
        st_off_idx=st_off_idx, st_off_x_m=st_off_x_m, st_off_y_m=st_off_y_m,
        st_off_z_m=st_off_z_m, st_on_total=st_on_total, st_off_total=st_off_total,
        # Row 94: the ACTUAL matched node (x,y,z), for the header location
        # stamp only -- st_off_x_m/y_m/z_m above (the REQUESTED coordinate)
        # stay what the file NAME is derived from; see offfault_filename vs
        # offfault_location_stamp (library_output.py) and this mission's
        # explicit ruling that the name stays request-derived.
        st_off_x_actual_m=st_off_actual_m[:, 0], st_off_y_actual_m=st_off_actual_m[:, 1],
        st_off_z_actual_m=st_off_actual_m[:, 2],
    )
    mesh = dict(meshCoor=meshCoor, nsmp=nsmp)
    if part is not None:
        n_local = (len(xline), len(yline), len(zline))
        mesh.update(
            n_local=n_local, offsets=offsets,
            n_global=tuple(len(a) for a in glines),
            flt_lists=meshgen.fault_boundary_lists(nsmp, *n_local),
            fault_box={k: params[k] for k in ('fxmin', 'fxmax', 'fymin', 'fymax',
                                               'fzmin', 'fzmax')},
            arn1=arn, mass1=nodalMassArr, fnms1=fnms,
            fault_census=meshgen.fault_census(*glines, params),
            equation_census=meshgen.equation_census(*glines, params, pmlb, bounds))
    return S, mesh


class Profile(dict):
    """Wall-clock per phase. Measured, never estimated.

    Why this exists: reasoning about where a run spends its time produced two
    wrong answers in a row (a threading claim built on a knob that does not
    work, and a memory-bandwidth claim contradicted by the achieved bandwidth).
    The phases below are the ones that can actually be confused for each other
    -- one-time setup vs one-time XLA compile vs the per-step loop -- and a
    single total hides all three.

    JAX dispatches asynchronously, so every phase boundary that follows device
    work calls block_until_ready; without it a phase records queue-submission
    latency and the next phase inherits the real cost.
    """

    def __init__(self, backend):
        super().__init__()
        self.backend = backend
        self.nelem = 0

    def _sync(self):
        if self.backend != 'jax':
            return
        try:
            import jax
        except ImportError:
            return
        for d in jax.devices():
            try:
                d.synchronize_all_activity()
            except AttributeError:
                pass

    @contextlib.contextmanager
    def phase(self, name):
        # t0 BEFORE the entry sync: the first _sync imports jax (~1 s), and a
        # t0 taken after it left that import outside every phase, i.e. in
        # unaccounted_s (3.1% of total_s on test.tpv8 x jax, this box).
        t0 = time.perf_counter()
        self._sync()
        try:
            yield
        finally:
            self._sync()
            self[name] = self.get(name, 0.0) + time.perf_counter() - t0

    def report(self, nsteps=None, nelem=None, stream=None):
        out = stream or sys.stderr
        total = sum(self.values())
        print('  --- profile (%s) ---' % self.backend, file=out)
        for k, v in self.items():
            share = 100.0*v/total if total else 0.0
            print('  %-22s %9.3f s  %5.1f%%' % (k, v, share), file=out)
        print('  %-22s %9.3f s' % ('TOTAL', total), file=out)
        step = self.get('solve')
        if step is not None and nsteps:
            per = step/nsteps
            line = '  %-22s %9.3f ms/step' % ('solve', per*1e3)
            if nelem:
                line += '   %8.0f ns/element/step' % (per/nelem*1e9)
            print(line, file=out)
        return self


def run_case_mpi(case_dir, comm, nsteps=None, verbose=True, profile=None):
    """run_case's one-process-per-rank twin: Fortran's decomposition, jax
    owning the local element kernel (driver.run_mpi, MPI4NodalQuant.py).

    Each rank builds only its own box (MPI4NodalQuant.DECOMP, Fortran's
    split) and writes `frt.txt<rank>` holding exactly the fault nodes it OWNS
    (the lowest rank holding each), which testsys/frt_canonical.py already
    globs for -- so an N-rank python run is compared against the SAME
    committed reference, through the same canonicalisation, with no new
    comparison path. (Fortran writes every fault node a rank holds, so its
    shared-plane nodes appear twice; this port writes each once.)

    Returns (path, report) -- the report carries this rank's element counts,
    halo size and ms/step, which every multi-rank measurement must print."""
    # total_s is an INDEPENDENT timer, not a bucket sum (profile_schema.py's
    # "the sum check, and why it is the one that matters" -- unaccounted_s
    # must be a real remainder against a total measured by its OWN clock,
    # exactly as compTimeInSeconds(9)/simuStartTime is in eqdyna3d.f90).
    #
    # EQDYNA_PROFILE read ONCE here, before setup, before the loop. When
    # off, run_t0 stays 0.0 and _emit() below is never called -- no
    # perf_counter() call this landing added runs, not just no file write.
    profile_on = _profile_emit.enabled()
    run_t0 = time.perf_counter() if profile_on else 0.0
    prof = profile if profile is not None else Profile('jax')
    from . import MPI4NodalQuant as MQ
    part = MQ.Partition.for_size(comm.Get_rank(), comm.Get_size())
    with prof.phase('setup (mesh+input)'):
        # THIS RANK'S BOX only (meshgen.f90's rank-local build), then the
        # setup exchanges Fortran does (MPI4arn, MPI4NodalQuant on mass/fnms).
        S, mesh = build_solver_state(case_dir, part=part)
        plan = MQ.setup_exchange(comm, part, S, mesh)
    # Row 114: python-jax-mpi does NOT port station output (build_solver_state
    # leaves S['st_on_idx']/['st_off_idx'] empty on this path; see driver.py's
    # run_mpi carry comment). Recorded clearly, once per job, rather than
    # silently writing no faultst*/body* files -- this must not be a hard
    # refusal: test.tpv8 (the one case opted into this mode) has stations in
    # its DEFAULT bStations.txt (scripts/defaultParameters.py's
    # st_coor_on_fault/st_coor_off_fault), so raising here would break the
    # existing gated cell rather than just leave a documented gap.
    if comm.Get_rank() == 0 and (S['st_on_total'] > 0 or S['st_off_total'] > 0):
        print('run_case_mpi: bStations.txt names %d on-fault / %d off-fault '
              'station(s), but python-jax-mpi does not port station output '
              '(row 114 scope: fortran and the serial python backends only) '
              '-- NO faultst*/body* files are written by this run.'
              % (S['st_on_total'], S['st_off_total']), file=sys.stderr, flush=True)
    with prof.phase('solve'):
        out = driver.run_mpi(S, comm, part, plan, nsteps=nsteps, verbose=verbose,
                             xp=_backend.array_module('jax'))

    rank, nranks = comm.Get_rank(), comm.Get_size()
    rep = out['report']

    def _emit(io_s):
        # See profile_emit.py's module docstring for why `fault`=0.0 (folded
        # into `element`, both Python backends) and why `wait`=0.0 here is a
        # REAL measurement (no barrier on the default path), not a gap.
        # setup = mesh+input build (Profile phase) PLUS driver.run_mpi's own
        # pre-loop invariants/to_device/jit-construction span (rep['setup_s'],
        # see driver.py's t_setup comment) -- both are real setup cost, and
        # omitting the latter is what left ~30% of total_s in
        # unaccounted_s on test.tpv8 x 4 ranks before this fix.
        buckets = dict(setup=prof.get('setup (mesh+input)', 0.0)
                            + rep.get('setup_s', 0.0),
                       element=rep.get('compute_s', 0.0), fault=0.0,
                       exchange=rep.get('mpi_s', 0.0),
                       wait=rep.get('wait_s', 0.0), io=io_s)
        total_s = time.perf_counter() - run_t0
        _profile_emit.write_profile(case_dir, 'python-jax-mpi', rank, nranks,
                                    rep['nsteps'], buckets,
                                    loop_s=rep['solve_s'], total_s=total_s)

    rows = out['fault_rows']
    sel = out['own_in_computed']
    n_own = int(rows.shape[0])
    if n_own == 0:
        # A rank that OWNS no fault node writes no frt.txt, and
        # library_output.write_frt refuses nftnd==0. Two shapes reach here:
        # a box that never meets the fault (Fortran writes nothing for it
        # either), and a box that holds fault nodes all of which a LOWER rank
        # also holds and therefore writes (e.g. the +y side of a y split lying
        # on the fault plane) -- Fortran would write duplicates there; this
        # port writes each fault row once so frt_canonical's duplicate branch
        # never fires. This cannot hide lost nodes: setup_exchange has already
        # checked the owned (count, key sum) against the global fault census,
        # so a missing file means "this rank owned none", never "these nodes
        # went missing". Which ranks own none depends on the decomposition --
        # test at the rank counts that change it, not at 2.
        # It still gets a profile.rank<r>.json -- the parity gate's byte-
        # identity check is about frt output, not about profile coverage,
        # and a rank that did real setup/compute/exchange work is not "no
        # data" just because it wrote no fault row.
        prof.nelem = rep['E']
        if profile_on:
            _emit(io_s=0.0)
        return None, rep
    fric_1idx = np.zeros((n_own + 1, 101))
    fric_1idx[1:, 1:101] = out['fric'][sel]
    fnft_1idx = np.zeros(n_own + 1)
    fnft_1idx[1:] = out['fnft'][sel]
    path = os.path.join(case_dir, 'frt.txt%d' % comm.Get_rank())
    with prof.phase('write frt'):
        library_output.write_frt(path, mesh['meshCoor'], mesh['nsmp'][rows],
                                 fnft_1idx, fric_1idx)
    prof.nelem = rep['E']
    if profile_on:
        _emit(io_s=prof.get('write frt', 0.0))
    return path, rep


def run_case(case_dir, nsteps=None, verbose=True, backend=DEFAULT_BACKEND,
             profile=None):
    """Builds S (zero pydump reads), runs driver.py's one time loop under the
    requested `backend` ('jax', the default, or 'numpy') with the friclaw
    dispatch inside faulting.py, and writes frt.txt0 via
    library_output.write_frt (byte-exact Fortran E18.7E4 format). Returns the
    path written.

    `profile` is an optional Profile; when given, each phase is timed
    separately so setup, solve and output cannot be confused for one another.
    """
    # EQDYNA_PROFILE read ONCE here, before setup, before the loop -- never
    # per step. When off, run_t0 stays 0.0 and the write_profile call below
    # is skipped entirely (no bucket dict built, no total_s taken) -- not
    # just gated at the file-write step inside write_profile itself.
    profile_on = _profile_emit.enabled()
    run_t0 = time.perf_counter() if profile_on else 0.0   # independent total_s timer, see run_case_mpi
    if backend == 'numpy':
        _narrow_numpy_affinity()
    prof = profile if profile is not None else Profile(backend)
    with prof.phase('setup (mesh+input)'):
        S, mesh = build_solver_state(case_dir)
    with prof.phase('resolve solver'):
        solver = _resolve_solver(S['friclaw'], backend)
    with prof.phase('solve'):
        out = solver(S, nsteps=nsteps, verbose=verbose)

    nftnd = S['nftnd']
    fric_1idx = np.zeros((nftnd + 1, 101))
    fric_1idx[1:, 1:101] = out['fric']
    fnft_1idx = np.zeros(nftnd + 1)
    fnft_1idx[1:] = out['fnft']

    # A run made under one of backend.py's timing-only sharding knobs computes
    # the wrong answer by construction, so it must not be able to land on the
    # path every comparison tool reads. It is still written (the measurement
    # wants a completed run, and a silently skipped write is its own defect),
    # under a name nothing gates on.
    name = ('frt.txt0.TIMING-ONLY-INVALID' if _backend.timing_only()
            else 'frt.txt0')
    frt_path = os.path.join(case_dir, name)
    with prof.phase('write frt'):
        library_output.write_frt(frt_path, mesh['meshCoor'], mesh['nsmp'],
                             fnft_1idx, fric_1idx)
    with prof.phase('write stations'):
        # Same TIMING-ONLY guard as frt_path above: a run made under one of
        # backend.py's timing-only sharding knobs computes the wrong answer
        # by construction and must not silently overwrite the real
        # faultst*/body* files (which carry no distinguishing suffix the way
        # frt.txt0.TIMING-ONLY-INVALID does).
        if _backend.timing_only():
            print('run_case: *** TIMING-ONLY RUN -- NOT writing faultst*/body* '
                  '(would overwrite valid station files under their real names) ***',
                  file=sys.stderr)
        else:
            library_output.write_onfault_stations(case_dir, S, out['on_st_hist'])
            library_output.write_offfault_stations(case_dir, S, out['off_st_hist'])
    prof.nelem = S.get('totalNumOfElements') or 0

    # ALWAYS-ON profile.rank0.json (nranks=1: this path is serial by
    # construction, checkInputConsistency/build_solver_state already refuse
    # npx/npy/npz>1). `exchange`/`wait` are genuinely 0.0: there is no MPI on
    # this path at all, not a folded or unmeasured cost.
    #
    # `fault`: split out of 'solve' along Fortran's own boundary (`fault` =
    # compTimeInSeconds(6), faulting.f90:28 -- rule 23) for python-numpy
    # ONLY. `driver.run` measures it with a plain perf_counter() pair around
    # the eagerly-executed FLT.faulting call (see driver.make_step_parts's
    # docstring) and returns it as `out['fault_s']`, always 0.0 on jax
    # because make_step_parts refuses to create the timer at all when
    # `B.is_jax(xp)` -- a timer inside a traced function would measure the
    # ONE-TIME trace, not the per-step cost (wrong, not just imprecise), so
    # python-jax and python-jax-mpi stay folded (`fault`=0.0, cost inside
    # `element`) exactly as before. See profile_emit.py's docstring.
    if profile_on:
        fault_s = out.get('fault_s', 0.0)
        solve_s = prof.get('solve', 0.0)
        total_s = time.perf_counter() - run_t0
        _profile_emit.write_profile(
            case_dir, 'python-%s' % backend, 0, 1, S['nstep'],
            serial_buckets(prof, fault_s),
            loop_s=solve_s, total_s=total_s)
    return frt_path


# Every Profile phase run_case records, and the schema bucket it belongs to.
# `io` is frt AND station output: Fortran's io figure (compTimeInSeconds(8),
# eqdyna3d.f90:182) spans the whole post-loop output stage, stations included
# (rule 23). The 'write stations' phase (row 114) was once left out of this
# map, so its time fell into unaccounted_s and tripped the schema's 5% sum
# check on master CI run 36082084976 (0.80 s of 15.40 s, test.tpv8 x jax).
_SERIAL_PHASE_BUCKET = {
    'setup (mesh+input)': 'setup',
    'resolve solver': 'setup',
    'solve': 'element',
    'write frt': 'io',
    'write stations': 'io',
}


def serial_buckets(prof, fault_s):
    """run_case's Profile phases -> the six profile_emit buckets. Refuses a
    phase with no bucket rather than letting its time vanish into
    unaccounted_s; `fault_s` (numpy only, else 0.0) is carved out of solve."""
    unmapped = sorted(set(prof) - set(_SERIAL_PHASE_BUCKET))
    if unmapped:
        raise ValueError('serial_buckets: profile phase(s) %s have no bucket in '
                         '_SERIAL_PHASE_BUCKET -- map each one, or its time is '
                         'silently counted as unaccounted_s' % unmapped)
    buckets = dict.fromkeys(_profile_emit.BUCKET_KEYS, 0.0)
    for name, secs in prof.items():
        buckets[_SERIAL_PHASE_BUCKET[name]] += secs
    buckets['element'] -= fault_s
    buckets['fault'] = fault_s
    return buckets


def _select_device(device):
    """Pin the JAX platform BEFORE jax is imported anywhere. No fallback.

    JAX_PLATFORMS is read by jax at import time, so this must run before the
    first `import jax` -- which is why this module never imports jax at module
    level and `active_device`/`Profile._sync` import it inside the function.

    `--device cuda` with no GPU is a hard failure: silently running on CPU would
    put a row labelled gpu into a backend comparison whose whole purpose is to
    tell cpu and gpu apart (rule 2).
    """
    os.environ['JAX_PLATFORMS'] = device
    import jax
    got = jax.devices()[0].platform
    want = 'gpu' if device == 'cuda' else 'cpu'
    if got != want:
        raise RuntimeError(
            '--device %s was requested but JAX resolved to %r (devices: %r). '
            'This does NOT fall back: a run reported as %s must have run on %s.'
            % (device, got, jax.devices(), device, device))


def _abort(exc, rank=0):
    """checkInputConsistency.InputConsistencyError -> the FATAL block
    errorCodes.f90:153-162 (abortRun) prints: on STDOUT, with the rank line
    (Fortran always has MPI up, so a serial run is rank 0), then
    SystemExit(exc.code) so a Python run and a Fortran run of the same bad
    config exit with the SAME number (rule 23). Every rank raises the same
    refusal before any collective (build_solver_state runs first on all
    ranks), so THIS path needs no MPI_Abort. Any OTHER exception raised on
    one rank goes through the MPI entry point's traceback + comm.Abort
    handler instead, so its partners cannot hang."""
    print(flush=True)
    print(' ==================== EQdyna: FATAL ====================')
    print('  rank      : ', rank)
    print('  exit code : ', exc.code)
    print('  reason    : ', exc)
    print('  See the "Exit codes" table in README.md for this code.')
    print(' =======================================================')
    print(flush=True)
    raise SystemExit(exc.code)


DEVICE_CHOICES = ('cpu', 'cuda')


def _device_arg(value):
    """--device's parser. `auto` was removed by owner ruling (board row 56,
    2026-09-24): it let JAX pick a GPU on a serial run while --mpi silently
    mapped it to cpu. Refused by name, never re-mapped. The GPU choice is
    spelled `cuda` (the owner's wording, and the spelling run_e2e/run.py
    already used); the old `gpu` spelling is refused with that hint."""
    if value not in DEVICE_CHOICES:
        hint = {'auto': " -- 'auto' was removed (board row 56); pass --device cuda to request a GPU",
                'gpu': " -- the GPU choice is spelled 'cuda' (board row 56)"}.get(value, '')
        raise argparse.ArgumentTypeError(
            '%r is not a device: choose cpu (the default) or cuda%s' % (value, hint))
    return value


def main():
    ap = argparse.ArgumentParser(prog='python3 -m eqdyna')
    ap.add_argument('case_dir')
    ap.add_argument('nsteps', nargs='?', type=int, default=None)
    ap.add_argument('--backend', choices=('jax', 'numpy'), default=DEFAULT_BACKEND,
                     help='solver backend (default: %(default)s). No fallback: '
                          'jax with jaxlib missing is an error, not a demotion '
                          'to numpy.')
    ap.add_argument('--device', type=_device_arg, default='cpu',
                     help='JAX platform: cpu (default) or cuda. GPU runs only '
                          'on an explicit --device cuda, serial and --mpi '
                          'alike. No fallback: gpu with no GPU is an error.')
    ap.add_argument('--profile', action='store_true',
                     help='print wall-clock per phase (setup / solve / write) '
                          'so one-time cost and per-step cost cannot be '
                          'confused. JAX phases are synchronised, so the '
                          'numbers are compute, not queue submission.')
    ap.add_argument('--mpi', action='store_true',
                     help='one process per rank, Fortran-style domain '
                          'decomposition, jax owning the local element kernel '
                          '(driver.run_mpi). Launch under mpirun. Writes '
                          'frt.txt<rank>. No fallback: --mpi with mpi4py '
                          'missing is an error, and --mpi outside mpirun is a '
                          'legal 1-rank run, not a silent serial demotion.')
    args = ap.parse_args()
    if args.mpi:
        if args.backend != 'jax':
            raise SystemExit('--mpi is implemented for --backend jax only '
                             '(numpy is out of scope for the MPI path)')
        _select_device(args.device)
        from mpi4py import MPI      # ImportError is deliberate, not caught
        comm = MPI.COMM_WORLD
        prof = Profile('jax')
        try:
            path, report = run_case_mpi(args.case_dir, comm, nsteps=args.nsteps,
                                        profile=prof)
        except checkInputConsistency.InputConsistencyError as exc:
            _abort(exc, rank=comm.Get_rank())
        except BaseException:
            # A failure on ONE rank -- its own box's mesh, its own faces --
            # while the others wait in a Sendrecv would otherwise hang the job
            # forever (observed: a fault-free box raised in
            # build_node_coordinates and three ranks sat 27 min in the setup
            # exchange). Fortran's abortRun calls MPI_Abort for the same
            # reason. The traceback is printed first, flushed, so the cause is
            # attributable to a rank and a line.
            import traceback
            traceback.print_exc()
            sys.stdout.flush(); sys.stderr.flush()
            print('rank %d/%d: aborting the MPI job (see traceback above)'
                  % (comm.Get_rank(), comm.Get_size()), file=sys.stderr, flush=True)
            comm.Abort(1)
        if args.profile:
            prof.report(nsteps=args.nsteps, nelem=prof.nelem or None,
                        stream=sys.stdout)
        # `path` is None for a rank that owns no fault node (see
        # run_case_mpi): it wrote nothing, exactly as Fortran does, and the
        # line says so rather than printing a bare None that reads like a bug.
        # backend=/device= on EVERY rank, the same line the serial path
        # prints, because active_device's own contract is that a timing row
        # carries what ACTUALLY ran -- and this branch was the one place that
        # printed no device at all. Measured consequence: a 1-rank tpv104
        # point launched with JAX_PLATFORMS=cuda recorded 812.78 ms/step with
        # zero bytes allocated on any GPU, because this branch's
        # `--device auto` mapped to cpu here and nothing downstream
        # said so (that silent auto->cpu map was removed with `auto` itself,
        # board row 56; the default is now an explicit cpu). One line per rank; note that under a per-rank
        # CUDA_VISIBLE_DEVICES every rank legitimately reports gpu:0, so this
        # line proves the PLATFORM, and it is the per-device memory poll in
        # testsys/perf/run_mpi_scaling.py that proves four distinct devices.
        print('backend=%s device=%s' % (args.backend, active_device(args.backend)))
        print('rank %d/%d wrote %s  %s'
              % (comm.Get_rank(), comm.Get_size(),
                 path or 'NO-FRT-OWNS-0-FAULT-NODES',
                 ' '.join('%s=%s' % (k, report[k]) for k in
                          ('Ei', 'Ep', 'halo_eqs', 'ms_per_step',
                           'mpi_ms_per_step', 'wait_ms_per_step', 'sync',
                           'device_peak_gb', 'threads', 'cpus_allowed',
                           # The per-rank working set, on the same line as the
                           # per-step cost it explains. A scaling number
                           # without it cannot distinguish "the work shrank"
                           # from "the problem shrank"; these two rose and
                           # fell together is the whole claim of the rank-local
                           # carry, so they are recorded together by tooling.
                           'N_local', 'NEQ_local', 'carry_bytes_total'))))
        return
    if args.backend == 'jax':
        _select_device(args.device)
    elif args.device == 'cuda':
        raise SystemExit('--device cuda is meaningless with --backend numpy')
    prof = Profile(args.backend)
    try:
        path = run_case(args.case_dir, nsteps=args.nsteps, backend=args.backend,
                        profile=prof)
    except checkInputConsistency.InputConsistencyError as exc:
        _abort(exc)
    print('backend=%s device=%s' % (args.backend, active_device(args.backend)))
    if args.profile:
        prof.report(nsteps=args.nsteps, nelem=prof.nelem or None,
                    stream=sys.stdout)
    print('wrote', path)


if __name__ == '__main__':
    main()


def run(S, nsteps=None, verbose=True, backend='numpy'):
    """Step the solve on `backend` ('numpy' or 'jax').

    This is eqdyna3d.f90's call to driver, with the ONE argument the Fortran
    does not have. There is no per-backend and no per-friclaw solver module
    to pick between any more: driver.py is the single time loop and
    faulting.py dispatches friclaw inside it, exactly as faulting.f90:21-22
    does.
    """
    return driver.run(S, nsteps=nsteps, verbose=verbose,
                      xp=_backend.array_module(backend))
