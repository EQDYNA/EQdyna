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
partial runs: ntotft==1, serial (npx==npy==npz==1), C_degen==0 (planar) OR
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

from . import assembleGlobalMass, driver, func_lib, library_output, meshgen, readInputFiles
from . import backend as _backend

# friclaw -> the NumPy solver module whose run(S, nsteps, verbose) consumes
# this module's S dict. All three modules were independently verified
# (by the removed parity tier) to need the EXACT SAME S-dict keys (the
# per-friclaw modules
# port_tp.py only add nucleation-parameter reads and, for TP, an internally
# -owned onFaultTPHist scan-carry -- neither needs a NEW key from S beyond
# what build_solver_state already provides for tpv8, confirmed by grepping
# each module's `S[...]` accesses before wiring this dispatch) -- so
# build_solver_state below is friclaw-agnostic; only the solver CALLED
# differs.
# The friction laws this port implements. EVERY one of them is served by the
# SAME code -- eqdyna/{driver,faulting,fric,assembleGlobalKU,backend}.py --
# with the friclaw dispatch inside faulting.py exactly where faulting.f90:21-22
# puts it, and with the backend as an argument rather than a second module.
#
# This used to be two tables of three modules each: driver.py
# time loop. They are deleted.
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
    """Returns the run()-providing module for `friclaw` under `backend`.

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


def build_solver_state(case_dir):
    """Builds the full `S` dict eqdyna/driver.py's `run()` expects,
    reading ONLY case-input files (bGlobal.txt/bModelGeometry.txt/
    bFaultGeometry.txt/bMaterial.txt/on_fault_vars_input.nc) via
    readInputFiles.py, then meshgen.py/assembleGlobalMass.py -- no pydump_*.txt
    anywhere. Returns (S, mesh) where mesh is a dict of the raw 1-indexed
    builder outputs (meshCoor, nsmp, eq_nums, ...) library_output needs later
    (S itself is converted to loading.load()'s 0-indexed convention).
    """
    params, g = readInputFiles.build_params(case_dir)
    if g['ntotft'] != 1:
        raise NotImplementedError('build_solver_state: only ntotft==1 is supported (got %d)'
                                   % g['ntotft'])
    if g['friclaw'] not in SUPPORTED_FRICLAW:
        raise NotImplementedError('build_solver_state: friclaw=%d is not implemented '
                                   '(implemented: %r)' % (g['friclaw'], list(SUPPORTED_FRICLAW)))
    if (g['npx'], g['npy'], g['npz']) != (1, 1, 1):
        raise NotImplementedError('build_solver_state: only serial (npx=npy=npz=1) is supported '
                                   '(got %r)' % ((g['npx'], g['npy'], g['npz']),))
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

    xline, yline, zline, pmlb, _ = meshgen.build_grid_lines(params)
    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)
    conn, elem_type, mat, elem_depth = meshgen.build_elements(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    num_dof, eq_start, eq_nums, total_eqs = meshgen.build_equation_numbers(
        xline, yline, zline, params, pmlb)
    un, us, ud, arn = meshgen.build_fault_geometry(xline, yline, zline, params, nsmp)

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
    )
    mesh = dict(meshCoor=meshCoor, nsmp=nsmp)
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
        self._sync()
        t0 = time.perf_counter()
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


def run_case(case_dir, nsteps=None, verbose=True, backend=DEFAULT_BACKEND,
             profile=None):
    """Builds S (zero pydump reads), dispatches to the friclaw-appropriate
    solver's run() under the requested `backend` ('jax', the default, or
    'numpy') -- port.py/port_jax.py friclaw==1, port_rsf.py/port_rsf_jax.py
    friclaw==4, port_tp.py/port_tp_jax.py friclaw==5 -- writes frt.txt0 via
    library_output.write_frt (byte-exact Fortran E18.7E4 format). Returns the
    path written.

    `profile` is an optional Profile; when given, each phase is timed
    separately so setup, solve and output cannot be confused for one another.
    """
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

    frt_path = os.path.join(case_dir, 'frt.txt0')
    with prof.phase('write frt'):
        library_output.write_frt(frt_path, mesh['meshCoor'], mesh['nsmp'],
                             fnft_1idx, fric_1idx)
    prof.nelem = S.get('totalNumOfElements') or 0
    return frt_path


def _select_device(device):
    """Pin the JAX platform BEFORE jax is imported anywhere. No fallback.

    JAX_PLATFORMS is read by jax at import time, so this must run before the
    first `import jax` -- which is why this module never imports jax at module
    level (see the note above _JAX_MODULE_BY_FRICLAW).

    `--device gpu` with no GPU is a hard failure: silently running on CPU would
    put a row labelled gpu into a backend comparison whose whole purpose is to
    tell cpu and gpu apart (rule 2).
    """
    if device == 'auto':
        return
    os.environ['JAX_PLATFORMS'] = {'cpu': 'cpu', 'gpu': 'cuda'}[device]
    import jax
    got = jax.devices()[0].platform
    want = 'gpu' if device == 'gpu' else 'cpu'
    if got != want:
        raise RuntimeError(
            '--device %s was requested but JAX resolved to %r (devices: %r). '
            'This does NOT fall back: a run reported as %s must have run on %s.'
            % (device, got, jax.devices(), device, device))


def main():
    ap = argparse.ArgumentParser(prog='python3 -m eqdyna')
    ap.add_argument('case_dir')
    ap.add_argument('nsteps', nargs='?', type=int, default=None)
    ap.add_argument('--backend', choices=('jax', 'numpy'), default=DEFAULT_BACKEND,
                     help='solver backend (default: %(default)s). No fallback: '
                          'jax with jaxlib missing is an error, not a demotion '
                          'to numpy.')
    ap.add_argument('--device', choices=('auto', 'cpu', 'gpu'), default='auto',
                     help='JAX platform (default: auto = whatever JAX picks). '
                          'Only meaningful with --backend jax. No fallback: '
                          'gpu with no GPU is an error.')
    ap.add_argument('--profile', action='store_true',
                     help='print wall-clock per phase (setup / solve / write) '
                          'so one-time cost and per-step cost cannot be '
                          'confused. JAX phases are synchronised, so the '
                          'numbers are compute, not queue submission.')
    args = ap.parse_args()
    if args.backend == 'jax':
        _select_device(args.device)
    elif args.device != 'auto':
        raise SystemExit('--device %s is meaningless with --backend numpy'
                         % args.device)
    prof = Profile(args.backend)
    path = run_case(args.case_dir, nsteps=args.nsteps, backend=args.backend,
                    profile=prof)
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
