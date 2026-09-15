"""
Milestone 8: the standalone entry point. Assembles the solver state dict
`S` ENTIRELY from python/eqdyna/standalone/{native_input,meshgen,
mass_assembly}.py -- ZERO pydump_*.txt reads (grep this file for
'pydump': there are none outside this docstring) -- and hands it to the
EXISTING, already-parity-verified NumPy time-stepping kernel
(python/eqdyna/port.py's `run()`, which wraps kernels_numpy.py's
velDispUpdate + assembleGlobalKU + hrglss + faulting), then writes
frt.txt via frt_writer.write_frt.

Scope: friclaw in {1, 4, 5} (slip-weakening/tpv8, rate-and-state/tpv104,
thermal-pressurization/tpv1053d -- dispatched to port.py/port_rsf.py/
port_tp.py respectively, see _NUMPY_SOLVER_BY_FRICLAW/_JAX_MODULE_BY_FRICLAW
below), single planar or dipping/rough fault (ntotft==1), C_degen==0,
npx==npy==npz==1 (serial). insertFaultType>0 (tpv10 dipping, drv.a6
fractal-rough -- Milestone 9, meshgen.py's insert_fault_interface +
build_fault_geometry's pfx/pfz un/us/ud branch) is verified combined with
friclaw==1 (tpv10) and friclaw==4 (drv.a6, faulting.f90's min_norm/max_norm
normal-stress clamp ported into port_rsf.py/port_rsf_jax.py); friclaw==5
(TP) with insertFaultType>0 still needs that same clamp ported into
port_tp.py/port_tp_jax.py, not done yet (see the guard below).
Raises loudly if a case violates any of these (no silent partial run).

S-dict field-by-field provenance (every key python/eqdyna/loading.py's
`load()` reads from a pydump_* file, this module instead computes from
the M1-M7.5 builders -- same shapes/conventions, so kernels_numpy.py/
port.py need ZERO changes to consume it):

  N, E, NEQ, nen, ned            -- meshgen.py M1-M3 (nen=8, ned=3: the
                                     Fortran globalvar.f90 hex-element/
                                     dof-per-node constants, never
                                     case-input, reproduced verbatim)
  dt, w, rdampk, rdampm, kapa_hg,
  R, nPML, vmaxPML, PMLb, grav,
  C_elastic, roumax, rhow, gamar,
  slipRateThres, xsource, ysource,
  zsource, nucR, nucRuptVel,
  nucdtau0, nucT, TPV, C_nuclea,
  nucfault, friclaw, nstep         -- native_input.read_bglobal (bGlobal.txt)
                                     + globalvar.f90 PARAMETER constants
                                     (w=8.0 calcLocalShapeFunc.f90; rdampm=0.0,
                                     kapa_hg=0.1, R=0.01, nPML=6, grav=9.8 --
                                     see native_input.build_params/this
                                     module for where each lives)
  meshCoor, conn, elemType, mat    -- meshgen.py M1/M2 (build_node_coordinates/
                                     build_elements), converted from this
                                     port's 1-indexed convention to
                                     loading.load()'s 0-indexed convention
  eledet, eleshp, ss, phi          -- mass_assembly.py M7.5
                                     (compute_element_shape/compute_hourglass)
  nodalMassArr, fnms               -- mass_assembly.py M6.5 (assemble_mass)
  ndof, eq_ids                     -- meshgen.py M3 (build_equation_numbers)
  nsmp1, nsmp2, un, us, ud, arn    -- meshgen.py M4 (build_fault_geometry)
  fric_init                        -- native_input.py M6 (read_on_fault_vars)
  v1_init                          -- mass_assembly.py M7.5 (init_vel) --
                                     NOTE (pre-existing, not introduced
                                     here): port.py's run() calls
                                     kernels_numpy.init_state, which
                                     zero-inits v1 and never reads
                                     S['v1_init'] -- the SAME gap exists in
                                     loading.load() (it computes v1_init
                                     but run() ignores it). For friclaw==1/
                                     mode==1 this is a documented no-op
                                     (fric(31:36) are exactly 0.0, verified
                                     by M7.5's init_vel parity check), so it
                                     does not affect tpv8's result -- flagged
                                     here, not silently worked around,
                                     because a future case where mode==2 or
                                     fric(31:36)!=0 WOULD need this wired in.

Milestone 10 (drv.a6, C_elastic==0, friclaw==4): wires
src/calcElemKU.f90:127-161's Drucker-Prager viscoplastic return-mapping
(deviatoric projection to a cohesion/friction-angle yield surface, tv-scale
viscoplastic relaxation) and src/assembleGlobalKU.f90:16's
`(1-C_elastic)*grav*(roumax-(gamar+1)*rhow)/roumax` gravity body-force term
into kernels_numpy.py/kernels_jax.py's shared elastic_step (gated behind a
static `if S['C_elastic']==0:` branch -- an exact, zero-cost no-op for every
other case: the gravity term's own `(1-C_elastic)` factor is EXACTLY 0.0 for
C_elastic==1, and the plasticity block is skipped entirely, not merely
masked, so tpv8/tpv104/tpv1053d/tpv10 parity is provably untouched -- see
kernels_numpy.py's build()/elastic_step docstrings for the derivation of why
the gravity term reduces to `-nodalMassArr[z-eq]*const` and why
src/meshgen.f90:103's setPlasticStress lithostatic pre-stress is
precomputed once here (`ccosphi`/`sinphi`/`tv`/`init_stress`) rather than
per-step). eleporep (pore pressure) is hardcoded 0.0 in the plasticity
kernel: confirmed by grep that src/meshgen.f90's setPlasticStress and
src/eqdyna3d.f90's zero-init are the ONLY writes to `eleporep` in all of
src/*.f90 -- the lithostatic pore-pressure formula is commented out in the
Fortran itself, so eleporep is provably always exactly 0.0, not an
assumption. pstrain (the plastic-strain accumulator) is explicitly NOT
ported: test.reference.results/test.drv.a6 has no pstr.txt* output to gate
against (only frt.txt1/frt.txt3), and pstrain is write-only (never read back
by any downstream physics), so computing it would be unverifiable, silent
scope creep -- flagged here, not silently included.

Verified (tpv8, friclaw==1): coordinate-aligned comparison vs the
committed 4-rank test.reference.results/test.tpv8 references (frt.txt0 +
frt.txt2, deduped by rounded coordinates and lexsorted against this
module's serial output, since the 4-rank references use a different
domain decomposition and therefore a different frt.txt node order than a
serial run -- a flat positional text diff cannot compare the two
directly). See the M7.5/M8 landing commit message for the exact
methodology and numbers; extended to tpv104 (friclaw==4) and tpv1053d
(friclaw==5) the same way, same-session follow-on.
"""
import argparse
import contextlib
import importlib
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from eqdyna.standalone import meshgen, native_input, mass_assembly, frt_writer  # noqa: E402
from eqdyna import port, port_rsf, port_tp  # noqa: E402

# friclaw -> the NumPy solver module whose run(S, nsteps, verbose) consumes
# this module's S dict. All three modules were independently verified
# (README-parity.md) to need the EXACT SAME S-dict keys (port_rsf.py/
# port_tp.py only add nucleation-parameter reads and, for TP, an internally
# -owned onFaultTPHist scan-carry -- neither needs a NEW key from S beyond
# what build_solver_state already provides for tpv8, confirmed by grepping
# each module's `S[...]` accesses before wiring this dispatch) -- so
# build_solver_state below is friclaw-agnostic; only the solver CALLED
# differs.
_NUMPY_SOLVER_BY_FRICLAW = {1: port, 4: port_rsf, 5: port_tp}

# JAX counterparts, same S-dict, same run(S, nsteps, verbose) signature,
# same parity-verified formulas (README-parity.md Updates 3/4/5) -- 2-2.6x
# faster than serial Fortran on this hardware (testsys/perf/baseline.json).
# Imported LAZILY (module names, not modules) so a numpy-only environment
# with no jaxlib installed never pays an ImportError just for choosing
# --backend numpy; each jax module itself does
# `jax.config.update("jax_enable_x64", True)` as the FIRST line after
# `import jax`, before `import jax.numpy` -- that ordering guarantee is
# preserved here because this file never imports jax/jax.numpy directly,
# only these modules, lazily, on first actual use.
_JAX_MODULE_BY_FRICLAW = {1: 'eqdyna.port_jax', 4: 'eqdyna.port_rsf_jax', 5: 'eqdyna.port_tp_jax'}

DEFAULT_BACKEND = 'jax'


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
    if friclaw not in _NUMPY_SOLVER_BY_FRICLAW:
        raise NotImplementedError('_resolve_solver: friclaw=%d has no wired standalone solver '
                                   '(wired: %r)' % (friclaw, sorted(_NUMPY_SOLVER_BY_FRICLAW)))
    if backend == 'numpy':
        return _NUMPY_SOLVER_BY_FRICLAW[friclaw]
    if backend != 'jax':
        raise ValueError("_resolve_solver: backend must be 'jax' or 'numpy' (got %r)" % backend)
    try:
        return importlib.import_module(_JAX_MODULE_BY_FRICLAW[friclaw])
    except ImportError as e:
        raise RuntimeError(
            "backend 'jax' was requested but %s is not importable (%s). "
            "Install jaxlib, or pass --backend numpy explicitly. This does NOT "
            "fall back: a run reported as jax must have been jax."
            % (_JAX_MODULE_BY_FRICLAW[friclaw], e))


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
    """Builds the full `S` dict python/eqdyna/port.py's `run()` expects,
    reading ONLY case-input files (bGlobal.txt/bModelGeometry.txt/
    bFaultGeometry.txt/bMaterial.txt/on_fault_vars_input.nc) via
    native_input.py, then meshgen.py/mass_assembly.py -- no pydump_*.txt
    anywhere. Returns (S, mesh) where mesh is a dict of the raw 1-indexed
    builder outputs (meshCoor, nsmp, eq_nums, ...) frt_writer needs later
    (S itself is converted to loading.load()'s 0-indexed convention).
    """
    params, g = native_input.build_params(case_dir)
    if g['ntotft'] != 1:
        raise NotImplementedError('build_solver_state: only ntotft==1 is supported (got %d)'
                                   % g['ntotft'])
    if g['friclaw'] not in _NUMPY_SOLVER_BY_FRICLAW:
        raise NotImplementedError('build_solver_state: friclaw=%d has no wired standalone '
                                   'solver (wired: %r)' % (g['friclaw'], sorted(_NUMPY_SOLVER_BY_FRICLAW)))
    if (g['npx'], g['npy'], g['npz']) != (1, 1, 1):
        raise NotImplementedError('build_solver_state: only serial (npx=npy=npz=1) is supported '
                                   '(got %r)' % ((g['npx'], g['npy'], g['npz']),))
    if g['insertFaultType'] != 0 and g['friclaw'] not in (1, 4):
        # Milestone 9 (dipping/rough fault, meshgen.py's insert_fault_interface
        # + build_fault_geometry's pfx/pfz un/us/ud branch) ported and
        # verified for friclaw==1 (tpv10). friclaw==4 (RSF, drv.a6) is now
        # also verified: faulting.f90's solveRSF has its own
        # min_norm/max_norm normal-stress clamp (line ~201, "if
        # (insertFaultType>0 .and. C_elastic==1)"), ported verbatim into
        # port_rsf.py/port_rsf_jax.py (see S['insertFaultType'] below).
        # friclaw==5 (TP, port_tp.py/port_tp_jax.py) hits the SAME clamp in
        # the Fortran but has NOT had it ported -- still guarded loudly.
        raise NotImplementedError(
            'build_solver_state: insertFaultType=%d combined with friclaw=%d is not yet '
            'ported -- verified combinations are friclaw==1 (tpv10) and friclaw==4 (drv.a6); '
            'friclaw==5 additionally needs faulting.f90\'s min_norm/max_norm clamp ported '
            'into port_tp.py/port_tp_jax.py' % (g['insertFaultType'], g['friclaw']))

    material = native_input.read_bmaterial(
        os.path.join(case_dir, 'bMaterial.txt'), g['nmat'], g['n2mat'])

    xline, yline, zline, pmlb, _ = meshgen.build_grid_lines(params)
    meshCoor, nftnd, nsmp = meshgen.build_node_coordinates(xline, yline, zline, params)
    conn, elem_type, mat, elem_depth = meshgen.build_elements(
        xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
    num_dof, eq_start, eq_nums, total_eqs = meshgen.build_equation_numbers(
        xline, yline, zline, params, pmlb)
    un, us, ud, arn = meshgen.build_fault_geometry(xline, yline, zline, params, nsmp)

    fric = native_input.read_on_fault_vars(
        os.path.join(case_dir, 'on_fault_vars_input.nc'), params['fxmin'], params['fzmin'],
        params['dx'], params['dz'], meshCoor, nsmp)

    xl = meshCoor[conn]
    det, eleshp3, xs = mass_assembly.compute_element_shape(xl)
    ss, phi48 = mass_assembly.compute_hourglass(xl, xs, mat, eleshp3)
    nodalMassArr, fnms = mass_assembly.assemble_mass(
        conn, mat, det, num_dof, eq_start, eq_nums, total_eqs, meshCoor.shape[0] - 1)
    v1 = mass_assembly.init_vel(nsmp, eq_nums, fric, total_eqs)

    N = meshCoor.shape[0] - 1  # drop the unused row-0
    E = conn.shape[0]

    # Milestone 10 (drv.a6, C_elastic==0): readInputFiles.f90's readmaterial
    # computes ccosphi/sinphi/tv AFTER bMaterial.txt is read (needs dz, read
    # earlier in readmodelgeometry) but BEFORE the time loop -- reproduced
    # here, verbatim formula, NUC_VS_FIXED=3464.0 (globalvar.f90's own
    # comment: "fixed shear-wave speed ... also readInputFiles.f90's tv
    # init"). g['bulk']/g['coheplas'] are read unconditionally by
    # native_input.read_bglobal regardless of C_elastic (same bGlobal.txt
    # line for every case) -- always computed here too, harmless for
    # C_elastic==1 cases since kernels_numpy/kernels_jax gate their USE on
    # C_elastic==0 (see those modules' build()).
    _NUC_VS_FIXED = 3464.0
    ccosphi = g['coheplas'] * np.cos(np.arctan(g['bulk']))
    sinphi = np.sin(np.arctan(g['bulk']))
    tv = 2.0 * params['dz'] / _NUC_VS_FIXED

    # meshgen.f90:103's setPlasticStress (called for EVERY element, both
    # interior and PML, only when C_elastic==0): lithostatic per-element
    # pre-stress, Voigt order [xx,yy,zz,yz,xz,xy] (calcB.f90's b(4,*)/
    # b(5,*)/b(6,*) confirm 4=yz,5=xz,6=xy) -- ALWAYS computed (cheap,
    # harmless when unused) so kernels_numpy.build/kernels_jax's build can
    # raise loudly (not silently default) if C_elastic==0 but this key is
    # somehow missing, rather than silently zero-initializing plastic runs.
    strVert = -(g['roumax'] - g['rhow'] * (g['gamar'] + 1.0)) * elem_depth * 9.8
    devStr = np.abs(strVert) * g['devStrToStrVertRatio']
    theta2 = 2.0 * g['str1ToFaultAngle']
    init_stress = np.zeros((E, 6))
    init_stress[:, 0] = strVert - devStr * np.cos(theta2)  # xx
    init_stress[:, 1] = strVert + devStr * np.cos(theta2)  # yy
    init_stress[:, 2] = strVert  # zz
    init_stress[:, 5] = devStr * np.sin(theta2)  # xy

    # ---- convert to loading.load()'s 0-indexed convention ----
    eq_ids = np.zeros((N, 12), dtype=np.int64)
    ndof0 = np.zeros(N, dtype=np.int64)
    for node in range(1, N + 1):
        nd = int(num_dof[node])
        ndof0[node - 1] = nd
        eqs = np.array(eq_nums[node], dtype=np.int64)
        eqs = np.where(eqs > 0, eqs, 0)  # -1 fixed-boundary sentinel -> sink 0
        eq_ids[node - 1, :nd] = eqs

    PMLb = np.array([pmlb['xmax0'], pmlb['xmin0'], pmlb['ymax0'], pmlb['ymin0'], pmlb['zmin0'],
                      pmlb['maxdx'], pmlb['maxdy'], pmlb['maxdz']])

    S = dict(
        N=N, E=E, NEQ=total_eqs, nen=8, ned=3, nftnd=nftnd, ntotft=1,
        nstep=g['nstep'], dt=g['dt'], w=mass_assembly._W, rdampk=g['rdampk'], rdampm=0.0,
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
    frt_writer.write_frt (byte-exact Fortran E18.7E4 format). Returns the
    path written.

    `profile` is an optional Profile; when given, each phase is timed
    separately so setup, solve and output cannot be confused for one another.
    """
    prof = profile if profile is not None else Profile(backend)
    with prof.phase('setup (mesh+input)'):
        S, mesh = build_solver_state(case_dir)
    with prof.phase('resolve solver'):
        solver = _resolve_solver(S['friclaw'], backend)
    with prof.phase('solve'):
        out = solver.run(S, nsteps=nsteps, verbose=verbose)

    nftnd = S['nftnd']
    fric_1idx = np.zeros((nftnd + 1, 101))
    fric_1idx[1:, 1:101] = out['fric']
    fnft_1idx = np.zeros(nftnd + 1)
    fnft_1idx[1:] = out['fnft']

    frt_path = os.path.join(case_dir, 'frt.txt0')
    with prof.phase('write frt'):
        frt_writer.write_frt(frt_path, mesh['meshCoor'], mesh['nsmp'],
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
    ap = argparse.ArgumentParser(prog='python3 -m eqdyna.standalone')
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
