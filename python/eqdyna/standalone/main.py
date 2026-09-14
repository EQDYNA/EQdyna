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
import importlib
import os
import sys

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
    """Returns the run()-providing module for `friclaw` under `backend`
    ('jax' or 'numpy'). backend='jax' (the default) falls back to numpy
    with an explicit printed notice -- never silently -- if jaxlib isn't
    importable in this environment."""
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
        print('NOTICE: --backend jax requested (the default) but %s is not importable (%s) -- '
              'falling back to the NumPy solver for friclaw=%d.'
              % (_JAX_MODULE_BY_FRICLAW[friclaw], e, friclaw), file=sys.stderr)
        return _NUMPY_SOLVER_BY_FRICLAW[friclaw]


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
    conn, elem_type, mat = meshgen.build_elements(xline, yline, zline, params, pmlb, nsmp, material, meshCoor)
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
    )
    mesh = dict(meshCoor=meshCoor, nsmp=nsmp)
    return S, mesh


def run_case(case_dir, nsteps=None, verbose=True, backend=DEFAULT_BACKEND):
    """Builds S (zero pydump reads), dispatches to the friclaw-appropriate
    solver's run() under the requested `backend` ('jax', the default, or
    'numpy') -- port.py/port_jax.py friclaw==1, port_rsf.py/port_rsf_jax.py
    friclaw==4, port_tp.py/port_tp_jax.py friclaw==5 -- writes frt.txt0 via
    frt_writer.write_frt (byte-exact Fortran E18.7E4 format). Returns the
    path written."""
    S, mesh = build_solver_state(case_dir)
    solver = _resolve_solver(S['friclaw'], backend)
    out = solver.run(S, nsteps=nsteps, verbose=verbose)

    nftnd = S['nftnd']
    fric_1idx = np.zeros((nftnd + 1, 101))
    fric_1idx[1:, 1:101] = out['fric']
    fnft_1idx = np.zeros(nftnd + 1)
    fnft_1idx[1:] = out['fnft']

    frt_path = os.path.join(case_dir, 'frt.txt0')
    frt_writer.write_frt(frt_path, mesh['meshCoor'], mesh['nsmp'], fnft_1idx, fric_1idx)
    return frt_path


def main():
    ap = argparse.ArgumentParser(prog='python3 -m eqdyna.standalone')
    ap.add_argument('case_dir')
    ap.add_argument('nsteps', nargs='?', type=int, default=None)
    ap.add_argument('--backend', choices=('jax', 'numpy'), default=DEFAULT_BACKEND,
                     help="solver backend (default: %(default)s; falls back to numpy "
                          "with a printed notice if jaxlib isn't importable)")
    args = ap.parse_args()
    path = run_case(args.case_dir, nsteps=args.nsteps, backend=args.backend)
    print('wrote', path)


if __name__ == '__main__':
    main()
