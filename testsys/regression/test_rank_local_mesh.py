#! /usr/bin/env python3
"""
Regression guard (pathway item 64): a python-jax-mpi rank's RANK-LOCAL mesh is
the serial mesh's box, bit for bit, under the analytic global<->local map.

WHY THIS EXISTS. Before item 64 every rank built the serial mesh and
RESTRICTED it, so every index array was the gated serial array under an
injective relabelling and a partition bug could only raise. A rank now
GENERATES its box (eqdyna3d.build_solver_state(case_dir, part=...), the
meshgen.f90 way), and three failure classes open up that the frt gate cannot
reliably see -- the worst being a wrong equation number or boundary
classification far from the fault, whose effect never reaches a fault node in
a short gated run (design 9d189d7, section 1.2, class C). This test observes
the mesh DIRECTLY: it builds the serial mesh once and every rank's box of
every MPI4NodalQuant.DECOMP decomposition in one process, maps each local node
to its serial id analytically (grid offsets for regular nodes, the slave's
serial fault index for masters), and requires

  1. coordinates                       -- bitwise equal
  2. element rows (conn under the map, elemType, mat, eledet, eleshp, ss,
     phi, init_stress), in SERIAL ORDER -- bitwise equal
  3. dof counts and the fixed-boundary pattern equal, and the local equation
     numbering an order-preserving injection into the serial one
  4. fault rows (nsmp under the map, un/us/ud, fric_init) -- bitwise equal
  5. the partial nodal mass / fnms / arn, completed by an INDEPENDENT
     in-process simulation of MPI4NodalQuant's x/y/z relay and MPI4arn's
     DIVIDE/DUPLICATE rule, equal to the serial values: bitwise at nodes on
     no shared plane, to 1e-14 relative on shared planes (a different, and
     Fortran's, summation order)
  6. the exchange faces of every neighbour pair name the same serial nodes
     and equations in the same order (what MPI4NodalQuant.handshake checks at
     run time), and the owned fault rows partition the serial fault exactly
  7. meshgen's global fault/equation censuses equal the serial counts
  8. the zero-fault-node contract: how many ranks write a frt file and how
     many hold no fault node, per decomposition, equal FRT_SHAPE's data

Cases: test.tpv8 at its gated dx (planar, the opted-in MPI case), test.tpv10
on a shrunken, y-symmetric domain (insertFaultType>0: the rough y-blend must
use the MODEL's ymin/ymax, and a y boundary ON the fault plane is the
DUPLICATE arn path) and test.tpv36
coarsened to dx=1000 m (C_degen>3 wedges; a y boundary crossing the dipping
fault is the DIVIDE path).
Rank counts are the smallest set that reaches every path and both
zero-fault-node shapes (CI cost, ci_shard.py); 2 and 32 ranks were run once,
green, when this guard was written (FRT_SHAPE keeps their data). The test prints, per case, which arn paths it actually
exercised and fails if a path it is responsible for was never reached (a green
result that tested nothing is this repo's recurring failure).

No MPI, no Fortran, no time stepping. Pure python, one process.
"""
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))

from eqdyna import MPI4NodalQuant as MQ   # noqa: E402
from eqdyna import eqdyna3d, meshgen      # noqa: E402
from eqdyna import readInputFiles         # noqa: E402

# test.tpv23 (ntotft=2, item 145/PR #78 follow-up) at its native 250 m dx
# costs ~4 min just to build its serial mesh (measured: 236 s at test.tpv22's
# 200 m, same O(N) shape) -- far too slow for this file's "run constantly"
# tier. Coarsened to dx=1000 m (same value test.tpv36 already coarsens to,
# below): commensurate with BOTH fault boxes' extents (30000/1000=30),
# fault #2's own 1000 m y-stepover (1000/1000=1) and the domain extents
# (64000/1000=64, 20000/1000=20, 27000/1000=27) -- every commensurability
# checkInputConsistency enforces still holds, so this is the same benchmark
# geometry at a coarser sample, not a different one. `on_fault_vars` (and
# nfx/nfz/fx/fz) must be rebuilt AT the override dx -- case.setup's netCDF
# writer shapes it to the actual fault node count, which is dx-derived; a
# plain `par.dx = 1000.0` override (as test.tpv36's below) left the native-
# resolution array in place and case.setup raised a netCDF shape mismatch
# (21,31) vs (81,121) when first tried. The two per-fault builders invoked
# below are tpv22_23_common's OWN (not reinvented): this rebuilds the same
# physics field at a coarser sample, the same thing build_params does
# internally for its native dx, not a second formula.
_TPV23_COARSE_DX = 1000.0
_TPV23_COARSE_RAW = '''
import numpy as np
import tpv22_23_common as common
par.dx = %r
par.dy = %r
par.dz = %r
par.dt = 0.5 * par.dx / par.vp
par.nfx = round((common.FAULT1_XRANGE[1] - common.FAULT1_XRANGE[0]) / par.dx + 1)
par.nfz = round((common.FZMAX - common.FZMIN) / par.dz + 1)
par.fx = np.linspace(common.FAULT1_XRANGE[0], common.FAULT1_XRANGE[1], par.nfx)
par.fz = np.linspace(common.FZMIN, common.FZMAX, par.nfz)
par.on_fault_vars = common.build_on_fault_vars(par.fx, par.fz, par.nfx, par.nfz)
nfx2 = round((common.FAULT2_XRANGE[1] - common.FAULT2_XRANGE[0]) / par.dx + 1)
fx2 = np.linspace(common.FAULT2_XRANGE[0], common.FAULT2_XRANGE[1], nfx2)
fault2_vars = common.build_on_fault_vars(fx2, par.fz, nfx2, par.nfz)
par.onFaultVarsPerFault = [par.on_fault_vars, fault2_vars]
''' % (_TPV23_COARSE_DX, _TPV23_COARSE_DX, _TPV23_COARSE_DX)

# (label, case, par overrides, decompositions by rank count, arn paths it must reach)
CASES = (
    ('tpv8', 'test.tpv8', {}, (4, 8, 16), ('xsplit', 'zsplit')),
    # tpv10 keeps its dx (its rough-geometry file is sampled at it) and
    # shrinks its DOMAIN, symmetric in y so npy=2 puts a boundary ON y=0.
    ('tpv10-ysym', 'test.tpv10', {'xmin': -20.0e3, 'xmax': 20.0e3, 'ymin': -15.0e3,
                                  'ymax': 15.0e3, 'zmin': -20.0e3}, (4,), ('xsplit', 'ydup')),
    # ...and with its own asymmetric y (-20/+40 km), so the mey=1 boxes never
    # meet the fault: the rough-fault builders' fault-free-box path (found by
    # the full-scale tpv10 4-rank run, not by this guard's first version).
    ('tpv10-yasym', 'test.tpv10', {'xmin': -20.0e3, 'xmax': 20.0e3, 'zmin': -20.0e3},
     (4,), ('xsplit',)),
    ('tpv36', 'test.tpv36', {'dx': 1000.0}, (4, 8), ('xsplit', 'ydiv')),
    # test.tpv23 (ntotft=2): the only multi-fault case this file covers
    # (victor-reyes PR #78 audit, item 145 follow-up). 2 ranks picks
    # MPI4NodalQuant.DECOMP[2] = (2, 1, 1) -- npy=1, so a rank's box spans
    # BOTH faults' y (fault #1 at y=0, fault #2 at y=+1000 m) and, wherever
    # its x-slab meets the two faults' 10 km x-overlap band, owns fault rows
    # from BOTH faults in SCAN order (interleaved: fault_id goes 0,1,0,1,...
    # down each shared (x,z) column) -- the exact shape build_faces' Row 17
    # comment describes, and the configuration MPI4NodalQuant.DECOMP[4] =
    # (2, 2, 1) (this case's REGISTERED python-jax-mpi rank count) does NOT
    # reach: measured directly (mira, this follow-up), at 4 ranks the y-split
    # falls exactly on fault #1's own y=0 plane, so each rank owns only ONE
    # fault's rows and the interleave never occurs. 4 ranks is included too,
    # to mesh-check the actually-registered decomposition, but 2 is the one
    # that is vacuous without the interleave check below.
    ('tpv23-2fault', 'test.tpv23', {'_raw': _TPV23_COARSE_RAW}, (2, 4), ('xsplit',)),
)
# THE ZERO-FAULT-NODE CONTRACT, per decomposition, as MEASURED DATA (first run
# of this guard, 2026-09-24): (ranks that write a frt file, ranks whose box
# holds no fault node at all). The 3D split moves both numbers with the rank
# count -- tpv8 has fault-free boxes from 4 ranks up, tpv36 only at 8, and
# tpv10's y-symmetric domain at 4 ranks is the other shape: all four boxes
# COMPUTE the fault (it lies on their shared y plane) and only the two lower
# ones write it. A change to either number is a change to what the e2e cell's
# PY_MPI_EXPECTED_FRT_FILES counts, so it must be made here deliberately.
FRT_SHAPE = {
    ('tpv8', 2): (2, 0), ('tpv8', 4): (2, 2), ('tpv8', 8): (4, 4),
    ('tpv8', 16): (8, 8), ('tpv8', 32): (8, 24),
    ('tpv10-ysym', 2): (2, 0), ('tpv10-ysym', 4): (2, 0), ('tpv10-ysym', 8): (4, 0),
    ('tpv10-yasym', 4): (2, 2),
    ('tpv36', 2): (2, 0), ('tpv36', 4): (4, 0), ('tpv36', 8): (6, 2),
    ('tpv23-2fault', 2): (2, 0), ('tpv23-2fault', 4): (4, 0),
}
ELEM_KEYS = ('elemType', 'mat', 'eledet', 'eleshp', 'ss', 'phi', 'init_stress')
FAULT_KEYS = ('un', 'us', 'ud', 'fric_init')
REL = 1e-14


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'), env.get('PATH', '')])
    return env


def _make_case(case, case_dir, overrides):
    """overrides is `par.key = repr(value)` pairs, EXCEPT the reserved key
    `_raw`, whose value is appended VERBATIM (not repr'd) -- for overrides a
    scalar repr can't express, e.g. test.tpv22/23's per-fault `on_fault_vars`
    ndarray, which must be rebuilt at the override dx (case.setup's netCDF
    writer requires on_fault_vars shaped to the ACTUAL fault node count at
    whatever par.dx is when it runs, not the case's native resolution)."""
    env = _env()
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, case], env=env, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('create.newcase %s failed: %s' % (case, (r.stderr or '')[-800:]))
    p = os.path.join(case_dir, 'user_defined_params.py')
    text = open(p).read().rstrip('\n')
    raw = overrides.get('_raw', '')
    scalars = {k: v for k, v in overrides.items() if k != '_raw'}
    extra = ''.join('par.%s = %r\n' % kv for kv in scalars.items())
    open(p, 'w').write(text + '\n\n# test_rank_local_mesh.py: serial case inputs\n'
                       'par.nx = 1\npar.ny = 1\npar.nz = 1\n' + extra + raw)
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError('case.setup %s failed: %s' % (case, (r.stderr or r.stdout)[-800:]))


def _eq(name, a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.array_equal(a, b):
        nbad = (int(np.count_nonzero(a != b)) if a.shape == b.shape else 'shape %r vs %r'
                % (a.shape, b.shape))
        raise AssertionError('%s: rank-local differs from serial (%s)' % (name, nbad))


def _g_of_l(mesh_l, mesh_s):
    """1-based serial node id of every local node (index 0 unused)."""
    nx, ny, nz = mesh_l['n_local']
    nxg, nyg, nzg = mesh_l['n_global']
    ox, oy, oz = mesh_l['offsets']
    n_reg = nx * ny * nz
    s = np.arange(n_reg, dtype=np.int64)
    g_reg = ((ox + s // (nz * ny)) * nzg * nyg + (oz + (s % (nz * ny)) // ny) * nyg
             + (oy + s % ny) + 1)
    slaves_g = g_reg[mesh_l['nsmp'][:, 0] - 1]
    serial_slaves = mesh_s['nsmp'][:, 0]
    # mesh_s['nsmp'] is GROUPED BY FAULT (Row 17: stable-sorted by fault id),
    # not globally sorted by slave node id -- a bare searchsorted on it is
    # only valid at ntotft==1, where grouping is a no-op and scan-order slave
    # ids are already ascending. At ntotft>1 two faults' slave ids interleave
    # in absolute value (build_node_coordinates's docstring), so sort once
    # explicitly and map back through the permutation -- correct at any
    # ntotft, including 1 (sorting an already-sorted array is a no-op).
    order = np.argsort(serial_slaves, kind='stable')
    sorted_slaves = serial_slaves[order]
    pos = np.searchsorted(sorted_slaves, slaves_g)
    pos = np.clip(pos, 0, max(sorted_slaves.size - 1, 0))
    frow = order[pos] if sorted_slaves.size else pos
    if frow.size and not np.array_equal(serial_slaves[frow], slaves_g):
        raise AssertionError('a local fault node is not a serial fault node')
    g_master_by_row = mesh_s['nsmp'][frow, 1]
    # meshCoor's OWN master-node block is written in LOCAL SCAN order
    # (build_node_coordinates fills it straight from slave_ids_scan), but
    # `frow`/`g_master_by_row` above are in mesh_l['nsmp']'s ROW order --
    # fault-GROUPED (Row 17's stable sort), the same order FAULT_KEYS/`S[k]`
    # and `mesh_s['nsmp']` itself are read in elsewhere in this file. The two
    # orders coincide only at ntotft==1 (no grouping reorder happens). Each
    # row's OWN master id tells us which scan position it is -- a fault row
    # whose master id is n_reg+1+p was written at meshCoor's scan position p
    # by construction -- so scatter this row's mapped serial master id there
    # to build g_master in meshCoor's own order, without assuming the two
    # orders match.
    n_reg_l = n_reg
    scan_pos = mesh_l['nsmp'][:, 1] - n_reg_l - 1
    nftnd = mesh_l['nsmp'].shape[0]
    g_master = np.empty(nftnd, dtype=np.int64)
    if nftnd:
        g_master[scan_pos] = g_master_by_row
    return np.concatenate(([0], g_reg, g_master)), frow


def _serial_elem_rows(S_s, mesh_s, n_glob, mesh_l):
    """Serial element rows whose brick (top corner = max grid index over its
    nodes; masters count as their slave's grid point) lies in the box, in
    SERIAL order. Works for bricks and both wedge halves.

    `grid[n_reg:]` must be indexed by RAW master node id (i.e. SCAN position,
    since master ids are n_reg+1+scan_position by construction), but
    mesh_s['nsmp']'s rows are fault-GROUPED (Row 17's stable sort) -- the two
    orders coincide only at ntotft==1. Each row's own master id tells us its
    scan position (n_reg+1+p), so scatter each row's slave id there rather
    than assuming row order is scan order (same fix as _g_of_l's g_master)."""
    nxg, nyg, nzg = n_glob
    n_reg = nxg * nyg * nzg
    grid = np.arange(S_s['N'], dtype=np.int64)
    if mesh_s['nsmp'].shape[0]:
        scan_pos = mesh_s['nsmp'][:, 1] - n_reg - 1
        grid[n_reg + scan_pos] = mesh_s['nsmp'][:, 0] - 1
    cg = grid[S_s['conn']]
    gx, gz, gy = cg // (nzg * nyg), (cg % (nzg * nyg)) // nyg, cg % nyg
    top = (gx.max(axis=1), gy.max(axis=1), gz.max(axis=1))
    sel = np.ones(cg.shape[0], dtype=bool)
    for d in range(3):
        o, n = mesh_l['offsets'][d], mesh_l['n_local'][d]
        sel &= (top[d] >= o + 1) & (top[d] <= o + n - 1)
    return np.nonzero(sel)[0]


def _simulate_relay(parts, planes, key, arrs, num_dof_label):
    """Independent oracle for MPI4NodalQuant.relay across all ranks at once.
    Within one dimension a rank's minus and plus faces are disjoint node
    sets, so every exchange of that dimension sees the values as they stood
    when the dimension began; dimensions run x, y, z in sequence (the relay).
    """
    for d in range(3):
        snap = [a.copy() for a in arrs]
        for r, faces in enumerate(planes):
            for f in faces:
                if f['d'] != d:
                    continue
                partner = [g for g in planes[f['nb']] if g['d'] == d and g['ib'] == 1 - f['ib']]
                if len(partner) != 1:
                    raise AssertionError('%s: rank %d face %d has %d partners'
                                         % (num_dof_label, r, f['k'], len(partner)))
                arrs[r][f[key]] += snap[f['nb']][partner[0][key]]


def _simulate_arn(parts, meshes, arns, fault_boxes):
    """Independent oracle for meshgen.mpi4arn. Row 17 (multi-fault): the
    divide/duplicate decision is now PER ROW (keyed by that row's own
    fault's extent, `fault_boxes[fid]`), not once per face -- mirrors
    mpi4arn's own generalization. At ntotft==1 (every CASES entry here) the
    per-row mask is constant across a face, so this reduces to the old
    per-face branch exactly."""
    fext = np.array([((b[0], b[1]), (b[2], b[3]), (b[4], b[5])) for b in fault_boxes])
    used = set()
    flt_mpi = [[False] * 6 for _ in parts]
    for d in range(3):
        snap = [a.copy() for a in arns]
        for r, part in enumerate(parts):
            for ib in (0, 1):
                k = 2 * d + ib
                if part.dims[d] <= 1 or part.at_model_edge(d, ib):
                    continue
                idx = meshes[r]['flt_lists'][k]
                if idx.size == 0:
                    continue
                flt_mpi[r][k] = True
                nb = part.neighbour(d, ib)
                pidx = meshes[nb]['flt_lists'][2 * d + 1 - ib]
                if pidx.size != idx.size:
                    raise AssertionError('rank %d face %d: %d fault nodes, neighbour %d: %d'
                                         % (r, k, idx.size, nb, pidx.size))
                fid = meshes[r]['nsmp'][idx - 1, 2]
                divide = fext[fid, d, 1] != fext[fid, d, 0]
                if divide.any():
                    arns[r][idx[divide]] += snap[nb][pidx[divide]]
                    used.add(('xsplit', 'ydiv', 'zsplit')[d])
                if (~divide).any():
                    used.add(('xdup', 'ydup', 'zdup')[d])
    return flt_mpi, used


def _close(name, loc, ser, shared):
    loc, ser = np.asarray(loc), np.asarray(ser)
    _eq(name + ' (nodes on no shared plane)', loc[~shared], ser[~shared])
    scale = np.maximum(np.abs(ser[shared]), np.finfo(float).tiny)
    worst = float(np.max(np.abs(loc[shared] - ser[shared]) / scale)) if shared.any() else 0.0
    if worst > REL:
        raise AssertionError('%s: shared-plane relative diff %.3e > %.0e' % (name, worst, REL))
    return worst


def check_case(label, case, overrides, rank_counts, must_reach, tmp):
    case_dir = os.path.join(tmp, label)
    _make_case(case, case_dir, overrides)
    t0 = time.time()
    S_s, mesh_s = eqdyna3d.build_solver_state(case_dir)
    params, _g = readInputFiles.build_params(case_dir)
    xg, yg, zg, pmlb, bounds = meshgen.build_grid_lines(params)
    n_glob = (len(xg), len(yg), len(zg))
    fnms_s1 = np.concatenate(([0.0], S_s['fnms']))
    mass_s1 = np.concatenate(([0.0], S_s['nodalMassArr']))
    print('  %s serial: N=%d E=%d NEQ=%d nftnd=%d (%.1f s)'
          % (label, S_s['N'], S_s['E'], S_s['NEQ'], S_s['nftnd'], time.time() - t0))
    if tuple(meshgen.fault_census(xg, yg, zg, params))[0] != S_s['nftnd']:
        raise AssertionError('fault_census count != serial nftnd')
    if meshgen.equation_census(xg, yg, zg, params, pmlb, bounds) != S_s['NEQ']:
        raise AssertionError('equation_census != serial NEQ')
    if mesh_s['nsmp'].shape[0] and np.unique(mesh_s['nsmp'][:, 2]).size > 1:
        # Multi-fault vacuousness guard (item 145 follow-up, Fix 2): a case
        # whose GLOBAL master-node ids already sit at n_reg+1+row_position
        # for every fault-grouped row would pass build_faces' old buggy
        # n_reg+owned_rows arithmetic too, and this test would be exercising
        # nothing new over the ntotft==1 cases above. This checks a DATA
        # property (is this case's own serial mesh actually interleaved?),
        # not the fix's formula -- the correctness oracle stays the existing
        # bit-identical `nsmp` comparison against this same serial mesh,
        # below, per rank.
        n_reg_s = S_s['N'] - mesh_s['nsmp'].shape[0]
        expect_seq = n_reg_s + 1 + np.arange(mesh_s['nsmp'].shape[0])
        if np.array_equal(mesh_s['nsmp'][:, 1], expect_seq):
            raise AssertionError(
                '%s: serial mesh has %d fault ids but every master-node id '
                'already sits at n_reg+1+row_position -- this case is '
                'VACUOUS for the Row 17 master-id regression (build_faces\' '
                'old n_reg+owned_rows arithmetic would also have passed)'
                % (label, np.unique(mesh_s['nsmp'][:, 2]).size))
    reached = set()
    for nranks in rank_counts:
        parts = [MQ.Partition.for_size(r, nranks) for r in range(nranks)]
        built = [eqdyna3d.build_solver_state(case_dir, part=p) for p in parts]
        Ss = [b[0] for b in built]
        ms = [b[1] for b in built]
        maps = [_g_of_l(m, mesh_s) for m in ms]
        owned = []
        worst = 0.0
        for r, (S, m) in enumerate(built):
            g, frow = maps[r]
            if tuple(m['n_global']) != n_glob:
                raise AssertionError('n_global mismatch')
            _eq('meshCoor', m['meshCoor'][1:], mesh_s['meshCoor'][g[1:]])
            rows = _serial_elem_rows(S_s, mesh_s, n_glob, m)
            _eq('conn', g[S['conn'] + 1], S_s['conn'][rows] + 1)
            for k in ELEM_KEYS:
                _eq(k, S[k], np.asarray(S_s[k])[rows])
            _eq('ndof', S['ndof'], S_s['ndof'][g[1:] - 1])
            e_l, e_s = S['eq_ids'], S_s['eq_ids'][g[1:] - 1]
            _eq('fixed/sink pattern', e_l > 0, e_s > 0)
            pairs = np.stack((e_l[e_l > 0], e_s[e_l > 0]), axis=1)
            pairs = pairs[np.argsort(pairs[:, 0], kind='stable')]
            if not (np.all(np.diff(pairs[:, 1]) > 0)
                    and np.array_equal(pairs[:, 0], np.arange(1, S['NEQ'] + 1))):
                raise AssertionError('local equations are not an order-preserving '
                                     'injection into the serial numbering')
            # nsmp's 3rd column is fault_id (Row 17), a small int tag, NOT a
            # node id -- only columns 0/1 (slave/master) go through `g`. At
            # ntotft==1 fault_id is always 0 and g[0]==0 trivially agreed, so
            # this never mattered until a case with a real fault_id column.
            nsmp_g = np.empty_like(m['nsmp'])
            nsmp_g[:, :2] = g[m['nsmp'][:, :2]]
            nsmp_g[:, 2] = m['nsmp'][:, 2]
            _eq('nsmp', nsmp_g, mesh_s['nsmp'][frow])
            for k in FAULT_KEYS:
                _eq(k, S[k], np.asarray(S_s[k])[frow])
        # --- setup exchanges, simulated across all ranks
        arns = [m['arn1'].copy() for m in ms]
        flt_mpi, used = _simulate_arn(parts, ms, arns, ms[0]['fault_boxes'])
        reached |= used
        planes = [MQ.build_faces(parts[r], ms[r]['n_local'], Ss[r]['ndof'], Ss[r]['eq_ids'],
                                 ms[r]['flt_lists'], flt_mpi[r], ms[r]['nsmp'])
                 for r in range(nranks)]
        mass = [m['mass1'].copy() for m in ms]
        fnms = [m['fnms1'].copy() for m in ms]
        _simulate_relay(parts, planes, 'eqs', mass, 'mass')
        _simulate_relay(parts, planes, 'nodes', fnms, 'fnms')
        for r in range(nranks):
            g, frow = maps[r]
            S, m = built[r]
            for f in planes[r]:
                pf = [h for h in planes[f['nb']] if h['d'] == f['d'] and h['ib'] == 1 - f['ib']][0]
                _eq('face %d nodes vs neighbour %d' % (f['k'], f['nb']),
                    g[f['nodes']], maps[f['nb']][0][pf['nodes']])
                _eq('face %d eq counts vs neighbour' % f['k'], f['eqs_per_node'], pf['eqs_per_node'])
            shared_node = np.zeros(S['N'] + 1, dtype=bool)
            for f in planes[r]:
                shared_node[f['nodes']] = True
            worst = max(worst, _close('fnms', fnms[r][1:], fnms_s1[g[1:]],
                                      shared_node[1:]))
            e = S['eq_ids']
            live = (np.arange(e.shape[1])[None, :] < S['ndof'][:, None]) & (e > 0)
            shared_eq = np.zeros(S['NEQ'] + 1, dtype=bool)
            shared_eq[e[shared_node[1:]][live[shared_node[1:]]]] = True
            s_eq = np.zeros(S['NEQ'] + 1, dtype=np.int64)
            s_eq[e[live]] = S_s['eq_ids'][g[1:] - 1][live]
            worst = max(worst, _close('nodalMassArr', mass[r][1:],
                                      mass_s1[s_eq[1:]], shared_eq[1:]))
            if frow.size:
                on_face = np.zeros(frow.size, dtype=bool)
                for k in range(6):
                    on_face[m['flt_lists'][k] - 1] = True
                worst = max(worst, _close('arn', arns[r][1:], np.asarray(S_s['arn'])[frow],
                                          on_face))
            own = MQ.owned_mask(parts[r], m['nsmp'][:, 0], m['n_local']) if frow.size else \
                np.zeros(0, dtype=bool)
            owned.append(frow[own])
            if tuple(m['fault_census']) != tuple(meshgen.fault_census(xg, yg, zg, params)):
                raise AssertionError('rank %d fault census differs' % r)
        allown = np.sort(np.concatenate(owned)) if owned else np.zeros(0)
        _eq('owned fault rows partition the serial fault', allown, np.arange(S_s['nftnd']))
        files = sum(1 for o in owned if o.size)
        empty = sum(1 for S in Ss if int(S['nftnd']) == 0)
        if (files, empty) != FRT_SHAPE[(label, nranks)]:
            raise AssertionError(
                '%s at %d ranks: %d rank(s) write frt and %d hold no fault node; '
                'FRT_SHAPE records %r' % (label, nranks, files, empty,
                                          FRT_SHAPE[(label, nranks)]))
        print('    %2d ranks %r: %d ranks write frt, fault computed per rank %s, '
              'shared-plane max rel diff %.2e, arn paths %s'
              % (nranks, parts[0].dims, files, [int(s['nftnd']) for s in Ss], worst,
                 sorted(used) or '-'))
    missing = [p for p in must_reach if p not in reached]
    if missing:
        raise AssertionError('%s: arn path(s) %s never exercised by %r -- this test would '
                             'be green without testing them' % (label, missing, rank_counts))


def main():
    print('Regression guard: rank-local box mesh == serial mesh under the analytic map '
          '(pathway item 64)')
    fails = []
    with tempfile.TemporaryDirectory(prefix='ranklocal.') as tmp:
        for label, case, ov, counts, must in CASES:
            try:
                check_case(label, case, ov, counts, must, tmp)
            except Exception as exc:          # noqa: BLE001 -- reported, then FAIL
                import traceback
                traceback.print_exc()
                fails.append('%s: %s' % (label, exc))
    if fails:
        print('FAIL test_rank_local_mesh')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_rank_local_mesh (%d cases)' % len(CASES))
    return 0


if __name__ == '__main__':
    sys.exit(main())
