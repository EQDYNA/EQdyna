"""Post-meshing screen (board row 167: wedge-mesh mis-wiring).

Counterpart: src/fortran/screenMesh.f90 (screenMesh): same two checks, same
counts. Called after the fault-area accumulation and BEFORE MPI4arn. Both
backends STOP the run (hard error, no switch, no fallback) when a check fails;
nothing here repairs a mesh. Each check has its own exit code
(src/fortran/errorCodes.f90); all are evaluated, the first failing check sets
the exit status.

  sides     a (cell,node) attachment to the MASTER node of a pair has the cell
            centre on the +un side of the pair, one to its SLAVE node on the -un
            side (signed cosine strictly > 0). Cells whose centre is beyond the
            fault patch in x or z are not tested.
  bridging  no cell holds both nodes of one pair; the two nodes of a pair
            coincide to the mesh tolerance tol.
"""
import numpy as np

ERR_MESH_SCREEN_SIDES = 81
ERR_MESH_SCREEN_BRIDGING = 82
SCREEN_CODES = {'sides': ERR_MESH_SCREEN_SIDES, 'bridging': ERR_MESH_SCREEN_BRIDGING}


class MeshScreenError(RuntimeError):
    """Raised by screen_mesh. .failures maps each failed check name to its
    message; .code is the exit status of the FIRST failed check in
    SCREEN_CODES order (the number screenMesh.f90 aborts with)."""

    def __init__(self, failures, measured=None):
        self.failures = dict(failures)
        self.measured = dict(measured or {})
        first = next(k for k in SCREEN_CODES if k in self.failures)
        self.code = SCREEN_CODES[first]
        super().__init__('mesh screening FAILED: ' + ' | '.join(
            '[%s] %s' % (k, self.failures[k]) for k in SCREEN_CODES if k in self.failures))


def _e(v):
    return '%.3E' % v


def _screen_measure(meshCoor, conn, elem_type, nsmp, un, tol):
    """Evaluate both checks. Returns (m, failures): m the measured scalars,
    failures {check: msg}. meshCoor (N+1,3) and un (nft+1,3) keep the 1-indexed
    row-0 convention of build_node_coordinates/build_fault_geometry; conn is
    (E,8) 1-indexed node ids; elem_type (E,); nsmp (nft,3) is [slave, master,
    0-indexed fault]."""
    X = np.asarray(meshCoor, dtype=float)
    conn = np.asarray(conn, dtype=np.int64)
    nsmp = np.asarray(nsmp, dtype=np.int64)
    nft = nsmp.shape[0]
    fail = {}
    m = {'elements': int(conn.shape[0]), 'nodes': int(X.shape[0] - 1), 'fault_pairs': int(nft)}

    role = np.zeros(X.shape[0], dtype=np.int8)     # 1 slave, 2 master
    pid = np.full(X.shape[0], -1, dtype=np.int64)  # 0-based pair row
    role[nsmp[:, 0]] = 1
    role[nsmp[:, 1]] = 2
    pid[nsmp[:, 0]] = np.arange(nft)
    pid[nsmp[:, 1]] = np.arange(nft)
    flt = nsmp[:, 2]
    # patch extent of each fault in strike (x) and depth (z), from the master nodes
    xm, zm = X[nsmp[:, 1], 0], X[nsmp[:, 1], 2]
    x_lo, x_hi, z_lo, z_hi = (np.empty(nft) for _ in range(4))
    for ift in np.unique(flt):
        s_ = flt == ift
        x_lo[s_], x_hi[s_], z_lo[s_], z_hi[s_] = xm[s_].min(), xm[s_].max(), zm[s_].min(), zm[s_].max()

    # distinct (cell, fault node) attachments, element order
    ce, sl = np.nonzero(role[conn] != 0)
    ce, nodes = np.unique(np.stack([ce, conn[ce, sl]], axis=1), axis=0).T
    pp, rl = pid[nodes], role[nodes]

    # ---- bridging: a cell holding both nodes of one pair
    key = ce * max(nft, 1) + pp
    br = np.unique(np.intersect1d(key[rl == 1], key[rl == 2]) // max(nft, 1))
    m['bridging_cells'] = int(br.size)
    bad = []
    if br.size:
        bad.append('%d cell(s) hold both nodes of a split pair, first element %d (type %d)' % (
            br.size, br[0] + 1, elem_type[br[0]]))
    sep = np.abs(X[nsmp[:, 0]] - X[nsmp[:, 1]]).max(axis=1) if nft else np.zeros(1)
    m['coincidence_max'] = float(sep.max())
    if not sep.max() <= tol:
        w = int(np.argmax(sep))
        bad.append('nodes %d/%d of a pair are %s m apart (tol %s)' % (nsmp[w, 0], nsmp[w, 1], _e(sep[w]), _e(tol)))
    if bad:
        fail['bridging'] = ' ; '.join(bad)

    # ---- sides: cell centre on the own side of the fault plane
    cen = X[conn[ce]].sum(axis=1) / 8.0
    inside = ((cen[:, 0] >= x_lo[pp]) & (cen[:, 0] <= x_hi[pp]) &
              (cen[:, 2] >= z_lo[pp]) & (cen[:, 2] <= z_hi[pp]))
    ce, nodes, pp, rl, cen = ce[inside], nodes[inside], pp[inside], rl[inside], cen[inside]
    dv = cen - X[nodes]
    nv = np.sqrt(np.einsum('ij,ij->i', dv, dv))
    cosd = np.einsum('ij,ij->i', dv, un[1:nft + 1][pp]) / np.where(nv > 0, nv, 1.0)
    margin = np.where(rl == 2, cosd, -cosd)
    wrong = np.nonzero(~(margin > 0.0))[0]
    mw, sw = int((rl[wrong] == 2).sum()), int((rl[wrong] == 1).sum())
    m.update(sides_attachments=int(margin.size),
             sides_min_margin=float(margin.min()) if margin.size else float('nan'),
             sides_wrong_master=mw, sides_wrong_slave=sw)
    if wrong.size:
        fail['sides'] = ('%d of %d (cell,node) attachments on the wrong side of the fault plane '
                         '(%d master-attached on the footwall side, %d slave-attached on the hanging-wall side), '
                         'first element %d (type %d)' % (wrong.size, margin.size, mw, sw,
                                                         ce[wrong[0]] + 1, elem_type[ce[wrong[0]]]))
    return m, fail


def screen_mesh(meshCoor, conn, elem_type, nsmp, un, tol):
    """Post-meshing screen; raises MeshScreenError if any check fails, returns
    the measured scalars otherwise (see the module docstring)."""
    m, fail = _screen_measure(meshCoor, conn, elem_type, nsmp, un, tol)
    if fail:
        raise MeshScreenError(fail, m)
    return m
