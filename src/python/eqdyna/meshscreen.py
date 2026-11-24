"""Post-meshing screening (owner decision 2026-11-23, board row "wedge-mesh
mis-wiring", deliverable 1).

Counterpart: src/fortran/screenMesh.f90 (screenMesh): same algorithm, same
thresholds, same message text up to float formatting. It is called from
meshgen.f90 after the fault-area accumulation and BEFORE MPI4arn, so arn is
this rank's partial sum exactly as meshgen.build_fault_geometry returns it.
Both backends STOP the run (hard error, no switch, no fallback) when a check
fails; nothing here repairs a mesh.

Five checks, each its own named failure with its own exit code (81-85,
src/fortran/errorCodes.f90). All five are evaluated, so one message lists
every check that failed; the first failing check's code is the exit status.

  sides        a (cell,node) attachment to a pair's MASTER node has the cell
               centre on the +un (hanging-wall) side of the pair, one to its
               SLAVE node on the -un side (signed cosine strictly > 0). Cells
               whose centre is beyond the fault patch in x or z touch the
               fault only through an edge node and are not tested.
  bridging     no cell holds both nodes of one pair; the two nodes of a pair
               coincide to the mesh tolerance tol.
  fault_area   each pair's arn equals a quarter of the vector areas of the
               meshed facets (element faces whose four nodes are all master
               ids; slave ids where no master-side cell exists) touching it,
               to SCREEN_AREA_REL; every pair has 1, 2 or 4 such facets.
  volume_mass  every cell's trilinear Jacobian is > SCREEN_DET_MIN of its mean
               at all 8 Gauss points, cell volumes sum to the rank box volume
               taken from the node layers, each interior fault node's volume
               share (the mass proxy: rho is piecewise-constant) is positive
               on both sides with master/slave ratio <= SCREEN_VOL_RATIO_MAX.
  normals      un, us, ud are an orthonormal triple and un is parallel to the
               area-weighted normal of the meshed facets at the node.

The limits are repeated in screenMesh.f90; test_mesh_screening.py fails if
the two files differ.
"""
import numpy as np

ERR_MESH_SCREEN_SIDES = 81
ERR_MESH_SCREEN_BRIDGING = 82
ERR_MESH_SCREEN_FAULT_AREA = 83
ERR_MESH_SCREEN_VOLUME_MASS = 84
ERR_MESH_SCREEN_NORMALS = 85
SCREEN_CODES = {'sides': ERR_MESH_SCREEN_SIDES, 'bridging': ERR_MESH_SCREEN_BRIDGING,
                'fault_area': ERR_MESH_SCREEN_FAULT_AREA,
                'volume_mass': ERR_MESH_SCREEN_VOLUME_MASS,
                'normals': ERR_MESH_SCREEN_NORMALS}

SCREEN_SIDE_MIN = 0.0            # signed cosine must be strictly greater
SCREEN_AREA_REL = 1.0e-12        # facet identity, relative
SCREEN_VOL_REL = 1.0e-10         # rank-box volume, relative
SCREEN_DET_MIN = 0.1             # min detJ at a Gauss point / (cell volume / 8)
SCREEN_VOL_RATIO_MAX = 2.0       # interior fault node master/slave volume share
SCREEN_NORMAL_TOL = 1.0e-9       # orthonormality of (un,us,ud), and planar facet sin
SCREEN_NORMAL_TOL_ROUGH = 0.2    # facet sin on a non-planar (insertFaultType>0) fault

_HEX_FACES = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))
_HEX_A = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                   [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)


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


def _hex_gauss_detj(xc):
    """det J of the trilinear brick map at the 2x2x2 Gauss points, (E, 8)."""
    g = 1.0 / np.sqrt(3.0)
    out = np.empty((xc.shape[0], 8))
    k = 0
    for sz in (-g, g):
        for sy in (-g, g):
            for sx in (-g, g):
                f = 1.0 + _HEX_A * np.array([sx, sy, sz])
                dn = np.empty((8, 3))
                dn[:, 0] = _HEX_A[:, 0] * f[:, 1] * f[:, 2] / 8.0
                dn[:, 1] = _HEX_A[:, 1] * f[:, 0] * f[:, 2] / 8.0
                dn[:, 2] = _HEX_A[:, 2] * f[:, 0] * f[:, 1] / 8.0
                jac = np.einsum('ij,eik->ejk', dn, xc)
                out[:, k] = np.linalg.det(jac)
                k += 1
    return out


def _e(v):
    return '%.3E' % v


def _screen_measure(meshCoor, conn, elem_type, nsmp, un, us, ud, arn, lines, tol, planar):
    """Evaluate the five checks. Returns (m, failures): m the measured scalars,
    failures {check: msg}. meshCoor (N+1,3) and un/us/ud/arn (nft+1,...) keep
    the 1-indexed row-0 convention of build_node_coordinates/
    build_fault_geometry; conn is (E,8) 1-indexed node ids; elem_type (E,);
    nsmp (nft,3) is [slave, master, 0-indexed fault]; lines is the rank's
    (xline, yline, zline); planar[ift] says fault ift has insertFaultType==0."""
    xline, yline, zline = lines
    nx, ny, nz = len(xline), len(yline), len(zline)
    X = np.asarray(meshCoor, dtype=float)
    conn = np.asarray(conn, dtype=np.int64)
    nsmp = np.asarray(nsmp, dtype=np.int64)
    nft = nsmp.shape[0]
    nnode = X.shape[0]
    fail = {}
    m = {'elements': int(conn.shape[0]), 'nodes': int(nnode - 1), 'fault_pairs': int(nft)}

    def bad(name, msg):
        fail.setdefault(name, []).append(msg)

    role = np.zeros(nnode, dtype=np.int8)     # 1 slave, 2 master
    pid = np.full(nnode, -1, dtype=np.int64)  # 0-based pair row
    role[nsmp[:, 0]] = 1
    role[nsmp[:, 1]] = 2
    pid[nsmp[:, 0]] = np.arange(nft)
    pid[nsmp[:, 1]] = np.arange(nft)
    flt = nsmp[:, 2]
    unp, usp, udp, arnp = un[1:nft + 1], us[1:nft + 1], ud[1:nft + 1], arn[1:nft + 1]

    # patch extent of each pair's fault in strike (x) and depth (z)
    xm_, zm_ = X[nsmp[:, 1], 0], X[nsmp[:, 1], 2]
    x_lo, x_hi, z_lo, z_hi = (np.empty(nft) for _ in range(4))
    for ift in np.unique(flt):
        s_ = flt == ift
        x_lo[s_], x_hi[s_] = xm_[s_].min(), xm_[s_].max()
        z_lo[s_], z_hi[s_] = zm_[s_].min(), zm_[s_].max()
    on_edge = ((np.abs(xm_ - x_lo) <= tol) | (np.abs(xm_ - x_hi) <= tol) |
               (np.abs(zm_ - z_lo) <= tol) | (np.abs(zm_ - z_hi) <= tol))

    rc = role[conn]
    touch = np.nonzero((rc != 0).any(axis=1))[0]
    ct = conn[touch]
    cen = X[ct].sum(axis=1) / 8.0

    # ---- sides ------------------------------------------------------------
    ci, cj = np.nonzero(rc[touch] != 0)
    nodes = ct[ci, cj]
    uq = np.unique(np.stack([ci, nodes], axis=1), axis=0)
    ci, nodes = uq[:, 0], uq[:, 1]
    pp = pid[nodes]
    dv = cen[ci] - X[nodes]
    dn = np.einsum('ij,ij->i', dv, unp[pp])
    nv = np.sqrt(np.einsum('ij,ij->i', dv, dv))
    cosd = dn / np.where(nv > 0, nv, 1.0)
    sgn = np.where(role[nodes] == 2, 1.0, -1.0)
    margin = sgn * cosd
    inside = ((cen[ci, 0] >= x_lo[pp]) & (cen[ci, 0] <= x_hi[pp]) &
              (cen[ci, 2] >= z_lo[pp]) & (cen[ci, 2] <= z_hi[pp]))
    ci, nodes, pp, margin = ci[inside], nodes[inside], pp[inside], margin[inside]
    m['sides_attachments'] = int(margin.size)
    m['sides_min_margin'] = float(margin.min()) if margin.size else float('nan')
    wrong = np.nonzero(~(margin > SCREEN_SIDE_MIN))[0]
    m['sides_wrong_master'] = 0
    m['sides_wrong_slave'] = 0
    if wrong.size:
        mw = int((role[nodes[wrong]] == 2).sum())
        sw = int((role[nodes[wrong]] == 1).sum())
        m['sides_wrong_master'] = mw
        m['sides_wrong_slave'] = sw
        # worst = smallest margin, ties to the lowest element then node id
        order = np.lexsort((nodes[wrong], touch[ci[wrong]], margin[wrong]))
        w = wrong[order[0]]
        bad('sides', '%d of %d (cell,node) attachments on the wrong side of the fault plane '
            '(%d master-attached on the footwall side, %d slave-attached on the hanging-wall side); '
            'worst: element %d (type %d) holds %s node %d of pair %d (fault %d), margin %s' % (
                wrong.size, margin.size, mw, sw, touch[ci[w]] + 1, elem_type[touch[ci[w]]],
                'master' if role[nodes[w]] == 2 else 'slave', nodes[w], pp[w] + 1,
                flt[pp[w]] + 1, _e(margin[w])))

    # ---- bridging / coincidence -------------------------------------------
    sp = np.where(rc[touch] == 1, pid[ct], -1)
    mp = np.where(rc[touch] == 2, pid[ct], -1)
    both = ((sp[:, :, None] == mp[:, None, :]) & (sp[:, :, None] >= 0)).any(axis=(1, 2))
    m['bridging_cells'] = int(both.sum())
    if both.any():
        f0 = touch[np.nonzero(both)[0][0]]
        bad('bridging', '%d cell(s) hold both nodes of a split pair, first element %d (type %d)' % (
            int(both.sum()), f0 + 1, elem_type[f0]))
    sep = np.abs(X[nsmp[:, 0]] - X[nsmp[:, 1]]).max(axis=1) if nft else np.zeros(0)
    m['coincidence_max'] = float(sep.max()) if nft else 0.0
    if nft and not (sep.max() <= tol):
        w = int(np.argmax(sep))
        bad('bridging', 'pair %d (fault %d) nodes %d/%d are %s m apart (tol %s)' % (
            w + 1, flt[w] + 1, nsmp[w, 0], nsmp[w, 1], _e(sep[w]), _e(tol)))

    # ---- facets: element faces whose four nodes are all master (else slave) ids
    fnodes = ct[:, list(_HEX_FACES)].reshape(-1, 4)
    fr = role[fnodes]
    distinct = np.ones(fnodes.shape[0], dtype=bool)
    for a in range(4):
        for b in range(a + 1, 4):
            distinct &= fnodes[:, a] != fnodes[:, b]
    fcnt = np.zeros(nft, dtype=np.int64)
    fexp = np.zeros(nft)
    nacc = np.zeros((nft, 3))
    for r in (2, 1):
        sel = np.nonzero(distinct & (fr == r).all(axis=1))[0]
        fn = fnodes[sel]
        va = 0.5 * np.cross(X[fn[:, 2]] - X[fn[:, 0]], X[fn[:, 3]] - X[fn[:, 1]])
        ar = np.sqrt(np.einsum('ij,ij->i', va, va))
        c_r = np.zeros(nft, dtype=np.int64)
        e_r = np.zeros(nft)
        n_r = np.zeros((nft, 3))
        for k in range(4):
            p = pid[fn[:, k]]
            np.add.at(c_r, p, 1)
            np.add.at(e_r, p, 0.25 * ar)
            sg = np.sign(np.einsum('ij,ij->i', va, unp[p]))
            sg[sg == 0] = 1.0
            np.add.at(n_r, p, va * sg[:, None])
        if r == 2:
            fcnt, fexp, nacc = c_r, e_r, n_r
        else:
            use = fcnt == 0
            fcnt = np.where(use, c_r, fcnt)
            fexp = np.where(use, e_r, fexp)
            nacc = np.where(use[:, None], n_r, nacc)

    # ---- fault_area -------------------------------------------------------
    nofacet = np.nonzero(fcnt == 0)[0]
    m['area_nofacet_pairs'] = int(nofacet.size)
    if nofacet.size:
        bad('fault_area', '%d pair(s) have no meshed facet on either side, first pair %d (fault %d)' % (
            nofacet.size, nofacet[0] + 1, flt[nofacet[0]] + 1))
    badcnt = np.nonzero((fcnt > 0) & (fcnt != 1) & (fcnt != 2) & (fcnt != 4))[0]
    m['area_badcount_pairs'] = int(badcnt.size)
    if badcnt.size:
        bad('fault_area', '%d pair(s) carry a facet count other than 1, 2 or 4, first pair %d has %d' % (
            badcnt.size, badcnt[0] + 1, fcnt[badcnt[0]]))
    rel = np.where(fcnt > 0, np.abs(arnp - fexp) / np.where(fexp > 0, fexp, 1.0), 0.0)
    m['area_facet_rel_max'] = float(rel.max()) if nft else 0.0
    for ift in np.unique(flt):
        sel = np.nonzero(flt == ift)[0]
        w = sel[np.argmax(rel[sel])]
        if not rel[w] <= SCREEN_AREA_REL:
            bad('fault_area', 'fault %d: node area of pair %d is %s but the incident facets give %s '
                '(rel %s > %s); %d of %d pairs off' % (
                    ift + 1, w + 1, _e(arnp[w]), _e(fexp[w]), _e(rel[w]), _e(SCREEN_AREA_REL),
                    int((rel[sel] > SCREEN_AREA_REL).sum()), sel.size))

    # ---- normals ----------------------------------------------------------
    def dot(a, b):
        return np.einsum('ij,ij->i', a, b)
    dev = np.maximum.reduce([np.abs(np.sqrt(dot(v, v)) - 1.0) for v in (unp, usp, udp)] +
                            [np.abs(dot(unp, usp)), np.abs(dot(unp, udp)), np.abs(dot(usp, udp))])
    nn = np.sqrt(dot(nacc, nacc))
    nf = nacc / np.where(nn > 0, nn, 1.0)[:, None]
    cr = np.cross(unp, nf)
    ang = np.where(fcnt > 0, np.sqrt(dot(cr, cr)), 0.0)
    m['normal_orthonormal_dev_max'] = float(dev.max()) if nft else 0.0
    m['normal_facet_sin_max'] = float(ang.max()) if nft else 0.0
    for ift in np.unique(flt):
        sel = np.nonzero(flt == ift)[0]
        lim = SCREEN_NORMAL_TOL if planar[ift] else SCREEN_NORMAL_TOL_ROUGH
        w = sel[np.argmax(dev[sel])]
        if not dev[w] <= SCREEN_NORMAL_TOL:
            bad('normals', 'fault %d: (un,us,ud) of pair %d is not orthonormal (worst deviation %s > %s)' % (
                ift + 1, w + 1, _e(dev[w]), _e(SCREEN_NORMAL_TOL)))
        w = sel[np.argmax(ang[sel])]
        if not ang[w] <= lim:
            bad('normals', 'fault %d: un of pair %d is not perpendicular to its meshed facets '
                '(sin angle %s > %s); %d of %d pairs off' % (
                    ift + 1, w + 1, _e(ang[w]), _e(lim), int((ang[sel] > lim).sum()), sel.size))

    # ---- volume_mass ------------------------------------------------------
    nel = conn.shape[0]
    vol = np.empty(nel)
    detmin = np.empty(nel)
    for s in range(0, nel, 200000):
        dj = _hex_gauss_detj(X[conn[s:s + 200000]])
        vol[s:s + 200000] = dj.sum(axis=1)
        detmin[s:s + 200000] = dj.min(axis=1)
    rel_det = detmin / np.where(vol > 0, vol / 8.0, 1.0)
    m['det_min_rel'] = float(rel_det.min()) if nel else float('nan')
    badj = np.nonzero(~(rel_det > SCREEN_DET_MIN))[0]
    m['det_bad_cells'] = int(badj.size)
    if badj.size:
        bad('volume_mass', '%d of %d cells have a Gauss-point Jacobian below %s of their mean, first element %d '
            '(type %d, volume %s)' % (badj.size, nel, _e(SCREEN_DET_MIN), badj[0] + 1,
                                      elem_type[badj[0]], _e(vol[badj[0]])))
    gid = (np.arange(nx)[:, None, None] * nz + np.arange(nz)[None, :, None]) * ny + \
        np.arange(ny)[None, None, :] + 1
    h = X[gid[:, :, ny - 1], 1] - X[gid[:, :, 0], 1]
    hc = 0.25 * (h[:-1, :-1] + h[1:, :-1] + h[:-1, 1:] + h[1:, 1:])
    vexp = float(np.sum(np.diff(np.asarray(xline, dtype=float))[:, None] *
                        np.diff(np.asarray(zline, dtype=float))[None, :] * hc))
    vsum = float(vol.sum())
    m['volume_rel'] = abs(vsum - vexp) / vexp
    if not m['volume_rel'] <= SCREEN_VOL_REL:
        bad('volume_mass', 'cell volumes sum to %s m^3 but the rank box holds %s m^3 (rel %s > %s)' % (
            _e(vsum), _e(vexp), _e(m['volume_rel']), _e(SCREEN_VOL_REL)))
    share = np.zeros(nnode)
    np.add.at(share, conn.ravel(), np.repeat(vol / 8.0, 8))
    vs, vm = share[nsmp[:, 0]], share[nsmp[:, 1]]
    iy_s = (nsmp[:, 0] - 1) % ny
    edge_y = (iy_s == 0) | (iy_s == ny - 1)
    empty = np.nonzero(((vs <= 0) | (vm <= 0)) & ~edge_y)[0]
    m['share_empty_pairs'] = int(empty.size)
    if empty.size:
        bad('volume_mass', '%d fault node(s) have no volume share on a side, first pair %d (fault %d)' % (
            empty.size, empty[0] + 1, flt[empty[0]] + 1))
    both_ok = (vs > 0) & (vm > 0) & ~on_edge
    ratio = np.where(both_ok, vm / np.where(vs > 0, vs, 1.0), 1.0)
    lr = np.abs(np.log(ratio))
    m['mass_ratio_max'] = float(np.exp(lr.max())) if nft else 1.0
    nout = int((lr > np.log(SCREEN_VOL_RATIO_MAX)).sum())
    m['mass_ratio_out_pairs'] = nout
    if nout:
        w = int(np.argmax(lr))
        bad('volume_mass', 'fault node pair %d (fault %d) has master/slave volume share %s '
            '(limit %s); %d pairs outside' % (w + 1, flt[w] + 1, _e(ratio[w]), _e(SCREEN_VOL_RATIO_MAX), nout))
    return m, {k: ' ; '.join(v) for k, v in fail.items()}


def screen_mesh(meshCoor, conn, elem_type, nsmp, un, us, ud, arn, lines, tol, planar):
    """Post-meshing screen; raises MeshScreenError if any check fails, returns
    the measured scalars otherwise (see the module docstring)."""
    m, fail = _screen_measure(meshCoor, conn, elem_type, nsmp, un, us, ud, arn, lines, tol, planar)
    if fail:
        raise MeshScreenError(fail, m)
    return m
