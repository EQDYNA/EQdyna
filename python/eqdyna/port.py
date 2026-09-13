"""
Vectorized NumPy port of EQdyna's explicit time loop for tpv8 (friclaw=1,
serial/np=1). NOT a from-scratch mesh generator: mesh/mass/shape-function
state is loaded from Fortran `pydump_*` files written by `src/pydump.f90`
(instrumentation added for this spike, called once before the time loop
starts). This ports only the per-step kernels: velDispUpdate, assembleGlobalKU
(interior calcElemKU + PML calcPMLElemKU), hrglss (C_hg==1 KF78), faulting
(getNsdSlipSliprateTraction + solveSWTW friclaw==1 + swtwNucleation +
storeRuptureTime). calcElemMass's contribution is analytically zero for this
case (rdampm=0.0d0, C_elastic=1) and is not computed at all -- see README.

UPDATE (post a844e92, master): `region_damp` now ports the FIXED PML
region cascade (src/func_lib.f90's `pmlRegionDistance`, post-refactor) --
region 14 correctly tests y against ymax0, not xmax0 (pathway_forward.md
item 8's bug, since fixed). Also reproduces `pmlRegionDistance`'s
`bound_inclusive` flag exactly: computePMLDampingVector.f90 (node-based,
called from velDispUpdate) always uses inclusive (>=/<=) comparisons;
assembleGlobalKU.f90's calcPMLElemKU (element-center-based) always uses
strict (>/<) -- a genuine, pre-existing difference between the two call
sites, not something this refactor (or this port) unifies. Literal
constant forms (3*vmaxPML/2/delta*log(1/R)*(.)**2) are preserved.
"""
import numpy as np
import os

def _read_scalars(path):
    vals = []
    with open(path) as f:
        for line in f:
            vals.extend(line.split())
    return vals

def load(case_dir):
    h = _read_scalars(os.path.join(case_dir, 'pydump_header.txt'))
    it = iter(h)
    N = int(next(it)); E = int(next(it)); NEQ = int(next(it))
    nen = int(next(it)); ned = int(next(it))
    nftnd = int(next(it)); ntotft = int(next(it)); nstep = int(next(it))
    dt = float(next(it)); w = float(next(it)); rdampk = float(next(it))
    rdampm = float(next(it)); kapa_hg = float(next(it)); R = float(next(it))
    nPML = int(next(it)); vmaxPML = float(next(it))
    PMLb = np.array([float(next(it)) for _ in range(8)])
    grav = float(next(it)); C_elastic = int(next(it))
    roumax = float(next(it)); rhow = float(next(it)); gamar = float(next(it))
    slipRateThres = float(next(it))
    xsource, ysource, zsource = float(next(it)), float(next(it)), float(next(it))
    nucR, nucRuptVel, nucdtau0, nucT = [float(next(it)) for _ in range(4)]
    TPV = int(next(it)); C_nuclea = int(next(it)); nucfault = int(next(it))
    friclaw = int(next(it)); timeElapsed0 = float(next(it))

    meshCoor = np.loadtxt(os.path.join(case_dir, 'pydump_meshCoor.txt'))  # (N,3)
    conn_raw = np.loadtxt(os.path.join(case_dir, 'pydump_conn.txt'))      # (E,14)
    conn = conn_raw[:, 0:8].astype(np.int64) - 1  # 0-indexed node ids
    elemType = conn_raw[:, 8].astype(np.int64)
    mat = conn_raw[:, 9:14]  # vp,vs,rho,lam,miu

    eg = np.loadtxt(os.path.join(case_dir, 'pydump_elemgeo.txt'))
    eledet = eg[:, 0]
    eleshp = eg[:, 1:25].reshape(E, 8, 3)          # (E,8,3): stored ((eleshp(j,k)),j=1,3),k=1,8) -> j fastest
    ss = eg[:, 25:31]                               # (E,6)
    phi = eg[:, 31:63].reshape(E, 4, 8)             # stored ((phi(j,k)),j=1,8),k=1,4) -> j fastest -> phi[e,k-1,j-1]

    nodalMassArr = np.loadtxt(os.path.join(case_dir, 'pydump_nodalmass.txt'))  # (NEQ,)
    fnms = np.loadtxt(os.path.join(case_dir, 'pydump_fnms.txt'))              # (N,)
    v1 = np.loadtxt(os.path.join(case_dir, 'pydump_v1.txt'))                  # (NEQ,)

    ndof = np.zeros(N, dtype=np.int64)
    eq_ids = np.zeros((N, 12), dtype=np.int64)  # 0 = "no equation" sink
    with open(os.path.join(case_dir, 'pydump_nodeinfo.txt')) as f:
        for i, line in enumerate(f):
            parts = line.split()
            nd = int(parts[0])
            ndof[i] = nd
            eqs = np.array(parts[2:2 + nd], dtype=np.int64)
            eqs = np.where(eqs > 0, eqs, 0)  # -1 sentinel (fixed bnd) -> sink 0
            eq_ids[i, :nd] = eqs

    fq = np.loadtxt(os.path.join(case_dir, 'pydump_fault.txt'))  # (nftnd, 2+3+3+3+1+100)
    nsmp1 = fq[:, 0].astype(np.int64) - 1
    nsmp2 = fq[:, 1].astype(np.int64) - 1
    un = fq[:, 2:5]; us = fq[:, 5:8]; ud = fq[:, 8:11]
    arn = fq[:, 11]
    fric = fq[:, 12:112]  # columns 1..100 -> index 0..99

    return dict(N=N, E=E, NEQ=NEQ, nen=nen, ned=ned, nftnd=nftnd, ntotft=ntotft,
                nstep=nstep, dt=dt, w=w, rdampk=rdampk, rdampm=rdampm,
                kapa_hg=kapa_hg, R=R, nPML=nPML, vmaxPML=vmaxPML, PMLb=PMLb,
                grav=grav, C_elastic=C_elastic, roumax=roumax, rhow=rhow,
                gamar=gamar, slipRateThres=slipRateThres, xsource=xsource,
                ysource=ysource, zsource=zsource, nucR=nucR, nucRuptVel=nucRuptVel,
                nucdtau0=nucdtau0, nucT=nucT, TPV=TPV, C_nuclea=C_nuclea,
                nucfault=nucfault, friclaw=friclaw, meshCoor=meshCoor, conn=conn,
                elemType=elemType, mat=mat, eledet=eledet, eleshp=eleshp, ss=ss,
                phi=phi, nodalMassArr=nodalMassArr, fnms=fnms, v1_init=v1,
                ndof=ndof, eq_ids=eq_ids, nsmp1=nsmp1, nsmp2=nsmp2, un=un, us=us,
                ud=ud, arn=arn, fric_init=fric)


def region_damp(x, y, z, PMLb, nPML, vmaxPML, R, bound_inclusive):
    """Vectorized port of src/func_lib.f90's `pmlRegionDistance` (post
    a844e92: region-14 now correctly tests y against ymax0, not xmax0;
    z<=zmin0 gap closed). `bound_inclusive` reproduces the Fortran's own
    pre-existing, NOT-unified difference between the two call sites:
    computePMLDampingVector.f90 (node-based velDispUpdate damping) always
    calls with boundInclusive=.true. (>=/<=); assembleGlobalKU.f90's
    calcPMLElemKU (element-center damping) always calls with
    boundInclusive=.false. (strict >/<). Get this flag right per caller,
    not just once -- it changes which points land exactly on a PML shell
    boundary, which is common on this code's structured grids."""
    xmax0, xmin0, ymax0, ymin0, zmin0, maxdx, maxdy, maxdz = PMLb
    n = x.shape[0]
    if bound_inclusive:
        xHi = x >= xmax0; xLo = x <= xmin0; yHi = y >= ymax0; yLo = y <= ymin0
    else:
        xHi = x > xmax0; xLo = x < xmin0; yHi = y > ymax0; yLo = y < ymin0

    d3 = np.where(z <= zmin0, np.abs(z - zmin0), 0.0)

    d1 = np.zeros(n); d2 = np.zeros(n)
    taken = np.zeros(n, dtype=bool)
    for cond, a, b in [
        (xHi & yHi, np.abs(x - xmax0), np.abs(y - ymax0)),          # region 11
        (xHi & yLo, np.abs(x - xmax0), np.abs(y - ymin0)),          # region 12
        (xLo & yLo, np.abs(x - xmin0), np.abs(y - ymin0)),          # region 13
        (xLo & yHi, np.abs(x - xmin0), np.abs(y - ymax0)),          # region 14 (fixed)
    ]:
        sel = cond & ~taken
        d1 = np.where(sel, a, d1); d2 = np.where(sel, b, d2); taken |= sel
    sel = ~taken & xHi & (y > ymin0) & (y < ymax0)                  # region 1_12
    d1 = np.where(sel, np.abs(x - xmax0), d1); d2 = np.where(sel, 0.0, d2); taken |= sel
    sel = ~taken & yLo & (x > xmin0) & (x < xmax0)                  # region 1_23
    d1 = np.where(sel, 0.0, d1); d2 = np.where(sel, np.abs(y - ymin0), d2); taken |= sel
    sel = ~taken & xLo & (y > ymin0) & (y < ymax0)                  # region 1_34
    d1 = np.where(sel, np.abs(x - xmin0), d1); d2 = np.where(sel, 0.0, d2); taken |= sel
    sel = ~taken & yHi & (x > xmin0) & (x < xmax0)                  # region 1_41
    d1 = np.where(sel, 0.0, d1); d2 = np.where(sel, np.abs(y - ymax0), d2); taken |= sel
    # else (middle area 9): d1=d2=0, already the default.

    out = []
    for d, delta in ((d1, nPML * maxdx), (d2, nPML * maxdy), (d3, nPML * maxdz)):
        out.append(3.0 * vmaxPML / 2.0 / delta * np.log(1.0 / R) * (d / delta) ** 2.0)
    return out  # damp1, damp2, damp3


def _c(dN, v):
    """Specialized replacement for np.einsum('ei,ei->e', dN, v): a plain
    broadcast-multiply + axis-sum. Same math as the einsum call it replaces,
    verified to roundoff against the pre-optimization einsum-based run in
    README-parity.md's optimization log."""
    return (dN * v).sum(axis=1)


def run(S, nsteps=None, verbose=True):
    N = S['N']; E = S['E']; dt = S['dt']; w = S['w']; rdampk = S['rdampk']
    conn = S['conn']; elemType = S['elemType']; mat = S['mat']; eledet = S['eledet']
    eleshp = S['eleshp']; ss = S['ss']; phi = S['phi']
    ndof = S['ndof']; eq_ids = S['eq_ids']; NEQ = S['NEQ']
    nsteps = nsteps or S['nstep']

    is_int = elemType == 1
    is_pml = elemType == 2
    E_int = np.nonzero(is_int)[0]
    E_pml = np.nonzero(is_pml)[0]

    velArr = np.zeros((N, 3)); dispArr = np.zeros((N, 3))
    v1 = np.zeros(NEQ + 1)  # index 0 = sink for -1/no-equation dof
    force = np.zeros(NEQ + 1)
    mass = np.concatenate(([1.0], S['nodalMassArr']))  # index0 unused (never divided)

    # ---- loop-invariant setup (hoisted OUT of the time loop; opt #1) ----
    int_nodes = ndof == 3
    idx3_v = eq_ids[int_nodes, 0:3]  # (Nint,3), static

    pml_nodes = np.nonzero(ndof == 12)[0]
    d1n, d2n, d3n = region_damp(S['meshCoor'][pml_nodes, 0], S['meshCoor'][pml_nodes, 1],
                                 S['meshCoor'][pml_nodes, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'], True)
    dampv_pml = np.zeros((pml_nodes.shape[0], 9))
    for k, dk in enumerate((d1n, d2n, d3n)):
        dampv_pml[:, k] = dk; dampv_pml[:, k + 3] = dk; dampv_pml[:, k + 6] = dk
    idx12_v = eq_ids[pml_nodes, :]  # (Npml,12), static
    a9 = 1.0 / dt - dampv_pml / 2.0; b9 = 1.0 / dt + dampv_pml / 2.0  # static per-node damping factors

    lam_i = mat[E_int, 3]; miu_i = mat[E_int, 4]
    dN_i = eleshp[E_int]
    dNx_i, dNy_i, dNz_i = dN_i[:, :, 0], dN_i[:, :, 1], dN_i[:, :, 2]
    constk_w_i = (-eledet[E_int]) * w
    stress_i = np.zeros((E_int.shape[0], 6))
    conn_i = conn[E_int]
    idxIx = eq_ids[conn_i, 0].ravel(); idxIy = eq_ids[conn_i, 1].ravel(); idxIz = eq_ids[conn_i, 2].ravel()

    lam_p = mat[E_pml, 3]; miu_p = mat[E_pml, 4]
    dN_p = eleshp[E_pml]
    dNx_p, dNy_p, dNz_p = dN_p[:, :, 0], dN_p[:, :, 1], dN_p[:, :, 2]
    det_w_p = eledet[E_pml] * w
    conn_p = conn[E_pml]
    xc_p = S['meshCoor'][conn_p].mean(axis=1)
    d1p, d2p, d3p = region_damp(xc_p[:, 0], xc_p[:, 1], xc_p[:, 2], S['PMLb'], S['nPML'],
                                 S['vmaxPML'], S['R'], False)
    a1 = 1.0 / dt - d1p / 2.0; b1 = 1.0 / dt + d1p / 2.0
    a2 = 1.0 / dt - d2p / 2.0; b2 = 1.0 / dt + d2p / 2.0
    a3 = 1.0 / dt - d3p / 2.0; b3 = 1.0 / dt + d3p / 2.0
    s_p = np.zeros((E_pml.shape[0], 15))  # s16..21 stay 0 forever (never written in Fortran either)
    nd_p = ndof[conn_p]; is12_p = nd_p == 12; is3_p = nd_p == 3
    # For nodes with ndof==12: scatter all 12 split components (masked to sink 0 where not is12_p).
    idxP12 = [np.where(is12_p, eq_ids[conn_p, j], 0).ravel() for j in range(12)]
    # For nodes with ndof==3 (mixed-dof PML element corner): scatter grouped sums at slots 0,1,2.
    idxP3 = [np.where(is3_p, eq_ids[conn_p, g], 0).ravel() for g in range(3)]

    # hrglss: node's own ndof decides slot 0-2 (interior) or 9-11 (PML); static, so index once.
    nd_all = ndof[conn]  # (E,8)
    slot0 = np.where(nd_all == 3, 0, 9); slot1 = np.where(nd_all == 3, 1, 10); slot2 = np.where(nd_all == 3, 2, 11)
    idxH0 = np.take_along_axis(eq_ids[conn], slot0[:, :, None], axis=2)[:, :, 0].ravel()
    idxH1 = np.take_along_axis(eq_ids[conn], slot1[:, :, None], axis=2)[:, :, 0].ravel()
    idxH2 = np.take_along_axis(eq_ids[conn], slot2[:, :, None], axis=2)[:, :, 0].ravel()

    nftnd = S['nftnd']; nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    un = S['un']; us = S['us']; ud = S['ud']; arn = S['arn']
    fric = S['fric_init'].copy()
    fnft = np.full(nftnd, 99999.0)
    timeElapsed = 0.0
    slipRateThres = S['slipRateThres']; C_elastic = S['C_elastic']
    idxF_s = [eq_ids[nsmp1, d] for d in range(3)]
    idxF_m = [eq_ids[nsmp2, d] for d in range(3)]

    NEQ1 = NEQ + 1

    for nt in range(1, nsteps + 1):
        timeElapsed += dt

        # ---- velDispUpdate ----
        accel3 = force[idx3_v]
        v1[idx3_v] = v1[idx3_v] + accel3 * dt
        velArr[int_nodes] = v1[idx3_v]
        dispArr[int_nodes] += v1[idx3_v] * dt

        if pml_nodes.size:
            f9 = force[idx12_v[:, 0:9]]
            vold9 = v1[idx12_v[:, 0:9]]
            vnew9 = (f9 + vold9 * a9) / b9
            v1[idx12_v[:, 0:9]] = vnew9  # index0 sink is safe to overwrite; reset below

            f3v9 = force[idx12_v[:, 9:12]]
            v1[idx12_v[:, 9:12]] = v1[idx12_v[:, 9:12]] + f3v9 * dt

            v1[0] = 0.0  # scrub the sink before the sums below re-read it
            vA = v1[idx12_v[:, 0]] + v1[idx12_v[:, 1]] + v1[idx12_v[:, 2]] + v1[idx12_v[:, 9]]
            vB = v1[idx12_v[:, 3]] + v1[idx12_v[:, 4]] + v1[idx12_v[:, 5]] + v1[idx12_v[:, 10]]
            vC = v1[idx12_v[:, 6]] + v1[idx12_v[:, 7]] + v1[idx12_v[:, 8]] + v1[idx12_v[:, 11]]
            has_eq = idx12_v[:, 0] > 0
            velArr[pml_nodes, 0] = np.where(has_eq, vA, 0.0)
            velArr[pml_nodes, 1] = np.where(has_eq, vB, 0.0)
            velArr[pml_nodes, 2] = np.where(has_eq, vC, 0.0)
            dispArr[pml_nodes, 0] = np.where(has_eq, dispArr[pml_nodes, 0] + vA * dt, 0.0)
            dispArr[pml_nodes, 1] = np.where(has_eq, dispArr[pml_nodes, 1] + vB * dt, 0.0)
            dispArr[pml_nodes, 2] = np.where(has_eq, dispArr[pml_nodes, 2] + vC * dt, 0.0)

        force[:] = 0.0
        scat_idx = []; scat_val = []  # opt #2: one combined bincount instead of many np.add.at calls

        # ---- assembleGlobalKU: interior elements ----
        vl = velArr[conn_i]
        sr1 = _c(dNx_i, vl[:, :, 0]); sr2 = _c(dNy_i, vl[:, :, 1]); sr3 = _c(dNz_i, vl[:, :, 2])
        sr4 = _c(dNz_i, vl[:, :, 1]) + _c(dNy_i, vl[:, :, 2])
        sr5 = _c(dNz_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 2])
        sr6 = _c(dNy_i, vl[:, :, 0]) + _c(dNx_i, vl[:, :, 1])
        vol_sr = sr1 + sr2 + sr3
        strr1 = lam_i * vol_sr + 2 * miu_i * sr1
        strr2 = lam_i * vol_sr + 2 * miu_i * sr2
        strr3 = lam_i * vol_sr + 2 * miu_i * sr3
        strr4 = miu_i * sr4; strr5 = miu_i * sr5; strr6 = miu_i * sr6
        stress_i[:, 0] += strr1 * dt; stress_i[:, 1] += strr2 * dt; stress_i[:, 2] += strr3 * dt
        stress_i[:, 3] += strr4 * dt; stress_i[:, 4] += strr5 * dt; stress_i[:, 5] += strr6 * dt
        st1 = constk_w_i * (stress_i[:, 0] + rdampk * strr1)
        st2 = constk_w_i * (stress_i[:, 1] + rdampk * strr2)
        st3 = constk_w_i * (stress_i[:, 2] + rdampk * strr3)
        st4 = constk_w_i * (stress_i[:, 3] + rdampk * strr4)
        st5 = constk_w_i * (stress_i[:, 4] + rdampk * strr5)
        st6 = constk_w_i * (stress_i[:, 5] + rdampk * strr6)
        Fx = dNx_i * st1[:, None] + dNz_i * st5[:, None] + dNy_i * st6[:, None]
        Fy = dNy_i * st2[:, None] + dNz_i * st4[:, None] + dNx_i * st6[:, None]
        Fz = dNz_i * st3[:, None] + dNy_i * st4[:, None] + dNx_i * st5[:, None]
        scat_idx += [idxIx, idxIy, idxIz]
        scat_val += [Fx.ravel(), Fy.ravel(), Fz.ravel()]

        # ---- assembleGlobalKU: PML elements ----
        if E_pml.shape[0]:
            vlp = velArr[conn_p]
            srp1 = _c(dNx_p, vlp[:, :, 0]); srp2 = _c(dNy_p, vlp[:, :, 1]); srp3 = _c(dNz_p, vlp[:, :, 2])
            srp4 = _c(dNz_p, vlp[:, :, 1]) + _c(dNy_p, vlp[:, :, 2])
            srp5 = _c(dNz_p, vlp[:, :, 0]) + _c(dNx_p, vlp[:, :, 2])
            srp6 = _c(dNy_p, vlp[:, :, 0]) + _c(dNx_p, vlp[:, :, 1])
            volp = srp1 + srp2 + srp3
            srate = [lam_p * volp + 2 * miu_p * srp1, lam_p * volp + 2 * miu_p * srp2,
                     lam_p * volp + 2 * miu_p * srp3, miu_p * srp4, miu_p * srp5, miu_p * srp6]

            Dxvx = _c(dNx_p, vlp[:, :, 0]); Dyvy = _c(dNy_p, vlp[:, :, 1]); Dzvz = _c(dNz_p, vlp[:, :, 2])
            Dxvy = _c(dNx_p, vlp[:, :, 1]); Dyvx = _c(dNy_p, vlp[:, :, 0])
            Dxvz = _c(dNx_p, vlp[:, :, 2]); Dzvx = _c(dNz_p, vlp[:, :, 0])
            Dyvz = _c(dNy_p, vlp[:, :, 2]); Dzvy = _c(dNz_p, vlp[:, :, 1])

            s_p[:, 0] = ((lam_p + 2 * miu_p) * Dxvx + a1 * s_p[:, 0]) / b1
            s_p[:, 1] = (lam_p * Dyvy + a2 * s_p[:, 1]) / b2
            s_p[:, 2] = (lam_p * Dzvz + a3 * s_p[:, 2]) / b3
            s_p[:, 3] = (lam_p * Dxvx + a1 * s_p[:, 3]) / b1
            s_p[:, 4] = ((lam_p + 2 * miu_p) * Dyvy + a2 * s_p[:, 4]) / b2
            s_p[:, 5] = (lam_p * Dzvz + a3 * s_p[:, 5]) / b3
            s_p[:, 6] = (lam_p * Dxvx + a1 * s_p[:, 6]) / b1
            s_p[:, 7] = (lam_p * Dyvy + a2 * s_p[:, 7]) / b2
            s_p[:, 8] = ((lam_p + 2 * miu_p) * Dzvz + a3 * s_p[:, 8]) / b3
            s_p[:, 9] = (miu_p * Dxvy + a1 * s_p[:, 9]) / b1
            s_p[:, 10] = (miu_p * Dyvx + a2 * s_p[:, 10]) / b2
            s_p[:, 11] = (miu_p * Dxvz + a1 * s_p[:, 11]) / b1
            s_p[:, 12] = (miu_p * Dzvx + a3 * s_p[:, 12]) / b3
            s_p[:, 13] = (miu_p * Dyvz + a2 * s_p[:, 13]) / b2
            s_p[:, 14] = (miu_p * Dzvy + a3 * s_p[:, 14]) / b3

            sxx = s_p[:, 0] + s_p[:, 1] + s_p[:, 2]
            syy = s_p[:, 3] + s_p[:, 4] + s_p[:, 5]
            szz = s_p[:, 6] + s_p[:, 7] + s_p[:, 8]
            sxy = s_p[:, 9] + s_p[:, 10]
            sxz = s_p[:, 11] + s_p[:, 12]
            syz = s_p[:, 13] + s_p[:, 14]
            s0 = [rdampk * srate[k] for k in range(6)]  # s16..21==0 forever (dead storage in Fortran)

            f1 = -det_w_p[:, None] * dNx_p * sxx[:, None]
            f2 = -det_w_p[:, None] * dNy_p * sxy[:, None]
            f3v = -det_w_p[:, None] * dNz_p * sxz[:, None]
            f4 = -det_w_p[:, None] * dNx_p * sxy[:, None]
            f5 = -det_w_p[:, None] * dNy_p * syy[:, None]
            f6 = -det_w_p[:, None] * dNz_p * syz[:, None]
            f7 = -det_w_p[:, None] * dNx_p * sxz[:, None]
            f8 = -det_w_p[:, None] * dNy_p * syz[:, None]
            f9v = -det_w_p[:, None] * dNz_p * szz[:, None]
            f10 = -det_w_p[:, None] * (dNx_p * s0[0][:, None] + dNz_p * s0[4][:, None] + dNy_p * s0[5][:, None])
            f11 = -det_w_p[:, None] * (dNy_p * s0[1][:, None] + dNz_p * s0[3][:, None] + dNx_p * s0[5][:, None])
            f12 = -det_w_p[:, None] * (dNz_p * s0[2][:, None] + dNy_p * s0[3][:, None] + dNx_p * s0[4][:, None])
            efPML12 = [f1, f2, f3v, f4, f5, f6, f7, f8, f9v, f10, f11, f12]

            for j in range(12):
                scat_idx.append(idxP12[j]); scat_val.append(efPML12[j].ravel())
            grp_sums = [efPML12[0] + efPML12[1] + efPML12[2] + efPML12[9],
                        efPML12[3] + efPML12[4] + efPML12[5] + efPML12[10],
                        efPML12[6] + efPML12[7] + efPML12[8] + efPML12[11]]
            for g in range(3):
                scat_idx.append(idxP3[g]); scat_val.append(grp_sums[g].ravel())

        # ---- hrglss (C_hg==1, all elements) ----
        dl_all = dispArr[conn] + rdampk * velArr[conn]
        for m in range(4):
            phid = (phi[:, m, :, None] * dl_all).sum(axis=1)  # (E,3); replaces einsum('ei,eik->ek', ...)
            r0 = ss[:, 0] * phid[:, 0] + ss[:, 1] * phid[:, 1] + ss[:, 2] * phid[:, 2]
            r1 = ss[:, 1] * phid[:, 0] + ss[:, 3] * phid[:, 1] + ss[:, 4] * phid[:, 2]
            r2 = ss[:, 2] * phid[:, 0] + ss[:, 4] * phid[:, 1] + ss[:, 5] * phid[:, 2]
            fh0 = -(phi[:, m, :] * r0[:, None]).ravel()
            fh1 = -(phi[:, m, :] * r1[:, None]).ravel()
            fh2 = -(phi[:, m, :] * r2[:, None]).ravel()
            scat_idx += [idxH0, idxH1, idxH2]
            scat_val += [fh0, fh1, fh2]

        all_idx = np.concatenate(scat_idx)
        all_val = np.concatenate(scat_val)
        force += np.bincount(all_idx, weights=all_val, minlength=NEQ1)
        force[0] = 0.0  # scrub sink again (bincount may have accumulated masked-out contributions there)

        # ---- faulting (friclaw==1: solveSWTW), C_nuclea==1 nucfault==ift==1 ----
        fx_s = force[idxF_s[0]]; fy_s = force[idxF_s[1]]; fz_s = force[idxF_s[2]]
        fx_m = force[idxF_m[0]]; fy_m = force[idxF_m[1]]; fz_m = force[idxF_m[2]]
        vx_s, vy_s, vz_s = velArr[nsmp1, 0], velArr[nsmp1, 1], velArr[nsmp1, 2]
        vx_m, vy_m, vz_m = velArr[nsmp2, 0], velArr[nsmp2, 1], velArr[nsmp2, 2]
        dx_s, dy_s, dz_s = dispArr[nsmp1, 0], dispArr[nsmp1, 1], dispArr[nsmp1, 2]
        dx_m, dy_m, dz_m = dispArr[nsmp2, 0], dispArr[nsmp2, 1], dispArr[nsmp2, 2]

        def rot(vx, vy, vz, dirvec):
            return vx * dirvec[:, 0] + vy * dirvec[:, 1] + vz * dirvec[:, 2]

        fN_s, fS_s, fD_s = rot(fx_s, fy_s, fz_s, un), rot(fx_s, fy_s, fz_s, us), rot(fx_s, fy_s, fz_s, ud)
        fN_m, fS_m, fD_m = rot(fx_m, fy_m, fz_m, un), rot(fx_m, fy_m, fz_m, us), rot(fx_m, fy_m, fz_m, ud)
        vN_s, vS_s, vD_s = rot(vx_s, vy_s, vz_s, un), rot(vx_s, vy_s, vz_s, us), rot(vx_s, vy_s, vz_s, ud)
        vN_m, vS_m, vD_m = rot(vx_m, vy_m, vz_m, un), rot(vx_m, vy_m, vz_m, us), rot(vx_m, vy_m, vz_m, ud)
        dN_s, dS_s, dD_s = rot(dx_s, dy_s, dz_s, un), rot(dx_s, dy_s, dz_s, us), rot(dx_s, dy_s, dz_s, ud)
        dN_m, dS_m, dD_m = rot(dx_m, dy_m, dz_m, un), rot(dx_m, dy_m, dz_m, us), rot(dx_m, dy_m, dz_m, ud)

        slipN = dN_m - dN_s; slipS = dS_m - dS_s; slipD = dD_m - dD_s
        srN = vN_m - vN_s; srS = vS_m - vS_s; srD = vD_m - vD_s
        srMag = np.sqrt(srN ** 2 + srS ** 2 + srD ** 2)

        fric[:, 70] = slipS; fric[:, 71] = slipD; fric[:, 72] = slipN
        fric[:, 73] = srS; fric[:, 74] = srD
        fric[:, 75] = np.maximum(fric[:, 75], srMag)
        fric[:, 76] = fric[:, 76] + srMag * dt

        massSlave = S['fnms'][nsmp1]; massMaster = S['fnms'][nsmp2]
        totalMass = (massSlave + massMaster) * arn
        Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
              + massSlave * fN_m - massMaster * fN_s) / totalMass + fric[:, 6] * C_elastic
        Ts = (massSlave * massMaster * (vS_m - vS_s) / dt + massSlave * fS_m - massMaster * fS_s) / totalMass \
             + fric[:, 7] * C_elastic
        Td = (massSlave * massMaster * (vD_m - vD_s) / dt + massSlave * fD_m - massMaster * fD_s) / totalMass \
             + fric[:, 48] * C_elastic
        Tmag = np.sqrt(Ts ** 2 + Td ** 2)

        fs = fric[:, 0]; fd = fric[:, 1]; D0 = fric[:, 2]
        slip = fric[:, 76]
        fricCoeff = np.where(np.abs(slip) < 1.0e-10, fs, fs - (fs - fd) * slip / D0)
        fricCoeff = np.where(slip >= D0, fd, fricCoeff)

        fricCoeff = np.minimum(fs, fricCoeff)  # swtwNucleation no-op for TPV==8 (tr default 1e9)

        effNorm = np.where((Tn + fric[:, 5]) > 0.0, 0.0, Tn + fric[:, 5])
        trialShear = fric[:, 3] - fricCoeff * effNorm
        over = Tmag > trialShear
        scale = np.where(over, trialShear / np.where(Tmag == 0, 1.0, Tmag), 1.0)
        Ts = Ts * scale; Td = Td * scale

        xN = Tn * un[:, 0] + Ts * us[:, 0] + Td * ud[:, 0]
        xN_y = Tn * un[:, 1] + Ts * us[:, 1] + Td * ud[:, 1]
        xN_z = Tn * un[:, 2] + Ts * us[:, 2] + Td * ud[:, 2]
        xTrac = np.stack([xN, xN_y, xN_z], axis=1) * arn[:, None]
        initN, initS, initD = fric[:, 6], fric[:, 7], fric[:, 48]
        xInit = np.stack([initN * un[:, 0] + initS * us[:, 0] + initD * ud[:, 0],
                           initN * un[:, 1] + initS * us[:, 1] + initD * ud[:, 1],
                           initN * un[:, 2] + initS * us[:, 2] + initD * ud[:, 2]], axis=1) * arn[:, None]
        delta_f = xTrac - xInit * C_elastic
        flt_idx = np.concatenate([idxF_s[0], idxF_m[0], idxF_s[1], idxF_m[1], idxF_s[2], idxF_m[2]])
        flt_val = np.concatenate([delta_f[:, 0], -delta_f[:, 0], delta_f[:, 1], -delta_f[:, 1],
                                   delta_f[:, 2], -delta_f[:, 2]])
        force += np.bincount(flt_idx, weights=flt_val, minlength=NEQ1)
        force[0] = 0.0

        fric[:, 77] = Tn; fric[:, 78] = Ts; fric[:, 79] = Td

        need = fnft > 5000.0
        fnft = np.where(need & (srMag >= slipRateThres), timeElapsed, fnft)

        force[0] = 0.0
        force[1:] = force[1:] / mass[1:]
        if verbose and nt % 20 == 0:
            print('step', nt, '/', nsteps, 'max|v|', np.abs(velArr).max())

    return dict(velArr=velArr, dispArr=dispArr, fnft=fnft, fric=fric, force=force)

