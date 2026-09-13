"""
Shared static-state loader + PML region-damping helper, used by every
eqdyna Python port (NumPy and JAX, all three friction families). Moved out
of port.py so port_rsf.py/port_tp.py/port_jax.py/port_rsf_jax.py/
port_tp_jax.py import ONE copy instead of `from port import load,
region_damp` (which worked, but made port.py -- the tpv8/friclaw=1 port --
an accidental shared-infrastructure module for every other friclaw).

`load()` reads the Fortran `pydump_*` files written by `src/pydump.f90`
(mesh/mass/shape-function/PML/fault-static state) -- see python/README-
parity.md for the full field-by-field provenance.

`region_damp()` is a vectorized NumPy port of src/func_lib.f90's
`pmlRegionDistance`. As of pathway_forward.md item 13 (2026-09-13),
`pmlRegionDistance` classifies BOTH its callers (node-based, from
computePMLDampingVector.f90/velDispUpdate, and element-center-based, from
assembleGlobalKU.f90's calcPMLElemKU) with INCLUSIVE (>=/<=) boundary
tests -- the two call sites are no longer different, and the subroutine no
longer takes a boundInclusive argument at all. This function matches that:
no bound_inclusive parameter, always inclusive. (Before this change, this
port had FOLLOWED an earlier version of the Fortran where the two call
sites genuinely differed -- element centers strict, nodes inclusive. That
distinction no longer exists in src/func_lib.f90 and carrying it here would
be reproducing a bug that was fixed upstream. Re-run parity after this
change, per the risk list -- element centers are additionally guaranteed
never to land exactly on a PML bound by checkPMLAlignment (called from
meshgen), so in practice this was already a no-op for the element-center
caller; parity confirms it.)
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


def region_damp(x, y, z, PMLb, nPML, vmaxPML, R):
    """Vectorized port of src/func_lib.f90's `pmlRegionDistance` (current
    master: region-14 tests y against ymax0, not xmax0 -- pathway_forward.md
    item 8's bug, fixed upstream in a844e92; boundary tests are INCLUSIVE
    (>=/<=) for both callers -- pathway_forward.md item 13, standardized
    2026-09-13. There is no longer a bound_inclusive distinction to
    reproduce."""
    xmax0, xmin0, ymax0, ymin0, zmin0, maxdx, maxdy, maxdz = PMLb
    n = x.shape[0]
    xHi = x >= xmax0; xLo = x <= xmin0; yHi = y >= ymax0; yLo = y <= ymin0

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
