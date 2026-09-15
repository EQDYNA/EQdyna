"""pydump.py <- src/pydump.f90 (the reader for what it writes).

src/pydump.f90 is Fortran-side INSTRUMENTATION: called once before the time
loop, it dumps mesh, mass, shape-function, PML and fault-static state as
pydump_*.txt. This module reads those files back.

It is used by ONE thing -- python/run_parity.py, the per-step parity
diagnostic that answers "at which step does the port start to diverge from
Fortran's own internal state?". The solver itself never touches it:
eqdyna3d.build_solver_state builds every one of these fields from the case
inputs instead (readInputFiles.py + meshgen.py + assembleGlobalMass.py), so a
normal run needs no Fortran in the loop at all.

Kept separate from readInputFiles.py for that reason: that module reads the
CASE INPUTS the Fortran also reads, this one reads the Fortran's debug dump.
Merging them would put a diagnostic-only dependency in the solver's path.
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
