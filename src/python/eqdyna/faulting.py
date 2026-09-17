"""faulting.py <- src/faulting.f90.

ONE faulting routine, and the friclaw dispatch is INSIDE it, exactly as
faulting.f90:17-18 does it:

    if (friclaw<=2) call solveSWTW(...)
    if (friclaw>=3) call solveRSF(...)

friclaw is static per run, so under jit that `if` is a TRACE-TIME branch on
a Python int and costs nothing at runtime -- the same reason kernels_*.py's
`if C_elastic == 0:` is free.

WHY THE DISPATCH BEING HERE MATTERS, CONCRETELY
The port this replaces hoisted the friclaw branch to MODULE level: port.py
(friclaw 1), port_rsf.py (4), port_tp.py (5), each duplicated again per
backend -- six copies of one loop. swtwNucleation's forced-rupture branch
was then fixed in port.py alone, so test.tpv29 passed on python-numpy and
failed on python-jax at the identical value: one bug, two fixes, one of
them never made. There is now one place for that fix to live.

nucleation.py, created as a separate module to share swtwNucleation between
the two backends, is folded in here instead: swtwNucleation lives in
faulting.f90, and a module-per-subroutine is the same mistake at finer
grain. (It was also never wired in -- port.py kept its own inline copy.)
"""
from . import backend as B
from . import fric as F
from . import globalvar as gv


# ---------------------------------------------------------------------------
# loop-invariant fault setup
# ---------------------------------------------------------------------------
def build(S):
    """Per-fault-node invariants. Host/numpy side, like assembleGlobalKU.build.

    `nuc_radius` is the distance from each fault node to the hypocentre,
    computed from the SLAVE node coordinate to match faulting.f90:388's
    meshCoor(:, nsmp(1,...)). It is geometry, so it is computed once here
    rather than every step as the scalar Fortran does.
    """
    import numpy as np

    eq_ids = S['eq_ids']
    nsmp1 = S['nsmp1']; nsmp2 = S['nsmp2']
    coords = S['meshCoor'][nsmp1]
    nuc_radius = np.sqrt((coords[:, 0] - S['xsource']) ** 2
                         + (coords[:, 1] - S['ysource']) ** 2
                         + (coords[:, 2] - S['zsource']) ** 2)
    massSlave = S['fnms'][nsmp1]
    massMaster = S['fnms'][nsmp2]
    return dict(
        nftnd=S['nftnd'], nsmp1=nsmp1, nsmp2=nsmp2,
        un=S['un'], us=S['us'], ud=S['ud'], arn=S['arn'],
        idxF_s=[eq_ids[nsmp1, d] for d in range(3)],
        idxF_m=[eq_ids[nsmp2, d] for d in range(3)],
        massSlave=massSlave, massMaster=massMaster,
        mr=massMaster * massSlave / (massMaster + massSlave),   # reduced mass
        nuc_radius=nuc_radius,
        # static per-run scalars the step branches on at trace time
        friclaw=S['friclaw'], TPV=S['TPV'], C_elastic=S['C_elastic'],
        C_nuclea=S.get('C_nuclea', 0), nucfault=S.get('nucfault', 1),
        nucR=S.get('nucR', 0.0), nucRuptVel=S.get('nucRuptVel', 0.0),
        nucT=S.get('nucT', 0.0), nucdtau0=S.get('nucdtau0', 0.0),
        insertFaultType=S.get('insertFaultType', 0),
        slipRateThres=S['slipRateThres'],
    )


def nucleation_enabled(finv):
    """Whether swtwNucleation/rsfNucleation does anything for this case.

    faulting.f90:125 and :151 both gate on `C_nuclea==1 .and. ift==nucfault`.
    ntotft==1 throughout this port, so ift is always 1.
    """
    return finv['C_nuclea'] == 1 and finv['nucfault'] == 1


# ---------------------------------------------------------------------------
# faulting.f90:382-417  swtwNucleation
# ---------------------------------------------------------------------------
def forced_rupture_time(xp, finv):
    """tr per fault node -- faulting.f90:392-404. Loop-invariant (geometry
    only), so the caller computes it once outside the step.

    THE FORMULA (TPV29 spec, TPV29_30_Description_v06 Part 6; TPV36/37/201
    reuse it):

        tr = (r + 0.081*nucR*(1/(1-(r/nucR)^2) - 1)) / (0.7 * 3464)

    with the three constants taken from globalvar.f90 VERBATIM, including
    the truncated 3464 m/s shear-wave speed -- this must reproduce the
    reference's arithmetic, not a more accurate version of it.

    Outside that TPV set tr stays at its 1.0e9 default and the whole thing
    degenerates to min(fs, fricCoeff). Implementing ONLY that degenerate
    branch is what made test.tpv29 nucleate 0 of 3321 fault nodes against a
    reference of 2974 -- not a tolerance question, nothing nucleated at all.

    At r == nucR exactly the taper divides by zero and Fortran yields +Inf,
    i.e. that node never forces. Reproduced rather than guarded, so the edge
    matches.
    """
    radius = finv['nuc_radius']; nucR = finv['nucR']; TPV = finv['TPV']
    never = xp.full(radius.shape, 1.0e9)
    if nucR <= 0.0:
        return never
    inside = radius <= nucR
    # TPV30 is TPV29 with off-fault Drucker-Prager viscoplasticity added
    # (spec p.11: "the material properties are the only difference"); p.16
    # states Parts 5/6 (friction + this nucleation formula) apply to
    # "Benchmarks TPV29 and TPV30" identically -- same as faulting.f90:405.
    if TPV in (29, 30, 36, 37, 201):
        ratio = xp.where(inside, radius / nucR, 0.0)
        taper = 1.0 / (1.0 - ratio ** 2) - 1.0
        tr = (radius + gv.NUC_TAPER_COEF * nucR * taper) / (
            gv.NUC_VR_TO_VS * gv.NUC_VS_FIXED)
        return xp.where(inside, tr, never)
    if TPV == 202:
        return xp.where(inside, radius / finv['nucRuptVel'], never)
    return never


def swtwNucleation(xp, fric, fricCoeff, tr, timeElapsed):
    """faulting.f90:406-413 -- drive mu from fs toward fd over t0 after tr.

    With tr at its 1e9 default this reduces EXACTLY to min(fs, fricCoeff):
    timeElapsed < tr, so tc == 0.0, so fs + (fd-fs)*0.0 == fs.

    The middle arm divides by TW_T0, which is 0.0 for cases that do not use
    time-weakening. xp.where discards that arm's value, but numpy still
    evaluates it and warns. The division is left unguarded deliberately:
    guarding it would mean substituting a value into the arm that IS taken
    for time-weakening cases, and that arm has to match the Fortran exactly.
    """
    fs = fric[:, gv.SW_FS]; fd = fric[:, gv.SW_FD]; t0 = fric[:, gv.TW_T0]
    tc = xp.where(timeElapsed < tr, 0.0,
                  xp.where(timeElapsed < tr + t0, (timeElapsed - tr) / t0, 1.0))
    return xp.minimum(fs + (fd - fs) * tc, fricCoeff)


# ---------------------------------------------------------------------------
# faulting.f90:54-130  getNsdSlipSliprateTraction
# ---------------------------------------------------------------------------
def _rot(vx, vy, vz, dirvec):
    """Project an xyz nodal vector onto one fault-local direction."""
    return vx * dirvec[:, 0] + vy * dirvec[:, 1] + vz * dirvec[:, 2]


def getNsdSlipSliprateTraction(xp, finv, fric, velArr, dispArr, force, dt):
    """faulting.f90:54-130 -- slip, slip rate and traction on every fault
    split-node pair, in fault-local (normal, strike, dip) coordinates.

    Returns (fric, comps) where comps carries everything solveSWTW/solveRSF
    need. The Fortran writes fric slots 71-77 here, BEFORE solveRSF redefines
    the slip-rate magnitude as strike+dip only, so those slots always hold
    the 3-component magnitude -- reproduced by writing them here and
    recomputing the 2-component magnitude inside solveRSF.
    """
    nsmp1, nsmp2 = finv['nsmp1'], finv['nsmp2']
    un, us, ud, arn = finv['un'], finv['us'], finv['ud'], finv['arn']
    idxF_s, idxF_m = finv['idxF_s'], finv['idxF_m']

    fx_s = force[idxF_s[0]]; fy_s = force[idxF_s[1]]; fz_s = force[idxF_s[2]]
    fx_m = force[idxF_m[0]]; fy_m = force[idxF_m[1]]; fz_m = force[idxF_m[2]]
    vx_s, vy_s, vz_s = velArr[nsmp1, 0], velArr[nsmp1, 1], velArr[nsmp1, 2]
    vx_m, vy_m, vz_m = velArr[nsmp2, 0], velArr[nsmp2, 1], velArr[nsmp2, 2]
    dx_s, dy_s, dz_s = dispArr[nsmp1, 0], dispArr[nsmp1, 1], dispArr[nsmp1, 2]
    dx_m, dy_m, dz_m = dispArr[nsmp2, 0], dispArr[nsmp2, 1], dispArr[nsmp2, 2]

    fN_s = _rot(fx_s, fy_s, fz_s, un); fS_s = _rot(fx_s, fy_s, fz_s, us)
    fD_s = _rot(fx_s, fy_s, fz_s, ud)
    fN_m = _rot(fx_m, fy_m, fz_m, un); fS_m = _rot(fx_m, fy_m, fz_m, us)
    fD_m = _rot(fx_m, fy_m, fz_m, ud)
    vN_s = _rot(vx_s, vy_s, vz_s, un); vS_s = _rot(vx_s, vy_s, vz_s, us)
    vD_s = _rot(vx_s, vy_s, vz_s, ud)
    vN_m = _rot(vx_m, vy_m, vz_m, un); vS_m = _rot(vx_m, vy_m, vz_m, us)
    vD_m = _rot(vx_m, vy_m, vz_m, ud)
    dN_s = _rot(dx_s, dy_s, dz_s, un); dS_s = _rot(dx_s, dy_s, dz_s, us)
    dD_s = _rot(dx_s, dy_s, dz_s, ud)
    dN_m = _rot(dx_m, dy_m, dz_m, un); dS_m = _rot(dx_m, dy_m, dz_m, us)
    dD_m = _rot(dx_m, dy_m, dz_m, ud)

    slipN = dN_m - dN_s; slipS = dS_m - dS_s; slipD = dD_m - dD_s
    srN = vN_m - vN_s; srS = vS_m - vS_s; srD = vD_m - vD_s
    srMag = xp.sqrt(srN ** 2 + srS ** 2 + srD ** 2)   # 3-component, faulting.f90:101

    fric = B.setat(xp, fric, (slice(None), gv.SLIP_STRIKE), slipS)
    fric = B.setat(xp, fric, (slice(None), gv.SLIP_DIP), slipD)
    fric = B.setat(xp, fric, (slice(None), gv.SLIP_NORM), slipN)
    fric = B.setat(xp, fric, (slice(None), gv.SLIPRATE_STRIKE), srS)
    fric = B.setat(xp, fric, (slice(None), gv.SLIPRATE_DIP), srD)
    fric = B.setat(xp, fric, (slice(None), gv.SLIPRATE_MAX),
                   xp.maximum(fric[:, gv.SLIPRATE_MAX], srMag))
    fric = B.setat(xp, fric, (slice(None), gv.CUM_SLIP),
                   fric[:, gv.CUM_SLIP] + srMag * dt)

    massSlave = finv['massSlave']; massMaster = finv['massMaster']
    totalMass = (massSlave + massMaster) * arn
    C_elastic = finv['C_elastic']
    Tn = (massSlave * massMaster * ((vN_m - vN_s) + (dN_m - dN_s) / dt) / dt
          + massSlave * fN_m - massMaster * fN_s) / totalMass \
        + fric[:, gv.INIT_NORM] * C_elastic
    Ts = (massSlave * massMaster * (vS_m - vS_s) / dt
          + massSlave * fS_m - massMaster * fS_s) / totalMass \
        + fric[:, gv.INIT_STRIKE_SHEAR] * C_elastic
    Td = (massSlave * massMaster * (vD_m - vD_s) / dt
          + massSlave * fD_m - massMaster * fD_s) / totalMass \
        + fric[:, gv.INIT_DIP_SHEAR] * C_elastic

    comps = dict(Tn=Tn, Ts=Ts, Td=Td, srMag=srMag,
                 slipN=slipN, slipS=slipS, slipD=slipD,
                 srN=srN, srS=srS, srD=srD,
                 fx_s=fx_s, fy_s=fy_s, fz_s=fz_s,
                 fx_m=fx_m, fy_m=fy_m, fz_m=fz_m,
                 vx_s=vx_s, vy_s=vy_s, vz_s=vz_s,
                 vx_m=vx_m, vy_m=vy_m, vz_m=vz_m)
    return fric, comps


# ---------------------------------------------------------------------------
# faulting.f90:132-180  solveSWTW   (friclaw 1 and 2)
# ---------------------------------------------------------------------------
def solveSWTW(xp, finv, fric, fnft, comps, force, timeElapsed, tr, friclaw):
    """faulting.f90:132-180 -- slip-weakening (friclaw 1) or time-weakening
    (friclaw 2) friction, then the fault-node force update.

    The friclaw 1-vs-2 choice is a TRACE-TIME Python branch on a static int
    (faulting.f90:144-149), so it costs nothing inside jit.
    """
    un, us, ud, arn = finv['un'], finv['us'], finv['ud'], finv['arn']
    Tn, Ts, Td = comps['Tn'], comps['Ts'], comps['Td']

    if friclaw == 1:
        fricCoeff = F.slip_weak(xp, fric[:, gv.CUM_SLIP], fric)
    elif friclaw == 2:
        # faulting.f90:147 -- trupt = timeElapsed - fnft. An unruptured node
        # still carries the 99999.0 sentinel, giving a large NEGATIVE trupt,
        # so time_weak's first arm returns fs and no special case is needed.
        fricCoeff = F.time_weak(xp, timeElapsed - fnft, fric)
    else:
        raise ValueError('solveSWTW: friclaw must be 1 or 2, got %r' % (friclaw,))

    if tr is not None:
        fricCoeff = swtwNucleation(xp, fric, fricCoeff, tr, timeElapsed)

    # faulting.f90:153-158. Slot 6 (the retired "surface pore pressure" field,
    # globalvar.f90:27) is NOT added: the Fortran removed those dead reads and
    # the column is identically 0.0 -- verified column-wise on the real case,
    # not assumed. The ports this replaces still carried `Tn + fric[:, 5]`.
    effNorm = xp.where(Tn > 0.0, 0.0, Tn)
    trialShear = fric[:, gv.COHESION] - fricCoeff * effNorm
    Tmag = xp.sqrt(Ts ** 2 + Td ** 2)
    over = Tmag > trialShear
    # DEVIATION FROM THE REFERENCE, pre-existing and deliberately PRESERVED
    # so that this restructure changes nothing on its own: faulting.f90:162
    # computes Ts*trialShear/Tmag, i.e. (Ts*trialShear)/Tmag, while this
    # computes Ts*(trialShear/Tmag). Different associations of the same
    # product, not bit-identical. Correcting it is a separate gated change,
    # not something to fold into a restructure whose claim is that it moved
    # nothing.
    scale = xp.where(over, trialShear / xp.where(Tmag == 0, 1.0, Tmag), 1.0)
    Ts = Ts * scale
    Td = Td * scale

    xT0 = Tn * un[:, 0] + Ts * us[:, 0] + Td * ud[:, 0]
    xT1 = Tn * un[:, 1] + Ts * us[:, 1] + Td * ud[:, 1]
    xT2 = Tn * un[:, 2] + Ts * us[:, 2] + Td * ud[:, 2]
    xTrac = xp.stack([xT0, xT1, xT2], axis=1) * arn[:, None]
    iN = fric[:, gv.INIT_NORM]; iS = fric[:, gv.INIT_STRIKE_SHEAR]
    iD = fric[:, gv.INIT_DIP_SHEAR]
    xInit = xp.stack([iN * un[:, 0] + iS * us[:, 0] + iD * ud[:, 0],
                      iN * un[:, 1] + iS * us[:, 1] + iD * ud[:, 1],
                      iN * un[:, 2] + iS * us[:, 2] + iD * ud[:, 2]],
                     axis=1) * arn[:, None]
    delta_f = xTrac - xInit * finv['C_elastic']

    # faulting.f90:174-177: slave += delta, master -= delta. For a split-node
    # fault each node belongs to exactly one pair, so these indices are unique
    # within each direction and one scatter over the concatenated blocks
    # accumulates exactly what six separate ones would.
    idxF_s, idxF_m = finv['idxF_s'], finv['idxF_m']
    flt_idx = xp.concatenate([idxF_s[0], idxF_m[0], idxF_s[1], idxF_m[1],
                              idxF_s[2], idxF_m[2]])
    flt_val = xp.concatenate([delta_f[:, 0], -delta_f[:, 0],
                              delta_f[:, 1], -delta_f[:, 1],
                              delta_f[:, 2], -delta_f[:, 2]])
    force = B.scatter_add(xp, force, flt_idx, flt_val)
    force = B.setat(xp, force, 0, 0.0)

    fric = B.setat(xp, fric, (slice(None), gv.TRACT_NORM), Tn)
    fric = B.setat(xp, fric, (slice(None), gv.TRACT_STRIKE), Ts)
    fric = B.setat(xp, fric, (slice(None), gv.TRACT_DIP), Td)
    return fric, force, comps['srMag']


# ---------------------------------------------------------------------------
# faulting.f90:302-314  storeRuptureTime
# ---------------------------------------------------------------------------
def storeRuptureTime(xp, finv, fnft, srMag, timeElapsed):
    """faulting.f90:302-314 -- first time a node reaches slipRateThres.
    `fnft > 5000.0` is the Fortran's own test against the 99999.0 sentinel."""
    need = fnft > gv.FNFT_UNRUPTURED_ABOVE
    return xp.where(need & (srMag >= finv['slipRateThres']), timeElapsed, fnft)


# ---------------------------------------------------------------------------
# faulting.f90:3-25  faulting  -- THE dispatch
# ---------------------------------------------------------------------------
def faulting(xp, finv, fric, fnft, velArr, dispArr, force, dt, timeElapsed,
             tr, nt):
    """faulting.f90:3-25, with the friclaw dispatch INSIDE it (lines 17-18).

    finv['friclaw'] is a static Python int, so under jit both branches are
    resolved at trace time and only one is ever compiled.
    """
    friclaw = finv['friclaw']
    fric, comps = getNsdSlipSliprateTraction(xp, finv, fric, velArr, dispArr,
                                             force, dt)
    if friclaw <= 2:
        fric, force, srMag = solveSWTW(xp, finv, fric, fnft, comps, force,
                                       timeElapsed, tr, friclaw)
    else:
        fric, force, srMag = solveRSF(xp, finv, fric, comps, force,
                                      timeElapsed, dt, nt, friclaw)
    fnft = storeRuptureTime(xp, finv, fnft, srMag, timeElapsed)
    return fric, fnft, force


# ---------------------------------------------------------------------------
# faulting.f90:340-380  rsfNucleation
# ---------------------------------------------------------------------------
def rsfNucleation(xp, finv, fric, Tn, Ts, Td, srS, srD, timeElapsed, nt):
    """faulting.f90:340-380 -- the smooth space-time stress perturbation that
    starts an RSF rupture. Returns (fric, Ts).

    F and G are the Fortran's own expressions verbatim, including that F is
    0 outside the patch and G is 1 after nucT, and including the singular
    forms at radius==nucR and timeElapsed==0 that the Fortran evaluates
    without guarding.
    """
    radius = finv['nuc_radius']; nucR = finv['nucR']; nucT = finv['nucT']
    TPV = finv['TPV']
    Fr = xp.where(radius < nucR, xp.exp(radius ** 2 / (radius ** 2 - nucR ** 2)), 0.0)
    G = xp.where(timeElapsed <= nucT,
                 xp.exp((timeElapsed - nucT) ** 2
                        / (timeElapsed * (timeElapsed - 2.0 * nucT))), 1.0)

    if TPV in (104, 105):
        dtau = finv['nucdtau0'] * Fr * G
    elif TPV == 2802:
        # drv.a6, C_elastic==0 (plastic). On the VERY FIRST step only,
        # faulting.f90:365-373 re-derives the RSF state variable and theta_pc
        # from the ACTUAL post-elastic-solve traction: on_fault_vars_input.nc's
        # initial STATE/THETA_PC are set up for the C_elastic==1 perturbation
        # convention and are not consistent with C_elastic==0's absolute-stress
        # one. Omitting this leaves the friction law at the wrong operating
        # point and rupture never nucleates -- found by direct Fortran-vs-python
        # cross-check: identical Tn/Ts/Td at step 1, then python stalls.
        ttao = xp.sqrt(Ts ** 2 + Td ** 2)
        backSliprate = xp.sqrt((srS + fric[:, gv.VINI_X]) ** 2
                               + (srD + fric[:, gv.VINI_Z]) ** 2)
        a0 = fric[:, gv.RSF_A]; v00 = fric[:, gv.RSF_V0]
        first = a0 * xp.log(2.0 * v00 / backSliprate
                            * xp.sinh(ttao / xp.abs(Tn) / a0))
        is1 = (nt == 1)
        fric = B.setat(xp, fric, (slice(None), gv.STATE),
                       xp.where(is1, first, fric[:, gv.STATE]))
        fric = B.setat(xp, fric, (slice(None), gv.THETA_PC),
                       xp.where(is1, xp.abs(Tn), fric[:, gv.THETA_PC]))
        fric = B.setat(xp, fric, (slice(None), gv.NUC_DTAU0),
                       xp.where(is1, finv['nucdtau0'], fric[:, gv.NUC_DTAU0]))
        dtau = fric[:, gv.NUC_DTAU0] * Fr * G
    else:
        dtau = 0.0 * Fr
    return fric, Ts + dtau


# ---------------------------------------------------------------------------
# faulting.f90:433-496  NewtonRaphson
# ---------------------------------------------------------------------------
def NewtonRaphson(xp, friclaw, fric, v_trial, state0, thetaPc0, tnrm,
                  trialMag, T_coeff, dt):
    """faulting.f90:433-496 -- solve for next-step slip rate.

    A fixed 20-iteration loop with a per-node `done` mask, which reproduces
    the scalar Fortran's `exit` EXACTLY rather than approximately. The
    argument: state/theta_pc/xmu/taoc/rsfeq/drsfeqdv are PURE functions of
    (v_trial, the FROZEN state0/thetaPc0 baseline, and fixed inputs) -- they
    are not accumulated across iterations -- so once a node's v_trial is
    frozen, every later iteration recomputes bit-identical values for it.
    Masking the v_trial update with the PRE-iteration flag (so a node that
    converges on this iteration does NOT receive this iteration's update) is
    what matches Fortran's `exit` landing before the update.

    A node that never converges in 20 iterations DOES receive the 20th
    update, because Fortran's `do` loop runs its full body on iv==ivmax when
    it does not exit. That is reproduced, not corrected.
    """
    done = xp.zeros(v_trial.shape, dtype=bool)
    state_out = state0
    thetaPc_out = thetaPc0
    taoc_new = xp.zeros(v_trial.shape)
    ivmax = 20

    for _ in range(ivmax):
        if friclaw == 3:
            xmu, dxmudv, state_this = F.rate_state_ageing_law(
                xp, v_trial, state0, fric, dt)
        else:
            xmu, dxmudv, state_this = F.rate_state_slip_law(
                xp, v_trial, state0, fric, dt)

        if friclaw < 5:
            # thetaPcTmp is re-seeded from the frozen baseline every iteration
            # (faulting.f90:460), not carried forward.
            thetaPc_this, _ = F.rate_state_normal_stress(
                xp, v_trial, thetaPc0, tnrm, fric, dt)
            taoc_this = xmu * thetaPc_this
            drsfeqdv = 1.0 + T_coeff * (dxmudv * thetaPc_this) * 0.5
        else:
            # faulting.f90:466-468. MIN(tnrm, 0) is written literally to match
            # the source's own guard, even though tnrm is already clamped <= 0
            # upstream, so it is a no-op here rather than a divergence.
            thetaPc_this = thetaPc0
            tn0 = xp.minimum(tnrm, 0.0)
            taoc_this = fric[:, gv.COHESION] - xmu * tn0
            drsfeqdv = 1.0 + T_coeff * (-dxmudv * tn0) * 0.5

        rsfeq = v_trial + T_coeff * (taoc_this * 0.5 - trialMag)
        # Both exit criteria, ANDed, exactly as faulting.f90:474.
        exit_now = ((xp.abs(rsfeq / drsfeqdv) < 1.0e-14 * xp.abs(v_trial))
                    & (xp.abs(rsfeq) < 1.0e-6 * xp.abs(v_trial)) & (~done))
        done_after = done | exit_now
        newSliprate = v_trial - rsfeq / drsfeqdv
        v_next = xp.where(newSliprate <= 0.0, v_trial / 2.0, newSliprate)
        v_trial = xp.where(done_after, v_trial, v_next)
        state_out = state_this
        thetaPc_out = thetaPc_this
        taoc_new = taoc_this
        done = done_after

    return v_trial, taoc_new, state_out, thetaPc_out


# ---------------------------------------------------------------------------
# faulting.f90:182-300  solveRSF   (friclaw 3, 4 and 5)
# ---------------------------------------------------------------------------
def solveRSF(xp, finv, fric, comps, force, timeElapsed, dt, nt, friclaw):
    """faulting.f90:182-300 -- rate-and-state friction.

    One function for friclaw 3 (ageing), 4 (slip law with strong rate
    weakening) and 5 (4 + thermal pressurization). Every friclaw branch below
    is a TRACE-TIME Python `if` on a static int, exactly like the Fortran's.

    Two subtleties this reproduces deliberately:

    1. solveRSF REDEFINES the slip-rate magnitude as strike+dip only
       (faulting.f90:218), dropping the normal component -- and it is THAT
       value, not the 3-component one, that storeRuptureTime then sees. The
       fric slots 71-77 written back in getNsdSlipSliprateTraction keep the
       3-component magnitude, because they were written before this point.

    2. The two "pre-loop" constitutive calls mutate fric(20)/fric(23) as a
       side effect in the Fortran, but those mutations are provably discarded:
       both slots are overwritten immediately after NewtonRaphson returns,
       from NewtonRaphson's own state, which was seeded from a snapshot taken
       BEFORE the pre-loop call. Only the pre-loop xmu survives, via
       taoc_old. theta_pc_dot (fric(24)) is fully dead and is not computed.
    """
    un, us, ud, arn = finv['un'], finv['us'], finv['ud'], finv['arn']
    Tn, Ts, Td = comps['Tn'], comps['Ts'], comps['Td']
    srN0, srS0, srD0 = comps['srN'], comps['srS'], comps['srD']
    slipN0 = comps['slipN']

    if nucleation_enabled(finv):
        fric, Ts = rsfNucleation(xp, finv, fric, Tn, Ts, Td, srS0, srD0,
                                 timeElapsed, nt)

    # --- normal traction: TP offset, geometry caps, positivity clamp ---
    if friclaw == 5:
        tnrm = Tn + fric[:, gv.TP_NORM_TP]        # faulting.f90:194-195
    else:
        tnrm = Tn
    if finv['insertFaultType'] > 0 and finv['C_elastic'] == 1:
        # faulting.f90:201-208. min_norm > max_norm, so the two arms are
        # mutually exclusive and values strictly between max_norm and
        # min_norm pass through UNTOUCHED. That gap is the Fortran's, kept.
        #
        # This clamp was present in the friclaw-4 port and MISSING from the
        # friclaw-5 one, which is why the pre-restructure entry point refused
        # insertFaultType>0 combined with friclaw==5. Unifying solveRSF
        # closes that gap by construction rather than by a second edit.
        tnrm = xp.where(tnrm >= gv.MIN_NORM, gv.MIN_NORM,
                        xp.where(tnrm <= gv.MAX_NORM, gv.MAX_NORM, tnrm))
    tnrm = xp.where(tnrm > 0.0, 0.0, tnrm)        # faulting.f90:211

    # --- background (creep) slip rate on top, faulting.f90:215-218 ---
    vini_n = fric[:, gv.VINI_N]; vini_s = fric[:, gv.VINI_X]
    vini_d = fric[:, gv.VINI_Z]
    slipN = slipN0 + vini_n * timeElapsed
    srN = srN0 + vini_n
    srS = srS0 + vini_s
    srD = srD0 + vini_d
    srMag_sd = xp.sqrt(srS ** 2 + srD ** 2)       # RSF-redefined magnitude

    v_trial = srMag_sd
    state0 = fric[:, gv.STATE]
    thetaPc0 = fric[:, gv.THETA_PC]

    # --- pre-loop xmu -> taoc_old (faulting.f90:237-251) ---
    if friclaw == 3:
        xmu_pre, _, _ = F.rate_state_ageing_law(xp, v_trial, state0, fric, dt)
    else:
        xmu_pre, _, _ = F.rate_state_slip_law(xp, v_trial, state0, fric, dt)
    if friclaw == 5:
        taoc_old = fric[:, gv.COHESION] - xmu_pre * tnrm
    else:
        taoc_old = xmu_pre * thetaPc0

    mr = finv['mr']
    T_coeff = arn * dt / mr
    trialS = Ts - taoc_old * 0.5 * (srS / v_trial) + vini_s / T_coeff
    trialD = Td - taoc_old * 0.5 * (srD / v_trial) + vini_d / T_coeff
    trialMag = xp.sqrt(trialS ** 2 + trialD ** 2)

    v_trial, taoc_new, state_out, thetaPc_out = NewtonRaphson(
        xp, friclaw, fric, v_trial, state0, thetaPc0, tnrm, trialMag,
        T_coeff, dt)

    # faulting.f90:491 -- creeping-rate floor. Only TPV==105 reaches it; for
    # every other TPV this is a documented no-op, evaluated as a trace-time
    # Python branch so it is not merely masked.
    if finv['TPV'] == 105:
        creep = fric[:, gv.CREEP_VMIN]
        v_trial = xp.where(v_trial < creep, creep, v_trial)

    fric = B.setat(xp, fric, (slice(None), gv.STATE), state_out)
    if friclaw < 5:
        # faulting.f90:264. friclaw==5 never evolves theta_pc (NewtonRaphson's
        # friclaw==5 branch does not call rate_state_normal_stress), so the
        # slot is left untouched rather than written back unchanged.
        fric = B.setat(xp, fric, (slice(None), gv.THETA_PC), thetaPc_out)

    Ts = taoc_old * 0.5 * (srS / srMag_sd) + taoc_new * 0.5 * (trialS / trialMag)
    Td = taoc_old * 0.5 * (srD / srMag_sd) + taoc_new * 0.5 * (trialD / trialMag)

    fric = B.setat(xp, fric, (slice(None), gv.TRACT_NORM), tnrm)
    fric = B.setat(xp, fric, (slice(None), gv.TRACT_STRIKE), Ts)
    fric = B.setat(xp, fric, (slice(None), gv.TRACT_DIP), Td)
    fric = B.setat(xp, fric, (slice(None), gv.PEAK_SLIPRATE), v_trial)
    fric = B.setat(xp, fric, (slice(None), gv.SHEAR_MAG),
                   xp.sqrt(Ts ** 2 + Td ** 2))

    # --- relative acceleration -> nodal forces (faulting.f90:277-298) ---
    accN = -srN / dt - slipN / dt / dt
    accS = (v_trial * (trialS / trialMag) - srS) / dt
    accD = (v_trial * (trialD / trialMag) - srD) / dt
    xAcc = accN * un[:, 0] + accS * us[:, 0] + accD * ud[:, 0]
    yAcc = accN * un[:, 1] + accS * us[:, 1] + accD * ud[:, 1]
    zAcc = accN * un[:, 2] + accS * us[:, 2] + accD * ud[:, 2]
    xR = comps['fx_s'] + comps['fx_m']
    yR = comps['fy_s'] + comps['fy_m']
    zR = comps['fz_s'] + comps['fz_m']

    massSlave = finv['massSlave']; massMaster = finv['massMaster']
    idxF_s, idxF_m = finv['idxF_s'], finv['idxF_m']
    # solveRSF REPLACES the fault-node force (it does not accumulate into it
    # the way solveSWTW does) -- faulting.f90:292-293.
    # The bracketed acceleration terms are shared by the force write
    # (faulting.f90:292-293, scaled by mr) and the stored master/slave
    # velocities (:296-297, scaled by dt). They are computed ONCE, here, and
    # used by both -- NOT recovered from the force by dividing mr back out,
    # which would be a multiply and a divide instead of the Fortran's single
    # multiply and would not be bit-identical.
    acc_s = (-xAcc + xR / massMaster, -yAcc + yR / massMaster,
             -zAcc + zR / massMaster)
    acc_m = (xAcc + xR / massSlave, yAcc + yR / massSlave,
             zAcc + zR / massSlave)
    for d in range(3):
        force = B.setat(xp, force, idxF_s[d], acc_s[d] * mr)
        force = B.setat(xp, force, idxF_m[d], acc_m[d] * mr)

    vm = (comps['vx_m'], comps['vy_m'], comps['vz_m'])
    vs = (comps['vx_s'], comps['vy_s'], comps['vz_s'])
    for d in range(3):
        fric = B.setat(xp, fric, (slice(None), gv.VEL_MASTER_X + d),
                       vm[d] + acc_m[d] * dt)
        fric = B.setat(xp, fric, (slice(None), gv.VEL_SLAVE_X + d),
                       vs[d] + acc_s[d] * dt)
    force = B.setat(xp, force, 0, 0.0)
    return fric, force, srMag_sd
